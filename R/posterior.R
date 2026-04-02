# Posterior Computation and EM Algorithm

#' @importFrom stats nlm
NULL

# Variance Updates =============================================================

.update_measurement_error_sd <- function(data_list, parameters,
                                         random_effects) {
  n_total_obs <- sum(vapply(data_list, function(d) {
    length(d$longitudinal$measurements)
  }, integer(1)))
  ode_solutions <- .solve_batch_joint(
    data_list = data_list, random_effects = random_effects,
    parameters = parameters
  )
  rss <- sum(vapply(seq_along(ode_solutions), function(k) {
    measurements <- data_list[[k]]$longitudinal$measurements
    if (length(measurements) == 0) return(0)
    obs_times <- data_list[[k]]$longitudinal$times
    match_idx <- match(obs_times, ode_solutions[[k]]$times)
    sum((measurements - ode_solutions[[k]]$biomarker[match_idx])^2)
  }, numeric(1)))
  sqrt(rss / n_total_obs)
}

#' @noRd
.update_random_effect_sigma <- function(random_effects, posteriors) {
  n <- nrow(random_effects)
  sigma <- crossprod(random_effects)
  for (i in seq_len(n)) sigma <- sigma + posteriors[[i]]$cov
  sigma / n
}

# Expected Objective ===========================================================

.compute_objective_expected <- function(
  params, data_list, random_effects, parameters,
  gradient = TRUE, hessian = FALSE,
  weights = NULL,
  parallel = FALSE, n_cores = 0
) {
  eval_chunk <- function(idx) {
    w <- if (!is.null(weights)) weights[idx] else NULL
    .compute_joint_objective(
      params = params, data_list = data_list[idx],
      random_effects = random_effects[idx, , drop = FALSE],
      parameters = parameters,
      weights = w,
      gradient = gradient, hessian = hessian
    )
  }

  n <- length(data_list)
  if (!parallel || n <= 1) return(eval_chunk(seq_len(n)))

  n_cores <- .resolve_cores(n_cores)
  chunks <- split(seq_len(n), cut(seq_len(n), min(n_cores, n),
    labels = FALSE
  ))

  results <- .parallel_apply(chunks, eval_chunk,
    parallel = TRUE, n_cores = n_cores, setup = FALSE
  )

  obj <- sum(vapply(results, as.numeric, numeric(1)))
  result <- structure(obj, names = NULL)
  if (gradient) {
    attr(result, "gradient") <- Reduce(`+`, lapply(
      results, function(r) as.vector(attr(r, "gradient"))
    ))
  }
  if (hessian) {
    attr(result, "hessian") <- Reduce(`+`, lapply(
      results, function(r) attr(r, "hessian")
    ))
  }
  result
}

# Laplace Approximation ========================================================

.compute_posterior_laplace <- function(data, random_effect, parameters) {
  objective <- function(theta) {
    result <- .compute_joint_logpost(
      random_effect = theta, data = data,
      parameters = parameters, gradient = TRUE, hessian = TRUE
    )
    value <- -as.numeric(result)
    attr(value, "gradient") <- -as.vector(attr(result, "gradient"))
    attr(value, "hessian") <- -as.matrix(attr(result, "hessian"))
    value
  }

  fit <- nlm(
    f = objective, p = random_effect,
    hessian = FALSE, check.analyticals = FALSE
  )

  result_at_mode <- .compute_joint_logpost(
    random_effect = fit$estimate, data = data,
    parameters = parameters, gradient = FALSE, hessian = TRUE
  )
  R <- .safe_chol(-attr(result_at_mode, "hessian"))

  list(mode = fit$estimate, cov = chol2inv(R))
}

.compute_posteriors <- function(
  data_list, parameters, random_effects,
  parallel = FALSE, n_cores = 0, setup = TRUE
) {
  .parallel_apply(
    seq_along(data_list),
    function(i) {
      .compute_posterior_laplace(
        data = data_list[[i]],
        random_effect = random_effects[i, ],
        parameters = parameters
      )
    },
    parallel = parallel, n_cores = n_cores, setup = setup
  )
}

# SEM Variance-Covariance =====================================================

#' @noRd
.compute_vcov_sem <- function(
  data_list, parameters, random_effects, control
) {
  sem_control <- control
  sem_control$mc_samples <- max(1L, as.integer(sem_control$mc_samples))

  theta <- .coef_to_vector(parameters)

  if (control$verbose > 0) cli::cli_alert_info("Computing SEM vcov...")

  obj <- .compute_objective_expected(
    theta, data_list, random_effects, parameters,
    gradient = FALSE, hessian = TRUE,
    parallel = sem_control$parallel, n_cores = sem_control$n_cores
  )
  hess <- attr(obj, "hessian")
  n_coef <- nrow(hess)
  info_complete_inv <- solve(hess)

  em_map <- function(theta_input) {
    # Fixed seed makes the MCEM mapping deterministic for numerical Jacobian.
    set.seed(20260402)
    params_input <- .vector_to_coef(parameters, theta_input)
    result <- .em_step(data_list, params_input, random_effects, sem_control)
    .coef_to_vector(result$parameters)
  }

  dm_matrix <- numDeriv::jacobian(func = em_map, x = theta,
    method.args = list(r = 2))

  info_observed_inv <- info_complete_inv %*% solve(diag(n_coef) - t(dm_matrix))

  asymmetry <- max(abs(info_observed_inv - t(info_observed_inv)))
  if (control$verbose > 0 && asymmetry > 1e-4) {
    cli::cli_alert_warning(sprintf(
      "Vcov matrix asymmetry: %.2e", asymmetry
    ))
  }
  vcov_sym <- (info_observed_inv + t(info_observed_inv)) / 2

  eigenvalues <- eigen(vcov_sym, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigenvalues <= 0) && control$verbose > 0) {
    cli::cli_alert_warning(sprintf(
      "%d non-positive eigenvalue(s) (min: %.2e)",
      sum(eigenvalues <= 0), min(eigenvalues)
    ))
  }

  diag_ratio <- diag(vcov_sym) / diag(info_complete_inv)
  if (any(diag_ratio < 1 - 0.01) && control$verbose > 0) {
    cli::cli_alert_warning(sprintf(
      "%d parameter(s) violate missing information principle (min ratio: %.3f)",
      sum(diag_ratio < 1 - 0.01), min(diag_ratio)
    ))
  }

  vcov_sym
}

# MCEM Step ====================================================================

#' @noRd
.em_step <- function(data_list, parameters, random_effects, control) {
  .em_step_mcem(data_list, parameters, random_effects, control)
}

# Monte Carlo EM ==============================================================

#' MCEM E-step: Independent Metropolis-Hastings sampling from posterior
#' @references Booth & Hobert (1999), Wei & Tanner (1990)
#' @noRd
.compute_posteriors_mcem <- function(data_list, parameters, random_effects,
                                     M, mc_burnin, parallel, n_cores) {
  .parallel_apply(seq_along(data_list), function(i) {
    # Laplace proposal distribution: N(mode, cov)
    post <- .compute_posterior_laplace(
      data_list[[i]], random_effects[i, ], parameters
    )
    n_re <- length(post$mode)
    R <- .safe_chol(post$cov)
    n_total <- mc_burnin + M

    # Pre-generate all proposals from N(mode, cov)
    z <- matrix(rnorm(n_total * n_re), n_total, n_re)
    proposals <- sweep(z %*% R, 2, post$mode, `+`)

    # Batch-evaluate log-posteriors for all proposals
    all_logp <- as.numeric(.compute_joint_logpost_batch(
      samples = proposals, data = data_list[[i]], parameters = parameters
    ))

    # Log-proposal densities: -0.5 * (b - mode)' Sigma^{-1} (b - mode)
    centered <- sweep(proposals, 2, post$mode, `-`)
    V <- forwardsolve(t(R), t(centered))
    all_logq <- -0.5 * colSums(V^2)

    # Initialize chain at posterior mode
    mode_logp <- as.numeric(.compute_joint_logpost_batch(
      samples = matrix(post$mode, nrow = 1),
      data = data_list[[i]], parameters = parameters
    ))
    current <- post$mode
    current_logp <- mode_logp
    current_logq <- 0.0  # (mode - mode)' Sigma^{-1} (mode - mode) = 0

    # MH accept/reject loop
    samples <- matrix(NA_real_, M, n_re)
    u <- log(runif(n_total))
    n_accept <- 0L

    for (t in seq_len(n_total)) {
      log_alpha <- all_logp[t] - current_logp + current_logq - all_logq[t]
      if (u[t] < log_alpha) {
        current <- proposals[t, ]
        current_logp <- all_logp[t]
        current_logq <- all_logq[t]
        if (t > mc_burnin) n_accept <- n_accept + 1L
      }
      if (t > mc_burnin) {
        samples[t - mc_burnin, ] <- current
      }
    }

    # Equal weights for MCMC samples
    weights <- rep(1.0 / M, M)

    list(samples = samples, weights = weights,
         mode = post$mode, cov = post$cov,
         accept_rate = n_accept / M)
  }, parallel = parallel, n_cores = n_cores, setup = FALSE)
}

#' @noRd
.em_step_mcem <- function(data_list, parameters, random_effects, control) {
  M <- control$mc_samples
  if (M < 1) {
    stop("mc_samples must be a positive integer for MCEM optimization")
  }
  n <- length(data_list)
  n_re <- ncol(random_effects)

  # E-step: Independent Metropolis-Hastings from posterior
  mc_post <- .compute_posteriors_mcem(
    data_list, parameters, random_effects, M,
    control$mc_burnin, control$parallel, control$n_cores
  )

  if (control$verbose >= 2) {
    rates <- vapply(mc_post, `[[`, numeric(1), "accept_rate")
    cli::cli_text(sprintf(
      "    MH accept rate: mean=%.1f%%, min=%.1f%%, max=%.1f%%",
      mean(rates) * 100, min(rates) * 100, max(rates) * 100
    ))
  }

  # MCMC posterior mean (equal weights)
  random_effects <- t(vapply(mc_post, function(p) {
    colMeans(p$samples)
  }, numeric(n_re)))

  # M-step: one-step Newton on MC-averaged complete-data log-likelihood
  idx_map <- rep(seq_len(n), each = M)
  data_list_rep <- data_list[idx_map]
  re_matrix <- do.call(rbind, lapply(mc_post, `[[`, "samples"))
  weights_vec <- unlist(lapply(mc_post, `[[`, "weights"))

  eval_objective <- function(theta_val, gradient = FALSE, hessian = FALSE) {
    .compute_objective_expected(
      theta_val, data_list_rep, re_matrix,
      .vector_to_coef(parameters, theta_val),
      gradient = gradient, hessian = hessian,
      weights = weights_vec,
      parallel = control$parallel, n_cores = control$n_cores
    )
  }

  theta <- .coef_to_vector(parameters)
  obj <- eval_objective(theta, gradient = TRUE, hessian = TRUE)
  g <- as.vector(attr(obj, "gradient"))
  H <- as.matrix(attr(obj, "hessian"))
  R_h <- .safe_chol(H)
  direction <- -backsolve(R_h, forwardsolve(t(R_h), g))

  # One-step Newton update
  theta_new <- theta + direction
  if (any(!is.finite(theta_new))) {
    stop("One-step Newton update produced non-finite parameters")
  }
  f_new <- as.numeric(eval_objective(theta_new, gradient = FALSE, hessian = FALSE))
  if (!is.finite(f_new)) {
    stop("One-step Newton update produced non-finite objective")
  }
  parameters <- .vector_to_coef(parameters, theta_new)
  loglik <- -f_new

  # MCMC Sigma_b = (1/n) sum_i (1/M) sum_m b_{im} b_{im}'
  sigma_b <- Reduce(`+`, lapply(seq_len(n), function(i) {
    crossprod(mc_post[[i]]$samples) / M
  })) / n
  parameters$coefficients$random_effect_sigma <- sigma_b

  # sigma_e at MCMC posterior mean
  parameters$coefficients$measurement_error_sd <-
    .update_measurement_error_sd(data_list, parameters, random_effects)

  list(
    parameters = parameters,
    random_effects = random_effects,
    loglik = loglik
  )
}
