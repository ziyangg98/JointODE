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
  parallel = FALSE, n_cores = 0
) {
  eval_chunk <- function(idx) {
    .compute_joint_objective(
      params = params, data_list = data_list[idx],
      random_effects = random_effects[idx, , drop = FALSE],
      parameters = parameters,
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
      parameters = parameters, gradient = TRUE, hessian = FALSE
    )
    value <- -as.numeric(result)
    attr(value, "gradient") <- -as.vector(attr(result, "gradient"))
    value
  }

  fit <- nlm(
    f = objective, p = random_effect,
    hessian = FALSE, check.analyticals = FALSE
  )
  if (fit$code >= 4) {
    warning(sprintf("Laplace nlm did not converge (code=%d)", fit$code),
            call. = FALSE)
  }

  # Hessian only at mode (not during nlm iterations)
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
  theta <- .coef_to_vector(parameters)

  if (control$verbose > 0) cli::cli_alert_info("Computing SEM vcov...")

  obj <- .compute_objective_expected(
    theta, data_list, random_effects, parameters,
    gradient = FALSE, hessian = TRUE,
    parallel = control$parallel, n_cores = control$n_cores
  )
  hess <- attr(obj, "hessian")
  n_coef <- nrow(hess)
  info_complete_inv <- solve(hess)

  em_map <- function(theta_input) {
    params_input <- .vector_to_coef(parameters, theta_input)
    result <- .em_step(data_list, params_input, random_effects, control)
    .coef_to_vector(result$parameters)
  }

  dm_matrix <- numDeriv::jacobian(func = em_map, x = theta, method = "simple")

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

# Laplace EM Step ==============================================================

#' @noRd
.em_step <- function(data_list, parameters, random_effects, control) {
  n <- length(data_list)
  n_re <- ncol(random_effects)

  # E-step: Laplace approximation — posterior mode + covariance
  posteriors <- .compute_posteriors(
    data_list, parameters, random_effects,
    parallel = control$parallel, n_cores = control$n_cores, setup = FALSE
  )
  random_effects <- t(vapply(posteriors, `[[`, numeric(n_re), "mode"))

  # M-step: one-step Newton on complete-data log-likelihood at posterior modes
  eval_objective <- function(theta_val, gradient = FALSE, hessian = FALSE) {
    .compute_objective_expected(
      theta_val, data_list, random_effects,
      .vector_to_coef(parameters, theta_val),
      gradient = gradient, hessian = hessian,
      parallel = control$parallel, n_cores = control$n_cores
    )
  }

  theta <- .coef_to_vector(parameters)
  obj <- eval_objective(theta, gradient = TRUE, hessian = TRUE)
  g <- as.vector(attr(obj, "gradient"))
  H <- as.matrix(attr(obj, "hessian"))
  R_h <- .safe_chol(H)
  direction <- -backsolve(R_h, forwardsolve(t(R_h), g))

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

  # Laplace-corrected variance updates
  # Sigma_b = (1/n) sum_i (b_i b_i^T + H_i^{-1})
  parameters$coefficients$random_effect_sigma <-
    .update_random_effect_sigma(random_effects, posteriors)

  # sigma_e at posterior modes
  parameters$coefficients$measurement_error_sd <-
    .update_measurement_error_sd(data_list, parameters, random_effects)

  list(
    parameters = parameters,
    random_effects = random_effects,
    loglik = loglik
  )
}
