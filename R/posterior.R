# Posterior Computation and EM Algorithm

#' @importFrom stats nlm
NULL

# Quadrature Cache =============================================================

.quad_cache <- new.env(parent = emptyenv())

.get_quad_grid <- function(dim, level) {
  if (dim >= 4 && level > 3) level <- 3
  key <- paste0(dim, "_", level)
  if (!exists(key, envir = .quad_cache)) {
    grid <- mvQuad::createNIGrid(
      dim = dim, type = "GHe",
      level = level, ndConstruction = "product"
    )
    assign(key, list(
      nodes = mvQuad::getNodes(grid),
      weights = as.numeric(mvQuad::getWeights(grid))
    ), envir = .quad_cache)
  }
  get(key, envir = .quad_cache)
}

# Posterior Expansion ==========================================================

.expand_posteriors <- function(data_list, posteriors) {
  n_subjects <- length(data_list)
  n_nodes_per_subject <- vapply(posteriors$nodes, nrow, integer(1))
  total_nodes <- sum(n_nodes_per_subject)

  all_nodes <- do.call(rbind, posteriors$nodes)
  all_weights <- unlist(posteriors$weights, use.names = FALSE)
  all_data <- vector("list", total_nodes)

  idx <- 1L
  for (i in seq_len(n_subjects)) {
    n <- n_nodes_per_subject[i]
    all_data[idx:(idx + n - 1L)] <- rep(list(data_list[[i]]), n)
    idx <- idx + n
  }

  list(
    nodes = all_nodes, data = all_data, weights = all_weights,
    node_to_subject = rep(seq_len(n_subjects), n_nodes_per_subject)
  )
}

# Variance Updates =============================================================

.update_measurement_error_sd <- function(data_list, parameters, posteriors) {
  n_total_obs <- sum(vapply(data_list, function(d) {
    length(d$longitudinal$measurements)
  }, integer(1)))

  expanded <- .expand_posteriors(data_list, posteriors)

  ode_solutions <- .solve_batch_joint(
    data_list = expanded$data,
    random_effects = expanded$nodes,
    parameters = parameters
  )

  sigma_e_squared <- sum(vapply(seq_along(ode_solutions), function(k) {
    i <- expanded$node_to_subject[k]
    measurements <- data_list[[i]]$longitudinal$measurements
    if (length(measurements) == 0) {
      return(0)
    }

    obs_times <- data_list[[i]]$longitudinal$times
    result_times <- ode_solutions[[k]]$times
    match_idx <- vapply(obs_times, function(t) {
      which.min(abs(result_times - t))
    }, integer(1))
    biomarker_pred <- ode_solutions[[k]]$biomarker[match_idx]

    expanded$weights[k] * sum((measurements - biomarker_pred)^2)
  }, numeric(1)))

  sqrt(sigma_e_squared / n_total_obs)
}

#' @noRd
.update_random_effect_sigma <- function(posterior_moments, n_subjects) {
  n_re <- length(posterior_moments[[1]]$mean)
  sigma_sum <- matrix(0, n_re, n_re)
  for (i in seq_len(n_subjects)) {
    sigma_sum <- sigma_sum + posterior_moments[[i]]$second_moment
  }
  sigma_sum / n_subjects
}

# Expected Objective ===========================================================

.compute_objective_expected <- function(
  params, data_list, posteriors, parameters,
  gradient = TRUE, hessian = FALSE,
  parallel = FALSE, n_cores = 0
) {
  n_subjects <- length(data_list)

  if (parallel && n_subjects > 1) {
    if (n_cores == 0) {
      n_cores <- max(1L, parallel::detectCores() - 1L)
    }
    n_chunks <- min(n_cores, n_subjects)
    chunks <- split(seq_len(n_subjects), cut(seq_len(n_subjects),
      n_chunks,
      labels = FALSE
    ))

    compute_chunk <- function(subject_indices) {
      chunk_data <- data_list[subject_indices]
      chunk_posteriors <- list(
        nodes = posteriors$nodes[subject_indices],
        weights = posteriors$weights[subject_indices]
      )
      expanded <- .expand_posteriors(chunk_data, chunk_posteriors)
      .compute_joint_objective(
        params = params,
        data_list = expanded$data,
        random_effects = expanded$nodes,
        parameters = parameters,
        weights = expanded$weights,
        gradient = gradient,
        hessian = hessian
      )
    }

    results <- .parallel_apply(
      chunks, compute_chunk,
      parallel = TRUE, n_cores = n_cores, setup = FALSE
    )

    total_obj <- sum(vapply(results, as.numeric, numeric(1)))
    result <- structure(total_obj, names = NULL)

    if (gradient) {
      total_grad <- Reduce(`+`, lapply(results, function(r) {
        as.vector(attr(r, "gradient"))
      }))
      attr(result, "gradient") <- total_grad
    }
    if (hessian) {
      total_hess <- Reduce(`+`, lapply(results, function(r) {
        attr(r, "hessian")
      }))
      attr(result, "hessian") <- total_hess
    }
    result
  } else {
    expanded <- .expand_posteriors(data_list, posteriors)
    .compute_joint_objective(
      params = params,
      data_list = expanded$data,
      random_effects = expanded$nodes,
      parameters = parameters,
      weights = expanded$weights,
      gradient = gradient,
      hessian = hessian
    )
  }
}

# Laplace Approximation + AGHQ ================================================

.compute_posterior_laplace <- function(
  data, random_effect, parameters, level = 3
) {
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
  hessian_neglogpost <- -attr(result_at_mode, "hessian")

  R <- .regularized_chol(hessian_neglogpost)
  chol_factor <- t(backsolve(R, diag(nrow(hessian_neglogpost))))

  quad <- .get_quad_grid(length(fit$estimate), level)
  n_nodes <- nrow(quad$nodes)
  n_re <- length(fit$estimate)

  aghq_nodes <- t(vapply(seq_len(n_nodes), function(k) {
    fit$estimate + tcrossprod(quad$nodes[k, ], chol_factor)
  }, numeric(n_re)))

  logpost_at_nodes <- vapply(seq_len(n_nodes), function(k) {
    as.numeric(.compute_joint_logpost(
      random_effect = aghq_nodes[k, ],
      data = data, parameters = parameters,
      gradient = FALSE, hessian = FALSE
    ))
  }, numeric(1))

  log_jacobian <- -sum(log(diag(R)))

  log_weights <- log(quad$weights) + log_jacobian + logpost_at_nodes
  max_log_weight <- max(log_weights)
  aghq_weights <- exp(log_weights - max_log_weight)
  aghq_weights <- aghq_weights / sum(aghq_weights)

  list(nodes = aghq_nodes, weights = aghq_weights)
}

.compute_posterior_moments <- function(posterior_result) {
  nodes <- posterior_result$nodes
  weights <- posterior_result$weights
  mean <- colSums(weights * nodes)
  weighted_nodes <- sqrt(weights) * nodes
  second_moment <- crossprod(weighted_nodes)
  list(mean = mean, second_moment = second_moment)
}

.compute_posteriors <- function(
  data_list, parameters, random_effects,
  parallel = FALSE, n_cores = 0, level = 3, setup = TRUE
) {
  n_subjects <- length(data_list)
  .get_quad_grid(ncol(random_effects), level)

  posterior_results <- .parallel_apply(
    seq_len(n_subjects),
    function(i) {
      .compute_posterior_laplace(
        data = data_list[[i]],
        random_effect = random_effects[i, ],
        parameters = parameters,
        level = level
      )
    },
    parallel = parallel, n_cores = n_cores, setup = setup
  )

  list(
    nodes = lapply(posterior_results, `[[`, "nodes"),
    weights = lapply(posterior_results, `[[`, "weights")
  )
}

# SEM Variance-Covariance =====================================================

#' @noRd
.m_step_map <- function(
  theta_input, data_list, posteriors, posterior_moments, parameters, control
) {
  params_input <- .vector_to_coef(parameters, theta_input)
  n_subjects <- length(data_list)

  random_effect_sigma <- .update_random_effect_sigma(
    posterior_moments, n_subjects
  )
  measurement_error_sd <- .update_measurement_error_sd(
    data_list, params_input, posteriors
  )
  params_input$coefficients$measurement_error_sd <- measurement_error_sd
  params_input$coefficients$random_effect_sigma <- random_effect_sigma

  theta <- .coef_to_vector(params_input)
  obj <- .compute_objective_expected(
    theta, data_list, posteriors, params_input,
    gradient = TRUE, hessian = TRUE,
    parallel = control$parallel, n_cores = control$n_cores
  )
  as.vector(theta - .regularized_solve(
    attr(obj, "hessian"), attr(obj, "gradient")
  ))
}

#' @noRd
.compute_vcov_sem <- function(
  data_list, posteriors, posterior_moments,
  parameters, random_effects, control
) {
  theta <- .coef_to_vector(parameters)
  n_coef <- length(theta)

  if (control$verbose > 0) cli::cli_alert_info("Computing SEM vcov...")

  obj <- .compute_objective_expected(
    theta, data_list, posteriors, parameters,
    gradient = FALSE, hessian = TRUE,
    parallel = control$parallel, n_cores = control$n_cores
  )
  hess <- attr(obj, "hessian")
  n_coef <- nrow(hess)
  r_hess <- try(.regularized_chol(hess), silent = TRUE)
  if (inherits(r_hess, "try-error")) {
    if (control$verbose > 0) {
      cli::cli_alert_warning("Hessian singular, vcov unavailable")
    }
    return(matrix(NA, n_coef, n_coef))
  }
  info_complete_inv <- chol2inv(r_hess)

  em_map <- function(theta_input) {
    .m_step_map(
      theta_input, data_list, posteriors, posterior_moments,
      parameters, control
    )
  }

  dm_matrix <- numDeriv::jacobian(
    func = em_map, x = theta, method = "simple"
  )

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
  if (any(diag_ratio < 1 - 1e-6) && control$verbose > 0) {
    cli::cli_alert_warning(sprintf(
      "%d parameter(s) violate missing information principle (min ratio: %.3f)",
      sum(diag_ratio < 1 - 1e-6), min(diag_ratio)
    ))
  }

  vcov_sym
}

# EM Step ======================================================================

#' @noRd
.em_step <- function(data_list, parameters, random_effects, control) {
  posteriors <- .compute_posteriors(
    data_list, parameters, random_effects,
    control$parallel, control$n_cores,
    control$quad_level,
    setup = FALSE
  )

  n_subjects <- length(data_list)
  posterior_moments <- lapply(seq_len(n_subjects), function(i) {
    .compute_posterior_moments(
      list(nodes = posteriors$nodes[[i]], weights = posteriors$weights[[i]])
    )
  })
  n_re <- length(posterior_moments[[1]]$mean)
  random_effects <- t(vapply(posterior_moments, `[[`, numeric(n_re), "mean"))

  random_effect_sigma <- .update_random_effect_sigma(
    posterior_moments, n_subjects
  )
  measurement_error_sd <- .update_measurement_error_sd(
    data_list, parameters, posteriors
  )
  parameters$coefficients$measurement_error_sd <- measurement_error_sd
  parameters$coefficients$random_effect_sigma <- random_effect_sigma

  # M-step: one Newton step for fixed effects
  theta <- .coef_to_vector(parameters)
  obj <- .compute_objective_expected(
    theta, data_list, posteriors, parameters,
    gradient = TRUE, hessian = TRUE,
    parallel = control$parallel, n_cores = control$n_cores
  )
  theta_new <- as.vector(theta - .regularized_solve(
    attr(obj, "hessian"), attr(obj, "gradient")
  ))
  loglik_value <- -as.numeric(obj)
  parameters <- .vector_to_coef(parameters, theta_new)

  list(
    parameters = parameters,
    random_effects = random_effects,
    loglik = loglik_value,
    posteriors = posteriors,
    posterior_moments = posterior_moments
  )
}
