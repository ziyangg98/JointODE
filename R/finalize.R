# Model Finalization (Joint + Marginal)

#' @noRd
.finalize_joint <- function(
  data_list,
  parameters,
  loglik,
  control,
  coef_names,
  converged,
  posteriors,
  posterior_moments
) {
  n_subjects <- length(data_list)
  n_re <- length(posterior_moments[[1]]$mean)
  random_effects <- t(vapply(posterior_moments, `[[`, numeric(n_re), "mean"))

  names(parameters$coefficients$baseline) <- coef_names$baseline
  names(parameters$coefficients$hazard) <- coef_names$hazard
  names(parameters$coefficients$longitudinal) <- coef_names$longitudinal
  names(parameters$coefficients$initial_state) <- coef_names$initial_state

  coef_names_expanded <- c(
    paste0("baseline:", coef_names$baseline),
    paste0("hazard:", coef_names$hazard),
    paste0("longitudinal:", coef_names$longitudinal),
    paste0("initial state:", coef_names$initial_state)
  )

  vcov_matrix <- .compute_vcov_sem(
    data_list,
    posteriors,
    posterior_moments,
    parameters,
    random_effects,
    control
  )
  dimnames(vcov_matrix) <- list(coef_names_expanded, coef_names_expanded)

  n_params <- .count_params(parameters)
  aic <- -2 * loglik + 2 * n_params
  bic <- -2 * loglik + n_params * log(n_subjects)

  ode_solutions <- .solve_batch_joint(data_list, random_effects, parameters)
  # nolint start: object_usage_linter
  risk_scores <- vapply(
    ode_solutions,
    function(x) tail(x$log_hazard, 1),
    numeric(1)
  )
  event_times <- vapply(data_list, `[[`, numeric(1), "time")
  event_status <- vapply(data_list, `[[`, numeric(1), "status")

  cindex <- survival::concordance(
    Surv(event_times, event_status) ~ risk_scores,
    reverse = TRUE
  )$concordance
  # nolint end

  if (control$verbose > 0) {
    cli::cli_alert_info(sprintf("C-index (concordance): %.3f", cindex))
  }

  list(
    random_effects = random_effects,
    vcov = vcov_matrix,
    loglik = loglik,
    aic = aic,
    bic = bic,
    cindex = cindex
  )
}

#' @noRd
.finalize_marginal <- function(
  theta, sse, data_list, biomarker_clamp,
  param_names, converged, n_iter, has_state,
  control, cl
) {
  hess <- attr(
    .compute_marginal_objective(
      theta, data_list, biomarker_clamp, FALSE, TRUE
    ),
    "hessian"
  )
  vcov_mat <- solve(hess)
  dimnames(vcov_mat) <- list(param_names, param_names)

  n_obs <- .n_obs(data_list)
  n_params <- length(theta)
  sigma_e <- sqrt(sse / (n_obs - n_params))
  ll <- -0.5 * n_obs * (log(2 * pi) + log(sse / n_obs) + 1)

  result <- structure(list(
    parameters = setNames(theta, param_names),
    measurement_error_sd = sigma_e,
    logLik = ll,
    AIC = 2 * n_params - 2 * ll,
    BIC = log(n_obs) * n_params - 2 * ll,
    convergence = list(
      converged = converged, iterations = n_iter,
      message = sprintf(
        "%s after %d iterations",
        if (converged) "Converged" else "Did not converge",
        n_iter
      )
    ),
    vcov = vcov_mat, data = data_list,
    control = control, call = cl
  ), class = "MarginalODE")

  if (!has_state) {
    mat <- do.call(
      rbind, lapply(data_list, `[[`, "initial_state")
    )
    dimnames(mat) <- list(names(data_list), c("m0", "v0"))
    result$initial_states <- mat
  }

  result
}
