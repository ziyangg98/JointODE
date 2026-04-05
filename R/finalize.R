# Model Finalization (Marginal only — Joint finalization is inline in JointODE.R)

#' @noRd
.finalize_marginal <- function(
  theta, sse, data_list, biomarker_clamp,
  param_names, converged, n_iter,
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

  mat <- do.call(
    rbind, lapply(data_list, `[[`, "initial_state")
  )
  dimnames(mat) <- list(names(data_list), c("biomarker", "velocity"))
  result$initial_states <- mat

  result
}
