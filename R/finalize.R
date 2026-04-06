# TMB Result Extraction

#' Extract results from fitted MarginalODE TMB object
#' @noRd
.finalize_marginal <- function(obj, opt, coef_names, n_re, n_subjects) {
  sdr <- TMB::sdreport(obj)
  reported <- obj$report()
  par <- obj$env$last.par.best
  pn <- names(par)

  longitudinal <- as.numeric(par[pn == "longitudinal"])
  names(longitudinal) <- coef_names$longitudinal
  initial_state <- as.numeric(par[pn == "initial_state"])
  names(initial_state) <- coef_names$initial_state
  sigma_e <- unname(exp(par[pn == "log_sigma_e"]))

  parameters <- c(longitudinal, initial_state)

  # Vcov of fixed effects
  n_fixed <- length(parameters)
  vcov_matrix <- if (!is.null(sdr$cov.fixed) && nrow(sdr$cov.fixed) >= n_fixed) {
    sdr$cov.fixed[seq_len(n_fixed), seq_len(n_fixed), drop = FALSE]
  } else {
    matrix(NA, n_fixed, n_fixed)
  }
  dimnames(vcov_matrix) <- list(names(parameters), names(parameters))

  # Variance component SEs
  sdr_report <- summary(sdr, "report")
  sdr_names <- rownames(sdr_report)
  sigma_e_se <- as.numeric(sdr_report[sdr_names == "sigma_e", "Std. Error"])
  sigma_b <- as.matrix(reported$Sigma_b)
  sigma_b_se <- matrix(
    sdr_report[sdr_names == "Sigma_b", "Std. Error"], n_re, n_re)

  # Random effects posterior modes
  random_effects <- matrix(par[pn == "random_effects"],
                           nrow = n_subjects, ncol = n_re)

  loglik <- -opt$objective
  n_total_params <- n_fixed + 1  # +1 for sigma_e
  converged <- opt$convergence == 0

  list(
    parameters = parameters,
    measurement_error_sd = sigma_e,
    measurement_error_sd_se = sigma_e_se,
    random_effect_sigma = sigma_b,
    random_effect_sigma_se = sigma_b_se,
    logLik = loglik,
    AIC = -2 * loglik + 2 * n_total_params,
    BIC = -2 * loglik + n_total_params * log(n_subjects),
    vcov = vcov_matrix,
    random_effects = random_effects,
    convergence = list(
      converged = converged,
      iterations = opt$iterations,
      message = sprintf("%s (%s)",
        if (converged) "Converged" else "Did not converge", opt$message))
  )
}

#' @noRd
.finalize_joint <- function(obj, opt, parameters, coef_names,
                          data_list, n_re, control) {
  sdr <- TMB::sdreport(obj)
  reported <- obj$report()
  par <- obj$env$last.par.best
  pn <- names(par)

  # Fixed effects
  cf <- parameters$coefficients
  for (nm in c("baseline", "hazard", "longitudinal", "initial_state"))
    cf[[nm]] <- setNames(as.numeric(par[pn == nm]), coef_names[[nm]])
  cf$measurement_error_sd <- unname(exp(par[pn == "log_sigma_e"]))
  cf$random_effect_sigma <- as.matrix(reported$Sigma_b)

  # Variance component SEs (delta method via ADREPORT)
  sdr_report <- summary(sdr, "report")
  sdr_names <- rownames(sdr_report)
  cf$measurement_error_sd_se <- as.numeric(
    sdr_report[sdr_names == "sigma_e", "Std. Error"])
  cf$random_effect_sigma_se <- matrix(
    sdr_report[sdr_names == "Sigma_b", "Std. Error"], n_re, n_re)

  parameters$coefficients <- cf
  parameters$random_effects_init <- NULL

  # Vcov of fixed effects
  coef_names_exp <- .prefixed_coef_names(coef_names)
  n_fixed <- length(coef_names_exp)
  vcov_matrix <- if (!is.null(sdr$cov.fixed) && nrow(sdr$cov.fixed) >= n_fixed) {
    sdr$cov.fixed[seq_len(n_fixed), seq_len(n_fixed), drop = FALSE]
  } else {
    matrix(NA, n_fixed, n_fixed)
  }
  dimnames(vcov_matrix) <- list(coef_names_exp, coef_names_exp)

  # C-index
  n_subjects <- length(data_list)
  event_t <- vapply(data_list, `[[`, numeric(1), "time")
  event_s <- vapply(data_list, `[[`, numeric(1), "status")
  cindex <- survival::concordance(
    Surv(event_t, event_s) ~ as.numeric(reported$log_hazard_at_event),
    reverse = TRUE
  )$concordance

  # Convergence reporting
  loglik <- -opt$objective
  converged <- opt$convergence == 0
  if (control$verbose > 0) {
    if (converged) cli::cli_alert_success(sprintf("Converged (%s)", opt$message))
    else cli::cli_alert_warning(sprintf("Did not converge: %s", opt$message))
    cli::cli_alert_info(sprintf("Log-likelihood: %.2f", loglik))
    cli::cli_alert_info(sprintf("C-index: %.3f", cindex))
  }

  n_params <- .count_params(parameters)
  list(
    parameters = parameters,
    logLik = loglik,
    AIC = -2 * loglik + 2 * n_params,
    BIC = -2 * loglik + n_params * log(n_subjects),
    cindex = cindex,
    convergence = list(
      converged = converged,
      iterations = opt$iterations,
      message = sprintf("%s (%s)",
        if (converged) "Converged" else "Did not converge", opt$message)),
    random_effects = matrix(par[pn == "random_effects"], nrow = n_subjects, ncol = n_re),
    vcov = vcov_matrix,
    tmb_report = reported
  )
}
