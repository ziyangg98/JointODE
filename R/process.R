# Data Processing ==============================================================

#' @importFrom stats na.pass
#' @noRd
.process_joint <- function(
  longitudinal_data,
  survival_data,
  parsed_long = NULL,
  parsed_surv = NULL,
  survival_formula = NULL,
  longitudinal_formula = NULL
) {
  # Allow backward-compatible calls with formula arguments
  if (is.null(parsed_long) && !is.null(longitudinal_formula)) {
    parsed_long <- .parse_longitudinal_formula(longitudinal_formula)
  }
  if (is.null(parsed_surv) && !is.null(survival_formula)) {
    parsed_surv <- .parse_survival_formula(survival_formula)
  }
  id <- parsed_long$grouping
  time <- parsed_surv$time_var

  # Build formulas
  fixed_formula <- .build_formula(
    parsed_long$fixed_terms,
    response = parsed_long$response
  )
  random_formula <- if (is.null(parsed_long$random_terms)) {
    NULL # No covariate random effects
  } else {
    .build_formula(parsed_long$random_terms, is_random = TRUE)
  }

  unique_ids <- unique(survival_data[[id]])
  n_subjects <- length(unique_ids)
  data_process <- vector("list", n_subjects)
  names(data_process) <- as.character(unique_ids)

  surv_frame <- model.frame(survival_formula, data = survival_data)
  surv_response_matrix <- model.response(surv_frame)
  surv_matrix <- model.matrix(survival_formula, surv_frame)
  surv_design <- if (ncol(surv_matrix) > 1) {
    surv_matrix[, -1, drop = FALSE]
  } else {
    NULL
  }

  survival_index_map <- match(unique_ids, survival_data[[id]])
  long_id_groups <- split(
    seq_len(nrow(longitudinal_data)),
    longitudinal_data[[id]]
  )

  for (i in seq_along(unique_ids)) {
    sid <- as.character(unique_ids[i])

    long_rows <- long_id_groups[[sid]]
    if (!is.null(long_rows) && length(long_rows) > 0) {
      long_subset <- longitudinal_data[long_rows, , drop = FALSE]
      long_subset <- long_subset[order(long_subset[[time]]), , drop = FALSE]

      long_frame <- model.frame(fixed_formula, data = long_subset)
      long_times <- long_subset[[time]]
      long_measurements <- model.response(long_frame)
      long_fixed_covariates <- model.matrix(fixed_formula, long_frame)
      long_random_covariates <- if (is.null(random_formula)) {
        matrix(nrow = nrow(long_subset), ncol = 0)
      } else {
        model.matrix(random_formula, long_subset)
      }
    } else {
      long_times <- numeric(0)
      long_measurements <- numeric(0)
      # Get dimensions from empty data frame
      empty_data <- longitudinal_data[0, , drop = FALSE]
      n_fixed <- ncol(model.matrix(fixed_formula, empty_data))
      long_fixed_covariates <- matrix(nrow = 0, ncol = n_fixed)
      long_random_covariates <- if (is.null(random_formula)) {
        matrix(nrow = 0, ncol = 0)
      } else {
        n_random <- ncol(model.matrix(random_formula, empty_data))
        matrix(nrow = 0, ncol = n_random)
      }
    }

    survival_row <- survival_index_map[i]
    event_time <- surv_response_matrix[survival_row, 1]
    event_status <- surv_response_matrix[survival_row, 2]
    covariates <- if (!is.null(surv_design)) {
      surv_design[survival_row, , drop = FALSE]
    } else {
      data.frame()
    }
    data_process[[i]] <- list(
      id = unique_ids[i],
      time = event_time,
      status = event_status,
      covariates = covariates,
      longitudinal = list(
        times = long_times,
        measurements = long_measurements,
        covariates = list(
          fixed = long_fixed_covariates,
          random = long_random_covariates
        )
      )
    )
  }
  data_process
}

#' Pack per-subject data into flat TMB input
#' @noRd
.pack_data <- function(data_list, parameters, control) {
  n_subjects <- length(data_list)
  configs <- parameters$configurations
  coefs <- parameters$coefficients

  # Concatenate per-subject longitudinal fields
  extract_long <- function(field) {
    lapply(data_list, function(d) d$longitudinal[[field]])
  }
  obs_times <- unlist(extract_long("times")) %||% numeric(0)
  obs_values <- unlist(extract_long("measurements")) %||% numeric(0)
  n_observations <- vapply(data_list, function(d) {
    length(d$longitudinal$measurements)
  }, integer(1))

  # Covariate counts (constant across subjects)
  n_fixed_covariates <- ncol(data_list[[1]]$longitudinal$covariates$fixed)
  n_random_covariates <- ncol(data_list[[1]]$longitudinal$covariates$random)

  # Flatten design matrices (row-major)
  flatten_design <- function(type) {
    result <- unlist(lapply(data_list, function(d)
      as.vector(t(d$longitudinal$covariates[[type]]))))
    result %||% numeric(0)
  }

  # Survival covariates
  first_covariates <- data_list[[1]]$covariates
  n_surv_covariates <- if (is.data.frame(first_covariates) ||
                           is.matrix(first_covariates)) {
    ncol(first_covariates)
  } else {
    0L
  }
  surv_covariates <- if (n_surv_covariates > 0) {
    unlist(lapply(data_list, function(d) as.numeric(d$covariates)))
  } else {
    numeric(0)
  }

  # ODE branch (same for all subjects at initialization)
  b1_init <- if (configs$biomarker$fixed) coefs$longitudinal[1] else 0
  b2_init <- if (configs$velocity$fixed) {
    coefs$longitudinal[1 + configs$biomarker$fixed]
  } else {
    0
  }
  ode_branch <- .classify_disc(b1_init, b2_init)

  baseline_config <- configs$baseline
  list(
    model_type = 0L,
    n_subjects = as.integer(n_subjects),
    n_random_effects = as.integer(ncol(parameters$random_effects_init)),
    n_observations = as.integer(n_observations),
    event_times = vapply(data_list, `[[`, numeric(1), "time"),
    event_status = as.integer(vapply(data_list, `[[`, numeric(1), "status")),
    obs_times = obs_times,
    obs_values = obs_values,
    n_fixed_covariates = as.integer(n_fixed_covariates),
    n_random_covariates = as.integer(n_random_covariates),
    long_fixed_covariates = flatten_design("fixed"),
    long_random_covariates = flatten_design("random"),
    surv_covariates = surv_covariates,
    n_surv_covariates = as.integer(n_surv_covariates),
    ode_branch = as.integer(ode_branch),
    clamp_value = configs$biomarker_clamp,
    hazard_quadrature = as.integer(control$hazard_quadrature),
    velocity_power = configs$gamma,
    baseline_degree = as.integer(baseline_config$degree),
    baseline_knots = if (length(baseline_config$knots) > 0) {
      baseline_config$knots
    } else {
      numeric(0)
    },
    baseline_boundary = baseline_config$boundary_knots,
    biomarker_fixed = as.integer(configs$biomarker$fixed),
    biomarker_random = as.integer(configs$biomarker$random),
    velocity_fixed = as.integer(configs$velocity$fixed),
    velocity_random = as.integer(configs$velocity$random)
  )
}

#' Convert R parameters to TMB parameterization
#' @noRd
.pack_params <- function(parameters) {
  coefs <- parameters$coefficients
  n_random_effects <- ncol(parameters$random_effects_init)
  sigma_matrix <- coefs$random_effect_sigma

  # Decompose Sigma_b -> log(SD) + correlation theta
  marginal_sds <- sqrt(diag(sigma_matrix))
  scale_inv <- diag(1 / pmax(marginal_sds, 1e-10))
  corr_matrix <- scale_inv %*% sigma_matrix %*% scale_inv

  # theta = rho / sqrt(1 - rho^2), lower-triangular off-diagonal
  n_corr <- n_random_effects * (n_random_effects - 1) / 2
  corr_theta <- numeric(n_corr)
  idx <- 1L
  for (col in seq_len(n_random_effects - 1)) {
    for (row in (col + 1):n_random_effects) {
      rho <- max(-0.99, min(0.99, corr_matrix[row, col]))
      corr_theta[idx] <- rho / sqrt(1 - rho^2)
      idx <- idx + 1L
    }
  }

  list(
    baseline = coefs$baseline,
    hazard = coefs$hazard,
    longitudinal = coefs$longitudinal,
    initial_state = coefs$initial_state,
    log_sigma_e = log(coefs$measurement_error_sd),
    log_sd_re = log(pmax(marginal_sds, 1e-10)),
    corr_par = corr_theta,
    random_effects = parameters$random_effects_init
  )
}

#' Classify ODE branch from discriminant (R version)
#' @noRd
.classify_disc <- function(b1, b2) {
  eps <- 1e-8
  disc_tol <- 1e-12
  if (abs(b1) < eps) {
    return(if (abs(b2) > eps) 1L else 4L) # FIRST_ORD or ZERO
  }
  D <- b2^2 + 4 * b1
  if (D > disc_tol) return(0L)   # REAL
  if (D < -disc_tol) return(2L)  # COMPLEX
  return(3L)                       # REPEATED
}

#' @importFrom stats model.frame model.matrix model.response
#' @noRd
#' @importFrom stats model.frame model.matrix model.response
#' @noRd
.process_marginal <- function(formula, data, time, id, parsed_long) {
  if (is.matrix(data)) data <- as.data.frame(data)

  fixed_formula <- .build_formula(
    parsed_long$fixed_terms, response = parsed_long$response
  )
  random_formula <- if (!is.null(parsed_long$random_terms)) {
    .build_formula(parsed_long$random_terms, is_random = TRUE)
  } else {
    NULL
  }

  unique_ids <- unique(data[[id]])
  n_subjects <- length(unique_ids)
  id_groups <- split(seq_len(nrow(data)), data[[id]])

  subject_data <- vector("list", n_subjects)
  names(subject_data) <- as.character(unique_ids)

  for (i in seq_along(unique_ids)) {
    sid <- as.character(unique_ids[i])
    rows <- id_groups[[sid]]

    if (!is.null(rows) && length(rows) > 0) {
      subset <- data[rows, , drop = FALSE]
      subset <- subset[order(subset[[time]]), , drop = FALSE]
      long_frame <- model.frame(fixed_formula, data = subset)
      long_times <- subset[[time]]
      long_measurements <- model.response(long_frame)
      long_fixed_covariates <- model.matrix(fixed_formula, long_frame)
      long_random_covariates <- if (is.null(random_formula)) {
        matrix(nrow = nrow(subset), ncol = 0)
      } else {
        model.matrix(random_formula, subset)
      }
    } else {
      long_times <- numeric(0)
      long_measurements <- numeric(0)
      empty <- data[0, , drop = FALSE]
      n_fixed <- ncol(model.matrix(fixed_formula, empty))
      long_fixed_covariates <- matrix(nrow = 0, ncol = n_fixed)
      long_random_covariates <- if (is.null(random_formula)) {
        matrix(nrow = 0, ncol = 0)
      } else {
        matrix(nrow = 0, ncol = ncol(model.matrix(random_formula, empty)))
      }
    }

    subject_data[[i]] <- list(
      longitudinal = list(
        times = long_times,
        measurements = long_measurements,
        covariates = list(
          fixed = long_fixed_covariates,
          random = long_random_covariates
        )
      )
    )
  }

  all_y <- unlist(lapply(subject_data, function(d) d$longitudinal$measurements))
  attr(subject_data, "biomarker_clamp") <- max(abs(all_y)) * 5
  subject_data
}

#' Pack marginal data into flat TMB input
#' @noRd
.pack_marginal_data <- function(data_list, parsed_long, n_re) {
  n_subjects <- length(data_list)
  clamp_value <- attr(data_list, "biomarker_clamp")

  n_observations <- vapply(data_list, function(d) {
    length(d$longitudinal$measurements)
  }, integer(1))

  obs_times <- unlist(lapply(data_list, function(d) d$longitudinal$times))
  obs_values <- unlist(lapply(data_list, function(d) d$longitudinal$measurements))

  n_fixed_covariates <- ncol(data_list[[1]]$longitudinal$covariates$fixed)
  n_random_covariates <- ncol(data_list[[1]]$longitudinal$covariates$random)

  flatten_design <- function(type) {
    result <- unlist(lapply(data_list, function(d)
      as.vector(t(d$longitudinal$covariates[[type]]))))
    result %||% numeric(0)
  }

  list(
    model_type = 1L,
    n_subjects = as.integer(n_subjects),
    n_random_effects = as.integer(n_re),
    n_observations = as.integer(n_observations),
    obs_times = obs_times %||% numeric(0),
    obs_values = obs_values %||% numeric(0),
    n_fixed_covariates = as.integer(n_fixed_covariates),
    n_random_covariates = as.integer(n_random_covariates),
    long_fixed_covariates = flatten_design("fixed"),
    long_random_covariates = flatten_design("random"),
    ode_branch = 0L,  # BR_REAL default
    clamp_value = clamp_value,
    biomarker_fixed = as.integer(parsed_long$biomarker$fixed),
    biomarker_random = as.integer(parsed_long$biomarker$random),
    velocity_fixed = as.integer(parsed_long$velocity$fixed),
    velocity_random = as.integer(parsed_long$velocity$random)
  )
}

#' Build marginal TMB parameter list
#' @noRd
.pack_marginal_params <- function(n_long_coef, n_re, n_subjects,
                                  mean_y, sd_y) {
  # Decompose identity Sigma_b
  marginal_sds <- rep(1, n_re)
  n_corr <- n_re * (n_re - 1) / 2

  list(
    longitudinal = rep(0, n_long_coef),
    initial_state = c(mean_y, 0),
    log_sigma_e = log(sd_y),
    log_sd_re = log(marginal_sds),
    corr_par = rep(0, n_corr),
    random_effects = matrix(0, nrow = n_subjects, ncol = n_re)
  )
}

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
  sigma_e <- exp(par[pn == "log_sigma_e"])

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

# TMB Result Extraction ========================================================

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
  cf$measurement_error_sd <- exp(par[pn == "log_sigma_e"])
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
