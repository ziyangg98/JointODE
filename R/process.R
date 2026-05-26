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

#' Decompose Sigma_b into log(SD) and correlation theta parameters
#' @noRd
.pack_correlation_theta <- function(sigma_matrix, n_re) {
  marginal_sds <- sqrt(diag(sigma_matrix))
  scale_inv <- diag(1 / pmax(marginal_sds, 1e-10))
  corr_matrix <- scale_inv %*% sigma_matrix %*% scale_inv

  corr_theta <- numeric(n_re * (n_re - 1) / 2)
  idx <- 1L
  for (row in 2:n_re) {
    for (col in seq_len(row - 1L)) {
      rho <- max(-0.99, min(0.99, corr_matrix[row, col]))
      corr_theta[idx] <- rho / sqrt(1 - rho^2)
      idx <- idx + 1L
    }
  }

  list(
    log_sd_re = log(pmax(marginal_sds, 1e-10)),
    corr_par = corr_theta
  )
}

#' Pack per-subject data into flat TMB input
#' @noRd
.pack_joint_data <- function(data_list, parameters, control) {
  n_subjects <- length(data_list)
  configs <- parameters$configurations
  coefs <- parameters$coefficients

  # Concatenate per-subject longitudinal fields
  extract_long <- function(field) {
    lapply(data_list, function(d) d$longitudinal[[field]])
  }
  obs_times <- unlist(extract_long("times"))
  obs_values <- unlist(extract_long("measurements"))
  n_observations <- vapply(data_list, function(d) {
    length(d$longitudinal$measurements)
  }, integer(1))

  # Covariate counts (constant across subjects)
  n_fixed_covariates <- ncol(data_list[[1]]$longitudinal$covariates$fixed)
  n_random_covariates <- ncol(data_list[[1]]$longitudinal$covariates$random)

  # Flatten design matrices (row-major)
  flatten_design <- function(type) {
    result <- unlist(lapply(data_list, function(d) {
      as.vector(t(d$longitudinal$covariates[[type]]))
    }))
    if (length(result) == 0) numeric(0) else result
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
    velocity_random = as.integer(configs$velocity$random),
    diagonal_re = as.integer(isTRUE(configs$covariance == "diagonal"))
  )
}

#' Convert R parameters to TMB parameterization
#' @noRd
.pack_joint_params <- function(parameters) {
  coefs <- parameters$coefficients
  n_random_effects <- ncol(parameters$random_effects_init)
  corr <- .pack_correlation_theta(coefs$random_effect_sigma, n_random_effects)

  list(
    baseline = coefs$baseline,
    hazard = coefs$hazard,
    longitudinal = coefs$longitudinal,
    initial_state = coefs$initial_state,
    log_sigma_e = log(coefs$measurement_error_sd),
    log_sd_re = corr$log_sd_re,
    corr_par = corr$corr_par,
    random_effects = parameters$random_effects_init
  )
}

#' @importFrom stats model.frame model.matrix model.response
#' @noRd
.process_marginal <- function(formula, data, time, id, parsed_long) {
  if (is.matrix(data)) data <- as.data.frame(data)

  fixed_formula <- .build_formula(
    parsed_long$fixed_terms,
    response = parsed_long$response
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

  subject_data
}

#' Pack marginal data into flat TMB input
#' @noRd
.pack_marginal_data <- function(data_list, parsed_long, n_re) {
  n_subjects <- length(data_list)

  n_observations <- vapply(data_list, function(d) {
    length(d$longitudinal$measurements)
  }, integer(1))

  obs_times <- unlist(lapply(data_list, function(d) d$longitudinal$times))
  obs_values <- unlist(lapply(data_list, function(d) d$longitudinal$measurements))

  n_fixed_covariates <- ncol(data_list[[1]]$longitudinal$covariates$fixed)
  n_random_covariates <- ncol(data_list[[1]]$longitudinal$covariates$random)

  flatten_design <- function(type) {
    result <- unlist(lapply(data_list, function(d) {
      as.vector(t(d$longitudinal$covariates[[type]]))
    }))
    if (length(result) == 0) numeric(0) else result
  }

  list(
    model_type = 1L,
    n_subjects = as.integer(n_subjects),
    n_random_effects = as.integer(n_re),
    n_observations = as.integer(n_observations),
    obs_times = if (length(obs_times) == 0) numeric(0) else obs_times,
    obs_values = if (length(obs_values) == 0) numeric(0) else obs_values,
    n_fixed_covariates = as.integer(n_fixed_covariates),
    n_random_covariates = as.integer(n_random_covariates),
    long_fixed_covariates = flatten_design("fixed"),
    long_random_covariates = flatten_design("random"),
    biomarker_fixed = as.integer(parsed_long$biomarker$fixed),
    biomarker_random = as.integer(parsed_long$biomarker$random),
    velocity_fixed = as.integer(parsed_long$velocity$fixed),
    velocity_random = as.integer(parsed_long$velocity$random),
    diagonal_re = as.integer(isTRUE(parsed_long$diagonal))
  )
}

#' Default parameters for MarginalODE
#' @noRd
.default_marginal_parameters <- function(model_config, parsed_long) {
  longitudinal <- rep(0, model_config$n_longitudinal_coef)

  list(
    coefficients = list(
      longitudinal = longitudinal,
      initial_state = c(0, 0),
      measurement_error_sd = 1,
      random_effect_sigma = diag(1, model_config$n_re)
    ),
    configurations = list(
      biomarker = list(
        fixed = parsed_long$biomarker$fixed,
        random = parsed_long$biomarker$random
      ),
      velocity = list(
        fixed = parsed_long$velocity$fixed,
        random = parsed_long$velocity$random
      ),
      covariance = if (isTRUE(parsed_long$diagonal)) "diagonal" else "full"
    )
  )
}

#' Convert marginal parameters to TMB parameterization
#' @noRd
.pack_marginal_params <- function(parameters) {
  coefs <- parameters$coefficients
  n_re <- ncol(parameters$random_effects_init)
  corr <- .pack_correlation_theta(coefs$random_effect_sigma, n_re)

  list(
    longitudinal = coefs$longitudinal,
    initial_state = coefs$initial_state,
    log_sigma_e = log(coefs$measurement_error_sd),
    log_sd_re = corr$log_sd_re,
    corr_par = corr$corr_par,
    random_effects = parameters$random_effects_init
  )
}
