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

#' Build TMB data list from processed data
#' @noRd
.build_tmb_data <- function(data_list, parameters, control) {
  n_subjects <- length(data_list)
  configs <- parameters$configurations
  coefs <- parameters$coefficients

  # Count observations per subject
  n_obs <- vapply(data_list, function(d) {
    length(d$longitudinal$measurements)
  }, integer(1))

  # Concatenate observation times and values
  obs_times_all <- unlist(lapply(data_list, function(d) d$longitudinal$times))
  obs_values_all <- unlist(lapply(data_list, function(d) d$longitudinal$measurements))
  if (is.null(obs_times_all)) obs_times_all <- numeric(0)
  if (is.null(obs_values_all)) obs_values_all <- numeric(0)

  # Covariate dimensions per subject
  cov_dims <- t(vapply(data_list, function(d) {
    c(ncol(d$longitudinal$covariates$fixed),
      ncol(d$longitudinal$covariates$random))
  }, integer(2)))

  # Concatenate covariates (row-major flat)
  X_fixed_all <- unlist(lapply(data_list, function(d) {
    as.vector(t(d$longitudinal$covariates$fixed))
  }))
  X_random_all <- unlist(lapply(data_list, function(d) {
    as.vector(t(d$longitudinal$covariates$random))
  }))
  if (is.null(X_fixed_all)) X_fixed_all <- numeric(0)
  if (is.null(X_random_all)) X_random_all <- numeric(0)

  # Event data
  event_times <- vapply(data_list, `[[`, numeric(1), "time")
  event_status <- vapply(data_list, `[[`, numeric(1), "status")

  # Survival covariates
  n_surv_cov <- if (is.data.frame(data_list[[1]]$covariates) &&
                    ncol(data_list[[1]]$covariates) > 0) {
    ncol(data_list[[1]]$covariates)
  } else if (is.matrix(data_list[[1]]$covariates)) {
    ncol(data_list[[1]]$covariates)
  } else {
    0L
  }

  W_surv_all <- if (n_surv_cov > 0) {
    unlist(lapply(data_list, function(d) as.numeric(d$covariates)))
  } else {
    numeric(0)
  }

  # Classify ODE branch from initial b1/b2
  b1_val <- if (configs$biomarker$fixed) coefs$longitudinal[1] else 0
  b2_val <- if (configs$velocity$fixed) {
    coefs$longitudinal[1 + configs$biomarker$fixed]
  } else {
    0
  }
  branches <- rep(.classify_disc(b1_val, b2_val), n_subjects)

  sbc <- configs$baseline

  list(
    n_subjects = as.integer(n_subjects),
    n_random_effects = as.integer(ncol(parameters$random_effects_init)),
    n_observations = as.integer(n_obs),
    event_times = event_times,
    event_status = as.integer(event_status),
    obs_times_all = obs_times_all,
    obs_values_all = obs_values_all,
    covariate_dims = matrix(as.integer(cov_dims), nrow = n_subjects, ncol = 2),
    X_fixed_all = X_fixed_all,
    X_random_all = X_random_all,
    W_survival_all = W_surv_all,
    n_survival_covariates = as.integer(n_surv_cov),
    ode_branch = as.integer(branches),
    biomarker_clamp = configs$biomarker_clamp,
    hazard_quadrature = as.integer(control$hazard_quadrature),
    gamma = configs$gamma,
    spline_degree = as.integer(sbc$degree),
    spline_knots = if (length(sbc$knots) > 0) sbc$knots else numeric(0),
    spline_boundary = sbc$boundary_knots,
    biomarker_fixed = as.integer(configs$biomarker$fixed),
    biomarker_random = as.integer(configs$biomarker$random),
    velocity_fixed = as.integer(configs$velocity$fixed),
    velocity_random = as.integer(configs$velocity$random)
  )
}

#' Build TMB parameter list
#' @noRd
.build_tmb_parameters <- function(parameters) {
  coefs <- parameters$coefficients
  n_re <- ncol(parameters$random_effects_init)
  Sigma <- coefs$random_effect_sigma

  # Decompose Sigma_b into SDs + correlation parameters
  sds <- sqrt(diag(Sigma))
  # Correlation matrix
  D_inv <- diag(1 / pmax(sds, 1e-10))
  R <- D_inv %*% Sigma %*% D_inv
  # Off-diagonal correlation -> theta parameterization: theta = rho / sqrt(1 - rho^2)
  corr_par <- numeric(n_re * (n_re - 1) / 2)
  idx <- 1L
  for (j in seq_len(n_re - 1)) {
    for (i_row in (j + 1):n_re) {
      rho <- R[i_row, j]
      rho <- max(-0.99, min(0.99, rho))
      corr_par[idx] <- rho / sqrt(1 - rho^2)
      idx <- idx + 1L
    }
  }

  list(
    baseline = coefs$baseline,
    hazard = coefs$hazard,
    longitudinal = coefs$longitudinal,
    initial_state = coefs$initial_state,
    log_sigma_e = log(coefs$measurement_error_sd),
    log_sd_re = log(pmax(sds, 1e-10)),
    corr_par = corr_par,
    b = parameters$random_effects_init
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
.process_marginal <- function(formula, data, time, id) {
  if (is.matrix(data)) data <- as.data.frame(data)
  stopifnot(
    "Data cannot be empty" = nrow(data) > 0,
    "Data must contain a time column" = time %in% names(data),
    "Data must contain an id column" = id %in% names(data)
  )

  mf <- model.frame(formula, data = data, na.action = na.omit)
  y <- model.response(mf)
  X <- model.matrix(formula, data = mf) # nolint: object_name_linter
  stopifnot("Formula must include a response" = !is.null(y))

  row_idx <- as.numeric(rownames(mf))
  times <- data[[time]][row_idx]
  ids <- data[[id]][row_idx]
  subjects <- unique(ids)

  subject_data <- lapply(seq_along(subjects), function(i) {
    idx <- which(ids == subjects[i])
    idx <- idx[order(times[idx])]
    t_subj <- times[idx]
    list(
      time = max(t_subj),
      initial_state = c(0, 0),
      longitudinal = list(
        times = t_subj,
        measurements = y[idx],
        covariates = list(
          fixed = X[idx, , drop = FALSE],
          random = matrix(nrow = length(idx), ncol = 0)
        )
      )
    )
  })

  names(subject_data) <- as.character(subjects)
  attr(subject_data, "n_covariates") <- ncol(X)
  attr(subject_data, "covariate_names") <- colnames(X)
  attr(subject_data, "biomarker_clamp") <- max(abs(y)) * 5
  subject_data
}
