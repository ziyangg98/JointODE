# Utility Functions for JointODE Package

#' @importFrom stats setNames model.frame model.response model.matrix quantile
#' @importFrom utils head
#' @importFrom parallel detectCores
#' @importFrom future plan multicore multisession sequential
#' @importFrom future.apply future_lapply
NULL

# Null coalescing operator
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Setup parallel plan for computation
.setup_parallel_plan <- function(n_cores = 0) {
  # Auto-detect cores if not specified
  if (n_cores == 0) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  # Set up future plan based on platform
  if (.Platform$OS.type == "unix") {
    future::plan(future::multicore, workers = n_cores)
  } else {
    future::plan(future::multisession, workers = n_cores)
  }
  # Return cleanup function
  function() future::plan(future::sequential)
}

# ===== SECTION 1: DATA VALIDATION AND PROCESSING =====
.validate <- function(
  longitudinal_formula,
  survival_formula,
  longitudinal_data,
  survival_data,
  state,
  id,
  time,
  spline_baseline,
  init
) {
  stopifnot(
    "longitudinal_formula must be a formula" = inherits(
      longitudinal_formula,
      "formula"
    ),
    "longitudinal_data must be a data frame" = is.data.frame(longitudinal_data),
    "survival_formula must be a formula" = inherits(
      survival_formula,
      "formula"
    ),
    "survival_data must be a data frame" = is.data.frame(survival_data),
    "id must be a single character string" = is.character(id) &&
      length(id) == 1,
    "time must be a single character string" = is.character(time) &&
      length(time) == 1,
    "spline_baseline must be a list" = is.list(spline_baseline)
  )

  stopifnot(
    "Longitudinal data has no rows" = nrow(longitudinal_data) > 0,
    "Survival data has no rows" = nrow(survival_data) > 0
  )

  if (!id %in% names(longitudinal_data)) {
    stop(sprintf("ID variable '%s' not found in longitudinal data", id))
  }
  if (!time %in% names(longitudinal_data)) {
    stop(sprintf("Time variable '%s' not found in longitudinal data", time))
  }

  if (!id %in% names(survival_data)) {
    stop(sprintf("ID variable '%s' not found in survival data", id))
  }

  surv_response <- survival_formula[[2]]
  if (
    !inherits(surv_response, "call") ||
      !identical(as.character(surv_response[[1]]), "Surv")
  ) {
    stop("Survival formula must have Surv() on the left-hand side")
  }

  surv_call_vars <- all.vars(surv_response)
  if (length(surv_call_vars) < 2) {
    stop("Surv() must have at least time and status arguments")
  }

  time_var <- surv_call_vars[1]
  status_var <- surv_call_vars[2]

  long_vars <- all.vars(longitudinal_formula)
  missing_long_vars <- setdiff(long_vars, names(longitudinal_data))
  if (length(missing_long_vars) > 0) {
    stop(sprintf(
      "Variables in longitudinal formula not found in data: %s",
      paste(missing_long_vars, collapse = ", ")
    ))
  }

  surv_vars <- all.vars(survival_formula)
  missing_surv_vars <- setdiff(surv_vars, names(survival_data))
  if (length(missing_surv_vars) > 0) {
    stop(sprintf(
      "Variables in survival formula not found in data: %s",
      paste(missing_surv_vars, collapse = ", ")
    ))
  }

  critical_cols <- list(
    "ID in longitudinal data" = longitudinal_data[[id]],
    "Time in longitudinal data" = longitudinal_data[[time]],
    "ID in survival data" = survival_data[[id]]
  )

  for (col_name in names(critical_cols)) {
    if (any(is.na(critical_cols[[col_name]]))) {
      stop(sprintf("Missing values found in %s", col_name))
    }
  }

  if (any(duplicated(survival_data[[id]]))) {
    stop(paste(
      "Duplicate IDs found in survival data -",
      "each subject should have one record"
    ))
  }

  long_ids <- unique(longitudinal_data[[id]])
  surv_ids <- unique(survival_data[[id]])

  orphaned_ids <- setdiff(long_ids, surv_ids)
  if (length(orphaned_ids) > 0) {
    stop(sprintf(
      "Subjects in longitudinal data not found in survival data: %s",
      paste(head(orphaned_ids, 5), collapse = ", ")
    ))
  }

  missing_long <- setdiff(surv_ids, long_ids)
  if (length(missing_long) > 0) {
    warning(sprintf(
      "Subjects in survival data without longitudinal data: %s",
      paste(head(missing_long, 5), collapse = ", ")
    ))
  }

  if (any(longitudinal_data[[time]] < 0, na.rm = TRUE)) {
    stop("Negative time values found in longitudinal data")
  }

  if (any(survival_data[[time_var]] <= 0, na.rm = TRUE)) {
    stop("Invalid observation times in survival data (must be positive)")
  }

  # Pre-compute grouped data for efficiency
  long_times_by_id <- split(longitudinal_data[[time]], longitudinal_data[[id]])
  surv_times_map <- setNames(survival_data[[time_var]], survival_data[[id]])

  subjects_with_late_obs <- character()
  for (subject_id in long_ids) {
    if (subject_id %in% surv_ids) {
      long_times <- long_times_by_id[[as.character(subject_id)]]
      surv_time <- surv_times_map[as.character(subject_id)]

      if (!is.na(surv_time) && any(long_times > surv_time + 1e-6)) {
        subjects_with_late_obs <- c(subjects_with_late_obs, subject_id)
      }
    }
  }

  if (length(subjects_with_late_obs) > 0) {
    warning(sprintf(
      "%d subjects have measurements after observation time: %s",
      length(subjects_with_late_obs),
      paste(head(subjects_with_late_obs, 3), collapse = ", ")
    ))
  }

  unique_status <- unique(survival_data[[status_var]])
  invalid_status <- setdiff(unique_status, c(0, 1, NA))

  if (length(invalid_status) > 0) {
    stop(sprintf(
      "Invalid status values found: %s. Must be 0 (censored) or 1 (event)",
      paste(invalid_status, collapse = ", ")
    ))
  }

  if (any(is.na(survival_data[[status_var]]))) {
    stop("Missing values in status variable")
  }

  observations_per_subject <- table(longitudinal_data[[id]])
  if (all(observations_per_subject == 1)) {
    warning(paste(
      "Each subject has only one longitudinal observation -",
      "joint modeling may not be appropriate"
    ))
  }

  valid_baseline_params <- c(
    "degree",
    "n_knots",
    "knot_placement",
    "boundary_knots"
  )
  invalid_baseline <- setdiff(names(spline_baseline), valid_baseline_params)
  if (length(invalid_baseline) > 0) {
    stop(sprintf(
      "Invalid parameters in spline_baseline: %s. Valid parameters are: %s",
      paste(invalid_baseline, collapse = ", "),
      paste(valid_baseline_params, collapse = ", ")
    ))
  }

  if ("degree" %in% names(spline_baseline)) {
    if (
      !is.numeric(spline_baseline$degree) ||
        length(spline_baseline$degree) != 1 ||
        spline_baseline$degree < 1 ||
        spline_baseline$degree > 5
    ) {
      stop("spline_baseline$degree must be a single integer between 1 and 5")
    }
  }

  if ("n_knots" %in% names(spline_baseline)) {
    if (
      !is.numeric(spline_baseline$n_knots) ||
        length(spline_baseline$n_knots) != 1 ||
        spline_baseline$n_knots < 0 ||
        spline_baseline$n_knots > 20
    ) {
      stop("spline_baseline$n_knots must be a single integer between 0 and 20")
    }
  }

  if ("knot_placement" %in% names(spline_baseline)) {
    valid_placements <- c("quantile", "equal")
    if (!spline_baseline$knot_placement %in% valid_placements) {
      stop(sprintf(
        "spline_baseline$knot_placement must be one of: %s",
        paste(valid_placements, collapse = ", ")
      ))
    }
  }

  if ("boundary_knots" %in% names(spline_baseline)) {
    if (!is.null(spline_baseline$boundary_knots)) {
      if (
        !is.numeric(spline_baseline$boundary_knots) ||
          length(spline_baseline$boundary_knots) != 2
      ) {
        stop(paste(
          "spline_baseline$boundary_knots must be NULL or",
          "a numeric vector of length 2"
        ))
      }
      if (
        spline_baseline$boundary_knots[1] >= spline_baseline$boundary_knots[2]
      ) {
        stop(paste(
          "spline_baseline$boundary_knots[1] must be less than",
          "boundary_knots[2]"
        ))
      }
    }
  }

  # Validate state parameter if provided
  if (!is.null(state)) {
    if (!is.matrix(state)) {
      stop("'state' must be a matrix")
    }

    # Check dimensions
    n_subjects_survival <- nrow(survival_data)
    if (nrow(state) != n_subjects_survival) {
      stop(sprintf(
        paste(
          "Invalid 'state': number of rows (%d)",
          "must match number of subjects (%d)"
        ),
        nrow(state),
        n_subjects_survival
      ))
    }

    if (ncol(state) != 2) {
      stop(sprintf(
        paste(
          "Invalid 'state': must have exactly 2 columns",
          "(biomarker and velocity), got %d"
        ),
        ncol(state)
      ))
    }

    # Check for non-finite values
    if (any(!is.finite(state))) {
      stop("Invalid 'state': all values must be finite")
    }
  }

  # Validate init parameter if provided
  if (!is.null(init)) {
    # Calculate dimensions needed for init validation
    # Get number of longitudinal covariates
    long_model_matrix <- model.matrix(longitudinal_formula, longitudinal_data)
    n_longitudinal_covariates <- ncol(long_model_matrix)

    # Get number of survival covariates
    surv_model_matrix <- model.matrix(survival_formula, survival_data)
    n_survival_covariates <- ncol(surv_model_matrix) - 1 # exclude intercept

    # Calculate spline basis dimension
    spline_config <- modifyList(
      list(degree = 3, n_knots = 5),
      spline_baseline
    )
    n_spline_basis <- spline_config$degree + spline_config$n_knots + 1

    if (!is.list(init)) {
      stop("Invalid 'init' parameter: must be a list")
    }

    # Check for unknown top-level components
    valid_components <- c("coefficients", "configurations")
    unknown_components <- setdiff(names(init), valid_components)
    if (length(unknown_components) > 0) {
      stop(sprintf(
        "Invalid 'init': unknown components '%s' (valid: %s)",
        paste(unknown_components, collapse = "', '"),
        paste(valid_components, collapse = ", ")
      ))
    }

    if (!is.null(init$coefficients)) {
      if (!is.list(init$coefficients)) {
        stop("Invalid 'init$coefficients': must be a list")
      }

      # Check for unknown coefficient types
      valid_coefs <- c(
        "baseline",
        "hazard",
        "acceleration",
        "measurement_error_sd",
        "random_effect_sd"
      )
      unknown_coefs <- setdiff(names(init$coefficients), valid_coefs)
      if (length(unknown_coefs) > 0) {
        stop(sprintf(
          "Invalid 'init$coefficients': unknown types '%s' (valid: %s)",
          paste(unknown_coefs, collapse = ", "),
          paste(valid_coefs, collapse = ", ")
        ))
      }

      # Validate baseline coefficients
      if (!is.null(init$coefficients$baseline)) {
        if (!is.numeric(init$coefficients$baseline)) {
          stop("Invalid 'init$coefficients$baseline': must be numeric")
        }
        if (any(!is.finite(init$coefficients$baseline))) {
          stop(paste(
            "Invalid 'init$coefficients$baseline':",
            "must contain finite values"
          ))
        }
        # Check length
        if (length(init$coefficients$baseline) != n_spline_basis) {
          stop(sprintf(
            paste(
              "Invalid 'init$coefficients$baseline':",
              "wrong length (expected %d, got %d)"
            ),
            n_spline_basis,
            length(init$coefficients$baseline)
          ))
        }
      }

      # Validate hazard coefficients
      if (!is.null(init$coefficients$hazard)) {
        if (!is.numeric(init$coefficients$hazard)) {
          stop("Invalid 'init$coefficients$hazard': must be numeric")
        }
        if (any(!is.finite(init$coefficients$hazard))) {
          stop("Invalid 'init$coefficients$hazard': must contain finite values")
        }
        if (length(init$coefficients$hazard) < 2) {
          stop(paste(
            "Invalid 'init$coefficients$hazard':",
            "must have at least 2 elements for association parameters"
          ))
        }
        # Check exact length
        expected_len <- n_survival_covariates + 2
        if (length(init$coefficients$hazard) != expected_len) {
          stop(sprintf(
            paste(
              "Invalid 'init$coefficients$hazard':",
              "wrong length (expected %d, got %d)"
            ),
            expected_len,
            length(init$coefficients$hazard)
          ))
        }
      }

      # Validate acceleration coefficients
      if (!is.null(init$coefficients$acceleration)) {
        if (!is.numeric(init$coefficients$acceleration)) {
          stop(paste(
            "Invalid 'init$coefficients$acceleration':",
            "must be numeric"
          ))
        }
        if (any(!is.finite(init$coefficients$acceleration))) {
          stop(paste(
            "Invalid 'init$coefficients$acceleration':",
            "must contain finite values"
          ))
        }
        if (length(init$coefficients$acceleration) < 2) {
          stop(paste(
            "Invalid 'init$coefficients$acceleration':",
            "must have at least 2 elements"
          ))
        }
        # Check exact length
        expected_len <- n_longitudinal_covariates + 2
        if (length(init$coefficients$acceleration) != expected_len) {
          stop(sprintf(
            paste(
              "Invalid 'init$coefficients$acceleration':",
              "wrong length (expected %d, got %d)"
            ),
            expected_len,
            length(init$coefficients$acceleration)
          ))
        }
      }

      # Validate measurement error SD
      if (!is.null(init$coefficients$measurement_error_sd)) {
        if (
          !is.numeric(init$coefficients$measurement_error_sd) ||
            length(init$coefficients$measurement_error_sd) != 1
        ) {
          stop(paste(
            "Invalid 'init$coefficients$measurement_error_sd':",
            "must be a single numeric value"
          ))
        }
        if (!is.finite(init$coefficients$measurement_error_sd)) {
          stop(paste(
            "Invalid 'init$coefficients$measurement_error_sd':",
            "must be finite"
          ))
        }
        if (init$coefficients$measurement_error_sd <= 0) {
          stop(paste(
            "Invalid 'init$coefficients$measurement_error_sd':",
            "must be positive"
          ))
        }
      }

      # Validate random effect SD
      if (!is.null(init$coefficients$random_effect_sd)) {
        if (
          !is.numeric(init$coefficients$random_effect_sd) ||
            length(init$coefficients$random_effect_sd) != 1
        ) {
          stop(paste(
            "Invalid 'init$coefficients$random_effect_sd':",
            "must be a single numeric value"
          ))
        }
        if (!is.finite(init$coefficients$random_effect_sd)) {
          stop("Invalid 'init$coefficients$random_effect_sd': must be finite")
        }
        if (init$coefficients$random_effect_sd <= 0) {
          stop("Invalid 'init$coefficients$random_effect_sd': must be positive")
        }
      }
    }

    # Validate configurations if provided
    if (!is.null(init$configurations)) {
      if (!is.list(init$configurations)) {
        stop("Invalid 'init$configurations': must be a list")
      }
      if (!is.null(init$configurations$baseline)) {
        if (!is.list(init$configurations$baseline)) {
          stop("Invalid 'init$configurations$baseline': must be a list")
        }
        # Could add more specific validation for spline configuration here
      }
    }
  }

  invisible(NULL)
}

.process <- function(
  longitudinal_formula,
  survival_formula,
  longitudinal_data,
  survival_data,
  state,
  id,
  time
) {
  unique_ids <- unique(survival_data[[id]])
  n_subjects <- length(unique_ids)
  data_process <- vector("list", n_subjects)
  names(data_process) <- as.character(unique_ids)

  surv_frame <- model.frame(survival_formula, data = survival_data)
  surv_response_matrix <- model.response(surv_frame)
  has_surv_covs <- length(all.vars(survival_formula[[3]])) > 0 &&
    survival_formula[[3]] != 1
  surv_design <- if (has_surv_covs) {
    model.matrix(survival_formula, surv_frame)[, -1, drop = FALSE]
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
      ord <- order(long_subset[[time]])
      if (!identical(ord, seq_along(ord))) {
        long_subset <- long_subset[ord, , drop = FALSE]
      }
      long_frame <- model.frame(longitudinal_formula, data = long_subset)
      long_times <- long_subset[[time]]
      long_measurements <- model.response(long_frame)
      long_covariates <- as.matrix(
        model.matrix(longitudinal_formula, long_frame)
      )
    } else {
      long_times <- numeric(0)
      long_measurements <- numeric(0)
      long_covariates <- matrix()
    }

    survival_row <- survival_index_map[i]
    event_time <- surv_response_matrix[survival_row, 1]
    event_status <- surv_response_matrix[survival_row, 2]
    covariates <- if (!is.null(surv_design)) {
      surv_design[survival_row, , drop = FALSE]
    } else {
      data.frame()
    }

    # Get initial state for this subject if provided
    initial_state <- if (!is.null(state)) {
      state[i, , drop = TRUE] # Extract row i as vector
    } else {
      c(0, 0) # Default initial state
    }

    data_process[[i]] <- list(
      id = unique_ids[i],
      time = event_time,
      status = event_status,
      covariates = covariates,
      initial_state = initial_state,
      longitudinal = list(
        times = long_times,
        measurements = long_measurements,
        covariates = long_covariates,
        n_obs = length(long_times)
      )
    )
  }

  statuses <- vapply(data_process, `[[`, numeric(1), "status")
  attr(data_process, "n_subjects") <- n_subjects
  attr(data_process, "n_observations") <- nrow(longitudinal_data)
  attr(data_process, "event_rate") <- mean(statuses == 1, na.rm = TRUE)

  data_process
}

.get_spline_config <- function(
  x,
  degree = 3,
  n_knots = 5,
  knot_placement = "quantile",
  boundary_knots = NULL
) {
  if (is.null(boundary_knots)) {
    boundary_knots <- range(x, na.rm = TRUE)
  } else {
    boundary_knots <- boundary_knots
  }

  if (knot_placement == "quantile") {
    probs <- seq(0, 1, length.out = n_knots + 2)[-c(1, n_knots + 2)]
    knots <- quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  } else if (knot_placement == "equal") {
    knots <- seq(
      boundary_knots[1],
      boundary_knots[2],
      length.out = n_knots + 2
    )[-c(1, n_knots + 2)]
  } else {
    stop("knot_placement must be 'quantile' or 'equal'")
  }

  list(
    degree = degree,
    knots = knots,
    boundary_knots = boundary_knots,
    df = length(knots) + degree + 1
  )
}

.compute_spline_basis <- function(x, config) {
  splines2::bSpline(
    x = x,
    knots = config$knots,
    degree = config$degree,
    Boundary.knots = config$boundary_knots,
    intercept = TRUE,
    warn.outside = FALSE
  )
}

.compute_spline_basis_deriv <- function(x, config) {
  splines2::dbs(
    x = x,
    knots = config$knots,
    degree = config$degree,
    Boundary.knots = config$boundary_knots,
    derivs = 1,
    intercept = TRUE,
    warn.outside = FALSE
  )
}


.get_longitudinal_covariates <- function(data, time) {
  long_cov <- data$longitudinal$covariates
  if (is.null(long_cov) || nrow(long_cov) == 0) {
    return(numeric(0))
  }
  row_index <- which.min(abs(data$longitudinal$times - time))
  as.numeric(long_cov[row_index, , drop = TRUE])
}

.compute_acceleration <- function(
  biomarker,
  velocity,
  time,
  data,
  parameters
) {
  long_cov <- .get_longitudinal_covariates(data, time)
  # Construct z vector
  z <- c(biomarker, velocity, long_cov)
  sum(parameters$coefficients$acceleration * z)
}

.compute_acceleration_deriv <- function(
  biomarker,
  velocity,
  time,
  data,
  parameters,
  dbiomarker = NULL,
  dvelocity = NULL
) {
  long_cov <- .get_longitudinal_covariates(data, time)

  # Construct Z vector
  z_vec <- c(biomarker, velocity, long_cov)
  beta <- parameters$coefficients$acceleration

  # For linear model: acceleration = beta' * Z
  # d(acceleration)/d(beta_j) = Z_j + beta' * dZ/d(beta_j)

  n_beta <- length(beta)
  dacceleration_dbeta <- numeric(n_beta)

  for (j in seq_len(n_beta)) {
    # Direct effect: Z_j
    direct_effect <- if (j <= length(z_vec)) z_vec[j] else 0

    # Indirect effect through state feedback
    # Only the first two components of Z (biomarker and velocity) depend on beta
    indirect_effect <- 0
    if (!is.null(dbiomarker) && !is.null(dvelocity)) {
      # All beta parameters affect the dynamics, so all have indirect effects
      indirect_effect <- beta[1] * dbiomarker[j] + beta[2] * dvelocity[j]
    }

    dacceleration_dbeta[j] <- direct_effect + indirect_effect
  }

  return(dacceleration_dbeta)
}

.compute_log_hazard <- function(
  time,
  biomarker,
  velocity,
  data,
  parameters
) {
  # Baseline hazard
  basis_lambda <- .compute_spline_basis(
    time,
    parameters$configurations$baseline
  )
  log_baseline <- sum(basis_lambda * parameters$coefficients$baseline)

  # Biomarker effects: m_i(t)' * α
  biomarker_vec <- c(biomarker, velocity)
  log_biomarker <- sum(biomarker_vec * parameters$coefficients$hazard[1:2])

  # Covariate effects: W_i' * φ
  n_hazard <- length(parameters$coefficients$hazard)
  log_covariate <- if (n_hazard > 2 && !is.null(data$covariates)) {
    w <- data$covariates
    if (nrow(w) > 0) {
      sum(parameters$coefficients$hazard[3:n_hazard] * w)
    } else {
      0
    }
  } else {
    0
  }
  log_baseline + log_biomarker + log_covariate
}

.prepare_initial_conditions <- function(
  biomarker_initial,
  parameters,
  sensitivity_type
) {
  # Basic state: [Λ(t), m(t), ṁ(t)]
  basic_state <- c(0, biomarker_initial)

  switch(
    sensitivity_type,
    basic = basic_state,
    forward = {
      # Augmented for parameters needing ODE: (η, α, β)
      # [Λ, m, ṁ, ∂Λ/∂η, ∂Λ/∂α, ∂m/∂β, ∂ṁ/∂β, ∂Λ/∂β]
      n_eta <- parameters$config$baseline$df
      n_alpha <- 2
      n_beta <- nrow(parameters$coef$acceleration)
      c(
        basic_state,
        rep(0, n_eta), # ∂Λ/∂η
        rep(0, n_alpha), # ∂Λ/∂α
        rep(0, n_beta), # ∂m/∂β
        rep(0, n_beta), # ∂ṁ/∂β
        rep(0, n_beta) # ∂Λ/∂β
      )
    },
    {
      stop("Invalid sensitivity_type: ", sensitivity_type)
    }
  )
}

.extract_ode_results <- function(
  solution,
  data,
  parameters,
  sensitivity_type,
  times = NULL
) {
  # Get final state
  final_state <- solution[nrow(solution), -1]
  event_time <- data$time

  # Compute final values
  biomarker_final <- final_state[2]
  velocity_final <- final_state[3]

  acceleration_final <- .compute_acceleration(
    biomarker_final,
    velocity_final,
    event_time,
    data,
    parameters
  )

  log_hazard_final <- .compute_log_hazard(
    event_time,
    biomarker_final,
    velocity_final,
    data,
    parameters
  )

  times_to_extract <- if (is.null(times)) {
    data$longitudinal$times
  } else {
    times
  }

  # Extract values at specified times
  if (length(times_to_extract) > 0) {
    obs_indices <- match(times_to_extract, solution[, 1])
    cumhazard_values <- solution[obs_indices, 2]
    biomarker_values <- solution[obs_indices, 3]
    velocity_values <- solution[obs_indices, 4]
    acceleration_values <- sapply(seq_along(obs_indices), function(i) {
      idx <- obs_indices[i]
      .compute_acceleration(
        solution[idx, 3],
        solution[idx, 4],
        solution[idx, 1],
        data,
        parameters
      )
    })
  } else {
    cumhazard_values <- NULL
    biomarker_values <- NULL
    velocity_values <- NULL
    acceleration_values <- NULL
  }

  # Build base result
  result <- list(
    cum_hazard_at_event = final_state[1],
    log_hazard_at_event = log_hazard_final,
    cum_hazard = cumhazard_values,
    biomarker = biomarker_values,
    velocity = velocity_values,
    acceleration = acceleration_values
  )

  # Add sensitivity-specific outputs
  if (sensitivity_type == "forward") {
    # Extract sensitivities for θ = (η, α, β)
    n_eta <- parameters$configurations$baseline$df
    n_alpha <- 2
    n_beta <- length(parameters$coefficients$acceleration)

    # Extract sensitivities from final state (matching ODE output order)
    # State order: [Λ, m, ṁ, ∂Λ/∂η, ∂Λ/∂α, ∂m/∂β, ∂ṁ/∂β, ∂Λ/∂β]
    idx <- 4
    dcumhazard_deta <- final_state[idx:(idx + n_eta - 1)]
    idx <- idx + n_eta
    dcumhazard_dalpha <- final_state[idx:(idx + n_alpha - 1)]
    idx <- idx + n_alpha
    dbiomarker_dbeta <- final_state[idx:(idx + n_beta - 1)]
    idx <- idx + n_beta
    dvelocity_dbeta <- final_state[idx:(idx + n_beta - 1)]
    idx <- idx + n_beta
    dcumhazard_dbeta <- final_state[idx:(idx + n_beta - 1)]

    # Extract β sensitivities at observation times
    dbiomarker_dbeta_values <- if (length(times_to_extract) > 0) {
      obs_indices <- match(times_to_extract, solution[, 1])
      # Column 1 is time, columns 2-4 are [Λ, m, ṁ]
      # Columns 5-(4+n_eta) are ∂Λ/∂η
      # Columns (5+n_eta)-(4+n_eta+n_alpha) are ∂Λ/∂α
      # Columns (5+n_eta+n_alpha)-(4+n_eta+n_alpha+n_beta) are ∂m/∂β
      col_start <- 5 + n_eta + n_alpha # Adjust for time column (column 1)
      col_end <- col_start + n_beta - 1
      solution[obs_indices, col_start:col_end, drop = FALSE]
    } else {
      NULL
    }

    # Compute acceleration sensitivity at event time
    dacceleration_dbeta <- .compute_acceleration_deriv(
      biomarker_final,
      velocity_final,
      event_time,
      data,
      parameters,
      dbiomarker = dbiomarker_dbeta,
      dvelocity = dvelocity_dbeta
    )

    result <- list(
      log_hazard_at_event = log_hazard_final,
      cum_hazard_at_event = final_state[1],
      biomarker_at_event = final_state[2],
      velocity_at_event = final_state[3],
      acceleration_at_event = acceleration_final,
      biomarker = biomarker_values,
      velocity = velocity_values,
      acceleration = acceleration_values,
      # Survival parameter sensitivities
      dcumhazard_deta_at_event = dcumhazard_deta,
      dcumhazard_dalpha_at_event = dcumhazard_dalpha,
      # Beta sensitivities
      dbiomarker_dbeta = dbiomarker_dbeta_values,
      dbiomarker_dbeta_at_event = dbiomarker_dbeta,
      dvelocity_dbeta_at_event = dvelocity_dbeta,
      dacceleration_dbeta_at_event = dacceleration_dbeta,
      dcumhazard_dbeta_at_event = dcumhazard_dbeta
    )
  }

  result
}

.solve_joint_ode <- function(
  data,
  parameters,
  sensitivity_type = "basic",
  times = NULL
) {
  # Define ODE derivatives function based on sensitivity type
  ode_derivatives <- function(time, state, parameters) {
    biomarker <- state[2]
    velocity <- state[3]

    # Compute acceleration
    acceleration <- .compute_acceleration(
      biomarker,
      velocity,
      time,
      parameters$data,
      parameters$parameters
    )

    # Compute hazard (without random effect b)
    log_hazard <- .compute_log_hazard(
      time,
      biomarker,
      velocity,
      parameters$data,
      parameters$parameters
    )
    # Cap log_hazard to prevent numerical overflow
    log_hazard_capped <- pmin(log_hazard, 50) # exp(50) ≈ 5e21
    hazard <- exp(log_hazard_capped)

    # Basic derivatives: [dΛ/dt, dm/dt, dṁ/dt]
    basic_derivs <- c(hazard, velocity, acceleration)

    if (parameters$sensitivity_type == "basic") {
      list(basic_derivs)
    } else if (parameters$sensitivity_type == "forward") {
      # Forward sensitivity for parameters needing ODE: (η, α, β)
      # Note: φ (survival covariates) don't need ODE integration
      n_eta <- parameters$parameters$configurations$baseline$df
      n_alpha <- 2
      n_beta <- length(parameters$parameters$coefficients$acceleration)

      # Compute direct hazard sensitivities
      basis_lambda <- .compute_spline_basis(
        time,
        parameters$parameters$configurations$baseline
      )
      m_vec <- c(biomarker, velocity)

      # Extract β sensitivity states from augmented state
      idx <- 4 + n_eta + n_alpha
      dbiomarker_dbeta <- state[idx:(idx + n_beta - 1)]
      idx <- idx + n_beta
      dvelocity_dbeta <- state[idx:(idx + n_beta - 1)]

      # Compute acceleration sensitivity (linear model)
      dacceleration_dbeta <- .compute_acceleration_deriv(
        biomarker,
        velocity,
        time,
        parameters$data,
        parameters$parameters,
        dbiomarker = dbiomarker_dbeta,
        dvelocity = dvelocity_dbeta
      )

      # Compute ∂λ/∂β using chain rule
      alpha <- parameters$parameters$coefficients$hazard[1:2]
      dm_vec_dbeta <- rbind(
        dbiomarker_dbeta,
        dvelocity_dbeta
      )
      # Fix: Use proper matrix multiplication
      dhazard_dbeta <- as.vector(
        t(alpha) %*% dm_vec_dbeta
      ) *
        hazard

      # Return augmented derivatives
      # State order: [Λ, m, ṁ, ∂Λ/∂η, ∂Λ/∂α, ∂m/∂β, ∂ṁ/∂β, ∂Λ/∂β]
      list(c(
        basic_derivs, # [dΛ/dt, dm/dt, dṁ/dt]
        as.numeric(basis_lambda) * hazard, # d(∂Λ/∂η)/dt = ∂λ/∂η
        m_vec * hazard, # d(∂Λ/∂α)/dt = ∂λ/∂α
        dvelocity_dbeta, # d(∂m/∂β)/dt = ∂ṁ/∂β
        dacceleration_dbeta, # d(∂ṁ/∂β)/dt = ∂m̈/∂β
        dhazard_dbeta # d(∂Λ/∂β)/dt = ∂λ/∂β
      ))
    } else {
      stop("Unknown sensitivity type: ", sensitivity_type)
    }
  }

  # Prepare initial conditions based on sensitivity type
  initial_extended <- .prepare_initial_conditions(
    data$initial_state,
    parameters,
    sensitivity_type
  )

  # If custom times provided, add them to ensure they're in the solution
  ode_times <- if (!is.null(times)) {
    sort(unique(c(0, times)))
  } else {
    sort(unique(c(0, data$longitudinal$times, data$time)))
  }

  n_groups <- parameters$configurations$n_groups
  ode_results <- vector("list", n_groups)
  for (k in seq_len(n_groups)) {
    parameters_subject <- list(
      coefficients = list(
        baseline = parameters$coefficients$baseline,
        hazard = parameters$coefficients$hazard,
        acceleration = parameters$coefficients$acceleration[, k],
        measurement_error_sd = parameters$coefficients$measurement_error_sd,
        random_effect_sd = parameters$coefficients$random_effect_sd
      ),
      configurations = parameters$configurations
    )
    ode_parameters <- list(
      data = data,
      parameters = parameters_subject,
      sensitivity_type = sensitivity_type
    )
    solution <- if (length(ode_times) < 2) {
      matrix(c(0, initial_extended), nrow = 1)
    } else {
      # Use deSolve to solve the ODEs
      deSolve::ode(
        y = initial_extended,
        times = ode_times,
        func = ode_derivatives,
        parms = ode_parameters,
        method = "ode45",
        atol = 1e-8,
        rtol = 1e-10
      )
    }
    ode_results[[k]] <- .extract_ode_results(
      solution,
      data,
      parameters_subject,
      sensitivity_type,
      times
    )
  }
  ode_results
}

# Helper functions for posterior computation

# Statistical Distributions

.compute_posteriors <- function(
  data_list,
  parameters,
  n_points = 7,
  parallel = TRUE,
  n_cores = 0,
  return_ode_solutions = FALSE
) {
  n_subjects <- length(data_list)
  init <- parameters$coefficients$random_effect

  # Function to compute posterior for a single subject
  compute_subject_posterior <- function(i) {
    n_groups <- parameters$configurations$n_groups
    ode_solution <- .solve_joint_ode(data_list[[i]], parameters)
    posterior <- vector("list", n_groups)
    for (k in seq_len(n_groups)) {
      posterior[[k]] <- .compute_posterior_aghq(
        ode_solution = ode_solution[[k]],
        data = data_list[[i]],
        b_hat_init = init[i, k],
        measurement_error_sd = parameters$coefficients$measurement_error_sd,
        random_effect_sd = parameters$coefficients$random_effect_sd,
        n_points = n_points
      )

      # Optionally include ODE solution for later reuse
      if (return_ode_solutions) {
        posterior[[k]]$ode_solution <- ode_solution[[k]]
      }
    }
    posterior
  }

  # Compute posteriors in parallel or sequentially
  if (parallel) {
    cleanup <- .setup_parallel_plan(n_cores)
    on.exit(cleanup(), add = TRUE)

    # Parallel computation
    posterior_results <- future.apply::future_lapply(
      seq_len(n_subjects),
      compute_subject_posterior,
      future.seed = TRUE
    )
  } else {
    # Sequential computation
    posterior_results <- lapply(seq_len(n_subjects), compute_subject_posterior)
  }

  # Aggregate results - extract matrices for each parameter
  marginal <- do.call(
    rbind,
    lapply(posterior_results, function(x) {
      vapply(x, `[[`, numeric(1), "marginal")
    })
  )
  # Compute posterior subgroup probabilities
  weighted_marginal <- sweep(
    marginal,
    2,
    parameters$coefficients$prob,
    "*"
  )
  w <- weighted_marginal / rowSums(weighted_marginal)

  posteriors <- list(
    b = do.call(
      rbind,
      lapply(posterior_results, function(x) sapply(x, `[[`, "b"))
    ),
    v = do.call(
      rbind,
      lapply(posterior_results, function(x) sapply(x, `[[`, "v"))
    ),
    e = do.call(
      rbind,
      lapply(posterior_results, function(x) sapply(x, `[[`, "e"))
    ),
    w = w
  )

  # Store ODE solutions if requested (each subject has 4 different solutions)
  if (return_ode_solutions) {
    posteriors$ode_solutions <- lapply(
      posterior_results,
      function(subject_results) {
        lapply(subject_results, function(x) {
          if (!is.null(x$ode_solution)) x$ode_solution else NULL
        })
      }
    )
  }

  posteriors
}

.compute_posterior_aghq <- function(
  ode_solution,
  data,
  b_hat_init,
  measurement_error_sd,
  random_effect_sd,
  n_points = 7
) {
  # Extract components
  cum_hazard_0 <- ode_solution$cum_hazard_at_event
  biomarker_0 <- ode_solution$biomarker

  status <- data$status
  measurements <- data$longitudinal$measurements
  n_obs <- data$longitudinal$n_obs

  # Compute residual sum
  s_i <- sum(measurements - biomarker_0)
  inv_measurement_error_sd2 <- 1 / measurement_error_sd^2
  inv_random_effect_sd2 <- 1 / random_effect_sd^2

  # Define log-posterior function (up to normalizing constant)
  logpost <- function(b) {
    -0.5 *
      b^2 *
      (n_obs * inv_measurement_error_sd2 + inv_random_effect_sd2) +
      b * (s_i * inv_measurement_error_sd2 + status) -
      exp(b) * cum_hazard_0
  }
  logpost_grad <- function(b) {
    -b *
      (n_obs * inv_measurement_error_sd2 + inv_random_effect_sd2) +
      (s_i * inv_measurement_error_sd2 + status) -
      exp(b) * cum_hazard_0
  }
  logpost_hess <- function(b) {
    -(n_obs * inv_measurement_error_sd2 + inv_random_effect_sd2) -
      exp(b) * cum_hazard_0
  }

  fit_aghq <- aghq::aghq(
    ff = list(
      fn = logpost,
      gr = logpost_grad,
      he = logpost_hess
    ),
    k = n_points,
    startingvalue = b_hat_init
  )

  b_hat <- aghq::compute_moment(
    fit_aghq$normalized_posterior,
    ff = function(x) x
  )
  b2_hat <- aghq::compute_moment(
    fit_aghq$normalized_posterior,
    ff = function(x) x^2
  )
  v_hat <- b2_hat - b_hat^2
  e <- aghq::compute_moment(
    fit_aghq$normalized_posterior,
    ff = function(x) exp(x)
  )

  lognormconst <- fit_aghq$normalized_posterior$lognormconst
  log_marginal <- pmin(lognormconst, 30)

  list(
    b = b_hat,
    v = v_hat,
    e = e,
    marginal = exp(log_marginal)
  )
}

.m_step_closed_form <- function(
  data_list,
  parameters,
  posteriors,
  ode_solutions
) {
  n_subjects <- length(data_list)
  n_groups <- parameters$configurations$n_groups
  n_observations <- sum(sapply(data_list, function(d) d$longitudinal$n_obs))

  prob <- colMeans(posteriors$w)

  # Update variance components using weighted contributions
  measurement_error_variance <- 0
  random_effect_variance <- 0

  for (i in seq_len(n_subjects)) {
    n_i <- data_list[[i]]$longitudinal$n_obs

    for (k in seq_len(n_groups)) {
      # Weight for this subject-group combination
      w_ik <- posteriors$w[i, k]

      if (n_i > 0) {
        # Longitudinal residuals
        residuals <- data_list[[i]]$longitudinal$measurements -
          ode_solutions[[i]][[k]]$biomarker -
          posteriors$b[i, k]

        measurement_error_variance <- measurement_error_variance +
          w_ik * (sum(residuals^2) + n_i * posteriors$v[i, k])
      }

      # Random effect contribution
      random_effect_variance <- random_effect_variance +
        w_ik * (posteriors$b[i, k]^2 + posteriors$v[i, k])
    }
  }

  measurement_error_sd <- sqrt(measurement_error_variance / n_observations)
  random_effect_sd <- sqrt(random_effect_variance / n_subjects)

  list(
    prob = prob, # Updated mixture probabilities
    measurement_error_sd = measurement_error_sd,
    random_effect_sd = random_effect_sd
  )
}


# Joint optimization functions for linear model

.compute_objective_joint <- function(
  params,
  data_list,
  posteriors,
  parameters,
  parallel = TRUE,
  n_cores = 0
) {
  # Compute negative expected complete-data log-likelihood
  # Input: params = c(eta, alpha, phi, beta) - parameter vector for optim
  # Returns: scalar value (negative log-likelihood) for minimization

  # Parse parameters
  n_baseline <- parameters$configurations$baseline$df
  n_hazard <- length(parameters$coefficients$hazard)
  n_acceleration <- length(parameters$coefficients$acceleration)

  idx <- 1
  parameters$coefficients$baseline <- params[idx:(idx + n_baseline - 1)]
  idx <- idx + n_baseline
  parameters$coefficients$hazard <- params[idx:(idx + n_hazard - 1)]
  idx <- idx + n_hazard
  parameters$coefficients$acceleration <- matrix(
    params[idx:(idx + n_acceleration - 1)],
    ncol = parameters$configurations$n_groups
  )

  n_subjects <- length(data_list)
  n_groups <- parameters$configurations$n_groups
  sigma_e <- parameters$coefficients$measurement_error_sd
  sigma_b <- parameters$coefficients$random_effect_sd
  inv_sigma_e2 <- 1 / (sigma_e^2)
  inv_sigma_b2 <- 1 / (sigma_b^2)

  # Function to compute objective for a single subject
  compute_subject_objective <- function(i) {
    subject_data <- data_list[[i]]

    # Solve ODE for this subject
    ode_sol <- .solve_joint_ode(
      subject_data,
      parameters,
      sensitivity_type = "basic"
    )

    # Initialize subject objective and observation count
    obj_i <- 0
    n_obs_i <- 0

    for (k in seq_len(n_groups)) {
      # Survival component
      if (subject_data$status == 1) {
        obj_i <- obj_i + posteriors$w[i, k] * ode_sol[[k]]$log_hazard_at_event
      }
      obj_i <- obj_i -
        posteriors$e[i, k] *
          posteriors$w[i, k] *
          ode_sol[[k]]$cum_hazard_at_event

      # Longitudinal component
      if (subject_data$longitudinal$n_obs > 0) {
        residuals <- subject_data$longitudinal$measurements -
          ode_sol[[k]]$biomarker -
          posteriors$b[i, k]
        obj_i <- obj_i -
          0.5 * posteriors$w[i, k] * sum(residuals^2) * inv_sigma_e2
        n_obs_i <- subject_data$longitudinal$n_obs
      }

      # Random effect component
      obj_i <- obj_i -
        0.5 *
          posteriors$w[i, k] *
          (posteriors$b[i, k]^2 + posteriors$v[i, k]) *
          inv_sigma_b2

      # Group mixture weight
      obj_i <- obj_i + posteriors$w[i, k] * log(parameters$coefficients$prob[k])
    }

    list(objective = obj_i, n_obs = n_obs_i)
  }

  # Compute objectives in parallel or sequentially
  if (parallel) {
    cleanup <- .setup_parallel_plan(n_cores)
    on.exit(cleanup(), add = TRUE)

    # Parallel computation
    obj_results <- future.apply::future_lapply(
      seq_len(n_subjects),
      compute_subject_objective,
      future.seed = TRUE
    )
  } else {
    # Sequential computation
    obj_results <- lapply(seq_len(n_subjects), compute_subject_objective)
  }

  # Aggregate results
  objective <- sum(vapply(obj_results, `[[`, numeric(1), "objective"))
  n_obs_total <- sum(vapply(obj_results, `[[`, numeric(1), "n_obs"))

  # Add constant terms
  objective <- objective - 0.5 * n_obs_total * log(2 * pi * sigma_e^2)
  objective <- objective - 0.5 * n_subjects * log(2 * pi * sigma_b^2)

  # Return negative for minimization
  -objective
}

.compute_gradient_joint <- function(
  params,
  data_list,
  posteriors,
  parameters,
  parallel = TRUE,
  n_cores = 0,
  return_individual = FALSE
) {
  n_baseline <- parameters$configurations$baseline$df
  n_hazard <- length(parameters$coefficients$hazard)
  n_acceleration <- length(parameters$coefficients$acceleration)
  n_groups <- parameters$configurations$n_groups

  idx <- 1
  parameters$coefficients$baseline <- params[idx:(idx + n_baseline - 1)]
  idx <- idx + n_baseline
  parameters$coefficients$hazard <- params[idx:(idx + n_hazard - 1)]
  idx <- idx + n_hazard
  parameters$coefficients$acceleration <- matrix(
    params[idx:(idx + n_acceleration - 1)],
    ncol = n_groups
  )

  # Extract variance parameters
  sigma_e <- parameters$coefficients$measurement_error_sd
  sigma_b <- parameters$coefficients$random_effect_sd

  # Use provided posteriors (passed as parameter)

  n_subjects <- length(data_list)
  inv_sigma_e2 <- 1 / (sigma_e^2)

  # Function to compute gradient for a single subject
  compute_subject_gradient <- function(i) {
    subject_data <- data_list[[i]]

    # Solve ODE for all groups at once
    ode_sol <- .solve_joint_ode(
      subject_data,
      parameters,
      sensitivity_type = "forward"
    )

    status_i <- subject_data$status

    # Initialize subject-specific gradients
    grad_baseline_i <- numeric(n_baseline)
    grad_hazard_i <- numeric(n_hazard)
    grad_acceleration_i <- matrix(
      0,
      nrow = nrow(parameters$coefficients$acceleration),
      ncol = n_groups
    )

    # Aggregate gradients across groups
    for (k in seq_len(n_groups)) {
      w_ik <- posteriors$w[i, k]
      b_ik <- posteriors$b[i, k]
      e_ik <- posteriors$e[i, k]

      # Survival gradients
      # η (baseline hazard)
      basis_lambda <- .compute_spline_basis(
        subject_data$time,
        parameters$configurations$baseline
      )
      grad_baseline_i <- grad_baseline_i +
        w_ik *
          (status_i *
            basis_lambda -
            e_ik * ode_sol[[k]]$dcumhazard_deta_at_event)

      # α (association parameters)
      m_vec <- c(
        ode_sol[[k]]$biomarker_at_event,
        ode_sol[[k]]$velocity_at_event
      )
      # First two components of hazard gradient (association parameters)
      grad_hazard_i[1:2] <- grad_hazard_i[1:2] +
        w_ik *
          (status_i *
            m_vec -
            e_ik * ode_sol[[k]]$dcumhazard_dalpha_at_event)

      # Remaining hazard components (survival covariates)
      if (n_hazard > 2 && !is.null(subject_data$covariates)) {
        w_vec <- as.numeric(subject_data$covariates)
        grad_hazard_i[3:n_hazard] <- grad_hazard_i[3:n_hazard] +
          w_ik *
            ((status_i -
              e_ik * ode_sol[[k]]$cum_hazard_at_event) *
              w_vec)
      }

      # β (acceleration coefficients for group k)
      # Longitudinal component
      if (subject_data$longitudinal$n_obs > 0) {
        residuals <- subject_data$longitudinal$measurements -
          ode_sol[[k]]$biomarker -
          b_ik

        if (!is.null(ode_sol[[k]]$dbiomarker_dbeta)) {
          grad_beta_long <- as.vector(
            crossprod(residuals, ode_sol[[k]]$dbiomarker_dbeta)
          )
          grad_acceleration_i[, k] <- grad_acceleration_i[, k] +
            w_ik * grad_beta_long * inv_sigma_e2
        }
      }

      # Survival component for β
      # Second term: direct contribution at event time
      if (status_i == 1 && !is.null(ode_sol[[k]]$dbiomarker_dbeta_at_event)) {
        alpha <- parameters$coefficients$hazard[1:2]
        dm_dbeta_at_event <- rbind(
          ode_sol[[k]]$dbiomarker_dbeta_at_event,
          ode_sol[[k]]$dvelocity_dbeta_at_event
        )
        grad_acceleration_i[, k] <- grad_acceleration_i[, k] +
          w_ik * as.vector(t(alpha) %*% dm_dbeta_at_event)
      }

      # Third term: integral contribution
      if (!is.null(ode_sol[[k]]$dcumhazard_dbeta_at_event)) {
        grad_acceleration_i[, k] <- grad_acceleration_i[, k] -
          w_ik * e_ik * ode_sol[[k]]$dcumhazard_dbeta_at_event
      }
    }

    # Return gradient components
    list(
      grad_baseline = grad_baseline_i,
      grad_hazard = grad_hazard_i,
      grad_acceleration = as.vector(grad_acceleration_i)
    )
  }

  # Compute gradients in parallel or sequentially
  if (parallel) {
    cleanup <- .setup_parallel_plan(n_cores)
    on.exit(cleanup(), add = TRUE)

    # Parallel computation
    grad_results <- future.apply::future_lapply(
      seq_len(n_subjects),
      compute_subject_gradient,
      future.seed = TRUE
    )
  } else {
    # Sequential computation
    grad_results <- lapply(seq_len(n_subjects), compute_subject_gradient)
  }

  # Check if individual gradients should be returned
  if (return_individual) {
    # Return matrix where each row is a subject's gradient
    grad_matrix <- do.call(
      rbind,
      lapply(grad_results, function(res) {
        -c(res$grad_baseline, res$grad_hazard, res$grad_acceleration)
      })
    )
    return(grad_matrix)
  }

  # Aggregate results (default behavior)
  grad_baseline <- Reduce(`+`, lapply(grad_results, `[[`, "grad_baseline"))
  grad_hazard <- Reduce(`+`, lapply(grad_results, `[[`, "grad_hazard"))
  grad_acceleration <- Reduce(
    `+`,
    lapply(grad_results, `[[`, "grad_acceleration")
  )

  # Combine all gradients (negative for minimization)
  grad_all <- -c(grad_baseline, grad_hazard, grad_acceleration)

  grad_all
}
