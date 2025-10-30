#' @importFrom stats predict reformulate
#' @noRd
.compute_initial <- function(
  longitudinal_formula,
  longitudinal_data,
  survival_formula,
  survival_data,
  spline_baseline,
  gamma,
  state,
  control_opts
) {
  verbose <- control_opts$verbose

  parsed_long <- .parse_longitudinal_formula(longitudinal_formula)
  parsed_surv <- .parse_survival_formula(survival_formula)
  spline_baseline <- modifyList(.default_spline, spline_baseline)

  id <- parsed_long$grouping
  time <- parsed_surv$time_var

  dims <- .compute_dimensions(
    longitudinal_formula,
    survival_formula,
    spline_baseline
  )

  default_init <- list(
    coefficients = list(
      baseline = rep(0, dims$n_spline_basis),
      hazard = rep(0, dims$n_survival_covariates + 2),
      longitudinal = rep(0, dims$n_longitudinal_coef),
      measurement_error_sd = 1,
      random_effect_sigma = diag(
        1.0,
        dims$n_random_effects,
        dims$n_random_effects
      )
    ),
    configurations = list(
      baseline = dims$spline_config,
      gamma = gamma,
      biomarker = list(
        fixed = parsed_long$biomarker$fixed,
        random = parsed_long$biomarker$random
      ),
      velocity = list(
        fixed = parsed_long$velocity$fixed,
        random = parsed_long$velocity$random
      )
    )
  )

  if (verbose > 0) {
    cli::cli_alert_info("Initializing with MarginalODE...")
  }

  # MarginalODE only supports fixed effects, so build a fixed-only formula
  fixed_formula <- .build_formula(
    parsed_long$fixed_terms,
    response = parsed_long$response
  )

  marginal_fit <- tryCatch(
    {
      MarginalODE(
        fixed_formula,
        longitudinal_data,
        time,
        id,
        as.matrix(state),
        control_opts
      )
    },
    error = function(e) {
      if (verbose > 0) {
        cli::cli_alert_warning(
          "MarginalODE failed with error: {e$message}"
        )
        cli::cli_alert_warning("Using default initial values")
      }
      NULL
    }
  )

  # Check if MarginalODE failed or didn't converge
  if (is.null(marginal_fit) || !marginal_fit$convergence$converged) {
    if (!is.null(marginal_fit) && verbose > 0) {
      cli::cli_alert_warning("MarginalODE failed to converge, using defaults")
    }
    return(default_init)
  }

  longitudinal <- marginal_fit$parameters

  surv_vars <- all.vars(survival_formula[[2]])
  surv_cov_names <- all.vars(survival_formula[[3]])

  pred_result <- predict(marginal_fit)
  names(pred_result)[names(pred_result) == "time"] <- "obstime"

  merged_data <- merge(pred_result, survival_data, by = id, all.x = TRUE)
  merged_data <- merged_data[order(merged_data[[id]], merged_data$obstime), ]

  merged_data$start <- merged_data$obstime
  subj_ids <- as.vector(merged_data[[id]])
  subj_idx <- c(which(!duplicated(subj_ids)), nrow(merged_data) + 1)
  ni <- diff(subj_idx)

  stop_vec <- numeric(nrow(merged_data))
  event_vec <- numeric(nrow(merged_data))
  write_idx <- 1

  for (i in seq_along(ni)) {
    idx_start <- subj_idx[i]
    idx_end <- subj_idx[i + 1] - 1
    n_write <- ni[i]

    if (ni[i] == 1) {
      stop_vec[write_idx] <- merged_data[[surv_vars[1]]][idx_start]
      event_vec[write_idx] <- merged_data[[surv_vars[2]]][idx_start]
    } else {
      stop_vec[write_idx:(write_idx + n_write - 1)] <- c(
        merged_data$obstime[(idx_start + 1):idx_end],
        merged_data[[surv_vars[1]]][idx_end]
      )
      event_vec[write_idx:(write_idx + n_write - 1)] <- c(
        rep(0, ni[i] - 1),
        merged_data[[surv_vars[2]]][idx_end]
      )
    }
    write_idx <- write_idx + n_write
  }

  merged_data$stop <- stop_vec
  merged_data$event <- event_vec

  cox_predictors <- c("biomarker", "velocity", surv_cov_names)
  cox_formula <- reformulate(
    cox_predictors,
    response = quote(Surv(start, stop, event))
  )
  cox_data <- merged_data[, c("start", "stop", "event", cox_predictors)]

  cox_fit <- tryCatch(
    survival::coxph(cox_formula, data = cox_data),
    error = function(e) {
      stop(
        "Cox model fit failed: ",
        e$message,
        "\nCheck for sufficient events and covariate variation.",
        call. = FALSE
      )
    }
  )

  required_coefs <- c("biomarker", "velocity", surv_cov_names)
  missing_coefs <- setdiff(required_coefs, names(coef(cox_fit)))
  if (length(missing_coefs) > 0) {
    stop(
      "Missing coefficients: ",
      paste(missing_coefs, collapse = ", "),
      "\nCheck for collinearity or insufficient variation.",
      call. = FALSE
    )
  }

  hazard <- coef(cox_fit)[required_coefs]
  if (any(is.na(hazard))) {
    stop("Cox model produced NA coefficients.", call. = FALSE)
  }

  if (verbose > 0) {
    cli::cli_alert_success("Hazard initialized")
  }

  spline_baseline_config <- .get_spline_config(
    x = survival_data[[surv_vars[1]]],
    degree = spline_baseline$degree,
    n_knots = spline_baseline$n_knots,
    knot_placement = spline_baseline$knot_placement,
    boundary_knots = spline_baseline$boundary_knots
  )
  spline_baseline_config$boundary_knots[1] <- 0

  baseline <- .init_baseline_spline(
    survival_data,
    surv_vars[1],
    surv_vars[2],
    spline_baseline_config,
    verbose
  )

  list(
    coefficients = list(
      baseline = baseline,
      hazard = hazard,
      longitudinal = longitudinal,
      measurement_error_sd = default_init$coefficients$measurement_error_sd,
      random_effect_sigma = default_init$coefficients$random_effect_sigma
    ),
    configurations = list(
      baseline = spline_baseline_config,
      gamma = gamma,
      biomarker = list(
        fixed = parsed_long$biomarker$fixed,
        random = parsed_long$biomarker$random
      ),
      velocity = list(
        fixed = parsed_long$velocity$fixed,
        random = parsed_long$velocity$random
      )
    )
  )
}

.init_baseline_spline <- function(
  survival_data,
  time_var,
  status_var,
  spline_config,
  verbose
) {
  weibull_formula <- as.formula(paste0(
    "survival::Surv(",
    time_var,
    ", ",
    status_var,
    ") ~ 1"
  ))

  weibull_fit <- tryCatch(
    survival::survreg(weibull_formula, data = survival_data, dist = "weibull"),
    error = function(e) {
      if (verbose > 0) {
        cli::cli_alert_warning("Weibull fit failed: {e$message}")
      }
      NULL
    }
  )

  if (is.null(weibull_fit)) {
    if (verbose > 0) {
      cli::cli_alert_warning("Using constant baseline hazard")
    }
    n_events <- sum(survival_data[[status_var]])
    total_time <- sum(survival_data[[time_var]])
    const_hazard <- log(pmax(n_events / total_time, 0.01))
    return(rep(const_hazard, spline_config$df))
  }

  # Weibull parameterization: log(T) = μ + σ*ε
  # Hazard: λ_0(t) = (γ/α) * (t/α)^(γ-1), α = exp(μ), γ = 1/σ
  mu <- coef(weibull_fit)[1]
  sigma <- weibull_fit$scale
  alpha <- exp(mu)
  gamma <- 1 / sigma

  if (verbose > 0) {
    cli::cli_alert_info(
      "Weibull baseline: shape={round(gamma, 3)}, scale={round(alpha, 3)}"
    )
  }

  event_times <- survival_data[[time_var]][survival_data[[status_var]] == 1]
  if (length(event_times) < 2) {
    if (verbose > 0) {
      cli::cli_alert_warning("Insufficient events, using constant baseline")
    }
    return(rep(log(0.01), spline_config$df))
  }

  max_time <- max(survival_data[[time_var]])
  min_time <- max(min(event_times), 1e-6)
  time_grid <- seq(min_time, max_time, length.out = 100)

  weibull_hazard <- (gamma / alpha) * (time_grid / alpha)^(gamma - 1)
  log_hazard <- log(pmax(weibull_hazard, 1e-10))

  basis <- splines2::bSpline(
    x = time_grid,
    knots = spline_config$knots,
    Boundary.knots = spline_config$boundary_knots,
    degree = spline_config$degree,
    intercept = TRUE
  )
  qr_decomp <- qr(basis)

  if (qr_decomp$rank < ncol(basis)) {
    if (verbose > 0) {
      cli::cli_alert_warning("Spline basis rank-deficient, using constant")
    }
    coef <- rep(mean(log_hazard), ncol(basis))
  } else {
    coef <- qr.coef(qr_decomp, log_hazard)
    if (anyNA(coef)) {
      if (verbose > 0) {
        cli::cli_alert_warning("QR produced NA, using mean log-hazard")
      }
      coef <- rep(mean(log_hazard), ncol(basis))
    }
  }

  if (verbose > 0) {
    cli::cli_alert_success("Baseline initialized via Weibull")
  }

  coef
}


#' @importFrom stats model.frame model.matrix model.response
#' @noRd
.initialize <- function(
  longitudinal_formula,
  survival_formula,
  longitudinal_data,
  survival_data,
  spline_baseline,
  gamma,
  state,
  init,
  control
) {
  # Set defaults for spline_baseline if not provided
  spline_baseline <- modifyList(.default_spline, spline_baseline)

  # Extract data dimensions from raw data
  n_subjects <- nrow(survival_data)

  # Parse formulas once
  parsed_long <- .parse_longitudinal_formula(longitudinal_formula)
  parsed_surv <- .parse_survival_formula(survival_formula)

  # Build formulas and matrices
  fixed_formula <- .build_formula(
    parsed_long$fixed_terms,
    response = parsed_long$response
  )
  long_fixed_frame <- model.frame(fixed_formula, longitudinal_data)
  long_fixed_matrix <- model.matrix(fixed_formula, long_fixed_frame)

  has_biomarker_velocity_random <- parsed_long$biomarker$random ||
    parsed_long$velocity$random

  random_terms <- if (
    is.null(parsed_long$random_terms) &&
      !has_biomarker_velocity_random
  ) {
    "(Intercept)"
  } else if (!is.null(parsed_long$random_terms)) {
    parsed_long$random_terms
  } else {
    character(0)
  }

  if (length(random_terms) > 0) {
    random_formula <- .build_formula(random_terms, is_random = TRUE)
    long_random_matrix <- model.matrix(random_formula, longitudinal_data)
  } else {
    long_random_matrix <- matrix(0, nrow = nrow(longitudinal_data), ncol = 0)
  }

  surv_frame <- model.frame(survival_formula, survival_data)
  surv_response <- model.response(surv_frame)
  event_times <- surv_response[, 1]

  # Extract dimensions
  n_longitudinal_fixed <- ncol(long_fixed_matrix)
  n_longitudinal_random <- ncol(long_random_matrix)

  has_surv_covs <- length(all.vars(survival_formula[[3]])) > 0 &&
    survival_formula[[3]] != 1
  surv_matrix <- if (has_surv_covs) {
    model.matrix(survival_formula, surv_frame)[, -1, drop = FALSE]
  } else {
    NULL
  }

  # Calculate biomarker/velocity fixed effects count
  n_biomarker_velocity_fixed <- sum(
    parsed_long$biomarker$fixed,
    parsed_long$velocity$fixed
  )

  # Configure splines
  spline_baseline_config <- .get_spline_config(
    x = event_times,
    degree = spline_baseline$degree,
    n_knots = spline_baseline$n_knots,
    knot_placement = spline_baseline$knot_placement,
    boundary_knots = spline_baseline$boundary_knots
  )
  spline_baseline_config$boundary_knots[1] <- 0

  # Extract variable names
  long_fixed_vars_names <- colnames(long_fixed_matrix)
  surv_vars_names <- if (!is.null(surv_matrix)) {
    colnames(surv_matrix)
  } else {
    character(0)
  }

  # Build longitudinal coefficient names
  # (biomarker/velocity first, then covariates)
  accel_names <- character(0)
  if (parsed_long$biomarker$fixed) {
    accel_names <- c(accel_names, "biomarker")
  }
  if (parsed_long$velocity$fixed) {
    accel_names <- c(accel_names, "velocity")
  }
  accel_names <- c(accel_names, long_fixed_vars_names)

  # Create coefficient names
  coef_names <- list(
    baseline = paste0("bs", seq_len(spline_baseline_config$df)),
    hazard = c("alpha_1", "alpha_2", surv_vars_names),
    longitudinal = accel_names
  )

  # Calculate random effects dimension for initialization
  n_biomarker_velocity_random <- sum(
    parsed_long$biomarker$random,
    parsed_long$velocity$random
  )
  n_random_effects <- n_longitudinal_random + n_biomarker_velocity_random
  random_effects <- matrix(0, n_subjects, n_random_effects)

  if (is.null(init)) {
    parameters <- .compute_initial(
      longitudinal_formula,
      longitudinal_data,
      survival_formula,
      survival_data,
      spline_baseline,
      gamma,
      state,
      control
    )
  } else {
    parameters <- init
  }

  # Set correct names after initialization
  names(parameters$coefficients$baseline) <- coef_names$baseline
  names(parameters$coefficients$hazard) <- coef_names$hazard
  names(parameters$coefficients$longitudinal) <- coef_names$longitudinal

  list(
    parameters = parameters,
    random_effects = random_effects,
    coef_names = coef_names
  )
}
