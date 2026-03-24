#' @noRd
.default_parameters <- function(
  dims, gamma, parsed_long, spline_baseline_config
) {
  list(
    coefficients = list(
      baseline = rep(0, dims$n_spline_basis),
      hazard = rep(0, dims$n_survival_covariates + 2),
      longitudinal = rep(0, dims$n_longitudinal_coef),
      measurement_error_sd = 1,
      random_effect_sigma = diag(1, dims$n_random_effects)
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

#' @importFrom stats predict reformulate var
#' @noRd
.compute_initial <- function(
  longitudinal_data,
  survival_data,
  gamma,
  state,
  control_opts,
  parsed_long,
  parsed_surv,
  spline_config,
  spline_baseline_config
) {
  verbose <- control_opts$verbose

  id <- parsed_long$grouping
  time <- parsed_surv$time_var

  dims <- .compute_dimensions(parsed_long, parsed_surv, spline_config)
  default_init <- .default_parameters(
    dims, gamma, parsed_long, spline_baseline_config
  )

  if (verbose > 0) {
    cli::cli_alert_info("Initializing with MarginalODE...")
  }

  # MarginalODE only supports fixed effects, so build a fixed-only formula
  fixed_formula <- .build_formula(
    parsed_long$fixed_terms,
    response = parsed_long$response
  )

  marginal_fit <- MarginalODE(
    fixed_formula, longitudinal_data, time, id, state, control_opts
  )
  if (!marginal_fit$convergence$converged) {
    stop("MarginalODE failed to converge", call. = FALSE)
  }

  longitudinal <- marginal_fit$parameters
  measurement_error_sd <- marginal_fit$measurement_error_sd

  # Estimate random effect covariance from individual-level predictions
  pred_result <- predict(marginal_fit)
  subj_means <- tapply(pred_result$biomarker, pred_result[[id]], mean)
  subj_vel <- tapply(pred_result$velocity, pred_result[[id]], mean)
  re_vars <- c(var(subj_means, na.rm = TRUE), var(subj_vel, na.rm = TRUE))
  re_vars[is.na(re_vars) | re_vars < 1e-4] <- 1e-2
  n_re <- dims$n_random_effects
  random_effect_sigma <- diag(
    rep_len(re_vars, n_re),
    n_re, n_re
  )

  surv_vars <- c(parsed_surv$time_var, parsed_surv$status_var)
  surv_cov_names <- if (is.null(parsed_surv$covariate_terms)) {
    character(0)
  } else {
    parsed_surv$covariate_terms
  }

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

  cox_fit <- survival::coxph(cox_formula, data = cox_data)

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
      measurement_error_sd = measurement_error_sd,
      random_effect_sigma = random_effect_sigma
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

  weibull_fit <- survival::survreg(
    weibull_formula, data = survival_data, dist = "weibull"
  )

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
    stop("Insufficient events for baseline hazard estimation", call. = FALSE)
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
  coef <- qr.coef(qr(basis), log_hazard)

  if (verbose > 0) {
    cli::cli_alert_success("Baseline initialized via Weibull")
  }

  coef
}


#' @importFrom stats model.frame model.matrix model.response
#' @noRd
.initialize <- function(
  longitudinal_data,
  survival_data,
  survival_formula,
  gamma,
  state,
  init,
  control,
  parsed_long,
  parsed_surv,
  spline_config
) {
  n_subjects <- nrow(survival_data)

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
  n_longitudinal_random <- ncol(long_random_matrix)

  has_surv_covs <- length(all.vars(survival_formula[[3]])) > 0 &&
    survival_formula[[3]] != 1
  surv_matrix <- if (has_surv_covs) {
    model.matrix(survival_formula, surv_frame)[, -1, drop = FALSE]
  } else {
    NULL
  }

  # Configure splines
  spline_baseline_config <- .get_spline_config(
    x = event_times,
    degree = spline_config$degree,
    n_knots = spline_config$n_knots,
    knot_placement = spline_config$knot_placement,
    boundary_knots = spline_config$boundary_knots
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

  if (is.character(init) && init == "marginal") {
    parameters <- .compute_initial(
      longitudinal_data,
      survival_data,
      gamma,
      state,
      control,
      parsed_long,
      parsed_surv,
      spline_config,
      spline_baseline_config
    )
  } else if (is.list(init)) {
    parameters <- init
  } else {
    dims <- .compute_dimensions(parsed_long, parsed_surv, spline_config)
    parameters <- .default_parameters(
      dims, gamma, parsed_long, spline_baseline_config
    )
  }

  # Ensure baseline config always uses the computed spline configuration
  # (.compute_initial() fallback may use raw config without knots/boundary)
  parameters$configurations$baseline <- spline_baseline_config

  # Data-adaptive biomarker clamp: 5x max observed value (standardized)
  response <- model.response(long_fixed_frame)
  parameters$configurations$biomarker_clamp <- max(abs(response)) * 5

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
