# Model Initialization

#' @noRd
.default_parameters <- function(
  dims, gamma, parsed_long, spline_baseline_config
) {
  list(
    coefficients = list(
      baseline = rep(0, dims$n_spline_basis),
      hazard = rep(0, dims$n_survival_covariates + 2),
      longitudinal = rep(0, dims$n_longitudinal_coef),
      initial_state = c(0, 0),
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
.initialize_from_marginal <- function(
  longitudinal_data, survival_data, gamma, control,
  parsed_long, parsed_surv, model_config
) {
  verbose <- control$verbose
  dims <- model_config$dims
  sbc <- model_config$spline_baseline_config

  # Start from defaults, then override with marginal estimates
  params <- .default_parameters(dims, gamma, parsed_long, sbc)

  if (verbose > 0) cli::cli_alert_info("Initializing with MarginalODE...")

  # Fit marginal model
  id <- parsed_long$grouping
  time <- parsed_surv$time_var
  marginal_fit <- MarginalODE(
    .build_formula(parsed_long$fixed_terms, response = parsed_long$response),
    longitudinal_data, time, id,
    control = MarginalODE.control(
      maxit = 200, verbose = control$verbose,
      parallel = control$parallel, n_cores = control$n_cores
    )
  )
  if (!marginal_fit$convergence$converged) {
    stop("MarginalODE failed to converge", call. = FALSE)
  }

  params$coefficients$longitudinal <- marginal_fit$parameters
  params$coefficients$measurement_error_sd <- marginal_fit$measurement_error_sd

  # Initial state: population mean from per-subject estimates
  pred <- predict(marginal_fit)
  subj_m0 <- tapply(pred$biomarker, pred[[id]], function(x) x[1])
  subj_v0 <- tapply(pred$velocity, pred[[id]], function(x) x[1])
  params$coefficients$initial_state <- c(
    mean(subj_m0, na.rm = TRUE),
    mean(subj_v0, na.rm = TRUE)
  )

  # Random effect covariance: [state(2), coef(...)]
  n_re <- dims$n_random_effects
  sigma <- diag(1e-2, n_re)

  # State block (2×2) from empirical variance of initial conditions
  state_mat <- cbind(subj_m0 - mean(subj_m0, na.rm = TRUE),
                     subj_v0 - mean(subj_v0, na.rm = TRUE))
  state_cov <- cov(state_mat, use = "complete.obs")
  state_cov[is.na(state_cov)] <- 1e-2
  diag(state_cov) <- pmax(diag(state_cov), 1e-4)
  sigma[1:2, 1:2] <- state_cov

  # Coef block from between-subject biomarker/velocity variance
  coef_vars <- c(
    var(tapply(pred$biomarker, pred[[id]], mean), na.rm = TRUE),
    var(tapply(pred$velocity, pred[[id]], mean), na.rm = TRUE)
  )
  coef_vars[is.na(coef_vars) | coef_vars < 1e-4] <- 1e-2
  for (k in 3:n_re) sigma[k, k] <- coef_vars[((k - 3) %% 2) + 1]
  params$coefficients$random_effect_sigma <- sigma

  # Hazard from time-varying Cox model
  surv_vars <- c(parsed_surv$time_var, parsed_surv$status_var)
  surv_cov_names <- if (is.null(parsed_surv$covariate_terms)) {
    character(0)
  } else {
    parsed_surv$covariate_terms
  }

  merged <- .build_counting_process(pred, survival_data, id, surv_vars)
  cox_predictors <- c("biomarker", "velocity", surv_cov_names)
  cox_fit <- survival::coxph(
    reformulate(cox_predictors, response = quote(Surv(start, stop, event))),
    data = merged[, c("start", "stop", "event", cox_predictors)]
  )

  hazard <- coef(cox_fit)[cox_predictors]
  if (any(is.na(hazard))) {
    stop("Cox model produced NA coefficients.", call. = FALSE)
  }
  params$coefficients$hazard <- hazard
  if (verbose > 0) cli::cli_alert_success("Hazard initialized")

  # Baseline hazard via Weibull fit → B-spline projection
  weibull_fit <- survival::survreg(
    as.formula(paste0(
      "survival::Surv(", surv_vars[1], ", ", surv_vars[2], ") ~ 1"
    )),
    data = survival_data, dist = "weibull"
  )
  wb_alpha <- exp(coef(weibull_fit)[1])
  wb_shape <- 1 / weibull_fit$scale
  if (verbose > 0) {
    cli::cli_alert_info(
      "Weibull baseline: shape={round(wb_shape, 3)}, scale={round(wb_alpha, 3)}"
    )
  }

  evt <- survival_data[[surv_vars[1]]][survival_data[[surv_vars[2]]] == 1]
  if (length(evt) < 2) {
    stop("Insufficient events for baseline hazard estimation", call. = FALSE)
  }
  tg <- seq(max(min(evt), 1e-6), max(survival_data[[surv_vars[1]]]),
    length.out = 100
  )
  log_h <- log(pmax(
    (wb_shape / wb_alpha) * (tg / wb_alpha)^(wb_shape - 1), 1e-10
  ))
  basis <- splines2::bSpline(
    x = tg, knots = sbc$knots, Boundary.knots = sbc$boundary_knots,
    degree = sbc$degree, intercept = TRUE
  )
  params$coefficients$baseline <- qr.coef(qr(basis), log_h)
  if (verbose > 0) cli::cli_alert_success("Baseline initialized via Weibull")

  params
}

#' Build counting process (start, stop, event) from predictions + survival
#' @noRd
.build_counting_process <- function(pred, survival_data, id, surv_vars) {
  names(pred)[names(pred) == "time"] <- "obstime"
  merged <- merge(pred, survival_data, by = id, all.x = TRUE)
  merged <- merged[order(merged[[id]], merged$obstime), ]

  merged$start <- merged$obstime
  subj_ids <- as.vector(merged[[id]])
  subj_idx <- c(which(!duplicated(subj_ids)), nrow(merged) + 1)
  ni <- diff(subj_idx)

  stop_vec <- event_vec <- numeric(nrow(merged))
  for (i in seq_along(ni)) {
    s <- subj_idx[i]
    e <- subj_idx[i + 1] - 1
    if (ni[i] == 1) {
      stop_vec[s] <- merged[[surv_vars[1]]][s]
      event_vec[s] <- merged[[surv_vars[2]]][s]
    } else {
      stop_vec[s:e] <- c(merged$obstime[(s + 1):e], merged[[surv_vars[1]]][e])
      event_vec[s:e] <- c(rep(0, ni[i] - 1), merged[[surv_vars[2]]][e])
    }
  }
  merged$stop <- stop_vec
  merged$event <- event_vec
  merged
}

# Model Setup ==================================================================

#' @importFrom stats model.frame model.matrix model.response
#' @noRd
.setup_model <- function(
  longitudinal_data, survival_data, survival_formula,
  gamma, parsed_long, parsed_surv, spline_config
) {
  # Fixed covariate names
  fixed_formula <- .build_formula(
    parsed_long$fixed_terms,
    response = parsed_long$response
  )
  long_fixed_names <- colnames(model.matrix(
    fixed_formula, model.frame(fixed_formula, longitudinal_data)
  ))

  # Random effects dimension
  has_re_covs <- !is.null(parsed_long$random_terms)
  has_re_bv <- parsed_long$biomarker$random || parsed_long$velocity$random
  random_terms <- if (!has_re_covs && !has_re_bv) {
    "(Intercept)"
  } else if (has_re_covs) {
    parsed_long$random_terms
  } else {
    character(0)
  }

  n_long_random <- if (length(random_terms) > 0) {
    ncol(model.matrix(
      .build_formula(random_terms, is_random = TRUE), longitudinal_data
    ))
  } else {
    0
  }
  n_re <- n_long_random + sum(
    parsed_long$biomarker$random, parsed_long$velocity$random
  )

  # Survival dimensions
  surv_frame <- model.frame(survival_formula, survival_data)
  event_times <- model.response(surv_frame)[, 1]
  has_surv_covs <- length(all.vars(survival_formula[[3]])) > 0 &&
    survival_formula[[3]] != 1
  surv_names <- if (has_surv_covs) {
    colnames(model.matrix(survival_formula, surv_frame)[, -1, drop = FALSE])
  } else {
    character(0)
  }

  # Spline configuration from data
  sbc <- .get_spline_config(
    x = event_times,
    degree = spline_config$degree,
    n_knots = spline_config$n_knots,
    knot_placement = spline_config$knot_placement,
    boundary_knots = spline_config$boundary_knots
  )
  sbc$boundary_knots[1] <- 0

  # Longitudinal coefficient names (biomarker/velocity first)
  long_names <- character(0)
  if (parsed_long$biomarker$fixed) long_names <- c(long_names, "biomarker")
  if (parsed_long$velocity$fixed) long_names <- c(long_names, "velocity")
  long_names <- c(long_names, long_fixed_names)

  # RE layout: [b_m0, b_v0, b_coef...] — front 2 are state RE
  n_re_total <- n_re + 2  # +2 for initial state RE
  random_effects <- matrix(0, nrow(survival_data), n_re_total)

  list(
    dims = .compute_dimensions(parsed_long, parsed_surv, spline_config),
    random_effects = random_effects,
    coef_names = list(
      baseline = paste0("bs", seq_len(sbc$df)),
      hazard = c("alpha_1", "alpha_2", surv_names),
      longitudinal = long_names,
      initial_state = c("m0", "v0")
    ),
    spline_baseline_config = sbc
  )
}
