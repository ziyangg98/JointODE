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
      initial_state = NULL,
      measurement_error_sd = NULL,
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

#' @importFrom stats cov reformulate var
#' @noRd
.initialize_from_marginal <- function(
  longitudinal_formula, longitudinal_data, survival_data,
  gamma, control, parsed_long, parsed_surv, model_config
) {
  verbose <- control$verbose
  dims <- model_config$dims
  sbc <- model_config$spline_baseline_config
  id <- parsed_long$grouping
  time <- parsed_surv$time_var
  surv_vars <- c(parsed_surv$time_var, parsed_surv$status_var)

  params <- .default_parameters(dims, gamma, parsed_long, sbc)

  if (verbose > 0) cli::cli_alert_info("Initializing with MarginalODE...")

  marginal_fit <- MarginalODE(
    longitudinal_formula, longitudinal_data, time, id,
    control = MarginalODE.control(
      maxit = 200, verbose = max(0, verbose - 1),
      parallel = control$parallel, n_cores = control$n_cores
    )
  )
  if (!marginal_fit$convergence$converged) {
    stop("MarginalODE failed to converge", call. = FALSE)
  }

  # Transfer longitudinal coefficients
  mpar <- marginal_fit$parameters
  params$coefficients$longitudinal <- mpar[
    !names(mpar) %in% .init_state_names
  ]
  params$coefficients$initial_state <- unname(mpar[.init_state_names])
  params$coefficients$measurement_error_sd <-
    marginal_fit$measurement_error_sd

  # Random effect covariance from marginal RE posterior
  n_re <- dims$n_random_effects
  sigma <- diag(1e-2, n_re)
  marginal_re <- marginal_fit$random_effects
  if (!is.null(marginal_re) && nrow(marginal_re) > 1) {
    n_shared <- min(ncol(marginal_re), n_re)
    re_cov <- cov(marginal_re[, seq_len(n_shared), drop = FALSE],
                  use = "complete.obs")
    re_cov[is.na(re_cov)] <- 1e-2
    diag(re_cov) <- pmax(diag(re_cov), 1e-4)
    sigma[seq_len(n_shared), seq_len(n_shared)] <- re_cov
  }
  params$coefficients$random_effect_sigma <- sigma
  params$random_effects_init <- marginal_re

  # --- Survival: time-dependent Cox with predicted trajectories ---
  surv_cov_names <- parsed_surv$covariate_terms %||% character(0)

  # Get time-varying biomarker/velocity from marginal fit
  pred <- predict(marginal_fit)

  # Build counting process dataset
  cp <- .build_counting_process(pred, survival_data, id, surv_vars)

  # Merge survival covariates
  cox_predictors <- c("biomarker", "velocity", surv_cov_names)
  cox_formula <- reformulate(
    cox_predictors,
    response = call(
      "Surv", as.name("start"), as.name("stop"), as.name("event")
    )
  )
  cox_fit <- survival::coxph(cox_formula, data = cp)

  hazard <- coef(cox_fit)[cox_predictors]
  if (any(is.na(hazard))) {
    stop("Cox model produced NA coefficients.", call. = FALSE)
  }
  params$coefficients$hazard <- hazard
  # Baseline hazard via Weibull → B-spline projection
  weibull_fit <- survival::survreg(
    reformulate("1", response = call(
      "Surv", as.name(surv_vars[1]), as.name(surv_vars[2])
    )),
    data = survival_data, dist = "weibull"
  )
  wb_scale <- exp(coef(weibull_fit)[1])
  wb_shape <- 1 / weibull_fit$scale

  evt <- survival_data[[surv_vars[1]]][
    survival_data[[surv_vars[2]]] == 1
  ]
  if (length(evt) < 2) {
    stop("Insufficient events for baseline hazard", call. = FALSE)
  }
  tg <- seq(
    max(min(evt), 1e-6),
    max(survival_data[[surv_vars[1]]]),
    length.out = 100
  )
  log_h <- log(pmax(
    (wb_shape / wb_scale) * (tg / wb_scale)^(wb_shape - 1),
    1e-10
  ))
  basis <- splines2::bSpline(
    x = tg, knots = sbc$knots,
    Boundary.knots = sbc$boundary_knots,
    degree = sbc$degree, intercept = TRUE
  )
  params$coefficients$baseline <- qr.coef(qr(basis), log_h)
  if (verbose > 0) {
    cli::cli_alert_success("Initialization complete")
    cf <- params$coefficients
    vfmt <- function(x) {
      paste0("[", paste(round(x, 3), collapse = ", "), "]")
    }
    cli::cli_h3("Initialized parameters")
    cat("Longitudinal:\n")
    cat("  Coefficients:", vfmt(cf$longitudinal), "\n")
    cat("  Initial state:",
        vfmt(cf$initial_state), "\n")
    cat("Survival:\n")
    cat("  Hazard:      ", vfmt(cf$hazard), "\n")
    cat("  Baseline:    ", vfmt(cf$baseline), "\n")
    cat("Variance:\n")
    cat("  sigma_e:     ",
        format(cf$measurement_error_sd, digits = 4), "\n")
    cat("  Random SD:   ",
        vfmt(sqrt(diag(cf$random_effect_sigma))), "\n")
  }

  params
}

#' Warm-start for MarginalODE with dynamics RE
#'
#' Fit a reduced model (no dynamics RE) via MarginalODE and return
#' parameters suitable as starting values for the full model.
#' @noRd
.warm_start_marginal <- function(
  formula, data, time, id, parsed_long, control
) {
  # Build reduced formula: keep fixed effects, drop dynamics RE
  fixed_terms <- c(
    if (parsed_long$biomarker$fixed) "biomarker",
    if (parsed_long$velocity$fixed) "velocity",
    parsed_long$fixed_terms
  )
  reduced_formula <- .build_formula(
    fixed_terms, response = parsed_long$response
  )

  if (control$verbose > 0)
    cli::cli_alert_info(
      "Phase 1: fitting reduced model (no dynamics RE)"
    )

  fit0 <- MarginalODE(
    formula = reduced_formula, data = data,
    time = time, id = id, control = control
  )

  if (control$verbose > 0)
    cli::cli_alert_success(sprintf(
      "Phase 1 done (logLik: %.2f)", fit0$logLik
    ))

  list(
    longitudinal = unname(fit0$parameters[
      !names(fit0$parameters) %in% .init_state_names
    ]),
    initial_state = unname(fit0$parameters[.init_state_names]),
    measurement_error_sd = fit0$measurement_error_sd,
    random_effects = fit0$random_effects
  )
}
