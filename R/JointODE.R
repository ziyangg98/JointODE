#' Joint Modeling of Longitudinal and Survival Data Using ODEs
#'
#' @description
#' Implements a unified framework for jointly modeling longitudinal biomarker
#' trajectories and time-to-event outcomes using ordinary differential
#' equations (ODEs). The model captures complex non-linear dynamics in
#' biomarker evolution while simultaneously quantifying their association with
#' survival risk through
#' shared random effects and flexible hazard specifications.
#'
#' @param longitudinal_formula A formula specifying the longitudinal submodel.
#'   The left-hand side defines the response variable, while the right-hand
#'   side specifies fixed effects including time-varying and baseline
#'   covariates
#'   (e.g., \code{biomarker ~ time + treatment + age}).
#' @param survival_formula A formula for the survival submodel using
#'   \code{Surv(time, status)} notation on the left-hand side. The right-hand
#'   side specifies baseline hazard covariates
#'   (e.g., \code{Surv(event_time, event) ~ treatment + age}).
#' @param longitudinal_data A data frame containing repeated measurements with
#'   one row per observation. Required columns include subject identifier,
#'   measurement times, response values, and any covariates specified in the
#'   formula.
#' @param survival_data A data frame with time-to-event information containing
#'   one row per subject. Must include event/censoring times, event indicators,
#'   and baseline covariates.
#' @param gamma Numeric scalar specifying the power parameter for the velocity
#'   effect in the hazard function. When \code{gamma = 0}, velocity has no
#'   effect; \code{gamma = 1} uses linear velocity; \code{gamma = 2} uses
#'   squared velocity. Default is 1 (default: 1).
#' @param spline_baseline A list controlling the B-spline representation of the
#'   baseline hazard function with the following components:
#'   \describe{
#'     \item{\code{degree}}{Polynomial degree of the B-spline basis functions
#'       (default: 2, quadratic splines)}
#'     \item{\code{n_knots}}{Number of interior knots for flexibility
#'       (default: 1)}
#'     \item{\code{knot_placement}}{Strategy for positioning knots:
#'       \code{"quantile"} places knots at quantiles of observed event times,
#'       \code{"equal"} uses equally-spaced knots (default: \code{"equal"})}
#'     \item{\code{boundary_knots}}{A numeric vector of length 2 specifying
#'       the boundary knot locations. If \code{NULL}, automatically set to the
#'       range of observed event times (default: \code{NULL})}
#'   }
#' @param init Initial values for model parameters. Can be:
#'   \itemize{
#'     \item \code{"default"} (default): Use zero/default initial values.
#'     \item \code{"marginal"}: Use \code{\link{MarginalODE}} to compute
#'       data-driven initial estimates.
#'     \item A list with the same structure as the fitted model's
#'       \code{parameters} component for full manual control.
#'   }
#' @param control A list of control parameters for optimization, or output from
#'   \code{\link{JointODE.control}}.
#' @param ... Additional arguments passed to internal optimization routines.
#'
#' @return An S3 object of class \code{"JointODE"} containing fitted model
#'   results.
#'
#' @importFrom utils modifyList
#' @importFrom survival Surv
#' @importFrom cli cli_h2 cli_text cli_alert_success
#' @importFrom cli cli_alert_warning cli_alert_info
#'
#' @examples
#' \dontrun{
#' fit <- JointODE(
#'   longitudinal_formula = observed ~ x1 + x2,
#'   survival_formula = Surv(time, status) ~ w1 + w2,
#'   longitudinal_data = sim$data$longitudinal_data,
#'   survival_data = sim$data$survival_data
#' )
#' summary(fit)
#' }
#'
#' @concept model-fitting
#' @export
# nolint next: object_name_linter
JointODE <- function(
  longitudinal_formula,
  survival_formula,
  longitudinal_data,
  survival_data,
  gamma = 1,
  spline_baseline = list(
    degree = 2,
    n_knots = 1,
    knot_placement = "equal",
    boundary_knots = NULL
  ),
  init = "default",
  control = list(),
  ...
) {
  cl <- match.call()

  # Parse formulas and config once
  parsed_long <- .parse_longitudinal_formula(longitudinal_formula)
  parsed_surv <- .parse_survival_formula(survival_formula)
  spline_config <- modifyList(.default_spline, spline_baseline)

  # Process control settings
  if (is.null(control)) {
    control <- JointODE.control()
  } else if (is.list(control)) {
    control <- JointODE.control(.list = control)
  } else {
    stop("control must be a list or NULL")
  }

  # Validate
  .validate_joint(
    longitudinal_formula = longitudinal_formula,
    survival_formula = survival_formula,
    longitudinal_data = longitudinal_data,
    survival_data = survival_data,
    gamma = gamma,
    spline_baseline = spline_baseline,
    init = init,
    parsed_long = parsed_long,
    parsed_surv = parsed_surv,
    spline_config = spline_config
  )

  # Process data
  data_list <- .process_joint(
    longitudinal_data = longitudinal_data,
    survival_formula = survival_formula,
    survival_data = survival_data,
    parsed_long = parsed_long,
    parsed_surv = parsed_surv
  )

  # Initialize parameters
  model_config <- .setup_model(
    longitudinal_data = longitudinal_data,
    survival_data = survival_data,
    survival_formula = survival_formula,
    gamma = gamma,
    parsed_long = parsed_long,
    parsed_surv = parsed_surv,
    spline_config = spline_config
  )

  if (is.character(init) && init == "marginal") {
    stop("init = 'marginal' is not yet available in the TMB version. ",
         "Use init = 'default' or provide a parameter list.", call. = FALSE)
  } else if (is.list(init)) {
    parameters <- init
  } else {
    parameters <- .default_parameters(
      model_config$dims, gamma, parsed_long,
      model_config$spline_baseline_config
    )
  }

  parameters$configurations$baseline <- model_config$spline_baseline_config
  parameters$configurations$biomarker_clamp <- max(abs(unlist(
    lapply(data_list, function(d) d$longitudinal$measurements)
  ))) * 5
  parameters$configurations$hazard_quadrature <- control$hazard_quadrature
  parameters$configurations$gamma <- gamma

  names(parameters$coefficients$baseline) <- model_config$coef_names$baseline
  names(parameters$coefficients$hazard) <- model_config$coef_names$hazard
  names(parameters$coefficients$longitudinal) <-
    model_config$coef_names$longitudinal
  names(parameters$coefficients$initial_state) <-
    model_config$coef_names$initial_state

  coef_names <- model_config$coef_names

  # Initialize random effects
  random_effects <- model_config$random_effects
  for (i in seq_along(data_list)) {
    obs <- data_list[[i]]$longitudinal
    if (length(obs$measurements) >= 1) {
      random_effects[i, 1] <- obs$measurements[1] -
        parameters$coefficients$initial_state[1]
    }
    if (length(obs$measurements) >= 2) {
      dt <- obs$times[2] - obs$times[1]
      if (dt > 0) {
        random_effects[i, 2] <- (obs$measurements[2] - obs$measurements[1]) /
          dt - parameters$coefficients$initial_state[2]
      }
    }
  }
  parameters$random_effects_init <- random_effects

  if (control$verbose > 0) {
    cli::cli_h2("Joint ODE Model Estimation (TMB)")
  }

  # Enable OpenMP if requested
  if (control$parallel && control$n_cores > 0) {
    TMB::openmp(control$n_cores)
  }

  # Build TMB data and parameter lists
  tmb_data <- .build_tmb_data(data_list, parameters, control)
  tmb_params <- .build_tmb_parameters(parameters)

  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = "b",
    DLL = "JointODE",
    silent = control$verbose < 2
  )

  opt <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(
      iter.max = control$maxit,
      eval.max = control$maxit * 10,
      rel.tol = control$tol
    )
  )

  # Extract results
  results <- .extract_tmb_results(obj, opt, parameters, coef_names,
                                  data_list, tmb_data$n_random_effects, control)

  structure(c(results, list(
    data = data_list, control = control, call = cl,
    tmb_obj = obj, tmb_opt = opt
  )), class = "JointODE")
}

#' Summary of JointODE Model
#'
#' @param object A JointODE object
#' @param ... Additional arguments
#' @return A summary.JointODE object with coefficients and test statistics
#'
#' @concept model-summary
#' @importFrom stats coef pnorm qnorm
#' @export
summary.JointODE <- function(object, ...) {
  coefs <- coef(object)
  se <- if (!is.null(object$vcov)) {
    sqrt(diag(object$vcov))
  } else {
    rep(NA_real_, length(coefs))
  }

  # Split coefficients by component
  n_baseline <- length(object$parameters$coefficients$baseline)
  n_hazard <- length(object$parameters$coefficients$hazard)
  n_longitudinal <- length(object$parameters$coefficients$longitudinal)
  n_init <- length(object$parameters$coefficients$initial_state)

  idx_baseline <- seq_len(n_baseline)
  idx_survival <- n_baseline + seq_len(n_hazard)
  idx_longitudinal <- (n_baseline + n_hazard) + seq_len(n_longitudinal)
  idx_initial <- (n_baseline + n_hazard + n_longitudinal) + seq_len(n_init)

  coef_baseline <- .coef_table(coefs[idx_baseline], se[idx_baseline])
  rownames(coef_baseline) <- gsub("baseline:", "", rownames(coef_baseline))

  coef_survival <- .coef_table(coefs[idx_survival], se[idx_survival])
  rownames(coef_survival) <- gsub("hazard:", "", rownames(coef_survival))

  coef_longitudinal <- .coef_table(
    coefs[idx_longitudinal], se[idx_longitudinal]
  )
  rownames(coef_longitudinal) <- gsub(
    "longitudinal:", "", rownames(coef_longitudinal)
  )

  coef_initial <- .coef_table(coefs[idx_initial], se[idx_initial])
  rownames(coef_initial) <- gsub(
    "initial state:", "", rownames(coef_initial)
  )

  # Delta method for derived parameters (period and xi)
  derived_params <- NULL
  if (!is.null(object$vcov) && n_longitudinal >= 2) {
    value_coef <- coefs[idx_longitudinal[1]]
    slope_coef <- coefs[idx_longitudinal[2]]
    var_value <- object$vcov[idx_longitudinal[1], idx_longitudinal[1]]
    var_slope <- object$vcov[idx_longitudinal[2], idx_longitudinal[2]]
    cov_value_slope <- object$vcov[idx_longitudinal[1], idx_longitudinal[2]]

    if (value_coef < 0) {
      omega_est <- sqrt(-value_coef)
      period_est <- 2 * pi / omega_est
      xi_est <- -slope_coef / (2 * omega_est)

      grad_period_value <- pi / ((-value_coef)^(3 / 2))
      var_period <- grad_period_value^2 * var_value
      se_period <- sqrt(var_period)

      grad_xi_value <- -slope_coef / (4 * (-value_coef)^(3 / 2))
      grad_xi_slope <- -1 / (2 * sqrt(-value_coef))
      var_xi <- grad_xi_value^2 * var_value +
        grad_xi_slope^2 * var_slope +
        2 * grad_xi_value * grad_xi_slope * cov_value_slope
      se_xi <- sqrt(var_xi)

      derived_params <- cbind(
        Estimate = c(period_est, xi_est),
        `Std. Error` = c(se_period, se_xi),
        `z value` = c(period_est / se_period, xi_est / se_xi),
        `Pr(>|z|)` = 2 * pnorm(-abs(c(
          period_est / se_period, xi_est / se_xi
        )))
      )
      rownames(derived_params) <- c("T (period)", "xi (damping ratio)")
    }
  }

  n_subjects <- length(object$data)
  n_observations <- .n_obs(object$data)
  n_events <- sum(vapply(object$data, `[[`, numeric(1), "status"))
  event_rate <- n_events / n_subjects

  structure(
    list(
      call = object$call,
      coef_baseline = coef_baseline,
      coef_survival = coef_survival,
      coef_longitudinal = coef_longitudinal,
      coef_initial = coef_initial,
      derived_params = derived_params,
      sigma = c(sigma_e = object$parameters$coefficients$measurement_error_sd),
      sigma_b_matrix = object$parameters$coefficients$random_effect_sigma,
      logLik = object$logLik,
      AIC = object$AIC,
      BIC = object$BIC,
      cindex = object$cindex,
      nobs = n_subjects,
      n_observations = n_observations,
      n_events = n_events,
      event_rate = event_rate,
      convergence = object$convergence
    ),
    class = "summary.JointODE"
  )
}

#' @concept model-summary
#' @importFrom stats printCoefmat
#' @export
print.summary.JointODE <- function(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nData Descriptives:\n")
  cat("Longitudinal Process            Survival Process\n")
  cat(sprintf(
    "Number of Observations: %-7d Number of Events: %d (%.0f%%)\n",
    x$n_observations, x$n_events, x$event_rate * 100
  ))
  cat(sprintf("Number of Subjects: %d\n", x$nobs))

  cat(sprintf("\n%10s %10s %10s\n", "AIC", "BIC", "logLik"))
  cat(sprintf("%10.3f %10.3f %10.3f\n", x$AIC, x$BIC, x$logLik))

  cat("\nCoefficients:\n")

  cat("Longitudinal Process: Second-Order ODE Model\n")
  if (!is.null(x$coef_longitudinal)) {
    .print_coefmat(x$coef_longitudinal, digits = digits,
                   signif.stars = signif.stars, ...)
  }

  if (!is.null(x$derived_params)) {
    cat("\nODE System Characteristics:\n")
    .print_coefmat(x$derived_params, digits = digits,
                   signif.stars = signif.stars, ...)
  }

  cat("\nSurvival Process: Proportional Hazards Model\n")
  if (!is.null(x$coef_survival)) {
    .print_coefmat(x$coef_survival, digits = digits,
                   signif.stars = signif.stars, ...)
  }

  if (!is.null(x$coef_baseline)) {
    cat("\nBaseline Hazard: B-spline with", nrow(x$coef_baseline),
        "basis functions\n")
    cat("(Coefficients range:",
        sprintf("[%.3f, %.3f]",
                min(x$coef_baseline[, "Estimate"]),
                max(x$coef_baseline[, "Estimate"])), ")\n")
  }

  if (!is.null(x$coef_initial)) {
    cat("\nInitial State: Population Mean\n")
    .print_coefmat(x$coef_initial, digits = digits,
                   signif.stars = signif.stars, ...)
  }

  cat("\nVariance Components:\n")
  cat(sprintf("Measurement Error SD: %.6f\n", x$sigma["sigma_e"]))
  if (!is.null(x$sigma_b_matrix)) {
    cat("Random Effect Covariance Matrix:\n")
    print(x$sigma_b_matrix, digits = 4)
  }

  cat("\nModel Diagnostics:\n")
  if (!is.na(x$cindex)) {
    cat(sprintf("C-index (Concordance): %.3f\n", x$cindex))
  }
  cat(sprintf("Convergence: %s\n", x$convergence$message))

  invisible(x)
}

#' Extract Model Coefficients
#'
#' @param object A JointODE object
#' @param ... Additional arguments
#' @return Named numeric vector of fixed effects coefficients
#' @concept model-inspection
#' @export
coef.JointODE <- function(object, ...) {
  cf <- object$parameters$coefficients
  coefs <- c(cf$baseline, cf$hazard, cf$longitudinal, cf$initial_state)
  names(coefs) <- .prefixed_coef_names(lapply(
    list(baseline = cf$baseline, hazard = cf$hazard,
         longitudinal = cf$longitudinal, initial_state = cf$initial_state),
    names
  ))
  coefs
}

#' Extract Variance-Covariance Matrix
#'
#' @param object A JointODE object
#' @param ... Additional arguments
#' @return Variance-covariance matrix of fixed effects
#' @concept model-inspection
#' @export
vcov.JointODE <- function(object, ...) {
  object$vcov
}

#' Extract Log-Likelihood
#'
#' @param object A JointODE object
#' @param ... Additional arguments
#' @return Log-likelihood with df and nobs attributes
#' @concept model-inspection
#' @importFrom stats coef
#' @export
logLik.JointODE <- function(object, ...) {
  structure(
    object$logLik,
    df = .count_params(object$parameters),
    nobs = length(object$data),
    class = "logLik"
  )
}

#' Print JointODE Model
#'
#' @param x A JointODE object
#' @param digits Number of digits for numeric output
#' @param ... Additional arguments
#' @return Invisibly returns the object
#' @concept model-display
#' @importFrom stats coef
#' @export
print.JointODE <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nJoint ODE Model (TMB)\n")
  cat("Call: ")
  print(x$call)

  n_params <- .count_params(x$parameters)
  cat(
    "\nLog-likelihood:",
    format(x$logLik, digits = digits),
    "on", n_params, "degrees of freedom\n"
  )
  cat("AIC:", format(x$AIC, digits = digits),
      "  BIC:", format(x$BIC, digits = digits), "\n")
  invisible(x)
}

#' Predict Method for JointODE Objects
#'
#' @description
#' Computes predicted biomarker trajectories, velocities, and accelerations
#' for subjects based on the fitted joint ODE model. Predictions are obtained
#' from the TMB REPORT output.
#'
#' @param object An object of class \code{JointODE}
#' @param newdata Optional data frame with new subjects. If NULL, uses the
#'   training data from the model fit.
#' @param times Optional time points for prediction (currently unused).
#' @param ... Additional arguments (currently unused)
#'
#' @return A data.frame with columns id, time, biomarker, velocity.
#'
#' @concept model-prediction
#' @export
predict.JointODE <- function(object, newdata = NULL, times = NULL, ...) {
  if (!is.null(newdata)) stop("newdata not yet supported")

  data_list <- object$data
  reported <- object$tmb_report
  cf <- object$parameters$coefficients

  # ODE coefficients for acceleration computation
  configs <- object$parameters$configurations
  fi <- 0L
  b1_pop <- if (configs$biomarker$fixed) cf$longitudinal[fi <- fi + 1L] else 0
  b2_pop <- if (configs$velocity$fixed) cf$longitudinal[fi <- fi + 1L] else 0

  obs_offset <- 0L
  results <- vector("list", length(data_list))
  for (i in seq_along(data_list)) {
    ni <- length(data_list[[i]]$longitudinal$measurements)
    idx <- obs_offset + seq_len(ni)
    m <- as.numeric(reported$fitted_biomarker[idx])
    v <- as.numeric(reported$fitted_velocity[idx])
    cum_h <- as.numeric(reported$cumulative_hazard[i])
    log_h <- as.numeric(reported$log_hazard_at_event[i])

    # Per-subject b1, b2 (population + random effects)
    re <- object$random_effects[i, ]
    ri <- 2L
    b1 <- b1_pop
    b2 <- b2_pop
    if (configs$biomarker$random) b1 <- b1 + re[ri <- ri + 1L]
    if (configs$velocity$random) b2 <- b2 + re[ri <- ri + 1L]

    results[[i]] <- data.frame(
      id = names(data_list)[i],
      time = data_list[[i]]$longitudinal$times,
      biomarker = m,
      velocity = v,
      acceleration = b1 * m + b2 * v,
      cumhaz = cum_h,
      survival = exp(-cum_h),
      log_hazard = log_h,
      stringsAsFactors = FALSE
    )
    obs_offset <- obs_offset + ni
  }

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  result_df
}
