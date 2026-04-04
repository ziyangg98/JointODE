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
#'   When providing a list, it should have elements:
#'   \describe{
#'     \item{\code{coefficients}}{A list containing:
#'       \itemize{
#'         \item \code{baseline}: Vector of B-spline coefficients for baseline
#'           hazard (length = number of spline basis functions)
#'         \item \code{hazard}: Vector of hazard parameters including
#'           association parameters (2) and survival covariates
#'         \item \code{longitudinal}: Vector of longitudinal fixed effects
#'           including intercept and covariates
#'         \item \code{measurement_error_sd}: Residual standard deviation
#'           (positive scalar)
#'         \item \code{random_effect_sigma}: Random effect covariance matrix
#'           (positive definite matrix)
#'       }}
#'     \item{\code{configurations}}{Optional; if not provided, will use
#'       spline configuration from \code{spline_baseline}}
#'   }
#' @param control A list of control parameters for optimization, or output from
#'   \code{\link{JointODE.control}}. Key parameters include:
#'   \describe{
#'     \item{\code{maxit}}{Maximum number of EM iterations (default: 200)}
#'     \item{\code{tol}}{Convergence tolerance on max absolute
#'       parameter change (default: 1e-4)}
#'     \item{\code{verbose}}{Verbosity level: FALSE/0 for silent, TRUE/1 for
#'       progress messages, 2 for detailed output (default: FALSE)}
#'     \item{\code{parallel}}{Logical flag enabling parallel computation
#'       (default: FALSE)}
#'     \item{\code{n_cores}}{Number of CPU cores for parallel processing.
#'       If 0, automatically detects available cores (default: 0)}
#'   }
#'   See \code{\link{JointODE.control}} for complete details and examples.
#' @param ... Additional arguments passed to internal optimization routines.
#'
#' @return An S3 object of class \code{"JointODE"} containing fitted model
#'   results:
#'   \describe{
#'     \item{\code{parameters}}{A list containing all estimated parameters:
#'       \itemize{
#'         \item \code{coefficients}: Named list with \code{baseline} (B-spline
#'           coefficients for baseline hazard), \code{hazard} (association and
#'           survival covariate effects), \code{longitudinal} (longitudinal
#'           fixed effects), \code{measurement_error_sd} (residual standard
#'           deviation), and \code{random_effect_sigma} (random effect
#'           covariance matrix)
#'         \item \code{configurations}: Model configuration including spline
#'           basis specifications
#'       }}
#'     \item{\code{logLik}}{Maximum log-likelihood value achieved at
#'       convergence}
#'     \item{\code{AIC}}{Akaike Information Criterion for model comparison}
#'     \item{\code{BIC}}{Bayesian Information Criterion adjusted for sample
#'       size}
#'     \item{\code{cindex}}{Concordance index (C-index) measuring the model's
#'       discrimination ability for survival prediction}
#'     \item{\code{convergence}}{List containing convergence diagnostics:
#'       \itemize{
#'         \item \code{converged}: Logical indicating convergence status
#'         \item \code{iterations}: Number of EM iterations performed
#'         \item \code{message}: Descriptive convergence message
#'       }}
#'     \item{\code{random_effects}}{Matrix of posterior mode random effects
#'       (n_subjects x n_re)}
#'     \item{\code{data}}{Processed data used for model fitting in internal
#'       format}
#'     \item{\code{control}}{List of control parameters used in optimization}
#'     \item{\code{call}}{The matched function call for reproducibility}
#'   }
#'
#' @details
#' The joint modeling framework integrates longitudinal and survival processes
#' through a shared random effects structure. The longitudinal biomarker
#' evolution is characterized by a system of ODEs that can accommodate
#' non-linear dynamics, feedback mechanisms, and complex temporal patterns.
#' The survival component employs a proportional hazards model where the
#' instantaneous risk depends on
#' features derived from the longitudinal trajectory.
#'
#' Two association structures are supported:
#' \itemize{
#'   \item Current value: hazard depends on the biomarker level at time t
#'   \item Rate of change: hazard depends on the biomarker's instantaneous
#'     slope
#' }
#'
#' Parameter estimation employs a Laplace EM (PQL) algorithm with:
#' \itemize{
#'   \item E-step: Laplace approximation for posterior mode and covariance
#'     of random effects
#'   \item M-step: one-step Newton update on the complete-data log-likelihood
#'     evaluated at posterior modes
#' }
#'
#' @importFrom utils modifyList
#' @importFrom survival Surv
#' @importFrom cli cli_h2 cli_text cli_alert_success
#' @importFrom cli cli_alert_warning cli_alert_info
#' @importFrom future.apply future_lapply
#' @importFrom numDeriv jacobian
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' sim <- simulate(JointODE, n_subjects = 100)
#'
#' # Fit with default control parameters
#' fit1 <- JointODE(
#'   longitudinal_formula = observed ~ x1 + x2,
#'   survival_formula = Surv(time, status) ~ w1 + w2,
#'   longitudinal_data = sim$data$longitudinal_data,
#'   survival_data = sim$data$survival_data
#' )
#'
#' # Fit with custom control parameters using JointODE.control()
#' control <- JointODE.control(
#'   maxit = 200, tol = 1e-4, verbose = TRUE
#' )
#' fit2 <- JointODE(
#'   longitudinal_formula = observed ~ x1 + x2,
#'   survival_formula = Surv(event_time, event) ~ x1 + x2,
#'   longitudinal_data = sim$longitudinal_data,
#'   survival_data = sim$survival_data,
#'   control = control
#' )
#'
#' # Fit with control parameters as a list
#' # By default, uses MarginalODE for initialization (init = NULL)
#' fit3 <- JointODE(
#'   longitudinal_formula = observed ~ x1 + x2,
#'   survival_formula = Surv(event_time, event) ~ x1 + x2,
#'   longitudinal_data = sim$longitudinal_data,
#'   survival_data = sim$survival_data,
#'   control = list(maxit = 50, verbose = TRUE)
#' )
#'
#' summary(fit1)
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
    parameters <- .initialize_from_marginal(
      longitudinal_data, survival_data, gamma, control,
      parsed_long, parsed_surv, model_config
    )
  } else if (is.list(init)) {
    parameters <- init
  } else {
    parameters <- .default_parameters(
      model_config$dims, gamma, parsed_long,
      model_config$spline_baseline_config
    )
  }

  parameters$configurations$baseline <- model_config$spline_baseline_config
  response <- stats::model.response(stats::model.frame(
    .build_formula(parsed_long$fixed_terms, response = parsed_long$response),
    longitudinal_data
  ))
  parameters$configurations$biomarker_clamp <- max(abs(response)) * 5
  parameters$configurations$hazard_quadrature <- control$hazard_quadrature

  names(parameters$coefficients$baseline) <- model_config$coef_names$baseline
  names(parameters$coefficients$hazard) <- model_config$coef_names$hazard
  names(parameters$coefficients$longitudinal) <-
    model_config$coef_names$longitudinal
  names(parameters$coefficients$initial_state) <-
    model_config$coef_names$initial_state

  random_effects <- model_config$random_effects
  coef_names <- model_config$coef_names

  # Initialize state random effects from first observations
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

  # Laplace EM Algorithm
  if (control$verbose > 0) {
    cli::cli_h2("Joint ODE Model Estimation")
    cli::cli_text(sprintf(
      "Convergence: tol=%.1e, max_iter=%d",
      control$tol,
      control$maxit
    ))
    cli::cli_text("")
  }

  # Set up parallel plan once before EM loop
  if (control$parallel) {
    parallel_cleanup <- .setup_parallel_plan(control$n_cores)
    on.exit(parallel_cleanup(), add = TRUE)
  }

  # EM loop
  curr <- list(
    parameters = parameters,
    random_effects = random_effects,
    loglik = -Inf
  )
  prev <- curr
  converged <- FALSE
  loglik_history <- rep(NA_real_, control$maxit)
  delta_theta_history <- rep(NA_real_, control$maxit)
  delta_loglik_history <- rep(NA_real_, control$maxit)

  for (em_iter in seq_len(control$maxit)) {
    curr <- .em_step(
      data_list, curr$parameters, curr$random_effects, control
    )

    status <- .track(em_iter, curr, prev, control)
    loglik_history[em_iter] <- curr$loglik
    delta_theta_history[em_iter] <- status$metrics$delta_theta
    delta_loglik_history[em_iter] <- status$metrics$delta_l

    if (status$converged) {
      converged <- TRUE
      break
    }
    prev <- curr
  }

  n_iter <- if (converged) em_iter else control$maxit

  # Finalize model
  final_results <- .finalize_joint(
    data_list = data_list,
    parameters = curr$parameters,
    loglik = curr$loglik,
    control = control,
    coef_names = coef_names,
    converged = converged,
    random_effects = curr$random_effects
  )

  # Return fitted model
  structure(
    list(
      parameters = final_results$parameters,
      logLik = final_results$loglik,
      AIC = final_results$aic,
      BIC = final_results$bic,
      cindex = final_results$cindex,
      convergence = list(
        converged = converged,
        iterations = n_iter,
        message = sprintf(
          "%s after %d iterations",
          if (converged) "Converged" else "Did not converge",
          n_iter
        ),
        loglik_history = loglik_history[seq_len(n_iter)],
        delta_theta_history = delta_theta_history[seq_len(n_iter)],
        delta_loglik_history = delta_loglik_history[seq_len(n_iter)]
      ),
      random_effects = final_results$random_effects,
      vcov = final_results$vcov,
      data = data_list,
      control = control,
      call = cl
    ),
    class = "JointODE"
  )
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
    # Extract coefficients from the longitudinal ODE equation:
    # d²y/dt² = β₁ * y + β₂ * dy/dt + ...
    # Comparing with damped harmonic oscillator: ẍ = -ω²x - 2ξωẋ + kω²f
    # We have: β₁ = -ω² and β₂ = -2ξω

    value_coef <- coefs[idx_longitudinal[1]] # β₁ = -ω²
    slope_coef <- coefs[idx_longitudinal[2]] # β₂ = -2ξω

    # Variances and covariance
    var_value <- object$vcov[idx_longitudinal[1], idx_longitudinal[1]]
    var_slope <- object$vcov[idx_longitudinal[2], idx_longitudinal[2]]
    cov_value_slope <- object$vcov[idx_longitudinal[1], idx_longitudinal[2]]

    # Check if value_coef is negative (as expected for -ω²)
    if (value_coef < 0) {
      # Calculate omega_n from β₁ = -ωₙ²
      # ωₙ = √(-β₁)
      omega_est <- sqrt(-value_coef)

      # Calculate period T = 2π/ωₙ
      period_est <- 2 * pi / omega_est

      # Calculate xi from β₂ = -2ξωₙ
      # ξ = -β₂ / (2ωₙ) = -β₂ / (2√(-β₁))
      xi_est <- -slope_coef / (2 * omega_est)

      # Delta method for period
      # T = 2π/ωₙ = 2π/√(-β₁)
      # ∂T/∂β₁ = π/((-β₁)^(3/2))
      grad_period_value <- pi / ((-value_coef)^(3 / 2))
      var_period <- grad_period_value^2 * var_value
      se_period <- sqrt(var_period)

      # Delta method for xi
      # ξ = -β₂/(2√(-β₁))
      # ∂ξ/∂β₁ = -β₂/(4*(-β₁)^(3/2))
      # ∂ξ/∂β₂ = -1/(2√(-β₁))
      grad_xi_value <- -slope_coef / (4 * (-value_coef)^(3 / 2))
      grad_xi_slope <- -1 / (2 * sqrt(-value_coef))

      # Variance of xi using Delta method with covariance
      var_xi <- grad_xi_value^2 *
        var_value +
        grad_xi_slope^2 * var_slope +
        2 * grad_xi_value * grad_xi_slope * cov_value_slope
      se_xi <- sqrt(var_xi)

      # Create coefficient matrix for derived parameters
      derived_params <- cbind(
        Estimate = c(period_est, xi_est),
        `Std. Error` = c(se_period, se_xi),
        `z value` = c(period_est / se_period, xi_est / se_xi),
        `Pr(>|z|)` = 2 * pnorm(-abs(c(period_est / se_period, xi_est / se_xi)))
      )
      rownames(derived_params) <- c("T (period)", "xi (damping ratio)")
    }
  }

  # Count longitudinal observations and events from data
  n_subjects <- length(object$data)
  n_observations <- sum(vapply(object$data, function(subject) {
    length(subject$longitudinal$measurements)
  }, integer(1)))
  n_events <- sum(vapply(object$data, function(subject) {
    subject$status
  }, numeric(1)))
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

  # Data Descriptives
  cat("\nData Descriptives:\n")
  cat("Longitudinal Process            Survival Process\n")
  cat(sprintf(
    "Number of Observations: %-7d Number of Events: %d (%.0f%%)\n",
    x$n_observations,
    x$n_events,
    x$event_rate * 100
  ))
  cat(sprintf("Number of Subjects: %d\n", x$nobs))

  # Model fit statistics
  cat(sprintf("\n%10s %10s %10s\n", "AIC", "BIC", "logLik"))
  cat(sprintf("%10.3f %10.3f %10.3f\n", x$AIC, x$BIC, x$logLik))

  cat("\nCoefficients:\n")

  # Longitudinal Process (ODE model)
  cat("Longitudinal Process: Second-Order ODE Model\n")
  if (!is.null(x$coef_longitudinal)) {
    .print_coefmat(
      x$coef_longitudinal,
      digits = digits,
      signif.stars = signif.stars,
      ...
    )
  }

  # Derived ODE characteristics
  if (!is.null(x$derived_params)) {
    cat("\nODE System Characteristics:\n")
    .print_coefmat(
      x$derived_params,
      digits = digits,
      signif.stars = signif.stars,
      ...
    )
  }

  # Survival Process
  cat("\nSurvival Process: Proportional Hazards Model\n")
  if (!is.null(x$coef_survival)) {
    .print_coefmat(
      x$coef_survival,
      digits = digits,
      signif.stars = signif.stars,
      ...
    )
  }

  # Baseline hazard (spline coefficients - optional, summarized)
  if (!is.null(x$coef_baseline)) {
    cat(
      "\nBaseline Hazard: B-spline with",
      nrow(x$coef_baseline),
      "basis functions\n"
    )
    cat(
      "(Coefficients range:",
      sprintf(
        "[%.3f, %.3f]",
        min(x$coef_baseline[, "Estimate"]),
        max(x$coef_baseline[, "Estimate"])
      ),
      ")\n"
    )
  }

  # Initial State
  if (!is.null(x$coef_initial)) {
    cat("\nInitial State: Population Mean\n")
    .print_coefmat(
      x$coef_initial,
      digits = digits, signif.stars = signif.stars, ...
    )
  }

  # Variance Components
  cat("\nVariance Components:\n")
  cat(sprintf("Measurement Error SD: %.6f\n", x$sigma["sigma_e"]))
  if (!is.null(x$sigma_b_matrix)) {
    cat("Random Effect Covariance Matrix:\n")
    print(x$sigma_b_matrix, digits = 4)
  }

  # Model diagnostics
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
  names(coefs) <- c(
    paste0("baseline:", names(cf$baseline)),
    paste0("hazard:", names(cf$hazard)),
    paste0("longitudinal:", names(cf$longitudinal)),
    paste0("initial state:", names(cf$initial_state))
  )
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
  cat("\nJoint ODE Model\n")
  cat("Call: ")
  print(x$call)

  n_params <- .count_params(x$parameters)
  cat(
    "\nLog-likelihood:",
    format(x$logLik, digits = digits),
    "on",
    n_params,
    "degrees of freedom\n"
  )
  cat(
    "AIC:",
    format(x$AIC, digits = digits),
    "  BIC:",
    format(x$BIC, digits = digits),
    "\n"
  )
  invisible(x)
}

#' Predict Method for JointODE Objects
#'
#' @description
#' Computes predicted biomarker trajectories, velocities, and accelerations
#' for subjects based on the fitted joint ODE model. Predictions incorporate
#' both fixed effects and subject-specific random effects.
#'
#' @param object An object of class \code{JointODE}
#' @param newdata Optional data frame with new subjects. If NULL, uses the
#'   training data from the model fit.
#' @param times Optional time points for prediction. Can be:
#'   \itemize{
#'     \item A numeric vector for the same times across all subjects
#'     \item A named list with subject-specific time vectors
#'     \item NULL to use observation times from the data (default)
#'   }
#' @param parallel Logical flag for parallel computation (default: FALSE)
#' @param n_cores Number of cores for parallel processing. If 0, automatically
#'   detects available cores (default: 0)
#' @param ... Additional arguments (currently unused)
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{\code{id}}{Subject identifier}
#'     \item{\code{time}}{Time points for predictions}
#'     \item{\code{cumhaz}}{Predicted cumulative hazard at each time point}
#'     \item{\code{survival}}{Predicted survival probability, computed as
#'       \eqn{S(t) = \exp(-\text{cumhaz})}}
#'     \item{\code{log_hazard}}{Log instantaneous hazard at each time point}
#'     \item{\code{biomarker}}{Predicted biomarker values}
#'     \item{\code{velocity}}{Predicted velocity (first derivative)}
#'     \item{\code{acceleration}}{Predicted acceleration (second derivative)}
#'   }
#'
#' @concept model-prediction
#' @export
predict.JointODE <- function(
  object,
  newdata = NULL,
  times = NULL,
  parallel = FALSE,
  n_cores = 0,
  ...
) {
  if (!is.null(newdata)) {
    stop("newdata not yet supported")
  }

  parameters <- object$parameters
  data_list <- object$data

  random_effects <- object$random_effects

  # Prepare prediction data for all subjects
  pred_data_list <- lapply(seq_along(data_list), function(i) {
    subj <- data_list[[i]]
    subj_id <- names(data_list)[i]

    pred_times <- if (!is.null(times)) {
      time_vec <- if (is.list(times)) {
        if (!is.null(times[[subj_id]])) {
          times[[subj_id]]
        } else {
          subj$longitudinal$times
        }
      } else {
        times
      }
      sort(unique(time_vec))
    } else {
      subj$longitudinal$times
    }

    if (length(pred_times) == 0) {
      return(NULL)
    }

    # Extend longitudinal covariates to match pred_times
    # Use last observation carried forward (LOCF) for times beyond observed data
    orig_times <- subj$longitudinal$times

    extended_fixed <- .extend_covariates(
      subj$longitudinal$covariates$fixed, orig_times, pred_times
    )
    extended_random <- .extend_covariates(
      subj$longitudinal$covariates$random, orig_times, pred_times
    )

    extended_time <- if (!is.null(times)) {
      max(pred_times)
    } else {
      subj$time
    }

    list(
      id = subj$id,
      time = extended_time,
      status = subj$status,
      covariates = subj$covariates,
      longitudinal = list(
        times = pred_times,
        measurements = rep(0, length(pred_times)),
        covariates = list(
          fixed = extended_fixed,
          random = extended_random
        )
      )
    )
  })

  # Remove NULL entries
  valid_idx <- !vapply(pred_data_list, is.null, logical(1))
  pred_data_list <- pred_data_list[valid_idx]
  valid_random_effects <- random_effects[valid_idx, , drop = FALSE]
  valid_ids <- names(data_list)[valid_idx]

  # Solve ODE for all subjects at once
  ode_solutions <- .solve_batch_joint(
    data_list = pred_data_list,
    random_effects = valid_random_effects,
    parameters = parameters
  )

  # Combine results
  results <- lapply(seq_along(ode_solutions), function(i) {
    data.frame(
      id = valid_ids[i],
      time = ode_solutions[[i]]$times,
      cumhaz = ode_solutions[[i]]$cum_hazard,
      survival = exp(-ode_solutions[[i]]$cum_hazard),
      log_hazard = ode_solutions[[i]]$log_hazard,
      biomarker = ode_solutions[[i]]$biomarker,
      velocity = ode_solutions[[i]]$velocity,
      acceleration = ode_solutions[[i]]$acceleration,
      stringsAsFactors = FALSE
    )
  })

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  result_df
}
