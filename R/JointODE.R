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
#' @param state A matrix specifying initial conditions for the ODE system with
#'   two columns: initial biomarker values and their first derivatives. Each row
#'   corresponds to one subject. If \code{NULL}, defaults to a zero matrix with
#'   appropriate dimensions (default: \code{NULL}).
#' @param spline_baseline A list controlling the B-spline representation of the
#'   baseline hazard function with the following components:
#'   \describe{
#'     \item{\code{degree}}{Polynomial degree of the B-spline basis functions
#'       (default: 2, quadratic splines)}
#'     \item{\code{n_knots}}{Number of interior knots for flexibility
#'       (default: 0, providing moderate flexibility)}
#'     \item{\code{knot_placement}}{Strategy for positioning knots:
#'       \code{"quantile"} places knots at quantiles of observed event times,
#'       \code{"equal"} uses equally-spaced knots (default: \code{"equal"})}
#'     \item{\code{boundary_knots}}{A numeric vector of length 2 specifying
#'       the boundary knot locations. If \code{NULL}, automatically set to the
#'       range of observed event times (default: \code{NULL})}
#'   }
#' @param init Optional initial values for model parameters. Can be either:
#'   \itemize{
#'     \item \code{NULL} (default): Automatically calls
#'       \code{\link{MarginalODE}} to compute data-driven initial estimates
#'       for longitudinal parameters and random effects. Recommended for
#'       most applications.
#'     \item A list with the same structure as the fitted model's
#'       \code{parameters} component for full manual control
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
#'     \item{\code{atol}}{Absolute tolerance for parameter convergence
#'       (default: 5e-4)}
#'     \item{\code{rtol}}{Relative tolerance for log-likelihood
#'       convergence (default: 1e-5)}
#'     \item{\code{verbose}}{Verbosity level: FALSE/0 for silent, TRUE/1 for
#'       progress messages, 2 for detailed output (default: FALSE)}
#'     \item{\code{parallel}}{Logical flag enabling parallel computation
#'       (default: FALSE)}
#'     \item{\code{n_cores}}{Number of CPU cores for parallel processing.
#'       If 0, automatically detects available cores (default: 0)}
#'     \item{\code{quad_level}}{Quadrature level for numerical integration
#'       (default: 4)}
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
#'         \item \code{em_iterations}: Number of EM iterations performed
#'         \item \code{message}: Descriptive convergence message
#'       }}
#'     \item{\code{random_effects}}{List containing random effects estimates:
#'       \itemize{
#'         \item \code{estimates}: Posterior means of subject-specific random
#'           effects
#'         \item \code{variances}: Posterior variances of random effects
#'       }}
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
#' Parameter estimation employs an Expectation-Maximization (EM) algorithm
#' with:
#' \itemize{
#'   \item E-step: Multivariate Gauss-Hermite quadrature (mvQuad) for
#'     posterior computation of random effects
#'   \item M-step: Optimization for fixed effects and mvQuad-based
#'     closed-form updates for variance parameters
#' }
#'
#' @importFrom stats optim
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
#' sim <- simulate(n_subjects = 100, n_times = 10)
#'
#' # Fit with default control parameters
#' fit1 <- JointODE(
#'   longitudinal_formula = observed ~ x1 + x2,
#'   survival_formula = Surv(event_time, event) ~ x1 + x2,
#'   longitudinal_data = sim$longitudinal_data,
#'   survival_data = sim$survival_data
#' )
#'
#' # Fit with custom control parameters using JointODE.control()
#' control <- JointODE.control(
#'   maxit = 200, atol = 1e-4, rtol = 1e-6, verbose = TRUE
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
  state = NULL,
  spline_baseline = list(
    degree = 2,
    n_knots = 0,
    knot_placement = "equal",
    boundary_knots = NULL
  ),
  init = NULL,
  control = list(),
  ...
) {
  cl <- match.call()

  # Validate and process data
  .validate(
    longitudinal_formula = longitudinal_formula,
    survival_formula = survival_formula,
    longitudinal_data = longitudinal_data,
    survival_data = survival_data,
    state = state,
    gamma = gamma,
    spline_baseline = spline_baseline,
    init = init
  )

  # Process data
  data_list <- .process(
    longitudinal_formula = longitudinal_formula,
    survival_formula = survival_formula,
    longitudinal_data = longitudinal_data,
    survival_data = survival_data,
    state = state
  )

  # Determine if we need to optimize initial states
  has_state <- !is.null(state)

  # Process control settings
  if (is.null(control)) {
    control <- JointODE.control()
  } else if (is.list(control)) {
    control <- JointODE.control(.list = control)
  } else {
    stop("control must be a list or NULL")
  }

  # Initialize parameters
  init_state <- .initialize(
    longitudinal_formula,
    survival_formula,
    longitudinal_data,
    survival_data,
    spline_baseline,
    gamma,
    state,
    init,
    control
  )
  parameters <- init_state$parameters
  random_effects <- init_state$random_effects
  coef_names <- init_state$coef_names

  # EM Algorithm
  if (control$verbose > 0) {
    cli::cli_h2("Joint ODE Model Estimation")
    cli::cli_text(sprintf(
      "Convergence: atol=%.1e, rtol=%.1e, max_iter=%d",
      control$atol,
      control$rtol,
      control$maxit
    ))
    cli::cli_text("")
  }

  # EM loop
  curr <- list(
    parameters = parameters,
    random_effects = random_effects,
    loglik = -Inf
  )
  prev <- curr
  converged <- FALSE

  for (em_iter in seq_len(control$maxit)) {
    if (!has_state) {
      opt_results <- .parallel_apply(
        seq_along(data_list),
        function(i) {
          .estimate_state(
            data_list[[i]]$initial_state,
            data_list[[i]],
            curr$random_effects[i, ],
            curr$parameters,
            tol = control$atol
          )
        },
        control$parallel,
        control$n_cores
      )
      for (i in seq_along(data_list)) {
        data_list[[i]]$initial_state <- opt_results[[i]]$state
      }

      if (control$verbose > 0) {
        obj_changes <- vapply(opt_results, `[[`, numeric(1), "obj_change")
        n_increased <- sum(obj_changes > 1e-6)

        base_msg <- sprintf(
          "    State: dL=%.2e [%.2e, %.2e]",
          mean(obj_changes),
          min(obj_changes),
          max(obj_changes)
        )

        if (n_increased > 0) {
          cli::cli_alert_warning(
            sprintf("%s | %d increased", base_msg, n_increased)
          )
        } else {
          cli::cli_text(base_msg)
        }
      }
    }

    # EM step
    curr <- .em_step(
      data_list,
      curr$parameters,
      curr$random_effects,
      control
    )

    # Track progress and check convergence
    status <- .track(em_iter, curr, prev, control)

    if (status$converged) {
      converged <- TRUE
      break
    }

    prev <- curr

    gc(verbose = FALSE, full = TRUE)
  }

  # Finalize model
  final_results <- .finalize(
    data_list = data_list,
    parameters = curr$parameters,
    random_effects = curr$random_effects,
    loglik = curr$loglik,
    control = control,
    coef_names = coef_names,
    converged = converged
  )

  # Return fitted model
  structure(
    list(
      parameters = curr$parameters,
      logLik = final_results$loglik,
      AIC = final_results$aic,
      BIC = final_results$bic,
      cindex = final_results$cindex,
      convergence = list(
        converged = converged,
        em_iterations = if (converged) em_iter else control$maxit,
        message = sprintf(
          "EM algorithm %s after %d iterations",
          if (converged) "converged" else "did not converge within",
          if (converged) em_iter else control$maxit
        )
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

  # Baseline hazard coefficients
  idx_baseline <- seq_len(n_baseline)
  coef_baseline <- cbind(
    Estimate = coefs[idx_baseline],
    `Std. Error` = se[idx_baseline],
    `z value` = coefs[idx_baseline] / se[idx_baseline],
    `Pr(>|z|)` = 2 * pnorm(-abs(coefs[idx_baseline] / se[idx_baseline]))
  )
  rownames(coef_baseline) <- gsub("baseline:", "", rownames(coef_baseline))

  # Survival process coefficients (hazard parameters)
  idx_survival <- n_baseline + seq_len(n_hazard)
  coef_survival <- cbind(
    Estimate = coefs[idx_survival],
    `Std. Error` = se[idx_survival],
    `z value` = coefs[idx_survival] / se[idx_survival],
    `Pr(>|z|)` = 2 * pnorm(-abs(coefs[idx_survival] / se[idx_survival]))
  )
  rownames(coef_survival) <- gsub("hazard:", "", rownames(coef_survival))

  # Longitudinal process coefficients (ODE parameters)
  idx_longitudinal <- (n_baseline + n_hazard) + seq_len(n_longitudinal)
  coef_longitudinal <- cbind(
    Estimate = coefs[idx_longitudinal],
    `Std. Error` = se[idx_longitudinal],
    `z value` = coefs[idx_longitudinal] / se[idx_longitudinal],
    `Pr(>|z|)` = 2 * pnorm(-abs(coefs[idx_longitudinal] / se[idx_longitudinal]))
  )
  rownames(coef_longitudinal) <- gsub(
    "longitudinal:",
    "",
    rownames(coef_longitudinal)
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
  n_observations <- sum(sapply(object$data, function(subject) {
    length(subject$longitudinal$measurements)
  }))
  n_events <- sum(sapply(object$data, function(subject) {
    subject$status
  }))
  event_rate <- n_events / n_subjects

  structure(
    list(
      call = object$call,
      coef_baseline = coef_baseline,
      coef_survival = coef_survival,
      coef_longitudinal = coef_longitudinal,
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
    printCoefmat(
      x$coef_longitudinal,
      digits = digits,
      signif.stars = signif.stars,
      ...
    )
  }

  # Derived ODE characteristics
  if (!is.null(x$derived_params)) {
    cat("\nODE System Characteristics:\n")
    printCoefmat(
      x$derived_params,
      digits = digits,
      signif.stars = signif.stars,
      ...
    )
  }

  # Survival Process
  cat("\nSurvival Process: Proportional Hazards Model\n")
  if (!is.null(x$coef_survival)) {
    printCoefmat(
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
  coefs <- c(cf$baseline, cf$hazard, cf$longitudinal)
  names(coefs) <- c(
    paste0("baseline:", names(cf$baseline)),
    paste0("hazard:", names(cf$hazard)),
    paste0("longitudinal:", names(cf$longitudinal))
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
  # Count parameters: fixed effects + sigma_e + Sigma_b (symmetric matrix)
  p_random <- nrow(object$parameters$coefficients$random_effect_sigma)
  n_sigma_b_params <- p_random * (p_random + 1) / 2
  n_params <- length(coef(object)) + 1 + n_sigma_b_params

  structure(
    object$logLik,
    df = n_params,
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

  # Count parameters: fixed effects + sigma_e + Sigma_b (symmetric matrix)
  p_random <- nrow(x$parameters$coefficients$random_effect_sigma)
  n_sigma_b_params <- p_random * (p_random + 1) / 2
  n_params <- length(coef(x)) + 1 + n_sigma_b_params

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

  # Handle both old (matrix) and new (list) random effects format
  random_effects <- if (
    is.list(object$random_effects) &&
      !is.null(object$random_effects$estimates)
  ) {
    object$random_effects$estimates
  } else {
    object$random_effects
  }

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

    extend_covariates_matrix <- function(cov_mat, orig_times, pred_times) {
      if (is.null(cov_mat) || length(cov_mat) == 0) {
        # Return empty matrix with correct dimensions
        if (is.matrix(cov_mat)) {
          return(matrix(numeric(0), nrow = length(pred_times), ncol = 0))
        } else {
          return(cov_mat) # Keep original structure if not a matrix
        }
      }

      # For each pred_time, find the corresponding row in cov_mat using LOCF
      indices <- sapply(pred_times, function(t) {
        if (t <= max(orig_times)) {
          # Find the closest time <= t
          max(which(orig_times <= t))
        } else {
          # Use last observation
          length(orig_times)
        }
      })

      # Extract rows based on indices
      if (is.matrix(cov_mat)) {
        cov_mat[indices, , drop = FALSE]
      } else {
        # If it's a vector, return as vector
        cov_mat[indices]
      }
    }

    extended_fixed <- extend_covariates_matrix(
      subj$longitudinal$covariates$fixed,
      orig_times,
      pred_times
    )

    extended_random <- extend_covariates_matrix(
      subj$longitudinal$covariates$random,
      orig_times,
      pred_times
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
      initial_state = subj$initial_state,
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
  valid_idx <- !sapply(pred_data_list, is.null)
  pred_data_list <- pred_data_list[valid_idx]
  valid_random_effects <- random_effects[valid_idx, , drop = FALSE]
  valid_ids <- names(data_list)[valid_idx]

  # Solve ODE for all subjects at once
  ode_solutions <- .solve_batch_ode_cppad(
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
