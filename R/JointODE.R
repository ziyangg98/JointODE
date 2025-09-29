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
#' @param state A matrix specifying initial conditions for the ODE system with
#'   two columns: initial biomarker values and their first derivatives. Each row
#'   corresponds to one subject. If \code{NULL}, defaults to a zero matrix with
#'   appropriate dimensions (default: \code{NULL}).
#' @param id Character string specifying the column name for subject
#'   identifiers. This variable must be present in both longitudinal and
#'   survival datasets
#'   to link observations (default: \code{"id"}).
#' @param time Character string specifying the column name for measurement
#'   times
#'   in the longitudinal dataset (default: \code{"time"}).
#' @param spline_baseline A list controlling the B-spline representation of the
#'   baseline hazard function with the following components:
#'   \describe{
#'     \item{\code{degree}}{Polynomial degree of the B-spline basis functions
#'       (default: 3, cubic splines)}
#'     \item{\code{n_knots}}{Number of interior knots for flexibility
#'       (default: 3, providing moderate flexibility)}
#'     \item{\code{knot_placement}}{Strategy for positioning knots:
#'       \code{"quantile"} places knots at quantiles of observed event times,
#'       \code{"equal"} uses equally-spaced knots (default: \code{"quantile"})}
#'     \item{\code{boundary_knots}}{A numeric vector of length 2 specifying
#'       the boundary knot locations. If \code{NULL}, automatically set to the
#'       range of observed event times (default: \code{NULL})}
#'   }
#' @param robust Logical flag indicating whether to use robust variance
#'   estimation (default: \code{FALSE}).
#' @param init Optional list providing initial values for model parameters.
#'   Should have the same structure as the fitted model's \code{parameters}
#'   component with elements:
#'   \describe{
#'     \item{\code{coefficients}}{A list containing:
#'       \itemize{
#'         \item \code{baseline}: Vector of B-spline coefficients for baseline
#'           hazard (length = number of spline basis functions)
#'         \item \code{hazard}: Vector of hazard parameters including
#'           association parameters (2) and survival covariates
#'         \item \code{acceleration}: Vector of longitudinal fixed effects
#'           including intercept and covariates
#'         \item \code{measurement_error_sd}: Residual standard deviation
#'           (positive scalar)
#'         \item \code{random_effect_sd}: Random effect standard deviation
#'           (positive scalar)
#'       }}
#'     \item{\code{configurations}}{Optional; if not provided, will use
#'       spline configuration from \code{spline_baseline}}
#'   }
#'   If \code{NULL}, default initial values are used (default: \code{NULL}).
#' @param control A list of control parameters for optimization, or output from
#'   \code{\link{JointODE.control}}. Key parameters include:
#'   \describe{
#'     \item{\code{maxit}}{Maximum number of EM iterations (default: 100)}
#'     \item{\code{tol}}{Convergence tolerance for EM algorithm (default: 1e-6)}
#'     \item{\code{verbose}}{Verbosity level: FALSE/0 for silent, TRUE/1 for
#'       progress messages, 2 for detailed output (default: FALSE)}
#'     \item{\code{optim.maxit}}{Maximum iterations for M-step optimization
#'       (default: 50)}
#'     \item{\code{optim.factr}}{Convergence factor for L-BFGS-B optimizer
#'       (default: 1e10)}
#'   }
#'   See \code{\link{JointODE.control}} for more details.
#' @param parallel Logical flag enabling parallel computation for
#'   computationally intensive operations including posterior calculations,
#'   gradient evaluations, and likelihood computations. Requires \pkg{future}
#'   and \pkg{future.apply} packages
#'   (default: \code{FALSE}).
#' @param n_cores Integer specifying the number of CPU cores for parallel
#'   processing. If \code{NULL}, automatically detects and uses all available
#'   cores minus one
#'   (default: \code{NULL}).
#' @param ... Additional arguments passed to internal optimization routines.
#'
#' @return An S3 object of class \code{"JointODE"} containing fitted model
#'   results:
#'   \describe{
#'     \item{\code{parameters}}{A list containing all estimated parameters:
#'       \itemize{
#'         \item \code{coefficients}: Named list with \code{baseline} (B-spline
#'           coefficients for baseline hazard), \code{hazard} (association and
#'           survival covariate effects), \code{acceleration} (longitudinal
#'           fixed effects), \code{measurement_error_sd} (residual standard
#'           deviation), and \code{random_effect_sd} (random effect standard
#'           deviation)
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
#'   \item E-step: Adaptive Gauss-Hermite quadrature for numerical integration
#'   \item M-step: Quasi-Newton optimization for parameter updates
#' }
#'
#' @note
#' \itemize{
#'   \item Data validation is performed automatically with informative error
#'     messages
#'   \item For high-dimensional problems, parallel computation is strongly
#'     recommended
#'   \item Convergence issues may arise with sparse event data or limited
#'     follow-up
#'   \item Initial values are computed using separate model fits when not
#'     provided
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
#'   survival_data = sim$survival_data,
#'   id = "id"
#' )
#'
#' # Fit with custom control parameters using JointODE.control()
#' control <- JointODE.control(maxit = 200, tol = 1e-8, verbose = TRUE)
#' fit2 <- JointODE(
#'   longitudinal_formula = observed ~ x1 + x2,
#'   survival_formula = Surv(event_time, event) ~ x1 + x2,
#'   longitudinal_data = sim$longitudinal_data,
#'   survival_data = sim$survival_data,
#'   id = "id",
#'   control = control
#' )
#'
#' # Fit with control parameters as a list
#' fit3 <- JointODE(
#'   longitudinal_formula = observed ~ x1 + x2,
#'   survival_formula = Surv(event_time, event) ~ x1 + x2,
#'   longitudinal_data = sim$longitudinal_data,
#'   survival_data = sim$survival_data,
#'   id = "id",
#'   control = list(maxit = 50, verbose = TRUE)
#' )
#'
#' summary(fit1)
#' }
#'
#' @concept model-fitting
#' @export
JointODE <- function(
  longitudinal_formula,
  survival_formula,
  longitudinal_data,
  survival_data,
  n_groups = 2,
  state = NULL,
  id = "id",
  time = "time",
  spline_baseline = list(
    degree = 3,
    n_knots = 3,
    knot_placement = "equal",
    boundary_knots = NULL
  ),
  robust = FALSE,
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
    id = id,
    time = time,
    spline_baseline = spline_baseline,
    init = init
  )

  data_list <- .process(
    longitudinal_formula = longitudinal_formula,
    survival_formula = survival_formula,
    longitudinal_data = longitudinal_data,
    survival_data = survival_data,
    state = state,
    id = id,
    time = time
  )

  # Extract data dimensions
  event_times <- vapply(data_list, `[[`, numeric(1), "time")
  n_subjects <- attr(data_list, "n_subjects")
  subjects_with_long <- Filter(
    function(s) s$longitudinal$n_obs > 0,
    data_list
  )
  n_longitudinal_covariates <- if (length(subjects_with_long) > 0) {
    ncol(subjects_with_long[[1]]$longitudinal$covariates)
  } else {
    0
  }
  n_survival_covariates <- ncol(subjects_with_long[[1]]$covariates)

  # Configure splines
  spline_baseline_config <- .get_spline_config(
    x = event_times,
    degree = spline_baseline$degree,
    n_knots = spline_baseline$n_knots,
    knot_placement = spline_baseline$knot_placement,
    boundary_knots = spline_baseline$boundary_knots
  )
  spline_baseline_config$boundary_knots[1] <- 0

  # Extract variable names from formulas
  long_vars_names <- if (n_longitudinal_covariates > 0) {
    colnames(subjects_with_long[[1]]$longitudinal$covariates)
  } else {
    character(0)
  }
  surv_vars_names <- if (n_survival_covariates > 0) {
    colnames(subjects_with_long[[1]]$covariates)
  } else {
    character(0)
  }

  # Use mathematical notation for parameter names
  hazard_names <- c("alpha_1", "alpha_2", surv_vars_names)
  accel_names <- c("-omega_n^2", "-2*xi*omega_n", long_vars_names)

  coef_names <- list(
    baseline = paste0("bs", seq_len(spline_baseline_config$df)),
    hazard = hazard_names,
    acceleration = accel_names
  )
  coef_names_expanded <- c(
    paste0("baseline:", coef_names$baseline),
    paste0("hazard:", coef_names$hazard),
    paste0("acceleration:", coef_names$acceleration)
  )

  # Initialize parameters
  parameters <- list(
    coefficients = list(
      baseline = rep(0, spline_baseline_config$df),
      hazard = rep(0, n_survival_covariates + 2),
      acceleration = matrix(
        rnorm((2 + n_longitudinal_covariates) * n_groups, mean = -1, sd = 0.1),
        nrow = 2 + n_longitudinal_covariates,
        ncol = n_groups
      ),
      membership = matrix(1 / n_groups, nrow = n_subjects, ncol = n_groups),
      prob = rep(1 / n_groups, n_groups),
      random_effect = matrix(0, nrow = n_subjects, ncol = n_groups),
      measurement_error_sd = 1,
      random_effect_sd = 1
    ),
    configurations = list(
      baseline = spline_baseline_config,
      n_groups = n_groups
    )
  )

  # Override with user-provided initial values if available
  if (!is.null(init)) {
    if (!is.null(init$coefficients)) {
      parameters$coefficients <- modifyList(
        parameters$coefficients,
        init$coefficients
      )
    }

    if (!is.null(init$configurations)) {
      parameters$configurations <- modifyList(
        parameters$configurations,
        init$configurations
      )
    }
  }

  # Control settings - use JointODE.control to handle defaults
  if (is.null(control)) {
    control <- JointODE.control()
  } else if (is.list(control)) {
    control <- JointODE.control(.list = control)
  } else {
    stop("control must be a list or NULL")
  }

  # EM Algorithm setup
  converged <- FALSE
  old_loglik <- -Inf
  n_baseline <- spline_baseline_config$df
  n_hazard <- n_survival_covariates + 2

  # Check verbose setting
  if (control$verbose > 0) {
    cli::cli_h2("EM Algorithm Optimization")
  }

  for (em_iter in seq_len(control$maxit)) {
    # E-step
    posteriors <- .compute_posteriors(
      data_list,
      parameters,
      parallel = control$parallel,
      n_cores = control$n_cores,
      return_ode_solutions = TRUE
    )
    ode_solutions <- posteriors$ode_solutions
    parameters$coefficients$random_effect <- posteriors$b
    parameters$coefficients$membership <- posteriors$w

    # M-step: Closed-form update for mixture and variance parameters
    m_step_results <- .m_step_closed_form(
      data_list,
      parameters,
      posteriors,
      ode_solutions
    )
    parameters$coefficients$prob <- m_step_results$prob
    parameters$coefficients$measurement_error_sd <-
      m_step_results$measurement_error_sd
    parameters$coefficients$random_effect_sd <-
      m_step_results$random_effect_sd

    # M-step
    res <- optim(
      par = c(
        parameters$coefficients$baseline,
        parameters$coefficients$hazard,
        as.vector(parameters$coefficients$acceleration)
      ),
      fn = .compute_objective_joint,
      gr = .compute_gradient_joint,
      method = "L-BFGS-B",
      control = list(
        factr = control$optim.factr,
        pgtol = control$tol,
        maxit = control$optim.maxit,
        trace = if (control$verbose >= 3) 3 else 0,
        REPORT = 1
      ),
      data_list = data_list,
      posteriors = posteriors,
      parameters = parameters,
      parallel = control$parallel,
      n_cores = control$n_cores
    )

    # Update parameters
    idx_end <- cumsum(c(
      n_baseline,
      n_hazard,
      length(res$par) - n_baseline - n_hazard
    ))
    idx_start <- c(1, idx_end[-length(idx_end)] + 1)
    parameters$coefficients$baseline <- res$par[idx_start[1]:idx_end[1]]
    parameters$coefficients$hazard <- res$par[idx_start[2]:idx_end[2]]
    parameters$coefficients$acceleration <- matrix(
      res$par[idx_start[3]:idx_end[3]],
      nrow = 2 + n_longitudinal_covariates,
      ncol = n_groups
    )

    # Track convergence
    new_loglik <- -res$value
    loglik_change <- new_loglik - old_loglik

    if (control$verbose > 0) {
      cli::cli_text(
        "[{sprintf('%3d', em_iter)}/{control$maxit}] ",
        "LogLik: {sprintf('%.4f', new_loglik)} ",
        "(Change: {sprintf('%.4f', loglik_change)}) ",
        "\u03c3_e={sprintf('%.3f',
                     parameters$coefficients$measurement_error_sd)} ",
        "\u03c3_b={sprintf('%.3f', parameters$coefficients$random_effect_sd)}"
      )

      # Detailed parameter output for verbose=2
      if (control$verbose >= 2) {
        cli::cli_text(
          "  Baseline hazard: [{paste(sprintf('%.3f',
            head(parameters$coefficients$baseline, 5)), collapse=', ')}",
          if (length(parameters$coefficients$baseline) > 5) "..." else "",
          "]"
        )
        cli::cli_text(
          "  Association params: [{paste(sprintf('%.3f',
            parameters$coefficients$hazard[1:2]), collapse=', ')}]"
        )
        if (n_hazard > 2) {
          cli::cli_text(
            "  Survival covariates: [{paste(sprintf('%.3f',
              parameters$coefficients$hazard[3:n_hazard]), collapse=', ')}]"
          )
        }
        cli::cli_text(
          "  ODE coefficients: [{paste(sprintf('%.3f',
            parameters$coefficients$acceleration), collapse=', ')}]"
        )
        cli::cli_text("")
      }
    }

    # Check convergence
    if (abs(loglik_change) < control$tol && em_iter > 1) {
      converged <- TRUE
      if (control$verbose > 0) {
        cli::cli_text("")
        cli::cli_alert_success(
          "Converged after {em_iter} iterations ",
          "(LogLik: {sprintf('%.4f', new_loglik)})"
        )
      }
      break
    }

    old_loglik <- new_loglik
  }

  if (!converged && control$verbose > 0) {
    cli::cli_text("")
    cli::cli_alert_warning(
      "Did not converge within {control$maxit} iterations ",
      "(LogLik: {sprintf('%.4f', new_loglik)})"
    )
    cli::cli_alert_info(
      "Try increasing em_maxit or adjusting initial values"
    )
  }

  # Final computations
  final_posteriors <- .compute_posteriors(
    data_list,
    parameters,
    parallel = control$parallel,
    n_cores = control$n_cores,
    return_ode_solutions = TRUE
  )

  names(parameters$coefficients$baseline) <- coef_names$baseline
  names(parameters$coefficients$hazard) <- coef_names$hazard
  names(parameters$coefficients$acceleration) <- coef_names$acceleration

  # Compute variance-covariance matrix
  vcov_matrix <- if (!converged) {
    if (control$verbose > 0) {
      cli::cli_alert_warning(
        "Skipping variance-covariance computation (EM did not converge)"
      )
    }
    NULL
  } else {
    tryCatch(
      {
        if (control$verbose > 0) {
          cli::cli_alert_info("Computing variance-covariance matrix...")
        }

        # Compute Hessian via numerical differentiation of gradient
        # H = d²(-logL)/dθ² ≈ J(∇(-logL))
        H <- numDeriv::jacobian(
          func = function(p) {
            .compute_gradient_joint(
              p,
              data_list = data_list,
              posteriors = final_posteriors,
              parameters = parameters,
              parallel = control$parallel,
              n_cores = control$n_cores
            )
          },
          x = res$par,
          method = "Richardson" # More accurate for second derivatives
        )

        # Ensure symmetry (numerical errors can break symmetry)
        H <- 0.5 * (H + t(H))

        # Check condition before inversion
        eigen_vals <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        if (min(eigen_vals) <= 0) {
          stop("Hessian is not positive definite")
        }

        if (robust) {
          # Robust sandwich estimator
          score_list <- lapply(seq_len(n_subjects), function(i) {
            .compute_gradient_joint(
              res$par,
              data_list = data_list[i],
              posteriors = final_posteriors[i],
              parameters = parameters,
              parallel = control$parallel,
              n_cores = control$n_cores,
              return_individual = TRUE
            )
          })
          S <- Reduce("+", lapply(score_list, function(s) tcrossprod(s)))
          s_inv <- solve(S)
          V <- s_inv %*% H %*% s_inv
        } else {
          V <- solve(H)
        }

        # Set parameter names
        dimnames(V) <- list(coef_names_expanded, coef_names_expanded)

        if (control$verbose > 0) {
          cli::cli_alert_success("Variance-covariance matrix computed")
        }
        V
      },
      error = function(e) {
        if (control$verbose > 0) {
          cli::cli_alert_warning(
            "Variance-covariance matrix computation failed: {e$message}"
          )
          cli::cli_alert_info(
            "Standard errors will not be available in summary output"
          )
        }
        NULL
      }
    )
  }

  # Compute model fit statistics
  log_lik <- -res$value
  n_params <- length(res$par) + 2
  aic <- -2 * log_lik + 2 * n_params
  bic <- -2 * log_lik + n_params * log(n_subjects)

  # Compute C-index for model discrimination using log hazard
  risk_scores <- sapply(seq_len(n_subjects), function(i) {
    ode_sol <- final_posteriors$ode_solutions[[i]]
    ode_sol$log_hazard_at_event + final_posteriors$b[i]
  })
  event_times <- vapply(data_list, `[[`, numeric(1), "time")
  event_status <- vapply(data_list, `[[`, numeric(1), "status")

  # Create a data frame for concordance calculation
  conc_data <- data.frame(
    time = event_times,
    status = event_status,
    risk = risk_scores
  )

  cindex <- survival::concordance(
    Surv(time, status) ~ risk,
    data = conc_data,
    reverse = TRUE # Higher risk should correspond to earlier events
  )$concordance

  if (control$verbose > 0) {
    cli::cli_alert_info(
      "C-index (concordance): {sprintf('%.3f', cindex)}"
    )
  }

  # Return fitted model
  structure(
    list(
      parameters = parameters,
      logLik = log_lik,
      AIC = aic,
      BIC = bic,
      cindex = cindex,
      convergence = list(
        converged = converged,
        em_iterations = if (converged) em_iter else control$maxit,
        message = sprintf(
          "EM algorithm %s after %d iterations",
          if (converged) "converged" else "did not converge within",
          if (converged) em_iter else control$maxit
        )
      ),
      random_effects = list(
        estimates = final_posteriors$b,
        variances = final_posteriors$v
      ),
      vcov = vcov_matrix,
      data = data_list,
      control = control,
      robust = robust,
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
  n_longitudinal <- length(object$parameters$coefficients$acceleration)

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
    "acceleration:",
    "",
    rownames(coef_longitudinal)
  )

  # Delta method for derived parameters (period and xi)
  derived_params <- NULL
  if (!is.null(object$vcov) && n_longitudinal >= 2) {
    # Extract coefficients from the acceleration equation:
    # acceleration = β₁ * biomarker + β₂ * velocity + ...
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

  # Count longitudinal observations and events
  n_observations <- attr(object$data, "n_observations")
  n_subjects <- attr(object$data, "n_subjects")
  event_rate <- attr(object$data, "event_rate")
  n_events <- n_subjects * event_rate

  structure(
    list(
      call = object$call,
      coef_baseline = coef_baseline,
      coef_survival = coef_survival,
      coef_longitudinal = coef_longitudinal,
      derived_params = derived_params,
      sigma = with(
        object$parameters$coefficients,
        c(sigma_e = measurement_error_sd, sigma_b = random_effect_sd)
      ),
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
  cat(sprintf("%20s\n", "StdDev"))
  cat(sprintf("Random Effect %14.6f\n", x$sigma["sigma_b"]))
  cat(sprintf("Residual %19.6f\n", x$sigma["sigma_e"]))

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
  coefs <- c(cf$baseline, cf$hazard, cf$acceleration)
  names(coefs) <- c(
    paste0("baseline:", names(cf$baseline)),
    paste0("hazard:", names(cf$hazard)),
    paste0("acceleration:", names(cf$acceleration))
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
    df = length(coef(object)) + 2,
    nobs = attr(object$data, "n_subjects"),
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
  cat(
    "\nLog-likelihood:",
    format(x$logLik, digits = digits),
    "on",
    length(coef(x)) + 2,
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

#' Predict Method for JointODE Models
#'
#' @description
#' Computes predictions for biomarker trajectories, velocities, accelerations,
#' and survival functions from a fitted JointODE model. This method uses the
#' fitted model parameters and random effects to solve the ODE system and
#' generate predictions for each subject.
#'
#' @param object A fitted JointODE model object.
#' @param times Optional numeric vector of time points at which to evaluate
#'   predictions. If NULL, uses the union of observation times and event time
#'   for each subject.
#' @param parallel Logical flag for parallel computation. Default is FALSE.
#' @param n_cores Number of CPU cores for parallel processing. If NULL,
#'   uses all available cores minus one.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing predictions for each subject, where each
#'   element contains:
#'   \describe{
#'     \item{\code{id}}{Subject identifier}
#'     \item{\code{times}}{Time points for predictions}
#'     \item{\code{biomarker}}{Predicted biomarker values with random effects}
#'     \item{\code{velocity}}{Predicted biomarker velocity}
#'     \item{\code{acceleration}}{Predicted biomarker acceleration}
#'     \item{\code{cumhazard}}{Cumulative hazard values}
#'     \item{\code{survival}}{Survival probability values}
#'   }
#'
#' @details
#' The prediction process involves:
#' \enumerate{
#'   \item Using the fitted model parameters and random effects
#'   \item Solving the ODE system for each subject at specified time points
#'   \item Computing biomarker trajectories with random effects included
#'   \item Calculating cumulative hazard and survival probabilities
#' }
#'
#' Random effects from the fitted model are always included in predictions
#' to provide subject-specific trajectories.
#' @concept model-prediction
#' @importFrom stats predict
#' @export
predict.JointODE <- function(
  object,
  times = NULL,
  parallel = FALSE,
  n_cores = 0,
  ...
) {
  # Extract model components
  parameters <- object$parameters
  data_list <- object$data
  random_effects <- object$random_effects$estimates

  # Function to compute predictions for a single subject
  compute_subject_predictions <- function(i) {
    subject_data <- data_list[[i]]
    b_i <- random_effects[i]

    # Determine time points for prediction
    if (!is.null(times)) {
      pred_times <- sort(unique(times))
    } else {
      # Use observation times (without event time for cleaner output)
      pred_times <- sort(unique(
        subject_data$longitudinal$times,
        subject_data$time
      ))
    }

    # Solve ODE with specified times
    ode_solution <- .solve_joint_ode(
      subject_data,
      parameters,
      sensitivity_type = "basic",
      times = pred_times
    )

    cumhazard_values <- ode_solution$cum_hazard * exp(b_i)
    survival_values <- exp(-cumhazard_values)

    predictions <- list(
      id = subject_data$id,
      times = pred_times,
      biomarker = ode_solution$biomarker + b_i,
      velocity = ode_solution$velocity,
      acceleration = ode_solution$acceleration,
      survival = survival_values
    )

    predictions
  }

  # Compute predictions for all subjects
  n_subjects <- length(data_list)

  if (parallel) {
    cleanup <- .setup_parallel_plan(n_cores)
    on.exit(cleanup(), add = TRUE)

    prediction_results <- future.apply::future_lapply(
      seq_len(n_subjects),
      compute_subject_predictions,
      future.seed = TRUE
    )
  } else {
    prediction_results <- lapply(
      seq_len(n_subjects),
      compute_subject_predictions
    )
  }

  # Return predictions
  prediction_results
}
