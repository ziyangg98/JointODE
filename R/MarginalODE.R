#' @importFrom stats na.omit nlm setNames predict
NULL

#' Marginal Second-Order ODE Parameter Estimation
#'
#' @description
#' Estimates population-level ODE parameters for longitudinal
#' biomarker trajectories:
#' \eqn{\ddot{m}(t) = \beta_1 m(t) + \beta_2 \dot{m}(t) + X\beta}
#'
#' @param formula Response and covariates
#'   (e.g., \code{biomarker ~ x1 + x2})
#' @param data Data frame with longitudinal measurements
#' @param time Time variable name (default: \code{"time"})
#' @param id Subject identifier name (default: \code{"id"})
#' @param control List of control parameters.
#'   See \code{\link{MarginalODE.control}}.
#'
#' @return S3 object of class \code{MarginalODE}
#'
#' @examples
#' \dontrun{
#' fit <- MarginalODE(
#'   formula = observed ~ x1 + x2,
#'   data = sim$data$longitudinal_data
#' )
#' }
#'
#' @concept model-fitting
#' @export
# nolint next: object_name_linter
MarginalODE <- function(
  formula, data,
  time = "time", id = "id",
  control = list()
) {
  cl <- match.call()
  .validate_marginal(formula, data, time, id)

  if (is.null(control)) {
    control <- MarginalODE.control()
  } else if (is.list(control)) {
    control <- MarginalODE.control(.list = control)
  } else {
    stop("control must be a list or NULL")
  }

  # Parse formula (same as JointODE, supports (... | id) RE syntax)
  parsed_long <- .parse_longitudinal_formula(formula)
  if (is.null(parsed_long$grouping)) parsed_long$grouping <- id

  # Process data
  data_list <- .process_marginal(formula, data, time, id, parsed_long)
  n_subjects <- length(data_list)

  # Model setup
  model_config <- .setup_marginal_model(data, parsed_long)
  coef_names <- model_config$coef_names
  n_re <- model_config$n_re

  # Build parameters (default or user-provided)
  parameters <- .default_marginal_parameters(model_config, parsed_long)

  # Initialize RE from observations
  all_y <- unlist(lapply(data_list, function(d) d$longitudinal$measurements))
  mean_y <- mean(all_y)
  parameters$coefficients$initial_state <- c(mean_y, 0)

  random_effects_init <- matrix(0, nrow = n_subjects, ncol = n_re)
  for (i in seq_along(data_list)) {
    obs <- data_list[[i]]$longitudinal
    if (length(obs$measurements) >= 1)
      random_effects_init[i, 1] <- obs$measurements[1] - mean_y
    if (length(obs$measurements) >= 2) {
      dt <- obs$times[2] - obs$times[1]
      if (dt > 0)
        random_effects_init[i, 2] <- (obs$measurements[2] - obs$measurements[1]) / dt
    }
  }
  re_sds <- pmax(apply(random_effects_init, 2, sd), 1e-4)
  parameters$coefficients$measurement_error_sd <- sd(all_y)
  parameters$coefficients$random_effect_sigma <- diag(re_sds^2, n_re)
  parameters$random_effects_init <- random_effects_init

  parameters$configurations$biomarker_clamp <- max(abs(all_y)) * 5

  # Pack TMB inputs
  tmb_data <- .pack_marginal_data(data_list, parsed_long, n_re)
  tmb_params <- .pack_marginal_params(parameters)

  # Configure OpenMP
  if (control$parallel && control$n_cores > 0) {
    TMB::openmp(control$n_cores)
  }

  if (control$verbose > 0) cli::cli_h2("Marginal ODE Model Estimation (TMB)")

  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = "random_effects",
    DLL = "JointODE",
    silent = control$verbose < 2
  )

  opt <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(iter.max = control$maxit, eval.max = control$maxit * 2,
                   rel.tol = control$tol)
  )

  results <- .finalize_marginal(obj, opt, coef_names, n_re, n_subjects)

  if (control$verbose > 0) {
    if (results$convergence$converged) {
      cli::cli_alert_success(results$convergence$message)
    } else {
      cli::cli_alert_warning(results$convergence$message)
    }
    cli::cli_alert_info(sprintf("Log-likelihood: %.2f", results$logLik))
  }

  structure(c(results, list(
    data = data_list, control = control, call = cl,
    tmb_obj = obj, tmb_opt = opt
  )), class = "MarginalODE")
}

# -- S3 methods ---------------------------------------------------------------

#' Print MarginalODE Model
#' @param x A MarginalODE object
#' @param digits Number of digits for numeric output
#' @param ... Additional arguments
#' @return Invisibly returns the object
#' @concept model-display
#' @export
print.MarginalODE <- function(
  x, digits = max(3L, getOption("digits") - 3L), ...
) {
  cat("\nMarginal Second-Order ODE Model\n")
  cat("Call: ")
  print(x$call)
  cat(
    "\nLog-likelihood:", format(x$logLik, digits = digits),
    "on", length(x$parameters), "degrees of freedom\n"
  )
  cat(
    "AIC:", format(x$AIC, digits = digits),
    "  BIC:", format(x$BIC, digits = digits), "\n"
  )
  invisible(x)
}

#' Extract Model Coefficients
#' @param object A MarginalODE object
#' @param ... Additional arguments
#' @return Named numeric vector of parameter estimates
#' @concept model-inspection
#' @export
coef.MarginalODE <- function(object, ...) object$parameters

#' Extract Variance-Covariance Matrix
#' @param object A MarginalODE object
#' @param ... Additional arguments
#' @return Variance-covariance matrix
#' @concept model-inspection
#' @export
vcov.MarginalODE <- function(object, ...) object$vcov

#' Extract Log-Likelihood
#' @param object A MarginalODE object
#' @param ... Additional arguments
#' @return Log-likelihood with df and nobs attributes
#' @concept model-inspection
#' @export
logLik.MarginalODE <- function(object, ...) {
  structure(
    object$logLik,
    df = length(object$parameters),
    nobs = .n_obs(object$data),
    class = "logLik"
  )
}

#' Summary of MarginalODE Model
#' @param object A MarginalODE object
#' @param ... Additional arguments
#' @return A summary.MarginalODE object
#' @concept model-summary
#' @export
summary.MarginalODE <- function(object, ...) {
  se <- if (!is.null(object$vcov) && !all(is.na(object$vcov))) {
    sqrt(diag(object$vcov))
  } else {
    rep(NA_real_, length(object$parameters))
  }

  structure(
    list(
      call = object$call,
      coefficients = .coef_table(object$parameters, se),
      sigma = c(sigma_e = object$measurement_error_sd),
      nobs = length(object$data),
      n_observations = .n_obs(object$data),
      AIC = object$AIC, BIC = object$BIC,
      logLik = object$logLik,
      convergence = object$convergence
    ),
    class = "summary.MarginalODE"
  )
}

#' Print Summary of MarginalODE Model
#' @param x A summary.MarginalODE object
#' @param digits Number of digits to display
#' @param signif.stars Logical; show significance stars
#' @param ... Additional arguments passed to \code{printCoefmat}
#' @return Invisibly returns the object
#' @concept model-summary
#' @importFrom stats printCoefmat
#' @export
print.summary.MarginalODE <- function(
  x, digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"), ...
) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nData Descriptives:\n")
  cat(sprintf("Number of Observations: %d\n", x$n_observations))
  cat(sprintf("Number of Subjects: %d\n", x$nobs))

  cat(sprintf("\n%10s %10s %10s\n", "AIC", "BIC", "logLik"))
  cat(sprintf("%10.3f %10.3f %10.3f\n", x$AIC, x$BIC, x$logLik))
  cat("\nCoefficients:\n")
  .print_coefmat(
    x$coefficients,
    digits = digits, signif.stars = signif.stars, ...
  )
  cat(sprintf(
    "\nMeasurement Error SD: %.6f\n", x$sigma["sigma_e"]
  ))
  cat(sprintf("Convergence: %s\n", x$convergence$message))
  invisible(x)
}

#' Predict Method for MarginalODE Objects
#' @param object A MarginalODE object
#' @param newdata Not yet supported
#' @param times Prediction times (NULL = observed)
#' @param parallel Logical for parallel computation
#' @param n_cores Number of cores (0 = auto)
#' @param ... Additional arguments
#' @return A data.frame with id, time, biomarker, velocity, acceleration
#' @concept model-prediction
#' @export
predict.MarginalODE <- function(object, newdata = NULL, times = NULL, ...) {
  stop("MarginalODE predict is not yet available in the TMB version.",
       call. = FALSE)
}
