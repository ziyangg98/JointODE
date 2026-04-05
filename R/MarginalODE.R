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

  stop("MarginalODE is not yet available in the TMB version. ",
       "Use init = 'default' in JointODE().", call. = FALSE)
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
