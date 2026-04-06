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

  all_y <- unlist(lapply(data_list, function(d) d$longitudinal$measurements))
  mean_y <- mean(all_y)
  has_dynamics_re <- parsed_long$biomarker$random ||
    parsed_long$velocity$random

  if (has_dynamics_re) {
    # Phase 1: fit reduced model to warm-start dynamics RE
    ws <- .warm_start_marginal(
      formula, data, time, id, parsed_long, control
    )
    parameters$coefficients$longitudinal <- ws$longitudinal
    parameters$coefficients$initial_state <- ws$initial_state
    parameters$coefficients$measurement_error_sd <- ws$measurement_error_sd

    # Expand RE matrix: cols 1-2 from warm start, rest = 0
    random_effects_init <- matrix(0, n_subjects, n_re)
    random_effects_init[, 1:2] <- ws$random_effects
  } else {
    parameters$coefficients$initial_state <- c(mean_y, 0)
    parameters$coefficients$measurement_error_sd <- sd(all_y)
    random_effects_init <- .init_re_from_observations(
      data_list, mean_y, 0, n_re
    )
  }

  re_sds <- pmax(apply(random_effects_init, 2, sd), sd(all_y) * 0.1)
  parameters$coefficients$random_effect_sigma <- diag(re_sds^2, n_re)
  parameters$random_effects_init <- random_effects_init
  parameters$configurations$biomarker_clamp <- max(abs(all_y)) * 5

  # Pack TMB inputs
  tmb_data <- .pack_marginal_data(data_list, parsed_long, n_re)
  tmb_params <- .pack_marginal_params(parameters)

  .setup_openmp(control)

  if (control$verbose > 0) {
    phase_label <- if (has_dynamics_re) "Phase 2: " else ""
    cli::cli_h2(paste0(phase_label, "Marginal ODE Model Estimation (TMB)"))
  }

  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = "random_effects",
    DLL = "JointODE",
    silent = control$verbose < 3,
    normalize = TRUE,
    inner.control = list(ustep = 1)
  )

  opt <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(
      iter.max = control$maxit, eval.max = control$maxit * 2,
      rel.tol = control$tol,
      trace = as.integer(control$verbose >= 2)
    )
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

  # Store configs needed for predict
  configs <- list(
    biomarker = list(
      fixed = parsed_long$biomarker$fixed,
      random = parsed_long$biomarker$random
    ),
    velocity = list(
      fixed = parsed_long$velocity$fixed,
      random = parsed_long$velocity$random
    ),
    biomarker_clamp = parameters$configurations$biomarker_clamp
  )
  # Store longitudinal coefficients in named form
  coefs <- list(
    longitudinal = setNames(
      as.numeric(obj$env$last.par.best[
        names(obj$env$last.par.best) == "longitudinal"]),
      coef_names$longitudinal
    ),
    initial_state = setNames(
      as.numeric(obj$env$last.par.best[
        names(obj$env$last.par.best) == "initial_state"]),
      coef_names$initial_state
    )
  )

  structure(c(results, list(
    data = data_list, control = control, call = cl,
    tmb_obj = obj, tmb_opt = opt,
    configs = configs, coefs = coefs
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
coef.MarginalODE <- function(object, ...) {
  p <- object$parameters
  is_init <- names(p) %in% .init_state_names
  names(p)[!is_init] <- paste0("longitudinal:", names(p)[!is_init])
  names(p)[is_init] <- paste0("initial:", names(p)[is_init])
  p
}

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

  # Derived ODE physical parameters (period and damping ratio)
  derived_params <- NULL
  nms <- names(object$parameters)
  has_bv <- "biomarker" %in% nms && "velocity" %in% nms
  if (!is.null(object$vcov) && !all(is.na(object$vcov)) && has_bv) {
    idx_b <- which(nms == "biomarker")
    idx_v <- which(nms == "velocity")
    value_coef <- object$parameters[idx_b]
    slope_coef <- object$parameters[idx_v]
    var_value <- object$vcov[idx_b, idx_b]
    var_slope <- object$vcov[idx_v, idx_v]
    cov_value_slope <- object$vcov[idx_b, idx_v]

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

  structure(
    list(
      call = object$call,
      coefficients = .coef_table(object$parameters, se),
      derived_params = derived_params,
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

  if (!is.null(x$derived_params)) {
    cat("\nODE System Characteristics:\n")
    .print_coefmat(x$derived_params, digits = digits,
                   signif.stars = signif.stars, ...)
  }
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
  if (!is.null(newdata)) stop("newdata not yet supported")

  data_list <- object$data

  if (is.null(times)) {
    all_t <- unlist(lapply(data_list, function(d) d$longitudinal$times))
    times <- sort(unique(c(0, all_t)))
  }

  .predict_marginal_trajectories(
    data_list, times, object$coefs, object$configs, object$random_effects
  )
}
