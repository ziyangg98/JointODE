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
#' @param state Optional \eqn{n \times 2} matrix of initial
#'   conditions \eqn{[m(0), \dot{m}(0)]}. If \code{NULL},
#'   estimated from data.
#' @param control List of control parameters.
#'   See \code{\link{MarginalODE.control}}.
#'
#' @return S3 object of class \code{MarginalODE}
#'
#' @examples
#' \dontrun{
#' fit <- MarginalODE(
#'   formula = observed ~ x1 + x2,
#'   data = sim$data$longitudinal_data,
#'   state = as.matrix(sim$data$state)
#' )
#' }
#'
#' @concept model-fitting
#' @export
# nolint next: object_name_linter
MarginalODE <- function(
  formula, data,
  time = "time", id = "id",
  state = NULL, control = list()
) {
  cl <- match.call()
  .validate_marginal(formula, data, time, id, state)

  if (is.null(control)) {
    control <- MarginalODE.control()
  } else if (is.list(control)) {
    control <- MarginalODE.control(.list = control)
  } else {
    stop("control must be a list or NULL")
  }

  data_list <- .process_marginal(formula, data, time, id, state)
  n_covariates <- attr(data_list, "n_covariates")
  covariate_names <- attr(data_list, "covariate_names")
  biomarker_clamp <- attr(data_list, "biomarker_clamp")
  n_params <- 2 + n_covariates
  param_names <- c("value", "slope", covariate_names)
  has_state <- !is.null(state)

  if (control$verbose > 0) {
    cli::cli_h2("Marginal ODE Model Estimation")
    cli::cli_text(sprintf(
      "Data: %d subjects, %d obs, %d params",
      length(data_list), .n_obs(data_list), n_params
    ))
  }

  if (control$parallel) {
    cleanup <- .setup_parallel_plan(control$n_cores)
    on.exit(cleanup(), add = TRUE)
  }

  # Optimization loop (1 iteration when state provided)
  theta <- rep(0, n_params)
  sse <- Inf

  for (iter in seq_len(if (has_state) 1L else control$maxit)) {
    if (!has_state) {
      opt <- .parallel_apply(
        seq_along(data_list),
        function(i) {
          .estimate_marginal_state(
            data_list[[i]]$initial_state, data_list[[i]],
            theta, biomarker_clamp
          )
        },
        control$parallel, control$n_cores,
        setup = FALSE
      )
      for (i in seq_along(data_list)) {
        data_list[[i]]$initial_state <- opt[[i]]
      }
    }

    fit <- nlm(
      function(th) {
        .compute_marginal_objective(
          th, data_list, biomarker_clamp, TRUE, TRUE
        )
      },
      theta,
      print.level = 0,
      hessian = FALSE, check.analyticals = FALSE
    )
    theta <- fit$estimate
    prev_sse <- sse
    sse <- fit$minimum

    if (control$verbose > 1) {
      cli::cli_alert_info(sprintf(
        "Iter %d: SSE=%.4f", iter, sse
      ))
    }

    rel <- abs(sse - prev_sse) / (abs(sse) + 1)
    if (iter > 1 && rel < control$tol) break
  }

  converged <- if (has_state) {
    fit$code <= 2
  } else {
    iter > 1 && rel < control$tol
  }
  n_iter <- if (has_state) fit$iterations else iter

  if (control$verbose > 0) {
    cli::cli_alert_info(sprintf(
      "%s after %d iterations",
      if (converged) "Converged" else "Did not converge",
      n_iter
    ))
  }

  .finalize_marginal(
    theta, sse, data_list, biomarker_clamp,
    param_names, converged, n_iter, has_state,
    control, cl
  )
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
predict.MarginalODE <- function(
  object, newdata = NULL, times = NULL,
  parallel = FALSE, n_cores = 0, ...
) {
  if (!is.null(newdata)) stop("newdata not yet supported")

  data_list <- object$data

  pred_data <- if (!is.null(times)) {
    lapply(seq_along(data_list), function(i) {
      subj <- data_list[[i]]
      obs_t <- subj$longitudinal$times
      pred_t <- if (is.list(times)) {
        tv <- times[[names(data_list)[i]]]
        if (!is.null(tv)) sort(unique(tv)) else obs_t
      } else {
        sort(unique(times))
      }
      subj$longitudinal$times <- pred_t
      subj$longitudinal$measurements <- rep(0, length(pred_t))
      subj$longitudinal$covariates$fixed <- .extend_covariates(
        subj$longitudinal$covariates$fixed, obs_t, pred_t
      )
      subj
    })
  } else {
    data_list
  }

  sols <- .solve_batch_marginal(
    pred_data, object$parameters,
    attr(object$data, "biomarker_clamp")
  )

  result <- do.call(rbind, lapply(
    seq_along(sols),
    function(i) {
      data.frame(
        id = names(data_list)[i],
        time = sols[[i]]$times,
        biomarker = sols[[i]]$biomarker,
        velocity = sols[[i]]$velocity,
        acceleration = sols[[i]]$acceleration,
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(result) <- NULL
  result
}
