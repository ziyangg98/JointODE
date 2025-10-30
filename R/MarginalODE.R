# Marginal ODE Utility Functions

#' @importFrom stats na.omit nlm setNames
NULL

#' Marginal Second-Order ODE Parameter Estimation
#'
#' @description
#' Estimates population-level (marginal) ODE parameters for longitudinal
#' biomarker trajectories using the second-order differential equation with
#' covariate support:
#' \deqn{\ddot{m}(t) = \text{value} \cdot m(t) +
#'   \text{slope} \cdot \dot{m}(t) + X\beta}
#'
#' @param formula A formula object specifying the response variable and
#'   covariates (e.g., \code{biomarker ~ x1 + x2} or \code{biomarker ~ 1}
#'   for intercept-only)
#' @param data A data frame containing the longitudinal measurements
#' @param time Character string specifying the time variable name
#'   (default: \code{"time"})
#' @param id Character string specifying the subject identifier variable name
#'   (default: \code{"id"})
#' @param state Optional matrix of initial conditions with two columns:
#'   \itemize{
#'     \item Column 1: Initial biomarker values \eqn{m(0)}
#'     \item Column 2: Initial velocities \eqn{\dot{m}(0)}
#'   }
#'   Each row corresponds to one subject. If \code{NULL}, initial values are
#'   estimated from the data.
#' @param control A list of control parameters for optimization, or output
#'   from \code{\link{JointODE.control}}. Key parameters include:
#'   \describe{
#'     \item{\code{verbose}}{Verbosity level (default: FALSE)}
#'     \item{\code{parallel}}{Logical flag enabling parallel computation
#'       (default: FALSE)}
#'     \item{\code{n_cores}}{Number of CPU cores for parallel processing.
#'       If 0, automatically detects available cores (default: 0)}
#'   }
#'   See \code{\link{JointODE.control}} for complete details.
#'
#' @concept model-fitting
#'
#' @return A list (S3 class \code{MarginalODE}) containing:
#'   \describe{
#'     \item{\code{parameters}}{Named numeric vector of parameter estimates
#'       (value, slope, and covariate coefficients)}
#'     \item{\code{measurement_error_sd}}{Residual standard deviation}
#'     \item{\code{logLik}}{Log-likelihood value}
#'     \item{\code{AIC}}{Akaike Information Criterion}
#'     \item{\code{BIC}}{Bayesian Information Criterion}
#'     \item{\code{convergence}}{List with \code{converged} (logical),
#'       \code{iterations}, and \code{message}}
#'     \item{\code{vcov}}{Variance-covariance matrix from inverse Hessian
#'       (may be NA matrix if Hessian is singular)}
#'     \item{\code{data}}{Processed data list used for fitting}
#'     \item{\code{control}}{Control parameters used}
#'     \item{\code{call}}{Matched function call}
#'   }
#'
#' @examples
#' \dontrun{
#' # Generate simulated data
#' fit <- MarginalODE(
#'   formula = observed ~ x1 + x2,
#'   data = sim$data$longitudinal_data,
#'   state = as.matrix(sim$data$state)
#' )
#' }
#'
#' @export
MarginalODE <- function(
  formula,
  data,
  time = "time",
  id = "id",
  state = NULL,
  control = list()
) {
  data_list <- .process_long(formula, data, time, id, state)

  control <- if (is.null(control)) {
    JointODE.control()
  } else if (is.list(control)) {
    JointODE.control(.list = modifyList(JointODE.control(), control))
  } else {
    stop("control must be a list or NULL")
  }

  n_subjects <- attr(data_list, "n_subjects")
  n_covariates <- attr(data_list, "n_covariates")
  covariate_names <- attr(data_list, "covariate_names")
  n_params <- 2 + n_covariates
  param_names <- c("value", "slope", covariate_names)

  if (control$verbose > 0) {
    cli::cli_h2("Second-Order ODE Model Estimation")
    n_obs <- sum(vapply(data_list, function(s) length(s$response), integer(1)))
    cli::cli_text(
      "Data: {n_subjects} subjects, {n_obs} observations, {n_params} parameters"
    )
    if (n_covariates > 0) {
      cli::cli_text("Covariates: {paste(covariate_names, collapse = ', ')}")
    }
    cli::cli_text("")
  }

  # Objective function with gradient and Hessian for nlm
  nlm_func <- function(theta) {
    .compute_marginal_objective_cppad(
      theta,
      data_list,
      gradient = TRUE,
      hessian = TRUE
    )
  }

  res <- nlm(
    nlm_func,
    rep(0, n_params),
    print.level = if (control$verbose > 1) 2 else 0,
    hessian = FALSE,
    check.analyticals = FALSE
  )

  if (control$verbose > 0) {
    if (res$code <= 2) {
      cli::cli_alert_success(
        "Optimization converged in {res$iterations} iterations"
      )
    } else {
      cli::cli_alert_warning(
        "Optimization did not converge (code {res$code})"
      )
    }
  }

  # Get Hessian at final parameters for vcov
  final_result <- .compute_marginal_objective_cppad(
    res$estimate,
    data_list,
    gradient = FALSE,
    hessian = TRUE
  )

  vcov_matrix <- tryCatch(
    {
      hess <- attr(final_result, "hessian")
      inv_h <- solve(hess)
      if (any(eigen(inv_h, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
        if (control$verbose > 0) {
          cli::cli_alert_warning("Hessian inverse not positive definite")
        }
        matrix(
          NA,
          n_params,
          n_params,
          dimnames = list(param_names, param_names)
        )
      } else {
        dimnames(inv_h) <- list(param_names, param_names)
        inv_h
      }
    },
    error = function(e) {
      if (control$verbose > 0) {
        cli::cli_alert_warning("Failed to compute vcov: {e$message}")
      }
      matrix(NA, n_params, n_params, dimnames = list(param_names, param_names))
    }
  )

  n_obs <- sum(vapply(data_list, function(s) length(s$response), integer(1)))
  residual_sd <- sqrt(res$minimum / (n_obs - n_params))
  loglik <- -0.5 * n_obs * (log(2 * pi) + log(res$minimum / n_obs) + 1)

  structure(
    list(
      parameters = setNames(res$estimate, param_names),
      measurement_error_sd = residual_sd,
      logLik = loglik,
      AIC = 2 * n_params - 2 * loglik,
      BIC = log(n_obs) * n_params - 2 * loglik,
      convergence = list(
        converged = res$code <= 2,
        iterations = res$iterations,
        message = if (res$code <= 2) "converged" else paste("code", res$code)
      ),
      vcov = vcov_matrix,
      data = data_list,
      control = control,
      call = match.call()
    ),
    class = "MarginalODE"
  )
}


#' Summary Method for MarginalODE Objects
#'
#' @param object An object of class \code{MarginalODE}
#' @param ... Additional arguments (currently unused)
#'
#' @return A list of class \code{summary.MarginalODE} containing:
#'   \item{call}{The model call}
#'   \item{coefficients}{Coefficient table for longitudinal ODE parameters}
#'   \item{sigma}{Named vector with measurement error SD}
#'   \item{nobs}{Number of subjects}
#'   \item{n_observations}{Total number of observations}
#'   \item{AIC}{Akaike Information Criterion}
#'   \item{BIC}{Bayesian Information Criterion}
#'   \item{logLik}{Log-likelihood}
#'   \item{convergence}{Convergence information}
#'
#' @concept model-summary
#' @export
summary.MarginalODE <- function(object, ...) {
  coef_vec <- object$parameters
  se_vec <- if (!is.null(object$vcov) && !all(is.na(object$vcov))) {
    sqrt(diag(object$vcov))
  } else {
    rep(NA_real_, length(coef_vec))
  }

  z_val <- coef_vec / se_vec
  p_val <- 2 * pnorm(-abs(z_val))

  structure(
    list(
      call = object$call,
      coefficients = cbind(
        Estimate = coef_vec,
        `Std. Error` = se_vec,
        `z value` = z_val,
        `Pr(>|z|)` = p_val
      ),
      sigma = c(sigma_e = object$measurement_error_sd),
      nobs = length(object$data),
      n_observations = sum(sapply(object$data, function(s) length(s$response))),
      AIC = object$AIC,
      BIC = object$BIC,
      logLik = object$logLik,
      convergence = object$convergence
    ),
    class = "summary.MarginalODE"
  )
}


#' Print Summary of MarginalODE Fit
#'
#' @param x An object of class \code{summary.MarginalODE}
#' @param digits Number of digits to display
#'   (default: max(3L, getOption("digits") - 3L))
#' @param signif.stars Logical; show significance stars
#'   (default: getOption("show.signif.stars"))
#' @param ... Additional arguments passed to \code{printCoefmat}
#'
#' @concept model-summary
#' @export
print.summary.MarginalODE <- function(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nData Descriptives:\n")
  cat(sprintf("Number of Observations: %d\n", x$n_observations))
  cat(sprintf("Number of Subjects: %d\n", x$nobs))

  cat(sprintf("\n%10s %10s %10s\n", "AIC", "BIC", "logLik"))
  cat(sprintf("%10.3f %10.3f %10.3f\n", x$AIC, x$BIC, x$logLik))

  cat("\nCoefficients:\n")
  cat("Second-Order ODE Model: Second-Order Dynamics\n")
  if (!is.null(x$coefficients)) {
    printCoefmat(
      x$coefficients,
      digits = digits,
      signif.stars = signif.stars,
      ...
    )
  }

  cat("\nVariance Components:\n")
  cat(sprintf("Measurement Error SD: %.6f\n", x$sigma["sigma_e"]))

  cat("\nConvergence:", ifelse(x$convergence$converged, "Yes", "No"))
  if (!x$convergence$converged) {
    cat("\n  Message:", x$convergence$message, "\n")
  } else {
    cat(sprintf(" (%d iterations)\n", x$convergence$iterations))
  }

  invisible(x)
}


#' Predict Method for MarginalODE Objects
#'
#' @description
#' Computes predicted biomarker trajectories, velocities, and accelerations
#' for subjects based on the fitted marginal ODE model.
#'
#' @param object An object of class \code{MarginalODE}
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
#'     \item{\code{biomarker}}{Predicted biomarker values}
#'     \item{\code{velocity}}{Predicted velocity (first derivative)}
#'     \item{\code{acceleration}}{Predicted acceleration (second derivative)}
#'   }
#'
#' @concept model-prediction
#' @export
predict.MarginalODE <- function(
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

  theta <- object$parameters
  control <- object$control
  data_list <- object$data

  compute_pred <- function(i) {
    subj <- data_list[[i]]
    subj_id <- names(data_list)[i]

    pred_times <- if (!is.null(times)) {
      time_vec <- if (is.list(times)) {
        if (!is.null(times[[subj_id]])) times[[subj_id]] else subj$time
      } else {
        times
      }
      sort(unique(time_vec))
    } else {
      subj$time
    }

    n_cov <- ncol(subj$covariates)
    pred_cov <- if (n_cov > 0) {
      cov_mat <- t(sapply(pred_times, function(t) {
        subj$covariates[which.min(abs(subj$time - t)), ]
      }))
      colnames(cov_mat) <- colnames(subj$covariates)
      cov_mat
    } else {
      matrix(0, length(pred_times), 0)
    }

    sol <- .solve_marginal_ode_cppad(
      theta,
      subj$initial,
      pred_times,
      pred_cov
    )

    data.frame(
      id = subj_id,
      time = pred_times,
      biomarker = sol$biomarker,
      velocity = sol$velocity,
      acceleration = sol$acceleration,
      stringsAsFactors = FALSE
    )
  }

  results <- .parallel_apply(
    seq_along(data_list),
    compute_pred,
    parallel = parallel,
    n_cores = n_cores
  )

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  result_df
}


# ===== SECTION 1: DATA PROCESSING =====

#' @importFrom stats model.frame model.matrix model.response
#' @noRd
.process_long <- function(formula, data, time, id, state) {
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  stopifnot(
    "Data cannot be empty" = nrow(data) > 0,
    "Data must contain a time column" = time %in% names(data),
    "Data must contain an id column" = id %in% names(data)
  )

  # Extract response and covariates using model.frame
  mf <- model.frame(formula, data = data, na.action = na.omit)
  y <- model.response(mf)
  X <- model.matrix(formula, data = mf)

  stopifnot(
    "Formula must include a response variable" = !is.null(y)
  )

  # Get corresponding time and id after na.omit
  row_idx <- as.numeric(rownames(mf))
  times <- data[[time]][row_idx]
  ids <- data[[id]][row_idx]
  subjects <- unique(ids)
  n_subjects <- length(subjects)

  if (!is.null(state)) {
    stopifnot(
      "state must be a matrix with 2 columns [m(0), m'(0)]" = is.matrix(
        state
      ) &&
        ncol(state) == 2,
      "state must have one row per subject" = nrow(state) == n_subjects
    )
  }

  # Process data for each subject
  subject_data <- lapply(seq_along(subjects), function(i) {
    idx <- which(ids == subjects[i])
    idx <- idx[order(times[idx])]

    if (length(idx) == 0) {
      stop(sprintf("No data for subject %s", subjects[i]))
    }

    y_subj <- y[idx]
    t_subj <- times[idx]
    x_subj <- X[idx, , drop = FALSE]

    if (length(t_subj) == 1 && t_subj[1] == 0) {
      warning(sprintf(
        "Subject %s: only one observation at t=0, skipping",
        subjects[i]
      ))
      return(NULL)
    }

    initial_state <- if (!is.null(state)) {
      c(state[i, 1], state[i, 2])
    } else {
      dt <- if (length(y_subj) > 1) t_subj[2] - t_subj[1] else 0
      c(
        y_subj[1],
        if (abs(dt) > .Machine$double.eps) (y_subj[2] - y_subj[1]) / dt else 0
      )
    }

    list(
      subject = subjects[i],
      response = y_subj,
      time = t_subj,
      covariates = x_subj,
      initial = initial_state
    )
  })

  # Remove NULL entries (subjects with insufficient data)
  subject_data <- Filter(Negate(is.null), subject_data)
  names(subject_data) <- as.character(vapply(
    subject_data,
    function(x) x$subject,
    FUN.VALUE = if (is.numeric(subjects[1])) numeric(1) else character(1)
  ))

  # Store metadata
  attr(subject_data, "n_subjects") <- length(subject_data)
  attr(subject_data, "n_covariates") <- ncol(X)
  attr(subject_data, "covariate_names") <- colnames(X)

  subject_data
}
