#' @importFrom stats na.omit nlm setNames predict
NULL

#' Marginal Second-Order ODE Parameter Estimation
#'
#' @description
#' Estimates the longitudinal part of a second-order ODE model:
#' \deqn{\ddot{m}_i(t) + 2\xi_i\omega_i\dot{m}_i(t) +
#'   \omega_i^2m_i(t) = f_i(t).}
#' The reserved formula terms \code{biomarker} and \code{velocity} represent
#' the latent log-coefficients \eqn{\log\omega_i^2} and
#' \eqn{\log(2\xi_i\omega_i)}, respectively.
#'
#' @param formula Longitudinal formula. The left-hand side is the observed
#'   response. On the right-hand side, \code{biomarker} and \code{velocity}
#'   activate the ODE dynamic parameters; other terms enter the forcing
#'   function. Random effects use lme-style syntax: \code{|} for full
#'   covariance and \code{||} for diagonal covariance.
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
#' data(sim)
#' fit <- MarginalODE(
#'   formula = observed ~ biomarker + velocity + x1 + x2 +
#'     (1 + biomarker + velocity || id),
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

  # Parse formula
  parsed_long <- .parse_longitudinal_formula(formula)
  if (is.null(parsed_long$grouping)) parsed_long$grouping <- id

  # Process data
  data_list <- .process_marginal(formula, data, time, id, parsed_long)
  scales <- list(time = .estimate_time_scale(data_list))
  data_fit <- lapply(data_list, function(d) {
    d$longitudinal$times <- d$longitudinal$times / scales$time
    d
  })
  n_subjects <- length(data_list)

  # Model setup
  model_config <- .setup_marginal_model(data, parsed_long)
  coef_names <- model_config$coef_names
  n_re <- model_config$n_re

  # Build parameters
  parameters <- .default_marginal_parameters(model_config, parsed_long)
  all_y <- unlist(lapply(data_list, function(d) d$longitudinal$measurements))
  mean_y <- mean(all_y)
  parameters$coefficients$initial_state <- setNames(c(mean_y, 0), .init_state_names)
  parameters$coefficients$measurement_error_sd <- sd(all_y)
  random_effects_init <- matrix(0, nrow = n_subjects, ncol = n_re)
  for (i in seq_along(data_list)) {
    obs <- data_list[[i]]$longitudinal
    if (length(obs$measurements) >= 1) {
      random_effects_init[i, 1] <- obs$measurements[1] - mean_y
    }
    if (length(obs$measurements) >= 2) {
      dt <- obs$times[2] - obs$times[1]
      if (is.finite(dt) && dt > 0) {
        random_effects_init[i, 2] <- (obs$measurements[2] - obs$measurements[1]) / dt
      }
    }
  }

  re_sds <- pmax(apply(random_effects_init, 2, sd), sd(all_y) * 0.1)
  parameters$coefficients$random_effect_sigma <- diag(re_sds^2, n_re)
  parameters$random_effects_init <- random_effects_init
  idx <- 1L
  if (parsed_long$biomarker$fixed) {
    parameters$coefficients$longitudinal[idx] <- -2 * log(scales$time)
    idx <- idx + 1L
  }
  if (parsed_long$velocity$fixed) {
    parameters$coefficients$longitudinal[idx] <- -log(scales$time)
  }
  parameters_fit <- parameters
  idx <- 1L
  if (parsed_long$biomarker$fixed) {
    parameters_fit$coefficients$longitudinal[idx] <-
      parameters_fit$coefficients$longitudinal[idx] + 2 * log(scales$time)
    idx <- idx + 1L
  }
  if (parsed_long$velocity$fixed) {
    parameters_fit$coefficients$longitudinal[idx] <-
      parameters_fit$coefficients$longitudinal[idx] + log(scales$time)
    idx <- idx + 1L
  }
  if (idx <= length(parameters_fit$coefficients$longitudinal)) {
    parameters_fit$coefficients$longitudinal[
      idx:length(parameters_fit$coefficients$longitudinal)
    ] <- parameters_fit$coefficients$longitudinal[
      idx:length(parameters_fit$coefficients$longitudinal)
    ] * scales$time^2
  }
  parameters_fit$coefficients$initial_state[2] <-
    parameters_fit$coefficients$initial_state[2] * scales$time
  re_scale <- rep(1, ncol(parameters_fit$random_effects_init))
  re_scale[2] <- scales$time
  forcing_start <- 3L + sum(
    parsed_long$biomarker$random, parsed_long$velocity$random
  )
  if (forcing_start <= length(re_scale)) {
    re_scale[forcing_start:length(re_scale)] <- scales$time^2
  }
  parameters_fit$random_effects_init <- sweep(
    parameters_fit$random_effects_init, 2, re_scale, `*`
  )
  parameters_fit$coefficients$random_effect_sigma <-
    parameters_fit$coefficients$random_effect_sigma *
      outer(re_scale, re_scale)

  # Pack TMB inputs
  tmb_data <- .pack_marginal_data(data_fit, parsed_long, n_re)
  tmb_params <- .pack_marginal_params(parameters_fit)

  .setup_openmp(control)

  if (control$verbose > 0) {
    cli::cli_h2("Marginal ODE Model Estimation (TMB)")
  }

  fit <- .fit_tmb(
    tmb_data = tmb_data,
    tmb_params = tmb_params,
    control = control,
    map = .correlation_map(parsed_long, tmb_params)
  )
  obj <- fit$obj
  opt <- fit$opt

  results <- .finalize_marginal(obj, opt, coef_names, n_re, n_subjects)
  re_nms <- model_config$re_names
  dimnames(results$random_effect_sigma) <- list(re_nms, re_nms)
  colnames(results$random_effects) <- re_nms
  n_init <- length(coef_names$initial_state)
  n_long <- length(results$parameters) - n_init
  idx <- 1L
  if (parsed_long$biomarker$fixed) {
    results$parameters[idx] <- results$parameters[idx] - 2 * log(scales$time)
    idx <- idx + 1L
  }
  if (parsed_long$velocity$fixed) {
    results$parameters[idx] <- results$parameters[idx] - log(scales$time)
    idx <- idx + 1L
  }
  if (idx <= n_long) {
    results$parameters[idx:n_long] <- results$parameters[idx:n_long] /
      scales$time^2
  }
  results$parameters["initial_velocity"] <-
    results$parameters["initial_velocity"] / scales$time
  par_scale <- rep(1 / scales$time^2, n_long)
  idx <- 1L
  if (parsed_long$biomarker$fixed) {
    par_scale[idx] <- 1
    idx <- idx + 1L
  }
  if (parsed_long$velocity$fixed) par_scale[idx] <- 1
  vcov_scale <- c(par_scale, 1, 1 / scales$time)
  results$vcov <- results$vcov * outer(vcov_scale, vcov_scale)

  re_scale <- rep(1, ncol(results$random_effects))
  re_scale[2] <- 1 / scales$time
  forcing_start <- 3L + sum(
    parsed_long$biomarker$random, parsed_long$velocity$random
  )
  if (forcing_start <= length(re_scale)) {
    re_scale[forcing_start:length(re_scale)] <- 1 / scales$time^2
  }
  results$random_effects <- sweep(results$random_effects, 2, re_scale, `*`)
  results$random_effect_sigma <- results$random_effect_sigma *
    outer(re_scale, re_scale)

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
    covariance = if (isTRUE(parsed_long$diagonal)) "diagonal" else "full"
  )
  coefs <- list(
    longitudinal = results$parameters[coef_names$longitudinal],
    initial_state = results$parameters[coef_names$initial_state]
  )

  structure(c(results, list(
    data = data_list, control = control, call = cl,
    tmb_obj = obj, tmb_opt = opt,
    configs = configs, coefs = coefs, scales = scales
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
#' @method summary MarginalODE
#' @export
summary.MarginalODE <- function(object, ...) {
  se <- if (!is.null(object$vcov) && !all(is.na(object$vcov))) {
    sqrt(pmax(diag(object$vcov), 0))
  } else {
    rep(NA_real_, length(object$parameters))
  }

  # Derived ODE physical parameters
  derived_params <- NULL
  nms <- names(object$parameters)
  has_dynamics <- all(c("log_omega2", "log_2xi_omega") %in% nms)
  if (!is.null(object$vcov) && !all(is.na(object$vcov)) && has_dynamics) {
    idx_o <- which(nms == "log_omega2")
    idx_d <- which(nms == "log_2xi_omega")
    log_omega2 <- object$parameters[idx_o]
    log_2xi_omega <- object$parameters[idx_d]
    var_o <- object$vcov[idx_o, idx_o]
    var_d <- object$vcov[idx_d, idx_d]
    cov_od <- object$vcov[idx_o, idx_d]

    omega_est <- sqrt(exp(log_omega2))
    xi_est <- exp(log_2xi_omega) / (2 * omega_est)

    vc_dyn <- matrix(c(var_o, cov_od, cov_od, var_d), 2, 2)
    grad_omega <- c(0.5 * omega_est, 0)
    grad_xi <- c(-0.5 * xi_est, xi_est)
    se_omega <- sqrt(pmax(drop(t(grad_omega) %*% vc_dyn %*% grad_omega), 0))
    se_xi <- sqrt(pmax(drop(t(grad_xi) %*% vc_dyn %*% grad_xi), 0))
    estimates <- c(omega_est, xi_est)
    ses <- c(se_omega, se_xi)
    z_values <- estimates / ses

    derived_params <- cbind(
      Estimate = estimates,
      `Std. Error` = ses,
      `z value` = z_values,
      `Pr(>|z|)` = 2 * pnorm(-abs(z_values))
    )
    rownames(derived_params) <- c(
      "omega (natural frequency)",
      "xi (damping ratio)"
    )
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
    .print_coefmat(x$derived_params,
      digits = digits,
      signif.stars = signif.stars, ...
    )
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
#' @param ... Additional arguments
#' @return A data.frame with id, time, biomarker, velocity, acceleration
#' @concept model-prediction
#' @export
predict.MarginalODE <- function(object, newdata = NULL, times = NULL, ...) {
  if (!is.null(newdata)) stop("newdata not yet supported")

  data_list <- object$data

  if (!is.null(times)) {
    .predict_marginal_trajectories(
      data_list, times, object$coefs, object$configs,
      object$random_effects
    )
  } else {
    # Per-subject at own observation times
    do.call(rbind, lapply(seq_along(data_list), function(i) {
      obs_t <- data_list[[i]]$longitudinal$times
      .predict_marginal_trajectories(
        data_list[i], obs_t, object$coefs, object$configs,
        object$random_effects[i, , drop = FALSE]
      )
    }))
  }
}
