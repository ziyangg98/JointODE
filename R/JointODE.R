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
#' data(sim)
#' fit <- JointODE(
#'   longitudinal_formula = observed ~
#'     biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
#'   survival_formula = Surv(time, status) ~ w1 + w2,
#'   longitudinal_data = sim$data$longitudinal_data,
#'   survival_data = sim$data$survival_data,
#'   init = sim$init
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
  control = list()
) {
  cl <- match.call()

  # Parse formulas and config once
  parsed_long <- .parse_longitudinal_formula(longitudinal_formula)
  if (is.null(parsed_long$grouping)) parsed_long$grouping <- "id"
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
  scales <- list(time = .estimate_time_scale(data_list))
  data_fit <- lapply(data_list, function(d) {
    d$longitudinal$times <- d$longitudinal$times / scales$time
    d$time <- d$time / scales$time
    d
  })

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

  default_init <- FALSE
  if (is.character(init) && init == "marginal") {
    parameters <- .initialize_from_marginal(
      longitudinal_formula, longitudinal_data, survival_data,
      gamma, control, parsed_long, parsed_surv, model_config
    )
  } else if (is.list(init)) {
    parameters <- init
  } else {
    default_init <- TRUE
    parameters <- .default_parameters(
      model_config$dims, gamma, parsed_long,
      model_config$spline_baseline_config
    )
  }

  all_y <- unlist(lapply(
    data_list, function(d) d$longitudinal$measurements
  ))

  parameters$configurations$baseline <- model_config$spline_baseline_config
  parameters$configurations$hazard_quadrature <- control$hazard_quadrature
  parameters$configurations$gamma <- gamma
  parameters$configurations$biomarker <- list(
    fixed = parsed_long$biomarker$fixed,
    random = parsed_long$biomarker$random
  )
  parameters$configurations$velocity <- list(
    fixed = parsed_long$velocity$fixed,
    random = parsed_long$velocity$random
  )
  parameters$configurations$covariance <-
    if (isTRUE(parsed_long$diagonal)) "diagonal" else "full"

  # Data-driven defaults for initial state and sigma_e
  # (prevents inner Newton divergence when init = "default")
  cf <- parameters$coefficients
  if (is.null(cf$initial_state)) {
    cf$initial_state <- mean(all_y)
  }
  if (is.null(cf$measurement_error_sd)) {
    cf$measurement_error_sd <- sd(all_y)
  }
  if (default_init) {
    idx <- 1L
    if (parsed_long$biomarker$fixed) {
      cf$longitudinal[idx] <- -2 * log(scales$time)
      idx <- idx + 1L
    }
    if (parsed_long$velocity$fixed) {
      cf$longitudinal[idx] <- -log(scales$time)
    }
  }
  parameters$coefficients <- cf

  names(parameters$coefficients$baseline) <- model_config$coef_names$baseline
  names(parameters$coefficients$hazard) <- model_config$coef_names$hazard
  names(parameters$coefficients$longitudinal) <-
    model_config$coef_names$longitudinal
  names(parameters$coefficients$initial_state) <-
    model_config$coef_names$initial_state

  coef_names <- model_config$coef_names

  # Initialize random effects (zero if not set by marginal init)
  n_re <- ncol(model_config$random_effects)
  if (is.null(parameters$random_effects_init)) {
    parameters$random_effects_init <- matrix(
      0,
      nrow = length(data_list), ncol = n_re
    )
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
  parameters_fit$coefficients$baseline <-
    parameters_fit$coefficients$baseline + log(scales$time)
  if (length(parameters_fit$coefficients$hazard) >= 2) {
    parameters_fit$coefficients$hazard[2] <-
      parameters_fit$coefficients$hazard[2] / scales$time^gamma
  }
  re_scale <- rep(1, ncol(parameters_fit$random_effects_init))
  forcing_start <- 2L + sum(
    parsed_long$biomarker$random, parsed_long$velocity$random
  )
  if (forcing_start <= length(re_scale)) {
    re_scale[forcing_start:length(re_scale)] <- scales$time^2
  }
  parameters_fit$random_effects_init <- sweep(
    parameters_fit$random_effects_init, 2, re_scale, `*`
  )
  parameters_fit$coefficients$random_effect_sigma <-
    parameters_fit$coefficients$random_effect_sigma * outer(re_scale, re_scale)
  bc <- parameters_fit$configurations$baseline
  bc$knots <- bc$knots / scales$time
  bc$boundary_knots <- bc$boundary_knots / scales$time
  parameters_fit$configurations$baseline <- bc

  if (control$verbose > 0) {
    cli::cli_h2("Joint ODE Model Estimation (TMB)")
  }

  .setup_openmp(control)

  # Build TMB data and parameter lists
  tmb_data <- .pack_joint_data(data_fit, parameters_fit, control)
  tmb_params <- .pack_joint_params(parameters_fit)

  fit <- .fit_tmb(
    tmb_data = tmb_data,
    tmb_params = tmb_params,
    control = control,
    map = .correlation_map(parsed_long, tmb_params)
  )
  obj <- fit$obj
  opt <- fit$opt

  # Extract results
  results <- .finalize_joint(
    obj, opt, parameters_fit, coef_names,
    data_list, tmb_data$n_random_effects, control
  )
  re_nms <- model_config$re_names
  sigma_b <- results$parameters$coefficients$random_effect_sigma
  dimnames(sigma_b) <- list(re_nms, re_nms)
  results$parameters$coefficients$random_effect_sigma <- sigma_b
  colnames(results$random_effects) <- re_nms
  cf <- results$parameters$coefficients
  cf$baseline <- cf$baseline - log(scales$time)
  if (length(cf$hazard) >= 2) cf$hazard[2] <- cf$hazard[2] * scales$time^gamma
  idx <- 1L
  if (parsed_long$biomarker$fixed) {
    cf$longitudinal[idx] <- cf$longitudinal[idx] - 2 * log(scales$time)
    idx <- idx + 1L
  }
  if (parsed_long$velocity$fixed) {
    cf$longitudinal[idx] <- cf$longitudinal[idx] - log(scales$time)
    idx <- idx + 1L
  }
  if (idx <= length(cf$longitudinal)) {
    cf$longitudinal[idx:length(cf$longitudinal)] <-
      cf$longitudinal[idx:length(cf$longitudinal)] / scales$time^2
  }
  re_scale <- rep(1, ncol(results$random_effects))
  forcing_start <- 2L + sum(
    parsed_long$biomarker$random, parsed_long$velocity$random
  )
  if (forcing_start <= length(re_scale)) {
    re_scale[forcing_start:length(re_scale)] <- 1 / scales$time^2
  }
  results$random_effects <- sweep(results$random_effects, 2, re_scale, `*`)
  cf$random_effect_sigma <- cf$random_effect_sigma * outer(re_scale, re_scale)
  n_baseline <- length(cf$baseline)
  n_hazard <- length(cf$hazard)
  n_long <- length(cf$longitudinal)
  long_scale <- rep(1 / scales$time^2, n_long)
  idx <- 1L
  if (parsed_long$biomarker$fixed) {
    long_scale[idx] <- 1
    idx <- idx + 1L
  }
  if (parsed_long$velocity$fixed) long_scale[idx] <- 1
  vcov_scale <- c(
    rep(1, n_baseline),
    c(1, scales$time^gamma, rep(1, max(0, n_hazard - 2L))),
    long_scale,
    1
  )
  results$vcov <- results$vcov * outer(vcov_scale, vcov_scale)
  bc <- results$parameters$configurations$baseline
  bc$knots <- bc$knots * scales$time
  bc$boundary_knots <- bc$boundary_knots * scales$time
  results$parameters$configurations$baseline <- bc
  results$parameters$coefficients <- cf

  structure(c(results, list(
    data = data_list, control = control, call = cl,
    tmb_obj = obj, tmb_opt = opt,
    parsed_long = parsed_long, parsed_surv = parsed_surv,
    survival_formula = survival_formula, scales = scales
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
    sqrt(pmax(diag(object$vcov), 0))
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

  # Derived ODE physical parameters
  derived_params <- NULL
  long_names <- names(object$parameters$coefficients$longitudinal)
  has_roots <- all(c("log_omega2", "log_2xi_omega") %in% long_names)
  if (!is.null(object$vcov) && has_roots) {
    idx_o <- idx_longitudinal[which(long_names == "log_omega2")]
    idx_x <- idx_longitudinal[which(long_names == "log_2xi_omega")]
    log_omega2 <- coefs[idx_o]
    log_2xi_omega <- coefs[idx_x]
    var_o <- object$vcov[idx_o, idx_o]
    var_x <- object$vcov[idx_x, idx_x]
    cov_ox <- object$vcov[idx_o, idx_x]

    omega_est <- sqrt(exp(log_omega2))
    xi_est <- exp(log_2xi_omega) / (2 * omega_est)

    se_omega <- sqrt(pmax(0.25 * omega_est^2 * var_o, 0))
    grad_xi <- c(-0.5 * xi_est, xi_est)
    vc_damping <- matrix(c(var_o, cov_ox, cov_ox, var_x), 2, 2)
    se_xi <- sqrt(pmax(drop(t(grad_xi) %*% vc_damping %*% grad_xi), 0))
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
    .print_coefmat(x$coef_longitudinal,
      digits = digits,
      signif.stars = signif.stars, ...
    )
  }

  if (!is.null(x$derived_params)) {
    cat("\nODE System Characteristics:\n")
    .print_coefmat(x$derived_params,
      digits = digits,
      signif.stars = signif.stars, ...
    )
  }

  cat("\nSurvival Process: Proportional Hazards Model\n")
  if (!is.null(x$coef_survival)) {
    .print_coefmat(x$coef_survival,
      digits = digits,
      signif.stars = signif.stars, ...
    )
  }

  if (!is.null(x$coef_baseline)) {
    cat(
      "\nBaseline Hazard: B-spline with", nrow(x$coef_baseline),
      "basis functions\n"
    )
    cat(
      "(Coefficients range:",
      sprintf(
        "[%.3f, %.3f]",
        min(x$coef_baseline[, "Estimate"]),
        max(x$coef_baseline[, "Estimate"])
      ), ")\n"
    )
  }

  if (!is.null(x$coef_initial)) {
    cat("\nInitial State: Population Mean\n")
    .print_coefmat(x$coef_initial,
      digits = digits,
      signif.stars = signif.stars, ...
    )
  }

  cat("\nVariance Components:\n")
  cat(sprintf(
    "Measurement Error SD: %.6f\n",
    x$sigma["sigma_e"]
  ))
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
    list(
      baseline = cf$baseline, hazard = cf$hazard,
      longitudinal = cf$longitudinal, initial_state = cf$initial_state
    ),
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
  cat(
    "AIC:", format(x$AIC, digits = digits),
    "  BIC:", format(x$BIC, digits = digits), "\n"
  )
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
#' @param times Optional time points for prediction. If NULL, uses observed
#'   time points for each subject.
#' @param ... Additional arguments (currently unused)
#'
#' @return A data.frame with columns id, time, biomarker, velocity.
#'
#' @concept model-prediction
#' @export
predict.JointODE <- function(object, newdata = NULL, times = NULL, ...) {
  cf <- object$parameters$coefficients
  configs <- object$parameters$configurations
  re <- object$random_effects

  if (!is.null(newdata)) {
    if (!is.list(newdata) || is.null(newdata$longitudinal_data) ||
      is.null(newdata$survival_data)) {
      stop("newdata must be a list with $longitudinal_data and $survival_data")
    }
    data_list <- .process_joint(
      longitudinal_data = newdata$longitudinal_data,
      survival_data = newdata$survival_data,
      parsed_long = object$parsed_long,
      parsed_surv = object$parsed_surv,
      survival_formula = object$survival_formula
    )
    # Match RE by subject ID; zero for unseen subjects
    orig_ids <- names(object$data)
    new_ids <- names(data_list)
    n_re <- ncol(re)
    re_new <- matrix(0, nrow = length(new_ids), ncol = n_re)
    matched <- match(new_ids, orig_ids)
    has_match <- !is.na(matched)
    if (any(has_match)) re_new[has_match, ] <- re[matched[has_match], , drop = FALSE]
    re <- re_new
  } else {
    data_list <- object$data
  }

  n_quad <- object$control$hazard_quadrature

  if (!is.null(times)) {
    .predict_trajectories(data_list, times, cf, configs, re, n_quad)
  } else {
    # Per-subject at own observation times
    do.call(rbind, lapply(seq_along(data_list), function(i) {
      obs_t <- data_list[[i]]$longitudinal$times
      .predict_trajectories(
        data_list[i], obs_t, cf, configs,
        re[i, , drop = FALSE], n_quad
      )
    }))
  }
}
