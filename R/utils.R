# Utility Functions for JointODE Package

# Package Constants ============================================================

.default_spline <- list(
  degree = 3,
  n_knots = 5,
  knot_placement = "quantile",
  boundary_knots = NULL
)

.reserved_words <- c("biomarker", "velocity")

# Formatting ===================================================================

#' @importFrom utils head
#' @noRd
.format_vector <- function(x, n = 4) {
  if (!length(x)) {
    return("[]")
  }
  shown <- head(x, n)
  rest <- length(x) - n
  paste0(
    "[",
    paste(sprintf("%.3f", shown), collapse = ", "),
    if (rest > 0) paste0(", ...+", rest) else "",
    "]"
  )
}

# Parallel Computing ===========================================================

#' @importFrom parallel detectCores
#' @importFrom future plan multicore multisession sequential
#' @noRd
.setup_parallel_plan <- function(n_cores = 0) {
  if (n_cores == 0) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  if (.Platform$OS.type == "unix") {
    future::plan(future::multicore, workers = n_cores)
  } else {
    future::plan(future::multisession, workers = n_cores)
  }
  function() future::plan(future::sequential)
}

#' @importFrom future.apply future_lapply
#' @noRd
.parallel_apply <- function(indices, fn, parallel = TRUE, n_cores = 0) {
  if (parallel) {
    cleanup <- .setup_parallel_plan(n_cores)
    on.exit(cleanup(), add = TRUE)
    future.apply::future_lapply(indices, fn, future.seed = TRUE)
  } else {
    lapply(indices, fn)
  }
}

# Formula Parsing ==============================================================

#' @importFrom stats as.formula
#' @noRd
.build_formula <- function(terms, response = NULL, is_random = FALSE) {
  covars <- setdiff(terms, "(Intercept)")
  has_intercept <- "(Intercept)" %in% terms

  terms_str <- if (length(covars) > 0) {
    if (is_random) {
      covars <- gsub("^\\(Intercept\\)$", "1", covars)
    }
    paste(covars, collapse = " + ")
  } else if (has_intercept) {
    "1"
  } else {
    "0"
  }

  if (!has_intercept && length(covars) > 0 && !is_random) {
    terms_str <- paste("0 +", terms_str)
  }

  if (is_random || is.null(response)) {
    as.formula(paste("~", terms_str))
  } else {
    as.formula(paste(response, "~", terms_str))
  }
}

#' @importFrom stats as.formula terms
#' @noRd
.parse_longitudinal_formula <- function(formula) {
  formula_str <- paste(deparse(formula, width.cutoff = 500L), collapse = "")
  re_pattern <- "\\(([^)]+\\|[^)]+)\\)"
  re_match <- regmatches(formula_str, regexpr(re_pattern, formula_str))

  if (length(re_match) == 0) {
    fixed_formula <- formula
    random_terms <- NULL
    grouping <- NULL
    biomarker_random <- FALSE
    velocity_random <- FALSE
  } else {
    inner <- substr(re_match, 2, nchar(re_match) - 1)
    parts <- strsplit(inner, "\\|")[[1]]
    grouping <- trimws(parts[2])
    re_terms_str <- trimws(parts[1])

    random_terms <- if (re_terms_str == "1") {
      "(Intercept)"
    } else {
      terms <- trimws(strsplit(re_terms_str, "\\+")[[1]])
      terms[terms == "1"] <- "(Intercept)"
      terms
    }

    fixed_str <- sub(re_match, "", formula_str, fixed = TRUE)
    fixed_str <- gsub("\\s+\\+\\s*$", "", fixed_str)
    fixed_str <- gsub("~\\s*\\+\\s*", "~ ", fixed_str)
    fixed_str <- gsub("\\+\\s*\\+", "+", fixed_str)
    fixed_formula <- as.formula(fixed_str)

    biomarker_random <- "biomarker" %in% random_terms
    velocity_random <- "velocity" %in% random_terms
    random_terms <- setdiff(random_terms, .reserved_words)
    if (length(random_terms) == 0) random_terms <- NULL
  }

  fixed_terms_obj <- terms(fixed_formula)
  fixed_terms_labels <- attr(fixed_terms_obj, "term.labels")
  has_intercept <- attr(fixed_terms_obj, "intercept") == 1

  biomarker_in_fixed <- "biomarker" %in% fixed_terms_labels
  velocity_in_fixed <- "velocity" %in% fixed_terms_labels
  fixed_covariates <- setdiff(fixed_terms_labels, .reserved_words)

  if (has_intercept) {
    fixed_covariates <- c("(Intercept)", fixed_covariates)
  }

  list(
    response = as.character(formula[[2]]),
    fixed_terms = fixed_covariates,
    random_terms = random_terms,
    biomarker = list(fixed = biomarker_in_fixed, random = biomarker_random),
    velocity = list(fixed = velocity_in_fixed, random = velocity_random),
    grouping = grouping
  )
}

#' @importFrom stats terms
#' @noRd
.parse_survival_formula <- function(formula) {
  surv_response <- formula[[2]]

  if (
    !inherits(surv_response, "call") ||
      !identical(as.character(surv_response[[1]]), "Surv")
  ) {
    stop(
      "Survival formula must have Surv() on the left-hand side",
      call. = FALSE
    )
  }

  surv_call_vars <- all.vars(surv_response)
  if (length(surv_call_vars) < 2) {
    stop("Surv() must have at least time and status arguments", call. = FALSE)
  }

  time_var <- surv_call_vars[1]
  status_var <- surv_call_vars[2]

  all_terms <- attr(terms(formula), "term.labels")
  covariate_terms <- if (length(all_terms) == 0) NULL else all_terms

  list(
    time_var = time_var,
    status_var = status_var,
    covariate_terms = covariate_terms
  )
}

# Model Configuration ==========================================================

#' @noRd
.compute_dimensions <- function(
  longitudinal_formula,
  survival_formula,
  spline_baseline
) {
  parsed_long <- .parse_longitudinal_formula(longitudinal_formula)
  parsed_surv <- .parse_survival_formula(survival_formula)

  n_longitudinal_fixed <- length(parsed_long$fixed_terms)
  n_longitudinal_random <- if (is.null(parsed_long$random_terms)) {
    0
  } else {
    length(parsed_long$random_terms)
  }
  n_survival_covariates <- if (is.null(parsed_surv$covariate_terms)) {
    0
  } else {
    length(parsed_surv$covariate_terms)
  }

  n_biomarker_velocity_fixed <- sum(
    parsed_long$biomarker$fixed,
    parsed_long$velocity$fixed
  )
  n_biomarker_velocity_random <- sum(
    parsed_long$biomarker$random,
    parsed_long$velocity$random
  )

  n_longitudinal_coef <- n_longitudinal_fixed + n_biomarker_velocity_fixed
  n_random_effects <- n_longitudinal_random + n_biomarker_velocity_random

  spline_config <- modifyList(.default_spline, spline_baseline)
  n_spline_basis <- spline_config$degree + spline_config$n_knots + 1

  list(
    n_longitudinal_coef = n_longitudinal_coef,
    n_random_effects = n_random_effects,
    n_survival_covariates = n_survival_covariates,
    n_spline_basis = n_spline_basis,
    spline_config = spline_config
  )
}

#' @importFrom stats quantile
#' @noRd
.get_spline_config <- function(
  x,
  degree = 3,
  n_knots = 5,
  knot_placement = "quantile",
  boundary_knots = NULL
) {
  if (is.null(boundary_knots)) {
    boundary_knots <- range(x, na.rm = TRUE)
  }

  if (knot_placement == "quantile") {
    probs <- seq(0, 1, length.out = n_knots + 2)[-c(1, n_knots + 2)]
    knots <- quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  } else if (knot_placement == "equal") {
    knots <- seq(
      boundary_knots[1],
      boundary_knots[2],
      length.out = n_knots + 2
    )[-c(1, n_knots + 2)]
  } else {
    stop("knot_placement must be 'quantile' or 'equal'")
  }

  list(
    degree = degree,
    knots = knots,
    boundary_knots = boundary_knots,
    df = length(knots) + degree + 1
  )
}

# Variance-Covariance ==========================================================

#' @noRd
.compute_vcov_sem <- function(
  data_list,
  posteriors,
  parameters,
  random_effects,
  control,
  verbose = 0
) {
  theta <- .coef_to_vector(parameters)
  n_coef <- length(theta)

  if (verbose > 0) {
    cli::cli_alert_info("Computing SEM vcov...")
  }

  obj <- .compute_objective_expected(
    theta,
    data_list,
    posteriors,
    parameters,
    gradient = FALSE,
    hessian = TRUE
  )
  info_complete_inv <- solve(attr(obj, "hessian"))

  em_map <- function(theta_input) {
    params_input <- .vector_to_coef(parameters, theta_input)
    em_result <- .em_step(
      data_list,
      params_input,
      random_effects,
      control
    )
    .coef_to_vector(em_result$parameters)
  }

  eps_weight <- sqrt(diag(info_complete_inv))
  eps_weight <- eps_weight / min(eps_weight) * sqrt(.Machine$double.eps)

  dm_matrix <- numDeriv::jacobian(
    func = em_map,
    x = theta,
    method = "simple",
    method.args = list(eps = eps_weight)
  )

  info_observed_inv <- info_complete_inv %*% solve(diag(n_coef) - t(dm_matrix))

  asymmetry <- max(abs(info_observed_inv - t(info_observed_inv)))
  if (verbose > 0 && asymmetry > 1e-4) {
    cli::cli_alert_warning(sprintf(
      "Vcov matrix asymmetry: %.2e (numerical instability detected)",
      asymmetry
    ))
  }
  vcov_sym <- (info_observed_inv + t(info_observed_inv)) / 2

  eigenvalues <- eigen(vcov_sym, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigenvalues <= 0) && verbose > 0) {
    n_neg <- sum(eigenvalues <= 0)
    min_eigen <- min(eigenvalues)
    cli::cli_alert_warning(sprintf(
      "%d non-positive eigenvalue%s (min: %.2e)",
      n_neg,
      if (n_neg > 1) "s" else "",
      min_eigen
    ))
  }

  diag_ratio <- diag(vcov_sym) / diag(info_complete_inv)
  if (any(diag_ratio < 1 - 1e-6) && verbose > 0) {
    n_violate <- sum(diag_ratio < 1 - 1e-6)
    min_ratio <- min(diag_ratio)
    cli::cli_alert_warning(sprintf(
      "%d parameter%s violate missing information principle (min ratio: %.3f)",
      n_violate,
      if (n_violate > 1) "s" else "",
      min_ratio
    ))
  }

  vcov_sym
}

# Parameter Vector Conversion =================================================

#' @noRd
.coef_to_vector <- function(parameters) {
  with(parameters$coefficients, c(baseline, hazard, longitudinal))
}

#' @noRd
.vector_to_coef <- function(parameters, theta) {
  n_base <- length(parameters$coefficients$baseline)
  n_haz <- length(parameters$coefficients$hazard)
  idx <- cumsum(c(n_base, n_haz, length(theta) - n_base - n_haz))

  parameters$coefficients$baseline <- theta[1:idx[1]]
  parameters$coefficients$hazard <- theta[(idx[1] + 1):idx[2]]
  parameters$coefficients$longitudinal <- theta[(idx[2] + 1):idx[3]]
  parameters
}

# State Optimization ===========================================================

#' Compute state log-likelihood and derivatives
#'
#' @noRd
.compute_state_objective <- function(
  initial_state,
  data,
  random_effect,
  parameters
) {
  result <- .compute_state_loglik_cppad(
    initial_state = initial_state,
    data = data,
    random_effect = random_effect,
    parameters = parameters,
    gradient = TRUE,
    hessian = TRUE
  )

  # Return negative for minimization
  value <- -as.numeric(result)
  attr(value, "gradient") <- -as.vector(attr(result, "gradient"))
  attr(value, "hessian") <- -as.matrix(attr(result, "hessian"))
  value
}

#' Estimate initial state
#'
#' @noRd
.estimate_state <- function(
  initial_guess,
  data,
  random_effect,
  parameters,
  max_iter = 50,
  tol = 1e-6
) {
  objective <- function(state) {
    .compute_state_objective(state, data, random_effect, parameters)
  }

  fit <- suppressWarnings(
    nlm(
      f = objective,
      p = initial_guess,
      iterlim = max_iter,
      gradtol = tol,
      hessian = FALSE,
      check.analyticals = FALSE
    )
  )

  list(
    state = fit$estimate,
    converged = fit$code <= 2,
    iterations = fit$iterations
  )
}

# EM Algorithm =================================================================

#' @noRd
.em_step <- function(data_list, parameters, random_effects, control) {
  # E-step: Compute posteriors
  posteriors <- .compute_posteriors(
    data_list,
    parameters,
    random_effects,
    control$parallel,
    control$n_cores,
    control$quad_level
  )

  # Update random effects from posterior
  n_subjects <- length(data_list)
  posterior_moments <- lapply(seq_len(n_subjects), function(i) {
    .compute_posterior_moments(
      list(nodes = posteriors$nodes[[i]], weights = posteriors$weights[[i]])
    )
  })
  random_effects <- t(sapply(posterior_moments, `[[`, "mean"))

  random_effect_sigma <- Reduce(
    `+`,
    lapply(posterior_moments, `[[`, "second_moment")
  )
  random_effect_sigma <- random_effect_sigma / n_subjects
  measurement_error_sd <- .update_measurement_error_sd(
    data_list,
    parameters,
    posteriors
  )
  parameters$coefficients$measurement_error_sd <- measurement_error_sd
  parameters$coefficients$random_effect_sigma <- random_effect_sigma

  theta <- .coef_to_vector(parameters)

  obj <- .compute_objective_expected(
    theta,
    data_list,
    posteriors,
    parameters,
    gradient = TRUE,
    hessian = TRUE
  )

  gradient <- attr(obj, "gradient")
  hessian <- attr(obj, "hessian")

  theta_new <- theta - solve(hessian, gradient)

  rm(obj, gradient, hessian)

  obj_new <- .compute_objective_expected(
    theta_new,
    data_list,
    posteriors,
    parameters,
    gradient = FALSE,
    hessian = FALSE
  )

  loglik_value <- -as.numeric(obj_new)
  rm(obj_new)

  parameters <- .vector_to_coef(parameters, theta_new)

  rm(posteriors, posterior_moments, theta, theta_new)
  gc(verbose = FALSE, full = FALSE)

  list(
    parameters = parameters,
    random_effects = random_effects,
    loglik = loglik_value
  )
}

#' @noRd
.compute_metrics <- function(curr, prev, iter) {
  if (iter > 1) {
    delta_l <- curr$loglik - prev$loglik
    rel_l <- abs(delta_l) / (abs(curr$loglik) + 1e-8)
    delta_theta <- with(curr$parameters$coefficients, {
      max(abs(c(
        baseline - prev$parameters$coefficients$baseline,
        hazard - prev$parameters$coefficients$hazard,
        longitudinal - prev$parameters$coefficients$longitudinal,
        measurement_error_sd -
          prev$parameters$coefficients$measurement_error_sd,
        random_effect_sigma - prev$parameters$coefficients$random_effect_sigma
      )))
    })
  } else {
    delta_l <- curr$loglik
    rel_l <- 1
    delta_theta <- Inf
  }

  list(
    delta_l = delta_l,
    rel_l = rel_l,
    delta_theta = delta_theta
  )
}

# Progress Tracking ============================================================

#' @noRd
.print_iteration <- function(iter, curr, metrics, control) {
  cf <- curr$parameters$coefficients

  cli::cli_text(sprintf(
    "[%3d/%3d] L=%10.2f | dL=%.3e (%.2e) | dtheta=%.2e",
    iter,
    control$maxit,
    curr$loglik,
    metrics$delta_l,
    metrics$rel_l,
    metrics$delta_theta
  ))

  if (control$verbose >= 2) {
    cli::cli_text(sprintf(
      "    variance: sigma_e=%.3f, tr(Sigma_b)=%.2f",
      cf$measurement_error_sd,
      sum(diag(cf$random_effect_sigma))
    ))
    cli::cli_text(sprintf("    baseline: %s", .format_vector(cf$baseline)))
    cli::cli_text(sprintf("    hazard: %s", .format_vector(cf$hazard)))
    cli::cli_text(sprintf(
      "    longitudinal: %s",
      .format_vector(cf$longitudinal, 6)
    ))
  }
}

#' @noRd
.track <- function(iter, curr, prev, control) {
  metrics <- .compute_metrics(curr, prev, iter)

  if (iter > 1 && metrics$delta_l < -1e-6 && control$verbose > 0) {
    cli::cli_alert_warning(sprintf(
      "Log-likelihood decreased by %.4e at iteration %d",
      abs(metrics$delta_l),
      iter
    ))
  }

  converged <- iter > 1 &&
    metrics$rel_l < control$rtol &&
    metrics$delta_theta < control$atol

  is_final <- iter == control$maxit

  if (control$verbose > 0 && !(is_final && !converged)) {
    .print_iteration(iter, curr, metrics, control)
  }

  if (control$verbose > 0) {
    if (converged) {
      cli::cli_text("")
      cli::cli_alert_success(sprintf(
        "Converged in %d iterations (rel dL=%.2e, dtheta=%.2e)",
        iter,
        metrics$rel_l,
        metrics$delta_theta
      ))
    } else if (is_final) {
      cli::cli_text("")
      cli::cli_alert_warning(sprintf(
        paste0(
          "Not converged after %d iterations ",
          "(rel dL=%.2e > %.2e, dtheta=%.2e > %.2e)"
        ),
        control$maxit,
        metrics$rel_l,
        control$rtol,
        metrics$delta_theta,
        control$atol
      ))
      cli::cli_alert_info(
        paste0(
          "Try increasing maxit, relaxing tolerances (rtol/atol), ",
          "or adjusting initial values"
        )
      )
    }
  }

  list(converged = converged)
}

# Model Finalization ===========================================================

#' @noRd
.finalize <- function(
  data_list,
  parameters,
  random_effects,
  loglik,
  control,
  coef_names,
  converged
) {
  posteriors <- .compute_posteriors(
    data_list,
    parameters,
    random_effects,
    control$parallel,
    control$n_cores,
    control$quad_level
  )

  n_subjects <- length(data_list)

  posterior_moments <- lapply(seq_len(n_subjects), function(i) {
    .compute_posterior_moments(
      list(nodes = posteriors$nodes[[i]], weights = posteriors$weights[[i]])
    )
  })

  random_effects <- t(sapply(posterior_moments, `[[`, "mean"))

  names(parameters$coefficients$baseline) <- coef_names$baseline
  names(parameters$coefficients$hazard) <- coef_names$hazard
  names(parameters$coefficients$longitudinal) <- coef_names$longitudinal

  coef_names_expanded <- c(
    paste0("baseline:", coef_names$baseline),
    paste0("hazard:", coef_names$hazard),
    paste0("longitudinal:", coef_names$longitudinal)
  )

  vcov_matrix <- .compute_vcov_sem(
    data_list,
    posteriors,
    parameters,
    random_effects,
    control,
    verbose = control$verbose
  )
  dimnames(vcov_matrix) <- list(coef_names_expanded, coef_names_expanded)

  n_coef <- with(
    parameters$coefficients,
    length(baseline) + length(hazard) + length(longitudinal)
  )
  n_sigma_b <- with(
    parameters$coefficients,
    length(longitudinal) * (length(longitudinal) + 1) / 2
  )
  n_params <- n_coef + 1 + n_sigma_b
  aic <- -2 * loglik + 2 * n_params
  bic <- -2 * loglik + n_params * log(n_subjects)

  ode_solutions <- .solve_batch_ode_cppad(data_list, random_effects, parameters)
  # nolint start: object_usage_linter
  risk_scores <- vapply(
    ode_solutions,
    function(x) tail(x$log_hazard, 1),
    numeric(1)
  )
  event_times <- vapply(data_list, `[[`, numeric(1), "time")
  event_status <- vapply(data_list, `[[`, numeric(1), "status")

  # Compute concordance index (C-index) for model discrimination
  cindex <- survival::concordance(
    Surv(event_times, event_status) ~ risk_scores,
    reverse = TRUE
  )$concordance
  # nolint end

  if (control$verbose > 0) {
    cli::cli_alert_info(sprintf("C-index (concordance): %.3f", cindex))
  }

  rm(posteriors, posterior_moments)
  gc(verbose = FALSE, full = FALSE)

  list(
    random_effects = random_effects,
    vcov = vcov_matrix,
    loglik = loglik,
    aic = aic,
    bic = bic,
    cindex = cindex
  )
}
