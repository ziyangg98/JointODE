# Utility Functions for JointODE Package

# Shared Helpers ===============================================================

#' @noRd
.logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#' @noRd
.n_obs <- function(data_list) {
  sum(vapply(
    data_list,
    function(s) length(s$longitudinal$measurements),
    integer(1)
  ))
}

#' @noRd
.coef_table <- function(estimates, std_errors) {
  z <- estimates / std_errors
  cbind(
    Estimate = estimates, `Std. Error` = std_errors,
    `z value` = z, `Pr(>|z|)` = 2 * pnorm(-abs(z))
  )
}

#' @noRd
.print_coefmat <- function(x, digits = 4, signif.stars = TRUE, ...) {
  x[, 1:2] <- round(x[, 1:2], digits)
  printCoefmat(x, digits = digits, signif.stars = signif.stars, ...)
}

#' @noRd
.extend_covariates <- function(cov_mat, orig_times, pred_times) {
  if (is.null(cov_mat)) {
    return(NULL)
  }
  if (is.matrix(cov_mat) && length(cov_mat) == 0) {
    return(matrix(numeric(0), nrow = length(pred_times), ncol = 0))
  }
  indices <- findInterval(pred_times, orig_times)
  indices[indices == 0L] <- 1L
  if (is.matrix(cov_mat)) cov_mat[indices, , drop = FALSE] else cov_mat[indices]
}

# Constants ====================================================================

.default_spline <- list(
  degree = 2,
  n_knots = 1,
  knot_placement = "equal",
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

#' @noRd
.format_matrix <- function(mat) {
  if (!is.matrix(mat) || !length(mat)) {
    return("[]")
  }
  rows <- apply(mat, 1, function(row) {
    paste(sprintf("%.3f", row), collapse = ", ")
  })
  paste0("[", paste(paste0("[", rows, "]"), collapse = ", "), "]")
}

# Parallel Computing ===========================================================

#' @importFrom parallel detectCores
#' @importFrom future plan multicore multisession sequential supportsMulticore
#' @noRd
.resolve_cores <- function(n_cores) {
  if (n_cores == 0L) {
    detected <- parallel::detectCores()
    if (is.na(detected)) 1L else max(1L, detected - 1L)
  } else {
    n_cores
  }
}

#' @noRd
.setup_parallel_plan <- function(n_cores = 0) {
  n_cores <- .resolve_cores(n_cores)
  strategy <- if (future::supportsMulticore()) {
    future::multicore
  } else {
    future::multisession
  }
  future::plan(strategy, workers = n_cores)
  function() future::plan(future::sequential)
}

#' @importFrom future.apply future_lapply
#' @noRd
.parallel_apply <- function(
  indices, fn, parallel = TRUE, n_cores = 0, setup = TRUE
) {
  if (!parallel) {
    return(lapply(indices, fn))
  }

  if (setup) {
    cleanup <- .setup_parallel_plan(n_cores)
    on.exit(cleanup(), add = TRUE)
  }

  future.apply::future_lapply(
    indices, fn,
    future.seed = TRUE,
    future.packages = "JointODE"
  )
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
.compute_dimensions <- function(parsed_long, parsed_surv, spline_config) {
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
    parsed_long$biomarker$fixed, parsed_long$velocity$fixed
  )
  n_biomarker_velocity_random <- sum(
    parsed_long$biomarker$random, parsed_long$velocity$random
  )

  list(
    n_longitudinal_coef = n_longitudinal_fixed + n_biomarker_velocity_fixed,
    n_random_effects = n_longitudinal_random + n_biomarker_velocity_random + 2,
    n_survival_covariates = n_survival_covariates,
    n_spline_basis = spline_config$degree + spline_config$n_knots + 1,
    spline_config = spline_config
  )
}

#' @importFrom stats quantile
#' @noRd
.get_spline_config <- function(
  x,
  degree = 2,
  n_knots = 1,
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

#' @noRd
.safe_chol <- function(H) {
  R <- try(chol(H), silent = TRUE)
  if (!inherits(R, "try-error")) {
    return(R)
  }
  tau <- 1e-4 * max(abs(diag(H)), 1)
  for (k in seq_len(10)) {
    R <- try(chol(H + diag(tau, nrow(H))), silent = TRUE)
    if (!inherits(R, "try-error")) {
      return(R)
    }
    tau <- tau * 10
  }
  stop("Hessian is not positive definite")
}


# Parameter Counting & Conversion =============================================

#' @noRd
.count_params <- function(parameters) {
  cf <- parameters$coefficients
  p <- nrow(cf$random_effect_sigma)
  length(cf$baseline) + length(cf$hazard) + length(cf$longitudinal) +
    length(cf$initial_state) + 1 + p * (p + 1) / 2
}

#' @noRd
.coef_to_vector <- function(parameters) {
  with(
    parameters$coefficients,
    c(baseline, hazard, longitudinal, initial_state)
  )
}

#' @noRd
.vector_to_coef <- function(parameters, theta) {
  cf <- parameters$coefficients
  n <- c(
    length(cf$baseline), length(cf$hazard),
    length(cf$longitudinal), length(cf$initial_state)
  )
  idx <- cumsum(n)

  parameters$coefficients$baseline <- theta[1:idx[1]]
  parameters$coefficients$hazard <- theta[(idx[1] + 1):idx[2]]
  parameters$coefficients$longitudinal <- theta[(idx[2] + 1):idx[3]]
  parameters$coefficients$initial_state <- theta[(idx[3] + 1):idx[4]]
  parameters
}

# EM Progress Tracking =========================================================

#' @noRd
.compute_metrics <- function(curr, prev, iter) {
  if (iter > 1) {
    delta_l <- curr$loglik - prev$loglik
    rel_l <- abs(delta_l) / (abs(curr$loglik) + 1)

    # Parameter change: max absolute change across fixed effects + variance
    theta_curr <- c(
      .coef_to_vector(curr$parameters),
      curr$parameters$coefficients$measurement_error_sd,
      as.vector(curr$parameters$coefficients$random_effect_sigma)
    )
    theta_prev <- c(
      .coef_to_vector(prev$parameters),
      prev$parameters$coefficients$measurement_error_sd,
      as.vector(prev$parameters$coefficients$random_effect_sigma)
    )
    delta_theta <- max(abs(theta_curr - theta_prev))
  } else {
    delta_l <- curr$loglik
    rel_l <- 1
    delta_theta <- 1
  }
  list(delta_l = delta_l, rel_l = rel_l, delta_theta = delta_theta)
}

#' @noRd
.print_iteration <- function(iter, curr, metrics, control) {
  cf <- curr$parameters$coefficients

  cli::cli_text(sprintf(
    "[%3d/%3d] L=%10.2f | dL=%+.2e  dTheta=%.2e",
    iter, control$maxit, curr$loglik,
    metrics$delta_l, metrics$delta_theta
  ))

  if (control$verbose >= 2) {
    cli::cli_text(sprintf("    sigma_e: %.5f", cf$measurement_error_sd))
    cli::cli_text(sprintf(
      "    Sigma_b diag: %s",
      paste(sprintf("%.4f", diag(cf$random_effect_sigma)), collapse = ", ")
    ))
    cli::cli_text(sprintf("    baseline: %s", .format_vector(cf$baseline)))
    cli::cli_text(sprintf("    hazard: %s", .format_vector(cf$hazard)))
    cli::cli_text(sprintf(
      "    longitudinal: %s", .format_vector(cf$longitudinal, 6)
    ))
    if (!is.null(cf$initial_state)) {
      cli::cli_text(sprintf(
        "    initial_state: %s",
        paste(sprintf("%.4f", cf$initial_state), collapse = ", ")
      ))
    }
  }
  if (control$verbose >= 3) {
    re <- curr$random_effects
    if (!is.null(re)) {
      for (k in seq_len(ncol(re))) {
        cli::cli_text(sprintf(
          "    Random effect[,%d] range: [%.3f, %.3f]",
          k, min(re[, k]), max(re[, k])
        ))
      }
    }
  }
}

#' @noRd
.track <- function(iter, curr, prev, control) {
  metrics <- .compute_metrics(curr, prev, iter)

  if (is.na(metrics$delta_l) || is.na(metrics$rel_l)) {
    if (control$verbose > 0) {
      cli::cli_alert_warning("Log-likelihood is NA at iteration {iter}")
    }
    return(list(converged = FALSE, metrics = metrics))
  }

  converged <- iter > 1 && metrics$delta_theta < control$tol
  is_final <- iter == control$maxit

  if (control$verbose > 0 && !(is_final && !converged)) {
    .print_iteration(iter, curr, metrics, control)
  }

  if (control$verbose > 0) {
    if (converged) {
      cli::cli_text("")
      cli::cli_alert_success(sprintf(
        "Converged in %d iterations (dTheta=%.2e)",
        iter, metrics$delta_theta
      ))
    } else if (is_final) {
      cli::cli_text("")
      cli::cli_alert_warning(sprintf(
        "Not converged after %d iterations (dTheta=%.2e > %.2e)",
        control$maxit, metrics$delta_theta, control$tol
      ))
      cli::cli_alert_info(
        "Try increasing maxit, relaxing tolerances, or adjusting initial values"
      )
    }
  }

  list(converged = converged, metrics = metrics)
}
