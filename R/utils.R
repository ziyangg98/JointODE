# Utility Functions for JointODE Package

# Shared Helpers ===============================================================

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
  if (is.null(cov_mat)) return(NULL)
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
  if (!length(x)) return("[]")
  shown <- head(x, n)
  rest <- length(x) - n
  paste0(
    "[", paste(sprintf("%.3f", shown), collapse = ", "),
    if (rest > 0) paste0(", ...+", rest) else "", "]"
  )
}

# Formula Parsing ==============================================================

#' @importFrom stats as.formula
#' @noRd
.build_formula <- function(terms, response = NULL, is_random = FALSE) {
  covars <- setdiff(terms, "(Intercept)")
  has_intercept <- "(Intercept)" %in% terms

  terms_str <- if (length(covars) > 0) {
    if (is_random) covars <- gsub("^\\(Intercept\\)$", "1", covars)
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
      boundary_knots[1], boundary_knots[2],
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

# TMB Result Extraction ========================================================

#' @noRd
.extract_tmb_results <- function(obj, opt, parameters, coef_names,
                                  data_list, n_re, control) {
  rep <- TMB::sdreport(obj)
  reported <- obj$report()
  par_full <- obj$env$last.par.best
  pn <- names(par_full)

  # Update coefficients from optimized values
  cf <- parameters$coefficients
  cf$baseline <- setNames(as.numeric(par_full[pn == "baseline"]),
                          coef_names$baseline)
  cf$hazard <- setNames(as.numeric(par_full[pn == "hazard"]),
                        coef_names$hazard)
  cf$longitudinal <- setNames(as.numeric(par_full[pn == "longitudinal"]),
                              coef_names$longitudinal)
  cf$initial_state <- setNames(as.numeric(par_full[pn == "initial_state"]),
                               coef_names$initial_state)
  cf$measurement_error_sd <- exp(par_full[pn == "log_sigma_e"])
  cf$random_effect_sigma <- as.matrix(reported$Sigma_b)
  parameters$coefficients <- cf
  parameters$random_effects_init <- NULL

  # Random effects posterior modes
  random_effects <- matrix(par_full[pn == "b"],
                           nrow = length(data_list), ncol = n_re)

  coef_names_exp <- .prefixed_coef_names(coef_names)
  n_fixed <- length(coef_names_exp)
  vcov_full <- rep$cov.fixed
  vcov_matrix <- if (!is.null(vcov_full) && nrow(vcov_full) >= n_fixed) {
    vcov_full[seq_len(n_fixed), seq_len(n_fixed), drop = FALSE]
  } else {
    matrix(NA, n_fixed, n_fixed)
  }
  dimnames(vcov_matrix) <- list(coef_names_exp, coef_names_exp)

  # Model fit statistics
  loglik <- -opt$objective
  n_params <- .count_params(parameters)
  n_subjects <- length(data_list)

  # C-index
  event_t <- vapply(data_list, `[[`, numeric(1), "time")
  event_s <- vapply(data_list, `[[`, numeric(1), "status")
  cindex <- survival::concordance(
    Surv(event_t, event_s) ~ as.numeric(reported$log_hazard_at_event),
    reverse = TRUE
  )$concordance

  converged <- opt$convergence == 0
  if (control$verbose > 0) {
    if (converged) cli::cli_alert_success(sprintf("Converged (%s)", opt$message))
    else cli::cli_alert_warning(sprintf("Did not converge: %s", opt$message))
    cli::cli_alert_info(sprintf("Log-likelihood: %.2f", loglik))
    cli::cli_alert_info(sprintf("C-index: %.3f", cindex))
  }

  list(
    parameters = parameters,
    logLik = loglik,
    AIC = -2 * loglik + 2 * n_params,
    BIC = -2 * loglik + n_params * log(n_subjects),
    cindex = cindex,
    convergence = list(
      converged = converged,
      iterations = opt$iterations,
      message = sprintf("%s (%s)",
        if (converged) "Converged" else "Did not converge", opt$message)
    ),
    random_effects = random_effects,
    vcov = vcov_matrix,
    tmb_report = reported
  )
}

# Parameter Helpers =============================================================

#' @noRd
.prefixed_coef_names <- function(coef_names) {
  c(
    paste0("baseline:", coef_names$baseline),
    paste0("hazard:", coef_names$hazard),
    paste0("longitudinal:", coef_names$longitudinal),
    paste0("initial state:", coef_names$initial_state)
  )
}

#' @noRd
.count_params <- function(parameters) {
  cf <- parameters$coefficients
  p <- nrow(cf$random_effect_sigma)
  length(cf$baseline) + length(cf$hazard) + length(cf$longitudinal) +
    length(cf$initial_state) + 1 + p * (p + 1) / 2
}
