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
.setup_openmp <- function(control) {
  if (control$parallel) {
    n_cores <- if (control$n_cores > 0) control$n_cores else parallel::detectCores()
    TMB::openmp(n_cores)
  }
}

#' Initialize RE columns 1-2 from first observations
#' @noRd
.init_re_from_observations <- function(data_list, m0, v0 = 0, n_re = 2L) {
  n_subjects <- length(data_list)
  re <- matrix(0, nrow = n_subjects, ncol = max(n_re, 2L))
  for (i in seq_along(data_list)) {
    obs <- data_list[[i]]$longitudinal
    if (length(obs$measurements) >= 1)
      re[i, 1] <- obs$measurements[1] - m0
    if (length(obs$measurements) >= 2) {
      dt <- obs$times[2] - obs$times[1]
      if (dt > 0)
        re[i, 2] <- (obs$measurements[2] - obs$measurements[1]) / dt - v0
    }
  }
  re
}

#' Build counting process (start, stop, event) from predictions + survival
#' @noRd
.build_counting_process <- function(pred, survival_data, id, surv_vars) {
  names(pred)[names(pred) == "time"] <- "obstime"
  merged <- merge(pred, survival_data, by = id, all.x = TRUE)
  merged <- merged[order(merged[[id]], merged$obstime), ]

  merged$start <- merged$obstime
  subj_ids <- as.vector(merged[[id]])
  subj_idx <- c(which(!duplicated(subj_ids)), nrow(merged) + 1)
  ni <- diff(subj_idx)

  stop_vec <- event_vec <- numeric(nrow(merged))
  for (i in seq_along(ni)) {
    s <- subj_idx[i]
    e <- subj_idx[i + 1] - 1
    if (ni[i] == 1) {
      stop_vec[s] <- merged[[surv_vars[1]]][s]
      event_vec[s] <- merged[[surv_vars[2]]][s]
    } else {
      stop_vec[s:e] <- c(merged$obstime[(s + 1):e], merged[[surv_vars[1]]][e])
      event_vec[s:e] <- c(rep(0, ni[i] - 1), merged[[surv_vars[2]]][e])
    }
  }
  merged$stop <- stop_vec
  merged$event <- event_vec
  merged
}

# Constants ====================================================================

.default_spline <- list(
  degree = 2,
  n_knots = 1,
  knot_placement = "equal",
  boundary_knots = NULL
)

.reserved_words <- c("biomarker", "velocity")

.init_state_names <- c("init_biomarker", "init_velocity")

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
