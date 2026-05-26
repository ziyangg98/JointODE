# Model Setup

#' @importFrom stats model.frame model.matrix model.response
#' @noRd
.setup_model <- function(
  longitudinal_data, survival_data, survival_formula,
  gamma, parsed_long, parsed_surv, spline_config
) {
  # Fixed covariate names
  fixed_formula <- .build_formula(
    parsed_long$fixed_terms,
    response = parsed_long$response
  )
  long_fixed_names <- colnames(model.matrix(
    fixed_formula, model.frame(fixed_formula, longitudinal_data)
  ))

  # Random effects dimension
  has_re_covs <- !is.null(parsed_long$random_terms)
  random_terms <- if (has_re_covs) parsed_long$random_terms else character(0)

  n_long_random <- if (length(random_terms) > 0) {
    ncol(model.matrix(
      .build_formula(random_terms, is_random = TRUE), longitudinal_data
    ))
  } else {
    0
  }
  n_re <- n_long_random + sum(
    parsed_long$lambda$random, parsed_long$tau$random
  )

  # Survival dimensions
  surv_frame <- model.frame(survival_formula, survival_data)
  event_times <- model.response(surv_frame)[, 1]
  has_surv_covs <- length(all.vars(survival_formula[[3]])) > 0 &&
    survival_formula[[3]] != 1
  surv_names <- if (has_surv_covs) {
    colnames(model.matrix(survival_formula, surv_frame)[, -1, drop = FALSE])
  } else {
    character(0)
  }

  # Spline configuration from data
  sbc <- .get_spline_config(
    x = event_times,
    degree = spline_config$degree,
    n_knots = spline_config$n_knots,
    knot_placement = spline_config$knot_placement,
    boundary_knots = spline_config$boundary_knots
  )
  sbc$boundary_knots[1] <- 0

  # Longitudinal coefficient names (dynamic parameters first)
  long_names <- character(0)
  if (parsed_long$lambda$fixed) long_names <- c(long_names, "lambda")
  if (parsed_long$tau$fixed) long_names <- c(long_names, "tau")
  long_names <- c(long_names, long_fixed_names)

  # Random effects layout: [initial_biomarker, initial_velocity, dyn_coefs...]
  n_re_total <- n_re + 2
  random_effects <- matrix(0, nrow(survival_data), n_re_total)

  re_names <- c("initial_biomarker", "initial_velocity")
  if (parsed_long$lambda$random) re_names <- c(re_names, "lambda")
  if (parsed_long$tau$random) re_names <- c(re_names, "tau")
  if (length(random_terms) > 0) {
    re_cov_names <- colnames(model.matrix(
      .build_formula(random_terms, is_random = TRUE), longitudinal_data
    ))
    re_names <- c(re_names, paste0("forcing_", re_cov_names))
  }

  list(
    dims = .compute_dimensions(parsed_long, parsed_surv, spline_config),
    random_effects = random_effects,
    re_names = re_names,
    coef_names = list(
      baseline = paste0("bs", seq_len(sbc$df)),
      hazard = c("value", "velocity", surv_names),
      longitudinal = long_names,
      initial_state = c("initial_biomarker", "initial_velocity")
    ),
    spline_baseline_config = sbc
  )
}

#' Setup marginal model dimensions and coefficient names
#' @noRd
.setup_marginal_model <- function(data, parsed_long) {
  fixed_formula <- .build_formula(parsed_long$fixed_terms,
    response = parsed_long$response
  )
  fixed_names <- colnames(model.matrix(
    fixed_formula, model.frame(fixed_formula, data)
  ))

  long_names <- character(0)
  if (parsed_long$lambda$fixed) long_names <- c(long_names, "lambda")
  if (parsed_long$tau$fixed) long_names <- c(long_names, "tau")
  long_names <- c(long_names, fixed_names)

  # RE dimension
  has_re_covs <- !is.null(parsed_long$random_terms)
  random_terms <- if (has_re_covs) parsed_long$random_terms else character(0)
  n_long_random <- if (length(random_terms) > 0) {
    ncol(model.matrix(.build_formula(random_terms, is_random = TRUE), data))
  } else {
    0L
  }
  n_re <- 2L + n_long_random +
    sum(parsed_long$lambda$random, parsed_long$tau$random)

  re_names <- c("initial_biomarker", "initial_velocity")
  if (parsed_long$lambda$random) re_names <- c(re_names, "lambda")
  if (parsed_long$tau$random) re_names <- c(re_names, "tau")
  if (length(random_terms) > 0) {
    re_cov_names <- colnames(model.matrix(
      .build_formula(random_terms, is_random = TRUE), data
    ))
    re_names <- c(re_names, paste0("forcing_", re_cov_names))
  }

  list(
    n_longitudinal_coef = length(long_names),
    n_re = n_re,
    re_names = re_names,
    coef_names = list(
      longitudinal = long_names,
      initial_state = .init_state_names
    )
  )
}
