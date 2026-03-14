# ==============================================================================
# Utility Functions Unit Tests
# ==============================================================================

# --- .count_params ---

test_that(".count_params counts correctly", {
  params <- list(
    coefficients = list(
      baseline = c(0.1, 0.2, 0.3),
      hazard = c(0.5, -0.3),
      longitudinal = c(-0.01, 0.02, 0.1),
      measurement_error_sd = 0.5,
      random_effect_sigma = diag(2)
    )
  )
  # 3 + 2 + 3 + 1 (sigma_e) + 2*(2+1)/2 (lower tri of 2x2) = 12
  expect_equal(.count_params(params), 12)
})

test_that(".count_params handles different dimensions", {
  params <- list(
    coefficients = list(
      baseline = numeric(5),
      hazard = numeric(3),
      longitudinal = numeric(4),
      measurement_error_sd = 0.1,
      random_effect_sigma = diag(3)
    )
  )
  # 5 + 3 + 4 + 1 + 3*(3+1)/2 = 19
  expect_equal(.count_params(params), 19)
})

# --- .coef_table ---

test_that(".coef_table returns correct structure", {
  est <- c(a = 1.5, b = -0.3, c = 0.0)
  se <- c(a = 0.5, b = 0.1, c = 0.2)
  tbl <- .coef_table(est, se)

  expect_true(is.matrix(tbl))
  expect_equal(ncol(tbl), 4)
  expect_equal(colnames(tbl), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_equal(nrow(tbl), 3)
})

test_that(".coef_table computes z and p correctly", {
  est <- c(2.0, -1.0)
  se <- c(0.5, 0.25)
  tbl <- .coef_table(est, se)

  expect_equal(tbl[, "z value"], est / se)
  expect_true(all(tbl[, "Pr(>|z|)"] >= 0 & tbl[, "Pr(>|z|)"] <= 1))
  expect_equal(tbl[, "Pr(>|z|)"], 2 * pnorm(-abs(est / se)))
})

test_that(".coef_table handles zero std error without error", {
  tbl <- .coef_table(c(1.0), c(0.0))
  expect_true(is.infinite(tbl[1, "z value"]))
})

# --- .coef_to_vector / .vector_to_coef ---

test_that(".coef_to_vector roundtrips with .vector_to_coef", {
  params <- sim$init
  theta <- .coef_to_vector(params)

  expect_true(is.numeric(theta))
  n_expected <- length(params$coefficients$baseline) +
    length(params$coefficients$hazard) +
    length(params$coefficients$longitudinal)
  expect_equal(length(theta), n_expected)

  params_recovered <- .vector_to_coef(params, theta)
  expect_equal(unname(params_recovered$coefficients$baseline),
               unname(params$coefficients$baseline))
  expect_equal(unname(params_recovered$coefficients$hazard),
               unname(params$coefficients$hazard))
  expect_equal(unname(params_recovered$coefficients$longitudinal),
               unname(params$coefficients$longitudinal))
})

test_that(".coef_to_vector preserves parameter ordering", {
  params <- sim$init
  theta <- .coef_to_vector(params)

  n_b <- length(params$coefficients$baseline)
  n_h <- length(params$coefficients$hazard)
  n_l <- length(params$coefficients$longitudinal)

  expect_equal(unname(theta[1:n_b]), unname(params$coefficients$baseline))
  expect_equal(unname(theta[(n_b + 1):(n_b + n_h)]),
               unname(params$coefficients$hazard))
  expect_equal(unname(theta[(n_b + n_h + 1):(n_b + n_h + n_l)]),
               unname(params$coefficients$longitudinal))
})

# --- .compute_dimensions ---

test_that(".compute_dimensions returns correct values", {
  parsed_long <- .parse_longitudinal_formula(observed ~ x1 + x2)
  parsed_surv <- .parse_survival_formula(Surv(time, status) ~ w1 + w2)
  spline_config <- list(degree = 2, n_knots = 0)

  dims <- .compute_dimensions(parsed_long, parsed_surv, spline_config)

  expect_type(dims, "list")
  expect_true("n_longitudinal_coef" %in% names(dims))
  expect_true("n_random_effects" %in% names(dims))
  expect_true("n_survival_covariates" %in% names(dims))
  expect_true("n_spline_basis" %in% names(dims))
  expect_equal(dims$n_survival_covariates, 2) # w1, w2
  expect_equal(dims$n_spline_basis, 3) # degree + n_knots + 1
})

test_that(".compute_dimensions handles no covariates", {
  parsed_long <- .parse_longitudinal_formula(y ~ 1)
  parsed_surv <- .parse_survival_formula(Surv(time, status) ~ 1)
  spline_config <- list(degree = 2, n_knots = 2)

  dims <- .compute_dimensions(parsed_long, parsed_surv, spline_config)
  expect_equal(dims$n_survival_covariates, 0)
  expect_equal(dims$n_spline_basis, 5) # 2 + 2 + 1
})

# --- .parse_longitudinal_formula ---

test_that(".parse_longitudinal_formula extracts components", {
  parsed <- .parse_longitudinal_formula(y ~ x1 + x2)

  expect_equal(parsed$response, "y")
  expect_true("x1" %in% parsed$fixed_terms)
  expect_true("x2" %in% parsed$fixed_terms)
  expect_null(parsed$random_terms)
  expect_null(parsed$grouping)
})

test_that(".parse_longitudinal_formula handles random effects", {
  parsed <- .parse_longitudinal_formula(
    y ~ biomarker + velocity + x1 + (biomarker + velocity | id)
  )

  expect_equal(parsed$response, "y")
  expect_true(parsed$biomarker$fixed)
  expect_true(parsed$biomarker$random)
  expect_true(parsed$velocity$fixed)
  expect_true(parsed$velocity$random)
  expect_equal(parsed$grouping, "id")
  # biomarker/velocity removed from random_terms (reserved words)
  expect_null(parsed$random_terms)
})

test_that(".parse_longitudinal_formula handles intercept only", {
  parsed <- .parse_longitudinal_formula(y ~ 1)

  expect_equal(parsed$response, "y")
  expect_true("(Intercept)" %in% parsed$fixed_terms)
  expect_null(parsed$random_terms)
})

# --- .parse_survival_formula ---

test_that(".parse_survival_formula extracts components", {
  parsed <- .parse_survival_formula(Surv(time, status) ~ w1 + w2)

  expect_equal(parsed$time_var, "time")
  expect_equal(parsed$status_var, "status")
  expect_equal(parsed$covariate_terms, c("w1", "w2"))
})

test_that(".parse_survival_formula handles no covariates", {
  parsed <- .parse_survival_formula(Surv(t, d) ~ 1)

  expect_equal(parsed$time_var, "t")
  expect_equal(parsed$status_var, "d")
  expect_null(parsed$covariate_terms)
})

test_that(".parse_survival_formula validates Surv()", {
  expect_error(
    .parse_survival_formula(y ~ x1),
    "Surv"
  )
})
