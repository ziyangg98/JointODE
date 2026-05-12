# ==============================================================================
# Input Module: Formula Parsing, Data Processing, Control, Parameters
# ==============================================================================

# --- Local helper for formula validation tests ---

.formula_test_data <- function(n = 10, n_subjects = 2) {
  data.frame(
    y = rnorm(n),
    id = rep(1:n_subjects, each = n / n_subjects),
    x = rnorm(n),
    time = rep(1:(n / n_subjects), n_subjects)
  )
}

# ==============================================================================
# Formula Parsing
# ==============================================================================

# --- .parse_longitudinal_formula ---

test_that("parse formula with random effects", {
  result <- .parse_longitudinal_formula(y ~ x1 + x2 + (1 + x1 | subject))
  expect_equal(result$response, "y")
  expect_equal(result$fixed_terms, c("(Intercept)", "x1", "x2"))
  expect_equal(result$random_terms, c("(Intercept)", "x1"))
  expect_equal(result$grouping, "subject")
})

test_that("parse formula without random effects", {
  result <- .parse_longitudinal_formula(y ~ x1 + x2)
  expect_equal(result$response, "y")
  expect_null(result$random_terms)
  expect_false(result$biomarker$random)
})

test_that("parse formula with reserved words", {
  result <- .parse_longitudinal_formula(y ~ x + (biomarker + velocity | id))
  expect_null(result$random_terms)
  expect_true(result$biomarker$random)
  expect_true(result$velocity$random)
  expect_false(result$biomarker$fixed)
})

test_that("parse formula with reserved words in fixed effects", {
  result <- .parse_longitudinal_formula(y ~ biomarker + velocity + x + (1 | id))
  expect_true(result$biomarker$fixed)
  expect_true(result$velocity$fixed)
  expect_equal(result$fixed_terms, c("(Intercept)", "x"))
})

test_that("parse formula with reserved words in both", {
  result <- .parse_longitudinal_formula(
    y ~ biomarker + (biomarker + velocity | id)
  )
  expect_true(result$biomarker$fixed)
  expect_true(result$biomarker$random)
  expect_true(result$velocity$random)
})

test_that("parse random intercept/slope variants", {
  expect_equal(
    .parse_longitudinal_formula(y ~ x + (1 | id))$random_terms,
    "(Intercept)"
  )
  expect_equal(
    .parse_longitudinal_formula(y ~ x + (x | id))$random_terms,
    "x"
  )
})

test_that("parse interaction terms", {
  result <- .parse_longitudinal_formula(y ~ x1 * x2 + (1 | id))
  expect_true("x1:x2" %in% result$fixed_terms)
})

test_that("variable in fixed only / random only / both", {
  r1 <- .parse_longitudinal_formula(y ~ time + (biomarker | id))
  expect_true("time" %in% r1$fixed_terms)

  r2 <- .parse_longitudinal_formula(y ~ x + (time | id))
  expect_true("time" %in% r2$random_terms)
  expect_false("time" %in% r2$fixed_terms)

  r3 <- .parse_longitudinal_formula(y ~ time + (time | id))
  expect_true("time" %in% r3$fixed_terms)
  expect_true("time" %in% r3$random_terms)
})

# --- .validate_longitudinal_formula ---

test_that("validate accepts valid formula", {
  expect_silent(.validate_longitudinal_formula(
    y ~ x + (1 | id), .formula_test_data()
  ))
  expect_silent(.validate_longitudinal_formula(
    y ~ x + (biomarker + velocity | id), .formula_test_data()
  ))
  # No RE term is valid (initial state RE are always included)
  expect_silent(.validate_longitudinal_formula(y ~ x, .formula_test_data()))
})

test_that("validate rejects invalid longitudinal formulas", {
  d <- .formula_test_data()
  expect_error(.validate_longitudinal_formula("bad", d), "formula object")
  expect_error(.validate_longitudinal_formula(~x, d), "two-sided")
  expect_error(
    .validate_longitudinal_formula(y ~ missing + (1 | id), d), "not found"
  )
  expect_error(
    .validate_longitudinal_formula(y ~ x + (1 | wrong), d), "Grouping"
  )
  expect_error(
    .validate_longitudinal_formula(y ~ x + (1 | id) + (1 | time), d),
    "Multiple random"
  )
})

test_that("validate rejects reserved words in data", {
  d <- data.frame(y = 1:10, biomarker = 1:10, id = rep(1:2, each = 5))
  expect_error(.validate_longitudinal_formula(y ~ 1 + (1 | id), d), "Reserved")
})

test_that("validate rejects nested random effects", {
  d <- data.frame(
    y = 1:20, x = 1:20, school = rep(1:2, 10), class = rep(1:4, 5)
  )
  expect_error(
    .validate_longitudinal_formula(y ~ x + (1 | school / class), d), "Nested"
  )
})

# --- .parse_survival_formula ---

test_that("parse survival formula", {
  r <- .parse_survival_formula(Surv(time, status) ~ age + treatment)
  expect_equal(r$time_var, "time")
  expect_equal(r$status_var, "status")
  expect_equal(r$covariate_terms, c("age", "treatment"))

  r2 <- .parse_survival_formula(Surv(survtime, event) ~ 1)
  expect_null(r2$covariate_terms)

  r3 <- .parse_survival_formula(Surv(t, d) ~ age * treatment)
  expect_true("age:treatment" %in% r3$covariate_terms)
})

test_that("parse survival formula rejects invalid LHS", {
  expect_error(.parse_survival_formula(time ~ age), "Surv\\(\\)")
  expect_error(.parse_survival_formula(Surv(time) ~ age), "at least")
})

# --- .validate_survival_formula ---

test_that("validate accepts valid survival formula", {
  d <- data.frame(
    time = c(10, 20, 30), status = c(1, 0, 1), age = c(50, 60, 70)
  )
  expect_silent(.validate_survival_formula(Surv(time, status) ~ age, d))
})

test_that("validate rejects invalid survival formulas", {
  d <- data.frame(time = 1:3, status = c(1, 0, 1))
  expect_error(.validate_survival_formula("bad", d), "formula object")
  expect_error(.validate_survival_formula(~age, d), "two-sided")
  expect_error(
    .validate_survival_formula(Surv(time, status) ~ missing, d),
    "not found"
  )
})

test_that("validate rejects invalid time/status values", {
  expect_error(
    .validate_survival_formula(
      Surv(time, status) ~ 1,
      data.frame(time = c(-1, 0, 10), status = c(1, 0, 1))
    ), "Invalid observation"
  )
  expect_error(
    .validate_survival_formula(
      Surv(time, status) ~ 1,
      data.frame(time = c(10, 20, 30), status = c(1, 2, 0))
    ), "Invalid status"
  )
  expect_error(
    .validate_survival_formula(
      Surv(time, status) ~ 1,
      data.frame(time = c(10, NA, 30), status = c(1, 0, 1))
    ), "Missing values.*time"
  )
  expect_error(
    .validate_survival_formula(
      Surv(time, status) ~ 1,
      data.frame(time = c(10, 20, 30), status = c(1, NA, 1))
    ), "Missing values.*status"
  )
})

# --- .build_formula ---

test_that(".build_formula constructs formulas correctly", {
  expect_equal(
    deparse(.build_formula(c("(Intercept)", "x1", "x2"), response = "y")),
    "y ~ x1 + x2"
  )
  expect_equal(
    deparse(.build_formula("(Intercept)", response = "y")),
    "y ~ 1"
  )
  expect_equal(
    deparse(.build_formula(c("x1", "x2"), response = "y")),
    "y ~ 0 + x1 + x2"
  )
  expect_true(inherits(
    .build_formula(c("(Intercept)", "x1"), is_random = TRUE), "formula"
  ))
})

# ==============================================================================
# Data Processing
# ==============================================================================

test_that(".process_joint structures data correctly", {
  long_data <- data.frame(
    id = rep(1:2, each = 3), obstime = rep(0:2, 2),
    v = rnorm(6), x1 = rnorm(6)
  )
  surv_data <- data.frame(
    id = 1:2, obstime = c(3, 4), status = c(1, 0), w1 = rnorm(2)
  )

  result <- .process_joint(
    longitudinal_formula = v ~ x1 + (1 | id),
    survival_formula = Surv(obstime, status) ~ w1,
    longitudinal_data = long_data, survival_data = surv_data
  )

  expect_equal(length(result), 2)
  expect_equal(names(result), c("1", "2"))

  s1 <- result[[1]]
  expect_equal(s1$id, 1)
  expect_equal(s1$time, 3)
  expect_equal(s1$status, 1)
  expect_equal(length(s1$longitudinal$times), 3)
  expect_equal(as.numeric(s1$longitudinal$measurements), long_data$v[1:3])
})

test_that(".process_joint handles missing longitudinal data", {
  result <- .process_joint(
    longitudinal_formula = v ~ 1 + (1 | id),
    survival_formula = Surv(obstime, status) ~ 1,
    longitudinal_data = data.frame(id = 1, obstime = 0, v = 1),
    survival_data = data.frame(id = 1:2, obstime = c(1, 2), status = c(1, 0))
  )
  expect_equal(length(result[[2]]$longitudinal$times), 0)
})

test_that(".process_joint sorts unordered time data", {
  long_data <- data.frame(
    id = c(1, 1, 1, 2, 2, 2),
    obstime = c(2, 0, 1, 1, 2, 0),
    v = c(3, 1, 2, 5, 6, 4)
  )
  result <- .process_joint(
    longitudinal_formula = v ~ 1 + (1 | id),
    survival_formula = Surv(obstime, status) ~ 1,
    longitudinal_data = long_data,
    survival_data = data.frame(id = 1:2, obstime = c(3, 3), status = c(1, 0))
  )
  expect_equal(result[[1]]$longitudinal$times, 0:2)
  expect_equal(as.numeric(result[[1]]$longitudinal$measurements), c(1, 2, 3))
})

# ==============================================================================
# Control Parameters
# ==============================================================================

test_that("JointODE.control returns complete defaults", {
  ctrl <- JointODE.control()
  expect_equal(ctrl$maxit, 200)
  expect_equal(ctrl$tol, 1e-4)
  expect_equal(ctrl$verbose, 0)
  expect_false(ctrl$parallel)
  expect_equal(ctrl$n_cores, 0)
  expect_equal(ctrl$hazard_quadrature, 1)
})

test_that("JointODE.control overrides and merges", {
  expect_equal(JointODE.control(maxit = 50)$maxit, 50)
  ctrl <- JointODE.control(.list = list(maxit = 100, verbose = TRUE))
  expect_equal(ctrl$maxit, 100)
  expect_equal(ctrl$verbose, 1)
})

test_that("JointODE.control validates inputs", {
  expect_error(JointODE.control(maxit = 0), "positive")
  expect_error(JointODE.control(maxit = -1), "positive")
  expect_error(JointODE.control(tol = 0), "positive")
  expect_error(JointODE.control(parallel = "yes"), "TRUE or FALSE")
  expect_error(JointODE.control(n_cores = -1), "non-negative")
  expect_error(JointODE.control(hazard_quadrature = 0), "positive integer")
  expect_error(JointODE.control(hazard_quadrature = 2.5), "positive integer")
  expect_error(JointODE.control(.list = "bad"), "\\.list must be a list")
})

test_that("JointODE.control verbose coercion", {
  expect_equal(JointODE.control(verbose = TRUE)$verbose, 1)
  expect_equal(JointODE.control(verbose = FALSE)$verbose, 0)
  expect_equal(JointODE.control(verbose = 2)$verbose, 2)
  expect_warning(JointODE.control(verbose = "yes"))
})

test_that("JointODE.control passes through extra params", {
  expect_equal(JointODE.control(custom_param = 42)$custom_param, 42)
})

# ==============================================================================
# Parameter Utilities
# ==============================================================================

# --- .count_params ---

test_that(".count_params counts correctly", {
  params <- list(
    coefficients = list(
      baseline = c(0.1, 0.2, 0.3), hazard = c(0.5, -0.3),
      longitudinal = c(-0.01, 0.02, 0.1), initial_state = c(0, 0),
      measurement_error_sd = 0.5, random_effect_sigma = diag(2)
    ),
    configurations = list(residual = "gaussian")
  )
  expect_equal(.count_params(params), 14)

  params2 <- list(
    coefficients = list(
      baseline = numeric(5), hazard = numeric(3),
      longitudinal = numeric(4), initial_state = c(0, 0),
      measurement_error_sd = 0.1, random_effect_sigma = diag(3)
    ),
    configurations = list(residual = "student_t")
  )
  expect_equal(.count_params(params2), 22)
})

# --- .coef_table ---

test_that(".coef_table structure and computation", {
  tbl <- .coef_table(c(a = 2.0, b = -1.0), c(a = 0.5, b = 0.25))
  expect_equal(ncol(tbl), 4)
  expect_equal(unname(tbl[, "z value"]), c(4.0, -4.0))
  expect_equal(unname(tbl[, "Pr(>|z|)"]), 2 * pnorm(-abs(c(4.0, -4.0))))

  tbl2 <- .coef_table(c(1.0), c(0.0))
  expect_true(is.infinite(tbl2[1, "z value"]))
})


# --- .compute_dimensions ---

test_that(".compute_dimensions returns correct values", {
  dims <- .compute_dimensions(
    .parse_longitudinal_formula(observed ~ x1 + x2),
    .parse_survival_formula(Surv(time, status) ~ w1 + w2),
    list(degree = 2, n_knots = 0)
  )
  expect_equal(dims$n_survival_covariates, 2)
  expect_equal(dims$n_spline_basis, 3)

  dims2 <- .compute_dimensions(
    .parse_longitudinal_formula(
      y ~ biomarker + velocity + (biomarker + velocity | id)
    ),
    .parse_survival_formula(Surv(time, status) ~ 1),
    list(degree = 2, n_knots = 0)
  )
  expect_equal(dims2$n_random_effects, 4)
})

# --- .get_spline_config ---

test_that(".get_spline_config variants", {
  x <- seq(0, 10, length.out = 100)
  c1 <- .get_spline_config(
    x,
    degree = 2, n_knots = 1, knot_placement = "quantile"
  )
  expect_equal(length(c1$knots), 1)
  expect_equal(c1$df, 4)

  c2 <- .get_spline_config(x, degree = 2, n_knots = 3, knot_placement = "equal")
  expect_equal(length(c2$knots), 3)

  c3 <- .get_spline_config(
    seq(1, 5, length.out = 50),
    degree = 2, n_knots = 1, boundary_knots = c(0, 10)
  )
  expect_equal(c3$boundary_knots, c(0, 10))

  expect_error(
    .get_spline_config(1:10, knot_placement = "invalid"),
    "knot_placement"
  )
})
