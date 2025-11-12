# Tests for Longitudinal Formula Parsing

# Helper function to create test data for formula parsing tests
formula_test_data <- function(n = 10, n_subjects = 2) {
  df <- data.frame(
    y = rnorm(n),
    id = rep(1:n_subjects, each = n / n_subjects)
  )
  df$x <- rnorm(n)
  df$time <- rep(1:(n / n_subjects), n_subjects)
  df
}


# ==============================================================================
# .parse_longitudinal_formula() tests
# ==============================================================================

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
  expect_equal(result$fixed_terms, c("(Intercept)", "x1", "x2"))
  expect_null(result$random_terms)
  expect_false(result$biomarker$random)
  expect_false(result$velocity$random)
})


test_that("parse formula with reserved words", {
  result <- .parse_longitudinal_formula(y ~ x + (biomarker + velocity | id))

  expect_null(result$random_terms) # No covariate random effects
  expect_true(result$biomarker$random)
  expect_true(result$velocity$random)
  expect_false(result$biomarker$fixed)
  expect_false(result$velocity$fixed)
})


test_that("parse formula with reserved words in fixed effects", {
  result <- .parse_longitudinal_formula(y ~ biomarker + velocity + x + (1 | id))

  expect_true(result$biomarker$fixed)
  expect_true(result$velocity$fixed)
  expect_false(result$biomarker$random)
  expect_false(result$velocity$random)
  # Reserved words excluded
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


test_that("parse random intercept only", {
  result <- .parse_longitudinal_formula(y ~ x + (1 | id))
  expect_equal(result$random_terms, "(Intercept)")
})


test_that("parse random slope without intercept", {
  result <- .parse_longitudinal_formula(y ~ x + (x | id))
  expect_equal(result$random_terms, "x")
})


test_that("parse interaction terms", {
  result <- .parse_longitudinal_formula(y ~ x1 * x2 + (1 | id))
  expect_true("x1:x2" %in% result$fixed_terms)
})


# ==============================================================================
# .validate_longitudinal_formula() tests
# ==============================================================================

test_that("validate accepts valid formula", {
  data <- formula_test_data()
  expect_silent(.validate_longitudinal_formula(y ~ x + (1 | id), data))
})


test_that("validate rejects formula without random effects", {
  data <- formula_test_data()
  expect_error(
    .validate_longitudinal_formula(y ~ x, data),
    "Random effects specification is required"
  )
})


test_that("validate rejects non-formula", {
  data <- formula_test_data()
  expect_error(
    .validate_longitudinal_formula("not formula", data),
    "Formula must be a formula object"
  )
})


test_that("validate rejects one-sided formula", {
  data <- formula_test_data()
  expect_error(
    .validate_longitudinal_formula(~x, data),
    "two-sided"
  )
})


test_that("validate rejects missing variable", {
  data <- formula_test_data()
  expect_error(
    .validate_longitudinal_formula(y ~ missing + (1 | id), data),
    "Variables not found"
  )
})


test_that("validate rejects missing grouping variable", {
  data <- formula_test_data()
  expect_error(
    .validate_longitudinal_formula(y ~ x + (1 | wrong), data),
    "Grouping variable"
  )
})


test_that("validate rejects reserved words in data", {
  data <- data.frame(
    y = 1:10,
    biomarker = 1:10, # Reserved word!
    id = rep(1:2, each = 5)
  )
  expect_error(
    .validate_longitudinal_formula(y ~ 1 + (1 | id), data),
    "Reserved words cannot be used as variable names"
  )
})


test_that("validate accepts reserved words in formula only", {
  data <- formula_test_data()
  expect_silent(
    .validate_longitudinal_formula(y ~ x + (biomarker + velocity | id), data)
  )
})


test_that("validate rejects multiple random groupings", {
  data <- formula_test_data()
  expect_error(
    .validate_longitudinal_formula(y ~ x + (1 | id) + (1 | time), data),
    "Multiple random effects groupings are not supported"
  )
})


test_that("validate rejects nested random effects", {
  data <- data.frame(
    y = 1:20,
    x = 1:20,
    school = rep(1:2, each = 10),
    class = rep(1:4, each = 5)
  )
  expect_error(
    .validate_longitudinal_formula(y ~ x + (1 | school / class), data),
    "Nested random effects are not supported"
  )
})


# ==============================================================================
# Variable modes: fixed only / random only / both
# ==============================================================================

test_that("variable in fixed effects only", {
  result <- .parse_longitudinal_formula(y ~ time + (biomarker | id))
  expect_true("time" %in% result$fixed_terms)
  expect_false("time" %in% result$random_terms)
})


test_that("variable in random effects only", {
  result <- .parse_longitudinal_formula(y ~ x + (time | id))
  expect_false("time" %in% result$fixed_terms)
  expect_true("time" %in% result$random_terms)
})


test_that("variable in both fixed and random", {
  result <- .parse_longitudinal_formula(y ~ time + (time | id))
  expect_true("time" %in% result$fixed_terms)
  expect_true("time" %in% result$random_terms)
})


# ==============================================================================
# .parse_survival_formula() tests
# ==============================================================================

test_that("parse survival formula with covariates", {
  result <- .parse_survival_formula(Surv(time, status) ~ age + treatment)

  expect_equal(result$time_var, "time")
  expect_equal(result$status_var, "status")
  expect_equal(result$covariate_terms, c("age", "treatment"))
})


test_that("parse survival formula without covariates", {
  result <- .parse_survival_formula(Surv(survtime, event) ~ 1)

  expect_equal(result$time_var, "survtime")
  expect_equal(result$status_var, "event")
  expect_null(result$covariate_terms)
})


test_that("parse survival formula with interaction", {
  result <- .parse_survival_formula(Surv(t, d) ~ age * treatment)

  expect_true("age:treatment" %in% result$covariate_terms)
})


test_that("parse survival formula rejects non-Surv LHS", {
  expect_error(
    .parse_survival_formula(time ~ age + treatment),
    "Survival formula must have Surv\\(\\) on the left-hand side"
  )
})


test_that("parse survival formula rejects Surv with insufficient args", {
  expect_error(
    .parse_survival_formula(Surv(time) ~ age),
    "Surv\\(\\) must have at least time and status arguments"
  )
})


# ==============================================================================
# .validate_survival_formula() tests
# ==============================================================================

test_that("validate accepts valid survival formula", {
  data <- data.frame(
    time = c(10, 20, 30),
    status = c(1, 0, 1),
    age = c(50, 60, 70)
  )
  expect_silent(.validate_survival_formula(Surv(time, status) ~ age, data))
})


test_that("validate rejects non-formula", {
  data <- data.frame(time = 1:3, status = c(1, 0, 1))
  expect_error(
    .validate_survival_formula("not formula", data),
    "Formula must be a formula object"
  )
})


test_that("validate rejects one-sided formula", {
  data <- data.frame(time = 1:3, status = c(1, 0, 1))
  expect_error(
    .validate_survival_formula(~age, data),
    "two-sided"
  )
})


test_that("validate rejects missing variables", {
  data <- data.frame(time = 1:3, status = c(1, 0, 1))
  expect_error(
    .validate_survival_formula(Surv(time, status) ~ missing_var, data),
    "Variables not found"
  )
})


test_that("validate rejects missing time variable", {
  data <- data.frame(status = c(1, 0, 1), age = c(50, 60, 70))
  expect_error(
    .validate_survival_formula(Surv(time, status) ~ age, data),
    "Variables not found"
  )
})


test_that("validate rejects non-positive times", {
  data <- data.frame(
    time = c(-1, 0, 10),
    status = c(1, 0, 1),
    age = c(50, 60, 70)
  )
  expect_error(
    .validate_survival_formula(Surv(time, status) ~ age, data),
    "Invalid observation times"
  )
})


test_that("validate rejects invalid status values", {
  data <- data.frame(
    time = c(10, 20, 30),
    status = c(1, 2, 0), # 2 is invalid
    age = c(50, 60, 70)
  )
  expect_error(
    .validate_survival_formula(Surv(time, status) ~ age, data),
    "Invalid status values"
  )
})


test_that("validate rejects missing time values", {
  data <- data.frame(
    time = c(10, NA, 30),
    status = c(1, 0, 1),
    age = c(50, 60, 70)
  )
  expect_error(
    .validate_survival_formula(Surv(time, status) ~ age, data),
    "Missing values found in time variable"
  )
})


test_that("validate rejects missing status values", {
  data <- data.frame(
    time = c(10, 20, 30),
    status = c(1, NA, 1),
    age = c(50, 60, 70)
  )
  expect_error(
    .validate_survival_formula(Surv(time, status) ~ age, data),
    "Missing values found in status variable"
  )
})
