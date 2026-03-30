# ==============================================================================
# JointODE S3 Methods Unit Tests (using mock object)
# ==============================================================================

mock <- create_mock_jointode()

# --- print.JointODE ---

test_that("print.JointODE works", {
  out <- capture.output(result <- print(mock))
  expect_identical(result, mock) # returns invisible(x)
  expect_true(any(grepl("Joint ODE Model", out)))
  expect_true(any(grepl("Log-likelihood", out)))
  expect_true(any(grepl("AIC", out)))
})

# --- summary.JointODE ---

test_that("summary.JointODE returns correct structure", {
  s <- summary(mock)
  expect_s3_class(s, "summary.JointODE")

  expect_true(!is.null(s$coef_baseline))
  expect_true(!is.null(s$coef_survival))
  expect_true(!is.null(s$coef_longitudinal))

  expect_equal(ncol(s$coef_baseline), 4)
  expect_equal(ncol(s$coef_survival), 4)
  expect_equal(ncol(s$coef_longitudinal), 4)

  expect_equal(
    nrow(s$coef_baseline),
    length(mock$parameters$coefficients$baseline)
  )
  expect_equal(
    nrow(s$coef_survival),
    length(mock$parameters$coefficients$hazard)
  )
  expect_equal(
    nrow(s$coef_longitudinal),
    length(mock$parameters$coefficients$longitudinal)
  )
})

test_that("summary.JointODE includes data descriptives", {
  s <- summary(mock)
  expect_true(s$nobs > 0)
  expect_true(s$n_observations > 0)
  expect_true(s$n_events >= 0)
  expect_true(s$event_rate >= 0 && s$event_rate <= 1)
})

test_that("summary.JointODE includes derived params", {
  s <- summary(mock)
  # longitudinal[1] is dyn_value which should be negative -> derived exists
  if (mock$parameters$coefficients$longitudinal[1] < 0) {
    expect_true(!is.null(s$derived_params))
    expect_equal(nrow(s$derived_params), 2) # period and xi
    expect_equal(ncol(s$derived_params), 4)
  }
})

test_that("print.summary.JointODE works", {
  s <- summary(mock)
  out <- capture.output(result <- print(s))
  expect_identical(result, s) # returns invisible(x)
  expect_true(any(grepl("Longitudinal Process", out)))
  expect_true(any(grepl("Survival Process", out)))
  expect_true(any(grepl("Variance Components", out)))
})

# --- coef.JointODE ---

test_that("coef.JointODE returns named vector", {
  cf <- coef(mock)
  expect_true(is.numeric(cf))
  expect_true(!is.null(names(cf)))

  n_expected <- length(mock$parameters$coefficients$baseline) +
    length(mock$parameters$coefficients$hazard) +
    length(mock$parameters$coefficients$longitudinal) +
    length(mock$parameters$coefficients$initial_state)
  expect_equal(length(cf), n_expected)

  # Check prefix naming
  pat <- "^(baseline|hazard|longitudinal|initial state):"
  expect_true(all(grepl(pat, names(cf))))
})

test_that("coef.JointODE preserves values", {
  cf <- coef(mock)
  n_b <- length(mock$parameters$coefficients$baseline)
  expect_equal(
    unname(cf[seq_len(n_b)]),
    unname(mock$parameters$coefficients$baseline)
  )
})

# --- vcov.JointODE ---

test_that("vcov.JointODE returns symmetric matrix", {
  v <- vcov(mock)
  expect_true(is.matrix(v))
  expect_equal(nrow(v), ncol(v))
  expect_equal(v, t(v)) # symmetric
})

# --- logLik.JointODE ---

test_that("logLik.JointODE returns logLik class", {
  ll <- logLik(mock)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), mock$logLik)
  expect_true(!is.null(attr(ll, "df")))
  expect_equal(attr(ll, "df"), .count_params(mock$parameters))
  expect_equal(attr(ll, "nobs"), length(mock$data))
})

# --- predict.JointODE ---

test_that("predict.JointODE returns data.frame", {
  pred <- predict(mock)
  expect_s3_class(pred, "data.frame")

  expected_cols <- c(
    "id", "time", "cumhaz", "survival",
    "biomarker", "velocity", "acceleration"
  )
  expect_true(all(expected_cols %in% names(pred)))
  expect_true(nrow(pred) > 0)
})

test_that("predict.JointODE errors on newdata", {
  expect_error(predict(mock, newdata = data.frame()), "not yet supported")
})

# --- summary.JointODE with NULL vcov ---

test_that("summary.JointODE handles NULL vcov", {
  mock_no_vcov <- mock
  mock_no_vcov$vcov <- NULL
  s <- summary(mock_no_vcov)
  expect_true(all(is.na(s$coef_baseline[, "Std. Error"])))
})
