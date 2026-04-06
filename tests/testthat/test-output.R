# ==============================================================================
# Output Module: S3 Methods, Plot, Predict Helpers
# ==============================================================================

mock <- .create_mock_jointode()

# --- print / summary ---

test_that("print.JointODE works", {
  out <- capture.output(result <- print(mock))
  expect_identical(result, mock)
  expect_true(any(grepl("Joint ODE Model", out)))
  expect_true(any(grepl("Log-likelihood", out)))
})

test_that("summary.JointODE returns correct structure", {
  s <- summary(mock)
  expect_s3_class(s, "summary.JointODE")
  expect_equal(ncol(s$coef_baseline), 4)
  expect_equal(ncol(s$coef_survival), 4)
  expect_equal(ncol(s$coef_longitudinal), 4)
  expect_equal(
    nrow(s$coef_baseline),
    length(mock$parameters$coefficients$baseline)
  )
  expect_true(s$nobs > 0)
  expect_true(s$event_rate >= 0 && s$event_rate <= 1)
})

test_that("summary.JointODE includes derived params", {
  s <- summary(mock)
  if (mock$parameters$coefficients$longitudinal[1] < 0) {
    expect_equal(nrow(s$derived_params), 2)
  }
})

test_that("print.summary.JointODE works", {
  out <- capture.output(result <- print(summary(mock)))
  expect_true(any(grepl("Longitudinal Process", out)))
  expect_true(any(grepl("Survival Process", out)))
})

test_that("summary.JointODE handles NULL vcov", {
  mock_nv <- mock
  mock_nv$vcov <- NULL
  expect_true(all(is.na(summary(mock_nv)$coef_baseline[, "Std. Error"])))
})

# --- coef / vcov / logLik ---

test_that("coef.JointODE returns named vector", {
  cf <- coef(mock)
  pat <- "^(baseline|hazard|longitudinal|initial state):"
  expect_true(all(grepl(pat, names(cf))))
  expect_equal(
    unname(cf[seq_along(mock$parameters$coefficients$baseline)]),
    unname(mock$parameters$coefficients$baseline)
  )
})

test_that("vcov.JointODE returns symmetric matrix", {
  v <- vcov(mock)
  expect_equal(nrow(v), ncol(v))
  expect_equal(v, t(v))
})

test_that("logLik.JointODE returns logLik class", {
  ll <- logLik(mock)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), mock$logLik)
  expect_equal(attr(ll, "df"), .count_params(mock$parameters))
})

# ==============================================================================
# predict
# ==============================================================================

# predict requires tmb_obj (not available in mock); tested in test-fit.R

test_that("predict.JointODE errors on newdata", {
  expect_error(predict(mock, newdata = data.frame()), "not yet supported")
})


# Plot and predict tests require TMB fit object — covered in test-fit.R
