# ==============================================================================
# JointODE End-to-End Integration Tests
# ==============================================================================

skip_on_cran()

# Use small dataset for speed
e2e_sim <- simulate(
  n_subjects = 30,
  longitudinal = list(
    xi = c(mean = 0.4, sd = 0.05),
    period = c(mean = 6, sd = 0.1),
    excitation = list(offset = 0, covariates = c(x1 = 0.5)),
    initial = list(
      offset = c(biomarker = -0.5, velocity = -0.1),
      covariates = list(
        biomarker = c(x1 = 0.08),
        velocity = c(x1 = 0.1)
      )
    ),
    error_sd = 0.1,
    n_measurements = 15
  ),
  survival = list(
    baseline = list(type = "weibull", shape = 2.0, scale = 15.0),
    value = 0.8,
    slope = 2.0,
    gamma = 1,
    covariates = c(w1 = 0.6)
  ),
  covariates = list(
    x1 = list(type = "normal", mean = 0, sd = 1),
    w1 = list(type = "normal", mean = 0, sd = 1)
  ),
  seed = 42
)

# Remove reserved columns from longitudinal data
e2e_sim$longitudinal_data$biomarker <- NULL
e2e_sim$longitudinal_data$velocity <- NULL
e2e_sim$longitudinal_data$acceleration <- NULL

test_that("JointODE fits and converges", {
  fit <- JointODE(
    longitudinal_formula = observed ~ biomarker + velocity + x1 +
      (biomarker + velocity | id),
    survival_formula = Surv(time, status) ~ w1,
    longitudinal_data = e2e_sim$longitudinal_data,
    survival_data = e2e_sim$survival_data,
    state = e2e_sim$state,
    control = list(maxit = 30, tol = 1e-3)
  )

  expect_s3_class(fit, "JointODE")
  expect_true(is.finite(fit$logLik))
  expect_true(!is.null(fit$vcov))
  expect_true(is.matrix(fit$vcov))
})

test_that("JointODE S3 methods work on fitted object", {
  fit <- JointODE(
    longitudinal_formula = observed ~ biomarker + velocity + x1 +
      (biomarker + velocity | id),
    survival_formula = Surv(time, status) ~ w1,
    longitudinal_data = e2e_sim$longitudinal_data,
    survival_data = e2e_sim$survival_data,
    state = e2e_sim$state,
    control = list(maxit = 30, tol = 1e-3)
  )

  # coef
  cf <- coef(fit)
  expect_true(is.numeric(cf))
  expect_true(length(cf) > 0)

  # logLik
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))

  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.JointODE")

  # print (no error)
  out <- capture.output(print(fit))
  expect_true(length(out) > 0)

  # predict
  pred <- predict(fit)
  expect_s3_class(pred, "data.frame")
  expect_true(nrow(pred) > 0)
})
