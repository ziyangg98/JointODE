# ==============================================================================
# End-to-end model fitting tests
# ==============================================================================

test_that("MarginalODE basic fit converges", {
  fit <- MarginalODE(
    formula = observed ~ x1 + x2,
    data = sim$data$longitudinal_data,
    control = list(maxit = 100, verbose = 0)
  )
  expect_true(fit$convergence$converged)
  expect_true(is.finite(fit$logLik))
  # (Intercept), x1, x2, init_biomarker, init_velocity
  expect_length(coef(fit), 5L)
})

test_that("MarginalODE with dynamics RE converges", {
  ids30 <- unique(sim$data$longitudinal_data$id)[1:30]
  sub_data <- sim$data$longitudinal_data[
    sim$data$longitudinal_data$id %in% ids30,
  ]
  fit <- MarginalODE(
    formula = observed ~ biomarker + velocity + x1 + x2 +
      (biomarker + velocity | id),
    data = sub_data,
    control = list(maxit = 200, verbose = 0)
  )
  expect_true(is.finite(fit$logLik))
  expect_length(coef(fit), 7L)
})

test_that("MarginalODE estimates Student-t residual df", {
  ids10 <- unique(sim$data$longitudinal_data$id)[1:10]
  sub_data <- sim$data$longitudinal_data[
    sim$data$longitudinal_data$id %in% ids10,
  ]
  fit <- MarginalODE(
    formula = observed ~ x1 + x2,
    data = sub_data,
    residual = "student_t",
    control = list(maxit = 100, verbose = 0)
  )
  expect_true(is.finite(fit$logLik))
  expect_equal(fit$residual$family, "student_t")
  expect_true(is.finite(fit$residual$nu))
  expect_gt(fit$residual$nu, 2)
  expect_true(is.finite(summary(fit)$residual$nu))
})

test_that("MarginalODE Student-t works with dynamics RE warm start", {
  ids8 <- unique(sim$data$longitudinal_data$id)[1:8]
  sub_data <- sim$data$longitudinal_data[
    sim$data$longitudinal_data$id %in% ids8,
  ]
  fit <- MarginalODE(
    formula = observed ~ biomarker + velocity + x1 + x2 + (biomarker | id),
    data = sub_data,
    residual = "student_t",
    control = list(maxit = 50, verbose = 0)
  )
  expect_true(is.finite(fit$logLik))
  expect_equal(fit$residual$family, "student_t")
  expect_true(is.finite(fit$residual$nu))
  expect_gt(attr(logLik(fit), "df"), length(fit$parameters))
})

test_that("JointODE with sim$init converges", {
  ld <- sim$data$longitudinal_data[
    , c("id", "time", "observed", "x1", "x2")
  ]
  ids20 <- unique(ld$id)[1:20]
  fit <- JointODE(
    longitudinal_formula = observed ~ biomarker + velocity +
      x1 + x2 + (biomarker + velocity | id),
    survival_formula = Surv(time, status) ~ w1 + w2,
    longitudinal_data = ld[ld$id %in% ids20, ],
    survival_data = sim$data$survival_data[
      sim$data$survival_data$id %in% ids20,
    ],
    init = sim$init,
    control = list(maxit = 30, verbose = 0)
  )
  expect_true(is.finite(fit$logLik))
  expect_true(fit$cindex > 0.4)
})

test_that("JointODE estimates Student-t residual df", {
  ld <- sim$data$longitudinal_data[
    , c("id", "time", "observed", "x1", "x2")
  ]
  ids <- unique(ld$id)[1:50]
  fit <- JointODE(
    longitudinal_formula = observed ~ biomarker + velocity +
      x1 + x2 + (biomarker | id),
    survival_formula = Surv(time, status) ~ w1 + w2,
    longitudinal_data = ld[ld$id %in% ids, ],
    survival_data = sim$data$survival_data[
      sim$data$survival_data$id %in% ids,
    ],
    init = "marginal",
    residual = "student_t",
    control = list(maxit = 20, verbose = 0)
  )
  expect_true(is.finite(fit$logLik))
  expect_equal(fit$residual$family, "student_t")
  expect_true(is.finite(fit$residual$nu))
  expect_gt(fit$residual$nu, 2)
})

test_that("coef returns named vector with correct length", {
  fit <- MarginalODE(
    formula = observed ~ x1 + x2,
    data = sim$data$longitudinal_data,
    control = list(maxit = 50, verbose = 0)
  )
  co <- coef(fit)
  expect_true(is.numeric(co))
  expect_true(!is.null(names(co)))
  expect_length(co, 5L)
})

test_that("predict returns expected columns", {
  fit <- MarginalODE(
    formula = observed ~ x1 + x2,
    data = sim$data$longitudinal_data,
    control = list(maxit = 50, verbose = 0)
  )
  pred <- predict(fit)
  expect_s3_class(pred, "data.frame")
  expect_true(all(
    c("id", "time", "biomarker", "velocity") %in% names(pred)
  ))
  expect_true(nrow(pred) > 0)
})
