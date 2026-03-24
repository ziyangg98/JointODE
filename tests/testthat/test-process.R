# Essential tests for .process function

test_that(".process structures data correctly", {
  # Create test data
  long_data <- data.frame(
    id = rep(1:2, each = 3),
    obstime = rep(0:2, 2),
    v = rnorm(6),
    x1 = rnorm(6)
  )

  surv_data <- data.frame(
    id = 1:2,
    obstime = c(3, 4),
    status = c(1, 0),
    w1 = rnorm(2)
  )

  result <- .process(
    longitudinal_formula = v ~ x1 + (1 | id),
    survival_formula = Surv(obstime, status) ~ w1,
    longitudinal_data = long_data,
    survival_data = surv_data,
    state = NULL
  )

  # Check structure
  expect_equal(length(result), 2)
  expect_equal(names(result), c("1", "2"))

  # Check subject 1
  s1 <- result[[1]]
  expect_equal(s1$id, 1)
  expect_equal(s1$time, 3)
  expect_equal(s1$status, 1)
  expect_equal(as.numeric(s1$covariates), as.numeric(surv_data[1, "w1"]))
  m0 <- long_data$v[1]
  v0 <- (long_data$v[2] - long_data$v[1]) / (1 - 0)
  expect_equal(as.numeric(s1$initial_state), c(m0, v0))
  expect_equal(length(s1$longitudinal$times), 3)
  expect_equal(s1$longitudinal$times, 0:2)
  expect_equal(as.numeric(s1$longitudinal$measurements), long_data$v[1:3])
  # Check fixed effects covariates (intercept + x1)
  expect_equal(
    as.matrix(s1$longitudinal$covariates$fixed, dimnames = NULL),
    cbind(1, long_data$x1[1:3]),
    ignore_attr = TRUE
  )
  # Check random effects covariates (only intercept for (1 | id))
  expect_equal(
    as.matrix(s1$longitudinal$covariates$random, dimnames = NULL),
    matrix(1, nrow = 3, ncol = 1),
    ignore_attr = TRUE
  )
})

test_that(".process handles state parameter", {
  long_data <- data.frame(
    id = 1,
    obstime = 0,
    v = 1
  )

  surv_data <- data.frame(
    id = 1,
    obstime = 1,
    status = 1
  )

  # With state matrix
  state <- matrix(c(2, 3), nrow = 1)
  result <- .process(
    longitudinal_formula = v ~ 1 + (1 | id),
    survival_formula = Surv(obstime, status) ~ 1,
    longitudinal_data = long_data,
    survival_data = surv_data,
    state = state
  )
  expect_equal(result[[1]]$initial_state, c(2, 3))
})

test_that(".process handles missing longitudinal data", {
  # Subject 2 has no longitudinal data
  long_data <- data.frame(
    id = 1,
    obstime = 0,
    v = 1
  )

  surv_data <- data.frame(
    id = 1:2,
    obstime = c(1, 2),
    status = c(1, 0)
  )

  result <- .process(
    longitudinal_formula = v ~ 1 + (1 | id),
    survival_formula = Surv(obstime, status) ~ 1,
    longitudinal_data = long_data,
    survival_data = surv_data,
    state = NULL
  )

  # Subject 2 should have no observations
  expect_equal(length(result[[2]]$longitudinal$times), 0)
})

test_that(".process handles unordered time data", {
  # Longitudinal data not in time order
  long_data <- data.frame(
    id = c(1, 1, 1, 2, 2, 2),
    obstime = c(2, 0, 1, 1, 2, 0), # Intentionally out of order
    v = c(3, 1, 2, 5, 6, 4),
    x1 = c(30, 10, 20, 50, 60, 40)
  )

  surv_data <- data.frame(
    id = 1:2,
    obstime = c(3, 3),
    status = c(1, 0)
  )

  result <- .process(
    longitudinal_formula = v ~ x1 + (1 | id),
    survival_formula = Surv(obstime, status) ~ 1,
    longitudinal_data = long_data,
    survival_data = surv_data,
    state = NULL
  )

  # Check that data is correctly ordered by time
  expect_equal(result[[1]]$longitudinal$times, 0:2)
  expect_equal(as.numeric(result[[1]]$longitudinal$measurements), c(1, 2, 3))
  expect_equal(result[[2]]$longitudinal$times, 0:2)
  expect_equal(as.numeric(result[[2]]$longitudinal$measurements), c(4, 5, 6))
})
