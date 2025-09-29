test_that("IndividualODE works with new interface", {
  # Generate test data
  set.seed(123)
  sim_data <- simulate(
    n_subjects = 5,
    longitudinal = list(
      xi = 0.707,
      period = 5,
      excitation = list(
        offset = 1.0,
        covariates = c(x1 = 0.5, x2 = -0.3)
      ),
      initial = list(
        offset = -2,
        covariates = c(x1 = 0.2, x2 = -0.1),
        random_coef = 0.1
      ),
      n_measurements = 20,
      error_sd = 0.05
    ),
    seed = 123
  )

  # Test single subject with covariates
  single_subject_data <- sim_data$longitudinal_data[
    sim_data$longitudinal_data$id == 1,
  ]
  fit1 <- IndividualODE(
    formula = observed ~ x1 + x2,
    data = single_subject_data
  )

  # Check structure
  expect_type(fit1, "list")
  expect_equal(fit1$subject_id, 1)
  expect_true("coefficients" %in% names(fit1))
  expect_true("initial_state" %in% names(fit1))
  expect_true("fitted_values" %in% names(fit1))
  expect_true("residuals" %in% names(fit1))
  expect_true("sse" %in% names(fit1))
  expect_true("convergence" %in% names(fit1))

  # Check coefficients
  expect_true("value" %in% names(fit1$coefficients))
  expect_true("slope" %in% names(fit1$coefficients))
  expect_true("covariates" %in% names(fit1$coefficients))
  expect_true(length(fit1$coefficients$covariates) == 3) # intercept, x1 and x2

  # Check fitted values match observations
  n_obs_subject1 <- sum(sim_data$longitudinal_data$id == 1)
  expect_equal(length(fit1$fitted_values), n_obs_subject1)
  expect_equal(length(fit1$residuals), n_obs_subject1)

  # Test without covariates (intercept only)
  single_subject_data2 <- sim_data$longitudinal_data[
    sim_data$longitudinal_data$id == 2,
  ]
  fit2 <- IndividualODE(
    formula = observed ~ 1,
    data = single_subject_data2
  )

  expect_equal(length(fit2$coefficients$covariates), 1) # Only intercept
  expect_true(fit2$convergence)

  # Test all subjects
  fits_all <- IndividualODE(
    formula = observed ~ x1 + x2,
    data = sim_data$longitudinal_data
  )

  expect_type(fits_all, "list")
  expect_equal(length(fits_all), 5)
  expect_equal(names(fits_all), as.character(1:5))
})


test_that("IndividualODE uses L-BFGS-B optimization", {
  # Simple test data
  set.seed(456)
  test_data <- data.frame(
    id = rep(1, 15),
    time = seq(0, 7, by = 0.5),
    biomarker = -2 * exp(-0.3 * seq(0, 7, by = 0.5)) + rnorm(15, 0, 0.1)
  )

  # Test default optimization
  fit <- IndividualODE(
    formula = biomarker ~ 1,
    data = test_data
  )
  expect_true(fit$convergence)
  expect_type(fit$coefficients$value, "double")
  expect_type(fit$coefficients$slope, "double")
})


test_that("IndividualODE handles matrix input", {
  # Create matrix data
  n_obs <- 30
  test_matrix <- matrix(
    c(
      rep(1:2, each = 15), # id
      rep(seq(0, 7, by = 0.5), 2), # time
      rnorm(30, -2, 0.5), # biomarker
      rnorm(30, 0, 1)
    ), # x1
    ncol = 4
  )
  colnames(test_matrix) <- c("id", "time", "biomarker", "x1")

  # Filter matrix to single subject
  test_matrix_single <- test_matrix[test_matrix[, "id"] == 1, ]
  fit <- IndividualODE(
    formula = biomarker ~ x1,
    data = test_matrix_single
  )

  expect_type(fit, "list")
  expect_true("value" %in% names(fit$coefficients))
  expect_true(length(fit$coefficients$covariates) == 2) # intercept and x1
  expect_true("x1" %in% names(fit$coefficients$covariates))
})


test_that("IndividualODE error handling", {
  test_data <- data.frame(
    id = 1:10,
    time = 1:10,
    y = rnorm(10)
  )

  # Missing response variable
  expect_error(
    IndividualODE(
      formula = biomarker ~ 1,
      data = test_data
    )
  )

  # Missing time column
  test_no_time <- data.frame(
    id = 1:10,
    y = rnorm(10)
  )
  expect_error(
    IndividualODE(
      formula = y ~ 1,
      data = test_no_time
    ),
    "Data must contain a time column"
  )

  # Missing id column
  test_no_id <- data.frame(
    time = 1:10,
    y = rnorm(10)
  )
  expect_error(
    IndividualODE(
      formula = y ~ 1,
      data = test_no_id
    ),
    "Data must contain an id column"
  )

  # Empty data
  empty_data <- data.frame(
    id = integer(0),
    time = numeric(0),
    y = numeric(0)
  )
  expect_error(
    IndividualODE(
      formula = y ~ 1,
      data = empty_data
    )
  )
})


test_that("IndividualODE with custom initial values", {
  # Generate simple data
  set.seed(789)
  test_data <- data.frame(
    id = rep(1, 20),
    time = seq(0, 10, length.out = 20),
    y = -exp(-0.5 * seq(0, 10, length.out = 20)) + rnorm(20, 0, 0.05)
  )

  # Fit with default initial values
  fit_default <- IndividualODE(
    formula = y ~ 1,
    data = test_data
  )

  # Fit with custom initial state
  state_matrix <- matrix(c(-1.0, 0.5), nrow = 1)
  fit_custom <- IndividualODE(
    formula = y ~ 1,
    data = test_data,
    state = state_matrix
  )

  # Both should converge
  expect_true(fit_default$convergence)
  expect_true(fit_custom$convergence)

  # Custom initial values might lead to different (but valid) solution
  expect_type(fit_custom$coefficients$value, "double")
  expect_type(fit_custom$coefficients$slope, "double")
})
