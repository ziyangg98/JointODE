test_that("simulate generates valid data structure", {
  sim <- JointODE::simulate(n_subjects = 10, seed = 123)

  # Check output structure
  expect_type(sim, "list")
  expect_named(sim, c("longitudinal_data", "survival_data", "state"))

  # Check data frames
  expect_s3_class(sim$longitudinal_data, "data.frame")
  expect_s3_class(sim$survival_data, "data.frame")
  expect_s3_class(sim$state, "data.frame")

  # Check required columns
  expect_true(all(
    c("id", "time", "observed", "biomarker", "velocity", "acceleration") %in%
      names(sim$longitudinal_data)
  ))
  expect_true(all(c("id", "time", "status", "b") %in% names(sim$survival_data)))
  expect_named(sim$state, c("biomarker", "velocity"))
})

test_that("simulate respects n_subjects parameter", {
  n <- 10
  sim <- JointODE::simulate(n_subjects = n, seed = 456)

  expect_equal(nrow(sim$survival_data), n)
  expect_equal(length(unique(sim$longitudinal_data$id)), n)
  expect_equal(nrow(sim$state), n)
})

test_that("simulate respects seed for reproducibility", {
  sim1 <- JointODE::simulate(n_subjects = 10, seed = 789)
  sim2 <- JointODE::simulate(n_subjects = 10, seed = 789)

  expect_identical(sim1$survival_data$time, sim2$survival_data$time)
  expect_identical(
    sim1$longitudinal_data$observed,
    sim2$longitudinal_data$observed
  )
})

test_that("simulate handles custom parameters correctly", {
  sim <- JointODE::simulate(
    n_subjects = 10,
    shared_sd = 1.2,
    longitudinal = list(
      xi = 0.5,
      period = 3,
      excitation = list(
        offset = 0.2,
        covariates = c(x1 = -0.2, x2 = 0.1)
      ),
      initial = list(
        offset = 1.2,
        covariates = c(x1 = -0.4, x2 = 0.2),
        random_coef = 1.0
      ),
      error_sd = 0.5,
      n_measurements = 10
    ),
    survival = list(
      baseline = list(type = "weibull", shape = 2, scale = 100),
      value = 0.01,
      slope = 0.5,
      covariates = c(w1 = -0.5, w2 = 0.3)
    ),
    covariates = list(
      x1 = list(type = "normal", mean = 0, sd = 1),
      x2 = list(type = "normal", mean = 0, sd = 1),
      w1 = list(type = "normal", mean = 0, sd = 1),
      w2 = list(type = "binary", prob = 0.6)
    ),
    maxt = 50,
    seed = 111
  )

  expect_equal(ncol(sim$longitudinal_data), 10)  # Now includes xi and period
  expect_true(all(sim$survival_data$time <= 50))
  expect_true(
    mean(sim$survival_data$w2) > 0.4 && mean(sim$survival_data$w2) < 0.8
  )
})

test_that("simulate validates input parameters", {
  # Test only key validation cases to speed up tests
  expect_error(JointODE::simulate(n_subjects = -5), "positive integer")
  expect_error(JointODE::simulate(n_subjects = 10, shared_sd = -1), "positive")
  expect_error(
    JointODE::simulate(n_subjects = 10, longitudinal = "not_a_list"),
    "list"
  )
  expect_error(
    JointODE::simulate(
      n_subjects = 10,
      survival = list(baseline = list(type = "unknown"))
    ),
    "weibull"
  )
  # Dimension mismatch
  expect_error(
    JointODE::simulate(
      n_subjects = 10,
      longitudinal = list(
        xi = 0.5,
        period = 5,
          excitation = list(
          offset = 0,
          covariates = c(x1 = 1, x2 = 2) # 2 covariates
        ),
        initial = list(
          offset = 1.0,
          covariates = c(x1 = 1), # Only 1 covariate - mismatch!
          random_coef = 0.5
        ),
        error_sd = 0.2,
        n_measurements = 10
      ),
      covariates = list(
        x1 = list(type = "normal", mean = 0, sd = 1)
        # Missing x2 definition
      )
    ),
    "Missing covariate"
  )
})

test_that("simulate generates reasonable biomarker trajectories", {
  sim <- JointODE::simulate(n_subjects = 50, seed = 222)

  # Check biomarker values are finite
  expect_true(all(is.finite(sim$longitudinal_data$biomarker)))
  expect_true(all(is.finite(sim$longitudinal_data$observed)))

  # Check velocities and accelerations are computed
  expect_true(all(is.finite(sim$longitudinal_data$velocity)))
  expect_true(all(is.finite(sim$longitudinal_data$acceleration)))

  # Check measurement error adds noise
  residuals <- sim$longitudinal_data$observed - sim$longitudinal_data$biomarker
  expect_true(sd(residuals) > 0)
})

test_that("simulate generates valid survival times", {
  sim <- JointODE::simulate(n_subjects = 50, maxt = 5, seed = 333)

  # All times should be positive and within maxt
  expect_true(all(sim$survival_data$time > 0))
  expect_true(all(sim$survival_data$time <= 5))

  # Status should be 0 or 1
  expect_true(all(sim$survival_data$status %in% c(0, 1)))

  # Should have at least some events (censoring is not guaranteed)
  expect_true(any(sim$survival_data$status == 1))

  # Test with parameters that ensure censoring
  sim_cens <- JointODE::simulate(
    n_subjects = 50,
    maxt = 2,
    survival = list(
      baseline = list(type = "weibull", shape = 1.5, scale = 20),
      value = 0.1,
      slope = 0.2,
      covariates = c(w1 = 0.1, w2 = -0.1)
    ),
    seed = 333
  )
  # With these parameters, we should have some censoring
  expect_true(any(sim_cens$survival_data$status == 0))
})

test_that("simulate handles patient-specific dynamics", {
  # Test multiple dynamics groups
  sim_groups <- JointODE::simulate(
    n_subjects = 100,
    longitudinal = list(
      xi = c(0.3, 0.707, 1.0, 2.0),
      period = c(2, 3, 5, 10),
      prob = c(0.4, 0.3, 0.2, 0.1),  # Different probabilities
      excitation = list(
        offset = 0,
        covariates = c(x1 = 0.5)
      ),
      initial = list(
        offset = 0,
        covariates = c(x1 = 0.2),
        random_coef = 0.5
      ),
      n_measurements = 10,
      error_sd = 0.1
    ),
    seed = 999
  )

  # Check that dynamics columns exist
  expect_true("xi" %in% names(sim_groups$longitudinal_data))
  expect_true("period" %in% names(sim_groups$longitudinal_data))

  # Extract unique dynamics per patient
  patient_dynamics <- sim_groups$longitudinal_data[
    !duplicated(sim_groups$longitudinal_data$id),
    c("id", "xi", "period")
  ]

  # Check that all expected xi values are present
  unique_xi <- unique(patient_dynamics$xi)
  expect_true(all(unique_xi %in% c(0.3, 0.707, 1.0, 2.0)))

  # Check that xi and period are paired correctly
  for(i in seq_len(nrow(patient_dynamics))) {
    xi_val <- patient_dynamics$xi[i]
    period_val <- patient_dynamics$period[i]
    if(xi_val == 0.3) expect_equal(period_val, 2)
    if(xi_val == 0.707) expect_equal(period_val, 3)
    if(xi_val == 1.0) expect_equal(period_val, 5)
    if(xi_val == 2.0) expect_equal(period_val, 10)
  }

  # Check probability distribution (with tolerance for sampling)
  xi_table <- table(patient_dynamics$xi)
  xi_props <- as.numeric(xi_table) / sum(xi_table)
  # With 100 subjects, we expect roughly 40%, 30%, 20%, 10%
  expect_true(abs(xi_props[names(xi_table) == "0.3"] - 0.4) < 0.15)
})

test_that("simulate handles single dynamics (backward compatibility)", {
  # Test with single xi and period values
  sim_single <- JointODE::simulate(
    n_subjects = 50,
    longitudinal = list(
      xi = 0.5,
      period = 4,
      excitation = list(
        offset = 0,
        covariates = numeric(0)
      ),
      initial = list(
        offset = -1,
        covariates = numeric(0),
        random_coef = 0.5
      ),
      n_measurements = 10,
      error_sd = 0.1
    ),
    seed = 888
  )

  # All patients should have the same dynamics
  patient_dynamics <- sim_single$longitudinal_data[
    !duplicated(sim_single$longitudinal_data$id),
    c("xi", "period")
  ]

  expect_true(all(patient_dynamics$xi == 0.5))
  expect_true(all(patient_dynamics$period == 4))
})

test_that("simulate validates dynamics parameters correctly", {
  # xi and period must have same length
  expect_error(
    JointODE::simulate(
      n_subjects = 10,
      longitudinal = list(
        xi = c(0.3, 0.7),
        period = c(2, 3, 4),  # Different length!
        excitation = list(offset = 0, covariates = numeric(0)),
        initial = list(offset = 0, covariates = numeric(0), random_coef = 0),
        n_measurements = 10,
        error_sd = 0.1
      )
    ),
    "same length"
  )

  # prob must match xi/period length if provided
  expect_error(
    JointODE::simulate(
      n_subjects = 10,
      longitudinal = list(
        xi = c(0.3, 0.7),
        period = c(2, 3),
        prob = c(0.5, 0.3, 0.2),  # Wrong length!
        excitation = list(offset = 0, covariates = numeric(0)),
        initial = list(offset = 0, covariates = numeric(0), random_coef = 0),
        n_measurements = 10,
        error_sd = 0.1
      )
    ),
    "prob must match"
  )

  # prob must sum to 1
  expect_error(
    JointODE::simulate(
      n_subjects = 10,
      longitudinal = list(
        xi = c(0.3, 0.7),
        period = c(2, 3),
        prob = c(0.3, 0.3),  # Doesn't sum to 1!
        excitation = list(offset = 0, covariates = numeric(0)),
        initial = list(offset = 0, covariates = numeric(0), random_coef = 0),
        n_measurements = 10,
        error_sd = 0.1
      )
    ),
    "prob must match"
  )
})

test_that("simulate handles edge cases", {
  # Small sample size
  sim_small <- JointODE::simulate(n_subjects = 1, seed = 444)
  expect_equal(nrow(sim_small$survival_data), 1)

  # Large error variance
  sim_noisy <- JointODE::simulate(
    n_subjects = 10,
    longitudinal = list(
      xi = 0.7,
      period = 5,
      excitation = list(
        offset = 0,
        covariates = numeric(0)
      ),
      initial = list(
        offset = 1.5,
        covariates = numeric(0),
        random_coef = 1.0
      ),
      error_sd = 10,
      n_measurements = 10
    ),
    survival = list(
      baseline = list(type = "weibull", shape = 1.5, scale = 8),
      value = 0.3,
      slope = 0.7,
      covariates = numeric(0)
    ),
    covariates = list(
      w1 = list(type = "normal", mean = 0, sd = 1),
      w2 = list(type = "binary", prob = 0.5)
    ),
    seed = 555
  )
  expect_true(sd(sim_noisy$longitudinal_data$observed) > 5)

  # No covariates
  sim_no_cov <- JointODE::simulate(
    n_subjects = 10,
    longitudinal = list(
      xi = 0.7,
      period = 5,
      excitation = list(
        offset = 0,
        covariates = numeric(0)
      ),
      initial = list(
        offset = 1.5,
        covariates = numeric(0),
        random_coef = 1.0
      ),
      error_sd = 0.4,
      n_measurements = 10
    ),
    survival = list(
      baseline = list(type = "weibull", shape = 1.6, scale = 150),
      value = 0.08,
      slope = 0.4,
      covariates = numeric(0)
    ),
    covariates = list(),
    seed = 666
  )
  expect_equal(ncol(sim_no_cov$longitudinal_data), 8) # Now includes xi and period
})

test_that("simulate creates consistent longitudinal observations", {
  skip_on_cran() # Skip on CRAN due to loop checks
  sim <- JointODE::simulate(n_subjects = 10, seed = 777)

  # Each subject should have observations up to their event/censoring time
  for (id in unique(sim$longitudinal_data$id)) {
    long_times <- sim$longitudinal_data$time[sim$longitudinal_data$id == id]
    surv_time <- sim$survival_data$time[sim$survival_data$id == id]

    expect_true(max(long_times) <= surv_time + 1e-6)
    expect_true(min(long_times) >= 0)

    # Times should be ordered
    expect_equal(long_times, sort(long_times))
  }
})

test_that("simulate respects covariate distributions", {
  skip_on_cran() # Skip on CRAN due to large sample size
  n <- 100 # Further reduced for faster tests
  sim <- JointODE::simulate(
    n_subjects = n,
    longitudinal = list(
      xi = 0.7,
      period = 5,
      excitation = list(
        offset = 0,
        covariates = c(x1 = -0.2, x2 = 0.1)
      ),
      initial = list(
        offset = 1.5,
        covariates = c(x1 = -0.5, x2 = 0.3),
        random_coef = 1.0
      ),
      error_sd = 0.4,
      n_measurements = 10
    ),
    survival = list(
      baseline = list(type = "weibull", shape = 1.6, scale = 150),
      value = 0.08,
      slope = 0.4,
      covariates = c(w1 = -0.7, w2 = 0.3)
    ),
    covariates = list(
      x1 = list(type = "normal", mean = 0, sd = 1),
      x2 = list(type = "normal", mean = 0, sd = 1),
      w1 = list(type = "normal", mean = 0, sd = 1),
      w2 = list(type = "binary", prob = 0.3)
    ),
    seed = 888
  )

  # Binary variable
  prop <- mean(sim$survival_data$w2)
  expect_true(abs(prop - 0.3) < 0.15) # Wider tolerance for smaller sample

  # Normal variable (w1)
  expect_true(abs(mean(sim$survival_data$w1) - 0) < 0.4)
  expect_true(abs(sd(sim$survival_data$w1) - 1) < 0.4)
})

# Tests for .create_example_data
test_that(".create_example_data works correctly", {
  skip_on_cran() # Skip on CRAN - internal function test
  # Single call to create example data
  example <- JointODE:::.create_example_data(n_subjects = 50, seed = 123)

  # Test structure
  expect_type(example, "list")
  expect_named(example, c("data", "init"))
  expect_type(example$data, "list")
  expect_type(example$init, "list")
  expect_named(example$init, c("coefficients", "configurations"))

  # Test data generation
  expect_true(all(
    c("longitudinal_data", "survival_data", "state") %in%
      names(example$data)
  ))
  expect_equal(nrow(example$data$survival_data), 50)
  expect_equal(nrow(example$data$state), 50)

  # Test coefficients
  coef <- example$init$coefficients
  expect_true(all(
    c(
      "baseline",
      "acceleration",
      "hazard",
      "measurement_error_sd",
      "random_effect_sd"
    ) %in%
      names(coef)
  ))
  expect_type(coef$baseline, "double")
  expect_type(coef$acceleration, "double")
  expect_type(coef$hazard, "double")
  expect_equal(coef$measurement_error_sd, 0.1)
  expect_equal(coef$random_effect_sd, 0.1)

  # Test configurations
  config <- example$init$configurations
  expect_named(config, c("baseline", "autonomous"))
  expect_type(config$baseline, "list")
  expect_true(config$autonomous)
  expect_true(all(
    c("degree", "knots", "boundary_knots", "df") %in% names(config$baseline)
  ))

  # Test reproducibility
  example2 <- JointODE:::.create_example_data(n_subjects = 50, seed = 123)
  expect_identical(
    example$data$survival_data$time,
    example2$data$survival_data$time
  )
  expect_identical(example$init$coefficients, example2$init$coefficients)
})
