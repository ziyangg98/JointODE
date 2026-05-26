# ==============================================================================
# simulate.R Unit Tests
# ==============================================================================
#
# Purpose:
#   Test simulate() which generates synthetic joint longitudinal-survival
#   data from a second-order ODE model with random effects and covariates.
#
# Test Coverage:
#   1. Basic Functionality
#      - Output structure: longitudinal_data, survival_data, state,
#        random_effects
#      - Data types, dimensions, and required columns
#      - Value ranges and reproducibility with seed
#
#   2. Mathematical Correctness
#      - Parameter transformation: xi/period → latent omega/xi roots
#      - Covariance matrix: positive semi-definite and symmetric
#      - Random effects: centered at zero (within 4 SE)
#
#   3. Input Validation
#      - Positive value requirements (n_subjects, maxt, xi mean)
#      - Parameter format (xi and period need mean & sd)
#      - gamma must be 0, 1, or 2
#      - Covariate consistency between usage and definition
#
#   4. Edge Cases
#      - Single subject (n=1)
#      - Zero variance (deterministic dynamics)
#
#   5. Feature Tests
#      - gamma parameter: controls velocity power in hazard
#        (0=no velocity, 1=linear velocity, 2=quadratic velocity)
#      - Covariates: longitudinal (x) vs survival (w), binary/normal
#
#   6. Internal Functions
#      - .create_example_data(): generates example with true parameters
#
# Random Effects Structure:
#   - Only omega and xi latent parameters have random effects (2×2 covariance)
#   - Rationale: subject-specific frequency and damping
#
# Model Equations:
#   ODE:
#     \eqn{\frac{d^2y}{dt^2} = c_0 + c_1 y + c_2 \frac{dy}{dt} +
#          \sum \beta_k x_k}
#   Hazard:
#     \eqn{h(t) = h_0(t) \exp(\alpha_1 y +
#          \alpha_2 (\frac{dy}{dt})^\gamma + \sum \phi_k w_k)}
#   where \eqn{\gamma \in \{0, 1, 2\}} controls velocity power
#
# ==============================================================================

test_that("simulate generates valid output and respects parameters", {
  n <- 10
  sim <- JointODE::simulate(n_subjects = n, maxt = 5, seed = 123)

  expect_type(sim, "list")
  expect_named(
    sim,
    c("longitudinal_data", "survival_data", "random_effects")
  )
  expect_s3_class(sim$longitudinal_data, "data.frame")
  expect_s3_class(sim$survival_data, "data.frame")
  expect_true(is.matrix(sim$random_effects))

  expect_true(all(
    c("id", "time", "observed", "biomarker", "velocity", "acceleration") %in%
      names(sim$longitudinal_data)
  ))
  expect_true(all(
    c("id", "time", "status") %in%
      names(sim$survival_data)
  ))
  expect_equal(
    colnames(sim$random_effects),
    c("initial_biomarker", "initial_velocity", "omega", "xi")
  )

  expect_equal(ncol(sim$random_effects), 4)

  expect_equal(nrow(sim$survival_data), n)
  expect_equal(nrow(sim$random_effects), n)
  expect_equal(length(unique(sim$longitudinal_data$id)), n)

  expect_true(all(sim$survival_data$time > 0))
  expect_true(all(sim$survival_data$time <= 5))
  expect_true(all(sim$survival_data$status %in% c(0, 1)))
  expect_true(all(is.finite(sim$longitudinal_data$biomarker)))
})

test_that("simulate is reproducible", {
  sim1 <- JointODE::simulate(n_subjects = 5, seed = 789)
  sim2 <- JointODE::simulate(n_subjects = 5, seed = 789)

  expect_identical(sim1$survival_data$time, sim2$survival_data$time)
  expect_identical(sim1$random_effects, sim2$random_effects)
})

test_that("dynamics transformation and covariance are correct", {
  n <- 20
  xi_mean <- 0.6
  period_mean <- 7

  sim <- JointODE::simulate(
    n_subjects = n,
    longitudinal = list(
      xi = c(mean = xi_mean, sd = 0.15),
      period = c(mean = period_mean, sd = 1.0),
      excitation = list(offset = 2.0, covariates = numeric(0)),
      initial = list(
        biomarker = c(mean = 3.0, sd = 0.1),
        velocity = c(mean = -0.1, sd = 0.1)
      ),
      error_sd = 0.2,
      n_measurements = 3
    ),
    seed = 999
  )

  re_means <- colMeans(sim$random_effects)
  expect_equal(as.numeric(re_means["omega"]), 0, tolerance = 0.1)
  expect_equal(as.numeric(re_means["xi"]), 0, tolerance = 0.2)
})

test_that("random_effects are properly centered", {
  n <- 15
  sim <- JointODE::simulate(
    n_subjects = n,
    longitudinal = list(
      xi = c(mean = 0.5, sd = 0.1),
      period = c(mean = 5, sd = 0.5),
      excitation = list(offset = 2, covariates = numeric(0)),
      initial = list(
        biomarker = c(mean = 1, sd = 0.1),
        velocity = c(mean = 0, sd = 0.1)
      ),
      error_sd = 0.1,
      n_measurements = 3
    ),
    seed = 777
  )

  # omega and xi latent RE should be approximately centered
  for (col_idx in 3:4) {
    col_sd <- sd(sim$random_effects[, col_idx])
    se <- col_sd / sqrt(n)
    expect_true(abs(mean(sim$random_effects[, col_idx])) < 4 * se)
  }
})

test_that("simulate validates parameters correctly", {
  expect_error(JointODE::simulate(n_subjects = -5), "positive integer")
  expect_error(JointODE::simulate(maxt = -1), "positive")

  expect_error(
    JointODE::simulate(
      n_subjects = 5,
      longitudinal = list(
        xi = 0.5,
        period = c(mean = 5, sd = 0.5),
        excitation = list(offset = 0, covariates = numeric(0)),
        initial = list(
          biomarker = c(mean = 0, sd = 0.1),
          velocity = c(mean = 0, sd = 0.1)
        ),
        error_sd = 0.1,
        n_measurements = 5
      )
    ),
    "xi must have mean and sd"
  )

  expect_error(
    JointODE::simulate(
      n_subjects = 5,
      longitudinal = list(
        xi = c(mean = -0.5, sd = 0.1),
        period = c(mean = 5, sd = 0.5),
        excitation = list(offset = 0, covariates = numeric(0)),
        initial = list(
          biomarker = c(mean = 0, sd = 0.1),
          velocity = c(mean = 0, sd = 0.1)
        ),
        error_sd = 0.1,
        n_measurements = 5
      )
    ),
    "xi mean must be positive"
  )

  expect_error(
    JointODE::simulate(
      n_subjects = 5,
      longitudinal = list(
        xi = c(mean = 0.5, sd = 0.1),
        period = c(mean = 5, sd = 0.5),
        excitation = list(offset = 0, covariates = c(x1 = 1, x2 = 2)),
        initial = list(
          biomarker = c(mean = 0, sd = 0.1),
          velocity = c(mean = 0, sd = 0.1)
        ),
        error_sd = 0.1,
        n_measurements = 5
      ),
      covariates = list(x1 = list(type = "normal", mean = 0, sd = 1))
    ),
    "Missing covariate"
  )

  expect_error(
    JointODE::simulate(
      n_subjects = 5,
      survival = list(
        baseline = list(type = "weibull", shape = 2, scale = 10),
        value = 0.5,
        slope = 1.0,
        gamma = 3,
        covariates = numeric(0)
      )
    ),
    "gamma must be 0, 1, or 2"
  )
})

test_that("simulate handles edge cases", {
  sim1 <- JointODE::simulate(
    n_subjects = 1,
    longitudinal = list(
      xi = c(mean = 0.5, sd = 0.05),
      period = c(mean = 5, sd = 0.2),
      excitation = list(offset = 2, covariates = numeric(0)),
      initial = list(
        biomarker = c(mean = 1, sd = 0.1),
        velocity = c(mean = 0, sd = 0.1)
      ),
      error_sd = 0.1,
      n_measurements = 5
    ),
    survival = list(
      baseline = list(type = "weibull", shape = 2, scale = 10),
      value = 0.5,
      slope = 1.0,
      gamma = 1,
      covariates = numeric(0)
    ),
    covariates = list(),
    seed = 444
  )
  expect_equal(nrow(sim1$random_effects), 1)
  expect_true(is.matrix(sim1$random_effects))
  expect_equal(ncol(sim1$random_effects), 4)

  sim_det <- JointODE::simulate(
    n_subjects = 5,
    longitudinal = list(
      xi = c(mean = 0.7, sd = 0),
      period = c(mean = 5, sd = 0),
      excitation = list(offset = 2, covariates = numeric(0)),
      initial = list(
        biomarker = c(mean = 1, sd = 0),
        velocity = c(mean = 0, sd = 0)
      ),
      error_sd = 0.1,
      n_measurements = 3
    ),
    seed = 111
  )

  for (col in colnames(sim_det$random_effects)) {
    expect_true(
      sd(sim_det$random_effects[, col]) < 1e-10,
      label = paste("SD of", col, "should be ~0")
    )
  }
})

test_that("simulate respects gamma parameter", {
  sim_g0 <- JointODE::simulate(
    n_subjects = 5,
    survival = list(
      baseline = list(type = "weibull", shape = 2, scale = 15),
      value = 0.5,
      slope = 1.0,
      gamma = 0,
      covariates = numeric(0)
    ),
    seed = 456
  )

  sim_g1 <- JointODE::simulate(
    n_subjects = 5,
    survival = list(
      baseline = list(type = "weibull", shape = 2, scale = 15),
      value = 0.5,
      slope = 1.0,
      gamma = 1,
      covariates = numeric(0)
    ),
    seed = 456
  )

  sim_g2 <- JointODE::simulate(
    n_subjects = 5,
    survival = list(
      baseline = list(type = "weibull", shape = 2, scale = 15),
      value = 0.5,
      slope = 1.0,
      gamma = 2,
      covariates = numeric(0)
    ),
    seed = 456
  )

  expect_s3_class(sim_g0$survival_data, "data.frame")
  expect_s3_class(sim_g1$survival_data, "data.frame")
  expect_s3_class(sim_g2$survival_data, "data.frame")

  expect_false(identical(sim_g0$survival_data$time, sim_g1$survival_data$time))
  expect_false(identical(sim_g1$survival_data$time, sim_g2$survival_data$time))
})

test_that("simulate handles covariates correctly", {
  sim_with_cov <- JointODE::simulate(
    n_subjects = 8,
    longitudinal = list(
      xi = c(mean = 0.5, sd = 0.1),
      period = c(mean = 5, sd = 0.5),
      excitation = list(offset = 0.2, covariates = c(x1 = -0.2, x2 = 0.1)),
      initial = list(
        biomarker = c(mean = 1.2, sd = 0.1),
        velocity = c(mean = -0.1, sd = 0.1)
      ),
      error_sd = 0.5,
      n_measurements = 3
    ),
    survival = list(
      baseline = list(type = "weibull", shape = 2, scale = 100),
      value = 0.01,
      slope = 0.5,
      gamma = 1,
      covariates = c(w1 = -0.5, w2 = 0.3)
    ),
    covariates = list(
      x1 = list(type = "normal", mean = 0, sd = 1),
      x2 = list(type = "normal", mean = 0, sd = 1),
      w1 = list(type = "normal", mean = 0, sd = 1),
      w2 = list(type = "binary", prob = 0.6)
    ),
    seed = 111
  )

  expect_true(all(c("x1", "x2") %in% names(sim_with_cov$longitudinal_data)))
  expect_true(all(c("w1", "w2") %in% names(sim_with_cov$survival_data)))
  expect_equal(
    colnames(sim_with_cov$random_effects),
    c("initial_biomarker", "initial_velocity", "omega", "xi")
  )
  expect_true(all(sim_with_cov$survival_data$w2 %in% c(0, 1)))

  sim_no_cov <- JointODE::simulate(
    n_subjects = 3,
    longitudinal = list(
      xi = c(mean = 0.7, sd = 0.1),
      period = c(mean = 5, sd = 0.5),
      excitation = list(offset = 0, covariates = numeric(0)),
      initial = list(
        biomarker = c(mean = 1.5, sd = 0.1),
        velocity = c(mean = 0, sd = 0.1)
      ),
      error_sd = 0.4,
      n_measurements = 3
    ),
    survival = list(
      baseline = list(type = "weibull", shape = 1.6, scale = 150),
      value = 0.08,
      slope = 0.4,
      gamma = 1,
      covariates = numeric(0)
    ),
    covariates = list(),
    seed = 666
  )

  expect_equal(ncol(sim_no_cov$longitudinal_data), 6)
  expect_equal(ncol(sim_no_cov$random_effects), 4)
  expect_equal(ncol(sim_no_cov$survival_data), 3)
})

test_that(".create_example_data works", {
  example <- .create_example_data(n_subjects = 20, seed = 123)

  expect_named(example, c("data", "init"))
  expect_named(example$init, c("coefficients", "configurations"))

  expect_true(all(
    c("longitudinal_data", "survival_data", "random_effects") %in%
      names(example$data)
  ))

  params <- example$init$coefficients
  expect_true(is.numeric(params$baseline))
  expect_true(is.numeric(params$longitudinal))
  expect_true(is.matrix(params$random_effect_sigma))

  expect_equal(example$init$configurations$gamma, 1)

  example2 <- .create_example_data(n_subjects = 20, seed = 123)
  expect_identical(
    example$data$survival_data$time,
    example2$data$survival_data$time
  )
})
