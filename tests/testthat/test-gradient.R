# ==============================================================================
# Test: Gradient Computation
# ==============================================================================
# Purpose: Verify analytical gradients and parallel computation
# Coverage:
#   - Gradient accuracy vs numerical gradients
#   - Edge cases (single subject)
#   - Parallel vs sequential consistency
#   - Dimension checks

test_that("gradient computation comprehensive tests", {
  skip_on_cran() # Skip on CRAN due to long runtime

  # Setup: Load test data once for multiple tests
  test_env <- load_test_data(n_subjects = 20)
  posteriors <- .compute_posteriors(test_env$data, test_env$parameters)

  coefficients <- test_env$parameters$coefficients
  configurations <- test_env$parameters$configurations

  params <- c(
    coefficients$baseline,
    coefficients$hazard,
    as.vector(coefficients$acceleration)
  )
  params <- params + rnorm(length(params), 0, 0.1)

  # Test 1: Analytical vs Numerical gradient
  grad_analytical <- .compute_gradient_joint(
    params = params,
    data_list = test_env$data,
    posteriors = posteriors,
    parameters = test_env$parameters
  )
  names(grad_analytical) <- NULL

  grad_numerical <- numDeriv::grad(
    func = .compute_objective_joint,
    x = params,
    data_list = test_env$data,
    posteriors = posteriors,
    parameters = test_env$parameters,
    method = "Richardson"
  )

  # Compare survival gradient components
  n_hazard <- length(c(coefficients$baseline, coefficients$hazard))
  compare_gradient_component(
    grad_analytical[1:n_hazard],
    grad_numerical[1:n_hazard],
    "Survival gradient",
    tolerance = 1e-2
  )

  # Test 2: Check dimensions
  n_baseline <- length(coefficients$baseline)
  n_hazard_only <- length(coefficients$hazard)
  n_acceleration <- length(coefficients$acceleration)

  expect_equal(
    length(grad_analytical),
    n_baseline + n_hazard_only + n_acceleration
  )
  expect_true(all(is.finite(grad_analytical)))

  # Test 3: Verify gradient structure
  grad_baseline <- grad_analytical[1:n_baseline]
  grad_hazard <- grad_analytical[(n_baseline + 1):(n_baseline + n_hazard_only)]
  grad_acceleration <- grad_analytical[
    (n_baseline + n_hazard_only + 1):length(grad_analytical)
  ]

  expect_equal(length(grad_baseline), n_baseline)
  expect_equal(length(grad_hazard), n_hazard_only)
  expect_equal(length(grad_acceleration), n_acceleration)
})

test_that("gradient handles edge cases", {
  skip_on_cran() # Skip on CRAN due to long runtime

  # Test with single subject
  test_env <- load_test_data(n_subjects = 1)
  posteriors <- .compute_posteriors(test_env$data, test_env$parameters)

  coefficients <- test_env$parameters$coefficients
  configurations <- test_env$parameters$configurations

  params <- c(
    coefficients$baseline,
    coefficients$hazard,
    coefficients$acceleration
  )
  fixed_params <- list(
    measurement_error_sd = coefficients$measurement_error_sd,
    random_effect_sd = coefficients$random_effect_sd
  )

  # Should not error with single subject
  grad <- .compute_gradient_joint(
    params = params,
    data_list = test_env$data,
    posteriors = posteriors,
    configurations = configurations,
    fixed_parameters = fixed_params
  )

  expect_true(all(is.finite(grad)))
  expect_equal(length(grad), length(params))
})

test_that("parallel computation consistency", {
  skip_on_cran() # Skip on CRAN due to long runtime

  # Load test data once
  test_env <- load_test_data(n_subjects = 20)
  posteriors <- .compute_posteriors(test_env$data, test_env$parameters)

  coefficients <- test_env$parameters$coefficients
  configurations <- test_env$parameters$configurations
  params <- c(
    coefficients$baseline,
    coefficients$hazard,
    coefficients$acceleration
  )
  fixed_params <- list(
    measurement_error_sd = coefficients$measurement_error_sd,
    random_effect_sd = coefficients$random_effect_sd
  )

  # Test both gradient and objective in parallel vs sequential

  # Gradient computation
  grad_sequential <- .compute_gradient_joint(
    params = params,
    data_list = test_env$data,
    posteriors = posteriors,
    configurations = configurations,
    fixed_parameters = fixed_params,
    parallel = FALSE
  )

  grad_parallel <- .compute_gradient_joint(
    params = params,
    data_list = test_env$data,
    posteriors = posteriors,
    configurations = configurations,
    fixed_parameters = fixed_params,
    parallel = TRUE,
    n_cores = 2
  )

  expect_equal(
    grad_sequential,
    grad_parallel,
    tolerance = 1e-10,
    label = "Parallel vs Sequential gradient"
  )

  # Objective computation
  obj_sequential <- .compute_objective_joint(
    params = params,
    data_list = test_env$data,
    posteriors = posteriors,
    configurations = configurations,
    fixed_parameters = fixed_params,
    parallel = FALSE
  )

  obj_parallel <- .compute_objective_joint(
    params = params,
    data_list = test_env$data,
    posteriors = posteriors,
    configurations = configurations,
    fixed_parameters = fixed_params,
    parallel = TRUE,
    n_cores = 2
  )

  expect_equal(
    obj_sequential,
    obj_parallel,
    tolerance = 1e-10,
    label = "Parallel vs Sequential objective"
  )

  # Test auto-detection of cores
  grad_parallel_auto <- .compute_gradient_joint(
    params = params,
    data_list = test_env$data,
    posteriors = posteriors,
    configurations = configurations,
    fixed_parameters = fixed_params,
    parallel = TRUE
  )

  expect_equal(
    grad_sequential,
    grad_parallel_auto,
    tolerance = 1e-10,
    label = "Parallel (auto cores) vs Sequential gradient"
  )
})
