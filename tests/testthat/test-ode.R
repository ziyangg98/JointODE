# Comprehensive tests for .solve_batch_ode function
# Tests basic mode, forward sensitivity, custom times, output structure,
# and accuracy against true simulation values

# ==========================================================================
# Setup: Process data once for all tests
# ==========================================================================
data_processed <- .process(
  longitudinal_data = sim$data$longitudinal_data[, c(
    "id",
    "time",
    "observed",
    "x1",
    "x2"
  )],
  longitudinal_formula = observed ~
    biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
  survival_data = sim$data$survival_data,
  survival_formula = Surv(time, status) ~ w1 + w2,
  state = as.matrix(sim$data$state)
)

parameters <- sim$init
random_effects <- sim$data$random_effects
n_subjects <- length(data_processed)

# ==========================================================================
# Test 1: CppAD basic mode - structure, sanity checks, and accuracy
# ==========================================================================

test_that(".solve_batch_ode_cppad - basic mode structure and accuracy", {
  # Compute batch ODE using CppAD
  batch_cppad <- .solve_batch_ode_cppad(
    data_list = data_processed,
    random_effects = random_effects,
    parameters = parameters
  )

  # Check that we got results for all subjects
  expect_equal(length(batch_cppad), n_subjects)

  # Check output structure for first subject
  result_basic <- batch_cppad[[1]]

  basic_fields <- c(
    "times",
    "cum_hazard",
    "biomarker",
    "velocity",
    "acceleration",
    "log_hazard"
  )

  for (field in basic_fields) {
    expect_true(
      field %in% names(result_basic),
      label = sprintf("Has field '%s'", field)
    )
  }

  # Check that times vector matches trajectory lengths
  expect_equal(
    length(result_basic$times),
    length(result_basic$biomarker),
    label = "times length matches biomarker length"
  )

  # Should NOT have sensitivity fields in basic mode
  expect_false(
    "dcumhazard_deta_at_event" %in% names(result_basic),
    label = "Should not have sensitivity fields"
  )

  # Check all results have finite values
  for (i in seq_len(n_subjects)) {
    expect_true(
      all(is.finite(batch_cppad[[i]]$cum_hazard)),
      label = sprintf("Subject %d: finite cum_hazard", i)
    )
    expect_true(
      all(is.finite(batch_cppad[[i]]$log_hazard)),
      label = sprintf("Subject %d: finite log_hazard", i)
    )
    expect_true(
      all(is.finite(batch_cppad[[i]]$biomarker)),
      label = sprintf("Subject %d: finite biomarker", i)
    )
  }

  # Test accuracy against true simulation values
  # Verify that ODE solutions match the true biomarker/velocity/acceleration
  # from the simulation. This confirms the ODE solver is correct.
  for (i in seq_len(length(data_processed))) {
    subject_id <- names(data_processed)[i]

    # Get true values from simulation at observation times
    idx <- sim$data$longitudinal_data$id == subject_id
    biomarker_data <- sim$data$longitudinal_data[idx, ]
    obs_times <- biomarker_data$time

    # Match ODE solution times with observation times
    result_times <- batch_cppad[[i]]$times
    match_idx <- match(obs_times, result_times)

    biomarker_true <- biomarker_data$biomarker
    biomarker_pred <- batch_cppad[[i]]$biomarker[match_idx]
    names(biomarker_pred) <- NULL

    velocity_true <- biomarker_data$velocity
    velocity_pred <- batch_cppad[[i]]$velocity[match_idx]
    names(velocity_pred) <- NULL

    acceleration_true <- biomarker_data$acceleration
    acceleration_pred <- batch_cppad[[i]]$acceleration[match_idx]
    names(acceleration_pred) <- NULL

    # Check accuracy
    expect_equal(
      biomarker_pred,
      biomarker_true,
      tolerance = 1e-5,
      label = sprintf("Subject %s: biomarker accuracy", subject_id)
    )
    expect_equal(
      velocity_pred,
      velocity_true,
      tolerance = 1e-5,
      label = sprintf("Subject %s: velocity accuracy", subject_id)
    )
    expect_equal(
      acceleration_pred,
      acceleration_true,
      tolerance = 1e-5,
      label = sprintf("Subject %s: acceleration accuracy", subject_id)
    )
  }
})
