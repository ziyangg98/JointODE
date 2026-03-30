# ODE Solver Tests

td <- .make_test_data(10)

test_that(".solve_batch_joint structure and accuracy", {
  batch_sol <- .solve_batch_joint(
    td$data_list, td$random_effects, td$parameters
  )

  n <- length(td$data_list)
  expect_equal(length(batch_sol), n)

  # Structure
  expected_fields <- c(
    "times", "cum_hazard", "biomarker",
    "velocity", "acceleration", "log_hazard"
  )
  for (field in expected_fields) {
    expect_true(field %in% names(batch_sol[[1]]))
  }
  expect_equal(
    length(batch_sol[[1]]$times),
    length(batch_sol[[1]]$biomarker)
  )
  expect_false("dcumhazard_deta_at_event" %in% names(batch_sol[[1]]))

  # Finite values
  for (i in seq_len(n)) {
    expect_true(all(is.finite(batch_sol[[i]]$cum_hazard)))
    expect_true(all(is.finite(batch_sol[[i]]$biomarker)))
  }

  # Accuracy vs true simulation values
  for (i in seq_along(td$data_list)) {
    subject_id <- names(td$data_list)[i]
    idx <- sim$data$longitudinal_data$id == subject_id
    true_data <- sim$data$longitudinal_data[idx, ]
    match_idx <- match(true_data$time, batch_sol[[i]]$times)

    expect_equal(
      unname(batch_sol[[i]]$biomarker[match_idx]),
      true_data$biomarker, tolerance = 1e-5,
      info = sprintf("subject %s biomarker", subject_id)
    )
    expect_equal(
      unname(batch_sol[[i]]$velocity[match_idx]),
      true_data$velocity, tolerance = 1e-5,
      info = sprintf("subject %s velocity", subject_id)
    )
    expect_equal(
      unname(batch_sol[[i]]$acceleration[match_idx]),
      true_data$acceleration, tolerance = 1e-5,
      info = sprintf("subject %s acceleration", subject_id)
    )
  }
})
