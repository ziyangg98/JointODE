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

test_that("hazard_quadrature Simpson quadrature convergence", {
  params_k1 <- td$parameters
  params_k1$configurations$hazard_quadrature <- 1L
  params_k5 <- td$parameters
  params_k5$configurations$hazard_quadrature <- 5L
  params_k50 <- td$parameters
  params_k50$configurations$hazard_quadrature <- 50L

  sol_k1 <- .solve_batch_joint(td$data_list, td$random_effects, params_k1)
  sol_k5 <- .solve_batch_joint(td$data_list, td$random_effects, params_k5)
  sol_k50 <- .solve_batch_joint(td$data_list, td$random_effects, params_k50)

  for (i in seq_along(td$data_list)) {
    # Same time grid
    expect_equal(sol_k5[[i]]$times, sol_k1[[i]]$times)

    # Biomarker/velocity identical (both use exact matexp)
    expect_equal(sol_k5[[i]]$biomarker, sol_k1[[i]]$biomarker, tolerance = 1e-10)
    expect_equal(sol_k5[[i]]$velocity, sol_k1[[i]]$velocity, tolerance = 1e-10)

    # Cumulative hazard converges: k=5 closer to k=50 than k=1
    H1 <- tail(sol_k1[[i]]$cum_hazard, 1)
    H5 <- tail(sol_k5[[i]]$cum_hazard, 1)
    H50 <- tail(sol_k50[[i]]$cum_hazard, 1)
    expect_lt(abs(H5 - H50), abs(H1 - H50) + 1e-15)
  }
})
