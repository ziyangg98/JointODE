# ==============================================================================
# Hazard Integration Stability Tests
# ==============================================================================

td_hazard <- .make_test_data(5)

test_that("log_hazard output is clamped under extreme coefficients", {
  params <- td_hazard$parameters
  params$coefficients$hazard[] <- 1e6
  params$coefficients$baseline[] <- 50

  re_zero <- matrix(0,
    nrow = nrow(td_hazard$random_effects),
    ncol = ncol(td_hazard$random_effects)
  )

  sol <- .solve_batch_joint(td_hazard$data_list, re_zero, params)

  for (i in seq_along(sol)) {
    expect_true(all(sol[[i]]$log_hazard <= 20 + 1e-12),
      info = sprintf("subject %d upper clamp", i)
    )
    expect_true(all(sol[[i]]$log_hazard >= -20 - 1e-12),
      info = sprintf("subject %d lower clamp", i)
    )
    expect_true(all(is.finite(sol[[i]]$cum_hazard)),
      info = sprintf("subject %d cum_hazard finite", i)
    )
  }
})

test_that("hazard quadrature remains finite under high curvature", {
  params_k1 <- td_hazard$parameters
  params_k1$coefficients$hazard[] <- c(6, -6, 4, -4)
  params_k1$configurations$hazard_quadrature <- 1L

  params_k20 <- params_k1
  params_k20$configurations$hazard_quadrature <- 20L

  sol_k1 <- .solve_batch_joint(
    td_hazard$data_list, td_hazard$random_effects, params_k1
  )
  sol_k20 <- .solve_batch_joint(
    td_hazard$data_list, td_hazard$random_effects, params_k20
  )

  for (i in seq_along(sol_k1)) {
    expect_true(all(is.finite(sol_k1[[i]]$cum_hazard)),
      info = sprintf("subject %d k=1 finite", i)
    )
    expect_true(all(is.finite(sol_k20[[i]]$cum_hazard)),
      info = sprintf("subject %d k=20 finite", i)
    )
    expect_true(min(diff(sol_k1[[i]]$cum_hazard)) >= -1e-10,
      info = sprintf("subject %d k=1 monotone", i)
    )
    expect_true(min(diff(sol_k20[[i]]$cum_hazard)) >= -1e-10,
      info = sprintf("subject %d k=20 monotone", i)
    )
  }
})
