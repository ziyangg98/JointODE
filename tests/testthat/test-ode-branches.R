# ==============================================================================
# ODE Branch Coverage Tests
# ==============================================================================

td_branch <- .make_test_data(3)

.branch_cases <- list(
  REAL = c(b1 = 0.5, b2 = 0.0),
  COMPLEX = c(b1 = -1.0, b2 = 0.0),
  REPEATED = c(b1 = -0.25, b2 = 1.0),
  FIRST_ORD = c(b1 = 0.0, b2 = 0.7),
  ZERO = c(b1 = 0.0, b2 = 0.0)
)

test_that(".solve_batch_joint is stable across ODE branch regimes", {
  for (case_name in names(.branch_cases)) {
    case <- .branch_cases[[case_name]]
    params <- td_branch$parameters

    # longitudinal[1:2] map to dyn_biomarker (b1) and dyn_velocity (b2)
    params$coefficients$longitudinal[1] <- unname(case["b1"])
    params$coefficients$longitudinal[2] <- unname(case["b2"])

    # Keep branch classification determined by fixed effects only.
    re_zero <- matrix(0, nrow = nrow(td_branch$random_effects),
      ncol = ncol(td_branch$random_effects)
    )

    sol <- JointODE:::.solve_batch_joint(td_branch$data_list, re_zero, params)

    expect_equal(length(sol), length(td_branch$data_list),
      info = sprintf("%s: number of subjects", case_name)
    )

    for (i in seq_along(sol)) {
      for (field in c("biomarker", "velocity", "cum_hazard")) {
        .expect_all_finite(sol[[i]][[field]])
      }

      expect_equal(sol[[i]]$times, sort(sol[[i]]$times),
        info = sprintf("%s: subject %d time monotonic", case_name, i)
      )
      expect_true(min(diff(sol[[i]]$cum_hazard)) >= -1e-10,
        info = sprintf("%s: subject %d cumulative hazard monotone", case_name, i)
      )
    }
  }
})
