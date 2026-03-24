test_that("State optimization with single random effect", {
  data("sim", package = "JointODE")

  n_test <- 50
  test_ids <- unique(sim$data$longitudinal_data$id)[seq_len(n_test)]
  data_list <- .process_joint(
    longitudinal_formula = observed ~ biomarker +
      velocity + x1 + x2 + (biomarker + velocity | id),
    survival_formula = Surv(time, status) ~ w1 + w2,
    longitudinal_data = sim$data$longitudinal_data[
      sim$data$longitudinal_data$id %in% test_ids,
    ],
    survival_data = sim$data$survival_data[
      sim$data$survival_data$id %in% test_ids,
    ]
  )

  n_subjects <- length(data_list)
  results <- matrix(NA, n_subjects, 2)

  for (i in 1:n_subjects) {
    subj_data <- data_list[[i]]
    params <- sim$init
    # Pass only covariate RE (cols 3-4), skip state RE (cols 1-2)
    random_effect <- sim$data$random_effects[i, 3:4]

    results[i, ] <- .estimate_joint_state(
      initial_guess = c(0, 0),
      data = subj_data,
      random_effect = random_effect,
      parameters = params
    )
  }

  errors <- results - sim$data$state[seq_len(n_test), ]
  bias <- colMeans(errors)
  rmse <- sqrt(colMeans(errors^2))

  expect_true(all(abs(bias) < 0.05),
    info = paste("bias too large:", paste(round(bias, 4), collapse = ", "))
  )
  expect_true(all(rmse < 0.1),
    info = paste("rmse too large:", paste(round(rmse, 4), collapse = ", "))
  )
})
