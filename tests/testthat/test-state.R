test_that("State optimization with single random effect", {
  data("sim", package = "JointODE")

  data_list <- .process(
    longitudinal_formula = observed ~ biomarker +
      velocity +
      x1 +
      x2 +
      (biomarker + velocity | id),
    survival_formula = Surv(time, status) ~ w1 + w2,
    longitudinal_data = sim$data$longitudinal_data,
    survival_data = sim$data$survival_data,
    state = as.matrix(sim$data$state)
  )

  n_subjects <- length(data_list)
  results <- matrix(NA, n_subjects, 2)

  for (i in 1:n_subjects) {
    subj_data <- data_list[[i]]
    params <- sim$init
    random_effect <- sim$data$random_effects[i, ]

    opt_result <- .estimate_state(
      initial_guess = c(0, 0),
      data = subj_data,
      random_effect = random_effect,
      parameters = params
    )

    # Single Newton step should decrease objective
    expect_true(opt_result$obj_change <= 1e-6,
      info = sprintf("subject %d: obj increased by %.4e", i,
                     opt_result$obj_change))

    results[i, ] <- opt_result$state
  }

  errors <- results - sim$data$state
  bias <- colMeans(errors)
  rmse <- sqrt(colMeans(errors^2))

  expect_true(all(abs(bias) < 0.5),
    info = paste("bias too large:", paste(round(bias, 4), collapse = ", ")))
  expect_true(all(rmse < 1.0),
    info = paste("rmse too large:", paste(round(rmse, 4), collapse = ", ")))
})
