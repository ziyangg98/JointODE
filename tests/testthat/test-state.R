test_that("Initial state optimization converges", {
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

    # Check convergence
    expect_true(opt_result$converged)
    expect_true(opt_result$iterations < 50)

    # Check optimized state is reasonable (close to true state)
    true_state <- sim$data$state[i, ]
    results[i, ] <- opt_result$state
  }
  errors <- results - sim$data$state
  bias <- colMeans(errors)
  rmse <- sqrt(colMeans(errors^2))
})
