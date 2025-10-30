test_that("AD gradient matches forward sensitivity gradient", {
  # Setup: Use full dataset
  data_list <- .process(
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
  random_effects <- sim$data$random_effects

  # Compute posteriors and measure time
  cat("\n=== Computing Posteriors ===\n")
  time_posterior <- system.time({
    posteriors <- .compute_posteriors(
      data_list = data_list,
      random_effects = random_effects,
      parameters = sim$init,
      parallel = TRUE,
      n_cores = 0,
      level = 3
    )
  })
  cat("Posterior computation time:", time_posterior["elapsed"], "seconds\n")

  coefficients <- sim$init$coefficients
  params <- c(
    coefficients$baseline,
    coefficients$hazard,
    as.vector(coefficients$longitudinal)
  )
  # Compute using AD and measure time
  cat("\n=== Computing Gradient & Hessian ===\n")
  time_ad <- system.time({
    result_ad <- .compute_objective_expected(
      params = params,
      data_list = data_list,
      posteriors = posteriors,
      parameters = sim$init,
      gradient = TRUE,
      hessian = TRUE
    )
  })
  cat("CppAD gradient time:", time_ad["elapsed"], "seconds\n")
})
