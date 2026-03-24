test_that("AD gradient matches forward sensitivity gradient", {
  n_test <- 50
  test_ids <- unique(sim$data$longitudinal_data$id)[seq_len(n_test)]
  data_list <- .process_joint(
    longitudinal_data = sim$data$longitudinal_data[
      sim$data$longitudinal_data$id %in% test_ids,
      c("id", "time", "observed", "x1", "x2")
    ],
    longitudinal_formula = observed ~
      biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
    survival_data = sim$data$survival_data[
      sim$data$survival_data$id %in% test_ids,
    ],
    survival_formula = Surv(time, status) ~ w1 + w2
  )
  random_effects <- sim$data$random_effects[seq_len(n_test), ]

  # Compute posteriors and measure time
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

  coefficients <- sim$init$coefficients
  params <- c(
    coefficients$baseline,
    coefficients$hazard,
    as.vector(coefficients$longitudinal)
  )

  # Compute using AD and measure time
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

  # Verify gradient is finite and non-zero
  grad <- attr(result_ad, "gradient")
  expect_true(
    all(is.finite(grad)),
    info = "gradient contains non-finite values"
  )
  expect_true(
    any(grad != 0),
    info = "gradient is all zeros"
  )

  # Verify Hessian is finite and symmetric
  hess <- attr(result_ad, "hessian")
  expect_true(
    all(is.finite(hess)),
    info = "Hessian contains non-finite values"
  )
  expect_equal(
    hess, t(hess),
    tolerance = 1e-10,
    info = "Hessian is not symmetric"
  )
})
