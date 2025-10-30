test_that("Posterior means are unbiased (AGHQ)", {
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

  posteriors <- .compute_posteriors(
    data_list = data_list,
    random_effects = random_effects,
    parameters = sim$init,
    parallel = TRUE,
    n_cores = 0,
    level = 3
  )

  # Test AGHQ structure: should return nodes and weights
  expect_true(is.list(posteriors$nodes))
  expect_true(is.list(posteriors$weights))
  n_subjects <- nrow(random_effects)
  expect_equal(length(posteriors$nodes), n_subjects)
  expect_equal(length(posteriors$weights), n_subjects)

  # Compute posterior means from AGHQ nodes and weights
  posterior_means <- t(sapply(seq_len(n_subjects), function(i) {
    colSums(sweep(posteriors$nodes[[i]], 1, posteriors$weights[[i]], "*"))
  }))

  # Test bias
  true_b <- random_effects
  bias <- colMeans(posterior_means - true_b)
  ese <- apply(posterior_means, 2, sd)
  expect_true(all(3 * abs(bias) < ese))

  # Test AGHQ nodes structure
  n_random <- ncol(sim$data$random_effects)
  for (i in seq_len(n_subjects)) {
    nodes_i <- posteriors$nodes[[i]]
    weights_i <- posteriors$weights[[i]]

    # Check nodes dimensions
    expect_true(is.matrix(nodes_i))
    expect_equal(ncol(nodes_i), n_random)

    # Check weights properties
    expect_true(is.numeric(weights_i))
    expect_equal(length(weights_i), nrow(nodes_i))
    expect_true(all(weights_i > 0))
    expect_equal(sum(weights_i), 1, tolerance = 1e-10) # Normalized
  }
})
