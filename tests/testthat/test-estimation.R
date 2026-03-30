# ==============================================================================
# Estimation Module Tests (E-step, M-step, gradient, parallel)
# ==============================================================================

# --- Shared data (computed once) ---

td <- .make_test_data(10)

posteriors <- .compute_posteriors(
  td$data_list, td$parameters, td$random_effects,
  parallel = FALSE, level = 3
)

# --- .expand_posteriors ---

test_that(".expand_posteriors expands correctly", {
  post_2 <- list(
    nodes = list(
      matrix(rnorm(6), nrow = 3, ncol = 2),
      matrix(rnorm(4), nrow = 2, ncol = 2)
    ),
    weights = list(c(0.3, 0.5, 0.2), c(0.6, 0.4))
  )
  data_2 <- list(list(id = 1, v = "a"), list(id = 2, v = "b"))
  expanded <- .expand_posteriors(data_2, post_2)

  expect_equal(nrow(expanded$nodes), 5)
  expect_equal(length(expanded$weights), 5)
  expect_equal(length(expanded$data), 5)
  expect_equal(expanded$node_to_subject, c(1, 1, 1, 2, 2))
  expect_equal(expanded$data[[1]]$v, "a")
  expect_equal(expanded$data[[4]]$v, "b")
})

test_that(".expand_posteriors preserves weights", {
  post_2 <- list(
    nodes = list(matrix(0, 2, 1), matrix(0, 3, 1)),
    weights = list(c(0.4, 0.6), c(0.2, 0.3, 0.5))
  )
  expanded <- .expand_posteriors(
    list(list(id = 1), list(id = 2)), post_2
  )
  expect_equal(expanded$weights, c(0.4, 0.6, 0.2, 0.3, 0.5))
})

# --- .compute_posteriors: structure + unbiasedness ---

test_that("Posterior AGHQ structure is valid", {
  n <- length(td$data_list)
  expect_equal(length(posteriors$nodes), n)
  expect_equal(length(posteriors$weights), n)

  for (i in seq_len(n)) {
    nodes_i <- posteriors$nodes[[i]]
    weights_i <- posteriors$weights[[i]]

    expect_true(is.matrix(nodes_i))
    expect_equal(ncol(nodes_i), ncol(td$random_effects))
    expect_equal(length(weights_i), nrow(nodes_i))
    expect_true(all(weights_i > 0))
    expect_equal(sum(weights_i), 1, tolerance = 1e-14)
  }
})

test_that("Posterior means are unbiased (AGHQ)", {
  n <- length(td$data_list)
  posterior_means <- t(sapply(seq_len(n), function(i) {
    colSums(sweep(posteriors$nodes[[i]], 1, posteriors$weights[[i]], "*"))
  }))

  bias <- colMeans(posterior_means - td$random_effects)
  ese <- apply(posterior_means, 2, sd)
  expect_true(all(3 * abs(bias) < ese))
})

# --- .compute_objective_expected: gradient + hessian ---

test_that("AD gradient and hessian are valid", {
  result <- .compute_objective_expected(
    td$params, td$data_list, posteriors, td$parameters,
    gradient = TRUE, hessian = TRUE
  )

  grad <- attr(result, "gradient")
  expect_true(all(is.finite(grad)), info = "gradient non-finite")
  expect_true(any(grad != 0), info = "gradient all zeros")

  hess <- attr(result, "hessian")
  expect_true(all(is.finite(hess)), info = "hessian non-finite")
  expect_equal(hess, t(hess), tolerance = 1e-10, info = "hessian asymmetric")
})

# --- .compute_metrics ---

test_that(".compute_metrics computes correctly", {
  curr <- list(
    loglik = -100,
    parameters = list(coefficients = list(
      baseline = c(0.1, 0.2), hazard = 0.5,
      longitudinal = -0.01, initial_state = c(0, 0),
      measurement_error_sd = 0.5, random_effect_sigma = diag(4)
    ))
  )
  prev <- curr
  prev$loglik <- -110

  m <- .compute_metrics(curr, prev, iter = 2)
  expect_equal(m$delta_l, 10)
  expect_true(m$rel_l > 0)

  m1 <- .compute_metrics(curr, NULL, iter = 1)
  expect_equal(m1$delta_l, -100)
})

# --- .update_random_effect_sigma ---

test_that(".update_random_effect_sigma computes mean of second moments", {
  moments <- list(
    list(mean = c(0, 0), second_moment = diag(2)),
    list(mean = c(0, 0), second_moment = diag(2) * 3)
  )
  expect_equal(.update_random_effect_sigma(moments, 2), diag(2) * 2)
})

# --- Parallel consistency (reuse td and posteriors) ---

test_that(".compute_posteriors parallel matches sequential", {
  result_par <- .compute_posteriors(
    td$data_list, td$parameters, td$random_effects,
    parallel = TRUE, n_cores = 2, level = 3
  )

  for (i in seq_along(td$data_list)) {
    expect_equal(result_par$nodes[[i]], posteriors$nodes[[i]],
      tolerance = 1e-10, info = sprintf("subject %d nodes", i)
    )
    expect_equal(result_par$weights[[i]], posteriors$weights[[i]],
      tolerance = 1e-10, info = sprintf("subject %d weights", i)
    )
  }
})

test_that(".compute_objective_expected parallel matches sequential", {
  result_seq <- .compute_objective_expected(
    td$params, td$data_list, posteriors, td$parameters,
    gradient = TRUE, hessian = TRUE, parallel = FALSE
  )
  result_par <- .compute_objective_expected(
    td$params, td$data_list, posteriors, td$parameters,
    gradient = TRUE, hessian = TRUE, parallel = TRUE, n_cores = 2
  )

  expect_equal(as.numeric(result_par), as.numeric(result_seq),
    tolerance = 1e-10, info = "objective value"
  )
  expect_equal(attr(result_par, "gradient"), attr(result_seq, "gradient"),
    tolerance = 1e-10, info = "gradient"
  )
  expect_equal(attr(result_par, "hessian"), attr(result_seq, "hessian"),
    tolerance = 1e-10, info = "hessian"
  )
})
