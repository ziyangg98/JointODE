# ==============================================================================
# Estimation Module Tests (E-step, M-step, gradient, parallel)
# ==============================================================================

# --- Shared data (computed once) ---

td <- .make_test_data(10)

posteriors <- .compute_posteriors(
  td$data_list, td$parameters, td$random_effects,
  parallel = FALSE
)

# --- .compute_posteriors: Laplace structure ---

test_that("Posterior Laplace structure is valid", {
  n <- length(td$data_list)
  expect_equal(length(posteriors), n)

  for (i in seq_len(n)) {
    expect_true(is.numeric(posteriors[[i]]$mode))
    expect_equal(length(posteriors[[i]]$mode), ncol(td$random_effects))
    expect_true(is.matrix(posteriors[[i]]$cov))
    expect_equal(nrow(posteriors[[i]]$cov), ncol(td$random_effects))
    # Covariance must be positive definite
    eigenvalues <- eigen(posteriors[[i]]$cov, symmetric = TRUE,
      only.values = TRUE
    )$values
    expect_true(all(eigenvalues > 0))
  }
})

test_that("Posterior modes are close to true random effects", {
  n <- length(td$data_list)
  posterior_modes <- t(vapply(posteriors, `[[`, numeric(ncol(td$random_effects)),
    "mode"
  ))

  bias <- colMeans(posterior_modes - td$random_effects)
  ese <- apply(posterior_modes, 2, sd)
  expect_true(all(3 * abs(bias) < ese))
})

# --- .compute_objective_expected: gradient + hessian ---

test_that("AD gradient and hessian are valid", {
  re <- t(vapply(posteriors, `[[`, numeric(ncol(td$random_effects)), "mode"))
  result <- .compute_objective_expected(
    td$params, td$data_list, re, td$parameters,
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

# --- Parallel consistency ---

test_that(".compute_posteriors parallel matches sequential", {
  result_par <- .compute_posteriors(
    td$data_list, td$parameters, td$random_effects,
    parallel = TRUE, n_cores = 2
  )

  for (i in seq_along(td$data_list)) {
    expect_equal(result_par[[i]]$mode, posteriors[[i]]$mode,
      tolerance = 1e-10, info = sprintf("subject %d mode", i)
    )
    expect_equal(result_par[[i]]$cov, posteriors[[i]]$cov,
      tolerance = 1e-10, info = sprintf("subject %d cov", i)
    )
  }
})

test_that(".compute_objective_expected parallel matches sequential", {
  re <- t(vapply(posteriors, `[[`, numeric(ncol(td$random_effects)), "mode"))
  result_seq <- .compute_objective_expected(
    td$params, td$data_list, re, td$parameters,
    gradient = TRUE, hessian = TRUE, parallel = FALSE
  )
  result_par <- .compute_objective_expected(
    td$params, td$data_list, re, td$parameters,
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
