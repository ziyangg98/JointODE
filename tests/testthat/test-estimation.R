# ==============================================================================
# Estimation Module Tests (E-step, M-step, gradient, parallel)
# ==============================================================================

# --- Shared data (computed once) ---

td <- .make_test_data(20)

posteriors <- JointODE:::.compute_posteriors(
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
    evals <- eigen(
      posteriors[[i]]$cov, symmetric = TRUE, only.values = TRUE
    )$values
    expect_true(all(evals > 0))
  }
})

test_that("Posterior modes are close to true random effects", {
  n <- length(td$data_list)
  posterior_modes <- t(vapply(
    posteriors, `[[`, numeric(ncol(td$random_effects)),
    "mode"
  ))

  bias <- colMeans(posterior_modes - td$random_effects)
  ese <- apply(posterior_modes, 2, sd)
  expect_true(all(3 * abs(bias) < ese))
})

# --- .compute_objective_expected: gradient + hessian ---

test_that("objective gradient matches finite differences", {
  re <- t(vapply(posteriors, `[[`, numeric(ncol(td$random_effects)), "mode"))

  obj_func <- function(theta) {
    params <- JointODE:::.vector_to_coef(td$parameters, theta)
    as.numeric(JointODE:::.compute_objective_expected(
      theta, td$data_list, re, params,
      gradient = FALSE, hessian = FALSE
    ))
  }

  result <- JointODE:::.compute_objective_expected(
    td$params, td$data_list, re, td$parameters,
    gradient = TRUE, hessian = TRUE
  )

  grad_ad <- as.vector(attr(result, "gradient"))
  grad_fd <- .finite_diff_gradient(obj_func, td$params)

  expect_lt(max(abs(grad_ad - grad_fd)), 0.01)
  expect_true(isSymmetric(attr(result, "hessian")))
})

test_that(".compute_joint_logpost gradient matches numDeriv", {
  data_i <- td$data_list[[1]]
  b0 <- posteriors[[1]]$mode

  result <- JointODE:::.compute_joint_logpost(
    random_effect = b0, data = data_i, parameters = td$parameters,
    gradient = TRUE, hessian = FALSE
  )
  grad_ad <- as.vector(attr(result, "gradient"))

  grad_num <- numDeriv::grad(
    func = function(b) {
      as.numeric(JointODE:::.compute_joint_logpost(
        random_effect = b, data = data_i, parameters = td$parameters,
        gradient = FALSE, hessian = FALSE
      ))
    },
    x = b0
  )

  expect_equal(grad_ad, grad_num, tolerance = 1e-4)
})

test_that(".compute_joint_logpost hessian is symmetric", {
  data_i <- td$data_list[[1]]
  b0 <- posteriors[[1]]$mode

  result <- JointODE:::.compute_joint_logpost(
    random_effect = b0, data = data_i, parameters = td$parameters,
    gradient = TRUE, hessian = TRUE
  )

  hess <- attr(result, "hessian")
  expect_true(is.matrix(hess))
  expect_equal(hess, t(hess), tolerance = 1e-10)
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

  m <- JointODE:::.compute_metrics(curr, prev, iter = 2)
  expect_equal(m$delta_l, 10)
  expect_true(m$rel_l > 0)

  m1 <- JointODE:::.compute_metrics(curr, NULL, iter = 1)
  expect_equal(m1$delta_l, -100)
})

# --- .update_random_effect_sigma ---

test_that(".update_random_effect_sigma computes Laplace-corrected estimate", {
  re <- matrix(c(1, 0, 0, 2), nrow = 2, byrow = TRUE)
  posteriors <- list(
    list(cov = diag(2) * 0.1),
    list(cov = diag(2) * 0.3)
  )
  expected <- (crossprod(re) + diag(2) * 0.1 + diag(2) * 0.3) / 2
  expect_equal(JointODE:::.update_random_effect_sigma(re, posteriors), expected)
})

# --- Parallel consistency ---

test_that(".compute_posteriors parallel matches sequential", {
  result_par <- JointODE:::.compute_posteriors(
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
  result_seq <- JointODE:::.compute_objective_expected(
    td$params, td$data_list, re, td$parameters,
    gradient = TRUE, hessian = TRUE, parallel = FALSE
  )
  result_par <- JointODE:::.compute_objective_expected(
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

test_that("non-PD random_effect_sigma fails fast", {
  td_bad <- .make_test_data(2)
  params_bad <- td_bad$parameters
  sigma_bad <- diag(nrow(params_bad$coefficients$random_effect_sigma))
  sigma_bad[1, 1] <- -1e-6
  params_bad$coefficients$random_effect_sigma <- sigma_bad
  theta_bad <- JointODE:::.coef_to_vector(params_bad)

  expect_error(
    JointODE:::.compute_objective_expected(
      theta_bad, td_bad$data_list, td_bad$random_effects, params_bad,
      gradient = FALSE, hessian = FALSE
    ),
    "random_effect_sigma must be positive definite"
  )

  expect_error(
    JointODE:::.compute_joint_logpost(
      random_effect = td_bad$random_effects[1, ],
      data = td_bad$data_list[[1]],
      parameters = params_bad,
      gradient = TRUE,
      hessian = TRUE
    ),
    "random_effect_sigma must be positive definite"
  )
})
