# ==============================================================================
# EM Diagnostic Logic Tests
# ==============================================================================

td_diag <- .make_test_data(3)

.make_control <- function(tol = 1e-3, maxit = 10, verbose = 0) {
  list(tol = tol, maxit = maxit, verbose = verbose)
}

test_that(".track converges when delta_theta below tolerance", {
  prev <- list(parameters = td_diag$parameters, loglik = -100)
  curr <- prev

  curr$parameters$coefficients$baseline[1] <-
    curr$parameters$coefficients$baseline[1] + 1e-8
  curr$loglik <- -99.5

  status <- .track(iter = 2, curr = curr, prev = prev,
    control = .make_control(tol = 1e-6)
  )

  expect_true(status$converged)
})

test_that(".track does not converge at first iteration", {
  prev <- list(parameters = td_diag$parameters, loglik = -Inf)
  curr <- list(parameters = td_diag$parameters, loglik = -100)

  status <- .track(iter = 1, curr = curr, prev = prev,
    control = .make_control(tol = 1e-6)
  )

  expect_false(status$converged)
})

test_that(".track returns non-converged when log-likelihood is NA", {
  prev <- list(parameters = td_diag$parameters, loglik = -100)
  curr <- list(parameters = td_diag$parameters, loglik = NA_real_)

  status <- .track(iter = 2, curr = curr, prev = prev,
    control = .make_control(tol = 1e-6)
  )

  expect_false(status$converged)
})
