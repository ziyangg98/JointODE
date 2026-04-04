# ==============================================================================
# End-to-End Integration Tests: JointODE + MarginalODE
# ==============================================================================

# ==============================================================================
# JointODE
# ==============================================================================

e2e_sim <- JointODE::simulate(
  n_subjects = 20,
  longitudinal = list(
    xi = c(mean = 0.4, sd = 0.05),
    period = c(mean = 6, sd = 0.1),
    excitation = list(offset = 0, covariates = c(x1 = 0.5)),
    initial = list(
      biomarker = c(mean = -0.5, sd = 0.1),
      velocity = c(mean = -0.1, sd = 0.1)
    ),
    error_sd = 0.1, n_measurements = 15
  ),
  survival = list(
    baseline = list(type = "weibull", shape = 2.0, scale = 15.0),
    value = 0.8, slope = 2.0, gamma = 1, covariates = c(w1 = 0.6)
  ),
  covariates = list(
    x1 = list(type = "normal", mean = 0, sd = 1),
    w1 = list(type = "normal", mean = 0, sd = 1)
  ),
  seed = 42
)

e2e_sim$longitudinal_data$biomarker <- NULL
e2e_sim$longitudinal_data$velocity <- NULL
e2e_sim$longitudinal_data$acceleration <- NULL

e2e_fit <- JointODE(
  longitudinal_formula = observed ~ biomarker + velocity + x1 +
    (biomarker + velocity | id),
  survival_formula = Surv(time, status) ~ w1,
  longitudinal_data = e2e_sim$longitudinal_data,
  survival_data = e2e_sim$survival_data,
  init = "marginal",
  control = list(maxit = 10, tol = 1e-2)
)

test_that("JointODE fits and converges", {
  expect_s3_class(e2e_fit, "JointODE")
  expect_true(is.finite(e2e_fit$logLik))
  expect_true(is.matrix(e2e_fit$vcov))
})

test_that("JointODE exposes EM history in convergence info", {
  conv <- e2e_fit$convergence

  expect_true(is.numeric(conv$loglik_history))
  expect_true(is.numeric(conv$delta_theta_history))
  expect_true(is.numeric(conv$delta_loglik_history))

  expect_equal(length(conv$loglik_history), conv$iterations)
  expect_equal(length(conv$delta_theta_history), conv$iterations)
  expect_equal(length(conv$delta_loglik_history), conv$iterations)

  expect_true(all(is.finite(conv$loglik_history)))
  expect_true(all(is.finite(conv$delta_theta_history)))
})

test_that("JointODE S3 methods work on fitted object", {
  expect_true(is.numeric(coef(e2e_fit)))
  expect_s3_class(logLik(e2e_fit), "logLik")
  expect_s3_class(summary(e2e_fit), "summary.JointODE")
  expect_true(length(capture.output(print(e2e_fit))) > 0)
  expect_s3_class(predict(e2e_fit), "data.frame")
})

# ==============================================================================
# MarginalODE
# ==============================================================================

marginal_sim <- JointODE::simulate(
  n_subjects = 10,
  longitudinal = list(
    xi = c(mean = 0.5, sd = 0.05),
    period = c(mean = 6, sd = 0.1),
    excitation = list(offset = 0, covariates = c(x1 = 0.5)),
    initial = list(
      biomarker = c(mean = -0.5, sd = 0.1),
      velocity = c(mean = -0.1, sd = 0.1)
    ),
    error_sd = 0.1, n_measurements = 15
  ),
  covariates = list(
    x1 = list(type = "normal", mean = 0, sd = 1),
    w1 = list(type = "normal", mean = 0, sd = 1),
    w2 = list(type = "binary", prob = 0.5)
  ),
  seed = 42
)

# --- Fitting ---

test_that("MarginalODE fits and estimates initial states", {
  fit <- MarginalODE(
    formula = observed ~ x1,
    data = marginal_sim$longitudinal_data,
    control = list(maxit = 20, tol = 1e-3)
  )
  expect_s3_class(fit, "MarginalODE")
  expect_equal(length(fit$parameters), 4)
  expect_true(fit$measurement_error_sd > 0)
  expect_equal(nrow(fit$initial_states), 10)
  expect_equal(ncol(fit$initial_states), 2)
})

# --- S3 methods ---

marginal_fit <- MarginalODE(
  formula = observed ~ x1,
  data = marginal_sim$longitudinal_data,
  control = list(maxit = 20, tol = 1e-3)
)

test_that("MarginalODE S3 methods work", {
  # print
  out <- capture.output(result <- print(marginal_fit))
  expect_identical(result, marginal_fit)
  expect_true(any(grepl("Marginal", out)))

  # summary
  s <- summary(marginal_fit)
  expect_s3_class(s, "summary.MarginalODE")
  expect_equal(ncol(s$coefficients), 4)
  out2 <- capture.output(print(s))
  expect_true(any(grepl("Coefficients", out2)))

  expect_equal(coef(marginal_fit), marginal_fit$parameters)
  expect_equal(vcov(marginal_fit), marginal_fit$vcov)
  ll <- logLik(marginal_fit)
  expect_s3_class(ll, "logLik")
  expect_true(attr(ll, "df") > 0)

  # predict
  pred <- predict(marginal_fit)
  expect_s3_class(pred, "data.frame")
  expect_true(all(c("id", "time", "biomarker", "velocity") %in% names(pred)))
  expect_error(predict(marginal_fit, newdata = data.frame()), "not yet")
})

test_that("MarginalODE parameters are named", {
  expect_true("value" %in% names(marginal_fit$parameters))
  expect_true("slope" %in% names(marginal_fit$parameters))
})

# --- C++ gradient verification ---

test_that(".compute_marginal_objective gradient is correct", {
  data_list <- JointODE:::.process_marginal(
    observed ~ x1, marginal_sim$longitudinal_data,
    "time", "id"
  )
  theta <- c(-0.5, -0.3, 0.2)

  result <- .compute_marginal_objective(
    theta, data_list,
    gradient = TRUE, hessian = FALSE
  )
  analytic <- attr(result, "gradient")
  numeric <- numDeriv::grad(function(x) {
    as.numeric(.compute_marginal_objective(
      x, data_list,
      gradient = FALSE, hessian = FALSE
    ))
  }, theta)
  expect_equal(analytic, numeric, tolerance = 1e-4)
})

test_that(".solve_batch_marginal returns correct structure", {
  data_list <- JointODE:::.process_marginal(
    observed ~ x1, marginal_sim$longitudinal_data,
    "time", "id"
  )
  sols <- .solve_batch_marginal(data_list, c(-0.5, -0.3, 0.2))
  expect_equal(length(sols), length(data_list))
  expect_named(sols[[1]], c("times", "biomarker", "velocity", "acceleration"))
  expect_true(all(is.finite(sols[[1]]$biomarker)))
})

test_that(".compute_marginal_state gradient is correct", {
  data_list <- JointODE:::.process_marginal(
    observed ~ x1, marginal_sim$longitudinal_data,
    "time", "id"
  )
  theta <- c(-0.5, -0.3, 0.2)
  state0 <- c(0.5, 0.1)

  result <- .compute_marginal_state(
    state0, data_list[[1]], theta,
    gradient = TRUE, hessian = FALSE
  )
  analytic <- attr(result, "gradient")
  numeric <- numDeriv::grad(function(x) {
    as.numeric(.compute_marginal_state(
      x, data_list[[1]], theta,
      gradient = FALSE, hessian = FALSE
    ))
  }, state0)
  expect_equal(analytic, numeric, tolerance = 1e-4)
})
