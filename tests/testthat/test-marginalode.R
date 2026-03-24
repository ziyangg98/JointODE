# ==============================================================================
# MarginalODE Unit Tests
# ==============================================================================

# Generate small test data for MarginalODE
marginal_sim <- JointODE::simulate(
  n_subjects = 20,
  longitudinal = list(
    xi = c(mean = 0.5, sd = 0.05),
    period = c(mean = 6, sd = 0.1),
    excitation = list(offset = 0, covariates = c(x1 = 0.5)),
    initial = list(
      offset = c(biomarker = -0.5, velocity = -0.1),
      covariates = list(
        biomarker = c(x1 = 0.08),
        velocity = c(x1 = 0.1)
      )
    ),
    error_sd = 0.1,
    n_measurements = 15
  ),
  covariates = list(
    x1 = list(type = "normal", mean = 0, sd = 1),
    w1 = list(type = "normal", mean = 0, sd = 1),
    w2 = list(type = "binary", prob = 0.5)
  ),
  seed = 42
)

# --- MarginalODE fitting ---

test_that("MarginalODE fits with provided state", {
  fit <- MarginalODE(
    formula = observed ~ x1,
    data = marginal_sim$longitudinal_data,
    state = marginal_sim$state
  )

  expect_s3_class(fit, "MarginalODE")
  expect_true(fit$convergence$converged)
  # value, slope, (Intercept), x1
  expect_equal(length(fit$parameters), 4)
  expect_true(is.numeric(fit$parameters))
  expect_true(is.finite(fit$logLik))
  expect_true(fit$measurement_error_sd > 0)
})

test_that("MarginalODE fits without state (iterative)", {
  fit <- MarginalODE(
    formula = observed ~ x1,
    data = marginal_sim$longitudinal_data,
    control = list(maxit = 20, tol = 1e-3)
  )

  expect_s3_class(fit, "MarginalODE")
  expect_true(is.numeric(fit$parameters))
  expect_true(!is.null(fit$initial_states))
  expect_equal(nrow(fit$initial_states), 20)
  expect_equal(ncol(fit$initial_states), 2)
})

# --- S3 methods ---

marginal_fit <- MarginalODE(
  formula = observed ~ x1,
  data = marginal_sim$longitudinal_data,
  state = marginal_sim$state
)

test_that("print.MarginalODE works", {
  out <- capture.output(result <- print(marginal_fit))
  expect_identical(result, marginal_fit)
  expect_true(any(grepl("Marginal", out)))
  expect_true(any(grepl("Log-likelihood", out)))
})

test_that("summary.MarginalODE returns correct structure", {
  s <- summary(marginal_fit)
  expect_s3_class(s, "summary.MarginalODE")

  expect_true(!is.null(s$coefficients))
  expect_equal(ncol(s$coefficients), 4)
  expect_equal(nrow(s$coefficients), length(marginal_fit$parameters))
  expect_true(s$nobs > 0)
  expect_true(s$n_observations > 0)
})

test_that("print.summary.MarginalODE works", {
  s <- summary(marginal_fit)
  out <- capture.output(result <- print(s))
  expect_identical(result, s)
  expect_true(any(grepl("Coefficients", out)))
})

test_that("MarginalODE parameters are named", {
  expect_true(is.numeric(marginal_fit$parameters))
  expect_true(!is.null(names(marginal_fit$parameters)))
  expect_true("value" %in% names(marginal_fit$parameters))
  expect_true("slope" %in% names(marginal_fit$parameters))
})

test_that("MarginalODE vcov is accessible", {
  expect_true(is.matrix(marginal_fit$vcov))
  n <- length(marginal_fit$parameters)
  expect_equal(dim(marginal_fit$vcov), c(n, n))
})

test_that("predict.MarginalODE returns data.frame", {
  pred <- predict(marginal_fit)
  expect_s3_class(pred, "data.frame")
  expect_true(
    all(c("id", "time", "biomarker", "velocity", "acceleration") %in%
      names(pred))
  )
  expect_true(nrow(pred) > 0)
})

test_that("predict.MarginalODE errors on newdata", {
  expect_error(
    predict(marginal_fit, newdata = data.frame()),
    "not yet supported"
  )
})

test_that("coef.MarginalODE works", {
  expect_equal(coef(marginal_fit), marginal_fit$parameters)
})

test_that("vcov.MarginalODE works", {
  expect_equal(vcov(marginal_fit), marginal_fit$vcov)
})

test_that("logLik.MarginalODE works", {
  ll <- logLik(marginal_fit)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), marginal_fit$logLik)
  expect_true(attr(ll, "df") > 0)
  expect_true(attr(ll, "nobs") > 0)
})

# --- C++ functions ---

test_that(".compute_marginal_objective gradient is correct", {
  data_list <- JointODE:::.process_marginal(
    observed ~ x1,
    marginal_sim$longitudinal_data,
    "time", "id",
    marginal_sim$state
  )

  theta <- c(-0.5, -0.3, 0.2)

  result <- .compute_marginal_objective(
    theta, data_list,
    gradient = TRUE, hessian = FALSE
  )
  analytic_grad <- attr(result, "gradient")

  numeric_grad <- numDeriv::grad(
    function(x) {
      as.numeric(.compute_marginal_objective(
        x, data_list,
        gradient = FALSE, hessian = FALSE
      ))
    },
    theta
  )

  expect_equal(analytic_grad, numeric_grad, tolerance = 1e-4)
})

test_that(".solve_batch_marginal returns correct structure", {
  data_list <- JointODE:::.process_marginal(
    observed ~ x1,
    marginal_sim$longitudinal_data,
    "time", "id",
    marginal_sim$state
  )
  theta <- c(-0.5, -0.3, 0.2)

  sols <- .solve_batch_marginal(data_list, theta)

  expect_type(sols, "list")
  expect_equal(length(sols), length(data_list))
  sol <- sols[[1]]
  expect_named(sol, c("times", "biomarker", "velocity", "acceleration"))
  expect_true(all(is.finite(sol$biomarker)))
})

test_that(".compute_marginal_state gradient is correct", {
  data_list <- JointODE:::.process_marginal(
    observed ~ x1,
    marginal_sim$longitudinal_data,
    "time", "id",
    marginal_sim$state
  )

  theta <- c(-0.5, -0.3, 0.2)
  initial_state <- c(0.5, 0.1)

  result <- .compute_marginal_state(
    initial_state, data_list[[1]], theta,
    gradient = TRUE, hessian = FALSE
  )
  analytic_grad <- attr(result, "gradient")

  numeric_grad <- numDeriv::grad(
    function(x) {
      as.numeric(.compute_marginal_state(
        x, data_list[[1]], theta,
        gradient = FALSE, hessian = FALSE
      ))
    },
    initial_state
  )

  expect_equal(analytic_grad, numeric_grad, tolerance = 1e-4)
})
