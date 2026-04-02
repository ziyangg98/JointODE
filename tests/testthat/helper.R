# ==============================================================================
# Shared Test Helpers
# ==============================================================================

# nolint start: object_usage_linter

# Load required packages
library(survival)

# Load bundled test data
data("sim", package = "JointODE", envir = environment())

#' Subset sim data into processed test data
.make_test_data <- function(n = 10) {
  test_ids <- unique(sim$data$longitudinal_data$id)[seq_len(n)]
  data_list <- JointODE:::.process_joint(
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
  list(
    data_list = data_list,
    random_effects = sim$data$random_effects[seq_len(n), , drop = FALSE],
    params = JointODE:::.coef_to_vector(sim$init),
    parameters = sim$init
  )
}

#' Create mock JointODE fitted object for S3 method tests
.create_mock_jointode <- function(n_subjects = 10L) {
  td <- .make_test_data(n_subjects)
  parameters <- td$parameters

  parameters$coefficients$baseline <- setNames(
    parameters$coefficients$baseline,
    paste0("bs", seq_along(parameters$coefficients$baseline))
  )
  parameters$coefficients$hazard <- setNames(
    parameters$coefficients$hazard,
    c("value", "slope", "w1", "w2")
  )

  n_params <- .count_params(parameters)
  n_re <- nrow(parameters$coefficients$random_effect_sigma)
  coef_names <- c(
    paste0("baseline:", names(parameters$coefficients$baseline)),
    paste0("hazard:", names(parameters$coefficients$hazard)),
    paste0("longitudinal:", names(parameters$coefficients$longitudinal)),
    paste0("initial state:", names(parameters$coefficients$initial_state))
  )
  vcov_matrix <- diag(0.01, length(coef_names))
  dimnames(vcov_matrix) <- list(coef_names, coef_names)

  structure(
    list(
      parameters = parameters,
      logLik = -500.0,
      AIC = 1000.0 + 2 * n_params,
      BIC = 1000.0 + n_params * log(n_subjects),
      cindex = 0.65,
      convergence = list(
        converged = TRUE, iterations = 10,
        message = "Converged after 10 iterations"
      ),
      random_effects = td$random_effects,
      vcov = vcov_matrix,
      data = td$data_list,
      control = JointODE.control(),
      call = quote(JointODE(
        longitudinal_formula = observed ~ biomarker + velocity +
          x1 + x2 + (biomarker + velocity | id),
        survival_formula = Surv(time, status) ~ w1 + w2
      ))
    ),
    class = "JointODE"
  )
}

#' Central finite difference gradient helper
.finite_diff_gradient <- function(func, x, eps = 1e-5) {
  grad <- numeric(length(x))
  for (i in seq_along(x)) {
    x_plus <- x
    x_minus <- x
    x_plus[i] <- x_plus[i] + eps
    x_minus[i] <- x_minus[i] - eps
    grad[i] <- (func(x_plus) - func(x_minus)) / (2 * eps)
  }
  grad
}

#' Zero random effects matrix with same shape as td$random_effects
.zero_random_effects <- function(td) {
  matrix(0, nrow = nrow(td$random_effects), ncol = ncol(td$random_effects))
}

#' Expect all values finite
.expect_all_finite <- function(x) {
  expect_true(all(is.finite(x)))
}

# nolint end
