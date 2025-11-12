# ==============================================================================
# Helper Functions for Gradient Tests
# ==============================================================================

#' Compare gradient components with detailed diagnostics
#' @param analytical Analytical gradient vector
#' @param numerical Numerical gradient vector
#' @param name Component name for error messages
#' @param tolerance Comparison tolerance (default: 1e-4)
compare_gradient_component <- function(
  analytical,
  numerical,
  name,
  tolerance = 1e-2
) {
  rel_error <- max(abs((analytical - numerical) / (abs(numerical) + 1e-10)))

  # expect_lt is from testthat package (loaded in test environment)
  # nolint next: object_usage_linter
  expect_lt(rel_error, tolerance, label = paste0(name, " gradient"))
}

#' Load and process test data from sim dataset
#' @param n_subjects Number of subjects to use (NULL = all)
#' @return List with data, parameters, and n_subjects
load_test_data <- function(n_subjects = NULL) {
  # sim is loaded from package data (data/sim.rda)
  # nolint start: object_usage_linter
  # Process data using sim dataset
  data_processed <- .process(
    longitudinal_formula = observed ~ x1 + x2,
    survival_formula = Surv(time, status) ~ w1 + w2,
    longitudinal_data = sim$data$longitudinal_data,
    survival_data = sim$data$survival_data,
    state = NULL
  )

  parameters <- sim$init

  # Subset if requested for faster testing
  if (!is.null(n_subjects)) {
    n_subjects <- min(n_subjects, length(data_processed))
    data_processed <- data_processed[1:n_subjects]
    random_effects <- sim$data$random_effects[1:n_subjects, , drop = FALSE]
  }
  # nolint end

  list(
    data = data_processed,
    parameters = parameters,
    random_effects = random_effects
  )
}
