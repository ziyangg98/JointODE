# ==============================================================================
# Helper Functions Unit Tests
# ==============================================================================

# --- .expand_posteriors ---

test_that(".expand_posteriors expands correctly", {
  # Construct minimal posteriors structure
  posteriors <- list(
    nodes = list(
      matrix(rnorm(6), nrow = 3, ncol = 2),
      matrix(rnorm(4), nrow = 2, ncol = 2)
    ),
    weights = list(
      c(0.3, 0.5, 0.2),
      c(0.6, 0.4)
    )
  )
  data_list <- list(
    list(id = 1, value = "a"),
    list(id = 2, value = "b")
  )

  expanded <- .expand_posteriors(data_list, posteriors)

  expect_equal(nrow(expanded$nodes), 5) # 3 + 2
  expect_equal(length(expanded$weights), 5)
  expect_equal(length(expanded$data), 5)
  expect_equal(expanded$node_to_subject, c(1, 1, 1, 2, 2))

  # Data replicated correctly
  expect_equal(expanded$data[[1]]$value, "a")
  expect_equal(expanded$data[[3]]$value, "a")
  expect_equal(expanded$data[[4]]$value, "b")
})

test_that(".expand_posteriors preserves weights", {
  posteriors <- list(
    nodes = list(
      matrix(0, nrow = 2, ncol = 1),
      matrix(0, nrow = 3, ncol = 1)
    ),
    weights = list(
      c(0.4, 0.6),
      c(0.2, 0.3, 0.5)
    )
  )
  data_list <- list(list(id = 1), list(id = 2))

  expanded <- .expand_posteriors(data_list, posteriors)
  expect_equal(
    expanded$weights,
    c(0.4, 0.6, 0.2, 0.3, 0.5)
  )
})

# --- .extend_covariates ---

test_that(".extend_covariates does LOCF interpolation", {
  cov_mat <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3)
  orig_times <- c(0, 1, 2)
  pred_times <- c(0, 0.5, 1, 1.5, 2, 2.5)

  result <- .extend_covariates(cov_mat, orig_times, pred_times)

  expect_equal(nrow(result), length(pred_times))
  expect_equal(ncol(result), 2)

  # At original times, values unchanged
  expect_equal(result[1, ], cov_mat[1, ])
  expect_equal(result[3, ], cov_mat[2, ])
  expect_equal(result[5, ], cov_mat[3, ])

  # Between times, uses LOCF (last observation carried forward)
  expect_equal(result[2, ], cov_mat[1, ]) # t=0.5 -> uses t=0
  expect_equal(result[4, ], cov_mat[2, ]) # t=1.5 -> uses t=1
  expect_equal(result[6, ], cov_mat[3, ]) # t=2.5 -> uses t=2
})

test_that(".extend_covariates handles NULL/empty input", {
  result <- .extend_covariates(NULL, c(0, 1), c(0, 0.5, 1))
  expect_null(result)

  empty <- matrix(numeric(0), nrow = 3, ncol = 0)
  result2 <- .extend_covariates(empty, c(0, 1, 2), c(0, 0.5))
  expect_equal(nrow(result2), 2)
  expect_equal(ncol(result2), 0)
})

# --- .format_vector ---

test_that(".format_vector formats correctly", {
  expect_equal(.format_vector(numeric(0)), "[]")
  expect_true(grepl("^\\[", .format_vector(c(1.0, 2.0))))
  expect_true(grepl("\\+", .format_vector(1:10, n = 3)))
})

# --- .format_matrix ---

test_that(".format_matrix formats correctly", {
  expect_equal(.format_matrix(NULL), "[]")
  result <- .format_matrix(diag(2))
  expect_true(grepl("1.000", result))
})

# --- .compute_metrics ---

test_that(".compute_metrics computes correctly at iter > 1", {
  curr <- list(
    loglik = -100,
    parameters = list(coefficients = list(
      baseline = c(0.1, 0.2),
      hazard = c(0.5),
      longitudinal = c(-0.01),
      initial_state = c(0, 0),
      measurement_error_sd = 0.5,
      random_effect_sigma = diag(4)
    ))
  )
  prev <- list(
    loglik = -110,
    parameters = list(coefficients = list(
      baseline = c(0.1, 0.2),
      hazard = c(0.5),
      longitudinal = c(-0.01),
      initial_state = c(0, 0),
      measurement_error_sd = 0.5,
      random_effect_sigma = diag(4)
    ))
  )

  m <- .compute_metrics(curr, prev, iter = 2)
  expect_equal(m$delta_l, 10)
  expect_true(m$rel_l > 0)
})

test_that(".compute_metrics at iter 1", {
  curr <- list(loglik = -100)
  m <- .compute_metrics(curr, NULL, iter = 1)
  expect_equal(m$delta_l, -100)
})
