# ==============================================================================
# Output Module: S3 Methods, Plot, Predict Helpers
# ==============================================================================

mock <- .create_mock_jointode()

# --- print / summary ---

test_that("print.JointODE works", {
  out <- capture.output(result <- print(mock))
  expect_identical(result, mock)
  expect_true(any(grepl("Joint ODE Model", out)))
  expect_true(any(grepl("Log-likelihood", out)))
})

test_that("summary.JointODE returns correct structure", {
  s <- summary(mock)
  expect_s3_class(s, "summary.JointODE")
  expect_equal(ncol(s$coef_baseline), 4)
  expect_equal(ncol(s$coef_survival), 4)
  expect_equal(ncol(s$coef_longitudinal), 4)
  expect_equal(
    nrow(s$coef_baseline),
    length(mock$parameters$coefficients$baseline)
  )
  expect_true(s$nobs > 0)
  expect_true(s$event_rate >= 0 && s$event_rate <= 1)
})

test_that("summary.JointODE includes derived params", {
  s <- summary(mock)
  if (mock$parameters$coefficients$longitudinal[1] < 0) {
    expect_equal(nrow(s$derived_params), 2)
  }
})

test_that("print.summary.JointODE works", {
  out <- capture.output(result <- print(summary(mock)))
  expect_true(any(grepl("Longitudinal Process", out)))
  expect_true(any(grepl("Survival Process", out)))
})

test_that("summary.JointODE handles NULL vcov", {
  mock_nv <- mock
  mock_nv$vcov <- NULL
  expect_true(all(is.na(summary(mock_nv)$coef_baseline[, "Std. Error"])))
})

# --- coef / vcov / logLik ---

test_that("coef.JointODE returns named vector", {
  cf <- coef(mock)
  pat <- "^(baseline|hazard|longitudinal|initial state):"
  expect_true(all(grepl(pat, names(cf))))
  expect_equal(
    unname(cf[seq_along(mock$parameters$coefficients$baseline)]),
    unname(mock$parameters$coefficients$baseline)
  )
})

test_that("vcov.JointODE returns symmetric matrix", {
  v <- vcov(mock)
  expect_equal(nrow(v), ncol(v))
  expect_equal(v, t(v))
})

test_that("logLik.JointODE returns logLik class", {
  ll <- logLik(mock)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), mock$logLik)
  expect_equal(attr(ll, "df"), .count_params(mock$parameters))
})

# ==============================================================================
# predict
# ==============================================================================

test_that("predict.JointODE returns data.frame", {
  pred <- predict(mock)
  expect_s3_class(pred, "data.frame")
  expected <- c("id", "time", "cumhaz", "survival", "biomarker")
  expect_true(all(expected %in% names(pred)))
})

test_that("predict.JointODE errors on newdata", {
  expect_error(predict(mock, newdata = data.frame()), "not yet supported")
})

# ==============================================================================
# Predict Helpers
# ==============================================================================

test_that(".extend_covariates does LOCF interpolation", {
  cov_mat <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3)
  result <- .extend_covariates(cov_mat, c(0, 1, 2), c(0, 0.5, 1, 1.5, 2, 2.5))
  expect_equal(nrow(result), 6)
  expect_equal(result[1, ], cov_mat[1, ])
  expect_equal(result[2, ], cov_mat[1, ]) # LOCF
  expect_equal(result[6, ], cov_mat[3, ]) # LOCF
})

test_that(".extend_covariates handles NULL input", {
  expect_null(.extend_covariates(NULL, c(0, 1), c(0, 0.5, 1)))
})

# ==============================================================================
# Format Helpers
# ==============================================================================

test_that(".format_vector formats correctly", {
  expect_equal(.format_vector(numeric(0)), "[]")
  expect_true(grepl("^\\[", .format_vector(c(1.0, 2.0))))
  expect_true(grepl("\\+", .format_vector(1:10, n = 3)))
})

test_that(".format_matrix formats correctly", {
  expect_equal(.format_matrix(NULL), "[]")
  expect_true(grepl("1.000", .format_matrix(diag(2))))
})

# ==============================================================================
# Plot
# ==============================================================================

.safe_plot <- function(...) {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  plot(...)
}

# --- All plot types ---

test_that("plot.JointODE all types produce gg", {
  types <- c(
    "trajectory_biomarker", "trajectory_velocity",
    "phase_biomarker_velocity", "phase_velocity_acceleration",
    "survival", "diagnostic_residuals", "diagnostic_residuals_time",
    "diagnostic_qq"
  )
  for (type in types) {
    expect_s3_class(.safe_plot(mock, type = type), "gg")
  }
})

test_that("plot.JointODE random_effects needs matrix conversion", {
  mock_re <- mock
  mock_re$random_effects <- mock$random_effects$estimates
  expect_s3_class(
    .safe_plot(mock_re, type = "diagnostic_random_effects"), "gg"
  )
})

test_that("plot.JointODE invalid type errors", {
  expect_error(plot(mock, type = "nonexistent"), "arg")
  expect_error(.safe_plot(mock, type = "diagnostic_association"), "not yet")
})

# --- subject_ids faceting ---

test_that("plot.JointODE subject_ids produces faceted plots", {
  ids <- names(mock$data)[1:2]
  types <- c(
    "trajectory_biomarker", "trajectory_velocity",
    "phase_biomarker_velocity", "phase_velocity_acceleration", "survival"
  )
  for (type in types) {
    expect_s3_class(
      .safe_plot(mock, type = type, subject_ids = ids), "gg"
    )
  }
})

# --- show_observed / show_individual ---

test_that("plot.JointODE show_observed=FALSE works", {
  ids <- names(mock$data)[1:2]
  p <- .safe_plot(
    mock, type = "trajectory_biomarker",
    subject_ids = ids, show_observed = FALSE
  )
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE show_individual=FALSE works", {
  types <- c(
    "trajectory_biomarker", "phase_biomarker_velocity",
    "phase_velocity_acceleration", "survival"
  )
  for (type in types) {
    p <- .safe_plot(mock, type = type, show_individual = FALSE)
    expect_s3_class(p, "gg")
  }
})

# --- overview ---

test_that("plot.JointODE overview works", {
  p <- .safe_plot(mock, type = "overview")
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

# --- by grouping ---

test_that("plot.JointODE by grouping works", {
  types <- c(
    "trajectory_biomarker", "trajectory_velocity",
    "phase_biomarker_velocity", "phase_velocity_acceleration", "survival"
  )
  for (type in types) {
    p <- .safe_plot(mock, type = type, by = "w1")
    expect_s3_class(p, "gg")
  }
  expect_error(
    .safe_plot(mock, type = "survival", by = "nonexistent"),
    "not found"
  )
})

# --- random effects edge cases ---

test_that("plot.JointODE random_effects edge cases", {
  mock_null <- mock
  mock_null$random_effects <- NULL
  expect_error(
    .safe_plot(mock_null, type = "diagnostic_random_effects"),
    "No random"
  )

  mock_vec <- mock
  mock_vec$random_effects <- rnorm(length(mock$data))
  expect_s3_class(
    .safe_plot(mock_vec, type = "diagnostic_random_effects"), "gg"
  )

  mock_df <- mock
  mock_df$random_effects <- as.data.frame(
    mock$random_effects$estimates
  )
  expect_s3_class(
    .safe_plot(mock_df, type = "diagnostic_random_effects"), "gg"
  )
})
