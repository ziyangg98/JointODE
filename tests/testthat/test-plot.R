# ==============================================================================
# plot.JointODE Unit Tests
# ==============================================================================

mock <- create_mock_jointode()

# Use pdf(NULL) to suppress actual plot rendering in tests
.safe_plot <- function(...) {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  plot(...)
}

test_that("plot.JointODE type=trajectory_biomarker works", {
  p <- .safe_plot(mock, type = "trajectory_biomarker")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE type=trajectory_velocity works", {
  p <- .safe_plot(mock, type = "trajectory_velocity")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE type=phase_biomarker_velocity works", {
  p <- .safe_plot(mock, type = "phase_biomarker_velocity")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE type=phase_velocity_acceleration works", {
  p <- .safe_plot(mock, type = "phase_velocity_acceleration")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE type=survival works", {
  p <- .safe_plot(mock, type = "survival")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE type=diagnostic_residuals works", {
  p <- .safe_plot(mock, type = "diagnostic_residuals")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE type=diagnostic_qq works", {
  p <- .safe_plot(mock, type = "diagnostic_qq")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE type=diagnostic_random_effects works", {
  mock_re <- mock
  mock_re$random_effects <- mock$random_effects$estimates
  p <- .safe_plot(mock_re, type = "diagnostic_random_effects")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE invalid type errors", {
  expect_error(
    plot(mock, type = "nonexistent"),
    "arg"
  )
})

test_that("plot.JointODE subject_ids filters subjects", {
  ids <- names(mock$data)[1:3]
  p <- .safe_plot(mock, type = "trajectory_biomarker", subject_ids = ids)
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE overview works", {
  p <- .safe_plot(mock, type = "overview")
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})
