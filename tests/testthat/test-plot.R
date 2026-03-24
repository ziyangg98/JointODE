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

# --- All plot types should produce gg objects ---

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

test_that("plot.JointODE type=diagnostic_residuals_time works", {
  p <- .safe_plot(mock, type = "diagnostic_residuals_time")
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

test_that("plot.JointODE type=diagnostic_association errors", {
  expect_error(
    .safe_plot(mock, type = "diagnostic_association"),
    "not yet implemented"
  )
})

test_that("plot.JointODE invalid type errors", {
  expect_error(plot(mock, type = "nonexistent"), "arg")
})

# --- subject_ids triggers faceted view ---

test_that("plot.JointODE subject_ids produces faceted biomarker plot", {
  ids <- names(mock$data)[1:3]
  p <- .safe_plot(mock, type = "trajectory_biomarker", subject_ids = ids)
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE subject_ids produces faceted velocity plot", {
  ids <- names(mock$data)[1:2]
  p <- .safe_plot(mock, type = "trajectory_velocity", subject_ids = ids)
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE subject_ids produces faceted phase plot", {
  ids <- names(mock$data)[1:2]
  p <- .safe_plot(mock, type = "phase_biomarker_velocity", subject_ids = ids)
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE subject_ids produces faceted vel_acc plot", {
  ids <- names(mock$data)[1:2]
  p <- .safe_plot(mock, type = "phase_velocity_acceleration", subject_ids = ids)
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE subject_ids produces faceted survival plot", {
  ids <- names(mock$data)[1:2]
  p <- .safe_plot(mock, type = "survival", subject_ids = ids)
  expect_s3_class(p, "gg")
})

# --- show_observed=FALSE ---

test_that("plot.JointODE show_observed=FALSE hides points in faceted view", {
  ids <- names(mock$data)[1:2]
  p <- .safe_plot(mock,
    type = "trajectory_biomarker",
    subject_ids = ids, show_observed = FALSE
  )
  expect_s3_class(p, "gg")
})

# --- show_individual=FALSE uses mean-only overlay ---

test_that("plot.JointODE show_individual=FALSE for biomarker", {
  p <- .safe_plot(mock, type = "trajectory_biomarker", show_individual = FALSE)
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE show_individual=FALSE for phase", {
  p <- .safe_plot(mock,
    type = "phase_biomarker_velocity",
    show_individual = FALSE
  )
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE show_individual=FALSE for vel_acc", {
  p <- .safe_plot(mock,
    type = "phase_velocity_acceleration",
    show_individual = FALSE
  )
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE show_individual=FALSE for survival", {
  p <- .safe_plot(mock, type = "survival", show_individual = FALSE)
  expect_s3_class(p, "gg")
})

# --- overview ---

test_that("plot.JointODE overview works", {
  p <- .safe_plot(mock, type = "overview")
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

# --- by grouping triggers grouped plot path ---

test_that("plot.JointODE trajectory_biomarker with by grouping", {
  p <- .safe_plot(mock, type = "trajectory_biomarker", by = "w1")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE trajectory_velocity with by grouping", {
  p <- .safe_plot(mock, type = "trajectory_velocity", by = "w1")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE phase with by grouping", {
  p <- .safe_plot(mock, type = "phase_biomarker_velocity", by = "w1")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE vel_acc with by grouping", {
  p <- .safe_plot(mock, type = "phase_velocity_acceleration", by = "w1")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE survival with by grouping", {
  p <- .safe_plot(mock, type = "survival", by = "w1")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE by invalid variable errors", {
  expect_error(
    .safe_plot(mock, type = "survival", by = "nonexistent_var"),
    "not found"
  )
})

# --- random effects edge cases ---

test_that("plot.JointODE random_effects NULL errors", {
  mock_null_re <- mock
  mock_null_re$random_effects <- NULL
  expect_error(
    .safe_plot(mock_null_re, type = "diagnostic_random_effects"),
    "No random effects"
  )
})

test_that("plot.JointODE random_effects as vector", {
  mock_vec_re <- mock
  mock_vec_re$random_effects <- rnorm(length(mock$data))
  p <- .safe_plot(mock_vec_re, type = "diagnostic_random_effects")
  expect_s3_class(p, "gg")
})

test_that("plot.JointODE random_effects as data.frame", {
  mock_df_re <- mock
  mock_df_re$random_effects <- as.data.frame(mock$random_effects$estimates)
  p <- .safe_plot(mock_df_re, type = "diagnostic_random_effects")
  expect_s3_class(p, "gg")
})
