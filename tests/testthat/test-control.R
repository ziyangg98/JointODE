# ==============================================================================
# JointODE.control() Unit Tests
# ==============================================================================

test_that("JointODE.control returns complete defaults", {
  ctrl <- JointODE.control()

  expect_type(ctrl, "list")
  expect_equal(ctrl$maxit, 200)
  expect_equal(ctrl$tol, 1e-4)
  expect_equal(ctrl$verbose, 0)
  expect_false(ctrl$parallel)
  expect_equal(ctrl$n_cores, 0)
  expect_equal(ctrl$quad_level, 4)
})

test_that("JointODE.control overrides individual parameters", {
  ctrl <- JointODE.control(maxit = 50)
  expect_equal(ctrl$maxit, 50)
  expect_equal(ctrl$tol, 1e-4) # other defaults unchanged

  ctrl2 <- JointODE.control(tol = 1e-6)
  expect_equal(ctrl2$tol, 1e-6)
})

test_that("JointODE.control .list merges with defaults", {
  ctrl <- JointODE.control(.list = list(maxit = 100, verbose = TRUE))
  expect_equal(ctrl$maxit, 100)
  expect_equal(ctrl$verbose, 1) # TRUE -> 1
  expect_equal(ctrl$tol, 1e-4) # default preserved
})

test_that("JointODE.control validates maxit", {
  expect_error(JointODE.control(maxit = 0), "positive")
  expect_error(JointODE.control(maxit = -1), "positive")
})

test_that("JointODE.control validates tol", {
  expect_error(JointODE.control(tol = 0), "positive")
  expect_error(JointODE.control(tol = -0.1), "positive")
})

test_that("JointODE.control validates verbose", {
  ctrl <- JointODE.control(verbose = TRUE)
  expect_equal(ctrl$verbose, 1)

  ctrl2 <- JointODE.control(verbose = FALSE)
  expect_equal(ctrl2$verbose, 0)

  ctrl3 <- JointODE.control(verbose = 2)
  expect_equal(ctrl3$verbose, 2)

  # verbose="yes" coerced to NA by as.numeric, which is accepted
  expect_warning(JointODE.control(verbose = "yes"))
})

test_that("JointODE.control validates parallel", {
  expect_error(JointODE.control(parallel = "yes"), "TRUE or FALSE")
})

test_that("JointODE.control validates n_cores", {
  expect_error(JointODE.control(n_cores = -1), "non-negative")
})

test_that("JointODE.control validates quad_level", {
  expect_error(JointODE.control(quad_level = 0), "positive")
})


test_that("JointODE.control validates .list type", {
  expect_error(JointODE.control(.list = "bad"), "\\.list must be a list")
})

test_that("JointODE.control passes through extra parameters via ...", {
  ctrl <- JointODE.control(custom_param = 42)
  expect_equal(ctrl$custom_param, 42)
})
