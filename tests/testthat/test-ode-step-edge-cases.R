test_that("R ODE step agrees with numerical integration in edge regimes", {
  de_step <- function(m, v, b1, b2, f, dt) {
    sol <- deSolve::ode(
      y = c(m = m, v = v),
      times = c(0, dt),
      func = function(t, y, parms) {
        list(c(y[["v"]], b1 * y[["m"]] + b2 * y[["v"]] + f))
      },
      parms = NULL,
      rtol = 1e-12,
      atol = 1e-12
    )
    unname(c(sol[2, "m"], sol[2, "v"]))
  }

  cases <- list(
    underdamped = c(b1 = -4, b2 = -0.4, dt = 1.2),
    repeated = c(b1 = -0.25, b2 = -1, dt = 2),
    near_repeated = c(b1 = -0.25000001, b2 = -1, dt = 2),
    weak_restoring = c(b1 = -1e-12, b2 = -0.5, dt = 3),
    first_order = c(b1 = 0, b2 = -0.8, dt = 3),
    zero_dynamics = c(b1 = 0, b2 = 0, dt = 1),
    strong_overdamped = c(b1 = -1e-4, b2 = -100, dt = 2)
  )

  for (case in cases) {
    got <- JointODE:::.ode_step_r(
      m = 0.7, v = -0.2,
      b1 = case[["b1"]], b2 = case[["b2"]],
      f = 0.15, dt = case[["dt"]]
    )
    ref <- de_step(
      m = 0.7, v = -0.2,
      b1 = case[["b1"]], b2 = case[["b2"]],
      f = 0.15, dt = case[["dt"]]
    )
    expect_true(all(is.finite(got)))
    expect_equal(got, ref, tolerance = 1e-6)
  }
})

test_that("R ODE step is continuous near removable singularities", {
  zero <- JointODE:::.ode_step_r(0.7, -0.2, 0, -0.5, 0.15, 3)
  weak <- JointODE:::.ode_step_r(0.7, -0.2, -1e-12, -0.5, 0.15, 3)
  expect_equal(weak, zero, tolerance = 1e-8)

  left <- JointODE:::.ode_step_r(0.7, -0.2, -0.25 * (1 - 1e-7), -1, 0.15, 2)
  right <- JointODE:::.ode_step_r(0.7, -0.2, -0.25 * (1 + 1e-7), -1, 0.15, 2)
  expect_equal(left, right, tolerance = 1e-6)

  below <- JointODE:::.ode_step_r(0.7, -0.2, -1e-4, -25, 0.15, 2)
  above <- JointODE:::.ode_step_r(0.7, -0.2, -1e-4, -25.000001, 0.15, 2)
  expect_equal(below, above, tolerance = 1e-6)
})
