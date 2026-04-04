# Unit tests: ODE analytical solutions for each branch

# Analytical solution for m''(t) = b1*m(t) + b2*m'(t) + f with constant f
.analytical_ode_step <- function(m0, v0, b1, b2, f, dt) {
  eps_b1 <- 1e-8
  eps_d <- 1e-12

  D <- b2^2 + 4 * b1

  if (abs(b1) < eps_b1) {
    if (abs(b2) < eps_b1) { # ZERO: constant acceleration m''=f
      m <- m0 + v0 * dt + (f / 2) * dt^2
      v <- v0 + f * dt
    } else { # FIRST_ORD: v' = b2*v + f
      eb <- exp(b2 * dt)
      m <- m0 + (v0 + f / b2) * (eb - 1) / b2 - f * dt / b2
      v <- (v0 + f / b2) * eb - f / b2
    }
  } else if (D > eps_d) { # REAL: two real roots
    s <- sqrt(D)
    h <- b2 / 2
    l1 <- h + s / 2
    l2 <- h - s / 2
    e1 <- exp(l1 * dt)
    e2 <- exp(l2 * dt)

    a0 <- (l1 * e2 - l2 * e1) / s
    a1 <- (e1 - e2) / s
    F1 <- (e1 - 1) / l1
    F2 <- (e2 - 1) / l2
    J1 <- (F1 - F2) / s
    J0 <- (l1 * F2 - l2 * F1) / s

    m <- a0 * m0 + a1 * v0 + f * J1
    v <- b1 * a1 * m0 + (a0 + b2 * a1) * v0 + f * (J0 + b2 * J1)
  } else if (D < -eps_d) { # COMPLEX: conjugate roots a ± iw
    h <- b2 / 2
    w <- sqrt(abs(D)) / 2
    r2 <- h^2 + w^2
    ehat <- exp(h * dt)
    cost <- cos(w * dt)
    sint <- sin(w * dt)

    a0 <- ehat * (cost - h * sint / w)
    a1 <- ehat * sint / w
    Ic <- (ehat * (h * cost + w * sint) - h) / r2
    Is <- (ehat * (h * sint - w * cost) + w) / r2
    J1 <- Is / w
    J0 <- Ic - (h / w) * Is

    m <- a0 * m0 + a1 * v0 + f * J1
    v <- b1 * a1 * m0 + (a0 + b2 * a1) * v0 + f * (J0 + b2 * J1)
  } else { # REPEATED: double root at b2/2
    h <- b2 / 2
    eh <- exp(h * dt)
    a0 <- eh * (1.0 - h * dt)
    a1 <- dt * eh
    Fh <- (eh - 1.0) / h
    J1 <- (dt * eh - Fh) / h
    J0 <- 2.0 * Fh - dt * eh

    m <- a0 * m0 + a1 * v0 + f * J1
    v <- b1 * a1 * m0 + (a0 + b2 * a1) * v0 + f * (J0 + b2 * J1)
  }

  list(m = m, v = v)
}

# Test: Analytical solutions for each ODE branch
test_that("analytical solutions are finite and directionally sane", {
  test_cases <- list(
    REAL = list(b1 = -2.0, b2 = 1.0),
    FIRST_ORD = list(b1 = 0, b2 = 0.5),
    COMPLEX = list(b1 = -2.0, b2 = 0.1),
    REPEATED = list(b1 = -0.249999, b2 = 1.0),
    ZERO = list(b1 = 0, b2 = 0)
  )

  m0 <- 1.0
  v0 <- 0.5
  f <- 0.1
  dt <- 0.05

  for (branch_name in names(test_cases)) {
    case <- test_cases[[branch_name]]
    sol <- .analytical_ode_step(m0, v0, case$b1, case$b2, f, dt)

    # Verify solution is finite
    expect_true(is.finite(sol$m) && is.finite(sol$v))

    # Verify direction of movement (sanity check)
    dm <- sol$m - m0
    expect_equal(dm >= 0, v0 >= 0)
  }
})

# Test: ZERO branch limiting case (constant acceleration)
test_that("ZERO branch recovers constant acceleration", {
  m0 <- 0.0
  v0 <- 1.0
  f <- 2.0
  dt <- 0.1
  b1 <- 0
  b2 <- 0

  sol <- .analytical_ode_step(m0, v0, b1, b2, f, dt)

  # Expected: m(t) = m0 + v0*t + (f/2)*t²
  m_expected <- m0 + v0 * dt + (f / 2) * dt^2
  # Expected velocity
  v_expected <- v0 + f * dt

  expect_equal(sol$m, m_expected, tolerance = 1e-10)
  expect_equal(sol$v, v_expected, tolerance = 1e-10)
})

# Test: FIRST_ORD limiting case (when b1 → 0)
test_that("FIRST_ORD branch matches first-order ODE limit", {
  m0 <- 1.0
  v0 <- 0.5
  f <- 0.2
  b2 <- 0.5
  dt <- 0.1
  b1_tiny <- 1e-10 # << 1e-8 threshold

  sol <- .analytical_ode_step(m0, v0, b1_tiny, b2, f, dt)

  # Expected velocity (first-order ODE): v(t) = (v0 + f/b2)*e^{b2*t} - f/b2
  v_expected <- (v0 + f / b2) * exp(b2 * dt) - f / b2

  expect_equal(sol$v, v_expected, tolerance = 1e-8)
})

# Test: Numerical stability with large parameters
test_that("analytical solver stays finite under extreme settings", {
  stress_cases <- list(
    list(m0 = 1.0, v0 = 0.5, b1 = -100.0, b2 = 10.0, f = 0.1, dt = 0.01),
    list(m0 = 1.0, v0 = 0.5, b1 = -0.5, b2 = 0.5, f = 0.1, dt = 1e-6),
    list(m0 = 1.0, v0 = 0.5, b1 = -0.25, b2 = 1.000001, f = 0.1, dt = 0.1)
  )

  for (case in stress_cases) {
    sol <- .analytical_ode_step(
      case$m0, case$v0, case$b1, case$b2, case$f, case$dt
    )
    expect_true(all(is.finite(c(sol$m, sol$v))))
  }
})
