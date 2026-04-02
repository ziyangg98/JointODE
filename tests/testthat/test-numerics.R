# ==============================================================================
# Numerical Utility Tests
# ==============================================================================

test_that(".safe_chol regularizes indefinite matrices", {
  H <- matrix(c(1, 2, 2, 1), 2, 2)

  R <- JointODE:::.safe_chol(H)
  H_reg <- t(R) %*% R

  ev <- eigen(H_reg, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))
  expect_true(all(is.finite(H_reg)))
})

test_that(".safe_chol returns standard chol for SPD matrices", {
  H <- matrix(c(2, 0.3, 0.3, 1.5), 2, 2)

  expect_equal(JointODE:::.safe_chol(H), chol(H), tolerance = 1e-12)
})
