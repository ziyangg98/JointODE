library(JointODE)
data(sim)

long_data <- sim$data$longitudinal_data[, c("id", "time", "observed", "x1", "x2")]

cat("=== MCEM E2E Test ===\n")
fit <- JointODE(
  longitudinal_formula = observed ~ biomarker + velocity + x1 + x2 +
    (biomarker + velocity | id),
  survival_formula = Surv(time, status) ~ w1 + w2,
  longitudinal_data = long_data,
  survival_data = sim$data$survival_data,
  init = sim$init,
  control = list(maxit = 10, mc_samples = 50, verbose = 2)
)
cat("\n=== Summary ===\n")
print(summary(fit))
