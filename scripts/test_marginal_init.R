devtools::load_all()
data("sim", package = "JointODE")

long_data <- sim$data$longitudinal_data[, c("id", "time", "observed", "x1", "x2")]
surv_data <- sim$data$survival_data

cat("=== Step 1: Marginal init ===\n")
t0 <- proc.time()
fit <- JointODE(
  observed ~ biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
  Surv(time, status) ~ w1 + w2,
  long_data, surv_data,
  init = "marginal",
  control = list(verbose = 2, parallel = TRUE)
)
elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("\nTotal time: %.1f s\n", elapsed))
cat(sprintf("logLik: %.2f | C-index: %.3f\n", fit$logLik, fit$cindex))
cat("\nCoefficients:\n")
print(round(coef(fit), 4))

cat("\n\n=== Step 2: Default init (for comparison) ===\n")
t0 <- proc.time()
fit2 <- JointODE(
  observed ~ biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
  Surv(time, status) ~ w1 + w2,
  long_data, surv_data,
  init = "default",
  control = list(verbose = 2, parallel = TRUE)
)
elapsed2 <- (proc.time() - t0)["elapsed"]
cat(sprintf("\nTotal time: %.1f s\n", elapsed2))
cat(sprintf("logLik: %.2f | C-index: %.3f\n", fit2$logLik, fit2$cindex))
cat("\nCoefficients:\n")
print(round(coef(fit2), 4))

cat("\n=== Summary ===\n")
cat(sprintf("Marginal init: %.1f s, logLik=%.2f, C=%.3f\n",
            elapsed, fit$logLik, fit$cindex))
cat(sprintf("Default init:  %.1f s, logLik=%.2f, C=%.3f\n",
            elapsed2, fit2$logLik, fit2$cindex))
