library(JointODE)
data(sim)
ld <- sim$data$longitudinal_data[, c("id", "time", "observed", "x1", "x2")]

cat("=== MarginalODE with sinhc/sinc ODE solver ===\n")
t0 <- proc.time()
fit <- MarginalODE(observed ~ x1 + x2, data = ld,
                   control = list(verbose = 1))
elapsed <- (proc.time() - t0)[3]
cat(sprintf("Elapsed: %.1f sec\n", elapsed))

cat("\n=== summary ===\n")
print(summary(fit))

cat("\n=== coef ===\n")
print(coef(fit))
