library(JointODE)

cat("=== Sim Dataset Test ===\n\n")

data(sim)
longitudinal_data <- sim$data$longitudinal_data[
  , c("id", "time", "observed", "x1", "x2")
]

cat(sprintf(
  "Subjects: %d | Obs: %d | Events: %d (%.0f%%)\n",
  nrow(sim$data$survival_data),
  nrow(longitudinal_data),
  sum(sim$data$survival_data$status),
  100 * mean(sim$data$survival_data$status)
))

cat("\n--- MarginalODE ---\n")
t0 <- proc.time()
fit_marginal <- MarginalODE(
  observed ~ x1 + x2,
  longitudinal_data,
  time = "time", id = "id",
  control = list(verbose = 2, parallel = TRUE)
)
cat(sprintf("Elapsed: %.1f s\n", (proc.time() - t0)["elapsed"]))
summary(fit_marginal)

cat("\n--- JointODE (init = marginal) ---\n")
t0 <- proc.time()
fit <- JointODE(
  longitudinal_formula = observed ~ biomarker + velocity + x1 + x2 +
    (biomarker + velocity | id),
  survival_formula = Surv(time, status) ~ w1 + w2,
  longitudinal_data = longitudinal_data,
  survival_data = sim$data$survival_data,
  init = "marginal",
  control = list(verbose = 2, maxit = 200, parallel = TRUE)
)
elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("\nTotal elapsed: %.1f s\n", elapsed))

cat("\n--- Summary ---\n")
summary(fit)

cat("\n--- Longitudinal SSE ---\n")
pred <- predict(fit)
long <- longitudinal_data
merged <- merge(pred, long, by.x = c("id", "time"), by.y = c("id", "time"))
sse <- sum((merged$biomarker - merged$observed)^2)
cat(sprintf("Total SSE: %.4f  (n=%d, MSE=%.6f, RMSE=%.4f)\n",
  sse, nrow(merged), sse / nrow(merged), sqrt(sse / nrow(merged))))

# Per-subject SSE
subj_sse <- tapply(
  (merged$biomarker - merged$observed)^2, merged$id, sum
)
cat(sprintf("Per-subject SSE: median=%.4f, mean=%.4f, max=%.4f\n",
  median(subj_sse), mean(subj_sse), max(subj_sse)))

cat("\n--- Save ---\n")
saveRDS(fit, "tests/fit_sim.rds")
cat("Saved to tests/fit_sim.rds\n")

cat("\nDone.\n")
