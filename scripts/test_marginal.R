devtools::load_all()
data("sim", package = "JointODE")

cat("=== Without biomarker/velocity fixed effects ===\n")
t1 <- system.time(fit1 <- MarginalODE(
  observed ~ x1 + x2,
  data = sim$data$longitudinal_data,
  control = list(verbose = 0)
))
cat(sprintf("Time: %.1f s | logLik: %.2f | sigma_e: %.4f\n",
            t1["elapsed"], fit1$logLik, fit1$measurement_error_sd))

cat("\n=== With biomarker/velocity fixed effects ===\n")
t2 <- system.time(fit2 <- MarginalODE(
  observed ~ biomarker + velocity + x1 + x2,
  data = sim$data$longitudinal_data,
  control = list(verbose = 0)
))
cat(sprintf("Time: %.1f s | logLik: %.2f | sigma_e: %.4f\n",
            t2["elapsed"], fit2$logLik, fit2$measurement_error_sd))

cat("\n=== With biomarker/velocity fixed + random effects ===\n")
t3 <- system.time(fit3 <- MarginalODE(
  observed ~ biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
  data = sim$data$longitudinal_data,
  control = list(verbose = 0)
))
cat(sprintf("Time: %.1f s | logLik: %.2f | sigma_e: %.4f\n",
            t3["elapsed"], fit3$logLik, fit3$measurement_error_sd))
