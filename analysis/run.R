library(JointODE)

dat <- read.csv("analysis/data/num_rdw.csv")
dat$time <- dat$time * 100

fit <- MarginalODE(
  observed ~ biomarker + velocity + x1 + x2 +
    (1 + biomarker + velocity | id),
  data = dat,
  control = MarginalODE.control(
    maxit = 200,
    verbose = 2,
    parallel = TRUE,
    n_cores = 4
  )
)

saveRDS(fit, "analysis/num_rdw_direct_fit.rds")

pred <- predict(fit)
rmse <- sqrt(mean((dat$observed - pred$biomarker)^2))

cat("converged:", fit$convergence$converged, "\n")
cat("message:", fit$convergence$message, "\n")
cat("logLik:", fit$logLik, "\n")
cat("AIC:", fit$AIC, "\n")
cat("BIC:", fit$BIC, "\n")
cat("RMSE:", rmse, "\n")
print(summary(fit))
