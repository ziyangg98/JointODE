# Test script: initial states as random effects
# Run on feature/state-as-re branch

devtools::load_all()
data(sim)

# True parameter values (from simulate defaults)
cat("=== True Values ===\n")
cat("initial_state:", sim$init$coefficients$initial_state, "\n")
cat("hazard:", sim$init$coefficients$hazard, "\n")
cat("sigma_e:", sim$init$coefficients$measurement_error_sd, "\n")
cat("sigma_b diag:", round(diag(sim$init$coefficients$random_effect_sigma), 4), "\n")
cat("RE dims:", dim(sim$data$random_effects), "\n\n")

# Prepare data (drop reserved columns)
longitudinal_data <- sim$data$longitudinal_data[
  , c("id", "time", "observed", "x1", "x2")
]

# Fit model
cat("=== Fitting JointODE (200 subjects, marginal init) ===\n")
t0 <- Sys.time()
fit <- JointODE(
  longitudinal_formula = observed ~ biomarker + velocity + x1 + x2 +
    (biomarker + velocity | id),
  survival_formula = Surv(time, status) ~ w1 + w2,
  longitudinal_data = longitudinal_data,
  survival_data = sim$data$survival_data,
  init = "marginal",
  control = list(maxit = 50, verbose = 1, parallel = TRUE)
)
elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 0)

# Results
cat(sprintf("\n=== Results (elapsed: %s sec) ===\n", elapsed))
cat("Converged:", fit$convergence$converged, "\n")
cat("Iterations:", fit$convergence$iterations, "\n")
cat("LogLik:", fit$logLik, "\n")
cat("C-index:", fit$cindex, "\n\n")

# Parameter comparison
cat("=== Parameter Recovery ===\n")
cf <- coef(fit)
true_hazard <- c(0.8, 2.0, 0.6, -0.8)
cat(sprintf("%-25s %10s %10s\n", "Parameter", "True", "Estimated"))
cat(sprintf("%-25s %10.4f %10.4f\n", "hazard:alpha_1", true_hazard[1], cf["hazard:alpha_1"]))
cat(sprintf("%-25s %10.4f %10.4f\n", "hazard:alpha_2", true_hazard[2], cf["hazard:alpha_2"]))
cat(sprintf("%-25s %10.4f %10.4f\n", "hazard:w1", true_hazard[3], cf["hazard:w1"]))
cat(sprintf("%-25s %10.4f %10.4f\n", "hazard:w2", true_hazard[4], cf["hazard:w2"]))
cat(sprintf("%-25s %10.4f %10.4f\n", "longitudinal:x1", 0.5, cf["longitudinal:x1"]))
cat(sprintf("%-25s %10.4f %10.4f\n", "longitudinal:x2", -0.45, cf["longitudinal:x2"]))
cat(sprintf("%-25s %10.4f %10.4f\n", "initial_state:m0", -0.5, cf["initial_state:m0"]))
cat(sprintf("%-25s %10.4f %10.4f\n", "initial_state:v0", -0.1, cf["initial_state:v0"]))

# Variance components
cat("\n=== Variance Components ===\n")
cat("sigma_e:", round(fit$parameters$coefficients$measurement_error_sd, 4),
    "(true: 0.1)\n")
cat("\nRandom Effect Sigma (estimated):\n")
print(round(fit$parameters$coefficients$random_effect_sigma, 4))
cat("\nRandom Effect Sigma (true):\n")
print(round(sim$init$coefficients$random_effect_sigma, 4))

# Derived ODE parameters
cat("\n=== Derived ODE Parameters ===\n")
b1 <- cf["longitudinal:biomarker"]
b2 <- cf["longitudinal:velocity"]
omega <- sqrt(-b1)
period <- 2 * pi / omega
xi <- -b2 / (2 * omega)
cat(sprintf("Period:  %.3f (true: 6.0)\n", period))
cat(sprintf("Damping: %.3f (true: 0.4)\n", xi))

# Summary
cat("\n=== Full Summary ===\n")
summary(fit)
