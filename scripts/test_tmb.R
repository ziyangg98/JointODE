library(JointODE)
data(sim)
TMB::openmp(8)
# Use a small subset for quick test
n_test <- 200
test_ids <- unique(sim$data$longitudinal_data$id)[seq_len(n_test)]

long_data <- sim$data$longitudinal_data[
  sim$data$longitudinal_data$id %in% test_ids,
  c("id", "time", "observed", "x1", "x2")
]
surv_data <- sim$data$survival_data[
  sim$data$survival_data$id %in% test_ids,
]

cat("Fitting JointODE (TMB) with", n_test, "subjects...\n")
fit <- JointODE(
  longitudinal_formula = observed ~ biomarker + velocity + x1 + x2 +
    (biomarker + velocity | id),
  survival_formula = Surv(time, status) ~ w1 + w2,
  longitudinal_data = long_data,
  survival_data = surv_data,
  init = sim$init,
  control = list(verbose = 3, maxit = 500)
)

cat("\n=== TMB Fitted Parameters ===\n")
cat("\nFixed effects (coef):\n")
print(coef(fit))

cat("\nsigma_e:", fit$parameters$coefficients$measurement_error_sd, "\n")
cat("\nSigma_b:\n")
print(round(fit$parameters$coefficients$random_effect_sigma, 4))

cat("\n=== Initial (sim$init) Parameters ===\n")
cat("\nFixed effects:\n")
init_coef <- with(sim$init$coefficients,
  c(baseline, hazard, longitudinal, initial_state))
names(init_coef) <- names(coef(fit))
print(init_coef)

cat("\nsigma_e:", sim$init$coefficients$measurement_error_sd, "\n")
cat("\nSigma_b:\n")
print(round(sim$init$coefficients$random_effect_sigma, 4))

cat("\n=== Comparison (TMB - init) ===\n")
diff <- coef(fit) - init_coef
print(round(diff, 4))

cat("\nRandom effects (first 5 subjects):\n")
print(round(fit$random_effects[1:5, ], 4))

cat("\nVcov diagonal (SE):\n")
print(round(sqrt(diag(fit$vcov)), 4))
