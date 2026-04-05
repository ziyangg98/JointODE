library(JointODE)
data(sim)

n_test <- 50
n_rep <- 10

# True values from sim$init
init_coef <- with(sim$init$coefficients,
  c(baseline, hazard, longitudinal, initial_state))

results <- vector("list", n_rep)

for (r in seq_len(n_rep)) {
  # Random subset of subjects
  all_ids <- unique(sim$data$longitudinal_data$id)
  set.seed(r * 42)
  test_ids <- sample(all_ids, n_test)

  long_data <- sim$data$longitudinal_data[
    sim$data$longitudinal_data$id %in% test_ids,
    c("id", "time", "observed", "x1", "x2")
  ]
  surv_data <- sim$data$survival_data[
    sim$data$survival_data$id %in% test_ids,
  ]

  fit <- tryCatch(
    JointODE(
      longitudinal_formula = observed ~ biomarker + velocity + x1 + x2 +
        (biomarker + velocity | id),
      survival_formula = Surv(time, status) ~ w1 + w2,
      longitudinal_data = long_data,
      survival_data = surv_data,
      init = sim$init,
      control = list(verbose = 0, maxit = 500)
    ),
    error = function(e) { cat("Rep", r, "FAILED:", e$message, "\n"); NULL }
  )

  if (!is.null(fit)) {
    results[[r]] <- list(
      coef = coef(fit),
      sigma_e = fit$parameters$coefficients$measurement_error_sd,
      cindex = fit$cindex,
      converged = fit$convergence$converged,
      loglik = fit$logLik
    )
    cat(sprintf("Rep %2d: converged=%s  logLik=%.1f  C=%.3f  sigma_e=%.4f\n",
                r, fit$convergence$converged, fit$logLik, fit$cindex,
                fit$parameters$coefficients$measurement_error_sd))
  }
}

# Summary table
ok <- !vapply(results, is.null, logical(1))
coef_mat <- do.call(rbind, lapply(results[ok], `[[`, "coef"))

cat("\n=== Summary across", sum(ok), "successful fits ===\n")
cat("\nTrue values (init):\n")
names(init_coef) <- colnames(coef_mat)
print(round(init_coef, 4))

cat("\nMean estimated:\n")
print(round(colMeans(coef_mat), 4))

cat("\nSD of estimates:\n")
print(round(apply(coef_mat, 2, sd), 4))

cat("\nMean bias (est - true):\n")
print(round(colMeans(coef_mat) - init_coef, 4))

cat("\nsigma_e: mean =", mean(vapply(results[ok], `[[`, numeric(1), "sigma_e")),
    " sd =", sd(vapply(results[ok], `[[`, numeric(1), "sigma_e")), "\n")
cat("C-index: mean =", mean(vapply(results[ok], `[[`, numeric(1), "cindex")),
    " sd =", sd(vapply(results[ok], `[[`, numeric(1), "cindex")), "\n")
cat("Converged:", sum(vapply(results[ok], `[[`, logical(1), "converged")),
    "/", sum(ok), "\n")
