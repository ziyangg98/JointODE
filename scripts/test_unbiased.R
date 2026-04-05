library(JointODE)
library(parallel)

n_rep <- 20
n_subjects <- 200
n_cores <- 4

# Get true parameters
ref <- JointODE:::.create_example_data(n_subjects = 100, seed = 0)
true_coef <- with(ref$init$coefficients,
  c(baseline, hazard, longitudinal, initial_state))
true_sigma_e <- ref$init$coefficients$measurement_error_sd

cat("True parameters:\n")
print(round(true_coef, 4))
cat("True sigma_e:", true_sigma_e, "\n\n")

# Single fit function
fit_one <- function(r) {
  sim_data <- JointODE:::.create_example_data(
    n_subjects = n_subjects, seed = r * 100
  )
  tryCatch({
    fit <- JointODE(
      observed ~ biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
      Surv(time, status) ~ w1 + w2,
      sim_data$data$longitudinal_data[, c("id", "time", "observed", "x1", "x2")],
      sim_data$data$survival_data,
      init = sim_data$init,
      control = list(verbose = 0, maxit = 500,
                     parallel = TRUE, n_cores = n_cores)
    )
    list(
      coef = coef(fit),
      sigma_e = fit$parameters$coefficients$measurement_error_sd,
      cindex = fit$cindex,
      converged = fit$convergence$converged
    )
  }, error = function(e) {
    cat("Rep", r, "FAILED:", e$message, "\n")
    NULL
  })
}

# Test OpenMP first with a single fit
cat("=== Testing OpenMP with", n_cores, "cores ===\n")
t0 <- proc.time()
test_fit <- fit_one(999)
t_openmp <- (proc.time() - t0)[3]
cat("Single fit (OpenMP):", round(t_openmp, 1), "sec\n")

if (is.null(test_fit)) stop("Test fit failed")
cat("Test fit converged:", test_fit$converged, "\n\n")

# Compare with single-thread
cat("=== Testing single-thread ===\n")
t0 <- proc.time()
sim_data <- JointODE:::.create_example_data(n_subjects = n_subjects, seed = 99900)
fit_single <- JointODE(
  observed ~ biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
  Surv(time, status) ~ w1 + w2,
  sim_data$data$longitudinal_data[, c("id", "time", "observed", "x1", "x2")],
  sim_data$data$survival_data,
  init = sim_data$init,
  control = list(verbose = 0, maxit = 500, parallel = FALSE)
)
t_single <- (proc.time() - t0)[3]
cat("Single fit (1 thread):", round(t_single, 1), "sec\n")
cat("Speedup:", round(t_single / t_openmp, 2), "x\n\n")

# Run all reps sequentially (each fit uses OpenMP internally)
cat("=== Running", n_rep, "replications ===\n")
results <- vector("list", n_rep)
for (r in seq_len(n_rep)) {
  results[[r]] <- fit_one(r)
  if (!is.null(results[[r]])) {
    cat(sprintf("Rep %2d: conv=%s  C=%.3f  sigma_e=%.4f\n",
                r, results[[r]]$converged, results[[r]]$cindex,
                results[[r]]$sigma_e))
  }
}

# Summary
ok <- !vapply(results, is.null, logical(1))
n_ok <- sum(ok)
cat(sprintf("\n========== SUMMARY: %d/%d successful ==========\n", n_ok, n_rep))

if (n_ok > 0) {
  coef_mat <- do.call(rbind, lapply(results[ok], `[[`, "coef"))
  sigma_es <- vapply(results[ok], `[[`, numeric(1), "sigma_e")
  cindexes <- vapply(results[ok], `[[`, numeric(1), "cindex")

  cat("\nConverged:", sum(vapply(results[ok], `[[`, logical(1), "converged")),
      "/", n_ok, "\n")
  cat("C-index: mean =", round(mean(cindexes), 3),
      " sd =", round(sd(cindexes), 3), "\n")
  cat("sigma_e: mean =", round(mean(sigma_es), 4),
      " sd =", round(sd(sigma_es), 4),
      " (true:", true_sigma_e, ")\n")

  bias <- colMeans(coef_mat) - true_coef
  se_est <- apply(coef_mat, 2, sd)

  tab <- data.frame(
    true = round(true_coef, 4),
    mean = round(colMeans(coef_mat), 4),
    bias = round(bias, 4),
    sd = round(se_est, 4),
    bias_se = round(bias / (se_est / sqrt(n_ok)), 2)
  )
  cat("\nParameter estimates (|bias_se| > 2 suggests bias):\n")
  print(tab)
}
