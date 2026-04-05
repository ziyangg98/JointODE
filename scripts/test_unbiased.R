library(JointODE)
library(parallel)

n_rep <- 20
n_subjects <- 200
n_cores_fit <- 2  # OpenMP threads per fit

# True parameters
ref <- JointODE:::.create_example_data(n_subjects = 100, seed = 0)
true_coef <- with(ref$init$coefficients,
  c(baseline, hazard, longitudinal, initial_state))
true_sigma_e <- ref$init$coefficients$measurement_error_sd

cat("True parameters:\n")
print(round(true_coef, 4))
cat("True sigma_e:", true_sigma_e, "\n\n")

fit_one <- function(r) {
  sim_data <- JointODE:::.create_example_data(n_subjects = n_subjects, seed = r * 100)
  tryCatch({
    fit <- JointODE(
      observed ~ biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
      Surv(time, status) ~ w1 + w2,
      sim_data$data$longitudinal_data[, c("id", "time", "observed", "x1", "x2")],
      sim_data$data$survival_data,
      init = sim_data$init,
      control = list(verbose = 0, maxit = 500,
                     parallel = n_cores_fit > 1, n_cores = n_cores_fit)
    )
    list(coef = coef(fit),
         sigma_e = fit$parameters$coefficients$measurement_error_sd,
         cindex = fit$cindex,
         converged = fit$convergence$converged)
  }, error = function(e) NULL)
}

# Run reps in parallel across R processes
cat(sprintf("Running %d reps x %d subjects (%d cores per fit, %d parallel R processes)\n",
            n_rep, n_subjects, n_cores_fit,
            max(1, detectCores() - 1) %/% max(n_cores_fit, 1)))

n_workers <- max(1, (detectCores() - 1) %/% max(n_cores_fit, 1))
cl <- makeCluster(n_workers)
clusterEvalQ(cl, library(JointODE))
clusterExport(cl, c("n_subjects", "n_cores_fit"))
t0 <- proc.time()
results <- parLapply(cl, seq_len(n_rep), fit_one)
stopCluster(cl)
elapsed <- (proc.time() - t0)[3]

cat(sprintf("Total time: %.0f sec (%.0f sec/rep)\n\n", elapsed, elapsed / n_rep))

# Summary
ok <- !vapply(results, is.null, logical(1))
n_ok <- sum(ok)
cat(sprintf("========== SUMMARY: %d/%d successful ==========\n", n_ok, n_rep))

if (n_ok > 0) {
  coef_mat <- do.call(rbind, lapply(results[ok], `[[`, "coef"))
  sigma_es <- vapply(results[ok], `[[`, numeric(1), "sigma_e")
  cindexes <- vapply(results[ok], `[[`, numeric(1), "cindex")

  cat("Converged:", sum(vapply(results[ok], `[[`, logical(1), "converged")),
      "/", n_ok, "\n")
  cat("C-index: mean =", round(mean(cindexes), 3),
      " sd =", round(sd(cindexes), 3), "\n")
  cat("sigma_e: mean =", round(mean(sigma_es), 4),
      " sd =", round(sd(sigma_es), 4),
      " (true:", true_sigma_e, ")\n\n")

  bias <- colMeans(coef_mat) - true_coef
  se_est <- apply(coef_mat, 2, sd)

  tab <- data.frame(
    true = round(true_coef, 4),
    mean = round(colMeans(coef_mat), 4),
    bias = round(bias, 4),
    sd = round(se_est, 4),
    bias_se = round(bias / (se_est / sqrt(n_ok)), 2)
  )
  cat("Parameter estimates (|bias_se| > 2 suggests bias):\n")
  print(tab)
}
