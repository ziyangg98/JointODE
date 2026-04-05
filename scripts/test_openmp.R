library(JointODE)
data(sim)

n_test <- 200
ids <- unique(sim$data$longitudinal_data$id)[seq_len(n_test)]
long_data <- sim$data$longitudinal_data[
  sim$data$longitudinal_data$id %in% ids, c("id", "time", "observed", "x1", "x2")]
surv_data <- sim$data$survival_data[sim$data$survival_data$id %in% ids, ]

run_fit <- function(n_threads) {
  TMB::openmp(n_threads)
  t0 <- proc.time()
  fit <- JointODE(
    observed ~ biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
    Surv(time, status) ~ w1 + w2,
    long_data, surv_data,
    init = sim$init,
    control = list(verbose = 0, maxit = 500)
  )
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("  %d thread(s): %.1f sec  (LL=%.1f, conv=%s)\n",
              n_threads, elapsed, fit$logLik, fit$convergence$converged))
  elapsed
}

cat(sprintf("Benchmark: %d subjects\n", n_test))
t1 <- run_fit(1)
t2 <- run_fit(2)
t4 <- run_fit(4)
t8 <- run_fit(8)

cat(sprintf("\nSpeedup: 2T=%.2fx  4T=%.2fx  8T=%.2fx\n", t1/t2, t1/t4, t1/t8))
