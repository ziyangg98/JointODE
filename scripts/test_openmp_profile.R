library(JointODE)
data(sim)

n_test <- 200
ids <- unique(sim$data$longitudinal_data$id)[seq_len(n_test)]
long_data <- sim$data$longitudinal_data[
  sim$data$longitudinal_data$id %in% ids, c("id", "time", "observed", "x1", "x2")]
surv_data <- sim$data$survival_data[sim$data$survival_data$id %in% ids, ]

run_fit <- function(nt, label) {
  TMB::openmp(nt)
  cat("  config nthreads before MakeADFun:", TMB::config()$nthreads, "\n")
  t0 <- proc.time()
  fit <- JointODE(
    observed ~ biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
    Surv(time, status) ~ w1 + w2, long_data, surv_data,
    init = sim$init, control = list(verbose = 0, maxit = 500,
                                     parallel = nt > 1, n_cores = nt))
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("  %s: %.1f sec (LL=%.1f)\n", label, elapsed, fit$logLik))
  elapsed
}

cat("=== OpenMP Benchmark (200 subjects) ===\n\n")
t1 <- run_fit(1, "1 thread")
t2 <- run_fit(2, "2 threads")
t4 <- run_fit(4, "4 threads")
cat(sprintf("\nSpeedup: 2T=%.2fx  4T=%.2fx\n", t1/t2, t1/t4))
