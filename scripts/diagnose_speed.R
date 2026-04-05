library(JointODE)
data(sim)

n_test <- 200
ids <- unique(sim$data$longitudinal_data$id)[seq_len(n_test)]
long_data <- sim$data$longitudinal_data[
  sim$data$longitudinal_data$id %in% ids, c("id", "time", "observed", "x1", "x2")]
surv_data <- sim$data$survival_data[sim$data$survival_data$id %in% ids, ]

# Fit once to get the TMB object
TMB::openmp(1)
fit <- JointODE(
  observed ~ biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
  Surv(time, status) ~ w1 + w2, long_data, surv_data,
  init = sim$init, control = list(verbose = 0, maxit = 500))
obj <- fit$tmb_obj
par0 <- fit$tmb_opt$par

# 1. Hessian sparsity
cat("=== Random effects Hessian structure ===\n")
n_re_total <- sum(names(obj$env$last.par) == "b")
cat("Total random effects:", n_re_total, "\n")
cat("Per subject:", n_re_total / n_test, "\n")
h <- obj$env$spHess(obj$env$last.par, random = TRUE)
cat("Hessian dim:", dim(h), "\n")
cat("Hessian nnz:", length(h@x), "\n")
cat("Density:", round(length(h@x) / prod(dim(h)) * 100, 2), "%\n")
cat("If block-diagonal (4x4):", n_test * 16, "nnz expected\n\n")

# 2. Profile fn/gr
cat("=== fn/gr timing (30 evals) ===\n")
n_eval <- 30
t_fn <- system.time(for (j in seq_len(n_eval)) obj$fn(par0))[3]
t_gr <- system.time(for (j in seq_len(n_eval)) obj$gr(par0))[3]
cat("fn:", round(t_fn / n_eval, 4), "s/eval\n")
cat("gr:", round(t_gr / n_eval, 4), "s/eval\n")
cat("Total per outer iter:", round((t_fn + t_gr) / n_eval, 4), "s\n")
cat("Estimated outer iters:", fit$tmb_opt$iterations[1], "\n")
cat("Estimated total:", round((t_fn + t_gr) / n_eval * fit$tmb_opt$iterations[1], 1), "s\n\n")

# 3. Inner Newton control tuning
cat("=== Inner Newton tuning ===\n")
for (maxit_inner in c(5, 10, 20, 50, 100)) {
  obj$env$inner.control <- list(maxit = maxit_inner)
  t_fn_i <- system.time(for (j in seq_len(10)) obj$fn(par0))[3]
  v_i <- obj$fn(par0)
  cat(sprintf("  inner_maxit=%3d: fn=%.4fs  nll=%.4f\n",
              maxit_inner, t_fn_i / 10, v_i))
}

# Reset
obj$env$inner.control <- list()

# 4. Compare full optimization with inner tuning
cat("\n=== Full optimization with inner_maxit=20 ===\n")
# Re-create to get fresh object
fit2_obj <- TMB::MakeADFun(
  data = obj$env$data, parameters = obj$env$parList(),
  random = "b", DLL = "JointODE", silent = TRUE)
fit2_obj$env$inner.control <- list(maxit = 20)
t0 <- proc.time()
opt2 <- nlminb(fit2_obj$par, fit2_obj$fn, fit2_obj$gr,
               control = list(iter.max = 500, eval.max = 5000, rel.tol = 1e-4))
t2 <- (proc.time() - t0)[3]
cat("Time:", round(t2, 1), "s  (vs", round(fit$tmb_opt$evaluations[1], 1),
    "fn evals in original)\n")
cat("LL:", -opt2$objective, " (original:", fit$logLik, ")\n")
cat("Converged:", opt2$convergence == 0, "\n")
