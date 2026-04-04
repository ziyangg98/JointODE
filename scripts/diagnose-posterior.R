#!/usr/bin/env Rscript
# Diagnose posterior shape at true parameter values
library(JointODE)
library(survival)

data("sim", package = "JointODE")
data_list <- JointODE:::.process_joint(
  longitudinal_data = sim$data$longitudinal_data[,
    c("id", "time", "observed", "x1", "x2")],
  longitudinal_formula = observed ~
    biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
  survival_data = sim$data$survival_data,
  survival_formula = Surv(time, status) ~ w1 + w2
)
parameters <- sim$init
true_re <- sim$data$random_effects

n <- length(data_list)
n_re <- ncol(true_re)
re_names <- colnames(true_re)
if (is.null(re_names)) re_names <- paste0("RE", seq_len(n_re))

# For each subject: compute Laplace mode, compare to true RE,
# and profile the posterior along each RE dimension
cat("=== Posterior shape at true values ===\n\n")

# Pick a few representative subjects
set.seed(42)
sample_ids <- sort(sample(n, min(5, n)))

for (i in sample_ids) {
  d <- data_list[[i]]

  # Laplace approximation starting from true RE
  post <- JointODE:::.compute_posterior_laplace(d, true_re[i, ], parameters)

  # Log-posterior at true RE, at Laplace mode, and at zero
  lp_true <- as.numeric(JointODE:::.compute_joint_logpost(
    random_effect = true_re[i, ], data = d,
    parameters = parameters, gradient = FALSE, hessian = FALSE
  ))
  lp_mode <- as.numeric(JointODE:::.compute_joint_logpost(
    random_effect = post$mode, data = d,
    parameters = parameters, gradient = FALSE, hessian = FALSE
  ))
  lp_zero <- as.numeric(JointODE:::.compute_joint_logpost(
    random_effect = rep(0, n_re), data = d,
    parameters = parameters, gradient = FALSE, hessian = FALSE
  ))

  cat(sprintf("Subject %d (n_obs=%d, event=%d):\n",
              i, length(d$longitudinal$measurements), d$survival$status))
  cat(sprintf("  true RE: %s\n", paste(sprintf("%8.4f", true_re[i, ]), collapse = "")))
  cat(sprintf("  mode:    %s\n", paste(sprintf("%8.4f", post$mode), collapse = "")))
  cat(sprintf("  diff:    %s\n", paste(sprintf("%8.4f", post$mode - true_re[i, ]), collapse = "")))
  cat(sprintf("  logpost: true=%.2f, mode=%.2f, zero=%.2f\n",
              lp_true, lp_mode, lp_zero))

  # Posterior covariance diagnostics
  eig <- eigen(post$cov, symmetric = TRUE, only.values = TRUE)$values
  cat(sprintf("  cov diag: %s\n", paste(sprintf("%.2e", diag(post$cov)), collapse = ", ")))
  cat(sprintf("  cov cond: %.1f\n", max(eig) / max(min(eig), .Machine$double.eps)))

  # Profile posterior: sweep each RE dimension, fix others at mode
  cat("  Posterior profiles (centered at mode, conditional SD):\n")
  prec <- solve(post$cov)
  for (d_re in seq_len(n_re)) {
    cond_sd <- 1 / sqrt(prec[d_re, d_re])
    grid <- post$mode[d_re] + cond_sd * seq(-3, 3, by = 0.5)
    lp_grid <- vapply(grid, function(val) {
      b <- post$mode
      b[d_re] <- val
      as.numeric(JointODE:::.compute_joint_logpost(
        random_effect = b, data = d,
        parameters = parameters, gradient = FALSE, hessian = FALSE
      ))
    }, numeric(1))

    # Compare to conditional Gaussian: -0.5 * prec[d,d] * (x - mode)^2
    lp_gauss <- lp_mode - 0.5 * prec[d_re, d_re] * (grid - post$mode[d_re])^2
    max_dev <- max(abs(lp_grid - lp_gauss))

    # Skewness: compare logpost at +2*cond_sd vs -2*cond_sd
    idx_neg2 <- which.min(abs(grid - (post$mode[d_re] - 2 * cond_sd)))
    idx_pos2 <- which.min(abs(grid - (post$mode[d_re] + 2 * cond_sd)))
    skew_indicator <- (lp_grid[idx_pos2] - lp_mode) - (lp_grid[idx_neg2] - lp_mode)

    cat(sprintf("    %s: cond_sd=%.4f, max_gauss_dev=%.2f, skew=%.2f, true_at=%.1f*sd\n",
                re_names[d_re], cond_sd, max_dev, skew_indicator,
                (true_re[i, d_re] - post$mode[d_re]) / cond_sd))
  }
  cat("\n")
}

# Summary across all subjects: mode vs true
cat("=== Summary: Laplace mode vs true RE (all subjects) ===\n")
modes <- t(vapply(seq_len(n), function(i) {
  post <- JointODE:::.compute_posterior_laplace(data_list[[i]], true_re[i, ], parameters)
  post$mode
}, numeric(n_re)))

bias <- colMeans(modes - true_re)
rmse <- sqrt(colMeans((modes - true_re)^2))
cor_diag <- vapply(seq_len(n_re), function(j) cor(modes[, j], true_re[, j]), numeric(1))

cat(sprintf("  %-12s %8s %8s %8s\n", "RE", "bias", "RMSE", "cor"))
for (j in seq_len(n_re)) {
  cat(sprintf("  %-12s %8.4f %8.4f %8.4f\n", re_names[j], bias[j], rmse[j], cor_diag[j]))
}
