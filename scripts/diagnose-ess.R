#!/usr/bin/env Rscript
# Diagnose low ESS subjects in MCEM importance sampling
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

# Run one E-step to get IS diagnostics
M <- 100
random_effects <- sim$data$random_effects

mc_post <- JointODE:::.compute_posteriors_mcem(
  data_list, parameters, random_effects, M,
  parallel = FALSE, n_cores = 0
)

# ESS per subject
ess <- vapply(mc_post, function(p) 1 / sum(p$weights^2), numeric(1))
cat(sprintf("ESS: mean=%.1f, sd=%.1f, min=%.1f, max=%.1f\n",
            mean(ess), sd(ess), min(ess), max(ess)))

# Find low ESS subjects (< 50)
low_idx <- which(ess < 50)
cat(sprintf("\n%d subjects with ESS < 50:\n", length(low_idx)))

for (i in low_idx) {
  d <- data_list[[i]]
  n_obs <- length(d$longitudinal$measurements)
  t_range <- range(d$longitudinal$times)
  event <- d$survival$status

  # Weight concentration: max weight / uniform weight
  max_w <- max(mc_post[[i]]$weights)
  w_sorted <- sort(mc_post[[i]]$weights, decreasing = TRUE)
  top5_mass <- sum(w_sorted[1:5])

  # Laplace mode vs warm-start
  mode <- mc_post[[i]]$mode
  init <- random_effects[i, ]
  mode_shift <- sqrt(sum((mode - init)^2))

  # Posterior Hessian condition number
  cov_mat <- mc_post[[i]]$cov
  eig <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
  cond <- max(eig) / max(min(eig), .Machine$double.eps)

  # Log-posterior at mode vs at samples — check for multimodality
  log_post_samples <- as.numeric(JointODE:::.compute_joint_logpost_batch(
    samples = mc_post[[i]]$samples, data = d, parameters = parameters
  ))
  log_post_mode <- as.numeric(JointODE:::.compute_joint_logpost(
    random_effect = mode, data = d,
    parameters = parameters, gradient = FALSE, hessian = FALSE
  ))
  # How many samples have log_post within 2 units of mode?
  near_mode <- sum(log_post_samples > log_post_mode - 2)

  cat(sprintf(
    "\n  Subject %d: ESS=%.1f, n_obs=%d, t=[%.1f,%.1f], event=%d\n",
    i, ess[i], n_obs, t_range[1], t_range[2], event
  ))
  cat(sprintf(
    "    max_weight=%.3f, top5_mass=%.3f, mode_shift=%.4f\n",
    max_w, top5_mass, mode_shift
  ))
  cat(sprintf(
    "    cov_cond=%.1f, eig_range=[%.2e, %.2e]\n",
    cond, min(eig), max(eig)
  ))
  cat(sprintf(
    "    near_mode_samples=%d/%d, logpost_range=[%.1f, %.1f]\n",
    near_mode, M, min(log_post_samples), max(log_post_samples)
  ))
  cat(sprintf("    mode: %s\n", paste(sprintf("%.4f", mode), collapse=", ")))
  cat(sprintf("    init: %s\n", paste(sprintf("%.4f", init), collapse=", ")))

  # Check which RE dimension has the most weight variation
  w <- mc_post[[i]]$weights
  for (d_re in seq_len(ncol(mc_post[[i]]$samples))) {
    weighted_mean <- sum(w * mc_post[[i]]$samples[, d_re])
    weighted_var <- sum(w * (mc_post[[i]]$samples[, d_re] - weighted_mean)^2)
    unweighted_var <- var(mc_post[[i]]$samples[, d_re])
    cat(sprintf("    RE[%d]: w_mean=%.4f, w_var=%.4f, q_var=%.4f, ratio=%.2f\n",
                d_re, weighted_mean, weighted_var, unweighted_var,
                weighted_var / unweighted_var))
  }
}

# Summary: correlate ESS with subject characteristics
n_obs_vec <- vapply(data_list, function(d) length(d$longitudinal$measurements), integer(1))
event_vec <- vapply(data_list, function(d) d$survival$status %||% 0L, integer(1))
t_max_vec <- vapply(data_list, function(d) d$survival$event_time %||% 0, numeric(1))
cov_cond_vec <- vapply(mc_post, function(p) {
  eig <- eigen(p$cov, symmetric = TRUE, only.values = TRUE)$values
  max(eig) / max(min(eig), .Machine$double.eps)
}, numeric(1))

cat("\n\nCorrelation with ESS:\n")
cat(sprintf("  n_obs: r=%.3f\n", cor(ess, n_obs_vec)))
cat(sprintf("  event: r=%.3f\n", cor(ess, event_vec)))
cat(sprintf("  t_max: r=%.3f\n", cor(ess, t_max_vec)))
cat(sprintf("  log(cov_cond): r=%.3f\n", cor(ess, log(cov_cond_vec))))
