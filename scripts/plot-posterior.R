#!/usr/bin/env Rscript
# Plot posterior profiles for representative subjects
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
n_re <- ncol(true_re)
re_names <- c("init_biomarker", "init_velocity", "biomarker", "velocity")

set.seed(42)
sample_ids <- sort(sample(length(data_list), 4))

pdf("scripts/posterior_profiles.pdf", width = 12, height = 12)
par(mfrow = c(4, 4), mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))

for (i in sample_ids) {
  d <- data_list[[i]]
  post <- JointODE:::.compute_posterior_laplace(d, true_re[i, ], parameters)
  lp_mode <- as.numeric(JointODE:::.compute_joint_logpost(
    random_effect = post$mode, data = d,
    parameters = parameters, gradient = FALSE, hessian = FALSE
  ))

  prec <- solve(post$cov)
  for (d_re in seq_len(n_re)) {
    cond_sd <- 1 / sqrt(prec[d_re, d_re])
    grid <- post$mode[d_re] + cond_sd * seq(-4, 4, length.out = 81)
    lp_grid <- vapply(grid, function(val) {
      b <- post$mode
      b[d_re] <- val
      as.numeric(JointODE:::.compute_joint_logpost(
        random_effect = b, data = d,
        parameters = parameters, gradient = FALSE, hessian = FALSE
      ))
    }, numeric(1))

    # Laplace conditional Gaussian: -0.5 * precision[d,d] * (x - mode)^2
    lp_gauss <- lp_mode - 0.5 * prec[d_re, d_re] * (grid - post$mode[d_re])^2

    # Plot
    ylim <- range(c(lp_grid, lp_gauss), na.rm = TRUE)
    plot(grid, lp_grid, type = "l", lwd = 2, col = "black",
         xlab = re_names[d_re], ylab = "log posterior",
         main = sprintf("Subject %d: %s", i, re_names[d_re]),
         ylim = ylim)
    lines(grid, lp_gauss, lwd = 2, col = "blue", lty = 2)
    abline(v = true_re[i, d_re], col = "red", lty = 3, lwd = 1.5)
    abline(v = post$mode[d_re], col = "darkgreen", lty = 3, lwd = 1.5)
    if (d_re == 1) {
      legend("topright", legend = c("Posterior", "Laplace", "True", "Mode"),
             col = c("black", "blue", "red", "darkgreen"),
             lty = c(1, 2, 3, 3), lwd = c(2, 2, 1.5, 1.5), cex = 0.7)
    }
  }
}

# Page 2: 2D contour for m0 vs v0
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (i in sample_ids) {
  d <- data_list[[i]]
  post <- JointODE:::.compute_posterior_laplace(d, true_re[i, ], parameters)

  sd1 <- sqrt(post$cov[1, 1])
  sd2 <- sqrt(post$cov[2, 2])
  g1 <- post$mode[1] + sd1 * seq(-3, 3, length.out = 41)
  g2 <- post$mode[2] + sd2 * seq(-3, 3, length.out = 41)

  lp_2d <- outer(g1, g2, Vectorize(function(x1, x2) {
    b <- post$mode
    b[1] <- x1; b[2] <- x2
    as.numeric(JointODE:::.compute_joint_logpost(
      random_effect = b, data = d,
      parameters = parameters, gradient = FALSE, hessian = FALSE
    ))
  }))

  lp_max <- max(lp_2d, na.rm = TRUE)
  levels <- lp_max - c(0.5, 2, 4.5, 8, 12.5)
  contour(g1, g2, lp_2d, levels = levels, xlab = "init_biomarker", ylab = "init_velocity",
          main = sprintf("Subject %d: init_bio vs init_vel", i),
          labcex = 0.7)
  points(post$mode[1], post$mode[2], pch = 3, col = "darkgreen", cex = 1.5, lwd = 2)
  points(true_re[i, 1], true_re[i, 2], pch = 4, col = "red", cex = 1.5, lwd = 2)
  legend("topright", legend = c("Mode", "True"),
         pch = c(3, 4), col = c("darkgreen", "red"), cex = 0.8)
}

# Page 3: 2D contour for dyn_biomarker vs dyn_velocity
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (i in sample_ids) {
  d <- data_list[[i]]
  post <- JointODE:::.compute_posterior_laplace(d, true_re[i, ], parameters)

  sd3 <- sqrt(post$cov[3, 3])
  sd4 <- sqrt(post$cov[4, 4])
  g3 <- post$mode[3] + sd3 * seq(-3, 3, length.out = 41)
  g4 <- post$mode[4] + sd4 * seq(-3, 3, length.out = 41)

  lp_2d <- outer(g3, g4, Vectorize(function(x3, x4) {
    b <- post$mode
    b[3] <- x3; b[4] <- x4
    as.numeric(JointODE:::.compute_joint_logpost(
      random_effect = b, data = d,
      parameters = parameters, gradient = FALSE, hessian = FALSE
    ))
  }))

  lp_max <- max(lp_2d, na.rm = TRUE)
  levels <- lp_max - c(0.5, 2, 4.5, 8, 12.5)
  contour(g3, g4, lp_2d, levels = levels, xlab = "biomarker", ylab = "velocity",
          main = sprintf("Subject %d: biomarker vs velocity", i),
          labcex = 0.7)
  points(post$mode[3], post$mode[4], pch = 3, col = "darkgreen", cex = 1.5, lwd = 2)
  points(true_re[i, 3], true_re[i, 4], pch = 4, col = "red", cex = 1.5, lwd = 2)
  legend("topright", legend = c("Mode", "True"),
         pch = c(3, 4), col = c("darkgreen", "red"), cex = 0.8)
}

dev.off()
cat("Saved to scripts/posterior_profiles.pdf\n")
