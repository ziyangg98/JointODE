library(JointODE)
library(ggplot2)
library(patchwork)

set.seed(1)

data_path <- "analysis/data/num_bili.csv"
trim_value <- 0
prefix <- "analysis/num_bili_x1x2_random_biomarker_velocity_trim0"
fit_path <- paste0(prefix, "_fit.rds")
summary_path <- paste0(prefix, "_summary.csv")
prediction_path <- paste0(prefix, "_predictions.csv")
diagnostic_plot_path <- paste0(prefix, "_prediction_diagnostics.png")
trajectory_plot_path <- paste0(prefix, "_prediction_trajectories_top12.png")
log_path <- paste0(prefix, ".log")

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

dat <- read.csv(data_path)
n_subjects <- length(unique(dat$id))
cat("data_path:", data_path, "\n")
cat("n_observations:", nrow(dat), "\n")
cat("n_subjects:", n_subjects, "\n")
cat("trim:", trim_value, "\n")
cat("formula: observed ~ biomarker + velocity + x1 + x2 +",
    "(1 + biomarker + velocity | id)\n\n")

t0 <- proc.time()
fit <- MarginalODE(
  observed ~ biomarker + velocity + x1 + x2 +
    (1 + biomarker + velocity | id),
  data = dat,
  control = MarginalODE.control(
    maxit = 200,
    outer.maxit = 100,
    trim = trim_value,
    verbose = 2,
    parallel = TRUE,
    n_cores = 4
  )
)
elapsed <- (proc.time() - t0)[["elapsed"]]
saveRDS(fit, fit_path)

obs <- dat[order(dat$id, dat$time, seq_len(nrow(dat))), ]
pred <- predict(fit)
pred <- pred[order(pred$id, pred$time, seq_len(nrow(pred))), ]
stopifnot(nrow(obs) == nrow(pred))

out <- data.frame(
  id = obs$id,
  time = obs$time,
  observed = obs$observed,
  x1 = obs$x1,
  x2 = obs$x2,
  fitted = pred$biomarker,
  velocity = pred$velocity
)
out$residual <- out$observed - out$fitted
write.csv(out, prediction_path, row.names = FALSE)

residual <- out$residual
fit_summary <- data.frame(
  model = "MarginalODE",
  trim = trim_value,
  elapsed_seconds = elapsed,
  n_observations = nrow(dat),
  n_subjects = n_subjects,
  converged = fit$convergence$converged,
  convergence_message = fit$convergence$message,
  logLik = fit$logLik,
  AIC = fit$AIC,
  BIC = fit$BIC,
  sigma_e = fit$measurement_error_sd,
  RMSE = sqrt(mean(residual^2)),
  MAE = mean(abs(residual)),
  bias = mean(residual),
  stringsAsFactors = FALSE
)
write.csv(fit_summary, summary_path, row.names = FALSE)

p_fit <- ggplot(out, aes(fitted, observed)) +
  geom_point(alpha = 0.18, size = 0.45, color = "#2f5f8f") +
  geom_abline(slope = 1, intercept = 0, color = "#b23a48", linewidth = 0.6) +
  coord_equal() +
  theme_bw(base_size = 11) +
  labs(
    title = "Observed vs fitted",
    subtitle = sprintf(
      "RMSE %.4f, MAE %.4f",
      fit_summary$RMSE,
      fit_summary$MAE
    ),
    x = "Fitted biomarker",
    y = "Observed biomarker"
  )

p_res_fit <- ggplot(out, aes(fitted, residual)) +
  geom_point(alpha = 0.15, size = 0.45, color = "#326b5f") +
  geom_hline(yintercept = 0, color = "#333333", linewidth = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "#b23a48", linewidth = 0.7) +
  theme_bw(base_size = 11) +
  labs(title = "Residual vs fitted", x = "Fitted biomarker", y = "Residual")

p_res_time <- ggplot(out, aes(time, residual)) +
  geom_point(alpha = 0.12, size = 0.4, color = "#6f5b9a") +
  geom_hline(yintercept = 0, color = "#333333", linewidth = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "#b23a48", linewidth = 0.7) +
  theme_bw(base_size = 11) +
  labs(title = "Residual vs time", x = "Time", y = "Residual")

p_hist <- ggplot(out, aes(residual)) +
  geom_histogram(bins = 80, fill = "#d8a24a", color = "white", linewidth = 0.1) +
  geom_vline(xintercept = 0, color = "#333333", linewidth = 0.4) +
  theme_bw(base_size = 11) +
  labs(title = "Residual distribution", x = "Residual", y = "Count")

ggsave(
  diagnostic_plot_path,
  (p_fit | p_res_fit) / (p_res_time | p_hist),
  width = 12,
  height = 8,
  dpi = 180
)

top_ids <- head(names(sort(table(out$id), decreasing = TRUE)), 12)
traj <- out[out$id %in% top_ids, ]
traj$id <- factor(traj$id, levels = top_ids)
p_traj <- ggplot(traj, aes(time)) +
  geom_point(aes(y = observed), alpha = 0.5, size = 0.75, color = "#2f5f8f") +
  geom_line(aes(y = fitted), color = "#b23a48", linewidth = 0.55) +
  facet_wrap(~id, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 10) +
  labs(title = "Top 12 subjects by observation count", x = "Time", y = "Biomarker")

ggsave(trajectory_plot_path, p_traj, width = 12, height = 8, dpi = 180)

cat("\nfit_path:", fit_path, "\n")
cat("summary_path:", summary_path, "\n")
cat("prediction_path:", prediction_path, "\n")
cat("diagnostic_plot_path:", diagnostic_plot_path, "\n")
cat("trajectory_plot_path:", trajectory_plot_path, "\n")
cat("log_path:", log_path, "\n\n")

print(fit_summary)
cat("\nResidual quantiles:\n")
print(quantile(residual, c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1)))
cat("\n")
print(fit)
print(summary(fit))
