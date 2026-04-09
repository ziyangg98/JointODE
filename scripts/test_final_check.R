library(JointODE)

# --- sim data: dynamics RE ---
cat("=== sim: Level 2c ===\n")
data(sim)
ld <- sim$data$longitudinal_data[, c("id", "time", "observed", "x1", "x2")]
fit_sim <- MarginalODE(
  observed ~ biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
  ld, control = list(verbose = 0))
cat("converged:", fit_sim$convergence$converged, "\n")
obj <- fit_sim$tmb_obj; par <- obj$env$last.par.best; pn <- names(par)
outer <- par[pn != "random_effects"]; grad <- obj$gr(outer)
max_ratio_err <- 0
for (k in seq_along(outer)) {
  eps <- max(abs(outer[k]) * 1e-4, 1e-6)
  p <- outer; p[k] <- p[k] + eps; m <- outer; m[k] <- m[k] - eps
  fd <- (obj$fn(p) - obj$fn(m)) / (2 * eps)
  if (abs(fd) > 1e-10) max_ratio_err <- max(max_ratio_err, abs(grad[k]/fd - 1))
}
cat(sprintf("max |AD/FD - 1| = %.2e\n", max_ratio_err))

# --- PBC: Level 2b ---
cat("\n=== PBC: Level 2b ===\n")
library(JMbayes2); library(dplyr)
sd <- pbc2.id |> transmute(id = as.integer(id), time = years, status = status2,
  drug = as.numeric(drug == "D-penicil"))
ld2 <- pbc2 |> transmute(id = as.integer(id), time = year,
  observed = scale(log(albumin))[,1]) |>
  left_join(sd |> select(id, drug), by = "id")
fit_pbc <- MarginalODE(
  observed ~ biomarker + velocity + drug + (velocity | id),
  ld2, control = list(verbose = 1))
cat("converged:", fit_pbc$convergence$converged, "\n")
cat("coefs:", coef(fit_pbc)[1:2], "\n")
