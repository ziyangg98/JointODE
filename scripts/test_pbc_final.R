library(JointODE)
library(JMbayes2)
library(dplyr)

survival_data <- pbc2.id |>
  transmute(
    id = as.integer(id), time = years, status = status2,
    drug = as.numeric(drug == "D-penicil"),
    age = scale(age)[, 1],
    sex = as.numeric(sex == "female")
  )

longitudinal_data <- pbc2 |>
  transmute(id = as.integer(id), time = year, observed = log(albumin)) |>
  left_join(survival_data |> select(id, drug, age, sex), by = "id")

alb_mean <- mean(longitudinal_data$observed)
alb_sd <- sd(longitudinal_data$observed)
longitudinal_data$observed <- (longitudinal_data$observed - alb_mean) / alb_sd

# Level 2b: random damping
cat("=== Level 2b ===\n")
t0 <- proc.time()
mfit2b <- MarginalODE(
  observed ~ biomarker + velocity + drug + (velocity | id),
  longitudinal_data,
  control = list(parallel = TRUE, verbose = 1)
)
cat(sprintf("Elapsed: %.1f sec\n", (proc.time() - t0)[3]))
print(summary(mfit2b))

# Gradient check
obj <- mfit2b$tmb_obj
par_best <- obj$env$last.par.best
pn <- names(par_best)
outer_par <- par_best[pn != "random_effects"]
grad <- obj$gr(outer_par)
cat("\n=== Gradient check ===\n")
for (k in seq_along(outer_par)) {
  eps <- max(abs(outer_par[k]) * 1e-4, 1e-6)
  p_plus <- outer_par; p_plus[k] <- p_plus[k] + eps
  p_minus <- outer_par; p_minus[k] <- p_minus[k] - eps
  nll_p <- tryCatch(as.numeric(obj$fn(p_plus)), error = function(e) NA)
  nll_m <- tryCatch(as.numeric(obj$fn(p_minus)), error = function(e) NA)
  fd <- if (!is.na(nll_p) && !is.na(nll_m)) (nll_p - nll_m) / (2 * eps) else NA
  ratio <- if (!is.na(fd) && abs(fd) > 1e-10) grad[k] / fd else NaN
  cat(sprintf("  par[%2d]: AD=%11.4e FD=%11.4e ratio=%.4f\n",
              k, grad[k], if (is.na(fd)) NaN else fd, ratio))
}

# Level 2c
cat("\n=== Level 2c ===\n")
t0 <- proc.time()
mfit2c <- MarginalODE(
  observed ~ biomarker + velocity + drug + (biomarker + velocity | id),
  longitudinal_data,
  control = list(parallel = TRUE, verbose = 1)
)
cat(sprintf("Elapsed: %.1f sec\n", (proc.time() - t0)[3]))
print(summary(mfit2c))
