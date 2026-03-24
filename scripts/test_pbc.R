library(JointODE)
library(JMbayes2)
library(survival)
library(dplyr)

cat("=== PBC Dataset Test ===\n\n")

# --- Data Preparation ---
survival_data <- pbc2.id |>
  transmute(
    id = as.integer(id), time = years, status = status2,
    drug = as.numeric(drug == "D-penicil"), age = age,
    sex = as.numeric(sex == "female")
  )

longitudinal_data <- pbc2 |>
  transmute(
    id = as.integer(id), time = year,
    observed = log(serBilir)
  ) |>
  left_join(
    survival_data |> select(id, drug, age, sex),
    by = "id"
  )

bili_mean <- mean(longitudinal_data$observed)
bili_sd <- sd(longitudinal_data$observed)
longitudinal_data$observed <-
  (longitudinal_data$observed - bili_mean) / bili_sd

cat(sprintf(
  "Patients: %d | Events: %d (%.0f%%) | Obs: %d\n",
  nrow(survival_data), sum(survival_data$status),
  100 * mean(survival_data$status),
  nrow(longitudinal_data)
))
cat(sprintf(
  "Log-bilirubin: mean=%.3f, sd=%.3f (standardized)\n",
  bili_mean, bili_sd
))

# --- MarginalODE ---
cat("\n--- MarginalODE ---\n")
t0 <- proc.time()
fit_marginal <- MarginalODE(
  observed ~ drug,
  longitudinal_data,
  time = "time", id = "id",
  control = list(verbose = 2, parallel = TRUE)
)
cat(sprintf("Elapsed: %.1f s\n", (proc.time() - t0)["elapsed"]))
summary(fit_marginal)

# --- JointODE ---
cat("\n--- JointODE ---\n")
t0 <- proc.time()
fit_ode <- JointODE(
  longitudinal_formula =
    observed ~ biomarker + velocity + drug +
      (biomarker + velocity | id),
  survival_formula = Surv(time, status) ~ drug + age + sex,
  longitudinal_data = longitudinal_data,
  survival_data = survival_data,
  init = "marginal",
  control = list(verbose = 2, maxit = 50, parallel = TRUE)
)
elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("\nTotal elapsed: %.1f s\n", elapsed))

cat("\n--- Summary ---\n")
summary(fit_ode)

# --- Predict ---
cat("\n--- Predict (first 3 subjects) ---\n")
pred_ode <- predict(fit_ode)
ids3 <- unique(pred_ode$id)[1:3]
for (sid in ids3) {
  p <- pred_ode[pred_ode$id == sid, ]
  cat(sprintf(
    "Subject %s: %d pts, biomarker [%.3f, %.3f], surv [%.3f, %.3f]\n",
    sid, nrow(p),
    min(p$biomarker), max(p$biomarker),
    min(p$survival), max(p$survival)
  ))
}

# --- Extended Cox for comparison ---
cat("\n--- Longitudinal SSE ---\n")
pred_ode <- predict(fit_ode)
merged <- merge(pred_ode, longitudinal_data,
  by.x = c("id", "time"), by.y = c("id", "time"))
sse <- sum((merged$biomarker - merged$observed)^2)
cat(sprintf("Total SSE: %.4f  (n=%d, MSE=%.6f, RMSE=%.4f)\n",
  sse, nrow(merged), sse / nrow(merged), sqrt(sse / nrow(merged))))
subj_sse <- tapply(
  (merged$biomarker - merged$observed)^2, merged$id, sum
)
cat(sprintf("Per-subject SSE: median=%.4f, mean=%.4f, max=%.4f\n",
  median(subj_sse), mean(subj_sse), max(subj_sse)))

cat("\n--- Save ---\n")
saveRDS(fit_ode, "tests/fit_pbc.rds")
cat("Saved to tests/fit_pbc.rds\n")

cat("\n--- Extended Cox ---\n")
td <- longitudinal_data |>
  arrange(id, time) |>
  left_join(
    survival_data |> select(id, time_event = time, status),
    by = "id"
  ) |>
  group_by(id) |>
  mutate(
    tstart = time,
    tstop = lead(time, default = first(time_event)),
    event = ifelse(row_number() == n(), first(status), 0L)
  ) |>
  ungroup() |>
  filter(tstop > tstart)

cox_td <- coxph(
  Surv(tstart, tstop, event) ~ observed + drug + age + sex,
  data = td, x = TRUE
)
cat("Extended Cox coefficients:\n")
print(summary(cox_td)$coefficients)

# --- Hazard Comparison ---
cat("\n--- Hazard Coefficient Comparison ---\n")
ode_c <- coef(fit_ode)
ode_se <- sqrt(diag(vcov(fit_ode)))
cox_s <- summary(cox_td)$coefficients

f <- function(est, se) sprintf("%.3f (%.3f)", est, se)
cat(sprintf("  %-12s %20s %20s\n", "Parameter", "JointODE", "Ext.Cox"))
cat(sprintf("  %-12s %20s %20s\n", "value",
  f(ode_c["hazard:alpha_1"], ode_se["hazard:alpha_1"]),
  f(cox_s["observed", 1], cox_s["observed", 3])))
cat(sprintf("  %-12s %20s %20s\n", "slope",
  f(ode_c["hazard:alpha_2"], ode_se["hazard:alpha_2"]),
  "--"))
cat(sprintf("  %-12s %20s %20s\n", "drug",
  f(ode_c["hazard:drug"], ode_se["hazard:drug"]),
  f(cox_s["drug", 1], cox_s["drug", 3])))
cat(sprintf("  %-12s %20s %20s\n", "age",
  f(ode_c["hazard:age"], ode_se["hazard:age"]),
  f(cox_s["age", 1], cox_s["age", 3])))
cat(sprintf("  %-12s %20s %20s\n", "sex",
  f(ode_c["hazard:sex"], ode_se["hazard:sex"]),
  f(cox_s["sex", 1], cox_s["sex", 3])))

cat("\nDone.\n")
