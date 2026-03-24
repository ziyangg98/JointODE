library(JointODE)
library(dplyr)
library(JMbayes2)

survival_data <- pbc2.id |> transmute(
  id = as.integer(id), time = years, status = status2,
  drug = as.numeric(drug == "D-penicil"), age = age,
  sex = as.numeric(sex == "female")
)

biomarkers <- c("serBilir", "albumin", "SGOT", "platelets", "alkaline")
results <- list()

for (bm in biomarkers) {
  cat(sprintf("\n========== %s ==========\n", bm))
  ld <- pbc2 |>
    transmute(
      id = as.integer(id), time = year,
      observed = log(.data[[bm]])
    ) |>
    left_join(survival_data |> select(id, drug, age, sex), by = "id") |>
    filter(is.finite(observed))
  m <- mean(ld$observed)
  s <- sd(ld$observed)
  ld$observed <- (ld$observed - m) / s

  t0 <- proc.time()
  fit <- tryCatch(
    JointODE(
      longitudinal_formula =
        observed ~ biomarker + velocity + drug +
        (biomarker + velocity | id),
      survival_formula = Surv(time, status) ~ drug + age + sex,
      longitudinal_data = ld,
      survival_data = survival_data,
      init = "marginal",
      control = list(parallel = TRUE, verbose = 1, maxit = 200)
    ),
    error = function(e) {
      cat(sprintf("FAILED: %s\n", e$message))
      NULL
    }
  )
  el <- (proc.time() - t0)["elapsed"]

  if (!is.null(fit)) {
    cat(sprintf("\nElapsed: %.1f s\n", el))
    print(summary(fit))
    results[[bm]] <- fit
  }
}

saveRDS(results, "pbc_results.rds")

cat("\n========== SUMMARY ==========\n")
for (bm in names(results)) {
  f <- results[[bm]]
  cat(sprintf(
    "%s: %d iter, LogLik=%.2f, C=%.3f\n",
    bm, f$convergence$iterations, f$logLik, f$cindex
  ))
}
