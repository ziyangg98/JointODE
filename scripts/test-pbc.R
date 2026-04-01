devtools::load_all()
library(JMbayes2)
library(dplyr)
library(survival)

survival_data <- pbc2.id |>
  transmute(
    id = as.integer(id), time = years, status = status2,
    drug = as.numeric(drug == "D-penicil"), age = age,
    sex = as.numeric(sex == "female")
  )

longitudinal_data <- pbc2 |>
  transmute(id = as.integer(id), time = year, observed = log(albumin)) |>
  left_join(survival_data |> select(id, drug, age, sex), by = "id")

m <- mean(longitudinal_data$observed)
s <- sd(longitudinal_data$observed)
longitudinal_data$observed <- (longitudinal_data$observed - m) / s

t0 <- proc.time()
fit <- JointODE(
  observed ~ biomarker + velocity + drug + (biomarker + velocity | id),
  Surv(time, status) ~ drug + age + sex,
  longitudinal_data, survival_data,
  init = "marginal",
  spline_baseline = list(degree = 1, n_knots = 1),
  control = list(parallel = TRUE, verbose = 3)
)
cat(sprintf("\nElapsed: %.1f s\n", (proc.time() - t0)["elapsed"]))
summary(fit)
