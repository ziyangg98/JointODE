devtools::load_all()
data("sim", package = "JointODE")

fit <- JointODE(
  longitudinal_formula = observed ~ biomarker + velocity + x1 + x2 +
    (biomarker + velocity | id),
  survival_formula = Surv(time, status) ~ w1 + w2,
  longitudinal_data = sim$data$longitudinal_data[
    , c("id", "time", "observed", "x1", "x2")
  ],
  survival_data = sim$data$survival_data,
  init = "marginal",
  control = list(verbose = 3, parallel = TRUE)
)
