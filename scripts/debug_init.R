devtools::load_all()
data("sim", package = "JointODE")

long_data <- sim$data$longitudinal_data[, c("id", "time", "observed", "x1", "x2")]
surv_data <- sim$data$survival_data

parsed_long <- .parse_longitudinal_formula(
  observed ~ biomarker + velocity + x1 + x2 + (biomarker + velocity | id)
)
parsed_surv <- .parse_survival_formula(Surv(time, status) ~ w1 + w2)

model_config <- .setup_model(
  long_data, surv_data, Surv(time, status) ~ w1 + w2,
  gamma = 1, parsed_long, parsed_surv,
  spline_config = list(degree = 2, n_knots = 1, knot_placement = "equal",
                       boundary_knots = NULL)
)

cat("=== Default init ===\n")
default_params <- .default_parameters(
  model_config$dims, 1, parsed_long, model_config$spline_baseline_config
)
cat("baseline:", round(default_params$coefficients$baseline, 4), "\n")
cat("hazard:", round(default_params$coefficients$hazard, 4), "\n")
cat("longitudinal:", round(default_params$coefficients$longitudinal, 4), "\n")
cat("initial_state:", round(default_params$coefficients$initial_state, 4), "\n")
cat("sigma_e:", default_params$coefficients$measurement_error_sd, "\n")
cat("Sigma_b diag:", round(diag(default_params$coefficients$random_effect_sigma), 4), "\n")

cat("\n=== Marginal init ===\n")
marginal_params <- .initialize_from_marginal(
  long_data, surv_data, 1,
  list(verbose = 1, parallel = FALSE, n_cores = 0),
  parsed_long, parsed_surv, model_config
)
cat("baseline:", round(marginal_params$coefficients$baseline, 4), "\n")
cat("hazard:", round(marginal_params$coefficients$hazard, 4), "\n")
cat("longitudinal:", round(marginal_params$coefficients$longitudinal, 4), "\n")
cat("initial_state:", round(marginal_params$coefficients$initial_state, 4), "\n")
cat("sigma_e:", round(marginal_params$coefficients$measurement_error_sd, 4), "\n")
cat("Sigma_b diag:", round(diag(marginal_params$coefficients$random_effect_sigma), 4), "\n")

if (!is.null(marginal_params$random_effects_init)) {
  cat("\nRE init (first 5):\n")
  print(round(marginal_params$random_effects_init[1:5, ], 4))
} else {
  cat("\nRE init: NULL\n")
}

cat("\n=== True values (from sim$init) ===\n")
true <- sim$init$coefficients
cat("baseline:", round(true$baseline, 4), "\n")
cat("hazard:", round(true$hazard, 4), "\n")
cat("longitudinal:", round(true$longitudinal, 4), "\n")
cat("initial_state:", round(true$initial_state, 4), "\n")
cat("sigma_e:", round(true$measurement_error_sd, 4), "\n")
cat("Sigma_b diag:", round(diag(true$random_effect_sigma), 4), "\n")
