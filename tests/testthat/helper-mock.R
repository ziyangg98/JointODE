# ==============================================================================
# Mock Object Factories for S3 Method Tests
# ==============================================================================

# Build processed data from sim dataset
.build_mock_data <- function() {
  long_formula <- observed ~ biomarker + velocity + x1 + x2 +
    (biomarker + velocity | id)
  surv_formula <- Surv(time, status) ~ w1 + w2

  .process(
    longitudinal_data = sim$data$longitudinal_data,
    survival_data = sim$data$survival_data,
    state = sim$data$state,
    survival_formula = surv_formula,
    longitudinal_formula = long_formula
  )
}

create_mock_jointode <- function(n_subjects = 20L) {
  data_list <- .build_mock_data()
  # Subset to n_subjects for faster tests (plot, S3 methods)
  n_subjects <- min(n_subjects, length(data_list))
  data_list <- data_list[seq_len(n_subjects)]

  parameters <- sim$init
  # Add names to coefficients
  parameters$coefficients$baseline <- setNames(
    parameters$coefficients$baseline,
    paste0("bs", seq_along(parameters$coefficients$baseline))
  )
  parameters$coefficients$hazard <- setNames(
    parameters$coefficients$hazard,
    c("value", "slope", "w1", "w2")
  )
  # longitudinal already named from sim$init

  n_params <- .count_params(parameters)
  n_re <- nrow(parameters$coefficients$random_effect_sigma)
  random_effects <- sim$data$random_effects[seq_len(n_subjects), , drop = FALSE]

  coef_names_expanded <- c(
    paste0("baseline:", names(parameters$coefficients$baseline)),
    paste0("hazard:", names(parameters$coefficients$hazard)),
    paste0(
      "longitudinal:",
      names(parameters$coefficients$longitudinal)
    )
  )
  n_coef <- length(coef_names_expanded)
  vcov_matrix <- diag(0.01, n_coef)
  dimnames(vcov_matrix) <- list(coef_names_expanded, coef_names_expanded)

  structure(
    list(
      parameters = parameters,
      logLik = -500.0,
      AIC = -500.0 * -2 + 2 * n_params,
      BIC = -500.0 * -2 + n_params * log(n_subjects),
      cindex = 0.65,
      convergence = list(
        converged = TRUE,
        em_iterations = 10,
        message = "EM algorithm converged after 10 iterations"
      ),
      random_effects = list(
        estimates = random_effects,
        variances = lapply(
          seq_len(n_subjects),
          function(i) diag(0.01, n_re)
        )
      ),
      vcov = vcov_matrix,
      data = data_list,
      control = JointODE.control(),
      call = quote(JointODE(
        longitudinal_formula = observed ~ biomarker + velocity +
          x1 + x2 + (biomarker + velocity | id),
        survival_formula = Surv(time, status) ~ w1 + w2,
        longitudinal_data = sim$data$longitudinal_data,
        survival_data = sim$data$survival_data
      ))
    ),
    class = "JointODE"
  )
}
