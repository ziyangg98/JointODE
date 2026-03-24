# State Optimization (Joint + Marginal)

# Negate C++ loglik result for nlm minimization
.negate_result <- function(result) {
  value <- -as.numeric(result)
  attr(value, "gradient") <- -as.vector(attr(result, "gradient"))
  attr(value, "hessian") <- -as.matrix(attr(result, "hessian"))
  value
}

#' @noRd
.estimate_joint_state <- function(
  initial_guess, data, random_effect, parameters
) {
  nlm(
    f = function(state) {
      .negate_result(.compute_joint_state(
        state, data, random_effect, parameters,
        gradient = TRUE, hessian = TRUE
      ))
    },
    p = initial_guess, hessian = FALSE,
    check.analyticals = FALSE
  )$estimate
}

#' @noRd
.estimate_marginal_state <- function(
  initial_guess, subject_data, parameters, biomarker_clamp
) {
  nlm(
    f = function(state) {
      .negate_result(.compute_marginal_state(
        state, subject_data, parameters,
        biomarker_clamp = biomarker_clamp,
        gradient = TRUE, hessian = TRUE
      ))
    },
    p = initial_guess, hessian = FALSE,
    check.analyticals = FALSE
  )$estimate
}
