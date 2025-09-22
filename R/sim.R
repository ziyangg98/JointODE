#' Simulated Dataset for Joint Ordinary Differential Equation Modeling
#'
#' @description
#' A pre-computed dataset comprising longitudinal biomarker trajectories and
#' time-to-event outcomes for 200 simulated subjects. This dataset was
#' generated using the internal function \code{.create_example_data()} which
#' calls \code{\link{simulate}} with default parameters and seed = 123.
#' The biomarker trajectories follow a damped harmonic oscillator model
#' with subject-specific initial conditions and external forcing.
#'
#' @format A list with two components:
#' \describe{
#'   \item{\code{data}}{A list containing the actual simulated data with:
#'     \describe{
#'       \item{\code{longitudinal_data}}{A data frame containing longitudinal
#'         observations with columns: id (subject identifier),
#'         time (measurement time), observed (biomarker with error),
#'         biomarker (true value), velocity, acceleration, x1, x2}
#'       \item{\code{survival_data}}{A data frame containing survival
#'         information for 200 subjects with columns: id, time (event/censor),
#'         status (1=event, 0=censored), w1, w2, b (random effect)}
#'       \item{\code{state}}{A 200 × 2 matrix containing initial conditions
#'         \eqn{[m_i(0), \dot{m}_i(0)]} for each subject}
#'     }
#'   }
#'   \item{\code{init}}{Initial parameter values for model fitting containing
#'     coefficients and configurations}
#' }
#'
#' @details
#' \subsection{Data Generation Process}{
#' Generated using \code{.create_example_data(n_subjects = 200, seed = 123)},
#' with follow-up period from 0 to 10 time units and shared random effect
#' \eqn{b_i \sim \mathcal{N}(0, 0.01)} (sd = 0.1).
#' }
#'
#' \subsection{Model Specification}{
#' \strong{Longitudinal Sub-model (Damped Harmonic Oscillator):}
#'
#' The biomarker trajectory follows:
#' \deqn{\ddot{m}_i(t) + 2\xi\omega\dot{m}_i(t) + \omega^2 m_i(t) =
#'       k\omega^2[f_0 + \mathbf{X}_i^T\boldsymbol{\beta}_{exc}]}
#'
#' With parameters:
#' \itemize{
#'   \item Damping ratio: \eqn{\xi = 0.707} (slightly underdamped)
#'   \item Natural period: \eqn{T = 5}, giving \eqn{\omega = 2\pi/5}
#'   \item Excitation amplitude: \eqn{k = 1.0}
#'   \item Excitation offset: \eqn{f_0 = 0.0}
#'   \item Excitation coefficients:
#'     \eqn{\boldsymbol{\beta}_{exc} = (0.8, -0.5)^T}
#' }
#'
#' Initial conditions:
#' \deqn{m_i(0) = -3.0 + \mathbf{X}_i^T\boldsymbol{\beta}_{init} + 0.2 b_i,
#'       \quad \dot{m}_i(0) = 0}
#' where \eqn{\boldsymbol{\beta}_{init} = (0.1, -0.1)^T} and \eqn{b_i} is the
#' random effect.
#'
#' Observed measurements:
#' \deqn{y_{ij} = m_i(t_{ij}) + b_i + \epsilon_{ij},
#'       \quad \epsilon_{ij} \sim \mathcal{N}(0, 0.01)}
#'
#' \strong{Survival Sub-model:}
#' \deqn{\lambda_i(t) = \frac{1}{8}
#'       \exp(0.3m_i(t) + 1.0\dot{m}_i(t) +
#'       \mathbf{W}_i^T\boldsymbol{\phi} + b_i)}
#'
#' where \eqn{\boldsymbol{\phi} = (0.4, -0.6)^T} for covariates w1 and w2.
#' Note: With shape parameter = 1.0, the Weibull baseline reduces to an
#' exponential distribution with constant hazard 1/8.
#' }
#'
#' @source
#' Generated via: \code{sim <- .create_example_data()}
#'
#' @examples
#' # Load the simulated dataset
#' data(sim)
#'
#' # Explore the data structure
#' str(sim$data, max.level = 2)
#'
#' # Summary statistics for longitudinal data
#' summary(sim$data$longitudinal_data[, c("biomarker", "velocity")])
#'
#' # Check event rate
#' mean(sim$data$survival_data$status)
#'
#' # Visualize biomarker trajectories for first 10 subjects
#' library(ggplot2)
#' library(dplyr)
#'
#' sim$data$longitudinal_data %>%
#'   filter(id <= 10) %>%
#'   ggplot(aes(x = time, y = biomarker, group = id)) +
#'   geom_line(alpha = 0.7) +
#'   geom_point(aes(y = observed), size = 0.5, alpha = 0.5) +
#'   facet_wrap(~id, scales = "free_y") +
#'   theme_minimal() +
#'   labs(title = "Biomarker Trajectories (True and Observed)",
#'        x = "Time", y = "Biomarker Value")
#'
#' # Survival analysis
#' library(survival)
#'
#' # Kaplan-Meier curve
#' km_fit <- survfit(Surv(time, status) ~ 1, data = sim$data$survival_data)
#' plot(km_fit, main = "Kaplan-Meier Survival Curve",
#'      xlab = "Time", ylab = "Survival Probability")
#'
#' # Examine relationship between initial biomarker and survival
#' # Note: Ensure proper alignment of subject IDs when extracting initial states
#' sim$data$survival_data$initial_biomarker <- sim$data$state[
#'   sim$data$survival_data$id, "biomarker"
#' ]
#' cox_fit <- coxph(Surv(time, status) ~ w1 + w2, data = sim$data$survival_data)
#' summary(cox_fit)
#'
#' # Verify ODE dynamics - compare numerical solution with data
#' # Extract parameters used in simulation
#' omega <- 2 * pi / 5  # Natural frequency
#' xi <- 0.707           # Damping ratio
#'
#' # Check phase space (velocity vs biomarker)
#' subject_5 <- sim$data$longitudinal_data[
#'   sim$data$longitudinal_data$id == 5,
#' ]
#' plot(subject_5$biomarker, subject_5$velocity,
#'      type = "l", main = "Phase Portrait (Subject 5)",
#'      xlab = "Biomarker", ylab = "Velocity")
#' points(subject_5$biomarker[1], subject_5$velocity[1],
#'        col = "red", pch = 16, cex = 1.5)
#' legend("topright", c("Trajectory", "Initial State"),
#'        col = c("black", "red"), pch = c(NA, 16), lty = c(1, NA))
#'
#' @seealso
#' \code{\link{simulate}} for generating customized datasets with different
#' parameter specifications (note: \code{sim} uses
#' \code{.create_example_data()} internally);
#' \code{\link{JointODE}} for fitting joint ODE models to longitudinal and
#' survival data
#'
#' @concept data-simulation
#' @keywords datasets
"sim"
