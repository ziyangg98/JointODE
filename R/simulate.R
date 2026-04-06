#' Simulate Data from a Joint Ordinary Differential Equation Model
#'
#' @description
#' Generates synthetic longitudinal and time-to-event data under a joint
#' modeling framework where the longitudinal biomarker trajectory follows a
#' damped harmonic oscillator model (second-order ODE), and the survival
#' process is associated with the biomarker dynamics through shared random
#' effects and trajectory-dependent hazard functions. The biomarker dynamics
#' are parameterized using physically interpretable parameters: damping ratio,
#' natural period, and excitation amplitude.
#'
#' @param n_subjects Integer specifying the number of subjects to simulate
#'   (default: 500)
#' @param longitudinal List specifying the longitudinal sub-model parameters:
#'   \describe{
#'     \item{xi}{Damping ratio \eqn{\xi} controlling the oscillation decay.
#'       Specified as c(mean = ..., sd = ...) for population mean and
#'       standard deviation. Values: \eqn{\xi < 1} (underdamped),
#'       \eqn{\xi = 1} (critically damped), \eqn{\xi > 1} (overdamped).
#'       Default: c(mean = 0.4, sd = 0.05)}
#'     \item{period}{Natural period \eqn{T} of oscillation in time units,
#'       related to natural frequency as \eqn{\omega = 2\pi/T}.
#'       Specified as c(mean = ..., sd = ...) for population mean and
#'       standard deviation. Default: c(mean = 6, sd = 1.0)}
#'     \item{excitation}{List specifying external forcing parameters:
#'       \describe{
#'         \item{offset}{Constant excitation term \eqn{f_0} (default: 0.0)}
#'         \item{covariates}{Named vector of covariate effects
#'           \eqn{\boldsymbol{\beta}_{exc}} on excitation
#'           (default: c(x1 = 0.5, x2 = -0.45))}
#'       }
#'     }
#'     \item{initial}{List specifying initial condition parameters:
#'       \describe{
#'         \item{offset}{Named vector with baseline initial values:
#'           \code{c(biomarker = ..., velocity = ...)}
#'           (default: c(biomarker = -0.5, velocity = -0.05))}
#'         \item{covariates}{List with two named vectors for covariate effects:
#'           \describe{
#'             \item{biomarker}{Covariate effects on initial biomarker
#'               (default: c(x1 = 0.08, x2 = -0.06))}
#'             \item{velocity}{Covariate effects on initial velocity
#'               (default: c(x1 = 0.25, x2 = -0.2))}
#'           }}
#'       }
#'     }
#'     \item{error_sd}{Standard deviation \eqn{\sigma_{\epsilon}} of the
#'       measurement error process (default: 0.1)}
#'     \item{n_measurements}{Number of longitudinal measurements per subject
#'       (default: 100)}
#'   }
#' @param survival List specifying the survival sub-model parameters:
#'   \describe{
#'     \item{baseline}{List defining the Weibull baseline hazard function:
#'       \describe{
#'         \item{type}{Character string specifying the baseline hazard type
#'           (currently only "weibull" is supported)}
#'         \item{shape}{Weibull shape parameter \eqn{\kappa > 0}
#'           (default: 2.0)}
#'         \item{scale}{Weibull scale parameter \eqn{\lambda > 0}
#'           (default: 100.0)}
#'       }
#'     }
#'     \item{value}{Association parameter \eqn{\alpha_1} linking current
#'       biomarker value to hazard (default: 0.8)}
#'     \item{slope}{Association parameter \eqn{\alpha_2} linking biomarker
#'       velocity (or its power) to hazard (default: 2.0)}
#'     \item{gamma}{Power parameter for velocity contribution, where
#'       \eqn{\gamma = 0} excludes velocity effect, \eqn{\gamma = 1} uses
#'       linear velocity \eqn{\alpha_2\dot{m}_i(t)}, and \eqn{\gamma = 2}
#'       uses quadratic velocity \eqn{\alpha_2[\dot{m}_i(t)]^2} (default: 1)}
#'     \item{covariates}{Named vector of regression coefficients
#'       \eqn{\boldsymbol{\phi}} for survival covariates
#'       (default: c(w1 = 0.6, w2 = -0.8))}
#'   }
#' @param covariates List defining the distributions of baseline covariates:
#'   \describe{
#'     \item{x1}{List with \code{type = "normal"}, \code{mean = 0},
#'       \code{sd = 1} for standardized continuous covariate (longitudinal)}
#'     \item{x2}{List with \code{type = "normal"}, \code{mean = 0},
#'       \code{sd = 1} for standardized continuous covariate (longitudinal)}
#'     \item{w1}{List with \code{type = "normal"}, \code{mean = 0},
#'       \code{sd = 1} for standardized continuous covariate (survival)}
#'     \item{w2}{List with \code{type = "binary"} and \code{prob = 0.5}
#'       for binary covariate (survival)}
#'   }
#' @param maxt Positive scalar specifying the maximum follow-up time in the
#'   study (default: 10 time units)
#' @param seed Integer seed for random number generation to ensure
#'   reproducibility (default: 42)
#'
#' @return A list containing four components:
#' \describe{
#'   \item{\code{longitudinal_data}}{Data frame comprising longitudinal
#'     observations with columns:
#'     \itemize{
#'       \item \code{id}: Subject identifier (integer)
#'       \item \code{time}: Observation time point (numeric)
#'       \item \code{observed}: Measured biomarker value including
#'         measurement error, \eqn{y_{ij}}
#'       \item \code{biomarker}: True underlying biomarker value,
#'         \eqn{m_i(t_{ij})}
#'       \item \code{velocity}: First derivative of the biomarker trajectory,
#'         \eqn{\dot{m}_i(t_{ij})}
#'       \item \code{acceleration}: Second derivative of the biomarker
#'         trajectory, \eqn{\ddot{m}_i(t_{ij})}
#'       \item \code{x1}, \code{x2}:
#'         Time-invariant covariates (if specified)
#'     }
#'   }
#'   \item{\code{survival_data}}{Data frame containing time-to-event data
#'     with columns:
#'     \itemize{
#'       \item \code{id}: Subject identifier
#'       \item \code{time}: Observed event or censoring time, \eqn{T_i}
#'       \item \code{status}: Event indicator, \eqn{\delta_i}
#'         (1 = event observed, 0 = censored)
#'       \item \code{w1}, \code{w2}:
#'         Baseline survival covariates (if specified)
#'     }
#'   }
#'   \item{\code{state}}{An \eqn{n \times 2} data frame containing
#'     initial states \eqn{[m_i(0), \dot{m}_i(0)]} for each subject}
#'   \item{\code{random_effects}}{An \eqn{n \times 2} matrix containing
#'     subject-specific random effects (centered at zero) for:
#'     \itemize{
#'       \item \code{dyn_value}: Random effect for biomarker level term
#'       \item \code{dyn_slope}: Random effect for velocity term
#'     }
#'     The matrix has attributes \code{mu} (population mean vector of length 2)
#'     and \code{sigma} (\eqn{2 \times 2} covariance matrix). These correspond
#'     to the random effects specified in the formula structure
#'     \code{(biomarker + velocity | id)}. Fixed effects (\code{dyn_const}
#'     and covariate terms) are not included as they have no between-subject
#'     variability.}
#' }
#'
#' @details
#' The simulation framework implements a joint model comprising longitudinal
#' and survival sub-models linked through shared random effects and
#' trajectory-dependent associations.
#'
#' \subsection{Default Parameter Design}{
#' The default parameters are calibrated to achieve the following properties:
#' \itemize{
#'   \item Approximately 60\% censoring rate (40\% event rate)
#'   \item Initial biomarker difference between risk groups: ~0.2
#'   \item Final biomarker difference between event and censored groups: ~0.7
#'   \item Damping ratio centered at 0.4 (underdamped oscillations)
#'   \item Significant velocity differences between patients with different
#'         risk profiles, enabling survival prediction based on trajectory
#'         dynamics
#' }
#' These settings produce realistic heterogeneity in biomarker trajectories
#' while maintaining numerical stability for joint model estimation.
#' }
#'
#' \subsection{Longitudinal Sub-model}{
#' The biomarker trajectory \eqn{m_i(t)} for subject \eqn{i} follows a damped
#' harmonic oscillator with external forcing:
#' \deqn{\ddot{m}_i(t) + 2\xi\omega\dot{m}_i(t) + \omega^2 m_i(t) =
#'   k \omega^2 [f_0 + \mathbf{X}_i^T\boldsymbol{\beta}_{exc}]}
#' where:
#' \itemize{
#'   \item \eqn{\omega = 2\pi/T} is the natural angular frequency
#'   \item \eqn{\xi} is the damping ratio determining oscillation behavior
#'   \item \eqn{k} scales the excitation amplitude
#'   \item \eqn{f_0} is the constant excitation term
#'   \item \eqn{\mathbf{X}_i} contains time-invariant covariates
#'   \item \eqn{\boldsymbol{\beta}_{exc}} represents covariate effects on
#'           excitation
#' }
#'
#' Initial conditions are specified as:
#' \itemize{
#'   \item \eqn{m_i(0) = \mu_{m,0} +
#'     \mathbf{X}_i^T\boldsymbol{\beta}_{m,init}}
#'   \item \eqn{\dot{m}_i(0) = \mu_{v,0} +
#'     \mathbf{X}_i^T\boldsymbol{\beta}_{v,init}}
#' }
#'
#' The observed longitudinal measurements incorporate additive Gaussian noise:
#' \deqn{y_{ij} = m_i(t_{ij}) + b_i + \epsilon_{ij}}
#' where \eqn{\epsilon_{ij} \sim \mathcal{N}(0, \sigma_\epsilon^2)}
#' represents independent measurement error.
#' }
#'
#' \subsection{Survival Sub-model}{
#' The instantaneous hazard function incorporates both current biomarker value
#' and velocity (or its power):
#' \deqn{\lambda_i(t) = \lambda_0(t)\exp(\alpha_1 m_i(t) +
#'   \alpha_2[\dot{m}_i(t)]^\gamma + \mathbf{W}_i^T\boldsymbol{\phi})}
#' where:
#' \itemize{
#'   \item \eqn{\lambda_0(t)} denotes the Weibull baseline hazard:
#'     \eqn{\lambda_0(t) = (\kappa/\lambda)(t/\lambda)^{\kappa-1}}
#'   \item \eqn{\alpha_1} quantifies the association with current biomarker
#'     value
#'   \item \eqn{\alpha_2} quantifies the association with biomarker velocity
#'     (or its power)
#'   \item \eqn{\gamma \in \{0, 1, 2\}} determines the power of velocity:
#'     0 (no velocity effect), 1 (linear), or 2 (quadratic)
#'   \item \eqn{\mathbf{W}_i} contains baseline covariates
#'   \item \eqn{\boldsymbol{\phi}} represents covariate effects on survival
#' }
#' }
#'
#' \subsection{Patient-Specific Dynamics}{
#' The function models patient heterogeneity through continuous distributions
#' of dynamics parameters \eqn{\xi_i} (damping ratio) and \eqn{T_i} (period):
#' \deqn{\xi_i \sim \mathcal{N}(\mu_\xi, \sigma_\xi^2), \quad
#'   T_i \sim \mathcal{N}(\mu_T, \sigma_T^2)}
#' These are transformed to ODE parameters via:
#' \deqn{\omega_i = 2\pi/T_i, \quad
#'   \beta_{1,i} = -\omega_i^2, \quad
#'   \beta_{2,i} = -2\xi_i\omega_i}
#' The transformation uses the Delta method to preserve the correct covariance
#' structure in the ODE parameter space.
#'
#' \strong{Random Effects Structure:} Only \eqn{\beta_{1,i}}
#' (\code{dyn_value}) and \eqn{\beta_{2,i}} (\code{dyn_slope}) have
#' subject-specific random effects. The constant excitation term
#' (\code{dyn_const}) and all covariate excitation coefficients
#' (\code{dyn_x1}, \code{dyn_x2}, etc.) are fixed effects that take
#' the same population mean value for all subjects, with no
#' between-subject variability. This structure corresponds to the
#' formula specification \code{observed ~ biomarker + velocity + x1 +
#' x2 + (biomarker + velocity | id)}.
#' }
#'
#' @examples
#' # Example 1: Simple simulation with default parameters
#' sim_basic <- simulate(n_subjects = 20, seed = 123)
#'
#' # Explore the output structure
#' names(sim_basic)
#' head(sim_basic$longitudinal_data)
#' head(sim_basic$survival_data)
#'
#' # Check patient-specific dynamics
#' # Each patient has unique dynamics drawn from population distribution
#' head(sim_basic$random_effects)
#'
#' # Population parameters (mean and covariance)
#' attr(sim_basic$random_effects, "mu") # Mean of dynamics distribution
#' attr(sim_basic$random_effects, "sigma") # Covariance matrix
#'
#' # Example 2: Heterogeneous patient dynamics
#' # Patients differ in damping ratio (xi) and oscillation period
#' sim_hetero <- simulate(
#'   n_subjects = 20,
#'   longitudinal = list(
#'     xi = c(mean = 0.5, sd = 0.2), # Mean damping ratio with variation
#'     period = c(mean = 8, sd = 2), # Mean period with variation
#'     excitation = list(
#'       offset = 4.0, # Constant external force
#'       covariates = c(x1 = 0.8, x2 = -0.5) # Covariate effects
#'     ),
#'     initial = list(
#'       offset = c(biomarker = 3.8, velocity = -0.1),
#'       covariates = list(
#'         biomarker = c(x1 = 0.1, x2 = -0.1),
#'         velocity = c(x1 = 0.1, x2 = -0.05)
#'       )
#'     ),
#'     error_sd = 0.1,
#'     n_measurements = 10
#'   ),
#'   covariates = list(
#'     x1 = list(type = "normal", mean = 0, sd = 1),
#'     x2 = list(type = "normal", mean = 0, sd = 1),
#'     w1 = list(type = "normal", mean = 0, sd = 1),
#'     w2 = list(type = "binary", prob = 0.5)
#'   ),
#'   seed = 456
#' )
#'
#' # Example 3: Nearly uniform dynamics (low variability)
#' # Useful for testing when patient heterogeneity should be minimal
#' sim_uniform <- simulate(
#'   n_subjects = 20,
#'   longitudinal = list(
#'     xi = c(mean = 0.707, sd = 0.01), # Very small variation
#'     period = c(mean = 5, sd = 0.1), # Very small variation
#'     excitation = list(offset = 0, covariates = numeric(0)),
#'     initial = list(
#'       offset = c(biomarker = -2, velocity = 0),
#'       covariates = list(biomarker = numeric(0), velocity = numeric(0))
#'     ),
#'     error_sd = 0.15,
#'     n_measurements = 20
#'   ),
#'   seed = 789
#' )
#'
#' # Verify low variability - all patients should have similar dynamics
#' apply(sim_uniform$random_effects, 2, sd)
#'
#' \donttest{
#' # Example 4: Comprehensive simulation with survival outcomes
#' sim_full <- simulate(
#'   n_subjects = 30,
#'   longitudinal = list(
#'     xi = c(mean = 0.4, sd = 0.05),
#'     period = c(mean = 8, sd = 0.3),
#'     excitation = list(
#'       offset = 4.0,
#'       covariates = c(x1 = 0.8, x2 = -0.5)
#'     ),
#'     initial = list(
#'       offset = c(biomarker = 3.8, velocity = -0.1),
#'       covariates = list(
#'         biomarker = c(x1 = 0.1, x2 = -0.1),
#'         velocity = c(x1 = 0.1, x2 = -0.05)
#'       )
#'     ),
#'     error_sd = 0.1,
#'     n_measurements = 20
#'   ),
#'   survival = list(
#'     baseline = list(type = "weibull", shape = 3.0, scale = 23),
#'     value = 0.4, # Effect of biomarker value on hazard
#'     slope = 1.5, # Effect of biomarker velocity on hazard
#'     gamma = 1, # Linear velocity effect
#'     covariates = c(w1 = 0.4, w2 = -0.6)
#'   ),
#'   covariates = list(
#'     x1 = list(type = "normal", mean = 0, sd = 1),
#'     x2 = list(type = "normal", mean = 0, sd = 1),
#'     w1 = list(type = "normal", mean = 0, sd = 1),
#'     w2 = list(type = "binary", prob = 0.5)
#'   ),
#'   maxt = 10,
#'   seed = 42
#' )
#'
#' # Visualize patient-specific dynamics distribution
#' library(ggplot2)
#'
#' # Recover original xi and period from dynamics parameters
#' mu <- attr(sim_full$random_effects, "mu")
#' dynamics_df <- as.data.frame(sim_full$random_effects)
#' dynamics_df$omega <- sqrt(-dynamics_df$dyn_value - mu["dyn_value"])
#' dynamics_df$xi_implied <- -(dynamics_df$dyn_slope + mu["dyn_slope"]) /
#'   (2 * dynamics_df$omega)
#' dynamics_df$period_implied <- 2 * pi / dynamics_df$omega
#'
#' # Distribution of damping ratios
#' ggplot(dynamics_df, aes(x = xi_implied)) +
#'   geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
#'   geom_vline(xintercept = 0.4, color = "red", linetype = "dashed") +
#'   labs(
#'     title = "Distribution of Patient Damping Ratios",
#'     x = "Damping Ratio (xi)", y = "Count"
#'   ) +
#'   theme_minimal()
#'
#' # Joint distribution of xi and period
#' ggplot(dynamics_df, aes(x = period_implied, y = xi_implied)) +
#'   geom_point(alpha = 0.5) +
#'   geom_density_2d(color = "steelblue") +
#'   labs(
#'     title = "Patient-Specific Dynamics",
#'     x = "Period", y = "Damping Ratio (xi)"
#'   ) +
#'   theme_minimal()
#'
#' # Compare trajectories for patients with different dynamics
#' # Select patients at 10th, 50th, 90th percentiles of xi
#' quantiles <- quantile(dynamics_df$xi_implied, probs = c(0.1, 0.5, 0.9))
#' example_ids <- sapply(quantiles, function(q) {
#'   which.min(abs(dynamics_df$xi_implied - q))
#' })
#'
#' plot_data <- sim_full$longitudinal_data[
#'   sim_full$longitudinal_data$id %in% example_ids,
#' ]
#' plot_data$xi_group <- factor(
#'   plot_data$id,
#'   levels = example_ids,
#'   labels = c("Low xi (10%)", "Median xi (50%)", "High xi (90%)")
#' )
#'
#' ggplot(plot_data, aes(x = time, y = biomarker, color = xi_group)) +
#'   geom_line(linewidth = 1) +
#'   labs(
#'     title = "Biomarker Trajectories by Damping Ratio",
#'     x = "Time", y = "Biomarker Value",
#'     color = "Patient Group"
#'   ) +
#'   theme_minimal()
#'
#' # Kaplan-Meier survival curve
#' library(survival)
#' km_fit <- survfit(Surv(time, status) ~ 1, data = sim_full$survival_data)
#' plot(km_fit,
#'   xlab = "Time", ylab = "Survival Probability",
#'   main = "Kaplan-Meier Survival Curve"
#' )
#' }
#' @concept data-simulation
#'
#' @importFrom stats rnorm rbinom
#' @importFrom utils tail
#'
#' @export
simulate <- function(
  n_subjects = 500,
  longitudinal = list(
    xi = c(mean = 0.4, sd = 0.1),
    period = c(mean = 6, sd = 0.1),
    excitation = list(
      offset = 0.0,
      covariates = c(x1 = 0.5, x2 = -0.45)
    ),
    initial = list(
      offset = c(biomarker = -0.5, velocity = -0.1),
      covariates = list(
        biomarker = c(x1 = 0.08, x2 = -0.06),
        velocity = c(x1 = 0.25, x2 = -0.2)
      )
    ),
    error_sd = 0.1,
    n_measurements = 100
  ),
  survival = list(
    baseline = list(
      type = "weibull",
      shape = 2.0,
      scale = 15.0
    ),
    value = 0.8,
    slope = 2.0,
    gamma = 1,
    covariates = c(w1 = 0.6, w2 = -0.8)
  ),
  covariates = list(
    x1 = list(type = "normal", mean = 0, sd = 1),
    x2 = list(type = "normal", mean = 0, sd = 1),
    w1 = list(type = "normal", mean = 0, sd = 1),
    w2 = list(type = "binary", prob = 0.5)
  ),
  maxt = 10,
  seed = 42
) {
  # Validate basic parameters
  stopifnot(
    "n_subjects must be a positive integer" = is.numeric(n_subjects) &&
      length(n_subjects) == 1 &&
      n_subjects > 0 &&
      n_subjects == round(n_subjects),
    "maxt must be positive" = is.numeric(maxt) && length(maxt) == 1 && maxt > 0
  )

  # Validate longitudinal structure first
  stopifnot(
    "longitudinal must be a list" = is.list(longitudinal)
  )

  # Set default for n_measurements if not provided
  if (is.null(longitudinal$n_measurements)) {
    longitudinal$n_measurements <- 10
  }

  # Validate remaining longitudinal parameters
  # Only support mean/sd format
  stopifnot(
    "xi must have mean and sd components" = !is.null(names(longitudinal$xi)) &&
      length(longitudinal$xi) == 2 &&
      all(c("mean", "sd") %in% names(longitudinal$xi)),
    "period must have mean and sd components" = !is.null(names(
      longitudinal$period
    )) &&
      length(longitudinal$period) == 2 &&
      all(c("mean", "sd") %in% names(longitudinal$period)),
    "xi mean must be positive" = longitudinal$xi["mean"] > 0,
    "xi sd must be non-negative" = longitudinal$xi["sd"] >= 0,
    "period mean must be positive" = longitudinal$period["mean"] > 0,
    "period sd must be non-negative" = longitudinal$period["sd"] >= 0
  )

  stopifnot(
    "longitudinal$error_sd must be positive" = is.numeric(
      longitudinal$error_sd
    ) &&
      longitudinal$error_sd > 0,
    "longitudinal$n_measurements must be a positive integer" = is.numeric(
      longitudinal$n_measurements
    ) &&
      longitudinal$n_measurements > 0 &&
      longitudinal$n_measurements == round(longitudinal$n_measurements),
    "longitudinal$excitation must be a list" = is.list(longitudinal$excitation),
    "excitation parameters must be numeric" = is.numeric(
      longitudinal$excitation$offset
    ),
    "excitation$covariates must be numeric vector" = is.numeric(
      longitudinal$excitation$covariates
    ) &&
      is.vector(longitudinal$excitation$covariates),
    "longitudinal$initial must be a list" = is.list(longitudinal$initial),
    "initial$offset: named vector with biomarker and velocity" = is.numeric(
      longitudinal$initial$offset
    ) &&
      length(longitudinal$initial$offset) == 2 &&
      all(c("biomarker", "velocity") %in% names(longitudinal$initial$offset)),
    "initial$covariates must be a list" = is.list(
      longitudinal$initial$covariates
    ),
    "initial$covariates must have biomarker and velocity components" = all(
      c("biomarker", "velocity") %in% names(longitudinal$initial$covariates)
    ),
    "initial$covariates$biomarker must be numeric vector" = is.numeric(
      longitudinal$initial$covariates$biomarker
    ) &&
      is.vector(longitudinal$initial$covariates$biomarker),
    "initial$covariates$velocity must be numeric vector" = is.numeric(
      longitudinal$initial$covariates$velocity
    ) &&
      is.vector(longitudinal$initial$covariates$velocity)
  )

  # Handle empty covariates (convert numeric(0) to named numeric vector)
  if (length(longitudinal$excitation$covariates) == 0) {
    longitudinal$excitation$covariates <- numeric(0)
    names(longitudinal$excitation$covariates) <- character(0)
  }
  if (length(longitudinal$initial$covariates$biomarker) == 0) {
    longitudinal$initial$covariates$biomarker <- numeric(0)
    names(longitudinal$initial$covariates$biomarker) <- character(0)
  }
  if (length(longitudinal$initial$covariates$velocity) == 0) {
    longitudinal$initial$covariates$velocity <- numeric(0)
    names(longitudinal$initial$covariates$velocity) <- character(0)
  }
  if (length(survival$covariates) == 0) {
    survival$covariates <- numeric(0)
    names(survival$covariates) <- character(0)
  }

  # Validate survival structure
  stopifnot(
    "survival must be a list" = is.list(survival),
    "survival$baseline must be a list" = is.list(survival$baseline),
    "Only weibull baseline hazard is supported" = survival$baseline$type ==
      "weibull",
    "Weibull parameters must be positive" = is.numeric(
      survival$baseline$shape
    ) &&
      survival$baseline$shape > 0 &&
      is.numeric(survival$baseline$scale) &&
      survival$baseline$scale > 0,
    "survival coefficients must be numeric" = is.numeric(survival$value) &&
      is.numeric(survival$slope),
    "gamma must be 0, 1, or 2" = survival$gamma %in% c(0, 1, 2),
    "survival$covariates must be numeric vector" = is.numeric(
      survival$covariates
    ) &&
      is.vector(survival$covariates)
  )

  # Validate covariates
  stopifnot("covariates must be a list" = is.list(covariates))
  for (cov_name in names(covariates)) {
    cov <- covariates[[cov_name]]
    if (!is.list(cov)) {
      stop(paste("Covariate", cov_name, "must be a list"))
    }
    if (!cov$type %in% c("binary", "normal")) {
      stop(paste("Covariate", cov_name, "type must be binary or normal"))
    }
    if (cov$type == "binary") {
      if (!is.numeric(cov$prob) || cov$prob < 0 || cov$prob > 1) {
        stop(paste("Covariate", cov_name, "prob must be in [0,1]"))
      }
    } else {
      if (!is.numeric(cov$mean)) {
        stop(paste("Covariate", cov_name, "mean must be numeric"))
      }
      if (!is.numeric(cov$sd) || cov$sd <= 0) {
        stop(paste("Covariate", cov_name, "sd must be positive"))
      }
    }
  }

  # Validate dimension consistency
  long_cov_names <- names(longitudinal$excitation$covariates)
  init_cov_names <- unique(c(
    names(longitudinal$initial$covariates$biomarker),
    names(longitudinal$initial$covariates$velocity)
  ))

  # Handle survival covariates - if user passes numeric(0), ignore defaults
  if (length(survival$covariates) == 0) {
    surv_cov_names <- character(0)
  } else {
    surv_cov_names <- names(survival$covariates)
  }

  # All mentioned covariates should exist in covariates list
  all_cov_names <- unique(c(long_cov_names, init_cov_names, surv_cov_names))
  if (length(all_cov_names) > 0) {
    missing_covs <- setdiff(all_cov_names, names(covariates))
    if (length(missing_covs) > 0) {
      stop(paste(
        "Missing covariate definitions:",
        paste(missing_covs, collapse = ", ")
      ))
    }
  }

  set.seed(seed)
  x <- .generate_covariates(n_subjects, covariates)

  # Sample patient-specific dynamics
  dynamics <- .generate_patient_dynamics(
    n_subjects,
    longitudinal$xi,
    longitudinal$period,
    longitudinal$excitation
  )
  x_init <- x[, c("id", init_cov_names), drop = FALSE]
  x_long <- x[, c("id", long_cov_names), drop = FALSE]
  x_surv <- x[, c("id", surv_cov_names), drop = FALSE]

  init <- .compute_initial_biomarker(
    x_init,
    longitudinal$initial
  )

  surv <- .generate_survival_data(
    x_surv,
    x_long,
    dynamics,
    init,
    longitudinal$excitation,
    survival,
    maxt
  )

  long <- .generate_longitudinal_data(
    x_long,
    dynamics,
    init,
    surv,
    longitudinal$n_measurements,
    longitudinal$error_sd
  )

  long_data <- merge(long, x_long, by = "id")

  surv_data <- merge(surv, x_surv, by = "id")
  names(surv_data)[names(surv_data) == "eventtime"] <- "time"

  init_final <- as.matrix(init[, c("biomarker", "velocity")])

  mu <- attr(dynamics, "mu")
  sigma <- attr(dynamics, "sigma")
  # ODE coefficient RE: deviations from population mean
  coef_re <- dynamics[, c("dyn_value", "dyn_slope"), drop = FALSE]
  coef_re <- sweep(coef_re, 2, mu[c("dyn_value", "dyn_slope")], "-")
  coef_re <- as.matrix(coef_re)

  # State RE: per-subject initial state deviations from population mean
  pop_m0 <- mean(init_final[, 1])
  pop_v0 <- mean(init_final[, 2])
  state_re <- cbind(
    init_final[, 1] - pop_m0,
    init_final[, 2] - pop_v0
  )

  # Full RE: [b_m0, b_v0, b_coef...]
  random_effects <- cbind(state_re, coef_re)
  colnames(random_effects) <- c("m0", "v0", "dyn_value", "dyn_slope")
  attr(random_effects, "mu") <- mu[c("dyn_value", "dyn_slope")]
  attr(random_effects, "sigma") <- sigma[1:2, 1:2]

  list(
    longitudinal_data = long_data,
    survival_data = surv_data,
    state = init_final,
    initial_state_mean = c(m0 = pop_m0, v0 = pop_v0),
    random_effects = random_effects
  )
}

.create_example_data <- function(
  n_subjects = 200,
  n_measurements = 100,
  seed = 42
) {
  coef_args <- formals(simulate)
  long_params <- eval(coef_args$longitudinal)
  long_params$n_measurements <- n_measurements

  data <- simulate(
    n_subjects = n_subjects,
    longitudinal = long_params,
    seed = seed
  )
  f0 <- function(t) {
    baseline_type <- eval(coef_args$survival$baseline$type)
    res <- switch(baseline_type,
      weibull = {
        shape <- eval(coef_args$survival$baseline$shape)
        scale <- eval(coef_args$survival$baseline$scale)
        (shape / scale) * (t / scale)^(shape - 1)
      },
      stop(paste(
        "Unsupported baseline hazard type:",
        baseline_type
      ))
    )
    log(res)
  }
  spline_config <- formals(JointODE)$spline_baseline
  spline_config <- .get_spline_config(
    x = data$survival_data$time,
    degree = spline_config$degree,
    n_knots = spline_config$n_knots,
    knot_placement = spline_config$knot_placement,
    boundary_knots = spline_config$boundary_knots
  )
  spline_config$boundary_knots[1] <- 0
  baseline_coef <- .estimate_bspline_coef(
    x = data$survival_data$time,
    f0 = f0,
    config = spline_config
  )
  dynamics_params <- .compute_dynamics_parameters(
    long_params$xi,
    long_params$period,
    long_params$excitation
  )
  longitudinal_coef <- dynamics_params$mu
  random_effect_sigma <- dynamics_params$sigma

  hazard_coef <- c(
    eval(coef_args$survival$value),
    eval(coef_args$survival$slope),
    eval(coef_args$survival$covariates)
  )
  # Full RE sigma: [state(2), coef(2)]
  # State RE variance from simulation parameters
  cov_defs <- eval(coef_args$covariates)
  init_covs <- long_params$initial$covariates
  # Compute Var(X_j) for each covariate
  cov_variance <- function(nm, cov_defs) {
    d <- cov_defs[[nm]]
    if (d$type == "normal") d$sd^2
    else if (d$type == "binary") d$prob * (1 - d$prob)
    else 1
  }
  # State RE covariance: Cov(m0, v0) = sum(beta_m_j * beta_v_j * Var(X_j))
  beta_m <- init_covs$biomarker
  beta_v <- init_covs$velocity
  shared_covs <- intersect(names(beta_m), names(beta_v))

  var_m0 <- sum(vapply(names(beta_m), function(nm) {
    beta_m[nm]^2 * cov_variance(nm, cov_defs)
  }, numeric(1)))
  var_v0 <- sum(vapply(names(beta_v), function(nm) {
    beta_v[nm]^2 * cov_variance(nm, cov_defs)
  }, numeric(1)))
  cov_m0_v0 <- sum(vapply(shared_covs, function(nm) {
    beta_m[nm] * beta_v[nm] * cov_variance(nm, cov_defs)
  }, numeric(1)))

  # Full RE sigma: [state(2), coef(2)]
  # State-coef cross-covariance = 0 (independent)
  n_coef_re <- nrow(random_effect_sigma)
  n_re_total <- n_coef_re + 2
  full_sigma <- matrix(0, n_re_total, n_re_total)
  full_sigma[1, 1] <- var_m0
  full_sigma[2, 2] <- var_v0
  full_sigma[1, 2] <- full_sigma[2, 1] <- cov_m0_v0
  full_sigma[3:n_re_total, 3:n_re_total] <- random_effect_sigma

  parameters <- list(
    coefficients = list(
      baseline = baseline_coef,
      hazard = hazard_coef,
      longitudinal = longitudinal_coef,
      initial_state = data$initial_state_mean,
      measurement_error_sd = eval(coef_args$longitudinal$error_sd),
      random_effect_sigma = full_sigma
    ),
    configurations = list(
      baseline = spline_config,
      gamma = eval(coef_args$survival$gamma),
      biomarker = list(fixed = TRUE, random = TRUE),
      velocity = list(fixed = TRUE, random = TRUE)
    )
  )
  list(data = data, init = parameters)
}

.estimate_bspline_coef <- function(x, f0, config) {
  basis_matrix <- splines2::bSpline(
    x,
    knots = config$knots,
    Boundary.knots = config$boundary_knots,
    degree = config$degree,
    intercept = TRUE
  )

  y_target <- f0(x)
  if (!is.matrix(y_target)) {
    y_target <- matrix(y_target, ncol = 1)
  }

  spline_coefficients <- solve(
    t(basis_matrix) %*% basis_matrix,
    t(basis_matrix) %*% y_target
  )
  as.vector(spline_coefficients)
}

.generate_covariates <- function(n_subjects, covariates) {
  covariate_data <- data.frame(id = seq_len(n_subjects))
  for (cov_name in names(covariates)) {
    cov_info <- covariates[[cov_name]]
    if (cov_info$type == "binary") {
      covariate_data[[cov_name]] <- rbinom(n_subjects, 1, cov_info$prob)
    } else if (cov_info$type == "normal") {
      covariate_data[[cov_name]] <- rnorm(
        n_subjects,
        cov_info$mean,
        cov_info$sd
      )
    } else {
      stop(paste("Unsupported covariate type:", cov_info$type))
    }
  }
  covariate_data
}

.compute_dynamics_parameters <- function(xi, period, excitation) {
  mu_period <- period["mean"]
  sigma_period <- period["sd"]
  mu_xi <- xi["mean"]
  sigma_xi <- xi["sd"]
  mu_omega <- 2 * pi / mu_period
  sigma_omega <- (2 * pi / mu_period^2) * sigma_period

  mu <- c(
    -mu_omega^2,
    -2 * mu_xi * mu_omega,
    mu_omega^2 * excitation$offset,
    mu_omega^2 * excitation$covariates
  )
  cov_names <- if (length(excitation$covariates) > 0) {
    paste0("dyn_", names(excitation$covariates))
  } else {
    character(0)
  }
  names(mu) <- c("dyn_value", "dyn_slope", "dyn_const", cov_names)

  jacobian <- matrix(0, nrow = 2, ncol = 2)
  jacobian[1, 1] <- -2 * mu_omega
  jacobian[2, 1] <- -2 * mu_xi
  jacobian[2, 2] <- -2 * mu_omega

  sigma_orig <- diag(c(sigma_omega^2, sigma_xi^2))
  sigma <- jacobian %*% sigma_orig %*% t(jacobian)

  list(mu = mu, sigma = sigma)
}

.generate_patient_dynamics <- function(n_subjects, xi, period, excitation) {
  params <- .compute_dynamics_parameters(xi, period, excitation)

  random_effects_2d <- MASS::mvrnorm(
    n_subjects,
    mu = params$mu[1:2],
    Sigma = params$sigma
  )

  if (n_subjects == 1) {
    random_effects_2d <- matrix(random_effects_2d, nrow = 1)
  }

  n_beta <- length(params$mu)
  random_effects <- matrix(0, nrow = n_subjects, ncol = n_beta)
  random_effects[, 1:2] <- random_effects_2d
  random_effects[, 3:n_beta] <- matrix(
    rep(params$mu[3:n_beta], each = n_subjects),
    nrow = n_subjects,
    byrow = FALSE
  )
  colnames(random_effects) <- names(params$mu)

  dynamics <- data.frame(id = seq_len(n_subjects), random_effects)
  attr(dynamics, "mu") <- params$mu
  attr(dynamics, "sigma") <- params$sigma
  dynamics
}

.compute_initial_biomarker <- function(x, initial) {
  biomarker_offset <- initial$offset["biomarker"]
  biomarker_covs <- initial$covariates$biomarker

  if (length(biomarker_covs) > 0) {
    x_biomarker <- x[, names(biomarker_covs), drop = FALSE]
    biomarker_value <- biomarker_offset +
      as.matrix(x_biomarker) %*% biomarker_covs
  } else {
    biomarker_value <- rep(biomarker_offset, nrow(x))
  }

  velocity_offset <- initial$offset["velocity"]
  velocity_covs <- initial$covariates$velocity

  if (length(velocity_covs) > 0) {
    x_velocity <- x[, names(velocity_covs), drop = FALSE]
    velocity_value <- velocity_offset +
      as.matrix(x_velocity) %*% velocity_covs
  } else {
    velocity_value <- rep(velocity_offset, nrow(x))
  }

  data.frame(id = x$id, biomarker = biomarker_value, velocity = velocity_value)
}

.solve_biomarker_ode <- function(times, x, init, dynamic) {
  .ode_deriv <- function(t, state, parms) {
    list(c(
      state[2],
      parms$offset + parms$value * state[1] + parms$slope * state[2]
    ))
  }

  x_vec <- as.numeric(x)
  dynamic_vec <- as.numeric(dynamic)

  ode_parms <- list(
    offset = dynamic_vec[3] + sum(dynamic_vec[-(1:3)] * x_vec),
    value = dynamic_vec[1],
    slope = dynamic_vec[2]
  )

  ode_solution <- deSolve::ode(
    y = c(init$biomarker, init$velocity),
    times = sort(c(0, times)),
    func = .ode_deriv,
    parms = ode_parms
  )
  idx <- match(times, ode_solution[, 1])
  biomarker <- ode_solution[idx, 2]
  velocity <- ode_solution[idx, 3]
  acceleration <- rep(ode_parms$offset, length(times)) +
    ode_parms$value * biomarker +
    ode_parms$slope * velocity
  data.frame(
    time = times,
    biomarker = biomarker,
    velocity = velocity,
    acceleration = acceleration
  )
}

.generate_survival_data <- function(
  x_surv,
  x_long,
  dynamics,
  init,
  excitation,
  survival,
  maxt
) {
  hazard_function <- function(t, x, betas, ...) {
    h0 <- switch(survival$baseline$type,
      weibull = {
        shape <- survival$baseline$shape
        scale <- survival$baseline$scale
        (shape / scale) * (t / scale)^(shape - 1)
      },
      stop(paste("Unsupported baseline hazard type:", survival$baseline$type))
    )

    x_surv <- x[names(survival$covariates)]
    x_long <- x[names(excitation$covariates)]
    init <- x[c("biomarker", "velocity")]
    dyn_cols <- c("dyn_value", "dyn_slope", "dyn_const")
    if (length(excitation$covariates) > 0) {
      dyn_cols <- c(dyn_cols, paste0("dyn_", names(excitation$covariates)))
    }
    dynamic <- x[, dyn_cols, drop = FALSE]

    biomarker <- .solve_biomarker_ode(t, x_long, init, dynamic)

    biomarker_term <- survival$value * biomarker$biomarker
    velocity_term <- if (survival$gamma == 0) {
      0
    } else if (survival$gamma == 1) {
      survival$slope * biomarker$velocity
    } else if (survival$gamma == 2) {
      survival$slope * biomarker$velocity^2
    } else {
      stop("gamma must be 0, 1, or 2")
    }
    covariate_term <- sum(x_surv * survival$covariates)

    linpred <- biomarker_term + velocity_term + covariate_term
    h0 * exp(linpred)
  }

  covs <- merge(x_surv, x_long, by = "id")
  covs <- merge(covs, init, by = "id")
  covs <- merge(covs, dynamics, by = "id")

  simsurv::simsurv(hazard = hazard_function, x = covs, maxt = maxt)
}

.generate_longitudinal_data <- function(
  x,
  dynamics,
  init,
  surv,
  n_measurements,
  error_sd
) {
  n <- nrow(surv)
  data_list <- vector("list", n)
  maxt <- max(surv$eventtime)

  for (i in seq_len(n)) {
    times <- seq(
      0,
      surv[i, "eventtime"],
      by = maxt / n_measurements
    )
    if (tail(times, 1) == surv[i, "eventtime"]) {
      times <- times[-length(times)]
    }

    patient_id <- surv[i, "id"]
    biomarkers <- .solve_biomarker_ode(
      times,
      x[x$id == patient_id, names(x) != "id"],
      init[init$id == patient_id, names(init) != "id"],
      dynamics[dynamics$id == patient_id, names(dynamics) != "id"]
    )
    biomarkers$id <- patient_id

    measurement_error <- rnorm(
      nrow(biomarkers),
      mean = 0,
      sd = error_sd
    )

    biomarkers$observed <- biomarkers$biomarker + measurement_error

    data_list[[i]] <- biomarkers
  }
  long <- do.call(rbind, data_list)
  long
}
