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
#'   (default: 200)
#' @param shared_sd Positive scalar defining the standard deviation of the
#'   shared random effects, \eqn{\sigma_b} (default: 0.1)
#' @param longitudinal List specifying the longitudinal sub-model parameters:
#'   \describe{
#'     \item{xi}{Damping ratio \eqn{\xi} controlling the oscillation decay.
#'       Values: \eqn{\xi < 1} (underdamped), \eqn{\xi = 1} (critically damped),
#'       \eqn{\xi > 1} (overdamped) (default: 0.707, slightly underdamped)}
#'     \item{period}{Natural period \eqn{T} of oscillation in time units,
#'       related to natural frequency as \eqn{\omega = 2\pi/T} (default: 5)}
#'     \item{k}{Excitation amplitude scaling factor \eqn{k} that modulates
#'       the strength of external forcing (default: 1.0)}
#'     \item{excitation}{List specifying external forcing parameters:
#'       \describe{
#'         \item{offset}{Constant excitation term \eqn{f_0} (default: 0.0)}
#'         \item{covariates}{Named vector of covariate effects
#'           \eqn{\boldsymbol{\beta}_{exc}} on excitation
#'           (default: c(x1 = 0.8, x2 = -0.5))}
#'       }
#'     }
#'     \item{initial}{List specifying initial condition parameters:
#'       \describe{
#'         \item{offset}{Baseline initial biomarker value \eqn{\mu_0}
#'           (default: -3.0)}
#'         \item{covariates}{Named vector of covariate effects
#'           \eqn{\boldsymbol{\beta}_{init}} on initial biomarker level
#'           (default: c(x1 = 0.1, x2 = -0.1))}
#'         \item{random_coef}{Scaling coefficient \eqn{\psi} for random effect
#'           influence on initial conditions (default: 0.2)}
#'       }
#'     }
#'     \item{error_sd}{Standard deviation \eqn{\sigma_{\epsilon}} of the
#'       measurement error process (default: 0.1)}
#'     \item{n_measurements}{Number of longitudinal measurements per subject
#'       (default: 10)}
#'   }
#' @param survival List specifying the survival sub-model parameters:
#'   \describe{
#'     \item{baseline}{List defining the Weibull baseline hazard function:
#'       \describe{
#'         \item{type}{Character string specifying the baseline hazard type
#'           (currently only "weibull" is supported)}
#'         \item{shape}{Weibull shape parameter \eqn{\kappa > 0}
#'           (default: 1.0)}
#'         \item{scale}{Weibull scale parameter \eqn{\lambda > 0}
#'           (default: 8.0)}
#'       }
#'     }
#'     \item{value}{Association parameter \eqn{\alpha_1} linking current
#'       biomarker value to hazard (default: 0.3)}
#'     \item{slope}{Association parameter \eqn{\alpha_2} linking biomarker
#'       velocity to hazard (default: 1.0)}
#'     \item{covariates}{Named vector of regression coefficients
#'       \eqn{\boldsymbol{\phi}} for survival covariates
#'       (default: c(w1 = 0.4, w2 = -0.6))}
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
#' @return A list containing three components:
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
#'         Time-invariant covariates
#'     }
#'   }
#'   \item{\code{survival_data}}{Data frame containing time-to-event data
#'     with columns:
#'     \itemize{
#'       \item \code{id}: Subject identifier
#'       \item \code{time}: Observed event or censoring time, \eqn{T_i}
#'       \item \code{status}: Event indicator, \eqn{\delta_i}
#'         (1 = event observed, 0 = censored)
#'       \item \code{b}: Realized subject-specific random effect, \eqn{b_i}
#'       \item \code{w1}, \code{w2}:
#'         Baseline survival covariates
#'     }
#'   }
#'   \item{\code{state}}{An \eqn{n \times 2} matrix containing initial states
#'     \eqn{[m_i(0), \dot{m}_i(0)]} for each subject}
#' }
#'
#' @details
#' The simulation framework implements a joint model comprising longitudinal
#' and survival sub-models linked through shared random effects and
#' trajectory-dependent associations.
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
#'   \item
#'     \eqn{m_i(0) = \mu_0 + \mathbf{X}_i^T\boldsymbol{\beta}_{init} + \psi b_i}
#'   \item \eqn{\dot{m}_i(0) = 0} (zero initial velocity)
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
#' and velocity:
#' \deqn{\lambda_i(t) = \lambda_0(t)\exp(\alpha_1 m_i(t) + \alpha_2\dot{m}_i(t)
#'  + \mathbf{W}_i^T\boldsymbol{\phi} + b_i)}
#' where:
#' \itemize{
#'   \item \eqn{\lambda_0(t)} denotes the Weibull baseline hazard:
#'     \eqn{\lambda_0(t) = (\kappa/\lambda)(t/\lambda)^{\kappa-1}}
#'   \item \eqn{\alpha_1} quantifies the association with current biomarker
#'     value
#'   \item \eqn{\alpha_2} quantifies the association with biomarker velocity
#'   \item \eqn{\mathbf{W}_i} contains baseline covariates
#'   \item \eqn{\boldsymbol{\phi}} represents covariate effects on survival
#' }
#' }
#'
#' \subsection{Random Effects Structure}{
#' Subject-specific random effects \eqn{b_i \sim \mathcal{N}(0, \sigma_b^2)}
#' induce correlation
#' between the longitudinal and survival processes through:
#' \itemize{
#'   \item Initial biomarker level via \eqn{m_i(0)}
#'   \item Hazard function via the frailty term
#'   \item Longitudinal observations via additive shift
#' }
#' }
#'
#' @examples
#' # Generate a small dataset with default parameters
#' sim_data <- simulate(n_subjects = 50, seed = 123)
#'
#' # Examine the structure of generated data
#' str(sim_data)
#'
#' \dontrun{
#' # Customized simulation with underdamped oscillator
#' sim_data <- simulate(
#'   n_subjects = 200,
#'   longitudinal = list(
#'     xi = 0.3,          # Low damping (more oscillatory)
#'     period = 3,        # Faster oscillations
#'     k = 1.5,           # Stronger excitation
#'     excitation = list(
#'       offset = 0.5,
#'       covariates = c(x1 = 1.0, x2 = -0.3)
#'     ),
#'     initial = list(
#'       offset = 0.0,    # Start at equilibrium
#'       covariates = c(x1 = 0.5, x2 = -0.2),
#'       random_coef = 0.5
#'     ),
#'     error_sd = 0.2     # Moderate measurement error
#'   )
#' )
#'
#' # Overdamped system (no oscillations)
#' sim_data_overdamped <- simulate(
#'   n_subjects = 100,
#'   longitudinal = list(
#'     xi = 2.0,          # Overdamped (exponential decay)
#'     period = 10,       # Slow dynamics
#'     k = 0.5,           # Weak excitation
#'     excitation = list(
#'       offset = -1.0,
#'       covariates = c(x1 = 0.2, x2 = 0.1)
#'     )
#'   )
#' )
#' }
#' @concept data-simulation
#'
#' @importFrom stats rnorm rbinom
#' @importFrom utils tail
#'
#' @export
simulate <- function(
  n_subjects = 200,
  shared_sd = 0.1,
  longitudinal = list(
    xi = 0.707,
    period = 5,
    k = 1.0,
    excitation = list(
      offset = 0.0,
      covariates = c(x1 = 0.8, x2 = -0.5)
    ),
    initial = list(
      offset = -3.0,
      covariates = c(x1 = 0.1, x2 = -0.1),
      random_coef = 0.2
    ),
    error_sd = 0.1,
    n_measurements = 10
  ),
  survival = list(
    baseline = list(
      type = "weibull",
      shape = 1.0,
      scale = 8.0
    ),
    value = 0.3,
    slope = 1.0,
    covariates = c(w1 = 0.4, w2 = -0.6)
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
    "shared_sd must be positive" = is.numeric(shared_sd) &&
      length(shared_sd) == 1 &&
      shared_sd > 0,
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
  stopifnot(
    "longitudinal parameters must be numeric" = is.numeric(
      longitudinal$xi
    ) &&
      is.numeric(longitudinal$period) &&
      is.numeric(longitudinal$k),
    "longitudinal$period must be positive" = longitudinal$period > 0,
    "longitudinal$k must be numeric" = is.numeric(longitudinal$k),
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
    "initial parameters must be numeric" = is.numeric(
      longitudinal$initial$offset
    ) &&
      is.numeric(longitudinal$initial$random_coef),
    "initial$covariates must be numeric vector" = is.numeric(
      longitudinal$initial$covariates
    ) &&
      is.vector(longitudinal$initial$covariates)
  )

  # Handle empty covariates (convert numeric(0) to named numeric vector)
  if (length(longitudinal$excitation$covariates) == 0) {
    longitudinal$excitation$covariates <- numeric(0)
    names(longitudinal$excitation$covariates) <- character(0)
  }
  if (length(longitudinal$initial$covariates) == 0) {
    longitudinal$initial$covariates <- numeric(0)
    names(longitudinal$initial$covariates) <- character(0)
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
  init_cov_names <- names(longitudinal$initial$covariates)

  # Handle survival covariates - if user passes numeric(0), ignore defaults
  if (
    length(survival$covariates) == 1 &&
      is.numeric(survival$covariates) &&
      length(survival$covariates) == 0
  ) {
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
  b <- .generate_shared_effects(n_subjects, shared_sd)
  x <- .generate_covariates(n_subjects, covariates)
  x_init <- x[, c("id", init_cov_names), drop = FALSE]
  x_long <- x[, c("id", long_cov_names), drop = FALSE]
  x_surv <- x[, c("id", surv_cov_names), drop = FALSE]
  init <- .compute_initial_biomarker(
    x_init,
    b,
    longitudinal$initial
  )
  surv <- .generate_survival_data(
    x_surv,
    x_long,
    b,
    init,
    longitudinal,
    survival,
    maxt
  )
  long <- .generate_longitudinal_data(
    x_long,
    b,
    init,
    surv,
    longitudinal
  )

  surv_data <- merge(surv, x_surv, by = "id")
  surv_data <- merge(surv_data, b, by = "id")
  names(surv_data)[names(surv_data) == "eventtime"] <- "time"

  long_data <- merge(long, x_long, by = "id")
  col_order <- c(
    "id",
    "time",
    "observed",
    "biomarker",
    "velocity",
    "acceleration",
    long_cov_names
  )
  long_data <- long_data[, col_order]

  init_final <- init[, c("biomarker", "velocity")]

  list(
    longitudinal_data = long_data,
    survival_data = surv_data,
    state = init_final
  )
}

.create_example_data <- function(n_subjects = 200, seed = 123) {
  data <- simulate(n_subjects = n_subjects, seed = seed)
  coef_args <- formals(simulate)
  f0 <- function(t) {
    baseline_type <- eval(coef_args$survival$baseline$type)
    res <- switch(
      baseline_type,
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

  omega <- 2 * pi / eval(coef_args$longitudinal$period)
  xi <- eval(coef_args$longitudinal$xi)
  k <- eval(coef_args$longitudinal$k)
  excitation_value <- -omega^2
  excitation_slope <- -2 * xi * omega
  excitation_offset <- k *
    omega^2 *
    eval(coef_args$longitudinal$excitation$offset)
  excitation_coefs <- k *
    omega^2 *
    eval(coef_args$longitudinal$excitation$covariates)

  acceleration_coef <- c(
    excitation_value,
    excitation_slope,
    excitation_offset,
    excitation_coefs
  )
  hazard_coef <- c(
    eval(coef_args$survival$value),
    eval(coef_args$survival$slope),
    eval(coef_args$survival$covariates)
  )
  parameters <- list(
    coefficients = list(
      baseline = baseline_coef,
      acceleration = acceleration_coef,
      hazard = hazard_coef,
      measurement_error_sd = eval(coef_args$longitudinal$error_sd),
      random_effect_sd = eval(coef_args$shared_sd)
    ),
    configurations = list(
      baseline = spline_config,
      autonomous = TRUE
    )
  )
  list(data = data, init = parameters)
}

.estimate_bspline_coef <- function(x, f0, config) {
  # Compute B-spline basis at grid points
  basis_matrix <- .compute_spline_basis(x, config)

  # Compute target values
  y_target <- f0(x)

  # Ensure y_target is a column vector
  if (!is.matrix(y_target)) {
    y_target <- matrix(y_target, ncol = 1)
  }

  # Fit coefficients using least squares
  spline_coefficients <- solve(
    t(basis_matrix) %*% basis_matrix,
    t(basis_matrix) %*% y_target
  )
  as.vector(spline_coefficients)
}

.generate_shared_effects <- function(n_subjects, shared_sd) {
  data.frame(
    id = seq_len(n_subjects),
    b = rnorm(n_subjects, mean = 0, sd = shared_sd)
  )
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

.compute_initial_biomarker <- function(x, b, initial) {
  value <- initial$offset +
    as.matrix(x[, names(x) != "id"]) %*% initial$covariates +
    initial$random_coef * b[, names(b) != "id"]
  slope <- rep(0, nrow(x))
  data.frame(id = x$id, biomarker = value, velocity = slope)
}

.solve_biomarker_ode <- function(times, x, init, longitudinal) {
  .ode_deriv <- function(t, state, parms) {
    biomarker <- state[1]
    velocity <- state[2]
    acceleration <- parms$offset +
      parms$value * biomarker +
      parms$slope * velocity
    list(c(velocity, acceleration))
  }
  omega <- 2 * pi / longitudinal$period
  offset <- longitudinal$k *
    omega^2 *
    (longitudinal$excitation$offset +
      sum(x * longitudinal$excitation$covariates))
  ode_parms <- list(
    offset = offset,
    value = -omega^2,
    slope = -2 * longitudinal$xi * omega
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
  acceleration <- rep(offset, length(times)) +
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
  b,
  init,
  longitudinal,
  survival,
  maxt
) {
  # Define hazard function
  hazard_function <- function(t, x, betas, ...) {
    # baseline hazard
    h0 <- switch(
      survival$baseline$type,
      weibull = {
        shape <- survival$baseline$shape
        scale <- survival$baseline$scale
        (shape / scale) * (t / scale)^(shape - 1)
      },
      stop(paste("Unsupported baseline hazard type:", survival$baseline$type))
    )
    x_surv <- x[names(survival$covariates)]
    x_long <- x[names(longitudinal$excitation$covariates)]
    init <- x[c("biomarker", "velocity")]

    biomarker <- .solve_biomarker_ode(
      t,
      x_long,
      init,
      longitudinal
    )

    linpred <- survival$value *
      biomarker$biomarker +
      survival$slope * biomarker$velocity +
      sum(x_surv * survival$covariates) +
      x$b

    h0 * exp(linpred)
  }

  covs <- merge(x_surv, x_long, by = "id")
  covs <- merge(covs, b, by = "id")
  covs <- merge(covs, init, by = "id")

  # Generate event times
  simsurv::simsurv(
    hazard = hazard_function,
    x = covs,
    maxt = maxt
  )
}

.generate_longitudinal_data <- function(x, b, init, surv, longitudinal) {
  n <- nrow(surv)
  data_list <- vector("list", n)
  maxt <- max(surv$eventtime)
  for (i in seq_len(n)) {
    times <- seq(
      0,
      surv[i, "eventtime"],
      by = maxt / longitudinal$n_measurements
    )
    if (tail(times, 1) == surv[i, "eventtime"]) {
      times <- times[-length(times)]
    }

    biomarkers <- .solve_biomarker_ode(
      times,
      x[x$id == surv[i, "id"], names(x) != "id"],
      init[init$id == surv[i, "id"], names(init) != "id"],
      longitudinal
    )
    biomarkers$id <- surv[i, "id"]
    random_effect <- b[b$id == surv[i, "id"], "b"]
    measurement_error <- rnorm(
      nrow(biomarkers),
      mean = 0,
      sd = longitudinal$error_sd
    )

    biomarkers$observed <- biomarkers$biomarker +
      random_effect +
      measurement_error

    data_list[[i]] <- biomarkers
  }
  long <- do.call(rbind, data_list)
  long
}
