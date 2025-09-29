# Individual ODE Utility Functions

#' Individual ODE Parameter Estimation
#'
#' @description
#' Estimates ODE parameters for individual longitudinal biomarker trajectories
#' using the second-order differential equation:
#' \deqn{\ddot{m}(t) = \text{value} \cdot m(t) + \text{slope} \cdot \dot{m}(t) +
#'   \mathbf{X}(t)^T \boldsymbol{\beta}_{\text{cov}}}
#'
#' @param formula A formula object specifying the longitudinal model
#'   (e.g., \code{biomarker ~ x1 + x2})
#' @param data A data frame containing the longitudinal measurements
#' @param time Character string specifying the time variable name
#'   (default: \code{"time"})
#' @param id Character string specifying the subject identifier variable name
#'   (default: \code{"id"})
#' @param state Optional matrix of initial conditions with two columns:
#'   \itemize{
#'     \item Column 1: Initial biomarker values \eqn{m(0)}
#'     \item Column 2: Initial velocities \eqn{\dot{m}(0)}
#'   }
#'   Each row corresponds to one subject. If \code{NULL}, initial values are
#'   estimated from the data.
#' @param control A list of control parameters passed to \code{optim}
#'
#' @concept model-fitting
#'
#' @return For a single subject, returns a list containing:
#'   \describe{
#'     \item{\code{subject_id}}{Subject identifier}
#'     \item{\code{coefficients}}{Named list with estimated parameters:
#'       \itemize{
#'         \item \code{value}: Biomarker value coefficient
#'         \item \code{slope}: Biomarker velocity/slope coefficient
#'         \item \code{covariates}: Named vector of covariate effects
#'       }}
#'     \item{\code{initial_state}}{Vector with initial values \code{c(m0, v0)}}
#'     \item{\code{fitted_values}}{Vector of fitted biomarker values}
#'     \item{\code{residuals}}{Vector of residuals (observed - fitted)}
#'     \item{\code{sse}}{Sum of squared errors (objective function value)}
#'     \item{\code{convergence}}{Logical indicating convergence status}
#'     \item{\code{message}}{Optimization message from \code{optim}}
#'   }
#'   For multiple subjects, returns a named list with one element per subject.
#'
#' @examples
#' \dontrun{
#' # Generate simulated data
#' sim <- simulate(n_subjects = 5)
#'
#' # Fit individual ODE models with covariates
#' fit <- IndividualODE(
#'   formula = observed ~ x1 + x2,
#'   data = sim$longitudinal_data,
#'   state = as.matrix(sim$state)
#' )
#'
#' # Fit without covariates (intercept only)
#' fit2 <- IndividualODE(
#'   formula = observed ~ 1,
#'   data = sim$longitudinal_data
#' )
#' }
#'
#' @export
IndividualODE <- function(
  formula,
  data,
  time = "time",
  id = "id",
  state = NULL,
  control = list()
) {
  # Capture call for potential debugging
  cl <- match.call()

  # Process formula and extract components
  processed <- .process_individual_longitudinal_data(
    formula = formula,
    data = data,
    time = time,
    id = id,
    state = state
  )

  # Set optimization control defaults
  control <- modifyList(
    list(maxit = 1000, factr = 1e7, pgtol = 1e-8),
    control
  )

  # Fit models for each subject
  results <- lapply(processed, function(subject_data) {
    if (is.null(subject_data)) {
      return(NULL)
    }
    .fit_individual_ode(
      data = subject_data,
      control = control
    )
  })

  names(results) <- names(processed)

  # Return single result if only one subject
  if (length(results) == 1) {
    results[[1]]
  } else {
    structure(
      results,
      class = "IndividualODE",
      call = cl
    )
  }
}


# ===== SECTION: INDIVIDUAL ODE FITTING =====

# Process longitudinal data for individual ODE fitting
.process_individual_longitudinal_data <- function(
  formula,
  data,
  time,
  id,
  state
) {
  # Convert matrix to data frame if necessary
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  # Validate inputs
  stopifnot(
    "Data cannot be empty" = nrow(data) > 0,
    "Data must contain a time column" = time %in% names(data),
    "Data must contain an id column" = id %in% names(data)
  )

  # Parse formula
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)

  stopifnot(
    "Formula must include a response variable" = !is.null(y)
  )

  X <- model.matrix(formula, data = data)

  # Extract time and id
  times <- data[[time]]
  ids <- data[[id]]

  # Get unique subjects
  subjects <- unique(ids)
  n_subjects <- length(subjects)

  # Validate state if provided
  if (!is.null(state)) {
    stopifnot(
      "state must be a matrix with 2 columns [m(0), m'(0)]" = is.matrix(
        state
      ) &&
        ncol(state) == 2,
      "state must have one row per subject" = nrow(state) == n_subjects
    )
  }

  # Process each subject's data
  processed <- lapply(seq_along(subjects), function(i) {
    subject_id <- subjects[i]
    idx <- which(ids == subject_id)

    # Validate subject has data
    if (length(idx) == 0) {
      stop(sprintf("No data found for subject %s", subject_id))
    }

    # Sort by time
    time_order <- order(times[idx])
    idx_sorted <- idx[time_order]

    # Extract sorted data
    y_subj <- y[idx_sorted]
    t_subj <- times[idx_sorted]
    x_subj <- X[idx_sorted, , drop = FALSE]

    # Skip subjects with only one observation at time 0
    if (length(t_subj) == 1 && t_subj[1] == 0) {
      warning(sprintf(
        "Subject %s has only one observation at t=0. Skipping ODE estimation.",
        subject_id
      ))
      return(NULL)
    }

    # Determine initial state
    if (!is.null(state)) {
      initial_state <- c(state[i, 1], state[i, 2])
    } else {
      # Estimate from data
      m0 <- y_subj[1]
      v0 <- if (length(y_subj) > 1) {
        dt <- t_subj[2] - t_subj[1]
        if (abs(dt) < .Machine$double.eps) {
          0 # Avoid division by near-zero
        } else {
          (y_subj[2] - y_subj[1]) / dt
        }
      } else {
        0
      }
      initial_state <- c(m0, v0)
    }

    # Return structured data for this subject
    list(
      subject = subject_id,
      response = y_subj,
      time = t_subj,
      covariates = x_subj,
      initial = initial_state
    )
  })

  names(processed) <- as.character(subjects)
  processed
}


# Fit ODE parameters for individual subject
.fit_individual_ode <- function(data, control) {
  n_cov <- ncol(data$covariates)

  # Initialize parameter vector: [value, slope, covariates]
  theta_init <- c(
    value = 0,
    slope = 0,
    rep(0, n_cov)
  )

  # Define objective function (sum of squared errors)
  objective <- function(theta) {
    # Extract parameters
    parameters <- list(
      value = theta[1],
      slope = theta[2],
      covariates = if (n_cov > 0) theta[3:(2 + n_cov)] else numeric(0)
    )

    # Solve ODE and compute predictions
    predictions <- .solve_individual_longitudinal_ode(
      data = data,
      parameters = parameters
    )

    # Return large value if ODE solution failed
    if (is.null(predictions)) {
      return(.Machine$double.xmax)
    }

    # Compute sum of squared errors
    sum((data$response - predictions)^2)
  }

  # Perform optimization
  result <- optim(
    par = theta_init,
    fn = objective,
    method = "L-BFGS-B",
    control = control
  )

  # Extract optimized parameters
  value <- result$par[1]
  slope <- result$par[2]
  covariate_effects <- if (n_cov > 0) {
    effects <- result$par[3:(2 + n_cov)]
    names(effects) <- colnames(data$covariates)
    effects
  } else {
    numeric(0)
  }

  # Prepare final parameter list
  final_parameters <- list(
    value = value,
    slope = slope,
    covariates = covariate_effects
  )

  # Compute fitted values and residuals
  fitted_values <- .solve_individual_longitudinal_ode(
    data = data,
    parameters = final_parameters
  )
  residuals <- data$response - fitted_values

  # Return results
  list(
    subject_id = data$subject,
    coefficients = list(
      value = value,
      slope = slope,
      covariates = covariate_effects
    ),
    initial_state = data$initial,
    fitted_values = fitted_values,
    residuals = residuals,
    sse = result$value,
    convergence = result$convergence == 0,
    message = result$message
  )
}


# Solve longitudinal ODE for individual subject
.solve_individual_longitudinal_ode <- function(data, parameters) {
  # Extract data components
  times <- data$time
  X <- data$covariates
  initial_state <- data$initial

  # Extract parameters
  value <- parameters$value
  slope <- parameters$slope
  covariate_effects <- parameters$covariates

  # Define ODE derivative function
  ode_deriv <- function(t, state, parms) {
    m <- state[1] # biomarker
    v <- state[2] # velocity

    # Find nearest time point for covariate interpolation
    idx <- findInterval(t, times)
    idx <- max(1, min(idx, length(times)))

    # Compute excitation term from covariates
    excitation <- if (length(covariate_effects) > 0) {
      sum(X[idx, ] * covariate_effects)
    } else {
      0
    }

    # Compute acceleration: ẍ = value*x + slope*ẋ + excitation
    acceleration <- value * m + slope * v + excitation

    # Return derivatives [dx/dt, dv/dt]
    list(c(v, acceleration))
  }

  # Solve ODE system
  solution <- tryCatch(
    {
      suppressWarnings(
        deSolve::ode(
          y = initial_state,
          times = times,
          func = ode_deriv,
          parms = NULL
        )
      )
    },
    error = function(e) NULL
  )

  # Return biomarker values (column 2 of solution)
  if (is.null(solution)) NULL else solution[, 2]
}
