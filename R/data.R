#' Example Dataset for Joint ODE Model
#'
#' A simulated dataset generated using the Joint ODE Model framework,
#' demonstrating heterogeneous patient dynamics with longitudinal
#' biomarker measurements and survival outcomes.
#'
#' @format A list with two components:
#' \describe{
#'   \item{data}{A list containing the simulated data:
#'     \describe{
#'       \item{longitudinal_data}{Data frame with longitudinal measurements
#'         (938 observations):
#'         \itemize{
#'           \item id: Patient identifier (1-200)
#'           \item time: Measurement time
#'           \item observed: Observed biomarker (with error)
#'           \item biomarker: True biomarker value
#'           \item velocity: First derivative
#'           \item acceleration: Second derivative
#'           \item x1, x2: Longitudinal covariates
#'           \item xi: Patient-specific damping ratio
#'           \item period: Patient-specific natural period
#'         }
#'       }
#'       \item{survival_data}{Data frame with survival outcomes (200 patients):
#'         \itemize{
#'           \item id: Patient identifier
#'           \item time: Event or censoring time
#'           \item status: Event indicator (1=event, 0=censored)
#'           \item b: Shared random effect
#'           \item w1, w2: Survival covariates
#'         }
#'       }
#'       \item{state}{Data frame with initial states (200 patients):
#'         \itemize{
#'           \item biomarker: Initial biomarker value
#'           \item velocity: Initial velocity (0)
#'         }
#'       }
#'     }
#'   }
#'   \item{init}{A list containing initial parameter values:
#'     \describe{
#'       \item{coefficients}{True parameter values used in simulation:
#'         \itemize{
#'           \item baseline: B-spline coefficients for baseline hazard
#'           \item acceleration: ODE acceleration parameters
#'           \item hazard: Hazard function coefficients
#'           \item measurement_error_sd: Measurement error SD (0.15)
#'           \item random_effect_sd: Random effect SD (0.1)
#'         }
#'       }
#'       \item{configurations}{Model configuration:
#'         \itemize{
#'           \item baseline: B-spline configuration
#'           \item autonomous: TRUE (autonomous ODE system)
#'         }
#'       }
#'     }
#'   }
#' }
#'
#' @details
#' The dataset contains 200 patients with heterogeneous dynamics.
#' The patients are divided into four groups based on their xi and
#' period values:
#' \itemize{
#'   \item Underdamped (xi=0.3, period=5): oscillatory behavior
#'   \item Critically damped (xi=0.707, period=5): optimal damping
#'   \item Critically damped (xi=1.0, period=3): faster dynamics
#'   \item Overdamped (xi=2.0, period=3): exponential decay
#' }
#'
#' The "true" parameters stored in init$coefficients correspond to
#' the first group's dynamics for model fitting demonstration.
#'
#' @source Generated using \code{.create_example_data()} with seed 12345
#' @seealso \code{\link{JointODE}} for model fitting,
#'   \code{\link{simulate}} for data generation
#' @concept data-simulation
#' @examples
#' # Load the data
#' data(sim)
#'
#' # Examine structure
#' str(sim, max.level = 2)
#'
#' # Access longitudinal data
#' head(sim$data$longitudinal_data)
#'
#' # Summary of survival outcomes
#' summary(sim$data$survival_data$time)
#' table(sim$data$survival_data$status)
#'
#' # Visualize different dynamics behaviors
#' library(ggplot2)
#' library(dplyr)
#'
#' # Patient dynamics distribution
#' dynamics_summary <- sim$data$longitudinal_data %>%
#'   group_by(xi, period) %>%
#'   summarise(n_patients = n_distinct(id), .groups = "drop") %>%
#'   mutate(
#'     group = paste0("(ξ=", xi, ", T=", period, ")"),
#'     dynamics_type = case_when(
#'       xi < 1 ~ "Underdamped",
#'       xi == 1 ~ "Critically damped",
#'       xi > 1 ~ "Overdamped"
#'     )
#'   ) %>%
#'   select(group, dynamics_type, n_patients)
#' print(dynamics_summary)
#'
#' # Select example patients from each dynamics group
#' # Choose patients with at least 3 observations
#' dynamics_groups <- sim$data$longitudinal_data %>%
#'   group_by(id) %>%
#'   mutate(n_obs = n()) %>%
#'   filter(n_obs >= 3) %>%
#'   group_by(xi, period) %>%
#'   summarise(example_id = first(id), .groups = "drop")
#'
#' example_ids <- dynamics_groups$example_id
#'
#' # Plot trajectories for different dynamics
#' plot_data <- sim$data$longitudinal_data %>%
#'   filter(id %in% example_ids) %>%
#'   mutate(
#'     dynamics_type = case_when(
#'       xi < 1 ~ "Underdamped (Oscillatory)",
#'       xi == 1 ~ "Critically Damped",
#'       xi > 1 ~ "Overdamped (Exponential)"
#'     ),
#'     label = paste0("xi=", xi, ", T=", period)
#'   )
#'
#' # True trajectories
#' p1 <- ggplot(plot_data, aes(x = time, y = biomarker)) +
#'   geom_line(aes(color = dynamics_type), linewidth = 1.2) +
#'   geom_point(aes(y = observed), alpha = 0.4, size = 0.8) +
#'   facet_wrap(~ label, scales = "free_y", ncol = 2) +
#'   labs(
#'     title = "Biomarker Trajectories for Different Dynamics Groups",
#'     subtitle = "Lines: true values, Points: observed with measurement error",
#'     x = "Time",
#'     y = "Biomarker Value",
#'     color = "Dynamics Type"
#'   ) +
#'   theme_minimal() +
#'   theme(legend.position = "bottom")
#'
#' print(p1)
#'
#' # Phase plane plot (velocity vs biomarker)
#' p2 <- ggplot(plot_data, aes(x = biomarker, y = velocity)) +
#'   geom_path(aes(color = dynamics_type), linewidth = 1) +
#'   geom_point(data = filter(plot_data, time == 0),
#'              size = 3, shape = 21, fill = "white") +
#'   facet_wrap(~ label, scales = "free") +
#'   labs(
#'     title = "Phase Plane Trajectories",
#'     subtitle = "White circles indicate initial states",
#'     x = "Biomarker",
#'     y = "Velocity",
#'     color = "Dynamics Type"
#'   ) +
#'   theme_minimal() +
#'   theme(legend.position = "bottom")
#'
#' print(p2)
#'
#' \dontrun{
#' # Fit a Joint ODE model using this data
#' fit <- JointODE(
#'   longitudinal_formula = observed ~ x1 + x2,
#'   survival_formula = Surv(time, status) ~ w1 + w2,
#'   longitudinal_data = sim$data$longitudinal_data,
#'   survival_data = sim$data$survival_data,
#'   state = as.matrix(sim$data$state)
#' )
#' summary(fit)
#' }
"sim"
