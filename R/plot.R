# Declare global variables to avoid R CMD check warnings
# These variables are used in NSE contexts (ggplot2 aes(), dplyr verbs)
utils::globalVariables(c(
  "id",
  "biomarker",
  "velocity",
  "acceleration",
  "time",
  "survival",
  "group",
  "parameter",
  "residual",
  "value"
))

#' Plot Method for JointODE Objects
#'
#' @description
#' Provides comprehensive visualization tools for fitted JointODE models,
#' including longitudinal trajectories, phase space plots, survival curves,
#' diagnostic plots, and association patterns. All plots are generated using
#' ggplot2 for modern, publication-ready graphics.
#'
#' @param x An object of class \code{JointODE}
#' @param type Character string specifying the plot type:
#'   \itemize{
#'     \item \code{"overview"}: Panel of 4 key plots (survival, biomarker,
#'       velocity, phase space)
#'     \item \code{"trajectory_biomarker"}: Longitudinal biomarker trajectories.
#'       Can show individual curves or group averages based on \code{by}
#'     \item \code{"trajectory_velocity"}: Longitudinal velocity trajectories.
#'       Can show individual curves or group averages based on \code{by}
#'     \item \code{"phase_biomarker_velocity"}: Phase space plots
#'       (biomarker vs velocity). Can show individual curves or group
#'       averages based on \code{by}
#'     \item \code{"phase_velocity_acceleration"}: Phase space plots
#'       (velocity vs acceleration). Can show individual curves or group
#'       averages based on \code{by}
#'     \item \code{"survival"}: Survival probability curves.
#'       Can be stratified by \code{by} variable
#'     \item \code{"diagnostic_residuals"}: Residuals vs fitted values
#'     \item \code{"diagnostic_residuals_time"}: Residuals vs time
#'     \item \code{"diagnostic_qq"}: Normal Q-Q plot of residuals
#'     \item \code{"diagnostic_random_effects"}: Random effects distribution
#'     \item \code{"diagnostic_association"}: Association between
#'       biomarker features and hazard
#'   }
#' @param subject_ids Character or numeric vector of subject IDs to plot.
#'   If \code{NULL}, displays all subjects.
#'   Ignored when \code{by} is specified (shows group averages instead).
#' @param show_observed Logical; whether to show observed data points in
#'   faceted view when specific subject IDs are provided (default: TRUE).
#'   Ignored when showing all subjects (overlay view)
#' @param show_individual Logical; whether to show individual subject curves
#'   in addition to group means (default: TRUE). Applies to overview,
#'   trajectories, phase space, and survival plots. When FALSE, only shows
#'   group mean curves. Ignored for residuals, residuals_time, qq, and
#'   random_effects plots.
#' @param by Character string naming a covariate in the data to use for
#'   grouping.
#'   When specified for trajectories/phase/survival plots, shows
#'   group-averaged curves instead of individual subjects. If \code{NULL},
#'   shows overall population average (for trajectories/phase) or
#'   individual curves (default behavior).
#' @param n_groups Integer; number of groups to create when \code{by} is a
#'   continuous variable (default: 4). Subjects are divided into groups using
#'   quantiles. Ignored when \code{by} is categorical or \code{NULL}.
#' @param cols Color palette for plots (default: NULL uses built-in colors)
#' @param span Numeric value controlling the degree of smoothing for loess
#'   curves (default: 0.75). Larger values produce smoother curves. Only
#'   applies to mean curves in overview, survival, biomarker, and velocity
#'   plots.
#' @param ... Additional parameters passed to ggplot2 functions
#'
#' @return A ggplot2 object (or patchwork object for multi-panel plots) that
#'   can be further customized
#'
#' @details
#' The function provides several visualization types to assess model fit and
#' understand the joint modeling results:
#'
#' \strong{Overview plots} provide a quick diagnostic panel with four key
#' plots: survival probability (with individual curves), biomarker trajectory
#' (with individual curves), velocity trajectory (with individual curves),
#' and phase space plots.
#'
#' \strong{Trajectory plots} show observed biomarker values overlaid with
#' model-fitted trajectories, allowing assessment of longitudinal model fit
#' at the individual level.
#'
#' \strong{Dynamics plots} display phase portraits (biomarker vs velocity)
#' revealing the ODE system's behavior and stability properties.
#'
#' \strong{Survival plots} show predicted survival probabilities over time,
#' optionally stratified by covariates.
#'
#' \strong{Diagnostic plots} include residual plots and QQ plots for
#' checking model assumptions.
#'
#' \strong{Association plots} visualize how biomarker features (value and
#' slope) relate to hazard ratios.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_smooth
#'   geom_ribbon geom_abline geom_hline facet_wrap scale_color_manual
#'   scale_fill_manual scale_color_viridis_c theme_bw theme element_text
#'   element_blank element_rect labs coord_cartesian xlim ylim stat_qq
#'   stat_qq_line geom_density alpha geom_path geom_contour_filled
#'   guides guide_legend
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom dplyr group_by summarise mutate filter ungroup .data left_join
#' @importFrom tidyr pivot_longer
#' @importFrom stats loess predict residuals
#' @importFrom magrittr %>%
#' @importFrom viridisLite viridis
#'
#' @examples
#' \dontrun{
#' library(JointODE)
#' library(survival)
#'
#' # Load example dataset
#' data(sim)
#'
#' # Prepare data
#' longitudinal_data <- sim$data$longitudinal_data[
#'   , c("id", "time", "observed", "x1", "x2")
#' ]
#'
#' # Fit joint ODE model
#' fit <- JointODE(
#'   longitudinal_formula = observed ~ biomarker + velocity + x1 + x2 +
#'     (biomarker + velocity | id),
#'   survival_formula = Surv(time, status) ~ w1 + w2,
#'   longitudinal_data = longitudinal_data,
#'   survival_data = sim$data$survival_data,
#'   state = as.matrix(sim$data$state)
#' )
#'
#' # Overview plot
#' plot(fit)
#' plot(fit, type = "overview")
#'
#' # Individual trajectory plots
#' plot(fit, type = "trajectory_biomarker")  # Shows all subjects
#' plot(fit, type = "trajectory_biomarker", subject_ids = c("1", "2", "3"))
#' plot(fit, type = "trajectory_velocity")  # Velocity trajectories
#'
#' # Phase space plots
#' plot(fit, type = "phase_biomarker_velocity")
#' plot(fit, type = "phase_velocity_acceleration")
#'
#' # Group-stratified plots by survival covariates
#' plot(fit, type = "survival", by = "w1")  # Continuous variable (auto-grouped)
#' plot(fit, type = "trajectory_biomarker", by = "w2")
#' plot(fit, type = "trajectory_velocity", by = "w2")
#' plot(fit, type = "phase_biomarker_velocity", by = "w1")
#'
#' # Stratify survival by biomarker/velocity
#' plot(fit, type = "survival", by = "biomarker")
#' plot(fit, type = "survival", by = "velocity")
#'
#' # Diagnostic plots
#' plot(fit, type = "diagnostic_residuals")
#' plot(fit, type = "diagnostic_residuals_time")
#' plot(fit, type = "diagnostic_qq")
#' plot(fit, type = "diagnostic_random_effects")
#'
#' # Customize colors
#' plot(fit, type = "survival", cols = c("red", "blue", "green"))
#'
#' # Adjust smoothing
#' plot(fit, type = "overview", span = 0.5)  # Less smooth
#' plot(fit, type = "survival", span = 1.0)  # More smooth
#' }
#'
#' @concept model-display
#' @export
plot.JointODE <- function(
  x,
  type = c(
    "overview",
    "trajectory_biomarker",
    "trajectory_velocity",
    "phase_biomarker_velocity",
    "phase_velocity_acceleration",
    "survival",
    "diagnostic_residuals",
    "diagnostic_residuals_time",
    "diagnostic_qq",
    "diagnostic_random_effects",
    "diagnostic_association"
  ),
  subject_ids = NULL,
  show_observed = TRUE,
  show_individual = TRUE,
  by = NULL,
  n_groups = 4,
  cols = NULL,
  span = 0.75,
  ...
) {
  type <- match.arg(type)

  # Route to appropriate plot function
  plot_obj <- switch(
    type,
    overview = .plot_overview(
      x,
      subject_ids = subject_ids,
      show_individual = show_individual,
      cols = cols,
      span = span,
      ...
    ),
    trajectory_biomarker = .plot_biomarker(
      x,
      subject_ids = subject_ids,
      show_observed = show_observed,
      show_individual = show_individual,
      by = by,
      times = NULL,
      cols = cols,
      span = span,
      ...
    ),
    trajectory_velocity = .plot_velocity(
      x,
      subject_ids = subject_ids,
      show_observed = FALSE,
      show_individual = show_individual,
      by = by,
      times = NULL,
      cols = cols,
      span = span,
      ...
    ),
    phase_biomarker_velocity = .plot_phase(
      x,
      subject_ids = subject_ids,
      show_individual = show_individual,
      by = by,
      times = NULL,
      cols = cols,
      span = span,
      ...
    ),
    phase_velocity_acceleration = .plot_vel_acc(
      x,
      subject_ids = subject_ids,
      show_individual = show_individual,
      by = by,
      times = NULL,
      cols = cols,
      span = span,
      ...
    ),
    survival = .plot_survival(
      x,
      by = by,
      n_groups = n_groups,
      show_individual = show_individual,
      times = NULL,
      cols = cols,
      span = span,
      ...
    ),
    diagnostic_residuals = .plot_residuals_vs_fitted(
      x,
      cols = cols,
      span = span,
      ...
    ),
    diagnostic_residuals_time = .plot_residuals_vs_time(
      x,
      cols = cols,
      span = span,
      ...
    ),
    diagnostic_qq = .plot_qq(
      x,
      cols = cols,
      ...
    ),
    diagnostic_random_effects = .plot_random_effects(
      x,
      cols = cols,
      ...
    ),
    diagnostic_association = stop(
      "Association plot not yet implemented. ",
      "This plot type is reserved for future functionality."
    )
  )

  print(plot_obj)
  invisible(plot_obj)
}


# ==============================================================================
# 2. OVERVIEW PANEL
# ==============================================================================

.plot_overview <- function(
  x,
  subject_ids = NULL,
  show_individual = TRUE,
  cols = NULL,
  span = 0.75,
  ...
) {
  # Plot 1: Survival Probability
  p1 <- .plot_survival(
    x = x,
    subject_ids = subject_ids,
    show_individual = show_individual,
    by = NULL,
    times = NULL,
    cols = cols,
    span = span,
    ...
  ) +
    labs(title = "Survival Probability") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom"
    )

  # Plot 2: Biomarker Trajectory
  p2 <- .plot_biomarker(
    x = x,
    subject_ids = subject_ids,
    show_observed = FALSE,
    show_individual = show_individual,
    by = NULL,
    times = NULL,
    cols = cols,
    span = span,
    ...
  ) +
    labs(title = "Biomarker") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom"
    )

  # Plot 3: Velocity Trajectory
  p3 <- .plot_velocity(
    x = x,
    subject_ids = subject_ids,
    show_observed = FALSE,
    show_individual = show_individual,
    by = NULL,
    times = NULL,
    cols = cols,
    span = span,
    ...
  ) +
    labs(title = "Velocity") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom"
    )

  # Plot 4: Phase Space
  p4 <- .plot_phase(
    x = x,
    subject_ids = subject_ids,
    show_individual = show_individual,
    by = NULL,
    times = NULL,
    cols = cols,
    span = span,
    ...
  ) +
    labs(title = "Phase Space") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom"
    )

  # Combine plots in 2x2 layout (without collecting guides)
  combined_plot <- (p1 + p2 + p3 + p4) +
    patchwork::plot_layout(ncol = 2)

  return(combined_plot)
}


# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

# Create time grid for predictions
.get_time_grid <- function(x, times = NULL, n_points = 100) {
  if (!is.null(times)) {
    return(times)
  }

  all_times <- unlist(lapply(x$data, function(subj) {
    subj$longitudinal$times
  }))

  seq(min(all_times), max(all_times), length.out = n_points)
}

# Standard theme for individual subject faceted plots
.theme_faceted <- function() {
  theme_bw() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9)
    )
}

# Standard theme for grouped plots
.theme_grouped <- function() {
  theme_bw() +
    theme(
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    )
}

# Standard theme for simple plots (no legend)
.theme_simple <- function() {
  theme_bw() +
    theme(
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9)
    )
}

# Build faceted plot for user-specified subjects
.build_faceted_trajectory_plot <- function(
  data,
  x_var,
  y_var,
  x_label,
  y_label,
  color = "#2E86AB"
) {
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], group = id)) +
    geom_line(color = color, linewidth = 0.8) +
    facet_wrap(~id, scales = "free") +
    labs(x = x_label, y = y_label) +
    .theme_faceted()
}

# Build overlay plot with individual curves and mean
.build_overlay_plot <- function(
  data,
  x_var,
  y_var,
  x_label,
  y_label,
  show_individual = TRUE,
  span = 0.75
) {
  p <- ggplot()

  if (show_individual) {
    p <- p +
      geom_line(
        data = data,
        aes(
          x = .data[[x_var]],
          y = .data[[y_var]],
          group = id,
          color = "Individual"
        ),
        alpha = 0.15,
        linewidth = 0.4
      ) +
      geom_smooth(
        data = data,
        aes(x = .data[[x_var]], y = .data[[y_var]], color = "Mean"),
        method = "loess",
        span = span,
        se = FALSE,
        linewidth = 1.2
      ) +
      scale_color_manual(
        name = NULL,
        values = c("Individual" = "gray60", "Mean" = "#2E86AB"),
        breaks = c("Mean", "Individual")
      )
  } else {
    p <- p +
      geom_smooth(
        data = data,
        aes(x = .data[[x_var]], y = .data[[y_var]]),
        method = "loess",
        span = span,
        se = FALSE,
        color = "#2E86AB",
        linewidth = 1.2
      )
  }

  p <- p + labs(x = x_label, y = y_label) + .theme_simple()

  if (show_individual) {
    p <- p + theme(legend.position = "bottom")
  }

  p
}

# Build grouped trajectory plot (for biomarker/velocity)
.build_grouped_trajectory_plot <- function(
  pred_data,
  group_means,
  x_var,
  y_var,
  x_label,
  y_label,
  group_colors,
  group_name,
  show_individual = TRUE,
  span = 0.75
) {
  p <- ggplot()

  if (show_individual) {
    p <- p +
      geom_line(
        data = pred_data,
        aes(x = .data[[x_var]], y = .data[[y_var]], group = id, color = group),
        alpha = 0.1,
        linewidth = 0.3
      )
  }

  p <- p +
    geom_smooth(
      data = group_means,
      aes(x = .data[[x_var]], y = .data[[y_var]], color = group),
      method = "loess",
      span = span,
      se = FALSE,
      linewidth = 1.5
    ) +
    scale_color_manual(values = group_colors, name = group_name) +
    labs(x = x_label, y = y_label) +
    .theme_grouped()

  return(p)
}

# Build grouped phase space plot (no smoothing, uses geom_path)
.build_grouped_phase_plot <- function(
  pred_data,
  group_means,
  group_colors,
  group_name,
  show_individual = TRUE
) {
  p <- ggplot()

  if (show_individual) {
    p <- p +
      geom_path(
        data = pred_data,
        aes(x = biomarker, y = velocity, group = id, color = group),
        alpha = 0.1,
        linewidth = 0.3
      )
  }

  p <- p +
    geom_path(
      data = group_means,
      aes(x = biomarker, y = velocity, color = group),
      linewidth = 1.5
    ) +
    scale_color_manual(values = group_colors, name = group_name) +
    labs(x = "Biomarker", y = "Velocity") +
    .theme_grouped()

  return(p)
}

# Build grouped velocity-acceleration plot (no smoothing, uses geom_path)
.build_grouped_vel_acc_plot <- function(
  pred_data,
  group_means,
  group_colors,
  group_name,
  show_individual = TRUE
) {
  p <- ggplot()

  if (show_individual) {
    p <- p +
      geom_path(
        data = pred_data,
        aes(x = velocity, y = acceleration, group = id, color = group),
        alpha = 0.1,
        linewidth = 0.3
      )
  }

  p <- p +
    geom_path(
      data = group_means,
      aes(x = velocity, y = acceleration, color = group),
      linewidth = 1.5
    ) +
    scale_color_manual(values = group_colors, name = group_name) +
    labs(x = "Velocity", y = "Acceleration") +
    .theme_grouped()

  return(p)
}

# Build overlay phase space plot with individual curves and mean
.build_overlay_phase_plot <- function(
  data,
  show_individual = TRUE
) {
  p <- ggplot()

  if (show_individual) {
    p <- p +
      geom_path(
        data = data,
        aes(x = biomarker, y = velocity, group = id, color = "Individual"),
        alpha = 0.15,
        linewidth = 0.4
      )

    # Calculate mean phase space trajectory
    mean_data <- data %>%
      group_by(time) %>%
      summarise(
        biomarker = mean(biomarker, na.rm = TRUE),
        velocity = mean(velocity, na.rm = TRUE),
        .groups = "drop"
      )

    p <- p +
      geom_path(
        data = mean_data,
        aes(x = biomarker, y = velocity, color = "Mean"),
        linewidth = 1.2
      ) +
      scale_color_manual(
        name = NULL,
        values = c("Individual" = "gray60", "Mean" = "#2E86AB"),
        breaks = c("Mean", "Individual")
      )
  } else {
    # Only show mean phase space trajectory
    mean_data <- data %>%
      group_by(time) %>%
      summarise(
        biomarker = mean(biomarker, na.rm = TRUE),
        velocity = mean(velocity, na.rm = TRUE),
        .groups = "drop"
      )

    p <- p +
      geom_path(
        data = mean_data,
        aes(x = biomarker, y = velocity),
        color = "#2E86AB",
        linewidth = 1.2
      )
  }

  p <- p + labs(x = "Biomarker", y = "Velocity") + .theme_simple()

  if (show_individual) {
    p <- p + theme(legend.position = "bottom")
  }

  p
}

# Build overlay velocity-acceleration plot with individual curves and mean
.build_overlay_vel_acc_plot <- function(
  data,
  show_individual = TRUE
) {
  p <- ggplot()

  if (show_individual) {
    p <- p +
      geom_path(
        data = data,
        aes(x = velocity, y = acceleration, group = id, color = "Individual"),
        alpha = 0.15,
        linewidth = 0.4
      )

    # Calculate mean velocity-acceleration trajectory
    mean_data <- data %>%
      group_by(time) %>%
      summarise(
        velocity = mean(velocity, na.rm = TRUE),
        acceleration = mean(acceleration, na.rm = TRUE),
        .groups = "drop"
      )

    p <- p +
      geom_path(
        data = mean_data,
        aes(x = velocity, y = acceleration, color = "Mean"),
        linewidth = 1.2
      ) +
      scale_color_manual(
        name = NULL,
        values = c("Individual" = "gray60", "Mean" = "#2E86AB"),
        breaks = c("Mean", "Individual")
      )
  } else {
    # Only show mean velocity-acceleration trajectory
    mean_data <- data %>%
      group_by(time) %>%
      summarise(
        velocity = mean(velocity, na.rm = TRUE),
        acceleration = mean(acceleration, na.rm = TRUE),
        .groups = "drop"
      )

    p <- p +
      geom_path(
        data = mean_data,
        aes(x = velocity, y = acceleration),
        color = "#2E86AB",
        linewidth = 1.2
      )
  }

  p <- p + labs(x = "Velocity", y = "Acceleration") + .theme_simple()

  if (show_individual) {
    p <- p + theme(legend.position = "bottom")
  }

  p
}

# Prepare prediction data with group information
.prepare_grouped_data <- function(x, by, times = NULL, n_groups = 4) {
  # Extract groups
  groups <- .extract_groups(x, by, n_groups = n_groups)
  group_colors <- .setup_group_colors(groups, cols = NULL)

  # Get predictions
  times <- .get_time_grid(x, times)
  pred_data <- predict(x, times = times)

  # Add group information
  pred_data$group <- groups[as.character(pred_data$id)]

  list(
    pred_data = pred_data,
    groups = groups,
    group_colors = group_colors
  )
}

# Build grouped survival plot
.build_grouped_survival_plot <- function(
  pred_data,
  group_colors,
  group_name,
  show_individual = TRUE,
  span = 0.75
) {
  p <- ggplot()

  if (show_individual) {
    p <- p +
      geom_line(
        data = pred_data,
        aes(
          x = time,
          y = survival,
          group = interaction(id, group),
          color = group
        ),
        alpha = 0.1,
        linewidth = 0.3
      )
  }

  p <- p +
    geom_smooth(
      data = pred_data,
      aes(x = time, y = survival, color = group, group = group),
      method = "loess",
      span = span,
      se = FALSE,
      linewidth = 1.5
    ) +
    scale_color_manual(values = group_colors, name = group_name) +
    labs(x = "Time", y = "Survival Probability") +
    .theme_grouped() +
    coord_cartesian(ylim = c(0, 1))

  return(p)
}

# Build residual plot (shared by residuals vs fitted and residuals vs time)
.build_residual_plot <- function(resid_data, x_var, x_label, span = 0.75) {
  ggplot(resid_data, aes(x = .data[[x_var]], y = residual)) +
    geom_point(color = "#2E86AB", alpha = 0.5, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(
      method = "loess",
      span = span,
      se = TRUE,
      color = "#A23B72",
      fill = "#A23B72",
      alpha = 0.2
    ) +
    labs(x = x_label, y = "Residuals") +
    .theme_simple()
}

.extract_groups <- function(x, by, n_groups = 4) {
  # Get survival covariate names from hazard coefficients
  # The first two are alpha_1 and alpha_2, rest are covariate names
  hazard_names <- names(x$parameters$coefficients$hazard)
  survival_cov_names <- hazard_names[-(1:2)] # Remove alpha_1 and alpha_2

  # Track if we're using individual means for time-varying covariates
  using_individual_means <- FALSE

  group_var <- vapply(x$data, function(subj) {
    # Check if 'by' is in survival covariates
    if (by %in% survival_cov_names) {
      idx <- which(survival_cov_names == by)
      if (idx <= length(subj$covariates)) {
        return(subj$covariates[idx])
      }
    }

    # Check if 'by' is in longitudinal covariates
    if (!is.null(subj$longitudinal$covariates)) {
      # Check fixed covariates
      if (!is.null(subj$longitudinal$covariates$fixed)) {
        fixed_cov_names <- colnames(subj$longitudinal$covariates$fixed)
        if (by %in% fixed_cov_names) {
          # Use mean value across all time points for time-varying covariates
          using_individual_means <<- TRUE
          return(mean(subj$longitudinal$covariates$fixed[, by], na.rm = TRUE))
        }
      }
      # Check random covariates
      if (!is.null(subj$longitudinal$covariates$random)) {
        random_cov_names <- colnames(subj$longitudinal$covariates$random)
        if (by %in% random_cov_names) {
          # Use mean value across all time points for time-varying covariates
          using_individual_means <<- TRUE
          return(mean(subj$longitudinal$covariates$random[, by], na.rm = TRUE))
        }
      }
    }

    # Check if 'by' is directly available in subject data
    if (!is.null(subj[[by]])) {
      return(subj[[by]])
    }

    NA
  }, numeric(1))

  # Inform user if grouping by individual means
  if (using_individual_means) {
    cli::cli_alert_info(
      "Grouping by individual means of time-varying covariate '{by}'."
    )
  }

  # Check if all values are NA
  if (all(is.na(group_var))) {
    long_cov_example <- "none"
    if (!is.null(x$data[[1]]$longitudinal$covariates)) {
      fixed_names <- if (!is.null(x$data[[1]]$longitudinal$covariates$fixed)) {
        colnames(x$data[[1]]$longitudinal$covariates$fixed)
      } else {
        character(0)
      }
      random_names <- if (
        !is.null(x$data[[1]]$longitudinal$covariates$random)
      ) {
        colnames(x$data[[1]]$longitudinal$covariates$random)
      } else {
        character(0)
      }
      all_long_names <- c(fixed_names, random_names)
      if (length(all_long_names) > 0) {
        long_cov_example <- paste(all_long_names, collapse = ", ")
      }
    }

    stop(
      sprintf(
        paste0(
          "Grouping variable '%s' not found in covariates.\n",
          "Available survival covariates: %s\n",
          "Available longitudinal covariates: %s"
        ),
        by,
        paste(survival_cov_names, collapse = ", "),
        long_cov_example
      )
    )
  }

  # Check if continuous variable - if so, discretize into groups
  unique_vals <- unique(group_var[!is.na(group_var)])
  n_unique <- length(unique_vals)

  # If more than 5 unique values, treat as continuous and create groups
  if (n_unique > 5) {
    # Use quantile-based grouping with n_groups
    probs <- seq(0, 1, length.out = n_groups + 1)
    quantile_breaks <- quantile(group_var, probs = probs, na.rm = TRUE)

    # Generate labels for groups
    group_labels <- vapply(seq_len(n_groups), function(i) {
      sprintf("G%d (%.0f-%.0f%%)", i, probs[i] * 100, probs[i + 1] * 100)
    }, character(1))

    # Assign each observation to a quantile group
    group_var_discrete <- cut(
      group_var,
      breaks = quantile_breaks,
      labels = group_labels,
      include.lowest = TRUE
    )
    group_var_discrete <- as.character(group_var_discrete)

    cli::cli_alert_info(
      paste0(
        "'{by}' is continuous (n={n_unique} unique values). ",
        "Created {n_groups} groups using quantiles."
      )
    )

    names(group_var_discrete) <- names(x$data)
    return(group_var_discrete)
  }

  result <- as.character(group_var)
  names(result) <- names(x$data)
  return(result)
}

.setup_group_colors <- function(groups, cols = NULL) {
  unique_groups <- unique(groups[!is.na(groups)])
  n_groups <- length(unique_groups)

  if (is.null(cols)) {
    # Use viridis color palette for any number of groups
    cols <- viridisLite::viridis(n_groups)
  } else {
    if (length(cols) < n_groups) {
      cols <- rep(cols, length.out = n_groups)
    }
  }

  names(cols) <- unique_groups
  return(cols)
}

.calculate_residuals <- function(x) {
  # Get predictions at observed times
  pred_data <- predict(x, times = NULL)

  residuals_list <- list()

  for (i in seq_along(x$data)) {
    subj <- x$data[[i]]
    sid <- names(x$data)[i]
    obs_times <- subj$longitudinal$times
    obs_values <- subj$longitudinal$measurements

    # Get fitted values from predictions
    subj_pred <- pred_data[pred_data$id == sid, ]

    # Match observed times to predicted times
    matched_idx <- match(obs_times, subj_pred$time)

    if (any(!is.na(matched_idx))) {
      valid_idx <- !is.na(matched_idx)
      fitted_vals <- subj_pred$biomarker[matched_idx[valid_idx]]
      observed_vals <- obs_values[valid_idx]
      times_vals <- obs_times[valid_idx]

      residuals_list[[i]] <- data.frame(
        id = sid,
        time = times_vals,
        observed = observed_vals,
        fitted = fitted_vals,
        residual = observed_vals - fitted_vals
      )
    }
  }

  do.call(rbind, residuals_list)
}

.select_subjects <- function(x, subject_ids = NULL) {
  all_ids <- names(x$data)

  if (is.null(subject_ids)) {
    subject_ids <- all_ids
  }

  return(subject_ids)
}

.calculate_group_means <- function(pred_data, groups, type = "trajectory") {
  # Add groups to prediction data
  pred_data$group <- rep(groups, each = length(unique(pred_data$time)))

  # Calculate group means
  if (type == "trajectory") {
    group_means <- pred_data %>%
      group_by(group, time) %>%
      summarise(
        biomarker = mean(biomarker, na.rm = TRUE),
        .groups = "drop"
      )
  } else if (type == "velocity") {
    group_means <- pred_data %>%
      group_by(group, time) %>%
      summarise(
        velocity = mean(velocity, na.rm = TRUE),
        .groups = "drop"
      )
  } else if (type == "phase") {
    group_means <- pred_data %>%
      group_by(group, time) %>%
      summarise(
        biomarker = mean(biomarker, na.rm = TRUE),
        velocity = mean(velocity, na.rm = TRUE),
        .groups = "drop"
      )
  } else if (type == "velocity_acceleration") {
    group_means <- pred_data %>%
      group_by(group, time) %>%
      summarise(
        velocity = mean(velocity, na.rm = TRUE),
        acceleration = mean(acceleration, na.rm = TRUE),
        .groups = "drop"
      )
  } else if (type == "survival") {
    group_means <- pred_data %>%
      group_by(group, time) %>%
      summarise(
        survival = mean(survival, na.rm = TRUE),
        .groups = "drop"
      )
  }

  return(group_means)
}


# ==============================================================================
# 4. TRAJECTORY PLOTS
# ==============================================================================

.plot_biomarker <- function(
  x,
  subject_ids = NULL,
  show_observed = TRUE,
  show_individual = TRUE,
  by = NULL,
  times = NULL,
  cols = NULL,
  span = 0.75,
  ...
) {
  if (!is.null(by)) {
    return(.plot_biomarker_grouped(
      x,
      by,
      show_individual,
      times,
      cols,
      span,
      ...
    ))
  }

  user_specified_subjects <- !is.null(subject_ids)
  times <- .get_time_grid(x, times)
  pred_data <- predict(x, times = times)

  if (!is.data.frame(pred_data) || nrow(pred_data) == 0) {
    stop("predict() returned invalid or empty data")
  }

  subject_ids <- .select_subjects(x, subject_ids)
  pred_subset <- pred_data[pred_data$id %in% subject_ids, ]

  if (user_specified_subjects) {
    p <- ggplot() +
      geom_line(
        data = pred_subset,
        aes(x = time, y = biomarker, group = id),
        color = "#2E86AB",
        linewidth = 0.8
      )

    if (show_observed) {
      obs_data_list <- lapply(subject_ids, function(sid) {
        subj <- x$data[[sid]]
        data.frame(
          id = sid,
          time = subj$longitudinal$times,
          biomarker = subj$longitudinal$measurements
        )
      })

      if (length(obs_data_list) > 0) {
        obs_data <- do.call(rbind, obs_data_list)
        p <- p +
          geom_point(
            data = obs_data,
            aes(x = time, y = biomarker, group = id),
            color = "#A23B72",
            size = 1.5,
            alpha = 0.6
          )
      }
    }

    p <- p +
      facet_wrap(~id, scales = "free") +
      labs(x = "Time", y = "Biomarker") +
      .theme_faceted()
  } else {
    # Overlay view: show_observed is ignored (only for faceted view)
    p <- .build_overlay_plot(
      data = pred_subset,
      x_var = "time",
      y_var = "biomarker",
      x_label = "Time",
      y_label = "Biomarker",
      show_individual = show_individual,
      span = span
    )
  }

  return(p)
}

.plot_biomarker_grouped <- function(
  x,
  by,
  show_individual = TRUE,
  times = NULL,
  cols = NULL,
  span = 0.75,
  ...
) {
  grouped <- .prepare_grouped_data(x, by, times)
  if (!is.null(cols)) {
    grouped$group_colors <- .setup_group_colors(grouped$groups, cols)
  }

  group_means <- .calculate_group_means(
    grouped$pred_data,
    grouped$groups,
    type = "trajectory"
  )

  .build_grouped_trajectory_plot(
    pred_data = grouped$pred_data,
    group_means = group_means,
    x_var = "time",
    y_var = "biomarker",
    x_label = "Time",
    y_label = "Biomarker",
    group_colors = grouped$group_colors,
    group_name = by,
    show_individual = show_individual,
    span = span
  )
}

.plot_velocity <- function(
  x,
  subject_ids = NULL,
  show_observed = FALSE,
  show_individual = TRUE,
  by = NULL,
  times = NULL,
  cols = NULL,
  span = 0.75,
  ...
) {
  if (!is.null(by)) {
    return(.plot_velocity_grouped(
      x,
      by,
      show_individual,
      times,
      cols,
      span,
      ...
    ))
  }

  user_specified_subjects <- !is.null(subject_ids)
  times <- .get_time_grid(x, times)
  pred_data <- predict(x, times = times)

  if (!is.data.frame(pred_data) || nrow(pred_data) == 0) {
    stop("predict() returned invalid or empty data")
  }

  subject_ids <- .select_subjects(x, subject_ids)
  pred_subset <- pred_data[pred_data$id %in% subject_ids, ]

  if (user_specified_subjects) {
    return(.build_faceted_trajectory_plot(
      data = pred_subset,
      x_var = "time",
      y_var = "velocity",
      x_label = "Time",
      y_label = "Velocity"
    ))
  } else {
    return(.build_overlay_plot(
      data = pred_subset,
      x_var = "time",
      y_var = "velocity",
      x_label = "Time",
      y_label = "Velocity",
      show_individual = show_individual,
      span = span
    ))
  }
}

.plot_velocity_grouped <- function(
  x,
  by,
  show_individual = TRUE,
  times = NULL,
  cols = NULL,
  span = 0.75,
  ...
) {
  # Prepare grouped data
  grouped <- .prepare_grouped_data(x, by, times)
  if (!is.null(cols)) {
    grouped$group_colors <- .setup_group_colors(grouped$groups, cols)
  }

  # Calculate group means
  group_means <- .calculate_group_means(
    grouped$pred_data,
    grouped$groups,
    type = "velocity"
  )

  # Build plot
  .build_grouped_trajectory_plot(
    pred_data = grouped$pred_data,
    group_means = group_means,
    x_var = "time",
    y_var = "velocity",
    x_label = "Time",
    y_label = "Velocity",
    group_colors = grouped$group_colors,
    group_name = by,
    show_individual = show_individual,
    span = span
  )
}


# ==============================================================================
# 5. PHASE SPACE PLOTS
# ==============================================================================

.plot_phase <- function(
  x,
  subject_ids = NULL,
  show_individual = TRUE,
  by = NULL,
  times = NULL,
  cols = NULL,
  span = 0.75,
  ...
) {
  if (!is.null(by)) {
    return(.plot_phase_grouped(
      x,
      by = by,
      show_individual = show_individual,
      times = times,
      cols = cols,
      span = span,
      ...
    ))
  }

  # Remember if user explicitly specified subject_ids
  user_specified_subjects <- !is.null(subject_ids)

  # Get predictions
  times <- .get_time_grid(x, times)
  pred_data <- predict(x, times = times)

  # Select subjects
  subject_ids <- .select_subjects(x, subject_ids)
  pred_data <- pred_data[pred_data$id %in% subject_ids, ]

  # If user specified subjects explicitly, use faceting for detailed view
  if (user_specified_subjects) {
    return(
      ggplot(pred_data, aes(x = biomarker, y = velocity, group = id)) +
        geom_path(color = "#2E86AB", linewidth = 0.8) +
        facet_wrap(~id, scales = "free") +
        labs(x = "Biomarker", y = "Velocity") +
        .theme_faceted()
    )
  }

  # If showing all subjects, use overlay with mean trajectory
  .build_overlay_phase_plot(data = pred_data, show_individual = show_individual)
}

.plot_phase_grouped <- function(
  x,
  by,
  show_individual = TRUE,
  times = NULL,
  cols = NULL,
  span = 0.75,
  ...
) {
  # Prepare grouped data
  grouped <- .prepare_grouped_data(x, by, times)
  if (!is.null(cols)) {
    grouped$group_colors <- .setup_group_colors(grouped$groups, cols)
  }

  # Calculate group means
  group_means <- .calculate_group_means(
    grouped$pred_data,
    grouped$groups,
    type = "phase"
  )

  # Build plot
  .build_grouped_phase_plot(
    pred_data = grouped$pred_data,
    group_means = group_means,
    group_colors = grouped$group_colors,
    group_name = by,
    show_individual = show_individual
  )
}


# ==============================================================================
# 6. VELOCITY-ACCELERATION PHASE SPACE PLOTS
# ==============================================================================

.plot_vel_acc <- function(
  x,
  subject_ids = NULL,
  show_individual = TRUE,
  by = NULL,
  times = NULL,
  cols = NULL,
  span = 0.75,
  ...
) {
  if (!is.null(by)) {
    return(.plot_vel_acc_grouped(
      x,
      by = by,
      show_individual = show_individual,
      times = times,
      cols = cols,
      span = span,
      ...
    ))
  }

  # Remember if user explicitly specified subject_ids
  user_specified_subjects <- !is.null(subject_ids)

  # Get predictions
  times <- .get_time_grid(x, times)
  pred_data <- predict(x, times = times)

  # Select subjects
  subject_ids <- .select_subjects(x, subject_ids)
  pred_data <- pred_data[pred_data$id %in% subject_ids, ]

  # If user specified subjects explicitly, use faceting for detailed view
  if (user_specified_subjects) {
    return(
      ggplot(pred_data, aes(x = velocity, y = acceleration, group = id)) +
        geom_path(color = "#2E86AB", linewidth = 0.8) +
        facet_wrap(~id, scales = "free") +
        labs(x = "Velocity", y = "Acceleration") +
        .theme_faceted()
    )
  }

  # If showing all subjects, use overlay with mean trajectory
  .build_overlay_vel_acc_plot(
    data = pred_data,
    show_individual = show_individual
  )
}

.plot_vel_acc_grouped <- function(
  x,
  by,
  show_individual = TRUE,
  times = NULL,
  cols = NULL,
  span = 0.75,
  ...
) {
  # Prepare grouped data
  grouped <- .prepare_grouped_data(x, by, times)
  if (!is.null(cols)) {
    grouped$group_colors <- .setup_group_colors(grouped$groups, cols)
  }

  # Calculate group means
  group_means <- .calculate_group_means(
    grouped$pred_data,
    grouped$groups,
    type = "velocity_acceleration"
  )

  # Build plot
  .build_grouped_vel_acc_plot(
    pred_data = grouped$pred_data,
    group_means = group_means,
    group_colors = grouped$group_colors,
    group_name = by,
    show_individual = show_individual
  )
}


# ==============================================================================
# 6. SURVIVAL PLOTS
# ==============================================================================

.plot_survival <- function(
  x,
  subject_ids = NULL,
  by = NULL,
  n_groups = 4,
  show_individual = TRUE,
  times = NULL,
  cols = NULL,
  span = 0.75,
  ...
) {
  # Remember if user explicitly specified subject_ids
  user_specified_subjects <- !is.null(subject_ids)

  # Get predictions
  times <- .get_time_grid(x, times)
  pred_data <- predict(x, times = times)

  # Select subjects
  subject_ids <- .select_subjects(x, subject_ids)
  pred_data <- pred_data[pred_data$id %in% subject_ids, ]

  # Case 1: No grouping - overall survival curve
  if (is.null(by)) {
    if (user_specified_subjects) {
      # Faceted view for specific subjects
      return(
        ggplot(pred_data, aes(x = time, y = survival, group = id)) +
          geom_line(color = "black", linewidth = 0.8) +
          facet_wrap(~id, scales = "free_y") +
          labs(x = "Time", y = "Survival Probability") +
          .theme_faceted() +
          theme(panel.grid.minor = element_blank()) +
          coord_cartesian(ylim = c(0, 1))
      )
    }

    # Overall survival with individual curves and mean
    p <- ggplot()
    if (show_individual) {
      p <- p +
        geom_line(
          data = pred_data,
          aes(x = time, y = survival, group = id, color = "Individual"),
          alpha = 0.15,
          linewidth = 0.4
        ) +
        geom_smooth(
          data = pred_data,
          aes(x = time, y = survival, color = "Mean"),
          method = "loess",
          span = span,
          se = FALSE,
          linewidth = 1.2
        ) +
        scale_color_manual(
          name = NULL,
          values = c("Individual" = "gray60", "Mean" = "#2E86AB"),
          breaks = c("Mean", "Individual")
        )
    } else {
      p <- p +
        geom_smooth(
          data = pred_data,
          aes(x = time, y = survival),
          method = "loess",
          span = span,
          se = FALSE,
          color = "#2E86AB",
          linewidth = 1.2
        )
    }

    return(
      p +
        labs(x = "Time", y = "Survival Probability") +
        .theme_simple() +
        theme(
          legend.position = "bottom",
          panel.grid.minor = element_blank()
        ) +
        coord_cartesian(ylim = c(0, 1))
    )
  }

  # Case 2: Grouping by biomarker or velocity (use individual means)
  if (tolower(by) %in% c("biomarker", "velocity")) {
    var_name <- tolower(by)
    cli::cli_alert_info("Grouping by individual means of '{by}'.")

    # Calculate per-individual means and create temporary object
    individual_means <- pred_data %>%
      group_by(id) %>%
      summarise(
        individual_mean = mean(.data[[var_name]], na.rm = TRUE),
        .groups = "drop"
      )

    x_temp <- x
    for (i in seq_along(x_temp$data)) {
      id <- names(x_temp$data)[i]
      mean_val <- individual_means$individual_mean[individual_means$id == id]
      x_temp$data[[i]][[var_name]] <- mean_val
    }

    groups <- .extract_groups(x_temp, by = var_name, n_groups = n_groups)
    group_colors <- .setup_group_colors(groups, cols)
    pred_data$group <- rep(groups, each = length(unique(pred_data$time)))

    return(.build_grouped_survival_plot(
      pred_data,
      group_colors,
      by,
      show_individual,
      span
    ))
  }

  # Case 3: Grouping by other covariates
  groups <- .extract_groups(x, by, n_groups = n_groups)
  group_colors <- .setup_group_colors(groups, cols)
  pred_data$group <- rep(groups, each = length(unique(pred_data$time)))

  .build_grouped_survival_plot(
    pred_data,
    group_colors,
    by,
    show_individual,
    span
  )
}


# ==============================================================================
# 7. DIAGNOSTIC PLOTS
# ==============================================================================

.plot_residuals_vs_fitted <- function(x, cols = NULL, span = 0.75, ...) {
  .build_residual_plot(.calculate_residuals(x), "fitted", "Fitted Values", span)
}

.plot_residuals_vs_time <- function(x, cols = NULL, span = 0.75, ...) {
  .build_residual_plot(.calculate_residuals(x), "time", "Time", span)
}

.plot_qq <- function(x, cols = NULL, ...) {
  resid_data <- .calculate_residuals(x)

  ggplot(resid_data, aes(sample = residual)) +
    stat_qq(color = "#2E86AB", alpha = 0.6, size = 2) +
    stat_qq_line(color = "red", linetype = "dashed") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    .theme_simple()
}

.plot_random_effects <- function(x, cols = NULL, ...) {
  random_effects <- x$random_effects

  if (is.null(random_effects)) {
    stop("No random effects found in the model")
  }

  # Convert to matrix
  if (!is.matrix(random_effects)) {
    if (is.data.frame(random_effects)) {
      random_effects <- as.matrix(random_effects)
    } else if (is.vector(random_effects)) {
      random_effects <- matrix(random_effects, ncol = 1)
      if (is.null(colnames(random_effects))) {
        colnames(random_effects) <- "Random Effect"
      }
    } else {
      stop("random_effects must be a matrix, data frame, or vector")
    }
  }

  # Validate dimensions
  if (nrow(random_effects) != length(x$data)) {
    stop(sprintf(
      paste0(
        "Mismatch between number of random effects (%d) ",
        "and number of subjects (%d)"
      ),
      nrow(random_effects),
      length(x$data)
    ))
  }

  # Ensure column names
  if (is.null(colnames(random_effects))) {
    colnames(random_effects) <- paste0("RE", seq_len(ncol(random_effects)))
  }

  # Prepare data and plot
  re_data <- data.frame(
    id = rep(names(x$data), each = ncol(random_effects)),
    parameter = rep(colnames(random_effects), length(x$data)),
    value = as.vector(t(random_effects))
  )

  ggplot(re_data, aes(x = value, fill = parameter)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~parameter, scales = "free") +
    labs(x = "Value", y = "Density") +
    .theme_faceted() +
    theme(legend.position = "none")
}
