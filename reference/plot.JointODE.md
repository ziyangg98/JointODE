# Plot Method for JointODE Objects

Provides comprehensive visualization tools for fitted JointODE models,
including longitudinal trajectories, phase space plots, survival curves,
diagnostic plots, and association patterns. All plots are generated
using ggplot2 for modern, publication-ready graphics.

## Usage

``` r
# S3 method for class 'JointODE'
plot(
  x,
  type = c("overview", "trajectory_biomarker", "trajectory_velocity",
    "phase_biomarker_velocity", "phase_velocity_acceleration", "survival",
    "diagnostic_residuals", "diagnostic_residuals_time", "diagnostic_qq",
    "diagnostic_random_effects", "diagnostic_association"),
  subject_ids = NULL,
  show_observed = TRUE,
  show_individual = TRUE,
  by = NULL,
  n_groups = 4,
  cols = NULL,
  span = 0.75,
  ...
)
```

## Arguments

- x:

  An object of class `JointODE`

- type:

  Character string specifying the plot type:

  - `"overview"`: Panel of 4 key plots (survival, biomarker, velocity,
    phase space)

  - `"trajectory_biomarker"`: Longitudinal biomarker trajectories. Can
    show individual curves or group averages based on `by`

  - `"trajectory_velocity"`: Longitudinal velocity trajectories. Can
    show individual curves or group averages based on `by`

  - `"phase_biomarker_velocity"`: Phase space plots (biomarker vs
    velocity). Can show individual curves or group averages based on
    `by`

  - `"phase_velocity_acceleration"`: Phase space plots (velocity vs
    acceleration). Can show individual curves or group averages based on
    `by`

  - `"survival"`: Survival probability curves. Can be stratified by `by`
    variable

  - `"diagnostic_residuals"`: Residuals vs fitted values

  - `"diagnostic_residuals_time"`: Residuals vs time

  - `"diagnostic_qq"`: Normal Q-Q plot of residuals

  - `"diagnostic_random_effects"`: Random effects distribution

  - `"diagnostic_association"`: Association between biomarker features
    and hazard (not yet implemented)

- subject_ids:

  Character or numeric vector of subject IDs to plot. If `NULL`,
  displays all subjects. Ignored when `by` is specified (shows group
  averages instead).

- show_observed:

  Logical; whether to show observed data points in faceted view when
  specific subject IDs are provided (default: TRUE). Ignored when
  showing all subjects (overlay view)

- show_individual:

  Logical; whether to show individual subject curves in addition to
  group means (default: TRUE). Applies to overview, trajectories, phase
  space, and survival plots. When FALSE, only shows group mean curves.
  Ignored for residuals, residuals_time, qq, and random_effects plots.

- by:

  Character string naming a covariate in the data to use for grouping.
  When specified for trajectories/phase/survival plots, shows
  group-averaged curves instead of individual subjects. If `NULL`, shows
  overall population average (for trajectories/phase) or individual
  curves (default behavior).

- n_groups:

  Integer; number of groups to create when `by` is a continuous variable
  (default: 4). Subjects are divided into groups using quantiles.
  Ignored when `by` is categorical or `NULL`.

- cols:

  Color palette for plots (default: NULL uses built-in colors)

- span:

  Numeric value controlling the degree of smoothing for loess curves
  (default: 0.75). Larger values produce smoother curves. Only applies
  to mean curves in overview, survival, biomarker, and velocity plots.

- ...:

  Additional parameters passed to ggplot2 functions

## Value

A ggplot2 object (or patchwork object for multi-panel plots) that can be
further customized

## Details

The function provides several visualization types to assess model fit
and understand the joint modeling results:

**Overview plots** provide a quick diagnostic panel with four key plots:
survival probability (with individual curves), biomarker trajectory
(with individual curves), velocity trajectory (with individual curves),
and phase space plots.

**Trajectory plots** show observed biomarker values overlaid with
model-fitted trajectories, allowing assessment of longitudinal model fit
at the individual level.

**Dynamics plots** display phase portraits (biomarker vs velocity)
revealing the ODE system's behavior and stability properties.

**Survival plots** show predicted survival probabilities over time,
optionally stratified by covariates.

**Diagnostic plots** include residual plots and QQ plots for checking
model assumptions.

**Association plots** visualize how biomarker features (value and slope)
relate to hazard ratios.

## Examples

``` r
if (FALSE) { # \dontrun{
library(JointODE)

# Load example dataset
data(sim)

# Fit joint ODE model
fit <- JointODE(
  longitudinal_formula = observed ~ x1 + x2 +
    (biomarker + velocity | id),
  survival_formula = Surv(time, status) ~ w1 + w2,
  longitudinal_data = sim$data$longitudinal_data,
  survival_data = sim$data$survival_data,
  control = list(maxit = 5)
)

# Overview plot
plot(fit)
plot(fit, type = "overview")

# Individual trajectory plots
plot(fit, type = "trajectory_biomarker") # Shows all subjects
plot(fit, type = "trajectory_biomarker", subject_ids = c("1", "2", "3"))
plot(fit, type = "trajectory_velocity") # Velocity trajectories

# Phase space plots
plot(fit, type = "phase_biomarker_velocity")
plot(fit, type = "phase_velocity_acceleration")

# Group-stratified plots by survival covariates
plot(fit, type = "survival", by = "w1") # Continuous variable (auto-grouped)
plot(fit, type = "trajectory_biomarker", by = "w2")
plot(fit, type = "trajectory_velocity", by = "w2")
plot(fit, type = "phase_biomarker_velocity", by = "w1")

# Stratify survival by biomarker/velocity
plot(fit, type = "survival", by = "biomarker")
plot(fit, type = "survival", by = "velocity")

# Diagnostic plots
plot(fit, type = "diagnostic_residuals")
plot(fit, type = "diagnostic_residuals_time")
plot(fit, type = "diagnostic_qq")
plot(fit, type = "diagnostic_random_effects")

# Customize colors
plot(fit, type = "survival", cols = c("red", "blue", "green"))

# Adjust smoothing
plot(fit, type = "overview", span = 0.5) # Less smooth
plot(fit, type = "survival", span = 1.0) # More smooth
} # }
```
