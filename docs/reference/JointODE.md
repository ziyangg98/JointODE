# Joint Modeling of Longitudinal and Survival Data Using ODEs

Implements a unified framework for jointly modeling longitudinal
biomarker trajectories and time-to-event outcomes using ordinary
differential equations (ODEs). The model captures complex non-linear
dynamics in biomarker evolution while simultaneously quantifying their
association with survival risk through shared random effects and
flexible hazard specifications.

## Usage

``` r
JointODE(
  longitudinal_formula,
  survival_formula,
  longitudinal_data,
  survival_data,
  gamma = 1,
  spline_baseline = list(degree = 2, n_knots = 1, knot_placement = "equal",
    boundary_knots = NULL),
  init = "default",
  control = list()
)
```

## Arguments

- longitudinal_formula:

  A formula specifying the longitudinal submodel. The left-hand side
  defines the observed biomarker. On the right-hand side, `biomarker`
  and `velocity` are reserved ODE terms for \\\log\omega_i^2\\ and
  \\\log(2\xi_i\omega_i)\\; all other fixed terms enter the forcing
  function. Random effects use lme-style syntax: `|` for full covariance
  and `||` for diagonal covariance.

- survival_formula:

  A formula for the survival submodel using `Surv(time, status)`
  notation on the left-hand side. The right-hand side specifies baseline
  hazard covariates (e.g., `Surv(event_time, event) ~ treatment + age`).

- longitudinal_data:

  A data frame containing repeated measurements with one row per
  observation. Required columns include subject identifier, measurement
  times, response values, and any covariates specified in the formula.

- survival_data:

  A data frame with time-to-event information containing one row per
  subject. Must include event/censoring times, event indicators, and
  baseline covariates.

- gamma:

  Numeric scalar specifying the power parameter for the velocity effect
  in the hazard function. When `gamma = 0`, velocity has no effect;
  `gamma = 1` uses linear velocity; `gamma = 2` uses squared velocity.
  Default is 1 (default: 1).

- spline_baseline:

  A list controlling the B-spline representation of the baseline hazard
  function with the following components:

  `degree`

  :   Polynomial degree of the B-spline basis functions (default: 2,
      quadratic splines)

  `n_knots`

  :   Number of interior knots for flexibility (default: 1)

  `knot_placement`

  :   Strategy for positioning knots: `"quantile"` places knots at
      quantiles of observed event times, `"equal"` uses equally-spaced
      knots (default: `"equal"`)

  `boundary_knots`

  :   A numeric vector of length 2 specifying the boundary knot
      locations. If `NULL`, automatically set to the range of observed
      event times (default: `NULL`)

- init:

  Initial values for model parameters. Can be:

  - `"default"` (default): Use zero/default initial values.

  - `"marginal"`: Use
    [`MarginalODE`](https://gongziyang.com/JointODE/reference/MarginalODE.md)
    to compute data-driven initial estimates.

  - A list with the same structure as the fitted model's `parameters`
    component for full manual control.

- control:

  A list of control parameters for optimization, or output from
  [`JointODE.control`](https://gongziyang.com/JointODE/reference/JointODE.control.md).

## Value

An S3 object of class `"JointODE"` containing fitted model results.

## Examples

``` r
if (FALSE) { # \dontrun{
data(sim)
fit <- JointODE(
  longitudinal_formula = observed ~
    biomarker + velocity + x1 + x2 + (1 + biomarker + velocity | id),
  survival_formula = Surv(time, status) ~ w1 + w2,
  longitudinal_data = sim$data$longitudinal_data,
  survival_data = sim$data$survival_data,
  init = sim$init
)
summary(fit)
} # }
```
