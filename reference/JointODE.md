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
  state = NULL,
  spline_baseline = list(degree = 2, n_knots = 0, knot_placement = "equal",
    boundary_knots = NULL),
  init = NULL,
  control = list(),
  ...
)
```

## Arguments

- longitudinal_formula:

  A formula specifying the longitudinal submodel. The left-hand side
  defines the response variable, while the right-hand side specifies
  fixed effects including time-varying and baseline covariates (e.g.,
  `biomarker ~ time + treatment + age`).

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

- state:

  A matrix specifying initial conditions for the ODE system with two
  columns: initial biomarker values and their first derivatives. Each
  row corresponds to one subject. If `NULL`, defaults to a zero matrix
  with appropriate dimensions (default: `NULL`).

- spline_baseline:

  A list controlling the B-spline representation of the baseline hazard
  function with the following components:

  `degree`

  :   Polynomial degree of the B-spline basis functions (default: 2,
      quadratic splines)

  `n_knots`

  :   Number of interior knots for flexibility (default: 0, providing
      moderate flexibility)

  `knot_placement`

  :   Strategy for positioning knots: `"quantile"` places knots at
      quantiles of observed event times, `"equal"` uses equally-spaced
      knots (default: `"equal"`)

  `boundary_knots`

  :   A numeric vector of length 2 specifying the boundary knot
      locations. If `NULL`, automatically set to the range of observed
      event times (default: `NULL`)

- init:

  Optional initial values for model parameters. Can be either:

  - `NULL` (default): Automatically calls
    [`MarginalODE`](http://gongziyang.com/JointODE/reference/MarginalODE.md)
    to compute data-driven initial estimates for longitudinal parameters
    and random effects. Recommended for most applications.

  - A list with the same structure as the fitted model's `parameters`
    component for full manual control

  When providing a list, it should have elements:

  `coefficients`

  :   A list containing:

      - `baseline`: Vector of B-spline coefficients for baseline hazard
        (length = number of spline basis functions)

      - `hazard`: Vector of hazard parameters including association
        parameters (2) and survival covariates

      - `longitudinal`: Vector of longitudinal fixed effects including
        intercept and covariates

      - `measurement_error_sd`: Residual standard deviation (positive
        scalar)

      - `random_effect_sigma`: Random effect covariance matrix (positive
        definite matrix)

  `configurations`

  :   Optional; if not provided, will use spline configuration from
      `spline_baseline`

- control:

  A list of control parameters for optimization, or output from
  [`JointODE.control`](http://gongziyang.com/JointODE/reference/JointODE.control.md).
  Key parameters include:

  `maxit`

  :   Maximum number of EM iterations (default: 200)

  `atol`

  :   Absolute tolerance for parameter convergence (default: 5e-4)

  `rtol`

  :   Relative tolerance for log-likelihood convergence (default: 1e-5)

  `verbose`

  :   Verbosity level: FALSE/0 for silent, TRUE/1 for progress messages,
      2 for detailed output (default: FALSE)

  `parallel`

  :   Logical flag enabling parallel computation (default: FALSE)

  `n_cores`

  :   Number of CPU cores for parallel processing. If 0, automatically
      detects available cores (default: 0)

  `quad_level`

  :   Quadrature level for numerical integration (default: 4)

  See
  [`JointODE.control`](http://gongziyang.com/JointODE/reference/JointODE.control.md)
  for complete details and examples.

- ...:

  Additional arguments passed to internal optimization routines.

## Value

An S3 object of class `"JointODE"` containing fitted model results:

- `parameters`:

  A list containing all estimated parameters:

  - `coefficients`: Named list with `baseline` (B-spline coefficients
    for baseline hazard), `hazard` (association and survival covariate
    effects), `longitudinal` (longitudinal fixed effects),
    `measurement_error_sd` (residual standard deviation), and
    `random_effect_sigma` (random effect covariance matrix)

  - `configurations`: Model configuration including spline basis
    specifications

- `logLik`:

  Maximum log-likelihood value achieved at convergence

- `AIC`:

  Akaike Information Criterion for model comparison

- `BIC`:

  Bayesian Information Criterion adjusted for sample size

- `cindex`:

  Concordance index (C-index) measuring the model's discrimination
  ability for survival prediction

- `convergence`:

  List containing convergence diagnostics:

  - `converged`: Logical indicating convergence status

  - `em_iterations`: Number of EM iterations performed

  - `message`: Descriptive convergence message

- `random_effects`:

  List containing random effects estimates:

  - `estimates`: Posterior means of subject-specific random effects

  - `variances`: Posterior variances of random effects

- `data`:

  Processed data used for model fitting in internal format

- `control`:

  List of control parameters used in optimization

- `call`:

  The matched function call for reproducibility

## Details

The joint modeling framework integrates longitudinal and survival
processes through a shared random effects structure. The longitudinal
biomarker evolution is characterized by a system of ODEs that can
accommodate non-linear dynamics, feedback mechanisms, and complex
temporal patterns. The survival component employs a proportional hazards
model where the instantaneous risk depends on features derived from the
longitudinal trajectory.

Two association structures are supported:

- Current value: hazard depends on the biomarker level at time t

- Rate of change: hazard depends on the biomarker's instantaneous slope

Parameter estimation employs an Expectation-Maximization (EM) algorithm
with:

- E-step: Multivariate Gauss-Hermite quadrature (mvQuad) for posterior
  computation of random effects

- M-step: Optimization for fixed effects and mvQuad-based closed-form
  updates for variance parameters

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate example data
sim <- simulate(n_subjects = 100, n_times = 10)

# Fit with default control parameters
fit1 <- JointODE(
  longitudinal_formula = observed ~ x1 + x2,
  survival_formula = Surv(event_time, event) ~ x1 + x2,
  longitudinal_data = sim$longitudinal_data,
  survival_data = sim$survival_data
)

# Fit with custom control parameters using JointODE.control()
control <- JointODE.control(
  maxit = 200, atol = 1e-4, rtol = 1e-6, verbose = TRUE
)
fit2 <- JointODE(
  longitudinal_formula = observed ~ x1 + x2,
  survival_formula = Surv(event_time, event) ~ x1 + x2,
  longitudinal_data = sim$longitudinal_data,
  survival_data = sim$survival_data,
  control = control
)

# Fit with control parameters as a list
# By default, uses MarginalODE for initialization (init = NULL)
fit3 <- JointODE(
  longitudinal_formula = observed ~ x1 + x2,
  survival_formula = Surv(event_time, event) ~ x1 + x2,
  longitudinal_data = sim$longitudinal_data,
  survival_data = sim$survival_data,
  control = list(maxit = 50, verbose = TRUE)
)

summary(fit1)
} # }
```
