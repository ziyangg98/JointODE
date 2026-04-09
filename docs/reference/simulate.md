# Simulate Data from a Joint Ordinary Differential Equation Model

Generates synthetic longitudinal and time-to-event data under a joint
modeling framework where the longitudinal biomarker trajectory follows a
damped harmonic oscillator model (second-order ODE), and the survival
process is associated with the biomarker dynamics through shared random
effects and trajectory-dependent hazard functions. The biomarker
dynamics are parameterized using physically interpretable parameters:
damping ratio, natural period, and excitation amplitude.

## Usage

``` r
simulate(
  n_subjects = 500,
  longitudinal = list(xi = c(mean = 0.4, sd = 0.1), period = c(mean = 6, sd = 0.1),
    excitation = list(offset = 0, covariates = c(x1 = 0.5, x2 = -0.45)), initial =
    list(biomarker = c(mean = -0.5, sd = 0.1), velocity = c(mean = -0.1, sd = 0.1)),
    error_sd = 0.1, n_measurements = 100),
  survival = list(baseline = list(type = "weibull", shape = 2, scale = 15), value = 0.8,
    slope = 2, gamma = 1, covariates = c(w1 = 0.6, w2 = -0.8)),
  covariates = list(x1 = list(type = "normal", mean = 0, sd = 1), x2 = list(type =
    "normal", mean = 0, sd = 1), w1 = list(type = "normal", mean = 0, sd = 1), w2 =
    list(type = "binary", prob = 0.5)),
  maxt = 10,
  seed = 42
)
```

## Arguments

- n_subjects:

  Integer specifying the number of subjects to simulate (default: 500)

- longitudinal:

  List specifying the longitudinal sub-model parameters:

  xi

  :   Damping ratio \\\xi\\ controlling the oscillation decay. Specified
      as c(mean = ..., sd = ...) for population mean and standard
      deviation. Values: \\\xi \< 1\\ (underdamped), \\\xi = 1\\
      (critically damped), \\\xi \> 1\\ (overdamped). Default: c(mean =
      0.4, sd = 0.1)

  period

  :   Natural period \\T\\ of oscillation in time units, related to
      natural frequency as \\\omega = 2\pi/T\\. Specified as c(mean =
      ..., sd = ...) for population mean and standard deviation.
      Default: c(mean = 6, sd = 0.1)

  excitation

  :   List specifying external forcing parameters:

      offset

      :   Constant excitation term \\f_0\\ (default: 0.0)

      covariates

      :   Named vector of covariate effects
          \\\boldsymbol{\beta}\_{exc}\\ on excitation (default: c(x1 =
          0.5, x2 = -0.45))

  initial

  :   List specifying initial condition parameters:

      biomarker

      :   Population distribution of initial biomarker value \\m_i(0)\\,
          specified as `c(mean = ..., sd = ...)`. (default: c(mean =
          -0.5, sd = 0.1))

      velocity

      :   Population distribution of initial velocity \\\dot{m}\_i(0)\\,
          specified as `c(mean = ..., sd = ...)`. (default: c(mean =
          -0.1, sd = 0.1))

  error_sd

  :   Standard deviation \\\sigma\_{\epsilon}\\ of the measurement error
      process (default: 0.1)

  n_measurements

  :   Number of longitudinal measurements per subject (default: 100)

- survival:

  List specifying the survival sub-model parameters:

  baseline

  :   List defining the Weibull baseline hazard function:

      type

      :   Character string specifying the baseline hazard type
          (currently only "weibull" is supported)

      shape

      :   Weibull shape parameter \\\kappa \> 0\\ (default: 2.0)

      scale

      :   Weibull scale parameter \\\lambda \> 0\\ (default: 15.0)

  value

  :   Association parameter \\\alpha_1\\ linking current biomarker value
      to hazard (default: 0.8)

  slope

  :   Association parameter \\\alpha_2\\ linking biomarker velocity (or
      its power) to hazard (default: 2.0)

  gamma

  :   Power parameter for velocity contribution, where \\\gamma = 0\\
      excludes velocity effect, \\\gamma = 1\\ uses linear velocity
      \\\alpha_2\dot{m}\_i(t)\\, and \\\gamma = 2\\ uses quadratic
      velocity \\\alpha_2\[\dot{m}\_i(t)\]^2\\ (default: 1)

  covariates

  :   Named vector of regression coefficients \\\boldsymbol{\phi}\\ for
      survival covariates (default: c(w1 = 0.6, w2 = -0.8))

- covariates:

  List defining the distributions of baseline covariates:

  x1

  :   List with `type = "normal"`, `mean = 0`, `sd = 1` for standardized
      continuous covariate (longitudinal)

  x2

  :   List with `type = "normal"`, `mean = 0`, `sd = 1` for standardized
      continuous covariate (longitudinal)

  w1

  :   List with `type = "normal"`, `mean = 0`, `sd = 1` for standardized
      continuous covariate (survival)

  w2

  :   List with `type = "binary"` and `prob = 0.5` for binary covariate
      (survival)

- maxt:

  Positive scalar specifying the maximum follow-up time in the study
  (default: 10 time units)

- seed:

  Integer seed for random number generation to ensure reproducibility
  (default: 42)

## Value

A list with two elements:

- `data`:

  A list containing three data components:

  `longitudinal_data`

  :   Data frame comprising longitudinal observations with columns:

      - `id`: Subject identifier (integer)

      - `time`: Observation time point (numeric)

      - `observed`: Measured biomarker value including measurement
        error, \\y\_{ij}\\

      - `biomarker`: True underlying biomarker value, \\m_i(t\_{ij})\\

      - `velocity`: First derivative of the biomarker trajectory,
        \\\dot{m}\_i(t\_{ij})\\

      - `acceleration`: Second derivative of the biomarker trajectory,
        \\\ddot{m}\_i(t\_{ij})\\

      - `x1`, `x2`: Time-invariant covariates (if specified)

  `survival_data`

  :   Data frame containing time-to-event data with columns:

      - `id`: Subject identifier

      - `time`: Observed event or censoring time, \\T_i\\

      - `status`: Event indicator, \\\delta_i\\ (1 = event observed, 0 =
        censored)

      - `w1`, `w2`: Baseline survival covariates (if specified)

  `random_effects`

  :   An \\n \times 4\\ matrix of subject-specific random effects
      (centered at zero). Columns `init_biomarker` and `init_velocity`
      capture initial state variability; `dyn_biomarker` and
      `dyn_velocity` capture ODE coefficient variability corresponding
      to the formula term `(biomarker + velocity | id)`.

- `init`:

  A list of initial parameter values suitable for passing to
  [`JointODE()`](https://gongziyang.com/JointODE/reference/JointODE.md),
  containing `$coefficients` and `$configurations`.

## Details

The simulation framework implements a joint model comprising
longitudinal and survival sub-models linked through shared random
effects and trajectory-dependent associations.

### Default Parameter Design

The default parameters are calibrated to achieve the following
properties:

- Approximately 60\\

- Initial biomarker difference between risk groups: ~0.2

- Final biomarker difference between event and censored groups: ~0.7

- Damping ratio centered at 0.4 (underdamped oscillations)

- Significant velocity differences between patients with different risk
  profiles, enabling survival prediction based on trajectory dynamics

These settings produce realistic heterogeneity in biomarker trajectories
while maintaining numerical stability for joint model estimation.

### Longitudinal Sub-model

The biomarker trajectory \\m_i(t)\\ for subject \\i\\ follows a damped
harmonic oscillator with external forcing: \$\$\ddot{m}\_i(t) +
2\xi\omega\dot{m}\_i(t) + \omega^2 m_i(t) = k \omega^2 \[f_0 +
\mathbf{X}\_i^T\boldsymbol{\beta}\_{exc}\]\$\$ where:

- \\\omega = 2\pi/T\\ is the natural angular frequency

- \\\xi\\ is the damping ratio determining oscillation behavior

- \\k\\ scales the excitation amplitude

- \\f_0\\ is the constant excitation term

- \\\mathbf{X}\_i\\ contains time-invariant covariates

- \\\boldsymbol{\beta}\_{exc}\\ represents covariate effects on
  excitation

Initial conditions are drawn from population distributions:

- \\m_i(0) \sim \mathcal{N}(\mu\_{m,0}, \sigma\_{m,0}^2)\\

- \\\dot{m}\_i(0) \sim \mathcal{N}(\mu\_{v,0}, \sigma\_{v,0}^2)\\

The observed longitudinal measurements incorporate additive Gaussian
noise: \$\$y\_{ij} = m_i(t\_{ij}) + b_i + \epsilon\_{ij}\$\$ where
\\\epsilon\_{ij} \sim \mathcal{N}(0, \sigma\_\epsilon^2)\\ represents
independent measurement error.

### Survival Sub-model

The instantaneous hazard function incorporates both current biomarker
value and velocity (or its power): \$\$\lambda_i(t) =
\lambda_0(t)\exp(\alpha_1 m_i(t) + \alpha_2\[\dot{m}\_i(t)\]^\gamma +
\mathbf{W}\_i^T\boldsymbol{\phi})\$\$ where:

- \\\lambda_0(t)\\ denotes the Weibull baseline hazard: \\\lambda_0(t) =
  (\kappa/\lambda)(t/\lambda)^{\kappa-1}\\

- \\\alpha_1\\ quantifies the association with current biomarker value

- \\\alpha_2\\ quantifies the association with biomarker velocity (or
  its power)

- \\\gamma \in \\0, 1, 2\\\\ determines the power of velocity: 0 (no
  velocity effect), 1 (linear), or 2 (quadratic)

- \\\mathbf{W}\_i\\ contains baseline covariates

- \\\boldsymbol{\phi}\\ represents covariate effects on survival

### Patient-Specific Dynamics

The function models patient heterogeneity through continuous
distributions of dynamics parameters \\\xi_i\\ (damping ratio) and
\\T_i\\ (period): \$\$\xi_i \sim \mathcal{N}(\mu\_\xi, \sigma\_\xi^2),
\quad T_i \sim \mathcal{N}(\mu_T, \sigma_T^2)\$\$ These are transformed
to ODE parameters via: \$\$\omega_i = 2\pi/T_i, \quad \beta\_{1,i} =
-\omega_i^2, \quad \beta\_{2,i} = -2\xi_i\omega_i\$\$ The transformation
uses the Delta method to preserve the correct covariance structure in
the ODE parameter space.

Only \\\beta\_{1,i}\\ and \\\beta\_{2,i}\\ vary across subjects; the
offset and covariate coefficients are shared fixed effects.

## Examples

``` r
# Example 1: Simple simulation with default parameters
sim_basic <- simulate(n_subjects = 20, seed = 123)

# Explore the output structure
names(sim_basic)
#> [1] "longitudinal_data" "survival_data"     "random_effects"
head(sim_basic$longitudinal_data)
#>   id time  biomarker   velocity acceleration   observed         x1        x2
#> 1  1  0.0 -0.5575347 0.04445509    0.7775724 -0.6043193 -0.5604756 -1.067824
#> 2  1  0.1 -0.5493074 0.11903932    0.7135863 -0.5364070 -0.5604756 -1.067824
#> 3  1  0.2 -0.5339453 0.18708254    0.6468857 -0.3797485 -0.5604756 -1.067824
#> 4  1  0.3 -0.5121163 0.24835639    0.5783474 -0.3657800 -0.5604756 -1.067824
#> 5  1  0.4 -0.4845043 0.30271866    0.5088018 -0.4791397 -0.5604756 -1.067824
#> 6  1  0.5 -0.4518045 0.35010838    0.4390307 -0.5044062 -0.5604756 -1.067824
head(sim_basic$survival_data)
#>   id      time status         w1 w2
#> 1  1 10.000000      0 -0.6947070  1
#> 2  2 10.000000      0 -0.2079173  0
#> 3  3  4.584052      1 -1.2653964  0
#> 4  4 10.000000      0  2.1689560  0
#> 5  5  9.142940      1  1.2079620  0
#> 6  6 10.000000      0 -1.1231086  1

# Check patient-specific dynamics
# Each patient has unique dynamics drawn from population distribution
head(sim_basic$random_effects)
#>      init_biomarker init_velocity dyn_biomarker dyn_velocity
#> [1,]   -0.057534696   0.144455086  -0.034149079   -0.1227326
#> [2,]    0.060796432   0.045150405  -0.023449893   -0.5784090
#> [3,]   -0.161788271   0.004123292  -0.005669442    0.2520641
#> [4,]   -0.005556197  -0.042249683   0.019281452   -0.1778817
#> [5,]    0.051940720  -0.205324722  -0.046803869   -0.1719762
#> [6,]    0.030115336   0.113133721   0.022277514    0.2567820

# Random effects structure
colnames(sim_basic$random_effects)
#> [1] "init_biomarker" "init_velocity"  "dyn_biomarker"  "dyn_velocity"
apply(sim_basic$random_effects, 2, sd)
#> init_biomarker  init_velocity  dyn_biomarker   dyn_velocity
#>     0.07598365     0.13030015     0.03326469     0.21422580

# Example 2: Custom dynamics and survival
# \donttest{
sim_custom <- simulate(
  n_subjects = 50,
  longitudinal = list(
    xi = c(mean = 0.5, sd = 0.1),
    period = c(mean = 8, sd = 1),
    excitation = list(offset = 4.0, covariates = c(x1 = 0.8, x2 = -0.5)),
    initial = list(
      biomarker = c(mean = 3.8, sd = 0.2),
      velocity = c(mean = -0.1, sd = 0.1)
    ),
    error_sd = 0.1,
    n_measurements = 20
  ),
  survival = list(
    baseline = list(type = "weibull", shape = 3.0, scale = 23),
    value = 0.4, slope = 1.5, gamma = 1,
    covariates = c(w1 = 0.4, w2 = -0.6)
  ),
  seed = 42
)
table(sim_custom$survival_data$status) # event vs censored
#>
#>  0  1
#> 45  5
# }
```
