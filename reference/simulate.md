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
  longitudinal = list(xi = c(mean = 0.4, sd = 0.1), period = c(mean = 6, sd = 0.8),
    excitation = list(offset = 0, covariates = c(x1 = 0.5, x2 = -0.45),
    random_intercept_sd = 0.2), initial = list(biomarker = c(mean = -0.5, sd = 0.1),
    velocity = c(mean = 0, sd = 0.1)), error_sd = 0.1, n_measurements = 100),
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
      Default: c(mean = 6, sd = 0.8)

  excitation

  :   List specifying external forcing parameters:

      offset

      :   Constant excitation term \\f_0\\ (default: 0.0)

      covariates

      :   Named vector of covariate effects
          \\\boldsymbol{\beta}\_{exc}\\ on excitation (default: c(x1 =
          0.5, x2 = -0.45))

      random_intercept_sd

      :   Standard deviation of subject-specific forcing intercepts
          (default: 0.2)

  initial

  :   List specifying initial condition parameters:

      biomarker

      :   Population distribution of initial biomarker value \\m_i(0)\\,
          specified as `c(mean = ..., sd = ...)`. (default: c(mean =
          -0.5, sd = 0.1))

      velocity

      :   Population distribution of initial velocity \\\dot{m}\_i(0)\\,
          specified as `c(mean = ..., sd = ...)`. (default: c(mean = 0,
          sd = 0.1))

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

  :   An \\n \times 5\\ matrix of subject-specific random effects
      (centered at zero). Columns `initial_biomarker` and
      `initial_velocity` capture initial state variability; `log_omega2`
      and `log_2xi_omega` capture latent ODE parameter variability;
      `forcing_(Intercept)` matches the lme-style random intercept in
      `(1 + biomarker + velocity | id)`.

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

The initial biomarker is drawn from a population distribution:

- \\m_i(0) \sim \mathcal{N}(\mu\_{m,0}, \sigma\_{m,0}^2)\\

- \\\dot{m}\_i(0) \sim \mathcal{N}(\mu\_{v,0}, \sigma\_{v,0}^2)\\

The observed longitudinal measurements incorporate additive Gaussian
noise: \$\$y\_{ij} = m_i(t\_{ij}) + \epsilon\_{ij}\$\$ where
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

The default simulation includes random effects on initial biomarker,
initial velocity, the two dynamic log-coefficients, and the forcing
intercept; forcing covariate coefficients are shared fixed effects.

## Examples

``` r
# Example 1: Simple simulation with default parameters
sim_basic <- simulate(n_subjects = 20, seed = 123)

# Explore the output structure
names(sim_basic)
#> [1] "longitudinal_data" "survival_data"     "random_effects"
head(sim_basic$longitudinal_data)
#>   id time  biomarker   velocity acceleration   observed         x1        x2
#> 1  1  0.0 -0.3555449 0.07877388    0.3465426 -0.4394285 -0.5604756 -1.067824
#> 2  1  0.1 -0.3459894 0.11182253    0.3143598 -0.5155927 -0.5604756 -1.067824
#> 3  1  0.2 -0.3332881 0.14163673    0.2818890 -0.3080600 -0.5604756 -1.067824
#> 4  1  0.3 -0.3177691 0.16820137    0.2494177 -0.3568550 -0.5604756 -1.067824
#> 5  1  0.4 -0.2997557 0.19152952    0.2172084 -0.4766036 -0.5604756 -1.067824
#> 6  1  0.5 -0.2795697 0.21165975    0.1855044 -0.2760971 -0.5604756 -1.067824
head(sim_basic$survival_data)
#>   id      time status         w1 w2
#> 1  1 10.000000      0 -0.6947070  1
#> 2  2 10.000000      0 -0.2079173  0
#> 3  3 10.000000      0 -1.2653964  0
#> 4  4  5.872357      1  2.1689560  0
#> 5  5  7.995048      1  1.2079620  0
#> 6  6  9.898148      1 -1.1231086  1

# Check patient-specific dynamics
# Each patient has unique dynamics drawn from population distribution
head(sim_basic$random_effects)
#>      initial_biomarker initial_velocity  log_omega2 log_2xi_omega
#> [1,]       0.144455086       0.07877388 -0.26376808    -0.1373538
#> [2,]       0.045150405       0.07690422 -0.14559542    -0.6459329
#> [3,]       0.004123292       0.03322026 -0.06338155     0.2813305
#> [4,]      -0.042249683      -0.10083766  0.16670452    -0.1983814
#> [5,]      -0.205324722      -0.01194526 -0.36124337    -0.1924534
#> [6,]       0.113133721      -0.02803953  0.15936457     0.2868782
#>      forcing_(Intercept)
#> [1,]         -0.11506939
#> [2,]          0.12159286
#> [3,]         -0.32357654
#> [4,]         -0.01111239
#> [5,]          0.10388144
#> [6,]          0.06023067

# Random effects structure
colnames(sim_basic$random_effects)
#> [1] "initial_biomarker"   "initial_velocity"    "log_omega2"
#> [4] "log_2xi_omega"       "forcing_(Intercept)"
apply(sim_basic$random_effects, 2, sd)
#>   initial_biomarker    initial_velocity          log_omega2       log_2xi_omega
#>           0.1303001           0.0992504           0.2600702           0.2392747
#> forcing_(Intercept)
#>           0.1519673

# Example 2: Custom dynamics and survival
# \donttest{
sim_custom <- simulate(
  n_subjects = 50,
  longitudinal = list(
    xi = c(mean = 0.5, sd = 0.1),
    period = c(mean = 8, sd = 1),
    excitation = list(
      offset = 4.0, covariates = c(x1 = 0.8, x2 = -0.5),
      random_intercept_sd = 0.2
    ),
    initial = list(
      biomarker = c(mean = 3.8, sd = 0.2)
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
#> 32 18
# }
```
