# Example Dataset for Joint ODE Model

A simulated dataset with 200 subjects generated using the Joint ODE
Model framework, demonstrating patient heterogeneity through continuous
distributions of dynamics parameters, with longitudinal biomarker
measurements and survival outcomes.

## Usage

``` r
sim
```

## Format

A list with two components:

- data:

  A list containing the simulated data:

  longitudinal_data

  :   Data frame with 17,347 longitudinal measurements from 200
      subjects:

      - id: Patient identifier (1-200)

      - time: Measurement time

      - biomarker: True biomarker value (ODE solution)

      - velocity: First derivative \\dy/dt\\

      - acceleration: Second derivative \\d^2y/dt^2\\

      - observed: Observed biomarker with measurement error

      - x1, x2: Longitudinal covariates (normal distribution)

  survival_data

  :   Data frame with survival outcomes (200 patients):

      - id: Patient identifier

      - time: Event or censoring time

      - status: Event indicator (1=event, 0=censored)

      - w1, w2: Survival covariates

  random_effects

  :   Matrix (200 x 4) with subject-specific random effects:

      - init_biomarker: Initial biomarker value

      - init_velocity: Initial velocity

      - dyn_biomarker: ODE biomarker coefficient (\\-\omega^2\\)

      - dyn_velocity: ODE velocity coefficient (\\-2\xi\omega\\)

- init:

  A list containing initial parameter values:

  coefficients

  :   True parameter values used in simulation:

      - baseline: B-spline coefficients for log baseline hazard

      - longitudinal: Fixed effects for ODE dynamics (dyn_offset,
        dyn_biomarker, dyn_velocity, covariate effects)

      - hazard: Association parameters (value, slope) and survival
        covariate effects

      - measurement_error_sd: Measurement error SD (0.1)

      - random_effect_sigma: 4x4 covariance matrix for random effects
        (init_biomarker, init_velocity, dyn_biomarker, dyn_velocity)

  configurations

  :   Model configuration:

      - baseline: B-spline configuration for baseline hazard

      - gamma: Power parameter for velocity effect (1)

## Source

Generated using `.create_example_data(n_subjects = 200, seed = 123)`

## Details

The dataset contains 200 subjects with heterogeneous ODE dynamics.
Subject-specific dynamics are characterized by random effects on
dyn_biomarker (\\-\omega^2\\) and dyn_velocity (\\-2\xi\omega\\)
parameters. Population means: damping ratio \\\xi \approx 0.4\\, period
\\T \approx 6\\.

## See also

[`JointODE`](https://gongziyang.com/JointODE/reference/JointODE.md) for
model fitting,
[`simulate`](https://gongziyang.com/JointODE/reference/simulate.md) for
data generation

## Examples

``` r
# Load the data
data(sim)

# Examine structure
str(sim, max.level = 2)
#> List of 2
#>  $ data:List of 3
#>   ..$ longitudinal_data:'data.frame':    17319 obs. of  8 variables:
#>   ..$ survival_data    :'data.frame':    200 obs. of  5 variables:
#>   ..$ random_effects   : num [1:200, 1:4] -0.065 -0.1003 -0.0535 -0.011 0.06 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>  $ init:List of 2
#>   ..$ coefficients  :List of 6
#>   ..$ configurations:List of 4

# Access longitudinal data
head(sim$data$longitudinal_data)
#>   id time  biomarker   velocity acceleration   observed       x1        x2
#> 1  1  0.0 -0.5650284 -0.2185547     2.537227 -0.5970006 1.370958 -2.000929
#> 2  1  0.1 -0.5745436  0.0247568     2.327828 -0.4739660 1.370958 -2.000929
#> 3  1  0.2 -0.5607854  0.2468262     2.112843 -0.5053470 1.370958 -2.000929
#> 4  1  0.3 -0.5259008  0.4472288     1.894937 -0.3926945 1.370958 -2.000929
#> 5  1  0.4 -0.4720671  0.6257977     1.676578 -0.3919803 1.370958 -2.000929
#> 6  1  0.5 -0.4014663  0.7826035     1.460021 -0.3512690 1.370958 -2.000929

# Summary of survival outcomes
summary(sim$data$survival_data$time)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#>   0.665   8.920  10.000   8.644  10.000  10.000
table(sim$data$survival_data$status)
#>
#>   0   1
#> 142  58

# Population parameters
sim$init$coefficients$longitudinal
#> dyn_biomarker  dyn_velocity    dyn_offset        dyn_x1        dyn_x2
#>    0.09223519   -0.17702595    0.00000000    0.54831136   -0.49348022
sim$init$coefficients$random_effect_sigma
#>      [,1] [,2]         [,3]         [,4]
#> [1,] 0.01 0.00 0.0000000000 0.0000000000
#> [2,] 0.00 0.01 0.0000000000 0.0000000000
#> [3,] 0.00 0.00 0.0011111111 0.0005555556
#> [4,] 0.00 0.00 0.0005555556 0.0627777778

# Random effects structure
colnames(sim$data$random_effects)
#> [1] "init_biomarker" "init_velocity"  "dyn_biomarker"  "dyn_velocity"
apply(sim$data$random_effects, 2, sd)
#> init_biomarker  init_velocity  dyn_biomarker   dyn_velocity
#>     0.09744381     0.09336398     0.03668010     0.24254412

if (FALSE) { # \dontrun{
# Fit a Joint ODE model using this data
fit <- JointODE(
  longitudinal_formula = observed ~
    biomarker + velocity + x1 + x2 + (biomarker + velocity | id),
  survival_formula = Surv(time, status) ~ w1 + w2,
  longitudinal_data = sim$data$longitudinal_data,
  survival_data = sim$data$survival_data,
  init = sim$init
)
summary(fit)
} # }
```
