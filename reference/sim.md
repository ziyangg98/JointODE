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

  :   Data frame with 17,222 longitudinal measurements from 200
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

  :   Matrix (200 x 5) with subject-specific random effects:

      - initial_biomarker: Initial biomarker value

      - initial_velocity: Initial biomarker velocity

      - log_omega2: Random effect on \\\log \omega_i^2\\

      - log_2xi_omega: Random effect on \\\log(2\xi_i\omega_i)\\

      - forcing\_(Intercept): Random forcing intercept

- init:

  A list containing initial parameter values:

  coefficients

  :   True parameter values used in simulation:

      - baseline: B-spline coefficients for log baseline hazard

      - longitudinal: Fixed effects for ODE dynamics (log_omega2,
        log_2xi_omega, forcing intercept, covariate effects)

      - hazard: Association parameters (value, slope) and survival
        covariate effects

      - measurement_error_sd: Measurement error SD (0.1)

      - random_effect_sigma: 5x5 covariance matrix for random effects
        (initial_biomarker, initial_velocity, log_omega2, log_2xi_omega,
        forcing intercept)

  configurations

  :   Model configuration:

      - baseline: B-spline configuration for baseline hazard

      - gamma: Power parameter for velocity effect (1)

## Source

Generated using `.create_example_data(n_subjects = 200, seed = 123)`

## Details

The dataset contains 200 subjects with heterogeneous ODE dynamics.
Subject-specific dynamics are characterized by random effects on latent
ODE log-coefficients. Population means: damping ratio \\\xi \approx
0.4\\, period \\T \approx 6\\.

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
#>   ..$ longitudinal_data:'data.frame':    17222 obs. of  8 variables:
#>   ..$ survival_data    :'data.frame':    200 obs. of  5 variables:
#>   ..$ random_effects   : num [1:200, 1:5] -0.075 -0.0322 -0.1148 0.0354 0.0425 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>  $ init:List of 2
#>   ..$ coefficients  :List of 6
#>   ..$ configurations:List of 4

# Access longitudinal data
head(sim$data$longitudinal_data)
#>   id time  biomarker    velocity acceleration   observed         x1      x2
#> 1  1  0.0 -0.5749726 -0.08209867   -0.3277411 -0.7508544 -0.5604756 2.19881
#> 2  1  0.1 -0.5847617 -0.11312631   -0.2924862 -0.5713544 -0.5604756 2.19881
#> 3  1  0.2 -0.5974767 -0.14053847   -0.2555021 -0.7562703 -0.5604756 2.19881
#> 4  1  0.3 -0.6127449 -0.16419177   -0.2174216 -0.6317456 -0.5604756 2.19881
#> 5  1  0.4 -0.6301871 -0.18400663   -0.1788435 -0.5694431 -0.5604756 2.19881
#> 6  1  0.5 -0.6494179 -0.19996219   -0.1403353 -0.6157630 -0.5604756 2.19881

# Summary of survival outcomes
summary(sim$data$survival_data$time)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#>  0.5283  7.7731 10.0000  8.5946 10.0000 10.0000
table(sim$data$survival_data$status)
#>
#>   0   1
#> 135  65

# Population parameters
sim$init$coefficients$longitudinal
#>    log_omega2 log_2xi_omega    dyn_offset        dyn_x1        dyn_x2
#>    0.09223519   -0.17702595    0.00000000    0.54831136   -0.49348022
sim$init$coefficients$random_effect_sigma
#>                     initial_biomarker initial_velocity log_omega2 log_2xi_omega
#> initial_biomarker                0.01             0.00 0.00000000    0.00000000
#> initial_velocity                 0.00             0.01 0.00000000    0.00000000
#> log_omega2                       0.00             0.00 0.07048641    0.00000000
#> log_2xi_omega                    0.00             0.00 0.00000000    0.07824622
#> forcing_(Intercept)              0.00             0.00 0.00000000    0.00000000
#>                     forcing_(Intercept)
#> initial_biomarker                  0.00
#> initial_velocity                   0.00
#> log_omega2                         0.00
#> log_2xi_omega                      0.00
#> forcing_(Intercept)                0.04

# Random effects structure
colnames(sim$data$random_effects)
#> [1] "initial_biomarker"   "initial_velocity"    "log_omega2"
#> [4] "log_2xi_omega"       "forcing_(Intercept)"
apply(sim$data$random_effects, 2, sd)
#>   initial_biomarker    initial_velocity          log_omega2       log_2xi_omega
#>          0.10005173          0.09961304          0.27385366          0.28770989
#> forcing_(Intercept)
#>          0.19221629

if (FALSE) { # \dontrun{
# Fit a Joint ODE model using this data
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
