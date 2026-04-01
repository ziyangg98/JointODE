# JointODE

The **JointODE** package provides a unified framework for joint modeling
of longitudinal biomarker measurements and time-to-event outcomes using
ordinary differential equations (ODEs).

## Overview

JointODE models biomarker evolution using **ordinary differential
equations** (ODEs) and links them to survival outcomes. Unlike
traditional approaches using splines or polynomials, ODEs naturally
capture the dynamic behavior of biomarker trajectories.

### Model Structure

**Longitudinal Model:** Biomarker trajectories evolve according to a
second-order ODE:

\ddot{m}\_i(t) + 2 \xi_i \omega_i \dot{m}\_i(t) + \omega_i^2 m_i(t) = k
\omega_i^2 \mu_i(t)

where m_i(t) is the biomarker value, \dot{m}\_i(t) is velocity (rate of
change), and \ddot{m}\_i(t) is acceleration. Subject-specific dynamics
are characterized by natural frequency \omega_i and damping ratio \xi_i,
with \mu_i(t) representing external forcing (e.g., treatment effects,
covariates). Individual heterogeneity is captured through random effects
on these ODE parameters.

**Survival Model:** The hazard function incorporates biomarker dynamics:

\lambda_i(t) = \lambda\_{0}(t)\exp\left\[\alpha_1 m_i(t) + \alpha_2
\dot{m}\_i^{\gamma}(t) + \mathbf{W}\_i^{\top}\boldsymbol{\phi}\right\]

where \gamma \in \\0, 1, 2\\ controls the velocity power (0: no
velocity, 1: linear, 2: quadratic).

### Key Features

- **Subject-specific dynamics**: Random effects on ODE parameters
  capture individual heterogeneity
- **Flexible hazard**: B-spline baseline hazard with biomarker value and
  velocity effects
- **Efficient computation**: C++ implementation with parallel processing
  support

For full mathematical details, see the [technical
documentation](http://gongziyang.com/JointODE/articles/technical-details.md).

## Installation

You can install the development version of JointODE from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("ziyangg98/JointODE")
```

## Quick Start

Here’s a basic example using the included simulated dataset:

``` r
library(JointODE)
#>
#> Attaching package: 'JointODE'
#> The following object is masked from 'package:stats':
#>
#>     simulate

# Load example dataset (200 subjects with longitudinal and survival data)
data(sim)

# Fit joint ODE model
longitudinal_data <- sim$data$longitudinal_data[
  , c("id", "time", "observed", "x1", "x2")
]
t0 <- proc.time()
fit <- JointODE(
  longitudinal_formula = observed ~ biomarker + velocity + x1 + x2 +
    (biomarker + velocity | id),
  survival_formula = Surv(time, status) ~ w1 + w2,
  longitudinal_data = longitudinal_data,
  survival_data = sim$data$survival_data,
  init = "marginal",
  control = list(parallel = TRUE)
)
cat(sprintf("Elapsed: %.1f s\n", (proc.time() - t0)["elapsed"]))
#> Elapsed: 175.1 s
```

``` r
# Model summary
summary(fit)
#>
#> Call:
#> JointODE(longitudinal_formula = observed ~ biomarker + velocity +
#>     x1 + x2 + (biomarker + velocity | id), survival_formula = Surv(time,
#>     status) ~ w1 + w2, longitudinal_data = longitudinal_data,
#>     survival_data = sim$data$survival_data, init = "marginal",
#>     control = list(parallel = TRUE))
#>
#> Data Descriptives:
#> Longitudinal Process            Survival Process
#> Number of Observations: 17350   Number of Events: 61 (30%)
#> Number of Subjects: 200
#>
#>        AIC        BIC     logLik
#> -32488.560 -32396.207  16272.280
#>
#> Coefficients:
#> Longitudinal Process: Second-Order ODE Model
#>             Estimate Std. Error  z value Pr(>|z|)
#> biomarker    -1.0844     0.0083 -129.995   <2e-16 ***
#> velocity     -0.8223     0.0162  -50.876   <2e-16 ***
#> (Intercept)   0.0003     0.0015    0.217    0.828
#> x1            0.5427     0.0042  128.237   <2e-16 ***
#> x2           -0.4885     0.0039 -124.359   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> ODE System Characteristics:
#>                    Estimate Std. Error z value Pr(>|z|)
#> T (period)           6.0336     0.0232  259.99   <2e-16 ***
#> xi (damping ratio)   0.3948     0.0075   52.64   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Survival Process: Proportional Hazards Model
#>         Estimate Std. Error z value Pr(>|z|)
#> alpha_1   0.6774     0.1873   3.617 0.000298 ***
#> alpha_2   1.9869     1.0253   1.938 0.052646 .
#> w1        0.6762     0.1240   5.453 4.95e-08 ***
#> w2       -1.3353     0.2857  -4.673 2.97e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Baseline Hazard: B-spline with 6 basis functions
#> (Coefficients range: [-6.506, -0.518] )
#>
#> Initial State: Population Mean
#>    Estimate Std. Error z value Pr(>|z|)
#> m0  -0.5061     0.0085  -59.76  < 2e-16 ***
#> v0  -0.0950     0.0228   -4.16 3.19e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Variance Components:
#> Measurement Error SD: 0.098143
#> Random Effect Covariance Matrix:
#>            [,1]       [,2]       [,3]      [,4]
#> [1,]  0.0119511  0.0307804 -0.0004738 -0.001493
#> [2,]  0.0307804  0.0947550 -0.0001958 -0.003552
#> [3,] -0.0004738 -0.0001958  0.0017180  0.001939
#> [4,] -0.0014930 -0.0035520  0.0019394  0.038128
#>
#> Model Diagnostics:
#> C-index (Concordance): 0.868
#> Convergence: Converged after 8 iterations

# Plot results
plot(fit)
```

![](reference/figures/README-output-1.png)

The formula specifies:

- **ODE terms**: `biomarker` and `velocity` are the state variables
  (value and slope) in the ODE
- **Covariates**: `x1` and `x2` are external variables affecting the
  dynamics
- **Random effects**: `(biomarker + velocity | id)` allows
  subject-specific coefficients on the ODE value and slope terms

## Learn More

- **Getting Started**: See `vignette("JointODE")` for a detailed
  tutorial
- **Technical Details**: See
  [`vignette("technical-details")`](https://gongziyang.com/JointODE/articles/technical-details.md)
  for mathematical formulations
- **Model Comparison**: See `vignette("comparison")` for comparisons
  with traditional joint models

## Code of Conduct

Please note that the JointODE project is released with a [Contributor
Code of Conduct](http://gongziyang.com/JointODE/CODE_OF_CONDUCT.md). By
contributing to this project, you agree to abide by its terms.
