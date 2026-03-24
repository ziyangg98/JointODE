
<!-- README.md is generated from README.Rmd. Please edit that file -->

# JointODE

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/ziyangg98/JointODE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ziyangg98/JointODE/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ziyangg98/JointODE/graph/badge.svg)](https://app.codecov.io/gh/ziyangg98/JointODE)

<!-- badges: end -->

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

$$\ddot{m}_i(t) + 2 \xi_i \omega_i \dot{m}_i(t) + \omega_i^2 m_i(t) = k \omega_i^2 \mu_i(t)$$

where $m_i(t)$ is the biomarker value, $\dot{m}_i(t)$ is velocity (rate
of change), and $\ddot{m}_i(t)$ is acceleration. Subject-specific
dynamics are characterized by natural frequency $\omega_i$ and damping
ratio $\xi_i$, with $\mu_i(t)$ representing external forcing (e.g.,
treatment effects, covariates). Individual heterogeneity is captured
through random effects on these ODE parameters.

**Survival Model:** The hazard function incorporates biomarker dynamics:

$$\lambda_i(t) = \lambda_{0}(t)\exp\left[\alpha_1 m_i(t) + \alpha_2 \dot{m}_i^{\gamma}(t) + \mathbf{W}_i^{\top}\boldsymbol{\phi}\right]$$

where $\gamma \in \{0, 1, 2\}$ controls the velocity power (0: no
velocity, 1: linear, 2: quadratic).

### Key Features

- **Subject-specific dynamics**: Random effects on ODE parameters
  capture individual heterogeneity
- **Flexible hazard**: B-spline baseline hazard with biomarker value and
  velocity effects
- **Efficient computation**: C++ implementation with parallel processing
  support

For full mathematical details, see the [technical
documentation](http://gongziyang.com/JointODE/articles/technical-details.html).

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
#> Elapsed: 2291.0 s
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
#> -38935.412 -38852.954  19492.706
#>
#> Coefficients:
#> Longitudinal Process: Second-Order ODE Model
#>             Estimate Std. Error z value Pr(>|z|)
#> biomarker    -1.0847     0.0059  -183.1   <2e-16 ***
#> velocity     -0.8122     0.0047  -174.2   <2e-16 ***
#> (Intercept)   0.0004     0.0009     0.4    0.689
#> x1            0.5418     0.0033   166.7   <2e-16 ***
#> x2           -0.4884     0.0029  -167.8   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> ODE System Characteristics:
#>                    Estimate Std. Error z value Pr(>|z|)
#> T (period)           6.0330     0.0165   366.2   <2e-16 ***
#> xi (damping ratio)   0.3899     0.0018   216.4   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Survival Process: Proportional Hazards Model
#>         Estimate Std. Error z value Pr(>|z|)
#> alpha_1   0.7226     0.1885   3.834 0.000126 ***
#> alpha_2   1.9644     0.7780   2.525 0.011567 *
#> w1        0.6876     0.1246   5.521 3.37e-08 ***
#> w2       -1.3422     0.2860  -4.693 2.69e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Baseline Hazard: B-spline with 3 basis functions
#> (Coefficients range: [-4.916, -1.955] )
#>
#> Initial State: Population Mean
#>    Estimate Std. Error z value Pr(>|z|)
#> m0  -0.5116     0.0030 -173.13   <2e-16 ***
#> v0  -0.1080     0.0051  -21.06   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Variance Components:
#> Measurement Error SD: 0.098762
#> Random Effect Covariance Matrix:
#>            [,1]       [,2]       [,3]      [,4]
#> [1,]  0.0108201  0.0314779 -0.0003538 -0.001884
#> [2,]  0.0314779  0.0915955 -0.0009834 -0.004753
#> [3,] -0.0003538 -0.0009834  0.0005906  0.001785
#> [4,] -0.0018843 -0.0047530  0.0017848  0.027303
#>
#> Model Diagnostics:
#> C-index (Concordance): 0.671
#> Convergence: Converged after 120 iterations

# Plot results
plot(fit)
#> `geom_smooth()` using formula = 'y ~ x'
#> `geom_smooth()` using formula = 'y ~ x'
#> `geom_smooth()` using formula = 'y ~ x'
```

<img src="man/figures/README-output-1.png" alt="" width="100%" />

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
- **Technical Details**: See `vignette("technical-details")` for
  mathematical formulations
- **Model Comparison**: See `vignette("comparison")` for comparisons
  with traditional joint models

## Code of Conduct

Please note that the JointODE project is released with a [Contributor
Code of Conduct](http://gongziyang.com/JointODE/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
