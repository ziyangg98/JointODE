
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
#> Elapsed: 279.8 s
```

``` r
# Model summary
summary(fit)
#> Warning in sqrt(diag(object$vcov)): NaNs produced
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
#> -31491.145 -31438.371  15761.572
#>
#> Coefficients:
#> Longitudinal Process: Second-Order ODE Model
#>             Estimate Std. Error z value Pr(>|z|)
#> biomarker    -1.0807     0.0044 -246.60   <2e-16 ***
#> velocity     -0.8312     0.0044 -188.22   <2e-16 ***
#> (Intercept)  -0.0005     0.0009   -0.57    0.568
#> x1            0.5406     0.0025  216.64   <2e-16 ***
#> x2           -0.4867     0.0023 -215.40   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> ODE System Characteristics:
#>                    Estimate Std. Error z value Pr(>|z|)
#> T (period)           6.0442     0.0123   493.2   <2e-16 ***
#> xi (damping ratio)   0.3998     0.0018   225.4   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Survival Process: Proportional Hazards Model
#>         Estimate Std. Error z value Pr(>|z|)
#> alpha_1   0.7243     0.1870   3.874 0.000107 ***
#> alpha_2   1.8728     0.6761   2.770 0.005606 **
#> w1        0.6862     0.1237   5.546 2.92e-08 ***
#> w2       -1.3423     0.2866  -4.684 2.81e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Baseline Hazard: B-spline with 3 basis functions
#> (Coefficients range: [-4.853, -1.996] )
#>
#> Variance Components:
#> Measurement Error SD: 0.098770
#> Random Effect Covariance Matrix:
#>          [,1]     [,2]
#> [1,] 0.001369 0.001522
#> [2,] 0.001522 0.030368
#>
#> Model Diagnostics:
#> C-index (Concordance): 0.668
#> Convergence: EM algorithm converged after 10 iterations

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
