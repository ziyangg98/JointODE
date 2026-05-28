
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

$$\ddot{m}_i(t) + 2 \xi_i \omega_i \dot{m}_i(t) + \omega_i^2 m_i(t) = f_i(t)$$

with initial conditions $m_i(0) = m_{0,i}$ and $\dot{m}_i(0)=v_{0,i}$.
Here $\omega_i$ is the natural frequency, $\xi_i$ is the damping ratio,
and $f_i(t)$ is covariate-driven forcing. The dynamic parameters are
estimated on log-coefficient scales: `biomarker` is $\log \omega_i^2$
and `velocity` is $\log(2\xi_i\omega_i)$. Individual heterogeneity is
captured through random effects on initial biomarker level, initial
velocity, ODE parameters, and optional forcing terms.

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
    (1 + biomarker + velocity | id),
  survival_formula = Surv(time, status) ~ w1 + w2,
  longitudinal_data = longitudinal_data,
  survival_data = sim$data$survival_data,
  init = "marginal"
)
cat(sprintf("Elapsed: %.1f s\n", (proc.time() - t0)["elapsed"]))
#> Elapsed: 87.4 s
```

``` r
# Model summary
summary(fit)
#>
#> Call:
#> JointODE(longitudinal_formula = observed ~ biomarker + velocity +
#>     x1 + x2 + (1 + biomarker + velocity | id), survival_formula = Surv(time,
#>     status) ~ w1 + w2, longitudinal_data = longitudinal_data,
#>     survival_data = sim$data$survival_data, init = "marginal")
#>
#> Data Descriptives:
#> Longitudinal Process            Survival Process
#> Number of Observations: 17222   Number of Events: 65 (32%)
#> Number of Subjects: 200
#>
#>        AIC        BIC     logLik
#> -27481.997 -27379.749  13771.999
#>
#> Coefficients:
#> Longitudinal Process: Second-Order ODE Model
#>               Estimate Std. Error z value Pr(>|z|)
#> log_omega2      0.0903     0.0228   3.965 7.34e-05 ***
#> log_2xi_omega  -0.1420     0.0253  -5.610 2.02e-08 ***
#> (Intercept)    -0.0091     0.0144  -0.632    0.527
#> x1              0.5441     0.0164  33.215  < 2e-16 ***
#> x2             -0.4821     0.0151 -31.942  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> ODE System Characteristics:
#>                           Estimate Std. Error z value Pr(>|z|)
#> omega (natural frequency)   1.0462     0.0119   87.84   <2e-16 ***
#> xi (damping ratio)          0.4147     0.0109   37.98   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Survival Process: Proportional Hazards Model
#>          Estimate Std. Error z value Pr(>|z|)
#> value      1.1387     0.2090   5.449 5.07e-08 ***
#> velocity   2.5307     0.7130   3.550 0.000386 ***
#> w1         0.9854     0.1621   6.078 1.22e-09 ***
#> w2        -1.0904     0.2726  -4.001 6.32e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Baseline Hazard: B-spline with 4 basis functions
#> (Coefficients range: [-5.532, -2.453] )
#>
#> Initial State: Population Mean
#>                   Estimate Std. Error z value Pr(>|z|)
#> initial_biomarker  -0.5069     0.0079 -63.785   <2e-16 ***
#> initial_velocity    0.0169     0.0110   1.534    0.125
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Variance Components:
#> Measurement Error SD: 0.099639
#> Random Effect Covariance Matrix:
#>                     initial_biomarker initial_velocity log_omega2 log_2xi_omega
#> initial_biomarker           0.0099902        0.0003266 -0.0017520    -0.0012379
#> initial_velocity            0.0003266        0.0108382 -0.0067544     0.0002496
#> log_omega2                 -0.0017520       -0.0067544  0.0703016     0.0002162
#> log_2xi_omega              -0.0012379        0.0002496  0.0002162     0.0824629
#> forcing_(Intercept)         0.0005011        0.0019067  0.0057967     0.0074494
#>                     forcing_(Intercept)
#> initial_biomarker             0.0005011
#> initial_velocity              0.0019067
#> log_omega2                    0.0057967
#> log_2xi_omega                 0.0074494
#> forcing_(Intercept)           0.0333415
#>
#> Model Diagnostics:
#> C-index (Concordance): 0.727
#> Convergence: Converged (relative convergence (4))

# Plot results
plot(fit)
```

<img src="man/figures/README-output-1.png" alt="" width="100%" />

The formula uses two reserved keywords (not data columns):

- **`biomarker`**: includes the latent $\log \omega_i^2$ parameter
- **`velocity`**: includes the latent $\log(2\xi_i\omega_i)$ parameter
- **`(1 + biomarker + velocity | id)`**: adds a subject-specific forcing
  intercept and random effects on these ODE parameters
- **`x1`, `x2`**: standard covariates driving the forcing function
  $f(t)$

## Learn More

- **PBC example workflow**: See `vignette("examples")`
- **Model illustration**: See `vignette("illustration")`
- **Technical details**: See `vignette("technical-details")`

## Code of Conduct

Please note that the JointODE project is released with a [Contributor
Code of Conduct](http://gongziyang.com/JointODE/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
