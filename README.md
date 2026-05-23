
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

with initial conditions $m_i(0) = m_{0,i}$ and $\dot{m}_i(0) = v_{0,i}$,
where $\omega_i$ is the natural frequency, $\xi_i$ is the damping ratio,
and $f_i(t)$ is covariate-driven forcing. Stability ($\omega > 0$,
$\xi > 0$) is enforced by estimating $\log(\omega^2)$ and
$\log(2\xi\omega)$. Individual heterogeneity is captured through random
effects on initial states $(m_{0,i}, v_{0,i})$ and multiplicative random
effects on ODE parameters.

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
  init = "marginal"
)
cat(sprintf("Elapsed: %.1f s\n", (proc.time() - t0)["elapsed"]))
#> Elapsed: 450.5 s
```

``` r
# Model summary
summary(fit)
#> 
#> Call:
#> JointODE(longitudinal_formula = observed ~ biomarker + velocity + 
#>     x1 + x2 + (biomarker + velocity | id), survival_formula = Surv(time, 
#>     status) ~ w1 + w2, longitudinal_data = longitudinal_data, 
#>     survival_data = sim$data$survival_data, init = "marginal")
#> 
#> Data Descriptives:
#> Longitudinal Process            Survival Process
#> Number of Observations: 17319   Number of Events: 58 (29%)
#> Number of Subjects: 200
#> 
#>        AIC        BIC     logLik
#> -28900.821 -28815.064  14476.410
#> 
#> Coefficients:
#> Longitudinal Process: Second-Order ODE Model
#>               Estimate Std. Error  z value Pr(>|z|)    
#> log_omega2      0.0872     0.0075   11.644   <2e-16 ***
#> log_2xi_omega  -0.1907     0.0202   -9.450   <2e-16 ***
#> (Intercept)    -0.0013     0.0015   -0.844    0.399    
#> x1              0.5477     0.0040  136.034   <2e-16 ***
#> x2             -0.4909     0.0037 -131.452   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> ODE System Characteristics:
#>                           Estimate Std. Error z value Pr(>|z|)    
#> omega (natural frequency)   1.0446     0.0039  267.04   <2e-16 ***
#> xi (damping ratio)          0.3955     0.0078   50.67   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Survival Process: Proportional Hazards Model
#>          Estimate Std. Error z value Pr(>|z|)    
#> value      0.9475     0.2138   4.432 9.35e-06 ***
#> velocity   2.1034     0.6429   3.272 0.001070 ** 
#> w1         0.7079     0.1358   5.211 1.88e-07 ***
#> w2        -0.9126     0.2759  -3.308 0.000941 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Baseline Hazard: B-spline with 4 basis functions
#> (Coefficients range: [-4.705, -2.446] )
#> 
#> Initial State: Population Mean
#>                   Estimate Std. Error z value Pr(>|z|)    
#> initial_biomarker  -0.4988     0.0073  -68.47   <2e-16 ***
#> initial_velocity   -0.1012     0.0093  -10.90   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Variance Components:
#> Measurement Error SD: 0.099369
#> Random Effect Covariance Matrix:
#>                   initial_biomarker initial_velocity log_omega2 log_2xi_omega
#> initial_biomarker         0.0083456        0.0008110  0.0002327     0.0008118
#> initial_velocity          0.0008110        0.0094564  0.0001412     0.0006467
#> log_omega2                0.0002327        0.0001412  0.0017059    -0.0004349
#> log_2xi_omega             0.0008118        0.0006467 -0.0004349     0.0514761
#> 
#> Model Diagnostics:
#> C-index (Concordance): 0.616
#> Convergence: Converged (outer parameter convergence; last M-step: relative convergence (4))

# Plot results
plot(fit)
```

<img src="man/figures/README-output-1.png" alt="" width="100%" />

The formula uses two reserved keywords (not data columns):

- **`biomarker`**: includes the frequency parameter $\omega^2$ in the
  ODE ($b_1 m(t)$ term)
- **`velocity`**: includes the damping parameter $2\xi\omega$ in the ODE
  ($b_2 \dot{m}(t)$ term)
- **`(biomarker + velocity | id)`**: adds subject-specific
  multiplicative random effects on these ODE parameters
- **`x1`, `x2`**: standard covariates driving the forcing function
  $f(t)$

## Learn More

- **Getting Started**: See `vignette("JointODE")` for a detailed
  tutorial
- **Technical Details**: See `vignette("technical-details")` for
  mathematical formulations
- **Model Comparison**: See `vignette("comparison")` for comparisons
  with traditional joint models

## Performance Baseline

Use the helper scripts in `scripts/` to generate reproducible runtime
baselines and compare two runs:

``` bash
# Generate baseline CSV (sequential + optional parallel case)
Rscript scripts/perf-baseline.R --n=20 --reps=3 --maxit=10 --tol=1e-2 --out=perf-base.csv

# Generate candidate CSV after code changes
Rscript scripts/perf-baseline.R --n=20 --reps=3 --maxit=10 --tol=1e-2 --out=perf-new.csv

# Compare elapsed_mean and fail if slowdown > 10%
Rscript scripts/perf-compare.R --base=perf-base.csv --new=perf-new.csv --metric=elapsed_mean --fail_pct=10
```

These scripts are intended for quick engineering checks, not
publication-grade benchmarking.

## Code of Conduct

Please note that the JointODE project is released with a [Contributor
Code of Conduct](http://gongziyang.com/JointODE/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
