
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

with $m_i(0) = m_{0,i}$ and a quasi-steady initial velocity
$\dot{m}_i(0)=\{f_i(0)-\omega_i^2m_{0,i}\}/(2\xi_i\omega_i)$. Here
$\omega_i$ is the natural frequency, $\xi_i$ is the damping ratio, and
$f_i(t)$ is covariate-driven forcing. The dynamic parameters are
estimated on log-coefficient scales: `biomarker` is $\log \omega_i^2$
and `velocity` is $\log(2\xi_i\omega_i)$. Individual heterogeneity is
captured through random effects on initial biomarker level, ODE
parameters, and optional forcing terms.

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
#> Elapsed: 87.2 s
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
#> Number of Observations: 16894   Number of Events: 67 (34%)
#> Number of Subjects: 200
#>
#>        AIC        BIC     logLik
#> -26764.725 -26682.267  13407.362
#>
#> Coefficients:
#> Longitudinal Process: Second-Order ODE Model
#>               Estimate Std. Error z value Pr(>|z|)
#> log_omega2      0.0949     0.0225   4.212 2.53e-05 ***
#> log_2xi_omega  -0.1369     0.0236  -5.793 6.93e-09 ***
#> (Intercept)    -0.0158     0.0146  -1.083    0.279
#> x1              0.5380     0.0165  32.652  < 2e-16 ***
#> x2             -0.4830     0.0152 -31.780  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> ODE System Characteristics:
#>                           Estimate Std. Error z value Pr(>|z|)
#> omega (natural frequency)   1.0486     0.0118   88.76   <2e-16 ***
#> xi (damping ratio)          0.4158     0.0100   41.55   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Survival Process: Proportional Hazards Model
#>          Estimate Std. Error z value Pr(>|z|)
#> value      0.7728     0.1808   4.273 1.93e-05 ***
#> velocity   1.7008     0.3848   4.420 9.85e-06 ***
#> w1         0.6337     0.1464   4.329 1.49e-05 ***
#> w2        -0.9578     0.2672  -3.585 0.000337 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Baseline Hazard: B-spline with 4 basis functions
#> (Coefficients range: [-5.253, -2.394] )
#>
#> Initial State: Population Mean
#>                   Estimate Std. Error z value Pr(>|z|)
#> initial_biomarker  -0.5016     0.0075  -66.98   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Variance Components:
#> Measurement Error SD: 0.099662
#> Random Effect Covariance Matrix:
#>                     initial_biomarker log_omega2 log_2xi_omega
#> initial_biomarker           0.0087971  -0.006464     0.0001618
#> log_omega2                 -0.0064640   0.075953     0.0070495
#> log_2xi_omega               0.0001618   0.007049     0.0862564
#> forcing_(Intercept)         0.0012301   0.001358     0.0043238
#>                     forcing_(Intercept)
#> initial_biomarker              0.001230
#> log_omega2                     0.001358
#> log_2xi_omega                  0.004324
#> forcing_(Intercept)            0.034458
#>
#> Model Diagnostics:
#> C-index (Concordance): 0.631
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
