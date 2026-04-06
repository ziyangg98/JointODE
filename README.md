
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
  init = "marginal"
)
#> Warning in stats::nlminb(start = obj$par, objective = obj$fn, gradient =
#> obj$gr, : NA/NaN function evaluation
cat(sprintf("Elapsed: %.1f s\n", (proc.time() - t0)["elapsed"]))
#> Elapsed: 110.8 s
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
#> Number of Observations: 17339   Number of Events: 59 (30%)
#> Number of Subjects: 200
#>
#>        AIC        BIC     logLik
#> -28960.371 -28874.614  14506.185
#>
#> Coefficients:
#> Longitudinal Process: Second-Order ODE Model
#>             Estimate Std. Error  z value Pr(>|z|)
#> biomarker    -1.0897     0.0078 -140.137   <2e-16 ***
#> velocity     -0.8449     0.0179  -47.094   <2e-16 ***
#> (Intercept)   0.0006     0.0014    0.392    0.695
#> x1            0.5437     0.0039  137.755   <2e-16 ***
#> x2           -0.4896     0.0036 -134.951   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> ODE System Characteristics:
#>                    Estimate Std. Error z value Pr(>|z|)
#> T (period)           6.0191     0.0215  280.27   <2e-16 ***
#> xi (damping ratio)   0.4047     0.0083   48.54   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Survival Process: Proportional Hazards Model
#>         Estimate Std. Error z value Pr(>|z|)
#> alpha_1   0.8963     0.2166   4.138 3.50e-05 ***
#> alpha_2   2.3209     0.6466   3.590 0.000331 ***
#> w1        0.6563     0.1336   4.912 9.02e-07 ***
#> w2       -0.8894     0.2761  -3.222 0.001275 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Baseline Hazard: B-spline with 4 basis functions
#> (Coefficients range: [-4.705, -2.450] )
#>
#> Initial State: Population Mean
#>           Estimate Std. Error z value Pr(>|z|)
#> biomarker  -0.4887     0.0073  -67.05   <2e-16 ***
#> velocity   -0.1146     0.0087  -13.15   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Variance Components:
#> Measurement Error SD: 0.099619 (SE: 0.000547)
#> Random Effect Covariance Matrix:
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  8.300e-03  2.218e-03  5.519e-04 -4.028e-05
#> [2,]  2.218e-03  7.116e-03 -5.438e-05  3.002e-04
#> [3,]  5.519e-04 -5.438e-05  1.338e-03  3.531e-04
#> [4,] -4.028e-05  3.002e-04  3.531e-04  4.206e-02
#>
#> Model Diagnostics:
#> C-index (Concordance): 0.613
#> Convergence: Converged (relative convergence (4))

# Plot results
plot(fit)
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
