
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
#> Elapsed: 124.0 s
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
#> -39111.219 -39018.866  19583.609
#>
#> Coefficients:
#> Longitudinal Process: Second-Order ODE Model
#>             Estimate Std. Error  z value Pr(>|z|)
#> biomarker    -1.0774     0.0059 -181.456   <2e-16 ***
#> velocity     -0.8139     0.0047 -174.335   <2e-16 ***
#> (Intercept)  -0.0001     0.0009   -0.147    0.883
#> x1            0.5386     0.0033  165.429   <2e-16 ***
#> x2           -0.4854     0.0029 -166.535   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> ODE System Characteristics:
#>                    Estimate Std. Error z value Pr(>|z|)
#> T (period)           6.0534     0.0167   362.9   <2e-16 ***
#> xi (damping ratio)   0.3921     0.0018   216.6   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Survival Process: Proportional Hazards Model
#>         Estimate Std. Error z value Pr(>|z|)
#> alpha_1   0.6792     0.1871   3.630 0.000284 ***
#> alpha_2   1.9311     1.0118   1.909 0.056309 .
#> w1        0.6767     0.1241   5.455 4.91e-08 ***
#> w2       -1.3368     0.2857  -4.680 2.87e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Baseline Hazard: B-spline with 6 basis functions
#> (Coefficients range: [-6.456, -0.531] )
#>
#> Initial State: Population Mean
#>    Estimate Std. Error z value Pr(>|z|)
#> m0  -0.5073     0.0030 -171.13   <2e-16 ***
#> v0  -0.0941     0.0051  -18.32   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Variance Components:
#> Measurement Error SD: 0.098863
#> Random Effect Covariance Matrix:
#>            [,1]       [,2]       [,3]      [,4]
#> [1,]  0.0105661  0.0315126 -0.0003927 -0.002762
#> [2,]  0.0315126  0.0949087 -0.0007375 -0.003633
#> [3,] -0.0003927 -0.0007375  0.0005088  0.001521
#> [4,] -0.0027615 -0.0036329  0.0015213  0.025505
#>
#> Model Diagnostics:
#> C-index (Concordance): 0.868
#> Convergence: Converged after 8 iterations

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
