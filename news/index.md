# Changelog

## JointODE 0.2.0

- Migrated estimation from EM algorithm to TMB Laplace approximation
  (~20x faster)
- Two-phase warm start for subject-specific ODE dynamics random effects
- [`predict()`](https://rdrr.io/r/stats/predict.html) supports custom
  time grids and new subjects
- Time-dependent Cox model for hazard initialization
- ODE period and damping ratio in
  [`summary()`](https://rdrr.io/r/base/summary.html)
- PBC vignette with real data analysis

## JointODE 0.1.2

### Major Features

- **Comprehensive Visualization System**: Implemented complete plotting
  framework with
  [`plot.JointODE()`](https://gongziyang.com/JointODE/reference/plot.JointODE.md)
  method
  - Model diagnostics: overview panels, biomarker/velocity trajectories,
    phase space diagrams
  - Survival analysis: Kaplan-Meier curves, hazard contributions
  - Residual diagnostics: standardized residuals, Q-Q plots
  - All plots built with ggplot2 for publication-ready graphics
- **Initial State Optimization**: Added iterative optimization for
  estimating initial biomarker and velocity values
  - Uses CppAD automatic differentiation for efficient gradient
    computation
  - Improves model convergence and parameter estimation accuracy
  - Available in both
    [`JointODE()`](https://gongziyang.com/JointODE/reference/JointODE.md)
    and
    [`MarginalODE()`](https://gongziyang.com/JointODE/reference/MarginalODE.md)
    functions
- **Enhanced Prediction**: Extended
  [`predict()`](https://rdrr.io/r/stats/predict.html) method with
  improved covariate handling via LOCF (Last Observation Carried
  Forward)

### Package Quality Improvements

- Fixed R CMD check warnings for CRAN compliance
  - Renamed `src/utils.hpp` to `src/utils.h` following R package
    conventions
  - Added `.covrignore` to exclude third-party libraries from test
    coverage
- Optimized documentation examples (72s → 1.5s, ~98% improvement)
- Fixed NEWS.md format for proper changelog display

## JointODE 0.1.1

- Implemented subject-specific random effects on ODE acceleration
  parameters to account for population heterogeneity
- Added Laplace approximation for efficient posterior computation of
  random effects
- Integrated AGHQ (Adaptive Gauss-Hermite Quadrature) for accurate
  numerical integration
- Implemented CppAD automatic differentiation for efficient gradient
  computation

## JointODE 0.1.0

- Joint modeling of longitudinal biomarkers and survival outcomes using
  ODEs
- Second-order differential equation formulation for biomarker dynamics
- EM algorithm for parameter estimation
- Parallel processing support for improved computational efficiency

------------------------------------------------------------------------

### Development Notes

#### Subject-level identifiability

Real-data experiments showed that full second-order heterogeneity can be
weakly identified for subjects with short, nearly linear, or abrupt
trajectories. In the current model this appears as flat directions among
initial velocity, dynamic log-coefficients, and forcing terms. These
cases are not necessarily data errors; they indicate that the observed
window may not contain enough curvature to estimate every
subject-specific ODE component.

The package therefore keeps the main model simple and explicit:

``` r

observed ~ biomarker + velocity + covariates +
  (1 + biomarker + velocity | id)
```

where `biomarker` is `log_omega2`, `velocity` is `log_2xi_omega`, and
`1` inside the random-effects term is the subject-specific forcing
intercept. Future work should focus on principled model-selection and
identifiability diagnostics for deciding when the full correlated
covariance structure is supported by a given biomarker.

### Version 0.2.0 (Planned - Q2 2025)

#### Subgroup Heterogeneity

- `JointODE_group()`: Latent subgroup modeling
- Group-specific ODE parameters (κ_g, γ_g)
- Model selection via ICL and entropy criteria
- K-means initialization with stability analysis

### Version 0.3.0 (Planned - Q4 2025)

#### Multiple Biomarkers

- `JointODE_multi()`: Multi-marker joint modeling
- `select_biomarkers()`: SIP-based variable selection
- Adaptive LASSO and group penalties
- Rcpp integration for performance

### Version 0.4.0 (Planned - Q4 2026)

#### Machine Learning

- `JointODE_nn()`: Neural ODE backend
- Python bridge via reticulate
- Cloud deployment tools
