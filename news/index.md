# Changelog

## JointODE 0.1.2

### Major Features

- **Comprehensive Visualization System**: Implemented complete plotting
  framework with
  [`plot.JointODE()`](http://gongziyang.com/JointODE/reference/plot.JointODE.md)
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
    [`JointODE()`](http://gongziyang.com/JointODE/reference/JointODE.md)
    and
    [`MarginalODE()`](http://gongziyang.com/JointODE/reference/MarginalODE.md)
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

### Development Roadmap

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
