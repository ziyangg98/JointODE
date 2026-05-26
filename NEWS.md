<!-- markdownlint-disable MD025 -->

# JointODE 0.2.0

* Migrated estimation from EM algorithm to TMB Laplace approximation (~20x faster)
* Two-phase warm start for subject-specific ODE dynamics random effects
* `predict()` supports custom time grids and new subjects
* Time-dependent Cox model for hazard initialization
* ODE period and damping ratio in `summary()`
* PBC vignette with real data analysis

# JointODE 0.1.2

## Major Features

* **Comprehensive Visualization System**: Implemented complete plotting framework with `plot.JointODE()` method
  * Model diagnostics: overview panels, biomarker/velocity trajectories, phase space diagrams
  * Survival analysis: Kaplan-Meier curves, hazard contributions
  * Residual diagnostics: standardized residuals, Q-Q plots
  * All plots built with ggplot2 for publication-ready graphics

* **Initial State Optimization**: Added iterative optimization for estimating initial biomarker and velocity values
  * Uses CppAD automatic differentiation for efficient gradient computation
  * Improves model convergence and parameter estimation accuracy
  * Available in both `JointODE()` and `MarginalODE()` functions

* **Enhanced Prediction**: Extended `predict()` method with improved covariate handling via LOCF (Last Observation Carried Forward)

## Package Quality Improvements

* Fixed R CMD check warnings for CRAN compliance
  * Renamed `src/utils.hpp` to `src/utils.h` following R package conventions
  * Added `.covrignore` to exclude third-party libraries from test coverage
* Optimized documentation examples (72s → 1.5s, ~98% improvement)
* Fixed NEWS.md format for proper changelog display

# JointODE 0.1.1

* Implemented subject-specific random effects on ODE acceleration parameters to account for population heterogeneity
* Added Laplace approximation for efficient posterior computation of random effects
* Integrated AGHQ (Adaptive Gauss-Hermite Quadrature) for accurate numerical integration
* Implemented CppAD automatic differentiation for efficient gradient computation

# JointODE 0.1.0

* Joint modeling of longitudinal biomarkers and survival outcomes using ODEs
* Second-order differential equation formulation for biomarker dynamics
* EM algorithm for parameter estimation
* Parallel processing support for improved computational efficiency

---

## Development Roadmap

### Future Modeling Considerations

#### Subject-level identifiability and future ODE parameterization

Real-data subject-level diagnostics showed that the original target-level
parameterization can create weakly identified directions when full
second-order individual heterogeneity is allowed. In the target form,
`tau_i m_i''(t) + m_i'(t) = lambda_i * (mu_i - m_i(t))`, the data often identify
the product `lambda_i * mu_i` more reliably than `lambda_i` and `mu_i`
separately. When `lambda_i` is small, `mu_i` can drift to unrealistic values
with little change in the fitted trajectory, producing nearly flat Hessian
directions.

The more identifiable parameterization is the forcing form
`tau_i m_i''(t) + m_i'(t) = eta_i(t) - lambda_i m_i(t)`, where
`eta_i(t) = lambda_i * mu_i(t)` is modeled directly and `mu_i(t)` is treated as
a derived quantity. This separates the statistically visible forcing term from
the less stable target-level interpretation.

After this reparameterization, the main remaining weak-identification pattern
in RDW and several comparison biomarkers was `lambda_i` near zero. These
subjects are better interpreted as weak-restoring or drift-like trajectories:
`tau_i m_i''(t) + m_i'(t) ~= eta_i(t)`. They are not necessarily data errors;
their observed time windows do not contain enough information to estimate a
subject-specific relaxation rate. Future model development should therefore
consider a zero-inflated or hurdle relaxation-rate structure, for example
`lambda_i ~ pi * delta_0 + (1 - pi) * LogNormal(...)`.

Subject-specific `tau_i` remains scientifically useful for full second-order
heterogeneity, but diagnostics also showed possible coupling among `lambda_i`,
`tau_i`, and initial velocity. This coupling arises because the equation can be
written with coefficients `1 / tau_i` and `lambda_i / tau_i`, while early-time
Taylor expansion depends on the combined acceleration term
`eta_i - lambda_i * m0_i - v0_i`. Future implementations should include
Hessian-based subject-level identifiability diagnostics before relying on
unrestricted full covariance structures for these random effects.

Finally, bilirubin-like biomarkers with abrupt spikes or non-smooth changes may
represent model mismatch rather than only boundary relaxation rates. These
cases may require time-varying or event-driven forcing terms, such as
low-dimensional spline forcing or clinically defined event indicators.

## Version 0.2.0 (Planned - Q2 2025)

### Subgroup Heterogeneity

* `JointODE_group()`: Latent subgroup modeling
* Group-specific ODE parameters (κ_g, γ_g)
* Model selection via ICL and entropy criteria
* K-means initialization with stability analysis

## Version 0.3.0 (Planned - Q4 2025)

### Multiple Biomarkers

* `JointODE_multi()`: Multi-marker joint modeling
* `select_biomarkers()`: SIP-based variable selection
* Adaptive LASSO and group penalties
* Rcpp integration for performance

## Version 0.4.0 (Planned - Q4 2026)

### Machine Learning

* `JointODE_nn()`: Neural ODE backend
* Python bridge via reticulate
* Cloud deployment tools
