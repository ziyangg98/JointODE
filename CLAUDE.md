# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is JointODE

An R package for joint modeling of longitudinal biomarker trajectories and time-to-event outcomes using ODEs. Biomarker dynamics follow a second-order linear ODE solved via matrix exponential (Cayley-Hamilton), linked to a Cox proportional hazards model with B-spline baseline hazard. Estimation uses TMB (Template Model Builder) for automatic differentiation and Laplace approximation of random effects.

## Build & Test Commands

```bash
# Full R CMD check
R CMD check .

# Rebuild C++ and install (fast iteration)
R CMD INSTALL --no-docs --no-multiarch .

# Clean rebuild (after C++ header changes)
rm -f src/*.o && R CMD INSTALL .

# Run all tests
Rscript -e 'testthat::test_local()'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test-fit.R")'

# Regenerate roxygen docs and NAMESPACE
Rscript -e 'roxygen2::roxygenise()'

# Build pkgdown site
Rscript -e 'pkgdown::build_site()'

# Build README from README.Rmd
Rscript -e 'devtools::build_readme()'
```

## Architecture

### TMB Estimation Flow

`JointODE()` in `R/JointODE.R` is the main entry point:

1. **Parse & validate** — Formula parsing (`.parse_longitudinal_formula`, `.parse_survival_formula` in `R/utils.R`), input validation (`R/validate.R`)
2. **Data processing** (`R/process.R`) — Build per-subject data lists, pack into flat TMB input vectors
3. **Model setup** (`R/setup.R`) — Compute dimensions, spline config, coefficient names, RE layout
4. **Initialization** (`R/initial.R`) — Either `init = "marginal"` (fits MarginalODE → time-dependent Cox → Weibull baseline) or user-provided parameters
5. **TMB optimization** — `MakeADFun(random = "random_effects", normalize = TRUE)` + `nlminb` with Laplace approximation
6. **Result extraction** (`R/finalize.R`) — `sdreport`, C-index, vcov, convergence

### MarginalODE Warm Start

When `MarginalODE()` has dynamics RE (`biomarker + velocity | id`), it uses two-phase fitting:
- **Phase 1**: Fit reduced model (no dynamics RE) — fast, inner Newton converges quickly
- **Phase 2**: Fit full model using Phase 1 estimates as starting values — inner Newton starts near mode

This is implemented in `.warm_start_marginal()` in `R/initial.R`. The same warm start is used by `.initialize_from_marginal()` for JointODE.

### C++ Layer (`src/`)

- **`src/include/utils.hpp`** — Shared: ODE solver (`ode_step` template, works with both AD and double via CondExp dispatch), B-spline basis (double only), hazard evaluation, time grid utilities, safe math (`safe_exp`, `clamp`, `safe_sqrt`)
- **`src/include/joint.hpp`** — Joint model TMB template: longitudinal + survival likelihood with Simpson's rule hazard integration
- **`src/include/marginal.hpp`** — Marginal (longitudinal-only) TMB template
- **`src/TMB/JointODE.hpp`** — Dispatch: `model_type == 1` → marginal, else → joint

### ODE Solver

Solved analytically via Cayley-Hamilton: exp(A·dt) decomposed using discriminant of characteristic equation. Uses `CppAD::CondExpGt` to select between REAL (D>0) and COMPLEX (D<0) branches without breaking the AD tape — both branches are always computed, selected via CondExp. This enables AD through branch changes when dynamics RE modify eigenvalues during optimization.

### R-side Prediction

`predict.JointODE()` and `predict.MarginalODE()` use pure R implementations (`.predict_trajectories`, `.predict_marginal_trajectories` in `R/utils.R`) with a 5-branch if/else ODE solver (`.ode_step_r`). No AD overhead. Supports custom time grids and per-subject observation times.

### Key R Modules

- `R/JointODE.R` — Main entry + S3 methods (summary, coef, vcov, logLik, predict, print)
- `R/MarginalODE.R` — Longitudinal-only model + S3 methods
- `R/process.R` — Data processing, TMB data/param packing, `.pack_correlation_theta()` helper
- `R/setup.R` — `.setup_model()`, `.setup_marginal_model()` for dimensions and config
- `R/initial.R` — `.default_parameters()`, `.initialize_from_marginal()`, `.warm_start_marginal()`
- `R/finalize.R` — `.finalize_joint()`, `.finalize_marginal()` for TMB result extraction
- `R/utils.R` — Formula parsing, constants, prediction helpers, `.build_counting_process()`
- `R/control.R` — `JointODE.control()` / `MarginalODE.control()`
- `R/validate.R` — Input validation
- `R/simulate.R` — `simulate.JointODE()` with physical parameterization (damping ratio, period)
- `R/plot.R` — `plot.JointODE()` diagnostic plots

### Test Data: `sim` Dataset

Tests use the bundled `sim` dataset (`data/sim.rda`):
- `sim$data$longitudinal_data` — columns: id, time, biomarker, velocity, observed, x1, x2
- `sim$data$survival_data` — columns: id, time, status, w1, w2
- `sim$data$random_effects` — matrix (n × 4): init_biomarker, init_velocity, dyn_biomarker, dyn_velocity
- `sim$init` — initial parameter list for direct use in `JointODE(init = sim$init)`

Note: `longitudinal_data` has reserved column names (`biomarker`, `velocity`). When passing to `JointODE()`, filter these out: `ld[, c("id", "time", "observed", "x1", "x2")]`

## Conventions

- Commit messages: English with prefix (e.g., `refactor:`, `fix:`, `feat:`)
- C++ uses TMB for AD; `LinkingTo: RcppEigen, TMB`
- Errors must `stop()` directly — never use `tryCatch` to silently fall back, `suppressWarnings`, `ginv`, or silent guards. Fix root causes in C++ when possible.
- testthat edition 3; test helpers in `tests/testthat/helper.R`
- Tests rely on the bundled `sim` dataset
- Ad-hoc test scripts go in `scripts/`, not `tests/`
- Verbose levels: 0 = silent, 1 = progress, 2 = outer iterations, 3 = inner Newton detail
