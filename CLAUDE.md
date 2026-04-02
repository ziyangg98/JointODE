# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is JointODE

An R package for joint modeling of longitudinal biomarker trajectories and time-to-event outcomes using ODEs. Biomarker dynamics follow a second-order linear ODE solved via matrix exponential (Cayley-Hamilton, not Runge-Kutta), linked to a Cox proportional hazards model with B-spline baseline hazard. Estimation uses an EM algorithm with CppAD automatic differentiation.

## Build & Test Commands

```bash
# Full R CMD check (includes build, tests, examples, vignettes)
R CMD check .

# Build and install
R CMD INSTALL .

# Run all tests
Rscript -e 'testthat::test_local()'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test-ode.R")'

# Regenerate Rcpp exports after changing C++ [[Rcpp::export]] functions
Rscript -e 'Rcpp::compileAttributes()'

# Regenerate roxygen docs and NAMESPACE
Rscript -e 'roxygen2::roxygenise()'

# Rebuild C++ only (faster iteration)
R CMD INSTALL --no-docs --no-multiarch .
```

## Architecture

### EM Algorithm Flow (R side)

`JointODE()` in `R/JointODE.R` is the main entry point. It runs a PQL-ECM loop:

1. **E-step** (`R/posterior.R`): For each subject, find posterior mode of random effects via Laplace approximation (`nlm`). Calls `.compute_joint_logpost()` (C++). Uses `.safe_chol()` (Gill-Murray) when the posterior Hessian is not positive definite.
2. **M-step** (`R/posterior.R`): Trust region optimization (`trust::trust`) of Laplace marginal likelihood. The C++ objective (`compute_joint_objective`) includes `-0.5·log|H_z(θ)|` via nested AD (`AD<AD<double>>`). Fixed effects θ = [baseline, hazard, longitudinal, initial_state]; σ_e and Σ_b are fixed during M-step.
3. **Variance updates** (`R/posterior.R`): Laplace-corrected closed form — Σ_b = (1/n)Σ(b̂b̂ᵀ + H_z⁻¹), σ_e = sqrt(RSS/N).
4. **SEM vcov** (`R/posterior.R`): Post-convergence variance-covariance via supplemented EM (Louis formula). Numerical Jacobian of EM map.

### C++ Layer (`src/`)

- **`solver.h`** — Header-only: all shared logic including `SubjectData`/`ODEParams` structs, ODE solver (`ode_step`, `ode_solve_joint`, `ode_solve_marginal`), joint/marginal log-likelihood, B-spline basis, AD tape evaluation (`eval_tape`), and log-linear hazard integration.
- **`exports.cpp`** — 7 `Rcpp::export` functions: `compute_joint_objective` (M-step, AD over theta), `compute_joint_logpost` (E-step, AD over b), `compute_joint_state` (state opt, AD over initial states), `solve_batch_joint` (prediction, no AD), `compute_marginal_objective`, `compute_marginal_state`, `solve_batch_marginal`.
- **`cppad_stub.cpp`** — Minimal CppAD compilation unit required for linking.

CppAD headers are vendored in `inst/include/cppad/` and included via `PKG_CPPFLAGS = -I../inst/include` in `src/Makevars`.

Three AD tapes share the same ODE solver and likelihood code; they differ only in **which variable is the AD independent variable**. A fourth nested AD tape (`AD<AD<double>>`) computes the Laplace correction: inner tape over b yields H_z as `AD<double>` values tracked by the outer theta tape, enabling exact gradient propagation of -0.5·log|H_z(θ)| through `CppAD::LuSolve`.

### ODE Solver: Matrix Exponential

The ODE is solved analytically via Cayley-Hamilton decomposition of exp(A*dt), with 5 branches based on discriminant of the characteristic equation (`REAL`, `COMPLEX`, `REPEATED`, `FIRST_ORD`, `ZERO`). Branch classification happens **before** `CppAD::Independent()` using double values — this is critical for correct AD taping. Biomarker/velocity values are clamped after each step to prevent overflow (data-adaptive clamp = 5x max observed value).

### Hazard Integration

Composite Simpson's rule over configurable sub-intervals (`hazard_quadrature` control parameter) for cumulative hazard integration within each ODE time step.

### Key R Modules

- `R/posterior.R` — E-step (Laplace), M-step (trust region + nested AD Laplace correction), SEM vcov
- `R/utils.R` — Formula parsing, `.safe_chol` (Gill-Murray), parameter conversion, convergence tracking
- `R/process.R` — Data preprocessing: parses formulas, builds per-subject data lists for C++
- `R/initial.R` — Default parameter construction and `MarginalODE`-based initial value computation
- `R/finalize.R` — Post-convergence: C-index, AIC/BIC, SEM vcov call
- `R/control.R` — `JointODE.control()` / `MarginalODE.control()` for algorithm settings
- `R/validate.R` — Input validation
- `R/MarginalODE.R` — Longitudinal-only model (no survival)
- `R/state.R` — Per-subject initial state optimization via `nlm`
- `R/simulate.R` — `simulate.JointODE()` data generation with physical parameterization (damping ratio, period)
- `R/plot.R` — `plot.JointODE()` diagnostic plots (trajectories, survival, residuals)

### Test Data: `sim` Dataset

Tests use the bundled `sim` dataset (`data/sim.rda`), structured as:
- `sim$data$longitudinal_data` — data frame with columns: id, time, biomarker, velocity, observed, x1, x2
- `sim$data$survival_data` — data frame with columns: id, time, status, w1, w2
- `sim$data$state` — matrix (n x 2): initial biomarker and velocity per subject
- `sim$data$random_effects` — matrix (n x 2): dyn_value, dyn_slope random effects
- `sim$init` — initial parameter list with `$coefficients` and `$configurations`, suitable for direct use in model fitting and test helpers

## Conventions

- Commit messages: English with emoji prefix (e.g., `refactor:`, `fix:`, `feat:`)
- C++ uses CppAD for AD and RcppArmadillo for linear algebra; `LinkingTo: Rcpp, RcppArmadillo`
- Errors must `stop()` directly — never use `tryCatch` to silently fall back, `suppressWarnings`, `ginv`, or silent guards. Fix root causes in C++ when possible.
- testthat edition 3; test helpers in `tests/testthat/helper-*.R`
- Tests rely on the bundled `sim` dataset (`data/sim.rda`) for processed data, initial parameters, and random effects
- Ad-hoc test scripts (PBC, sim benchmarks) go in `scripts/`, not `tests/`

## Performance Baseline

```bash
# Generate baseline CSV (sequential + optional parallel)
Rscript scripts/perf-baseline.R --n=20 --reps=3 --maxit=10 --tol=1e-2 --out=perf-base.csv

# Generate candidate CSV after code changes
Rscript scripts/perf-baseline.R --n=20 --reps=3 --maxit=10 --tol=1e-2 --out=perf-new.csv

# Compare elapsed_mean; fail if slowdown > 10%
Rscript scripts/perf-compare.R --base=perf-base.csv --new=perf-new.csv --metric=elapsed_mean --fail_pct=10
```

## Formula Syntax

`longitudinal_formula`: `observed ~ biomarker + velocity + x1 + x2 + (biomarker + velocity | id)`
- `biomarker` and `velocity` are reserved names for ODE state variables (value and slope)
- Other terms (`x1`, `x2`) are external covariates affecting the forcing function
- `(biomarker + velocity | id)` specifies subject-specific random effects on ODE coefficients

`survival_formula`: `Surv(time, status) ~ w1 + w2` — standard Cox-style formula
