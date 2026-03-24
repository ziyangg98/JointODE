# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is JointODE

An R package for joint modeling of longitudinal biomarker trajectories and time-to-event outcomes using ODEs. Biomarker dynamics follow a second-order linear ODE solved via matrix exponential (not Runge-Kutta), linked to a Cox proportional hazards model with B-spline baseline hazard. Estimation uses an EM algorithm with CppAD automatic differentiation.

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

`JointODE()` in `R/JointODE.R` is the main entry point. It runs an EM loop:

1. **E-step** (`R/posterior.R`): For each subject, find posterior mode of random effects via Laplace approximation. Calls `compute_logpost_cppad()` (C++).
2. **M-step** (`R/JointODE.R`): Update fixed effects by minimizing negative log-likelihood. Calls `compute_objective_cppad()` (C++). Variance components updated in closed form.
3. **State optimization**: Periodically re-optimizes initial ODE states `[m(0), v(0)]` via `compute_state_loglik_cppad()` (C++).

### C++ Layer (`src/`)

- **`solver.h`** — Header-only: all shared logic including `SubjectData`/`ODEParams` structs, ODE solver (`ode_step`, `ode_solve`), joint log-likelihood (`compute_joint_loglik`), B-spline basis, AD tape evaluation, hazard integration via RK45.
- **`exports.cpp`** — 4 `Rcpp::export` functions for the joint model. Each builds a CppAD tape with different AD independent variables (fixed effects theta, random effects b, or initial states).
- **`marginal.cpp`** — Standalone module for `MarginalODE` (longitudinal-only model, no survival component).

Three AD tapes share the same ODE solver and likelihood code; they differ only in **which variable is the AD independent variable**. See `src/DESIGN.md` for full mathematical details.

### ODE Solver: Matrix Exponential

The ODE is solved analytically via Cayley-Hamilton decomposition of exp(A*dt), with 5 branches based on discriminant of the characteristic equation (REAL, COMPLEX, REPEATED, FIRST_ORD, ZERO). Branch classification happens **before** `CppAD::Independent()` using double values — this is critical for correct AD taping.

### Key R Modules

- `R/process.R` — Data preprocessing: parses formulas, builds per-subject data lists for C++
- `R/initial.R` — Default parameter construction and initial value computation
- `R/control.R` — `JointODE.control()` for EM algorithm settings
- `R/utils.R` — Formula parsing, shared utilities
- `R/validate.R` — Input validation
- `R/MarginalODE.R` — Longitudinal-only model (no survival)

## Conventions

- Commit messages: English with emoji prefix (e.g., `refactor:`, `fix:`, `feat:`)
- C++ uses CppAD for AD and RcppArmadillo for linear algebra; `LinkingTo: Rcpp, RcppArmadillo`
- Errors must `stop()` directly — never use `tryCatch` to silently fall back
- testthat edition 3; test helpers in `tests/testthat/helper-*.R`
