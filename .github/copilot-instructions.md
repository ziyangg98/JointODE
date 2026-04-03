# Project Guidelines

JointODE: longitudinal + survival joint model via ODE (matrix exponential) + Cox PH. MCEM estimation with MCMCpack (random-walk Metropolis, Laplace-calibrated proposal). See [CLAUDE.md](../CLAUDE.md) for architecture.

## Build and Test

```bash
R CMD INSTALL --no-docs --no-multiarch .   # Fast C++ rebuild
Rscript -e 'testthat::test_local()'        # Run all tests
Rscript -e 'Rcpp::compileAttributes()'     # After changing [[Rcpp::export]]
```

## Conventions

- 2-space indent, 80-char limit (`.lintr`). snake_case functions, CamelCase classes.
- `stop()` directly — never silent fallback (`tryCatch`, `suppressWarnings`, `ginv`).
- C++ AD branch classification **before** `CppAD::Independent()` using `double` values.
- `R/RcppExports.R`, `src/RcppExports.cpp` auto-generated — don't edit.
- Test data: bundled `sim` dataset. Ad-hoc scripts go in `scripts/`.
- 如无必要勿增实体 — no speculative abstractions, helpers, or error handling.
