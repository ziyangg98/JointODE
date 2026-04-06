## R CMD check results

0 errors | 0 warnings | 1 note

- NOTE: `std::cout`, `puts`, `putchar` found in compiled code.
  These are from the vendored CppAD library (`inst/include/cppad/`),
  not from package code. CppAD's error handler uses these internally.

## Test environments

- local macOS (aarch64-apple-darwin), R 4.4.x
- GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release, oldrel-1)

## Downstream dependencies

This is a new package with no downstream dependencies.
