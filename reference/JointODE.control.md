# Control Parameters for JointODE

Construct control parameters for the JointODE optimization algorithm.
This function can be called with no arguments to get defaults, or can
process a list to fill in missing values with defaults.

## Usage

``` r
JointODE.control(
  maxit = 200,
  tol = 0.001,
  verbose = FALSE,
  parallel = FALSE,
  n_cores = 0,
  quad_level = 3,
  hazard_quadrature = 1,
  .list = NULL,
  ...
)
```

## Arguments

- maxit:

  Maximum number of EM iterations (default: 200)

- tol:

  Convergence tolerance. The EM algorithm converges when
  max\|theta_new - theta_old\| \< tol (default: 1e-3)

- verbose:

  Logical or numeric; controls verbosity level. FALSE/0 for silent,
  TRUE/1 for basic progress, 2 for detailed output (default: FALSE)

- parallel:

  Logical; whether to use parallel computation (default: FALSE)

- n_cores:

  Integer; number of cores to use for parallel computation. If 0, uses
  all available cores (default: 0)

- quad_level:

  Integer; quadrature level for numerical integration (default: 3)

- hazard_quadrature:

  Integer; number of Simpson sub-intervals per observation interval for
  hazard integration (default: 1)

- .list:

  Optional list of control parameters to process

- ...:

  Additional control parameters

## Value

A list of control parameters with all defaults filled in

## See also

[`JointODE`](https://gongziyang.com/JointODE/reference/JointODE.md)

## Examples

``` r
# Default settings
control <- JointODE.control()

# Custom settings for faster exploration
control <- JointODE.control(maxit = 30, tol = 1e-3)

# Verbose output for debugging
control <- JointODE.control(verbose = TRUE)

# Parallel computation
control <- JointODE.control(parallel = TRUE, n_cores = 4)

# Process an existing list
my_list <- list(maxit = 200)
control <- JointODE.control(.list = my_list)
```
