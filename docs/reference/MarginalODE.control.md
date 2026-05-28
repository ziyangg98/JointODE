# Control Parameters for MarginalODE

Construct control parameters for the MarginalODE optimization.

## Usage

``` r
MarginalODE.control(
  maxit = 200,
  tol = 1e-04,
  verbose = FALSE,
  parallel = FALSE,
  n_cores = 0,
  .list = NULL
)
```

## Arguments

- maxit:

  Maximum number of `nlminb` iterations (default: 200)

- tol:

  Relative convergence tolerance passed to `nlminb` (default: 1e-4)

- verbose:

  Logical or numeric; FALSE/0 for silent, TRUE/1 for basic progress, 2
  for detailed output (default: FALSE)

- parallel:

  Logical; whether to use parallel computation (default: FALSE)

- n_cores:

  Integer; number of cores (0 = auto) (default: 0)

- .list:

  Optional list of control parameters to process

## Value

A list of control parameters with all defaults filled in

## See also

[`MarginalODE`](https://gongziyang.com/JointODE/reference/MarginalODE.md)

## Examples

``` r
control <- MarginalODE.control()
control <- MarginalODE.control(maxit = 50, tol = 1e-4)
```
