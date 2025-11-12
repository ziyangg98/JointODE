# Summary Method for MarginalODE Objects

Summary Method for MarginalODE Objects

## Usage

``` r
# S3 method for class 'MarginalODE'
summary(object, ...)
```

## Arguments

- object:

  An object of class `MarginalODE`

- ...:

  Additional arguments (currently unused)

## Value

A list of class `summary.MarginalODE` containing:

- call:

  The model call

- coefficients:

  Coefficient table for longitudinal ODE parameters

- sigma:

  Named vector with measurement error SD

- nobs:

  Number of subjects

- n_observations:

  Total number of observations

- AIC:

  Akaike Information Criterion

- BIC:

  Bayesian Information Criterion

- logLik:

  Log-likelihood

- convergence:

  Convergence information
