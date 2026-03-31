# Predict Method for MarginalODE Objects

Predict Method for MarginalODE Objects

## Usage

``` r
# S3 method for class 'MarginalODE'
predict(
  object,
  newdata = NULL,
  times = NULL,
  parallel = FALSE,
  n_cores = 0,
  ...
)
```

## Arguments

- object:

  A MarginalODE object

- newdata:

  Not yet supported

- times:

  Prediction times (NULL = observed)

- parallel:

  Logical for parallel computation

- n_cores:

  Number of cores (0 = auto)

- ...:

  Additional arguments

## Value

A data.frame with id, time, biomarker, velocity, acceleration
