# Predict Method for MarginalODE Objects

Predict Method for MarginalODE Objects

## Usage

``` r
# S3 method for class 'MarginalODE'
predict(object, newdata = NULL, times = NULL, ...)
```

## Arguments

- object:

  A MarginalODE object

- newdata:

  Not yet supported

- times:

  Prediction times (NULL = observed)

- ...:

  Additional arguments

## Value

A data.frame with id, time, biomarker, velocity, acceleration
