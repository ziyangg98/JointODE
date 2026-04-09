# Predict Method for JointODE Objects

Computes predicted biomarker trajectories, velocities, and accelerations
for subjects based on the fitted joint ODE model. Predictions are
obtained from the TMB REPORT output.

## Usage

``` r
# S3 method for class 'JointODE'
predict(object, newdata = NULL, times = NULL, ...)
```

## Arguments

- object:

  An object of class `JointODE`

- newdata:

  Optional data frame with new subjects. If NULL, uses the training data
  from the model fit.

- times:

  Optional time points for prediction. If NULL, uses observed time
  points for each subject.

- ...:

  Additional arguments (currently unused)

## Value

A data.frame with columns id, time, biomarker, velocity.
