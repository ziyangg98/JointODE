# Predict Method for MarginalODE Objects

Computes predicted biomarker trajectories, velocities, and accelerations
for subjects based on the fitted marginal ODE model.

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

  An object of class `MarginalODE`

- newdata:

  Optional data frame with new subjects. If NULL, uses the training data
  from the model fit.

- times:

  Optional time points for prediction. Can be:

  - A numeric vector for the same times across all subjects

  - A named list with subject-specific time vectors

  - NULL to use observation times from the data (default)

- parallel:

  Logical flag for parallel computation (default: FALSE)

- n_cores:

  Number of cores for parallel processing. If 0, automatically detects
  available cores (default: 0)

- ...:

  Additional arguments (currently unused)

## Value

A data.frame with columns:

- `id`:

  Subject identifier

- `time`:

  Time points for predictions

- `biomarker`:

  Predicted biomarker values

- `velocity`:

  Predicted velocity (first derivative)

- `acceleration`:

  Predicted acceleration (second derivative)
