# Marginal Second-Order ODE Parameter Estimation

Estimates population-level ODE parameters for longitudinal biomarker
trajectories: \\\ddot{m}(t) = \beta_1 m(t) + \beta_2 \dot{m}(t) +
X\beta\\

## Usage

``` r
MarginalODE(
  formula,
  data,
  time = "time",
  id = "id",
  state = NULL,
  control = list()
)
```

## Arguments

- formula:

  Response and covariates (e.g., `biomarker ~ x1 + x2`)

- data:

  Data frame with longitudinal measurements

- time:

  Time variable name (default: `"time"`)

- id:

  Subject identifier name (default: `"id"`)

- state:

  Optional \\n \times 2\\ matrix of initial conditions \\\[m(0),
  \dot{m}(0)\]\\. If `NULL`, estimated from data.

- control:

  List of control parameters. See
  [`MarginalODE.control`](https://gongziyang.com/JointODE/reference/MarginalODE.control.md).

## Value

S3 object of class `MarginalODE`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- MarginalODE(
  formula = observed ~ x1 + x2,
  data = sim$data$longitudinal_data,
  state = as.matrix(sim$data$state)
)
} # }
```
