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
  residual = c("gaussian", "student_t"),
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

- residual:

  Residual distribution, either `"gaussian"` or `"student_t"`.
  Student-t residuals estimate the degrees of freedom.

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
  data = sim$data$longitudinal_data
)
} # }
```
