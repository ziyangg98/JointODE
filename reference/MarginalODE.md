# Marginal Second-Order ODE Parameter Estimation

Estimates the longitudinal part of a second-order ODE model:
\$\$\ddot{m}\_i(t) + 2\xi_i\omega_i\dot{m}\_i(t) + \omega_i^2m_i(t) =
f_i(t).\$\$ The reserved formula terms `biomarker` and `velocity`
represent the latent log-coefficients \\\log\omega_i^2\\ and
\\\log(2\xi_i\omega_i)\\, respectively.

## Usage

``` r
MarginalODE(formula, data, time = "time", id = "id", control = list())
```

## Arguments

- formula:

  Longitudinal formula. The left-hand side is the observed response. On
  the right-hand side, `biomarker` and `velocity` activate the ODE
  dynamic parameters; other terms enter the forcing function. Random
  effects use lme-style syntax: `|` for full covariance and `||` for
  diagonal covariance.

- data:

  Data frame with longitudinal measurements

- time:

  Time variable name (default: `"time"`)

- id:

  Subject identifier name (default: `"id"`)

- control:

  List of control parameters. See
  [`MarginalODE.control`](https://gongziyang.com/JointODE/reference/MarginalODE.control.md).

## Value

S3 object of class `MarginalODE`

## Examples

``` r
if (FALSE) { # \dontrun{
data(sim)
fit <- MarginalODE(
  formula = observed ~ biomarker + velocity + x1 + x2 +
    (1 + biomarker + velocity || id),
  data = sim$data$longitudinal_data
)
} # }
```
