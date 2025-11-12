# Marginal Second-Order ODE Parameter Estimation

Estimates population-level (marginal) ODE parameters for longitudinal
biomarker trajectories using the second-order differential equation with
covariate support: \$\$\ddot{m}(t) = \text{value} \cdot m(t) +
\text{slope} \cdot \dot{m}(t) + X\beta\$\$

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

  A formula object specifying the response variable and covariates
  (e.g., `biomarker ~ x1 + x2` or `biomarker ~ 1` for intercept-only)

- data:

  A data frame containing the longitudinal measurements

- time:

  Character string specifying the time variable name (default: `"time"`)

- id:

  Character string specifying the subject identifier variable name
  (default: `"id"`)

- state:

  Optional matrix of initial conditions with two columns:

  - Column 1: Initial biomarker values \\m(0)\\

  - Column 2: Initial velocities \\\dot{m}(0)\\

  Each row corresponds to one subject. If `NULL`, initial values are
  estimated from the data.

- control:

  A list of control parameters for optimization, or output from
  [`JointODE.control`](http://gongziyang.com/JointODE/reference/JointODE.control.md).
  Key parameters include:

  `verbose`

  :   Verbosity level (default: FALSE)

  `parallel`

  :   Logical flag enabling parallel computation (default: FALSE)

  `n_cores`

  :   Number of CPU cores for parallel processing. If 0, automatically
      detects available cores (default: 0)

  See
  [`JointODE.control`](http://gongziyang.com/JointODE/reference/JointODE.control.md)
  for complete details.

## Value

A list (S3 class `MarginalODE`) containing:

- `parameters`:

  Named numeric vector of parameter estimates (value, slope, and
  covariate coefficients)

- `measurement_error_sd`:

  Residual standard deviation

- `logLik`:

  Log-likelihood value

- `AIC`:

  Akaike Information Criterion

- `BIC`:

  Bayesian Information Criterion

- `convergence`:

  List with `converged` (logical), `iterations`, and `message`

- `vcov`:

  Variance-covariance matrix from inverse Hessian (may be NA matrix if
  Hessian is singular)

- `data`:

  Processed data list used for fitting

- `control`:

  Control parameters used

- `call`:

  Matched function call

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate simulated data
fit <- MarginalODE(
  formula = observed ~ x1 + x2,
  data = sim$data$longitudinal_data,
  state = as.matrix(sim$data$state)
)
} # }
```
