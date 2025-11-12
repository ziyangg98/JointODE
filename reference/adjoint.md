# Adjoint Sensitivity Analysis for ODE Systems

Computes gradients of scalar-valued objective functions with respect to
ODE parameters using the adjoint sensitivity method. This implementation
requires analytical derivatives for optimal performance and accuracy.

## Usage

``` r
adjoint(
  ode_func,
  jacobian_func,
  x0,
  params,
  times,
  data = NULL,
  objective_func = NULL,
  objective_grad = NULL,
  running_cost = NULL,
  running_cost_grad = NULL,
  rtol = 1e-08,
  atol = 1e-10,
  method = "lsoda",
  save_trajectory = FALSE
)
```

## Arguments

- ode_func:

  Function defining the ODE system: function(t, x, params, data)
  returning list with element 'dx' containing the state derivatives
  dx/dt

- jacobian_func:

  Required function providing analytical derivatives: function(t, x,
  params, data) returning list with:

  - df_dx: State Jacobian \\\partial f/\partial x\\ (n_states × n_states
    matrix)

  - df_dtheta: Parameter Jacobian \\\partial f/\partial \theta\\
    (n_states × n_params matrix)

- x0:

  Initial state vector at time t0

- params:

  Parameter vector \\\theta\\ to compute sensitivities for

- times:

  Time grid for ODE integration (must include t0 and T)

- data:

  Optional list containing constant auxiliary data

- objective_func:

  Terminal cost function g(x(T), data) Returns scalar objective value at
  final time

- objective_grad:

  Required when objective_func is provided: function(x_final, data)
  returning gradient \\\partial g/\partial x\\ at final time

- running_cost:

  Integrand function L(t, x, data) for running cost Returns scalar cost
  rate at time t

- running_cost_grad:

  Required when running_cost is provided: function(t, x, data) returning
  gradient \\\partial L/\partial x\\

- rtol:

  Relative error tolerance for ODE solver (default: 1e-8)

- atol:

  Absolute error tolerance for ODE solver (default: 1e-10)

- method:

  ODE solver algorithm (default: "lsoda" - adaptive solver)

- save_trajectory:

  Whether to return the full state trajectory

## Value

Object of class "adjoint" containing:

- objective:

  Scalar value of the objective function J

- gradient:

  Gradient vector \\dJ/d\theta\\ with respect to parameters

- final_state:

  System state at final time x(T)

- sensitivity_final:

  Sensitivity matrix \\\partial x(T)/\partial \theta\\ at final time

- trajectory:

  Full state trajectory (if save_trajectory = TRUE)

- n_states:

  Number of state variables

- n_params:

  Number of parameters

- times:

  Time grid used for integration

## Details

For a dynamical system described by ordinary differential equations:
\$\$dx/dt = f(t, x, \theta, data)\$\$ \$\$x(t_0) = x_0\$\$

And an objective functional: \$\$J = g(x(T), data) + \int\_{t_0}^{T}
L(t, x(t), data) dt\$\$

The adjoint method efficiently computes the gradient \\dJ/d\theta\\ by
solving:

1.  Forward ODE with sensitivity equations (forward pass)

2.  Adjoint ODE backward in time (if running cost present)

This is particularly efficient when the number of parameters exceeds the
number of objective functions.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example: Parameter estimation for exponential decay

# Define the ODE system: dx/dt = -\theta x
ode_system <- function(t, x, params, data) {
  list(dx = -params[1] * x)
}

# Provide analytical Jacobians
jacobians <- function(t, x, params, data) {
  list(
    df_dx = matrix(-params[1], 1, 1),      # \partial f/\partial x
    df_dtheta = matrix(-x, 1, 1)           # \partial f/\partial \theta
  )
}

# Define objective: squared error from target
target_value <- 0.5
objective <- function(x_final, data) {
  (x_final - target_value)^2
}

# Gradient of objective
objective_gradient <- function(x_final, data) {
  2 * (x_final - target_value)
}

# Compute sensitivity
result <- adjoint(
  ode_func = ode_system,
  jacobian_func = jacobians,
  x0 = 1,
  params = 0.5,
  times = seq(0, 2, length.out = 21),
  objective_func = objective,
  objective_grad = objective_gradient
)

print(result)
} # }
```
