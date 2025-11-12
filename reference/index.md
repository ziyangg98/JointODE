# Package index

## Model Fitting

Core functionality for fitting joint ODE models combining survival and
longitudinal data

- [`JointODE()`](http://gongziyang.com/JointODE/reference/JointODE.md) :
  Joint Modeling of Longitudinal and Survival Data Using ODEs
- [`MarginalODE()`](http://gongziyang.com/JointODE/reference/MarginalODE.md)
  : Marginal Second-Order ODE Parameter Estimation

## Model Inspection

Extract and examine fitted model components and parameters

- [`coef(`*`<JointODE>`*`)`](http://gongziyang.com/JointODE/reference/coef.JointODE.md)
  : Extract Model Coefficients
- [`logLik(`*`<JointODE>`*`)`](http://gongziyang.com/JointODE/reference/logLik.JointODE.md)
  : Extract Log-Likelihood
- [`vcov(`*`<JointODE>`*`)`](http://gongziyang.com/JointODE/reference/vcov.JointODE.md)
  : Extract Variance-Covariance Matrix
- [`predict(`*`<JointODE>`*`)`](http://gongziyang.com/JointODE/reference/predict.JointODE.md)
  : Predict Method for JointODE Objects
- [`predict(`*`<MarginalODE>`*`)`](http://gongziyang.com/JointODE/reference/predict.MarginalODE.md)
  : Predict Method for MarginalODE Objects

## Model Summary & Display

Summarize and display model results with statistical inference

- [`print(`*`<summary.MarginalODE>`*`)`](http://gongziyang.com/JointODE/reference/print.summary.MarginalODE.md)
  : Print Summary of MarginalODE Fit
- [`summary(`*`<JointODE>`*`)`](http://gongziyang.com/JointODE/reference/summary.JointODE.md)
  : Summary of JointODE Model
- [`summary(`*`<MarginalODE>`*`)`](http://gongziyang.com/JointODE/reference/summary.MarginalODE.md)
  : Summary Method for MarginalODE Objects
- [`plot(`*`<JointODE>`*`)`](http://gongziyang.com/JointODE/reference/plot.JointODE.md)
  : Plot Method for JointODE Objects
- [`print(`*`<JointODE>`*`)`](http://gongziyang.com/JointODE/reference/print.JointODE.md)
  : Print JointODE Model

## Advanced Methods

Computational methods for sensitivity analysis and optimization

- [`adjoint()`](http://gongziyang.com/JointODE/reference/adjoint.md) :
  Adjoint Sensitivity Analysis for ODE Systems
- [`print(`*`<adjoint>`*`)`](http://gongziyang.com/JointODE/reference/print.adjoint.md)
  : Print Method for Adjoint Results

## Data & Simulation

Tools for simulating joint ODE data and example datasets

- [`sim`](http://gongziyang.com/JointODE/reference/sim.md) : Example
  Dataset for Joint ODE Model
- [`simulate()`](http://gongziyang.com/JointODE/reference/simulate.md) :
  Simulate Data from a Joint Ordinary Differential Equation Model
