# Package index

## Model Fitting

Core functionality for fitting joint ODE models combining survival and
longitudinal data

- [`JointODE()`](https://gongziyang.com/JointODE/reference/JointODE.md)
  : Joint Modeling of Longitudinal and Survival Data Using ODEs
- [`MarginalODE()`](https://gongziyang.com/JointODE/reference/MarginalODE.md)
  : Marginal Second-Order ODE Parameter Estimation

## Model Inspection

Extract and examine fitted model components and parameters

- [`coef(`*`<JointODE>`*`)`](https://gongziyang.com/JointODE/reference/coef.JointODE.md)
  : Extract Model Coefficients
- [`coef(`*`<MarginalODE>`*`)`](https://gongziyang.com/JointODE/reference/coef.MarginalODE.md)
  : Extract Model Coefficients
- [`logLik(`*`<JointODE>`*`)`](https://gongziyang.com/JointODE/reference/logLik.JointODE.md)
  : Extract Log-Likelihood
- [`logLik(`*`<MarginalODE>`*`)`](https://gongziyang.com/JointODE/reference/logLik.MarginalODE.md)
  : Extract Log-Likelihood
- [`vcov(`*`<JointODE>`*`)`](https://gongziyang.com/JointODE/reference/vcov.JointODE.md)
  : Extract Variance-Covariance Matrix
- [`vcov(`*`<MarginalODE>`*`)`](https://gongziyang.com/JointODE/reference/vcov.MarginalODE.md)
  : Extract Variance-Covariance Matrix
- [`predict(`*`<JointODE>`*`)`](https://gongziyang.com/JointODE/reference/predict.JointODE.md)
  : Predict Method for JointODE Objects
- [`predict(`*`<MarginalODE>`*`)`](https://gongziyang.com/JointODE/reference/predict.MarginalODE.md)
  : Predict Method for MarginalODE Objects

## Model Summary & Display

Summarize and display model results with statistical inference

- [`print(`*`<summary.MarginalODE>`*`)`](https://gongziyang.com/JointODE/reference/print.summary.MarginalODE.md)
  : Print Summary of MarginalODE Model
- [`summary(`*`<JointODE>`*`)`](https://gongziyang.com/JointODE/reference/summary.JointODE.md)
  : Summary of JointODE Model
- [`summary(`*`<MarginalODE>`*`)`](https://gongziyang.com/JointODE/reference/summary.MarginalODE.md)
  : Summary of MarginalODE Model
- [`plot(`*`<JointODE>`*`)`](https://gongziyang.com/JointODE/reference/plot.JointODE.md)
  : Plot Method for JointODE Objects
- [`print(`*`<JointODE>`*`)`](https://gongziyang.com/JointODE/reference/print.JointODE.md)
  : Print JointODE Model
- [`print(`*`<MarginalODE>`*`)`](https://gongziyang.com/JointODE/reference/print.MarginalODE.md)
  : Print MarginalODE Model

## Advanced Methods

Computational methods for sensitivity analysis and optimization

## Data & Simulation

Tools for simulating joint ODE data and example datasets

- [`sim`](https://gongziyang.com/JointODE/reference/sim.md) : Example
  Dataset for Joint ODE Model
- [`simulate()`](https://gongziyang.com/JointODE/reference/simulate.md)
  : Simulate Data from a Joint Ordinary Differential Equation Model

## Re-exports

- [`Surv`](https://gongziyang.com/JointODE/reference/Surv.md) :
  Re-exported from survival
