#include <RcppArmadillo.h>
#include <cppad/cppad.hpp>
#include <algorithm>
#include "solver.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ============================================================================
// Marginal ODE Data Structures (独立于 utils.hpp)
// ============================================================================

template <typename Scalar>
struct MarginalODEParams {
  std::vector<double> times;
  std::vector<std::vector<double>> covariates;
  std::vector<Scalar> theta;  // [value, slope, cov1, cov2, ...]
  int n_cov;
  MatExpBranch branch;

  MarginalODEParams() : n_cov(0), branch(MatExpBranch::COMPLEX) {}
};

// ============================================================================
// ODE Solver using Matrix Exponential (exact analytical solution)
// ============================================================================

// Extract forcing term b at time index
template <typename Scalar>
inline Scalar marginal_forcing(const MarginalODEParams<Scalar>& params,
                               size_t time_idx) {
  Scalar b = Scalar(0.0);
  for (int c = 0; c < params.n_cov; c++) {
    b += Scalar(params.covariates[time_idx][c]) * params.theta[c + 2];
  }
  return b;
}

// Find time index for covariate lookup
inline size_t find_cov_idx(double t, const std::vector<double>& times) {
  if (t >= times.back()) return times.size() - 1;
  if (t <= times[0]) return 0;
  auto it = std::upper_bound(times.begin(), times.end(), t);
  return std::distance(times.begin(), it) - 1;
}

template <typename Scalar>
inline std::vector<std::vector<Scalar>> solve_marginal_ode(
    const std::vector<Scalar>& y0,
    const std::vector<double>& times,
    const MarginalODEParams<Scalar>& params) {

  const int n_times = times.size();
  std::vector<std::vector<Scalar>> solution(n_times, std::vector<Scalar>(2));

  Scalar m = y0[0];
  Scalar v = y0[1];
  solution[0] = {m, v};

  const Scalar beta1 = params.theta[0];
  const Scalar beta2 = params.theta[1];
  const MatExpBranch branch = params.branch;

  for (int idx = 1; idx < n_times; idx++) {
    const double dt = times[idx] - times[idx - 1];
    const size_t cov_idx = find_cov_idx(times[idx - 1], params.times);
    const Scalar b = marginal_forcing(params, cov_idx);

    ode_step(m, v, beta1, beta2, b, Scalar(dt), branch);
    solution[idx] = {m, v};
  }

  return solution;
}

// ============================================================================
// Helper Function (参考 objective.cpp 的 load_subject)
// ============================================================================

template <typename Scalar>
inline void fill_marginal_subject_data(
    MarginalODEParams<Scalar>& params,
    const List& subject_data) {

  NumericVector time_vec = subject_data["time"];
  NumericMatrix cov_mat = subject_data["covariates"];

  const int n_times = time_vec.size();
  const int n_cov = cov_mat.ncol();

  params.times.resize(n_times);
  params.covariates.resize(n_times);
  params.n_cov = n_cov;

  // Copy times
  for (int t = 0; t < n_times; t++) {
    params.times[t] = time_vec[t];
  }

  // Copy covariates
  for (int t = 0; t < n_times; t++) {
    params.covariates[t].resize(n_cov);
    for (int c = 0; c < n_cov; c++) {
      params.covariates[t][c] = cov_mat(t, c);
    }
  }
}

// ============================================================================
// Main Objective Function (参考 objective.cpp 结构)
// ============================================================================

// [[Rcpp::export(.compute_marginal_objective_cppad)]]
NumericVector compute_marginal_objective_cppad(
    const NumericVector& params,
    const List& data_list,
    bool gradient = true,
    bool hessian = false) {

  const int n_subjects = data_list.size();
  const int n_params = params.size();

  // Classify branch BEFORE Independent() using double values
  const MatExpBranch branch = classify_disc(params[0], params[1]);

  // Start recording for automatic differentiation
  ADvector ad_params(n_params);
  std::copy(params.begin(), params.end(), ad_params.begin());
  CppAD::Independent(ad_params);

  // Prepare theta vector
  std::vector<ADdouble> theta(n_params);
  for (int i = 0; i < n_params; i++) {
    theta[i] = ad_params[i];
  }

  ADdouble total_sse(0.0);

  // Process each subject
  for (int i = 0; i < n_subjects; i++) {
    List subject_data = data_list[i];

    NumericVector response = subject_data["response"];
    NumericVector initial = subject_data["initial"];
    const int n_obs = response.size();

    if (n_obs == 0) continue;

    // Setup parameters for this subject
    MarginalODEParams<ADdouble> ode_params;
    ode_params.theta = theta;
    ode_params.branch = branch;
    fill_marginal_subject_data(ode_params, subject_data);

    // Initial conditions [biomarker, velocity]
    std::vector<ADdouble> y0(2);
    y0[0] = ADdouble(initial[0]);
    y0[1] = ADdouble(initial[1]);

    // Solve ODE (同 objective.cpp: const auto solution = solve_ode(...))
    const std::vector<std::vector<ADdouble>> solution =
        solve_marginal_ode(y0, ode_params.times, ode_params);

    // Compute residuals
    for (int t = 0; t < n_obs; t++) {
      ADdouble residual = ADdouble(response[t]) - solution[t][0];
      total_sse += residual * residual;
    }
  }

  // Create tape (同 objective.cpp)
  CppAD::ADFun<double> tape;
  tape.Dependent(ad_params, ADvector{total_sse});

  // Optimize the tape to reduce memory and improve performance
  tape.optimize();

  // Evaluate at the parameter values (同 objective.cpp)
  tape.check_for_nan(false);
  const std::vector<double> x(params.begin(), params.end());
  const std::vector<double> y = tape.Forward(0, x);

  NumericVector result(1);
  result[0] = y[0];

  // Compute gradient if requested (同 objective.cpp)
  if (gradient) {
    const std::vector<double> grad = tape.Reverse(1, std::vector<double>(1, 1.0));
    result.attr("gradient") = wrap(grad);
  }

  // Compute Hessian if requested (同 objective.cpp)
  if (hessian) {
    const std::vector<double> hess = tape.Hessian(x, 0);
    arma::mat hess_mat(n_params, n_params);
    for (int i = 0; i < n_params; i++) {
      for (int j = 0; j < n_params; j++) {
        hess_mat(i, j) = hess[i * n_params + j];
      }
    }
    result.attr("hessian") = wrap(hess_mat);
  }

  // Free memory: clear Taylor coefficients
  // This releases memory used by Forward/Reverse/Hessian computations
  tape.capacity_order(0);

  return result;
}

// ============================================================================
// ODE Solver for Prediction (不需要梯度)
// ============================================================================

// [[Rcpp::export(.solve_marginal_ode_cppad)]]
List solve_marginal_ode_cppad(
    const NumericVector& theta,
    const NumericVector& initial,
    const NumericVector& times,
    const NumericMatrix& covariates) {

  const int n_times = times.size();
  const int n_cov = covariates.ncol();

  if (n_times == 0) {
    stop("times cannot be empty");
  }

  // Setup parameters
  MarginalODEParams<double> ode_params;
  ode_params.theta.resize(theta.size());
  for (int i = 0; i < theta.size(); i++) {
    ode_params.theta[i] = theta[i];
  }
  ode_params.n_cov = n_cov;
  ode_params.branch = classify_disc(theta[0], theta[1]);

  // Copy times
  ode_params.times.resize(n_times);
  for (int t = 0; t < n_times; t++) {
    ode_params.times[t] = times[t];
  }

  // Copy covariates
  ode_params.covariates.resize(n_times);
  for (int t = 0; t < n_times; t++) {
    ode_params.covariates[t].resize(n_cov);
    for (int c = 0; c < n_cov; c++) {
      ode_params.covariates[t][c] = covariates(t, c);
    }
  }

  // Initial conditions
  std::vector<double> y0(2);
  y0[0] = initial[0];
  y0[1] = initial[1];

  // Solve ODE using our local solve_marginal_ode function
  const std::vector<std::vector<double>> solution =
      solve_marginal_ode(y0, ode_params.times, ode_params);

  // Extract results
  NumericVector biomarker(n_times);
  NumericVector velocity(n_times);
  NumericVector acceleration(n_times);

  for (int t = 0; t < n_times; t++) {
    biomarker[t] = solution[t][0];
    velocity[t] = solution[t][1];

    // Compute acceleration: d²m/dt² = theta[0]*m + theta[1]*dm/dt + X*beta
    acceleration[t] = ode_params.theta[0] * solution[t][0] +
                      ode_params.theta[1] * solution[t][1];
    if (n_cov > 0) {
      for (int c = 0; c < n_cov; c++) {
        acceleration[t] += ode_params.covariates[t][c] * ode_params.theta[c + 2];
      }
    }
  }

  return List::create(
    Named("biomarker") = biomarker,
    Named("velocity") = velocity,
    Named("acceleration") = acceleration
  );
}

// ============================================================================
// Initial State Optimization (模仿 initial.cpp)
// ============================================================================

// [[Rcpp::export(.compute_marginal_state_loglik)]]
NumericVector compute_marginal_state_loglik(
    const NumericVector& initial_state,
    const List& subject_data,
    const NumericVector& theta,
    bool gradient = true,
    bool hessian = false) {

  if (initial_state.size() != 2) {
    stop("initial_state must have length 2");
  }

  // Setup AD tape for initial_state
  ADvector ad_initial_state(2);
  std::copy(initial_state.begin(), initial_state.end(), ad_initial_state.begin());
  CppAD::Independent(ad_initial_state);

  // Setup ODE parameters (fixed theta)
  MarginalODEParams<ADdouble> ode_params;
  ode_params.theta.resize(theta.size());
  for (int i = 0; i < theta.size(); i++) {
    ode_params.theta[i] = ADdouble(theta[i]);
  }
  ode_params.branch = classify_disc(theta[0], theta[1]);
  fill_marginal_subject_data(ode_params, subject_data);

  NumericVector response = subject_data["response"];
  const int n_obs = response.size();

  if (n_obs == 0) {
    stop("subject must have at least one observation");
  }

  // Initial conditions from AD variable
  std::vector<ADdouble> y0(2);
  y0[0] = ad_initial_state[0];
  y0[1] = ad_initial_state[1];

  // Solve ODE with AD initial state
  const std::vector<std::vector<ADdouble>> solution =
      solve_marginal_ode(y0, ode_params.times, ode_params);

  // Compute SSE (sum of squared errors)
  ADdouble sse(0.0);
  for (int t = 0; t < n_obs; t++) {
    ADdouble residual = ADdouble(response[t]) - solution[t][0];
    sse += residual * residual;
  }

  // Create tape
  CppAD::ADFun<double> tape;
  tape.Dependent(ad_initial_state, ADvector{sse});
  tape.optimize();

  // Evaluate
  tape.check_for_nan(false);
  std::vector<double> x(initial_state.begin(), initial_state.end());
  std::vector<double> y = tape.Forward(0, x);

  NumericVector result(1);
  result[0] = y[0];

  if (gradient || hessian) {
    std::vector<double> grad = tape.Reverse(1, std::vector<double>(1, 1.0));
    if (gradient) {
      result.attr("gradient") = wrap(grad);
    }
    if (hessian) {
      result.attr("hessian") = wrap(tape.Hessian(x, 0));
    }
  }

  tape.capacity_order(0);
  return result;
}
