#include <RcppArmadillo.h>
#include <cppad/cppad.hpp>
#include <algorithm>
#include "utils.h"

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
  std::vector<Scalar> theta;  // [value, slope, beta1, beta2, ...]
  int n_cov;

  MarginalODEParams() : n_cov(0) {}
};

// ============================================================================
// Marginal ODE Functor
// ============================================================================

template <typename Scalar>
class MarginalODEFunctor {
 private:
  const MarginalODEParams<Scalar>& params_;

 public:
  explicit MarginalODEFunctor(const MarginalODEParams<Scalar>& params)
      : params_(params) {}

  void Ode(const Scalar& t, const CppAD::vector<Scalar>& y,
           CppAD::vector<Scalar>& dydt) const {
    const Scalar& m = y[0];      // biomarker
    const Scalar& dmdt = y[1];   // velocity

    // Find time index for covariate lookup (same logic as utils.hpp)
    double t_val = get_value(t);
    size_t time_idx;
    if (t_val >= params_.times.back()) {
      time_idx = params_.times.size() - 1;
    } else if (t_val <= params_.times[0]) {
      time_idx = 0;
    } else {
      // Binary search for the interval containing t
      auto it = std::upper_bound(params_.times.begin(), params_.times.end(), t_val);
      time_idx = std::distance(params_.times.begin(), it) - 1;
    }

    // Compute acceleration: d²m/dt² = theta[0]*m + theta[1]*dm/dt + X*beta
    Scalar d2mdt2 = params_.theta[0] * m + params_.theta[1] * dmdt;

    if (params_.n_cov > 0) {
      for (int c = 0; c < params_.n_cov; c++) {
        d2mdt2 += Scalar(params_.covariates[time_idx][c]) * params_.theta[c + 2];
      }
    }

    dydt[0] = dmdt;
    dydt[1] = d2mdt2;
  }
};

// ============================================================================
// ODE Solver (参考 utils.hpp 的 solve_ode 模式)
// ============================================================================

template <typename Scalar>
inline std::vector<std::vector<Scalar>> solve_marginal_ode(
    const std::vector<Scalar>& y0,
    const std::vector<double>& times,
    const MarginalODEParams<Scalar>& params) {

  const int n_times = times.size();

  MarginalODEFunctor<Scalar> ode_fun(params);

  std::vector<std::vector<Scalar>> solution(n_times, std::vector<Scalar>(2));
  CppAD::vector<Scalar> yi(2);

  // Initialize with initial conditions
  std::copy(y0.begin(), y0.end(), yi.begin());
  std::copy(y0.begin(), y0.end(), solution[0].begin());

  Scalar ti = Scalar(times[0]);
  Scalar tf, dt;

  // Solve step by step
  for (int idx = 1; idx < n_times; idx++) {
    tf = Scalar(times[idx]);
    dt = tf - ti;

    // Compute ODE steps for this interval (same as utils.hpp)
    const size_t M = std::max(
        MIN_ODE_STEPS,
        static_cast<size_t>(std::ceil(get_value(dt) * STEPS_PER_UNIT_TIME)));

    yi = CppAD::Runge45(ode_fun, M, ti, tf, yi);

    for (int i = 0; i < 2; ++i) {
      solution[idx][i] = yi[i];
    }

    ti = tf;
  }

  return solution;
}

// ============================================================================
// Helper Function (参考 objective.cpp 的 fill_subject_data)
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

  // Start recording for automatic differentiation (同 objective.cpp)
  ADvector ad_params(n_params);
  std::copy(params.begin(), params.end(), ad_params.begin());
  CppAD::Independent(ad_params);

  // Prepare theta vector
  std::vector<ADdouble> theta(n_params);
  for (int i = 0; i < n_params; i++) {
    theta[i] = ad_params[i];
  }

  ADdouble total_sse(0.0);

  // Process each subject (同 objective.cpp 循环模式)
  for (int i = 0; i < n_subjects; i++) {
    List subject_data = data_list[i];

    NumericVector response = subject_data["response"];
    NumericVector initial = subject_data["initial"];
    const int n_obs = response.size();

    if (n_obs == 0) continue;

    // Setup parameters for this subject
    MarginalODEParams<ADdouble> ode_params;
    ode_params.theta = theta;
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

  return result;
}
