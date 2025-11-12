#include <RcppArmadillo.h>
#include <cppad/cppad.hpp>
#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export(.compute_state_loglik_cppad)]]
NumericVector compute_state_loglik_cppad(
    const NumericVector& initial_state,
    const List& data,
    const NumericVector& random_effect,
    const List& parameters,
    bool gradient = true,
    bool hessian = false
) {
  if (initial_state.size() != 2) {
    stop("initial_state must have length 2");
  }

  // Setup AD tape for initial_state
  ADvector ad_initial_state(2);
  std::copy(initial_state.begin(), initial_state.end(), ad_initial_state.begin());
  CppAD::Independent(ad_initial_state);

  // Setup ODE parameters (fixed)
  ODEParams<ADdouble> params;
  std::vector<ADdouble> subject_random_effects(random_effect.begin(), random_effect.end());

  fill_subject_data(params.subject, data, subject_random_effects);
  fill_ode_parameters(params, parameters);

  const double sigma_e = params.measurement_error_sd;
  const double inv_sigma_e2 = 1.0 / (sigma_e * sigma_e);
  const double log_2pi_sigma_e2 = std::log(2.0 * M_PI * sigma_e * sigma_e);

  // Initialize B-spline workspace
  BSplineWorkspace workspace;
  compute_bspline_basis(0.0, params.spline_degree, params.spline_knots,
                        params.spline_boundary, workspace.basis,
                        workspace.knots, workspace.work1, workspace.work2, false);

  // Build times and solve ODE with AD initial state
  const std::vector<double> times = build_times(
      params.subject.longitudinal_times, params.subject.event_time);
  const std::vector<ADdouble> y0 = {ADdouble(0.0), ad_initial_state[0], ad_initial_state[1]};
  const auto solution = solve_ode(y0, times, params);

  // Compute log-likelihood (longitudinal only)
  ADdouble log_lik(0.0);
  const ADdouble half_inv_sigma_e2(0.5 * inv_sigma_e2);
  const int n_obs = params.subject.longitudinal_times.size();

  for (int i = 0; i < n_obs; ++i) {
    const int idx = find_time_idx(times, params.subject.longitudinal_times[i]);
    if (idx >= 0) {
      const ADdouble residual =
          ADdouble(params.subject.longitudinal_measurements[i]) - solution[idx][1];
      log_lik -= half_inv_sigma_e2 * residual * residual;
    }
  }

  if (n_obs > 0) {
    log_lik -= ADdouble(0.5 * n_obs * log_2pi_sigma_e2);
  }

  // Create tape
  CppAD::ADFun<double> tape;
  tape.Dependent(ad_initial_state, ADvector{log_lik});
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
