#include <RcppArmadillo.h>

#include <cppad/cppad.hpp>

#include "utils.hpp"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export(.compute_logpost_cppad)]]
NumericVector compute_logpost_cppad(const NumericVector& random_effect,
                                    const List& data, const List& parameters,
                                    bool gradient = true,
                                    bool hessian = false) {
  const int n_random_effects = random_effect.size();

  ADvector ad_random_effects(n_random_effects);
  std::copy(random_effect.begin(), random_effect.end(),
            ad_random_effects.begin());
  CppAD::Independent(ad_random_effects);

  ODEParams<ADdouble> params;
  fill_subject_data(params.subject, data,
                    std::vector<ADdouble>(ad_random_effects.begin(),
                                          ad_random_effects.end()));
  fill_ode_parameters(params, parameters);

  const double sigma_e = params.measurement_error_sd;
  const double inv_sigma_e2 = 1.0 / (sigma_e * sigma_e);
  const double log_2pi_sigma_e2 = std::log(2.0 * M_PI * sigma_e * sigma_e);

  mat inv_sigma_b;
  if (!inv(inv_sigma_b, params.random_effect_sigma)) {
    stop("Failed to invert random effect covariance matrix");
  }
  double log_det_sb, sign;
  log_det(log_det_sb, sign, params.random_effect_sigma);
  const double const_re_term =
      0.5 * (n_random_effects * std::log(2.0 * M_PI) + log_det_sb);

  BSplineWorkspace workspace;
  compute_bspline_basis(0.0, params.spline_degree, params.spline_knots,
                        params.spline_boundary, workspace.basis,
                        workspace.knots, workspace.work1, workspace.work2,
                        false);

  const std::vector<double> times =
      build_times(params.subject.longitudinal_times, params.subject.event_time);
  const std::vector<ADdouble> y0 = {ADdouble(0.0),
                                    ADdouble(params.subject.initial_state[0]),
                                    ADdouble(params.subject.initial_state[1])};
  const auto solution = solve_ode(y0, times, params);

  ADdouble log_posterior(0.0);

  const int n_obs = params.subject.longitudinal_times.size();
  const ADdouble half_inv_sigma_e2(0.5 * inv_sigma_e2);
  for (int i = 0; i < n_obs; i++) {
    const int idx = find_time_idx(times, params.subject.longitudinal_times[i]);
    if (idx >= 0) {
      const ADdouble residual =
          ADdouble(params.subject.longitudinal_measurements[i]) -
          solution[idx][1];
      log_posterior -= half_inv_sigma_e2 * residual * residual;
    }
  }
  if (n_obs > 0) {
    log_posterior -= ADdouble(0.5 * n_obs * log_2pi_sigma_e2);
  }

  // Find the state at event time (not the last time point!)
  const int event_idx = find_time_idx(times, params.subject.event_time);
  if (event_idx < 0) {
    stop("Event time not found in ODE solution");
  }
  const auto& event_state = solution[event_idx];
  log_posterior -= event_state[0];

  if (params.subject.status == 1) {
    compute_bspline_basis(params.subject.event_time, params.spline_degree,
                          params.spline_knots, params.spline_boundary,
                          workspace.basis, workspace.knots, workspace.work1,
                          workspace.work2, true);
    log_posterior += compute_log_hazard(
        event_state[1], event_state[2], workspace.basis, params.baseline_coefs,
        params.hazard_coefs, params.subject.survival_covariates, params.gamma);
  }

  std::vector<ADdouble> inv_sigma_b_times_b(n_random_effects, ADdouble(0.0));
  for (int i = 0; i < n_random_effects; i++) {
    for (int j = 0; j < n_random_effects; j++) {
      inv_sigma_b_times_b[i] +=
          ADdouble(inv_sigma_b(i, j)) * ad_random_effects[j];
    }
  }
  ADdouble quadratic_form(0.0);
  for (int i = 0; i < n_random_effects; i++) {
    quadratic_form += ad_random_effects[i] * inv_sigma_b_times_b[i];
  }
  log_posterior -= ADdouble(0.5) * quadratic_form;
  log_posterior -= ADdouble(const_re_term);

  CppAD::ADFun<double> tape;
  tape.Dependent(ad_random_effects, ADvector{log_posterior});

  // Optimize the tape to reduce memory and improve performance
  tape.optimize();

  const std::vector<double> x(random_effect.begin(), random_effect.end());
  const std::vector<double> y = tape.Forward(0, x);

  NumericVector result(1);
  result[0] = y[0];

  if (gradient || hessian) {
    const std::vector<double> grad =
        tape.Reverse(1, std::vector<double>(1, 1.0));
    if (gradient) {
      result.attr("gradient") = wrap(grad);
    }
    if (hessian) {
      const std::vector<double> hess_vec = tape.Hessian(x, 0);
      NumericMatrix hess_matrix(n_random_effects, n_random_effects);
      for (int i = 0; i < n_random_effects; i++) {
        for (int j = 0; j < n_random_effects; j++) {
          hess_matrix(i, j) = hess_vec[i * n_random_effects + j];
        }
      }
      result.attr("hessian") = hess_matrix;
    }
  }

  // Free memory: clear Taylor coefficients
  // This releases memory used by Forward/Reverse/Hessian computations
  tape.capacity_order(0);

  return result;
}
