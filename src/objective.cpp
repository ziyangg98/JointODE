#include <RcppArmadillo.h>

#include <cppad/cppad.hpp>

#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Compute objective for a single subject, building its own small tape.
// Returns {value, gradient (if requested), hessian (if requested)}.
static void compute_subject_objective(
    const std::vector<double>& params_vec,
    const List& subject_data,
    const NumericVector& b_i,
    const List& parameters,
    const ODEParams<double>& params_template,
    double inv_sigma_e2,
    double log_2pi_sigma_e2,
    const mat& inv_sigma_b,
    double const_re_term,
    double weight,
    bool need_gradient,
    bool need_hessian,
    double& obj_out,
    std::vector<double>& grad_out,
    std::vector<double>& hess_out) {
  const int n_baseline = params_template.baseline_coefs.size();
  const int n_hazard = params_template.hazard_coefs.size();
  const int n_longitudinal = params_template.longitudinal_coefs.size();
  const int n_params = n_baseline + n_hazard + n_longitudinal;
  const int n_random_effects = inv_sigma_b.n_rows;

  ADvector ad_params(n_params);
  std::copy(params_vec.begin(), params_vec.end(), ad_params.begin());
  CppAD::Independent(ad_params);

  ODEParams<ADdouble> ode_params;
  ode_params.baseline_coefs.resize(n_baseline);
  ode_params.hazard_coefs.resize(n_hazard);
  ode_params.longitudinal_coefs.resize(n_longitudinal);

  std::copy_n(ad_params.begin(), n_baseline, ode_params.baseline_coefs.begin());
  std::copy_n(ad_params.begin() + n_baseline, n_hazard,
              ode_params.hazard_coefs.begin());
  std::copy_n(ad_params.begin() + n_baseline + n_hazard, n_longitudinal,
              ode_params.longitudinal_coefs.begin());

  fill_ode_config(ode_params, parameters);

  const std::vector<ADdouble> subject_random_effects(b_i.begin(), b_i.end());
  fill_subject_data(ode_params.subject, subject_data, subject_random_effects);

  const ADdouble half_inv_sigma_e2(0.5 * inv_sigma_e2);

  BSplineWorkspace workspace;
  compute_bspline_basis(0.0, ode_params.spline_degree, ode_params.spline_knots,
                        ode_params.spline_boundary, workspace.basis,
                        workspace.knots, workspace.work1, workspace.work2,
                        false);

  const std::vector<ADdouble> y0 = {
      ADdouble(0.0), ADdouble(ode_params.subject.initial_state[0]),
      ADdouble(ode_params.subject.initial_state[1])};
  const std::vector<double> times = build_times(
      ode_params.subject.longitudinal_times, ode_params.subject.event_time);
  const auto solution = solve_ode(y0, times, ode_params);

  ADdouble subject_ll(0.0);

  const int event_idx = find_time_idx(times, ode_params.subject.event_time);
  if (event_idx < 0) {
    stop("Subject: Event time not found in ODE solution");
  }

  const auto& event_state = solution[event_idx];
  subject_ll -= event_state[0];

  if (ode_params.subject.status == 1) {
    compute_bspline_basis(ode_params.subject.event_time,
                          ode_params.spline_degree, ode_params.spline_knots,
                          ode_params.spline_boundary, workspace.basis,
                          workspace.knots, workspace.work1, workspace.work2,
                          true);
    subject_ll += compute_log_hazard(
        event_state[1], event_state[2], workspace.basis,
        ode_params.baseline_coefs, ode_params.hazard_coefs,
        ode_params.subject.survival_covariates, ode_params.gamma);
  }

  const int n_obs = ode_params.subject.longitudinal_times.size();
  if (n_obs > 0) {
    for (int t = 0; t < n_obs; t++) {
      const int idx =
          find_time_idx(times, ode_params.subject.longitudinal_times[t]);
      if (idx < 0) {
        stop("Subject: Observation time not found in ODE solution");
      }
      const ADdouble residual =
          ADdouble(ode_params.subject.longitudinal_measurements[t]) -
          solution[idx][1];
      subject_ll -= half_inv_sigma_e2 * residual * residual;
    }
    subject_ll -= ADdouble(0.5 * n_obs * log_2pi_sigma_e2);
  }

  double quadratic_form = 0.0;
  for (int j = 0; j < n_random_effects; j++) {
    double temp = 0.0;
    for (int k = 0; k < n_random_effects; k++) {
      temp += inv_sigma_b(j, k) * b_i[k];
    }
    quadratic_form += b_i[j] * temp;
  }
  subject_ll -= ADdouble(0.5 * quadratic_form + const_re_term);

  // Negate (we minimize negative log-likelihood)
  ADdouble neg_ll = -subject_ll;

  CppAD::ADFun<double> tape;
  tape.Dependent(ad_params, ADvector{neg_ll});
  tape.optimize();

  const std::vector<double> y = tape.Forward(0, params_vec);
  obj_out = weight * y[0];

  if (need_gradient || need_hessian) {
    const std::vector<double> grad =
        tape.Reverse(1, std::vector<double>(1, 1.0));
    if (need_gradient) {
      for (int j = 0; j < n_params; j++) {
        grad_out[j] += weight * grad[j];
      }
    }
    if (need_hessian) {
      const std::vector<double> hess_vec = tape.Hessian(params_vec, 0);
      for (int j = 0; j < n_params * n_params; j++) {
        hess_out[j] += weight * hess_vec[j];
      }
    }
  }

  tape.capacity_order(0);
}

// [[Rcpp::export(.compute_objective_cppad)]]
NumericVector compute_objective_cppad(
    const NumericVector& params, const List& data_list,
    const NumericMatrix& random_effects, const List& parameters,
    Nullable<NumericVector> weights = R_NilValue, bool gradient = true,
    bool hessian = false) {
  const int n_subjects = data_list.size();

  // Setup weights: default to 1.0 for all subjects
  std::vector<double> subject_weights(n_subjects, 1.0);
  if (weights.isNotNull()) {
    NumericVector w(weights);
    if (w.size() != n_subjects) {
      stop("weights vector size (%d) does not match number of subjects (%d)",
           w.size(), n_subjects);
    }
    subject_weights = as<std::vector<double>>(w);
  }

  ODEParams<double> params_template;
  fill_ode_parameters(params_template, parameters);

  const int n_baseline = params_template.baseline_coefs.size();
  const int n_hazard = params_template.hazard_coefs.size();
  const int n_longitudinal = params_template.longitudinal_coefs.size();
  const int n_params = n_baseline + n_hazard + n_longitudinal;

  const double sigma_e = params_template.measurement_error_sd;
  const mat sigma_b = params_template.random_effect_sigma;

  if (params.size() != n_params) {
    stop("Parameter vector size (%d) does not match expected size (%d)",
         params.size(), n_params);
  }

  const double inv_sigma_e2 = 1.0 / (sigma_e * sigma_e);
  const double log_2pi_sigma_e2 = std::log(2.0 * M_PI * sigma_e * sigma_e);

  mat inv_sigma_b;
  if (!inv(inv_sigma_b, sigma_b)) {
    stop("Failed to invert random effect covariance matrix");
  }
  const int n_random_effects = sigma_b.n_rows;
  double log_det_sb, sign;
  log_det(log_det_sb, sign, sigma_b);
  const double const_re_term =
      0.5 * (n_random_effects * std::log(2.0 * M_PI) + log_det_sb);

  const std::vector<double> params_vec(params.begin(), params.end());

  double total_obj = 0.0;
  std::vector<double> total_grad(n_params, 0.0);
  std::vector<double> total_hess(n_params * n_params, 0.0);

  for (int i = 0; i < n_subjects; i++) {
    double subj_obj = 0.0;
    compute_subject_objective(
        params_vec, data_list[i], random_effects(i, _), parameters,
        params_template, inv_sigma_e2, log_2pi_sigma_e2, inv_sigma_b,
        const_re_term, subject_weights[i], gradient, hessian, subj_obj,
        total_grad, total_hess);
    total_obj += subj_obj;
  }

  NumericVector result(1);
  result[0] = total_obj;

  if (gradient) {
    result.attr("gradient") = wrap(total_grad);
  }
  if (hessian) {
    NumericMatrix hess_matrix(n_params, n_params);
    for (int i = 0; i < n_params; i++) {
      for (int j = 0; j < n_params; j++) {
        hess_matrix(i, j) = total_hess[i * n_params + j];
      }
    }
    result.attr("hessian") = hess_matrix;
  }

  return result;
}
