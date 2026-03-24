#include <RcppArmadillo.h>
#include <cppad/cppad.hpp>
#include "solver.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ============================================================================
// Objective: M-step negative log-likelihood w.r.t. fixed effects theta
// ============================================================================

// Per-subject AD tape: builds tape for theta, evaluates and accumulates
static void eval_subject_nll(
    const std::vector<double>& theta_vec,
    const List& subject_data,
    const NumericVector& random_effect,
    const List& parameters,
    const ODEParams<double>& params_template,
    double inv_sigma_e2,
    double log_2pi_sigma_e2,
    const mat& inv_sigma_b,
    double re_const,
    double weight,
    bool need_grad,
    bool need_hess,
    double& obj_out,
    std::vector<double>& grad_out,
    std::vector<double>& hess_out) {
  const int n_baseline = params_template.baseline_coefs.size();
  const int n_hazard = params_template.hazard_coefs.size();
  const int n_long = params_template.longitudinal_coefs.size();
  const int n_params = n_baseline + n_hazard + n_long;
  const int n_re = inv_sigma_b.n_rows;

  // Update branch and sub-steps with random effects (before Independent)
  NumericVector long_coefs = as<List>(parameters["coefficients"])["longitudinal"];
  MatExpBranch branch = params_template.branch;
  update_branch(branch, long_coefs, random_effect,
                params_template.biomarker_random,
                params_template.velocity_random);

  // Build AD tape: theta as independent variable
  ADvector ad_theta(n_params);
  std::copy(theta_vec.begin(), theta_vec.end(), ad_theta.begin());
  CppAD::Independent(ad_theta);

  ODEParams<ADdouble> ode;
  ode.baseline_coefs.resize(n_baseline);
  ode.hazard_coefs.resize(n_hazard);
  ode.longitudinal_coefs.resize(n_long);
  std::copy_n(ad_theta.begin(), n_baseline, ode.baseline_coefs.begin());
  std::copy_n(ad_theta.begin() + n_baseline, n_hazard,
              ode.hazard_coefs.begin());
  std::copy_n(ad_theta.begin() + n_baseline + n_hazard, n_long,
              ode.longitudinal_coefs.begin());
  load_config(ode, parameters);
  ode.branch = branch;
  load_subject(ode.subject, subject_data,
      std::vector<ADdouble>(random_effect.begin(), random_effect.end()));

  // Solve ODE and compute log-likelihood
  const auto times = time_grid(ode.subject.longitudinal_times,
                                  ode.subject.event_time);
  const std::vector<ADdouble> y0 = {ADdouble(0.0),
      ADdouble(ode.subject.initial_state[0]),
      ADdouble(ode.subject.initial_state[1])};
  const auto sol = ode_solve(y0, times, ode);

  BSplineWorkspace ws;
  bspline_basis(0.0, ode.spline_degree, ode.spline_knots,
      ode.spline_boundary, ws.basis, ws.knots, ws.work1, ws.work2, false);
  ADdouble ll = joint_loglik(sol, times, ode, ws,
                                      inv_sigma_e2, log_2pi_sigma_e2);

  // Add random effect prior (constant w.r.t. theta)
  double qf = 0.0;
  for (int j = 0; j < n_re; j++) {
    double s = 0.0;
    for (int k = 0; k < n_re; k++) s += inv_sigma_b(j, k) * random_effect[k];
    qf += random_effect[j] * s;
  }
  ll -= ADdouble(0.5 * qf + re_const);

  // Evaluate tape
  CppAD::ADFun<double> tape;
  tape.Dependent(ad_theta, ADvector{-ll});
  tape.optimize();
  tape.check_for_nan(false);
  obj_out = weight * tape.Forward(0, theta_vec)[0];

  if (need_grad || need_hess) {
    auto grad = tape.Reverse(1, std::vector<double>(1, 1.0));
    if (need_grad)
      for (int j = 0; j < n_params; j++) grad_out[j] += weight * grad[j];
    if (need_hess) {
      auto hess = tape.Hessian(theta_vec, 0);
      for (int j = 0; j < n_params * n_params; j++)
        hess_out[j] += weight * hess[j];
    }
  }
  tape.capacity_order(0);
}

// [[Rcpp::export(.compute_objective_cppad)]]
NumericVector compute_objective_cppad(
    const NumericVector& params, const List& data_list,
    const NumericMatrix& random_effects, const List& parameters,
    Nullable<NumericVector> weights = R_NilValue,
    bool gradient = true, bool hessian = false) {
  const int n_subjects = data_list.size();

  std::vector<double> subject_weights(n_subjects, 1.0);
  if (weights.isNotNull())
    subject_weights = as<std::vector<double>>(NumericVector(weights));

  ODEParams<double> params_template;
  load_params(params_template, parameters);

  const int n_baseline = params_template.baseline_coefs.size();
  const int n_hazard = params_template.hazard_coefs.size();
  const int n_long = params_template.longitudinal_coefs.size();
  const int n_params = n_baseline + n_hazard + n_long;

  const double sigma_e = params_template.measurement_error_sd;
  const double inv_sigma_e2 = 1.0 / (sigma_e * sigma_e);
  const double log_2pi_sigma_e2 = std::log(2.0 * M_PI * sigma_e * sigma_e);

  mat inv_sigma_b;
  if (!inv(inv_sigma_b, params_template.random_effect_sigma))
    stop("Failed to invert random effect covariance matrix");
  double log_det_val, sign;
  log_det(log_det_val, sign, params_template.random_effect_sigma);
  const double re_const =
      0.5 * (inv_sigma_b.n_rows * std::log(2.0 * M_PI) + log_det_val);

  const std::vector<double> theta_vec(params.begin(), params.end());
  double total_obj = 0.0;
  std::vector<double> total_grad(n_params, 0.0);
  std::vector<double> total_hess(n_params * n_params, 0.0);

  for (int i = 0; i < n_subjects; i++) {
    double subj_obj = 0.0;
    eval_subject_nll(theta_vec, data_list[i], random_effects(i, _),
        parameters, params_template, inv_sigma_e2, log_2pi_sigma_e2,
        inv_sigma_b, re_const, subject_weights[i],
        gradient, hessian, subj_obj, total_grad, total_hess);
    total_obj += subj_obj;
  }

  NumericVector result(1);
  result[0] = total_obj;
  if (gradient) result.attr("gradient") = wrap(total_grad);
  if (hessian) {
    NumericMatrix hess_mat(n_params, n_params);
    for (int i = 0; i < n_params; i++)
      for (int j = 0; j < n_params; j++)
        hess_mat(i, j) = total_hess[i * n_params + j];
    result.attr("hessian") = hess_mat;
  }
  return result;
}

// ============================================================================
// Posterior: E-step log-posterior w.r.t. random effects b
// ============================================================================

// [[Rcpp::export(.compute_logpost_cppad)]]
NumericVector compute_logpost_cppad(
    const NumericVector& random_effect,
    const List& data, const List& parameters,
    bool gradient = true, bool hessian = false) {
  const int n_re = random_effect.size();

  // Build AD tape: random effects as independent variable
  ADvector ad_re(n_re);
  std::copy(random_effect.begin(), random_effect.end(), ad_re.begin());
  CppAD::Independent(ad_re);

  ODEParams<ADdouble> ode;
  load_subject(ode.subject, data,
      std::vector<ADdouble>(ad_re.begin(), ad_re.end()));
  load_params(ode, parameters);
  NumericVector long_coefs = as<List>(parameters["coefficients"])["longitudinal"];
  update_branch(ode.branch, long_coefs, random_effect,
                ode.biomarker_random, ode.velocity_random);

  const double sigma_e = ode.measurement_error_sd;
  const double inv_sigma_e2 = 1.0 / (sigma_e * sigma_e);
  const double log_2pi_sigma_e2 = std::log(2.0 * M_PI * sigma_e * sigma_e);

  mat inv_sigma_b;
  if (!inv(inv_sigma_b, ode.random_effect_sigma))
    stop("Failed to invert random effect covariance matrix");
  double log_det_val, sign;
  log_det(log_det_val, sign, ode.random_effect_sigma);
  const double re_const = 0.5 * (n_re * std::log(2.0 * M_PI) + log_det_val);

  // Solve ODE and compute log-likelihood
  const auto times = time_grid(ode.subject.longitudinal_times,
                                  ode.subject.event_time);
  const std::vector<ADdouble> y0 = {ADdouble(0.0),
      ADdouble(ode.subject.initial_state[0]),
      ADdouble(ode.subject.initial_state[1])};
  const auto sol = ode_solve(y0, times, ode);

  BSplineWorkspace ws;
  bspline_basis(0.0, ode.spline_degree, ode.spline_knots,
      ode.spline_boundary, ws.basis, ws.knots, ws.work1, ws.work2, false);
  ADdouble log_post = joint_loglik(sol, times, ode, ws,
                                            inv_sigma_e2, log_2pi_sigma_e2);

  // Add random effect prior: -0.5 * b' * inv(Sigma_b) * b - const
  std::vector<ADdouble> Sb(n_re, ADdouble(0.0));
  for (int i = 0; i < n_re; i++)
    for (int j = 0; j < n_re; j++)
      Sb[i] += ADdouble(inv_sigma_b(i, j)) * ad_re[j];
  ADdouble qf(0.0);
  for (int i = 0; i < n_re; i++) qf += ad_re[i] * Sb[i];
  log_post -= ADdouble(0.5) * qf + ADdouble(re_const);

  CppAD::ADFun<double> tape;
  tape.Dependent(ad_re, ADvector{log_post});
  return eval_tape(tape,
      std::vector<double>(random_effect.begin(), random_effect.end()),
      n_re, gradient, hessian);
}

// ============================================================================
// State: log-likelihood w.r.t. initial state [m(0), v(0)]
// ============================================================================

// [[Rcpp::export(.compute_state_loglik_cppad)]]
NumericVector compute_state_loglik_cppad(
    const NumericVector& initial_state,
    const List& data,
    const NumericVector& random_effect,
    const List& parameters,
    bool gradient = true, bool hessian = false) {

  // Build AD tape: initial state as independent variable
  ADvector ad_state(2);
  std::copy(initial_state.begin(), initial_state.end(), ad_state.begin());
  CppAD::Independent(ad_state);

  ODEParams<ADdouble> ode;
  load_subject(ode.subject, data,
      std::vector<ADdouble>(random_effect.begin(), random_effect.end()));
  load_params(ode, parameters);
  NumericVector long_coefs = as<List>(parameters["coefficients"])["longitudinal"];
  update_branch(ode.branch, long_coefs, random_effect,
                ode.biomarker_random, ode.velocity_random);

  const double sigma_e = ode.measurement_error_sd;
  const double inv_sigma_e2 = 1.0 / (sigma_e * sigma_e);
  const double log_2pi_sigma_e2 = std::log(2.0 * M_PI * sigma_e * sigma_e);

  // Solve ODE and compute log-likelihood
  const auto times = time_grid(ode.subject.longitudinal_times,
                                  ode.subject.event_time);
  const std::vector<ADdouble> y0 = {ADdouble(0.0), ad_state[0], ad_state[1]};
  const auto sol = ode_solve(y0, times, ode);

  BSplineWorkspace ws;
  bspline_basis(0.0, ode.spline_degree, ode.spline_knots,
      ode.spline_boundary, ws.basis, ws.knots, ws.work1, ws.work2, false);
  ADdouble ll = joint_loglik(sol, times, ode, ws,
                                      inv_sigma_e2, log_2pi_sigma_e2);

  CppAD::ADFun<double> tape;
  tape.Dependent(ad_state, ADvector{ll});
  return eval_tape(tape,
      std::vector<double>(initial_state.begin(), initial_state.end()),
      2, gradient, hessian);
}

// ============================================================================
// Batch ODE solver: trajectories for diagnostics (no AD)
// ============================================================================

// [[Rcpp::export(.solve_batch_ode_cppad)]]
List solve_batch_ode_cppad(const List& data_list,
                           const NumericMatrix& random_effects,
                           const List& parameters) {
  const int n_subjects = data_list.size();

  ODEParams<double> params;
  load_params(params, parameters);

  BSplineWorkspace ws;
  bspline_basis(0.0, params.spline_degree, params.spline_knots,
      params.spline_boundary, ws.basis, ws.knots,
      ws.work1, ws.work2, false);

  List results(n_subjects);

  for (int i = 0; i < n_subjects; i++) {
    NumericVector re = random_effects(i, _);
    load_subject(params.subject, data_list[i],
        std::vector<double>(re.begin(), re.end()));
    NumericVector long_coefs = as<List>(parameters["coefficients"])["longitudinal"];
    update_branch(params.branch, long_coefs, re,
                  params.biomarker_random, params.velocity_random);

    const std::vector<double> y0 = {0.0, params.subject.initial_state[0],
                                    params.subject.initial_state[1]};
    const auto times = time_grid(params.subject.longitudinal_times,
                                    params.subject.event_time);
    const auto sol = ode_solve(y0, times, params);

    const int n_times = times.size();
    NumericVector t_out(n_times), m_out(n_times), v_out(n_times);
    NumericVector a_out(n_times), H_out(n_times), lh_out(n_times);

    for (int t = 0; t < n_times; t++) {
      t_out[t] = times[t];
      H_out[t] = sol[t][0];
      m_out[t] = sol[t][1];
      v_out[t] = sol[t][2];

      a_out[t] = acceleration(sol[t][1], sol[t][2], times[t],
          params.subject, params.longitudinal_coefs,
          params.biomarker_fixed, params.biomarker_random,
          params.velocity_fixed, params.velocity_random);

      bspline_basis(times[t], params.spline_degree, params.spline_knots,
          params.spline_boundary, ws.basis, ws.knots,
          ws.work1, ws.work2, true);
      double lh = log_hazard(sol[t][1], sol[t][2], ws.basis,
          params.baseline_coefs, params.hazard_coefs,
          params.subject.survival_covariates, params.gamma);
      lh_out[t] = std::max(HAZARD_CLAMP_MIN, std::min(lh, HAZARD_CLAMP_MAX));
    }

    results[i] = List::create(
        Named("times") = t_out, Named("biomarker") = m_out,
        Named("velocity") = v_out, Named("acceleration") = a_out,
        Named("cum_hazard") = H_out, Named("log_hazard") = lh_out);
  }
  return results;
}
