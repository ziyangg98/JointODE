#include <RcppArmadillo.h>
#include <cppad/cppad.hpp>
#include "solver.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ============================================================================
// Joint: M-step objective w.r.t. fixed effects theta
// ============================================================================

// Load coefficients from flat AD vector into ODEParams
template <typename DstScalar, typename SrcVec>
static void fill_coefs(ODEParams<DstScalar>& p, const SrcVec& theta,
                        int n_bl, int n_hz, int n_lg, int n_in) {
  int off = 0;
  p.baseline_coefs.resize(n_bl);
  for (int j = 0; j < n_bl; j++) p.baseline_coefs[j] = DstScalar(theta[off++]);
  p.hazard_coefs.resize(n_hz);
  for (int j = 0; j < n_hz; j++) p.hazard_coefs[j] = DstScalar(theta[off++]);
  p.longitudinal_coefs.resize(n_lg);
  for (int j = 0; j < n_lg; j++) p.longitudinal_coefs[j] = DstScalar(theta[off++]);
  p.initial_state_coefs.resize(n_in);
  for (int j = 0; j < n_in; j++) p.initial_state_coefs[j] = DstScalar(theta[off++]);
}

static const char* SPD_ERR_SIGMA =
    "random_effect_sigma must be positive definite";

static void assert_spd(const mat& M, const char* msg) {
  mat R;
  if (!chol(R, M)) stop(msg);
}

static void eval_subject_nll(
    const std::vector<double>& theta_vec,
    const List& subject_data,
    const NumericVector& random_effect,
    const List& parameters,
    const ODEParams<double>& params_template,
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
  const int n_init = params_template.initial_state_coefs.size();
  const int n_params = n_baseline + n_hazard + n_long + n_init;
  const int n_re = inv_sigma_b.n_rows;

  NumericVector coef_re(random_effect.begin() + 2, random_effect.end());
  NumericVector proposed_long(n_long);
  for (int j = 0; j < n_long; j++)
    proposed_long[j] = theta_vec[n_baseline + n_hazard + j];
  MatExpBranch branch = params_template.branch;
  update_branch(branch, proposed_long, coef_re,
                params_template.biomarker_random,
                params_template.velocity_random);

  // Outer tape: AD over theta = [baseline, hazard, longitudinal, initial_state]
  ADvector ad_theta(n_params);
  std::copy(theta_vec.begin(), theta_vec.end(), ad_theta.begin());
  CppAD::Independent(ad_theta);

  ODEParams<ADdouble> ode;
  fill_coefs<ADdouble>(ode, ad_theta, n_baseline, n_hazard, n_long, n_init);
  load_config(ode, parameters);
  ode.branch = branch;
  load_subject(ode.subject, subject_data,
      std::vector<ADdouble>(random_effect.begin() + 2, random_effect.end()));

  ode.log_sigma_e = ADdouble(params_template.log_sigma_e);

  const auto times = time_grid(ode.subject.longitudinal_times,
                                  ode.subject.event_time);
  const std::vector<ADdouble> y0 = {ADdouble(0.0),
      ode.initial_state_coefs[0] + ADdouble(random_effect[0]),
      ode.initial_state_coefs[1] + ADdouble(random_effect[1])};
  const auto sol = ode_solve_joint(y0, times, ode);

  BSplineWorkspace ws;
  bspline_basis(0.0, ode.spline_degree, ode.spline_knots,
      ode.spline_boundary, ws.basis, ws.knots, ws.work1, ws.work2, false);
  ADdouble ll = joint_loglik(sol, times, ode, ws);

  double qf = 0.0;
  for (int j = 0; j < n_re; j++) {
    double s = 0.0;
    for (int k = 0; k < n_re; k++) s += inv_sigma_b(j, k) * random_effect[k];
    qf += random_effect[j] * s;
  }
  ll -= ADdouble(0.5 * qf + re_const);

  // Finish outer tape: NLL = -ll
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

// [[Rcpp::export(.compute_joint_objective)]]
NumericVector compute_joint_objective(
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
  const int n_init = params_template.initial_state_coefs.size();
  const int n_params = n_baseline + n_hazard + n_long + n_init;

  mat inv_sigma_b;
  if (!inv(inv_sigma_b, params_template.random_effect_sigma))
    stop("Failed to invert random effect covariance matrix");
  assert_spd(params_template.random_effect_sigma, SPD_ERR_SIGMA);
  double log_det_val, sign;
  log_det(log_det_val, sign, params_template.random_effect_sigma);
  const double re_const =
      0.5 * (inv_sigma_b.n_rows * std::log(2.0 * M_PI) + log_det_val);

  const std::vector<double> theta_vec(params.begin(), params.end());
  double total_obj = 0.0;
  std::vector<double> total_grad(n_params, 0.0);
  std::vector<double> total_hess(n_params * n_params, 0.0);

  // Value-only fast path: pure double, no AD tape
  if (!gradient && !hessian) {
    const int n_re = random_effects.ncol();
    ODEParams<double> ode_val;
    load_params(ode_val, parameters);

    for (int i = 0; i < n_subjects; i++) {
      NumericVector re = random_effects(i, _);
      NumericVector cr(re.begin() + 2, re.end());
      load_subject(ode_val.subject, data_list[i],
          std::vector<double>(cr.begin(), cr.end()));
      NumericVector lc = as<List>(parameters["coefficients"])["longitudinal"];
      update_branch(ode_val.branch, lc, cr,
                    ode_val.biomarker_random, ode_val.velocity_random);

      std::vector<double> re_vec(re.begin(), re.end());
      double nll_i = -eval_logpost(re_vec, ode_val, inv_sigma_b, re_const);

      total_obj += subject_weights[i] * nll_i;
    }

    NumericVector result(1);
    result[0] = total_obj;
    return result;
  }

  for (int i = 0; i < n_subjects; i++) {
    double subj_obj = 0.0;
    eval_subject_nll(theta_vec, data_list[i], random_effects(i, _),
        parameters, params_template,
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
// Joint: E-step log-posterior w.r.t. random effects b
// ============================================================================

// Shared setup: load params, compute sigma constants, return loaded params
static ODEParams<double> setup_logpost(
    const List& parameters, const NumericVector& random_effect,
    mat& inv_sigma_b, double& re_const) {
  ODEParams<double> params;
  load_params(params, parameters);
  if (!inv(inv_sigma_b, params.random_effect_sigma))
    stop("Failed to invert random effect covariance matrix");
  assert_spd(params.random_effect_sigma, SPD_ERR_SIGMA);
  double ld, s;
  log_det(ld, s, params.random_effect_sigma);
  re_const = 0.5 * (random_effect.size() * std::log(2.0 * M_PI) + ld);
  NumericVector lc = as<List>(parameters["coefficients"])["longitudinal"];
  NumericVector cr(random_effect.begin() + 2, random_effect.end());
  update_branch(params.branch, lc, cr,
    params.biomarker_random, params.velocity_random);
  return params;
}

// [[Rcpp::export(.compute_joint_logpost)]]
NumericVector compute_joint_logpost(
    const NumericVector& random_effect,
    const List& data, const List& parameters,
    bool gradient = true, bool hessian = false) {
  const int n_re = random_effect.size();
  const std::vector<double> re_vec(random_effect.begin(), random_effect.end());

  mat inv_sigma_b;
  double re_const;
  ODEParams<double> base_params = setup_logpost(parameters, random_effect,
    inv_sigma_b, re_const);

  // Value-only: pure double, no AD tape
  if (!gradient && !hessian) {
    load_subject(base_params.subject, data,
        std::vector<double>(re_vec.begin() + 2, re_vec.end()));
    NumericVector result(1);
    result[0] = eval_logpost(re_vec, base_params, inv_sigma_b, re_const);
    return result;
  }

  // AD path
  ADvector ad_re(n_re);
  std::copy(re_vec.begin(), re_vec.end(), ad_re.begin());
  CppAD::Independent(ad_re);

  ODEParams<ADdouble> ode;
  load_subject(ode.subject, data,
      std::vector<ADdouble>(ad_re.begin() + 2, ad_re.end()));
  load_params(ode, parameters);
  ode.branch = base_params.branch;
  std::vector<ADdouble> ad_re_vec(ad_re.begin(), ad_re.end());
  ADdouble lp = eval_logpost(ad_re_vec, ode, inv_sigma_b, re_const);

  CppAD::ADFun<double> tape;
  tape.Dependent(ad_re, ADvector{lp});
  return eval_tape(tape, re_vec, n_re, gradient, hessian);
}

// ============================================================================
// Joint: state optimization w.r.t. initial state [m(0), v(0)]
// ============================================================================

// [[Rcpp::export(.compute_joint_state)]]
NumericVector compute_joint_state(
    const NumericVector& initial_state,
    const List& data,
    const NumericVector& random_effect,
    const List& parameters,
    bool gradient = true, bool hessian = false) {

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

  const auto times = time_grid(ode.subject.longitudinal_times,
                                  ode.subject.event_time);
  const std::vector<ADdouble> y0 = {ADdouble(0.0), ad_state[0], ad_state[1]};
  const auto sol = ode_solve_joint(y0, times, ode);

  BSplineWorkspace ws;
  bspline_basis(0.0, ode.spline_degree, ode.spline_knots,
      ode.spline_boundary, ws.basis, ws.knots, ws.work1, ws.work2, false);
  ADdouble ll = joint_loglik(sol, times, ode, ws);

  CppAD::ADFun<double> tape;
  tape.Dependent(ad_state, ADvector{ll});
  return eval_tape(tape,
      std::vector<double>(initial_state.begin(), initial_state.end()),
      2, gradient, hessian);
}

// ============================================================================
// Joint: batch solver for prediction (no AD)
// ============================================================================

// [[Rcpp::export(.solve_batch_joint)]]
List solve_batch_joint(const List& data_list,
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
    // RE layout: [b_m0, b_v0, b_coef...]
    NumericVector coef_re(re.begin() + 2, re.end());
    load_subject(params.subject, data_list[i],
        std::vector<double>(coef_re.begin(), coef_re.end()));
    NumericVector long_coefs = as<List>(parameters["coefficients"])["longitudinal"];
    update_branch(params.branch, long_coefs, coef_re,
                  params.biomarker_random, params.velocity_random);

    // y0 = [H=0, initial_state[0] + b_m0, initial_state[1] + b_v0]
    const std::vector<double> y0 = {0.0,
        params.initial_state_coefs[0] + re[0],
        params.initial_state_coefs[1] + re[1]};
    const auto times = time_grid(params.subject.longitudinal_times,
                                    params.subject.event_time);
    const auto sol = ode_solve_joint(y0, times, params);

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

// ============================================================================
// Marginal: objective w.r.t. theta (AD over theta)
// ============================================================================

// [[Rcpp::export(.compute_marginal_objective)]]
NumericVector compute_marginal_objective(
    const NumericVector& params,
    const List& data_list,
    double biomarker_clamp = 50.0,
    bool gradient = true,
    bool hessian = false) {

  const int n_params = params.size();
  const MatExpBranch branch = classify_disc(params[0], params[1]);

  ADvector ad_theta(n_params);
  std::copy(params.begin(), params.end(), ad_theta.begin());
  CppAD::Independent(ad_theta);

  std::vector<ADdouble> theta_vec(ad_theta.begin(), ad_theta.end());
  MarginalParams<ADdouble> mp;
  load_marginal_params(mp, theta_vec, branch, biomarker_clamp);

  ADdouble total(0.0);
  const int n_subjects = data_list.size();

  for (int i = 0; i < n_subjects; i++) {
    load_subject(mp.subject, data_list[i], std::vector<ADdouble>{});
    if (mp.subject.longitudinal_measurements.empty()) continue;

    std::vector<ADdouble> y0 = {
      ADdouble(mp.subject.initial_state[0]),
      ADdouble(mp.subject.initial_state[1])
    };
    const auto sol = ode_solve_marginal(
        y0, mp.subject.longitudinal_times, mp);
    total -= marginal_loglik(sol, mp);
  }

  CppAD::ADFun<double> tape;
  tape.Dependent(ad_theta, ADvector{total});
  return eval_tape(tape, std::vector<double>(params.begin(), params.end()),
                   n_params, gradient, hessian);
}

// ============================================================================
// Marginal: state optimization (AD over [m(0), v(0)])
// ============================================================================

// [[Rcpp::export(.compute_marginal_state)]]
NumericVector compute_marginal_state(
    const NumericVector& initial_state,
    const List& subject_data,
    const NumericVector& parameters,
    double biomarker_clamp = 50.0,
    bool gradient = true,
    bool hessian = false) {

  if (initial_state.size() != 2) stop("initial_state must have length 2");

  const MatExpBranch branch = classify_disc(parameters[0], parameters[1]);

  MarginalParams<ADdouble> mp;
  std::vector<ADdouble> theta_ad(parameters.begin(), parameters.end());
  load_marginal_params(mp, theta_ad, branch, biomarker_clamp);

  ADvector ad_state(2);
  std::copy(initial_state.begin(), initial_state.end(), ad_state.begin());
  CppAD::Independent(ad_state);

  load_subject(mp.subject, subject_data, std::vector<ADdouble>{});
  if (mp.subject.longitudinal_measurements.empty())
    stop("subject must have at least one observation");

  std::vector<ADdouble> y0 = {ad_state[0], ad_state[1]};
  const auto sol = ode_solve_marginal(
      y0, mp.subject.longitudinal_times, mp);
  ADdouble ll = marginal_loglik(sol, mp);

  CppAD::ADFun<double> tape;
  tape.Dependent(ad_state, ADvector{ll});
  return eval_tape(tape,
      std::vector<double>(initial_state.begin(), initial_state.end()),
      2, gradient, hessian);
}

// ============================================================================
// Marginal: batch solver for prediction (no AD)
// ============================================================================

// [[Rcpp::export(.solve_batch_marginal)]]
List solve_batch_marginal(
    const List& data_list,
    const NumericVector& parameters,
    double biomarker_clamp = 50.0) {

  const int n_subjects = data_list.size();
  std::vector<double> pvec(parameters.begin(), parameters.end());
  const MatExpBranch branch = classify_disc(pvec[0], pvec[1]);
  MarginalParams<double> mp;
  load_marginal_params(mp, pvec, branch, biomarker_clamp);
  const double b1 = pvec[0], b2 = pvec[1];
  std::vector<double> fc(pvec.begin() + 2, pvec.end());

  List results(n_subjects);

  for (int i = 0; i < n_subjects; i++) {
    load_subject(mp.subject, data_list[i], std::vector<double>{});

    const auto& times = mp.subject.longitudinal_times;
    const int nt = times.size();
    std::vector<double> y0 = {
      mp.subject.initial_state[0], mp.subject.initial_state[1]
    };
    const auto sol = ode_solve_marginal(y0, times, mp);

    NumericVector biomarker(nt), velocity(nt), accel(nt);
    for (int t = 0; t < nt; t++) {
      biomarker[t] = sol[t][0];
      velocity[t] = sol[t][1];

      std::vector<double> cf;
      covariates_at(times[t], times,
                    mp.subject.longitudinal_covariates_fixed, cf);
      accel[t] = b1 * sol[t][0] + b2 * sol[t][1];
      for (size_t c = 0; c < cf.size(); c++)
        accel[t] += cf[c] * fc[c];
    }

    results[i] = List::create(
        Named("times") = wrap(times),
        Named("biomarker") = biomarker,
        Named("velocity") = velocity,
        Named("acceleration") = accel);
  }
  return results;
}
