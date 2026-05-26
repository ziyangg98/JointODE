// joint.hpp — Joint model: longitudinal ODE + Cox PH survival with shared RE
#ifndef JOINT_HPP
#define JOINT_HPP

#include "utils.hpp"

template<class Type>
Type joint_ode_nll(objective_function<Type>* obj) {
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

  // --- DATA ---
  DATA_INTEGER(n_subjects);
  DATA_INTEGER(n_random_effects);
  DATA_IVECTOR(n_observations);        // per subject
  DATA_VECTOR(event_times);
  DATA_IVECTOR(event_status);
  DATA_VECTOR(obs_times);              // concatenated
  DATA_VECTOR(obs_values);             // concatenated
  DATA_INTEGER(n_fixed_covariates);
  DATA_INTEGER(n_random_covariates);
  DATA_VECTOR(long_fixed_covariates);  // row-major flat
  DATA_VECTOR(long_random_covariates);
  DATA_VECTOR(surv_covariates);
  DATA_INTEGER(n_surv_covariates);
  DATA_INTEGER(hazard_quadrature);
  DATA_SCALAR(velocity_power);
  DATA_INTEGER(baseline_degree);
  DATA_VECTOR(baseline_knots);
  DATA_VECTOR(baseline_boundary);
  DATA_INTEGER(omega_fixed);
  DATA_INTEGER(omega_random);
  DATA_INTEGER(xi_fixed);
  DATA_INTEGER(xi_random);
  DATA_INTEGER(diagonal_re);

  // Convert DATA_SCALAR/VECTOR to double for non-AD helper functions
  double velocity_power_val = asDouble(velocity_power);
  int n_quadrature = std::max(hazard_quadrature, 1);
  vector<double> knots_d(baseline_knots.size()), boundary_d(baseline_boundary.size());
  for (size_t i = 0; i < (size_t)knots_d.size(); i++) knots_d(i) = asDouble(baseline_knots(i));
  for (size_t i = 0; i < (size_t)boundary_d.size(); i++) boundary_d(i) = asDouble(baseline_boundary(i));
  vector<double> full_knots = build_knot_vector(baseline_degree, knots_d, boundary_d);

  // Pre-compute per-subject observation offsets for parallel access
  vector<int> obs_offset(n_subjects), fixed_offset(n_subjects);
  vector<int> random_offset(n_subjects), surv_offset(n_subjects);
  obs_offset(0) = fixed_offset(0) = random_offset(0) = surv_offset(0) = 0;
  for (int i = 1; i < n_subjects; i++) {
    int ni = n_observations(i - 1);
    obs_offset(i) = obs_offset(i - 1) + ni;
    fixed_offset(i) = fixed_offset(i - 1) + ni * n_fixed_covariates;
    random_offset(i) = random_offset(i - 1) + ni * n_random_covariates;
    surv_offset(i) = surv_offset(i - 1) + n_surv_covariates;
  }

  // --- PARAMETERS ---
  PARAMETER_VECTOR(baseline);
  PARAMETER_VECTOR(hazard);
  PARAMETER_VECTOR(longitudinal);
  PARAMETER_VECTOR(initial_state);
  PARAMETER(log_sigma_e);
  PARAMETER_VECTOR(log_sd_re);
  PARAMETER_VECTOR(corr_par);
  PARAMETER_MATRIX(random_effects);  // n_subjects x n_re

  // Covariance structure: VECSCALE(UNSTRUCTURED_CORR, exp(log_sd))
  using namespace density;
  UNSTRUCTURED_CORR_t<Type> corr_structure(corr_par);
  Type sigma_e = exp(log_sigma_e);

  // Output vectors for REPORT (filled per-subject, safe for parallel)
  int n_total_obs = obs_times.size();
  vector<Type> cumulative_hazard(n_subjects), log_hazard_at_event(n_subjects);
  vector<Type> fitted_biomarker(n_total_obs), fitted_velocity(n_total_obs);

  // Parallel accumulator for negative log-likelihood
  parallel_accumulator<Type> nll(obj);

  for (int i = 0; i < n_subjects; i++) {
    int ni = n_observations(i);
    int oi = obs_offset(i), fi = fixed_offset(i);
    int ri = random_offset(i), si = surv_offset(i);

    vector<Type> bi = random_effects.row(i);

    // --- Random effect prior: MVN(0, Sigma_b) ---
    if (diagonal_re) {
      for (int j = 0; j < n_random_effects; j++)
        nll -= dnorm(bi(j), Type(0), exp(log_sd_re(j)), true);
    } else {
      nll += VECSCALE(corr_structure, exp(log_sd_re))(bi);
    }

    // --- Unpack subject data as double ---
    vector<double> obs_t(ni), obs_y(ni);
    for (int j = 0; j < ni; j++) {
      obs_t(j) = asDouble(obs_times(oi + j));
      obs_y(j) = asDouble(obs_values(oi + j));
    }
    matrix<double> long_fixed_covariates_i(ni, n_fixed_covariates);
    for (int j = 0; j < ni * n_fixed_covariates; j++)
      long_fixed_covariates_i(j / n_fixed_covariates, j % n_fixed_covariates) = asDouble(long_fixed_covariates(fi + j));
    matrix<double> long_random_covariates_i(ni, n_random_covariates);
    for (int j = 0; j < ni * n_random_covariates; j++)
      long_random_covariates_i(j / n_random_covariates, j % n_random_covariates) = asDouble(long_random_covariates(ri + j));
    vector<double> surv_cov_i(n_surv_covariates);
    for (int k = 0; k < n_surv_covariates; k++)
      surv_cov_i(k) = asDouble(surv_covariates(si + k));

    double event_time = asDouble(event_times(i));
    int status = event_status(i);

    // --- ODE solve with Simpson hazard integration ---
    vector<double> time_pts = build_time_grid(obs_t, event_time);
    int n_times = time_pts.size();
    Type m = initial_state(0) + bi(0);
    Type v = initial_state(1) + bi(1);
    Type cum_haz(0);
    vector<Type> sol_m(n_times), sol_v(n_times), sol_H(n_times);
    sol_m(0) = m;  sol_v(0) = v;  sol_H(0) = cum_haz;

    // Dynamic parameters are constant across time steps.
    Type log_omega(0), log_xi(0);
    int fixed_idx = 0, re_idx = 2;
    bool omega_active = omega_fixed || omega_random;
    bool xi_active = xi_fixed || xi_random;
    if (omega_fixed)  log_omega += longitudinal(fixed_idx++);
    if (omega_random) log_omega += bi(re_idx++);
    if (xi_fixed)     log_xi += longitudinal(fixed_idx++);
    if (xi_random)    log_xi += bi(re_idx++);
    Type b1(0), b2(0), omega2(0);
    if (omega_active) {
      Type omega = exp(log_omega);
      omega2 = omega * omega;
      b1 = -omega2;
      if (xi_active) b2 = -Type(2) * exp(log_xi) * omega;
    }
    int forcing_fixed_start = fixed_idx;
    int forcing_re_start = re_idx;

    for (int ti = 1; ti < n_times; ti++) {
      double t0 = time_pts(ti - 1);
      double dt = time_pts(ti) - t0;
      double sub_dt = dt / n_quadrature;
      int cov_idx = covariate_index_at(t0, obs_t);

      // The forcing model is a target level mu(t); the ODE uses omega^2 * mu(t).
      Type mu(0);
      for (int k = 0; k < n_fixed_covariates; k++)
        mu += longitudinal(forcing_fixed_start + k) * Type(long_fixed_covariates_i(cov_idx, k));
      for (int k = 0; k < n_random_covariates; k++)
        mu += bi(forcing_re_start + k) * Type(long_random_covariates_i(cov_idx, k));
      Type forcing = omega2 * mu;

      // Composite Simpson's rule for cumulative hazard
      Type h_left = eval_hazard(m, v, t0, baseline_degree, full_knots,
                                baseline, hazard, surv_cov_i, velocity_power_val);
      for (int q = 0; q < n_quadrature; q++) {
        double t_start = t0 + q * sub_dt;
        // Midpoint
        Type m_mid = m, v_mid = v;
        ode_step(m_mid, v_mid, b1, b2, forcing, Type(sub_dt * 0.5));
        Type h_mid = eval_hazard(m_mid, v_mid, t_start + sub_dt * 0.5,
                                 baseline_degree, full_knots,
                                 baseline, hazard, surv_cov_i, velocity_power_val);
        // Right endpoint
        ode_step(m_mid, v_mid, b1, b2, forcing, Type(sub_dt * 0.5));
        Type h_right = eval_hazard(m_mid, v_mid, t_start + sub_dt,
                                   baseline_degree, full_knots,
                                   baseline, hazard, surv_cov_i, velocity_power_val);
        cum_haz += Type(sub_dt / 6.0) * (h_left + Type(4) * h_mid + h_right);
        m = m_mid;  v = v_mid;  h_left = h_right;
      }
      sol_m(ti) = m;  sol_v(ti) = v;  sol_H(ti) = cum_haz;
    }

    // --- Longitudinal likelihood: Gaussian via dnorm ---
    for (int j = 0; j < ni; j++) {
      int time_idx = find_time_index(time_pts, obs_t(j));
      if (time_idx >= 0) {
        nll -= dnorm(Type(obs_y(j)), sol_m(time_idx), sigma_e, true);
        fitted_biomarker(oi + j) = sol_m(time_idx);
        fitted_velocity(oi + j) = sol_v(time_idx);
      }
    }

    // --- Survival likelihood: Cox PH ---
    int event_idx = find_time_index(time_pts, event_time);
    if (event_idx >= 0) {
      nll += sol_H(event_idx);
      Type h_event = eval_hazard(sol_m(event_idx), sol_v(event_idx),
                                 event_time, baseline_degree, full_knots,
                                 baseline, hazard, surv_cov_i, velocity_power_val);
      Type log_h_event = log(h_event);
      if (status == 1) nll -= log_h_event;
      cumulative_hazard(i) = sol_H(event_idx);
      log_hazard_at_event(i) = log_h_event;
    }
  }

  // Reconstruct Sigma_b for reporting
  matrix<Type> Sigma_b(n_random_effects, n_random_effects);
  Sigma_b.setZero();
  if (diagonal_re) {
    for (int i = 0; i < n_random_effects; i++)
      Sigma_b(i, i) = exp(log_sd_re(i)) * exp(log_sd_re(i));
  } else {
    Sigma_b = corr_structure.cov();
    for (int i = 0; i < n_random_effects; i++)
      for (int j = 0; j < n_random_effects; j++)
        Sigma_b(i, j) *= exp(log_sd_re(i)) * exp(log_sd_re(j));
  }

  REPORT(fitted_biomarker);
  REPORT(fitted_velocity);
  REPORT(cumulative_hazard);
  REPORT(log_hazard_at_event);
  REPORT(Sigma_b);
  ADREPORT(Sigma_b);
  ADREPORT(sigma_e);

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
  return nll;
}

#endif
