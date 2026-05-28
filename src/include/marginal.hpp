// MarginalODE.hpp — TMB objective for longitudinal-only ODE model
// Supports random effects (Laplace) but no survival component
#ifndef MARGINAL_ODE_HPP
#define MARGINAL_ODE_HPP

#include "utils.hpp"

template<class Type>
Type marginal_ode_nll(objective_function<Type>* obj) {
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

  // --- DATA ---
  DATA_INTEGER(n_subjects);
  DATA_INTEGER(n_random_effects);
  DATA_IVECTOR(n_observations);
  DATA_VECTOR(obs_times);
  DATA_VECTOR(obs_values);
  DATA_INTEGER(n_fixed_covariates);
  DATA_INTEGER(n_random_covariates);
  DATA_VECTOR(long_fixed_covariates);
  DATA_VECTOR(long_random_covariates);
  DATA_INTEGER(biomarker_fixed);
  DATA_INTEGER(biomarker_random);
  DATA_INTEGER(velocity_fixed);
  DATA_INTEGER(velocity_random);
  DATA_INTEGER(diagonal_re);

  // Pre-compute per-subject offsets
  vector<int> obs_offset(n_subjects), fixed_offset(n_subjects);
  vector<int> random_offset(n_subjects);
  obs_offset(0) = fixed_offset(0) = random_offset(0) = 0;
  for (int i = 1; i < n_subjects; i++) {
    int ni = n_observations(i - 1);
    obs_offset(i) = obs_offset(i - 1) + ni;
    fixed_offset(i) = fixed_offset(i - 1) + ni * n_fixed_covariates;
    random_offset(i) = random_offset(i - 1) + ni * n_random_covariates;
  }

  // --- PARAMETERS ---
  PARAMETER_VECTOR(longitudinal);
  PARAMETER_VECTOR(initial_state);
  PARAMETER(log_sigma_e);
  PARAMETER_VECTOR(log_sd_re);
  PARAMETER_VECTOR(corr_par);
  PARAMETER_MATRIX(random_effects);  // n_subjects x n_re

  using namespace density;
  UNSTRUCTURED_CORR_t<Type> corr_structure(corr_par);
  Type sigma_e = exp(log_sigma_e);

  int n_total_obs = obs_times.size();
  vector<Type> fitted_biomarker(n_total_obs), fitted_velocity(n_total_obs);

  parallel_accumulator<Type> nll(obj);

  for (int i = 0; i < n_subjects; i++) {
    int ni = n_observations(i);
    int oi = obs_offset(i), fi = fixed_offset(i);
    int ri = random_offset(i);

    vector<Type> bi = random_effects.row(i);

    // --- Random effect prior: MVN(0, Sigma_b) ---
    if (diagonal_re) {
      for (int j = 0; j < n_random_effects; j++)
        nll -= dnorm(bi(j), Type(0), exp(log_sd_re(j)), true);
    } else {
      nll += VECSCALE(corr_structure, exp(log_sd_re))(bi);
    }
    // Unpack subject data as double
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

    // Build time grid: {0, obs_times}
    double max_t = obs_t(ni - 1);
    vector<double> time_pts = build_time_grid(obs_t, max_t);
    int n_times = time_pts.size();

    // Dynamic parameters are constant across time steps.
    Type m = initial_state(0) + bi(0);
    Type v = initial_state(1) + bi(1);
    Type log_omega2(0), log_2xi_omega(0);
    int fixed_idx = 0, re_idx = 2;
    if (biomarker_fixed) log_omega2 += longitudinal(fixed_idx++);
    if (biomarker_random) log_omega2 += bi(re_idx++);
    if (velocity_fixed) log_2xi_omega += longitudinal(fixed_idx++);
    if (velocity_random) log_2xi_omega += bi(re_idx++);
    Type b1 = -exp(log_omega2);
    Type b2 = -exp(log_2xi_omega);
    int forcing_fixed_start = fixed_idx;
    int forcing_re_start = re_idx;

    vector<Type> sol_m(n_times), sol_v(n_times);
    sol_m(0) = m;  sol_v(0) = v;

    for (int ti = 1; ti < n_times; ti++) {
      double t0 = time_pts(ti - 1);
      Type dt = Type(time_pts(ti) - t0);
      int cov_idx = covariate_index_at(t0, obs_t);

      Type eta(0);
      for (int k = 0; k < n_fixed_covariates; k++)
        eta += longitudinal(forcing_fixed_start + k) * Type(long_fixed_covariates_i(cov_idx, k));
      for (int k = 0; k < n_random_covariates; k++)
        eta += bi(forcing_re_start + k) * Type(long_random_covariates_i(cov_idx, k));

      ode_step(m, v, b1, b2, eta, dt);
      sol_m(ti) = m;  sol_v(ti) = v;
    }

    // Gaussian likelihood
    for (int j = 0; j < ni; j++) {
      int time_idx = find_time_index(time_pts, obs_t(j));
      if (time_idx >= 0) {
        nll -= dnorm(Type(obs_y(j)), sol_m(time_idx), sigma_e, true);
        fitted_biomarker(oi + j) = sol_m(time_idx);
        fitted_velocity(oi + j) = sol_v(time_idx);
      }
    }
  }

  REPORT(fitted_biomarker);
  REPORT(fitted_velocity);

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

  REPORT(Sigma_b);
  ADREPORT(Sigma_b);
  ADREPORT(sigma_e);

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
  return nll;
}

#endif
