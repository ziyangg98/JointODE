#include <RcppArmadillo.h>

#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export(.solve_batch_ode_cppad)]]
List solve_batch_ode_cppad(const List& data_list,
                           const NumericMatrix& random_effects,
                           const List& parameters) {
  const int n_subjects = data_list.size();

  ODEParams<double> params;
  fill_ode_parameters(params, parameters);

  BSplineWorkspace workspace;
  compute_bspline_basis(0.0, params.spline_degree, params.spline_knots,
                        params.spline_boundary, workspace.basis, workspace.knots,
                        workspace.work1, workspace.work2, false);

  List results(n_subjects);

  for (int i = 0; i < n_subjects; i++) {
    NumericVector re = random_effects(i, _);
    fill_subject_data(params.subject, data_list[i],
                      std::vector<double>(re.begin(), re.end()));

    const std::vector<double> y0 = {0.0, params.subject.initial_state[0],
                                    params.subject.initial_state[1]};
    const std::vector<double> subject_times = build_times(
        params.subject.longitudinal_times, params.subject.event_time);

    const std::vector<std::vector<double>> subject_solution =
        solve_ode(y0, subject_times, params);

    const int n_times = subject_times.size();
    NumericVector times(n_times);
    NumericVector biomarker_traj(n_times);
    NumericVector velocity_traj(n_times);
    NumericVector acceleration_traj(n_times);
    NumericVector cumhazard_traj(n_times);
    NumericVector log_hazard_traj(n_times);

    for (int t_idx = 0; t_idx < n_times; t_idx++) {
      const double t = subject_times[t_idx];
      const auto& state = subject_solution[t_idx];
      const double biomarker = state[1];
      const double velocity = state[2];

      times[t_idx] = t;
      cumhazard_traj[t_idx] = state[0];
      biomarker_traj[t_idx] = biomarker;
      velocity_traj[t_idx] = velocity;

      acceleration_traj[t_idx] = compute_acceleration(
          biomarker, velocity, t, params.subject, params.longitudinal_coefs,
          params.biomarker_fixed, params.biomarker_random,
          params.velocity_fixed, params.velocity_random);

      compute_bspline_basis(t, params.spline_degree, params.spline_knots,
                            params.spline_boundary, workspace.basis,
                            workspace.knots, workspace.work1, workspace.work2, true);

      const double log_hazard_raw =
          compute_log_hazard(biomarker, velocity, workspace.basis,
                             params.baseline_coefs, params.hazard_coefs,
                             params.subject.survival_covariates, params.gamma);
      log_hazard_traj[t_idx] = std::max(
          HAZARD_CLAMP_MIN, std::min(log_hazard_raw, HAZARD_CLAMP_MAX));
    }

    results[i] = List::create(Named("times") = times,
                              Named("biomarker") = biomarker_traj,
                              Named("velocity") = velocity_traj,
                              Named("acceleration") = acceleration_traj,
                              Named("cum_hazard") = cumhazard_traj,
                              Named("log_hazard") = log_hazard_traj);
  }

  return results;
}
