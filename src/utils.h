#ifndef UTILS_HPP
#define UTILS_HPP

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <cppad/cppad.hpp>
#include <vector>

using namespace Rcpp;
using namespace arma;
using CppAD::AD;

typedef AD<double> ADdouble;
typedef CppAD::vector<ADdouble> ADvector;

// ============================================================================
// Constants
// ============================================================================

const double HAZARD_CLAMP_MAX = 20.0;
const double HAZARD_CLAMP_MIN = -20.0;
const double SPLINE_DENOM_TOL = 1e-8;
const double TIME_MATCH_TOL = 1e-10;
const size_t MIN_ODE_STEPS = 2;
const double STEPS_PER_UNIT_TIME = 5.0;

// ============================================================================
// Basic utility functions
// ============================================================================

// Extract scalar value from AD or regular type
template <typename Scalar>
inline double get_value(const Scalar& x) {
  return CppAD::Value(x);
}

template <>
inline double get_value<double>(const double& x) {
  return x;
}

// Clamp value using CppAD conditional expressions
template <typename Scalar>
inline Scalar clamp_value(const Scalar& value, const Scalar& min_val,
                          const Scalar& max_val) {
  return CppAD::CondExpGt(value, max_val, max_val,
                          CppAD::CondExpLt(value, min_val, min_val, value));
}

// Convert Rcpp vector to std::vector<double>
template <typename T>
inline std::vector<double> to_double_vec(const T& v) {
  return std::vector<double>(v.begin(), v.end());
}

// ============================================================================
// Data structures
// ============================================================================

// Subject-level data
template <typename Scalar = double>
struct SubjectData {
  double event_time;
  int status;
  std::vector<double> initial_state;
  std::vector<Scalar> random_effect;
  std::vector<double> longitudinal_times;
  std::vector<double> longitudinal_measurements;
  std::vector<std::vector<double>> longitudinal_covariates_fixed;
  std::vector<std::vector<double>> longitudinal_covariates_random;
  std::vector<double> survival_covariates;
};

// ODE parameters combining subject data and model coefficients
template <typename Scalar>
struct ODEParams {
  SubjectData<Scalar> subject;

  // Coefficients
  std::vector<Scalar> longitudinal_coefs;
  std::vector<Scalar> baseline_coefs;
  std::vector<Scalar> hazard_coefs;
  double measurement_error_sd;
  arma::mat random_effect_sigma;

  // Configurations
  int spline_degree;
  std::vector<double> spline_knots;
  std::vector<double> spline_boundary;
  double gamma;
  bool biomarker_fixed;
  bool biomarker_random;
  bool velocity_fixed;
  bool velocity_random;
};

// Workspace for B-spline computation
class BSplineWorkspace {
 public:
  std::vector<double> basis;
  std::vector<double> knots;
  std::vector<double> work1;
  std::vector<double> work2;
  BSplineWorkspace() {
    basis.reserve(20);
    knots.reserve(50);
    work1.reserve(20);
    work2.reserve(20);
  }

  // Explicitly clear memory when workspace is no longer needed
  void clear() {
    std::vector<double>().swap(basis);
    std::vector<double>().swap(knots);
    std::vector<double>().swap(work1);
    std::vector<double>().swap(work2);
  }
};

// ============================================================================
// Data extraction from R Lists
// ============================================================================

// Fill subject data from R List
template <typename Scalar>
inline void fill_subject_data(SubjectData<Scalar>& subj, const List& data,
                              const std::vector<Scalar>& random_effects) {
  subj = SubjectData<Scalar>();
  subj.event_time = data["time"];
  subj.status =
      data.containsElementNamed("status") ? as<int>(data["status"]) : 0;
  subj.random_effect = random_effects;

  if (data.containsElementNamed("initial_state")) {
    NumericVector init_state = data["initial_state"];
    subj.initial_state = to_double_vec(init_state);
  }

  if (data.containsElementNamed("longitudinal")) {
    List longitudinal = data["longitudinal"];

    if (longitudinal.containsElementNamed("times")) {
      NumericVector times = longitudinal["times"];
      subj.longitudinal_times = to_double_vec(times);
    }

    if (longitudinal.containsElementNamed("measurements")) {
      NumericVector measurements = longitudinal["measurements"];
      subj.longitudinal_measurements = to_double_vec(measurements);
    }

    if (longitudinal.containsElementNamed("covariates")) {
      List covs = longitudinal["covariates"];

      if (covs.containsElementNamed("fixed")) {
        NumericMatrix covs_fixed = covs["fixed"];
        subj.longitudinal_covariates_fixed.resize(covs_fixed.nrow());
        for (int j = 0; j < covs_fixed.nrow(); j++) {
          subj.longitudinal_covariates_fixed[j].resize(covs_fixed.ncol());
          for (int k = 0; k < covs_fixed.ncol(); k++) {
            subj.longitudinal_covariates_fixed[j][k] = covs_fixed(j, k);
          }
        }
      }

      if (covs.containsElementNamed("random")) {
        NumericMatrix covs_random = covs["random"];
        subj.longitudinal_covariates_random.resize(covs_random.nrow());
        for (int j = 0; j < covs_random.nrow(); j++) {
          subj.longitudinal_covariates_random[j].resize(covs_random.ncol());
          for (int k = 0; k < covs_random.ncol(); k++) {
            subj.longitudinal_covariates_random[j][k] = covs_random(j, k);
          }
        }
      }
    }
  }

  if (data.containsElementNamed("covariates") &&
      !Rf_isNull(data["covariates"])) {
    if (Rf_isMatrix(data["covariates"])) {
      NumericMatrix cov_mat = data["covariates"];
      if (cov_mat.nrow() > 0) {
        // Extract first row as survival covariates
        NumericVector first_row = cov_mat(0, _);
        subj.survival_covariates = to_double_vec(first_row);
      }
    } else {
      NumericVector covs = data["covariates"];
      subj.survival_covariates = to_double_vec(covs);
    }
  }
}

// Fill ODE configuration from R List
template <typename Scalar>
inline void fill_ode_config(ODEParams<Scalar>& params, const List& parameters) {
  List configurations = parameters["configurations"];
  List baseline_config = configurations["baseline"];

  params.spline_degree = as<int>(baseline_config["degree"]);

  NumericVector knots = baseline_config.containsElementNamed("knots")
                            ? baseline_config["knots"]
                            : NumericVector();
  NumericVector boundary = baseline_config.containsElementNamed("boundary")
                               ? baseline_config["boundary"]
                               : baseline_config["boundary_knots"];

  params.spline_knots = to_double_vec(knots);
  params.spline_boundary = to_double_vec(boundary);
  params.gamma = as<double>(configurations["gamma"]);

  List biomarker_config = configurations["biomarker"];
  List velocity_config = configurations["velocity"];
  params.biomarker_fixed = as<bool>(biomarker_config["fixed"]);
  params.biomarker_random = as<bool>(biomarker_config["random"]);
  params.velocity_fixed = as<bool>(velocity_config["fixed"]);
  params.velocity_random = as<bool>(velocity_config["random"]);
}

// Fill all ODE parameters from R List
template <typename Scalar>
inline void fill_ode_parameters(ODEParams<Scalar>& params,
                                const List& parameters) {
  List coefficients = parameters["coefficients"];
  NumericVector baseline = coefficients["baseline"];
  NumericVector hazard = coefficients["hazard"];
  NumericVector longitudinal = coefficients["longitudinal"];

  params.baseline_coefs.resize(baseline.size());
  params.hazard_coefs.resize(hazard.size());
  params.longitudinal_coefs.resize(longitudinal.size());

  for (int i = 0; i < baseline.size(); i++)
    params.baseline_coefs[i] = Scalar(baseline[i]);
  for (int i = 0; i < hazard.size(); i++)
    params.hazard_coefs[i] = Scalar(hazard[i]);
  for (int i = 0; i < longitudinal.size(); i++)
    params.longitudinal_coefs[i] = Scalar(longitudinal[i]);

  params.measurement_error_sd =
      as<double>(coefficients["measurement_error_sd"]);
  params.random_effect_sigma =
      as<arma::mat>(coefficients["random_effect_sigma"]);

  fill_ode_config(params, parameters);
}

// ============================================================================
// Time utilities
// ============================================================================

// Build time points for ODE solving
inline std::vector<double> build_times(const std::vector<double>& obs_times,
                                       double event_time) {
  std::vector<double> times;
  times.reserve(obs_times.size() + 2);
  times.push_back(0.0);

  for (double t : obs_times) {
    times.push_back(std::max(0.0, std::min(event_time, t)));
  }
  times.push_back(event_time);

  std::sort(times.begin(), times.end());
  times.erase(std::unique(times.begin(), times.end()), times.end());

  return times;
}

// Find index of target time in sorted time vector
inline int find_time_idx(const std::vector<double>& times, double target) {
  auto it = std::lower_bound(times.begin(), times.end(), target);
  return (it != times.end() && std::abs(*it - target) < TIME_MATCH_TOL)
             ? std::distance(times.begin(), it)
             : -1;
}

// Extract covariates at a specific time
inline void extract_covariates_at_time(
    double t, const std::vector<double>& times,
    const std::vector<std::vector<double>>& covariates,
    std::vector<double>& result) {
  size_t idx;
  if (t >= times.back()) {
    idx = times.size() - 1;
  } else if (t <= times[0]) {
    idx = 0;
  } else {
    // Binary search for the interval containing t
    auto it = std::upper_bound(times.begin(), times.end(), t);
    idx = std::distance(times.begin(), it) - 1;
  }
  result = covariates[idx];
}

// ============================================================================
// Mathematical functions
// ============================================================================

// Compute B-spline basis functions
inline void compute_bspline_basis(double t, int degree,
                                  const std::vector<double>& interior_knots,
                                  const std::vector<double>& boundary,
                                  std::vector<double>& basis,
                                  std::vector<double>& knots,
                                  std::vector<double>& basis_curr,
                                  std::vector<double>& basis_prev,
                                  bool skip_knots_build = false) {
  const int n_interior_knots = interior_knots.size();
  const int n_total_knots = n_interior_knots + 2 * (degree + 1);

  if (!skip_knots_build) {
    if (static_cast<int>(knots.size()) != n_total_knots) {
      knots.resize(n_total_knots);
    }

    std::fill_n(knots.begin(), degree + 1, boundary[0]);
    std::copy(interior_knots.begin(), interior_knots.end(),
              knots.begin() + degree + 1);
    std::fill_n(knots.begin() + degree + 1 + n_interior_knots, degree + 1,
                boundary[1]);
  }

  const int n_basis_functions = n_total_knots - degree - 1;
  if (static_cast<int>(basis.size()) != n_basis_functions) {
    basis.resize(n_basis_functions);
  }

  t = std::max(boundary[0], std::min(boundary[1], t));

  int knot_span;
  if (t >= knots[n_total_knots - degree - 1]) {
    knot_span = n_total_knots - degree - 2;
  } else if (t <= knots[degree]) {
    knot_span = degree;
  } else {
    knot_span = degree;
    for (int i = degree; i < n_total_knots - degree - 1; ++i) {
      if (t < knots[i + 1]) {
        knot_span = i;
        break;
      }
    }
  }

  if (static_cast<int>(basis_curr.size()) != n_basis_functions) {
    basis_curr.resize(n_basis_functions);
    basis_prev.resize(n_basis_functions);
  }

  std::fill(basis_curr.begin(), basis_curr.end(), 0.0);
  if (knot_span < n_basis_functions) {
    basis_curr[knot_span] = 1.0;
  }
  std::fill(basis_prev.begin(), basis_prev.end(), 0.0);

  for (int order = 1; order <= degree; order++) {
    std::swap(basis_curr, basis_prev);
    std::fill(basis_curr.begin(), basis_curr.end(), 0.0);

    for (int i = 0; i < n_basis_functions; i++) {
      if (i < n_total_knots - order) {
        const double denominator = knots[i + order] - knots[i];
        if (denominator > SPLINE_DENOM_TOL) {
          basis_curr[i] += basis_prev[i] * (t - knots[i]) / denominator;
        }
      }

      if (i + 1 < n_total_knots - order) {
        const double denominator = knots[i + order + 1] - knots[i + 1];
        if (denominator > SPLINE_DENOM_TOL) {
          basis_curr[i] +=
              basis_prev[i + 1] * (knots[i + order + 1] - t) / denominator;
        }
      }
    }
  }

  std::copy(basis_curr.begin(), basis_curr.end(), basis.begin());
}

// Compute log hazard function
template <typename Scalar>
inline Scalar compute_log_hazard(Scalar biomarker, Scalar velocity,
                                 const std::vector<double>& basis_lambda,
                                 const std::vector<Scalar>& baseline_coefs,
                                 const std::vector<Scalar>& hazard_coefs,
                                 const std::vector<double>& surv_covs,
                                 double gamma) {
  Scalar log_baseline = Scalar(0.0);
  for (size_t i = 0; i < basis_lambda.size(); ++i) {
    log_baseline += Scalar(basis_lambda[i]) * baseline_coefs[i];
  }

  Scalar v_power;
  if (gamma == 0) {
    v_power = Scalar(0.0);
  } else if (gamma == 1) {
    v_power = velocity;
  } else {
    v_power = velocity * velocity;
  }

  Scalar log_assoc = hazard_coefs[0] * biomarker + hazard_coefs[1] * v_power;

  for (size_t i = 0; i < surv_covs.size(); i++) {
    log_assoc += hazard_coefs[i + 2] * Scalar(surv_covs[i]);
  }

  return log_baseline + log_assoc;
}

// Compute acceleration (biomarker second derivative)
template <typename Scalar>
inline Scalar compute_acceleration(
    Scalar biomarker, Scalar velocity, double t,
    const SubjectData<Scalar>& subj,
    const std::vector<Scalar>& longitudinal_coefs, bool biomarker_fixed,
    bool biomarker_random, bool velocity_fixed, bool velocity_random) {
  std::vector<double> covs_fixed, covs_random;
  extract_covariates_at_time(t, subj.longitudinal_times,
                             subj.longitudinal_covariates_fixed, covs_fixed);
  extract_covariates_at_time(t, subj.longitudinal_times,
                             subj.longitudinal_covariates_random, covs_random);

  Scalar acceleration = Scalar(0.0);
  size_t f_idx = 0, r_idx = 0;

  Scalar biomarker_coef = Scalar(0.0);
  if (biomarker_fixed) biomarker_coef = longitudinal_coefs[f_idx++];
  if (biomarker_random) biomarker_coef += subj.random_effect[r_idx++];
  if (biomarker_fixed || biomarker_random)
    acceleration += biomarker_coef * biomarker;

  Scalar velocity_coef = Scalar(0.0);
  if (velocity_fixed) velocity_coef = longitudinal_coefs[f_idx++];
  if (velocity_random) velocity_coef += subj.random_effect[r_idx++];
  if (velocity_fixed || velocity_random)
    acceleration += velocity_coef * velocity;

  for (size_t i = 0; i < covs_fixed.size(); i++)
    acceleration += longitudinal_coefs[f_idx++] * Scalar(covs_fixed[i]);

  for (size_t i = 0; i < covs_random.size(); i++)
    acceleration += subj.random_effect[r_idx++] * Scalar(covs_random[i]);

  return acceleration;
}

// ============================================================================
// ODE solving
// ============================================================================

// ODE functor for Runge-Kutta solver
template <typename Scalar>
class ODEFunctor {
 private:
  const ODEParams<Scalar>& params_;
  BSplineWorkspace workspace_;

 public:
  explicit ODEFunctor(const ODEParams<Scalar>& params)
      : params_(params), workspace_() {
    // Pre-compute knots once (they are constant during ODE solving)
    compute_bspline_basis(0.0, params_.spline_degree, params_.spline_knots,
                          params_.spline_boundary, workspace_.basis,
                          workspace_.knots, workspace_.work1, workspace_.work2,
                          false);  // Build knots
  }

  void Ode(const Scalar& t, const CppAD::vector<Scalar>& y,
           CppAD::vector<Scalar>& dydt) {
    const Scalar biomarker = y[1];
    const Scalar velocity = y[2];

    double t_val = get_value(t);
    const Scalar acceleration = compute_acceleration(
        biomarker, velocity, t_val, params_.subject, params_.longitudinal_coefs,
        params_.biomarker_fixed, params_.biomarker_random,
        params_.velocity_fixed, params_.velocity_random);
    compute_bspline_basis(t_val, params_.spline_degree, params_.spline_knots,
                          params_.spline_boundary, workspace_.basis,
                          workspace_.knots, workspace_.work1, workspace_.work2,
                          true);  // Skip knots rebuild (already initialized)

    const Scalar log_hazard =
        compute_log_hazard(biomarker, velocity, workspace_.basis,
                           params_.baseline_coefs, params_.hazard_coefs,
                           params_.subject.survival_covariates, params_.gamma);

    const Scalar clamped_log_hazard = clamp_value(
        log_hazard, Scalar(HAZARD_CLAMP_MIN), Scalar(HAZARD_CLAMP_MAX));

    const Scalar hazard = CppAD::exp(clamped_log_hazard);

    dydt[0] = hazard;
    dydt[1] = velocity;
    dydt[2] = acceleration;
  }
};

// Solve ODE system using Runge-Kutta method
template <typename Scalar>
inline std::vector<std::vector<Scalar>> solve_ode(
    const std::vector<Scalar>& y0, const std::vector<double>& times,
    const ODEParams<Scalar>& params) {
  const int n_state = y0.size();
  const int n_times = times.size();

  ODEFunctor<Scalar> ode_fun(params);

  std::vector<std::vector<Scalar>> solution(n_times,
                                            std::vector<Scalar>(n_state));
  CppAD::vector<Scalar> yi(n_state);

  // Initialize with initial conditions - copy once
  std::copy(y0.begin(), y0.end(), yi.begin());
  std::copy(y0.begin(), y0.end(), solution[0].begin());

  Scalar ti = Scalar(times[0]);
  Scalar tf, dt;

  for (int idx = 1; idx < n_times; idx++) {
    tf = Scalar(times[idx]);
    dt = tf - ti;

    const size_t M = std::max(
        MIN_ODE_STEPS,
        static_cast<size_t>(std::ceil(get_value(dt) * STEPS_PER_UNIT_TIME)));

    yi = CppAD::Runge45(ode_fun, M, ti, tf, yi);

    for (int i = 0; i < n_state; ++i) {
      solution[idx][i] = yi[i];
    }

    ti = tf;
  }

  // Clear temporary vector memory explicitly
  CppAD::vector<Scalar>().swap(yi);

  return solution;
}

#endif  // UTILS_HPP
