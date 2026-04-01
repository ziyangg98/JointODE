#ifndef SOLVER_H
#define SOLVER_H

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <cppad/cppad.hpp>
#include <vector>

using namespace Rcpp;
using namespace arma;
using CppAD::AD;

typedef AD<double> ADdouble;
typedef AD<ADdouble> AD2double;
typedef CppAD::vector<ADdouble> ADvector;

// ============================================================================
// Constants
// ============================================================================

const double HAZARD_CLAMP_MAX = 20.0;
const double HAZARD_CLAMP_MIN = -20.0;
const double BIOMARKER_CLAMP_DEFAULT = 50.0;  // fallback if not set from data
const double SPLINE_DENOM_TOL = 1e-8;
const double TIME_MATCH_TOL = 1e-10;
const double DISC_TOL = 1e-12;

// Branch type for matrix exponential
// ODE root types for d2m/dt2 = b1*m + b2*m' + f
// Roots of r^2 - b2*r - b1 = 0, discriminant D = b2^2 + 4*b1
//   REAL:       D > 0, |b1| > eps (both roots well-separated from 0)
//   FIRST_ORD:  |b1| ~ 0 (degenerates to first-order ODE v' = b2*v + f)
//   COMPLEX:    D < 0
//   REPEATED:   D ~ 0, |b2| > eps (double root at b2/2)
//   ZERO:       D ~ 0, |b2| ~ 0 (trivial: m'' ≈ f)
enum class MatExpBranch { REAL, FIRST_ORD, COMPLEX, REPEATED, ZERO };

inline MatExpBranch classify_disc(double b1, double b2) {
  if (std::abs(b1) < 1e-8) {
    return std::abs(b2) > 1e-8 ? MatExpBranch::FIRST_ORD : MatExpBranch::ZERO;
  }
  double D = b2 * b2 + 4.0 * b1;
  if (D > DISC_TOL) return MatExpBranch::REAL;
  if (D < -DISC_TOL) return MatExpBranch::COMPLEX;
  return MatExpBranch::REPEATED;
}


// ============================================================================
// Basic utility functions
// ============================================================================

// Extract scalar value from AD or regular type (recursive unwrap)
inline double get_value(double x) { return x; }
template <typename T>
inline double get_value(const CppAD::AD<T>& x) { return get_value(CppAD::Value(x)); }

// Clamp value using CppAD conditional expressions
template <typename Scalar>
inline Scalar clamp(const Scalar& value, const Scalar& min_val,
                          const Scalar& max_val) {
  return CppAD::CondExpGt(value, max_val, max_val,
                          CppAD::CondExpLt(value, min_val, min_val, value));
}

// Safe exp: linear extension beyond threshold preserves gradient for AD
template <typename Scalar>
inline Scalar safe_exp(const Scalar& x) {
  const Scalar M(500.0);
  const Scalar capped = CppAD::CondExpGt(x, M, M, x);
  const Scalar excess = CppAD::CondExpGt(x, M, x - M, Scalar(0));
  return CppAD::exp(capped) * (Scalar(1.0) + excess);
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
  std::vector<Scalar> initial_state_coefs;  // population mean [m(0), v(0)]
  Scalar log_sigma_e;                       // log measurement error SD
  double measurement_error_sd;              // double copy for non-AD paths
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
  MatExpBranch branch;
  double biomarker_clamp;
  int hazard_quadrature;  // Simpson sub-intervals for hazard integration
};

// Marginal ODE parameters (longitudinal only, no survival)
template <typename Scalar>
struct MarginalParams {
  SubjectData<Scalar> subject;
  std::vector<Scalar> parameters;  // [b1, b2, cov_coefs...]
  MatExpBranch branch;
  double biomarker_clamp;
};

// Workspace for B-spline computation
class BSplineWorkspace {
 public:
  std::vector<double> basis, knots, work1, work2;
};

// ============================================================================
// Data extraction from R Lists
// ============================================================================

// Fill subject data from R List
template <typename Scalar>
inline void load_subject(SubjectData<Scalar>& subj, const List& data,
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
inline void load_config(ODEParams<Scalar>& params, const List& parameters) {
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

  params.biomarker_clamp = configurations.containsElementNamed("biomarker_clamp")
    ? as<double>(configurations["biomarker_clamp"])
    : BIOMARKER_CLAMP_DEFAULT;
  params.hazard_quadrature = configurations.containsElementNamed("hazard_quadrature")
    ? as<int>(configurations["hazard_quadrature"])
    : 1;
}

// Fill all ODE parameters from R List
template <typename Scalar>
inline void load_params(ODEParams<Scalar>& params,
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

  NumericVector initial_state = coefficients["initial_state"];
  params.initial_state_coefs.resize(initial_state.size());
  for (int i = 0; i < initial_state.size(); i++)
    params.initial_state_coefs[i] = Scalar(initial_state[i]);

  params.measurement_error_sd =
      as<double>(coefficients["measurement_error_sd"]);
  params.log_sigma_e = Scalar(std::log(params.measurement_error_sd));
  params.random_effect_sigma =
      as<arma::mat>(coefficients["random_effect_sigma"]);

  load_config(params, parameters);

  // Branch from fixed effects (longitudinal[0]=b1, longitudinal[1]=b2)
  params.branch = classify_disc(longitudinal[0], longitudinal[1]);
}

// Fill marginal ODE parameters from flat theta vector
// Branch must be pre-classified with classify_disc() using double values
template <typename Scalar>
inline void load_marginal_params(
    MarginalParams<Scalar>& params,
    const std::vector<Scalar>& theta,
    MatExpBranch branch,
    double biomarker_clamp = BIOMARKER_CLAMP_DEFAULT) {
  params.parameters = theta;
  params.branch = branch;
  params.biomarker_clamp = biomarker_clamp;
}

// Update branch when random effects shift b1 or b2
inline void update_branch(MatExpBranch& branch,
                           const NumericVector& longitudinal,
                           const NumericVector& re,
                           bool biomarker_random, bool velocity_random) {
  double b1 = longitudinal[0], b2 = longitudinal[1];
  int ri = 0;
  if (biomarker_random) b1 += re[ri++];
  if (velocity_random) b2 += re[ri++];
  branch = classify_disc(b1, b2);
}

// ============================================================================
// Time utilities
// ============================================================================

// Build time points for ODE solving
inline std::vector<double> time_grid(const std::vector<double>& obs_times,
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
inline int time_index(const std::vector<double>& times, double target) {
  auto it = std::lower_bound(times.begin(), times.end(), target);
  return (it != times.end() && std::abs(*it - target) < TIME_MATCH_TOL)
             ? std::distance(times.begin(), it)
             : -1;
}

// Extract covariates at a specific time
inline void covariates_at(
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
inline void bspline_basis(double t, int degree,
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
inline Scalar log_hazard(Scalar biomarker, Scalar velocity,
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
inline Scalar acceleration(
    Scalar biomarker, Scalar velocity, double t,
    const SubjectData<Scalar>& subj,
    const std::vector<Scalar>& longitudinal_coefs, bool biomarker_fixed,
    bool biomarker_random, bool velocity_fixed, bool velocity_random) {
  std::vector<double> covs_fixed, covs_random;
  covariates_at(t, subj.longitudinal_times,
                             subj.longitudinal_covariates_fixed, covs_fixed);
  covariates_at(t, subj.longitudinal_times,
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
// Matrix Exponential ODE Solver
// ============================================================================

// Advance [m, v] by dt. Solves d2m/dt2 = b1*m + b2*m' + f exactly via
// matrix exponential (Cayley-Hamilton) + integral formulation for forcing.
// Branch pre-classified from doubles — no CondExp, minimal AD tape.
template <typename Scalar>
inline void ode_step(Scalar& m, Scalar& v,
                          Scalar b1, Scalar b2, Scalar f,
                          Scalar dt, MatExpBranch br) {
  const Scalar h = b2 * Scalar(0.5);
  Scalar a0, a1, J0, J1;

  if (br == MatExpBranch::REAL) {
    // D > 0, |b1| > eps: both roots well-separated from zero
    const Scalar s = CppAD::sqrt(b2 * b2 + Scalar(4.0) * b1);
    const Scalar l1 = h + s * Scalar(0.5), l2 = h - s * Scalar(0.5);
    const Scalar e1 = safe_exp(l1 * dt), e2 = safe_exp(l2 * dt);
    a0 = (l1 * e2 - l2 * e1) / s;
    a1 = (e1 - e2) / s;
    const Scalar F1 = (e1 - Scalar(1.0)) / l1;
    const Scalar F2 = (e2 - Scalar(1.0)) / l2;
    J1 = (F1 - F2) / s;
    J0 = (l1 * F2 - l2 * F1) / s;

  } else if (br == MatExpBranch::FIRST_ORD) {
    // b1 ~ 0: m'' = b2*m' + f (first-order in v = m')
    // v(t) = (v0 + f/b2)*e^{b2t} - f/b2
    // m(t) = m0 + (v0 + f/b2)*(e^{b2t}-1)/b2 - f*t/b2
    const Scalar eb = safe_exp(b2 * dt);
    const Scalar Fb = (eb - Scalar(1.0)) / b2;  // b2 != 0 (guaranteed)
    const Scalar vf = v + f / b2;  // v0 + f/b2
    m = m + vf * Fb - f * dt / b2;
    v = vf * eb - f / b2;
    return;  // skip the general formula below

  } else if (br == MatExpBranch::COMPLEX) {
    // D < 0: conjugate roots a ± iw, with a²+w² = -b1 > 0
    const Scalar w = CppAD::sqrt(-(b2 * b2 + Scalar(4.0) * b1)) * Scalar(0.5);
    const Scalar r2 = h * h + w * w;
    const Scalar ea = safe_exp(h * dt);
    const Scalar c = CppAD::cos(w * dt), s = CppAD::sin(w * dt);
    a0 = ea * (c - h * s / w);
    a1 = ea * s / w;
    const Scalar Ic = (ea * (h * c + w * s) - h) / r2;
    const Scalar Is = (ea * (h * s - w * c) + w) / r2;
    J1 = Is / w;
    J0 = Ic - (h / w) * Is;

  } else if (br == MatExpBranch::REPEATED) {
    // D ~ 0, |h| > eps: double root at h
    const Scalar e = safe_exp(h * dt);
    a0 = e * (Scalar(1.0) - h * dt);
    a1 = dt * e;
    const Scalar F = (e - Scalar(1.0)) / h;
    J1 = (dt * e - F) / h;
    J0 = Scalar(2.0) * F - dt * e;

  } else {  // ZERO: D ~ 0, h ~ 0
    a0 = Scalar(1.0);  a1 = dt;
    J1 = dt * dt * Scalar(0.5);  J0 = dt;
  }

  const Scalar mn = a0 * m + a1 * v + f * J1;
  const Scalar vn = b1 * a1 * m + (a0 + b2 * a1) * v + f * (J0 + b2 * J1);
  m = mn;  v = vn;
}

// Compute clamped hazard value from log-hazard
template <typename Scalar>
inline Scalar eval_hazard(Scalar biomarker, Scalar velocity,
                          const ODEParams<Scalar>& params,
                          BSplineWorkspace& ws, double t) {
  bspline_basis(t, params.spline_degree, params.spline_knots,
                params.spline_boundary, ws.basis, ws.knots,
                ws.work1, ws.work2, true);
  Scalar lh = log_hazard(biomarker, velocity, ws.basis,
                          params.baseline_coefs, params.hazard_coefs,
                          params.subject.survival_covariates, params.gamma);
  return clamp(lh, Scalar(HAZARD_CLAMP_MIN), Scalar(HAZARD_CLAMP_MAX));
}

// Solve joint ODE: matexp for biomarker, composite Simpson's rule for hazard
template <typename Scalar>
inline std::vector<std::vector<Scalar>> ode_solve_joint(
    const std::vector<Scalar>& y0, const std::vector<double>& times,
    const ODEParams<Scalar>& params) {
  const int nt = times.size();
  const int N = std::max(1, params.hazard_quadrature);
  BSplineWorkspace ws;
  bspline_basis(0.0, params.spline_degree, params.spline_knots,
                params.spline_boundary, ws.basis, ws.knots,
                ws.work1, ws.work2, false);

  std::vector<std::vector<Scalar>> sol(nt, std::vector<Scalar>(3));
  Scalar H = y0[0], m = y0[1], v = y0[2];
  const Scalar bc(params.biomarker_clamp);
  sol[0] = {H, m, v};

  for (int i = 1; i < nt; i++) {
    const double t0 = times[i - 1], t1 = times[i], dt = t1 - t0;

    // Extract beta1, beta2, forcing (piecewise constant at t0)
    Scalar b1(0), b2(0), f(0);
    std::vector<double> cf, cr;
    covariates_at(t0, params.subject.longitudinal_times,
                  params.subject.longitudinal_covariates_fixed, cf);
    covariates_at(t0, params.subject.longitudinal_times,
                  params.subject.longitudinal_covariates_random, cr);
    size_t fi = 0, ri = 0;
    if (params.biomarker_fixed) b1 = params.longitudinal_coefs[fi++];
    if (params.biomarker_random) b1 += params.subject.random_effect[ri++];
    if (params.velocity_fixed) b2 = params.longitudinal_coefs[fi++];
    if (params.velocity_random) b2 += params.subject.random_effect[ri++];
    for (size_t j = 0; j < cf.size(); j++)
      f += params.longitudinal_coefs[fi++] * Scalar(cf[j]);
    for (size_t j = 0; j < cr.size(); j++)
      f += params.subject.random_effect[ri++] * Scalar(cr[j]);

    // Composite Simpson's rule for cumulative hazard over N sub-intervals
    const double sub_dt = dt / N;
    Scalar h_left = safe_exp(eval_hazard(m, v, params, ws, t0));

    for (int k = 0; k < N; k++) {
      const double ta = t0 + k * sub_dt;

      // Advance to midpoint, evaluate hazard
      Scalar m_sub = m, v_sub = v;
      ode_step(m_sub, v_sub, b1, b2, f, Scalar(sub_dt * 0.5), params.branch);
      m_sub = clamp(m_sub, -bc, bc);
      v_sub = clamp(v_sub, -bc, bc);
      Scalar h_mid = safe_exp(
          eval_hazard(m_sub, v_sub, params, ws, ta + sub_dt * 0.5));

      // Advance midpoint to right boundary, evaluate hazard
      ode_step(m_sub, v_sub, b1, b2, f, Scalar(sub_dt * 0.5), params.branch);
      m_sub = clamp(m_sub, -bc, bc);
      v_sub = clamp(v_sub, -bc, bc);
      Scalar h_right = safe_exp(
          eval_hazard(m_sub, v_sub, params, ws, ta + sub_dt));

      H += Scalar(sub_dt / 6.0) *
           (h_left + Scalar(4.0) * h_mid + h_right);

      m = m_sub;
      v = v_sub;
      h_left = h_right;
    }

    sol[i] = {H, m, v};
  }
  return sol;
}

// ============================================================================
// Solve marginal ODE: [m, v] only, no hazard integration
// ============================================================================

// Mirrors ode_solve_joint(y0, times, params)
template <typename Scalar>
inline std::vector<std::vector<Scalar>> ode_solve_marginal(
    const std::vector<Scalar>& y0,
    const std::vector<double>& times,
    const MarginalParams<Scalar>& params) {
  const Scalar b1 = params.parameters[0], b2 = params.parameters[1];
  const double bc = params.biomarker_clamp;
  const int nt = times.size();

  std::vector<std::vector<Scalar>> sol(nt, std::vector<Scalar>(2));
  Scalar m = y0[0], v = y0[1];
  sol[0] = {m, v};

  for (int i = 1; i < nt; i++) {
    const double dt = times[i] - times[i - 1];

    std::vector<double> cf;
    covariates_at(times[i - 1], params.subject.longitudinal_times,
                  params.subject.longitudinal_covariates_fixed, cf);
    Scalar f(0);
    for (size_t j = 0; j < cf.size(); j++)
      f += params.parameters[j + 2] * Scalar(cf[j]);

    ode_step(m, v, b1, b2, f, Scalar(dt), params.branch);
    m = clamp(m, Scalar(-bc), Scalar(bc));
    v = clamp(v, Scalar(-bc), Scalar(bc));
    sol[i] = {m, v};
  }
  return sol;
}

// ============================================================================
// Shared objective functions
// ============================================================================

// Marginal log-likelihood: longitudinal only (-SSE)
// Mirrors joint_loglik(sol, times, params, ...)
template <typename Scalar>
inline Scalar marginal_loglik(
    const std::vector<std::vector<Scalar>>& sol,
    const MarginalParams<Scalar>& params) {
  Scalar ll(0);
  const int n_obs = params.subject.longitudinal_measurements.size();
  for (int i = 0; i < n_obs; i++) {
    Scalar r = Scalar(params.subject.longitudinal_measurements[i]) - sol[i][0];
    ll -= r * r;
  }
  return ll;
}

// Joint log-likelihood: longitudinal (Gaussian) + survival (Cox)
// Used by objective.cpp, posterior.cpp, state.cpp
template <typename Scalar>
inline Scalar joint_loglik(
    const std::vector<std::vector<Scalar>>& sol,
    const std::vector<double>& times,
    const ODEParams<Scalar>& params,
    BSplineWorkspace& ws) {
  Scalar ll(0);
  const int n_obs = params.subject.longitudinal_times.size();
  const Scalar inv_sigma_e2 = CppAD::exp(Scalar(-2) * params.log_sigma_e);
  const Scalar half_inv = Scalar(0.5) * inv_sigma_e2;

  // Longitudinal: -0.5 * sum((y - m)^2 / sigma_e^2)
  for (int i = 0; i < n_obs; i++) {
    const int idx = time_index(times, params.subject.longitudinal_times[i]);
    if (idx >= 0) {
      Scalar r = Scalar(params.subject.longitudinal_measurements[i]) - sol[idx][1];
      ll -= half_inv * r * r;
    }
  }
  if (n_obs > 0)
    ll -= Scalar(0.5) * Scalar(n_obs) *
          (Scalar(std::log(2.0 * M_PI)) + Scalar(2) * params.log_sigma_e);

  // Survival: -H(T) + delta * log h(T)
  const int ei = time_index(times, params.subject.event_time);
  if (ei >= 0) {
    ll -= sol[ei][0];
    if (params.subject.status == 1) {
      bspline_basis(params.subject.event_time, params.spline_degree,
                            params.spline_knots, params.spline_boundary,
                            ws.basis, ws.knots, ws.work1, ws.work2, true);
      ll += log_hazard(sol[ei][1], sol[ei][2], ws.basis,
                                params.baseline_coefs, params.hazard_coefs,
                                params.subject.survival_covariates, params.gamma);
    }
  }
  return ll;
}

// Evaluate CppAD tape: optimize, forward, reverse, hessian, cleanup
inline NumericVector eval_tape(CppAD::ADFun<double>& tape,
                                    const std::vector<double>& x,
                                    int n, bool grad, bool hess) {
  tape.optimize();
  tape.check_for_nan(false);
  std::vector<double> y = tape.Forward(0, x);

  NumericVector result(1);
  result[0] = y[0];

  if (grad || hess) {
    std::vector<double> g = tape.Reverse(1, std::vector<double>(1, 1.0));
    if (grad) result.attr("gradient") = wrap(g);
    if (hess) {
      std::vector<double> h = tape.Hessian(x, 0);
      NumericMatrix hm(n, n);
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          hm(i, j) = h[i * n + j];
      result.attr("hessian") = hm;
    }
  }
  tape.capacity_order(0);
  return result;
}

// Compute log-posterior for a single subject (shared by AD and non-AD paths)
template <typename Scalar>
inline Scalar eval_logpost(
    const std::vector<Scalar>& re,
    ODEParams<Scalar>& ode,
    const arma::mat& inv_sigma_b, double re_const) {
  const auto times = time_grid(ode.subject.longitudinal_times,
                                  ode.subject.event_time);
  const std::vector<Scalar> y0 = {Scalar(0.0),
      ode.initial_state_coefs[0] + re[0],
      ode.initial_state_coefs[1] + re[1]};
  const auto sol = ode_solve_joint(y0, times, ode);

  BSplineWorkspace ws;
  bspline_basis(0.0, ode.spline_degree, ode.spline_knots,
      ode.spline_boundary, ws.basis, ws.knots, ws.work1, ws.work2, false);
  Scalar lp = joint_loglik(sol, times, ode, ws);

  const int n_re = re.size();
  Scalar qf(0.0);
  for (int i = 0; i < n_re; i++)
    for (int j = 0; j < n_re; j++)
      qf += re[i] * Scalar(inv_sigma_b(i, j)) * re[j];
  lp -= Scalar(0.5) * qf + Scalar(re_const);
  return lp;
}

#endif  // SOLVER_H
