// utils.hpp — Shared utilities for JointODE TMB templates
// ODE solver, B-spline basis, hazard evaluation, time grid helpers
#ifndef JOINTODE_UTILS_HPP
#define JOINTODE_UTILS_HPP

// ============================================================================
// Core math utilities
// ============================================================================

// Safe exp with linear extension beyond threshold to prevent overflow
template <class Type>
Type safe_exp(Type x) {
  Type threshold(500.0);
  Type capped = CppAD::CondExpGt(x, threshold, threshold, x);
  Type excess  = CppAD::CondExpGt(x, threshold, x - threshold, Type(0));
  return exp(capped) * (Type(1) + excess);
}

// AD-safe clamping
template <class Type>
Type clamp(Type x, Type lower, Type upper) {
  return CppAD::CondExpLt(x, lower, lower,
           CppAD::CondExpGt(x, upper, upper, x));
}

// AD-safe sqrt: returns sqrt(max(x, eps)) to avoid NaN on negative input
template <class Type>
Type safe_sqrt(Type x) {
  Type eps(1e-20);
  return sqrt(CppAD::CondExpGt(x, eps, x, eps));
}

// ============================================================================
// Matrix exponential ODE solver
// Solves d²m/dt² = b1·m + b2·dm/dt + forcing exactly via Cayley-Hamilton
// Uses CondExp to select between REAL (D>0) and COMPLEX (D<0) branches,
// adapting to current parameter values during optimization.
// ============================================================================

template <class Type>
void ode_step(Type& biomarker, Type& velocity,
              Type b1, Type b2, Type forcing, Type dt) {
  Type half_b2 = b2 * Type(0.5);
  Type D = b2 * b2 + Type(4) * b1;  // discriminant

  // --- REAL branch (D > 0): two distinct real roots ---
  Type D_real = safe_sqrt(D);  // safe when D <= 0
  Type root1 = half_b2 + D_real * Type(0.5);
  Type root2 = half_b2 - D_real * Type(0.5);
  // Guard roots near zero to avoid division by zero
  Type root1_safe = CppAD::CondExpGt(root1 * root1, Type(1e-20), root1, Type(1e-10));
  Type root2_safe = CppAD::CondExpGt(root2 * root2, Type(1e-20), root2, Type(-1e-10));
  Type D_real_safe = CppAD::CondExpGt(D_real * D_real, Type(1e-20), D_real, Type(1e-10));
  Type exp1 = safe_exp(root1 * dt), exp2 = safe_exp(root2 * dt);
  Type a0_r = (root1 * exp2 - root2 * exp1) / D_real_safe;
  Type a1_r = (exp1 - exp2) / D_real_safe;
  Type int1 = (exp1 - Type(1)) / root1_safe;
  Type int2 = (exp2 - Type(1)) / root2_safe;
  Type J1_r = (int1 - int2) / D_real_safe;
  Type J0_r = (root1 * int2 - root2 * int1) / D_real_safe;

  // --- COMPLEX branch (D < 0): conjugate complex roots ---
  Type omega = safe_sqrt(-D) * Type(0.5);  // safe when D >= 0
  Type omega_safe = CppAD::CondExpGt(omega * omega, Type(1e-20), omega, Type(1e-10));
  Type mod_sq = half_b2 * half_b2 + omega * omega;
  Type mod_sq_safe = CppAD::CondExpGt(mod_sq, Type(1e-20), mod_sq, Type(1e-10));
  Type exp_real = safe_exp(half_b2 * dt);
  Type cos_w = cos(omega * dt), sin_w = sin(omega * dt);
  Type a0_c = exp_real * (cos_w - half_b2 * sin_w / omega_safe);
  Type a1_c = exp_real * sin_w / omega_safe;
  Type int_cos = (exp_real * (half_b2 * cos_w + omega * sin_w) - half_b2) / mod_sq_safe;
  Type int_sin = (exp_real * (half_b2 * sin_w - omega * cos_w) + omega) / mod_sq_safe;
  Type J1_c = int_sin / omega_safe;
  Type J0_c = int_cos - (half_b2 / omega_safe) * int_sin;

  // --- Select branch via CondExp based on discriminant sign ---
  Type zero(0);
  Type a0 = CppAD::CondExpGt(D, zero, a0_r, a0_c);
  Type a1 = CppAD::CondExpGt(D, zero, a1_r, a1_c);
  Type J0 = CppAD::CondExpGt(D, zero, J0_r, J0_c);
  Type J1 = CppAD::CondExpGt(D, zero, J1_r, J1_c);

  Type new_m = a0 * biomarker + a1 * velocity + forcing * J1;
  Type new_v = b1 * a1 * biomarker + (a0 + b2 * a1) * velocity +
               forcing * (J0 + b2 * J1);
  biomarker = new_m;
  velocity = new_v;
}

// ============================================================================
// B-spline basis evaluation (double only — basis values are data, not AD)
// De Boor algorithm with repeated boundary knots
// ============================================================================

// Build full knot vector once (called before subject loop)
vector<double> build_knot_vector(int degree,
                                 const vector<double>& interior_knots,
                                 const vector<double>& boundary_knots) {
  int n_interior = interior_knots.size();
  int n_knots = n_interior + 2 * (degree + 1);
  vector<double> knots(n_knots);
  for (int i = 0; i <= degree; i++) knots(i) = boundary_knots(0);
  for (int i = 0; i < n_interior; i++) knots(degree + 1 + i) = interior_knots(i);
  for (int i = 0; i <= degree; i++) knots(degree + 1 + n_interior + i) = boundary_knots(1);
  return knots;
}

// Evaluate B-spline basis at t using pre-built knot vector
vector<double> bspline_basis(double t, int degree,
                             const vector<double>& knots) {
  int n_knots = knots.size();
  int n_basis = n_knots - degree - 1;

  t = std::max(knots(0), std::min(knots(n_knots - 1), t));

  int span = degree;
  if (t >= knots(n_knots - degree - 1)) span = n_knots - degree - 2;
  else for (int i = degree; i < n_knots - degree - 1; ++i)
    if (t < knots(i + 1)) { span = i; break; }

  // De Boor recursion — only iterate over [span-degree, span] range
  vector<double> current(n_basis), previous(n_basis);
  current.setZero();
  if (span < n_basis) current(span) = 1.0;

  for (int p = 1; p <= degree; p++) {
    previous = current;
    current.setZero();
    int lo_p = std::max(0, span - p), hi_p = std::min(n_basis - 1, span);
    for (int i = lo_p; i <= hi_p; i++) {
      if (i < n_knots - p) {
        double d = knots(i + p) - knots(i);
        if (d > 1e-8) current(i) += previous(i) * (t - knots(i)) / d;
      }
      if (i + 1 < n_knots - p) {
        double d = knots(i + p + 1) - knots(i + 1);
        if (d > 1e-8) current(i) += previous(i + 1) * (knots(i + p + 1) - t) / d;
      }
    }
  }
  return current;
}

// ============================================================================
// Log-linear hazard evaluation
// log h(t) = B(t)·baseline + alpha_m·m + alpha_v·g(v) + W·beta_surv
// Returns exp(clamped log-hazard) for Simpson integration
// ============================================================================

template <class Type>
Type eval_hazard(Type biomarker, Type velocity, double t,
                 int spline_degree,
                 const vector<double>& knot_vector,
                 const vector<Type>& baseline_coefs,
                 const vector<Type>& hazard_coefs,
                 const vector<double>& survival_covariates,
                 double gamma_power) {
  vector<double> basis = bspline_basis(t, spline_degree, knot_vector);
  Type log_h(0);
  for (size_t k = 0; k < (size_t)baseline_coefs.size(); k++)
    log_h += Type(basis(k)) * baseline_coefs(k);
  log_h += hazard_coefs(0) * biomarker;
  if (gamma_power == 1)      log_h += hazard_coefs(1) * velocity;
  else if (gamma_power == 2) log_h += hazard_coefs(1) * velocity * velocity;
  for (size_t k = 0; k < (size_t)survival_covariates.size(); k++)
    log_h += hazard_coefs(k + 2) * Type(survival_covariates(k));
  return safe_exp(clamp(log_h, Type(-20.0), Type(20.0)));
}

// ============================================================================
// Time grid utilities
// ============================================================================

// Build sorted unique time grid: {0, obs_times ∩ [0, event_time], event_time}
vector<double> build_time_grid(const vector<double>& obs_times,
                               double event_time) {
  std::vector<double> grid;
  grid.reserve(obs_times.size() + 2);
  grid.push_back(0.0);
  for (size_t i = 0; i < (size_t)obs_times.size(); i++)
    grid.push_back(std::max(0.0, std::min(event_time, obs_times(i))));
  grid.push_back(event_time);
  std::sort(grid.begin(), grid.end());
  grid.erase(std::unique(grid.begin(), grid.end()), grid.end());
  vector<double> result(grid.size());
  for (size_t i = 0; i < grid.size(); i++) result(i) = grid[i];
  return result;
}

// Find index of target time in sorted grid (binary search, exact match)
int find_time_index(const vector<double>& times, double target) {
  int lo = 0, hi = times.size() - 1;
  while (lo <= hi) {
    int mid = (lo + hi) / 2;
    if (std::abs(times(mid) - target) < 1e-10) return mid;
    if (times(mid) < target) lo = mid + 1; else hi = mid - 1;
  }
  return -1;
}

// LOCF covariate index: find last observation time <= t
int covariate_index_at(double t, const vector<double>& obs_times) {
  int n = obs_times.size();
  if (n == 0) return 0;
  for (int i = 0; i < n - 1; i++)
    if (t < obs_times(i + 1)) return i;
  return n - 1;
}

#endif
