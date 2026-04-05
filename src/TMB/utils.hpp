// utils.hpp — Shared utilities for JointODE TMB templates
// ODE solver, B-spline basis, hazard evaluation, time grid helpers
#ifndef JOINTODE_UTILS_HPP
#define JOINTODE_UTILS_HPP

// ODE branch classification (BR_ prefix avoids R macro REAL/COMPLEX clash)
enum {
  BR_REAL      = 0,  // D > 0, |b1| > eps: two distinct real roots
  BR_FIRST_ORD = 1,  // |b1| ~ 0: degenerates to first-order ODE
  BR_COMPLEX   = 2,  // D < 0: conjugate complex roots
  BR_REPEATED  = 3,  // D ~ 0, |b2| > eps: double root
  BR_ZERO      = 4   // D ~ 0, |b2| ~ 0: trivial case
};

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

// ============================================================================
// Matrix exponential ODE solver
// Solves d²m/dt² = b1·m + b2·dm/dt + forcing exactly via Cayley-Hamilton
// Branch pre-classified from doubles → no CondExp overhead in AD tape
// ============================================================================

template <class Type>
void ode_step(Type& biomarker, Type& velocity,
              Type b1, Type b2, Type forcing, Type dt, int branch) {
  Type half_b2 = b2 * Type(0.5);
  Type a0, a1, J0, J1;

  if (branch == BR_REAL) {
    Type discriminant = sqrt(b2 * b2 + Type(4) * b1);
    Type root1 = half_b2 + discriminant * Type(0.5);
    Type root2 = half_b2 - discriminant * Type(0.5);
    Type exp1 = safe_exp(root1 * dt), exp2 = safe_exp(root2 * dt);
    a0 = (root1 * exp2 - root2 * exp1) / discriminant;
    a1 = (exp1 - exp2) / discriminant;
    Type int1 = (exp1 - Type(1)) / root1;
    Type int2 = (exp2 - Type(1)) / root2;
    J1 = (int1 - int2) / discriminant;
    J0 = (root1 * int2 - root2 * int1) / discriminant;

  } else if (branch == BR_FIRST_ORD) {
    Type exp_b2 = safe_exp(b2 * dt);
    Type int_b2 = (exp_b2 - Type(1)) / b2;
    Type shifted_vel = velocity + forcing / b2;
    biomarker += shifted_vel * int_b2 - forcing * dt / b2;
    velocity = shifted_vel * exp_b2 - forcing / b2;
    return;

  } else if (branch == BR_COMPLEX) {
    Type omega = sqrt(-(b2 * b2 + Type(4) * b1)) * Type(0.5);
    Type mod_sq = half_b2 * half_b2 + omega * omega;
    Type exp_real = safe_exp(half_b2 * dt);
    Type cos_w = cos(omega * dt), sin_w = sin(omega * dt);
    a0 = exp_real * (cos_w - half_b2 * sin_w / omega);
    a1 = exp_real * sin_w / omega;
    Type int_cos = (exp_real * (half_b2 * cos_w + omega * sin_w) - half_b2) / mod_sq;
    Type int_sin = (exp_real * (half_b2 * sin_w - omega * cos_w) + omega) / mod_sq;
    J1 = int_sin / omega;
    J0 = int_cos - (half_b2 / omega) * int_sin;

  } else if (branch == BR_REPEATED) {
    Type exp_h = safe_exp(half_b2 * dt);
    Type int_h = (exp_h - Type(1)) / half_b2;
    a0 = exp_h * (Type(1) - half_b2 * dt);
    a1 = dt * exp_h;
    J1 = (dt * exp_h - int_h) / half_b2;
    J0 = Type(2) * int_h - dt * exp_h;

  } else { // BR_ZERO
    a0 = Type(1);  a1 = dt;
    J1 = dt * dt * Type(0.5);  J0 = dt;
  }

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
  int lo = std::max(0, span - degree), hi = std::min(n_basis - 1, span);
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
  for (int k = 0; k < baseline_coefs.size(); k++)
    log_h += Type(basis(k)) * baseline_coefs(k);
  log_h += hazard_coefs(0) * biomarker;
  if (gamma_power == 1)      log_h += hazard_coefs(1) * velocity;
  else if (gamma_power == 2) log_h += hazard_coefs(1) * velocity * velocity;
  for (int k = 0; k < survival_covariates.size(); k++)
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
  for (int i = 0; i < obs_times.size(); i++)
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
