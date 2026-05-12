// utils.hpp — Shared utilities for JointODE TMB templates
// ODE solver, B-spline basis, hazard evaluation, time grid helpers
#ifndef JOINTODE_UTILS_HPP
#define JOINTODE_UTILS_HPP

// ============================================================================
// Smooth special functions for Cayley-Hamilton ODE solver.
// Both branches of every CondExp produce finite values and gradients
// for all inputs, so AD tape is never corrupted.
//
// Guard pattern: in the "direct" branch (large |x|), the denominator is
// replaced by x_safe = CondExp(|x|<eps, 1, x).  When the Taylor branch
// is selected, the direct branch still evaluates but divides by 1 instead
// of 0, preventing 0/0 gradient leakage through the AD tape.
// ============================================================================

// sinhc(x) = sinh(x)/x, smooth at x=0 (limit = 1)
// Direct branch guards denominator to avoid 0/0 AD gradient leak.
template <class Type>
Type sinhc(Type x) {
  Type x2 = x * x;
  Type x_safe = CppAD::CondExpLt(x2, Type(1e-8), Type(1), x);
  return CppAD::CondExpLt(x2, Type(1e-8),
    Type(1) + x2 / Type(6) + x2 * x2 / Type(120),
    sinh(x) / x_safe);
}

// sinc(x) = sin(x)/x, smooth at x=0 (limit = 1)
template <class Type>
Type sinc(Type x) {
  Type x2 = x * x;
  Type x_safe = CppAD::CondExpLt(x2, Type(1e-8), Type(1), x);
  return CppAD::CondExpLt(x2, Type(1e-8),
    Type(1) - x2 / Type(6) + x2 * x2 / Type(120),
    sin(x) / x_safe);
}

// expm1c2(x) = (exp(x) - 1 - x) / x^2, smooth at x=0 (limit = 1/2)
// Direct branch guards denominator to avoid 0/0 AD gradient leak.
template <class Type>
Type expm1c2(Type x) {
  Type x2 = x * x;
  Type x2_safe = CppAD::CondExpLt(x2, Type(1e-8), Type(1), x2);
  return CppAD::CondExpLt(x2, Type(1e-8),
    Type(0.5) + x / Type(6) + x2 / Type(24),
    (exp(x) - Type(1) - x) / x2_safe);
}

// ============================================================================
// Cayley-Hamilton ODE solver (sinhc/sinc formulation)
// Solves d²m/dt² = b1·m + b2·dm/dt + forcing exactly.
//
// Transition coefficients expressed as:
//   a0 = exp(m) * (C - m·Sc)     [state propagation]
//   a1 = exp(m) * dt * Sc         [velocity-to-state coupling]
// where m = b2·dt/2, and C/Sc = cosh/sinhc (D>=0) or cos/sinc (D<0).
//
// Near D=0, both branches share the same Taylor series in s2 = D·dt²/4
// (because sinh(ix)/ix = sin(x)/x), so a single Taylor expansion
// covers both cases with no branch needed.
//
// Three-way CondExp: Taylor (|s2|<eps) > real (s2>0) > complex (s2<0).
// The Taylor branch bypasses sqrt entirely, avoiding d(sqrt)/dx = Inf
// at s2=0.  In the trig/hyp branches, sqrt argument is clamped >= eps.
// ============================================================================

template <class Type>
Type observation_logdensity(Type y, Type mean, Type sigma,
                            int residual_family, Type nu) {
  if (residual_family == 0) return dnorm(y, mean, sigma, true);
  Type r = (y - mean) / sigma;
  return lgamma((nu + Type(1)) / Type(2)) -
         lgamma(nu / Type(2)) -
         Type(0.5) * (log(nu) + log(Type(3.14159265358979323846))) -
         log(sigma) -
         ((nu + Type(1)) / Type(2)) * log(Type(1) + r * r / nu);
}

template <class Type>
void ode_step(Type& biomarker, Type& velocity,
              Type b1, Type b2, Type forcing, Type dt) {
  Type m = b2 * dt * Type(0.5);
  Type D = b2 * b2 + Type(4) * b1;
  Type s2 = D * dt * dt * Type(0.25);  // delta^2 or -omega^2

  // Taylor expansion in s2 (valid for both D>0 and D<0 near D=0)
  Type s4 = s2 * s2;
  Type sc_taylor = Type(1) + s2 / Type(6) + s4 / Type(120);
  Type cc_taylor = Type(1) + s2 / Type(2) + s4 / Type(24);

  // For |s2| >= eps, compute via sqrt + trig/hyp
  Type eps_s2(1e-6);
  Type abs_s2 = CppAD::CondExpGt(s2, Type(0), s2, -s2);
  Type safe_abs_s2 = CppAD::CondExpGt(abs_s2, eps_s2, abs_s2, eps_s2);
  Type r = sqrt(safe_abs_s2);  // sqrt is safe since safe_abs_s2 >= eps

  // Real branch: sinhc(r), cosh(r)
  Type sc_real = sinhc(r);
  Type cc_real = cosh(r);
  // Complex branch: sinc(r), cos(r) — note r = sqrt(-s2) = omega
  Type sc_comp = sinc(r);
  Type cc_comp = cos(r);

  // Three-way selection: Taylor (|s2|<eps) > real (s2>0) > complex (s2<0)
  Type sc_rc = CppAD::CondExpGt(s2, Type(0), sc_real, sc_comp);
  Type cc_rc = CppAD::CondExpGt(s2, Type(0), cc_real, cc_comp);
  Type sc = CppAD::CondExpLt(abs_s2, eps_s2, sc_taylor, sc_rc);
  Type cc = CppAD::CondExpLt(abs_s2, eps_s2, cc_taylor, cc_rc);

  Type em = exp(m);
  Type a0 = em * (cc - m * sc);
  Type a1 = em * dt * sc;

  // Forcing integral for m: J1 = (a0 - 1) / b1
  // Removable singularity at b1=0: limit = dt^2 * expm1c2(b2*dt)
  Type b1_safe = CppAD::CondExpGt(b1 * b1, Type(1e-20), b1, Type(1));
  Type J1_direct = (a0 - Type(1)) / b1_safe;
  Type J1_taylor = dt * dt * expm1c2(b2 * dt);
  Type J1 = CppAD::CondExpGt(b1 * b1, Type(1e-20), J1_direct, J1_taylor);

  // Update: v_forcing = f*a1 (exact, from augmented 3x3 Cayley-Hamilton)
  Type new_m = a0 * biomarker + a1 * velocity + forcing * J1;
  Type new_v = b1 * a1 * biomarker + (a0 + b2 * a1) * velocity +
               forcing * a1;
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
// Returns exp(log-hazard) for Simpson integration
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
  return exp(log_h);
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
