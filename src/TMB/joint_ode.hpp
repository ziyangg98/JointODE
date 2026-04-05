// joint_ode.hpp — TMB template for JointODE
// Joint model: longitudinal ODE + Cox PH survival with shared random effects
#ifndef JOINT_ODE_HPP
#define JOINT_ODE_HPP

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

vector<double> bspline_basis(double t, int degree,
                             const vector<double>& interior_knots,
                             const vector<double>& boundary_knots) {
  int n_interior = interior_knots.size();
  int n_knots = n_interior + 2 * (degree + 1);
  int n_basis = n_knots - degree - 1;

  // Build full knot vector with repeated boundaries
  vector<double> knots(n_knots);
  for (int i = 0; i <= degree; i++) knots(i) = boundary_knots(0);
  for (int i = 0; i < n_interior; i++) knots(degree + 1 + i) = interior_knots(i);
  for (int i = 0; i <= degree; i++) knots(degree + 1 + n_interior + i) = boundary_knots(1);

  t = std::max(boundary_knots(0), std::min(boundary_knots(1), t));

  // Find knot span
  int span = degree;
  if (t >= knots(n_knots - degree - 1)) span = n_knots - degree - 2;
  else for (int i = degree; i < n_knots - degree - 1; ++i)
    if (t < knots(i + 1)) { span = i; break; }

  // De Boor recursion
  vector<double> current(n_basis), previous(n_basis);
  current.setZero();
  if (span < n_basis) current(span) = 1.0;

  for (int p = 1; p <= degree; p++) {
    previous = current;
    current.setZero();
    for (int i = 0; i < n_basis; i++) {
      if (i < n_knots - p) {
        double denom = knots(i + p) - knots(i);
        if (denom > 1e-8) current(i) += previous(i) * (t - knots(i)) / denom;
      }
      if (i + 1 < n_knots - p) {
        double denom = knots(i + p + 1) - knots(i + 1);
        if (denom > 1e-8) current(i) += previous(i + 1) * (knots(i + p + 1) - t) / denom;
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
                 const vector<double>& spline_interior_knots,
                 const vector<double>& spline_boundary_knots,
                 const vector<Type>& baseline_coefs,
                 const vector<Type>& hazard_coefs,
                 const vector<double>& survival_covariates,
                 double gamma_power) {
  vector<double> basis = bspline_basis(t, spline_degree,
                                        spline_interior_knots,
                                        spline_boundary_knots);
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

// ============================================================================
// TMB objective function (supports OpenMP via parallel_accumulator)
// ============================================================================

template<class Type>
Type objective_function<Type>::operator()() {

  // --- DATA ---
  DATA_INTEGER(n_subjects);
  DATA_INTEGER(n_random_effects);
  DATA_IVECTOR(n_observations);        // per subject
  DATA_VECTOR(event_times);
  DATA_IVECTOR(event_status);
  DATA_VECTOR(obs_times_all);          // concatenated
  DATA_VECTOR(obs_values_all);         // concatenated
  DATA_IMATRIX(covariate_dims);        // n_subjects x 2: [n_fixed, n_random]
  DATA_VECTOR(X_fixed_all);            // row-major flat
  DATA_VECTOR(X_random_all);
  DATA_VECTOR(W_survival_all);
  DATA_INTEGER(n_survival_covariates);
  DATA_IVECTOR(ode_branch);
  DATA_SCALAR(biomarker_clamp);
  DATA_INTEGER(hazard_quadrature);
  DATA_SCALAR(gamma);
  DATA_INTEGER(spline_degree);
  DATA_VECTOR(spline_knots);
  DATA_VECTOR(spline_boundary);
  DATA_INTEGER(biomarker_fixed);
  DATA_INTEGER(biomarker_random);
  DATA_INTEGER(velocity_fixed);
  DATA_INTEGER(velocity_random);

  // Convert DATA_SCALAR/VECTOR to double for non-AD helper functions
  double gamma_val = asDouble(gamma);
  double clamp_val = asDouble(biomarker_clamp);
  int n_quadrature = std::max(hazard_quadrature, 1);
  vector<double> knots_d(spline_knots.size()), boundary_d(spline_boundary.size());
  for (int i = 0; i < knots_d.size(); i++) knots_d(i) = asDouble(spline_knots(i));
  for (int i = 0; i < boundary_d.size(); i++) boundary_d(i) = asDouble(spline_boundary(i));

  // Pre-compute per-subject observation offsets for parallel access
  vector<int> obs_offset(n_subjects), fixed_offset(n_subjects);
  vector<int> random_offset(n_subjects), surv_offset(n_subjects);
  obs_offset(0) = fixed_offset(0) = random_offset(0) = surv_offset(0) = 0;
  for (int i = 1; i < n_subjects; i++) {
    int ni = n_observations(i - 1);
    int nf = covariate_dims(i - 1, 0), nr = covariate_dims(i - 1, 1);
    obs_offset(i) = obs_offset(i - 1) + ni;
    fixed_offset(i) = fixed_offset(i - 1) + ni * nf;
    random_offset(i) = random_offset(i - 1) + ni * nr;
    surv_offset(i) = surv_offset(i - 1) + n_survival_covariates;
  }

  // --- PARAMETERS ---
  PARAMETER_VECTOR(baseline);
  PARAMETER_VECTOR(hazard);
  PARAMETER_VECTOR(longitudinal);
  PARAMETER_VECTOR(initial_state);
  PARAMETER(log_sigma_e);
  PARAMETER_VECTOR(log_sd_re);
  PARAMETER_VECTOR(corr_par);
  PARAMETER_MATRIX(b);  // n_subjects x n_re (random effects)

  // Covariance structure: VECSCALE(UNSTRUCTURED_CORR, exp(log_sd))
  using namespace density;
  UNSTRUCTURED_CORR_t<Type> corr_structure(corr_par);
  Type sigma_e = exp(log_sigma_e);
  Type BC = Type(clamp_val);

  // Output vectors for REPORT (filled per-subject, safe for parallel)
  int n_total_obs = obs_values_all.size();
  vector<Type> fitted_biomarker(n_total_obs), fitted_velocity(n_total_obs);
  vector<Type> cumulative_hazard(n_subjects), log_hazard_at_event(n_subjects);

  // Parallel accumulator for negative log-likelihood
  parallel_accumulator<Type> nll(this);

  for (int i = 0; i < n_subjects; i++) {
    int ni = n_observations(i);
    int n_fixed_cov = covariate_dims(i, 0);
    int n_random_cov = covariate_dims(i, 1);
    int oi = obs_offset(i), fi = fixed_offset(i);
    int ri = random_offset(i), si = surv_offset(i);

    vector<Type> bi = b.row(i);

    // --- Random effect prior: MVN(0, Sigma_b) ---
    nll += VECSCALE(corr_structure, exp(log_sd_re))(bi);

    // --- Unpack subject data as double ---
    vector<double> obs_t(ni), obs_y(ni);
    for (int j = 0; j < ni; j++) {
      obs_t(j) = asDouble(obs_times_all(oi + j));
      obs_y(j) = asDouble(obs_values_all(oi + j));
    }
    matrix<double> X_fixed(ni, n_fixed_cov);
    for (int j = 0; j < ni * n_fixed_cov; j++)
      X_fixed(j / n_fixed_cov, j % n_fixed_cov) = asDouble(X_fixed_all(fi + j));
    matrix<double> X_random(ni, n_random_cov);
    for (int j = 0; j < ni * n_random_cov; j++)
      X_random(j / n_random_cov, j % n_random_cov) = asDouble(X_random_all(ri + j));
    vector<double> surv_cov(n_survival_covariates);
    for (int k = 0; k < n_survival_covariates; k++)
      surv_cov(k) = asDouble(W_survival_all(si + k));

    double event_time = asDouble(event_times(i));
    int status = event_status(i), branch = ode_branch(i);

    // --- ODE solve with Simpson hazard integration ---
    vector<double> time_pts = build_time_grid(obs_t, event_time);
    int n_times = time_pts.size();
    Type m = initial_state(0) + bi(0);
    Type v = initial_state(1) + bi(1);
    Type cum_haz(0);
    vector<Type> sol_m(n_times), sol_v(n_times), sol_H(n_times);
    sol_m(0) = m;  sol_v(0) = v;  sol_H(0) = cum_haz;

    // b1, b2 are constant across time steps (depend on fixed + RE only)
    Type b1(0), b2(0);
    int fixed_idx = 0, re_idx = 2;
    if (biomarker_fixed)  b1  = longitudinal(fixed_idx++);
    if (biomarker_random) b1 += bi(re_idx++);
    if (velocity_fixed)   b2  = longitudinal(fixed_idx++);
    if (velocity_random)  b2 += bi(re_idx++);
    int forcing_fixed_start = fixed_idx;
    int forcing_re_start = re_idx;

    for (int ti = 1; ti < n_times; ti++) {
      double t0 = time_pts(ti - 1);
      double dt = time_pts(ti) - t0;
      double sub_dt = dt / n_quadrature;
      int cov_idx = covariate_index_at(t0, obs_t);

      // Forcing varies with covariates at t0
      Type forcing(0);
      for (int k = 0; k < n_fixed_cov; k++)
        forcing += longitudinal(forcing_fixed_start + k) * Type(X_fixed(cov_idx, k));
      for (int k = 0; k < n_random_cov; k++)
        forcing += bi(forcing_re_start + k) * Type(X_random(cov_idx, k));

      // Composite Simpson's rule for cumulative hazard
      Type h_left = eval_hazard(m, v, t0, spline_degree, knots_d, boundary_d,
                                baseline, hazard, surv_cov, gamma_val);
      for (int q = 0; q < n_quadrature; q++) {
        double t_start = t0 + q * sub_dt;
        // Midpoint
        Type m_mid = m, v_mid = v;
        ode_step(m_mid, v_mid, b1, b2, forcing, Type(sub_dt * 0.5), branch);
        m_mid = clamp(m_mid, -BC, BC);
        v_mid = clamp(v_mid, -BC, BC);
        Type h_mid = eval_hazard(m_mid, v_mid, t_start + sub_dt * 0.5,
                                 spline_degree, knots_d, boundary_d,
                                 baseline, hazard, surv_cov, gamma_val);
        // Right endpoint
        ode_step(m_mid, v_mid, b1, b2, forcing, Type(sub_dt * 0.5), branch);
        m_mid = clamp(m_mid, -BC, BC);
        v_mid = clamp(v_mid, -BC, BC);
        Type h_right = eval_hazard(m_mid, v_mid, t_start + sub_dt,
                                   spline_degree, knots_d, boundary_d,
                                   baseline, hazard, surv_cov, gamma_val);
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
        fitted_velocity(oi + j)  = sol_v(time_idx);
      }
    }

    // --- Survival likelihood: Cox PH ---
    int event_idx = find_time_index(time_pts, event_time);
    if (event_idx >= 0) {
      nll += sol_H(event_idx);  // cumulative hazard contributes to nll
      Type h_event = eval_hazard(sol_m(event_idx), sol_v(event_idx),
                                 event_time, spline_degree, knots_d, boundary_d,
                                 baseline, hazard, surv_cov, gamma_val);
      Type log_h_event = log(h_event);
      if (status == 1) nll -= log_h_event;
      cumulative_hazard(i) = sol_H(event_idx);
      log_hazard_at_event(i) = log_h_event;
    }
  }

  // Reconstruct Sigma_b for reporting
  matrix<Type> Sigma_b = corr_structure.cov();
  for (int i = 0; i < n_random_effects; i++)
    for (int j = 0; j < n_random_effects; j++)
      Sigma_b(i, j) *= exp(log_sd_re(i)) * exp(log_sd_re(j));

  REPORT(fitted_biomarker);
  REPORT(fitted_velocity);
  REPORT(cumulative_hazard);
  REPORT(log_hazard_at_event);
  REPORT(Sigma_b);
  ADREPORT(Sigma_b);
  return nll;
}

#endif
