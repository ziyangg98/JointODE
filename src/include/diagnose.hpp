// Subject-level second-order ODE diagnostic objective
#ifndef DIAGNOSE_ODE_HPP
#define DIAGNOSE_ODE_HPP

#include "utils.hpp"

template<class Type>
Type diagnose_ode_nll(objective_function<Type>* obj) {
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

  DATA_VECTOR(time);
  DATA_VECTOR(y);
  PARAMETER_VECTOR(theta);

  Type m = theta(0);
  Type v = theta(1);
  Type b1 = -exp(theta(2));
  Type b2 = -exp(theta(3));
  Type forcing = theta(4);

  int n = y.size();
  Type nll(0);
  vector<Type> fitted(n);
  Type previous_time(0);

  for (int i = 0; i < n; ++i) {
    Type dt = time(i) - previous_time;
    ode_step(m, v, b1, b2, forcing, dt);
    fitted(i) = m;
    Type residual = y(i) - m;
    nll += residual * residual;
    previous_time = time(i);
  }

  REPORT(fitted);

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
  return nll;
}

#endif
