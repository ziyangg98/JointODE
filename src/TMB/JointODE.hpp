// JointODE.hpp — TMB entry point with model dispatch
#ifndef JOINT_ODE_HPP
#define JOINT_ODE_HPP

#include "../include/joint.hpp"
#include "../include/marginal.hpp"

template<class Type>
Type objective_function<Type>::operator()() {
  DATA_INTEGER(model_type);
  if (model_type == 1) return marginal_ode_nll(this);
  return joint_ode_nll(this);
}

#endif
