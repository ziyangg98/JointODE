#!/usr/bin/env Rscript
# Diagnose 6x performance regression between old (deee3aa) and new version

library(JointODE)
data(sim)
TMB::openmp(8)

long_data <- sim$data$longitudinal_data[, c("id", "time", "observed", "x1", "x2")]
surv_data <- sim$data$survival_data

cat("========== Step 1: Component profiling (NEW version) ==========\n\n")

# Replicate JointODE() internals to profile each step
parsed_long <- JointODE:::.parse_longitudinal_formula(
  observed ~ biomarker + velocity + x1 + x2 + (biomarker + velocity | id))
parsed_surv <- JointODE:::.parse_survival_formula(Surv(time, status) ~ w1 + w2)

data_list <- JointODE:::.process_joint(
  long_data, survival_formula = Surv(time, status) ~ w1 + w2,
  survival_data = surv_data, parsed_long = parsed_long, parsed_surv = parsed_surv)

model_config <- JointODE:::.setup_model(
  long_data, surv_data, Surv(time, status) ~ w1 + w2,
  gamma = 1, parsed_long, parsed_surv,
  spline_config = list(degree = 2, n_knots = 1, knot_placement = "equal",
                       boundary_knots = NULL))

parameters <- sim$init
parameters$configurations$baseline <- model_config$spline_baseline_config
parameters$configurations$biomarker_clamp <- max(abs(unlist(
  lapply(data_list, function(d) d$longitudinal$measurements)))) * 5
parameters$configurations$hazard_quadrature <- 1L
parameters$configurations$gamma <- 1

n_re <- ncol(model_config$random_effects)
parameters$random_effects_init <- JointODE:::.init_re_from_observations(
  data_list, parameters$coefficients$initial_state[1],
  parameters$coefficients$initial_state[2], n_re)

control <- JointODE::JointODE.control(parallel = TRUE)

t_pack <- system.time({
  tmb_data <- JointODE:::.pack_joint_data(data_list, parameters, control)
  tmb_params <- JointODE:::.pack_joint_params(parameters)
})
cat(sprintf("Pack data+params: %.2f s\n", t_pack["elapsed"]))

t_make <- system.time({
  obj <- TMB::MakeADFun(
    data = tmb_data, parameters = tmb_params,
    random = "random_effects", DLL = "JointODE", silent = TRUE)
})
cat(sprintf("MakeADFun:        %.2f s\n", t_make["elapsed"]))

# Single fn/gr evaluation
t_fn <- system.time(for (i in 1:10) obj$fn(obj$par))
t_gr <- system.time(for (i in 1:10) obj$gr(obj$par))
cat(sprintf("fn (x10):         %.2f s  (%.3f s/call)\n", t_fn["elapsed"], t_fn["elapsed"]/10))
cat(sprintf("gr (x10):         %.2f s  (%.3f s/call)\n", t_gr["elapsed"], t_gr["elapsed"]/10))

t_opt <- system.time({
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
    control = list(iter.max = 300, eval.max = 3000, rel.tol = 1e-8))
})
cat(sprintf("nlminb:           %.2f s  (%d iter, %s)\n",
            t_opt["elapsed"], opt$iterations, opt$message))

t_sdr <- system.time({ sdr <- TMB::sdreport(obj) })
cat(sprintf("sdreport:         %.2f s\n", t_sdr["elapsed"]))

cat(sprintf("\nTOTAL:            %.2f s\n",
  t_pack["elapsed"] + t_make["elapsed"] + t_opt["elapsed"] + t_sdr["elapsed"]))

cat("\n========== Step 2: Data structure comparison ==========\n\n")

cat("NEW tmb_data structure:\n")
for (nm in names(tmb_data)) {
  v <- tmb_data[[nm]]
  cat(sprintf("  %-25s class=%-10s length=%d\n", nm, class(v)[1], length(v)))
}

cat("\nNEW tmb_params structure:\n")
for (nm in names(tmb_params)) {
  v <- tmb_params[[nm]]
  if (is.matrix(v)) {
    cat(sprintf("  %-25s matrix %dx%d\n", nm, nrow(v), ncol(v)))
  } else {
    cat(sprintf("  %-25s length=%d\n", nm, length(v)))
  }
}
