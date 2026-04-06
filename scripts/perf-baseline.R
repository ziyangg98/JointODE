#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name, "=")
  hit <- grep(paste0("^", key), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(key, "", hit[1], fixed = TRUE)
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

subset_sim_data <- function(n_subjects) {
  env <- environment()
  utils::data(list = "sim", package = "JointODE", envir = env)
  sim_data <- get("sim", envir = env)
  ids <- unique(sim_data$data$longitudinal_data$id)[seq_len(n_subjects)]

  longitudinal_data <- sim_data$data$longitudinal_data[
    sim_data$data$longitudinal_data$id %in% ids,
    c("id", "time", "observed", "x1", "x2")
  ]
  survival_data <- sim_data$data$survival_data[
    sim_data$data$survival_data$id %in% ids,
  ]

  list(longitudinal_data = longitudinal_data, survival_data = survival_data)
}

run_case <- function(case_name, dat, control, reps) {
  elapsed <- numeric(reps)
  iterations <- integer(reps)
  converged <- logical(reps)
  loglik <- numeric(reps)

  for (i in seq_len(reps)) {
    gc(FALSE)
    t0 <- proc.time()["elapsed"]

    fit <- JointODE(
      longitudinal_formula = observed ~ biomarker + velocity + x1 + x2 +
        (biomarker + velocity | id),
      survival_formula = Surv(time, status) ~ w1 + w2,
      longitudinal_data = dat$longitudinal_data,
      survival_data = dat$survival_data,
      init = "marginal",
      control = control
    )

    elapsed[i] <- proc.time()["elapsed"] - t0
    iterations[i] <- fit$convergence$iterations
    converged[i] <- isTRUE(fit$convergence$converged)
    loglik[i] <- fit$logLik
  }

  data.frame(
    case = case_name,
    reps = reps,
    elapsed_mean = mean(elapsed),
    elapsed_median = median(elapsed),
    elapsed_p90 = as.numeric(stats::quantile(elapsed, 0.9, names = FALSE)),
    iter_mean = mean(iterations),
    converged_rate = mean(converged),
    logLik_mean = mean(loglik),
    stringsAsFactors = FALSE
  )
}

n_subjects <- as.integer(get_arg("n", "20"))
reps <- as.integer(get_arg("reps", "3"))
maxit <- as.integer(get_arg("maxit", "10"))
tol <- as.numeric(get_arg("tol", "1e-2"))
hazard_quadrature <- as.integer(get_arg("hazard_quadrature", "1"))
out <- get_arg("out", NULL)
use_parallel_case <- as_bool(get_arg("parallel_case", "true"), TRUE)

if (is.na(n_subjects) || n_subjects < 2) stop("--n must be >= 2")
if (is.na(reps) || reps < 1) stop("--reps must be >= 1")
if (is.na(maxit) || maxit < 1) stop("--maxit must be >= 1")
if (is.na(tol) || tol <= 0) stop("--tol must be > 0")
if (is.na(hazard_quadrature) || hazard_quadrature < 1) {
  stop("--hazard_quadrature must be >= 1")
}

dat <- subset_sim_data(n_subjects)
available_cores <- parallel::detectCores(logical = TRUE)

control_seq <- list(
  maxit = maxit,
  tol = tol,
  verbose = FALSE,
  parallel = FALSE,
  n_cores = 0,
  hazard_quadrature = hazard_quadrature
)

results <- run_case("sequential", dat, control_seq, reps)

if (use_parallel_case && !is.na(available_cores) && available_cores > 1) {
  control_par <- control_seq
  control_par$parallel <- TRUE
  control_par$n_cores <- min(2L, available_cores)
  results <- rbind(results, run_case("parallel_2cores", dat, control_par, reps))
}

cat("\nJointODE performance baseline\n")
cat(sprintf("subjects=%d, reps=%d, maxit=%d, tol=%g, hazard_quadrature=%d\n\n",
  n_subjects, reps, maxit, tol, hazard_quadrature
))
print(results, row.names = FALSE, digits = 4)

if (!is.null(out) && nzchar(out)) {
  utils::write.csv(results, out, row.names = FALSE)
  cat(sprintf("\nSaved: %s\n", out))
}
