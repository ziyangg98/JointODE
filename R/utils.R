# Utility Functions for JointODE Package

# Shared Helpers ===============================================================

#' @noRd
.n_obs <- function(data_list) {
  sum(vapply(
    data_list,
    function(s) length(s$longitudinal$measurements),
    integer(1)
  ))
}

#' @noRd
.coef_table <- function(estimates, std_errors) {
  z <- estimates / std_errors
  cbind(
    Estimate = estimates, `Std. Error` = std_errors,
    `z value` = z, `Pr(>|z|)` = 2 * pnorm(-abs(z))
  )
}

#' @noRd
.print_coefmat <- function(x, digits = 4, signif.stars = TRUE, ...) {
  x[, 1:2] <- round(x[, 1:2], digits)
  printCoefmat(x, digits = digits, signif.stars = signif.stars, ...)
}

#' @noRd
.setup_openmp <- function(control) {
  if (control$parallel) {
    n_cores <- if (control$n_cores > 0) control$n_cores else parallel::detectCores()
    TMB::openmp(n_cores)
  }
}

#' @noRd
.fit_tmb <- function(tmb_data, tmb_params, control, map = NULL) {
  opt_control <- list(
    iter.max = control$maxit,
    eval.max = control$maxit,
    rel.tol = control$tol,
    trace = as.integer(control$verbose >= 2)
  )

  args <- list(
    data = tmb_data,
    parameters = tmb_params,
    random = "random_effects",
    DLL = "JointODE",
    silent = control$verbose < 3
  )
  if (!is.null(map)) {
    if (!is.null(map$corr_par)) {
      tmb_params$corr_par[is.na(map$corr_par)] <- 0
      args$parameters <- tmb_params
    }
    args$map <- map
  }

  obj <- do.call(TMB::MakeADFun, args)
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = opt_control)

  list(obj = obj, opt = opt)
}

#' @noRd
.fit_tmb_em <- function(tmb_data, tmb_params, control, map = NULL) {
  opt_control <- list(
    iter.max = control$maxit,
    eval.max = control$maxit,
    rel.tol = control$tol,
    trace = as.integer(control$verbose >= 2)
  )
  variance_map <- list(
    log_sd_re = factor(rep(NA, length(tmb_params$log_sd_re))),
    corr_par = factor(rep(NA, length(tmb_params$corr_par)))
  )
  if (!is.null(map) && !is.null(map$corr_par)) {
    tmb_params$corr_par[is.na(map$corr_par)] <- 0
    variance_map$corr_par <- map$corr_par
  }

  previous_fixed <- NULL
  fixed_change <- Inf
  for (iter in seq_len(control$maxit)) {
    if (control$verbose >= 1) {
      cli::cli_alert_info("Laplace-EM iteration {iter}")
    }
    obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = tmb_params,
      random = "random_effects",
      map = variance_map,
      DLL = "JointODE",
      silent = control$verbose < 3
    )
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = opt_control)
    obj$fn(opt$par)
    par <- obj$env$parList()
    fixed <- c(par$longitudinal, par$initial_state, par$log_sigma_e)

    random_effects <- matrix(
      par$random_effects,
      nrow = tmb_data$n_subjects,
      ncol = tmb_data$n_random_effects
    )
    sigma_b <- .em_random_effect_sigma(obj, random_effects, tmb_data$n_random_effects)
    corr <- .pack_correlation_theta(sigma_b, tmb_data$n_random_effects)

    tmb_params$longitudinal <- par$longitudinal
    tmb_params$initial_state <- par$initial_state
    tmb_params$log_sigma_e <- par$log_sigma_e
    tmb_params$random_effects <- random_effects
    tmb_params$log_sd_re <- corr$log_sd_re
    tmb_params$corr_par <- corr$corr_par
    if (!is.null(map) && !is.null(map$corr_par)) {
      tmb_params$corr_par[is.na(map$corr_par)] <- 0
    }

    fixed_change <- if (is.null(previous_fixed)) {
      Inf
    } else {
      max(abs(fixed - previous_fixed))
    }
    if (control$verbose >= 1) {
      cli::cli_alert_info(sprintf(
        "objective = %.4f; fixed change = %.3g",
        opt$objective, fixed_change
      ))
    }
    if (is.finite(fixed_change) && fixed_change < control$tol) break
    previous_fixed <- fixed
  }

  if (is.finite(fixed_change) && fixed_change < control$tol) {
    opt$convergence <- 0L
    opt$message <- "EM fixed parameters converged"
  } else {
    opt$convergence <- 1L
    opt$message <- "EM reached the iteration limit"
  }

  list(obj = obj, opt = opt)
}

#' @noRd
.em_random_effect_sigma <- function(obj, random_effects, n_re) {
  H <- obj$env$spHess(random = TRUE)
  n_subjects <- nrow(random_effects)
  sigma_b <- crossprod(random_effects)
  for (i in seq_len(n_subjects)) {
    idx <- i + (seq_len(n_re) - 1L) * n_subjects
    Hi <- as.matrix(H[idx, idx, drop = FALSE])
    Hi <- (Hi + t(Hi)) / 2
    sigma_b <- sigma_b + solve(Hi)
  }
  sigma_b <- sigma_b / n_subjects
  sigma_b <- (sigma_b + t(sigma_b)) / 2
  ev <- eigen(sigma_b, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0 || any(!is.finite(ev))) {
    sigma_b <- sigma_b + diag(abs(min(ev, na.rm = TRUE)) + 1e-6, n_re)
  }
  sigma_b
}

#' @noRd
.correlation_map <- function(parsed_long, tmb_params) {
  if (!isTRUE(parsed_long$diagonal)) return(NULL)
  list(corr_par = factor(rep(NA, length(tmb_params$corr_par))))
}

#' @noRd
.estimate_time_scale <- function(data_list) {
  y <- unlist(lapply(data_list, function(d) d$longitudinal$measurements))
  slopes <- unlist(lapply(data_list, function(d) {
    tt <- d$longitudinal$times
    yy <- d$longitudinal$measurements
    if (length(tt) < 2) return(numeric(0))
    dt <- diff(tt)
    dy <- diff(yy)
    dy[dt > 0] / dt[dt > 0]
  }))
  slope <- median(abs(slopes[is.finite(slopes) & slopes != 0]), na.rm = TRUE)
  scale <- stats::sd(y, na.rm = TRUE) / slope
  if (is.finite(scale) && scale > 0) scale else 1
}

#' @noRd
.scale_data_time <- function(data_list, scale) {
  if (!is.finite(scale) || scale <= 0 || scale == 1) return(data_list)
  lapply(data_list, function(d) {
    d$longitudinal$times <- d$longitudinal$times / scale
    if (!is.null(d$time)) d$time <- d$time / scale
    d
  })
}

#' @noRd
.adjust_longitudinal_time_scale <- function(x, parsed_long, scale, to_internal) {
  if (!is.finite(scale) || scale <= 0 || scale == 1) return(x)
  direction <- if (to_internal) 1 else -1
  forcing_mult <- if (to_internal) scale^2 else 1 / scale^2
  idx <- 1L
  if (parsed_long$biomarker$fixed) {
    x[idx] <- x[idx] + direction * 2 * log(scale)
    idx <- idx + 1L
  }
  if (parsed_long$velocity$fixed) {
    x[idx] <- x[idx] + direction * log(scale)
    idx <- idx + 1L
  }
  if (idx <= length(x)) x[idx:length(x)] <- x[idx:length(x)] * forcing_mult
  x
}

#' @noRd
.random_effect_time_scale <- function(parsed_long, n_re, scale, to_internal) {
  mult <- rep(1, n_re)
  if (n_re >= 2) mult[2] <- scale
  forcing_start <- 3L + sum(
    parsed_long$biomarker$random, parsed_long$velocity$random
  )
  if (forcing_start <= n_re) mult[forcing_start:n_re] <- scale^2
  if (to_internal) mult else 1 / mult
}

#' @noRd
.longitudinal_vcov_time_scale <- function(parsed_long, n_long, scale) {
  if (!is.finite(scale) || scale <= 0 || scale == 1) return(rep(1, n_long))
  mult <- rep(1 / scale^2, n_long)
  idx <- 1L
  if (parsed_long$biomarker$fixed) {
    mult[idx] <- 1
    idx <- idx + 1L
  }
  if (parsed_long$velocity$fixed) {
    mult[idx] <- 1
    idx <- idx + 1L
  }
  mult
}

#' @noRd
.prepare_marginal_time_scale <- function(parameters, parsed_long, scale) {
  if (!is.finite(scale) || scale <= 0 || scale == 1) return(parameters)
  parameters$coefficients$longitudinal <- .adjust_longitudinal_time_scale(
    parameters$coefficients$longitudinal, parsed_long, scale, TRUE
  )
  parameters$coefficients$initial_state[2] <-
    parameters$coefficients$initial_state[2] * scale

  re_mult <- .random_effect_time_scale(
    parsed_long, ncol(parameters$random_effects_init), scale, TRUE
  )
  parameters$random_effects_init <- sweep(
    parameters$random_effects_init, 2, re_mult, `*`
  )
  parameters$coefficients$random_effect_sigma <-
    parameters$coefficients$random_effect_sigma * outer(re_mult, re_mult)
  parameters
}

#' @noRd
.restore_marginal_time_scale <- function(results, parsed_long, scale) {
  if (!is.finite(scale) || scale <= 0 || scale == 1) return(results)
  n_long <- length(results$parameters) - 2L
  long_idx <- seq_len(n_long)
  init_idx <- n_long + seq_len(2L)

  results$parameters[long_idx] <- .adjust_longitudinal_time_scale(
    results$parameters[long_idx], parsed_long, scale, FALSE
  )
  results$parameters[init_idx[2]] <- results$parameters[init_idx[2]] / scale

  vscale <- c(
    .longitudinal_vcov_time_scale(parsed_long, n_long, scale),
    1, 1 / scale
  )
  results$vcov <- results$vcov * outer(vscale, vscale)

  re_mult <- .random_effect_time_scale(
    parsed_long, ncol(results$random_effects), scale, FALSE
  )
  results$random_effects <- sweep(results$random_effects, 2, re_mult, `*`)
  results$random_effect_sigma <-
    results$random_effect_sigma * outer(re_mult, re_mult)
  results
}

#' @noRd
.prepare_joint_time_scale <- function(parameters, parsed_long, scale, gamma) {
  parameters <- .prepare_marginal_time_scale(parameters, parsed_long, scale)
  if (!is.finite(scale) || scale <= 0 || scale == 1) return(parameters)
  parameters$coefficients$baseline <-
    parameters$coefficients$baseline + log(scale)
  if (length(parameters$coefficients$hazard) >= 2) {
    parameters$coefficients$hazard[2] <-
      parameters$coefficients$hazard[2] / scale^gamma
  }
  bc <- parameters$configurations$baseline
  bc$knots <- bc$knots / scale
  bc$boundary_knots <- bc$boundary_knots / scale
  parameters$configurations$baseline <- bc
  parameters
}

#' @noRd
.restore_joint_time_scale <- function(results, parsed_long, scale, gamma) {
  if (!is.finite(scale) || scale <= 0 || scale == 1) return(results)
  cf <- results$parameters$coefficients
  cf$baseline <- cf$baseline - log(scale)
  if (length(cf$hazard) >= 2) cf$hazard[2] <- cf$hazard[2] * scale^gamma
  cf$longitudinal <- .adjust_longitudinal_time_scale(
    cf$longitudinal, parsed_long, scale, FALSE
  )
  cf$initial_state[2] <- cf$initial_state[2] / scale

  re_mult <- .random_effect_time_scale(
    parsed_long, ncol(results$random_effects), scale, FALSE
  )
  cf$random_effect_sigma <- cf$random_effect_sigma * outer(re_mult, re_mult)
  results$random_effects <- sweep(results$random_effects, 2, re_mult, `*`)

  bc <- results$parameters$configurations$baseline
  bc$knots <- bc$knots * scale
  bc$boundary_knots <- bc$boundary_knots * scale
  results$parameters$configurations$baseline <- bc

  n_baseline <- length(cf$baseline)
  n_hazard <- length(cf$hazard)
  n_long <- length(cf$longitudinal)
  vscale <- c(
    rep(1, n_baseline),
    c(1, scale^gamma, rep(1, max(0, n_hazard - 2L))),
    .longitudinal_vcov_time_scale(parsed_long, n_long, scale),
    1, 1 / scale
  )
  results$vcov <- results$vcov * outer(vscale, vscale)
  results$parameters$coefficients <- cf
  results
}

#' Initialize RE columns 1-2 from first observations
#' @noRd
.init_re_from_observations <- function(data_list, m0, v0 = 0, n_re = 2L) {
  n_subjects <- length(data_list)
  re <- matrix(0, nrow = n_subjects, ncol = max(n_re, 2L))
  for (i in seq_along(data_list)) {
    obs <- data_list[[i]]$longitudinal
    if (length(obs$measurements) >= 1) {
      re[i, 1] <- obs$measurements[1] - m0
    }
    if (length(obs$measurements) >= 2) {
      dt <- obs$times[2] - obs$times[1]
      if (dt > 0) {
        re[i, 2] <- (obs$measurements[2] - obs$measurements[1]) / dt - v0
      }
    }
  }
  re
}

#' Build counting process (start, stop, event) from predictions + survival
#' @noRd
.build_counting_process <- function(pred, survival_data, id, surv_vars) {
  names(pred)[names(pred) == "time"] <- "obstime"
  merged <- merge(pred, survival_data, by = id, all.x = TRUE)
  merged <- merged[order(merged[[id]], merged$obstime), ]

  merged$start <- merged$obstime
  subj_ids <- as.vector(merged[[id]])
  subj_idx <- c(which(!duplicated(subj_ids)), nrow(merged) + 1)
  ni <- diff(subj_idx)

  stop_vec <- event_vec <- numeric(nrow(merged))
  for (i in seq_along(ni)) {
    s <- subj_idx[i]
    e <- subj_idx[i + 1] - 1
    if (ni[i] == 1) {
      stop_vec[s] <- merged[[surv_vars[1]]][s]
      event_vec[s] <- merged[[surv_vars[2]]][s]
    } else {
      stop_vec[s:e] <- c(merged$obstime[(s + 1):e], merged[[surv_vars[1]]][e])
      event_vec[s:e] <- c(rep(0, ni[i] - 1), merged[[surv_vars[2]]][e])
    }
  }
  merged$stop <- stop_vec
  merged$event <- event_vec
  merged[merged$stop > merged$start, , drop = FALSE]
}

# Constants ====================================================================

.default_spline <- list(
  degree = 2,
  n_knots = 1,
  knot_placement = "equal",
  boundary_knots = NULL
)

.reserved_words <- c("biomarker", "velocity")

.init_state_names <- c("initial_biomarker", "initial_velocity")

# Formula Parsing ==============================================================

#' @importFrom stats as.formula
#' @noRd
.build_formula <- function(terms, response = NULL, is_random = FALSE) {
  covars <- setdiff(terms, "(Intercept)")
  has_intercept <- "(Intercept)" %in% terms

  terms_str <- if (length(covars) > 0) {
    if (is_random) covars <- gsub("^\\(Intercept\\)$", "1", covars)
    paste(covars, collapse = " + ")
  } else if (has_intercept) {
    "1"
  } else {
    "0"
  }

  if (!has_intercept && length(covars) > 0 && !is_random) {
    terms_str <- paste("0 +", terms_str)
  }

  if (is_random || is.null(response)) {
    as.formula(paste("~", terms_str))
  } else {
    as.formula(paste(response, "~", terms_str))
  }
}

#' @importFrom stats as.formula terms
#' @noRd
.parse_longitudinal_formula <- function(formula) {
  formula_str <- paste(deparse(formula, width.cutoff = 500L), collapse = "")
  re_pattern <- "\\(([^)]+\\|\\|?[^)]+)\\)"
  re_match <- regmatches(formula_str, regexpr(re_pattern, formula_str))

  if (length(re_match) == 0) {
    fixed_formula <- formula
    random_terms <- NULL
    grouping <- NULL
    biomarker_random <- FALSE
    velocity_random <- FALSE
    diagonal <- FALSE
  } else {
    inner <- substr(re_match, 2, nchar(re_match) - 1)
    diagonal <- grepl("\\|\\|", inner)
    parts <- strsplit(inner, if (diagonal) "\\|\\|" else "\\|")[[1]]
    grouping <- trimws(parts[2])
    re_terms_str <- trimws(parts[1])

    terms <- trimws(strsplit(re_terms_str, "\\+")[[1]])
    no_intercept <- any(terms %in% c("0", "-1"))
    has_intercept <- !no_intercept || any(terms == "1")
    terms <- terms[!terms %in% c("0", "1", "-1")]
    random_terms <- c(if (has_intercept) "(Intercept)", terms)

    fixed_str <- sub(re_match, "", formula_str, fixed = TRUE)
    fixed_str <- gsub("\\s+\\+\\s*$", "", fixed_str)
    fixed_str <- gsub("~\\s*\\+\\s*", "~ ", fixed_str)
    fixed_str <- gsub("\\+\\s*\\+", "+", fixed_str)
    fixed_formula <- as.formula(fixed_str)

    biomarker_random <- "biomarker" %in% random_terms
    velocity_random <- "velocity" %in% random_terms
    random_terms <- setdiff(random_terms, .reserved_words)
    if (length(random_terms) == 0) random_terms <- NULL
  }

  fixed_terms_obj <- terms(fixed_formula)
  fixed_terms_labels <- attr(fixed_terms_obj, "term.labels")
  has_intercept <- attr(fixed_terms_obj, "intercept") == 1

  biomarker_in_fixed <- "biomarker" %in% fixed_terms_labels
  velocity_in_fixed <- "velocity" %in% fixed_terms_labels
  fixed_covariates <- setdiff(fixed_terms_labels, .reserved_words)

  if (has_intercept) {
    fixed_covariates <- c("(Intercept)", fixed_covariates)
  }

  list(
    response = as.character(formula[[2]]),
    fixed_terms = fixed_covariates,
    random_terms = random_terms,
    biomarker = list(fixed = biomarker_in_fixed, random = biomarker_random),
    velocity = list(fixed = velocity_in_fixed, random = velocity_random),
    grouping = grouping,
    diagonal = diagonal
  )
}

#' @importFrom stats terms
#' @noRd
.parse_survival_formula <- function(formula) {
  surv_response <- formula[[2]]

  if (
    !inherits(surv_response, "call") ||
      !identical(as.character(surv_response[[1]]), "Surv")
  ) {
    stop(
      "Survival formula must have Surv() on the left-hand side",
      call. = FALSE
    )
  }

  surv_call_vars <- all.vars(surv_response)
  if (length(surv_call_vars) < 2) {
    stop("Surv() must have at least time and status arguments", call. = FALSE)
  }

  time_var <- surv_call_vars[1]
  status_var <- surv_call_vars[2]

  all_terms <- attr(terms(formula), "term.labels")
  covariate_terms <- if (length(all_terms) == 0) NULL else all_terms

  list(
    time_var = time_var,
    status_var = status_var,
    covariate_terms = covariate_terms
  )
}

# Model Configuration ==========================================================

#' @noRd
.compute_dimensions <- function(parsed_long, parsed_surv, spline_config) {
  n_longitudinal_fixed <- length(parsed_long$fixed_terms)
  n_longitudinal_random <- if (is.null(parsed_long$random_terms)) {
    0
  } else {
    length(parsed_long$random_terms)
  }
  n_survival_covariates <- if (is.null(parsed_surv$covariate_terms)) {
    0
  } else {
    length(parsed_surv$covariate_terms)
  }

  n_ode_fixed <- sum(parsed_long$biomarker$fixed, parsed_long$velocity$fixed)
  n_ode_random <- sum(parsed_long$biomarker$random, parsed_long$velocity$random)

  list(
    n_longitudinal_coef = n_longitudinal_fixed + n_ode_fixed,
    n_random_effects = n_longitudinal_random + n_ode_random + 2,
    n_survival_covariates = n_survival_covariates,
    n_spline_basis = spline_config$degree + spline_config$n_knots + 1,
    spline_config = spline_config
  )
}

#' @importFrom stats quantile
#' @noRd
.get_spline_config <- function(
  x,
  degree = 2,
  n_knots = 1,
  knot_placement = "quantile",
  boundary_knots = NULL
) {
  if (is.null(boundary_knots)) {
    boundary_knots <- range(x, na.rm = TRUE)
  }

  if (knot_placement == "quantile") {
    probs <- seq(0, 1, length.out = n_knots + 2)[-c(1, n_knots + 2)]
    knots <- quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  } else if (knot_placement == "equal") {
    knots <- seq(
      boundary_knots[1], boundary_knots[2],
      length.out = n_knots + 2
    )[-c(1, n_knots + 2)]
  } else {
    stop("knot_placement must be 'quantile' or 'equal'")
  }

  list(
    degree = degree,
    knots = knots,
    boundary_knots = boundary_knots,
    df = length(knots) + degree + 1
  )
}

# Parameter Helpers =============================================================

#' @noRd
.prefixed_coef_names <- function(coef_names) {
  c(
    paste0("baseline:", coef_names$baseline),
    paste0("hazard:", coef_names$hazard),
    paste0("longitudinal:", coef_names$longitudinal),
    paste0("initial state:", coef_names$initial_state)
  )
}

#' @noRd
.count_params <- function(parameters) {
  cf <- parameters$coefficients
  p <- nrow(cf$random_effect_sigma)
  length(cf$baseline) + length(cf$hazard) + length(cf$longitudinal) +
    length(cf$initial_state) + 1 + p * (p + 1) / 2
}

# Prediction Utilities =========================================================

#' ODE step via Cayley-Hamilton (sinhc/sinc formulation)
#' @noRd
.ode_step_r <- function(m, v, b1, b2, f, dt) {
  h <- b2 * dt / 2
  D <- b2^2 + 4 * b1
  s2 <- D * dt^2 / 4

  if (D >= 0) {
    delta <- sqrt(max(s2, 0))
    sc <- if (delta^2 < 1e-8) 1 + delta^2 / 6 else sinh(delta) / delta
    cc <- cosh(delta)
  } else {
    omega <- sqrt(max(-s2, 0))
    sc <- if (omega^2 < 1e-8) 1 - omega^2 / 6 else sin(omega) / omega
    cc <- cos(omega)
  }

  em <- exp(h)
  a0 <- em * (cc - h * sc)
  a1 <- em * dt * sc

  # Forcing integral: J1 = (a0-1)/b1, Taylor at b1~0
  if (abs(b1) > 1e-10) {
    J1 <- (a0 - 1) / b1
  } else {
    x <- b2 * dt
    J1 <- dt^2 * if (x^2 < 1e-8) 0.5 + x / 6 else (exp(x) - 1 - x) / x^2
  }

  c(
    a0 * m + a1 * v + f * J1,
    b1 * a1 * m + (a0 + b2 * a1) * v + f * a1
  )
}

#' Predict marginal trajectories for all subjects at given times
#' @noRd
.predict_marginal_trajectories <- function(data_list, times, cf, configs, re) {
  lc <- unname(cf$longitudinal)

  ids <- names(data_list)
  n_times <- length(times)
  out <- vector("list", length(ids))

  for (i in seq_along(data_list)) {
    subj <- data_list[[i]]
    obs_t <- subj$longitudinal$times
    fcov <- subj$longitudinal$covariates$fixed
    rcov <- subj$longitudinal$covariates$random

    bi <- re[i, ]
    m <- cf$initial_state[1] + bi[1]
    v <- cf$initial_state[2] + bi[2]

    log_omega2 <- 0
    log_2xi_omega <- 0
    fi <- 1
    ri <- 3
    if (configs$biomarker$fixed) {
      log_omega2 <- log_omega2 + lc[fi]
      fi <- fi + 1
    }
    if (configs$biomarker$random) {
      log_omega2 <- log_omega2 + bi[ri]
      ri <- ri + 1
    }
    if (configs$velocity$fixed) {
      log_2xi_omega <- log_2xi_omega + lc[fi]
      fi <- fi + 1
    }
    if (configs$velocity$random) {
      log_2xi_omega <- log_2xi_omega + bi[ri]
      ri <- ri + 1
    }
    b1 <- -exp(log_omega2)
    b2 <- -exp(log_2xi_omega)
    ffs <- fi
    rfs <- ri

    bio <- vel <- numeric(n_times)
    prev_t <- 0

    for (ti in seq_len(n_times)) {
      tgt <- times[ti]
      dt <- tgt - prev_t
      if (dt > 0 && length(obs_t) > 0) {
        ci <- findInterval(prev_t, obs_t)
        ci <- max(1, min(ci, nrow(fcov)))
        eta <- sum(lc[ffs:(ffs + ncol(fcov) - 1)] * fcov[ci, ])
        if (ncol(rcov) > 0) {
          eta <- eta + sum(bi[rfs:(rfs + ncol(rcov) - 1)] * rcov[ci, ])
        }

        mv <- .ode_step_r(m, v, b1, b2, eta, dt)
        m <- mv[1]
        v <- mv[2]
      }
      bio[ti] <- m
      vel[ti] <- v
      prev_t <- tgt
    }
    out[[i]] <- data.frame(
      id = ids[i], time = times,
      biomarker = bio, velocity = vel,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

#' Hazard at time t (double)
#' @noRd
.eval_hazard_r <- function(m, v, t, baseline, hazard, surv_cov,
                           spline_config, gamma) {
  basis <- splines2::bSpline(
    t,
    knots = spline_config$knots,
    Boundary.knots = spline_config$boundary_knots,
    degree = spline_config$degree, intercept = TRUE
  )
  log_h <- sum(basis * baseline) + hazard[1] * m
  if (gamma == 1) {
    log_h <- log_h + hazard[2] * v
  } else if (gamma == 2) log_h <- log_h + hazard[2] * v^2
  if (length(surv_cov) > 0) {
    log_h <- log_h + sum(hazard[-(1:2)] * surv_cov)
  }
  exp(log_h)
}

#' Predict trajectories for all subjects at given times
#' @noRd
.predict_trajectories <- function(
  data_list, times, cf, configs, re, n_quad
) {
  if (is.null(n_quad)) n_quad <- 1L
  n_quad <- max(n_quad, 1L)
  bl <- unname(cf$baseline)
  hz <- unname(cf$hazard)
  lc <- unname(cf$longitudinal)
  sbc <- configs$baseline

  ids <- names(data_list)
  n_subj <- length(ids)
  n_times <- length(times)

  out <- vector("list", n_subj)
  for (i in seq_along(data_list)) {
    subj <- data_list[[i]]
    obs_t <- subj$longitudinal$times
    fcov <- subj$longitudinal$covariates$fixed
    rcov <- subj$longitudinal$covariates$random
    surv_cov <- if (!is.null(subj$covariates)) {
      as.numeric(subj$covariates)
    } else {
      numeric(0)
    }

    bi <- re[i, ]
    m <- cf$initial_state[1] + bi[1]
    v <- cf$initial_state[2] + bi[2]

    log_omega2 <- 0
    log_2xi_omega <- 0
    fi <- 1
    ri <- 3
    if (configs$biomarker$fixed) {
      log_omega2 <- log_omega2 + lc[fi]
      fi <- fi + 1
    }
    if (configs$biomarker$random) {
      log_omega2 <- log_omega2 + bi[ri]
      ri <- ri + 1
    }
    if (configs$velocity$fixed) {
      log_2xi_omega <- log_2xi_omega + lc[fi]
      fi <- fi + 1
    }
    if (configs$velocity$random) {
      log_2xi_omega <- log_2xi_omega + bi[ri]
      ri <- ri + 1
    }
    b1 <- -exp(log_omega2)
    b2 <- -exp(log_2xi_omega)
    ffs <- fi
    rfs <- ri

    bio <- vel <- cumhaz <- numeric(n_times)
    ch <- 0
    prev_t <- 0

    for (ti in seq_len(n_times)) {
      tgt <- times[ti]
      dt <- tgt - prev_t
      if (dt > 0 && length(obs_t) > 0) {
        ci <- findInterval(prev_t, obs_t)
        ci <- max(1, min(ci, nrow(fcov)))
        forcing <- sum(lc[ffs:(ffs + ncol(fcov) - 1)] * fcov[ci, ])
        if (ncol(rcov) > 0) {
          forcing <- forcing +
            sum(bi[rfs:(rfs + ncol(rcov) - 1)] * rcov[ci, ])
        }
        sub_dt <- dt / n_quad
        hl <- .eval_hazard_r(m, v, prev_t, bl, hz, surv_cov, sbc, configs$gamma)
        for (q in seq_len(n_quad)) {
          ts <- prev_t + (q - 1) * sub_dt
          # Midpoint
          mv <- .ode_step_r(m, v, b1, b2, forcing, sub_dt / 2)
          mm <- mv[1]
          vv <- mv[2]
          hm <- .eval_hazard_r(mm, vv, ts + sub_dt / 2, bl, hz, surv_cov, sbc, configs$gamma)
          # Right
          mv <- .ode_step_r(mm, vv, b1, b2, forcing, sub_dt / 2)
          mm <- mv[1]
          vv <- mv[2]
          hr <- .eval_hazard_r(mm, vv, ts + sub_dt, bl, hz, surv_cov, sbc, configs$gamma)
          ch <- ch + (sub_dt / 6) * (hl + 4 * hm + hr)
          m <- mm
          v <- vv
          hl <- hr
        }
      }
      bio[ti] <- m
      vel[ti] <- v
      cumhaz[ti] <- ch
      prev_t <- tgt
    }
    out[[i]] <- data.frame(
      id = ids[i], time = times,
      biomarker = bio, velocity = vel,
      cumhaz = cumhaz, survival = exp(-cumhaz),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}
