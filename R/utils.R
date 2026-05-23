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
.fit_variance_tmb <- function(tmb_data, tmb_params, control, fixed,
                              eval_factor = 2, direct_trust = FALSE) {
  mstep_maxit <- min(control$maxit, 50L)
  opt_control <- list(
    iter.max = mstep_maxit,
    eval.max = mstep_maxit * eval_factor,
    rel.tol = control$tol,
    trace = as.integer(control$verbose >= 3)
  )

  if (isTRUE(direct_trust)) {
    obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = tmb_params,
      random = "random_effects",
      DLL = "JointODE",
      silent = control$verbose < 3
    )

    par <- obj$par
    objective <- obj$fn(par)
    step_norm <- Inf
    objective_change <- Inf
    converged <- FALSE
    last_opt <- NULL
    last_message <- NULL
    trust_radius <- 1
    opt_control$iter.max <- min(mstep_maxit, 20L)
    opt_control$eval.max <- opt_control$iter.max * eval_factor

    for (i in seq_len(control$maxit)) {
      if (control$verbose > 0) {
        cli::cli_alert_info(sprintf("Trust-region iteration %d", i))
      }

      grad0 <- tryCatch(obj$gr(par), error = function(e) NA_real_)
      grad0_max <- if (all(is.finite(grad0))) max(abs(grad0)) else Inf
      par_scale <- rep(1, length(par))
      if (is.finite(grad0_max) && grad0_max > 1000) {
        par_scale <- 1 / sqrt(pmax(abs(grad0), 1))
        par_scale <- pmin(1, pmax(1e-3, par_scale))
      }
      fn_scaled <- function(delta) obj$fn(par + par_scale * delta)
      gr_scaled <- function(delta) obj$gr(par + par_scale * delta) * par_scale

      opt <- tryCatch(
        stats::nlminb(
          rep(0, length(par)), fn_scaled, gr_scaled,
          control = c(opt_control, list(step.max = trust_radius))
        ),
        error = function(e) {
          list(
            par = rep(0, length(par)), objective = Inf, convergence = 1L,
            iterations = 0L, message = conditionMessage(e)
          )
        }
      )
      last_opt <- opt
      last_message <- opt$message

      if (any(!is.finite(opt$par))) {
        last_message <- "optimizer proposed non-finite parameters"
        break
      }

      step <- par_scale * opt$par
      alpha <- 1
      accepted <- FALSE
      trial_objective <- Inf
      objective_tol <- 1e-8 * max(1, abs(objective))
      while (alpha >= 1e-6) {
        trial_par <- par + alpha * step
        trial_objective <- tryCatch(
          obj$fn(trial_par),
          error = function(e) Inf
        )
        trial_grad <- tryCatch(obj$gr(trial_par), error = function(e) NA_real_)
        trial_grad_max <- if (all(is.finite(trial_grad))) {
          max(abs(trial_grad))
        } else {
          Inf
        }
        grad_limit <- if (is.finite(grad0_max)) {
          max(1000, 25 * grad0_max)
        } else {
          Inf
        }
        if (
          is.finite(trial_objective) &&
            trial_objective <= objective + objective_tol &&
            trial_grad_max <= grad_limit
        ) {
          accepted <- TRUE
          break
        }
        alpha <- alpha / 2
      }

      if (!accepted) {
        last_message <- "trust-region step could not improve objective"
        break
      }

      step_norm <- max(abs(alpha * step))
      objective_change <- abs(objective - trial_objective)
      par <- trial_par
      objective <- trial_objective
      trust_radius <- if (alpha == 1) {
        min(2 * trust_radius, 10)
      } else {
        max(alpha * trust_radius, 1e-4)
      }

      if (control$verbose > 0) {
        grad <- tryCatch(obj$gr(par), error = function(e) NA_real_)
        grad_max <- if (all(is.finite(grad))) max(abs(grad)) else Inf
        cli::cli_alert_info(sprintf(
          paste0(
            "convergence = %d; objective = %.4f; ",
            "step = %.3g; alpha = %.3g; trust radius = %.3g; ",
            "max gradient = %.3g"
          ),
          opt$convergence, objective, step_norm, alpha, trust_radius,
          grad_max
        ))
      }

      grad <- tryCatch(obj$gr(par), error = function(e) NA_real_)
      grad_max <- if (all(is.finite(grad))) max(abs(grad)) else Inf
      if (
        step_norm < control$tol &&
          objective_change < control$tol * max(1, abs(objective))
      ) {
        converged <- opt$convergence == 0 || grad_max < sqrt(control$tol)
        if (!converged) {
          last_message <- sprintf(
            "trust-region no progress; max gradient %.3g; last optimizer: %s",
            grad_max, opt$message
          )
        }
        break
      }
    }

    opt <- list(
      par = par,
      objective = objective,
      convergence = if (converged) 0L else 1L,
      iterations = if (is.null(last_opt)) 0L else last_opt$iterations,
      message = sprintf(
        "%s; last optimizer: %s",
        if (converged) {
          "trust-region convergence"
        } else {
          "trust-region iterations stopped before convergence"
        },
        last_message
      ),
      outer_delta = step_norm,
      objective_change = objective_change
    )
    return(list(obj = obj, opt = opt))
  }

  variance_map <- list(
    log_sd_re = factor(rep(NA, length(tmb_params$log_sd_re))),
    corr_par = factor(rep(NA, length(tmb_params$corr_par)))
  )

  outer_delta <- Inf
  outer_converged <- FALSE
  final_params <- tmb_params
  last_opt <- NULL
  last_message <- NULL
  previous_objective <- Inf
  damped_em <- FALSE

  for (i in seq_len(control$maxit)) {
    if (control$verbose > 0) {
      cli::cli_alert_info(sprintf("Variance update iteration %d", i))
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
    last_opt <- opt
    last_message <- opt$message

    if (control$verbose > 0) {
      msg <- sprintf(
        "M-step %d: convergence = %d; objective = %.4f; %s",
        i, opt$convergence, opt$objective, opt$message
      )
      if (opt$convergence == 0) {
        cli::cli_alert_success(msg)
      } else {
        cli::cli_alert_warning(msg)
      }
    }

    if (!is.finite(opt$objective) || any(!is.finite(opt$par))) {
      last_message <- "M-step produced non-finite objective or parameters"
      break
    }

    best <- obj$env$parList(opt$par)
    fixed_delta <- max(abs(c(
      unlist(best[fixed]) - unlist(tmb_params[fixed]),
      best$log_sigma_e - tmb_params$log_sigma_e
    )))
    objective_tol <- 1e-6 * max(1, abs(previous_objective))
    if (opt$objective > previous_objective + objective_tol) damped_em <- TRUE
    damping_rho <- if (damped_em) 0.5 else 1
    next_params <- .update_variance_parameters(
      tmb_data, obj, opt$par, damping_rho
    )
    previous_objective <- opt$objective
    variance_delta <- max(abs(c(
      next_params$log_sd_re - tmb_params$log_sd_re,
      next_params$corr_par - tmb_params$corr_par
    )))
    outer_delta <- max(fixed_delta, variance_delta)
    final_params <- next_params

    if (control$verbose > 0) {
      cli::cli_alert_info(sprintf(
        paste0(
          "sigma_e = %.4g; fixed change = %.3g; ",
          "variance change = %.3g; damping rho = %.2g"
        ),
        exp(next_params$log_sigma_e), fixed_delta, variance_delta,
        damping_rho
      ))
      cli::cli_alert_info(sprintf("outer change = %.3g", outer_delta))
    }

    if (outer_delta < control$tol) {
      outer_converged <- TRUE
      break
    }
    tmb_params <- next_params
  }

  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = final_params,
    random = "random_effects",
    map = variance_map,
    DLL = "JointODE",
    silent = control$verbose < 3
  )
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = opt_control)
  opt$outer_delta <- outer_delta
  opt$m_step_convergence <- last_opt$convergence
  opt$m_step_message <- last_message
  opt$convergence <- if (outer_converged) 0L else 1L
  opt$message <- sprintf(
    "%s; last M-step: %s",
    if (outer_converged) {
      "outer parameter convergence"
    } else {
      "outer iterations stopped before EM convergence"
    },
    last_message
  )

  list(obj = obj, opt = opt)
}

#' @noRd
.update_variance_parameters <- function(tmb_data, obj, par, damping_rho = 1) {
  obj$fn(par)
  p <- obj$env$parList(par)
  b <- as.matrix(p$random_effects)
  n_subjects <- tmb_data$n_subjects
  n_re <- tmb_data$n_random_effects
  h <- obj$env$spHess(random = TRUE)

  sigma <- crossprod(b)
  for (i in seq_len(n_subjects)) {
    idx <- i + (seq_len(n_re) - 1L) * n_subjects
    hi <- as.matrix(h[idx, idx, drop = FALSE])
    vi <- tryCatch(
      solve(hi),
      error = function(e) MASS::ginv(hi)
    )
    sigma <- sigma + vi
  }
  sigma <- (sigma + t(sigma)) / (2 * n_subjects)

  ev <- eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) {
    sigma <- sigma + diag(abs(min(ev)) + 1e-8, n_re)
  }
  if (damping_rho < 1) {
    old_sigma <- as.matrix(obj$report()$Sigma_b)
    sigma <- (1 - damping_rho) * old_sigma + damping_rho * sigma
  }

  re <- .pack_correlation_theta(sigma, n_re)
  p$log_sd_re <- re$log_sd_re
  p$corr_par <- re$corr_par
  attr(p, "damping_rho") <- damping_rho
  p
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
  idx <- 1L
  if (parsed_long$biomarker$fixed) {
    x[idx] <- x[idx] + direction * 2 * log(scale)
    idx <- idx + 1L
  }
  if (parsed_long$velocity$fixed) {
    x[idx] <- x[idx] + direction * log(scale)
    idx <- idx + 1L
  }
  if (idx <= length(x)) {
    x[idx:length(x)] <- x[idx:length(x)] *
      if (to_internal) scale^2 else 1 / scale^2
  }
  x
}

#' @noRd
.random_effect_time_scale <- function(parsed_long, n_re, scale, to_internal) {
  mult <- rep(1, n_re)
  if (n_re >= 2) mult[2] <- scale
  idx <- 3L
  if (parsed_long$biomarker$random) idx <- idx + 1L
  if (parsed_long$velocity$random) idx <- idx + 1L
  if (idx <= n_re) mult[idx:n_re] <- scale^2
  if (to_internal) mult else 1 / mult
}

#' @noRd
.longitudinal_vcov_time_scale <- function(parsed_long, n_long, scale) {
  mult <- rep(1, n_long)
  idx <- 1L
  if (parsed_long$biomarker$fixed) idx <- idx + 1L
  if (parsed_long$velocity$fixed) idx <- idx + 1L
  if (idx <= n_long) mult[idx:n_long] <- 1 / scale^2
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
  re_pattern <- "\\(([^)]+\\|[^)]+)\\)"
  re_match <- regmatches(formula_str, regexpr(re_pattern, formula_str))

  if (length(re_match) == 0) {
    fixed_formula <- formula
    random_terms <- NULL
    grouping <- NULL
    biomarker_random <- FALSE
    velocity_random <- FALSE
  } else {
    inner <- substr(re_match, 2, nchar(re_match) - 1)
    parts <- strsplit(inner, "\\|")[[1]]
    grouping <- trimws(parts[2])
    re_terms_str <- trimws(parts[1])

    random_terms <- if (re_terms_str == "1") {
      "(Intercept)"
    } else {
      terms <- trimws(strsplit(re_terms_str, "\\+")[[1]])
      terms[terms == "1"] <- "(Intercept)"
      terms
    }

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
    grouping = grouping
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

  n_biomarker_velocity_fixed <- sum(
    parsed_long$biomarker$fixed, parsed_long$velocity$fixed
  )
  n_biomarker_velocity_random <- sum(
    parsed_long$biomarker$random, parsed_long$velocity$random
  )

  list(
    n_longitudinal_coef = n_longitudinal_fixed + n_biomarker_velocity_fixed,
    n_random_effects = n_longitudinal_random + n_biomarker_velocity_random + 2,
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

    # b1, b2 from longitudinal + RE
    b1 <- 0
    b2 <- 0
    fi <- 1
    ri <- 3
    if (configs$biomarker$fixed) {
      b1 <- -exp(lc[fi])
      fi <- fi + 1
    }
    if (configs$biomarker$random) {
      b1 <- b1 * exp(bi[ri])
      ri <- ri + 1
    }
    if (configs$velocity$fixed) {
      b2 <- -exp(lc[fi])
      fi <- fi + 1
    }
    if (configs$velocity$random) {
      b2 <- b2 * exp(bi[ri])
      ri <- ri + 1
    }
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
        forcing <- sum(lc[ffs:(ffs + ncol(fcov) - 1)] * fcov[ci, ])
        if (ncol(rcov) > 0) {
          forcing <- forcing + sum(bi[rfs:(rfs + ncol(rcov) - 1)] * rcov[ci, ])
        }

        mv <- .ode_step_r(m, v, b1, b2, forcing, dt)
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
  n_quad <- max(n_quad %||% 1L, 1L)
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
    surv_cov <- if (!is.null(subj$survival$covariates)) {
      as.numeric(subj$survival$covariates)
    } else {
      numeric(0)
    }

    bi <- re[i, ]
    m <- cf$initial_state[1] + bi[1]
    v <- cf$initial_state[2] + bi[2]

    # b1, b2 from longitudinal + RE
    b1 <- 0
    b2 <- 0
    fi <- 1
    ri <- 3
    if (configs$biomarker$fixed) {
      b1 <- -exp(lc[fi])
      fi <- fi + 1
    }
    if (configs$biomarker$random) {
      b1 <- b1 * exp(bi[ri])
      ri <- ri + 1
    }
    if (configs$velocity$fixed) {
      b2 <- -exp(lc[fi])
      fi <- fi + 1
    }
    if (configs$velocity$random) {
      b2 <- b2 * exp(bi[ri])
      ri <- ri + 1
    }
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
          forcing <- forcing + sum(bi[rfs:(rfs + ncol(rcov) - 1)] * rcov[ci, ])
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
