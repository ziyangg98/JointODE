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

#' Initialize RE columns 1-2 from first observations
#' @noRd
.init_re_from_observations <- function(data_list, m0, v0 = 0, n_re = 2L) {
  n_subjects <- length(data_list)
  re <- matrix(0, nrow = n_subjects, ncol = max(n_re, 2L))
  for (i in seq_along(data_list)) {
    obs <- data_list[[i]]$longitudinal
    if (length(obs$measurements) >= 1)
      re[i, 1] <- obs$measurements[1] - m0
    if (length(obs$measurements) >= 2) {
      dt <- obs$times[2] - obs$times[1]
      if (dt > 0)
        re[i, 2] <- (obs$measurements[2] - obs$measurements[1]) / dt - v0
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
  merged
}

# Constants ====================================================================

.default_spline <- list(
  degree = 2,
  n_knots = 1,
  knot_placement = "equal",
  boundary_knots = NULL
)

.reserved_words <- c("biomarker", "velocity")

.init_state_names <- c("init_biomarker", "init_velocity")

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

#' ODE step (double, 5-branch)
#' @noRd
.ode_step_r <- function(m, v, b1, b2, f, dt) {
  h <- b2 / 2
  D <- b2^2 + 4 * b1
  eps <- 1e-10

  if (abs(b1) < 1e-12 && abs(b2) < 1e-12) {
    a0 <- 1; a1 <- dt; J0 <- dt; J1 <- dt^2 / 2
  } else if (abs(b1) < 1e-12) {
    eb <- exp(b2 * dt); a0 <- 1
    a1 <- if (abs(b2) > eps) (eb - 1) / b2 else dt
    J1 <- if (abs(b2) > eps) (eb - 1) / b2^2 - dt / b2 else dt^2 / 2
    J0 <- if (abs(b2) > eps) (eb - 1) / b2 else dt
  } else if (abs(D) < 1e-10) {
    ert <- exp(h * dt)
    a0 <- ert * (1 - h * dt); a1 <- ert * dt
    if (abs(h) > eps) {
      ir <- (ert - 1) / h
      isr <- (ert * (h * dt - 1) + 1) / h^2
      J1 <- isr; J0 <- ir - h * isr
    } else { J0 <- dt; J1 <- dt^2 / 2 }
  } else if (D > 0) {
    sD <- sqrt(D)
    r1 <- h + sD / 2; r2 <- h - sD / 2
    r1s <- if (r1^2 > 1e-20) r1 else eps
    r2s <- if (r2^2 > 1e-20) r2 else -eps
    sDs <- max(sD, eps)
    e1 <- exp(r1 * dt); e2 <- exp(r2 * dt)
    a0 <- (r1 * e2 - r2 * e1) / sDs
    a1 <- (e1 - e2) / sDs
    i1 <- (e1 - 1) / r1s; i2 <- (e2 - 1) / r2s
    J1 <- (i1 - i2) / sDs; J0 <- (r1 * i2 - r2 * i1) / sDs
  } else {
    w <- sqrt(-D) / 2; ws <- max(w, eps)
    ms <- max(h^2 + w^2, 1e-20)
    er <- exp(h * dt)
    cw <- cos(w * dt); sw <- sin(w * dt)
    a0 <- er * (cw - h * sw / ws)
    a1 <- er * sw / ws
    ic <- (er * (h * cw + w * sw) - h) / ms
    is_ <- (er * (h * sw - w * cw) + w) / ms
    J1 <- is_ / ws; J0 <- ic - (h / ws) * is_
  }
  c(a0 * m + a1 * v + f * J1,
    b1 * a1 * m + (a0 + b2 * a1) * v + f * (J0 + b2 * J1))
}

#' Predict marginal trajectories for all subjects at given times
#' @noRd
.predict_marginal_trajectories <- function(data_list, times, cf, configs, re) {
  clamp_val <- configs$biomarker_clamp
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
    b1 <- 0; b2 <- 0; fi <- 1; ri <- 3
    if (configs$biomarker$fixed)  { b1 <- lc[fi]; fi <- fi + 1 }
    if (configs$biomarker$random) { b1 <- b1 + bi[ri]; ri <- ri + 1 }
    if (configs$velocity$fixed)   { b2 <- lc[fi]; fi <- fi + 1 }
    if (configs$velocity$random)  { b2 <- b2 + bi[ri]; ri <- ri + 1 }
    ffs <- fi; rfs <- ri

    bio <- vel <- numeric(n_times)
    prev_t <- 0

    for (ti in seq_len(n_times)) {
      tgt <- times[ti]; dt <- tgt - prev_t
      if (dt > 0 && length(obs_t) > 0) {
        ci <- findInterval(prev_t, obs_t)
        ci <- max(1, min(ci, nrow(fcov)))
        forcing <- sum(lc[ffs:(ffs + ncol(fcov) - 1)] * fcov[ci, ])
        if (ncol(rcov) > 0)
          forcing <- forcing + sum(bi[rfs:(rfs + ncol(rcov) - 1)] * rcov[ci, ])

        mv <- .ode_step_r(m, v, b1, b2, forcing, dt)
        m <- max(-clamp_val, min(clamp_val, mv[1]))
        v <- max(-clamp_val, min(clamp_val, mv[2]))
      }
      bio[ti] <- m; vel[ti] <- v; prev_t <- tgt
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
    t, knots = spline_config$knots,
    Boundary.knots = spline_config$boundary_knots,
    degree = spline_config$degree, intercept = TRUE
  )
  log_h <- sum(basis * baseline) + hazard[1] * m
  if (gamma == 1) log_h <- log_h + hazard[2] * v
  else if (gamma == 2) log_h <- log_h + hazard[2] * v^2
  if (length(surv_cov) > 0)
    log_h <- log_h + sum(hazard[-(1:2)] * surv_cov)
  exp(max(-20, min(20, log_h)))
}

#' Predict trajectories for all subjects at given times
#' @noRd
.predict_trajectories <- function(
  data_list, times, cf, configs, re, n_quad
) {
  n_quad <- max(n_quad %||% 1L, 1L)
  clamp_val <- configs$biomarker_clamp
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
    b1 <- 0; b2 <- 0; fi <- 1; ri <- 3
    if (configs$biomarker$fixed)  { b1 <- lc[fi]; fi <- fi + 1 }
    if (configs$biomarker$random) { b1 <- b1 + bi[ri]; ri <- ri + 1 }
    if (configs$velocity$fixed)   { b2 <- lc[fi]; fi <- fi + 1 }
    if (configs$velocity$random)  { b2 <- b2 + bi[ri]; ri <- ri + 1 }
    ffs <- fi; rfs <- ri

    bio <- vel <- cumhaz <- numeric(n_times)
    ch <- 0; prev_t <- 0

    for (ti in seq_len(n_times)) {
      tgt <- times[ti]; dt <- tgt - prev_t
      if (dt > 0 && length(obs_t) > 0) {
        ci <- findInterval(prev_t, obs_t)
        ci <- max(1, min(ci, nrow(fcov)))
        forcing <- sum(lc[ffs:(ffs + ncol(fcov) - 1)] * fcov[ci, ])
        if (ncol(rcov) > 0)
          forcing <- forcing + sum(bi[rfs:(rfs + ncol(rcov) - 1)] * rcov[ci, ])

        sub_dt <- dt / n_quad
        hl <- .eval_hazard_r(m, v, prev_t, bl, hz, surv_cov, sbc, configs$gamma)
        for (q in seq_len(n_quad)) {
          ts <- prev_t + (q - 1) * sub_dt
          # Midpoint
          mv <- .ode_step_r(m, v, b1, b2, forcing, sub_dt / 2)
          mm <- max(-clamp_val, min(clamp_val, mv[1]))
          vv <- max(-clamp_val, min(clamp_val, mv[2]))
          hm <- .eval_hazard_r(mm, vv, ts + sub_dt / 2, bl, hz, surv_cov, sbc, configs$gamma)
          # Right
          mv <- .ode_step_r(mm, vv, b1, b2, forcing, sub_dt / 2)
          mm <- max(-clamp_val, min(clamp_val, mv[1]))
          vv <- max(-clamp_val, min(clamp_val, mv[2]))
          hr <- .eval_hazard_r(mm, vv, ts + sub_dt, bl, hz, surv_cov, sbc, configs$gamma)
          ch <- ch + (sub_dt / 6) * (hl + 4 * hm + hr)
          m <- mm; v <- vv; hl <- hr
        }
      }
      bio[ti] <- m; vel[ti] <- v
      cumhaz[ti] <- ch; prev_t <- tgt
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
