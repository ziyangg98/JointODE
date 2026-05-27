#' Diagnose Subject-Level ODE Identifiability
#'
#' @description
#' Fits the second-order ODE separately for each subject and diagnoses
#' identifiability from the subject-level dynamic information ratio.
#'
#' @param formula Formula whose left-hand side is the longitudinal response.
#' @param data Data frame containing the response, time, and subject ID.
#' @param time Name of the time variable.
#' @param id Name of the subject ID variable.
#' @param q_threshold Optional minimum subject-level dynamic information ratio.
#'   Use a positive number or \code{"elbow"}. Filtering uses only
#'   \code{q_i >= q_threshold}; Hessian diagnostics are still reported but do
#'   not determine \code{keep_ids}.
#' @param verbose Whether to print progress.
#'
#' @return A \code{JointODEDiagnose} object.
#' @export
diagnose <- function(
  formula,
  data,
  time = "time",
  id = "id",
  q_threshold = NULL,
  verbose = TRUE
) {
  cl <- match.call()
  if (!inherits(formula, "formula") || length(formula) < 3) {
    stop("formula must have a response on the left-hand side", call. = FALSE)
  }
  valid_threshold <- is.null(q_threshold) ||
    identical(q_threshold, "elbow") ||
    (is.numeric(q_threshold) && length(q_threshold) == 1 &&
      is.finite(q_threshold) && q_threshold > 0)
  if (!valid_threshold) {
    stop("q_threshold must be a positive finite number or 'elbow'", call. = FALSE)
  }
  response <- as.character(formula[[2]])
  missing <- setdiff(c(response, time, id), names(data))
  if (length(missing) > 0) {
    stop(
      sprintf("Missing required columns: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  data <- data[!is.na(data[[id]]) &
    is.finite(data[[time]]) &
    is.finite(data[[response]]), , drop = FALSE]
  parsed_long <- .parse_longitudinal_formula(formula)
  if (is.null(parsed_long$grouping)) parsed_long$grouping <- id
  time_scale <- .estimate_time_scale(
    .process_marginal(formula, data, time, id, parsed_long)
  )
  ids <- unique(data[[id]])

  subject <- do.call(rbind, lapply(seq_along(ids), function(i) {
    if (isTRUE(verbose) && (i == 1 || i %% 50 == 0 || i == length(ids))) {
      cli::cli_alert_info("Diagnosing subject {i}/{length(ids)}")
    }
    .diagnose_subject(
      data[data[[id]] == ids[[i]], , drop = FALSE],
      response = response,
      time = time,
      id = id,
      time_scale = time_scale
    )
  }))
  q_threshold_method <- if (identical(q_threshold, "elbow")) "elbow" else "fixed"
  if (identical(q_threshold, "elbow")) {
    log_q <- sort(log10(subject$q_i[is.finite(subject$q_i) & subject$q_i > 0]))
    if (length(log_q) < 3 || diff(range(log_q)) == 0) {
      stop("Cannot estimate elbow q_threshold from fewer than 3 finite positive q_i values.", call. = FALSE)
    }
    x <- seq(0, 1, length.out = length(log_q))
    y <- (log_q - min(log_q)) / diff(range(log_q))
    q_threshold <- 10^log_q[which.max(abs(y - x))]
  }
  if (is.null(q_threshold)) {
    subject$status <- "identifiable"
    subject$reason <- "q_threshold_not_supplied"
    subject$mechanism <- "not_filtered"
    subject$explanation <- "No q_threshold was supplied; all subjects are retained."
    drop_ids <- subject[[id]][FALSE]
  } else {
    keep_q <- is.finite(subject$q_i) & subject$q_i >= q_threshold
    subject$status <- ifelse(keep_q, "identifiable", "not_identifiable")
    subject$reason <- ifelse(keep_q, "q_threshold_passed", "low_dynamic_information")
    subject$mechanism <- "q_threshold_passed"
    subject$mechanism[!keep_q & !is.finite(subject$q_i)] <- "optimization_failed"
    subject$mechanism[!keep_q & is.finite(subject$q_i)] <- "weak_dynamic_information"
    subject$mechanism[!keep_q & is.finite(subject$omega) & subject$omega < 0.05] <-
      "omega_near_zero"
    xi_extreme <- !keep_q & is.finite(subject$xi) &
      (subject$xi < 0.05 | subject$xi > 5)
    subject$mechanism[xi_extreme] <- "xi_extreme"
    subject$mechanism[!keep_q & is.finite(subject$omega) & subject$omega < 0.05 &
      xi_extreme] <- "omega_near_zero_and_xi_extreme"
    subject$explanation <- ifelse(
      keep_q,
      "Subject retained: dynamic information ratio is above q_threshold.",
      paste0(
        "Subject removed: q_i is below q_threshold; single-subject estimates ",
        "indicate weak or unstable second-order dynamic information."
      )
    )
    subject$explanation[subject$mechanism == "optimization_failed"] <-
      "Subject removed: single-subject diagnostic did not produce a finite q_i."
    subject$explanation[subject$mechanism == "omega_near_zero"] <-
      "Subject removed: estimated omega is close to zero, so restoring-force dynamics are weakly identified."
    subject$explanation[subject$mechanism == "xi_extreme"] <-
      "Subject removed: estimated damping ratio is near zero or extremely large, making damping dynamics weakly identified."
    subject$explanation[subject$mechanism == "omega_near_zero_and_xi_extreme"] <-
      "Subject removed: estimated omega is close to zero and damping is extreme, so second-order dynamics are weakly identified."
    drop_ids <- subject[[id]][!keep_q]
  }
  keep_ids <- setdiff(ids, drop_ids)
  summary <- data.frame(
    n_subjects = length(ids),
    n_identifiable = sum(subject$status == "identifiable"),
    n_not_identifiable = sum(subject$status == "not_identifiable"),
    pct_identifiable = mean(subject$status == "identifiable"),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      subject = subject,
      summary = summary,
      filtered_data = data[data[[id]] %in% keep_ids, , drop = FALSE],
      keep_ids = keep_ids,
      drop_ids = drop_ids,
      q_threshold = q_threshold,
      q_threshold_method = q_threshold_method,
      scales = list(time = time_scale),
      call = cl
    ),
    class = "JointODEDiagnose"
  )
}

#' @noRd
.diagnose_subject <- function(data, response, time, id, time_scale) {
  data <- data[order(data[[time]]), , drop = FALSE]
  t <- data[[time]] / time_scale
  y <- data[[response]]
  parameter_names <- c("m0", "log_omega2", "log_2xi_omega", "forcing")
  theta0 <- c(y[[1]], 0, log(0.8), mean(y, na.rm = TRUE))

  obj <- TMB::MakeADFun(
    data = list(model_type = 2L, time = t, y = y),
    parameters = list(theta = theta0),
    DLL = "JointODE",
    silent = TRUE
  )

  opt <- stats::nlminb(obj$par, obj$fn, obj$gr)

  gradient <- obj$gr(opt$par)
  hessian <- obj$he(opt$par)
  hessian <- (hessian + t(hessian)) / 2
  if (any(!is.finite(hessian))) {
    values <- rep(NA_real_, length(parameter_names))
    loading <- rep(NA_real_, length(parameter_names))
  } else {
    eig <- eigen(hessian, symmetric = TRUE)
    values <- eig$values
    loading <- eig$vectors[, which.min(values)]
  }
  min_eigen <- suppressWarnings(min(values, na.rm = TRUE))
  condition_number <- if (all(is.finite(values)) && min(abs(values)) > 0) {
    max(abs(values)) / min(abs(values))
  } else {
    Inf
  }
  weakest_direction <- if (all(is.na(loading))) {
    NA_character_
  } else {
    parameter_names[which.max(abs(loading))]
  }
  hessian_reason <- if (!is.finite(min_eigen)) {
    "hessian_nonfinite"
  } else if (min_eigen <= 0) {
    "hessian_not_positive_definite"
  } else {
    "hessian_positive_definite"
  }

  fitted <- as.numeric(obj$report(opt$par)$fitted)
  rmse <- sqrt(mean((y - fitted)^2, na.rm = TRUE))
  q_i <- NA_real_
  if (length(t) >= length(theta0) && all(is.finite(opt$par))) {
    fitted_at <- function(par) {
      m <- par[[1]]
      b1 <- -exp(par[[2]])
      b2 <- -exp(par[[3]])
      forcing <- par[[4]]
      v <- (forcing + b1 * m) / (-b2)
      out <- numeric(length(t))
      previous <- 0
      for (j in seq_along(t)) {
        z <- .ode_step_r(m, v, b1, b2, forcing, t[[j]] - previous)
        m <- z[[1]]
        v <- z[[2]]
        out[[j]] <- m
        previous <- t[[j]]
      }
      out
    }
    jac <- matrix(NA_real_, length(t), length(opt$par))
    for (j in seq_along(opt$par)) {
      h <- 1e-4 * max(1, abs(opt$par[[j]]))
      plus <- minus <- opt$par
      plus[[j]] <- plus[[j]] + h
      minus[[j]] <- minus[[j]] - h
      jac[, j] <- (fitted_at(plus) - fitted_at(minus)) / (2 * h)
    }
    eig <- eigen(crossprod(jac), symmetric = TRUE, only.values = TRUE)$values
    tr <- sum(eig)
    if (is.finite(tr) && tr > 0) q_i <- min(eig) / tr
  }

  par_out <- opt$par
  par_out[[2]] <- par_out[[2]] - 2 * log(time_scale)
  par_out[[3]] <- par_out[[3]] - log(time_scale)
  par_out[[4]] <- par_out[[4]] / time_scale^2

  out <- data.frame(
    n_obs = nrow(data),
    convergence = opt$convergence,
    objective = opt$objective,
    rmse = rmse,
    max_abs_gradient = suppressWarnings(max(abs(gradient), na.rm = TRUE)),
    q_i = q_i,
    m0 = par_out[[1]],
    log_omega2 = par_out[[2]],
    log_2xi_omega = par_out[[3]],
    forcing = par_out[[4]],
    initial_velocity = (par_out[[4]] - exp(par_out[[2]]) * par_out[[1]]) /
      exp(par_out[[3]]),
    omega = sqrt(exp(par_out[[2]])),
    xi = exp(par_out[[3]]) / (2 * sqrt(exp(par_out[[2]]))),
    min_eigen = min_eigen,
    condition_number = condition_number,
    weakest_direction = weakest_direction,
    loading_m0 = loading[[1]],
    loading_log_omega2 = loading[[2]],
    loading_log_2xi_omega = loading[[3]],
    loading_forcing = loading[[4]],
    status = "not_identifiable",
    reason = "hessian_unavailable",
    hessian_reason = hessian_reason,
    stringsAsFactors = FALSE
  )
  out[[id]] <- data[[id]][[1]]
  out[c(id, setdiff(names(out), id))]
}

#' @export
print.JointODEDiagnose <- function(
  x, digits = max(3L, getOption("digits") - 3L), ...
) {
  cat("\nJointODE subject identifiability diagnostic\n")
  cat("Call: ")
  print(x$call)
  cat("\n")
  print(x$summary, digits = digits, row.names = FALSE)
  if (!is.null(x$q_threshold)) {
    cat(
      "\nFiltering: q_i >= ",
      format(x$q_threshold, digits = digits),
      " (",
      x$q_threshold_method,
      ")",
      "\n",
      sep = ""
    )
  }
  reasons <- sort(table(x$subject$reason), decreasing = TRUE)
  cat("\nFiltering reasons:\n")
  print(reasons)
  hessian_reasons <- sort(table(x$subject$hessian_reason), decreasing = TRUE)
  cat("\nHessian diagnostics:\n")
  print(hessian_reasons)
  if (!is.null(x$q_threshold) && length(x$drop_ids) > 0) {
    mechanisms <- sort(table(x$subject$mechanism[x$subject$status == "not_identifiable"]),
      decreasing = TRUE
    )
    cat("\nRemoval mechanisms:\n")
    print(mechanisms)
    examples <- x$subject[x$subject$status == "not_identifiable", , drop = FALSE]
    examples <- examples[order(examples$q_i, na.last = TRUE), , drop = FALSE]
    examples <- head(examples, 5)
    id_name <- names(x$subject)[1]
    cat("\nRemoved subject examples:\n")
    print(
      examples[, c(id_name, "n_obs", "q_i", "omega", "xi", "mechanism", "explanation"),
        drop = FALSE
      ],
      row.names = FALSE
    )
  }
  invisible(x)
}
