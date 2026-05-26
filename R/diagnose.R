#' Diagnose Subject-Level ODE Identifiability
#'
#' @description
#' Fits the second-order ODE separately for each subject and diagnoses
#' identifiability from the mixed posterior Hessian.
#'
#' @param formula Formula whose left-hand side is the longitudinal response.
#' @param data Data frame containing the response, time, and subject ID.
#' @param time Name of the time variable.
#' @param id Name of the subject ID variable.
#' @param hessian_condition Maximum allowed mixed posterior Hessian condition
#'   number. Defaults to \code{1e10}.
#' @param verbose Whether to print progress.
#'
#' @return A \code{JointODEDiagnose} object.
#' @export
diagnose <- function(
  formula,
  data,
  time = "time",
  id = "id",
  hessian_condition = 1e10,
  verbose = TRUE
) {
  cl <- match.call()
  if (!inherits(formula, "formula") || length(formula) < 3) {
    stop("formula must have a response on the left-hand side", call. = FALSE)
  }
  if (!is.numeric(hessian_condition) || length(hessian_condition) != 1 ||
    !is.finite(hessian_condition) || hessian_condition <= 0) {
    stop("hessian_condition must be a positive finite number", call. = FALSE)
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
  subject <- .diagnose_mixed_hessian(subject, id, hessian_condition)

  drop_ids <- subject[[id]][subject$status == "not_identifiable"]
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
      hessian_condition = hessian_condition,
      scales = list(time = time_scale),
      sigma_e = attr(subject, "sigma_e"),
      sigma_b = attr(subject, "sigma_b"),
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
  theta0 <- .diagnose_start(t, y)

  obj <- tryCatch(
    TMB::MakeADFun(
      data = list(model_type = 2L, time = t, y = y),
      parameters = list(theta = theta0),
      DLL = "JointODE",
      silent = TRUE
    ),
    error = identity
  )
  if (inherits(obj, "error")) {
    return(.diagnose_row(data[[id]][[1]], id, nrow(data), theta0))
  }

  opt <- tryCatch(stats::nlminb(obj$par, obj$fn, obj$gr), error = identity)
  if (inherits(opt, "error")) {
    return(.diagnose_row(data[[id]][[1]], id, nrow(data), theta0))
  }

  gradient <- tryCatch(obj$gr(opt$par), error = function(e) rep(NA_real_, 5))
  hessian <- tryCatch(obj$he(opt$par), error = function(e) NULL)
  fitted <- tryCatch(as.numeric(obj$report(opt$par)$fitted), error = function(e) {
    rep(NA_real_, length(y))
  })
  rmse <- sqrt(mean((y - fitted)^2, na.rm = TRUE))

  .diagnose_row(
    data[[id]][[1]], id, nrow(data), opt$par,
    convergence = opt$convergence,
    objective = opt$objective,
    gradient = gradient,
    hessian = hessian,
    rmse = rmse
  )
}

#' @noRd
.diagnose_start <- function(t, y) {
  dt <- diff(t)
  dy <- diff(y)
  slope <- if (any(dt > 0)) dy[which(dt > 0)[1]] / dt[which(dt > 0)[1]] else 0
  c(
    y[[1]],
    if (is.finite(slope)) slope else 0,
    0,
    log(0.8),
    mean(y, na.rm = TRUE)
  )
}

#' @noRd
.diagnose_row <- function(
  subject_id, id, n_obs, par,
  convergence = NA_integer_, objective = NA_real_,
  gradient = rep(NA_real_, 5), hessian = NULL, rmse = NA_real_
) {
  parameter_names <- .diagnose_parameter_names()
  hessian <- if (is.null(hessian)) matrix(NA_real_, 5, 5) else hessian
  hessian <- (hessian + t(hessian)) / 2
  metrics <- .diagnose_hessian_metrics(hessian, parameter_names)

  out <- data.frame(
    n_obs = n_obs,
    convergence = convergence,
    objective = objective,
    rmse = rmse,
    max_abs_gradient = suppressWarnings(max(abs(gradient), na.rm = TRUE)),
    m0 = par[[1]],
    v0 = par[[2]],
    log_omega2 = par[[3]],
    log_2xi_omega = par[[4]],
    forcing = par[[5]],
    omega = sqrt(exp(par[[3]])),
    xi = exp(par[[4]]) / (2 * sqrt(exp(par[[3]]))),
    raw_min_eigen = metrics$min_eigen,
    raw_condition_number = metrics$condition_number,
    raw_weakest_direction = metrics$weakest_direction,
    raw_loading_m0 = metrics$loading[[1]],
    raw_loading_v0 = metrics$loading[[2]],
    raw_loading_log_omega2 = metrics$loading[[3]],
    raw_loading_log_2xi_omega = metrics$loading[[4]],
    raw_loading_forcing = metrics$loading[[5]],
    min_eigen = metrics$min_eigen,
    condition_number = metrics$condition_number,
    weakest_direction = metrics$weakest_direction,
    loading_m0 = metrics$loading[[1]],
    loading_v0 = metrics$loading[[2]],
    loading_log_omega2 = metrics$loading[[3]],
    loading_log_2xi_omega = metrics$loading[[4]],
    loading_forcing = metrics$loading[[5]],
    status = "not_identifiable",
    reason = "hessian_unavailable",
    stringsAsFactors = FALSE
  )
  out[.diagnose_hessian_names(parameter_names)] <-
    as.list(as.numeric(hessian[lower.tri(hessian, diag = TRUE)]))
  out[[id]] <- subject_id
  out[c(id, setdiff(names(out), id))]
}

#' @noRd
.diagnose_mixed_hessian <- function(subject, id, hessian_condition) {
  parameter_names <- .diagnose_parameter_names()
  ok <- is.finite(subject$objective) &
    is.finite(subject$max_abs_gradient) &
    Reduce(`&`, lapply(parameter_names, function(nm) {
      is.finite(subject[[nm]])
    }))

  if (sum(ok) < length(parameter_names) + 2L) {
    attr(subject, "sigma_e") <- NA_real_
    attr(subject, "sigma_b") <- matrix(NA_real_, 5, 5)
    return(subject)
  }

  theta <- as.matrix(subject[ok, parameter_names, drop = FALSE])
  sigma_b <- stats::cov(theta) + diag(1e-6, length(parameter_names))
  precision_b <- solve(sigma_b)
  sigma_e <- sqrt(sum(subject$objective[ok]) / sum(subject$n_obs[ok]))

  for (i in seq_len(nrow(subject))) {
    if (!isTRUE(ok[[i]])) next
    raw_hessian <- .diagnose_subject_hessian(subject[i, ], parameter_names)
    if (any(!is.finite(raw_hessian))) next

    hessian <- raw_hessian / (2 * sigma_e^2) + precision_b
    metrics <- .diagnose_hessian_metrics(hessian, parameter_names)
    reason <- .diagnose_hessian_reason(
      metrics$min_eigen, metrics$condition_number, hessian_condition
    )

    subject$min_eigen[[i]] <- metrics$min_eigen
    subject$condition_number[[i]] <- metrics$condition_number
    subject$weakest_direction[[i]] <- metrics$weakest_direction
    subject$loading_m0[[i]] <- metrics$loading[[1]]
    subject$loading_v0[[i]] <- metrics$loading[[2]]
    subject$loading_log_omega2[[i]] <- metrics$loading[[3]]
    subject$loading_log_2xi_omega[[i]] <- metrics$loading[[4]]
    subject$loading_forcing[[i]] <- metrics$loading[[5]]
    subject$reason[[i]] <- reason
    subject$status[[i]] <- if (reason == "hessian_positive_definite") {
      "identifiable"
    } else {
      "not_identifiable"
    }
  }

  attr(subject, "sigma_e") <- sigma_e
  attr(subject, "sigma_b") <- sigma_b
  subject
}

#' @noRd
.diagnose_parameter_names <- function() {
  c("m0", "v0", "log_omega2", "log_2xi_omega", "forcing")
}

#' @noRd
.diagnose_hessian_names <- function(parameter_names) {
  hessian_names <- outer(parameter_names, parameter_names, paste, sep = "__")
  paste0("raw_hessian_", hessian_names[lower.tri(hessian_names, diag = TRUE)])
}

#' @noRd
.diagnose_subject_hessian <- function(row, parameter_names) {
  values <- unlist(row[.diagnose_hessian_names(parameter_names)], use.names = FALSE)
  if (any(!is.finite(values))) return(matrix(NA_real_, 5, 5))
  hessian <- matrix(0, 5, 5, dimnames = list(parameter_names, parameter_names))
  hessian[lower.tri(hessian, diag = TRUE)] <- values
  hessian[upper.tri(hessian)] <- t(hessian)[upper.tri(hessian)]
  hessian
}

#' @noRd
.diagnose_hessian_metrics <- function(hessian, parameter_names) {
  eig <- tryCatch(eigen(hessian, symmetric = TRUE), error = function(e) NULL)
  values <- if (is.null(eig)) rep(NA_real_, length(parameter_names)) else eig$values
  min_eigen <- suppressWarnings(min(values, na.rm = TRUE))
  condition <- if (all(is.finite(values)) && min(abs(values)) > 0) {
    max(abs(values)) / min(abs(values))
  } else {
    Inf
  }
  loading <- if (is.null(eig) || any(!is.finite(values))) {
    rep(NA_real_, length(parameter_names))
  } else {
    eig$vectors[, which.min(values)]
  }
  weakest_direction <- if (all(is.na(loading))) {
    NA_character_
  } else {
    parameter_names[which.max(abs(loading))]
  }
  list(
    min_eigen = min_eigen,
    condition_number = condition,
    weakest_direction = weakest_direction,
    loading = loading
  )
}

#' @noRd
.diagnose_hessian_reason <- function(min_eigen, condition, hessian_condition) {
  if (!is.finite(min_eigen) || !is.finite(condition)) {
    return("hessian_nonfinite")
  }
  if (min_eigen <= 0) {
    return("hessian_not_positive_definite")
  }
  if (condition > hessian_condition) {
    return("hessian_ill_conditioned")
  }
  "hessian_positive_definite"
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
  reasons <- sort(table(x$subject$reason), decreasing = TRUE)
  cat("\nHessian reasons:\n")
  print(reasons)
  invisible(x)
}
