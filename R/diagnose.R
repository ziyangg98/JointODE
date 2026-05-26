#' Diagnose Subject-Level ODE Identifiability
#'
#' @description
#' Fits a tau-lambda ODE separately for each subject and diagnoses
#' identifiability from the subject-level AD Hessian.
#'
#' @param formula Formula whose left-hand side is the longitudinal response.
#'   The right-hand side is ignored by this diagnostic.
#' @param data Data frame containing the response, time, and subject ID.
#' @param time Name of the time variable.
#' @param id Name of the subject ID variable.
#' @param verbose Whether to print progress.
#'
#' @return A \code{JointODEDiagnose} object.
#' @export
diagnose <- function(
  formula,
  data,
  time = "time",
  id = "id",
  verbose = TRUE
) {
  cl <- match.call()
  hessian_condition <- 1e8
  if (!inherits(formula, "formula") || length(formula) < 3) {
    stop("formula must have a response on the left-hand side", call. = FALSE)
  }
  response <- as.character(formula[[2]])
  required <- c(response, time, id)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      sprintf("Missing required columns: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  data <- data[!is.na(data[[id]]) &
    is.finite(data[[time]]) &
    is.finite(data[[response]]), , drop = FALSE]
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
      hessian_condition = hessian_condition
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
      sigma_e = attr(subject, "sigma_e"),
      sigma_b = attr(subject, "sigma_b"),
      call = cl
    ),
    class = "JointODEDiagnose"
  )
}

#' @noRd
.diagnose_subject <- function(data, response, time, id, hessian_condition) {
  data <- data[order(data[[time]]), , drop = FALSE]
  t <- data[[time]] - min(data[[time]])
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
    return(.diagnose_row(
      data[[id]][[1]], id, nrow(data), rep(NA_real_, 5),
      convergence = NA_integer_, objective = NA_real_,
      gradient = rep(NA_real_, 5), hessian = NULL,
      reason = "hessian_unavailable", hessian_condition = hessian_condition
    ))
  }

  opt <- tryCatch(
    stats::nlminb(obj$par, obj$fn, obj$gr),
    error = identity
  )
  if (inherits(opt, "error")) {
    return(.diagnose_row(
      data[[id]][[1]], id, nrow(data), rep(NA_real_, 5),
      convergence = NA_integer_, objective = NA_real_,
      gradient = rep(NA_real_, 5), hessian = NULL,
      reason = "hessian_unavailable", hessian_condition = hessian_condition
    ))
  }

  gradient <- tryCatch(obj$gr(opt$par), error = function(e) rep(NA_real_, 5))
  hessian <- tryCatch(obj$he(opt$par), error = function(e) NULL)
  .diagnose_row(
    data[[id]][[1]], id, nrow(data), opt$par,
    convergence = opt$convergence, objective = opt$objective,
    gradient = gradient, hessian = hessian,
    hessian_condition = hessian_condition
  )
}

#' @noRd
.diagnose_start <- function(t, y) {
  dt <- diff(t)
  dy <- diff(y)
  slope <- if (any(dt > 0)) dy[which(dt > 0)[1]] / dt[which(dt > 0)[1]] else 0
  lambda <- 0.05
  c(
    y[[1]],
    if (is.finite(slope)) slope else 0,
    log(lambda),
    log(5),
    lambda * mean(y, na.rm = TRUE)
  )
}

#' @noRd
.diagnose_row <- function(
  subject_id, id, n_obs, par, convergence, objective,
  gradient, hessian, reason = NULL, hessian_condition
) {
  hessian <- if (is.null(hessian)) matrix(NA_real_, 5, 5) else hessian
  hessian <- (hessian + t(hessian)) / 2
  eig <- tryCatch(eigen(hessian, symmetric = TRUE), error = function(e) NULL)
  values <- if (is.null(eig)) rep(NA_real_, 5) else eig$values
  min_eigen <- suppressWarnings(min(values, na.rm = TRUE))
  condition <- if (all(is.finite(values)) && min(abs(values)) > 0) {
    max(abs(values)) / min(abs(values))
  } else {
    Inf
  }
  weakest <- if (is.null(eig) || any(!is.finite(values))) {
    rep(NA_real_, 5)
  } else {
    eig$vectors[, which.min(values)]
  }
  parameter_names <- c("m0", "v0", "log_lambda", "log_tau", "eta")
  weakest_direction <- if (all(is.na(weakest))) {
    NA_character_
  } else {
    parameter_names[which.max(abs(weakest))]
  }

  if (is.null(reason)) {
    reason <- .diagnose_hessian_reason(
      min_eigen, condition, hessian_condition
    )
  }
  status <- if (reason == "hessian_positive_definite") {
    "identifiable"
  } else {
    "not_identifiable"
  }

  out <- data.frame(
    n_obs = n_obs,
    convergence = convergence,
    objective = objective,
    max_abs_gradient = suppressWarnings(max(abs(gradient), na.rm = TRUE)),
    m0 = par[[1]],
    v0 = par[[2]],
    lambda = exp(par[[3]]),
    tau = exp(par[[4]]),
    eta = par[[5]],
    raw_min_eigen = min_eigen,
    raw_condition_number = condition,
    raw_weakest_direction = weakest_direction,
    raw_loading_m0 = weakest[[1]],
    raw_loading_v0 = weakest[[2]],
    raw_loading_log_lambda = weakest[[3]],
    raw_loading_log_tau = weakest[[4]],
    raw_loading_eta = weakest[[5]],
    min_eigen = min_eigen,
    condition_number = condition,
    weakest_direction = weakest_direction,
    loading_m0 = weakest[[1]],
    loading_v0 = weakest[[2]],
    loading_log_lambda = weakest[[3]],
    loading_log_tau = weakest[[4]],
    loading_eta = weakest[[5]],
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
  hessian_names <- outer(parameter_names, parameter_names, paste, sep = "__")
  hessian_names <- paste0("raw_hessian_", hessian_names[lower.tri(hessian_names, diag = TRUE)])
  hessian_values <- hessian[lower.tri(hessian, diag = TRUE)]
  out[hessian_names] <- as.list(as.numeric(hessian_values))
  out[[id]] <- subject_id
  out[c(id, setdiff(names(out), id))]
}

#' @noRd
.diagnose_mixed_hessian <- function(subject, id, hessian_condition) {
  parameter_names <- c("m0", "v0", "log_lambda", "log_tau", "eta")
  ok <- subject$convergence == 0 &
    is.finite(subject$objective) &
    Reduce(`&`, lapply(
      c("m0", "v0", "lambda", "tau", "eta"),
      function(nm) is.finite(subject[[nm]])
    ))

  if (sum(ok) < length(parameter_names) + 2L) {
    attr(subject, "sigma_e") <- NA_real_
    attr(subject, "sigma_b") <- matrix(NA_real_, 5, 5)
    return(subject)
  }

  theta <- cbind(
    m0 = subject$m0,
    v0 = subject$v0,
    log_lambda = log(subject$lambda),
    log_tau = log(subject$tau),
    eta = subject$eta
  )
  sigma_b <- stats::cov(theta[ok, , drop = FALSE])
  sigma_b <- sigma_b + diag(1e-6, ncol(sigma_b))
  precision_b <- solve(sigma_b)
  sigma_e <- sqrt(sum(subject$objective[ok]) / sum(subject$n_obs[ok]))

  for (i in seq_len(nrow(subject))) {
    if (!isTRUE(ok[[i]])) next
    raw_hessian <- .diagnose_subject_hessian(subject[i, ])
    if (any(!is.finite(raw_hessian))) next

    hessian <- raw_hessian / (2 * sigma_e^2) + precision_b
    eig <- tryCatch(eigen(hessian, symmetric = TRUE), error = function(e) NULL)
    if (is.null(eig)) next

    values <- eig$values
    min_eigen <- suppressWarnings(min(values, na.rm = TRUE))
    condition <- if (all(is.finite(values)) && min(abs(values)) > 0) {
      max(abs(values)) / min(abs(values))
    } else {
      Inf
    }
    weakest <- if (any(!is.finite(values))) {
      rep(NA_real_, 5)
    } else {
      eig$vectors[, which.min(values)]
    }
    weakest_direction <- if (all(is.na(weakest))) {
      NA_character_
    } else {
      parameter_names[which.max(abs(weakest))]
    }

    reason <- .diagnose_hessian_reason(
      min_eigen, condition, hessian_condition
    )
    subject$min_eigen[[i]] <- min_eigen
    subject$condition_number[[i]] <- condition
    subject$weakest_direction[[i]] <- weakest_direction
    subject$loading_m0[[i]] <- weakest[[1]]
    subject$loading_v0[[i]] <- weakest[[2]]
    subject$loading_log_lambda[[i]] <- weakest[[3]]
    subject$loading_log_tau[[i]] <- weakest[[4]]
    subject$loading_eta[[i]] <- weakest[[5]]
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
.diagnose_subject_hessian <- function(row) {
  parameter_names <- c("m0", "v0", "log_lambda", "log_tau", "eta")
  hessian_names <- outer(parameter_names, parameter_names, paste, sep = "__")
  hessian_names <- paste0("raw_hessian_", hessian_names[lower.tri(hessian_names, diag = TRUE)])
  values <- unlist(row[hessian_names], use.names = FALSE)
  if (any(!is.finite(values))) return(matrix(NA_real_, 5, 5))

  hessian <- matrix(0, 5, 5, dimnames = list(parameter_names, parameter_names))
  hessian[lower.tri(hessian, diag = TRUE)] <- values
  hessian[upper.tri(hessian)] <- t(hessian)[upper.tri(hessian)]
  hessian
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
