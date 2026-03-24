#' @importFrom stats nlm
NULL

# Package-level cache for quadrature grids, keyed by "dim_level"
.quad_cache <- new.env(parent = emptyenv())

.get_quad_grid <- function(dim, level) {
  key <- paste0(dim, "_", level)
  if (!exists(key, envir = .quad_cache)) {
    grid <- mvQuad::createNIGrid(
      dim = dim, type = "GHe",
      level = level, ndConstruction = "product"
    )
    assign(key, list(
      nodes = mvQuad::getNodes(grid),
      weights = as.numeric(mvQuad::getWeights(grid))
    ), envir = .quad_cache)
  }
  get(key, envir = .quad_cache)
}

.expand_posteriors <- function(data_list, posteriors) {
  n_subjects <- length(data_list)
  n_nodes_per_subject <- vapply(posteriors$nodes, nrow, integer(1))
  total_nodes <- sum(n_nodes_per_subject)

  all_nodes <- do.call(rbind, posteriors$nodes)
  all_weights <- unlist(posteriors$weights, use.names = FALSE)
  all_data <- vector("list", total_nodes)

  idx <- 1L
  for (i in seq_len(n_subjects)) {
    n <- n_nodes_per_subject[i]
    all_data[idx:(idx + n - 1L)] <- rep(list(data_list[[i]]), n)
    idx <- idx + n
  }

  list(
    nodes = all_nodes, data = all_data, weights = all_weights,
    node_to_subject = rep(seq_len(n_subjects), n_nodes_per_subject)
  )
}

.update_measurement_error_sd <- function(data_list, parameters, posteriors) {
  n_total_obs <- sum(vapply(data_list, function(d) {
    length(d$longitudinal$measurements)
  }, integer(1)))

  expanded <- .expand_posteriors(data_list, posteriors)

  ode_solutions <- .solve_batch_ode_cppad(
    data_list = expanded$data,
    random_effects = expanded$nodes,
    parameters = parameters
  )

  sigma_e_squared <- sum(vapply(seq_along(ode_solutions), function(k) {
    i <- expanded$node_to_subject[k]
    measurements <- data_list[[i]]$longitudinal$measurements
    if (length(measurements) == 0) return(0)

    obs_times <- data_list[[i]]$longitudinal$times
    result_times <- ode_solutions[[k]]$times
    match_idx <- vapply(obs_times, function(t) {
      which.min(abs(result_times - t))
    }, integer(1))
    biomarker_pred <- ode_solutions[[k]]$biomarker[match_idx]

    expanded$weights[k] * sum((measurements - biomarker_pred)^2)
  }, numeric(1)))

  sqrt(sigma_e_squared / n_total_obs)
}


.compute_objective_expected <- function(
  params, data_list, posteriors, parameters,
  gradient = TRUE, hessian = FALSE,
  parallel = FALSE, n_cores = 0
) {
  n_subjects <- length(data_list)

  if (parallel && n_subjects > 1) {
    if (n_cores == 0) {
      n_cores <- max(1L, parallel::detectCores() - 1L)
    }
    n_chunks <- min(n_cores, n_subjects)
    chunks <- split(seq_len(n_subjects), cut(seq_len(n_subjects),
      n_chunks, labels = FALSE))

    compute_chunk <- function(subject_indices) {
      chunk_data <- data_list[subject_indices]
      chunk_posteriors <- list(
        nodes = posteriors$nodes[subject_indices],
        weights = posteriors$weights[subject_indices]
      )
      expanded <- .expand_posteriors(chunk_data, chunk_posteriors)
      .compute_objective_cppad(
        params = params,
        data_list = expanded$data,
        random_effects = expanded$nodes,
        parameters = parameters,
        weights = expanded$weights,
        gradient = gradient,
        hessian = hessian
      )
    }

    results <- .parallel_apply(
      chunks, compute_chunk,
      parallel = TRUE, n_cores = n_cores, setup = FALSE
    )

    total_obj <- sum(vapply(results, as.numeric, numeric(1)))
    result <- structure(total_obj, names = NULL)

    if (gradient) {
      total_grad <- Reduce(`+`, lapply(results, function(r) {
        as.vector(attr(r, "gradient"))
      }))
      attr(result, "gradient") <- total_grad
    }
    if (hessian) {
      total_hess <- Reduce(`+`, lapply(results, function(r) {
        attr(r, "hessian")
      }))
      attr(result, "hessian") <- total_hess
    }
    result
  } else {
    expanded <- .expand_posteriors(data_list, posteriors)
    .compute_objective_cppad(
      params = params,
      data_list = expanded$data,
      random_effects = expanded$nodes,
      parameters = parameters,
      weights = expanded$weights,
      gradient = gradient,
      hessian = hessian
    )
  }
}


.compute_posterior_laplace <- function(
  data,
  random_effect,
  parameters,
  level = 3
) {
  objective <- function(theta) {
    result <- .compute_logpost_cppad(
      random_effect = theta, data = data,
      parameters = parameters, gradient = TRUE, hessian = TRUE
    )
    value <- -as.numeric(result)
    attr(value, "gradient") <- -as.vector(attr(result, "gradient"))
    attr(value, "hessian") <- -as.matrix(attr(result, "hessian"))
    value
  }

  fit <- nlm(
    f = objective, p = random_effect,
    hessian = FALSE, check.analyticals = FALSE
  )

  result_at_mode <- .compute_logpost_cppad(
    random_effect = fit$estimate, data = data,
    parameters = parameters, gradient = FALSE, hessian = TRUE
  )
  hessian_neglogpost <- -attr(result_at_mode, "hessian")

  # Gill-Murray: find minimal tau s.t. H + tau*I is positive definite
  n <- nrow(hessian_neglogpost)
  tau <- 0
  ds <- max(abs(diag(hessian_neglogpost)), 1, na.rm = TRUE)
  if (!is.finite(ds)) ds <- 1
  chol_factor <- NULL
  for (k in seq_len(20)) {
    R <- try(chol(hessian_neglogpost + diag(tau, n)), silent = TRUE)
    if (!inherits(R, "try-error")) {
      chol_factor <- t(backsolve(R, diag(n)))
      break
    }
    tau <- if (tau == 0) 1e-4 * ds else tau * 10
  }
  if (is.null(chol_factor)) chol_factor <- diag(1e-4, n)

  quad <- .get_quad_grid(length(fit$estimate), level)
  quad_nodes <- quad$nodes
  quad_weights <- quad$weights
  n_nodes <- nrow(quad_nodes)

  n_re <- length(fit$estimate)
  aghq_nodes <- t(vapply(seq_len(n_nodes), function(k) {
    fit$estimate + tcrossprod(quad_nodes[k, ], chol_factor)
  }, numeric(n_re)))

  logpost_at_nodes <- vapply(seq_len(n_nodes), function(k) {
    result <- .compute_logpost_cppad(
      random_effect = aghq_nodes[k, ],
      data = data,
      parameters = parameters,
      gradient = FALSE,
      hessian = FALSE
    )
    as.numeric(result)
  }, numeric(1))

  # AGHQ change-of-variables: |det(L)| Jacobian factor
  det_result <- determinant(hessian_neglogpost, logarithm = TRUE)
  log_jacobian <- -0.5 * as.numeric(det_result$modulus)

  log_weights <- log(quad_weights) + log_jacobian + logpost_at_nodes
  max_log_weight <- max(log_weights)
  aghq_weights <- exp(log_weights - max_log_weight)
  aghq_weights <- aghq_weights / sum(aghq_weights)

  list(
    nodes = aghq_nodes,
    weights = aghq_weights
  )
}

.compute_posterior_moments <- function(posterior_result) {
  nodes <- posterior_result$nodes
  weights <- posterior_result$weights

  mean <- colSums(weights * nodes)
  weighted_nodes <- sqrt(weights) * nodes
  second_moment <- crossprod(weighted_nodes)

  list(mean = mean, second_moment = second_moment)
}

.compute_posteriors <- function(
  data_list,
  parameters,
  random_effects,
  parallel = FALSE,
  n_cores = 0,
  level = 3,
  setup = TRUE
) {
  n_subjects <- length(data_list)

  # Warm the cache once before parallel dispatch
  .get_quad_grid(ncol(random_effects), level)

  compute_subject_posterior <- function(i) {
    .compute_posterior_laplace(
      data = data_list[[i]],
      random_effect = random_effects[i, ],
      parameters = parameters,
      level = level
    )
  }

  posterior_results <- .parallel_apply(
    seq_len(n_subjects),
    compute_subject_posterior,
    parallel = parallel,
    n_cores = n_cores,
    setup = setup
  )

  list(
    nodes = lapply(posterior_results, `[[`, "nodes"),
    weights = lapply(posterior_results, `[[`, "weights")
  )
}
