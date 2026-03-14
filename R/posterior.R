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
  gradient = TRUE, hessian = FALSE
) {
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


.compute_posterior_laplace <- function(
  data,
  random_effect,
  parameters,
  level = 3
) {
  objective <- function(theta) {
    result <- .compute_logpost_cppad(
      random_effect = theta,
      data = data,
      parameters = parameters,
      gradient = TRUE,
      hessian = TRUE
    )
    value <- -as.numeric(result)
    attr(value, "gradient") <- -as.vector(attr(result, "gradient"))
    attr(value, "hessian") <- -as.matrix(attr(result, "hessian"))
    value
  }

  fit <- tryCatch(
    suppressWarnings(
      nlm(
        f = objective,
        p = random_effect,
        hessian = FALSE,
        check.analyticals = FALSE
      )
    ),
    error = function(e) {
      # nlm failed (e.g., non-finite value); fall back to initial estimate
      list(estimate = random_effect, minimum = NA_real_, code = 5)
    }
  )

  result_at_mode <- .compute_logpost_cppad(
    random_effect = fit$estimate,
    data = data,
    parameters = parameters,
    gradient = FALSE,
    hessian = TRUE
  )
  hessian_neglogpost <- -attr(result_at_mode, "hessian")

  # Compute Cholesky factor for AGHQ
  # If direct solve fails, try regularization to ensure positive definiteness
  chol_factor <- tryCatch(
    {
      t(chol(solve(hessian_neglogpost)))
    },
    error = function(e) {
      # Fallback 1: Add small diagonal regularization
      tryCatch(
        {
          diag_reg <- mean(diag(hessian_neglogpost)) * 1e-6
          regularized_hess <- hessian_neglogpost +
            diag(diag_reg, nrow(hessian_neglogpost))
          t(chol(solve(regularized_hess)))
        },
        error = function(e2) {
          # Fallback 2: Use eigenvalue decomposition for more stable inverse
          tryCatch(
            {
              eig <- eigen(hessian_neglogpost, symmetric = TRUE)
              # Threshold small eigenvalues
              eig$values[eig$values < max(eig$values) * 1e-10] <- max(
                eig$values
              ) *
                1e-10
              hess_inv <- eig$vectors %*%
                diag(1 / eig$values) %*%
                t(eig$vectors)
              t(chol(hess_inv))
            },
            error = function(e3) {
              # Last resort: use identity scaled by Hessian trace
              warning(
                "Cholesky decomposition failed, using diagonal approximation"
              )
              diag(1 / sqrt(diag(hessian_neglogpost)))
            }
          )
        }
      )
    }
  )

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

  # Compute Jacobian factor for coordinate transformation
  # When transforming from standard GHN points z to posterior space
  # theta = m + L*z, where L = chol(solve(H)), the weights need to
  # be multiplied by |det(L)|
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
