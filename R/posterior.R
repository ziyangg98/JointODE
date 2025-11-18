#' @importFrom stats nlm
NULL

.update_measurement_error_sd <- function(
  data_list,
  parameters,
  posteriors
) {
  n_subjects <- length(data_list)
  n_total_obs <- sum(sapply(data_list, function(d) {
    length(d$longitudinal$measurements)
  }))

  # Build data structures for batch ODE solving
  n_nodes_per_subject <- sapply(posteriors$nodes, nrow)
  total_nodes <- sum(n_nodes_per_subject)

  # Pre-allocate
  all_nodes <- do.call(rbind, posteriors$nodes)
  all_weights <- unlist(posteriors$weights, use.names = FALSE)
  all_data <- vector("list", total_nodes)

  idx <- 1
  for (i in seq_len(n_subjects)) {
    n_nodes <- n_nodes_per_subject[i]
    idx_range <- idx:(idx + n_nodes - 1)
    all_data[idx_range] <- rep(list(data_list[[i]]), n_nodes)
    idx <- idx + n_nodes
  }

  node_to_subject <- rep(seq_len(n_subjects), n_nodes_per_subject)

  ode_solutions <- .solve_batch_ode_cppad(
    data_list = all_data,
    random_effects = all_nodes,
    parameters = parameters
  )

  # Compute weighted residual sum of squares
  sigma_e_squared <- sum(sapply(seq_along(ode_solutions), function(k) {
    i <- node_to_subject[k]
    data_i <- data_list[[i]]
    measurements <- data_i$longitudinal$measurements

    if (length(measurements) == 0) {
      return(0)
    }

    obs_times <- data_i$longitudinal$times
    result_times <- ode_solutions[[k]]$times
    match_idx <- match(obs_times, result_times)
    biomarker_pred <- ode_solutions[[k]]$biomarker[match_idx]
    residuals <- measurements - biomarker_pred

    all_weights[k] * sum(residuals^2)
  }))

  # Clean up large objects
  rm(all_nodes, all_data, all_weights, ode_solutions, node_to_subject)

  sqrt(sigma_e_squared / n_total_obs)
}


.compute_objective_expected <- function(
  params,
  data_list,
  posteriors,
  parameters,
  gradient = TRUE,
  hessian = FALSE
) {
  n_subjects <- length(data_list)

  # Pre-allocate for efficiency
  total_nodes <- sum(sapply(posteriors$nodes, nrow))
  n_random_effects <- ncol(posteriors$nodes[[1]])

  all_nodes <- matrix(0, nrow = total_nodes, ncol = n_random_effects)
  all_data <- vector("list", total_nodes)
  all_weights <- numeric(total_nodes)

  idx <- 1
  for (i in seq_len(n_subjects)) {
    n_nodes <- nrow(posteriors$nodes[[i]])
    idx_range <- idx:(idx + n_nodes - 1)

    all_nodes[idx_range, ] <- posteriors$nodes[[i]]
    all_data[idx_range] <- rep(list(data_list[[i]]), n_nodes)
    all_weights[idx_range] <- posteriors$weights[[i]]

    idx <- idx + n_nodes
  }

  result <- .compute_objective_cppad(
    params = params,
    data_list = all_data,
    random_effects = all_nodes,
    parameters = parameters,
    weights = all_weights,
    gradient = gradient,
    hessian = hessian
  )

  # Explicitly clean up large temporary objects
  rm(all_nodes, all_data, all_weights)

  result
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

  fit <- suppressWarnings(
    nlm(
      f = objective,
      p = random_effect,
      hessian = FALSE,
      check.analyticals = FALSE
    )
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

  n_random_effects <- length(fit$estimate)
  quad_grid <- mvQuad::createNIGrid(
    dim = n_random_effects,
    type = "GHe",
    level = level,
    ndConstruction = "product"
  )

  quad_nodes <- mvQuad::getNodes(quad_grid)
  quad_weights <- mvQuad::getWeights(quad_grid)
  n_nodes <- nrow(quad_nodes)

  aghq_nodes <- t(sapply(seq_len(n_nodes), function(k) {
    fit$estimate + tcrossprod(quad_nodes[k, ], chol_factor)
  }))

  logpost_at_nodes <- sapply(seq_len(n_nodes), function(k) {
    result <- .compute_logpost_cppad(
      random_effect = aghq_nodes[k, ],
      data = data,
      parameters = parameters,
      gradient = FALSE,
      hessian = FALSE
    )
    as.numeric(result)
  })

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

  n_random_effects <- ncol(nodes)

  # Compute mean
  mean <- colSums(sweep(nodes, 1, weights, "*"))

  # Compute second moment (uncentered)
  second_moment <- matrix(0, nrow = n_random_effects, ncol = n_random_effects)
  for (k in seq_len(nrow(nodes))) {
    second_moment <- second_moment + weights[k] * tcrossprod(nodes[k, ])
  }

  list(mean = mean, second_moment = second_moment)
}

.compute_posteriors <- function(
  data_list,
  parameters,
  random_effects,
  parallel = FALSE,
  n_cores = 0,
  level = 3
) {
  n_subjects <- length(data_list)

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
    n_cores = n_cores
  )

  list(
    nodes = lapply(posterior_results, `[[`, "nodes"),
    weights = lapply(posterior_results, `[[`, "weights")
  )
}
