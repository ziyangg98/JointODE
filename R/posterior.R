# Posterior Computation and EM Algorithm

#' @importFrom stats nlm
NULL

# Variance Updates =============================================================

#' @noRd
.update_random_effect_sigma <- function(posterior_moments, n_subjects) {
  n_re <- length(posterior_moments[[1]]$mean)
  sigma_sum <- matrix(0, n_re, n_re)
  for (i in seq_len(n_subjects)) {
    sigma_sum <- sigma_sum + posterior_moments[[i]]$second_moment
  }
  sigma_sum / n_subjects
}

# Expected Objective ===========================================================

.compute_objective_expected <- function(
  params, data_list, random_effects, parameters,
  gradient = TRUE, hessian = FALSE,
  parallel = FALSE, n_cores = 0
) {
  eval_chunk <- function(idx) {
    .compute_joint_objective(
      params = params, data_list = data_list[idx],
      random_effects = random_effects[idx, , drop = FALSE],
      parameters = parameters,
      gradient = gradient, hessian = hessian
    )
  }

  n <- length(data_list)
  if (!parallel || n <= 1) return(eval_chunk(seq_len(n)))

  n_cores <- .resolve_cores(n_cores)
  chunks <- split(seq_len(n), cut(seq_len(n), min(n_cores, n),
    labels = FALSE
  ))

  results <- .parallel_apply(chunks, eval_chunk,
    parallel = TRUE, n_cores = n_cores, setup = FALSE
  )

  obj <- sum(vapply(results, as.numeric, numeric(1)))
  result <- structure(obj, names = NULL)
  if (gradient) {
    attr(result, "gradient") <- Reduce(`+`, lapply(
      results, function(r) as.vector(attr(r, "gradient"))
    ))
  }
  if (hessian) {
    attr(result, "hessian") <- Reduce(`+`, lapply(
      results, function(r) attr(r, "hessian")
    ))
  }
  result
}

# Laplace Approximation ========================================================

.compute_posterior_laplace <- function(data, random_effect, parameters) {
  objective <- function(theta) {
    result <- .compute_joint_logpost(
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

  result_at_mode <- .compute_joint_logpost(
    random_effect = fit$estimate, data = data,
    parameters = parameters, gradient = FALSE, hessian = TRUE
  )
  R <- .safe_chol(-attr(result_at_mode, "hessian"))

  list(mode = fit$estimate, cov = chol2inv(R))
}

.compute_posteriors <- function(
  data_list, parameters, random_effects,
  parallel = FALSE, n_cores = 0, setup = TRUE
) {
  .parallel_apply(
    seq_along(data_list),
    function(i) {
      .compute_posterior_laplace(
        data = data_list[[i]],
        random_effect = random_effects[i, ],
        parameters = parameters
      )
    },
    parallel = parallel, n_cores = n_cores, setup = setup
  )
}

# SEM Variance-Covariance =====================================================

#' @noRd
.compute_vcov_sem <- function(
  data_list, parameters, random_effects, control
) {
  theta <- .coef_to_vector(parameters)

  if (control$verbose > 0) cli::cli_alert_info("Computing SEM vcov...")

  obj <- .compute_objective_expected(
    theta, data_list, random_effects, parameters,
    gradient = FALSE, hessian = TRUE,
    parallel = control$parallel, n_cores = control$n_cores
  )
  hess <- attr(obj, "hessian")
  n_coef <- nrow(hess)
  info_complete_inv <- solve(hess)

  em_map <- function(theta_input) {
    params_input <- .vector_to_coef(parameters, theta_input)
    result <- .em_step(data_list, params_input, random_effects, control)
    .coef_to_vector(result$parameters)
  }

  eps_weight <- sqrt(abs(diag(info_complete_inv)))
  eps_weight <- eps_weight / min(eps_weight) * sqrt(.Machine$double.eps)

  dm_matrix <- numDeriv::jacobian(
    func = em_map, x = theta,
    method = "simple", method.args = list(eps = eps_weight)
  )

  info_observed_inv <- info_complete_inv %*% solve(diag(n_coef) - t(dm_matrix))

  asymmetry <- max(abs(info_observed_inv - t(info_observed_inv)))
  if (control$verbose > 0 && asymmetry > 1e-4) {
    cli::cli_alert_warning(sprintf(
      "Vcov matrix asymmetry: %.2e", asymmetry
    ))
  }
  vcov_sym <- (info_observed_inv + t(info_observed_inv)) / 2

  eigenvalues <- eigen(vcov_sym, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigenvalues <= 0) && control$verbose > 0) {
    cli::cli_alert_warning(sprintf(
      "%d non-positive eigenvalue(s) (min: %.2e)",
      sum(eigenvalues <= 0), min(eigenvalues)
    ))
  }

  diag_ratio <- diag(vcov_sym) / diag(info_complete_inv)
  if (any(diag_ratio < 1 - 0.01) && control$verbose > 0) {
    cli::cli_alert_warning(sprintf(
      "%d parameter(s) violate missing information principle (min ratio: %.3f)",
      sum(diag_ratio < 1 - 0.01), min(diag_ratio)
    ))
  }

  vcov_sym
}

# EM Step ======================================================================

#' @noRd
.em_step <- function(data_list, parameters, random_effects, control) {
  n_subjects <- length(data_list)
  n_re <- ncol(random_effects)

  # Cache last accepted E-step results to avoid redundant recomputation
  last <- list(posteriors = NULL, moments = NULL, re = random_effects)

  # Profile Laplace: re-solve b̂(θ) at each evaluation
  theta <- .coef_to_vector(parameters)
  fit_m <- trust::trust(
    objfun = function(th) {
      params_th <- .vector_to_coef(parameters, th)
      posts <- .compute_posteriors(
        data_list, params_th, last$re,
        control$parallel, control$n_cores, setup = FALSE
      )
      re_new <- t(vapply(posts, function(p) p$mode, numeric(n_re)))
      moments <- lapply(posts, function(p) {
        list(mean = p$mode, second_moment = tcrossprod(p$mode) + p$cov)
      })
      params_th$coefficients$random_effect_sigma <-
        .update_random_effect_sigma(moments, n_subjects)
      last <<- list(posteriors = posts, moments = moments, re = re_new)

      obj <- .compute_objective_expected(
        th, data_list, re_new, params_th,
        gradient = TRUE, hessian = TRUE,
        parallel = control$parallel, n_cores = control$n_cores
      )
      list(value = as.numeric(obj),
           gradient = as.vector(attr(obj, "gradient")),
           hessian = as.matrix(attr(obj, "hessian")))
    },
    parinit = theta, rinit = 1, rmax = 10, iterlim = 1
  )

  parameters <- .vector_to_coef(parameters, fit_m$argument)
  random_effects <- last$re
  parameters$coefficients$random_effect_sigma <-
    .update_random_effect_sigma(last$moments, n_subjects)

  list(
    parameters = parameters,
    random_effects = random_effects,
    loglik = -fit_m$value,
    posteriors = last$posteriors,
    posterior_moments = last$moments
  )
}
