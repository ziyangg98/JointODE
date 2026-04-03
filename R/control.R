#' Control Parameters for JointODE
#'
#' @description
#' Construct control parameters for the JointODE optimization algorithm.
#' This function can be called with no arguments to get defaults, or can
#' process a list to fill in missing values with defaults.
#'
#' @param maxit Maximum number of MCEM iterations (default: 200)
#' @param tol Convergence tolerance. The MCEM algorithm converges
#'   when max|theta_new - theta_old| < tol (default: 1e-4)
#' @param verbose Logical or numeric; controls verbosity level. FALSE/0 for
#'   silent, TRUE/1 for basic progress, 2 for detailed output (default: FALSE)
#' @param parallel Logical; whether to use parallel computation (default: FALSE)
#' @param n_cores Integer; number of cores to use for parallel computation.
#'   If 0, uses all available cores (default: 0)
#' @param hazard_quadrature Integer; number of Simpson sub-intervals per
#'   observation interval for hazard integration (default: 1)
#' @param mc_samples Integer; number of MCMC samples per subject for
#'   Monte Carlo EM (MCEM). Must be a positive integer.
#'   Recommended range: 50-500 for MCEM (default: 200).
#' @param mc_burnin Integer; number of burn-in iterations for the
#'   random-walk Metropolis sampler in the MCEM E-step (default: 200).
#' @param .list Optional list of control parameters to process
#' @param ... Additional control parameters
#'
#' @return A list of control parameters with all defaults filled in
#'
#' @concept model-fitting
#' @keywords internal
#' @export
#'
#' @examples
#' # Default settings
#' control <- JointODE.control()
#'
#' # Custom settings for faster exploration
#' control <- JointODE.control(maxit = 30, tol = 1e-4)
#'
#' # Verbose output for debugging
#' control <- JointODE.control(verbose = TRUE)
#'
#' # Parallel computation
#' control <- JointODE.control(parallel = TRUE, n_cores = 4)
#'
#' # Process an existing list
#' my_list <- list(maxit = 200)
#' control <- JointODE.control(.list = my_list)
#'
#' @seealso \code{\link{JointODE}}
# nolint next: object_name_linter
JointODE.control <- function(
  maxit = 200,
  tol = 1e-4,
  verbose = FALSE,
  parallel = FALSE,
  n_cores = 0,
  hazard_quadrature = 1,
  mc_samples = 200,
  mc_burnin = 200,
  .list = NULL,
  ...
) {
  defaults <- list(
    maxit = maxit, tol = tol, verbose = verbose,
    parallel = parallel, n_cores = n_cores,
    hazard_quadrature = hazard_quadrature,
    mc_samples = mc_samples,
    mc_burnin = mc_burnin
  )

  if (!is.null(.list)) {
    if (!is.list(.list)) stop(".list must be a list or NULL")
    control <- defaults
    for (name in names(.list)) control[[name]] <- .list[[name]]
  } else {
    control <- defaults
  }

  dots <- list(...)
  for (name in names(dots)) control[[name]] <- dots[[name]]

  control$verbose <- as.numeric(control$verbose)

  if (control$maxit <= 0) stop("maxit must be positive")
  if (control$tol <= 0) stop("tol must be positive")
  if (!is.logical(control$parallel)) stop("parallel must be TRUE or FALSE")
  if (!is.numeric(control$n_cores) || control$n_cores < 0) {
    stop("n_cores must be a non-negative integer")
  }
  if (!is.numeric(control$hazard_quadrature) || control$hazard_quadrature < 1 ||
      control$hazard_quadrature != as.integer(control$hazard_quadrature)) {
    stop("hazard_quadrature must be a positive integer")
  }
  if (!is.numeric(control$mc_samples) || control$mc_samples < 1 ||
      control$mc_samples != as.integer(control$mc_samples)) {
    stop("mc_samples must be a positive integer")
  }
  if (!is.numeric(control$mc_burnin) || control$mc_burnin < 0 ||
      control$mc_burnin != as.integer(control$mc_burnin)) {
    stop("mc_burnin must be a non-negative integer")
  }

  control
}

#' Control Parameters for MarginalODE
#'
#' @description
#' Construct control parameters for the MarginalODE optimization.
#'
#' @param maxit Maximum number of alternating optimization iterations
#'   (default: 200)
#' @param tol Convergence tolerance on max absolute parameter
#'   change (default: 1e-4)
#' @param verbose Logical or numeric; FALSE/0 for silent, TRUE/1 for basic
#'   progress, 2 for detailed output (default: FALSE)
#' @param parallel Logical; whether to use parallel computation (default: FALSE)
#' @param n_cores Integer; number of cores (0 = auto) (default: 0)
#' @param .list Optional list of control parameters to process
#' @param ... Additional control parameters
#'
#' @return A list of control parameters with all defaults filled in
#'
#' @concept model-fitting
#' @keywords internal
#' @export
#'
#' @examples
#' control <- MarginalODE.control()
#' control <- MarginalODE.control(maxit = 50, tol = 1e-4)
#'
#' @seealso \code{\link{MarginalODE}}
# nolint next: object_name_linter
MarginalODE.control <- function(
  maxit = 200,
  tol = 1e-4,
  verbose = FALSE,
  parallel = FALSE,
  n_cores = 0,
  .list = NULL,
  ...
) {
  defaults <- list(
    maxit = maxit, tol = tol, verbose = verbose,
    parallel = parallel, n_cores = n_cores
  )

  if (!is.null(.list)) {
    if (!is.list(.list)) stop(".list must be a list or NULL")
    control <- defaults
    for (name in names(.list)) control[[name]] <- .list[[name]]
  } else {
    control <- defaults
  }

  dots <- list(...)
  for (name in names(dots)) control[[name]] <- dots[[name]]

  control$verbose <- as.numeric(control$verbose)

  if (control$maxit <= 0) stop("maxit must be positive")
  if (control$tol <= 0) stop("tol must be positive")
  if (!is.logical(control$parallel)) stop("parallel must be TRUE or FALSE")
  if (!is.numeric(control$n_cores) || control$n_cores < 0) {
    stop("n_cores must be a non-negative integer")
  }

  control
}
