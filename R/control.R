#' Control Parameters for JointODE
#'
#' @description
#' Construct control parameters for the JointODE optimization algorithm.
#' This function can be called with no arguments to get defaults, or can
#' process a list to fill in missing values with defaults.
#'
#' @param maxit Maximum number of EM iterations (default: 200)
#' @param atol Absolute tolerance for parameter convergence. The EM algorithm
#'   converges when the maximum absolute change in any parameter (including
#'   variance parameters) is less than this value:
#'   max(|theta_new - theta_old|) < atol (default: 5e-4)
#' @param rtol Relative tolerance for log-likelihood convergence. The EM
#'   algorithm converges when the relative change in log-likelihood is less
#'   than this value: |L_new - L_old| / (|L_new| + epsilon) < rtol, where
#'   epsilon = 1e-8 prevents division by zero (default: 1e-5)
#' @param verbose Logical or numeric; controls verbosity level. FALSE/0 for
#'   silent, TRUE/1 for basic progress, 2 for detailed output (default: FALSE)
#' @param parallel Logical; whether to use parallel computation (default: FALSE)
#' @param n_cores Integer; number of cores to use for parallel computation.
#'   If 0, uses all available cores (default: 0)
#' @param quad_level Integer; quadrature level for numerical integration
#'   (default: 4)
#' @param .list Optional list of control parameters to process
#' @param ... Additional control parameters
#'
#' @return A list of control parameters with all defaults filled in
#'
#' @keywords internal
#' @export
#'
#' @examples
#' # Default settings
#' control <- JointODE.control()
#'
#' # Custom settings for faster exploration
#' control <- JointODE.control(maxit = 30, atol = 1e-3, rtol = 1e-3)
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
JointODE.control <- function(
  maxit = 200,
  atol = 1e-3,
  rtol = 1e-5,
  verbose = FALSE,
  parallel = FALSE,
  n_cores = 0,
  quad_level = 4,
  .list = NULL,
  ...
) {
  # Define defaults (use function arguments as source of truth)
  defaults <- list(
    maxit = maxit,
    atol = atol,
    rtol = rtol,
    verbose = verbose,
    parallel = parallel,
    n_cores = n_cores,
    quad_level = quad_level
  )

  # If a list is provided, merge with defaults
  if (!is.null(.list)) {
    if (!is.list(.list)) {
      stop(".list must be a list or NULL")
    }
    control <- defaults
    for (name in names(.list)) {
      control[[name]] <- .list[[name]]
    }
  } else {
    control <- defaults
  }

  # Add any additional parameters
  dots <- list(...)
  if (length(dots) > 0) {
    for (name in names(dots)) {
      control[[name]] <- dots[[name]]
    }
  }

  # Process verbose parameter - ensure it's numeric
  if (!is.null(control$verbose)) {
    if (is.logical(control$verbose)) {
      control$verbose <- as.numeric(control$verbose)
    } else if (!is.numeric(control$verbose)) {
      stop("verbose must be logical (TRUE/FALSE) or numeric (0, 1, 2)")
    }
  } else {
    control$verbose <- 0
  }

  # Validate parameters
  if (control$maxit <= 0) {
    stop("maxit must be positive")
  }
  if (control$atol <= 0) {
    stop("atol must be positive")
  }
  if (control$rtol <= 0) {
    stop("rtol must be positive")
  }
  if (!is.logical(control$parallel)) {
    stop("parallel must be TRUE or FALSE")
  }
  if (!is.numeric(control$n_cores) || control$n_cores < 0) {
    stop("n_cores must be a non-negative integer")
  }
  if (!is.numeric(control$quad_level) || control$quad_level < 1) {
    stop("quad_level must be a positive integer")
  }

  control
}
