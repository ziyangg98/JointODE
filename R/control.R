#' Control Parameters for JointODE
#'
#' @description
#' Construct control parameters for the JointODE optimization algorithm.
#' This function can be called with no arguments to get defaults, or can
#' process a list to fill in missing values with defaults.
#'
#' @param maxit Maximum number of EM iterations (default: 100)
#' @param tol Convergence tolerance for EM algorithm (default: 1e-6)
#' @param verbose Logical or numeric; controls verbosity level. FALSE/0 for
#'   silent, TRUE/1 for basic progress, 2 for detailed output (default: FALSE)
#' @param optim.maxit Maximum iterations for M-step optimization (default: 50)
#' @param optim.factr Convergence factor for L-BFGS-B (default: 1e10)
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
#' control <- JointODE.control(maxit = 50, tol = 1e-4)
#'
#' # Verbose output for debugging
#' control <- JointODE.control(verbose = TRUE)
#'
#' # Process an existing list
#' my_list <- list(maxit = 200)
#' control <- JointODE.control(.list = my_list)
#'
#' @seealso \code{\link{JointODE}}
JointODE.control <- function(
  maxit = 100,
  tol = 1e-6,
  verbose = FALSE,
  optim.maxit = 50,
  optim.factr = 1e10,
  parallel = FALSE,
  n_cores = 0,
  .list = NULL,
  ...
) {
  # If a list is provided, use it as the base
  if (!is.null(.list)) {
    if (!is.list(.list)) {
      stop(".list must be a list or NULL")
    }
    control <- .list
  } else {
    # Otherwise create from arguments
    control <- list(
      maxit = maxit,
      tol = tol,
      verbose = verbose,
      optim.maxit = optim.maxit,
      optim.factr = optim.factr,
      parallel = parallel,
      n_cores = n_cores
    )

    # Add any additional parameters
    dots <- list(...)
    if (length(dots) > 0) {
      control <- c(control, dots)
    }
  }

  # Define defaults
  defaults <- list(
    maxit = 100,
    tol = 1e-6,
    verbose = FALSE,
    optim.maxit = 50,
    optim.factr = 1e10,
    parallel = FALSE,
    n_cores = 0
  )

  # Fill in any missing values with defaults
  for (name in names(defaults)) {
    if (is.null(control[[name]])) {
      control[[name]] <- defaults[[name]]
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
  if (control$tol <= 0) {
    stop("tol must be positive")
  }
  if (control$optim.maxit <= 0) {
    stop("optim.maxit must be positive")
  }
  if (control$optim.factr <= 0) {
    stop("optim.factr must be positive")
  }
  if (!is.logical(control$parallel)) {
    stop("parallel must be TRUE or FALSE")
  }
  if (!is.numeric(control$n_cores) || control$n_cores < 0) {
    stop("n_cores must be a non-negative integer")
  }

  control
}
