#' @importFrom stats na.pass
#' @noRd
.process <- function(
  longitudinal_data,
  survival_data,
  state,
  parsed_long = NULL,
  parsed_surv = NULL,
  survival_formula = NULL,
  longitudinal_formula = NULL
) {
  # Allow backward-compatible calls with formula arguments
  if (is.null(parsed_long) && !is.null(longitudinal_formula)) {
    parsed_long <- .parse_longitudinal_formula(longitudinal_formula)
  }
  if (is.null(parsed_surv) && !is.null(survival_formula)) {
    parsed_surv <- .parse_survival_formula(survival_formula)
  }
  id <- parsed_long$grouping
  time <- parsed_surv$time_var

  # Build formulas
  fixed_formula <- .build_formula(
    parsed_long$fixed_terms,
    response = parsed_long$response
  )
  random_formula <- if (is.null(parsed_long$random_terms)) {
    NULL # No covariate random effects
  } else {
    .build_formula(parsed_long$random_terms, is_random = TRUE)
  }

  unique_ids <- unique(survival_data[[id]])
  n_subjects <- length(unique_ids)
  data_process <- vector("list", n_subjects)
  names(data_process) <- as.character(unique_ids)

  surv_frame <- model.frame(survival_formula, data = survival_data)
  surv_response_matrix <- model.response(surv_frame)
  surv_matrix <- model.matrix(survival_formula, surv_frame)
  surv_design <- if (ncol(surv_matrix) > 1) {
    surv_matrix[, -1, drop = FALSE]
  } else {
    NULL
  }

  survival_index_map <- match(unique_ids, survival_data[[id]])
  long_id_groups <- split(
    seq_len(nrow(longitudinal_data)),
    longitudinal_data[[id]]
  )

  for (i in seq_along(unique_ids)) {
    sid <- as.character(unique_ids[i])

    long_rows <- long_id_groups[[sid]]
    if (!is.null(long_rows) && length(long_rows) > 0) {
      long_subset <- longitudinal_data[long_rows, , drop = FALSE]
      long_subset <- long_subset[order(long_subset[[time]]), , drop = FALSE]

      long_frame <- model.frame(fixed_formula, data = long_subset)
      long_times <- long_subset[[time]]
      long_measurements <- model.response(long_frame)
      long_fixed_covariates <- model.matrix(fixed_formula, long_frame)
      long_random_covariates <- if (is.null(random_formula)) {
        matrix(nrow = nrow(long_subset), ncol = 0)
      } else {
        model.matrix(random_formula, long_subset)
      }
    } else {
      long_times <- numeric(0)
      long_measurements <- numeric(0)
      # Get dimensions from empty data frame
      empty_data <- longitudinal_data[0, , drop = FALSE]
      n_fixed <- ncol(model.matrix(fixed_formula, empty_data))
      long_fixed_covariates <- matrix(nrow = 0, ncol = n_fixed)
      long_random_covariates <- if (is.null(random_formula)) {
        matrix(nrow = 0, ncol = 0)
      } else {
        n_random <- ncol(model.matrix(random_formula, empty_data))
        matrix(nrow = 0, ncol = n_random)
      }
    }

    survival_row <- survival_index_map[i]
    event_time <- surv_response_matrix[survival_row, 1]
    event_status <- surv_response_matrix[survival_row, 2]
    covariates <- if (!is.null(surv_design)) {
      surv_design[survival_row, , drop = FALSE]
    } else {
      data.frame()
    }
    initial_state <- if (!is.null(state)) {
      state[i, , drop = TRUE]
    } else {
      m0 <- long_measurements[1]
      v0 <- if (length(long_measurements) >= 2) {
        (long_measurements[2] - long_measurements[1]) /
          (long_times[2] - long_times[1])
      } else {
        0
      }
      c(m0, v0)
    }

    data_process[[i]] <- list(
      id = unique_ids[i],
      time = event_time,
      status = event_status,
      covariates = covariates,
      initial_state = initial_state,
      longitudinal = list(
        times = long_times,
        measurements = long_measurements,
        covariates = list(
          fixed = long_fixed_covariates,
          random = long_random_covariates
        )
      )
    )
  }
  data_process
}
