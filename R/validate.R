.validate_joint <- function(
  longitudinal_formula,
  survival_formula,
  longitudinal_data,
  survival_data,
  gamma,
  spline_baseline,
  init,
  parsed_long = .parse_longitudinal_formula(longitudinal_formula),
  parsed_surv = .parse_survival_formula(survival_formula),
  spline_config = modifyList(.default_spline, spline_baseline)
) {
  .validate_longitudinal_formula(
    formula = longitudinal_formula,
    data = longitudinal_data,
    parsed = parsed_long
  )

  .validate_survival_formula(
    formula = survival_formula,
    data = survival_data,
    parsed = parsed_surv
  )

  .validate_data(longitudinal_data, survival_data, parsed_long, parsed_surv)
  .validate_gamma(gamma)
  .validate_spline(spline_baseline)

  .validate_initial(
    init = init,
    parsed_long = parsed_long,
    parsed_surv = parsed_surv,
    spline_config = spline_config
  )

  invisible(NULL)
}


.validate_gamma <- function(gamma) {
  if (!is.numeric(gamma) || length(gamma) != 1) {
    stop("gamma must be a single numeric value", call. = FALSE)
  }
  if (!gamma %in% c(0, 1, 2)) {
    stop("gamma must be 0, 1, or 2", call. = FALSE)
  }

  invisible(TRUE)
}


#' @importFrom stats setNames
#' @noRd
.validate_data <- function(
  longitudinal_data,
  survival_data,
  parsed_long,
  parsed_surv
) {
  id <- parsed_long$grouping
  time <- parsed_surv$time_var

  if (!time %in% names(longitudinal_data)) {
    stop(
      sprintf("Time variable '%s' not in longitudinal data", time),
      call. = FALSE
    )
  }
  if (!id %in% names(survival_data)) {
    stop(
      sprintf("ID variable '%s' not in survival data", id),
      call. = FALSE
    )
  }

  critical_cols <- list(
    "ID in longitudinal data" = longitudinal_data[[id]],
    "Time in longitudinal data" = longitudinal_data[[time]],
    "ID in survival data" = survival_data[[id]]
  )

  for (col_name in names(critical_cols)) {
    if (any(is.na(critical_cols[[col_name]]))) {
      stop(sprintf("Missing values in %s", col_name), call. = FALSE)
    }
  }

  if (any(duplicated(survival_data[[id]]))) {
    stop("Duplicate IDs in survival data", call. = FALSE)
  }

  long_ids <- unique(longitudinal_data[[id]])
  surv_ids <- unique(survival_data[[id]])

  orphaned_ids <- setdiff(long_ids, surv_ids)
  if (length(orphaned_ids) > 0) {
    stop(
      sprintf(
        "Subjects in longitudinal but not survival data: %s",
        paste(head(orphaned_ids, 5), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  missing_long <- setdiff(surv_ids, long_ids)
  if (length(missing_long) > 0) {
    stop(
      sprintf(
        "Subjects in survival but not longitudinal data: %s",
        paste(head(missing_long, 5), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (any(longitudinal_data[[time]] < 0, na.rm = TRUE)) {
    stop("Negative time values in longitudinal data", call. = FALSE)
  }

  long_times_by_id <- split(longitudinal_data[[time]], longitudinal_data[[id]])
  surv_times_map <- setNames(survival_data[[time]], survival_data[[id]])

  subjects_with_late_obs <- character()
  for (subject_id in long_ids) {
    if (subject_id %in% surv_ids) {
      long_times <- long_times_by_id[[as.character(subject_id)]]
      surv_time <- surv_times_map[as.character(subject_id)]

      if (!is.na(surv_time) && any(long_times > surv_time + 1e-6)) {
        subjects_with_late_obs <- c(subjects_with_late_obs, subject_id)
      }
    }
  }

  if (length(subjects_with_late_obs) > 0) {
    warning(sprintf(
      "%d subjects have measurements after observation time: %s",
      length(subjects_with_late_obs),
      paste(head(subjects_with_late_obs, 3), collapse = ", ")
    ))
  }

  observations_per_subject <- table(longitudinal_data[[id]])
  if (all(observations_per_subject == 1)) {
    warning(
      paste(
        "Each subject has only one longitudinal observation -",
        "joint modeling may not be appropriate"
      )
    )
  }

  invisible(TRUE)
}

.validate_longitudinal_formula <- function(
  formula,
  data,
  parsed = .parse_longitudinal_formula(formula)
) {
  if (!inherits(formula, "formula")) {
    stop("Formula must be a formula object", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("Data must be a data frame", call. = FALSE)
  }
  if (nrow(data) == 0) {
    stop("Longitudinal data has no rows", call. = FALSE)
  }
  if (length(formula) != 3) {
    stop("Formula must be two-sided: y ~ x", call. = FALSE)
  }

  formula_str <- deparse(formula, width.cutoff = 500)

  pipe_matches <- gregexpr("\\|", formula_str)[[1]]
  pipe_count <- if (pipe_matches[1] == -1L) 0L else length(pipe_matches)
  if (pipe_count > 1) {
    stop(
      "Multiple random effects groupings are not supported. ",
      "Use a single grouping: y ~ x + (terms | id)",
      call. = FALSE
    )
  }

  if (grepl("\\|[^)]*\\/", formula_str)) {
    stop(
      "Nested random effects are not supported. ",
      "Use simple grouping: y ~ x + (terms | id)",
      call. = FALSE
    )
  }

  parts <- parsed

  if (is.null(parts$grouping)) {
    stop(
      "Random effects specification is required in longitudinal formula.\n",
      "Use syntax: response ~ fixed_terms + (random_terms | grouping)\n",
      "Example: y ~ x + (biomarker + velocity | id)",
      call. = FALSE
    )
  }

  if (!(parts$grouping %in% names(data))) {
    stop(
      "Grouping variable '",
      parts$grouping,
      "' not found in data",
      call. = FALSE
    )
  }

  reserved_in_data <- intersect(.reserved_words, names(data))
  if (length(reserved_in_data) > 0) {
    stop(
      "Reserved words cannot be used as variable names in data: ",
      paste(reserved_in_data, collapse = ", "),
      call. = FALSE
    )
  }

  all_vars <- all.vars(formula)
  reserved_in_formula <- intersect(.reserved_words, all_vars)
  vars_to_check <- setdiff(all_vars, reserved_in_formula)
  missing_vars <- setdiff(vars_to_check, names(data))

  if (length(missing_vars) > 0) {
    stop(
      "Variables not found: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.validate_survival_formula <- function(
  formula,
  data,
  parsed = .parse_survival_formula(formula)
) {
  if (!inherits(formula, "formula")) {
    stop("Formula must be a formula object", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("Data must be a data frame", call. = FALSE)
  }
  if (nrow(data) == 0) {
    stop("Survival data has no rows", call. = FALSE)
  }
  if (length(formula) != 3) {
    stop("Formula must be two-sided: Surv(...) ~ x", call. = FALSE)
  }

  parts <- parsed

  all_vars <- c(parts$time_var, parts$status_var, parts$covariate_terms)
  missing_vars <- setdiff(all_vars, names(data))

  if (length(missing_vars) > 0) {
    stop(
      "Variables not found in data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  if (any(data[[parts$time_var]] <= 0, na.rm = TRUE)) {
    stop(
      "Invalid observation times in '",
      parts$time_var,
      "' (must be positive)",
      call. = FALSE
    )
  }

  unique_status <- unique(data[[parts$status_var]])
  invalid_status <- setdiff(unique_status, c(0, 1, NA))

  if (length(invalid_status) > 0) {
    stop(
      "Invalid status values in '",
      parts$status_var,
      "': ",
      paste(invalid_status, collapse = ", "),
      ". Must be 0 (censored) or 1 (event)",
      call. = FALSE
    )
  }

  if (any(is.na(data[[parts$time_var]]))) {
    stop(
      "Missing values found in time variable '",
      parts$time_var,
      "'",
      call. = FALSE
    )
  }

  if (any(is.na(data[[parts$status_var]]))) {
    stop(
      "Missing values found in status variable '",
      parts$status_var,
      "'",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.validate_spline <- function(spline_baseline) {
  if (!is.list(spline_baseline)) {
    stop("spline_baseline must be a list", call. = FALSE)
  }

  valid_baseline_params <- c(
    "degree",
    "n_knots",
    "knot_placement",
    "boundary_knots"
  )
  invalid_baseline <- setdiff(names(spline_baseline), valid_baseline_params)
  if (length(invalid_baseline) > 0) {
    stop(sprintf(
      "Invalid parameters in spline_baseline: %s. Valid parameters are: %s",
      paste(invalid_baseline, collapse = ", "),
      paste(valid_baseline_params, collapse = ", ")
    ))
  }

  if ("degree" %in% names(spline_baseline)) {
    if (
      !is.numeric(spline_baseline$degree) ||
        length(spline_baseline$degree) != 1 ||
        spline_baseline$degree < 1 ||
        spline_baseline$degree > 5
    ) {
      stop("spline_baseline$degree must be a single integer between 1 and 5")
    }
  }

  if ("n_knots" %in% names(spline_baseline)) {
    if (
      !is.numeric(spline_baseline$n_knots) ||
        length(spline_baseline$n_knots) != 1 ||
        spline_baseline$n_knots < 0 ||
        spline_baseline$n_knots > 20
    ) {
      stop("spline_baseline$n_knots must be a single integer between 0 and 20")
    }
  }

  if ("knot_placement" %in% names(spline_baseline)) {
    valid_placements <- c("quantile", "equal")
    if (!spline_baseline$knot_placement %in% valid_placements) {
      stop(sprintf(
        "spline_baseline$knot_placement must be one of: %s",
        paste(valid_placements, collapse = ", ")
      ))
    }
  }

  if ("boundary_knots" %in% names(spline_baseline)) {
    if (!is.null(spline_baseline$boundary_knots)) {
      if (
        !is.numeric(spline_baseline$boundary_knots) ||
          length(spline_baseline$boundary_knots) != 2
      ) {
        stop(paste(
          "spline_baseline$boundary_knots must be NULL or",
          "a numeric vector of length 2"
        ))
      }
      if (
        spline_baseline$boundary_knots[1] >= spline_baseline$boundary_knots[2]
      ) {
        stop(paste(
          "spline_baseline$boundary_knots[1] must be less than",
          "boundary_knots[2]"
        ))
      }
    }
  }

  invisible(TRUE)
}


.validate_initial <- function(
  init,
  parsed_long,
  parsed_surv,
  spline_config
) {
  if (is.character(init)) {
    if (!init %in% c("default", "marginal")) {
      stop("init: must be 'default', 'marginal', or a list", call. = FALSE)
    }
    return(invisible(TRUE))
  }

  dims <- .compute_dimensions(parsed_long, parsed_surv, spline_config)

  if (!is.list(init)) {
    stop("init: must be 'default', 'marginal', or a list", call. = FALSE)
  }

  valid_components <- c("coefficients", "configurations")
  unknown_components <- setdiff(names(init), valid_components)
  if (length(unknown_components) > 0) {
    stop(
      sprintf(
        "init: unknown components '%s'",
        paste(unknown_components, collapse = "', '")
      ),
      call. = FALSE
    )
  }

  if (!is.null(init$coefficients)) {
    if (!is.list(init$coefficients)) {
      stop("init$coefficients: must be a list", call. = FALSE)
    }

    valid_coefs <- c(
      "baseline",
      "hazard",
      "longitudinal",
      "initial_state",
      "measurement_error_sd",
      "random_effect_sigma"
    )
    unknown_coefs <- setdiff(names(init$coefficients), valid_coefs)
    if (length(unknown_coefs) > 0) {
      stop(
        sprintf(
          "init$coefficients: unknown types '%s'",
          paste(unknown_coefs, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    if (!is.null(init$coefficients$baseline)) {
      if (!is.numeric(init$coefficients$baseline)) {
        stop("baseline: must be numeric", call. = FALSE)
      }
      if (any(!is.finite(init$coefficients$baseline))) {
        stop("baseline: must contain finite values", call. = FALSE)
      }
      if (length(init$coefficients$baseline) != dims$n_spline_basis) {
        stop(
          sprintf(
            "baseline: wrong length (expected %d, got %d)",
            dims$n_spline_basis,
            length(init$coefficients$baseline)
          ),
          call. = FALSE
        )
      }
    }

    if (!is.null(init$coefficients$hazard)) {
      if (!is.numeric(init$coefficients$hazard)) {
        stop("hazard: must be numeric", call. = FALSE)
      }
      if (any(!is.finite(init$coefficients$hazard))) {
        stop("hazard: must contain finite values", call. = FALSE)
      }
      if (length(init$coefficients$hazard) < 2) {
        stop("hazard: need at least 2 elements", call. = FALSE)
      }
      expected_len <- dims$n_survival_covariates + 2
      if (length(init$coefficients$hazard) != expected_len) {
        stop(
          sprintf(
            "hazard: wrong length (expected %d, got %d)",
            expected_len,
            length(init$coefficients$hazard)
          ),
          call. = FALSE
        )
      }
    }

    if (!is.null(init$coefficients$longitudinal)) {
      if (!is.numeric(init$coefficients$longitudinal)) {
        stop("longitudinal: must be numeric", call. = FALSE)
      }
      if (any(!is.finite(init$coefficients$longitudinal))) {
        stop("longitudinal: must contain finite values", call. = FALSE)
      }
      if (length(init$coefficients$longitudinal) != dims$n_longitudinal_coef) {
        stop(
          sprintf(
            "longitudinal: wrong length (expected %d, got %d)",
            dims$n_longitudinal_coef,
            length(init$coefficients$longitudinal)
          ),
          call. = FALSE
        )
      }
    }

    if (!is.null(init$coefficients$measurement_error_sd)) {
      if (
        !is.numeric(init$coefficients$measurement_error_sd) ||
          length(init$coefficients$measurement_error_sd) != 1
      ) {
        stop("measurement_error_sd: must be scalar", call. = FALSE)
      }
      if (!is.finite(init$coefficients$measurement_error_sd)) {
        stop("measurement_error_sd: must be finite", call. = FALSE)
      }
      if (init$coefficients$measurement_error_sd <= 0) {
        stop("measurement_error_sd: must be positive", call. = FALSE)
      }
    }

    if (!is.null(init$coefficients$random_effect_sigma)) {
      if (!is.matrix(init$coefficients$random_effect_sigma)) {
        stop("random_effect_sigma: must be matrix", call. = FALSE)
      }
      if (
        nrow(init$coefficients$random_effect_sigma) !=
          ncol(init$coefficients$random_effect_sigma)
      ) {
        stop("random_effect_sigma: must be square", call. = FALSE)
      }
      if (
        nrow(init$coefficients$random_effect_sigma) != dims$n_random_effects
      ) {
        stop(
          sprintf(
            "random_effect_sigma: wrong dim (expected %dx%d, got %dx%d)",
            dims$n_random_effects,
            dims$n_random_effects,
            nrow(init$coefficients$random_effect_sigma),
            ncol(init$coefficients$random_effect_sigma)
          ),
          call. = FALSE
        )
      }
      if (!all(is.finite(init$coefficients$random_effect_sigma))) {
        stop("random_effect_sigma: must contain finite values", call. = FALSE)
      }
      eig_vals <- eigen(
        init$coefficients$random_effect_sigma,
        symmetric = TRUE,
        only.values = TRUE
      )$values
      if (any(eig_vals <= 0)) {
        stop("random_effect_sigma: must be positive definite", call. = FALSE)
      }
    }
  }

  if (!is.null(init$configurations)) {
    if (!is.list(init$configurations)) {
      stop("init$configurations: must be a list", call. = FALSE)
    }
    if (!is.null(init$configurations$baseline)) {
      if (!is.list(init$configurations$baseline)) {
        stop("init$configurations$baseline: must be a list", call. = FALSE)
      }
    }
    if (!is.null(init$configurations$biomarker)) {
      if (!is.list(init$configurations$biomarker)) {
        stop("init$configurations$biomarker: must be a list", call. = FALSE)
      }
      if (
        !all(c("fixed", "random") %in% names(init$configurations$biomarker))
      ) {
        stop(
          paste0(
            "init$configurations$biomarker: ",
            "must have 'fixed' and 'random' fields"
          ),
          call. = FALSE
        )
      }
      if (
        !is.logical(init$configurations$biomarker$fixed) ||
          !is.logical(init$configurations$biomarker$random)
      ) {
        stop(
          "init$configurations$biomarker: 'fixed' and 'random' must be logical",
          call. = FALSE
        )
      }
    }
    if (!is.null(init$configurations$velocity)) {
      if (!is.list(init$configurations$velocity)) {
        stop("init$configurations$velocity: must be a list", call. = FALSE)
      }
      if (!all(c("fixed", "random") %in% names(init$configurations$velocity))) {
        stop(
          "init$configurations$velocity: must have 'fixed' and 'random' fields",
          call. = FALSE
        )
      }
      if (
        !is.logical(init$configurations$velocity$fixed) ||
          !is.logical(init$configurations$velocity$random)
      ) {
        stop(
          "init$configurations$velocity: 'fixed' and 'random' must be logical",
          call. = FALSE
        )
      }
    }
  }

  invisible(TRUE)
}

#' @noRd
.validate_marginal <- function(formula, data, time, id, state) {
  stopifnot(
    "formula must be a formula" = inherits(formula, "formula"),
    "data must be a data.frame or matrix" =
      is.data.frame(data) || is.matrix(data),
    "time must be a character string" = is.character(time),
    "id must be a character string" = is.character(id)
  )
  if (is.data.frame(data) || is.matrix(data)) {
    stopifnot(
      "time column not found in data" = time %in% names(data),
      "id column not found in data" = id %in% names(data)
    )
  }
  if (!is.null(state)) {
    stopifnot(
      "state must be a matrix with 2 columns" =
        is.matrix(state) && ncol(state) == 2
    )
  }
}
