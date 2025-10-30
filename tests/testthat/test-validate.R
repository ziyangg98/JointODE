# Comprehensive validation tests

# ==============================================================================
# Formula and Data Type Validation
# ==============================================================================

test_that(".validate checks formula types", {
  # Non-formula longitudinal_formula
  expect_error(
    .validate(
      longitudinal_formula = "v ~ 1 + (1 | id)",
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Formula must be a formula object"
  )

  # Non-formula survival_formula
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = "Surv(time, status) ~ 1",
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Formula must be a formula object"
  )

  # Surv() with insufficient arguments
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Surv\\(\\) must have at least time and status arguments"
  )
})

test_that(".validate checks data frame types", {
  # Non-data.frame longitudinal_data
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = list(id = 1, time = 0, v = 1),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Data must be a data frame"
  )

  # Non-data.frame survival_data
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = list(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Data must be a data frame"
  )
})

test_that(".validate checks for empty data", {
  # Empty longitudinal data
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = numeric(0),
        time = numeric(0),
        v = numeric(0)
      ),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Longitudinal data has no rows"
  )

  # Empty survival data
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(
        id = numeric(0),
        time = numeric(0),
        status = numeric(0)
      ),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Survival data has no rows"
  )
})

# ==============================================================================
# Required Columns Validation
# ==============================================================================

test_that(".validate checks required columns", {
  # Missing id column in longitudinal data
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(time = 0, v = 1),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Grouping variable 'id' not found"
  )

  # Missing time column in longitudinal data
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = 1, v = 1),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Time variable 'time' not in longitudinal data"
  )

  # Missing id column in survival data
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "ID variable 'id' not in survival data"
  )
})

# ==============================================================================
# Formula Variable Validation
# ==============================================================================

test_that(".validate checks formula variables exist", {
  # Missing longitudinal response variable
  expect_error(
    .validate(
      longitudinal_formula = missing_var ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Variables not found: missing_var"
  )

  # Missing longitudinal predictor variable
  expect_error(
    .validate(
      longitudinal_formula = v ~ missing_pred + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Variables not found: missing_pred"
  )

  # Missing survival predictor variable
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ missing_surv,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Variables not found in data: missing_surv"
  )
})

# ==============================================================================
# Surv Formula Validation
# ==============================================================================

test_that(".validate checks Surv formula structure", {
  # Invalid Surv formula (no Surv on LHS)
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = time ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Survival formula must have Surv\\(\\) on the left-hand side"
  )
})

# ==============================================================================
# Missing Values Validation
# ==============================================================================

test_that(".validate checks for missing values", {
  # Missing values in longitudinal data
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = c(1, 1),
        time = c(NA, 1),
        v = c(1, 2)
      ),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Missing values in Time in longitudinal data"
  )

  # Missing values in survival status
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = NA),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Missing values found in status variable 'status'"
  )
})

# ==============================================================================
# ID Consistency Validation
# ==============================================================================

test_that(".validate checks ID consistency", {
  # Duplicate IDs in survival data
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = rep(1, 2), time = 0:1, v = 1:2),
      survival_data = data.frame(
        id = c(1, 1),
        time = c(1, 2),
        status = c(1, 0)
      ),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Duplicate IDs in survival data"
  )

  # Subjects in longitudinal not in survival
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = rep(1:3, each = 2),
        time = rep(c(0, 1), 3),
        v = 1:6
      ),
      survival_data = data.frame(id = 1:2, time = c(1, 2), status = c(1, 0)),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Subjects in longitudinal but not survival data: 3"
  )

  # Subjects in survival without longitudinal data (error)
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = rep(1:2, each = 2),
        time = rep(c(0, 1), 2),
        v = 1:4
      ),
      survival_data = data.frame(
        id = 1:4,
        time = c(1, 2, 3, 4),
        status = c(1, 0, 1, 1)
      ),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Subjects in survival but not longitudinal data: 3, 4"
  )
})

# ==============================================================================
# Time Values Validation
# ==============================================================================

test_that(".validate checks time values", {
  # Negative time in longitudinal data
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = c(1, 1),
        time = c(-1, 0),
        v = c(1, 2)
      ),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Negative time values in longitudinal data"
  )

  # Non-positive time in survival data
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 0, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Invalid observation times in 'time' \\(must be positive\\)"
  )

  # Warning: measurements after observation time
  expect_warning(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = rep(1:2, each = 3),
        time = rep(c(0, 1, 2), 2),
        v = 1:6
      ),
      survival_data = data.frame(
        id = 1:2,
        time = c(1.5, 2.5),
        status = c(1, 1)
      ),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "1 subjects have measurements after observation time: 1"
  )

  # Warning: single longitudinal observation
  expect_warning(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = 1:2,
        time = c(0, 0),
        v = c(1, 2)
      ),
      survival_data = data.frame(id = 1:2, time = c(1, 2), status = c(1, 1)),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Each subject has only one longitudinal observation"
  )
})

# ==============================================================================
# Status Values Validation
# ==============================================================================

test_that(".validate checks status values", {
  # Invalid status values
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 2),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "Invalid status values in 'status': 2\\. Must be 0 \\(censored\\) or 1 \\(event\\)"
  )
})

# ==============================================================================
# State Matrix Validation
# ==============================================================================

test_that(".validate checks state matrix", {
  # Non-matrix state
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = c(1, 2),
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "state.*must be a matrix"
  )

  # Wrong number of rows in state
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = rep(1:2, each = 2),
        time = rep(0:1, 2),
        v = 1:4
      ),
      survival_data = data.frame(id = 1:2, time = c(1, 2), status = c(1, 0)),
      state = matrix(c(1, 2), nrow = 1),
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "state.*wrong rows"
  )

  # Wrong number of columns in state
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = matrix(c(1, 2, 3), nrow = 1),
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "state.*must have 2 columns"
  )

  # Non-finite values in state
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = matrix(c(1, NA), nrow = 1),
      gamma = 1,
      spline_baseline = list(),
      init = NULL
    ),
    "state.*must contain finite values"
  )
})

# ==============================================================================
# Spline Baseline Validation
# ==============================================================================

test_that(".validate checks spline_baseline parameters", {
  # Invalid parameter names
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(invalid_param = 1),
      init = NULL
    ),
    "Invalid parameters in spline_baseline: invalid_param"
  )

  # Invalid degree
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      init = NULL,
      spline_baseline = list(degree = 10)
    ),
    "spline_baseline\\$degree must be a single integer between 1 and 5"
  )

  # Invalid n_knots
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      init = NULL,
      spline_baseline = list(n_knots = 50)
    ),
    "spline_baseline\\$n_knots must be a single integer between 0 and 20"
  )

  # Invalid boundary_knots
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      init = NULL,
      spline_baseline = list(boundary_knots = c(0, 1, 2))
    ),
    "spline_baseline\\$boundary_knots must be NULL or.*numeric vector"
  )

  # Invalid knot_placement
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      init = NULL,
      spline_baseline = list(knot_placement = "invalid")
    ),
    "spline_baseline\\$knot_placement must be one of: quantile, equal"
  )

  # Invalid boundary_knots
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      init = NULL,
      spline_baseline = list(boundary_knots = c(2, 1))
    ),
    "boundary_knots\\[1\\] must be less than"
  )
})

# ==============================================================================
# Init Parameter Validation
# ==============================================================================

test_that(".validate checks init parameter structure", {
  # Non-list init
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = "not a list"
    ),
    "init.*must be a list"
  )

  # Unknown init components
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(unknown_component = 1)
    ),
    "init.*unknown components 'unknown_component'"
  )
})

test_that(".validate checks init$coefficients structure", {
  # Non-list coefficients
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = "not a list")
    ),
    "init\\$coefficients.*must be a list"
  )

  # Unknown coefficient types
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(unknown_type = 1))
    ),
    "init\\$coefficients.*unknown types"
  )
})

test_that(".validate checks init$coefficients$baseline", {
  # Non-numeric baseline
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(baseline = "not numeric"))
    ),
    "baseline.*must be numeric"
  )

  # Non-finite baseline values
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(baseline = c(1, NA, 3)))
    ),
    "baseline.*must contain finite values"
  )

  # Wrong baseline length (requires spline_baseline setup)
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = rep(1:2, each = 2),
        time = rep(0:1, 2),
        v = 1:4
      ),
      survival_data = data.frame(id = 1:2, time = c(2, 3), status = c(1, 0)),
      state = NULL,
      gamma = 1,
      spline_baseline = list(degree = 3, n_knots = 2),
      init = list(coefficients = list(baseline = c(1, 2))) # Wrong length
    ),
    "baseline.*wrong length"
  )
})

test_that(".validate checks init$coefficients$hazard", {
  # Non-numeric hazard
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(hazard = "not numeric"))
    ),
    "hazard.*must be numeric"
  )

  # Non-finite hazard values
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(hazard = c(1, Inf)))
    ),
    "hazard.*must contain finite values"
  )

  # Wrong hazard length (too short)
  expect_error(
    .validate(
      longitudinal_formula = v ~ x1 + (1 | id),
      survival_formula = Surv(time, status) ~ w1,
      longitudinal_data = data.frame(
        id = c(1, 1),
        time = c(0, 1),
        v = c(1, 2),
        x1 = c(1, 2)
      ),
      survival_data = data.frame(id = 1, time = 1, status = 1, w1 = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(hazard = c(1)))
      # Should be length 3: biomarker, velocity, w1
    ),
    "hazard.*at least 2"
  )

  # Wrong hazard length (exact check)
  expect_error(
    .validate(
      longitudinal_formula = v ~ x1 + (1 | id),
      survival_formula = Surv(time, status) ~ w1,
      longitudinal_data = data.frame(
        id = c(1, 1),
        time = c(0, 1),
        v = c(1, 2),
        x1 = c(1, 2)
      ),
      survival_data = data.frame(id = 1, time = 1, status = 1, w1 = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(hazard = c(1, 2, 3, 4)))
      # Should be 3, not 4
    ),
    "hazard.*wrong length"
  )
})

test_that(".validate checks init$coefficients$longitudinal", {
  # Non-numeric longitudinal
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(longitudinal = "not numeric"))
    ),
    "longitudinal.*must be numeric"
  )

  # Non-finite longitudinal values
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(longitudinal = c(1, NA)))
    ),
    "longitudinal.*must contain finite values"
  )

  # Wrong longitudinal length
  expect_error(
    .validate(
      longitudinal_formula = v ~ x1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = c(1, 1),
        time = c(0, 1),
        v = c(1, 2),
        x1 = c(1, 2)
      ),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(longitudinal = c(1)))
      # Should be 2: fixed effects (intercept + x1), no ODE fixed effects
    ),
    "longitudinal.*wrong length.*expected 2"
  )
})

test_that(".validate checks init$coefficients$measurement_error_sd", {
  # Non-numeric measurement_error_sd
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(measurement_error_sd = "not numeric"))
    ),
    "measurement_error_sd.*must be"
  )

  # Non-finite measurement_error_sd
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(measurement_error_sd = Inf))
    ),
    "measurement_error_sd.*must be finite"
  )

  # Non-positive measurement_error_sd
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(measurement_error_sd = -1))
    ),
    "measurement_error_sd.*must be positive"
  )
})

test_that(".validate checks init$coefficients$random_effect_sigma", {
  # Non-numeric random_effect_sigma
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(random_effect_sigma = "not numeric"))
    ),
    "random_effect_sigma.*must be"
  )

  # Wrong dimension random_effect_sigma
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(random_effect_sigma = diag(4)))
      # Should be 1x1: only 1 random intercept (no ODE random effects specified)
    ),
    "random_effect_sigma.*wrong dim.*expected 1x1"
  )

  # Non-positive random_effect_sigma
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(coefficients = list(random_effect_sigma = matrix(0, 1, 1)))
      # 1x1 for 1 random intercept only (no ODE random effects)
    ),
    "random_effect_sigma.*must be positive"
  )
})

test_that(".validate checks init$configurations", {
  # Non-list configurations
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(configurations = "not a list")
    ),
    "init\\$configurations.*must be a list"
  )

  # Non-list baseline configuration
  expect_error(
    .validate(
      longitudinal_formula = v ~ 1 + (1 | id),
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      gamma = 1,
      spline_baseline = list(),
      init = list(configurations = list(baseline = "not a list"))
    ),
    "init\\$configurations\\$baseline.*must be a list"
  )
})

# ==============================================================================
# Valid Input Test
# ==============================================================================

test_that(".validate accepts valid inputs", {
  # Test passes without error or warnings
  result <- .validate(
    longitudinal_formula = v ~ 1 + (1 | id),
    survival_formula = Surv(time, status) ~ 1,
    longitudinal_data = data.frame(
      id = rep(1:3, each = 2),
      time = rep(0:1, 3),
      v = 1:6
    ),
    survival_data = data.frame(
      id = 1:3,
      time = c(2, 3, 4),
      status = c(0, 1, 1)
    ),
    state = NULL,
    gamma = 1,
    spline_baseline = list(),
    init = NULL
  )
  expect_true(is.null(result))
})
