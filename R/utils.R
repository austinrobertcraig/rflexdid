# Helpers for cohort construction, sample marking, and validation.
# Mirror the logic in stata/flexdid.ado lines 56-153.

#' @keywords internal
#' @noRd
validate_binary <- function(x, name) {
  vals <- sort(unique(x[!is.na(x)]))
  if (!identical(vals, c(0, 1)) && !identical(vals, c(0L, 1L))) {
    stop(
      sprintf(
        "Invalid treatment variable '%s' - tx must be binary with 0 for control observations and 1 for treated observations.",
        name
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
check_time_gaps <- function(time, has_usercohort) {
  u <- sort(unique(time))
  if (length(u) <= 1) return(invisible(TRUE))
  diffs <- diff(u)
  step <- min(diffs)
  has_gap <- any(diffs != step)
  if (has_gap && !has_usercohort) {
    stop(
      "Time variable has gaps. Specify usercohort= to correctly assign cohorts to groups that were first treated at times coincident with gaps in the data.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Compute the first-treatment cohort for each group.
#'
#' Replicates Stata: `egen long _Cohort = min(time/tx) by group;
#' replace _Cohort = 0 if _Cohort==.`. In Stata, `time/tx` evaluates to
#' `time` when `tx == 1` and missing when `tx == 0`. We achieve the same by
#' taking the minimum of `time` over rows where `tx == 1`, by group.
#'
#' @keywords internal
#' @noRd
make_cohort <- function(group, time, tx) {
  treated_time <- ifelse(tx == 1L, time, NA_real_)
  # tapply returns NaN if all NA; we replace with 0 to flag never-treated.
  cohort_per_group <- tapply(
    treated_time,
    group,
    function(z) {
      m <- suppressWarnings(min(z, na.rm = TRUE))
      if (is.infinite(m)) 0 else m
    }
  )
  out <- as.numeric(cohort_per_group[as.character(group)])
  out
}

#' Always-treated and never-treated checks.
#'
#' If the smallest positive cohort is at or before `min(time)`, the design
#' contains always-treated units; we abort. If the smallest cohort is > 0,
#' the design has no never-treated units; we recode the largest cohort to
#' "never-treated" and exclude all rows in time periods >= that cohort's
#' first-treatment year.
#'
#' Returns a list with the (possibly updated) cohort vector and a logical
#' `keep` mask aligned with the input.
#'
#' @keywords internal
#' @noRd
handle_treatment_status <- function(cohort, time, verbose = TRUE) {
  pos <- cohort > 0
  cmin <- if (any(pos)) min(cohort[pos]) else 0
  cmax <- if (any(pos)) max(cohort[pos]) else 0
  tmin <- min(time)

  if (cmin <= tmin && any(pos)) {
    stop(
      sprintf(
        "The first cohort is treated in or before the first time period (%s) observed in the data. This implies there are always-treated units. Remove always-treated units before using flexdid().",
        format(tmin)
      ),
      call. = FALSE
    )
  }

  keep <- rep(TRUE, length(cohort))
  if (length(unique(cohort)) > 0 && all(cohort > 0)) {
    if (verbose) {
      message(sprintf(
        "Note: There are no never-treated units. flexdid will define the last cohort (%s) as the never-treated group after dropping observations in all time periods in which the last cohort was treated.",
        format(cmax)
      ))
    }
    keep <- time < cmax
    cohort[cohort == cmax] <- 0
  }

  list(cohort = cohort, keep = keep)
}

#' Construct the lags-and-leads treatment indicator `_Tx`.
#'
#' Stata's recipe: start from `tx`, then for each treated cohort `c` set
#' `_Tx = 1` for observations with that cohort and `time < c - 1`. The base
#' period is `t = c - 1`, where `_Tx == 0`. The same recipe is used for both
#' specifications; the specification governs which interaction columns enter
#' the regression.
#'
#' @keywords internal
#' @noRd
make_tx_indicator <- function(tx, cohort, time) {
  out <- as.integer(tx)
  treated_cohorts <- sort(unique(cohort[cohort > 0]))
  for (c in treated_cohorts) {
    sel <- cohort == c & time < (c - 1)
    out[sel] <- 1L
  }
  out
}

#' Event time: time - cohort for treated, -1 for controls.
#' @keywords internal
#' @noRd
make_eventtime <- function(cohort, time) {
  ifelse(cohort > 0, time - cohort, -1L)
}

#' Resolve a `subset` argument given as a (possibly NSE) expression.
#'
#' Accepts a logical vector, an expression, or NULL. Always returns a logical
#' vector of length nrow(data).
#'
#' @keywords internal
#' @noRd
resolve_subset <- function(subset_expr, data, env) {
  if (is.null(subset_expr)) return(rep(TRUE, nrow(data)))
  if (is.logical(subset_expr) && length(subset_expr) == nrow(data)) {
    return(subset_expr & !is.na(subset_expr))
  }
  if (is.call(subset_expr) || is.name(subset_expr) || is.language(subset_expr)) {
    val <- eval(subset_expr, envir = data, enclos = env)
    if (!is.logical(val) || length(val) != nrow(data)) {
      stop("`subset` must evaluate to a logical vector of length nrow(data).",
           call. = FALSE)
    }
    return(val & !is.na(val))
  }
  stop("`subset` must be NULL, a logical vector, or an expression.",
       call. = FALSE)
}

#' Pull a column from `data` by name, with a friendly error.
#'
#' @keywords internal
#' @noRd
get_col <- function(data, name, role) {
  if (is.null(name)) return(NULL)
  if (!name %in% names(data)) {
    stop(sprintf("Column '%s' (role: %s) not found in data.", name, role),
         call. = FALSE)
  }
  data[[name]]
}
