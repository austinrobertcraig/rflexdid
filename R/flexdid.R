#' Flexible difference-in-differences with staggered treatment timing
#'
#' R port of the Stata `flexdid` command (Deb, Norton, Wooldridge & Zabel
#' 2025). Estimates a flexible OLS regression with cohort-by-group-by-time
#' treatment interactions, then aggregates cell-level effects to ATETs via
#' [atet()].
#'
#' @param formula Formula; two-sided, e.g. `y ~ x1 + x2`. The left-hand
#'   side is the outcome and the right-hand side lists the covariates that
#'   should be interacted with treatment, group, and time indicators. Pass
#'   `y ~ 1` for no interacted covariates.
#' @param data Data frame.
#' @param tx Character; column name of the binary treatment indicator
#'   (1 = treated this period, 0 = otherwise). Required.
#' @param group Character; column name of the group variable used for
#'   group fixed effects and (by default) the level at which ATETs are
#'   estimated. Required.
#' @param time Character; column name of the integer time variable.
#'   Required. Must be equally spaced unless `usercohort` is supplied.
#' @param specification Character; one of `"lagsonly"` (default) or
#'   `"lagsandleads"`.
#' @param xnotinteracted Optional additive controls. Character vector of
#'   column names, or a one-sided formula. Must be disjoint from
#'   `formula`'s covariates.
#' @param usercohort Character; optional column name with a user-supplied
#'   cohort variable. Overrides the internal cohort calculation. Useful
#'   when some time periods are missing.
#' @param weights Character; optional column name of observation weights
#'   (treated as pweights, matching the Stata default for clustered
#'   designs).
#' @param vcov Character; one of `"cluster"` (default) or `"robust"`.
#' @param cluster Character; optional column name for the cluster
#'   variable. Defaults to `group` when `vcov = "cluster"`.
#' @param subset Optional logical vector or unevaluated expression used to
#'   restrict the sample (analogous to Stata's `if`).
#' @param verbose Logical; if `TRUE`, print the underlying regression
#'   coefficient table after fitting.
#'
#' @return An object of class `flexdid` (a list); see [print.flexdid()] and
#'   [summary.flexdid()] for human-readable output and [atet()] for ATET
#'   aggregation.
#'
#' @references
#' Deb, P., Norton, E. C., Wooldridge, J. M., Zabel, J. E. (2025), "A
#' Flexible, Heterogeneous Treatment Effects Difference-in-Differences
#' Estimator for Repeated Cross-Sections", NBER Working Paper No. 33026.
#'
#' @export
flexdid <- function(formula,
                    data,
                    tx,
                    group,
                    time,
                    specification = c("lagsonly", "lagsandleads"),
                    xnotinteracted = NULL,
                    usercohort = NULL,
                    weights = NULL,
                    vcov = c("cluster", "robust"),
                    cluster = NULL,
                    subset = NULL,
                    verbose = FALSE) {

  cl <- match.call()
  specification <- match.arg(specification)
  vcov <- match.arg(vcov)

  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula such as `y ~ x1 + x2`.",
         call. = FALSE)
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  for (arg_nm in c("tx", "group", "time")) {
    val <- get(arg_nm)
    if (!is.character(val) || length(val) != 1L) {
      stop(sprintf("`%s` must be a single column name (character).", arg_nm),
           call. = FALSE)
    }
    if (!val %in% names(data)) {
      stop(sprintf("Column '%s' (role: %s) not found in data.", val, arg_nm),
           call. = FALSE)
    }
  }

  # --- Resolve subset
  subset_expr <- substitute(subset)
  if (is.null(subset_expr) || identical(subset_expr, quote(NULL))) {
    use <- rep(TRUE, nrow(data))
  } else {
    parent <- parent.frame()
    use <- resolve_subset(eval(subset_expr, parent), data, parent)
  }
  if (!any(use)) {
    stop("`subset` excludes all observations.", call. = FALSE)
  }

  # --- Pull required columns
  tx_raw <- data[[tx]]
  group_vec <- data[[group]]
  time_vec  <- data[[time]]

  tx_vec <- coerce_tx(tx_raw, tx)

  if (anyNA(tx_vec[use]) || anyNA(group_vec[use]) || anyNA(time_vec[use])) {
    use <- use & !is.na(tx_vec) & !is.na(group_vec) & !is.na(time_vec)
  }

  validate_binary(tx_vec[use], tx)

  if (!is.numeric(time_vec) && !is.integer(time_vec)) {
    stop(sprintf("Column '%s' (time) must be numeric/integer.", time), call. = FALSE)
  }

  has_user_cohort <- !is.null(usercohort)
  check_time_gaps(time_vec[use], has_user_cohort)

  # --- Weights
  w <- NULL
  if (!is.null(weights)) {
    if (!is.character(weights) || length(weights) != 1L) {
      stop("`weights` must be NULL or a single column name.", call. = FALSE)
    }
    w_vec <- data[[weights]]
    if (is.null(w_vec)) {
      stop(sprintf("Weights column '%s' not found in data.", weights),
           call. = FALSE)
    }
    if (anyNA(w_vec[use])) {
      use <- use & !is.na(w_vec)
    }
    w <- as.numeric(w_vec)
    if (any(w[use] < 0)) stop("Weights must be non-negative.", call. = FALSE)
  }

  # --- Cluster column
  if (vcov == "cluster") {
    if (is.null(cluster)) cluster <- group
    if (!is.character(cluster) || length(cluster) != 1L) {
      stop("`cluster` must be NULL or a single column name.", call. = FALSE)
    }
    if (!cluster %in% names(data)) {
      stop(sprintf("Cluster column '%s' not found in data.", cluster),
           call. = FALSE)
    }
    cluster_vec <- data[[cluster]]
    if (anyNA(cluster_vec[use])) {
      use <- use & !is.na(cluster_vec)
    }
  } else {
    cluster_vec <- NULL
  }

  # --- Covariate NA handling
  # model.matrix() in build_design() silently drops rows with NA covariates,
  # which would misalign row indices used to slice xvars_mm. Pre-filter here
  # so data_in already excludes those rows.
  covar_cols <- all.vars(formula[[3L]])
  if (!is.null(xnotinteracted)) {
    xni_vars <- if (is.character(xnotinteracted)) {
      xnotinteracted
    } else if (inherits(xnotinteracted, "formula")) {
      all.vars(xnotinteracted)
    } else {
      character(0)
    }
    covar_cols <- unique(c(covar_cols, xni_vars))
  }
  covar_cols <- intersect(covar_cols, names(data))
  if (length(covar_cols) > 0L) {
    covar_na <- Reduce(`|`, lapply(data[covar_cols], is.na))
    n_drop <- sum(use & covar_na)
    if (n_drop > 0L) {
      message(sprintf(
        "Note: dropped %d row(s) with NAs in covariates (%s).",
        n_drop, paste(covar_cols, collapse = ", ")
      ))
      use <- use & !covar_na
      if (!any(use)) {
        stop("All observations have NA covariates; nothing left to fit.",
             call. = FALSE)
      }
    }
  }

  # --- Outcome NA handling
  # Match the auto-filter behavior used for tx/group/time/cluster/weights
  # and covariates: drop rows with NA outcome and inform the user.
  yvar_name <- as.character(formula[[2L]])
  if (!yvar_name %in% names(data)) {
    stop(sprintf("Outcome variable '%s' not found in data.", yvar_name),
         call. = FALSE)
  }
  y_na <- is.na(data[[yvar_name]])
  n_drop_y <- sum(use & y_na)
  if (n_drop_y > 0L) {
    message(sprintf(
      "Note: dropped %d row(s) with NAs in outcome variable '%s'.",
      n_drop_y, yvar_name
    ))
    use <- use & !y_na
    if (!any(use)) {
      stop(sprintf("All observations have NA outcome ('%s'); nothing left to fit.",
                   yvar_name), call. = FALSE)
    }
  }

  # --- Cohort construction
  if (has_user_cohort) {
    if (!usercohort %in% names(data)) {
      stop(sprintf("Column '%s' (usercohort) not found in data.", usercohort),
           call. = FALSE)
    }
    cohort_full <- as.numeric(data[[usercohort]])
  } else {
    cohort_full <- rep(NA_real_, nrow(data))
    cohort_full[use] <- make_cohort(group_vec[use], time_vec[use], tx_vec[use])
  }

  # always-/never-treated handling, restricted to the in-sample rows
  cohort_in <- cohort_full[use]
  time_in   <- time_vec[use]
  ts <- handle_treatment_status(cohort_in, time_in, verbose = TRUE)
  cohort_in <- ts$cohort
  if (!all(ts$keep)) {
    rows_in <- which(use)[ts$keep]
    use2 <- rep(FALSE, nrow(data))
    use2[rows_in] <- TRUE
    use <- use2
    cohort_full <- rep(NA_real_, nrow(data))
    cohort_full[use] <- cohort_in
    time_in <- time_vec[use]
  } else {
    cohort_full[use] <- cohort_in
  }

  # --- Build _Tx (lags-and-leads indicator) and event time
  tx_indicator_in <- make_tx_indicator(tx_vec[use], cohort_in, time_vec[use])
  eventtime_in <- make_eventtime(cohort_in, time_vec[use])

  tx_indicator_full <- rep(NA_integer_, nrow(data))
  tx_indicator_full[use] <- tx_indicator_in
  eventtime_full <- rep(NA_integer_, nrow(data))
  eventtime_full[use] <- eventtime_in

  # Sanity warning matching the Stata "treatment varies within group/time" note
  cn <- sum(cohort_in > 0 & eventtime_in >= 0)
  tn <- sum(tx_indicator_in == 1 & eventtime_in >= 0)
  if (cn != tn) {
    message("Note: treatment varies within group and time. This is outside the formal scope of the standard model specification; interpret results accordingly.")
  }

  # --- Build the design matrix on the in-sample rows
  data_in <- data[use, , drop = FALSE]
  rownames(data_in) <- NULL
  y_in <- as.numeric(data_in[[yvar_name]])

  des <- build_design(
    data           = data_in,
    formula        = formula,
    group          = group_vec[use],
    time           = time_vec[use],
    cohort         = cohort_in,
    tx_indicator   = tx_indicator_in,
    xnotinteracted = xnotinteracted,
    specification  = specification
  )

  X_sparse <- des$X
  n <- nrow(X_sparse)

  # --- Fit on the sparse design (sparse Cholesky with rank-deficient fallback
  # to dense lm.fit; see R/sparse_ols.R).
  w_in <- if (is.null(w)) rep(1, n) else as.numeric(w[use])
  fit <- sparse_ols(X_sparse, y_in, w_in)

  beta <- fit$beta
  pivot_keep <- fit$pivot_keep
  bread_kept <- fit$bread
  rank <- sum(pivot_keep)
  k_total <- length(beta)
  # Stata / lm.fit convention: dropped columns get coefficient 0 (not NA) so
  # downstream uses (fitted values, ATET aggregation) sum cleanly. pivot_keep
  # still records which columns the OLS solution actually identifies.
  beta[!pivot_keep] <- 0
  fitted <- as.numeric(X_sparse %*% beta)
  residuals <- y_in - fitted

  rss <- sum(w_in * residuals^2)
  ybar <- if (is.null(w)) mean(y_in) else stats::weighted.mean(y_in, w_in)
  tss <- sum(w_in * (y_in - ybar)^2)
  mss <- tss - rss
  df_residual <- n - rank
  df_model <- rank - 1L  # subtract intercept
  r2 <- if (tss > 0) 1 - rss / tss else NA_real_
  r2_adj <- if (df_residual > 0 && tss > 0) {
    1 - (rss / df_residual) / (tss / (n - 1))
  } else NA_real_

  # --- VCE for beta. Pass the bread we already computed in sparse_ols() so
  # compute_vcov() does not refactor (X' W X).
  V <- compute_vcov(
    X = X_sparse,
    residuals = residuals,
    weights = w_in,
    pivot_keep = pivot_keep,
    rank = rank,
    vcov_type = vcov,
    cluster = if (vcov == "cluster") cluster_vec[use] else NULL,
    bread = bread_kept
  )
  # Pad V with zeros for pivoted-out columns so dim(V) = k x k.
  V_full <- matrix(0, k_total, k_total)
  V_full[pivot_keep, pivot_keep] <- V
  rownames(V_full) <- colnames(V_full) <- des$colnames

  # --- F statistic on all non-intercept terms (mirrors Stata's testparm).
  # The constraint matrix R is a row selector, so R %*% beta = beta[non_int]
  # and R %*% V %*% t(R) = V[non_int, non_int]; doing the matrix-multiplies
  # explicitly is wasteful (p x p dense) for large designs.
  non_int <- setdiff(seq_along(beta), des$block_index$intercept)
  non_int <- non_int[pivot_keep[non_int]]
  if (length(non_int) > 0L) {
    Rb <- beta[non_int]
    RVR <- V_full[non_int, non_int, drop = FALSE]
    inv <- tryCatch(solve(RVR), error = function(e) NULL)
    F_stat <- if (!is.null(inv)) {
      as.numeric(t(Rb) %*% inv %*% Rb) / length(non_int)
    } else NA_real_
  } else {
    F_stat <- NA_real_
  }

  names(beta) <- des$colnames

  out <- list(
    coefficients   = beta,
    pivot_keep     = pivot_keep,
    vcov           = V_full,
    vcov_type      = vcov,
    cluster_var    = if (vcov == "cluster") cluster else NULL,
    cluster_vec    = if (vcov == "cluster") cluster_vec[use] else NULL,
    n_clusters     = if (vcov == "cluster") length(unique(cluster_vec[use])) else NA_integer_,
    fitted         = fitted,
    residuals      = residuals,
    weights        = w_in,
    has_weights    = !is.null(w),
    rank           = rank,
    df_residual    = df_residual,
    df_model       = df_model,
    nobs           = n,
    rss            = rss,
    mss            = mss,
    r2             = r2,
    r2_adj         = r2_adj,
    f              = F_stat,
    f_df1          = length(non_int),
    f_df2          = df_residual,
    call           = cl,
    formula        = formula,
    specification  = specification,
    tx_var         = tx,
    group_var      = group,
    time_var       = time,
    usercohort     = usercohort,
    yvar           = yvar_name,
    # Sample-aligned vectors (length = n)
    sample_index   = which(use),
    cohort         = cohort_in,
    tx_indicator   = tx_indicator_in,
    eventtime      = eventtime_in,
    group_in       = group_vec[use],
    time_in        = time_vec[use],
    y_in           = y_in,
    # Design metadata
    design         = des,
    X              = X_sparse,
    bread          = bread_kept,
    data_in        = data_in
  )
  class(out) <- "flexdid"

  if (isTRUE(verbose)) {
    cat("Estimating", specification, "regression parameters\n")
    print(summary(out))
  }

  out
}
