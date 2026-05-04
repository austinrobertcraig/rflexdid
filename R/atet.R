# Post-estimation aggregation of cell-level effects to ATETs.
# Mirrors the Stata postestimation suite in stata/flexdid_atet.ado:
#   _ATET_Overall, _ATET_Byget, _ATET_Byexposure, _ATET_Bycalendar,
#   _ATET_Bycohort, _ATET_Bygroup.
#
# The core identity: for an observation i with treatment cell (g, t, kind),
#   TE_i = beta_{cell} + sum_k beta_{cell,k} * x_{i,k}.
# So ATET on a subpopulation S is (weighted) mean of TE_i over S, which is
# linear in beta:  ATET = c'beta with c the (weighted) mean of the per-obs
# gradient. For multiple ATETs (e.g., byexposure) we stack the c's into a
# matrix C and write ATET = C beta.
#
# Standard errors follow Stata's vce(unconditional). The influence function
# for obs i, ATET-level l is
#
#   IF_i^l = bar(c)_l' (X'WX)^{-1} w_i X_i' u_i
#            + (w_i / W_{S_l}) * 1[i in S_l] * (TE_i - ATET_l)
#
# and the variance is the sandwich on IF (cluster-summed if vcov_type =
# "cluster"), with the same small-sample factors as for the regression VCE.

#' Aggregate cell-level effects to ATETs.
#'
#' Postestimation function for [flexdid()], analogous to Stata's
#' `estat atet`.
#'
#' @param model A `flexdid` object.
#' @param type One of `"overall"`, `"byget"`, `"byexposure"`, `"bycalendar"`,
#'   `"bycohort"`, or `"bygroup"`.
#' @param values Optional numeric vector restricting which levels are
#'   reported. Semantics differ by `type` (see Details).
#' @param for_expr Optional one-sided formula restricting the subpopulation,
#'   e.g. `~ female == 1`. Multiple conditions can be combined with `&` or `|`
#'   (e.g. `~ female == 1 & age > 10`). A `quote()`-style unevaluated
#'   expression is also accepted. Evaluated in the modeling data.
#' @param dydx Logical; included for API parity with Stata. The R
#'   implementation uses the same gradient-based contrast in either case,
#'   so this flag is ignored.
#' @param aggregationweight Either `"obslevel"` (default) or `"grouplevel"`.
#'   Group-level weighting reweights observations by `(n_g/T)/n_gt`, then
#'   renormalizes the weights to average 1 on the subpopulation.
#' @param test Optional Wald test to append. `"zero"` tests H0: all ATETs = 0;
#'   `"equal"` tests H0: all ATETs are equal; `"pretrends"` tests H0: all
#'   pre-treatment ATETs = 0 (parallel-trends pre-test, only valid for
#'   `type = "byexposure"`).
#' @param level Numeric in (0,1). Confidence level for printed CIs.
#' @return An object of class `flexdid_atet`.
#' @export
atet <- function(model,
                 type = c("overall", "byexposure", "bycalendar",
                          "bycohort", "bygroup", "byget"),
                 values = NULL,
                 for_expr = NULL,
                 dydx = FALSE,
                 aggregationweight = c("obslevel", "grouplevel"),
                 test = NULL,
                 level = 0.95) {
  if (!inherits(model, "flexdid")) {
    stop("`model` must be a flexdid object.", call. = FALSE)
  }
  type <- match.arg(type)
  aggregationweight <- match.arg(aggregationweight)
  if (!is.null(test)) {
    if (!test %in% c("zero", "equal", "pretrends")) {
      stop('`test` must be "zero", "equal", or "pretrends".', call. = FALSE)
    }
  }

  n <- model$nobs
  cohort <- model$cohort
  time <- model$time_in
  group <- model$group_in
  eventtime <- model$eventtime
  tx_indicator <- model$tx_indicator
  ttx <- tx_indicator
  # Stata's `replace ttx = 1 if _Cohort > 0 & eventtime == -1`. Adds the base
  # period of treated cohorts to the subpopulation, per Stata's overall /
  # byexposure / byget specifications.
  ttx[cohort > 0 & eventtime == -1L] <- 1L

  # for_expr restriction
  for_mask <- rep(TRUE, n)
  if (!is.null(for_expr)) {
    expr <- if (inherits(for_expr, "formula")) for_expr[[2L]] else for_expr
    val <- eval(expr, envir = model$data_in)
    if (!is.logical(val) || length(val) != n) {
      stop("`for_expr` must evaluate to a logical vector of length nobs.",
           call. = FALSE)
    }
    for_mask <- val & !is.na(val)
  }

  # ---- Determine subpopulation and per-obs level
  build <- build_levels(type, values, n, ttx, tx_indicator, cohort, time,
                        group, eventtime, for_mask)
  subpop <- build$subpop
  level_idx <- build$level_idx        # integer in 1..L, NA outside subpop
  level_labels <- build$level_labels  # character, length L
  level_values <- build$level_values  # list, length L (for column names of result)

  L <- length(level_labels)
  if (L == 0L) {
    stop("No observations match the requested ATET specification.", call. = FALSE)
  }

  n_sub <- sum(subpop)
  n_full <- model$nobs

  # ---- Per-observation TE and gradient row.
  # Each treated, non-base obs maps to exactly one cell in the design; the
  # gradient there is 1 at col_indicator and x_{i,k} at col_x_offset+(k-1).
  cells <- model$design$cells
  k <- model$design$n_xvars
  X <- model$X
  beta <- model$coefficients

  obs_cell <- match_cell(group, time, cohort, tx_indicator, cells)
  TE <- numeric(n)
  has_cell <- !is.na(obs_cell)
  if (any(has_cell)) {
    rows_h <- which(has_cell)
    col_ind <- cells$col_indicator[obs_cell[rows_h]]
    TE[rows_h] <- beta[col_ind]
    if (k > 0L) {
      offsets <- cells$col_x_offset[obs_cell[rows_h]]
      xvars_mm <- model$design$xvars_matrix
      contrib <- numeric(length(rows_h))
      for (kk in seq_len(k)) {
        cols_kk <- offsets + (kk - 1L)
        contrib <- contrib + beta[cols_kk] * xvars_mm[rows_h, kk]
      }
      TE[rows_h] <- TE[rows_h] + contrib
    }
  }

  # ---- Aggregation weights (per the Stata `aggregationweight` option)
  w_obs <- model$weights
  if (aggregationweight == "grouplevel") {
    aw <- compute_grouplevel_weights(group, time, cohort)
    # Normalize on the subpop so weights average 1 there
    norm_subpop <- subpop & !is.na(aw)
    aw_mean <- mean(aw[norm_subpop])
    if (aw_mean == 0 || !is.finite(aw_mean)) aw_mean <- 1
    aw <- aw / aw_mean
    aw[is.na(aw)] <- 0
    w_eff <- w_obs * aw
  } else {
    w_eff <- w_obs
  }

  # ---- Per-level mean weight, ATET, and bar(c).
  W_per_level <- numeric(L)
  ATET <- numeric(L)
  for (l in seq_len(L)) {
    in_l <- !is.na(level_idx) & level_idx == l & subpop
    W_per_level[l] <- sum(w_eff[in_l])
    if (W_per_level[l] > 0) {
      ATET[l] <- sum(w_eff[in_l] * TE[in_l]) / W_per_level[l]
    } else {
      ATET[l] <- NA_real_
    }
  }

  # ---- Build C: L x p_kept (only kept columns; pivoted-out cols are zeroed)
  pivot_keep <- model$pivot_keep
  p_total <- length(beta)
  p_kept <- sum(pivot_keep)
  C <- matrix(0, nrow = L, ncol = p_kept)
  # Mapping from full column index to kept column index
  kept_index <- integer(p_total)
  kept_index[pivot_keep] <- seq_len(p_kept)

  for (l in seq_len(L)) {
    in_l <- !is.na(level_idx) & level_idx == l & subpop & has_cell
    if (W_per_level[l] <= 0 || !any(in_l)) next
    rows_l <- which(in_l)
    w_l <- w_eff[rows_l]
    cell_l <- obs_cell[rows_l]
    col_ind_full <- cells$col_indicator[cell_l]
    # contribution: w_i / W_{S_l} at col_indicator
    by_col <- tapply(w_l, col_ind_full, sum)
    cols <- as.integer(names(by_col))
    cols_k <- kept_index[cols]
    nz <- cols_k > 0
    C[l, cols_k[nz]] <- C[l, cols_k[nz]] + as.numeric(by_col)[nz] / W_per_level[l]

    if (k > 0L) {
      offsets <- cells$col_x_offset[cell_l]
      xvars_mm <- model$design$xvars_matrix
      for (kk in seq_len(k)) {
        cols_kk_full <- offsets + (kk - 1L)
        x_kk <- xvars_mm[rows_l, kk]
        # sum w_i * x_i,kk by column (in case multiple obs share a column,
        # which they do whenever multiple obs land in the same cell)
        contrib_by_col <- tapply(w_l * x_kk, cols_kk_full, sum)
        cols_kk_unique <- as.integer(names(contrib_by_col))
        cols_kk_kept <- kept_index[cols_kk_unique]
        nz <- cols_kk_kept > 0
        C[l, cols_kk_kept[nz]] <- C[l, cols_kk_kept[nz]] +
          as.numeric(contrib_by_col)[nz] / W_per_level[l]
      }
    }
  }

  # ---- Influence functions.
  # Score (per observation, on kept cols): w_i X_i u_i. X_kept stays sparse;
  # only the final IF_param product is densified.
  X_kept <- X[, pivot_keep, drop = FALSE]
  resid <- model$residuals
  score <- X_kept * (resid * w_obs)  # n x p_kept; rows are w_i X_i u_i

  # Bread = (X' W X)^{-1}, computed once in flexdid() and stored on the model.
  # Older flexdid objects (or callers that hand-construct one) may lack it; in
  # that case recompute here as a fallback.
  bread <- model$bread
  if (is.null(bread)) {
    if (all(w_obs == 1)) {
      bread_inv <- crossprod(as.matrix(X_kept))
    } else {
      bread_inv <- crossprod(as.matrix(X_kept) * sqrt(w_obs))
    }
    bread <- chol2inv(chol(bread_inv))
  }

  # IF parameter piece: (n x L) = score (n x p_kept) %*% bread %*% t(C)
  IF_param <- as.matrix(score %*% (bread %*% t(C)))

  # IF subpop piece: (n x L) where IF[i, l] = (w_i / W_{S_l}) * 1[i in S_l] * (TE_i - ATET_l)
  IF_subpop <- matrix(0, nrow = n, ncol = L)
  for (l in seq_len(L)) {
    in_l <- !is.na(level_idx) & level_idx == l & subpop
    if (W_per_level[l] <= 0) next
    IF_subpop[in_l, l] <- (w_eff[in_l] / W_per_level[l]) * (TE[in_l] - ATET[l])
  }

  IF <- IF_param + IF_subpop

  # ---- Variance via cluster/robust sandwich on IF.
  rank <- model$rank
  if (model$vcov_type == "cluster") {
    cluster_f <- as.factor(model$cluster_vec)
    G <- nlevels(cluster_f)
    Psi_g <- rowsum(IF, group = cluster_f, reorder = FALSE)
    meat <- crossprod(Psi_g)
    adj <- (G / (G - 1)) * ((n - 1) / (n - rank))
  } else {
    meat <- crossprod(IF)
    adj <- n / (n - rank)
  }
  V_atet <- adj * meat
  V_atet <- (V_atet + t(V_atet)) / 2

  rownames(V_atet) <- colnames(V_atet) <- level_labels
  names(ATET) <- level_labels

  # ---- Build the printed/coef table
  diag_v <- diag(V_atet)
  diag_v[diag_v < 0 & abs(diag_v) < 1e-12] <- 0  # tidy tiny negative numerical noise
  se <- sqrt(diag_v)
  df <- model$df_residual
  # Levels with zero SE (e.g. base period of byexposure in lagsonly) get
  # t / p / CI = NA so the table prints cleanly.
  # Degenerate SE: technically > 0 but below sqrt(.Machine$double.eps) ~ 1.49e-8.
  # Indicates a near-singular variance matrix (e.g., singleton or near-singleton
  # clusters in a given exposure period). Treated as NA, same as the base period.
  degenerate_se <- is.finite(se) & se > 0 & se < sqrt(.Machine$double.eps)
  zero_se <- !is.finite(se) | se < sqrt(.Machine$double.eps)
  degenerate_note <- if (any(degenerate_se)) {
    sprintf(
      "Note: %d level(s) have near-zero standard errors (SE < sqrt(.Machine$double.eps) ≈ 1.49e-8), likely due to singleton or near-singleton clusters. t, p, and CIs set to NA for: %s.",
      sum(degenerate_se),
      paste(level_labels[degenerate_se], collapse = ", ")
    )
  } else NULL
  tval <- ifelse(zero_se, NA_real_, ATET / se)
  pval <- ifelse(zero_se, NA_real_, 2 * stats::pt(-abs(tval), df = df))
  qcrit <- stats::qt(1 - (1 - level) / 2, df)
  ci_lo <- ifelse(zero_se, NA_real_, ATET - qcrit * se)
  ci_hi <- ifelse(zero_se, NA_real_, ATET + qcrit * se)
  tidy_table <- cbind(
    Estimate = ATET,
    `Std. Error` = se,
    `t value` = tval,
    `Pr(>|t|)` = pval,
    `[CI lo]` = ci_lo,
    `[CI hi]` = ci_hi
  )
  rownames(tidy_table) <- level_labels

  # ---- Optional Wald test (skip zero-SE rows; e.g. eventtime=-1 base period)
  test_result <- NULL
  if (!is.null(test)) {
    effective_test <- test
    if (identical(test, "pretrends")) {
      if (!identical(type, "byexposure")) {
        warning("`test = 'pretrends'` is only applicable for type = 'byexposure'; ",
                "falling back to `test = 'zero'`.")
        effective_test <- "zero"
      }
    }
    keep_l <- !zero_se
    if (identical(effective_test, "pretrends")) {
      keep_l_pre <- keep_l & (as.numeric(level_labels) < 0)
      if (sum(keep_l_pre) >= 1L) {
        test_result <- atet_wald(ATET[keep_l_pre],
                                 V_atet[keep_l_pre, keep_l_pre, drop = FALSE],
                                 df_resid = df, type_of_test = "zero")
        test_result$title <- "Pre-test of parallel trends assumption"
        test_result$h0    <- "All pre-treatment effects are equal to zero"
      }
    } else if (sum(keep_l) >= 1L) {
      test_result <- atet_wald(ATET[keep_l], V_atet[keep_l, keep_l, drop = FALSE],
                               df_resid = df, type_of_test = effective_test)
    }
  }

  out <- list(
    estimate           = ATET,
    vcov               = V_atet,
    df_residual        = df,
    n_full             = n_full,
    n_sub              = n_sub,
    type               = type,
    values             = values,
    for_expr           = for_expr,
    aggregation_weight = aggregationweight,
    test_result        = test_result,
    degenerate_note    = degenerate_note,
    level              = level,
    tidy_table         = tidy_table,
    level_labels       = level_labels,
    level_values       = level_values,
    yvar               = model$yvar,
    C                  = C,
    IF                 = IF
  )
  class(out) <- "flexdid_atet"
  out
}


# ---- Internal helpers ------------------------------------------------------

#' @keywords internal
#' @noRd
build_levels <- function(type, values, n, ttx, tx_indicator, cohort, time,
                         group, eventtime, for_mask) {
  # Returns subpop (logical n), level_idx (integer or NA, length n),
  # level_labels (character L), level_values (list of typed values, length L).

  if (type == "overall") {
    if (is.null(values)) {
      base_subpop <- ttx == 1L & eventtime >= 0L
      vals <- sort(unique(eventtime[base_subpop]))
    } else {
      vals <- values
      base_subpop <- ttx == 1L & eventtime %in% vals
    }
    subpop <- base_subpop & for_mask
    level_idx <- ifelse(subpop, 1L, NA_integer_)
    level_labels <- "Overall"
    level_values <- list(eventtimes = vals)
  } else if (type == "byexposure") {
    base_subpop <- ttx == 1L
    if (is.null(values)) {
      vals <- sort(unique(eventtime[base_subpop & for_mask]))
    } else {
      vals <- values
      base_subpop <- base_subpop & eventtime %in% vals
    }
    subpop <- base_subpop & for_mask
    level_idx <- match(eventtime, vals)
    level_idx[!subpop] <- NA_integer_
    level_labels <- as.character(vals)
    level_values <- as.list(vals)
  } else if (type == "bycalendar") {
    base_subpop <- tx_indicator == 1L
    if (is.null(values)) {
      vals <- sort(unique(time[base_subpop & for_mask]))
    } else {
      vals <- values
      base_subpop <- base_subpop & time %in% vals
    }
    subpop <- base_subpop & for_mask
    level_idx <- match(time, vals)
    level_idx[!subpop] <- NA_integer_
    level_labels <- as.character(vals)
    level_values <- as.list(vals)
  } else if (type == "bycohort") {
    base_subpop <- tx_indicator == 1L & eventtime >= 0L
    if (is.null(values)) {
      vals <- sort(unique(cohort[base_subpop & for_mask]))
    } else {
      vals <- values
      base_subpop <- base_subpop & cohort %in% vals
    }
    subpop <- base_subpop & for_mask
    level_idx <- match(cohort, vals)
    level_idx[!subpop] <- NA_integer_
    level_labels <- as.character(vals)
    level_values <- as.list(vals)
  } else if (type == "bygroup") {
    base_subpop <- tx_indicator == 1L & eventtime >= 0L
    if (is.null(values)) {
      vals <- sort(unique(group[base_subpop & for_mask]))
    } else {
      vals <- values
      base_subpop <- base_subpop & group %in% vals
    }
    subpop <- base_subpop & for_mask
    level_idx <- match(group, vals)
    level_idx[!subpop] <- NA_integer_
    level_labels <- as.character(vals)
    level_values <- as.list(vals)
  } else if (type == "byget") {
    base_subpop <- ttx == 1L
    if (is.null(values)) {
      sub <- base_subpop & for_mask
      pairs <- unique(data.frame(g = group[sub], et = eventtime[sub]))
      pairs <- pairs[order(pairs$g, pairs$et), , drop = FALSE]
    } else {
      base_subpop <- base_subpop & eventtime %in% values
      sub <- base_subpop & for_mask
      pairs <- unique(data.frame(g = group[sub], et = eventtime[sub]))
      pairs <- pairs[order(pairs$g, pairs$et), , drop = FALSE]
    }
    subpop <- base_subpop & for_mask
    pair_key_obs <- paste(group, eventtime, sep = "|")
    pair_key_lev <- paste(pairs$g, pairs$et, sep = "|")
    level_idx <- match(pair_key_obs, pair_key_lev)
    level_idx[!subpop] <- NA_integer_
    level_labels <- sprintf("g%s|et%s", pairs$g, pairs$et)
    level_values <- as.list(seq_len(nrow(pairs)))
  }

  list(
    subpop = subpop,
    level_idx = level_idx,
    level_labels = level_labels,
    level_values = level_values
  )
}

#' Match each observation to a row of `design$cells` (NA if none).
#'
#' Cell membership requires `tx_indicator == 1` (excluding base period
#' `t = c-1`) and the obs to be in a treated cohort. The cell kind is "lag"
#' for `t >= cohort` and "lead" for `t <= cohort - 2`. We compose a key
#' "kind|g|t" and match against the precomputed cells.
#' @keywords internal
#' @noRd
match_cell <- function(group, time, cohort, tx_indicator, cells) {
  n <- length(group)
  out <- rep(NA_integer_, n)
  active <- tx_indicator == 1L & cohort > 0L
  if (!any(active)) return(out)
  kind <- ifelse(time[active] >= cohort[active], "lag", "lead")
  key <- paste(kind, group[active], time[active], sep = "|")
  cell_key <- paste(cells$kind, cells$g, cells$t, sep = "|")
  out[active] <- match(key, cell_key)
  out
}

#' Group-level aggregation weights (mirrors the Stata `aggwt` calculation).
#' @keywords internal
#' @noRd
compute_grouplevel_weights <- function(group, time, cohort) {
  n <- length(group)
  T_total <- length(unique(time))
  treated <- cohort > 0L
  ngt <- numeric(n)
  ng <- numeric(n)
  if (any(treated)) {
    gt_key <- paste(group[treated], time[treated], sep = "|")
    g_key  <- as.character(group[treated])
    gt_count <- table(gt_key)
    g_count  <- table(g_key)
    ngt[treated] <- as.numeric(gt_count[gt_key])
    ng[treated]  <- as.numeric(g_count[g_key])
  }
  aw <- (ng / T_total) / ngt
  aw[!is.finite(aw)] <- 0
  aw
}

#' Wald test for ATETs.
#' @keywords internal
#' @noRd
atet_wald <- function(estimate, vcov_mat, df_resid, type_of_test) {
  L <- length(estimate)
  if (type_of_test == "zero") {
    R <- diag(L)
    rhs <- numeric(L)
    h0 <- "All effects are equal to zero"
  } else if (type_of_test == "equal") {
    if (L < 2L) {
      return(list(F = NA_real_, df1 = NA_integer_, df2 = df_resid,
                  p = NA_real_, h0 = "Effects are equal to each other",
                  title = "Test of equal ATETs"))
    }
    R <- cbind(matrix(-1, L - 1L, 1), diag(L - 1L))
    rhs <- numeric(L - 1L)
    h0 <- "Effects are equal to each other"
  }
  Rb  <- as.numeric(R %*% estimate - rhs)
  RVR <- R %*% vcov_mat %*% t(R)
  inv <- tryCatch(solve(RVR), error = function(e) NULL)
  pinv_note <- NULL
  if (is.null(inv)) {
    # RVR is rank-deficient (e.g. more ATET levels than clusters). Fall back to
    # SVD pseudoinverse; df1 is set to the numerical rank of RVR.
    s   <- svd(RVR)
    tol <- max(dim(RVR)) * max(s$d) * .Machine$double.eps
    r   <- sum(s$d > tol)
    if (r == 0L) {
      return(list(F = NA_real_, df1 = NA_integer_, df2 = df_resid,
                  p = NA_real_, h0 = h0, title = "Wald test", pinv_note = NULL))
    }
    inv <- s$v[, seq_len(r), drop = FALSE] %*%
           diag(1 / s$d[seq_len(r)], nrow = r) %*%
           t(s$u[, seq_len(r), drop = FALSE])
    df1 <- r
    pinv_note <- sprintf(
      "Note: VCov matrix is rank-deficient (rank %d of %d); pseudoinverse used for Wald test.",
      r, nrow(RVR)
    )
  } else {
    df1 <- nrow(R)
  }
  W     <- as.numeric(t(Rb) %*% inv %*% Rb)
  Fstat <- W / df1
  pval  <- stats::pf(Fstat, df1, df_resid, lower.tail = FALSE)
  list(F = Fstat, df1 = df1, df2 = df_resid, p = pval,
       h0 = h0, title = sprintf("Test of %s ATETs",
                                 if (type_of_test == "zero") "zero" else "equal"),
       pinv_note = pinv_note)
}
