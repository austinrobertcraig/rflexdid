# Design matrix construction for the FLEX regression.
# Mirrors the Stata regression in stata/flexdid.ado lines 156-222:
#
#   regress y TxGlags [TxGleads] TxGlagsXX [TxGleadsXX]
#                xvars i.group i.time i.group#c.xvars i.time#c.xvars
#                xnotinteracted [if touse] [pw=w], vce(...)
#
# Where, for each treated (g, t) cell with t >= cohort_g (lag) or
# t <= cohort_g - 2 (lead):
#
#   TxGlags column = 1[group==g & time==t & .Tx==1]
#   TxGlagsXX columns = the indicator times each xvar
#
# The base period t = cohort_g - 1 is never a column (.Tx == 0 there).
# Note: c.Cohort#g.group#t.time collapses to g.group when group is in cohort
# c (each group belongs to exactly one cohort), so we identify cells by
# (g, t) alone. See the plan for the full justification.

# All matrices are returned as Matrix::dgCMatrix so we can build a single
# sparse column-stacked design and pass either the sparse or a coerced
# dense form to lm.fit.

#' Build the FLEX design matrix.
#'
#' @return A list with:
#'   - `X`: sparse model matrix (n x p)
#'   - `colnames`: character vector of length p
#'   - `cells`: data.frame(g, t, kind, col_indicator, col_x_offset) with one
#'     row per (g, t, kind) cell. `col_indicator` is the column index of the
#'     cell's bare indicator; `col_x_offset` is the column where its first
#'     covariate interaction lives (followed by k-1 more contiguous columns,
#'     in the order of `xvars_matrix` columns).
#'   - `xvars_matrix`: the n x k matrix of interacted covariates as built
#'     by model.matrix (kept for ATET gradient construction).
#'   - `xvars_names`: character vector of length k.
#'   - `n_xvars`: integer k.
#'   - `block_index`: list of integer column ranges for each block, useful
#'     for diagnostics and for the F-test that mirrors `testparm`.
#' @keywords internal
#' @noRd
build_design <- function(data, formula, group, time, cohort, tx_indicator,
                          xnotinteracted, specification) {
  n <- nrow(data)
  if (length(group) != n || length(time) != n ||
      length(cohort) != n || length(tx_indicator) != n) {
    stop("Internal: design vectors not aligned to nrow(data).", call. = FALSE)
  }

  # 1. Interacted covariates xvars (drop intercept; allow factors via formula)
  xvar_form <- formula[c(1, 3)]  # ~ rhs
  if (length(all.vars(xvar_form)) == 0L) {
    xvars_mm <- matrix(0, n, 0)
    xvars_names <- character(0)
  } else {
    xvars_mm <- model.matrix(xvar_form, data = data)
    intercept_col <- which(colnames(xvars_mm) == "(Intercept)")
    if (length(intercept_col)) {
      xvars_mm <- xvars_mm[, -intercept_col, drop = FALSE]
    }
    xvars_names <- colnames(xvars_mm)
  }
  k <- ncol(xvars_mm)

  # 2. xnotinteracted (additive only, can include factors)
  if (is.null(xnotinteracted)) {
    xni_mm <- matrix(0, n, 0)
    xni_names <- character(0)
  } else {
    if (is.character(xnotinteracted)) {
      xni_form <- stats::as.formula(paste("~", paste(xnotinteracted, collapse = "+")))
    } else if (inherits(xnotinteracted, "formula")) {
      xni_form <- xnotinteracted
      if (length(xni_form) == 3L) xni_form <- xni_form[c(1, 3)]
    } else {
      stop("`xnotinteracted` must be NULL, a character vector of column names, or a formula.",
           call. = FALSE)
    }
    xni_mm <- model.matrix(xni_form, data = data)
    intercept_col <- which(colnames(xni_mm) == "(Intercept)")
    if (length(intercept_col)) {
      xni_mm <- xni_mm[, -intercept_col, drop = FALSE]
    }
    xni_names <- colnames(xni_mm)
  }

  # Disjointness check (the Stata "xvars and xnotinteracted must be distinct" rule)
  inter <- intersect(xvars_names, xni_names)
  if (length(inter) > 0L) {
    stop(sprintf(
      "Variables present in both formula and xnotinteracted: %s. They must be distinct.",
      paste(inter, collapse = ", ")
    ), call. = FALSE)
  }

  # 3. Determine cells -- unique (g, t, kind) combos for which we emit columns.
  is_lag  <- tx_indicator == 1L & cohort > 0 & time >= cohort
  is_lead <- tx_indicator == 1L & cohort > 0 & time <= (cohort - 2)
  if (specification == "lagsonly") {
    use <- is_lag
  } else {
    use <- is_lag | is_lead
  }
  cell_rows <- which(use)
  if (length(cell_rows) == 0L) {
    stop("No treated post-treatment observations available; design has no treatment cells.",
         call. = FALSE)
  }
  cell_g <- group[cell_rows]
  cell_t <- time[cell_rows]
  cell_kind <- ifelse(is_lag[cell_rows], "lag", "lead")

  cell_keys <- paste(cell_kind, cell_g, cell_t, sep = "|")
  uniq_first <- !duplicated(cell_keys)
  uniq_g <- cell_g[uniq_first]
  uniq_t <- cell_t[uniq_first]
  uniq_kind <- cell_kind[uniq_first]
  uniq_keys <- cell_keys[uniq_first]

  # Order: lag cells first (so column order matches Stata's TxGlags ... TxGleads ...)
  ord <- order(uniq_kind != "lag", uniq_g, uniq_t)
  uniq_g <- uniq_g[ord]
  uniq_t <- uniq_t[ord]
  uniq_kind <- uniq_kind[ord]
  uniq_keys <- uniq_keys[ord]
  num_cells <- length(uniq_keys)
  cell_id_for_row <- match(cell_keys, uniq_keys)

  cell_indicator_names <- sprintf("Tx_%s:g%s:t%s", uniq_kind, uniq_g, uniq_t)

  # 4. Build the indicator block: i = cell_rows, j = cell_id_for_row, x = 1
  block_indicators <- Matrix::sparseMatrix(
    i = cell_rows,
    j = cell_id_for_row,
    x = rep(1, length(cell_rows)),
    dims = c(n, num_cells)
  )
  colnames(block_indicators) <- cell_indicator_names

  # 5. Build the indicator x covariate block.
  #   For each of the |cell_rows| treated obs and each of k xvars, emit one
  #   triplet at (cell_rows[r], (cell_id_for_row[r] - 1)*k + xvar_index).
  if (k > 0L) {
    nrep <- length(cell_rows)
    i_xx <- rep(cell_rows, each = k)
    j_xx <- rep((cell_id_for_row - 1L) * k, each = k) +
            rep(seq_len(k), times = nrep)
    # x_xx must be the row of xvars_mm for each cell_row, repeated k entries
    # in the same order as the j offsets
    x_xx <- as.vector(t(xvars_mm[cell_rows, , drop = FALSE]))
    block_interactions <- Matrix::sparseMatrix(
      i = i_xx, j = j_xx, x = x_xx,
      dims = c(n, num_cells * k)
    )
    int_colnames <- character(num_cells * k)
    for (j in seq_len(num_cells)) {
      base <- (j - 1L) * k
      int_colnames[(base + 1L):(base + k)] <-
        paste0(cell_indicator_names[j], ":", xvars_names)
    }
    colnames(block_interactions) <- int_colnames
  } else {
    block_interactions <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(n, 0)
    )
    int_colnames <- character(0)
  }

  # 6. xvars additive block (sparse-from-dense)
  block_xvars <- if (k > 0L) {
    m <- methods::as(xvars_mm, "CsparseMatrix")
    colnames(m) <- xvars_names
    m
  } else {
    Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                         dims = c(n, 0))
  }

  # 7. factor(group) dummies (Stata-style: omit first level)
  group_f <- factor(group)
  group_levels <- levels(group_f)
  group_dummies <- diag(length(group_levels))[as.integer(group_f), , drop = FALSE]
  # drop first level
  group_dummies <- group_dummies[, -1, drop = FALSE]
  group_dum_names <- paste0("group:", group_levels[-1])
  block_group <- methods::as(group_dummies, "CsparseMatrix")
  colnames(block_group) <- group_dum_names

  # 8. factor(time) dummies
  time_f <- factor(time)
  time_levels <- levels(time_f)
  time_dummies <- diag(length(time_levels))[as.integer(time_f), , drop = FALSE]
  time_dummies <- time_dummies[, -1, drop = FALSE]
  time_dum_names <- paste0("time:", time_levels[-1])
  block_time <- methods::as(time_dummies, "CsparseMatrix")
  colnames(block_time) <- time_dum_names

  # 9. group x xvars
  if (k > 0L && length(group_levels) > 1L) {
    blk <- group_dummies
    cols <- list()
    cnames <- character(0)
    for (xj in seq_len(k)) {
      m <- blk * xvars_mm[, xj]
      cols[[xj]] <- m
      cnames <- c(cnames, paste0(group_dum_names, ":", xvars_names[xj]))
    }
    block_group_x <- do.call(cbind, cols)
    block_group_x <- methods::as(block_group_x, "CsparseMatrix")
    colnames(block_group_x) <- cnames
  } else {
    block_group_x <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                           x = numeric(0), dims = c(n, 0))
  }

  # 10. time x xvars
  if (k > 0L && length(time_levels) > 1L) {
    blk <- time_dummies
    cols <- list()
    cnames <- character(0)
    for (xj in seq_len(k)) {
      m <- blk * xvars_mm[, xj]
      cols[[xj]] <- m
      cnames <- c(cnames, paste0(time_dum_names, ":", xvars_names[xj]))
    }
    block_time_x <- do.call(cbind, cols)
    block_time_x <- methods::as(block_time_x, "CsparseMatrix")
    colnames(block_time_x) <- cnames
  } else {
    block_time_x <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                          x = numeric(0), dims = c(n, 0))
  }

  # 11. xnotinteracted
  block_xni <- if (length(xni_names) > 0L) {
    m <- methods::as(xni_mm, "CsparseMatrix")
    colnames(m) <- xni_names
    m
  } else {
    Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                         dims = c(n, 0))
  }

  # 12. Intercept (single column of 1s; we add explicitly to mirror lm()
  # default and so that group/time FEs are identified relative to base levels)
  block_int <- Matrix::sparseMatrix(
    i = seq_len(n), j = rep(1L, n), x = rep(1, n),
    dims = c(n, 1)
  )
  colnames(block_int) <- "(Intercept)"

  # cbind all blocks. Order: (Intercept), TxGlags*, TxGleads*, TxGlagsXX*,
  # TxGleadsXX*, xvars, group dummies, time dummies, group:xvars,
  # time:xvars, xnotinteracted. (Stata order is similar; we keep TxGlags* and
  # TxGleads* together for clarity but split the indicators from the
  # interactions for naming traceability. Order doesn't affect the OLS fit.)
  X <- cbind(
    block_int,
    block_indicators,
    block_interactions,
    block_xvars,
    block_group,
    block_time,
    block_group_x,
    block_time_x,
    block_xni
  )

  # 13. Cell -> column index lookup. Indicator columns sit right after the
  # intercept; interaction columns follow, in (cell, xvar) row-major order.
  base_indicator <- 1L  # column 1 is intercept; indicators start at 2
  base_interaction <- 1L + num_cells
  cells <- data.frame(
    g = uniq_g,
    t = uniq_t,
    kind = uniq_kind,
    col_indicator = base_indicator + seq_len(num_cells),
    col_x_offset = if (k > 0L) base_interaction + (seq_len(num_cells) - 1L) * k + 1L
                   else rep(NA_integer_, num_cells),
    stringsAsFactors = FALSE
  )

  # Block index ranges for diagnostic / F-test use
  p <- ncol(X)
  block_index <- list(
    intercept = 1L,
    indicators = base_indicator + seq_len(num_cells),
    interactions = if (k > 0L) base_interaction + seq_len(num_cells * k) else integer(0),
    xvars = if (k > 0L) {
      start <- 1L + num_cells + num_cells * k + 1L
      seq.int(start, length.out = k)
    } else integer(0),
    group = if (length(group_dum_names) > 0L) {
      start <- 1L + num_cells + num_cells * k + k + 1L
      seq.int(start, length.out = length(group_dum_names))
    } else integer(0),
    time = if (length(time_dum_names) > 0L) {
      start <- 1L + num_cells + num_cells * k + k + length(group_dum_names) + 1L
      seq.int(start, length.out = length(time_dum_names))
    } else integer(0),
    group_x = if (ncol(block_group_x) > 0L) {
      start <- 1L + num_cells + num_cells * k + k + length(group_dum_names) +
               length(time_dum_names) + 1L
      seq.int(start, length.out = ncol(block_group_x))
    } else integer(0),
    time_x = if (ncol(block_time_x) > 0L) {
      start <- 1L + num_cells + num_cells * k + k + length(group_dum_names) +
               length(time_dum_names) + ncol(block_group_x) + 1L
      seq.int(start, length.out = ncol(block_time_x))
    } else integer(0),
    xni = if (length(xni_names) > 0L) {
      start <- 1L + num_cells + num_cells * k + k + length(group_dum_names) +
               length(time_dum_names) + ncol(block_group_x) + ncol(block_time_x) + 1L
      seq.int(start, length.out = length(xni_names))
    } else integer(0)
  )

  list(
    X = X,
    colnames = colnames(X),
    cells = cells,
    xvars_matrix = xvars_mm,
    xvars_names = xvars_names,
    n_xvars = k,
    xnotinteracted_names = xni_names,
    block_index = block_index,
    group_levels = group_levels,
    time_levels = time_levels
  )
}
