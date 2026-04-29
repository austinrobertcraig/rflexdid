# Sparse OLS solver used by flexdid().
#
# Solves min || sqrt(w) (y - X beta) ||^2 by forming the weighted normal
# equations and factorizing X' diag(w) X with a sparse Cholesky.
#
# Rank deficiency. flexdid's design matrix is regularly rank-deficient -- it
# emits one cell coefficient per treated (group, time) so the group X xvars,
# time X xvars, and (cell X xvars) blocks overlap structurally. We detect
# rank-deficient columns from the LDL factor's D diagonal (a near-zero pivot
# means that column is linearly dependent on prior ones), drop those columns,
# and refactor on the kept submatrix. The dropped columns are filled with NA
# in the returned beta to match `lm.fit` semantics. Cell-indicator columns
# are well-identified by construction (one column per non-empty (kind, g, t))
# so they are never the columns dropped; the redundancy lives in the
# interaction blocks. That keeps Stata-comparable cell coefficients intact.
#
# Returns:
#   beta       numeric, length p, with NA in dropped positions
#   pivot_keep logical, length p
#   bread      p_kept x p_kept dense matrix, equal to (X' diag(w) X)^{-1}
#              restricted to the kept columns. Reused by compute_vcov() and
#              atet() to avoid refactoring.

#' @keywords internal
#' @noRd
sparse_ols <- function(X, y, w, tol = 1e-7) {
  if (!methods::is(X, "CsparseMatrix")) {
    X <- methods::as(X, "CsparseMatrix")
  }
  p <- ncol(X)
  cn <- colnames(X)

  unweighted <- isTRUE(all(w == 1))
  build_normal <- function(X) {
    if (unweighted) {
      list(XtX = Matrix::crossprod(X),
           Xty = as.numeric(Matrix::crossprod(X, y)))
    } else {
      sw <- sqrt(w)
      Xw <- sw * X
      list(XtX = Matrix::crossprod(Xw),
           Xty = as.numeric(Matrix::crossprod(Xw, sw * y)))
    }
  }

  ne <- build_normal(X)
  # CHOLMOD's Cholesky bails out on an *exact* zero leading principal minor
  # even with LDL=TRUE -- the README example with three covariates trips it.
  # A tiny relative ridge perturbs those exact zeros into detectably-small
  # pivots that we then flag and drop in the rank-revealing step. The ridge is
  # small enough (~1e-12 * max diag) not to move kept-coefficient values.
  diag_max <- max(Matrix::diag(ne$XtX))
  ridge <- if (diag_max > 0) 1e-12 * diag_max else 1e-12
  XtX_reg <- ne$XtX + ridge * Matrix::Diagonal(p)
  ch <- Matrix::Cholesky(XtX_reg, perm = TRUE, LDL = TRUE)
  Dvec <- as.numeric(Matrix::diag(Matrix::expand1(ch, "D")))
  # tol scales by the largest pivot. Anything below counts as a rank-deficient
  # column relative to the overall design magnitude.
  Dthr <- tol * max(abs(Dvec))
  zero_pos <- which(Dvec <= Dthr)

  if (length(zero_pos) == 0L) {
    # Solve and invert through the sparse factorization. Refactor without the
    # ridge first so beta and bread reflect the unperturbed system.
    ch_clean <- Matrix::Cholesky(ne$XtX, perm = TRUE, LDL = TRUE)
    beta <- as.numeric(Matrix::solve(ch_clean, ne$Xty))
    bread <- as.matrix(Matrix::solve(ch_clean, Matrix::Diagonal(p)))
    bread <- (bread + t(bread)) / 2
    names(beta) <- cn
    return(list(beta = beta, pivot_keep = rep(TRUE, p), bread = bread))
  }

  # Map zero pivots back to original column indices via the factor's
  # permutation, drop those columns, and refactor on the kept set.
  perm <- ch@perm + 1L           # CHOLMOD permutation, 1-based
  drop_cols <- perm[zero_pos]
  pivot_keep <- rep(TRUE, p)
  pivot_keep[drop_cols] <- FALSE

  X_kept <- X[, pivot_keep, drop = FALSE]
  ne2 <- build_normal(X_kept)
  # Same ridge guard for the kept submatrix; it should be cleanly PD now but
  # in pathological cases the first-pass column drop may not be enough.
  diag2_max <- max(Matrix::diag(ne2$XtX))
  ridge2 <- if (diag2_max > 0) 1e-12 * diag2_max else 1e-12
  XtX2_reg <- ne2$XtX + ridge2 * Matrix::Diagonal(sum(pivot_keep))
  ch2_chk <- Matrix::Cholesky(XtX2_reg, perm = TRUE, LDL = TRUE)
  D2 <- as.numeric(Matrix::diag(Matrix::expand1(ch2_chk, "D")))
  if (any(D2 <= tol * max(abs(D2)))) {
    return(dense_fallback(X, y, w, unweighted))
  }
  ch2 <- Matrix::Cholesky(ne2$XtX, perm = TRUE, LDL = TRUE)
  beta_kept <- as.numeric(Matrix::solve(ch2, ne2$Xty))
  bread <- as.matrix(Matrix::solve(ch2, Matrix::Diagonal(sum(pivot_keep))))
  bread <- (bread + t(bread)) / 2

  beta <- rep(NA_real_, p)
  beta[pivot_keep] <- beta_kept
  names(beta) <- cn
  list(beta = beta, pivot_keep = pivot_keep, bread = bread)
}

#' @keywords internal
#' @noRd
dense_fallback <- function(X, y, w, unweighted) {
  X_dense <- as.matrix(X)
  fit <- if (unweighted) stats::lm.fit(X_dense, y)
         else            stats::lm.wfit(X_dense, y, w)
  beta <- fit$coefficients
  names(beta) <- colnames(X)
  pivot_keep <- !is.na(beta)
  X_kept <- X_dense[, pivot_keep, drop = FALSE]
  bread <- if (unweighted) chol2inv(chol(crossprod(X_kept)))
           else            chol2inv(chol(crossprod(X_kept * sqrt(w))))
  bread <- (bread + t(bread)) / 2
  list(beta = beta, pivot_keep = pivot_keep, bread = bread)
}
