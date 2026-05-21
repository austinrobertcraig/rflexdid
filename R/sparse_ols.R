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
# and refactor on the kept submatrix. The drop pass is iterated up to
# `max_iter` times because a single pass through CHOLMOD's permutation may
# not reveal all linear dependencies -- the residual rank deficiency that
# remains after one drop becomes visible in the next factorization. If
# iteration does not converge we fall back to a pivoted-Cholesky pass on the
# dense p_kept x p_kept Gram matrix, which is rank-revealing without ever
# materializing the n x p design as a dense R matrix. The previous lm.fit
# fallback (which called as.matrix(X)) is removed because at n=15M and
# p~2500 that allocation is ~290 GB.
#
# The dropped columns are filled with NA in the returned beta to match
# `lm.fit` semantics. Cell-indicator columns are well-identified by
# construction (one column per non-empty (kind, g, t)) so they are never the
# columns dropped; the redundancy lives in the interaction blocks. That
# keeps Stata-comparable cell coefficients intact.
#
# Returns:
#   beta       numeric, length p, with NA in dropped positions
#   pivot_keep logical, length p
#   bread      p_kept x p_kept dense matrix, equal to (X' diag(w) X)^{-1}
#              restricted to the kept columns. Reused by compute_vcov() and
#              atet() to avoid refactoring.

#' @keywords internal
#' @noRd
sparse_ols <- function(X, y, w, tol = 1e-7, max_iter = 6L) {
  if (!methods::is(X, "CsparseMatrix")) {
    X <- methods::as(X, "CsparseMatrix")
  }
  p <- ncol(X)
  cn <- colnames(X)
  unweighted <- isTRUE(all(w == 1))

  build_normal <- function(X_sub) {
    if (unweighted) {
      list(XtX = Matrix::crossprod(X_sub),
           Xty = as.numeric(Matrix::crossprod(X_sub, y)))
    } else {
      sw <- sqrt(w)
      Xw <- sw * X_sub
      list(XtX = Matrix::crossprod(Xw),
           Xty = as.numeric(Matrix::crossprod(Xw, sw * y)))
    }
  }

  finish <- function(pivot_keep, ne) {
    # ne is the normal-equations object on the kept submatrix and is known
    # to be (numerically) PD here. Factor without the ridge so beta and
    # bread reflect the unperturbed system.
    p_kept <- sum(pivot_keep)
    ch <- Matrix::Cholesky(ne$XtX, perm = TRUE, LDL = TRUE)
    beta_kept <- as.numeric(Matrix::solve(ch, ne$Xty))
    bread <- as.matrix(Matrix::solve(ch, Matrix::Diagonal(p_kept)))
    bread <- (bread + t(bread)) / 2
    beta <- rep(NA_real_, p)
    beta[pivot_keep] <- beta_kept
    names(beta) <- cn
    list(beta = beta, pivot_keep = pivot_keep, bread = bread)
  }

  pivot_keep <- rep(TRUE, p)
  for (iter in seq_len(max_iter)) {
    p_kept <- sum(pivot_keep)
    X_kept <- if (p_kept == p) X else X[, pivot_keep, drop = FALSE]
    ne <- build_normal(X_kept)

    # CHOLMOD's Cholesky bails out on an *exact* zero leading principal minor
    # even with LDL=TRUE -- the README example with three covariates trips it.
    # A tiny relative ridge perturbs those exact zeros into detectably-small
    # pivots that we then flag and drop in the rank-revealing step. The ridge
    # is small enough (~1e-12 * max diag) not to move kept-coefficient values.
    diag_max <- max(Matrix::diag(ne$XtX))
    ridge <- if (diag_max > 0) 1e-12 * diag_max else 1e-12
    XtX_reg <- ne$XtX + ridge * Matrix::Diagonal(p_kept)
    ch <- Matrix::Cholesky(XtX_reg, perm = TRUE, LDL = TRUE)
    Dvec <- as.numeric(Matrix::diag(Matrix::expand1(ch, "D")))
    # tol scales by the largest pivot. Anything below counts as a
    # rank-deficient column relative to the overall design magnitude.
    Dthr <- tol * max(abs(Dvec))
    zero_pos <- which(Dvec <= Dthr)

    if (length(zero_pos) == 0L) {
      return(finish(pivot_keep, ne))
    }

    # Map zero pivot positions through CHOLMOD's permutation to columns of
    # the current kept submatrix, then translate to original X indices.
    perm <- ch@perm + 1L
    drop_in_kept <- perm[zero_pos]
    drop_in_orig <- which(pivot_keep)[drop_in_kept]
    pivot_keep[drop_in_orig] <- FALSE
  }

  # Iteration did not converge. Fall back to a pivoted-Cholesky pass on the
  # dense p_kept x p_kept Gram matrix. This is rank-revealing and lets us
  # identify and drop the remaining dependent columns without ever forming
  # the n x p_kept dense design.
  p_kept <- sum(pivot_keep)
  X_kept <- X[, pivot_keep, drop = FALSE]
  ne <- build_normal(X_kept)
  Gram <- as.matrix(ne$XtX)
  Gram <- (Gram + t(Gram)) / 2
  R_piv <- suppressWarnings(chol(Gram, pivot = TRUE))
  r <- attr(R_piv, "rank")
  piv <- attr(R_piv, "pivot")
  if (r < p_kept) {
    drop_in_kept <- piv[seq.int(r + 1L, p_kept)]
    drop_in_orig <- which(pivot_keep)[drop_in_kept]
    pivot_keep[drop_in_orig] <- FALSE
  }

  # Refactor sparsely on the (now full-rank) kept set and finish.
  p_kept2 <- sum(pivot_keep)
  X_kept2 <- X[, pivot_keep, drop = FALSE]
  ne2 <- build_normal(X_kept2)
  finish(pivot_keep, ne2)
}
