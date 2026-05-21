# Variance-covariance estimators for the FLEX regression.
#
# Stata's `regress, vce(robust)` is HC1 with multiplier n/(n-k).
# Stata's `regress, vce(cluster X)` is CR1 with multiplier
#   (G/(G-1)) * ((n-1)/(n-k)).
#
# For weighted OLS we treat the weights as pweights (Stata's typical default
# for clustered DiD designs). Under pweights, beta solves
#   sum_i w_i X_i' (y_i - X_i' beta) = 0,
# so the estimating function is psi_i = w_i X_i' u_i, the bread is
#   (sum_i w_i X_i X_i')^{-1} = (X' W X)^{-1},
# and the sandwich is
#   V_robust  = (X' W X)^{-1} (sum_i w_i^2 u_i^2 X_i X_i') (X' W X)^{-1}
#   V_cluster = (X' W X)^{-1} (sum_g psi_g psi_g') (X' W X)^{-1},
#                with psi_g = sum_{i in g} w_i X_i' u_i.
#
# We tested the unweighted output against sandwich::vcovHC(..., type="HC1")
# and sandwich::vcovCL(..., type="HC1", cadjust=TRUE); the formulas below
# reproduce those values. We hand-roll because sandwich's lm methods treat
# weights as aweights (different score function), which would not match
# Stata for pweighted models.
#
# Implementation note: we never materialize the n x k_kept per-observation
# score as a dense R matrix. On large n that allocation is what blows up
# memory (an n=15M, k_kept~2500 design needs ~290 GB for the dense score).
#
# For cluster vcov we aggregate to G x k_kept via a sparse G x n indicator
# matrix, then take V = crossprod(S_g %*% bread) -- PSD by construction.
#
# For robust vcov we cannot aggregate, so we compute meat = X' diag((u w)^2) X
# as a sparse k_kept x k_kept crossprod and then form V = bread meat bread.
# To guarantee numerical PSD (so summary()'s sqrt(diag(V)) never produces NaN),
# we factor meat = R' R via Cholesky and compute V = crossprod(R %*% bread).

#' @keywords internal
#' @noRd
compute_vcov <- function(X, residuals, weights, pivot_keep, rank,
                         vcov_type = c("cluster", "robust"),
                         cluster = NULL, bread = NULL) {
  vcov_type <- match.arg(vcov_type)
  n <- nrow(X)
  k <- rank
  X_kept <- X[, pivot_keep, drop = FALSE]
  u <- as.numeric(residuals)
  w <- as.numeric(weights)
  unweighted <- all(w == 1)

  # Row-scale X by the score weight; stays sparse if X_kept is sparse.
  score_sparse <- if (unweighted) X_kept * u else X_kept * (u * w)

  # Bread: (X' W X)^{-1}. flexdid() always supplies this; the recompute
  # branch is only for callers that build a vcov by hand.
  if (is.null(bread)) {
    bread_inv <- if (unweighted) Matrix::crossprod(X_kept) else
                                 Matrix::crossprod(X_kept * sqrt(w))
    bread <- chol2inv(chol(as.matrix(bread_inv)))
  }

  if (vcov_type == "robust") {
    # meat = sum_i (w_i u_i)^2 X_i X_i' as a sparse k_kept x k_kept crossprod.
    # Factor meat = R' R so V = bread meat bread = crossprod(R %*% bread),
    # which is PSD by construction at floating-point precision.
    meat <- as.matrix(Matrix::crossprod(score_sparse))
    meat <- (meat + t(meat)) / 2
    # Tiny ridge guards Cholesky against rank deficiency in meat (e.g. when
    # u has many exact zeros). Scaled to meat's diagonal so it never moves
    # the PD case at any meaningful precision.
    diag_meat <- diag(meat)
    ridge <- 1e-14 * max(c(abs(diag_meat), 1))
    meat_chol <- chol(meat + ridge * diag(nrow(meat)))
    M <- meat_chol %*% bread
    adj <- n / (n - k)
    V <- adj * crossprod(M)
  } else {
    if (is.null(cluster)) {
      stop("Internal: cluster vector required for vcov_type = 'cluster'.",
           call. = FALSE)
    }
    cluster_f <- as.factor(cluster)
    G <- nlevels(cluster_f)
    if (G < 2L) {
      stop("Cluster variable has fewer than 2 unique values.", call. = FALSE)
    }
    # Sparse cluster aggregator: row g of cluster_agg selects obs in cluster g.
    # cluster_agg %*% score_sparse forms per-cluster score sums directly,
    # avoiding the n x k_kept dense intermediate the rowsum() path required.
    cluster_agg <- Matrix::sparseMatrix(
      i = as.integer(cluster_f),
      j = seq_len(n),
      x = rep(1, n),
      dims = c(G, n)
    )
    S_mat <- as.matrix(cluster_agg %*% score_sparse)  # G x k_kept dense
    M <- S_mat %*% bread
    adj <- (G / (G - 1)) * ((n - 1) / (n - k))
    V <- adj * crossprod(M)
  }
  # Symmetrize against tiny numerical asymmetry from the matrix products.
  V <- (V + t(V)) / 2
  V
}
