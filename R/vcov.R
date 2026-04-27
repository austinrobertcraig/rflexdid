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

#' @keywords internal
#' @noRd
compute_vcov <- function(X, residuals, weights, pivot_keep, rank,
                         vcov_type = c("cluster", "robust"),
                         cluster = NULL) {
  vcov_type <- match.arg(vcov_type)
  n <- nrow(X)
  k <- rank
  X_kept <- X[, pivot_keep, drop = FALSE]
  u <- as.numeric(residuals)
  w <- as.numeric(weights)

  # Bread: (X' W X)^{-1}
  if (all(w == 1)) {
    bread_inv <- crossprod(X_kept)
    score_i <- X_kept * u  # n x k_kept; rows are X_i u_i
  } else {
    bread_inv <- crossprod(X_kept * sqrt(w))  # = sum_i w_i X_i X_i'
    score_i <- X_kept * (u * w)               # rows are w_i X_i u_i
  }
  bread <- chol2inv(chol(bread_inv))

  if (vcov_type == "robust") {
    meat <- crossprod(score_i)                # sum_i (w_i X_i u_i)(w_i X_i u_i)'
    adj <- n / (n - k)
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
    cluster_score <- rowsum(score_i, group = cluster_f, reorder = FALSE)
    meat <- crossprod(cluster_score)
    adj <- (G / (G - 1)) * ((n - 1) / (n - k))
  }

  V <- adj * (bread %*% meat %*% bread)
  # Symmetrize against tiny numerical asymmetry from the matrix products.
  V <- (V + t(V)) / 2
  V
}
