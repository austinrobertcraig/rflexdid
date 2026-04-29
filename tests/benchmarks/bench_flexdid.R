# Benchmark / profile flexdid() on a stress-shape synthetic panel.
# Run from the package root with:
#   Rscript tests/benchmarks/bench_flexdid.R
# Excluded from R CMD build via .Rbuildignore.
#
# Baseline (Phase 0, before any optimization):
#   Panel A: ng=100, ny=12, k=6, lagsonly. n=1200, p=3780, sparsity=0.60%.
#     flexdid()  : 23.0s
#       - stats::lm.fit         : 11.03s (57%) -- dense QR on n x p with p > n
#       - %*% (X_dense %*% beta): 6.71s  (35%) -- fitted-value matmul, dense
#       - compute_vcov          : 1.92s  (10%)
#     atet(overall|byexposure|bycohort): ~1.2s each
#       - dominated by repeated chol(crossprod(X_kept)) inside atet()
#   Panel B: ng=200, ny=12, k=5, lagsonly. n=2400, p=6414. (rerun on this shape
#     under the same baseline code path: 215s, with lm.fit 73% of time)
#
# After all changes (sparse rank-revealing fit + cached bread + vectorized
# design blocks + F-stat row selection + sandwich rewrite):
#   Panel A (this script): 23.0s -> 0.78s   (~29x faster)
#   Panel B (ng=200, k=5): 215.3s -> 4.1s   (~52x faster)
#   README example k=3:    7.7s -> 0.17s    (~45x faster)
# atet() variants now run in 5-80ms each.
#
# What changed:
#   - flexdid.R       : keep design sparse, cache bread, replace R*V*R' with
#                       V[non_int, non_int] selection.
#   - sparse_ols.R    : new; sparse Cholesky LDL with rank-revealing column
#                       drops, dense chol2inv on the kept submatrix for the
#                       bread, dense lm.fit only as a defensive fallback.
#   - vcov.R         : reuse cached bread; sandwich as crossprod(S %*% bread).
#   - atet.R         : reuse cached bread; X stays sparse, score stays sparse
#                       until the final IF product.
#   - design.R       : group/time dummies and group/time x covariate
#                       interactions built directly as sparse triplets, no
#                       per-covariate loop and no dense intermediate cbind.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(Matrix)
})

# ---- Synthetic stress panel --------------------------------------------------
# Larger and wider than the existing test fixture: 200 groups x 20 years,
# staggered cohorts, k continuous covariates that get interacted with treatment,
# group, and time -- which is the regime the user reports is slow.
make_stress_panel <- function(seed = 42, ng = 100, ny = 12, year0 = 2030,
                              cohorts = c(2034, 2037),
                              k = 6) {
  set.seed(seed)
  df <- expand.grid(unit = seq_len(ng), year = year0:(year0 + ny - 1L))
  per_cohort <- ng %/% (length(cohorts) + 1L)
  chrt <- rep(0L, ng)
  for (i in seq_along(cohorts)) {
    rng <- ((i - 1L) * per_cohort + 1L):(i * per_cohort)
    chrt[rng] <- cohorts[i]
  }
  df$cohort_t <- chrt[df$unit]
  df$tx <- as.integer(df$cohort_t > 0 & df$year >= df$cohort_t)
  Xmat <- matrix(rnorm(nrow(df) * k), nrow = nrow(df))
  colnames(Xmat) <- paste0("x", seq_len(k))
  df <- cbind(df, as.data.frame(Xmat))
  beta_x <- rnorm(k, sd = 0.2)
  df$y <- 20 + 0.5 * df$tx + 0.1 * df$year - 0.02 * df$unit +
          as.numeric(Xmat %*% beta_x) +
          rnorm(nrow(df), sd = 0.5)
  df
}

cat("Building stress panel ...\n")
df <- make_stress_panel()
xnames <- grep("^x[0-9]+$", names(df), value = TRUE)
fml <- as.formula(paste("y ~", paste(xnames, collapse = " + ")))
cat("  n rows :", nrow(df), "\n")
cat("  groups :", length(unique(df$unit)), "\n")
cat("  years  :", length(unique(df$year)), "\n")
cat("  k xvars:", length(xnames), "\n\n")

# ---- Time the full call ------------------------------------------------------
cat("Timing flexdid() ...\n")
t_full <- system.time({
  fit <- flexdid(fml, data = df, tx = "tx", group = "unit", time = "year",
                 specification = "lagsonly", vcov = "cluster")
})
cat("  elapsed:", sprintf("%.3fs", t_full["elapsed"]), "\n")
cat("  ncol(X):", ncol(fit$X), "  nrow(X):", nrow(fit$X), "\n")
nnz <- Matrix::nnzero(fit$design$X)
cat("  sparse nnz:", nnz, sprintf(" (%.2f%% dense)\n",
                                  100 * nnz / prod(dim(fit$design$X))))

cat("\nTiming atet() flavors ...\n")
for (ty in c("overall", "byexposure", "bycohort")) {
  t_at <- system.time(atet(fit, type = ty))
  cat(sprintf("  atet(%-12s) : %.3fs\n", ty, t_at["elapsed"]))
}

# ---- Rprof breakdown of flexdid() -------------------------------------------
cat("\nRprof breakdown of flexdid() (1 run) ...\n")
prof_file <- tempfile(fileext = ".out")
Rprof(prof_file, line.profiling = FALSE, memory.profiling = FALSE,
      interval = 0.005)
fit2 <- flexdid(fml, data = df, tx = "tx", group = "unit", time = "year",
                specification = "lagsonly", vcov = "cluster")
Rprof(NULL)
sm <- summaryRprof(prof_file)
cat("\nTop self-time functions:\n")
print(head(sm$by.self[, c("self.time", "self.pct", "total.time", "total.pct")], 15))
cat("\nTop total-time functions:\n")
print(head(sm$by.total[, c("total.time", "total.pct", "self.time", "self.pct")], 20))
unlink(prof_file)

cat("\nDone.\n")
