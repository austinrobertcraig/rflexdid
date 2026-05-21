# Sanity check that compute_vcov() no longer allocates an n x k_kept dense
# score matrix on a large panel. Mirrors the user's call shape:
#   formula = y ~ 1, specification = "lagsandleads", vcov = "cluster".
# Run with:
#   Rscript tests/benchmarks/bench_large_cluster.R

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})

set.seed(1)
# Replicate panel via repeated observations per (school, year) cell so n is
# much larger than the design column count. This mirrors the user's setup
# (15M rows but only ~2500 kept design columns).
ng <- 1000L
ny <- 15L
reps <- 80L   # observations per (school, year) cell
year0 <- 2010L
cohorts <- c(2014L, 2017L, 2020L)
base <- expand.grid(school = seq_len(ng), year = year0:(year0 + ny - 1L))
per_cohort <- ng %/% (length(cohorts) + 1L)
chrt <- rep(0L, ng)
for (i in seq_along(cohorts)) {
  rng <- ((i - 1L) * per_cohort + 1L):(i * per_cohort)
  chrt[rng] <- cohorts[i]
}
base$cohort_t <- chrt[base$school]
base$tx <- as.integer(base$cohort_t > 0L & base$year >= base$cohort_t)
df <- base[rep(seq_len(nrow(base)), each = reps), , drop = FALSE]
df$y <- 0.3 * df$tx + 0.05 * df$year + rnorm(nrow(df), sd = 0.5)

cat("Panel:\n")
cat("  n rows :", nrow(df), "\n")
cat("  groups :", ng, "\n")
cat("  years  :", ny, "\n")
cat("  cohorts:", length(cohorts), "treated +", "1 control\n\n")

gc(reset = TRUE)
t1 <- system.time({
  fit <- flexdid(y ~ 1, data = df, tx = "tx", group = "school", time = "year",
                 specification = "lagsandleads", vcov = "cluster")
})
g1 <- gc()

cat(sprintf("flexdid() elapsed : %.2fs\n", t1["elapsed"]))
cat(sprintf("  ncol(X)  : %d   nrow(X) : %d\n", ncol(fit$X), nrow(fit$X)))
cat(sprintf("  rank     : %d   df_res  : %d\n", fit$rank, fit$df_residual))
cat(sprintf("  peak Mb  : %.1f Mb (max used since gc(reset=TRUE))\n",
            sum(g1[, "max used"] * c(8, 1024) / (1024 * 1024))))

# Old code path would have allocated nrow * rank * 8 bytes for score_i.
old_bytes <- as.numeric(nrow(fit$X)) * as.numeric(fit$rank) * 8
cat(sprintf("  old dense-score allocation would have been: %.2f Gb\n",
            old_bytes / 1024^3))

# Cross-check: the V from compute_vcov should have PSD-ish diagonal (no NaNs).
se <- sqrt(diag(fit$vcov)[fit$pivot_keep])
cat(sprintf("  any NaN SE : %s   max SE: %.3f\n", any(is.na(se)), max(se)))

cat("\n-- atet() byexposure (cluster) --\n")
gc(reset = TRUE)
t2 <- system.time({
  a_be <- atet(fit, type = "byexposure")
})
g2 <- gc()
cat(sprintf("atet(byexposure) elapsed : %.2fs\n", t2["elapsed"]))
cat(sprintf("  L levels : %d\n", length(a_be$estimate)))
cat(sprintf("  peak Mb  : %.1f Mb (max used since gc(reset=TRUE))\n",
            sum(g2[, "max used"] * c(8, 1024) / (1024 * 1024))))
old_if_bytes <- as.numeric(nrow(fit$X)) * as.numeric(length(a_be$estimate)) * 8
cat(sprintf("  old IF (n x L) allocation would have been: %.2f Gb\n",
            old_if_bytes / 1024^3))
cat(sprintf("  any NaN SE in ATET : %s\n",
            any(is.na(a_be$tidy_table[, "Std. Error"]))))

cat("\n-- atet() overall (cluster) --\n")
t3 <- system.time({ a_ov <- atet(fit, type = "overall") })
cat(sprintf("atet(overall)   elapsed : %.2fs   estimate=%.4f  se=%.4f\n",
            t3["elapsed"], a_ov$estimate, sqrt(a_ov$vcov[1, 1])))

cat("\nDone.\n")
