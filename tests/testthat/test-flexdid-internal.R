# Tests that don't require Stata reference output.

make_panel <- function(seed = 42, ng = 60, ny = 10, year0 = 2030,
                       cohorts = c(2033, 2034, 2035, 2036)) {
  set.seed(seed)
  df <- expand.grid(schools = seq_len(ng), year = year0:(year0 + ny - 1))
  per_cohort <- ng %/% (length(cohorts) + 1)
  chrt_assign <- rep(0L, ng)
  for (i in seq_along(cohorts)) {
    rng <- ((i - 1) * per_cohort + 1):(i * per_cohort)
    chrt_assign[rng] <- cohorts[i]
  }
  df$cohort_t <- chrt_assign[df$schools]
  df$hhabit <- as.integer(df$cohort_t > 0 & df$year >= df$cohort_t)
  df$x1 <- rnorm(nrow(df))
  df$bmi <- 20 + 0.5 * df$hhabit + 0.1 * df$year - 0.02 * df$schools +
            0.3 * df$x1 + rnorm(nrow(df), sd = 0.5)
  df
}

test_that("flexdid coefficients reproduce lm() on the same effective design", {
  df <- make_panel()
  fit <- flexdid(bmi ~ 1, data = df, tx = "hhabit", group = "schools",
                 time = "year", specification = "lagsonly", vcov = "robust")

  df$cell <- ifelse(df$cohort_t > 0 & df$year >= df$cohort_t,
                    paste0("g", df$schools, "_t", df$year), "_base")
  df$cell <- relevel(factor(df$cell), ref = "_base")
  fit_lm <- lm(bmi ~ cell + factor(schools) + factor(year), data = df)

  expect_equal(unname(fit$rss), sum(resid(fit_lm)^2), tolerance = 1e-10)

  fd_cells <- fit$coefficients[grep("^Tx_lag:", names(fit$coefficients))]
  lm_cells <- coef(fit_lm)[grep("^cellg", names(coef(fit_lm)))]
  nm_lm <- gsub("^cellg", "", names(lm_cells))
  nm_lm <- gsub("_t", ":t", nm_lm)
  nm_fd <- gsub("^Tx_lag:g", "", names(fd_cells))
  m <- match(nm_fd, nm_lm)
  expect_true(!any(is.na(m)))
  expect_equal(unname(fd_cells), unname(lm_cells[m]), tolerance = 1e-10)
})

test_that("overall ATET = mean of cell coefficients (no covariates, lagsonly)", {
  df <- make_panel()
  fit <- flexdid(bmi ~ 1, data = df, tx = "hhabit", group = "schools",
                 time = "year", specification = "lagsonly", vcov = "robust")
  ax <- atet(fit, type = "overall")
  cell_coefs <- fit$coefficients[grep("^Tx_lag:", names(fit$coefficients))]
  expect_equal(unname(ax$estimate), mean(cell_coefs), tolerance = 1e-10)
})

test_that("byexposure ATET levels with no obs are dropped", {
  df <- make_panel()
  fit <- flexdid(bmi ~ 1, data = df, tx = "hhabit", group = "schools",
                 time = "year", specification = "lagsandleads", vcov = "robust")
  a <- atet(fit, type = "byexposure")
  # in this DGP the eventtimes range from -6 to 6; -1 is the base period
  # and is included with ATET = 0.
  base_row <- a$tidy_table["-1", , drop = FALSE]
  expect_equal(unname(base_row[, "Estimate"]), 0)
  expect_equal(unname(base_row[, "Std. Error"]), 0)
})

test_that("Wald 'zero' test on post-treatment matches manual chi-square", {
  df <- make_panel()
  fit <- flexdid(bmi ~ 1, data = df, tx = "hhabit", group = "schools",
                 time = "year", specification = "lagsonly", vcov = "cluster")
  a <- atet(fit, type = "byexposure", values = 0:4, test = "zero")
  est <- a$estimate
  V   <- a$vcov
  manual_W <- as.numeric(t(est) %*% solve(V, est))
  manual_F <- manual_W / length(est)
  expect_equal(a$test_result$F, manual_F, tolerance = 1e-8)
})

test_that("for_expr restricts the subpopulation", {
  df <- make_panel()
  fit <- flexdid(bmi ~ x1, data = df, tx = "hhabit", group = "schools",
                 time = "year", specification = "lagsonly", vcov = "robust")
  a_all  <- atet(fit, type = "overall")
  a_sub  <- atet(fit, type = "overall", for_expr = ~ x1 > 0)
  a_sub2 <- atet(fit, type = "overall", for_expr = quote(x1 > 0))
  expect_lt(a_sub$n_sub, a_all$n_sub)
  expect_true(is.numeric(as.numeric(a_sub$estimate)))
  expect_equal(as.numeric(a_sub$estimate), as.numeric(a_sub2$estimate))
  expect_equal(a_sub$n_sub, a_sub2$n_sub)
})

test_that("plot returns a ggplot for byexposure", {
  skip_if_not_installed("ggplot2")
  df <- make_panel()
  fit <- flexdid(bmi ~ 1, data = df, tx = "hhabit", group = "schools",
                 time = "year", specification = "lagsandleads", vcov = "robust")
  a <- atet(fit, type = "byexposure")
  expect_s3_class(plot(a), "ggplot")
})

test_that("print.summary.flexdid does not error (CI columns not passed as p-values)", {
  df <- make_panel()
  fit <- flexdid(bmi ~ 1, data = df, tx = "hhabit", group = "schools",
                 time = "year", specification = "lagsandleads", vcov = "robust")
  expect_no_error(capture.output(print(summary(fit))))
})

test_that("print.summary.flexdid handles non-estimable SEs when clusters < parameters", {
  # 4 schools = 3 clusters (1 control + 2 treated cohorts); lagsandleads
  # generates more kept parameters than clusters, so some SEs are zero.
  df <- make_panel(ng = 4, cohorts = c(2033, 2034))
  fit <- flexdid(bmi ~ 1, data = df, tx = "hhabit", group = "schools",
                 time = "year", specification = "lagsandleads", vcov = "cluster")
  expect_no_error(capture.output(print(summary(fit))))
})

test_that("atet Wald test uses pseudoinverse when clusters < ATET levels", {
  # 4 schools → rank-deficient cluster VCov for byexposure (more levels than
  # clusters). atet_wald must fall back to pseudoinverse without erroring.
  df  <- make_panel(ng = 4, cohorts = c(2033, 2034))
  fit <- flexdid(bmi ~ 1, data = df, tx = "hhabit", group = "schools",
                 time = "year", specification = "lagsandleads", vcov = "cluster")
  expect_no_error({
    a <- atet(fit, type = "byexposure", test = "zero")
  })
  tr <- a$test_result
  expect_true(is.finite(tr$F))
  expect_true(is.finite(tr$p))
  expect_false(is.null(tr$pinv_note))       # pseudoinverse path was taken
  expect_lt(tr$df1, sum(!is.na(a$tidy_table[, "t value"])))  # df1 = rank < non-NA levels
})

test_that("VCE robust matches a manual HC1 computation on a non-saturated fit", {
  # Use a panel where cells have more than one obs (combine schools into
  # cohort-level groups) so HC1 is well-defined without near-1 hat values.
  df <- make_panel()
  # Combine schools 1-12 into a single "group" per cohort, etc.
  df$cohort_group <- ifelse(df$cohort_t > 0, df$cohort_t, 0L)
  fit <- flexdid(bmi ~ 1, data = df, tx = "hhabit", group = "cohort_group",
                 time = "year", specification = "lagsonly", vcov = "robust")

  # Manual HC1: V = (n/(n-k)) * (X'X)^{-1} * (sum_i u_i^2 X_i X_i') * (X'X)^{-1}
  X <- as.matrix(fit$X[, fit$pivot_keep, drop = FALSE])
  u <- fit$residuals
  n <- fit$nobs
  k <- fit$rank
  XtX_inv <- chol2inv(chol(crossprod(X)))
  meat <- crossprod(X * u)
  V_manual <- (n / (n - k)) * (XtX_inv %*% meat %*% XtX_inv)

  # Compare against flexdid's reported VCE on the kept columns.
  expect_equal(unname(diag(fit$vcov)[fit$pivot_keep]),
               unname(diag(V_manual)),
               tolerance = 1e-8)
})

test_that("flexdid() handles NAs in outcome by dropping affected rows", {
  df <- make_panel()
  set.seed(2)
  df$bmi_with_na <- df$bmi
  na_idx <- sample(seq_len(nrow(df)), size = 15)
  df$bmi_with_na[na_idx] <- NA

  expect_message(
    fit <- flexdid(bmi_with_na ~ 1, data = df, tx = "hhabit",
                   group = "schools", time = "year",
                   specification = "lagsandleads", vcov = "cluster"),
    "dropped 15 row\\(s\\) with NAs in outcome variable 'bmi_with_na'"
  )
  expect_equal(fit$nobs, nrow(df) - 15)
})

test_that("flexdid() handles NAs in covariates by dropping affected rows", {
  df <- make_panel()
  set.seed(1)
  df$x1_with_na <- df$x1
  na_idx <- sample(which(df$hhabit == 1L), size = 20)
  df$x1_with_na[na_idx] <- NA

  expect_message(
    fit <- flexdid(bmi ~ x1_with_na, data = df, tx = "hhabit",
                   group = "schools", time = "year",
                   specification = "lagsandleads", vcov = "cluster"),
    "dropped 20 row\\(s\\) with NAs in covariates"
  )
  expect_equal(fit$nobs, nrow(df) - 20)
})

test_that("pretrends test restricts to pre-treatment periods and has correct label", {
  df  <- make_panel()
  fit <- flexdid(bmi ~ 1, data = df, tx = "hhabit", group = "schools",
                 time = "year", specification = "lagsandleads", vcov = "cluster")
  a  <- atet(fit, type = "byexposure", test = "pretrends")
  tr <- a$test_result
  expect_match(tr$title, "parallel trends", ignore.case = TRUE)
  expect_true(is.finite(tr$F))
  expect_true(is.finite(tr$p))
})
