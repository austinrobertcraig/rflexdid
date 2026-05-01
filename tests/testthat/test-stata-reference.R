# Tests that compare R output to the Stata reference CSVs produced by
# tests/testthat/stata_reference/make-stata-reference.do. Both R and Stata read the same simulated
# dataset (inst/extdata/example_data.csv, regenerated deterministically via
# simulate_flexdid_data()), so these tests verify R/Stata agreement to
# numerical tolerance with no quirk-handling for export-format artifacts.
#
# Tests skip when reference CSVs are absent (e.g. on machines without Stata).

ref_dir <- function() {
  test_path("stata_reference")
}

ref_path <- function(name) {
  file.path(ref_dir(), name)
}

read_ref <- function(name) {
  read.csv(ref_path(name), stringsAsFactors = FALSE)
}

sim_data <- function() {
  path <- system.file("extdata", "example_data.csv", package = "rflexdid")
  if (!nzchar(path)) skip("inst/extdata/example_data.csv not present.")
  read.csv(path)
}

fit_lagsandleads <- function(d) {
  flexdid(ssb_oz ~ 1, data = d, tx = "treated", group = "cohort",
          time = "year", specification = "lagsandleads",
          vcov = "cluster", cluster = "county")
}

fit_lagsonly <- function(d) {
  flexdid(ssb_oz ~ age + female, data = d, tx = "treated",
          group = "county", time = "year",
          specification = "lagsonly", vcov = "cluster", cluster = "county")
}

# ---------- coefficient tests ----------

test_that("lagsandleads coefficients match Stata reference", {
  if (!file.exists(ref_path("coefs_lagsandleads.csv"))) {
    skip("Stata reference CSV not present.")
  }
  ref <- read_ref("coefs_lagsandleads.csv")
  fit <- fit_lagsandleads(sim_data())

  ref_key <- paste(ref$g, ref$t, ref$kind, sep = "|")
  fd_key  <- paste(fit$design$cells$g, fit$design$cells$t,
                   fit$design$cells$kind, sep = "|")
  m  <- match(ref_key, fd_key)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  fd_b  <- unname(fit$coefficients[fit$design$cells$col_indicator[m[ok]]])
  fd_se <- unname(sqrt(diag(fit$vcov)[fit$design$cells$col_indicator[m[ok]]]))
  expect_equal(fd_b,  ref$b[ok],  tolerance = 1e-6)
  expect_equal(fd_se, ref$se[ok], tolerance = 1e-5)
})

test_that("lagsonly cell coefficients match Stata reference", {
  if (!file.exists(ref_path("coefs_lagsonly.csv"))) {
    skip("Stata reference CSV not present.")
  }
  ref <- read_ref("coefs_lagsonly.csv")
  fit <- fit_lagsonly(sim_data())

  ref_key <- paste(ref$g, ref$t, ref$kind, sep = "|")
  fd_key  <- paste(fit$design$cells$g, fit$design$cells$t,
                   fit$design$cells$kind, sep = "|")
  m  <- match(ref_key, fd_key)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  fd_b  <- unname(fit$coefficients[fit$design$cells$col_indicator[m[ok]]])
  fd_se <- unname(sqrt(diag(fit$vcov)[fit$design$cells$col_indicator[m[ok]]]))
  expect_equal(fd_b,  ref$b[ok],  tolerance = 1e-5)
  expect_equal(fd_se, ref$se[ok], tolerance = 1e-5)
})

# ---------- ATET tests, lagsandleads spec ----------

test_that("Overall ATET (lagsandleads) matches Stata reference", {
  if (!file.exists(ref_path("atet_overall.csv"))) {
    skip("Stata reference CSV not present.")
  }
  ref <- read_ref("atet_overall.csv")
  a <- atet(fit_lagsandleads(sim_data()), type = "overall")
  expect_equal(unname(as.numeric(a$estimate)), ref$b[1], tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov)[1])),  ref$se[1], tolerance = 1e-5)
})

test_that("byexposure ATET (lagsandleads) matches Stata reference", {
  if (!file.exists(ref_path("atet_byexposure.csv"))) {
    skip("Stata reference CSV not present.")
  }
  ref <- read_ref("atet_byexposure.csv")
  a <- atet(fit_lagsandleads(sim_data()), type = "byexposure")

  fd_lvl <- as.integer(rownames(a$tidy_table))
  m  <- match(ref$eventtime, fd_lvl)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  expect_equal(unname(a$estimate[m[ok]]),         ref$b[ok],  tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov))[m[ok]]), ref$se[ok], tolerance = 1e-5)
})

test_that("bycalendar ATET (lagsandleads) matches Stata reference", {
  if (!file.exists(ref_path("atet_bycalendar.csv"))) {
    skip("Stata reference CSV not present.")
  }
  ref <- read_ref("atet_bycalendar.csv")
  a <- atet(fit_lagsandleads(sim_data()), type = "bycalendar")

  fd_lvl <- as.numeric(rownames(a$tidy_table))
  m  <- match(ref$t, fd_lvl)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  expect_equal(unname(a$estimate[m[ok]]),         ref$b[ok],  tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov))[m[ok]]), ref$se[ok], tolerance = 1e-5)
})

test_that("bycohort ATET (lagsandleads) matches Stata reference", {
  if (!file.exists(ref_path("atet_bycohort.csv"))) {
    skip("Stata reference CSV not present.")
  }
  ref <- read_ref("atet_bycohort.csv")
  a <- atet(fit_lagsandleads(sim_data()), type = "bycohort")

  fd_lvl <- as.numeric(rownames(a$tidy_table))
  m  <- match(ref$g, fd_lvl)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  expect_equal(unname(a$estimate[m[ok]]),         ref$b[ok],  tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov))[m[ok]]), ref$se[ok], tolerance = 1e-5)
})

test_that("bygroup ATET (lagsandleads) matches Stata reference", {
  if (!file.exists(ref_path("atet_bygroup.csv"))) {
    skip("Stata reference CSV not present.")
  }
  ref <- read_ref("atet_bygroup.csv")
  a <- atet(fit_lagsandleads(sim_data()), type = "bygroup")

  fd_lvl <- as.numeric(rownames(a$tidy_table))
  m  <- match(ref$g, fd_lvl)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  expect_equal(unname(a$estimate[m[ok]]),         ref$b[ok],  tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov))[m[ok]]), ref$se[ok], tolerance = 1e-5)
})

test_that("byget ATET (lagsandleads) matches Stata reference", {
  if (!file.exists(ref_path("atet_byget.csv"))) {
    skip("Stata reference CSV not present.")
  }
  ref <- read_ref("atet_byget.csv")
  a <- atet(fit_lagsandleads(sim_data()), type = "byget")

  # R's byget rownames are formatted as "g<g>|et<et>" (see build_levels()).
  rn <- rownames(a$tidy_table)
  fd_g  <- as.numeric(sub("^g(-?[0-9]+)\\|et(-?[0-9]+)$", "\\1", rn))
  fd_et <- as.numeric(sub("^g(-?[0-9]+)\\|et(-?[0-9]+)$", "\\2", rn))
  fd_key  <- paste(fd_g,   fd_et,         sep = "|")
  ref_key <- paste(ref$g,  ref$eventtime, sep = "|")
  m  <- match(ref_key, fd_key)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  expect_equal(unname(a$estimate[m[ok]]),         ref$b[ok],  tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov))[m[ok]]), ref$se[ok], tolerance = 1e-5)
})

# ---------- ATET tests, lagsonly spec ----------

test_that("Overall ATET (lagsonly) matches Stata reference", {
  if (!file.exists(ref_path("atet_overall_lagsonly.csv"))) {
    skip("Stata reference CSV not present.")
  }
  ref <- read_ref("atet_overall_lagsonly.csv")
  a <- atet(fit_lagsonly(sim_data()), type = "overall")
  expect_equal(unname(as.numeric(a$estimate)), ref$b[1], tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov)[1])),  ref$se[1], tolerance = 1e-5)
})

test_that("byexposure ATET (lagsonly) matches Stata reference", {
  if (!file.exists(ref_path("atet_byexposure_lagsonly.csv"))) {
    skip("Stata reference CSV not present.")
  }
  ref <- read_ref("atet_byexposure_lagsonly.csv")
  a <- atet(fit_lagsonly(sim_data()), type = "byexposure")

  fd_lvl <- as.integer(rownames(a$tidy_table))
  m  <- match(ref$eventtime, fd_lvl)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  expect_equal(unname(a$estimate[m[ok]]),         ref$b[ok],  tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov))[m[ok]]), ref$se[ok], tolerance = 1e-5)
})
