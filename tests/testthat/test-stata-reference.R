# Tests that compare R output to the Stata reference CSVs produced by
# data-raw/make-stata-reference.do. They skip when the CSVs are not present
# (e.g. on machines without Stata).

ref_dir <- function() {
  test_path("stata_reference")
}

read_atet_ref <- function(name) {
  path <- file.path(ref_dir(), name)
  read.csv(path, stringsAsFactors = FALSE)
}

read_data <- function() {
  path <- system.file("extdata", "hhabits.csv", package = "flexdid")
  if (path == "") path <- test_path("..", "..", "inst", "extdata", "hhabits.csv")
  if (!file.exists(path)) {
    skip("hhabits.csv not present; run data-raw/make-stata-reference.do first.")
  }
  d <- read.csv(path, stringsAsFactors = FALSE)
  if (!"chrt" %in% names(d)) {
    d$chrt <- ave(ifelse(d$hhabit == 1, d$year, NA), d$schools,
                  FUN = function(z) suppressWarnings(min(z, na.rm = TRUE)))
    d$chrt[!is.finite(d$chrt)] <- 0
  }
  d
}

test_that("lagsandleads coefficients match Stata reference (group=chrt, cluster=schools)", {
  if (!file.exists(file.path(ref_dir(), "coefs_lagsandleads.csv"))) {
    skip("Stata reference CSV not present.")
  }
  d <- read_data()
  ref <- read.csv(file.path(ref_dir(), "coefs_lagsandleads.csv"),
                  stringsAsFactors = FALSE)

  fit <- flexdid(bmi ~ 1, data = d, tx = "hhabit", group = "chrt",
                 time = "year", specification = "lagsandleads",
                 vcov = "cluster", cluster = "schools")

  # Keep only the treatment-cell coefficients in the reference: those names
  # contain "_Cohort#" in Stata. Map them onto our naming.
  is_cell <- grepl("_Cohort#", ref$name)
  ref_cell <- ref[is_cell, , drop = FALSE]
  expect_gt(nrow(ref_cell), 0L)
  # Mapping is not 1:1 in name string but is 1:1 in (cohort, group, year);
  # do a robust match via parsing.
  parse_stata_cell <- function(nm) {
    # e.g. "2033bn._Cohort#1.chrt#2030.year#1._Tx" -> list(cohort=2033, g=1, t=2030)
    m <- regmatches(nm, regexec(
      "([0-9]+)[a-z]*\\._Cohort#([0-9]+)\\.[A-Za-z]+#([0-9]+)\\.year",
      nm
    ))[[1]]
    if (length(m) < 4L) return(c(NA, NA, NA))
    as.numeric(m[c(2, 3, 4)])
  }
  ref_keys <- t(vapply(ref_cell$name, parse_stata_cell, numeric(3)))
  fd_keys <- with(fit$design$cells, cbind(cohort_for_g = NA, g = g, t = t))
  # cohort_for_g lookup
  cohort_lookup <- tapply(d$chrt, d$schools, function(z) z[1])
  fd_keys[, "cohort_for_g"] <- as.numeric(cohort_lookup[as.character(fit$design$cells$g)])

  ref_key_str <- paste(ref_keys[, 1], ref_keys[, 2], ref_keys[, 3], sep = "|")
  fd_key_str  <- paste(fd_keys[, "cohort_for_g"], fd_keys[, "g"], fd_keys[, "t"],
                       sep = "|")
  m <- match(ref_key_str, fd_key_str)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  fd_b <- fit$coefficients[fit$design$cells$col_indicator[m[ok]]]
  fd_se <- sqrt(diag(fit$vcov)[fit$design$cells$col_indicator[m[ok]]])
  expect_equal(unname(fd_b), ref_cell$b[ok], tolerance = 1e-6)
  expect_equal(unname(fd_se), ref_cell$se[ok], tolerance = 1e-5)
})

test_that("Overall ATET (lagsandleads) matches Stata reference", {
  if (!file.exists(file.path(ref_dir(), "atet_overall.csv"))) {
    skip("Stata reference CSV not present.")
  }
  d <- read_data()
  ref <- read_atet_ref("atet_overall.csv")
  fit <- flexdid(bmi ~ 1, data = d, tx = "hhabit", group = "chrt",
                 time = "year", specification = "lagsandleads",
                 vcov = "cluster", cluster = "schools")
  a <- atet(fit, type = "overall")
  expect_equal(as.numeric(a$estimate), ref$b[1], tolerance = 1e-5)
  expect_equal(sqrt(diag(a$vcov)[1]), ref$se[1], tolerance = 1e-5)
})

test_that("byexposure ATET (lagsandleads) matches Stata reference", {
  if (!file.exists(file.path(ref_dir(), "atet_byexposure.csv"))) {
    skip("Stata reference CSV not present.")
  }
  d <- read_data()
  ref <- read_atet_ref("atet_byexposure.csv")
  fit <- flexdid(bmi ~ 1, data = d, tx = "hhabit", group = "chrt",
                 time = "year", specification = "lagsandleads",
                 vcov = "cluster", cluster = "schools")
  a <- atet(fit, type = "byexposure")

  # ref$label is like "0.bn", "1.bn", ... or numeric event-time labels
  # depending on Stata version. Parse the leading integer.
  ref_lvl <- suppressWarnings(as.integer(sub("^([\\-0-9]+).*", "\\1", ref$label)))
  fd_lvl  <- suppressWarnings(as.integer(rownames(a$tidy_table)))
  m <- match(ref_lvl, fd_lvl)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  expect_equal(unname(a$estimate[m[ok]]), ref$b[ok], tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov))[m[ok]]), ref$se[ok], tolerance = 1e-5)
})
