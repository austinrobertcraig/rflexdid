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
  path <- system.file("extdata", "hhabits.csv", package = "rflexdid")
  if (path == "") path <- test_path("..", "..", "inst", "extdata", "hhabits.csv")
  if (!file.exists(path)) {
    skip("hhabits.csv not present; run data-raw/make-stata-reference.do first.")
  }
  d <- read.csv(path, stringsAsFactors = FALSE)

  # If the CSV was exported without `nolabel`, the labeled variables come
  # through as character "Yes"/"No". Map back to integers so flexdid() can
  # consume the data unchanged.
  yesno <- function(z) {
    if (is.character(z) && all(z[!is.na(z)] %in% c("Yes", "No"))) {
      as.integer(z == "Yes")
    } else z
  }
  for (nm in c("hhabit", "girl", "sports")) {
    if (nm %in% names(d)) d[[nm]] <- yesno(d[[nm]])
  }

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
  # Stata's interaction names look like:
  #   "2034bn._Cohort#2034.chrt#2032b.year#1._Tx"
  # The factor-modifier suffixes (`bn`, `b`, `o`) make a structured regex
  # finicky; just pull the first three numeric tokens, which encode
  # (cohort, group, year) in that order.
  parse_stata_cell <- function(nm) {
    nums <- regmatches(nm, gregexpr("[0-9]+", nm))[[1]]
    if (length(nums) < 3L) return(c(NA, NA, NA))
    as.numeric(nums[1:3])
  }
  ref_keys <- t(vapply(ref_cell$name, parse_stata_cell, numeric(3)))
  # In this spec group=chrt and cohort=chrt, so columns 1 (cohort) and 2 (group)
  # are identical in the Stata names. Match on (group, year) which is what
  # uniquely identifies a flexdid cell anyway.
  ref_key_str <- paste(ref_keys[, 2], ref_keys[, 3], sep = "|")
  fd_key_str  <- paste(fit$design$cells$g, fit$design$cells$t, sep = "|")

  m <- match(ref_key_str, fd_key_str)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  fd_b  <- unname(fit$coefficients[fit$design$cells$col_indicator[m[ok]]])
  fd_se <- unname(sqrt(diag(fit$vcov)[fit$design$cells$col_indicator[m[ok]]]))
  expect_equal(fd_b, ref_cell$b[ok], tolerance = 1e-6)
  expect_equal(fd_se, ref_cell$se[ok], tolerance = 1e-5)
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
  expect_equal(unname(as.numeric(a$estimate)), ref$b[1], tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov)[1])), ref$se[1], tolerance = 1e-5)
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

  # ref$label is the event-time number (possibly negative), maybe with a
  # trailing "bn" / "o" Stata factor modifier. Pull the leading signed int.
  ref_lvl <- suppressWarnings(as.integer(sub("^(-?[0-9]+).*$", "\\1",
                                              as.character(ref$label))))
  fd_lvl  <- suppressWarnings(as.integer(rownames(a$tidy_table)))
  m <- match(ref_lvl, fd_lvl)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  expect_equal(unname(a$estimate[m[ok]]), ref$b[ok], tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov))[m[ok]]), ref$se[ok], tolerance = 1e-5)
})

test_that("bycalendar ATET (lagsandleads) matches Stata reference", {
  if (!file.exists(file.path(ref_dir(), "atet_bycalendar.csv"))) {
    skip("Stata reference CSV not present.")
  }
  d <- read_data()
  ref <- read_atet_ref("atet_bycalendar.csv")
  fit <- flexdid(bmi ~ 1, data = d, tx = "hhabit", group = "chrt",
                 time = "year", specification = "lagsandleads",
                 vcov = "cluster", cluster = "schools")
  a <- atet(fit, type = "bycalendar")

  ref_lvl <- suppressWarnings(as.numeric(ref$label))
  fd_lvl  <- suppressWarnings(as.numeric(rownames(a$tidy_table)))
  m <- match(ref_lvl, fd_lvl)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  expect_equal(unname(a$estimate[m[ok]]), ref$b[ok], tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov))[m[ok]]), ref$se[ok], tolerance = 1e-5)
})

test_that("bycohort ATET (lagsandleads) matches Stata reference", {
  if (!file.exists(file.path(ref_dir(), "atet_bycohort.csv"))) {
    skip("Stata reference CSV not present.")
  }
  d <- read_data()
  ref <- read_atet_ref("atet_bycohort.csv")
  fit <- flexdid(bmi ~ 1, data = d, tx = "hhabit", group = "chrt",
                 time = "year", specification = "lagsandleads",
                 vcov = "cluster", cluster = "schools")
  a <- atet(fit, type = "bycohort")

  # Stata labels cohorts as ordinal integers (1, 2, 3); R uses the actual
  # cohort values sorted ascending. Match positionally after excluding any
  # zero-SE rows (Stata marks base/omitted cohorts with se=0).
  ok_ref <- ref$se > 0
  ok_fd  <- unname(sqrt(diag(a$vcov))) > 0
  expect_equal(sum(ok_ref), sum(ok_fd))

  expect_equal(unname(a$estimate[ok_fd]),  ref$b[ok_ref],  tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov))[ok_fd]), ref$se[ok_ref], tolerance = 1e-5)
})

test_that("bygroup ATET (lagsandleads) matches Stata reference", {
  if (!file.exists(file.path(ref_dir(), "atet_bygroup.csv"))) {
    skip("Stata reference CSV not present.")
  }
  d <- read_data()
  ref <- read_atet_ref("atet_bygroup.csv")
  fit <- flexdid(bmi ~ 1, data = d, tx = "hhabit", group = "chrt",
                 time = "year", specification = "lagsandleads",
                 vcov = "cluster", cluster = "schools")
  a <- atet(fit, type = "bygroup")

  ref_lvl <- suppressWarnings(as.numeric(ref$label))
  fd_lvl  <- suppressWarnings(as.numeric(rownames(a$tidy_table)))
  m <- match(ref_lvl, fd_lvl)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  expect_equal(unname(a$estimate[m[ok]]), ref$b[ok], tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov))[m[ok]]), ref$se[ok], tolerance = 1e-5)
})

test_that("byget ATET (lagsandleads) matches Stata reference", {
  if (!file.exists(file.path(ref_dir(), "atet_byget.csv"))) {
    skip("Stata reference CSV not present.")
  }
  d <- read_data()
  ref <- read_atet_ref("atet_byget.csv")
  fit <- flexdid(bmi ~ 1, data = d, tx = "hhabit", group = "chrt",
                 time = "year", specification = "lagsandleads",
                 vcov = "cluster", cluster = "schools")
  a <- atet(fit, type = "byget")

  # Stata's dump_atet writes nrow(ref) rows where row i carries b = B[1,i]
  # (i-th cell ATET) and label = word i of `colnames r(b)`.  Because each
  # (group, eventtime) label is a two-word string ("2034 -2"), the word
  # positions are offset from the cell positions: the labels are wrong
  # identifiers, but the b values are correct and in cell order.
  # R's `atet(, type="byget")` sorts cells by (group, eventtime) ascending and
  # includes the base period (et=-1, ATET=0, SE=0).  Remove those base-period
  # rows and compare sequentially.
  fd_se <- unname(sqrt(diag(a$vcov)))
  keep  <- fd_se > 0
  expect_equal(sum(keep), nrow(ref))

  expect_equal(unname(a$estimate[keep]), ref$b, tolerance = 1e-5)
  expect_equal(fd_se[keep], ref$se, tolerance = 1e-5)
})

test_that("lagsonly cell coefficients match Stata reference (group=schools)", {
  if (!file.exists(file.path(ref_dir(), "coefs_lagsonly.csv"))) {
    skip("Stata reference CSV not present.")
  }
  d <- read_data()
  ref <- read.csv(file.path(ref_dir(), "coefs_lagsonly.csv"),
                  stringsAsFactors = FALSE)

  fit <- flexdid(bmi ~ girl + medu, data = d, tx = "hhabit",
                 group = "schools", time = "year",
                 specification = "lagsonly", vcov = "cluster")

  # Restrict to the indicator-only cells (exclude covariate-interacted cells
  # whose Stata names contain "_Tx#").  With covariates, "_Cohort#" also
  # matches x-interacted cells (e.g. "..._Tx#1.girl"), which would produce
  # duplicate (group, year) keys and misalign the comparison.
  is_cell <- grepl("_Cohort#", ref$name) & !grepl("_Tx#", ref$name)
  ref_cell <- ref[is_cell, , drop = FALSE]
  expect_gt(nrow(ref_cell), 0L)

  parse_stata_cell <- function(nm) {
    nums <- regmatches(nm, gregexpr("[0-9]+", nm))[[1]]
    if (length(nums) < 3L) return(c(NA, NA, NA))
    as.numeric(nums[1:3])
  }
  ref_keys <- t(vapply(ref_cell$name, parse_stata_cell, numeric(3)))
  ref_key_str <- paste(ref_keys[, 2], ref_keys[, 3], sep = "|")
  fd_key_str  <- paste(fit$design$cells$g, fit$design$cells$t, sep = "|")

  m <- match(ref_key_str, fd_key_str)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  fd_b  <- unname(fit$coefficients[fit$design$cells$col_indicator[m[ok]]])
  fd_se <- unname(sqrt(diag(fit$vcov)[fit$design$cells$col_indicator[m[ok]]]))
  expect_equal(fd_b, ref_cell$b[ok], tolerance = 1e-5)
  expect_equal(fd_se, ref_cell$se[ok], tolerance = 1e-5)
})

test_that("Overall ATET (lagsonly) matches Stata reference", {
  if (!file.exists(file.path(ref_dir(), "atet_overall_lagsonly.csv"))) {
    skip("Stata reference CSV not present.")
  }
  d <- read_data()
  ref <- read_atet_ref("atet_overall_lagsonly.csv")
  fit <- flexdid(bmi ~ girl + medu, data = d, tx = "hhabit",
                 group = "schools", time = "year",
                 specification = "lagsonly", vcov = "cluster")
  a <- atet(fit, type = "overall")
  expect_equal(unname(as.numeric(a$estimate)), ref$b[1], tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov)[1])), ref$se[1], tolerance = 1e-5)
})

test_that("byexposure ATET (lagsonly) matches Stata reference", {
  if (!file.exists(file.path(ref_dir(), "atet_byexposure_lagsonly.csv"))) {
    skip("Stata reference CSV not present.")
  }
  d <- read_data()
  ref <- read_atet_ref("atet_byexposure_lagsonly.csv")
  fit <- flexdid(bmi ~ girl + medu, data = d, tx = "hhabit",
                 group = "schools", time = "year",
                 specification = "lagsonly", vcov = "cluster")
  a <- atet(fit, type = "byexposure")

  # Stata shows pre-treatment rows with se=0 (omitted); filter them out.
  ok_ref <- ref$se > 0
  ref_lvl <- suppressWarnings(as.integer(sub("^(-?[0-9]+).*$", "\\1",
                                             as.character(ref$label[ok_ref]))))
  fd_lvl  <- suppressWarnings(as.integer(rownames(a$tidy_table)))
  m <- match(ref_lvl, fd_lvl)
  ok <- !is.na(m)
  expect_gt(sum(ok), 0L)

  expect_equal(unname(a$estimate[m[ok]]), ref$b[ok_ref][ok], tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(a$vcov))[m[ok]]), ref$se[ok_ref][ok], tolerance = 1e-5)
})
