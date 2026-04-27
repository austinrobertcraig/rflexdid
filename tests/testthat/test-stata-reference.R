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
