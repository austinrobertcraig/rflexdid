test_that("cohort = min(year) when treated, 0 otherwise", {
  group <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  time  <- c(2030, 2031, 2032, 2030, 2031, 2032, 2030, 2031, 2032)
  tx    <- c(0, 1, 1, 0, 0, 1, 0, 0, 0)

  ch <- flexdid:::make_cohort(group, time, tx)
  expect_equal(ch[group == 1], rep(2031, 3))
  expect_equal(ch[group == 2], rep(2032, 3))
  expect_equal(ch[group == 3], rep(0, 3))   # never-treated
})

test_that("validate_binary catches non-binary tx", {
  expect_error(flexdid:::validate_binary(c(0, 1, 2), "tx"), "binary")
  expect_error(flexdid:::validate_binary(c(1, 1, 1), "tx"), "binary")
  expect_invisible(flexdid:::validate_binary(c(0L, 1L, 0L, 1L), "tx"))
})

test_that("check_time_gaps requires usercohort when gaps present", {
  expect_invisible(flexdid:::check_time_gaps(2030:2034, FALSE))
  expect_error(flexdid:::check_time_gaps(c(2030, 2031, 2033), FALSE),
               "gaps")
  # OK if usercohort supplied
  expect_invisible(flexdid:::check_time_gaps(c(2030, 2031, 2033), TRUE))
})

test_that("make_tx_indicator sets _Tx=1 in pre-period for treated cohorts", {
  cohort <- c(2032, 2032, 2032, 2032, 2032)
  time   <- c(2029, 2030, 2031, 2032, 2033)
  tx     <- c(0,    0,    0,    1,    1)   # base = 2031
  out <- flexdid:::make_tx_indicator(tx, cohort, time)
  # pre-base periods (year < 2031) -> 1
  # base year (2031) -> 0
  # post (>=2032) -> 1
  expect_equal(out, c(1L, 1L, 0L, 1L, 1L))
})

test_that("handle_treatment_status drops last cohort when no never-treated", {
  cohort <- c(2031, 2031, 2032, 2032, 2033, 2033)
  time   <- c(2030, 2031, 2030, 2032, 2030, 2033)
  res <- flexdid:::handle_treatment_status(cohort, time, verbose = FALSE)
  # last cohort 2033 should be reset to 0; obs at year >= 2033 dropped
  expect_equal(res$cohort[5:6], c(0, 0))
  expect_equal(res$keep, c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))
})

test_that("handle_treatment_status errors on always-treated", {
  cohort <- c(2030, 2030, 2031, 2031)
  time   <- c(2030, 2031, 2030, 2031)
  expect_error(flexdid:::handle_treatment_status(cohort, time, verbose = FALSE),
               "always-treated")
})
