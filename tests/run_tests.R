# Run the testthat suite and write a detailed log to tests/test_results.txt
# (overwritten on each invocation). Usage from the repo root:
#
#     Rscript tests/run_tests.R
#
# Exits non-zero if any test fails or errors.

suppressMessages({
  if (!requireNamespace("devtools",  quietly = TRUE)) install.packages("devtools")
  if (!requireNamespace("testthat",  quietly = TRUE)) install.packages("testthat")
  if (!requireNamespace("here",      quietly = TRUE)) install.packages("here")
})

repo_root <- here::here()
out_path  <- file.path(repo_root, "tests", "test_results.txt")

devtools::load_all(repo_root, quiet = TRUE)

# Capture all output (stdout + messages from the SummaryReporter) into the
# log file while also echoing to the console so the user sees live progress.
log_con <- file(out_path, open = "wt")
on.exit({ try(close(log_con), silent = TRUE) }, add = TRUE)

sink(log_con, split = TRUE)
sink(log_con, type = "message", append = TRUE)

result <- testthat::test_dir(
  file.path(repo_root, "tests", "testthat"),
  reporter = testthat::SummaryReporter$new(),
  stop_on_failure = FALSE
)

sink(type = "message")
sink()
close(log_con)
on.exit()

message("Detailed results written to ", out_path)

df <- as.data.frame(result)
n_failed <- sum(df$failed) + sum(df$error)
if (n_failed > 0L) quit(status = 1L)
