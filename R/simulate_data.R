# Simulated dataset for the rflexdid demo vignette and Stata-reference tests.
#
# Setting: counties enact a sugar-sweetened-beverage (SSB) excise tax in
# different years. Each year an annual BRFSS-style survey draws a fresh
# cross section of adult residents and asks about daily SSB consumption.
# Treatment is assigned at the county level; observations are at the
# individual level (no individual is followed across years).
#
# The DGP has strict parallel pre-trends and known cohort-heterogeneous,
# dynamic treatment effects (early-cohort counties enacted larger taxes,
# and consumption falls smoothly with exposure time). See the docstring
# of simulate_flexdid_data() for the exact functional form.
#
# Usage (from the package root):
#   devtools::load_all()
#   simulate_flexdid_data()   # writes inst/extdata/example_data.csv

#' @noRd
simulate_flexdid_data <- function(n_counties = 40,
                                  years      = 2010:2019,
                                  cohorts    = c(2013, 2015, 2017, 0),
                                  n_per_cell = 50,
                                  seed       = 20260430) {
  set.seed(seed)

  # 1. Assign each county to a cohort (roughly equal shares; 0 = never-treated).
  cohort_assign <- sample(rep(cohorts, length.out = n_counties))

  # 2. Build the county x year x respondent grid.
  grid <- expand.grid(county = seq_len(n_counties),
                      year   = years,
                      person = seq_len(n_per_cell))
  grid$cohort  <- cohort_assign[grid$county]
  grid$treated <- as.integer(grid$cohort > 0 & grid$year >= grid$cohort)

  # 3. Individual-level covariates (drawn fresh per row).
  N <- nrow(grid)
  grid$age    <- pmin(80, pmax(18, stats::rnorm(N, mean = 40, sd = 15)))
  grid$female <- stats::rbinom(N, 1, 0.5)

  # 4. County-level region (constant within county).
  region_assign <- sample(1:4, n_counties, replace = TRUE)
  grid$region   <- factor(region_assign[grid$county],
                          levels = 1:4,
                          labels = c("Northeast", "Midwest", "South", "West"))
  region_effect <- c(Northeast = -0.5, Midwest = 0.0, South = 0.5, West = 1.0)

  # 5. Treatment effect by cohort: alpha_g - 0.20 * e - 0.02 * e^2.
  alpha_g <- c(`2013` = -1.5, `2015` = -1.0, `2017` = -0.5, `0` = 0)
  e   <- grid$year - grid$cohort
  tau <- ifelse(grid$treated == 1L,
                alpha_g[as.character(grid$cohort)] - 0.20 * e - 0.02 * e^2,
                0)

  # 6. Outcome: daily SSB consumption in fluid ounces.
  grid$ssb_oz <- 12.0 +
    -0.10 * (grid$year - 2010) +
     0.20 * (grid$county - 20) / 20 +
     0.50 * region_effect[as.character(grid$region)] +
    -0.04 * (grid$age - 40) +
    -1.50 * grid$female +
     tau +
     stats::rnorm(N, sd = 3.0)

  # Drop the helper person index (repeated cross section, no individual id).
  grid$person <- NULL

  # Reorder columns so the panel "skeleton" is up front.
  df <- grid[, c("county", "year", "cohort", "treated",
                 "age", "female", "region", "ssb_oz")]

  out_path <- "inst/extdata/example_data.csv"
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, out_path, row.names = FALSE)
  message("wrote ", nrow(df), " rows to ", out_path)
  invisible(NULL)
}
