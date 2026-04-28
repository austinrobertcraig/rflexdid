# rflexdid

> **Disclaimer:** This is a personal research tool, validated only against my own use cases and maintained only as needed for my own work. If you use it, verify your results independently.

R port of the Stata [`flexdid`](https://ideas.repec.org/c/boc/bocode/s459302.html)
command (Deb, Norton, Wooldridge & Zabel 2025; SSC). Estimates a flexible OLS
regression with cohort × group × time × treatment interactions for
staggered-timing difference-in-differences designs, then aggregates the
resulting cell-level effects to ATETs (overall, by exposure time, by calendar
time, by cohort, by group, or by group × exposure time). Standard errors
match Stata's `vce(unconditional)` via influence functions.

## Install

This package is not on CRAN. Install once from this GitHub repo:

```r
# install.packages("remotes")
remotes::install_github("austinrobertcraig/rflexdid")
```

Or, while iterating locally:

```r
remotes::install_local("~/GitHub/rflexdid")
```

`rflexdid` depends on `Matrix`, `sandwich`, and base R `stats`. `ggplot2` is
suggested (used only in `plot()`).

## Usage

```r
library(rflexdid)

# `hhabits` here is the standard Stata example dataset. Export it once via
# `webuse hhabits` and the data-raw/make-stata-reference.do script ships it
# to inst/extdata. From R:
hh <- read.csv(system.file("extdata", "hhabits.csv", package = "rflexdid"))

# Build the cohort variable used in the help-file examples.
hh$chrt <- with(hh, ave(ifelse(hhabit == 1, year, NA), schools,
                        FUN = function(z) suppressWarnings(min(z, na.rm = TRUE))))
hh$chrt[!is.finite(hh$chrt)] <- 0

# FLEX regression: lags-and-leads spec, cluster-robust SEs.
fit <- flexdid(bmi ~ 1, data = hh, tx = "hhabit",
               group = "chrt", time = "year",
               specification = "lagsandleads",
               vcov = "cluster", cluster = "schools")

summary(fit)

# Postestimation aggregations (analogues of `estat atet`).
atet(fit, type = "overall")
atet(fit, type = "byexposure", test = "zero")
atet(fit, type = "bycalendar")
atet(fit, type = "bycohort")
atet(fit, type = "bygroup")
atet(fit, type = "byget", values = -3:3)

# Event-study plot.
library(ggplot2)
plot(atet(fit, type = "byexposure"))

# Subgroup ATET (Stata's `for(...)` option).
atet(fit, type = "overall", for_expr = quote(girl == 1))
```

## Mapping from Stata to R

| Stata                                               | R                                                        |
|-----------------------------------------------------|----------------------------------------------------------|
| `flexdid bmi, tx(hhabit) group(s) time(year)`       | `flexdid(bmi ~ 1, data, tx="hhabit", group="s", time="year")` |
| `flexdid bmi medu girl, ...`                        | `flexdid(bmi ~ medu + girl, ...)`                        |
| `specification(lagsandleads)`                       | `specification = "lagsandleads"`                         |
| `xnotinteracted(medu i.s)`                          | `xnotinteracted = ~ medu + factor(s)`                    |
| `usercohort(chrt)`                                  | `usercohort = "chrt"`                                    |
| `vce(cluster s)`                                    | `vcov = "cluster", cluster = "s"`                        |
| `[pweight=w]`                                       | `weights = "w"`                                          |
| `estat atet, overall`                               | `atet(fit, type = "overall")`                            |
| `estat atet, byexposure(0/3) test(zero)`            | `atet(fit, type = "byexposure", values = 0:3, test = "zero")` |
| `estat atet, byexposure for(girl==1)`               | `atet(fit, type = "byexposure", for_expr = quote(girl==1))` |
| `estat atet, byget`                                 | `atet(fit, type = "byget")`                              |
| `aggregationweight(grouplevel)`                     | `aggregationweight = "grouplevel"`                       |

The R object returned by `flexdid()` carries the full design matrix and
sample-aligned vectors needed by `atet()`. `summary()`, `coef()`, `vcov()`,
and `print()` work as expected.

## Standard errors

Standard errors on the regression coefficients use the same HC1 / CR1
sandwich formulas as Stata's `regress, vce(robust)` and
`regress, vce(cluster)`. Standard errors on the aggregated ATETs use the
influence-function form that mirrors Stata's `margins, vce(unconditional)`
— see [R/atet.R](R/atet.R) for the explicit formula.

## Validating against Stata

To regenerate the reference CSVs that the testthat suite compares against,
run once from a Stata installation that has `flexdid` (>= v2.0):

```bash
cd ~/GitHub/rflexdid
stata-mp -b do data-raw/make-stata-reference.do
```

Then, in R:

```r
devtools::test()
```

Without the CSVs the Stata-reference tests skip; the internal numerical
tests (24 of them) run regardless.

## References

Deb, P., Norton, E. C., Wooldridge, J. M., Zabel, J. E. (2025), "A Flexible,
Heterogeneous Treatment Effects Difference-in-Differences Estimator for
Repeated Cross-Sections", NBER Working Paper No. 33026.
