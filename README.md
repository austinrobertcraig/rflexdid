# rflexdid

> **Disclaimer:** This is a personal research tool maintained only as needed for my own work. If you use it, verify your results independently.

`rflexdid` is an R port of the Stata [`flexdid`](https://ideas.repec.org/c/boc/bocode/s459517.html) command (Deb, Norton, Wooldridge & Zabel, 2025).
It implements the FLEX estimator, an OLS-based approach for staggered-timing difference-in-differences in panel and repeated cross-section data that allows treatment effects to vary by group, time, and covariates. The procedure estimates treatment effects at the treatment-cell level (by group and time) and then aggregates those estimates into interpretable average treatment-on-the-treated (ATET) summaries (overall, by exposure/event-time, by cohort, calendar time, group, or group × exposure).

[Deb et al. (2026)](https://www.nber.org/papers/w33026) refer to this approach as FLEX: **F**lexible **L**inear model **E**stimated by OLS with Covariates (**X**).

Identifying assumptions: (i) there must be at least one period in
which all units are untreated (no always‑treated units); (ii) conditional
parallel trends; (iii) for repeated
cross‑sections, any compositional changes over time must be unrelated to
treatment timing conditional on covariates. See Deb et al. (2026) for formal
conditions.

For example usage, see the [rflexdid demo vignette](vignettes/rflexdid_demo.md).

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

## Function reference

### `flexdid()`

Estimate the FLEX difference-in-differences regression.

**Usage**

```r
flexdid(
  formula,
  data,
  tx,
  group,
  time,
  specification  = c("lagsonly", "lagsandleads"),
  xnotinteracted = NULL,
  usercohort     = NULL,
  weights        = NULL,
  vcov           = c("cluster", "robust"),
  cluster        = NULL,
  subset         = NULL,
  verbose        = FALSE
)
```

**Arguments**

| Argument | Description |
|----------|-------------|
| `formula` | Formula. Two-sided, e.g. `y ~ x1 + x2`. The left-hand side is the outcome variable; the right-hand side lists covariates to interact with treatment, group, and time indicators. Use `y ~ 1` for no interacted covariates. |
| `data` | Data frame. |
| `tx` | Character. Column name of the binary treatment indicator (`1` = treated this period, `0` = otherwise). |
| `group` | Character. Column name of the group variable used for group fixed effects and, by default, the level at which ATETs are estimated. |
| `time` | Character. Column name of the integer time variable. Must be equally spaced unless `usercohort` is supplied. |
| `specification` | Character. One of `"lagsonly"` (default) or `"lagsandleads"`. Controls whether pre-treatment event-time indicators are included in the regression. |
| `xnotinteracted` | Optional additive controls that are not interacted with treatment cells. Character vector of column names, or a one-sided formula (e.g. `~ medu + factor(schools)`). Must be disjoint from the covariates in `formula`. |
| `usercohort` | Character. Optional column name of a user-supplied cohort variable. Overrides the internal cohort calculation. Useful when time periods are not equally spaced. |
| `weights` | Character. Optional column name of observation weights, treated as probability weights (matching Stata's default for clustered designs). |
| `vcov` | Character. One of `"cluster"` (default) or `"robust"`. |
| `cluster` | Character. Optional column name for the cluster variable. Defaults to `group` when `vcov = "cluster"`. |
| `subset` | Optional logical vector or unevaluated expression restricting the sample (analogous to Stata's `if`). |
| `verbose` | Logical. If `TRUE`, print the regression coefficient table after fitting. |

**Value**

An S3 object of class `flexdid`. Key components:

| Component | Description |
|-----------|-------------|
| `coefficients` | Named numeric vector of OLS estimates (pivoted-out columns set to zero). |
| `vcov` | Full `p × p` variance–covariance matrix for `coefficients`. |
| `vcov_type` | Character; the VCE method used (`"cluster"` or `"robust"`). |
| `n_clusters` | Integer; number of clusters (only when `vcov = "cluster"`). |
| `fitted` | Numeric vector of fitted values (length = number of in-sample observations). |
| `residuals` | Numeric vector of OLS residuals. |
| `df_residual` | Residual degrees of freedom (`n − rank`). |
| `nobs` | Number of observations used in estimation. |
| `r2` / `r2_adj` | R-squared and adjusted R-squared. |
| `f` | Overall F-statistic. |
| `specification` | The specification used (`"lagsonly"` or `"lagsandleads"`). |
| `cohort` | Integer vector of cohort assignments for each in-sample observation (0 for never-treated). |
| `tx_indicator` | Integer vector of the `_Tx` indicator (the treatment-cell indicator, not the raw treatment variable). |
| `eventtime` | Integer vector of event times (`time − cohort` for treated observations, `−1` for controls). |
| `design` | List of design metadata used by `atet()`, including the sparse model matrix and cell-index lookup. |

The object supports `print()`, `summary()`, `coef()`, and `vcov()`. Pass it to `atet()` to compute ATET aggregations.

---

### `atet()`

Aggregate cell-level effects to average treatment effects on the treated.

```r
atet(
  model,
  type = c("overall", "byexposure", "bycalendar",
    "bycohort", "bygroup", "byget"),
  values = NULL,
  for_expr = NULL,
  aggregationweight = c("obslevel", "grouplevel"),
  test = NULL,
  level = 0.95
)
```

| Argument | Description |
|----------|-------------|
| `model` | A `flexdid` object returned by `flexdid()`. |
| `type` | Character. Aggregation type; one of `"overall"`, `"byexposure"`, `"bycalendar"`, `"bycohort"`, `"bygroup"`, or `"byget"`. `"overall"` collapses all post-treatment cells to a single ATET; `"byexposure"` reports one ATET per event time; `"bycalendar"` by calendar period; `"bycohort"` by cohort; `"bygroup"` by group; `"byget"` by group × event time. |
| `values` | Numeric vector restricting which event times (for `"overall"`, `"byexposure"`, `"byget"`) or calendar periods (for `"bycalendar"`) to include. Defaults to all observed post-treatment values (all event times for `"lagsandleads"`). |
| `for_expr` | One-sided formula whose right-hand side is any R expression returning a logical vector; the ATET is averaged only over observations where the expression is `TRUE`. Variables are looked up in the modeling data. Compound conditions use the standard R logical operators `&` (AND), `\|` (OR), and `!` (NOT) — e.g. `~ female == 1`, `~ female == 1 & age > 10`. A `quote()`-style unevaluated expression is also accepted. Equivalent to Stata's `for(...)` option. |
| `aggregationweight` | Character. One of `"obslevel"` (default) or `"grouplevel"`. Group-level weighting reweights cells by group size rather than individual observations, matching Stata's `aggregationweight(grouplevel)`. |
| `test` | Character. Optional Wald F-test to append; one of `"zero"`, `"equal"`, or `"pretrends"`. `"zero"` tests H₀: all ATETs = 0; `"equal"` tests H₀: all ATETs are equal; `"pretrends"` tests H₀: all pre-treatment ATETs = 0 (parallel-trends pre-test, only valid for `type = "byexposure"`). |
| `level` | Numeric in (0, 1). Confidence level for the tidy table (default `0.95`). |

Returns a `flexdid_atet` object with components `estimate`, `vcov`, `tidy_table`, `test_result`, and others. Supports `print()`, `as.data.frame()`, and `plot()` (the last requires `ggplot2`).

## Mapping from Stata to R

| Stata                                                    | R                                                             |
|----------------------------------------------------------|---------------------------------------------------------------|
| `flexdid ssb_oz, tx(treated) group(county) time(year)`  | `flexdid(ssb_oz ~ 1, data, tx="treated", group="county", time="year")` |
| `flexdid ssb_oz age female, ...`                        | `flexdid(ssb_oz ~ age + female, ...)`                         |
| `specification(lagsandleads)`                            | `specification = "lagsandleads"`                              |
| `xnotinteracted(region)`                                 | `xnotinteracted = ~ region`                                   |
| `usercohort(cohort)`                                     | `usercohort = "cohort"`                                       |
| `vce(cluster county)`                                    | `vcov = "cluster", cluster = "county"`                        |
| `[pweight=w]`                                            | `weights = "w"`                                               |
| `estat atet, overall`                                    | `atet(fit, type = "overall")`                                 |
| `estat atet, byexposure(0/3) test(zero)`                | `atet(fit, type = "byexposure", values = 0:3, test = "zero")` |
| `estat atet, byexposure for(female==1)`                 | `atet(fit, type = "byexposure", for_expr = ~ female == 1)`    |
| `estat atet, byget`                                      | `atet(fit, type = "byget")`                                   |
| `aggregationweight(grouplevel)`                          | `aggregationweight = "grouplevel"`                            |

The R object returned by `flexdid()` carries the full design matrix and
sample-aligned vectors needed by `atet()`. `summary()`, `coef()`, `vcov()`,
and `print()` work as expected.

## Standard errors

Standard errors on the regression coefficients use the same HC1 / CR1
sandwich formulas as Stata's `regress, vce(robust)` and
`regress, vce(cluster)`. Standard errors on the aggregated ATETs use the
influence-function form that mirrors Stata's `margins, vce(unconditional)`
— see [R/atet.R](R/atet.R) for the explicit formula.

## Testing

The package has two test suites:

- **Internal tests** ([`test-cohort.R`](tests/testthat/test-cohort.R),
  [`test-flexdid-internal.R`](tests/testthat/test-flexdid-internal.R)) —
  verify numerical correctness on synthetic panels without requiring Stata.
  Checks include: OLS coefficients matching `lm()` to machine precision, ATET
  aggregation identities, influence-function SE formulas, Wald test
  statistics, and `for_expr` subgroup subsetting.

- **Stata-reference tests**
  ([`test-stata-reference.R`](tests/testthat/test-stata-reference.R)) —
  compare point estimates and standard errors against output produced by the
  original Stata `flexdid` command run on the simulated dataset in
  [`inst/extdata/example_data.csv`](inst/extdata/example_data.csv). Coverage spans
  all six ATET aggregation types (`overall`, `byexposure`, `bycalendar`,
  `bycohort`, `bygroup`, `byget`) and both specifications (`lagsonly` and
  `lagsandleads`), as well as the underlying regression coefficients. These
  tests skip automatically if the reference CSVs are not present.

The recommended invocation is the wrapper script, which mirrors progress to
the console and writes a log to [`tests/test_results.txt`](tests/test_results.txt):

```bash
Rscript tests/run_tests.R
```

`devtools::test()` from R also still works.

**Regenerating Stata reference CSVs.** On a machine with Stata:

```stata
do path/to/r-flexdid/tests/testthat/stata_reference/make-stata-reference.do
```

The do-file is working-directory independent — it resolves the repo root from
its own path. If [`inst/extdata/example_data.csv`](inst/extdata/example_data.csv)
is missing, regenerate it first:

```r
devtools::load_all()
simulate_flexdid_data()
```

## References

Deb, P., Norton, E. C., Wooldridge, J. M., & Zabel, J. E. (2026). A Flexible,
Heterogeneous Treatment Effects Difference-in-Differences Estimator for
Repeated Cross-Sections. NBER Working Paper No. 33026. https://www.nber.org/papers/w33026
