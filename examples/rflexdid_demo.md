<!--
Demo of rflexdid on simulated data. To regenerate rflexdid_demo.md:

    cd examples && Rscript -e 'knitr::knit("rflexdid_demo.Rmd")'

That produces rflexdid_demo.md plus a figure/ subdirectory of plot PNGs.
GitHub renders the .md natively (code, output, and embedded plots inline).
-->



# Demonstrating rflexdid on simulated data

## The simulated dataset

The script `simulate_data.R` builds a repeated cross section of 20,000
adults surveyed across 40 counties over 10 years (2010–2019). Each year
brings a fresh sample of 50 individuals per county; no individual is
followed across time. Counties enact a sugar-sweetened-beverage (SSB)
excise tax in different years (the staggered cohorts: 2013, 2015, 2017,
plus a never-treated cohort labeled `0`). The outcome is daily SSB
consumption in fluid ounces.

The data-generating process imposes strict parallel pre-trends and known
treatment effects: an immediate jump `alpha_g` (-1.5 oz/day for the 2013
cohort, -1.0 for 2015, -0.5 for 2017) plus a smooth dynamic decline
`-0.20 e - 0.02 e^2` in event time `e`. So early-cohort counties
enacted larger taxes, and consumption keeps falling with longer exposure.


``` r
library(rflexdid)
library(ggplot2)
library(here)
source(here("examples", "simulate_data.R"))
df <- simulate_flexdid_data()
head(df)
```

```
##   county year cohort treated      age female    region    ssb_oz
## 1      1 2010      0       0 25.60111      0      West 10.234163
## 2      2 2010      0       0 18.00000      1 Northeast  8.330653
## 3      3 2010   2013       0 38.73399      1 Northeast 14.632594
## 4      4 2010   2017       0 41.52287      1 Northeast  4.069476
## 5      5 2010   2017       0 21.34340      1     South 12.341625
## 6      6 2010      0       0 56.52613      0      West 16.259764
```

``` r
table(df$cohort)
```

```
## 
##    0 2013 2015 2017 
## 5000 5000 5000 5000
```

## Estimation

`flexdid()` runs a single OLS regression of `ssb_oz` on a saturated set of
cohort × group × time × treatment interactions, with `age` interacted in
on the right-hand side. County and year fixed effects are absorbed
automatically. `specification = "lagsandleads"` adds pre-treatment lead
indicators so the parallel-trends assumption is visually testable.
Standard errors cluster at the county level.


``` r
fit <- flexdid(ssb_oz ~ age,
               data = df,
               tx = "treated", group = "county", time = "year",
               specification = "lagsandleads",
               vcov = "cluster", cluster = "county")
fit
```

```
## Flexible difference-in-differences regression
##   Specification: lagsandleads 
##   Outcome:       ssb_oz 
##   Treatment:     treated  | Group: county  | Time: year 
##   Observations:  20000 
##   VCE:           cluster (county; 40 clusters)
##   R-squared:     0.1827   (adj. 0.1558)
```

## Overall ATET

A single number that averages the cell-level effects across all
post-treatment observations. The DGP's true mean treatment effect over
the treated subset is roughly -1.76 oz/day; the recovered estimate
should be close.


``` r
atet(fit, type = "overall")
```

```
## Overall ATET
## Observations: 20000 | Subpopulation observations: 7500
## Aggregation weight: obslevel
## 
##         Estimate Std. Error t value Pr(>|t|) [CI lo] [CI hi]
## Overall -1.5979   0.1863    -8.5779  0.0000  -1.9630 -1.2328
## 
## Note: Linearization (Stata's vce(unconditional)) is used for the standard errors.
```

## Event study (`byexposure`)

The most useful diagnostic. Pre-period estimates (event time `e < 0`)
should be statistically indistinguishable from zero; post-period
estimates should grow more negative with exposure time, tracing out the
dynamic profile of the tax's effect. The base period `e = -1` is
mechanically zero with SE = 0.


``` r
ax <- atet(fit, type = "byexposure", test = "zero")
ax
```

```
## ATET by exposure time
## Observations: 20000 | Subpopulation observations: 15000
## Aggregation weight: obslevel
## 
##    Estimate Std. Error t value  Pr(>|t|) [CI lo]  [CI hi] 
## -7  0.3339   0.2560     1.3040   0.1923  -0.1680   0.8357 
## -6  0.3460   0.3603     0.9603   0.3369  -0.3602   1.0522 
## -5 -0.0314   0.1660    -0.1892   0.8500  -0.3567   0.2939 
## -4  0.2874   0.2220     1.2944   0.1955  -0.1478   0.7226 
## -3  0.0592   0.1520     0.3896   0.6968  -0.2387   0.3571 
## -2  0.2756   0.1805     1.5265   0.1269  -0.0783   0.6294 
## -1  0.0000   0.0000          NA       NA       NA       NA
## 0  -0.9029   0.2057    -4.3898   0.0000  -1.3061  -0.4998 
## 1  -1.1250   0.1696    -6.6334   0.0000  -1.4574  -0.7926 
## 2  -1.2817   0.1771    -7.2355   0.0000  -1.6289  -0.9345 
## 3  -1.8229   0.2190    -8.3239   0.0000  -2.2521  -1.3936 
## 4  -2.2944   0.2416    -9.4984   0.0000  -2.7679  -1.8209 
## 5  -2.8562   0.3142    -9.0912   0.0000  -3.4720  -2.2404 
## 6  -2.9485   0.2777    -10.6166  0.0000  -3.4928  -2.4041 
## 
## Test of zero ATETs
##   H0: All effects are equal to zero
##   F(13, 19362) = 26.690   Prob > F = 0.0000
## 
## Note: Linearization (Stata's vce(unconditional)) is used for the standard errors.
```

``` r
plot(ax)
```

![plot of chunk byexposure](figure/byexposure-1.png)

The Wald test rejects the joint hypothesis that all event-time effects
are zero.

## By treated cohort

Each treated cohort gets one ATET, averaging that cohort's
post-treatment cell effects. The estimates should rank
2013 (most negative) > 2015 > 2017 (least negative), mirroring the tax
magnitudes baked into the DGP.


``` r
acoh <- atet(fit, type = "bycohort")
acoh
```

```
## ATET by treated cohort
## Observations: 20000 | Subpopulation observations: 7500
## Aggregation weight: obslevel
## 
##      Estimate Std. Error t value  Pr(>|t|) [CI lo] [CI hi]
## 2013 -2.2065   0.2095    -10.5326  0.0000  -2.6171 -1.7959
## 2015 -1.4897   0.2313    -6.4405   0.0000  -1.9431 -1.0363
## 2017 -0.3581   0.2459    -1.4564   0.1453  -0.8400  0.1238
## 
## Note: Linearization (Stata's vce(unconditional)) is used for the standard errors.
```

``` r
plot(acoh)
```

![plot of chunk bycohort](figure/bycohort-1.png)

## By calendar time

A separate ATET for each calendar year, averaging across whichever
cohorts were already post-treatment in that year. Effects start near
zero in 2010–2012 (no cohort treated yet), turn negative in 2013
(2013 cohort enters), and grow more negative as the 2015 and 2017
cohorts roll in.


``` r
acal <- atet(fit, type = "bycalendar")
acal
```

```
## ATET by calendar time
## Observations: 20000 | Subpopulation observations: 13500
## Aggregation weight: obslevel
## 
##      Estimate Std. Error t value Pr(>|t|) [CI lo] [CI hi]
## 2010  0.2167   0.1707     1.2696  0.2042  -0.1179  0.5513
## 2011  0.4084   0.2318     1.7621  0.0781  -0.0459  0.8626
## 2012 -0.2471   0.1896    -1.3029  0.1926  -0.6187  0.1246
## 2013 -0.4458   0.2794    -1.5954  0.1106  -0.9934  0.1019
## 2014 -0.6607   0.3007    -2.1974  0.0280  -1.2500 -0.0714
## 2015 -0.7429   0.2290    -3.2449  0.0012  -1.1917 -0.2942
## 2016 -1.7976   0.1991    -9.0294  0.0000  -2.1879 -1.4074
## 2017 -1.5574   0.2725    -5.7156  0.0000  -2.0914 -1.0233
## 2018 -1.5559   0.2803    -5.5509  0.0000  -2.1053 -1.0065
## 2019 -1.7769   0.2737    -6.4915  0.0000  -2.3134 -1.2404
## 
## Note: Linearization (Stata's vce(unconditional)) is used for the standard errors.
```

``` r
plot(acal)
```

![plot of chunk bycalendar](figure/bycalendar-1.png)

## By treated group (county)

A separate ATET for each treated county. Useful for spotting outlier
counties or for inspecting the within-cohort dispersion.


``` r
agrp <- atet(fit, type = "bygroup")
plot(agrp)
```

![plot of chunk bygroup](figure/bygroup-1.png)

## Group × event time (`byget`)

The most disaggregated view: one ATET per (county, event time) pair.
Restricting `values = 0:3` keeps the table readable.


``` r
aget <- atet(fit, type = "byget", values = 0:3)
head(as.data.frame(aget), 12)
```

```
##     label    Estimate Std. Error     t value     Pr(>|t|)     [CI lo]
## 1  g3|et0 -1.19173887  0.1960159  -6.0798062 1.225823e-09 -1.57594706
## 2  g3|et1 -1.77779029  0.1772682 -10.0288150 1.300401e-23 -2.12525136
## 3  g3|et2 -1.48958725  0.1384743 -10.7571350 6.534993e-27 -1.76100895
## 4  g3|et3 -2.65807660  0.1634677 -16.2605653 4.630006e-59 -2.97848736
## 5  g4|et0 -0.08695303  0.1753483  -0.4958876 6.199794e-01 -0.43065082
## 6  g4|et1 -0.46742651  0.1400854  -3.3367258 8.493108e-04 -0.74200599
## 7  g4|et2  0.28527450  0.2148237   1.3279472 1.842112e-01 -0.13579848
## 8  g5|et0  1.21994585  0.1801225   6.7728671 1.298788e-11  0.86689014
## 9  g5|et1  0.90231119  0.1480671   6.0939336 1.122524e-09  0.61208683
## 10 g5|et2  0.53800196  0.2249698   2.3914406 1.679188e-02  0.09704166
## 11 g7|et0 -1.03105848  0.1843624  -5.5925633 2.267631e-08 -1.39242477
## 12 g7|et1 -1.39257831  0.1941193  -7.1738286 7.556084e-13 -1.77306885
##       [CI hi]
## 1  -0.8075307
## 2  -1.4303292
## 3  -1.2181655
## 4  -2.3376658
## 5   0.2567448
## 6  -0.1928470
## 7   0.7063475
## 8   1.5730016
## 9   1.1925355
## 10  0.9789623
## 11 -0.6696922
## 12 -1.0120878
```

## Subgroup ATET via `for_expr`

Restrict the population over which the ATET is averaged. This is the
analogue of Stata's `for(female == 1)` option. The subgroup estimate
should be close to the overall estimate, since the DGP doesn't make
treatment effects depend on sex.


``` r
atet(fit, type = "byexposure", for_expr = quote(female == 1))
```

```
## ATET by exposure time for `female == 1`
## Observations: 20000 | Subpopulation observations: 7454
## Aggregation weight: obslevel
## 
##    Estimate Std. Error t value  Pr(>|t|) [CI lo]  [CI hi] 
## -7  0.2739   0.2528     1.0833   0.2787  -0.2217   0.7694 
## -6  0.3288   0.3554     0.9251   0.3549  -0.3678   1.0253 
## -5 -0.0683   0.1637    -0.4174   0.6764  -0.3892   0.2525 
## -4  0.2435   0.2208     1.1024   0.2703  -0.1894   0.6763 
## -3  0.0478   0.1530     0.3123   0.7548  -0.2520   0.3476 
## -2  0.2281   0.1886     1.2093   0.2265  -0.1416   0.5977 
## -1  0.0000   0.0000          NA       NA       NA       NA
## 0  -0.9111   0.2068    -4.4060   0.0000  -1.3164  -0.5058 
## 1  -1.1871   0.1735    -6.8433   0.0000  -1.5271  -0.8471 
## 2  -1.2850   0.1822    -7.0508   0.0000  -1.6422  -0.9278 
## 3  -1.8696   0.2116    -8.8354   0.0000  -2.2844  -1.4548 
## 4  -2.3580   0.2481    -9.5044   0.0000  -2.8443  -1.8717 
## 5  -2.8804   0.3211    -8.9698   0.0000  -3.5098  -2.2510 
## 6  -3.0591   0.2929    -10.4444  0.0000  -3.6332  -2.4850 
## 
## Note: Linearization (Stata's vce(unconditional)) is used for the standard errors.
```

## Non-interacted controls (`xnotinteracted`)

`xnotinteracted` adds additive controls that are not interacted with the
treatment-cell indicators — useful for absorbing variance from
group-level covariates without inflating the parameter count. Region is
constant within county, so it slots in cleanly here.


``` r
fit2 <- flexdid(ssb_oz ~ age,
                data = df,
                tx = "treated", group = "county", time = "year",
                specification = "lagsandleads",
                xnotinteracted = ~ region,
                vcov = "cluster", cluster = "county")
atet(fit2, type = "overall")
```

```
## Overall ATET
## Observations: 20000 | Subpopulation observations: 7500
## Aggregation weight: obslevel
## 
##         Estimate Std. Error t value Pr(>|t|) [CI lo] [CI hi]
## Overall -1.5979   0.1863    -8.5779  0.0000  -1.9630 -1.2328
## 
## Note: Linearization (Stata's vce(unconditional)) is used for the standard errors.
```

## Group-level aggregation weights

`aggregationweight = "grouplevel"` reweights cells so that each county
contributes equally to the ATET, regardless of cell size. In this
simulation every county has exactly 500 observations, so the result
matches the default `obslevel` weighting; with unbalanced cell sizes
the two would differ.


``` r
atet(fit, type = "overall", aggregationweight = "grouplevel")
```

```
## Overall ATET
## Observations: 20000 | Subpopulation observations: 7500
## Aggregation weight: grouplevel
## 
##         Estimate Std. Error t value Pr(>|t|) [CI lo] [CI hi]
## Overall -1.5979   0.1863    -8.5779  0.0000  -1.9630 -1.2328
## 
## Note: Linearization (Stata's vce(unconditional)) is used for the standard errors.
```

## Wald test of homogeneity

`test = "equal"` tests the null hypothesis that all event-time effects
are equal — i.e., that the dynamic profile is flat. The DGP makes
effects strictly increase in exposure, so the test should reject.


``` r
atet(fit, type = "byexposure", test = "equal")
```

```
## ATET by exposure time
## Observations: 20000 | Subpopulation observations: 15000
## Aggregation weight: obslevel
## 
##    Estimate Std. Error t value  Pr(>|t|) [CI lo]  [CI hi] 
## -7  0.3339   0.2560     1.3040   0.1923  -0.1680   0.8357 
## -6  0.3460   0.3603     0.9603   0.3369  -0.3602   1.0522 
## -5 -0.0314   0.1660    -0.1892   0.8500  -0.3567   0.2939 
## -4  0.2874   0.2220     1.2944   0.1955  -0.1478   0.7226 
## -3  0.0592   0.1520     0.3896   0.6968  -0.2387   0.3571 
## -2  0.2756   0.1805     1.5265   0.1269  -0.0783   0.6294 
## -1  0.0000   0.0000          NA       NA       NA       NA
## 0  -0.9029   0.2057    -4.3898   0.0000  -1.3061  -0.4998 
## 1  -1.1250   0.1696    -6.6334   0.0000  -1.4574  -0.7926 
## 2  -1.2817   0.1771    -7.2355   0.0000  -1.6289  -0.9345 
## 3  -1.8229   0.2190    -8.3239   0.0000  -2.2521  -1.3936 
## 4  -2.2944   0.2416    -9.4984   0.0000  -2.7679  -1.8209 
## 5  -2.8562   0.3142    -9.0912   0.0000  -3.4720  -2.2404 
## 6  -2.9485   0.2777    -10.6166  0.0000  -3.4928  -2.4041 
## 
## Test of equal ATETs
##   H0: Effects are equal to each other
##   F(12, 19362) = 25.138   Prob > F = 0.0000
## 
## Note: Linearization (Stata's vce(unconditional)) is used for the standard errors.
```
