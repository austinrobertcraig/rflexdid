* make-stata-reference.do
*
* Run this once on a machine with Stata >= 16 and the `flexdid` command (>= v2.0)
* to produce the reference CSVs that the testthat suite compares against.
* The script resolves the repo root without requiring any machine-specific paths.
* It first tries c(do_file) (set when you type `do file.do` in the Command window),
* then walks upward from c(pwd) looking for the R-package DESCRIPTION file.
* The walk-up approach works from the do-file editor (Do button / Ctrl+D) as long
* as Stata's working directory is anywhere inside the repo.
*
* Inputs:
*   inst/extdata/example_data.csv  (regenerate with: devtools::load_all(); simulate_flexdid_data())
*
* Outputs (overwritten):
*   tests/testthat/stata_reference/coefs_lagsandleads.csv     (g, t, kind, b, se)
*   tests/testthat/stata_reference/coefs_lagsonly.csv         (g, t, kind, b, se)
*   tests/testthat/stata_reference/atet_overall.csv           (b, se)
*   tests/testthat/stata_reference/atet_byexposure.csv        (eventtime, b, se)
*   tests/testthat/stata_reference/atet_bycalendar.csv        (t, b, se)
*   tests/testthat/stata_reference/atet_bycohort.csv          (g, b, se)
*   tests/testthat/stata_reference/atet_bygroup.csv           (g, b, se)
*   tests/testthat/stata_reference/atet_byget.csv             (g, eventtime, b, se)
*   tests/testthat/stata_reference/atet_overall_lagsonly.csv  (b, se)
*   tests/testthat/stata_reference/atet_byexposure_lagsonly.csv (eventtime, b, se)
*
* Usage (no machine-specific setup needed):
*   - Open the file in the do-file editor and click Do / Ctrl+D, OR
*   - Type in the Command window: do path/to/make-stata-reference.do
*
* If the walk-up fails (Stata's working directory is outside the repo), cd to
* the repo root first:
*     cd path/to/r-flexdid

clear all
set more off

* ---------- Resolve repo root ----------
* Strategy 1: c(do_file) — reliable when typed as `do file.do` in the Command window.
* Strategy 2: walk upward from c(pwd) looking for the R-package DESCRIPTION file.
*   This handles the do-file editor (Do button / Ctrl+D), which does not set c(do_file).
local repo_root ""

local do_path `"`c(do_file)'"'
local do_path = subinstr(`"`do_path'"', "\", "/", .)
if `"`do_path'"' != "" {
    local repo_root = subinstr(`"`do_path'"', "/tests/testthat/stata_reference/make-stata-reference.do", "", 1)
}

if `"`repo_root'"' == "" | `"`repo_root'"' == `"`do_path'"' {
    * c(do_file) was empty or the suffix wasn't found; walk up from c(pwd).
    local searchdir `"`c(pwd)'"'
    local searchdir = subinstr(`"`searchdir'"', "\", "/", .)
    forvalues _lvl = 1/20 {
        capture confirm file `"`searchdir'/DESCRIPTION"'
        if !_rc {
            local repo_root `"`searchdir'"'
            continue, break
        }
        local newdir = regexr(`"`searchdir'"', "/[^/]+$", "")
        if `"`newdir'"' == `"`searchdir'"' | `"`newdir'"' == "" {
            continue, break
        }
        local searchdir `"`newdir'"'
    }
}

if `"`repo_root'"' == "" {
    display as error "ERROR: could not locate the repo root."
    display as error "  Stata's working directory must be somewhere inside the r-flexdid repo."
    display as error "  Fix: cd to the repo root, then click Do again."
    exit 198
}

cd `"`repo_root'"'
display as text "Repo root: " `"`repo_root'"'

capture which flexdid
if _rc ssc install flexdid

capture confirm file "inst/extdata/example_data.csv"
if _rc {
    display as error "ERROR: inst/extdata/example_data.csv not found."
    display as error "  Regenerate it from R first: devtools::load_all(); simulate_flexdid_data()"
    exit 601
}

local OUT_REF "tests/testthat/stata_reference"
capture mkdir "tests"
capture mkdir "tests/testthat"
capture mkdir "`OUT_REF'"

* ---------- Load the simulated data ----------
import delimited "inst/extdata/example_data.csv", clear stringcols(_all)
* All columns import as strings because of the factor `region`. Coerce the
* numeric ones back to numeric.
foreach v in county year cohort treated age female ssb_oz {
    capture confirm variable `v'
    if !_rc destring `v', replace
}

* ---------- Helpers ----------
* dump_coefs <fname> <groupvar> <timevar>: write (g, t, kind, b, se) for each
* treatment-cell column of e(b). Stata's factor-variable expansion drops the
* level marker for the lowest level of each component variable, so the colname
* "2013._Cohort#2013.cohort#2013.year#1._Tx" can appear in any of these forms:
*   2013._Cohort#2013.cohort#2013.year#1._Tx   (full)
*    _Cohort#2013.cohort#2013.year#1._Tx       (cohort=base treated cohort)
*   2013._Cohort#2013.cohort#1._Tx             (year=base year)
*    _Cohort#2013.cohort#1._Tx                 (both base)
* Position-based numeric-token counting can't handle the missing tokens, so
* we extract by KEYWORD anchor instead: pull the integer preceding `._Cohort`,
* `.<groupvar>`, and `.<timevar>` (a `bn|b|o` factor modifier may attach to
* any of them and is captured-and-discarded). Missing tokens fall back to the
* corresponding base level pulled from the data.
capture program drop dump_coefs
program define dump_coefs
    args fname groupvar timevar
    matrix B = e(b)
    matrix V = e(V)
    local k = colsof(B)
    local names : colnames B

    * Base levels Stata may omit from colnames.
    quietly summarize _Cohort if _Cohort > 0
    local base_c = r(min)
    quietly summarize `timevar'
    local base_t = r(min)

    preserve
    clear
    set obs `k'
    gen str120 _name = ""
    gen double g    = .
    gen double t    = .
    gen str4   kind = ""
    gen double b    = .
    gen double se   = .
    forvalues i = 1/`k' {
        local nm : word `i' of `names'
        quietly replace _name = "`nm'" in `i'
        quietly replace b     = B[1, `i'] in `i'
        quietly replace se    = sqrt(V[`i', `i']) in `i'
    }
    * Keep indicator-only treatment cells: contain "_Cohort#" and don't have a
    * covariate interaction past _Tx (which would show up as "_Tx#").
    keep if regexm(_name, "_Cohort#") & !regexm(_name, "_Tx#")

    * Pull the leading integer for each of _Cohort, the user's group var, and
    * the time var (whichever appear). regexs(1) returns the first capture
    * group of the most recent regexm() call evaluated for that row.
    gen str20 _str_c = ""
    gen str20 _str_g = ""
    gen str20 _str_t = ""
    quietly replace _str_c = regexs(1) if regexm(_name, "(-?[0-9]+)(bn|b|o)?\._Cohort")
    quietly replace _str_g = regexs(1) if regexm(_name, "(-?[0-9]+)(bn|b|o)?\.`groupvar'")
    quietly replace _str_t = regexs(1) if regexm(_name, "(-?[0-9]+)(bn|b|o)?\.`timevar'")
    gen double _cohort = real(_str_c)
    gen double _gnum   = real(_str_g)
    gen double _ynum   = real(_str_t)
    replace _cohort = `base_c' if missing(_cohort)
    replace _ynum   = `base_t' if missing(_ynum)

    replace g = _gnum
    replace t = _ynum
    replace kind = cond(t >= _cohort, "lag", "lead")
    keep g t kind b se
    order g t kind b se
    sort g t kind
    export delimited using "`fname'", replace
    restore
end

* dump_atet_overall <fname>: estat atet, overall must already have been run.
capture program drop dump_atet_overall
program define dump_atet_overall
    args fname
    matrix B = r(b)
    matrix V = r(V)
    preserve
    clear
    set obs 1
    gen double b  = B[1, 1]
    gen double se = sqrt(V[1, 1])
    export delimited using "`fname'", replace
    restore
end

* dump_atet_1d <fname> <colname>: write (<colname>, b, se), one row per
* nonzero-SE column of r(b). The label/key is the leading integer token of
* the column name (event time, calendar year, group value, etc.).
capture program drop dump_atet_1d
program define dump_atet_1d
    args fname colname
    matrix B = r(b)
    matrix V = r(V)
    local k = colsof(B)
    local names : colnames B
    preserve
    clear
    set obs `k'
    gen str40  _name = ""
    gen double `colname' = .
    gen double b  = .
    gen double se = .
    forvalues i = 1/`k' {
        local nm : word `i' of `names'
        quietly replace _name = "`nm'" in `i'
        quietly replace b  = B[1, `i'] in `i'
        quietly replace se = sqrt(V[`i', `i']) in `i'
    }
    * Pull the first signed-integer token out of each Stata colname.
    gen str20 _key = regexs(0) if regexm(_name, "-?[0-9]+")
    replace `colname' = real(_key)
    drop _name _key
    * Drop omitted/base levels (Stata sets se=0 for those).
    keep if se > 0 & !missing(se) & !missing(`colname')
    order `colname' b se
    sort `colname'
    export delimited using "`fname'", replace
    restore
end

* dump_atet_bycohort <fname> <groupvar>: ordinal Stata cohort labels (1, 2, 3)
* are converted back to actual cohort values using levelsof on `groupvar'.
* Assumes `groupvar' is the cohort-defining variable in the active dataset
* (e.g. for spec 1, groupvar = "cohort"; the cohort levels are the positive
* values of that variable).
capture program drop dump_atet_bycohort
program define dump_atet_bycohort
    args fname groupvar
    matrix B = r(b)
    matrix V = r(V)
    local k = colsof(B)
    local names : colnames B
    * Build a list of actual treated-cohort values (sorted ascending), excluding
    * the never-treated 0.
    quietly levelsof `groupvar' if `groupvar' > 0 & `groupvar' < ., local(cohort_vals)
    preserve
    clear
    set obs `k'
    gen double g  = .
    gen double b  = .
    gen double se = .
    forvalues i = 1/`k' {
        local nm : word `i' of `names'
        * Stata's bycohort label is an ordinal integer: 1, 2, 3, ...
        if regexm("`nm'", "-?[0-9]+") {
            local ord = real(regexs(0))
        }
        else {
            local ord ""
        }
        if "`ord'" != "" & `ord' >= 1 {
            local val : word `ord' of `cohort_vals'
            if "`val'" != "" {
                quietly replace g = real("`val'") in `i'
            }
        }
        quietly replace b  = B[1, `i'] in `i'
        quietly replace se = sqrt(V[`i', `i']) in `i'
    }
    keep if se > 0 & !missing(se) & !missing(g)
    order g b se
    sort g
    export delimited using "`fname'", replace
    restore
end

* dump_atet_byget <fname> <groupvar> <timevar>: write (g, eventtime, b, se).
*
* Why we don't parse colnames: flexdid_atet.ado tries to assign
*   matrix colnames `beta' = `lofgt'
* where lofgt has compound-quoted "<g> <eventtime>" pairs. Stata's
* `matrix colnames` requires names without spaces, so this assignment is
* unreliable — the override may silently fail and r(b) keeps the
* margins-generated "<gt-level>.<tempvar>" names, which only carry one
* numeric token per column.
*
* Instead we replicate flexdid_atet's enumeration order in the do-file using
* the same `levelsof` walk on `groupvar` and event time. egen group() and
* margins, over() both visit levels in the same sorted order, so column i of
* r(b) corresponds to the i-th (g, eventtime) pair in our enumeration.
capture program drop dump_atet_byget
program define dump_atet_byget
    args fname groupvar timevar
    matrix B = r(b)
    matrix V = r(V)
    local k = colsof(B)

    * Mirror flexdid_atet.ado's lofgt construction (see byget block):
    *   levelsof `group' if _Cohort > 0   then per group, levelsof of event time
    *   among _Tx==1 cells. Sorted ascending, group-major, event-time-minor.
    tempvar et
    quietly gen double `et' = `timevar' - _Cohort if _Cohort > 0
    quietly levelsof `groupvar' if _Cohort > 0, local(lofg)
    local i = 0
    foreach gval of local lofg {
        quietly levelsof `et' if `groupvar'==`gval' & _Tx==1, local(loft)
        foreach tval of local loft {
            local ++i
            local g_`i' = `gval'
            local t_`i' = `tval'
        }
    }
    local n_pairs = `i'
    if `n_pairs' != `k' {
        display as error ///
          "dump_atet_byget: r(b) has `k' columns but enumerated `n_pairs' (g, eventtime) pairs."
        display as error ///
          "  Output rows may not correspond to the right cells. Check that flexdid_atet's"
        display as error ///
          "  byget enumeration matches: levelsof `groupvar' if _Cohort>0, then per-group"
        display as error ///
          "  levelsof of (timevar - _Cohort) if _Tx==1."
    }

    preserve
    clear
    set obs `k'
    gen double g         = .
    gen double eventtime = .
    gen double b  = .
    gen double se = .
    forvalues j = 1/`k' {
        if `j' <= `n_pairs' {
            quietly replace g         = `g_`j'' in `j'
            quietly replace eventtime = `t_`j'' in `j'
        }
        quietly replace b  = B[1, `j'] in `j'
        quietly replace se = sqrt(V[`j', `j']) in `j'
    }
    keep if se > 0 & !missing(se) & !missing(g) & !missing(eventtime)
    keep g eventtime b se
    order g eventtime b se
    sort g eventtime
    export delimited using "`fname'", replace
    restore
end

* ============================================================================
* Spec 1: lagsandleads, group=cohort, vce(cluster county)
* ============================================================================
flexdid ssb_oz, tx(treated) group(cohort) time(year) ///
        specification(lagsandleads) vce(cluster county)

dump_coefs "`OUT_REF'/coefs_lagsandleads.csv" "cohort" "year"

estat atet, overall
dump_atet_overall "`OUT_REF'/atet_overall.csv"

estat atet, byexposure nograph
dump_atet_1d "`OUT_REF'/atet_byexposure.csv" "eventtime"

estat atet, bycalendar nograph
dump_atet_1d "`OUT_REF'/atet_bycalendar.csv" "t"

estat atet, bycohort nograph
dump_atet_bycohort "`OUT_REF'/atet_bycohort.csv" "cohort"

estat atet, bygroup nograph
dump_atet_1d "`OUT_REF'/atet_bygroup.csv" "g"

estat atet, byget
dump_atet_byget "`OUT_REF'/atet_byget.csv" "cohort" "year"

* ============================================================================
* Spec 2: lagsonly with covariates (age, female), group=county, vce(cluster county)
* ============================================================================
flexdid ssb_oz age female, tx(treated) group(county) time(year) ///
        specification(lagsonly) vce(cluster county)

dump_coefs "`OUT_REF'/coefs_lagsonly.csv" "county" "year"

estat atet, overall
dump_atet_overall "`OUT_REF'/atet_overall_lagsonly.csv"

estat atet, byexposure nograph
dump_atet_1d "`OUT_REF'/atet_byexposure_lagsonly.csv" "eventtime"

display as text "Done. Wrote reference CSVs to `OUT_REF'/."
