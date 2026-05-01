* make-stata-reference.do
*
* Run this once on a machine with Stata >= 16 and the `flexdid` command (>= v2.0)
* to produce the reference CSVs that the testthat suite compares against.
* The script is robust to the working directory: it resolves the repo root
* from `c(do_file)` and `cd`s there before reading or writing anything.
*
* Inputs:
*   sim_data/example_data.csv  (regenerate with: Rscript sim_data/simulate_data.R)
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
* Usage (from anywhere inside Stata):
*     do path/to/r-flexdid/tests/testthat/stata_reference/make-stata-reference.do

clear all
set more off

* ---------- Resolve repo root from this do-file's path ----------
local do_path = `"`c(do_file)'"'
if `"`do_path'"' == "" {
    display as error "ERROR: c(do_file) is empty (Stata < 16). Upgrade Stata or run from the repo root and adjust paths manually."
    exit 198
}
* Normalize Windows backslashes to forward slashes for portable string ops.
local do_path = subinstr(`"`do_path'"', "\", "/", .)
* Strip "/tests/testthat/stata_reference/make-stata-reference.do" off the end
* to get the repo root.
local repo_root = subinstr(`"`do_path'"', "/tests/testthat/stata_reference/make-stata-reference.do", "", 1)
cd `"`repo_root'"'
display as text "Repo root: " `"`repo_root'"'

capture which flexdid
if _rc ssc install flexdid

capture confirm file "sim_data/example_data.csv"
if _rc {
    display as error "ERROR: sim_data/example_data.csv not found."
    display as error "  Regenerate it from R first: Rscript sim_data/simulate_data.R"
    exit 601
}

local OUT_REF "tests/testthat/stata_reference"
capture mkdir "tests"
capture mkdir "tests/testthat"
capture mkdir "`OUT_REF'"

* ---------- Load the simulated data ----------
import delimited "sim_data/example_data.csv", clear stringcols(_all)
* All columns import as strings because of the factor `region`. Coerce the
* numeric ones back to numeric.
foreach v in county year cohort treated age female ssb_oz {
    capture confirm variable `v'
    if !_rc destring `v', replace
}

* ---------- Helpers ----------
* dump_coefs <fname>: write (g, t, kind, b, se) for each treatment-cell column.
* Treatment cells are the indicator-only columns of e(b) whose names match
* "<cohort>...<groupvar>#<year>...._Tx" (no further "#<xvar>" suffix).
capture program drop dump_coefs
program define dump_coefs
    args fname groupvar
    matrix B = e(b)
    matrix V = e(V)
    local k = colsof(B)
    local names : colnames B
    preserve
    clear
    set obs `k'
    gen str40 _name = ""
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
    * Keep indicator-only treatment cells:
    *   - contain "_Cohort#"
    *   - end with "._Tx" (no covariate-interaction "#" after _Tx)
    keep if regexm(_name, "_Cohort#") & !regexm(_name, "_Tx#")
    * Extract the three leading numeric tokens: cohort, group, year.
    gen str20 _nums = _name
    quietly replace _nums = subinstr(_nums, ".", " ", .)
    quietly replace _nums = subinstr(_nums, "#", " ", .)
    * Drop Stata factor modifiers (bn, b, o) attached to numbers.
    quietly replace _nums = regexr(_nums, "([0-9])(bn|b|o)", "$1")
    quietly replace _nums = regexr(_nums, "([0-9])(bn|b|o)", "$1")
    quietly replace _nums = regexr(_nums, "([0-9])(bn|b|o)", "$1")
    * Now _nums looks like: "<cohort>  _Cohort <group> `groupvar' <year> year 1 _Tx"
    * Replace any non-numeric tokens with spaces, then split.
    gen double _cohort = real(word(_nums, 1))
    * The group-value token is the third word once non-numeric words are skipped;
    * easier: enumerate words and pull the first 3 numeric ones.
    quietly {
        gen double _gnum = .
        gen double _ynum = .
        local maxw = 12
        forvalues r = 1/`=_N' {
            local seenc = 0
            local seeng = 0
            local seent = 0
            forvalues w = 1/`maxw' {
                local tok = word(_nums[`r'], `w')
                if "`tok'" == "" continue
                capture confirm number `tok'
                if !_rc {
                    if `seenc' == 0 {
                        local seenc = 1
                    }
                    else if `seeng' == 0 {
                        local seeng = 1
                        replace _gnum = real("`tok'") in `r'
                    }
                    else if `seent' == 0 {
                        local seent = 1
                        replace _ynum = real("`tok'") in `r'
                    }
                }
            }
        }
    }
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
        local ord = real(regexs(0)) if regexm("`nm'", "-?[0-9]+")
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

* dump_atet_byget <fname>: write (g, eventtime, b, se). Each r(b) colname for
* byget encodes a (group, eventtime) pair in factor-interaction form, e.g.
* "2013.cohort#0.eventtime" or with bn/b/o modifiers attached. Pull the two
* leading signed-integer tokens.
capture program drop dump_atet_byget
program define dump_atet_byget
    args fname
    matrix B = r(b)
    matrix V = r(V)
    local k = colsof(B)
    local names : colnames B
    preserve
    clear
    set obs `k'
    gen str60  _name = ""
    gen double g         = .
    gen double eventtime = .
    gen double b  = .
    gen double se = .
    forvalues i = 1/`k' {
        local nm : word `i' of `names'
        quietly replace _name = "`nm'" in `i'
        quietly replace b  = B[1, `i'] in `i'
        quietly replace se = sqrt(V[`i', `i']) in `i'
    }
    * Strip Stata factor modifiers, then collect the first two signed integers.
    gen str60 _scrub = _name
    quietly replace _scrub = regexr(_scrub, "([0-9])(bn|b|o)", "$1")
    quietly replace _scrub = regexr(_scrub, "([0-9])(bn|b|o)", "$1")
    quietly replace _scrub = subinstr(_scrub, ".", " ", .)
    quietly replace _scrub = subinstr(_scrub, "#", " ", .)
    quietly {
        forvalues r = 1/`=_N' {
            local seen = 0
            local maxw = 12
            forvalues w = 1/`maxw' {
                local tok = word(_scrub[`r'], `w')
                if "`tok'" == "" continue
                capture confirm number `tok'
                if !_rc {
                    if `seen' == 0 {
                        replace g = real("`tok'") in `r'
                        local seen = 1
                    }
                    else {
                        replace eventtime = real("`tok'") in `r'
                        continue, break
                    }
                }
            }
        }
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

dump_coefs "`OUT_REF'/coefs_lagsandleads.csv" "cohort"

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
dump_atet_byget "`OUT_REF'/atet_byget.csv"

* ============================================================================
* Spec 2: lagsonly with covariates (age, female), group=county, vce(cluster county)
* ============================================================================
flexdid ssb_oz age female, tx(treated) group(county) time(year) ///
        specification(lagsonly) vce(cluster county)

dump_coefs "`OUT_REF'/coefs_lagsonly.csv" "county"

estat atet, overall
dump_atet_overall "`OUT_REF'/atet_overall_lagsonly.csv"

estat atet, byexposure nograph
dump_atet_1d "`OUT_REF'/atet_byexposure_lagsonly.csv" "eventtime"

display as text "Done. Wrote reference CSVs to `OUT_REF'/."
