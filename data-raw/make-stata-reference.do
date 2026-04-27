* make-stata-reference.do
*
* Run this once in Stata (>= 17) to produce the reference CSVs that the
* testthat suite compares against. After it finishes, the only files that
* should change in the repo are:
*
*   inst/extdata/hhabits.csv
*   tests/testthat/stata_reference/coefs_lagsonly.csv
*   tests/testthat/stata_reference/coefs_lagsandleads.csv
*   tests/testthat/stata_reference/atet_overall.csv
*   tests/testthat/stata_reference/atet_byexposure.csv
*   tests/testthat/stata_reference/atet_bycalendar.csv
*   tests/testthat/stata_reference/atet_bycohort.csv
*   tests/testthat/stata_reference/atet_bygroup.csv
*   tests/testthat/stata_reference/atet_byget.csv
*
* Usage from the repo root:
*   stata-mp -b do data-raw/make-stata-reference.do
*
* If the user has an older flexdid installed, it pulls the latest from SSC.

clear all
set more off
capture which flexdid
if _rc ssc install flexdid

* Repo root is the working directory for this run.
local OUT_DATA   "inst/extdata/hhabits.csv"
local OUT_REF    "tests/testthat/stata_reference"
capture mkdir "tests"
capture mkdir "tests/testthat"
capture mkdir "tests/testthat/stata_reference"

* ---------- Data ----------
webuse hhabits, clear

* Build the manual cohort variable used in the help-file examples.
egen chrt = min(year/hhabit), by(schools)
replace chrt = 0 if chrt==.

export delimited using "`OUT_DATA'", replace

* ---------- Spec 1: lagsandleads, group=chrt, cluster=schools ----------
* Matches the spec used throughout the help-file examples.
flexdid bmi, tx(hhabit) group(chrt) time(year) ///
        specification(lagsandleads) vce(cluster schools)

* Save the regression coefficient table.
matrix B = e(b)
matrix V = e(V)
local nb = colsof(B)
* Build a small dataset with name, b, se for every coefficient.
preserve
clear
set obs `nb'
gen str30 name = ""
gen double b   = .
gen double se  = .
forvalues i = 1/`nb' {
    local nm : word `i' of `: colnames B'
    quietly replace name = "`nm'" in `i'
    quietly replace b    = B[1, `i'] in `i'
    quietly replace se   = sqrt(V[`i', `i']) in `i'
}
export delimited using "`OUT_REF'/coefs_lagsandleads.csv", replace
restore

* ---------- ATETs from spec 1 ----------
* A small helper: dump r(b), r(V) to CSV with one row per ATET level.
capture program drop dump_atet
program define dump_atet
    args fname
    matrix B  = r(b)
    matrix V  = r(V)
    local k = colsof(B)
    preserve
    clear
    set obs `k'
    gen str40 label = ""
    gen double b    = .
    gen double se   = .
    forvalues i = 1/`k' {
        local nm : word `i' of `: colnames B'
        quietly replace label = "`nm'" in `i'
        quietly replace b     = B[1, `i'] in `i'
        quietly replace se    = sqrt(V[`i', `i']) in `i'
    }
    export delimited using "`fname'", replace
    restore
end

estat atet, overall
dump_atet "`OUT_REF'/atet_overall.csv"

estat atet, byexposure nograph
dump_atet "`OUT_REF'/atet_byexposure.csv"

estat atet, bycalendar nograph
dump_atet "`OUT_REF'/atet_bycalendar.csv"

estat atet, bycohort nograph
dump_atet "`OUT_REF'/atet_bycohort.csv"

estat atet, bygroup nograph
dump_atet "`OUT_REF'/atet_bygroup.csv"

estat atet, byget
dump_atet "`OUT_REF'/atet_byget.csv"

* ---------- Spec 2: lagsonly with covariates, cluster=schools ----------
* Matches the canonical "with covariates" example.
flexdid bmi girl medu, tx(hhabit) group(schools) time(year)

matrix B = e(b)
matrix V = e(V)
local nb = colsof(B)
preserve
clear
set obs `nb'
gen str30 name = ""
gen double b   = .
gen double se  = .
forvalues i = 1/`nb' {
    local nm : word `i' of `: colnames B'
    quietly replace name = "`nm'" in `i'
    quietly replace b    = B[1, `i'] in `i'
    quietly replace se   = sqrt(V[`i', `i']) in `i'
}
export delimited using "`OUT_REF'/coefs_lagsonly.csv", replace
restore

* ---------- ATETs from spec 2 ----------
estat atet, overall
dump_atet "`OUT_REF'/atet_overall_lagsonly.csv"

estat atet, byexposure nograph
dump_atet "`OUT_REF'/atet_byexposure_lagsonly.csv"

display "Done. Wrote reference CSVs to `OUT_REF'/."
