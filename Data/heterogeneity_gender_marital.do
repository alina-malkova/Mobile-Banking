* ==============================================================================
* Gender and Marital Status Heterogeneity Analysis
* Tests whether mobile banking → self-employment differs by gender/marital status
* ==============================================================================

clear all
set more off

use "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/_Research/Mobile_Money_Banking/Mobile banking USA/Data/analysis_dataset_with_se.dta", clear

* Sample restrictions
keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if banking_mode != .

gen se = (self_employed == 1) if self_employed != .
gen age2 = age * age

di "=============================================="
di "HETEROGENEITY ANALYSIS: GENDER & MARITAL STATUS"
di "=============================================="

* --------------------------------------------------------------------------
* 1. SPLIT-SAMPLE BY GENDER
* --------------------------------------------------------------------------
di _n "=== SPLIT SAMPLE BY GENDER ==="

* Female subsample
di _n "--- Female subsample ---"
reghdfe se mobile_user ///
    age age2 ///
    i.praceeth3 i.peducgrp i.hhincome ///
    married has_children ///
    [pw=hsupwgtk] if female == 1, absorb(cbsa year) cluster(cbsa)
eststo fem
local b_fem = _b[mobile_user]
local se_fem = _se[mobile_user]
local n_fem = e(N)

* Male subsample
di _n "--- Male subsample ---"
reghdfe se mobile_user ///
    age age2 ///
    i.praceeth3 i.peducgrp i.hhincome ///
    married has_children ///
    [pw=hsupwgtk] if female == 0, absorb(cbsa year) cluster(cbsa)
eststo male
local b_male = _b[mobile_user]
local se_male = _se[mobile_user]
local n_male = e(N)

di _n "Gender heterogeneity:"
di "  Female: b = " %7.4f `b_fem' " (SE = " %7.4f `se_fem' "), N = " `n_fem'
di "  Male:   b = " %7.4f `b_male' " (SE = " %7.4f `se_male' "), N = " `n_male'

* Difference test
local diff = `b_fem' - `b_male'
local se_diff = sqrt(`se_fem'^2 + `se_male'^2)
local z_diff = `diff' / `se_diff'
di "  Difference (F-M): " %7.4f `diff' " (SE = " %7.4f `se_diff' ")"
di "  Z-statistic: " %7.3f `z_diff'

* --------------------------------------------------------------------------
* 2. SPLIT-SAMPLE BY MARITAL STATUS
* --------------------------------------------------------------------------
di _n "=== SPLIT SAMPLE BY MARITAL STATUS ==="

* Married subsample
di _n "--- Married subsample ---"
reghdfe se mobile_user ///
    age age2 ///
    i.praceeth3 i.peducgrp i.hhincome ///
    female has_children ///
    [pw=hsupwgtk] if married == 1, absorb(cbsa year) cluster(cbsa)
eststo marr
local b_marr = _b[mobile_user]
local se_marr = _se[mobile_user]
local n_marr = e(N)

* Unmarried subsample
di _n "--- Unmarried subsample ---"
reghdfe se mobile_user ///
    age age2 ///
    i.praceeth3 i.peducgrp i.hhincome ///
    female has_children ///
    [pw=hsupwgtk] if married == 0, absorb(cbsa year) cluster(cbsa)
eststo unmarr
local b_unmarr = _b[mobile_user]
local se_unmarr = _se[mobile_user]
local n_unmarr = e(N)

di _n "Marital status heterogeneity:"
di "  Married:   b = " %7.4f `b_marr' " (SE = " %7.4f `se_marr' "), N = " `n_marr'
di "  Unmarried: b = " %7.4f `b_unmarr' " (SE = " %7.4f `se_unmarr' "), N = " `n_unmarr'

local diff2 = `b_marr' - `b_unmarr'
local se_diff2 = sqrt(`se_marr'^2 + `se_unmarr'^2)
local z_diff2 = `diff2' / `se_diff2'
di "  Difference (M-U): " %7.4f `diff2' " (SE = " %7.4f `se_diff2' ")"
di "  Z-statistic: " %7.3f `z_diff2'

* --------------------------------------------------------------------------
* 3. INTERACTION MODEL (FULL SAMPLE)
* --------------------------------------------------------------------------
di _n "=== INTERACTION MODEL ==="

gen mobile_female = mobile_user * female
gen mobile_married = mobile_user * married

reghdfe se mobile_user female married mobile_female mobile_married ///
    age age2 has_children ///
    i.praceeth3 i.peducgrp i.hhincome ///
    [pw=hsupwgtk], absorb(cbsa year) cluster(cbsa)
eststo interact

di _n "Interaction results:"
di "  mobile_user (base):         " %7.4f _b[mobile_user] " (" %7.4f _se[mobile_user] ")"
di "  mobile × female:            " %7.4f _b[mobile_female] " (" %7.4f _se[mobile_female] ")"
di "  mobile × married:           " %7.4f _b[mobile_married] " (" %7.4f _se[mobile_married] ")"
di "  female (main):              " %7.4f _b[female] " (" %7.4f _se[female] ")"
di "  married (main):             " %7.4f _b[married] " (" %7.4f _se[married] ")"

* Test joint significance of interactions
test mobile_female mobile_married
di "  Joint F-test p-value: " %7.4f r(p)

* --------------------------------------------------------------------------
* 4. EXPORT RESULTS
* --------------------------------------------------------------------------
di _n "=== EXPORT ==="

esttab fem male marr unmarr interact, ///
    keep(mobile_user mobile_female mobile_married female married) ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2_a) ///
    mtitles("Female" "Male" "Married" "Unmarried" "Interactions") ///
    title("Heterogeneity by Gender and Marital Status")

esttab fem male marr unmarr interact using ///
    "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/_Research/Mobile_Money_Banking/Mobile banking USA/Data/output/table_het_gender_marital.csv", ///
    replace csv ///
    keep(mobile_user mobile_female mobile_married female married) ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2_a) ///
    mtitles("Female" "Male" "Married" "Unmarried" "Interactions")

di _n "=============================================="
di "HETEROGENEITY ANALYSIS COMPLETE"
di "=============================================="
