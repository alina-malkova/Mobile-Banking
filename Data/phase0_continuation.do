/*******************************************************************************
* Phase 0 Continuation: IV and Heterogeneity Analysis
* Fix: Use state FE instead of CBSA FE for IV (broadband varies by CBSA, not within)
*******************************************************************************/

clear all
set more off

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

log using "$output/phase0_continuation.log", replace

* Load analysis sample
use "$datadir/analysis_dataset_with_se.dta", clear

* Apply same restrictions
keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .

di "Sample size: " _N

********************************************************************************
* PART 6 (Corrected): IV REGRESSION - State FE approach
********************************************************************************

di "=============================================="
di "IV REGRESSION (State FE instead of CBSA FE)"
di "=============================================="

* First stage: Mobile banking on broadband penetration (with state FE)
reghdfe mobile_user pct_broadband ///
    age c.age#c.age ///
    i.praceeth3 i.peducgrp i.hhincome ///
    [pw=weight], absorb(statefips year) cluster(statefips)
eststo fs_state

* Check first stage F-stat
test pct_broadband
scalar fs_F = r(F)
di "First-stage F-statistic: " fs_F

* Reduced form: Self-employment on broadband
reghdfe self_employed pct_broadband ///
    age c.age#c.age ///
    i.praceeth3 i.peducgrp i.hhincome ///
    [pw=weight], absorb(statefips year) cluster(statefips)
eststo rf_state

* IV: Self-employment on instrumented mobile banking (state FE)
ivregress 2sls self_employed ///
    (mobile_user = pct_broadband) ///
    age c.age#c.age ///
    i.praceeth3 i.peducgrp i.hhincome ///
    i.statefips i.year ///
    [pw=weight], cluster(statefips) first
eststo iv_state

* Export IV results
esttab fs_state rf_state iv_state using "$output/table4_iv_results.csv", replace ///
    keep(mobile_user pct_broadband age) ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2) ///
    title("IV Regression: Broadband as Instrument (State FE)") ///
    mtitles("First Stage" "Reduced Form" "IV-2SLS")

********************************************************************************
* PART 7: HETEROGENEITY ANALYSIS
********************************************************************************

di "=============================================="
di "HETEROGENEITY ANALYSIS"
di "=============================================="

* By race/ethnicity
eststo clear
foreach race in 1 2 6 {
    qui reghdfe self_employed mobile_user ///
        age c.age#c.age i.peducgrp i.hhincome ///
        [pw=weight] if praceeth3 == `race', absorb(cbsa year) cluster(cbsa)
    eststo het_race`race'
}

esttab het_race1 het_race2 het_race6 using "$output/table5_het_race.csv", replace ///
    keep(mobile_user) se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2) title("Heterogeneity by Race/Ethnicity") ///
    mtitles("Black" "Hispanic" "White")

* By education
eststo clear
forvalues edu = 1/4 {
    qui reghdfe self_employed mobile_user ///
        age c.age#c.age i.praceeth3 i.hhincome ///
        [pw=weight] if peducgrp == `edu', absorb(cbsa year) cluster(cbsa)
    eststo het_edu`edu'
}

esttab het_edu* using "$output/table5_het_education.csv", replace ///
    keep(mobile_user) se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2) title("Heterogeneity by Education") ///
    mtitles("No HS" "HS Diploma" "Some College" "College+")

* By income
eststo clear
forvalues inc = 1/5 {
    qui reghdfe self_employed mobile_user ///
        age c.age#c.age i.praceeth3 i.peducgrp ///
        [pw=weight] if hhincome == `inc', absorb(cbsa year) cluster(cbsa)
    eststo het_inc`inc'
}

esttab het_inc* using "$output/table5_het_income.csv", replace ///
    keep(mobile_user) se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2) title("Heterogeneity by Income") ///
    mtitles("<15K" "15-30K" "30-50K" "50-75K" "75K+")

* By metro/rural
eststo clear
qui reghdfe self_employed mobile_user ///
    age c.age#c.age i.praceeth3 i.peducgrp i.hhincome ///
    [pw=weight] if metro == 1, absorb(cbsa year) cluster(cbsa)
eststo het_metro

qui reghdfe self_employed mobile_user ///
    age c.age#c.age i.praceeth3 i.peducgrp i.hhincome ///
    [pw=weight] if rural == 1, absorb(statefips year) cluster(statefips)
eststo het_rural

esttab het_metro het_rural using "$output/table5_het_geography.csv", replace ///
    keep(mobile_user) se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2) title("Heterogeneity by Geography") ///
    mtitles("Metropolitan" "Rural")

********************************************************************************
* PART 8: CBSA-LEVEL ANALYSIS
********************************************************************************

di "=============================================="
di "CBSA-LEVEL PANEL ANALYSIS"
di "=============================================="

preserve
collapse (mean) self_employed mobile_user banked ///
         pct_broadband unemployment_rate ///
         (rawsum) pop=weight [aw=weight], by(cbsa year)

xtset cbsa year
xtreg self_employed mobile_user i.year, fe cluster(cbsa)
eststo cbsa_fe

xtreg self_employed mobile_user pct_broadband i.year, re cluster(cbsa)
eststo cbsa_re

esttab cbsa_fe cbsa_re using "$output/table6_cbsa_level.csv", replace ///
    keep(mobile_user pct_broadband) ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2_w r2_b) title("CBSA-Level Panel Regression") ///
    mtitles("CBSA FE" "Random Effects")
restore

********************************************************************************
* PART 9: CCPs FOR STRUCTURAL MODEL
********************************************************************************

di "=============================================="
di "CONDITIONAL CHOICE PROBABILITIES"
di "=============================================="

* Create age categories
gen age_cat = .
replace age_cat = 1 if age >= 18 & age <= 29
replace age_cat = 2 if age >= 30 & age <= 44
replace age_cat = 3 if age >= 45 & age <= 64

* Create joint choice variable
gen joint_choice = .
* Unbanked
replace joint_choice = 1 if banking_mode == 1 & wage_worker == 1
replace joint_choice = 2 if banking_mode == 1 & self_employed == 1
replace joint_choice = 3 if banking_mode == 1 & employed == 0
* Mobile/Online
replace joint_choice = 4 if banking_mode == 2 & wage_worker == 1
replace joint_choice = 5 if banking_mode == 2 & self_employed == 1
replace joint_choice = 6 if banking_mode == 2 & employed == 0
* Branch
replace joint_choice = 7 if banking_mode == 3 & wage_worker == 1
replace joint_choice = 8 if banking_mode == 3 & self_employed == 1
replace joint_choice = 9 if banking_mode == 3 & employed == 0

* Collapse to CBSA × year × demographic cell
preserve
keep if joint_choice != .

* Compute CCPs within cells
collapse (mean) prob_unbanked_wage = (joint_choice==1) ///
                prob_unbanked_se = (joint_choice==2) ///
                prob_unbanked_notwork = (joint_choice==3) ///
                prob_mobile_wage = (joint_choice==4) ///
                prob_mobile_se = (joint_choice==5) ///
                prob_mobile_notwork = (joint_choice==6) ///
                prob_branch_wage = (joint_choice==7) ///
                prob_branch_se = (joint_choice==8) ///
                prob_branch_notwork = (joint_choice==9) ///
                pct_broadband ///
         (count) n=joint_choice [aw=weight], ///
         by(cbsa year age_cat peducgrp)

* Keep cells with sufficient observations
keep if n >= 30

* Save CCPs
save "$output/ccps_for_structural.dta", replace
export delimited "$output/ccps_for_structural.csv", replace

di "CCPs saved: " _N " cells"
summarize prob_* n
restore

********************************************************************************
* SUMMARY STATISTICS
********************************************************************************

di "=============================================="
di "KEY FINDINGS SUMMARY"
di "=============================================="

di "1. Mobile banking adoption by year:"
table year [aw=weight], stat(mean mobile_user) nformat(%9.3f)

di "2. Self-employment by banking mode:"
table banking_mode [aw=weight] if banking_mode != ., stat(mean self_employed) nformat(%9.3f)

di "3. Sample composition:"
di "   Total observations: " _N
tab year

********************************************************************************
* SAVE FINAL ANALYSIS SAMPLE
********************************************************************************

compress
save "$datadir/analysis_sample_phase0.dta", replace

di "=============================================="
di "PHASE 0 COMPLETE"
di "Output saved to: $output"
di "=============================================="

log close
