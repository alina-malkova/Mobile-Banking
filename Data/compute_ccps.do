/*******************************************************************************
* Compute CCPs for Structural Model
*******************************************************************************/

clear all
set more off

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

use "$datadir/analysis_dataset_with_se.dta", clear

* Apply restrictions
keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .

* Create age categories
gen age_cat = .
replace age_cat = 1 if age >= 18 & age <= 29
replace age_cat = 2 if age >= 30 & age <= 44
replace age_cat = 3 if age >= 45 & age <= 64

* Create joint choice dummies
gen choice_unbanked_wage = (banking_mode == 1 & wage_worker == 1) if banking_mode != .
gen choice_unbanked_se = (banking_mode == 1 & self_employed == 1) if banking_mode != .
gen choice_unbanked_notwork = (banking_mode == 1 & employed == 0) if banking_mode != .
gen choice_mobile_wage = (banking_mode == 2 & wage_worker == 1) if banking_mode != .
gen choice_mobile_se = (banking_mode == 2 & self_employed == 1) if banking_mode != .
gen choice_mobile_notwork = (banking_mode == 2 & employed == 0) if banking_mode != .
gen choice_branch_wage = (banking_mode == 3 & wage_worker == 1) if banking_mode != .
gen choice_branch_se = (banking_mode == 3 & self_employed == 1) if banking_mode != .
gen choice_branch_notwork = (banking_mode == 3 & employed == 0) if banking_mode != .

* Keep observations with valid joint choices
keep if banking_mode != .

* Collapse to CBSA × year × demographic cell
collapse (mean) prob_unbanked_wage=choice_unbanked_wage ///
                prob_unbanked_se=choice_unbanked_se ///
                prob_unbanked_notwork=choice_unbanked_notwork ///
                prob_mobile_wage=choice_mobile_wage ///
                prob_mobile_se=choice_mobile_se ///
                prob_mobile_notwork=choice_mobile_notwork ///
                prob_branch_wage=choice_branch_wage ///
                prob_branch_se=choice_branch_se ///
                prob_branch_notwork=choice_branch_notwork ///
                pct_broadband ///
         (count) n=choice_branch_wage [aw=weight], ///
         by(cbsa year age_cat peducgrp)

* Keep cells with sufficient observations
keep if n >= 30

* Verify probabilities sum to 1
gen prob_sum = prob_unbanked_wage + prob_unbanked_se + prob_unbanked_notwork + ///
               prob_mobile_wage + prob_mobile_se + prob_mobile_notwork + ///
               prob_branch_wage + prob_branch_se + prob_branch_notwork
summarize prob_sum

* Save CCPs
save "$output/ccps_for_structural.dta", replace
export delimited "$output/ccps_for_structural.csv", replace

di "CCPs saved: " _N " cells"
summarize prob_* n pct_broadband

* Summary by education
table peducgrp, stat(mean prob_branch_se prob_mobile_se prob_unbanked_se) nformat(%9.4f)

* Summary by age
table age_cat, stat(mean prob_branch_se prob_mobile_se prob_unbanked_se) nformat(%9.4f)

di "CCP computation complete!"
