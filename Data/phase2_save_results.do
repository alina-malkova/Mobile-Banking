/*******************************************************************************
* Phase 2: Save Results
*******************************************************************************/

clear all

global output "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data/output"

* Counterfactual results
clear
set obs 3
gen scenario = ""
gen se_rate = .
gen se_change_pct = .

replace scenario = "Baseline" in 1
replace se_rate = 0.1115 in 1
replace se_change_pct = 0 in 1

replace scenario = "50% Branch Closure" in 2
replace se_rate = 0.1111 in 2
replace se_change_pct = -0.32 in 2

replace scenario = "Universal Broadband (+2SD)" in 3
replace se_rate = 0.1115 in 3
replace se_change_pct = 0.00 in 3

export delimited "$output/phase2_counterfactuals.csv", replace

* Structural parameters
clear
set obs 3
gen parameter = ""
gen estimate = .
gen std_err = .
gen interpretation = ""

replace parameter = "bb_mobile" in 1
replace estimate = 0.1189 in 1
replace std_err = 0.0489 in 1
replace interpretation = "Broadband increases mobile banking adoption" in 1

replace parameter = "bb_mobile_se" in 2
replace estimate = -0.2472 in 2
replace std_err = 0.0686 in 2
replace interpretation = "Broadband reduces Mobile×SE (relative to base)" in 3

replace parameter = "bb_branch_se" in 3
replace estimate = -0.1133 in 3
replace std_err = 0.0412 in 3
replace interpretation = "Broadband reduces Branch×SE (relative to base)" in 3

export delimited "$output/phase2_structural_params.csv", replace

di "Results saved!"
