/*******************************************************************************
* Phase 1 Completion: Extract Results and Compute Self-Employment Rates
*******************************************************************************/

clear all
set more off

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase1_completion.log", replace

use "$datadir/analysis_dataset_with_se.dta", clear

* Apply sample restrictions
keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .
keep if banking_mode != .

* Create employment status
capture drop emp_status
gen emp_status = .
replace emp_status = 1 if wage_worker == 1
replace emp_status = 2 if self_employed == 1
replace emp_status = 3 if employed == 0

* Create joint choice
capture drop joint_choice
gen joint_choice = (banking_mode - 1) * 3 + emp_status
keep if joint_choice >= 1 & joint_choice <= 9

* Create covariates
capture drop age_cat
gen age_cat = .
replace age_cat = 1 if age >= 18 & age <= 29
replace age_cat = 2 if age >= 30 & age <= 44
replace age_cat = 3 if age >= 45 & age <= 64

capture drop educ_hs educ_somecoll educ_college
gen educ_hs = (peducgrp == 2)
gen educ_somecoll = (peducgrp == 3)
gen educ_college = (peducgrp == 4)

capture drop race_black race_hispanic race_other
gen race_black = (praceeth == 2) if praceeth != .
gen race_hispanic = (praceeth == 3 | hispanic == 1) if praceeth != .
gen race_other = (praceeth == 4 | other_race == 1) if praceeth != .

capture drop broadband_std
sum pct_broadband
gen broadband_std = (pct_broadband - r(mean)) / r(sd)

tab year, gen(yr_)

/*******************************************************************************
* Re-estimate MNL (Base Model)
*******************************************************************************/

di _n "=== Estimating Base Multinomial Logit ==="

mlogit joint_choice i.age_cat educ_hs educ_somecoll educ_college ///
    race_black race_hispanic ///
    broadband_std ///
    yr_2 yr_3 yr_4 yr_5 yr_6 ///
    [pw=weight], base(7) vce(cluster cbsa)

estimates store mlogit_base

/*******************************************************************************
* Self-Employment Rates by Banking Mode
*******************************************************************************/

di _n "=== Computing Predicted Probabilities ==="

* Get predicted probabilities at means
margins, atmeans predict(pr outcome(#1)) predict(pr outcome(#2)) predict(pr outcome(#3)) ///
                predict(pr outcome(#4)) predict(pr outcome(#5)) predict(pr outcome(#6)) ///
                predict(pr outcome(#7)) predict(pr outcome(#8)) predict(pr outcome(#9)) post

matrix P = e(b)

di _n "=== Predicted Choice Probabilities (at means) ==="
di "P(Unbanked × Wage):    " %8.5f P[1,1]
di "P(Unbanked × SE):      " %8.5f P[1,2]
di "P(Unbanked × NotWork): " %8.5f P[1,3]
di "P(Mobile × Wage):      " %8.5f P[1,4]
di "P(Mobile × SE):        " %8.5f P[1,5]
di "P(Mobile × NotWork):   " %8.5f P[1,6]
di "P(Branch × Wage):      " %8.5f P[1,7]
di "P(Branch × SE):        " %8.5f P[1,8]
di "P(Branch × NotWork):   " %8.5f P[1,9]

* Compute conditional SE rates
scalar p_unbanked = P[1,1] + P[1,2] + P[1,3]
scalar p_mobile = P[1,4] + P[1,5] + P[1,6]
scalar p_branch = P[1,7] + P[1,8] + P[1,9]

scalar se_rate_unbanked = P[1,2] / p_unbanked
scalar se_rate_mobile = P[1,5] / p_mobile
scalar se_rate_branch = P[1,8] / p_branch

di _n "=== Banking Mode Shares ==="
di "P(Unbanked): " %8.5f p_unbanked
di "P(Mobile):   " %8.5f p_mobile
di "P(Branch):   " %8.5f p_branch

di _n "=== Self-Employment Rates Conditional on Banking Mode ==="
di "SE rate | Unbanked: " %8.5f se_rate_unbanked
di "SE rate | Mobile:   " %8.5f se_rate_mobile
di "SE rate | Branch:   " %8.5f se_rate_branch

/*******************************************************************************
* Counterfactual: Broadband Expansion Effect
*******************************************************************************/

di _n "=== Counterfactual: +1 SD Broadband Increase ==="

estimates restore mlogit_base

* Baseline predictions
margins, atmeans predict(pr outcome(#4)) predict(pr outcome(#5)) ///
                 predict(pr outcome(#7)) predict(pr outcome(#8))
matrix baseline = r(b)

* +1 SD broadband
margins, at(broadband_std=1) atmeans predict(pr outcome(#4)) predict(pr outcome(#5)) ///
                                      predict(pr outcome(#7)) predict(pr outcome(#8))
matrix high_bb = r(b)

di _n "=== Effect of +1 SD Broadband on Choice Probabilities ==="
di "                        Baseline    +1SD BB    Change"
di "P(Mobile × Wage):      " %8.5f baseline[1,1] "   " %8.5f high_bb[1,1] "   " %8.5f (high_bb[1,1] - baseline[1,1])
di "P(Mobile × SE):        " %8.5f baseline[1,2] "   " %8.5f high_bb[1,2] "   " %8.5f (high_bb[1,2] - baseline[1,2])
di "P(Branch × Wage):      " %8.5f baseline[1,3] "   " %8.5f high_bb[1,3] "   " %8.5f (high_bb[1,3] - baseline[1,3])
di "P(Branch × SE):        " %8.5f baseline[1,4] "   " %8.5f high_bb[1,4] "   " %8.5f (high_bb[1,4] - baseline[1,4])

/*******************************************************************************
* Summary Table for Paper
*******************************************************************************/

di _n "=== PHASE 1 SUMMARY FOR PAPER ==="
di "Sample size: " _N
di ""
di "Key Findings:"
di "1. Self-employment rate among branch users (" %5.3f se_rate_branch ") is higher than"
di "   mobile users (" %5.3f se_rate_mobile ") even after controlling for demographics."
di ""
di "2. Broadband expansion (+1 SD) increases:"
di "   - P(Mobile × SE) by " %6.4f (high_bb[1,2] - baseline[1,2])
di "   - P(Branch × SE) by " %6.4f (high_bb[1,4] - baseline[1,4])
di ""
di "3. Age is strongly associated with self-employment across all banking modes"
di "   (older workers more likely to be self-employed)"
di ""
di "4. Education reduces unbanked status but has mixed effects on self-employment"

/*******************************************************************************
* Save Results to CSV
*******************************************************************************/

* Create results matrix
matrix results = J(9, 3, .)
forvalues j = 1/9 {
    matrix results[`j', 1] = `j'
    matrix results[`j', 2] = P[1,`j']
}
matrix results[1,3] = se_rate_unbanked
matrix results[2,3] = .
matrix results[3,3] = .
matrix results[4,3] = se_rate_mobile
matrix results[5,3] = .
matrix results[6,3] = .
matrix results[7,3] = se_rate_branch
matrix results[8,3] = .
matrix results[9,3] = .

matrix colnames results = choice prob cond_se_rate

clear
svmat results, names(col)

gen choice_label = ""
replace choice_label = "Unbanked × Wage" if choice == 1
replace choice_label = "Unbanked × SE" if choice == 2
replace choice_label = "Unbanked × NotWork" if choice == 3
replace choice_label = "Mobile × Wage" if choice == 4
replace choice_label = "Mobile × SE" if choice == 5
replace choice_label = "Mobile × NotWork" if choice == 6
replace choice_label = "Branch × Wage" if choice == 7
replace choice_label = "Branch × SE" if choice == 8
replace choice_label = "Branch × NotWork" if choice == 9

export delimited "$output/phase1_choice_probs.csv", replace

di _n "Results saved to $output/phase1_choice_probs.csv"

log close
