/*******************************************************************************
* Phase 1: Static Multinomial Logit for Joint Banking Mode × Employment Choice
*
* Model: 9-alternative multinomial logit
*   Banking mode: {1=unbanked, 2=mobile, 3=branch}
*   Employment:   {1=wage, 2=self-employed, 3=not working}
*   Joint choice: 9 alternatives
*
* Key structural parameters:
*   - Credit access interaction: self-employment × banking mode × infrastructure
*   - Banking cost parameters: how infrastructure affects banking mode choice
*******************************************************************************/

clear all
set more off
set matsize 11000

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase1_structural.log", replace

/*******************************************************************************
* 1. Load and Prepare Data
*******************************************************************************/

use "$datadir/analysis_dataset_with_se.dta", clear

* Apply sample restrictions
keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .
keep if banking_mode != .

di "Sample size after restrictions: " _N

* Create employment status variable
capture drop emp_status
gen emp_status = .
replace emp_status = 1 if wage_worker == 1
replace emp_status = 2 if self_employed == 1
replace emp_status = 3 if employed == 0  // not working (unemployed in labor force)
label define emp_lbl 1 "Wage worker" 2 "Self-employed" 3 "Not working"
label values emp_status emp_lbl

* Create joint choice variable (9 alternatives)
capture drop joint_choice
* Encoding: joint_choice = (banking_mode - 1) * 3 + emp_status
gen joint_choice = (banking_mode - 1) * 3 + emp_status
label define joint_lbl ///
    1 "Unbanked × Wage" ///
    2 "Unbanked × SE" ///
    3 "Unbanked × NotWork" ///
    4 "Mobile × Wage" ///
    5 "Mobile × SE" ///
    6 "Mobile × NotWork" ///
    7 "Branch × Wage" ///
    8 "Branch × SE" ///
    9 "Branch × NotWork"
label values joint_choice joint_lbl

tab joint_choice [aw=weight], missing

* Keep only valid joint choices
keep if joint_choice >= 1 & joint_choice <= 9
di "Final estimation sample: " _N

/*******************************************************************************
* 2. Create Covariates
*******************************************************************************/

* Age categories
capture drop age_cat
gen age_cat = .
replace age_cat = 1 if age >= 18 & age <= 29
replace age_cat = 2 if age >= 30 & age <= 44
replace age_cat = 3 if age >= 45 & age <= 64

* Education dummies (base = less than high school)
capture drop educ_hs educ_somecoll educ_college
gen educ_hs = (peducgrp == 2)
gen educ_somecoll = (peducgrp == 3)
gen educ_college = (peducgrp == 4)

* Race dummies (base = white non-Hispanic)
* Using praceeth: 1=White, 2=Black, 3=Hispanic, 4=Other
capture drop race_black race_hispanic race_other
gen race_black = (praceeth == 2) if praceeth != .
gen race_hispanic = (praceeth == 3 | hispanic == 1) if praceeth != .
gen race_other = (praceeth == 4 | other_race == 1) if praceeth != .

* Income categories
capture drop inc_cat inc_30_50 inc_50_75 inc_75plus
gen inc_cat = .
replace inc_cat = 1 if hhincome < 30000
replace inc_cat = 2 if hhincome >= 30000 & hhincome < 50000
replace inc_cat = 3 if hhincome >= 50000 & hhincome < 75000
replace inc_cat = 4 if hhincome >= 75000 & hhincome != .

gen inc_30_50 = (inc_cat == 2)
gen inc_50_75 = (inc_cat == 3)
gen inc_75plus = (inc_cat == 4)

* Standardize broadband for easier interpretation
capture drop broadband_std
sum pct_broadband
gen broadband_std = (pct_broadband - r(mean)) / r(sd)

* Year dummies
tab year, gen(yr_)

/*******************************************************************************
* 3. Summary Statistics by Joint Choice
*******************************************************************************/

di _n "=== Joint Choice Distribution ==="
tab joint_choice [aw=weight]

di _n "=== Demographics by Joint Choice ==="
table joint_choice, stat(mean age educ_college inc_75plus pct_broadband) nformat(%9.3f)

/*******************************************************************************
* 4. Estimate Multinomial Logit - Base Specification
*******************************************************************************/

di _n "=== Model 1: Base Multinomial Logit ==="

* Base category: Branch × Wage (choice 7) - most common choice
mlogit joint_choice i.age_cat educ_hs educ_somecoll educ_college ///
    race_black race_hispanic ///
    inc_30_50 inc_50_75 inc_75plus ///
    broadband_std ///
    yr_2 yr_3 yr_4 yr_5 yr_6 ///
    [pw=weight], base(7) vce(cluster cbsa)

estimates store mlogit_base

* Save estimation results
esttab mlogit_base using "$output/phase1_mlogit_base.csv", ///
    cells(b(fmt(4)) se(fmt(4) par)) ///
    stats(N ll aic bic, fmt(0 2 2 2)) ///
    title("Multinomial Logit: Joint Banking-Employment Choice") ///
    replace

/*******************************************************************************
* 5. Compute Relative Risk Ratios
*******************************************************************************/

di _n "=== Relative Risk Ratios ==="
mlogit, rrr

/*******************************************************************************
* 6. Marginal Effects - Key Parameters
*******************************************************************************/

di _n "=== Average Marginal Effects ==="

* AME of college education on each choice
margins, dydx(educ_college) predict(outcome(#1)) predict(outcome(#2)) ///
    predict(outcome(#4)) predict(outcome(#5)) predict(outcome(#7)) predict(outcome(#8))

* AME of broadband on each choice
margins, dydx(broadband_std) predict(outcome(#1)) predict(outcome(#2)) ///
    predict(outcome(#4)) predict(outcome(#5)) predict(outcome(#7)) predict(outcome(#8))

/*******************************************************************************
* 7. Predicted Probabilities by Demographics
*******************************************************************************/

di _n "=== Predicted Choice Probabilities by Education ==="

* Predictions at different education levels
margins, at(educ_hs=0 educ_somecoll=0 educ_college=0) ///
         at(educ_hs=1 educ_somecoll=0 educ_college=0) ///
         at(educ_hs=0 educ_somecoll=1 educ_college=0) ///
         at(educ_hs=0 educ_somecoll=0 educ_college=1) ///
         predict(outcome(#5)) predict(outcome(#8))  // Mobile×SE and Branch×SE

marginsplot, name(se_by_educ, replace) ///
    title("Pr(Self-Employment) by Banking Mode and Education")
graph export "$output/phase1_se_by_education.png", replace

/*******************************************************************************
* 8. Model with Interactions - Credit Access Channel
*******************************************************************************/

di _n "=== Model 2: With Credit Access Interactions ==="

* Create interaction terms for credit access channel
* The key structural hypothesis: self-employment benefits more from branch banking
* when branch density is higher (soft information/relationship lending)

* First, we need to estimate separate effects by banking mode
* This is done by examining the coefficients on demographics across outcomes

* Alternative specification: nested logit structure
* Level 1: Choose banking mode
* Level 2: Choose employment status conditional on banking mode

di _n "=== Testing Credit Access Channel ==="

* Test if broadband effect differs across self-employment outcomes
mlogit joint_choice i.age_cat educ_hs educ_somecoll educ_college ///
    race_black race_hispanic ///
    inc_30_50 inc_50_75 inc_75plus ///
    c.broadband_std##i.age_cat ///
    yr_2 yr_3 yr_4 yr_5 yr_6 ///
    [pw=weight], base(7) vce(cluster cbsa)

estimates store mlogit_interact

* LR test for interaction significance
lrtest mlogit_base mlogit_interact

/*******************************************************************************
* 9. Extract Structural Parameters
*******************************************************************************/

di _n "=== Structural Parameter Interpretation ==="

* Key parameters of interest:
* 1. Banking mode effects on self-employment (comparing outcomes 2,5,8 vs 1,4,7)
* 2. How broadband shifts mobile banking adoption
* 3. Demographic selection into banking modes

* Restore base model
estimates restore mlogit_base

* Calculate implied probabilities at mean values
margins, atmeans predict(pr outcome(#1)) predict(pr outcome(#2)) predict(pr outcome(#3)) ///
                predict(pr outcome(#4)) predict(pr outcome(#5)) predict(pr outcome(#6)) ///
                predict(pr outcome(#7)) predict(pr outcome(#8)) predict(pr outcome(#9))

* Store predicted probabilities
matrix P = r(b)
matrix list P

* Calculate self-employment rate by banking mode
* P(SE | Unbanked) = P(2) / (P(1) + P(2) + P(3))
* P(SE | Mobile)   = P(5) / (P(4) + P(5) + P(6))
* P(SE | Branch)   = P(8) / (P(7) + P(8) + P(9))

scalar p_unbanked_total = P[1,1] + P[1,2] + P[1,3]
scalar p_mobile_total = P[1,4] + P[1,5] + P[1,6]
scalar p_branch_total = P[1,7] + P[1,8] + P[1,9]

scalar se_rate_unbanked = P[1,2] / p_unbanked_total
scalar se_rate_mobile = P[1,5] / p_mobile_total
scalar se_rate_branch = P[1,8] / p_branch_total

di _n "=== Self-Employment Rates by Banking Mode (at means) ==="
di "Unbanked: " %6.4f se_rate_unbanked
di "Mobile:   " %6.4f se_rate_mobile
di "Branch:   " %6.4f se_rate_branch

/*******************************************************************************
* 10. Counterfactual: Effect of Broadband Expansion
*******************************************************************************/

di _n "=== Counterfactual: +1 SD Broadband Increase ==="

* Current predictions
margins, atmeans predict(pr outcome(#5)) predict(pr outcome(#8)) post
matrix current = r(b)

* Restore and predict at higher broadband
estimates restore mlogit_base
margins, at(broadband_std=1) atmeans predict(pr outcome(#5)) predict(pr outcome(#8)) post
matrix higher_bb = r(b)

di "Current P(Mobile × SE):  " %8.5f current[1,1]
di "Current P(Branch × SE):  " %8.5f current[1,2]
di "+1SD BB P(Mobile × SE): " %8.5f higher_bb[1,1]
di "+1SD BB P(Branch × SE): " %8.5f higher_bb[1,2]
di "Change in Mobile × SE:   " %8.5f (higher_bb[1,1] - current[1,1])
di "Change in Branch × SE:   " %8.5f (higher_bb[1,2] - current[1,2])

/*******************************************************************************
* 11. Alternative: Conditional Logit with Alternative-Specific Variables
*******************************************************************************/

di _n "=== Reshaping for Conditional Logit ==="

* Reshape data to long format for conditional logit
preserve

* Keep necessary variables
keep joint_choice age_cat educ_* race_* inc_* broadband_std yr_* weight cbsa year ///
    pct_broadband

* Create individual ID
gen id = _n

* Reshape to long
reshape long, i(id) j(alt)

* The reshape didn't work as expected - need different approach
* For conditional logit, manually expand

restore

* Instead, estimate using asclogit (alternative-specific conditional logit)
* This requires data in long format

preserve

* Create ID
gen id = _n

* Save wide version
save "$datadir/temp_wide.dta", replace

* Create long version for conditional logit
expand 9
bysort id: gen alt = _n

* Create choice indicator
gen chosen = (joint_choice == alt)

* Create alternative-specific dummies
gen is_unbanked = inlist(alt, 1, 2, 3)
gen is_mobile = inlist(alt, 4, 5, 6)
gen is_branch = inlist(alt, 7, 8, 9)

gen is_se = inlist(alt, 2, 5, 8)
gen is_notwork = inlist(alt, 3, 6, 9)

* Interaction: self-employment × broadband × banking mode
gen se_broadband_mobile = is_se * is_mobile * broadband_std
gen se_broadband_branch = is_se * is_branch * broadband_std

* Alternative-specific constants (base = branch × wage = alt 7)
forvalues a = 1/9 {
    gen alt_`a' = (alt == `a')
}

di _n "=== Conditional Logit with Alternative-Specific Variables ==="

* Note: This is computationally intensive with large sample
* Subsample for initial estimation

sample 20, count(id)

clogit chosen alt_1 alt_2 alt_3 alt_4 alt_5 alt_6 alt_8 alt_9 ///
    c.educ_college#i.alt c.broadband_std#i.alt ///
    se_broadband_mobile se_broadband_branch ///
    , group(id) vce(cluster cbsa)

estimates store clogit_model

restore

/*******************************************************************************
* 12. Save Results Summary
*******************************************************************************/

di _n "=== Phase 1 Estimation Complete ==="

* Create summary table
estimates restore mlogit_base

* Export coefficient table for key outcomes (SE alternatives)
esttab mlogit_base using "$output/phase1_structural_results.csv", ///
    keep(2.age_cat 3.age_cat educ_college broadband_std) ///
    cells(b(fmt(4)) se(fmt(4) par) p(fmt(3))) ///
    eqlabels("Unbanked×Wage" "Unbanked×SE" "Unbanked×NotWork" ///
             "Mobile×Wage" "Mobile×SE" "Mobile×NotWork" ///
             "Branch×Wage" "Branch×SE" "Branch×NotWork", none) ///
    stats(N ll, fmt(0 2)) ///
    title("Phase 1: Structural MNL Estimates") ///
    replace

* Save predicted probabilities
estimates restore mlogit_base
margins, atmeans predict(pr) post
matrix probs = r(b)

putexcel set "$output/phase1_predicted_probs.xlsx", replace
putexcel A1 = "Choice"
putexcel B1 = "Predicted Probability"
putexcel A2 = "Unbanked × Wage"
putexcel A3 = "Unbanked × SE"
putexcel A4 = "Unbanked × NotWork"
putexcel A5 = "Mobile × Wage"
putexcel A6 = "Mobile × SE"
putexcel A7 = "Mobile × NotWork"
putexcel A8 = "Branch × Wage"
putexcel A9 = "Branch × SE"
putexcel A10 = "Branch × NotWork"

forvalues j = 1/9 {
    local row = `j' + 1
    putexcel B`row' = probs[1,`j']
}

di _n "Phase 1 output saved to $output/"
di "- phase1_mlogit_base.csv"
di "- phase1_structural_results.csv"
di "- phase1_predicted_probs.xlsx"
di "- phase1_se_by_education.png"

log close
