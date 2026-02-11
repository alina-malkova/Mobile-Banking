/*******************************************************************************
* Phase 2: Dynamic Discrete Choice via Arcidiacono-Miller CCP Method
*******************************************************************************/

clear all
set more off
set matsize 11000

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase2_dynamic_ccp.log", replace

di _n "============================================================"
di "PHASE 2: DYNAMIC DISCRETE CHOICE VIA CCP METHOD"
di "============================================================"

/*******************************************************************************
* 1. Load Individual-Level Data and Construct State Variables
*******************************************************************************/

use "$datadir/analysis_dataset_with_se.dta", clear

* Sample restrictions
keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .
keep if banking_mode != .

* Employment status
capture drop emp_status
gen emp_status = .
replace emp_status = 1 if wage_worker == 1
replace emp_status = 2 if self_employed == 1
replace emp_status = 3 if employed != 1

tab emp_status, missing

* Banking mode
capture drop bank_mode
gen bank_mode = banking_mode
tab bank_mode, missing

* Joint choice variable (3 banking × 3 employment = 9 alternatives)
capture drop joint_choice
gen joint_choice = (bank_mode - 1) * 3 + emp_status
tab joint_choice

* State variables for dynamic model
gen age_cat = 1 if age >= 18 & age < 30
replace age_cat = 2 if age >= 30 & age < 45
replace age_cat = 3 if age >= 45 & age <= 64

* Standardize broadband (using pct_broadband)
sum pct_broadband [aw=hsupwgtk]
scalar bb_mean = r(mean)
scalar bb_sd = r(sd)
gen bb_std = (pct_broadband - bb_mean) / bb_sd

* Education categories
gen educ = .
replace educ = 1 if no_hs == 1
replace educ = 2 if hs_diploma == 1
replace educ = 3 if some_college == 1
replace educ = 4 if college_degree == 1

/*******************************************************************************
* 2. Estimate CCPs by Cell
*******************************************************************************/

di _n "=== Step 1: CCP Estimation ==="

* Create dummies for collapse
forvalues j = 1/9 {
    gen choice`j' = (joint_choice == `j')
}

* Collapse to cell level
collapse (mean) p1=choice1 p2=choice2 p3=choice3 p4=choice4 p5=choice5 ///
         p6=choice6 p7=choice7 p8=choice8 p9=choice9 ///
         bb_std ///
         (count) n=joint_choice ///
         [aw=hsupwgtk], by(cbsa year age_cat educ)

* Drop small cells
drop if n < 25

di "Cell-level observations: " _N

/*******************************************************************************
* 3. Store CCPs
*******************************************************************************/

save "$output/ccps_dynamic.dta", replace

/*******************************************************************************
* 4. Normalize and compute log-CCPs
*******************************************************************************/

* Add small constant for numerical stability
scalar eps = 0.001
forvalues j = 1/9 {
    replace p`j' = max(p`j', eps)
}

* Normalize to sum to 1
gen psum = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9
forvalues j = 1/9 {
    replace p`j' = p`j' / psum
}
drop psum

* Compute log-CCPs
forvalues j = 1/9 {
    gen lnp`j' = ln(p`j')
}

* Inclusive value (logsum)
gen inclusive = ln(exp(lnp1) + exp(lnp2) + exp(lnp3) + exp(lnp4) + ///
                  exp(lnp5) + exp(lnp6) + exp(lnp7) + exp(lnp8) + exp(lnp9))

/*******************************************************************************
* 5. Construct Y for Berry Inversion
*******************************************************************************/

di _n "=== Step 2: Construct Berry-style Dependent Variable ==="

* Normalize to alternative 7 (Branch × Wage = renewal action)
forvalues j = 1/9 {
    gen y`j' = lnp`j' - lnp7
}

/*******************************************************************************
* 6. Dynamic Structural Estimation
*******************************************************************************/

di _n "=== Step 3: Dynamic Structural Estimation ==="

scalar beta = 0.90

* Stack alternatives
preserve
gen obs = _n
expand 9
bysort obs: gen alt = _n
gen y = .
forvalues j = 1/9 {
    replace y = y`j' if alt == `j'
}

* Alternative dummies
forvalues j = 1/9 {
    gen a`j' = (alt == `j')
}

* Banking indicators
gen is_mobile = inlist(alt, 4, 5, 6)
gen is_branch = inlist(alt, 7, 8, 9)
gen is_unbanked = inlist(alt, 1, 2, 3)

* Employment indicators
gen is_se = inlist(alt, 2, 5, 8)
gen is_wage = inlist(alt, 1, 4, 7)
gen is_notwork = inlist(alt, 3, 6, 9)

* State interactions
gen age2 = (age_cat == 2)
gen age3 = (age_cat == 3)
gen educ4 = (educ == 4)

gen age2_se = age2 * is_se
gen age3_se = age3 * is_se
gen educ4_se = educ4 * is_se
gen bb_mobile = bb_std * is_mobile
gen bb_se = bb_std * is_se

* Dynamic terms
gen dynamic_se = is_se * age_cat * beta
gen switch_cost = (is_mobile + is_unbanked)

di _n "=== Dynamic Structural Estimation ==="

reg y a1 a2 a3 a4 a5 a6 a8 a9 ///
    age2_se age3_se educ4_se ///
    bb_mobile bb_se ///
    dynamic_se switch_cost ///
    [aw=n], vce(cluster cbsa)

estimates store dynamic_model

* Extract key parameters
scalar gamma_age2_se = _b[age2_se]
scalar gamma_age3_se = _b[age3_se]
scalar gamma_bb_mobile = _b[bb_mobile]
scalar gamma_bb_se = _b[bb_se]
scalar gamma_dynamic = _b[dynamic_se]
scalar kappa_switch = _b[switch_cost]

di _n "=== Structural Parameter Estimates ==="
di "Age 30-44 × SE:     " %8.4f gamma_age2_se
di "Age 45-64 × SE:     " %8.4f gamma_age3_se
di "Broadband × Mobile: " %8.4f gamma_bb_mobile
di "Broadband × SE:     " %8.4f gamma_bb_se
di "Dynamic SE return:  " %8.4f gamma_dynamic
di "Switching cost:     " %8.4f kappa_switch

restore

/*******************************************************************************
* 7. Compute Implied Self-Employment Returns
*******************************************************************************/

di _n "=== Implied Returns to Self-Employment ==="

scalar exp_return = gamma_dynamic / beta
di "Implied experience return (per period): " %8.4f exp_return

scalar pv_exp = exp_return * beta / (1 - beta)
di "Present value of SE career (vs wage): " %8.4f pv_exp

di "Credit access: Branch vs Mobile = " %8.4f (gamma_bb_se - gamma_bb_mobile)

/*******************************************************************************
* 8. Counterfactuals with Dynamic Model
*******************************************************************************/

di _n "=== Dynamic Counterfactual: Branch Closure ==="

* Load CCPs
use "$output/ccps_dynamic.dta", clear

* Baseline SE rate
egen se_base = rowtotal(p2 p5 p8)
sum se_base [aw=n]
scalar baseline_se = r(mean)

preserve

* Shift from branch to mobile (80%) and unbanked (20%)
gen branch_total = p7 + p8 + p9
gen displaced = 0.50 * branch_total

gen branch_se_sh = p8 / (branch_total + 0.001)
gen branch_wage_sh = p7 / (branch_total + 0.001)
gen branch_notwork_sh = p9 / (branch_total + 0.001)

scalar se_adj_static = exp(gamma_bb_se - gamma_bb_mobile)
scalar se_adj_dynamic = exp(gamma_dynamic * (-0.5))
scalar se_adj_total = se_adj_static * se_adj_dynamic

di "SE adjustment factor (static):  " %6.4f se_adj_static
di "SE adjustment factor (dynamic): " %6.4f se_adj_dynamic
di "SE adjustment factor (total):   " %6.4f se_adj_total

gen p7_cf = p7 * 0.50
gen p8_cf = p8 * 0.50
gen p9_cf = p9 * 0.50

gen p4_cf = p4 + displaced * 0.80 * branch_wage_sh
gen p5_cf = p5 + displaced * 0.80 * branch_se_sh * se_adj_total
gen p6_cf = p6 + displaced * 0.80 * branch_notwork_sh

gen p1_cf = p1 + displaced * 0.20 * branch_wage_sh
gen p2_cf = p2 + displaced * 0.20 * branch_se_sh * 0.3
gen p3_cf = p3 + displaced * 0.20 * branch_notwork_sh

gen psum_cf = p1_cf + p2_cf + p3_cf + p4_cf + p5_cf + p6_cf + p7_cf + p8_cf + p9_cf
forvalues j = 1/9 {
    replace p`j'_cf = p`j'_cf / psum_cf
}

gen se_cf = p2_cf + p5_cf + p8_cf
sum se_cf [aw=n]
scalar cf_se_dynamic = r(mean)

local change = cf_se_dynamic - baseline_se
local pct = 100 * `change' / baseline_se

di _n "Dynamic Counterfactual Results:"
di "Baseline SE rate:        " %6.4f baseline_se
di "After 50% closure:       " %6.4f cf_se_dynamic
di "Change:                  " %7.4f `change' " (" %5.1f `pct' "%)"

restore

/*******************************************************************************
* 9. Summary
*******************************************************************************/

di _n "============================================================"
di "PHASE 2 SUMMARY: DYNAMIC STRUCTURAL ESTIMATION"
di "============================================================"

di _n "Model: Dynamic discrete choice with finite dependence"
di "Renewal action: Branch × Wage employment"
di "Discount factor: beta = " beta

di _n "Key Structural Parameters:"
di "  Experience return (SE):    " %8.4f gamma_dynamic
di "  Switching cost:            " %8.4f kappa_switch
di "  Broadband -> Mobile:       " %8.4f gamma_bb_mobile
di "  Broadband -> SE:           " %8.4f gamma_bb_se

di _n "Counterfactual (50% branch closure):"
di "  Baseline SE:    " %6.4f baseline_se
di "  Dynamic effect: " %6.4f cf_se_dynamic
di "  % Change:       " %5.1f 100*(cf_se_dynamic - baseline_se)/baseline_se "%"

/*******************************************************************************
* 10. Save Results
*******************************************************************************/

clear
set obs 6
gen parameter = ""
gen value = .

replace parameter = "gamma_dynamic" in 1
replace value = gamma_dynamic in 1

replace parameter = "kappa_switch" in 2
replace value = kappa_switch in 2

replace parameter = "gamma_bb_mobile" in 3
replace value = gamma_bb_mobile in 3

replace parameter = "gamma_bb_se" in 4
replace value = gamma_bb_se in 4

replace parameter = "baseline_se" in 5
replace value = baseline_se in 5

replace parameter = "cf_se_dynamic" in 6
replace value = cf_se_dynamic in 6

export delimited using "$output/phase2_dynamic_params.csv", replace

log close

di _n "Phase 2 complete. Results saved to:"
di "  $output/phase2_dynamic_params.csv"
di "  $output/phase2_dynamic_ccp.log"
