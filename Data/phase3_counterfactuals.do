/*******************************************************************************
* Phase 3: Policy Counterfactual Simulations
*
* Using structural estimates from Phase 2 to simulate:
*   1. Projected branch closure scenarios (25%, 50%, 75%)
*   2. Mobile banking subsidies (increase adoption)
*   3. Broadband infrastructure investment (+1 SD, +2 SD, universal)
*   4. Combined policies
*
* Key structural parameters (from Phase 2):
*   - bb_mobile = 0.1189 (broadband → mobile adoption)
*   - bb_mobile_se = -0.2472 (broadband × mobile SE interaction)
*   - bb_branch_se = -0.1133 (broadband × branch SE interaction)
*******************************************************************************/

clear all
set more off
set matsize 5000

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase3_counterfactuals.log", replace

/*******************************************************************************
* 1. Load Baseline Data and Structural Parameters
*******************************************************************************/

di _n "============================================================"
di "PHASE 3: POLICY COUNTERFACTUAL SIMULATIONS"
di "============================================================"

* Structural parameters from Phase 2
scalar bb_mobile = 0.1189
scalar bb_mobile_se = -0.2472
scalar bb_branch_se = -0.1133
scalar age2_se = 0.3718
scalar age3_se = 0.9616
scalar educ4_se = 0.7075

di _n "Structural Parameters:"
di "  Broadband → Mobile:    " %6.4f bb_mobile
di "  Broadband × Mobile SE: " %6.4f bb_mobile_se
di "  Broadband × Branch SE: " %6.4f bb_branch_se

* Load CCP data
import delimited "$output/ccps_for_structural.csv", clear

* Rename probabilities
rename prob_unbanked_wage p1
rename prob_unbanked_se p2
rename prob_unbanked_notwork p3
rename prob_mobile_wage p4
rename prob_mobile_se p5
rename prob_mobile_notwork p6
rename prob_branch_wage p7
rename prob_branch_se p8
rename prob_branch_notwork p9

* Normalize to sum to 1
egen psum = rowtotal(p1-p9)
forvalues j = 1/9 {
    replace p`j' = p`j' / psum
}
drop psum prob_sum

* Calculate baseline outcomes
egen se_baseline = rowtotal(p2 p5 p8)
egen mobile_baseline = rowtotal(p4 p5 p6)
egen branch_baseline = rowtotal(p7 p8 p9)
egen unbanked_baseline = rowtotal(p1 p2 p3)

* Standardize broadband
sum pct_broadband
scalar bb_mean = r(mean)
scalar bb_sd = r(sd)
gen bb_std = (pct_broadband - bb_mean) / bb_sd

di _n "Baseline Summary (weighted by cell size):"
sum se_baseline mobile_baseline branch_baseline unbanked_baseline [aw=n]

* Store baseline values
sum se_baseline [aw=n]
scalar baseline_se = r(mean)
sum mobile_baseline [aw=n]
scalar baseline_mobile = r(mean)
sum branch_baseline [aw=n]
scalar baseline_branch = r(mean)
sum unbanked_baseline [aw=n]
scalar baseline_unbanked = r(mean)

di _n "Baseline Rates:"
di "  Self-employment:  " %6.4f baseline_se
di "  Mobile banking:   " %6.4f baseline_mobile
di "  Branch banking:   " %6.4f baseline_branch
di "  Unbanked:         " %6.4f baseline_unbanked

/*******************************************************************************
* 2. Counterfactual 1: Branch Closure Scenarios
*
* Simulate 25%, 50%, 75% branch closures
* Assumption: Displaced branch users reallocate to mobile (80%) and unbanked (20%)
* SE rate among displaced users adjusts based on new banking mode
*******************************************************************************/

di _n "============================================================"
di "COUNTERFACTUAL 1: BRANCH CLOSURE SCENARIOS"
di "============================================================"

* Function to simulate branch closure
* closure_rate: fraction of branch users who lose access
* mobile_switch: fraction of displaced who switch to mobile (rest go unbanked)

foreach closure in 25 50 75 {

    preserve

    local closure_frac = `closure' / 100
    local mobile_switch = 0.80  // 80% switch to mobile
    local unbank_switch = 0.20  // 20% become unbanked

    di _n "--- `closure'% Branch Closure ---"

    * Calculate displaced branch users
    gen displaced = `closure_frac' * (p7 + p8 + p9)

    * Reallocate displaced users
    * Key insight: SE rate changes when banking mode changes
    * Mobile SE rate is lower than branch SE rate (from structural estimates)

    * Original SE share among branch users
    gen branch_se_share = p8 / (p7 + p8 + p9 + 0.0001)
    gen branch_wage_share = p7 / (p7 + p8 + p9 + 0.0001)
    gen branch_notwork_share = p9 / (p7 + p8 + p9 + 0.0001)

    * SE rate adjustment factor when switching to mobile
    * From structural estimates: mobile SE is relatively lower
    scalar se_adjustment_mobile = exp(bb_mobile_se - bb_branch_se)
    di "SE adjustment factor (branch→mobile): " %6.4f se_adjustment_mobile

    * New probabilities after closure
    * Remaining branch users
    gen p7_new = p7 * (1 - `closure_frac')
    gen p8_new = p8 * (1 - `closure_frac')
    gen p9_new = p9 * (1 - `closure_frac')

    * Displaced users switching to mobile
    gen mobile_inflow = displaced * `mobile_switch'
    gen p4_new = p4 + mobile_inflow * branch_wage_share
    gen p5_new = p5 + mobile_inflow * branch_se_share * se_adjustment_mobile
    gen p6_new = p6 + mobile_inflow * branch_notwork_share

    * Displaced users becoming unbanked
    gen unbank_inflow = displaced * `unbank_switch'
    * Unbanked have even lower SE rate
    scalar se_adjustment_unbank = 0.5  // Assume 50% SE rate relative to branch
    gen p1_new = p1 + unbank_inflow * branch_wage_share
    gen p2_new = p2 + unbank_inflow * branch_se_share * se_adjustment_unbank
    gen p3_new = p3 + unbank_inflow * branch_notwork_share

    * Normalize
    egen psum_new = rowtotal(p1_new p2_new p3_new p4_new p5_new p6_new p7_new p8_new p9_new)
    forvalues j = 1/9 {
        replace p`j'_new = p`j'_new / psum_new
    }

    * Calculate new outcomes
    egen se_cf`closure' = rowtotal(p2_new p5_new p8_new)
    egen mobile_cf`closure' = rowtotal(p4_new p5_new p6_new)
    egen branch_cf`closure' = rowtotal(p7_new p8_new p9_new)
    egen unbanked_cf`closure' = rowtotal(p1_new p2_new p3_new)

    sum se_cf`closure' [aw=n]
    scalar se_`closure' = r(mean)
    sum mobile_cf`closure' [aw=n]
    scalar mobile_`closure' = r(mean)
    sum branch_cf`closure' [aw=n]
    scalar branch_`closure' = r(mean)
    sum unbanked_cf`closure' [aw=n]
    scalar unbanked_`closure' = r(mean)

    di "Results after `closure'% branch closure:"
    di "  SE rate:        " %6.4f se_`closure' " (change: " %+6.4f (se_`closure' - baseline_se) ")"
    di "  Mobile rate:    " %6.4f mobile_`closure' " (change: " %+6.4f (mobile_`closure' - baseline_mobile) ")"
    di "  Branch rate:    " %6.4f branch_`closure' " (change: " %+6.4f (branch_`closure' - baseline_branch) ")"
    di "  Unbanked rate:  " %6.4f unbanked_`closure' " (change: " %+6.4f (unbanked_`closure' - baseline_unbanked) ")"

    restore
}

/*******************************************************************************
* 3. Counterfactual 2: Mobile Banking Subsidies
*
* Policy: Subsidize mobile banking to increase adoption
* Simulate 10%, 25%, 50% increase in mobile banking take-up
*******************************************************************************/

di _n "============================================================"
di "COUNTERFACTUAL 2: MOBILE BANKING SUBSIDIES"
di "============================================================"

foreach subsidy in 10 25 50 {

    preserve

    local increase_frac = `subsidy' / 100

    di _n "--- `subsidy'% Mobile Banking Subsidy ---"

    * Increase mobile adoption by drawing from branch and unbanked
    * Assume subsidy primarily affects branch users (they have smartphones)

    gen mobile_current = p4 + p5 + p6
    gen branch_current = p7 + p8 + p9

    * New mobile users come from branch (70%) and unbanked (30%)
    gen new_mobile = mobile_current * `increase_frac'
    gen from_branch = new_mobile * 0.70
    gen from_unbanked = new_mobile * 0.30

    * Proportional shift maintaining employment composition
    gen branch_se_share = p8 / (branch_current + 0.0001)
    gen unbank_se_share = p2 / (p1 + p2 + p3 + 0.0001)

    * New probabilities
    gen p4_new = p4 + from_branch * (p7/(branch_current+0.0001)) + from_unbanked * (p1/(p1+p2+p3+0.0001))
    gen p5_new = p5 + from_branch * branch_se_share * se_adjustment_mobile + from_unbanked * unbank_se_share * 1.5
    gen p6_new = p6 + from_branch * (p9/(branch_current+0.0001)) + from_unbanked * (p3/(p1+p2+p3+0.0001))

    gen p7_new = p7 - from_branch * (p7/(branch_current+0.0001))
    gen p8_new = p8 - from_branch * branch_se_share
    gen p9_new = p9 - from_branch * (p9/(branch_current+0.0001))

    gen p1_new = p1 - from_unbanked * (p1/(p1+p2+p3+0.0001))
    gen p2_new = p2 - from_unbanked * unbank_se_share
    gen p3_new = p3 - from_unbanked * (p3/(p1+p2+p3+0.0001))

    * Ensure non-negative
    forvalues j = 1/9 {
        replace p`j'_new = max(p`j'_new, 0.001)
    }

    * Normalize
    egen psum_new = rowtotal(p1_new p2_new p3_new p4_new p5_new p6_new p7_new p8_new p9_new)
    forvalues j = 1/9 {
        replace p`j'_new = p`j'_new / psum_new
    }

    * Calculate outcomes
    egen se_ms`subsidy' = rowtotal(p2_new p5_new p8_new)
    egen mobile_ms`subsidy' = rowtotal(p4_new p5_new p6_new)

    sum se_ms`subsidy' [aw=n]
    scalar se_ms_`subsidy' = r(mean)
    sum mobile_ms`subsidy' [aw=n]
    scalar mobile_ms_`subsidy' = r(mean)

    di "Results after `subsidy'% mobile subsidy:"
    di "  SE rate:     " %6.4f se_ms_`subsidy' " (change: " %+6.4f (se_ms_`subsidy' - baseline_se) ")"
    di "  Mobile rate: " %6.4f mobile_ms_`subsidy' " (change: " %+6.4f (mobile_ms_`subsidy' - baseline_mobile) ")"

    restore
}

/*******************************************************************************
* 4. Counterfactual 3: Broadband Infrastructure Investment
*
* Policy: Increase broadband penetration
* Simulate +1 SD, +2 SD, and universal (95th percentile) broadband
*******************************************************************************/

di _n "============================================================"
di "COUNTERFACTUAL 3: BROADBAND INFRASTRUCTURE INVESTMENT"
di "============================================================"

foreach bb_increase in 1 2 {

    preserve

    di _n "--- +`bb_increase' SD Broadband Expansion ---"

    * New broadband level
    gen bb_new = bb_std + `bb_increase'

    * Effect on choice probabilities via structural parameters
    * log(P_mobile/P_branch) increases by bb_mobile * delta_bb

    scalar delta_log_odds_mobile = bb_mobile * `bb_increase'
    scalar multiplier_mobile = exp(delta_log_odds_mobile)

    di "Log-odds shift (mobile vs branch): " %6.4f delta_log_odds_mobile
    di "Mobile probability multiplier:     " %6.4f multiplier_mobile

    * Shift probabilities
    * Mobile alternatives get multiplied, branch alternatives divided
    gen p4_new = p4 * multiplier_mobile
    gen p5_new = p5 * multiplier_mobile * exp(bb_mobile_se * `bb_increase')
    gen p6_new = p6 * multiplier_mobile

    gen p7_new = p7
    gen p8_new = p8 * exp(bb_branch_se * `bb_increase')
    gen p9_new = p9

    gen p1_new = p1
    gen p2_new = p2
    gen p3_new = p3

    * Normalize
    egen psum_new = rowtotal(p1_new p2_new p3_new p4_new p5_new p6_new p7_new p8_new p9_new)
    forvalues j = 1/9 {
        replace p`j'_new = p`j'_new / psum_new
    }

    * Calculate outcomes
    egen se_bb`bb_increase' = rowtotal(p2_new p5_new p8_new)
    egen mobile_bb`bb_increase' = rowtotal(p4_new p5_new p6_new)

    sum se_bb`bb_increase' [aw=n]
    scalar se_bb_`bb_increase' = r(mean)
    sum mobile_bb`bb_increase' [aw=n]
    scalar mobile_bb_`bb_increase' = r(mean)

    di "Results after +`bb_increase' SD broadband:"
    di "  SE rate:     " %6.4f se_bb_`bb_increase' " (change: " %+6.4f (se_bb_`bb_increase' - baseline_se) ")"
    di "  Mobile rate: " %6.4f mobile_bb_`bb_increase' " (change: " %+6.4f (mobile_bb_`bb_increase' - baseline_mobile) ")"

    restore
}

* Universal broadband (bring all areas to 95th percentile)
preserve

di _n "--- Universal Broadband (95th percentile) ---"

sum pct_broadband, detail
scalar bb_95 = r(p95)
di "95th percentile broadband: " %5.1f bb_95 "%"

* Only increase for areas below 95th percentile
gen bb_increase_needed = max(0, (bb_95 - pct_broadband) / bb_sd)
sum bb_increase_needed [aw=n]
di "Average increase needed (in SD): " %5.3f r(mean)

* Apply increase
gen p4_new = p4 * exp(bb_mobile * bb_increase_needed)
gen p5_new = p5 * exp(bb_mobile * bb_increase_needed) * exp(bb_mobile_se * bb_increase_needed)
gen p6_new = p6 * exp(bb_mobile * bb_increase_needed)

gen p7_new = p7
gen p8_new = p8 * exp(bb_branch_se * bb_increase_needed)
gen p9_new = p9

gen p1_new = p1
gen p2_new = p2
gen p3_new = p3

egen psum_new = rowtotal(p1_new p2_new p3_new p4_new p5_new p6_new p7_new p8_new p9_new)
forvalues j = 1/9 {
    replace p`j'_new = p`j'_new / psum_new
}

egen se_bb_univ = rowtotal(p2_new p5_new p8_new)
egen mobile_bb_univ = rowtotal(p4_new p5_new p6_new)

sum se_bb_univ [aw=n]
scalar se_bb_universal = r(mean)
sum mobile_bb_univ [aw=n]
scalar mobile_bb_universal = r(mean)

di "Results after universal broadband:"
di "  SE rate:     " %6.4f se_bb_universal " (change: " %+6.4f (se_bb_universal - baseline_se) ")"
di "  Mobile rate: " %6.4f mobile_bb_universal " (change: " %+6.4f (mobile_bb_universal - baseline_mobile) ")"

restore

/*******************************************************************************
* 5. Counterfactual 4: Combined Policy - Branch Closure + Broadband
*
* Realistic scenario: Branches close but broadband improves
*******************************************************************************/

di _n "============================================================"
di "COUNTERFACTUAL 4: COMBINED POLICY SCENARIOS"
di "============================================================"

preserve

di _n "--- 50% Branch Closure + Universal Broadband ---"

* First apply branch closure
local closure_frac = 0.50
local mobile_switch = 0.80
local unbank_switch = 0.20

gen displaced = `closure_frac' * (p7 + p8 + p9)
gen branch_se_share = p8 / (p7 + p8 + p9 + 0.0001)
gen branch_wage_share = p7 / (p7 + p8 + p9 + 0.0001)
gen branch_notwork_share = p9 / (p7 + p8 + p9 + 0.0001)

gen p7_temp = p7 * (1 - `closure_frac')
gen p8_temp = p8 * (1 - `closure_frac')
gen p9_temp = p9 * (1 - `closure_frac')

gen mobile_inflow = displaced * `mobile_switch'
gen p4_temp = p4 + mobile_inflow * branch_wage_share
gen p5_temp = p5 + mobile_inflow * branch_se_share * se_adjustment_mobile
gen p6_temp = p6 + mobile_inflow * branch_notwork_share

gen unbank_inflow = displaced * `unbank_switch'
gen p1_temp = p1 + unbank_inflow * branch_wage_share
gen p2_temp = p2 + unbank_inflow * branch_se_share * 0.5
gen p3_temp = p3 + unbank_inflow * branch_notwork_share

* Then apply broadband improvement
sum pct_broadband, detail
gen bb_increase_needed = max(0, (r(p95) - pct_broadband) / bb_sd)

gen p4_new = p4_temp * exp(bb_mobile * bb_increase_needed)
gen p5_new = p5_temp * exp(bb_mobile * bb_increase_needed) * exp(bb_mobile_se * bb_increase_needed)
gen p6_new = p6_temp * exp(bb_mobile * bb_increase_needed)

gen p7_new = p7_temp
gen p8_new = p8_temp * exp(bb_branch_se * bb_increase_needed)
gen p9_new = p9_temp

gen p1_new = p1_temp
gen p2_new = p2_temp
gen p3_new = p3_temp

egen psum_new = rowtotal(p1_new p2_new p3_new p4_new p5_new p6_new p7_new p8_new p9_new)
forvalues j = 1/9 {
    replace p`j'_new = p`j'_new / psum_new
}

egen se_combined = rowtotal(p2_new p5_new p8_new)
egen mobile_combined = rowtotal(p4_new p5_new p6_new)
egen unbanked_combined = rowtotal(p1_new p2_new p3_new)

sum se_combined [aw=n]
scalar se_combined_50bb = r(mean)
sum mobile_combined [aw=n]
scalar mobile_combined_50bb = r(mean)
sum unbanked_combined [aw=n]
scalar unbanked_combined_50bb = r(mean)

di "Results (50% closure + universal broadband):"
di "  SE rate:       " %6.4f se_combined_50bb " (change: " %+6.4f (se_combined_50bb - baseline_se) ")"
di "  Mobile rate:   " %6.4f mobile_combined_50bb " (change: " %+6.4f (mobile_combined_50bb - baseline_mobile) ")"
di "  Unbanked rate: " %6.4f unbanked_combined_50bb " (change: " %+6.4f (unbanked_combined_50bb - baseline_unbanked) ")"

* Compare to branch closure alone
di _n "Comparison:"
di "  50% closure alone:      SE = " %6.4f se_50
di "  50% closure + broadband: SE = " %6.4f se_combined_50bb
di "  Broadband mitigation:    " %+6.4f (se_combined_50bb - se_50)

restore

/*******************************************************************************
* 6. Heterogeneity Analysis: Effects by Demographics
*******************************************************************************/

di _n "============================================================"
di "HETEROGENEITY: COUNTERFACTUAL EFFECTS BY DEMOGRAPHICS"
di "============================================================"

* By education
di _n "--- Effects by Education Level ---"

forvalues educ = 1/4 {
    preserve
    keep if peducgrp == `educ'

    * Baseline SE
    egen se_base_e`educ' = rowtotal(p2 p5 p8)
    sum se_base_e`educ' [aw=n]
    local base_e`educ' = r(mean)

    * 50% branch closure
    local closure_frac = 0.50
    gen displaced = `closure_frac' * (p7 + p8 + p9)
    gen branch_se_share = p8 / (p7 + p8 + p9 + 0.0001)

    gen p8_new = p8 * 0.50
    gen p5_new = p5 + displaced * 0.80 * branch_se_share * se_adjustment_mobile
    gen p2_new = p2 + displaced * 0.20 * branch_se_share * 0.5

    egen se_cf_e`educ' = rowtotal(p2_new p5_new p8_new)
    sum se_cf_e`educ' [aw=n]
    local cf_e`educ' = r(mean)

    local change_e`educ' = `cf_e`educ'' - `base_e`educ''

    restore
}

di "Education     Baseline SE    After 50% Closure    Change"
di "Less than HS  " %8.4f `base_e1' "       " %8.4f `cf_e1' "          " %+8.4f `change_e1'
di "High School   " %8.4f `base_e2' "       " %8.4f `cf_e2' "          " %+8.4f `change_e2'
di "Some College  " %8.4f `base_e3' "       " %8.4f `cf_e3' "          " %+8.4f `change_e3'
di "College+      " %8.4f `base_e4' "       " %8.4f `cf_e4' "          " %+8.4f `change_e4'

* By age
di _n "--- Effects by Age Group ---"

forvalues age = 1/3 {
    preserve
    keep if age_cat == `age'

    egen se_base_a`age' = rowtotal(p2 p5 p8)
    sum se_base_a`age' [aw=n]
    local base_a`age' = r(mean)

    local closure_frac = 0.50
    gen displaced = `closure_frac' * (p7 + p8 + p9)
    gen branch_se_share = p8 / (p7 + p8 + p9 + 0.0001)

    gen p8_new = p8 * 0.50
    gen p5_new = p5 + displaced * 0.80 * branch_se_share * se_adjustment_mobile
    gen p2_new = p2 + displaced * 0.20 * branch_se_share * 0.5

    egen se_cf_a`age' = rowtotal(p2_new p5_new p8_new)
    sum se_cf_a`age' [aw=n]
    local cf_a`age' = r(mean)

    local change_a`age' = `cf_a`age'' - `base_a`age''

    restore
}

di "Age Group     Baseline SE    After 50% Closure    Change"
di "18-29         " %8.4f `base_a1' "       " %8.4f `cf_a1' "          " %+8.4f `change_a1'
di "30-44         " %8.4f `base_a2' "       " %8.4f `cf_a2' "          " %+8.4f `change_a2'
di "45-64         " %8.4f `base_a3' "       " %8.4f `cf_a3' "          " %+8.4f `change_a3'

/*******************************************************************************
* 7. Summary Table
*******************************************************************************/

di _n "============================================================"
di "PHASE 3 SUMMARY: COUNTERFACTUAL RESULTS"
di "============================================================"

di _n "SCENARIO                                SE Rate    Change    % Change"
di "--------------------------------------------------------------------"
di "Baseline                                " %6.4f baseline_se "      --        --"
di ""
di "BRANCH CLOSURES:"
di "  25% closure                           " %6.4f se_25 "   " %+7.4f (se_25 - baseline_se) "   " %+5.2f 100*(se_25 - baseline_se)/baseline_se "%"
di "  50% closure                           " %6.4f se_50 "   " %+7.4f (se_50 - baseline_se) "   " %+5.2f 100*(se_50 - baseline_se)/baseline_se "%"
di "  75% closure                           " %6.4f se_75 "   " %+7.4f (se_75 - baseline_se) "   " %+5.2f 100*(se_75 - baseline_se)/baseline_se "%"
di ""
di "MOBILE SUBSIDIES:"
di "  10% increase                          " %6.4f se_ms_10 "   " %+7.4f (se_ms_10 - baseline_se) "   " %+5.2f 100*(se_ms_10 - baseline_se)/baseline_se "%"
di "  25% increase                          " %6.4f se_ms_25 "   " %+7.4f (se_ms_25 - baseline_se) "   " %+5.2f 100*(se_ms_25 - baseline_se)/baseline_se "%"
di "  50% increase                          " %6.4f se_ms_50 "   " %+7.4f (se_ms_50 - baseline_se) "   " %+5.2f 100*(se_ms_50 - baseline_se)/baseline_se "%"
di ""
di "BROADBAND INVESTMENT:"
di "  +1 SD                                 " %6.4f se_bb_1 "   " %+7.4f (se_bb_1 - baseline_se) "   " %+5.2f 100*(se_bb_1 - baseline_se)/baseline_se "%"
di "  +2 SD                                 " %6.4f se_bb_2 "   " %+7.4f (se_bb_2 - baseline_se) "   " %+5.2f 100*(se_bb_2 - baseline_se)/baseline_se "%"
di "  Universal (95th pctile)               " %6.4f se_bb_universal "   " %+7.4f (se_bb_universal - baseline_se) "   " %+5.2f 100*(se_bb_universal - baseline_se)/baseline_se "%"
di ""
di "COMBINED POLICY:"
di "  50% closure + universal broadband     " %6.4f se_combined_50bb "   " %+7.4f (se_combined_50bb - baseline_se) "   " %+5.2f 100*(se_combined_50bb - baseline_se)/baseline_se "%"
di "--------------------------------------------------------------------"

/*******************************************************************************
* 8. Save Results
*******************************************************************************/

* Create comprehensive results dataset
clear
set obs 12

gen scenario = ""
gen se_rate = .
gen se_change = .
gen se_pct_change = .
gen mobile_rate = .
gen policy_type = ""

* Baseline
replace scenario = "Baseline" in 1
replace se_rate = baseline_se in 1
replace se_change = 0 in 1
replace se_pct_change = 0 in 1
replace mobile_rate = baseline_mobile in 1
replace policy_type = "Baseline" in 1

* Branch closures
replace scenario = "25% Branch Closure" in 2
replace se_rate = se_25 in 2
replace se_change = se_25 - baseline_se in 2
replace se_pct_change = 100*(se_25 - baseline_se)/baseline_se in 2
replace policy_type = "Branch Closure" in 2

replace scenario = "50% Branch Closure" in 3
replace se_rate = se_50 in 3
replace se_change = se_50 - baseline_se in 3
replace se_pct_change = 100*(se_50 - baseline_se)/baseline_se in 3
replace policy_type = "Branch Closure" in 3

replace scenario = "75% Branch Closure" in 4
replace se_rate = se_75 in 4
replace se_change = se_75 - baseline_se in 4
replace se_pct_change = 100*(se_75 - baseline_se)/baseline_se in 4
replace policy_type = "Branch Closure" in 4

* Mobile subsidies
replace scenario = "10% Mobile Subsidy" in 5
replace se_rate = se_ms_10 in 5
replace se_change = se_ms_10 - baseline_se in 5
replace se_pct_change = 100*(se_ms_10 - baseline_se)/baseline_se in 5
replace policy_type = "Mobile Subsidy" in 5

replace scenario = "25% Mobile Subsidy" in 6
replace se_rate = se_ms_25 in 6
replace se_change = se_ms_25 - baseline_se in 6
replace se_pct_change = 100*(se_ms_25 - baseline_se)/baseline_se in 6
replace policy_type = "Mobile Subsidy" in 6

replace scenario = "50% Mobile Subsidy" in 7
replace se_rate = se_ms_50 in 7
replace se_change = se_ms_50 - baseline_se in 7
replace se_pct_change = 100*(se_ms_50 - baseline_se)/baseline_se in 7
replace policy_type = "Mobile Subsidy" in 7

* Broadband
replace scenario = "+1 SD Broadband" in 8
replace se_rate = se_bb_1 in 8
replace se_change = se_bb_1 - baseline_se in 8
replace se_pct_change = 100*(se_bb_1 - baseline_se)/baseline_se in 8
replace policy_type = "Broadband" in 8

replace scenario = "+2 SD Broadband" in 9
replace se_rate = se_bb_2 in 9
replace se_change = se_bb_2 - baseline_se in 9
replace se_pct_change = 100*(se_bb_2 - baseline_se)/baseline_se in 9
replace policy_type = "Broadband" in 9

replace scenario = "Universal Broadband" in 10
replace se_rate = se_bb_universal in 10
replace se_change = se_bb_universal - baseline_se in 10
replace se_pct_change = 100*(se_bb_universal - baseline_se)/baseline_se in 10
replace policy_type = "Broadband" in 10

* Combined
replace scenario = "50% Closure + Universal BB" in 11
replace se_rate = se_combined_50bb in 11
replace se_change = se_combined_50bb - baseline_se in 11
replace se_pct_change = 100*(se_combined_50bb - baseline_se)/baseline_se in 11
replace mobile_rate = mobile_combined_50bb in 11
replace policy_type = "Combined" in 11

export delimited "$output/phase3_counterfactual_results.csv", replace

di _n "Results saved to $output/phase3_counterfactual_results.csv"

log close

di _n "Phase 3 counterfactual simulations complete!"
