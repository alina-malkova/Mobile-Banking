/*******************************************************************************
* Phase 3: Simplified Counterfactual Simulations
*******************************************************************************/

clear all
set more off

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase3_simple.log", replace

* Structural parameters
scalar bb_mobile = 0.1189
scalar bb_mobile_se = -0.2472
scalar bb_branch_se = -0.1133
scalar se_adj_mobile = exp(bb_mobile_se - bb_branch_se)

di "SE adjustment (branch→mobile): " se_adj_mobile

* Load CCPs
import delimited "$output/ccps_for_structural.csv", clear

rename prob_unbanked_wage p1
rename prob_unbanked_se p2
rename prob_unbanked_notwork p3
rename prob_mobile_wage p4
rename prob_mobile_se p5
rename prob_mobile_notwork p6
rename prob_branch_wage p7
rename prob_branch_se p8
rename prob_branch_notwork p9

egen psum = rowtotal(p1-p9)
forvalues j = 1/9 {
    replace p`j' = p`j' / psum
}
drop psum prob_sum

* Baseline
egen se_base = rowtotal(p2 p5 p8)
egen mobile_base = rowtotal(p4 p5 p6)
sum se_base [aw=n]
scalar baseline_se = r(mean)
sum mobile_base [aw=n]
scalar baseline_mobile = r(mean)

sum pct_broadband
scalar bb_mean = r(mean)
scalar bb_sd = r(sd)
gen bb_std = (pct_broadband - bb_mean) / bb_sd

di _n "BASELINE:"
di "  SE rate: " baseline_se
di "  Mobile rate: " baseline_mobile

/*******************************************************************************
* Branch Closure Simulations
*******************************************************************************/

di _n "=== BRANCH CLOSURES ==="

foreach closure in 25 50 75 {
    preserve

    local cf = `closure'/100

    gen displaced = `cf' * (p7 + p8 + p9)
    gen branch_se_sh = p8 / (p7 + p8 + p9 + 0.0001)
    gen branch_wage_sh = p7 / (p7 + p8 + p9 + 0.0001)
    gen branch_notwork_sh = p9 / (p7 + p8 + p9 + 0.0001)

    * 80% switch to mobile, 20% unbanked
    gen p7_n = p7 * (1 - `cf')
    gen p8_n = p8 * (1 - `cf')
    gen p9_n = p9 * (1 - `cf')

    gen p4_n = p4 + displaced * 0.80 * branch_wage_sh
    gen p5_n = p5 + displaced * 0.80 * branch_se_sh * se_adj_mobile
    gen p6_n = p6 + displaced * 0.80 * branch_notwork_sh

    gen p1_n = p1 + displaced * 0.20 * branch_wage_sh
    gen p2_n = p2 + displaced * 0.20 * branch_se_sh * 0.5
    gen p3_n = p3 + displaced * 0.20 * branch_notwork_sh

    egen ps = rowtotal(p1_n p2_n p3_n p4_n p5_n p6_n p7_n p8_n p9_n)
    forvalues j = 1/9 {
        replace p`j'_n = p`j'_n / ps
    }

    egen se_new = rowtotal(p2_n p5_n p8_n)
    sum se_new [aw=n]
    scalar se_`closure' = r(mean)

    local change = se_`closure' - baseline_se
    local pct = 100 * `change' / baseline_se

    di "`closure'% closure: SE = " %6.4f se_`closure' ", change = " %7.4f `change' " (" %5.2f `pct' "%)"

    restore
}

/*******************************************************************************
* Mobile Subsidy Simulations
*******************************************************************************/

di _n "=== MOBILE SUBSIDIES ==="

foreach sub in 10 25 50 {
    preserve

    local sf = `sub'/100

    gen mobile_curr = p4 + p5 + p6
    gen new_mobile = mobile_curr * `sf'

    gen branch_se_sh = p8 / (p7 + p8 + p9 + 0.0001)

    gen p4_n = p4 + new_mobile * 0.7 * 0.85
    gen p5_n = p5 + new_mobile * 0.7 * 0.10 * se_adj_mobile
    gen p6_n = p6 + new_mobile * 0.7 * 0.05

    gen p7_n = p7 - new_mobile * 0.7 * 0.85
    gen p8_n = max(p8 - new_mobile * 0.7 * 0.10, 0.001)
    gen p9_n = p9 - new_mobile * 0.7 * 0.05

    gen p1_n = max(p1 - new_mobile * 0.3 * 0.7, 0.001)
    gen p2_n = p2 + new_mobile * 0.3 * 0.2
    gen p3_n = max(p3 - new_mobile * 0.3 * 0.1, 0.001)

    egen ps = rowtotal(p1_n p2_n p3_n p4_n p5_n p6_n p7_n p8_n p9_n)
    forvalues j = 1/9 {
        replace p`j'_n = p`j'_n / ps
    }

    egen se_new = rowtotal(p2_n p5_n p8_n)
    sum se_new [aw=n]
    scalar se_ms_`sub' = r(mean)

    local change = se_ms_`sub' - baseline_se
    local pct = 100 * `change' / baseline_se

    di "`sub'% subsidy: SE = " %6.4f se_ms_`sub' ", change = " %7.4f `change' " (" %5.2f `pct' "%)"

    restore
}

/*******************************************************************************
* Broadband Investment Simulations
*******************************************************************************/

di _n "=== BROADBAND INVESTMENT ==="

foreach sd in 1 2 {
    preserve

    scalar mult = exp(bb_mobile * `sd')

    gen p4_n = p4 * mult
    gen p5_n = p5 * mult * exp(bb_mobile_se * `sd')
    gen p6_n = p6 * mult
    gen p7_n = p7
    gen p8_n = p8 * exp(bb_branch_se * `sd')
    gen p9_n = p9
    gen p1_n = p1
    gen p2_n = p2
    gen p3_n = p3

    egen ps = rowtotal(p1_n p2_n p3_n p4_n p5_n p6_n p7_n p8_n p9_n)
    forvalues j = 1/9 {
        replace p`j'_n = p`j'_n / ps
    }

    egen se_new = rowtotal(p2_n p5_n p8_n)
    egen mobile_new = rowtotal(p4_n p5_n p6_n)
    sum se_new [aw=n]
    scalar se_bb_`sd' = r(mean)
    sum mobile_new [aw=n]
    scalar mobile_bb_`sd' = r(mean)

    local change = se_bb_`sd' - baseline_se
    local pct = 100 * `change' / baseline_se
    local m_change = mobile_bb_`sd' - baseline_mobile

    di "+`sd' SD broadband: SE = " %6.4f se_bb_`sd' ", change = " %7.4f `change' " (" %5.2f `pct' "%)"
    di "                    Mobile = " %6.4f mobile_bb_`sd' ", change = " %7.4f `m_change'

    restore
}

/*******************************************************************************
* Combined: 50% Closure + Broadband
*******************************************************************************/

di _n "=== COMBINED: 50% CLOSURE + BROADBAND ==="

preserve

* First: 50% closure
gen displaced = 0.50 * (p7 + p8 + p9)
gen branch_se_sh = p8 / (p7 + p8 + p9 + 0.0001)
gen branch_wage_sh = p7 / (p7 + p8 + p9 + 0.0001)
gen branch_notwork_sh = p9 / (p7 + p8 + p9 + 0.0001)

gen p7_t = p7 * 0.50
gen p8_t = p8 * 0.50
gen p9_t = p9 * 0.50
gen p4_t = p4 + displaced * 0.80 * branch_wage_sh
gen p5_t = p5 + displaced * 0.80 * branch_se_sh * se_adj_mobile
gen p6_t = p6 + displaced * 0.80 * branch_notwork_sh
gen p1_t = p1 + displaced * 0.20 * branch_wage_sh
gen p2_t = p2 + displaced * 0.20 * branch_se_sh * 0.5
gen p3_t = p3 + displaced * 0.20 * branch_notwork_sh

* Then: +2 SD broadband
scalar mult = exp(bb_mobile * 2)
gen p4_n = p4_t * mult
gen p5_n = p5_t * mult * exp(bb_mobile_se * 2)
gen p6_n = p6_t * mult
gen p7_n = p7_t
gen p8_n = p8_t * exp(bb_branch_se * 2)
gen p9_n = p9_t
gen p1_n = p1_t
gen p2_n = p2_t
gen p3_n = p3_t

egen ps = rowtotal(p1_n p2_n p3_n p4_n p5_n p6_n p7_n p8_n p9_n)
forvalues j = 1/9 {
    replace p`j'_n = p`j'_n / ps
}

egen se_comb = rowtotal(p2_n p5_n p8_n)
egen mobile_comb = rowtotal(p4_n p5_n p6_n)
sum se_comb [aw=n]
scalar se_combined = r(mean)
sum mobile_comb [aw=n]
scalar mobile_combined = r(mean)

local change = se_combined - baseline_se
local pct = 100 * `change' / baseline_se

di "50% closure + 2SD broadband: SE = " %6.4f se_combined ", change = " %7.4f `change' " (" %5.2f `pct' "%)"
di "                             Mobile = " %6.4f mobile_combined

restore

/*******************************************************************************
* Summary Table
*******************************************************************************/

di _n "============================================================"
di "          PHASE 3 COUNTERFACTUAL SUMMARY"
di "============================================================"
di ""
di "Scenario                          SE Rate    Change    %Change"
di "------------------------------------------------------------"
di "Baseline                          " %6.4f baseline_se "       --        --"
di ""
di "Branch Closures:"
di "  25% closure                     " %6.4f se_25 "   " %7.4f (se_25-baseline_se) "   " %5.2f 100*(se_25-baseline_se)/baseline_se
di "  50% closure                     " %6.4f se_50 "   " %7.4f (se_50-baseline_se) "   " %5.2f 100*(se_50-baseline_se)/baseline_se
di "  75% closure                     " %6.4f se_75 "   " %7.4f (se_75-baseline_se) "   " %5.2f 100*(se_75-baseline_se)/baseline_se
di ""
di "Mobile Subsidies:"
di "  10% increase                    " %6.4f se_ms_10 "   " %7.4f (se_ms_10-baseline_se) "   " %5.2f 100*(se_ms_10-baseline_se)/baseline_se
di "  25% increase                    " %6.4f se_ms_25 "   " %7.4f (se_ms_25-baseline_se) "   " %5.2f 100*(se_ms_25-baseline_se)/baseline_se
di "  50% increase                    " %6.4f se_ms_50 "   " %7.4f (se_ms_50-baseline_se) "   " %5.2f 100*(se_ms_50-baseline_se)/baseline_se
di ""
di "Broadband Investment:"
di "  +1 SD                           " %6.4f se_bb_1 "   " %7.4f (se_bb_1-baseline_se) "   " %5.2f 100*(se_bb_1-baseline_se)/baseline_se
di "  +2 SD                           " %6.4f se_bb_2 "   " %7.4f (se_bb_2-baseline_se) "   " %5.2f 100*(se_bb_2-baseline_se)/baseline_se
di ""
di "Combined Policy:"
di "  50% closure + 2SD broadband     " %6.4f se_combined "   " %7.4f (se_combined-baseline_se) "   " %5.2f 100*(se_combined-baseline_se)/baseline_se
di "------------------------------------------------------------"

/*******************************************************************************
* Save Results
*******************************************************************************/

clear
set obs 10
gen scenario = ""
gen se_rate = .
gen se_change = .
gen pct_change = .

replace scenario = "Baseline" in 1
replace se_rate = baseline_se in 1
replace se_change = 0 in 1
replace pct_change = 0 in 1

replace scenario = "25% Branch Closure" in 2
replace se_rate = se_25 in 2
replace se_change = se_25 - baseline_se in 2
replace pct_change = 100*(se_25-baseline_se)/baseline_se in 2

replace scenario = "50% Branch Closure" in 3
replace se_rate = se_50 in 3
replace se_change = se_50 - baseline_se in 3
replace pct_change = 100*(se_50-baseline_se)/baseline_se in 3

replace scenario = "75% Branch Closure" in 4
replace se_rate = se_75 in 4
replace se_change = se_75 - baseline_se in 4
replace pct_change = 100*(se_75-baseline_se)/baseline_se in 4

replace scenario = "10% Mobile Subsidy" in 5
replace se_rate = se_ms_10 in 5
replace se_change = se_ms_10 - baseline_se in 5
replace pct_change = 100*(se_ms_10-baseline_se)/baseline_se in 5

replace scenario = "25% Mobile Subsidy" in 6
replace se_rate = se_ms_25 in 6
replace se_change = se_ms_25 - baseline_se in 6
replace pct_change = 100*(se_ms_25-baseline_se)/baseline_se in 6

replace scenario = "+1 SD Broadband" in 7
replace se_rate = se_bb_1 in 7
replace se_change = se_bb_1 - baseline_se in 7
replace pct_change = 100*(se_bb_1-baseline_se)/baseline_se in 7

replace scenario = "+2 SD Broadband" in 8
replace se_rate = se_bb_2 in 8
replace se_change = se_bb_2 - baseline_se in 8
replace pct_change = 100*(se_bb_2-baseline_se)/baseline_se in 8

replace scenario = "50% Closure + 2SD BB" in 9
replace se_rate = se_combined in 9
replace se_change = se_combined - baseline_se in 9
replace pct_change = 100*(se_combined-baseline_se)/baseline_se in 9

export delimited "$output/phase3_counterfactual_results.csv", replace

di _n "Results saved!"

log close
