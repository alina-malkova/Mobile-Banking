/*******************************************************************************
* Phase 2: Dynamic Discrete Choice with Unobserved Heterogeneity (v3)
* 
* Arcidiacono-Miller (2011) CCP Estimation with Finite Mixture
*******************************************************************************/

clear all
set more off
set matsize 11000
set seed 20260211

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase2_unobs_het.log", replace

di _n "============================================================"
di "PHASE 2: DYNAMIC DDC WITH UNOBSERVED HETEROGENEITY"
di "Arcidiacono-Miller (2011) Finite Mixture Approach"
di "============================================================"

/*******************************************************************************
* 1. Load and Prepare Data
*******************************************************************************/

use "$datadir/analysis_dataset_with_se.dta", clear

keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .
keep if banking_mode != .

* Employment and banking
capture drop emp_status
gen emp_status = .
replace emp_status = 1 if wage_worker == 1
replace emp_status = 2 if self_employed == 1
replace emp_status = 3 if employed != 1

capture drop bank_mode
gen bank_mode = banking_mode

capture drop joint_choice
gen joint_choice = (bank_mode - 1) * 3 + emp_status

* State variables
gen age_cat = 1 if age >= 18 & age < 30
replace age_cat = 2 if age >= 30 & age < 45
replace age_cat = 3 if age >= 45 & age <= 64

sum pct_broadband [aw=hsupwgtk]
gen bb_std = (pct_broadband - r(mean)) / r(sd)

gen educ = .
replace educ = 1 if no_hs == 1
replace educ = 2 if hs_diploma == 1
replace educ = 3 if some_college == 1
replace educ = 4 if college_degree == 1

/*******************************************************************************
* 2. Create Type Assignment Based on Observables
*******************************************************************************/

di _n "=== Initializing Types Based on Observables ==="

* Type 1: Entrepreneurial - older, educated, already SE
* Type 2: Risk-Averse - young, lower education  
* Type 3: Credit-Constrained - unbanked or mobile-only

gen type_score = 0
replace type_score = type_score + 2 if age_cat == 3
replace type_score = type_score + 1 if age_cat == 2
replace type_score = type_score - 1 if age_cat == 1
replace type_score = type_score + 2 if educ == 4
replace type_score = type_score + 1 if educ == 3
replace type_score = type_score - 1 if educ <= 2
replace type_score = type_score + 1 if bank_mode == 3
replace type_score = type_score - 1 if bank_mode == 1
replace type_score = type_score + 3 if emp_status == 2

egen type_rank = rank(type_score), unique
sum type_rank
gen init_type = 1 if type_rank <= r(max)/3
replace init_type = 2 if type_rank > r(max)/3 & type_rank <= 2*r(max)/3
replace init_type = 3 if type_rank > 2*r(max)/3

tab init_type

/*******************************************************************************
* 3. Estimate CCPs by Cell
*******************************************************************************/

di _n "=== Step 1: CCP Estimation ==="

forvalues j = 1/9 {
    gen choice`j' = (joint_choice == `j')
}

forvalues k = 1/3 {
    gen type`k' = (init_type == `k')
}

collapse (mean) p1=choice1 p2=choice2 p3=choice3 p4=choice4 p5=choice5 ///
         p6=choice6 p7=choice7 p8=choice8 p9=choice9 ///
         tau1_init=type1 tau2_init=type2 tau3_init=type3 ///
         bb_std ///
         (count) n=joint_choice ///
         [aw=hsupwgtk], by(cbsa year age_cat educ)

drop if n < 25
di "Cell-level observations: " _N

save "$output/ccps_for_em.dta", replace

/*******************************************************************************
* 4. Normalize CCPs and Compute Log-Odds
*******************************************************************************/

scalar eps = 0.001
forvalues j = 1/9 {
    replace p`j' = max(p`j', eps)
}

gen psum = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9
forvalues j = 1/9 {
    replace p`j' = p`j' / psum
}
drop psum

* First create all log-CCPs
forvalues j = 1/9 {
    gen lnp`j' = ln(p`j')
}

* Then create Berry inversions
forvalues j = 1/9 {
    gen y`j' = lnp`j' - lnp7
}

/*******************************************************************************
* 5. Stack Data for Estimation
*******************************************************************************/

gen obs = _n
expand 9
bysort obs: gen alt = _n

gen y = .
forvalues j = 1/9 {
    replace y = y`j' if alt == `j'
}

forvalues j = 1/9 {
    gen a`j' = (alt == `j')
}

gen is_mobile = inlist(alt, 4, 5, 6)
gen is_branch = inlist(alt, 7, 8, 9)
gen is_unbanked = inlist(alt, 1, 2, 3)
gen is_se = inlist(alt, 2, 5, 8)
gen is_wage = inlist(alt, 1, 4, 7)
gen is_notwork = inlist(alt, 3, 6, 9)

gen age2 = (age_cat == 2)
gen age3 = (age_cat == 3)
gen educ4 = (educ == 4)
gen age2_se = age2 * is_se
gen age3_se = age3 * is_se
gen educ4_se = educ4 * is_se
gen bb_mobile = bb_std * is_mobile
gen bb_se = bb_std * is_se

scalar beta = 0.90
gen dynamic_se = is_se * age_cat * beta
gen switch_cost = (is_mobile + is_unbanked)

/*******************************************************************************
* 6. Estimate Type-Specific Models Using Interactions
*******************************************************************************/

di _n "============================================================"
di "ESTIMATING TYPE-SPECIFIC MODELS"
di "============================================================"

* Normalize type shares
gen tau_sum = tau1_init + tau2_init + tau3_init
gen tau1 = tau1_init / tau_sum
gen tau2 = tau2_init / tau_sum
gen tau3 = tau3_init / tau_sum
drop tau_sum

di _n "Type shares in estimation sample:"
sum tau1 tau2 tau3 [aw=n]

* Create type-weighted interactions
gen dyn_t1 = dynamic_se * tau1
gen dyn_t2 = dynamic_se * tau2
gen dyn_t3 = dynamic_se * tau3

gen sw_t1 = switch_cost * tau1
gen sw_t2 = switch_cost * tau2
gen sw_t3 = switch_cost * tau3

gen bbse_t1 = bb_se * tau1
gen bbse_t2 = bb_se * tau2
gen bbse_t3 = bb_se * tau3

gen bbmob_t1 = bb_mobile * tau1
gen bbmob_t2 = bb_mobile * tau2
gen bbmob_t3 = bb_mobile * tau3

* SE preference by type
gen se_t1 = is_se * tau1
gen se_t2 = is_se * tau2
gen se_t3 = is_se * tau3

di _n "=== Pooled Model with Type-Specific Parameters ==="

reg y a1 a2 a3 a4 a5 a6 a8 a9 ///
    age2_se age3_se educ4_se ///
    dyn_t1 dyn_t2 dyn_t3 ///
    sw_t1 sw_t2 sw_t3 ///
    bbse_t1 bbse_t2 bbse_t3 ///
    bbmob_t1 bbmob_t2 bbmob_t3 ///
    [aw=n], vce(cluster cbsa)

estimates store types_pooled

/*******************************************************************************
* 7. Extract Type-Specific Parameters
*******************************************************************************/

di _n "============================================================"
di "TYPE-SPECIFIC PARAMETER ESTIMATES"
di "============================================================"

scalar type1_dynamic = _b[dyn_t1]
scalar type1_switch = _b[sw_t1]
scalar type1_bb_se = _b[bbse_t1]
scalar type1_bb_mobile = _b[bbmob_t1]

di _n "=== Type 1: Entrepreneurial ==="
di "  Dynamic SE return: " %8.4f type1_dynamic
di "  Switching cost:    " %8.4f type1_switch
di "  BB x SE:           " %8.4f type1_bb_se
di "  BB x Mobile:       " %8.4f type1_bb_mobile

scalar type2_dynamic = _b[dyn_t2]
scalar type2_switch = _b[sw_t2]
scalar type2_bb_se = _b[bbse_t2]
scalar type2_bb_mobile = _b[bbmob_t2]

di _n "=== Type 2: Risk-Averse ==="
di "  Dynamic SE return: " %8.4f type2_dynamic
di "  Switching cost:    " %8.4f type2_switch
di "  BB x SE:           " %8.4f type2_bb_se
di "  BB x Mobile:       " %8.4f type2_bb_mobile

scalar type3_dynamic = _b[dyn_t3]
scalar type3_switch = _b[sw_t3]
scalar type3_bb_se = _b[bbse_t3]
scalar type3_bb_mobile = _b[bbmob_t3]

di _n "=== Type 3: Credit-Constrained ==="
di "  Dynamic SE return: " %8.4f type3_dynamic
di "  Switching cost:    " %8.4f type3_switch
di "  BB x SE:           " %8.4f type3_bb_se
di "  BB x Mobile:       " %8.4f type3_bb_mobile

/*******************************************************************************
* 8. Compute Type Shares
*******************************************************************************/

preserve
collapse (mean) tau1 tau2 tau3 [aw=n]
scalar pi1 = tau1[1]
scalar pi2 = tau2[1]
scalar pi3 = tau3[1]
restore

di _n "Type Shares:"
di "  Type 1 (Entrepreneurial):    " %5.1f 100*pi1 "%"
di "  Type 2 (Risk-Averse):        " %5.1f 100*pi2 "%"
di "  Type 3 (Credit-Constrained): " %5.1f 100*pi3 "%"

/*******************************************************************************
* 9. Weighted Average Parameters
*******************************************************************************/

di _n "============================================================"
di "WEIGHTED AVERAGE PARAMETERS"
di "============================================================"

scalar avg_dynamic = pi1 * type1_dynamic + pi2 * type2_dynamic + pi3 * type3_dynamic
scalar avg_switch = pi1 * type1_switch + pi2 * type2_switch + pi3 * type3_switch
scalar avg_bb_se = pi1 * type1_bb_se + pi2 * type2_bb_se + pi3 * type3_bb_se
scalar avg_bb_mobile = pi1 * type1_bb_mobile + pi2 * type2_bb_mobile + pi3 * type3_bb_mobile

di "Dynamic SE return (weighted): " %8.4f avg_dynamic
di "Switching cost (weighted):    " %8.4f avg_switch
di "BB x SE (weighted):           " %8.4f avg_bb_se
di "BB x Mobile (weighted):       " %8.4f avg_bb_mobile

/*******************************************************************************
* 10. Type-Specific Counterfactuals
*******************************************************************************/

di _n "============================================================"
di "TYPE-SPECIFIC COUNTERFACTUALS: 50% BRANCH CLOSURE"
di "============================================================"

use "$output/ccps_for_em.dta", clear

egen se_base = rowtotal(p2 p5 p8)
sum se_base [aw=n]
scalar baseline_se = r(mean)
di "Baseline SE rate: " %6.4f baseline_se

foreach type in 1 2 3 {
    preserve
    
    if `type' == 1 {
        scalar t_bb_se = type1_bb_se
        scalar t_bb_mobile = type1_bb_mobile
        scalar t_dynamic = type1_dynamic
        scalar t_share = pi1
        local t_name "Entrepreneurial"
    }
    else if `type' == 2 {
        scalar t_bb_se = type2_bb_se
        scalar t_bb_mobile = type2_bb_mobile
        scalar t_dynamic = type2_dynamic
        scalar t_share = pi2
        local t_name "Risk-Averse"
    }
    else {
        scalar t_bb_se = type3_bb_se
        scalar t_bb_mobile = type3_bb_mobile
        scalar t_dynamic = type3_dynamic
        scalar t_share = pi3
        local t_name "Credit-Constrained"
    }
    
    scalar se_adj_static = exp(t_bb_se - t_bb_mobile)
    scalar se_adj_dynamic = exp(t_dynamic * (-0.5))
    scalar se_adj_total = se_adj_static * se_adj_dynamic
    
    gen branch_total = p7 + p8 + p9
    gen displaced = 0.50 * branch_total
    gen branch_se_sh = p8 / (branch_total + 0.001)
    gen branch_wage_sh = p7 / (branch_total + 0.001)
    gen branch_notwork_sh = p9 / (branch_total + 0.001)
    
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
    scalar cf_se_`type' = r(mean)
    
    local pct_change = 100 * (cf_se_`type' - baseline_se) / baseline_se
    
    di "Type `type' (`t_name', " %4.1f 100*t_share "%): SE " %6.4f cf_se_`type' " (" %5.1f `pct_change' "%)"
    
    restore
}

scalar cf_se_weighted = pi1 * cf_se_1 + pi2 * cf_se_2 + pi3 * cf_se_3
local pct_weighted = 100 * (cf_se_weighted - baseline_se) / baseline_se

di _n "Weighted average CF SE: " %6.4f cf_se_weighted " (" %5.1f `pct_weighted' "%)"

/*******************************************************************************
* 11. Summary
*******************************************************************************/

di _n "============================================================"
di "SUMMARY: UNOBSERVED HETEROGENEITY RESULTS"
di "============================================================"

di _n "Type Shares:"
di "  Type 1 (Entrepreneurial):     " %5.1f 100*pi1 "%"
di "  Type 2 (Risk-Averse):         " %5.1f 100*pi2 "%"
di "  Type 3 (Credit-Constrained):  " %5.1f 100*pi3 "%"

di _n "Key Parameters by Type:"
di "                        Type 1    Type 2    Type 3    Weighted"
di "  Dynamic SE return:   " %7.3f type1_dynamic "   " %7.3f type2_dynamic "   " %7.3f type3_dynamic "   " %7.3f avg_dynamic
di "  Switching cost:      " %7.3f type1_switch "   " %7.3f type2_switch "   " %7.3f type3_switch "   " %7.3f avg_switch
di "  BB x SE:             " %7.3f type1_bb_se "   " %7.3f type2_bb_se "   " %7.3f type3_bb_se "   " %7.3f avg_bb_se
di "  BB x Mobile:         " %7.3f type1_bb_mobile "   " %7.3f type2_bb_mobile "   " %7.3f type3_bb_mobile "   " %7.3f avg_bb_mobile

di _n "Counterfactual (50% Branch Closure):"
di "  Baseline SE:    " %6.4f baseline_se
di "  Type 1 CF SE:   " %6.4f cf_se_1 " (" %5.1f 100*(cf_se_1-baseline_se)/baseline_se "%)"
di "  Type 2 CF SE:   " %6.4f cf_se_2 " (" %5.1f 100*(cf_se_2-baseline_se)/baseline_se "%)"
di "  Type 3 CF SE:   " %6.4f cf_se_3 " (" %5.1f 100*(cf_se_3-baseline_se)/baseline_se "%)"
di "  Weighted CF SE: " %6.4f cf_se_weighted " (" %5.1f 100*(cf_se_weighted-baseline_se)/baseline_se "%)"

/*******************************************************************************
* 12. Save Results
*******************************************************************************/

clear
set obs 20

gen str20 parameter = ""
gen value = .
gen str15 type = ""

local row = 1

foreach t in 1 2 3 {
    replace parameter = "pi" in `row'
    replace type = "Type `t'" in `row'
    replace value = pi`t' in `row'
    local row = `row' + 1
}

foreach t in 1 2 3 {
    replace parameter = "gamma_dynamic" in `row'
    replace type = "Type `t'" in `row'
    replace value = type`t'_dynamic in `row'
    local row = `row' + 1
}

foreach t in 1 2 3 {
    replace parameter = "kappa_switch" in `row'
    replace type = "Type `t'" in `row'
    replace value = type`t'_switch in `row'
    local row = `row' + 1
}

foreach t in 1 2 3 {
    replace parameter = "gamma_bb_se" in `row'
    replace type = "Type `t'" in `row'
    replace value = type`t'_bb_se in `row'
    local row = `row' + 1
}

foreach t in 1 2 3 {
    replace parameter = "cf_se" in `row'
    replace type = "Type `t'" in `row'
    replace value = cf_se_`t' in `row'
    local row = `row' + 1
}

replace parameter = "baseline_se" in `row'
replace value = baseline_se in `row'
replace type = "All" in `row'
local row = `row' + 1

replace parameter = "cf_se_weighted" in `row'
replace value = cf_se_weighted in `row'
replace type = "All" in `row'

drop if parameter == ""

export delimited using "$output/phase2_unobs_het_params.csv", replace

log close

di _n "Phase 2 with unobserved heterogeneity complete."
