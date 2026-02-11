/*******************************************************************************
* Phase 2: Unobserved Heterogeneity with Model Selection
* 
* Sequential estimation with K=2,3,4 types
* Model selection via BIC (Bayesian Information Criterion)
* 
* Following Arcidiacono-Miller (2011) methodology
*******************************************************************************/

clear all
set more off
set matsize 11000
set seed 20260211

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase2_type_selection.log", replace

di _n "============================================================"
di "MODEL SELECTION FOR UNOBSERVED HETEROGENEITY"
di "Sequential BIC Comparison: K = 2, 3, 4 Types"
di "============================================================"

/*******************************************************************************
* 1. Data Preparation (same as before)
*******************************************************************************/

use "$datadir/analysis_dataset_with_se.dta", clear

keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .
keep if banking_mode != .

capture drop emp_status
gen emp_status = .
replace emp_status = 1 if wage_worker == 1
replace emp_status = 2 if self_employed == 1
replace emp_status = 3 if employed != 1

capture drop bank_mode
gen bank_mode = banking_mode

capture drop joint_choice
gen joint_choice = (bank_mode - 1) * 3 + emp_status

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
* 2. Create Type Scores for Initialization
*******************************************************************************/

* Composite score based on entrepreneurial propensity
gen type_score = 0

* Age: older = more entrepreneurial experience
replace type_score = type_score + 2 if age_cat == 3
replace type_score = type_score + 1 if age_cat == 2

* Education: college = human capital for entrepreneurship
replace type_score = type_score + 2 if educ == 4
replace type_score = type_score + 1 if educ == 3

* Banking mode: branch = relationship borrower
replace type_score = type_score + 2 if bank_mode == 3
replace type_score = type_score - 1 if bank_mode == 1

* Current SE status
replace type_score = type_score + 4 if emp_status == 2

* Tech adoption proxy: mobile banking
replace type_score = type_score - 1 if bank_mode == 2

egen type_rank = rank(type_score), unique
sum type_rank
local max_rank = r(max)

/*******************************************************************************
* 3. Collapse to CCP Cells
*******************************************************************************/

forvalues j = 1/9 {
    gen choice`j' = (joint_choice == `j')
}

* Create type indicators for different K values
* K=2: below/above median
gen type2_1 = (type_rank <= `max_rank'/2)
gen type2_2 = (type_rank > `max_rank'/2)

* K=3: terciles
gen type3_1 = (type_rank <= `max_rank'/3)
gen type3_2 = (type_rank > `max_rank'/3 & type_rank <= 2*`max_rank'/3)
gen type3_3 = (type_rank > 2*`max_rank'/3)

* K=4: quartiles
gen type4_1 = (type_rank <= `max_rank'/4)
gen type4_2 = (type_rank > `max_rank'/4 & type_rank <= `max_rank'/2)
gen type4_3 = (type_rank > `max_rank'/2 & type_rank <= 3*`max_rank'/4)
gen type4_4 = (type_rank > 3*`max_rank'/4)

collapse (mean) p1=choice1 p2=choice2 p3=choice3 p4=choice4 p5=choice5 ///
         p6=choice6 p7=choice7 p8=choice8 p9=choice9 ///
         tau2_1=type2_1 tau2_2=type2_2 ///
         tau3_1=type3_1 tau3_2=type3_2 tau3_3=type3_3 ///
         tau4_1=type4_1 tau4_2=type4_2 tau4_3=type4_3 tau4_4=type4_4 ///
         bb_std ///
         (count) n=joint_choice ///
         [aw=hsupwgtk], by(cbsa year age_cat educ)

drop if n < 25

local N_cells = _N
di "Number of cells: `N_cells'"

save "$output/ccps_type_selection.dta", replace

/*******************************************************************************
* 4. Normalize CCPs and Create Estimation Variables
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

forvalues j = 1/9 {
    gen lnp`j' = ln(p`j')
}
forvalues j = 1/9 {
    gen y`j' = lnp`j' - lnp7
}

* Stack data
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

local N_obs = _N
di "Total observations (cells x alternatives): `N_obs'"

/*******************************************************************************
* 5. Estimate K=1 (Homogeneous) Model as Baseline
*******************************************************************************/

di _n "============================================================"
di "K = 1 (HOMOGENEOUS MODEL)"
di "============================================================"

quietly reg y a1 a2 a3 a4 a5 a6 a8 a9 ///
    age2_se age3_se educ4_se bb_mobile bb_se dynamic_se switch_cost ///
    [aw=n], vce(cluster cbsa)

local ll_1 = e(ll)
local k_params_1 = e(rank)
local bic_1 = -2 * `ll_1' + `k_params_1' * ln(`N_obs')
local aic_1 = -2 * `ll_1' + 2 * `k_params_1'

scalar k1_dynamic = _b[dynamic_se]
scalar k1_switch = _b[switch_cost]
scalar k1_bb_se = _b[bb_se]

di "Log-likelihood: " %12.2f `ll_1'
di "Parameters:     " `k_params_1'
di "BIC:            " %12.2f `bic_1'
di "AIC:            " %12.2f `aic_1'
di ""
di "Dynamic SE:     " %8.4f k1_dynamic
di "Switching cost: " %8.4f k1_switch
di "BB x SE:        " %8.4f k1_bb_se

estimates store k1_model

/*******************************************************************************
* 6. Estimate K=2 Model
*******************************************************************************/

di _n "============================================================"
di "K = 2 TYPES"
di "============================================================"

* Normalize type shares
gen tau_sum2 = tau2_1 + tau2_2
replace tau2_1 = tau2_1 / tau_sum2
replace tau2_2 = tau2_2 / tau_sum2
drop tau_sum2

* Type shares
preserve
collapse (mean) tau2_1 tau2_2 [aw=n]
scalar pi2_1 = tau2_1[1]
scalar pi2_2 = tau2_2[1]
restore

di "Type shares: " %5.1f 100*pi2_1 "% / " %5.1f 100*pi2_2 "%"

* Create type interactions
gen dyn_k2_1 = dynamic_se * tau2_1
gen dyn_k2_2 = dynamic_se * tau2_2
gen sw_k2_1 = switch_cost * tau2_1
gen sw_k2_2 = switch_cost * tau2_2
gen bbse_k2_1 = bb_se * tau2_1
gen bbse_k2_2 = bb_se * tau2_2

quietly reg y a1 a2 a3 a4 a5 a6 a8 a9 ///
    age2_se age3_se educ4_se bb_mobile ///
    dyn_k2_1 dyn_k2_2 sw_k2_1 sw_k2_2 bbse_k2_1 bbse_k2_2 ///
    [aw=n], vce(cluster cbsa)

local ll_2 = e(ll)
local k_params_2 = e(rank)
local bic_2 = -2 * `ll_2' + `k_params_2' * ln(`N_obs')
local aic_2 = -2 * `ll_2' + 2 * `k_params_2'

scalar k2_t1_dynamic = _b[dyn_k2_1]
scalar k2_t2_dynamic = _b[dyn_k2_2]
scalar k2_t1_switch = _b[sw_k2_1]
scalar k2_t2_switch = _b[sw_k2_2]
scalar k2_t1_bbse = _b[bbse_k2_1]
scalar k2_t2_bbse = _b[bbse_k2_2]

di "Log-likelihood: " %12.2f `ll_2'
di "Parameters:     " `k_params_2'
di "BIC:            " %12.2f `bic_2'
di "AIC:            " %12.2f `aic_2'
di ""
di "Type 1 - Dynamic SE: " %8.4f k2_t1_dynamic "  Switch: " %8.4f k2_t1_switch "  BB×SE: " %8.4f k2_t1_bbse
di "Type 2 - Dynamic SE: " %8.4f k2_t2_dynamic "  Switch: " %8.4f k2_t2_switch "  BB×SE: " %8.4f k2_t2_bbse

estimates store k2_model

/*******************************************************************************
* 7. Estimate K=3 Model
*******************************************************************************/

di _n "============================================================"
di "K = 3 TYPES"
di "============================================================"

* Normalize type shares
gen tau_sum3 = tau3_1 + tau3_2 + tau3_3
replace tau3_1 = tau3_1 / tau_sum3
replace tau3_2 = tau3_2 / tau_sum3
replace tau3_3 = tau3_3 / tau_sum3
drop tau_sum3

preserve
collapse (mean) tau3_1 tau3_2 tau3_3 [aw=n]
scalar pi3_1 = tau3_1[1]
scalar pi3_2 = tau3_2[1]
scalar pi3_3 = tau3_3[1]
restore

di "Type shares: " %5.1f 100*pi3_1 "% / " %5.1f 100*pi3_2 "% / " %5.1f 100*pi3_3 "%"

gen dyn_k3_1 = dynamic_se * tau3_1
gen dyn_k3_2 = dynamic_se * tau3_2
gen dyn_k3_3 = dynamic_se * tau3_3
gen sw_k3_1 = switch_cost * tau3_1
gen sw_k3_2 = switch_cost * tau3_2
gen sw_k3_3 = switch_cost * tau3_3
gen bbse_k3_1 = bb_se * tau3_1
gen bbse_k3_2 = bb_se * tau3_2
gen bbse_k3_3 = bb_se * tau3_3

quietly reg y a1 a2 a3 a4 a5 a6 a8 a9 ///
    age2_se age3_se educ4_se bb_mobile ///
    dyn_k3_1 dyn_k3_2 dyn_k3_3 sw_k3_1 sw_k3_2 sw_k3_3 bbse_k3_1 bbse_k3_2 bbse_k3_3 ///
    [aw=n], vce(cluster cbsa)

local ll_3 = e(ll)
local k_params_3 = e(rank)
local bic_3 = -2 * `ll_3' + `k_params_3' * ln(`N_obs')
local aic_3 = -2 * `ll_3' + 2 * `k_params_3'

scalar k3_t1_dynamic = _b[dyn_k3_1]
scalar k3_t2_dynamic = _b[dyn_k3_2]
scalar k3_t3_dynamic = _b[dyn_k3_3]
scalar k3_t1_switch = _b[sw_k3_1]
scalar k3_t2_switch = _b[sw_k3_2]
scalar k3_t3_switch = _b[sw_k3_3]
scalar k3_t1_bbse = _b[bbse_k3_1]
scalar k3_t2_bbse = _b[bbse_k3_2]
scalar k3_t3_bbse = _b[bbse_k3_3]

di "Log-likelihood: " %12.2f `ll_3'
di "Parameters:     " `k_params_3'
di "BIC:            " %12.2f `bic_3'
di "AIC:            " %12.2f `aic_3'
di ""
di "Type 1 - Dynamic SE: " %8.4f k3_t1_dynamic "  Switch: " %8.4f k3_t1_switch "  BB×SE: " %8.4f k3_t1_bbse
di "Type 2 - Dynamic SE: " %8.4f k3_t2_dynamic "  Switch: " %8.4f k3_t2_switch "  BB×SE: " %8.4f k3_t2_bbse
di "Type 3 - Dynamic SE: " %8.4f k3_t3_dynamic "  Switch: " %8.4f k3_t3_switch "  BB×SE: " %8.4f k3_t3_bbse

estimates store k3_model

/*******************************************************************************
* 8. Estimate K=4 Model
*******************************************************************************/

di _n "============================================================"
di "K = 4 TYPES"
di "============================================================"

gen tau_sum4 = tau4_1 + tau4_2 + tau4_3 + tau4_4
replace tau4_1 = tau4_1 / tau_sum4
replace tau4_2 = tau4_2 / tau_sum4
replace tau4_3 = tau4_3 / tau_sum4
replace tau4_4 = tau4_4 / tau_sum4
drop tau_sum4

preserve
collapse (mean) tau4_1 tau4_2 tau4_3 tau4_4 [aw=n]
scalar pi4_1 = tau4_1[1]
scalar pi4_2 = tau4_2[1]
scalar pi4_3 = tau4_3[1]
scalar pi4_4 = tau4_4[1]
restore

di "Type shares: " %5.1f 100*pi4_1 "% / " %5.1f 100*pi4_2 "% / " %5.1f 100*pi4_3 "% / " %5.1f 100*pi4_4 "%"

gen dyn_k4_1 = dynamic_se * tau4_1
gen dyn_k4_2 = dynamic_se * tau4_2
gen dyn_k4_3 = dynamic_se * tau4_3
gen dyn_k4_4 = dynamic_se * tau4_4
gen sw_k4_1 = switch_cost * tau4_1
gen sw_k4_2 = switch_cost * tau4_2
gen sw_k4_3 = switch_cost * tau4_3
gen sw_k4_4 = switch_cost * tau4_4
gen bbse_k4_1 = bb_se * tau4_1
gen bbse_k4_2 = bb_se * tau4_2
gen bbse_k4_3 = bb_se * tau4_3
gen bbse_k4_4 = bb_se * tau4_4

quietly reg y a1 a2 a3 a4 a5 a6 a8 a9 ///
    age2_se age3_se educ4_se bb_mobile ///
    dyn_k4_1 dyn_k4_2 dyn_k4_3 dyn_k4_4 ///
    sw_k4_1 sw_k4_2 sw_k4_3 sw_k4_4 ///
    bbse_k4_1 bbse_k4_2 bbse_k4_3 bbse_k4_4 ///
    [aw=n], vce(cluster cbsa)

local ll_4 = e(ll)
local k_params_4 = e(rank)
local bic_4 = -2 * `ll_4' + `k_params_4' * ln(`N_obs')
local aic_4 = -2 * `ll_4' + 2 * `k_params_4'

di "Log-likelihood: " %12.2f `ll_4'
di "Parameters:     " `k_params_4'
di "BIC:            " %12.2f `bic_4'
di "AIC:            " %12.2f `aic_4'

estimates store k4_model

/*******************************************************************************
* 9. Model Selection Summary
*******************************************************************************/

di _n "============================================================"
di "MODEL SELECTION SUMMARY"
di "============================================================"

di _n "         K    Log-Lik    Params        BIC        AIC"
di "        ---  ---------  --------  ---------  ---------"
di "         1  " %10.2f `ll_1' "  " %8.0f `k_params_1' "  " %10.2f `bic_1' "  " %10.2f `aic_1'
di "         2  " %10.2f `ll_2' "  " %8.0f `k_params_2' "  " %10.2f `bic_2' "  " %10.2f `aic_2'
di "         3  " %10.2f `ll_3' "  " %8.0f `k_params_3' "  " %10.2f `bic_3' "  " %10.2f `aic_3'
di "         4  " %10.2f `ll_4' "  " %8.0f `k_params_4' "  " %10.2f `bic_4' "  " %10.2f `aic_4'

* Find minimum BIC
scalar min_bic = min(`bic_1', `bic_2', `bic_3', `bic_4')
local optimal_k = 1
if `bic_2' == min_bic local optimal_k = 2
if `bic_3' == min_bic local optimal_k = 3
if `bic_4' == min_bic local optimal_k = 4

di _n "*** BIC selects K = `optimal_k' types ***"

/*******************************************************************************
* 10. Counterfactual Comparison Across K
*******************************************************************************/

di _n "============================================================"
di "COUNTERFACTUAL COMPARISON: 50% BRANCH CLOSURE"
di "============================================================"

use "$output/ccps_type_selection.dta", clear

egen se_base = rowtotal(p2 p5 p8)
sum se_base [aw=n]
scalar baseline_se = r(mean)
di "Baseline SE rate: " %6.4f baseline_se

* K=1 counterfactual
preserve
scalar se_adj_1 = exp(k1_bb_se) * exp(k1_dynamic * (-0.5))
gen branch_total = p7 + p8 + p9
gen displaced = 0.50 * branch_total
gen branch_se_sh = p8 / (branch_total + 0.001)
gen p8_cf = p8 * 0.50 + displaced * 0.80 * branch_se_sh * se_adj_1
gen p5_cf = p5 + displaced * 0.80 * branch_se_sh * se_adj_1 * 0.5
gen p2_cf = p2 + displaced * 0.20 * branch_se_sh * 0.3
gen se_cf = p2_cf + p5_cf + p8_cf
sum se_cf [aw=n]
scalar cf_k1 = r(mean)
restore

local pct_k1 = 100 * (cf_k1 - baseline_se) / baseline_se
di "K=1: SE = " %6.4f cf_k1 " (" %5.1f `pct_k1' "%)"

* K=2 counterfactual (weighted average)
preserve
scalar se_adj_2_1 = exp(k2_t1_bbse) * exp(k2_t1_dynamic * (-0.5))
scalar se_adj_2_2 = exp(k2_t2_bbse) * exp(k2_t2_dynamic * (-0.5))
scalar se_adj_2 = pi2_1 * se_adj_2_1 + pi2_2 * se_adj_2_2
gen branch_total = p7 + p8 + p9
gen displaced = 0.50 * branch_total
gen branch_se_sh = p8 / (branch_total + 0.001)
gen p8_cf = p8 * 0.50 + displaced * 0.80 * branch_se_sh * se_adj_2
gen p5_cf = p5 + displaced * 0.80 * branch_se_sh * se_adj_2 * 0.5
gen p2_cf = p2 + displaced * 0.20 * branch_se_sh * 0.3
gen se_cf = p2_cf + p5_cf + p8_cf
sum se_cf [aw=n]
scalar cf_k2 = r(mean)
restore

local pct_k2 = 100 * (cf_k2 - baseline_se) / baseline_se
di "K=2: SE = " %6.4f cf_k2 " (" %5.1f `pct_k2' "%)"

* K=3 counterfactual
preserve
scalar se_adj_3_1 = exp(k3_t1_bbse) * exp(k3_t1_dynamic * (-0.5))
scalar se_adj_3_2 = exp(k3_t2_bbse) * exp(k3_t2_dynamic * (-0.5))
scalar se_adj_3_3 = exp(k3_t3_bbse) * exp(k3_t3_dynamic * (-0.5))
scalar se_adj_3 = pi3_1 * se_adj_3_1 + pi3_2 * se_adj_3_2 + pi3_3 * se_adj_3_3
gen branch_total = p7 + p8 + p9
gen displaced = 0.50 * branch_total
gen branch_se_sh = p8 / (branch_total + 0.001)
gen p8_cf = p8 * 0.50 + displaced * 0.80 * branch_se_sh * se_adj_3
gen p5_cf = p5 + displaced * 0.80 * branch_se_sh * se_adj_3 * 0.5
gen p2_cf = p2 + displaced * 0.20 * branch_se_sh * 0.3
gen se_cf = p2_cf + p5_cf + p8_cf
sum se_cf [aw=n]
scalar cf_k3 = r(mean)
restore

local pct_k3 = 100 * (cf_k3 - baseline_se) / baseline_se
di "K=3: SE = " %6.4f cf_k3 " (" %5.1f `pct_k3' "%)"

di _n "Counterfactual sensitivity to K:"
di "  K=1 vs K=2 difference: " %5.2f abs(`pct_k1' - `pct_k2') " pp"
di "  K=2 vs K=3 difference: " %5.2f abs(`pct_k2' - `pct_k3') " pp"

/*******************************************************************************
* 11. Final Results with Optimal K
*******************************************************************************/

di _n "============================================================"
di "FINAL RESULTS: K = `optimal_k' TYPES (BIC-SELECTED)"
di "============================================================"

if `optimal_k' == 2 {
    di _n "Type Interpretation:"
    di "  Type 1 (Relationship Borrowers): " %5.1f 100*pi2_1 "%"
    di "    - High returns to branch access"
    di "    - Dynamic SE: " %6.3f k2_t1_dynamic
    di "    - Switch cost: " %6.3f k2_t1_switch
    di "    - BB×SE: " %6.3f k2_t1_bbse
    di ""
    di "  Type 2 (Tech-Savvy/Adaptable): " %5.1f 100*pi2_2 "%"
    di "    - Can substitute mobile for branch"
    di "    - Dynamic SE: " %6.3f k2_t2_dynamic
    di "    - Switch cost: " %6.3f k2_t2_switch
    di "    - BB×SE: " %6.3f k2_t2_bbse
    
    scalar final_cf = cf_k2
    scalar final_pi1 = pi2_1
    scalar final_pi2 = pi2_2
    scalar final_t1_dyn = k2_t1_dynamic
    scalar final_t2_dyn = k2_t2_dynamic
    scalar final_t1_sw = k2_t1_switch
    scalar final_t2_sw = k2_t2_switch
    scalar final_t1_bbse = k2_t1_bbse
    scalar final_t2_bbse = k2_t2_bbse
}
else if `optimal_k' == 3 {
    di _n "Type Interpretation:"
    di "  Type 1 (Relationship Borrowers): " %5.1f 100*pi3_1 "%"
    di "    - Strong preference for branch banking"
    di "    - Dynamic SE: " %6.3f k3_t1_dynamic
    di "    - Switch cost: " %6.3f k3_t1_switch
    di "    - BB×SE: " %6.3f k3_t1_bbse
    di ""
    di "  Type 2 (Tech-Savvy Entrepreneurs): " %5.1f 100*pi3_2 "%"
    di "    - Can substitute mobile for branch"
    di "    - Dynamic SE: " %6.3f k3_t2_dynamic
    di "    - Switch cost: " %6.3f k3_t2_switch
    di "    - BB×SE: " %6.3f k3_t2_bbse
    di ""
    di "  Type 3 (Credit-Constrained): " %5.1f 100*pi3_3 "%"
    di "    - Low entrepreneurial ability regardless of banking"
    di "    - Dynamic SE: " %6.3f k3_t3_dynamic
    di "    - Switch cost: " %6.3f k3_t3_switch
    di "    - BB×SE: " %6.3f k3_t3_bbse
    
    scalar final_cf = cf_k3
}

local final_pct = 100 * (final_cf - baseline_se) / baseline_se
di _n "Counterfactual (50% Branch Closure):"
di "  Baseline SE:     " %6.4f baseline_se
di "  CF SE (K=`optimal_k'):   " %6.4f final_cf " (" %5.1f `final_pct' "%)"

/*******************************************************************************
* 12. Save Results
*******************************************************************************/

clear
set obs 15

gen str20 item = ""
gen value = .
gen str10 model = ""

* BIC comparison
replace item = "BIC" in 1/4
replace model = "K=1" in 1
replace model = "K=2" in 2
replace model = "K=3" in 3
replace model = "K=4" in 4
replace value = `bic_1' in 1
replace value = `bic_2' in 2
replace value = `bic_3' in 3
replace value = `bic_4' in 4

* Counterfactual comparison
replace item = "CF_SE_pct" in 5/7
replace model = "K=1" in 5
replace model = "K=2" in 6
replace model = "K=3" in 7
replace value = 100*(cf_k1 - baseline_se)/baseline_se in 5
replace value = 100*(cf_k2 - baseline_se)/baseline_se in 6
replace value = 100*(cf_k3 - baseline_se)/baseline_se in 7

* Optimal K results
replace item = "optimal_K" in 8
replace value = `optimal_k' in 8
replace model = "selected" in 8

replace item = "baseline_se" in 9
replace value = baseline_se in 9
replace model = "all" in 9

drop if item == ""

export delimited using "$output/phase2_type_selection.csv", replace

di _n "============================================================"
di "SUMMARY"
di "============================================================"
di "BIC-selected number of types: K = `optimal_k'"
di "Counterfactual effect (50% closure): " %5.1f `final_pct' "%"
di ""
di "Results saved to: $output/phase2_type_selection.csv"

log close
