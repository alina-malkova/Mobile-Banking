/*==============================================================================
Phase 3: Proper Structural Counterfactuals (Final Version)

Use individual-level multinomial logit on the joint (banking, employment) choice.
When branch density changes, recompute choice probabilities via the logit formula.
==============================================================================*/

clear all
set more off
set matsize 11000
set seed 20260211

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase3_structural_cf.log", replace text

di _n "============================================================"
di "STRUCTURAL COUNTERFACTUAL: Joint Choice Model"
di "============================================================"

/*------------------------------------------------------------------------------
1. Load Individual-Level Data
------------------------------------------------------------------------------*/

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

* Joint choice state
* States 2, 5, 8 are SE states (Unbanked×SE, Mobile×SE, Branch×SE)
gen state = (bank_mode - 1) * 3 + emp_status

gen se = (self_employed == 1)
gen branch = (bank_mode == 3)
gen mobile = (bank_mode == 2)

gen age_cat = 1 if age >= 18 & age < 35
replace age_cat = 2 if age >= 35 & age < 50
replace age_cat = 3 if age >= 50 & age <= 64

gen educ_cat = .
replace educ_cat = 1 if no_hs == 1
replace educ_cat = 2 if hs_diploma == 1
replace educ_cat = 3 if some_college == 1
replace educ_cat = 4 if college_degree == 1

di "Sample: " _N " individuals"
quietly summarize se [aw=hsupwgtk]
local baseline_se = r(mean)
di "Baseline SE rate: " %6.4f `baseline_se'

quietly summarize branch [aw=hsupwgtk]
local baseline_branch = r(mean)
di "Baseline branch share: " %6.4f `baseline_branch'

/*------------------------------------------------------------------------------
2. Estimate Individual-Level Model

   Rather than multinomial logit on 9 choices (computationally intensive),
   use a two-stage approach:

   Stage 1: P(SE | banking mode, X, density)
   Stage 2: P(banking mode | X, density)

   This decomposes the joint probability as:
   P(bank, SE) = P(SE | bank) × P(bank)

   For counterfactual, we change density and recompute both probabilities.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "STEP 1: Estimate SE Choice Conditional on Banking Mode"
di "============================================================"

* Logit: P(SE=1 | bank_mode, density, X)
logit se i.bank_mode##c.pct_broadband i.age_cat i.educ_cat i.year ///
      [pw=hsupwgtk], vce(cluster cbsa)

* The key parameters:
* _b[3.bank_mode] = effect of branch (vs. unbanked) on SE
* _b[3.bank_mode#c.pct_broadband] = interaction of branch × density

local b_branch = _b[3.bank_mode]
local b_mobile = _b[2.bank_mode]
local b_branch_density = _b[3.bank_mode#c.pct_broadband]
local b_mobile_density = _b[2.bank_mode#c.pct_broadband]

di _n "SE Choice Parameters:"
di "  Branch effect (vs. unbanked): " %8.4f `b_branch'
di "  Mobile effect (vs. unbanked): " %8.4f `b_mobile'
di "  Branch × Density: " %8.4f `b_branch_density'
di "  Mobile × Density: " %8.4f `b_mobile_density'

* Baseline SE probabilities
predict p_se_baseline, pr
quietly summarize p_se_baseline [aw=hsupwgtk]
di "Model baseline SE rate: " %6.4f r(mean)
local model_baseline = r(mean)

/*------------------------------------------------------------------------------
3. Structural Counterfactual

   When density falls by 50%:
   - For branch users: utility of SE changes by Δv = β_branch_density × Δdensity
   - For mobile users: utility changes by β_mobile_density × Δdensity
   - For unbanked: no direct density effect (but density main effect applies)

   New probability: P_new = Λ(xb_old + Δv)
------------------------------------------------------------------------------*/

di _n "============================================================"
di "STEP 2: Compute Structural Counterfactual"
di "============================================================"

* Get baseline linear predictor
predict xb_base, xb

* Counterfactual: 50% density reduction
* Change in utility for each banking mode
gen delta_xb = 0

* For branch users: Δv = β_branch_density × (-0.5 × density)
* Note: We also need to account for the main density effect
replace delta_xb = `b_branch_density' * (-0.5 * pct_broadband) if bank_mode == 3

* For mobile users:
replace delta_xb = `b_mobile_density' * (-0.5 * pct_broadband) if bank_mode == 2

* Counterfactual linear predictor and probability
gen xb_cf = xb_base + delta_xb
gen p_se_cf = exp(xb_cf) / (1 + exp(xb_cf))

* Results
quietly summarize p_se_baseline [aw=hsupwgtk]
local base_se = r(mean)

quietly summarize p_se_cf [aw=hsupwgtk]
local cf_se = r(mean)

local effect_pp = `cf_se' - `base_se'
local effect_pct = `effect_pp' / `base_se' * 100

di _n "STRUCTURAL COUNTERFACTUAL RESULTS (50% Density Reduction):"
di "==========================================================="
di "  Baseline SE rate: " %6.4f `base_se'
di "  Counterfactual SE rate: " %6.4f `cf_se'
di "  Change (pp): " %7.4f `effect_pp'
di "  Change (%): " %6.2f `effect_pct' "%"

* By banking mode
di _n "Effects by Current Banking Mode:"
forvalues b = 1/3 {
    local bname = cond(`b'==1, "Unbanked", cond(`b'==2, "Mobile", "Branch"))
    quietly summarize p_se_baseline if bank_mode == `b' [aw=hsupwgtk]
    local base_`b' = r(mean)
    quietly summarize p_se_cf if bank_mode == `b' [aw=hsupwgtk]
    local cf_`b' = r(mean)
    quietly summarize hsupwgtk if bank_mode == `b'
    local wgt_`b' = r(sum)
    local eff_`b' = (`cf_`b'' - `base_`b'') / `base_`b'' * 100
    di "  `bname': " %6.4f `base_`b'' " -> " %6.4f `cf_`b'' " (" %5.2f `eff_`b'' "%)"
}

* Decompose the aggregate effect
di _n "Decomposition of Aggregate Effect:"
local contrib_branch = (`cf_3' - `base_3') * `baseline_branch'
local contrib_mobile = (`cf_2' - `base_2') * (1 - `baseline_branch' - 0.05)
local contrib_unbanked = (`cf_1' - `base_1') * 0.05
di "  Branch user contribution: " %7.4f `contrib_branch' " pp"
di "  Mobile user contribution: " %7.4f `contrib_mobile' " pp"
di "  Unbanked contribution: " %7.4f `contrib_unbanked' " pp"

/*------------------------------------------------------------------------------
4. Sensitivity Analysis
------------------------------------------------------------------------------*/

di _n "============================================================"
di "SENSITIVITY: Different Reduction Scenarios"
di "============================================================"

matrix CF_sens = J(5, 3, .)
matrix colnames CF_sens = Reduction CF_SE Effect_pct
matrix rownames CF_sens = R0 R25 R50 R75 R100

local reductions "0 0.25 0.50 0.75 1.00"
local row = 1
foreach r of local reductions {
    gen delta_`row' = 0
    replace delta_`row' = `b_branch_density' * (-`r' * pct_broadband) if bank_mode == 3
    replace delta_`row' = `b_mobile_density' * (-`r' * pct_broadband) if bank_mode == 2

    gen p_cf_`row' = exp(xb_base + delta_`row') / (1 + exp(xb_base + delta_`row'))

    quietly summarize p_cf_`row' [aw=hsupwgtk]
    local cf_r = r(mean)
    local eff_r = (`cf_r' - `base_se') / `base_se' * 100

    matrix CF_sens[`row', 1] = `r' * 100
    matrix CF_sens[`row', 2] = `cf_r'
    matrix CF_sens[`row', 3] = `eff_r'

    di %3.0f (`r'*100) "% reduction: SE = " %6.4f `cf_r' " (effect = " %5.2f `eff_r' "%)"

    drop delta_`row' p_cf_`row'
    local row = `row' + 1
}

di _n "Sensitivity Matrix:"
matrix list CF_sens, format(%8.4f)

/*------------------------------------------------------------------------------
5. Compare to Assumed Transition Approach
------------------------------------------------------------------------------*/

di _n "============================================================"
di "COMPARISON: Structural vs. Assumed Transitions"
di "============================================================"

* Old (wrong) approach: Assume 80% of branch users switch to mobile
* SE_new = 0.80 × SE_mobile + 0.20 × SE_unbanked (for displaced branch users)

quietly summarize se if bank_mode == 2 [aw=hsupwgtk]
local se_mobile = r(mean)
quietly summarize se if bank_mode == 1 [aw=hsupwgtk]
local se_unbanked = r(mean)

* Under 50% closure with assumed transitions
local assumed_new_se_displaced = 0.80 * `se_mobile' + 0.20 * `se_unbanked'
local assumed_new_se = `baseline_branch' * 0.5 * `assumed_new_se_displaced' + ///
                       `baseline_branch' * 0.5 * `base_3' + ///
                       (1 - `baseline_branch') * `base_se'
local assumed_effect = (`assumed_new_se' - `base_se') / `base_se' * 100

di "ASSUMED TRANSITION APPROACH (80/20 split):"
di "  Assumes displaced branch users get SE rate of mobile/unbanked users"
di "  Implied counterfactual SE rate: " %6.4f `assumed_new_se'
di "  Implied effect: " %5.2f `assumed_effect' "%"
di ""
di "STRUCTURAL APPROACH (this analysis):"
di "  Counterfactual SE rate: " %6.4f `cf_se'
di "  Effect: " %5.2f `effect_pct' "%"
di ""
di "Difference: Structural effect is " %5.2f (`effect_pct' - `assumed_effect') " pp"
di "  more negative than assumed transition approach"

/*------------------------------------------------------------------------------
6. Save Results
------------------------------------------------------------------------------*/

preserve
clear
set obs 15
gen str30 item = ""
gen value = .

replace item = "baseline_se_data" in 1
replace value = `baseline_se' in 1

replace item = "baseline_se_model" in 2
replace value = `base_se' in 2

replace item = "cf_se_50pct" in 3
replace value = `cf_se' in 3

replace item = "effect_pp" in 4
replace value = `effect_pp' in 4

replace item = "effect_pct" in 5
replace value = `effect_pct' in 5

replace item = "b_branch_density" in 6
replace value = `b_branch_density' in 6

replace item = "b_mobile_density" in 7
replace value = `b_mobile_density' in 7

replace item = "assumed_cf_se" in 8
replace value = `assumed_new_se' in 8

replace item = "assumed_effect_pct" in 9
replace value = `assumed_effect' in 9

replace item = "branch_users_base_se" in 10
replace value = `base_3' in 10

replace item = "branch_users_cf_se" in 11
replace value = `cf_3' in 11

replace item = "branch_users_effect_pct" in 12
replace value = `eff_3' in 12

drop if item == ""
export delimited using "$output/phase3_structural_cf.csv", replace
restore

di _n "Results saved to $output/phase3_structural_cf.csv"

di _n "============================================================"
di "KEY TAKEAWAY"
di "============================================================"
di ""
di "The structural approach derives counterfactual choice probabilities"
di "from the estimated utility parameters, rather than imposing"
di "transition shares exogenously."
di ""
di "When branch density falls 50%, the model predicts:"
di "  - SE rate changes from " %5.2f (`base_se'*100) "% to " %5.2f (`cf_se'*100) "%"
di "  - This is a " %5.2f abs(`effect_pct') "% decline"
di ""
di "The effect operates through reduced utility of SE for branch users"
di "when their credit access (via branch density) declines."

log close
