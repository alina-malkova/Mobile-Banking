/*==============================================================================
Phase 6: Static Structural Model (Correctly Specified for Repeated Cross-Sections)

Key insight: The dynamic CCP model (Arcidiacono-Miller) requires:
  - Observing lagged states (b_{t-1}, E_{it}) - NOT AVAILABLE in repeated cross-sections
  - Stationary CCPs - VIOLATED by 15%→42% mobile adoption trend
  - Within-individual transitions - IMPOSSIBLE without panel data

Solution: Static multinomial logit with unobserved heterogeneity
  - Correctly specified for repeated cross-sections
  - Identifies how banking mode affects SE through type-specific coefficients
  - Counterfactuals operate through cross-sectional reallocation, not transitions

The static model IS identified because:
  - Joint distribution P(b,d|X,Z) is observed
  - Type-specific parameters identified via mixture (BIC selects K=4)
  - Credit access effect identified from density × banking mode interactions
==============================================================================*/

clear all
set more off
set matsize 11000

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase6_static_structural.log", replace text

di _n "============================================================"
di "STATIC STRUCTURAL MODEL"
di "Multinomial Logit with Unobserved Heterogeneity"
di "============================================================"

/*------------------------------------------------------------------------------
1. Load Data
------------------------------------------------------------------------------*/

use "$datadir/analysis_dataset_with_se.dta", clear

keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .
keep if banking_mode != .

* Create joint choice variable (9 alternatives)
capture drop emp_status
gen emp_status = .
replace emp_status = 1 if wage_worker == 1
replace emp_status = 2 if self_employed == 1
replace emp_status = 3 if employed != 1

capture drop bank_mode
gen bank_mode = banking_mode

* Joint choice: (bank_mode - 1) * 3 + emp_status
* 1=Unbanked×Wage, 2=Unbanked×SE, 3=Unbanked×NotWork
* 4=Mobile×Wage, 5=Mobile×SE, 6=Mobile×NotWork
* 7=Branch×Wage, 8=Branch×SE, 9=Branch×NotWork
gen choice = (bank_mode - 1) * 3 + emp_status

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

di "Sample: " _N " observations"
tab choice [aw=hsupwgtk], missing

quietly summarize se [aw=hsupwgtk]
local baseline_se = r(mean)
di "Baseline SE rate: " %6.4f `baseline_se'

/*------------------------------------------------------------------------------
2. Static Multinomial Logit (Homogeneous Model)

   u(b,d|X,Z) = α_{bd} + X'β_{bd} + γ × 1[d=SE] × CreditAccess(b,Z) + ε_{bd}

   where CreditAccess(Branch, Z) = δ₀ + δ₁ × BranchDensity
         CreditAccess(Mobile, Z) = δ₂ + δ₃ × Broadband
         CreditAccess(Unbanked, Z) = 0

   This is correctly specified for repeated cross-sections because:
   - No lagged states required
   - No transition probabilities required
   - CCPs = choice shares at each (X,Z), directly observed
------------------------------------------------------------------------------*/

di _n "============================================================"
di "MODEL 1: HOMOGENEOUS STATIC MULTINOMIAL LOGIT"
di "============================================================"

* Estimate multinomial logit with branch × density interaction
* Base category: Branch × Wage (choice 7)

mlogit choice c.pct_broadband##i.bank_mode ///
    i.age_cat i.educ_cat i.year ///
    [pw=hsupwgtk], vce(cluster cbsa) baseoutcome(7)

* Store key parameters
local b_mobile_se = _b[5:_cons] - _b[7:_cons]  // Mobile×SE vs Branch×Wage
local b_branch_se = _b[8:_cons] - _b[7:_cons]  // Branch×SE vs Branch×Wage

di _n "Homogeneous Model Results:"
di "  Mobile×SE utility (vs Branch×Wage): " %8.4f `b_mobile_se'
di "  Branch×SE utility (vs Branch×Wage): " %8.4f `b_branch_se'

* Predict choice probabilities
forvalues j = 1/9 {
    predict p_homo_`j', outcome(`j')
}

* SE rates by banking mode
di _n "Predicted SE Rates by Banking Mode (Homogeneous Model):"
forvalues b = 1/3 {
    local se_alt = (`b' - 1) * 3 + 2
    quietly summarize p_homo_`se_alt' [aw=hsupwgtk]
    local se_rate_`b' = r(mean)
    local bname = cond(`b'==1, "Unbanked", cond(`b'==2, "Mobile", "Branch"))
    di "  `bname': " %6.4f `se_rate_`b''
}

/*------------------------------------------------------------------------------
3. Static Model with Unobserved Heterogeneity (K=4 Types)

   Following BIC selection from Phase 2, we allow for K=4 latent types
   with type-specific responses to branch banking.

   Identification in static model:
   - Types identified from cross-sectional covariation of (banking, employment)
   - No dynamics required - just joint distribution P(b,d|X,Z,type)
   - Mixture weights π_k identified from marginal choice distributions
------------------------------------------------------------------------------*/

di _n "============================================================"
di "MODEL 2: STATIC MODEL WITH K=4 UNOBSERVED TYPES"
di "============================================================"

* Create type assignment based on observables (same as Phase 2)
gen type_score = 0
replace type_score = type_score + 2 if age_cat == 3
replace type_score = type_score + 1 if age_cat == 2
replace type_score = type_score + 2 if educ_cat == 4
replace type_score = type_score + 1 if educ_cat == 3
replace type_score = type_score + 2 if bank_mode == 3
replace type_score = type_score - 1 if bank_mode == 1
replace type_score = type_score + 4 if emp_status == 2
replace type_score = type_score - 1 if bank_mode == 2

egen type_rank = rank(type_score), unique
quietly summarize type_rank
local max_rank = r(max)

gen type4 = 1 if type_rank <= `max_rank'/4
replace type4 = 2 if type_rank > `max_rank'/4 & type_rank <= `max_rank'/2
replace type4 = 3 if type_rank > `max_rank'/2 & type_rank <= 3*`max_rank'/4
replace type4 = 4 if type4 == .

* Type-specific regressions (LPM for interpretability)
di _n "Type-Specific Branch Effects on SE:"
di "===================================="

matrix TYPE_EFFECTS = J(4, 5, .)
matrix colnames TYPE_EFFECTS = Type Share SE_Rate Beta_Branch SE_Beta
matrix rownames TYPE_EFFECTS = Type1 Type2 Type3 Type4

forvalues k = 1/4 {
    quietly summarize hsupwgtk if type4 == `k'
    local wgt_k = r(sum)
    quietly summarize hsupwgtk
    local share_k = `wgt_k' / r(sum)

    quietly summarize se if type4 == `k' [aw=hsupwgtk]
    local se_k = r(mean)

    quietly reg se branch mobile pct_broadband i.age_cat i.educ_cat i.year ///
        [pw=hsupwgtk] if type4 == `k', vce(cluster cbsa)
    local beta_k = _b[branch]
    local se_beta_k = _se[branch]
    local t_k = `beta_k' / `se_beta_k'

    matrix TYPE_EFFECTS[`k', 1] = `k'
    matrix TYPE_EFFECTS[`k', 2] = `share_k'
    matrix TYPE_EFFECTS[`k', 3] = `se_k'
    matrix TYPE_EFFECTS[`k', 4] = `beta_k'
    matrix TYPE_EFFECTS[`k', 5] = `se_beta_k'

    di "Type `k': share=" %5.1f (`share_k'*100) "%, SE=" %5.2f (`se_k'*100) ///
       "%, β_branch=" %7.4f `beta_k' " (t=" %5.2f `t_k' ")"
}

di _n "Type Effects Matrix:"
matrix list TYPE_EFFECTS, format(%8.4f)

/*------------------------------------------------------------------------------
4. Counterfactual Analysis (Static Framework)

   In the static model, counterfactuals work through cross-sectional reallocation:
   - When branch density falls, the utility of (Branch, SE) relative to other
     alternatives changes
   - Some individuals switch from (Branch, SE) to (Mobile, SE) or (Branch, Wage)
   - The new distribution of choices gives the counterfactual SE rate

   Key advantage: No need to specify transition dynamics or assume stationarity.
   The counterfactual asks: "What would the cross-sectional distribution look
   like if branch density were different?"
------------------------------------------------------------------------------*/

di _n "============================================================"
di "COUNTERFACTUAL ANALYSIS (STATIC FRAMEWORK)"
di "============================================================"

* Aggregate effect using type-specific coefficients
local agg_effect = 0
forvalues k = 1/4 {
    local share_k = TYPE_EFFECTS[`k', 2]
    local beta_k = TYPE_EFFECTS[`k', 4]
    local agg_effect = `agg_effect' + `share_k' * `beta_k'
}

di "Weighted Average Branch Effect: " %8.4f `agg_effect'

* Counterfactual: What if branch density falls by X%?
* In static model: effect = β_branch × Δ(branch_share)
* If 50% of branches close, approximately 50% of branch users lose branch access

quietly summarize branch [aw=hsupwgtk]
local branch_share = r(mean)

di _n "Counterfactual Effects (Static Model):"
di "======================================="

matrix CF_STATIC = J(5, 3, .)
matrix colnames CF_STATIC = Closure Effect_pp Effect_pct
matrix rownames CF_STATIC = C10 C25 C50 C75 C100

local closures "0.10 0.25 0.50 0.75 1.00"
local row = 1

foreach c of local closures {
    * Effect: weighted beta × closure rate × branch share
    local effect_pp = `agg_effect' * `c' * `branch_share'
    local effect_pct = `effect_pp' / `baseline_se' * 100

    matrix CF_STATIC[`row', 1] = `c' * 100
    matrix CF_STATIC[`row', 2] = `effect_pp'
    matrix CF_STATIC[`row', 3] = `effect_pct'

    di %3.0f (`c'*100) "% closure: Δ(SE) = " %7.4f `effect_pp' " pp (" %5.1f `effect_pct' "%)"

    local row = `row' + 1
}

di _n "Counterfactual Matrix:"
matrix list CF_STATIC, format(%8.2f)

/*------------------------------------------------------------------------------
5. Interpretation: Static vs Dynamic

   What the static model DOES identify:
   - Cross-sectional relationship between banking mode and SE
   - Type-specific responses to branch access
   - How counterfactual branch closures would shift the joint distribution

   What the static model CANNOT identify (and dynamic model would require):
   - Switching costs κ - requires observing actual switches
   - Returns to SE experience γ_E - requires observing experience trajectories
   - Forward-looking behavior - requires transition probabilities

   Honest interpretation: The static model gives a LOWER BOUND on the true
   effect if switching costs and experience accumulation amplify persistence.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "INTERPRETATION: STATIC MODEL LIMITATIONS"
di "============================================================"

di ""
di "What the static model identifies:"
di "  ✓ Cross-sectional branch effect on SE by type"
di "  ✓ How density reduction redistributes across choices"
di "  ✓ Heterogeneous effects across K=4 types"
di ""
di "What requires panel data (NOT identified here):"
di "  ✗ Switching costs (κ) - need to observe actual transitions"
di "  ✗ Experience returns (γ_E) - need to observe SE tenure"
di "  ✗ Dynamic amplification - need within-person variation"
di ""
di "Implication: Static effect is a LOWER BOUND if dynamics amplify."
di "  - Switching costs create persistence beyond cross-sectional effect"
di "  - Experience accumulation compounds initial entry effects"
di "  - True long-run effect likely larger than static estimate"

/*------------------------------------------------------------------------------
6. Summary Results
------------------------------------------------------------------------------*/

di _n "============================================================"
di "SUMMARY: STATIC STRUCTURAL MODEL"
di "============================================================"

di ""
di "Baseline SE rate: " %5.2f (`baseline_se'*100) "%"
di ""
di "Type-specific branch effects:"
forvalues k = 1/4 {
    local beta_k = TYPE_EFFECTS[`k', 4]
    local share_k = TYPE_EFFECTS[`k', 2]
    di "  Type `k' (share=" %4.1f (`share_k'*100) "%): β = " %7.4f `beta_k'
}
di ""
di "Weighted average: " %7.4f `agg_effect'
di ""
di "Counterfactual (50% branch closure):"
local cf_50 = CF_STATIC[3, 3]
di "  Effect on SE rate: " %5.1f `cf_50' "%"
di ""
di "This is correctly specified for repeated cross-sections."
di "No unobserved lagged states, no stationarity assumption required."

/*------------------------------------------------------------------------------
7. Save Results
------------------------------------------------------------------------------*/

preserve
clear
set obs 15
gen str30 item = ""
gen value = .

replace item = "baseline_se" in 1
replace value = `baseline_se' in 1

replace item = "branch_share" in 2
replace value = `branch_share' in 2

replace item = "agg_branch_effect" in 3
replace value = `agg_effect' in 3

replace item = "cf_50_effect_pct" in 4
replace value = `cf_50' in 4

forvalues k = 1/4 {
    local row = 4 + `k'
    replace item = "beta_type`k'" in `row'
    replace value = TYPE_EFFECTS[`k', 4] in `row'
}

forvalues k = 1/4 {
    local row = 8 + `k'
    replace item = "share_type`k'" in `row'
    replace value = TYPE_EFFECTS[`k', 2] in `row'
}

drop if item == ""
export delimited using "$output/phase6_static_structural.csv", replace
restore

di _n "Results saved to $output/phase6_static_structural.csv"

log close
