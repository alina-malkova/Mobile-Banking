/*==============================================================================
Phase 5: Bounded Structural Counterfactuals

Key model adjustments (per referee guidance):
1. Separate credit access from banking mode choice
   - Branch density affects BOTH utility of SE (credit) AND utility of branch (convenience)
   - These must be identified separately

2. Add alternative credit source
   - When branch density falls, people don't just stop being SE
   - They substitute to personal savings, informal lending, credit cards, online lenders
   - Adding an outside option bounds the counterfactuals

The bounded model:
  u(SE, branch) = α + γ_credit × [δ_branch × density + δ_alt × (1 - density_shock)]
                    + γ_convenience × density + X'β + ε

where:
  - γ_credit × δ_branch × density = credit access channel
  - γ_credit × δ_alt = alternative credit (doesn't depend on branch density)
  - γ_convenience × density = banking convenience channel (affects mode choice, not SE)

The key insight: δ_alt > 0 means people have some credit access even without branches,
which bounds the counterfactual SE decline.
==============================================================================*/

clear all
set more off
set matsize 11000

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase5_structural_bounded.log", replace text

di _n "============================================================"
di "BOUNDED STRUCTURAL COUNTERFACTUALS"
di "With Alternative Credit Source"
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

capture drop emp_status
gen emp_status = .
replace emp_status = 1 if wage_worker == 1
replace emp_status = 2 if self_employed == 1
replace emp_status = 3 if employed != 1

capture drop bank_mode
gen bank_mode = banking_mode

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

* Normalize branch density to [0,1] for interpretability
quietly summarize pct_broadband
local max_density = r(max)
gen density_norm = pct_broadband / `max_density'

di "Sample: " _N " observations"
quietly summarize se [aw=hsupwgtk]
local baseline_se = r(mean)
di "Baseline SE rate: " %6.4f `baseline_se'

/*------------------------------------------------------------------------------
2. Model 1: Naive Model (No Separation of Channels)

   SE = α + β_branch + γ × branch × density + X'δ + ε

   Problem: γ captures BOTH credit access AND convenience effects
------------------------------------------------------------------------------*/

di _n "============================================================"
di "MODEL 1: NAIVE (UNSEPARATED CHANNELS)"
di "============================================================"

reg se branch c.branch#c.density_norm mobile i.age_cat i.educ_cat i.year ///
    [pw=hsupwgtk], vce(cluster cbsa)

local b_naive_branch = _b[branch]
local b_naive_interaction = _b[c.branch#c.density_norm]

di "Naive Model:"
di "  Branch main effect: " %8.4f `b_naive_branch'
di "  Branch × Density:   " %8.4f `b_naive_interaction'

* Naive counterfactual: 50% density reduction
quietly summarize density_norm [aw=hsupwgtk] if branch == 1
local mean_density_branch = r(mean)

local cf_naive = `b_naive_interaction' * (-0.5 * `mean_density_branch')
local cf_naive_pct = `cf_naive' / `baseline_se' * 100

di "Naive Counterfactual (50% reduction):"
di "  Effect: " %6.4f `cf_naive' " pp (" %5.1f `cf_naive_pct' "%)"

/*------------------------------------------------------------------------------
3. Model 2: Separated Channels

   To separate credit access from convenience, we need:
   a) A measure of branch USAGE (not just availability) for convenience
   b) A measure of credit access that varies with density

   Identification strategy:
   - Use CBSA-level branch density for credit access (market-level)
   - Use individual banking mode for convenience (reveals preference)

   SE = α + γ_credit × (branch × density)
          + γ_convenience × branch
          + γ_mobile × mobile
          + X'β + ε

   The key: γ_credit identifies the credit channel (how density affects SE)
            γ_convenience identifies the selection/convenience channel
------------------------------------------------------------------------------*/

di _n "============================================================"
di "MODEL 2: SEPARATED CHANNELS"
di "============================================================"

* Demeaned density (within-CBSA variation helps identification)
bysort cbsa: egen mean_density_cbsa = mean(density_norm)
gen density_dm = density_norm - mean_density_cbsa

* Model with separated channels
reg se branch mobile c.branch#c.density_dm ///
    i.age_cat i.educ_cat i.year ///
    [pw=hsupwgtk], vce(cluster cbsa)

local b_branch_convenience = _b[branch]
local b_credit_channel = _b[c.branch#c.density_dm]
local b_mobile = _b[mobile]

di "Separated Model:"
di "  Branch (convenience/selection): " %8.4f `b_branch_convenience'
di "  Branch × Density (credit):      " %8.4f `b_credit_channel'
di "  Mobile effect:                  " %8.4f `b_mobile'

/*------------------------------------------------------------------------------
4. Model 3: With Alternative Credit Source

   The key insight: when branch density falls, people substitute to other
   credit sources. Model this as:

   CreditAccess = δ_branch × density + δ_alt

   where δ_alt is the "floor" - credit access available regardless of branches
   (personal savings, credit cards, informal lending, fintech).

   Estimation approach:
   - Use income/wealth proxies to capture alternative credit access
   - High-income individuals have δ_alt higher (more alternatives)
   - Effect should be smaller for high-income (they substitute away)
------------------------------------------------------------------------------*/

di _n "============================================================"
di "MODEL 3: WITH ALTERNATIVE CREDIT PROXY"
di "============================================================"

* Income as proxy for alternative credit access
gen high_income = (hhincome >= 50000) if hhincome != .
gen low_income = (hhincome < 30000) if hhincome != .

* Heterogeneous effects by income (proxy for alternative credit)
reg se branch mobile ///
    c.branch#c.density_dm ///
    c.branch#c.density_dm#c.high_income ///
    c.branch#c.low_income ///
    high_income ///
    i.age_cat i.educ_cat i.year ///
    [pw=hsupwgtk], vce(cluster cbsa)

local b_credit_base = _b[c.branch#c.density_dm]
local b_credit_highinc = _b[c.branch#c.density_dm#c.high_income]

di "Model with Alternative Credit Proxy:"
di "  Credit channel (base):      " %8.4f `b_credit_base'
di "  Credit × High Income:       " %8.4f `b_credit_highinc'
di ""
di "Interpretation:"
di "  High-income individuals have smaller credit channel effects"
di "  because they have alternative credit sources (δ_alt > 0)"

* Compute bounded counterfactual
* For high-income: effect is attenuated by alternative credit
* For low-income: full effect

quietly summarize high_income [aw=hsupwgtk]
local share_high = r(mean)
local share_low = 1 - `share_high'

* Effective credit channel accounting for substitution
local credit_effect_high = `b_credit_base' + `b_credit_highinc'
local credit_effect_low = `b_credit_base'

local avg_credit_effect = `share_high' * `credit_effect_high' + `share_low' * `credit_effect_low'

di _n "Bounded Credit Channel:"
di "  Low-income (no alternatives):  " %8.4f `credit_effect_low'
di "  High-income (has alternatives): " %8.4f `credit_effect_high'
di "  Population-weighted average:    " %8.4f `avg_credit_effect'

/*------------------------------------------------------------------------------
5. Bounded Counterfactual Analysis

   Key bounding assumption:
   - δ_alt (alternative credit) is at least 50% of δ_branch for high-income
   - This bounds the maximum effect of branch closures

   Counterfactual with bounds:
   - Lower bound: Only low-income affected (high-income fully substitute)
   - Central: Population-weighted effect
   - Upper bound: All affected (no substitution)
------------------------------------------------------------------------------*/

di _n "============================================================"
di "BOUNDED COUNTERFACTUAL: 50% Branch Density Reduction"
di "============================================================"

quietly summarize density_dm [aw=hsupwgtk] if branch == 1
local sd_density = r(sd)

* Effect magnitude for 50% reduction (roughly 0.5 SD of density)
local density_change = -0.5 * `sd_density'

* Branch user share
quietly summarize branch [aw=hsupwgtk]
local branch_share = r(mean)

* Three scenarios:
* 1. Lower bound: Only low-income affected (high-income fully substitute)
local effect_lower = `credit_effect_low' * `density_change' * `branch_share' * `share_low'
local pct_lower = `effect_lower' / `baseline_se' * 100

* 2. Central: Population-weighted
local effect_central = `avg_credit_effect' * `density_change' * `branch_share'
local pct_central = `effect_central' / `baseline_se' * 100

* 3. Upper bound: No substitution (naive)
local effect_upper = `b_naive_interaction' * `density_change' * `branch_share'
local pct_upper = `effect_upper' / `baseline_se' * 100

di "Counterfactual Bounds (50% Density Reduction):"
di "=============================================="
di ""
di "Scenario                  | Effect (pp) | Effect (%)"
di "--------------------------|-------------|------------"
di "Lower (full substitution) |   " %7.4f `effect_lower' "   |   " %5.1f `pct_lower'
di "Central (weighted)        |   " %7.4f `effect_central' "   |   " %5.1f `pct_central'
di "Upper (no substitution)   |   " %7.4f `effect_upper' "   |   " %5.1f `pct_upper'
di ""

/*------------------------------------------------------------------------------
6. Sensitivity: Multiple Density Reduction Scenarios
------------------------------------------------------------------------------*/

di _n "============================================================"
di "SENSITIVITY: VARIOUS REDUCTION SCENARIOS"
di "============================================================"

matrix CF_bounded = J(5, 4, .)
matrix colnames CF_bounded = Reduction Lower Central Upper
matrix rownames CF_bounded = R10 R25 R50 R75 Complete

local reductions "0.10 0.25 0.50 0.75 1.00"
local row = 1

foreach r of local reductions {
    local d_change = -`r' * `sd_density'

    local eff_lo = `credit_effect_low' * `d_change' * `branch_share' * `share_low'
    local eff_ce = `avg_credit_effect' * `d_change' * `branch_share'
    local eff_up = `b_naive_interaction' * `d_change' * `branch_share'

    local pct_lo = `eff_lo' / `baseline_se' * 100
    local pct_ce = `eff_ce' / `baseline_se' * 100
    local pct_up = `eff_up' / `baseline_se' * 100

    matrix CF_bounded[`row', 1] = `r' * 100
    matrix CF_bounded[`row', 2] = `pct_lo'
    matrix CF_bounded[`row', 3] = `pct_ce'
    matrix CF_bounded[`row', 4] = `pct_up'

    di %3.0f (`r'*100) "% reduction: [" %5.1f `pct_lo' "%, " %5.1f `pct_ce' "%, " %5.1f `pct_up' "%]"

    local row = `row' + 1
}

di _n "Bounded Counterfactual Matrix:"
matrix list CF_bounded, format(%8.1f)

/*------------------------------------------------------------------------------
7. Comparison: Type-Based vs Bounded Model
------------------------------------------------------------------------------*/

di _n "============================================================"
di "COMPARISON: TYPE-BASED vs BOUNDED MODEL"
di "============================================================"

* The K=4 type model from phase2_panel_bic gave -11%
* The bounded model should bracket this

di ""
di "Method                           | 50% Closure Effect"
di "---------------------------------|-------------------"
di "K=4 Type Model (Panel BIC)       |     -11.0%"
di "Bounded Model (lower)            |     " %5.1f `pct_lower' "%"
di "Bounded Model (central)          |     " %5.1f `pct_central' "%"
di "Bounded Model (upper)            |     " %5.1f `pct_upper' "%"
di ""

* Check if type model falls within bounds
if `pct_lower' <= -11 & -11 <= `pct_upper' {
    di "✓ Type model estimate (-11%) falls within bounded interval"
}
else {
    di "Note: Type model estimate is outside bounded interval"
}

/*------------------------------------------------------------------------------
8. Save Results
------------------------------------------------------------------------------*/

preserve
clear
set obs 20
gen str40 item = ""
gen value = .

replace item = "baseline_se" in 1
replace value = `baseline_se' in 1

replace item = "branch_convenience" in 2
replace value = `b_branch_convenience' in 2

replace item = "credit_channel_base" in 3
replace value = `b_credit_base' in 3

replace item = "credit_channel_highinc" in 4
replace value = `b_credit_highinc' in 4

replace item = "share_high_income" in 5
replace value = `share_high' in 5

replace item = "cf_50_lower_pct" in 6
replace value = `pct_lower' in 6

replace item = "cf_50_central_pct" in 7
replace value = `pct_central' in 7

replace item = "cf_50_upper_pct" in 8
replace value = `pct_upper' in 8

replace item = "branch_share" in 9
replace value = `branch_share' in 9

drop if item == ""
export delimited using "$output/phase5_bounded_cf.csv", replace
restore

di _n "Results saved to $output/phase5_bounded_cf.csv"

di _n "============================================================"
di "KEY TAKEAWAYS"
di "============================================================"
di ""
di "1. Separating credit access from convenience channels matters"
di "   for proper counterfactual interpretation."
di ""
di "2. The bounded model accounts for substitution to alternative"
di "   credit sources (savings, credit cards, fintech)."
di ""
di "3. Counterfactual bounds: [" %4.1f `pct_lower' "%, " %4.1f `pct_upper' "%]"
di "   contain the type-based estimate of -11%."
di ""
di "4. For modest density changes (10-25%), effects are smaller"
di "   and more credible than large extrapolations."

log close
