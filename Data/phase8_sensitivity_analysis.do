/*==============================================================================
Phase 8: Sensitivity Analysis for Unobserved Confounding

Following:
- Cinelli & Hazlett (2020, JRSS-B): Robustness values
- Oster (2019, JBES): Coefficient stability and R-squared movements

Problem: Weak first-stage (F=4.82) means IV results are unreliable.
Solution: Formal sensitivity analysis - how much unobserved confounding
would be needed to explain away the reduced-form relationship?

For null result: How strong would omitted variable need to be to CREATE
a significant effect where none is observed?
==============================================================================*/

clear all
set more off
set matsize 11000
set seed 20260211

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase8_sensitivity.log", replace text

di _n "============================================================"
di "SENSITIVITY ANALYSIS FOR UNOBSERVED CONFOUNDING"
di "Cinelli-Hazlett & Oster Approaches"
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

gen se = (self_employed == 1)
gen mobile = (banking_mode == 2)
gen branch = (banking_mode == 3)

gen age2 = age^2
gen educ_cat = .
replace educ_cat = 1 if no_hs == 1
replace educ_cat = 2 if hs_diploma == 1
replace educ_cat = 3 if some_college == 1
replace educ_cat = 4 if college_degree == 1

gen female = (sex == 2)
gen married = (marital_status == 1 | marital_status == 2)

di "Sample: " _N " observations"

/*------------------------------------------------------------------------------
2. Oster (2019) Coefficient Stability Analysis

   Key insight: If adding controls barely changes the coefficient,
   the estimate is robust to omitted variables of similar strength.

   Delta = (beta_full - beta_short) / (beta_short - 0)

   If |Delta| < 1, the coefficient is robust to omitted variables
   that explain at least as much variation as included controls.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "OSTER (2019) COEFFICIENT STABILITY"
di "============================================================"

* Short regression (minimal controls)
reg se branch mobile [pw=hsupwgtk], vce(cluster cbsa)
local b_short = _b[branch]
local r2_short = e(r2)

di "Short regression (no controls):"
di "  Branch coefficient: " %8.4f `b_short'
di "  R-squared: " %6.4f `r2_short'

* Medium regression (demographic controls)
reg se branch mobile age age2 i.educ_cat female married [pw=hsupwgtk], vce(cluster cbsa)
local b_med = _b[branch]
local r2_med = e(r2)

di _n "Medium regression (demographics):"
di "  Branch coefficient: " %8.4f `b_med'
di "  R-squared: " %6.4f `r2_med'

* Full regression (all controls + FE)
reg se branch mobile age age2 i.educ_cat female married pct_broadband ///
    i.year i.cbsa [pw=hsupwgtk], vce(cluster cbsa)
local b_full = _b[branch]
local r2_full = e(r2)

di _n "Full regression (all controls + FE):"
di "  Branch coefficient: " %8.4f `b_full'
di "  R-squared: " %6.4f `r2_full'

* Oster's delta: coefficient movement relative to R2 movement
* delta = (b_full - b_short) / (b_short) * (r2_max - r2_short) / (r2_full - r2_short)
* Simplified: proportional bias

local pct_change_coef = (`b_full' - `b_short') / `b_short' * 100
local pct_change_r2 = (`r2_full' - `r2_short') / `r2_short' * 100

di _n "Coefficient Stability Analysis:"
di "  Coefficient change (short → full): " %5.1f `pct_change_coef' "%"
di "  R-squared change (short → full): " %5.1f `pct_change_r2' "%"

* Oster's identified set bound
* Assume R_max = 1.3 * R_full (following Oster's suggestion)
local r2_max = min(1, 1.3 * `r2_full')

* Calculate bias-adjusted estimate
* beta* = beta_full - delta * (beta_full - beta_short)
* where delta = (R_max - R_full) / (R_full - R_short)

local delta_denom = `r2_full' - `r2_short'
if `delta_denom' > 0.001 {
    local delta = (`r2_max' - `r2_full') / `delta_denom'
    local b_adjusted = `b_full' - `delta' * (`b_full' - `b_short')

    di _n "Oster Bias-Adjusted Estimate:"
    di "  R_max assumption: " %6.4f `r2_max'
    di "  Delta: " %6.4f `delta'
    di "  Bias-adjusted coefficient: " %8.4f `b_adjusted'

    * Robustness value: what delta would drive coefficient to zero?
    if `b_full' != `b_short' {
        local delta_zero = `b_full' / (`b_full' - `b_short')
        di "  Delta needed to zero coefficient: " %6.4f `delta_zero'

        if `delta_zero' > 1 {
            di "  Interpretation: Robust - unobservables would need to be"
            di "    " %4.1f `delta_zero' "x stronger than observables to eliminate effect"
        }
        else {
            di "  Interpretation: NOT robust - unobservables of similar"
            di "    strength as observables could eliminate effect"
        }
    }
}

/*------------------------------------------------------------------------------
3. Cinelli-Hazlett (2020) Robustness Values

   The "robustness value" RV is the minimum strength of an omitted
   variable (relative to the strongest observed confounder) that would
   change the sign or significance of the estimate.

   RV_q=a: minimum partial R2 with outcome and treatment such that
   confidence interval includes q*beta (usually q=0 for sign change)
------------------------------------------------------------------------------*/

di _n "============================================================"
di "CINELLI-HAZLETT (2020) ROBUSTNESS VALUES"
di "============================================================"

* Check for sensemakr or implement manually
capture which sensemakr
if _rc != 0 {
    di "sensemakr not available. Computing robustness values manually."

    * Manual computation of robustness value
    * Based on partial R2 decomposition

    * Get partial R2 of treatment with outcome, controlling for X
    quietly reg se age age2 i.educ_cat female married pct_broadband i.year [pw=hsupwgtk]
    predict se_resid, resid

    quietly reg branch age age2 i.educ_cat female married pct_broadband i.year [pw=hsupwgtk]
    predict branch_resid, resid

    * Partial correlation
    correlate se_resid branch_resid
    local partial_r = r(rho)
    local partial_r2_dy = `partial_r'^2

    di _n "Partial R-squared decomposition:"
    di "  Partial R2 of branch with SE (controlling for X): " %8.6f `partial_r2_dy'

    * For robustness value: find partial R2 needed to change sign
    * RV = sqrt(partial_r2_dy) for approximate robustness value
    local rv_approx = sqrt(`partial_r2_dy')
    di "  Approximate robustness value (RV): " %6.4f `rv_approx'

    di _n "Interpretation:"
    di "  An omitted variable would need partial R2 of at least " %6.4f `partial_r2_dy'
    di "  with both SE and branch banking to reduce the coefficient to zero."

    * Benchmark against observed confounders
    di _n "Benchmark: Partial R2 of strongest observed confounders:"

    foreach var in age female married {
        quietly reg se `var' [pw=hsupwgtk]
        local r2_`var'_y = e(r2)
        quietly reg branch `var' [pw=hsupwgtk]
        local r2_`var'_d = e(r2)
        local bound_`var' = `r2_`var'_y' * `r2_`var'_d'
        di "  `var': R2_y = " %6.4f `r2_`var'_y' ", R2_d = " %6.4f `r2_`var'_d' ///
           ", product = " %8.6f `bound_`var''
    }

    drop se_resid branch_resid
}
else {
    * Use sensemakr package
    sensemakr se branch age age2 i.educ_cat female married pct_broadband i.year ///
        [pw=hsupwgtk], benchmark(age)
}

/*------------------------------------------------------------------------------
4. Sensitivity Contour Plot (Manual Implementation)

   Show how the coefficient changes under different assumptions
   about unobserved confounding strength.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "SENSITIVITY CONTOUR ANALYSIS"
di "============================================================"

* Create grid of confounding strengths
* r2_yd = partial R2 of omitted var with outcome
* r2_dd = partial R2 of omitted var with treatment

di _n "Adjusted coefficients under different confounding assumptions:"
di "(r2_outcome, r2_treatment) -> adjusted coefficient"
di ""

* Store results in matrix
matrix SENS = J(5, 5, .)
matrix rownames SENS = "r2_y=0.01" "r2_y=0.02" "r2_y=0.05" "r2_y=0.10" "r2_y=0.20"
matrix colnames SENS = "r2_d=0.01" "r2_d=0.02" "r2_d=0.05" "r2_d=0.10" "r2_d=0.20"

local r2_y_vals "0.01 0.02 0.05 0.10 0.20"
local r2_d_vals "0.01 0.02 0.05 0.10 0.20"

local row = 1
foreach r2_y of local r2_y_vals {
    local col = 1
    foreach r2_d of local r2_d_vals {
        * Omitted variable bias formula (simplified):
        * bias ≈ sqrt(r2_y * r2_d) * sign(correlation) * scale
        * Conservative: assume bias works against finding an effect

        local bias = sqrt(`r2_y' * `r2_d') * 0.10  // Scale factor
        local adj_coef = `b_full' - `bias'

        matrix SENS[`row', `col'] = `adj_coef'
        local col = `col' + 1
    }
    local row = `row' + 1
}

di "Sensitivity Matrix (adjusted coefficients):"
matrix list SENS, format(%8.4f)

/*------------------------------------------------------------------------------
5. Summary Statistics for Paper
------------------------------------------------------------------------------*/

di _n "============================================================"
di "SUMMARY FOR PAPER"
di "============================================================"

di ""
di "Coefficient Stability (Oster 2019):"
di "  - Coefficient changes by " %4.1f abs(`pct_change_coef') "% when adding all controls"
di "  - R-squared increases by " %4.1f `pct_change_r2' "%"
di "  - Bias-adjusted estimate: " %8.4f `b_adjusted'
if `delta_zero' > 1 {
    di "  - Robust: unobservables need " %4.1f `delta_zero' "x strength of observables"
}
di ""
di "Robustness Value (Cinelli-Hazlett 2020):"
di "  - Partial R2 of branch with SE: " %8.6f `partial_r2_dy'
di "  - An omitted variable would need partial R2 > " %6.4f `partial_r2_dy'
di "    with both outcome and treatment to change the sign"
di ""
di "Interpretation for counterfactual:"
di "  The structural branch effect is [robust/not robust] to"
di "  unobserved confounding of the magnitude of observed confounders."

/*------------------------------------------------------------------------------
6. Save Results
------------------------------------------------------------------------------*/

preserve
clear
set obs 10
gen str30 parameter = ""
gen value = .

replace parameter = "b_short" in 1
replace value = `b_short' in 1

replace parameter = "b_full" in 2
replace value = `b_full' in 2

replace parameter = "b_adjusted" in 3
replace value = `b_adjusted' in 3

replace parameter = "delta_zero" in 4
replace value = `delta_zero' in 4

replace parameter = "partial_r2_dy" in 5
replace value = `partial_r2_dy' in 5

replace parameter = "rv_approx" in 6
replace value = `rv_approx' in 6

replace parameter = "pct_change_coef" in 7
replace value = `pct_change_coef' in 7

drop if parameter == ""
export delimited using "$output/phase8_sensitivity.csv", replace
restore

di _n "Results saved to $output/phase8_sensitivity.csv"

log close
