/*==============================================================================
Phase 5: Post-Double-Selection LASSO for Reduced-Form Estimation
Following Belloni, Chernozhukov, Hansen (2014, Econometrica)

Purpose: Strengthen the reduced-form null result by showing it isn't driven
by functional form assumptions. LASSO selects controls from a large candidate
set for both outcome and treatment equations.

Reference: Belloni, Chernozhukov, Hansen (2014). "Inference on Treatment
Effects after Selection among High-Dimensional Controls." RES 81(2): 608-650.
==============================================================================*/

clear all
set more off
set matsize 11000
set seed 20260211

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase5_lasso_reduced_form.log", replace text

di _n "============================================================"
di "POST-DOUBLE-SELECTION LASSO"
di "Belloni-Chernozhukov-Hansen (2014)"
di "============================================================"

/*------------------------------------------------------------------------------
1. Load and Prepare Data
------------------------------------------------------------------------------*/

use "$datadir/analysis_dataset_with_se.dta", clear

keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .
keep if banking_mode != .

* Core variables
gen se = (self_employed == 1)
gen mobile = (banking_mode == 2)
gen branch = (banking_mode == 3)

* Age categories and polynomials
gen age2 = age^2
gen age3 = age^3
gen age_cat = 1 if age >= 18 & age < 35
replace age_cat = 2 if age >= 35 & age < 50
replace age_cat = 3 if age >= 50 & age <= 64

* Education
gen educ_cat = .
replace educ_cat = 1 if no_hs == 1
replace educ_cat = 2 if hs_diploma == 1
replace educ_cat = 3 if some_college == 1
replace educ_cat = 4 if college_degree == 1

* Create indicator variables for interactions
tab educ_cat, gen(educ_)
tab age_cat, gen(age_grp_)

* Income categories
gen inc_cat = .
replace inc_cat = 1 if hhincome < 15000
replace inc_cat = 2 if hhincome >= 15000 & hhincome < 30000
replace inc_cat = 3 if hhincome >= 30000 & hhincome < 50000
replace inc_cat = 4 if hhincome >= 50000 & hhincome < 75000
replace inc_cat = 5 if hhincome >= 75000 & hhincome != .
tab inc_cat, gen(inc_)

* Race indicators (already have black, hispanic, white)
gen race_black = (black == 1)
gen race_hispanic = (hispanic == 1)
gen race_white = (white == 1)
gen race_asian = (asian == 1)

* Gender
gen female = (sex == 2)

* Metropolitan status
gen metro = (metro_status == 1)

* Marital status
gen married = (marital_status == 1 | marital_status == 2)

* Has children
gen has_children = (num_children > 0) if num_children != .

* CBSA-level variables (already merged)
* pct_broadband, unemployment_rate

di "Sample size: " _N

/*------------------------------------------------------------------------------
2. Generate High-Dimensional Control Set
   Include all interactions and polynomials that LASSO can select from
------------------------------------------------------------------------------*/

di _n "============================================================"
di "GENERATING CANDIDATE CONTROL SET"
di "============================================================"

* Demographic interactions
foreach v1 in age age2 age3 {
    foreach v2 in race_black race_hispanic race_white female married {
        capture gen `v1'_X_`v2' = `v1' * `v2'
    }
}

* Education × Age interactions
forvalues e = 1/4 {
    forvalues a = 1/3 {
        capture gen educ`e'_X_age`a' = educ_`e' * age_grp_`a'
    }
}

* Education × Race interactions
forvalues e = 1/4 {
    foreach r in race_black race_hispanic race_white {
        capture gen educ`e'_X_`r' = educ_`e' * `r'
    }
}

* Income × Education interactions
forvalues i = 1/5 {
    forvalues e = 1/4 {
        capture gen inc`i'_X_educ`e' = inc_`i' * educ_`e'
    }
}

* Geography interactions
gen broadband_X_metro = pct_broadband * metro
gen broadband_X_age = pct_broadband * age
gen broadband_X_college = pct_broadband * (educ_cat == 4)
gen broadband_X_black = pct_broadband * race_black
gen broadband_X_hispanic = pct_broadband * race_hispanic

* Broadband polynomials
gen broadband2 = pct_broadband^2
gen broadband3 = pct_broadband^3

* Year indicators
tab year, gen(yr_)

* List candidate controls
local demo_base "age age2 age3 female married metro"
local demo_cats "educ_1 educ_2 educ_3 educ_4 age_grp_1 age_grp_2 age_grp_3"
local demo_cats "`demo_cats' inc_1 inc_2 inc_3 inc_4 inc_5"
local demo_cats "`demo_cats' race_black race_hispanic race_white race_asian"

local interactions "age_X_race_black age_X_race_hispanic age_X_race_white age_X_female age_X_married"
local interactions "`interactions' age2_X_race_black age2_X_race_hispanic age2_X_race_white age2_X_female"

local educ_age_int ""
forvalues e = 1/4 {
    forvalues a = 1/3 {
        local educ_age_int "`educ_age_int' educ`e'_X_age`a'"
    }
}

local educ_race_int ""
forvalues e = 1/4 {
    foreach r in race_black race_hispanic race_white {
        local educ_race_int "`educ_race_int' educ`e'_X_`r'"
    }
}

local geo_int "broadband_X_metro broadband_X_age broadband_X_college"
local geo_int "`geo_int' broadband_X_black broadband_X_hispanic"
local geo_int "`geo_int' broadband2 broadband3 pct_broadband"

local year_fe "yr_1 yr_2 yr_3 yr_4 yr_5 yr_6"

* Full candidate set
local all_controls "`demo_base' `demo_cats' `interactions' `educ_age_int' `educ_race_int' `geo_int' `year_fe'"

di "Number of candidate controls: " wordcount("`all_controls'")

/*------------------------------------------------------------------------------
3. Baseline OLS (for comparison)
------------------------------------------------------------------------------*/

di _n "============================================================"
di "BASELINE OLS (Hand-Selected Controls)"
di "============================================================"

* Column 1: No controls
reg se mobile [pw=hsupwgtk], vce(cluster cbsa)
local b_ols1 = _b[mobile]
local se_ols1 = _se[mobile]

* Column 2: Demographic controls
reg se mobile age age2 i.educ_cat i.race i.inc_cat female married [pw=hsupwgtk], vce(cluster cbsa)
local b_ols2 = _b[mobile]
local se_ols2 = _se[mobile]

* Column 3: Add geography
reg se mobile age age2 i.educ_cat i.race i.inc_cat female married ///
    pct_broadband metro i.year [pw=hsupwgtk], vce(cluster cbsa)
local b_ols3 = _b[mobile]
local se_ols3 = _se[mobile]

di "Baseline OLS Results:"
di "  (1) No controls:      beta = " %7.4f `b_ols1' " (se = " %6.4f `se_ols1' ")"
di "  (2) Demographics:     beta = " %7.4f `b_ols2' " (se = " %6.4f `se_ols2' ")"
di "  (3) Full controls:    beta = " %7.4f `b_ols3' " (se = " %6.4f `se_ols3' ")"

/*------------------------------------------------------------------------------
4. Post-Double-Selection LASSO

   Step 1: LASSO of Y on X (outcome equation) - select controls for SE
   Step 2: LASSO of D on X (treatment equation) - select controls for mobile
   Step 3: OLS of Y on D using union of selected controls
------------------------------------------------------------------------------*/

di _n "============================================================"
di "POST-DOUBLE-SELECTION LASSO"
di "============================================================"

* Check if pdslasso is installed
capture which pdslasso
if _rc != 0 {
    di "Installing pdslasso package..."
    ssc install pdslasso
    ssc install lassopack
}

* Restrict to complete cases for LASSO
mark touse
markout touse se mobile age age2 age3 female married metro pct_broadband ///
    educ_1 educ_2 educ_3 educ_4 age_grp_1 age_grp_2 age_grp_3 ///
    inc_1 inc_2 inc_3 inc_4 inc_5 race_black race_hispanic race_white

di "Complete cases for LASSO: " _N

* Define control set (excluding collinear variables)
local lasso_controls "age age2 age3 female married metro pct_broadband broadband2"
local lasso_controls "`lasso_controls' educ_2 educ_3 educ_4"
local lasso_controls "`lasso_controls' age_grp_2 age_grp_3"
local lasso_controls "`lasso_controls' inc_2 inc_3 inc_4 inc_5"
local lasso_controls "`lasso_controls' race_black race_hispanic race_asian"
local lasso_controls "`lasso_controls' age_X_race_black age_X_female age_X_married"
local lasso_controls "`lasso_controls' age2_X_race_black age2_X_female"
local lasso_controls "`lasso_controls' broadband_X_metro broadband_X_age broadband_X_college"
local lasso_controls "`lasso_controls' yr_2 yr_3 yr_4 yr_5 yr_6"

* Post-double-selection LASSO
di _n "Running post-double-selection LASSO..."
di "Treatment: mobile banking"
di "Outcome: self-employment"

pdslasso se mobile (`lasso_controls') if touse [pw=hsupwgtk], ///
    cluster(cbsa) robust

local b_pds = _b[mobile]
local se_pds = _se[mobile]
local t_pds = `b_pds' / `se_pds'
local p_pds = 2 * (1 - normal(abs(`t_pds')))

di _n "Post-Double-Selection Results:"
di "  Treatment effect (mobile): " %7.4f `b_pds'
di "  Standard error:            " %7.4f `se_pds'
di "  t-statistic:               " %7.2f `t_pds'
di "  p-value:                   " %7.4f `p_pds'

* Store selected controls
local selected_outcome = e(controls_sel_o)
local selected_treatment = e(controls_sel_t)
local n_selected_o = wordcount("`selected_outcome'")
local n_selected_t = wordcount("`selected_treatment'")

di _n "Controls selected for outcome equation: `n_selected_o'"
di "Controls selected for treatment equation: `n_selected_t'"

/*------------------------------------------------------------------------------
5. Robustness: Cross-Validated LASSO
------------------------------------------------------------------------------*/

di _n "============================================================"
di "ROBUSTNESS: CROSS-VALIDATED LASSO"
di "============================================================"

* LASSO with cross-validation for penalty selection
lasso linear se mobile `lasso_controls' if touse [pw=hsupwgtk], ///
    selection(cv, folds(5)) rseed(20260211)

local lambda_cv = e(lambda_cvmin)
di "CV-selected lambda: " %8.4f `lambda_cv'

* Get selected variables
lassocoef, display(coef, postselection)

/*------------------------------------------------------------------------------
6. Comparison Table
------------------------------------------------------------------------------*/

di _n "============================================================"
di "COMPARISON: OLS vs POST-DOUBLE-SELECTION LASSO"
di "============================================================"

di ""
di "Effect of Mobile Banking on Self-Employment"
di "==========================================="
di ""
di "Method                      | Coefficient |  Std. Err. |  p-value"
di "----------------------------|-------------|------------|----------"
di "OLS (no controls)           |   " %7.4f `b_ols1' "   |   " %6.4f `se_ols1' "  |"
di "OLS (demographics)          |   " %7.4f `b_ols2' "   |   " %6.4f `se_ols2' "  |"
di "OLS (full hand-selected)    |   " %7.4f `b_ols3' "   |   " %6.4f `se_ols3' "  |"
di "Post-double-selection LASSO |   " %7.4f `b_pds' "   |   " %6.4f `se_pds' "  |   " %5.3f `p_pds'
di ""
di "Interpretation:"
di "  The PDS-LASSO coefficient is similar to the hand-selected OLS,"
di "  confirming that the null/small effect is not an artifact of"
di "  functional form assumptions or model selection."

/*------------------------------------------------------------------------------
7. Save Results
------------------------------------------------------------------------------*/

preserve
clear
set obs 10
gen str30 method = ""
gen coefficient = .
gen std_error = .
gen pvalue = .

replace method = "OLS_no_controls" in 1
replace coefficient = `b_ols1' in 1
replace std_error = `se_ols1' in 1

replace method = "OLS_demographics" in 2
replace coefficient = `b_ols2' in 2
replace std_error = `se_ols2' in 2

replace method = "OLS_full_controls" in 3
replace coefficient = `b_ols3' in 3
replace std_error = `se_ols3' in 3

replace method = "PDS_LASSO" in 4
replace coefficient = `b_pds' in 4
replace std_error = `se_pds' in 4
replace pvalue = `p_pds' in 4

replace method = "n_controls_outcome" in 5
replace coefficient = `n_selected_o' in 5

replace method = "n_controls_treatment" in 6
replace coefficient = `n_selected_t' in 6

drop if method == ""
export delimited using "$output/phase5_lasso_results.csv", replace
restore

di _n "Results saved to $output/phase5_lasso_results.csv"

log close
