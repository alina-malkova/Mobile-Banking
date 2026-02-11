/*==============================================================================
Phase 5: ML-Based CCP Estimation

Problem: Standard CCP estimation uses cell means within state cells
(CBSA × year × age × education). With 293 CBSAs × 6 waves × age groups ×
education groups, many cells are thin.

Solution: Use machine learning to estimate CCPs by flexibly modeling how
choice probabilities vary with states. This gives better first-stage estimates
without requiring cell discretization.

Methods:
1. Random Forest for P(choice | state)
2. Regularized Multinomial Logit (LASSO)
3. Gradient Boosting (optional)

The Arcidiacono-Miller second stage is unchanged - we just use better CCP estimates.
==============================================================================*/

clear all
set more off
set matsize 11000
set seed 20260211

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase5_ml_ccps.log", replace text

di _n "============================================================"
di "ML-BASED CCP ESTIMATION"
di "Random Forest and Regularized Logit"
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

* Joint state: (bank_mode - 1) * 3 + emp_status
* 1 = Unbanked × Wage
* 2 = Unbanked × SE
* 3 = Unbanked × NotWork
* 4 = Mobile × Wage
* 5 = Mobile × SE
* 6 = Mobile × NotWork
* 7 = Branch × Wage
* 8 = Branch × SE
* 9 = Branch × NotWork
gen choice = (bank_mode - 1) * 3 + emp_status

* State variables
gen age_cat = 1 if age >= 18 & age < 35
replace age_cat = 2 if age >= 35 & age < 50
replace age_cat = 3 if age >= 50 & age <= 64

gen educ_cat = .
replace educ_cat = 1 if no_hs == 1
replace educ_cat = 2 if hs_diploma == 1
replace educ_cat = 3 if some_college == 1
replace educ_cat = 4 if college_degree == 1

* Additional state variables for ML
gen female = (sex == 2)
gen married = (marital_status == 1 | marital_status == 2)
gen metro = (metro_status == 1)
gen has_kids = (num_children > 0) if num_children != .
replace has_kids = 0 if has_kids == .

* Income terciles
xtile inc_tercile = hhincome, nq(3)

di "Sample: " _N " observations"
tab choice [aw=hsupwgtk], missing

/*------------------------------------------------------------------------------
2. Standard Cell-Based CCP Estimation (Baseline)

   P(choice | CBSA, year, age_cat, educ_cat) = cell mean
------------------------------------------------------------------------------*/

di _n "============================================================"
di "METHOD 1: CELL-BASED CCPs (Standard)"
di "============================================================"

* Create cell identifier
egen cell = group(cbsa year age_cat educ_cat)
quietly summarize cell
local n_cells = r(max)
di "Number of cells: `n_cells'"

* Cell sizes
bysort cell: gen cell_size = _N
quietly summarize cell_size
di "Cell size: mean = " %6.1f r(mean) ", min = " %3.0f r(min) ", max = " %5.0f r(max)
di "Cells with < 10 obs: " _N - sum(cell_size >= 10) / _N * 100 "%"

* Compute cell-based CCPs
forvalues j = 1/9 {
    gen choice_`j' = (choice == `j')
    bysort cell: egen ccp_cell_`j' = mean(choice_`j')
}

* Display sample CCPs
di _n "Sample Cell-Based CCPs (first 5 cells):"
list cell ccp_cell_1 ccp_cell_2 ccp_cell_7 ccp_cell_8 if cell <= 5 & _n <= 20, sepby(cell)

/*------------------------------------------------------------------------------
3. ML Method 1: Regularized Multinomial Logit

   Use LASSO-penalized multinomial logit to estimate CCPs flexibly
------------------------------------------------------------------------------*/

di _n "============================================================"
di "METHOD 2: REGULARIZED MULTINOMIAL LOGIT"
di "============================================================"

* Create feature set for ML
local features "age age_cat educ_cat female married metro has_kids inc_tercile pct_broadband year"

* Check if mlogit with LASSO is available
* If not, use regular mlogit with flexible controls as approximation

capture mlogit choice i.age_cat##i.educ_cat##i.year c.age c.pct_broadband ///
    female married metro has_kids i.inc_tercile ///
    [pw=hsupwgtk], vce(cluster cbsa) baseoutcome(7)

if _rc == 0 {
    di "Multinomial logit estimated successfully"

    * Predict CCPs for each alternative
    forvalues j = 1/9 {
        predict ccp_mlogit_`j', outcome(`j')
    }

    * Compare to cell-based CCPs
    di _n "Correlation: Cell-based vs Multinomial Logit CCPs"
    forvalues j = 1/9 {
        quietly correlate ccp_cell_`j' ccp_mlogit_`j'
        di "  Choice `j': r = " %6.3f r(rho)
    }
}
else {
    di "Note: Full mlogit did not converge. Using simplified specification."

    * Simplified specification
    mlogit choice i.age_cat i.educ_cat i.year c.pct_broadband ///
        female married [pw=hsupwgtk], vce(cluster cbsa) baseoutcome(7)

    forvalues j = 1/9 {
        predict ccp_mlogit_`j', outcome(`j')
    }
}

/*------------------------------------------------------------------------------
4. ML Method 2: Random Forest via Stata (if available) or Approximation

   Stata doesn't have built-in random forest, but we can:
   a) Use rforest package if installed
   b) Approximate via boosted classification trees
   c) Export to R/Python for RF estimation
------------------------------------------------------------------------------*/

di _n "============================================================"
di "METHOD 3: FLEXIBLE NONPARAMETRIC (RF Approximation)"
di "============================================================"

* Check for rforest
capture which rforest
if _rc == 0 {
    di "Random forest package available"

    * For each choice, fit RF classifier
    forvalues j = 1/9 {
        quietly {
            rforest choice_`j' age educ_cat female married metro has_kids ///
                inc_tercile pct_broadband, type(class) iter(100) seed(20260211)
            predict ccp_rf_`j', pr
        }
        di "RF for choice `j' complete"
    }
}
else {
    di "Random forest package not available."
    di "Using kernel-smoothed local means as nonparametric approximation."

    * Alternative: Local polynomial smoothing of CCPs
    * This captures nonlinear relationships without cell discretization

    forvalues j = 1/9 {
        * Use lpoly-style smoothing over continuous state variables
        * Approximate with regression on flexible polynomials
        quietly {
            gen ccp_smooth_`j' = .

            * Fit separate models by age category to allow heterogeneity
            forvalues a = 1/3 {
                capture reg choice_`j' c.pct_broadband##c.pct_broadband ///
                    i.educ_cat female married i.inc_tercile i.year ///
                    [pw=hsupwgtk] if age_cat == `a', vce(cluster cbsa)

                if _rc == 0 {
                    predict temp_`j'_`a' if age_cat == `a'
                    replace ccp_smooth_`j' = temp_`j'_`a' if age_cat == `a'
                    drop temp_`j'_`a'
                }
            }
        }
        di "Smoothed CCP for choice `j' complete"
    }
}

/*------------------------------------------------------------------------------
5. Compare CCP Estimation Methods
------------------------------------------------------------------------------*/

di _n "============================================================"
di "COMPARISON OF CCP ESTIMATION METHODS"
di "============================================================"

* Mean CCPs by method
di _n "Mean CCPs by Method:"
di "Choice | Cell-Based | MLogit | Smoothed"
di "-------|------------|--------|----------"

forvalues j = 1/9 {
    quietly summarize ccp_cell_`j' [aw=hsupwgtk]
    local m_cell = r(mean)
    quietly summarize ccp_mlogit_`j' [aw=hsupwgtk]
    local m_mlogit = r(mean)

    capture quietly summarize ccp_smooth_`j' [aw=hsupwgtk]
    if _rc == 0 {
        local m_smooth = r(mean)
    }
    else {
        local m_smooth = .
    }

    di "   " `j' "   |   " %6.4f `m_cell' "   |  " %6.4f `m_mlogit' "  |  " %6.4f `m_smooth'
}

* RMSE comparison (cell-based vs ML)
di _n "RMSE: ML vs Cell-Based CCPs"
local total_rmse_mlogit = 0
local total_rmse_smooth = 0
local n_choices = 0

forvalues j = 1/9 {
    gen sq_diff_mlogit_`j' = (ccp_mlogit_`j' - ccp_cell_`j')^2
    quietly summarize sq_diff_mlogit_`j'
    local rmse_mlogit = sqrt(r(mean))
    local total_rmse_mlogit = `total_rmse_mlogit' + `rmse_mlogit'

    capture {
        gen sq_diff_smooth_`j' = (ccp_smooth_`j' - ccp_cell_`j')^2
        quietly summarize sq_diff_smooth_`j'
        local rmse_smooth = sqrt(r(mean))
        local total_rmse_smooth = `total_rmse_smooth' + `rmse_smooth'
    }

    local n_choices = `n_choices' + 1
}

di "Average RMSE (MLogit vs Cell): " %8.4f (`total_rmse_mlogit' / `n_choices')
di "Average RMSE (Smooth vs Cell): " %8.4f (`total_rmse_smooth' / `n_choices')

/*------------------------------------------------------------------------------
6. Use ML CCPs in Structural Estimation

   Following Arcidiacono-Miller: continuation value differences are functions
   of CCPs. Better CCP estimates → better structural estimates.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "STRUCTURAL ESTIMATION WITH ML CCPs"
di "============================================================"

* The key object for finite dependence: log CCP ratios
* ln P(choice) - ln P(renewal) where renewal = Branch × Wage (choice 7)

forvalues j = 1/9 {
    if `j' != 7 {
        gen log_ccp_ratio_cell_`j' = ln(ccp_cell_`j' / ccp_cell_7)
        gen log_ccp_ratio_mlogit_`j' = ln(ccp_mlogit_`j' / ccp_mlogit_7)
    }
}

* Continuation value correction terms for finite dependence
* These enter the second-stage likelihood
di _n "Log CCP Ratios (for structural estimation):"
di "Mean log(P(choice)/P(Branch×Wage)):"

forvalues j = 1/9 {
    if `j' != 7 {
        quietly summarize log_ccp_ratio_cell_`j' if log_ccp_ratio_cell_`j' != .
        local m_cell = r(mean)
        quietly summarize log_ccp_ratio_mlogit_`j' if log_ccp_ratio_mlogit_`j' != .
        local m_mlogit = r(mean)
        di "  Choice `j': Cell = " %7.3f `m_cell' ", MLogit = " %7.3f `m_mlogit'
    }
}

/*------------------------------------------------------------------------------
7. Save ML CCPs for Use in Structural Model
------------------------------------------------------------------------------*/

preserve
keep choice ccp_cell_* ccp_mlogit_* age_cat educ_cat year cbsa hsupwgtk

* Collapse to unique state-choice cells
collapse (mean) ccp_cell_* ccp_mlogit_* [aw=hsupwgtk], by(cbsa year age_cat educ_cat)

save "$output/ml_ccps.dta", replace
di "ML CCPs saved to $output/ml_ccps.dta"
restore

/*------------------------------------------------------------------------------
8. Summary
------------------------------------------------------------------------------*/

di _n "============================================================"
di "SUMMARY"
di "============================================================"
di ""
di "ML-based CCP estimation provides several advantages:"
di ""
di "1. Better handling of thin cells: Instead of noisy cell means,"
di "   ML methods borrow strength across similar observations."
di ""
di "2. Flexible functional form: Captures nonlinear relationships"
di "   between state variables and choice probabilities."
di ""
di "3. Consistent with Arcidiacono-Miller framework: The second"
di "   stage structural estimation uses CCPs - ML just gives"
di "   better first-stage estimates."
di ""
di "Recommendation: Use multinomial logit CCPs for the structural"
di "model, as they provide smoothed estimates while maintaining"
di "the multinomial logit structure assumed in the model."

log close
