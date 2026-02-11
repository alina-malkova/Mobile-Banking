/*==============================================================================
Phase 2: Panel BIC Model Selection for Unobserved Heterogeneity
Following Hao & Kasahara (2025) for panel finite mixture models

Three-pronged approach:
1. Hao-Kasahara (2025) panel BIC - primary method
2. Bonhomme-Lamadon-Manresa (2022) counterfactual stability
3. Budanova (2025) OSCE penalized MLE - robustness check
==============================================================================*/

clear all
set more off
set matsize 11000
set seed 20260211

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase2_panel_bic.log", replace text

di _n "============================================================"
di "PANEL BIC MODEL SELECTION (HAO-KASAHARA 2025)"
di "Three-Pronged Approach to Type Selection"
di "============================================================"

/*------------------------------------------------------------------------------
1. Data Preparation
------------------------------------------------------------------------------*/

use "$datadir/analysis_dataset_with_se.dta", clear

keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .
keep if banking_mode != .

* Employment status
capture drop emp_status
gen emp_status = .
replace emp_status = 1 if wage_worker == 1
replace emp_status = 2 if self_employed == 1
replace emp_status = 3 if employed != 1

* Banking mode
capture drop bank_mode
gen bank_mode = banking_mode

* Age categories
gen age_cat = 1 if age >= 18 & age < 35
replace age_cat = 2 if age >= 35 & age < 50
replace age_cat = 3 if age >= 50 & age <= 64

* Education categories
gen educ_cat = .
replace educ_cat = 1 if no_hs == 1
replace educ_cat = 2 if hs_diploma == 1
replace educ_cat = 3 if some_college == 1
replace educ_cat = 4 if college_degree == 1

* SE indicator
gen se = (self_employed == 1)

* Branch user indicator
gen branch = (bank_mode == 3)

* Mobile user indicator
gen mobile = (bank_mode == 2)

* Count panels
egen cbsa_id = group(cbsa)
quietly summarize cbsa_id
local N_panels = r(max)
di "Number of panels (CBSAs): `N_panels'"

quietly count
local N_obs = r(N)
di "Total observations: `N_obs'"

* Panel structure
local T_periods = 6  // Survey waves
local c_T = 1 + 1/`T_periods'  // Hao-Kasahara correction factor
di "Panel correction factor c(T) = " %5.3f `c_T'

/*------------------------------------------------------------------------------
2. Construct Variables for Estimation
------------------------------------------------------------------------------*/

* Lagged SE rate by CBSA
bysort cbsa year: egen cbsa_se_rate = mean(se)
bysort cbsa (year): gen lag_se_rate = cbsa_se_rate[_n-1]
replace lag_se_rate = 0.11 if lag_se_rate == .

* Dynamic SE return interaction
gen dynamic_se = se * lag_se_rate

* Switching cost proxy
gen switch_cost = (bank_mode != 3) * branch

* Observable-based type scores for initialization
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
sum type_rank
local max_rank = r(max)

/*------------------------------------------------------------------------------
METHOD 1: Hao-Kasahara (2025) Panel BIC
Using Linear Probability Model for computational stability
Panel BIC = -2*ln(L) + p * ln(N_panels) * c(T)
==============================================================================*/

di _n "=============================================================="
di "METHOD 1: Hao-Kasahara (2025) Panel BIC"
di "=============================================================="

* Baseline SE rate
quietly summarize se [aw=hsupwgtk]
local baseline_se = r(mean)
di "Baseline self-employment rate: " %6.4f `baseline_se'

* Store results
matrix HK_results = J(4, 6, .)
matrix colnames HK_results = K logL params N_obs BIC_standard BIC_panel
matrix rownames HK_results = K1 K2 K3 K4

forvalues K = 1/4 {
    di _n "--- Estimating K = `K' types ---"

    preserve

    if `K' == 1 {
        * Homogeneous model (LPM)
        reg se i.age_cat i.educ_cat branch mobile ///
            dynamic_se switch_cost i.year [aw=hsupwgtk], vce(cluster cbsa)

        local ll = e(ll)
        local params = e(df_m) + 1
    }
    else {
        * Assign types based on type_rank percentiles
        if `K' == 2 {
            gen type = 1 if type_rank <= `max_rank'/2
            replace type = 2 if type == .
        }
        else if `K' == 3 {
            gen type = 1 if type_rank <= `max_rank'/3
            replace type = 2 if type_rank > `max_rank'/3 & type_rank <= 2*`max_rank'/3
            replace type = 3 if type == .
        }
        else if `K' == 4 {
            gen type = 1 if type_rank <= `max_rank'/4
            replace type = 2 if type_rank > `max_rank'/4 & type_rank <= `max_rank'/2
            replace type = 3 if type_rank > `max_rank'/2 & type_rank <= 3*`max_rank'/4
            replace type = 4 if type == .
        }

        * Create type-specific variables
        forvalues k = 1/`K' {
            gen tau_`k' = (type == `k')
            gen branch_`k' = branch * tau_`k'
            gen dyn_`k' = dynamic_se * tau_`k'
        }

        * Type-specific estimation (LPM)
        local br_vars ""
        local dyn_vars ""
        forvalues k = 1/`K' {
            local br_vars "`br_vars' branch_`k'"
            local dyn_vars "`dyn_vars' dyn_`k'"
        }

        reg se i.age_cat i.educ_cat mobile ///
            `br_vars' `dyn_vars' switch_cost i.year [aw=hsupwgtk], vce(cluster cbsa)

        local ll = e(ll)
        local params = e(df_m) + 1 + (`K' - 1)
    }

    restore

    * Compute BIC variants
    * Standard BIC uses total observations
    local BIC_std = -2 * `ll' + `params' * ln(`N_obs')
    * Panel BIC (Hao-Kasahara) uses number of panels with correction
    local BIC_panel = -2 * `ll' + `params' * ln(`N_panels') * `c_T'

    matrix HK_results[`K', 1] = `K'
    matrix HK_results[`K', 2] = `ll'
    matrix HK_results[`K', 3] = `params'
    matrix HK_results[`K', 4] = `N_obs'
    matrix HK_results[`K', 5] = `BIC_std'
    matrix HK_results[`K', 6] = `BIC_panel'

    di "K = `K': Log-L = " %12.2f `ll' ", params = " %3.0f `params'
    di "        Standard BIC = " %12.2f `BIC_std'
    di "        Panel BIC (Hao-Kasahara) = " %12.2f `BIC_panel'
}

di _n "Hao-Kasahara (2025) Panel BIC Summary:"
di "========================================"
matrix list HK_results, format(%12.2f)

* Find optimal K
local min_std_bic = HK_results[1, 5]
local min_panel_bic = HK_results[1, 6]
local optimal_K_std = 1
local optimal_K_panel = 1

forvalues K = 2/4 {
    if HK_results[`K', 5] < `min_std_bic' {
        local min_std_bic = HK_results[`K', 5]
        local optimal_K_std = `K'
    }
    if HK_results[`K', 6] < `min_panel_bic' {
        local min_panel_bic = HK_results[`K', 6]
        local optimal_K_panel = `K'
    }
}

di _n "Standard BIC selects: K = `optimal_K_std'"
di "Hao-Kasahara Panel BIC selects: K = `optimal_K_panel'"

/*------------------------------------------------------------------------------
METHOD 2: Bonhomme-Lamadon-Manresa (2022) Counterfactual Stability
==============================================================================*/

di _n "=============================================================="
di "METHOD 2: Bonhomme-Lamadon-Manresa (2022) Counterfactual Stability"
di "=============================================================="

matrix CF_results = J(4, 4, .)
matrix colnames CF_results = K CF_50pct CF_25pct CF_75pct
matrix rownames CF_results = K1 K2 K3 K4

forvalues K = 1/4 {
    di _n "--- Counterfactual for K = `K' ---"

    preserve

    if `K' == 1 {
        quietly reg se i.age_cat i.educ_cat branch mobile ///
            dynamic_se switch_cost i.year [aw=hsupwgtk], vce(cluster cbsa)

        local b_branch = _b[branch]
        local cf_50 = `b_branch' * 0.5 / `baseline_se' * 100
        local cf_25 = `b_branch' * 0.25 / `baseline_se' * 100
        local cf_75 = `b_branch' * 0.75 / `baseline_se' * 100

        di "Branch effect (beta) = " %8.5f `b_branch'
    }
    else {
        * Assign types
        if `K' == 2 {
            gen type = 1 if type_rank <= `max_rank'/2
            replace type = 2 if type == .
        }
        else if `K' == 3 {
            gen type = 1 if type_rank <= `max_rank'/3
            replace type = 2 if type_rank > `max_rank'/3 & type_rank <= 2*`max_rank'/3
            replace type = 3 if type == .
        }
        else if `K' == 4 {
            gen type = 1 if type_rank <= `max_rank'/4
            replace type = 2 if type_rank > `max_rank'/4 & type_rank <= `max_rank'/2
            replace type = 3 if type_rank > `max_rank'/2 & type_rank <= 3*`max_rank'/4
            replace type = 4 if type == .
        }

        * Type-specific variables
        forvalues k = 1/`K' {
            gen tau_`k' = (type == `k')
            gen branch_`k' = branch * tau_`k'
            gen dyn_`k' = dynamic_se * tau_`k'
        }

        local br_vars ""
        local dyn_vars ""
        forvalues k = 1/`K' {
            local br_vars "`br_vars' branch_`k'"
            local dyn_vars "`dyn_vars' dyn_`k'"
        }

        quietly reg se i.age_cat i.educ_cat mobile ///
            `br_vars' `dyn_vars' switch_cost i.year [aw=hsupwgtk], vce(cluster cbsa)

        * Compute weighted counterfactual
        local cf_50 = 0
        local cf_25 = 0
        local cf_75 = 0

        forvalues k = 1/`K' {
            quietly summarize tau_`k' [aw=hsupwgtk]
            local share_k = r(mean)
            local b_k = _b[branch_`k']

            local cf_50 = `cf_50' + `share_k' * `b_k' * 0.5 / `baseline_se' * 100
            local cf_25 = `cf_25' + `share_k' * `b_k' * 0.25 / `baseline_se' * 100
            local cf_75 = `cf_75' + `share_k' * `b_k' * 0.75 / `baseline_se' * 100

            di "Type `k': share = " %5.3f `share_k' ", beta = " %8.5f `b_k'
        }
    }

    restore

    matrix CF_results[`K', 1] = `K'
    matrix CF_results[`K', 2] = `cf_50'
    matrix CF_results[`K', 3] = `cf_25'
    matrix CF_results[`K', 4] = `cf_75'

    di "K = `K': CF effect (50% closure) = " %6.2f `cf_50' "%"
}

di _n "Counterfactual Stability Analysis:"
di "==================================="
matrix list CF_results, format(%8.2f)

* Check stability
di _n "Stability check (|change| < 2pp suggests stability):"
forvalues K = 2/4 {
    local prev = `K' - 1
    local change = CF_results[`K', 2] - CF_results[`prev', 2]
    local stable = "STABLE"
    if abs(`change') > 2 {
        local stable = "UNSTABLE"
    }
    di "K=`prev' to K=`K': change = " %6.2f `change' " pp  [`stable']"
}

/*------------------------------------------------------------------------------
METHOD 3: Budanova (2025) OSCE Approximation
Start with K=5, identify active types via significance
==============================================================================*/

di _n "=============================================================="
di "METHOD 3: Budanova (2025) OSCE Approximation"
di "=============================================================="

preserve

* Overspecified model with K=5
local K_max = 5
gen type = 1 if type_rank <= `max_rank'/5
replace type = 2 if type_rank > `max_rank'/5 & type_rank <= 2*`max_rank'/5
replace type = 3 if type_rank > 2*`max_rank'/5 & type_rank <= 3*`max_rank'/5
replace type = 4 if type_rank > 3*`max_rank'/5 & type_rank <= 4*`max_rank'/5
replace type = 5 if type == .

forvalues k = 1/`K_max' {
    gen tau_`k' = (type == `k')
    gen branch_`k' = branch * tau_`k'
    gen dyn_`k' = dynamic_se * tau_`k'
}

local br_vars ""
local dyn_vars ""
forvalues k = 1/`K_max' {
    local br_vars "`br_vars' branch_`k'"
    local dyn_vars "`dyn_vars' dyn_`k'"
}

quietly reg se i.age_cat i.educ_cat mobile ///
    `br_vars' `dyn_vars' switch_cost i.year [aw=hsupwgtk], vce(cluster cbsa)

di "Overspecified model (K=5) - Type-specific branch effects:"
di "=========================================================="

local active_types = 0
forvalues k = 1/`K_max' {
    quietly summarize tau_`k' [aw=hsupwgtk]
    local share = r(mean)
    local b = _b[branch_`k']
    local se_coef = _se[branch_`k']
    local t = `b' / `se_coef'
    local p = 2 * (1 - normal(abs(`t')))

    local status = "shrink"
    if abs(`t') > 1.5 {
        local status = "ACTIVE"
        local active_types = `active_types' + 1
    }

    di "Type `k': share=" %5.3f `share' ", beta=" %7.4f `b' ", t=" %6.2f `t' ", p=" %5.3f `p' "  [`status']"
}

di _n "OSCE approximation suggests K = `active_types' active types"

restore

/*------------------------------------------------------------------------------
SUMMARY AND CONSENSUS
==============================================================================*/

di _n "=============================================================="
di "SUMMARY: Three-Pronged Model Selection"
di "=============================================================="

di _n "1. Hao-Kasahara (2025) Panel BIC selects: K = `optimal_K_panel'"
di "   (Standard BIC selects: K = `optimal_K_std')"
di ""
di "2. Bonhomme-Lamadon-Manresa (2022) Counterfactual Stability:"
di "   See stability check above"
di ""
di "3. Budanova (2025) OSCE approximation: K = `active_types'"
di ""

* Determine consensus
local consensus = `optimal_K_panel'
if `optimal_K_panel' == `optimal_K_std' {
    di "CONSENSUS: K = `consensus' types (Panel BIC and Standard BIC agree)"
}
else {
    di "NOTE: Panel BIC suggests K=`optimal_K_panel', Standard BIC suggests K=`optimal_K_std'"
    di "      Following Hao-Kasahara (2025), prefer Panel BIC for panel data"
}

/*------------------------------------------------------------------------------
Save Results
==============================================================================*/

preserve
clear
set obs 20
gen str20 method = ""
gen str20 item = ""
gen value = .
gen K = .

local row = 1

* Panel BIC results
forvalues k = 1/4 {
    replace method = "Hao-Kasahara" in `row'
    replace item = "Panel_BIC" in `row'
    replace value = HK_results[`k', 6] in `row'
    replace K = `k' in `row'
    local row = `row' + 1
}

* Standard BIC for comparison
forvalues k = 1/4 {
    replace method = "Standard" in `row'
    replace item = "BIC" in `row'
    replace value = HK_results[`k', 5] in `row'
    replace K = `k' in `row'
    local row = `row' + 1
}

* Counterfactuals
forvalues k = 1/4 {
    replace method = "BLM" in `row'
    replace item = "CF_50pct" in `row'
    replace value = CF_results[`k', 2] in `row'
    replace K = `k' in `row'
    local row = `row' + 1
}

* Selection results
replace method = "Selected" in `row'
replace item = "optimal_K_panel" in `row'
replace value = `optimal_K_panel' in `row'
local row = `row' + 1

replace method = "Selected" in `row'
replace item = "optimal_K_std" in `row'
replace value = `optimal_K_std' in `row'
local row = `row' + 1

replace method = "OSCE" in `row'
replace item = "active_types" in `row'
replace value = `active_types' in `row'
local row = `row' + 1

replace method = "Baseline" in `row'
replace item = "se_rate" in `row'
replace value = `baseline_se' in `row'

drop if method == ""
export delimited using "$output/phase2_panel_bic.csv", replace

di _n "Results saved to $output/phase2_panel_bic.csv"
restore

log close
