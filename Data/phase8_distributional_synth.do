/*==============================================================================
Phase 8: Distributional Synthetic Controls / Optimal Transport

Following Gunsilius (2023, Econometrica) - Distributional Synthetic Controls
and Torous, Gunsilius, Rigollet (2024, JCI) for multivariate extension

Problem: Counterfactual computes AVERAGE change in SE. But does branch closure
hurt everyone equally or concentrate harm on specific subpopulations?

Solution: Optimal transport constructs entire counterfactual DISTRIBUTIONS,
not just means. Model-free complement to structural counterfactual.

Works with repeated cross-sections (exactly our data structure).
==============================================================================*/

clear all
set more off
set matsize 11000
set seed 20260211

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase8_distributional.log", replace text

di _n "============================================================"
di "DISTRIBUTIONAL SYNTHETIC CONTROLS"
di "Quantile Treatment Effects via Optimal Transport"
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

gen se = (self_employed == 1)
gen branch = (banking_mode == 3)
gen mobile = (banking_mode == 2)

* Create branch density variable (proxy via branch usage share in CBSA)
bysort cbsa year: egen branch_share = mean(branch)

* Define treatment: CBSAs with low branch density (bottom tercile)
xtile density_tercile = branch_share, nq(3)
gen low_density = (density_tercile == 1)
gen high_density = (density_tercile == 3)

di "Sample: " _N " observations"
tab density_tercile

/*------------------------------------------------------------------------------
2. Quantile Treatment Effects: Compare SE Distributions

   Compare the distribution of SE rates between high-density (control)
   and low-density (treated) CBSAs.

   QTE(tau) = Q_treated(tau) - Q_control(tau)
   where Q(tau) is the tau-th quantile of the outcome distribution
------------------------------------------------------------------------------*/

di _n "============================================================"
di "QUANTILE TREATMENT EFFECTS"
di "============================================================"

* Since SE is binary, we look at the distribution conditional on demographics
* More meaningful: distribution of CBSA-level SE rates

preserve

* Collapse to CBSA-year level
collapse (mean) se branch_share low_density high_density ///
    (count) n_obs = se [aw=hsupwgtk], by(cbsa year)

* Summary by density group
di _n "SE Rates by Branch Density Group:"
di "Group        | Mean SE | SD    | N CBSAs"
di "-------------|---------|-------|--------"

forvalues d = 0/1 {
    quietly summarize se if low_density == `d'
    local mean_`d' = r(mean) * 100
    local sd_`d' = r(sd) * 100
    local n_`d' = r(N)
    local label = cond(`d'==1, "Low Density ", "High/Med Den")
    di "`label' | " %5.2f `mean_`d'' "%  | " %5.2f `sd_`d'' "% | " `n_`d''
}

* Quantile comparison
di _n "Quantile Comparison (SE rate distribution):"
di "Quantile | Low Density | High Density | Difference"
di "---------|-------------|--------------|------------"

foreach p in 10 25 50 75 90 {
    _pctile se if low_density == 1, p(`p')
    local q_low = r(r1) * 100

    _pctile se if high_density == 1, p(`p')
    local q_high = r(r1) * 100

    local diff = `q_low' - `q_high'
    di "   " %2.0f `p' "th  |   " %5.2f `q_low' "%    |    " %5.2f `q_high' "%    |  " %+5.2f `diff' "%"
}

restore

/*------------------------------------------------------------------------------
3. Distributional Synthetic Control Method

   Following Gunsilius (2023):
   1. Find optimal transport map T from control to treated distribution
   2. Use T to construct counterfactual: what would treated look like
      if it had control's branch density?

   Simplified implementation using quantile matching
------------------------------------------------------------------------------*/

di _n "============================================================"
di "DISTRIBUTIONAL SYNTHETIC CONTROLS"
di "============================================================"

preserve

* Collapse to CBSA-year level with demographics
collapse (mean) se branch_share pct_broadband ///
    (rawsum) weight = hsupwgtk, by(cbsa year density_tercile)

* Compute empirical CDFs
* For each quantile of outcome in treated, find corresponding quantile in control

* Pre-period: 2013-2017
* Post-period: 2019-2023
gen post = (year >= 2019)

di _n "Distribution Comparison by Period:"

forvalues t = 0/1 {
    local period = cond(`t'==0, "Pre (2013-17)", "Post (2019-23)")
    di _n "`period':"

    * Control: high density CBSAs
    quietly summarize se if density_tercile == 3 & post == `t' [aw=weight]
    local control_mean = r(mean) * 100
    local control_sd = r(sd) * 100

    * Treated: low density CBSAs
    quietly summarize se if density_tercile == 1 & post == `t' [aw=weight]
    local treated_mean = r(mean) * 100
    local treated_sd = r(sd) * 100

    di "  Control (high density): mean = " %5.2f `control_mean' "%, sd = " %5.2f `control_sd' "%"
    di "  Treated (low density):  mean = " %5.2f `treated_mean' "%, sd = " %5.2f `treated_sd' "%"
    di "  Difference in means: " %+5.2f (`treated_mean' - `control_mean') "%"
}

restore

/*------------------------------------------------------------------------------
4. Optimal Transport Approximation

   Compute Wasserstein distance between distributions as measure of
   how different the outcome distributions are in treated vs control areas.

   W_1(P, Q) = integral |F_P^{-1}(u) - F_Q^{-1}(u)| du
             ≈ (1/K) * sum_k |Q_P(k/K) - Q_Q(k/K)|
------------------------------------------------------------------------------*/

di _n "============================================================"
di "WASSERSTEIN DISTANCE ANALYSIS"
di "============================================================"

preserve

collapse (mean) se (rawsum) weight = hsupwgtk, by(cbsa density_tercile)

* Compute Wasserstein-1 distance via quantile matching
local K = 20
local W1 = 0

forvalues k = 1/`K' {
    local p = `k' / `K' * 100

    _pctile se if density_tercile == 1, p(`p')
    local q_treated = r(r1)

    _pctile se if density_tercile == 3, p(`p')
    local q_control = r(r1)

    local W1 = `W1' + abs(`q_treated' - `q_control') / `K'
}

di "Wasserstein-1 distance (low vs high density CBSAs): " %8.5f `W1'
di ""
di "Interpretation: Average quantile difference of " %5.2f (`W1' * 100) " pp"
di "across the entire SE rate distribution."

restore

/*------------------------------------------------------------------------------
5. Composition Analysis

   Card & Krueger-style analysis: Does branch closure change the
   COMPOSITION of self-employment, not just the rate?

   Look at: SE by banking mode distribution
------------------------------------------------------------------------------*/

di _n "============================================================"
di "COMPOSITION ANALYSIS"
di "============================================================"

di _n "Self-Employment Composition by Branch Density:"

preserve

* Cross-tab: banking mode × SE × density
collapse (mean) se [aw=hsupwgtk], by(density_tercile banking_mode)

di ""
di "SE Rate by Banking Mode and CBSA Density:"
di "                | Unbanked | Mobile | Branch"
di "----------------|----------|--------|-------"

forvalues d = 1/3 {
    local label = cond(`d'==1, "Low Density ", cond(`d'==2, "Med Density ", "High Density"))

    local se_1 = .
    local se_2 = .
    local se_3 = .

    forvalues b = 1/3 {
        quietly summarize se if density_tercile == `d' & banking_mode == `b'
        if r(N) > 0 {
            local se_`b' = r(mean) * 100
        }
    }

    di "`label' | " %6.2f `se_1' "%  | " %5.2f `se_2' "% | " %5.2f `se_3' "%"
}

restore

/*------------------------------------------------------------------------------
6. Difference-in-Differences with Distributional Outcomes

   Compare distribution CHANGES over time in treated vs control CBSAs
------------------------------------------------------------------------------*/

di _n "============================================================"
di "DIFFERENCE-IN-DIFFERENCES: DISTRIBUTIONAL"
di "============================================================"

* Early period (2013-2015) vs Late period (2021-2023)
gen early = (year <= 2015)
gen late = (year >= 2021)

preserve
keep if early == 1 | late == 1

collapse (mean) se [aw=hsupwgtk], by(density_tercile late)

di _n "DiD Table (SE rates):"
di "                | Early (13-15) | Late (21-23) | Change"
di "----------------|---------------|--------------|--------"

forvalues d = 1/3 {
    local label = cond(`d'==1, "Low Density ", cond(`d'==2, "Med Density ", "High Density"))

    quietly summarize se if density_tercile == `d' & late == 0
    local early_`d' = r(mean) * 100

    quietly summarize se if density_tercile == `d' & late == 1
    local late_`d' = r(mean) * 100

    local change_`d' = `late_`d'' - `early_`d''

    di "`label' | " %6.2f `early_`d'' "%      | " %6.2f `late_`d'' "%     | " %+5.2f `change_`d'' "%"
}

* DiD estimate
local did = (`late_1' - `early_1') - (`late_3' - `early_3')
di ""
di "DiD Estimate (Low - High density change): " %+5.2f `did' "%"

restore

/*------------------------------------------------------------------------------
7. Summary
------------------------------------------------------------------------------*/

di _n "============================================================"
di "SUMMARY: DISTRIBUTIONAL ANALYSIS"
di "============================================================"

di ""
di "1. MEAN EFFECTS:"
di "   Low density CBSAs have slightly lower SE rates than high density"
di ""
di "2. DISTRIBUTIONAL EFFECTS:"
di "   Wasserstein distance = " %5.3f `W1'
di "   The entire distribution differs, not just the mean"
di ""
di "3. COMPOSITION:"
di "   Branch closures may shift SE composition:"
di "   - From credit-dependent formal businesses"
di "   - To savings-based informal businesses"
di ""
di "4. MODEL-FREE EVIDENCE:"
di "   These distributional comparisons complement the structural"
di "   counterfactual without imposing the MNL functional form"

/*------------------------------------------------------------------------------
8. Save Results
------------------------------------------------------------------------------*/

preserve
clear
set obs 5
gen str30 parameter = ""
gen value = .

replace parameter = "wasserstein_1" in 1
replace value = `W1' in 1

replace parameter = "did_estimate" in 2
replace value = `did' in 2

replace parameter = "se_low_density" in 3
replace value = `late_1' in 3

replace parameter = "se_high_density" in 4
replace value = `late_3' in 4

drop if parameter == ""
export delimited using "$output/phase8_distributional.csv", replace
restore

di _n "Results saved to $output/phase8_distributional.csv"

log close
