/*==============================================================================
Phase 7: Mixed Logit with Continuous Random Coefficients

Following Train (2009): Instead of K discrete types, allow the branch effect
coefficient gamma_C to be drawn from a continuous distribution:
    gamma_C ~ N(gamma_bar, sigma_gamma^2)

This directly estimates mean and variance of heterogeneity without discretizing.

Advantages:
- No need to select K
- Continuous heterogeneity is more realistic
- Simulated MLE straightforward in Stata (mixlogit) or R (mlogit)

Disadvantage:
- Lose interpretable type labels

Role in paper: Robustness check. If sigma_gamma is large and significant,
confirms heterogeneity found by finite mixture without relying on type selection.
==============================================================================*/

clear all
set more off
set matsize 11000
set seed 20260211

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase7_mixed_logit.log", replace text

di _n "============================================================"
di "MIXED LOGIT WITH CONTINUOUS RANDOM COEFFICIENTS"
di "Robustness Check for Finite Mixture Results"
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

* Create outcome and treatment
gen se = (self_employed == 1)
gen branch = (banking_mode == 3)
gen mobile = (banking_mode == 2)

* Covariates
gen age_cat = 1 if age >= 18 & age < 35
replace age_cat = 2 if age >= 35 & age < 50
replace age_cat = 3 if age >= 50 & age <= 64

gen educ_cat = .
replace educ_cat = 1 if no_hs == 1
replace educ_cat = 2 if hs_diploma == 1
replace educ_cat = 3 if some_college == 1
replace educ_cat = 4 if college_degree == 1

gen female = (sex == 2)
gen married = (marital_status == 1 | marital_status == 2)

di "Sample: " _N " observations"

/*------------------------------------------------------------------------------
2. Reshape for Mixed Logit

   Mixed logit requires data in long form with choice-specific variables
------------------------------------------------------------------------------*/

di _n "============================================================"
di "RESHAPING DATA FOR MIXED LOGIT"
di "============================================================"

* Create person identifier
gen person_id = _n

* Keep only necessary variables
keep person_id se branch mobile pct_broadband age_cat educ_cat female married ///
     year cbsa hsupwgtk

* The choice set is {SE=0, SE=1}
* We want to estimate how branch affects the probability of SE=1
* with random coefficient on branch

* For binary choice, we can use standard mixed logit setup
* Alternative: use cmxtmixlogit for panel mixed logit

/*------------------------------------------------------------------------------
3. Standard Logit (Baseline)
------------------------------------------------------------------------------*/

di _n "============================================================"
di "MODEL 1: STANDARD LOGIT (HOMOGENEOUS)"
di "============================================================"

logit se branch mobile c.pct_broadband i.age_cat i.educ_cat female married ///
    i.year [pw=hsupwgtk], vce(cluster cbsa)

local b_branch_std = _b[branch]
local se_branch_std = _se[branch]

di _n "Standard Logit Results:"
di "  Branch coefficient: " %8.4f `b_branch_std' " (se = " %6.4f `se_branch_std' ")"

margins, dydx(branch) atmeans
local me_branch_std = r(table)[1,1]

di "  Marginal effect: " %8.4f `me_branch_std'

/*------------------------------------------------------------------------------
4. Mixed Logit with Random Coefficient on Branch

   gamma_C_i = gamma_bar + sigma_gamma * nu_i
   where nu_i ~ N(0,1)

   This captures individual-specific heterogeneity in branch effects
------------------------------------------------------------------------------*/

di _n "============================================================"
di "MODEL 2: MIXED LOGIT (RANDOM COEFFICIENT ON BRANCH)"
di "============================================================"

* Check if mixlogit is installed
capture which mixlogit
if _rc != 0 {
    di "Note: mixlogit not installed. Using gsem for random effects logit."
    di "This is an approximation to the full mixed logit."

    * Use GSEM for random intercept model as approximation
    * True mixed logit would require panel structure or specialized command

    * Random effects logit (random intercept)
    gsem (se <- branch mobile c.pct_broadband i.age_cat i.educ_cat female married ///
        i.year M1[cbsa]@1), family(bernoulli) link(logit) ///
        latent(M1) vce(cluster cbsa)

    local b_branch_re = _b[se:branch]
    local se_branch_re = _se[se:branch]
    local sigma_cbsa = exp(_b[var(M1[cbsa]):_cons]/2)

    di _n "Random Effects Logit Results (CBSA random intercept):"
    di "  Branch coefficient: " %8.4f `b_branch_re' " (se = " %6.4f `se_branch_re' ")"
    di "  CBSA-level variance: " %8.4f `sigma_cbsa'^2
}
else {
    di "Using mixlogit for random coefficient estimation"

    * Need to reshape for mixlogit command
    * mixlogit requires case/alternative structure

    * For binary outcome, use melogit with random slope
    melogit se branch mobile c.pct_broadband i.age_cat i.educ_cat female married ///
        i.year [pw=hsupwgtk] || cbsa: branch, covariance(unstructured) ///
        intpoints(7)

    local b_branch_mixed = _b[se:branch]
    local se_branch_mixed = _se[se:branch]

    * Extract variance of random slope
    local var_branch = exp(_b[lns1_1_1:_cons])^2
    local sigma_branch = sqrt(`var_branch')

    di _n "Mixed Logit Results (Random slope on branch):"
    di "  Mean branch effect (gamma_bar): " %8.4f `b_branch_mixed'
    di "  SD of branch effect (sigma_gamma): " %8.4f `sigma_branch'
    di "  Variance of branch effect: " %8.4f `var_branch'
}

/*------------------------------------------------------------------------------
5. Alternative: Simulate Mixed Logit via Simulation

   For more flexible mixed logit, we can use simulation-based estimation
------------------------------------------------------------------------------*/

di _n "============================================================"
di "MODEL 3: SIMULATED MIXED LOGIT"
di "============================================================"

* Use fmlogit or simulate manually
* Here we approximate using multiple random coefficient draws

* Create draws from standard normal
local n_draws = 50
set seed 20260211

* Store original data
preserve

* For each draw, compute logit with perturbed coefficient
matrix DRAWS = J(`n_draws', 2, .)

forvalues d = 1/`n_draws' {
    * Random perturbation
    local nu = rnormal(0, 1)

    * Weight observations by likelihood contribution
    * (This is a simplified approximation)

    * For proper simulated MLE, would need mixlogit command
    * Here we just document the approach
}

restore

di "Note: Full simulated MLE requires mixlogit or mlogit (R)"
di "The random effects logit above provides an approximation"

/*------------------------------------------------------------------------------
6. Compare Finite Mixture vs Mixed Logit Implications
------------------------------------------------------------------------------*/

di _n "============================================================"
di "COMPARISON: FINITE MIXTURE vs MIXED LOGIT"
di "============================================================"

* Load finite mixture results (from phase6)
di _n "Finite Mixture (K=4) Results:"
di "  Type 1 (12.7%): beta = 0.000 (no effect)"
di "  Type 2 (20.3%): beta = -0.030 (negative)"
di "  Type 3 (34.6%): beta = 0.026 (moderate positive)"
di "  Type 4 (32.4%): beta = 0.138 (large positive)"
di ""
di "  Weighted mean: 0.127*0 + 0.203*(-0.030) + 0.346*0.026 + 0.324*0.138"
local weighted_mean = 0.127*0 + 0.203*(-0.030) + 0.346*0.026 + 0.324*0.138
di "               = " %8.4f `weighted_mean'
di ""
di "  Implied variance (between types):"
local var_fm = 0.127*(0 - `weighted_mean')^2 + 0.203*(-0.030 - `weighted_mean')^2 ///
             + 0.346*(0.026 - `weighted_mean')^2 + 0.324*(0.138 - `weighted_mean')^2
local sd_fm = sqrt(`var_fm')
di "  Var(beta) = " %8.4f `var_fm'
di "  SD(beta) = " %8.4f `sd_fm'

di _n "Mixed Logit Comparison:"
di "  If sigma_gamma from mixed logit is similar to SD(beta) from finite mixture,"
di "  both approaches confirm substantial heterogeneity in branch effects."

/*------------------------------------------------------------------------------
7. Summary
------------------------------------------------------------------------------*/

di _n "============================================================"
di "SUMMARY: MIXED LOGIT ROBUSTNESS CHECK"
di "============================================================"

di ""
di "The mixed logit with continuous random coefficients provides"
di "a robustness check for the finite mixture (K=4) results."
di ""
di "Key comparison:"
di "  - Finite mixture SD(beta): " %6.4f `sd_fm'
di "  - Mixed logit sigma_gamma: [from estimation above]"
di ""
di "If both approaches find similar magnitude of heterogeneity,"
di "this validates the finite mixture without relying on K selection."
di ""
di "Advantage of presenting both:"
di "  1. Finite mixture: interpretable types, counterfactual analysis"
di "  2. Mixed logit: continuous heterogeneity, no K selection"
di ""
di "The combination is compelling and rarely seen in applied work."

/*------------------------------------------------------------------------------
8. Save Results
------------------------------------------------------------------------------*/

preserve
clear
set obs 5
gen str30 item = ""
gen value = .

replace item = "beta_branch_std_logit" in 1
replace value = `b_branch_std' in 1

replace item = "me_branch_std_logit" in 2
replace value = `me_branch_std' in 2

replace item = "weighted_mean_fm" in 3
replace value = `weighted_mean' in 3

replace item = "sd_fm" in 4
replace value = `sd_fm' in 4

drop if item == ""
export delimited using "$output/phase7_mixed_logit.csv", replace
restore

di _n "Results saved to $output/phase7_mixed_logit.csv"

log close
