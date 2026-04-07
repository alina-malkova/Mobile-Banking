/*==============================================================================
Phase 9: Type-Specific Average Marginal Effects on the Probability Scale

Comment 9 — Converting log-odds coefficients to probability-scale AMEs.

Key idea:
  The logit model yields type-specific coefficients gamma_k on the
  log-odds scale. To interpret them as effects on P(SE), we compute
  average marginal effects (AMEs) within each type using the margins
  command after binary logit.

  AME_k = (1/N_k) * sum_i [ dP(SE_i)/d(branch_i) ]   for i in type k

  The weighted average AME = sum_k pi_k * AME_k is then compared to the
  reduced-form OLS coefficient (0.003).

Connection to Proposition 1 (Jensen's inequality):
  Because the logit link is nonlinear, the weighted average of
  type-specific AMEs != AME evaluated at the weighted-average coefficients.
  This is precisely the nonlinearity that makes the number of types K
  matter for counterfactual magnitudes.

==============================================================================*/

clear all
set more off
set matsize 11000

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/_Research/Mobile_Money_Banking/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase9_marginal_effects.log", replace text

di _n "============================================================"
di "PHASE 9: TYPE-SPECIFIC AVERAGE MARGINAL EFFECTS"
di "Converting log-odds to probability scale (Comment 9)"
di "============================================================"
di "Date: $S_DATE  Time: $S_TIME"

/*------------------------------------------------------------------------------
1. Load Data and Apply Sample Restrictions
   (Identical to phase3_counterfactual_final.do)
------------------------------------------------------------------------------*/

use "$datadir/analysis_dataset_with_se.dta", clear

keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .
keep if banking_mode != .

di _n "Sample size after restrictions: " _N

* Create employment status
capture drop emp_status
gen emp_status = .
replace emp_status = 1 if wage_worker == 1
replace emp_status = 2 if self_employed == 1
replace emp_status = 3 if employed != 1

* Banking mode
capture drop bank_mode
gen bank_mode = banking_mode

* Binary outcomes and treatments
gen se = (self_employed == 1)
gen branch = (bank_mode == 3)
gen mobile = (bank_mode == 2)

* Age categories (same as phase3)
gen age_cat = 1 if age >= 18 & age < 35
replace age_cat = 2 if age >= 35 & age < 50
replace age_cat = 3 if age >= 50 & age <= 64

* Education categories (same as phase3)
gen educ_cat = .
replace educ_cat = 1 if no_hs == 1
replace educ_cat = 2 if hs_diploma == 1
replace educ_cat = 3 if some_college == 1
replace educ_cat = 4 if college_degree == 1

* Lagged SE rate by CBSA (dynamic term)
bysort cbsa year: egen cbsa_se_rate = mean(se)
bysort cbsa (year): gen lag_se_rate = cbsa_se_rate[_n-1]
replace lag_se_rate = 0.11 if lag_se_rate == .
gen dynamic_se = se * lag_se_rate

/*------------------------------------------------------------------------------
2. Create Type Assignments (Observable-Based Scoring)
   (Identical logic to phase2_type_selection.do / phase3_counterfactual_final.do)
------------------------------------------------------------------------------*/

di _n "============================================================"
di "TYPE ASSIGNMENT (Observable-Based Score, K=4)"
di "============================================================"

gen type_score = 0

* Age: older = more entrepreneurial experience
replace type_score = type_score + 2 if age_cat == 3
replace type_score = type_score + 1 if age_cat == 2

* Education: college = human capital for entrepreneurship
replace type_score = type_score + 2 if educ_cat == 4
replace type_score = type_score + 1 if educ_cat == 3

* Banking mode: branch = relationship borrower
replace type_score = type_score + 2 if bank_mode == 3
replace type_score = type_score - 1 if bank_mode == 1

* Current SE status
replace type_score = type_score + 4 if emp_status == 2

* Tech adoption proxy: mobile banking
replace type_score = type_score - 1 if bank_mode == 2

* Rank and assign to quartiles (K=4)
egen type_rank = rank(type_score), unique
quietly summarize type_rank
local max_rank = r(max)

gen type4 = 1 if type_rank <= `max_rank'/4
replace type4 = 2 if type_rank > `max_rank'/4 & type_rank <= `max_rank'/2
replace type4 = 3 if type_rank > `max_rank'/2 & type_rank <= 3*`max_rank'/4
replace type4 = 4 if type4 == .

* Type shares and baseline SE rates
di _n "Type Descriptives:"
di "==================="
forvalues k = 1/4 {
    quietly summarize type4 if type4 == `k' [aw=hsupwgtk]
    local n_k = r(N)
    quietly count if type4 == `k'
    local count_k = r(N)
    quietly summarize se if type4 == `k' [aw=hsupwgtk]
    local se_rate_k = r(mean)
    quietly summarize branch if type4 == `k' [aw=hsupwgtk]
    local branch_rate_k = r(mean)
    di "Type `k': N=" %6.0f `count_k' ///
       ", SE rate=" %5.2f (`se_rate_k'*100) "%" ///
       ", Branch share=" %5.2f (`branch_rate_k'*100) "%"
}

* Compute weighted type shares
forvalues k = 1/4 {
    gen tau_`k' = (type4 == `k')
}
quietly summarize tau_1 [aw=hsupwgtk]
local share_1 = r(mean)
quietly summarize tau_2 [aw=hsupwgtk]
local share_2 = r(mean)
quietly summarize tau_3 [aw=hsupwgtk]
local share_3 = r(mean)
quietly summarize tau_4 [aw=hsupwgtk]
local share_4 = r(mean)

di _n "Weighted type shares:"
forvalues k = 1/4 {
    di "  pi_`k' = " %6.4f `share_`k'' " (" %5.1f (`share_`k''*100) "%)"
}

/*------------------------------------------------------------------------------
3. Reduced-Form OLS Benchmark
   Estimate the pooled (homogeneous) OLS coefficient for comparison.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "REDUCED-FORM OLS BENCHMARK"
di "============================================================"

reg se branch i.age_cat i.educ_cat ///
    female married has_children ///
    mobile dynamic_se i.year [aw=hsupwgtk], vce(cluster cbsa)

local ols_coef = _b[branch]
local ols_se   = _se[branch]
local ols_t    = `ols_coef' / `ols_se'

di _n "OLS coefficient on branch: " %8.5f `ols_coef' ///
     " (SE=" %8.5f `ols_se' ", t=" %5.2f `ols_t' ")"
di "This is the LPM effect: dP(SE)/d(branch), pooled across all types."

/*------------------------------------------------------------------------------
4. Type-Specific Binary Logit Estimation and AMEs
   For each type k, estimate logit of SE on branch + controls,
   then compute AME of branch on P(SE) via margins.

   NOTE: We use binary logit (SE vs not-SE) instead of 3-category mlogit
   because mlogit suffers from perfect prediction/separation in small type
   subsamples, yielding nonsensical coefficients. Binary logit directly
   gives the probability-scale AME of branch on P(SE) that we need.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "TYPE-SPECIFIC BINARY LOGIT & MARGINAL EFFECTS"
di "============================================================"

* Store results
matrix AME = J(4, 5, .)
matrix colnames AME = share gamma_logodds ame_prob_scale ame_se ame_z

forvalues k = 1/4 {

    di _n "------------------------------------------------------------"
    di "TYPE `k'"
    di "------------------------------------------------------------"

    * Estimate binary logit for type k: P(SE=1 | branch, controls)
    * Use simpler spec to avoid perfect prediction
    logit se branch i.age_cat i.educ_cat ///
        female married ///
        mobile i.year ///
        if type4 == `k' [pw=hsupwgtk], vce(cluster cbsa) iter(50)

    local gamma_`k' = _b[branch]
    local gamma_se_`k' = _se[branch]

    di "  gamma_`k' (log-odds) = " %8.5f `gamma_`k'' ///
       " (SE=" %8.5f `gamma_se_`k'' ")"

    * Compute average marginal effect of branch on P(SE)
    margins, dydx(branch) post

    local ame_`k'    = _b[branch]
    local ame_se_`k' = _se[branch]
    local ame_z_`k'  = `ame_`k'' / `ame_se_`k''

    di "  AME_`k' (prob scale) = " %8.5f `ame_`k'' ///
       " (SE=" %8.5f `ame_se_`k'' ", z=" %5.2f `ame_z_`k'' ")"
    if `gamma_`k'' != 0 {
        di "  Ratio AME/gamma = " %6.4f (`ame_`k'' / `gamma_`k'') ///
           "  (logit compression factor)"
    }

    * Store in matrix
    matrix AME[`k', 1] = `share_`k''
    matrix AME[`k', 2] = `gamma_`k''
    matrix AME[`k', 3] = `ame_`k''
    matrix AME[`k', 4] = `ame_se_`k''
    matrix AME[`k', 5] = `ame_z_`k''
}

di _n "============================================================"
di "SUMMARY: TYPE-SPECIFIC RESULTS"
di "============================================================"
di _n "  Type  |  Share  | gamma(log-odds) |   AME(prob)   |  AME SE   |    z"
di "  ------+--------+-----------------+---------------+-----------+------"
forvalues k = 1/4 {
    di "    `k'   | " %5.3f `share_`k'' " |   " %10.5f `gamma_`k'' ///
       "  | " %12.5f `ame_`k'' " | " %9.5f `ame_se_`k'' ///
       " | " %5.2f `ame_z_`k''
}

/*------------------------------------------------------------------------------
5. Weighted Average AME and Comparison to OLS
------------------------------------------------------------------------------*/

di _n "============================================================"
di "WEIGHTED AVERAGE AME vs. OLS"
di "============================================================"

* Weighted average AME: sum_k pi_k * AME_k
local wtd_ame = 0
forvalues k = 1/4 {
    local wtd_ame = `wtd_ame' + `share_`k'' * `ame_`k''
}

* Standard error of weighted average via delta method
* Var(wtd_ame) approx= sum_k pi_k^2 * Var(AME_k)
* (treating shares as fixed)
local var_wtd_ame = 0
forvalues k = 1/4 {
    local var_wtd_ame = `var_wtd_ame' + (`share_`k'')^2 * (`ame_se_`k'')^2
}
local se_wtd_ame = sqrt(`var_wtd_ame')

* Ratio
local ratio = `wtd_ame' / `ols_coef'

di "Weighted average AME (prob scale): " %10.5f `wtd_ame' ///
   " (SE=" %10.5f `se_wtd_ame' ")"
di "Reduced-form OLS coefficient:      " %10.5f `ols_coef' ///
   " (SE=" %10.5f `ols_se' ")"
di "Ratio (weighted AME / OLS):         " %8.4f `ratio'
di ""

* Decompose the contributions
di "Contribution of each type to weighted AME:"
di "-------------------------------------------"
forvalues k = 1/4 {
    local contrib_k = `share_`k'' * `ame_`k''
    local pct_k = `contrib_k' / `wtd_ame' * 100
    di "  Type `k': pi*AME = " %10.5f `contrib_k' ///
       " (" %5.1f `pct_k' "% of total)"
}

/*------------------------------------------------------------------------------
6. Jensen's Inequality Demonstration (Proposition 1)
   The nonlinearity of the logit link means:
     E[g(gamma_k)] != g(E[gamma_k])
   where g() maps log-odds to probability-scale AME.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "JENSEN'S INEQUALITY / PROPOSITION 1"
di "============================================================"

* Weighted average of log-odds coefficients
local wtd_gamma = 0
forvalues k = 1/4 {
    local wtd_gamma = `wtd_gamma' + `share_`k'' * `gamma_`k''
}
di "Weighted average gamma (log-odds):  " %10.5f `wtd_gamma'
di "Weighted average AME (prob scale):  " %10.5f `wtd_ame'

* If the logit were linear, AME at avg gamma would equal avg of AMEs.
* The difference quantifies the nonlinearity.
* Approximate AME at weighted-average gamma using the pooled sample:
* AME_homog approx = p_bar * (1 - p_bar) * gamma_bar
* where p_bar is the average P(SE) in the sample.
quietly summarize se [aw=hsupwgtk]
local p_bar = r(mean)
local ame_at_avg_gamma = `p_bar' * (1 - `p_bar') * `wtd_gamma'

di ""
di "Approximate AME at weighted-avg gamma: " %10.5f `ame_at_avg_gamma'
di "  (using p_bar*(1-p_bar)*gamma_bar, p_bar=" %6.4f `p_bar' ")"
di ""
di "Weighted avg of type-specific AMEs:    " %10.5f `wtd_ame'
di "AME at weighted-avg coefficients:      " %10.5f `ame_at_avg_gamma'
di "Difference (Jensen's gap):             " %10.5f (`wtd_ame' - `ame_at_avg_gamma')
di ""
di "Interpretation:"
di "  Because the logit link is nonlinear (concave/convex depending on"
di "  the region), the weighted average of type-specific AMEs differs"
di "  from the AME evaluated at the weighted-average coefficients."
di "  This is exactly the nonlinearity that makes K matter:"
di "  misspecifying the number of types changes the curvature-weighted"
di "  average and therefore the counterfactual magnitude."

/*------------------------------------------------------------------------------
7. Output: phase9_marginal_effects.csv
   Columns: type, share, gamma_logodds, ame_prob_scale, ame_se
   Plus summary row for weighted average.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "SAVING RESULTS"
di "============================================================"

preserve
clear
set obs 6

gen str25 type = ""
gen share = .
gen gamma_logodds = .
gen ame_prob_scale = .
gen ame_se = .

* Type-specific rows
replace type = "type_1" in 1
replace share = `share_1' in 1
replace gamma_logodds = `gamma_1' in 1
replace ame_prob_scale = `ame_1' in 1
replace ame_se = `ame_se_1' in 1

replace type = "type_2" in 2
replace share = `share_2' in 2
replace gamma_logodds = `gamma_2' in 2
replace ame_prob_scale = `ame_2' in 2
replace ame_se = `ame_se_2' in 2

replace type = "type_3" in 3
replace share = `share_3' in 3
replace gamma_logodds = `gamma_3' in 3
replace ame_prob_scale = `ame_3' in 3
replace ame_se = `ame_se_3' in 3

replace type = "type_4" in 4
replace share = `share_4' in 4
replace gamma_logodds = `gamma_4' in 4
replace ame_prob_scale = `ame_4' in 4
replace ame_se = `ame_se_4' in 4

* Summary row: weighted average AME
replace type = "weighted_avg_ame" in 5
replace share = 1 in 5
replace gamma_logodds = `wtd_gamma' in 5
replace ame_prob_scale = `wtd_ame' in 5
replace ame_se = `se_wtd_ame' in 5

* Summary row: OLS coefficient for comparison
replace type = "ols_coefficient" in 6
replace share = 1 in 6
replace gamma_logodds = . in 6
replace ame_prob_scale = `ols_coef' in 6
replace ame_se = `ols_se' in 6

export delimited using "$output/phase9_marginal_effects.csv", replace
restore

di "Results saved to: $output/phase9_marginal_effects.csv"

/*------------------------------------------------------------------------------
8. Summary Table for Paper
------------------------------------------------------------------------------*/

di _n "============================================================"
di "SUMMARY FOR PAPER (Comment 9)"
di "============================================================"
di ""
di "Table: Type-Specific Average Marginal Effects of Branch Access"
di "       on P(Self-Employment), K=4 Model"
di ""
di "---------------------------------------------------------------"
di "  Type  |  pi_k  | gamma_k     | AME_k        | AME SE       "
di "        | (share)| (log-odds)  | (prob scale) |              "
di "-------+--------+-------------+--------------+--------------"
forvalues k = 1/4 {
    di "    `k'   | " %5.3f `share_`k'' " | " %11.5f `gamma_`k'' ///
       " | " %12.5f `ame_`k'' " | " %12.5f `ame_se_`k''
}
di "-------+--------+-------------+--------------+--------------"
di "  Wtd  |  1.000 | " %11.5f `wtd_gamma' ///
   " | " %12.5f `wtd_ame' " | " %12.5f `se_wtd_ame'
di "  OLS  |        |             | " %12.5f `ols_coef' " | " %12.5f `ols_se'
di "---------------------------------------------------------------"
di ""
di "  Weighted AME / OLS ratio: " %8.4f `ratio'
di ""
di "Note: gamma_k is the binary logit coefficient on branch for SE"
di "  (log-odds scale). AME_k is the average marginal effect"
di "  dP(SE)/d(branch) averaged over individuals in type k"
di "  (probability scale). Weighted average uses type shares pi_k."
di "  Standard errors clustered at CBSA level."
di ""
di "Jensen's inequality gap: " %10.5f (`wtd_ame' - `ame_at_avg_gamma')
di "  (weighted avg of AMEs) - (AME at weighted avg coefficients)"
di "  This gap = 0 iff the logit link is linear (i.e., never)."

di _n "============================================================"
di "PHASE 9 COMPLETE"
di "============================================================"

log close
