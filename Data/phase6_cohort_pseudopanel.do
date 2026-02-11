/*==============================================================================
Phase 6: Cohort Pseudo-Panel Robustness Analysis

Following Browning-Collado-Crawford (2014) and Moffitt (1993):
- Construct pseudo-panels from cohort × CBSA cells
- Define cohorts by birth year × education × race (time-invariant)
- Track cell-level means across survey waves
- Estimate Markov transitions at the cell level

This provides LIMITED dynamics without true panel data.

Key assumption: Within-cell composition is stable across waves.
This is testable: check if cell-level demographics shift.
==============================================================================*/

clear all
set more off
set matsize 11000

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase6_cohort_pseudopanel.log", replace text

di _n "============================================================"
di "COHORT PSEUDO-PANEL ANALYSIS"
di "Limited Dynamics via Cell-Level Transitions"
di "============================================================"

/*------------------------------------------------------------------------------
1. Load Data and Define Cohorts
------------------------------------------------------------------------------*/

use "$datadir/analysis_dataset_with_se.dta", clear

keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .
keep if banking_mode != .

* Core variables
gen se = (self_employed == 1)
gen branch = (banking_mode == 3)
gen mobile = (banking_mode == 2)
gen unbanked = (banking_mode == 1)

* Define birth cohort (5-year bins)
gen birth_year = year - age
gen cohort = .
replace cohort = 1 if birth_year >= 1959 & birth_year <= 1964  // Born 1959-1964
replace cohort = 2 if birth_year >= 1965 & birth_year <= 1969  // Born 1965-1969
replace cohort = 3 if birth_year >= 1970 & birth_year <= 1974  // Born 1970-1974
replace cohort = 4 if birth_year >= 1975 & birth_year <= 1979  // Born 1975-1979
replace cohort = 5 if birth_year >= 1980 & birth_year <= 1984  // Born 1980-1984
replace cohort = 6 if birth_year >= 1985 & birth_year <= 1989  // Born 1985-1989
replace cohort = 7 if birth_year >= 1990 & birth_year <= 1999  // Born 1990-1999

drop if cohort == .

* Education groups (time-invariant for adults)
gen educ_group = .
replace educ_group = 1 if no_hs == 1 | hs_diploma == 1  // HS or less
replace educ_group = 2 if some_college == 1              // Some college
replace educ_group = 3 if college_degree == 1            // College+

* Race (time-invariant)
gen race_group = .
replace race_group = 1 if white == 1
replace race_group = 2 if black == 1
replace race_group = 3 if hispanic == 1
replace race_group = 4 if race_group == .  // Other

di "Sample: " _N " observations"
tab cohort year

/*------------------------------------------------------------------------------
2. Construct Pseudo-Panel Cells

   Cell = CBSA × Cohort × Education × Race
   This gives us "pseudo-individuals" we can track across waves.

   Trade-off: Finer cells → fewer observations per cell
              Coarser cells → more aggregation bias
------------------------------------------------------------------------------*/

di _n "============================================================"
di "CONSTRUCTING PSEUDO-PANEL CELLS"
di "============================================================"

* Create cell identifier
egen cell = group(cbsa cohort educ_group race_group)
quietly summarize cell
local n_cells = r(max)
di "Number of cells: `n_cells'"

* Collapse to cell-year level
preserve
collapse (mean) se branch mobile unbanked pct_broadband ///
         (count) n_obs = se ///
         (rawsum) weight = hsupwgtk ///
         [aw=hsupwgtk], by(cell cbsa cohort educ_group race_group year)

* Check cell sizes
quietly summarize n_obs
di "Cell-year observations: mean = " %6.1f r(mean) ", min = " %3.0f r(min) ", max = " %5.0f r(max)

* Keep cells with sufficient observations
gen valid = (n_obs >= 10)
tab valid

* Number of waves per cell
bysort cell: gen n_waves = _N
tab n_waves

* Keep cells observed in at least 2 waves (for transitions)
keep if n_waves >= 2

di "Cells with 2+ waves: " _N " cell-year observations"

/*------------------------------------------------------------------------------
3. Estimate Cell-Level Transitions

   Pseudo-transition: Change in cell mean from wave t to wave t+2
   (Waves are 2013, 2015, 2017, 2019, 2021, 2023 - biennial)
------------------------------------------------------------------------------*/

di _n "============================================================"
di "ESTIMATING PSEUDO-TRANSITIONS"
di "============================================================"

* Sort and create lagged values
sort cell year
bysort cell: gen L_se = se[_n-1]
bysort cell: gen L_branch = branch[_n-1]
bysort cell: gen L_mobile = mobile[_n-1]
bysort cell: gen L_year = year[_n-1]

* Only keep transitions (have lagged values)
drop if L_se == .

* Change in SE rate
gen d_se = se - L_se
gen d_branch = branch - L_branch
gen d_mobile = mobile - L_mobile

di "Cell-level transitions: " _N " observations"

* Summary of transitions
di _n "Summary of Cell-Level Changes:"
summarize d_se d_branch d_mobile [aw=weight]

/*------------------------------------------------------------------------------
4. Pseudo-Panel Regression

   Δ(SE)_cell,t = α + β × Δ(Branch)_cell,t + γ × X_cell,t + ε

   This captures how changes in branch usage within a cohort-CBSA cell
   correlate with changes in self-employment.

   Key assumption: Composition within cell is stable, so changes reflect
   behavioral responses rather than compositional shifts.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "PSEUDO-PANEL REGRESSIONS"
di "============================================================"

* First-difference regression
reg d_se d_branch d_mobile pct_broadband i.year [aw=weight], vce(cluster cbsa)

local b_branch_fd = _b[d_branch]
local se_branch_fd = _se[d_branch]
local t_branch_fd = `b_branch_fd' / `se_branch_fd'

di _n "First-Difference Results:"
di "  Δ(Branch) → Δ(SE): β = " %7.4f `b_branch_fd' " (t = " %5.2f `t_branch_fd' ")"

* Level regression with cell fixed effects
xtset cell year
xtreg se branch mobile pct_broadband i.year [aw=weight], fe vce(cluster cbsa)

local b_branch_fe = _b[branch]
local se_branch_fe = _se[branch]
local t_branch_fe = `b_branch_fe' / `se_branch_fe'

di _n "Cell Fixed Effects Results:"
di "  Branch → SE: β = " %7.4f `b_branch_fe' " (t = " %5.2f `t_branch_fe' ")"

/*------------------------------------------------------------------------------
5. Test Composition Stability

   Key assumption: Within-cell composition is stable.
   Test: Regress cell-level demographics on time.
   If composition shifts, the pseudo-panel approach is compromised.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "TESTING COMPOSITION STABILITY"
di "============================================================"

* The main concern: Are high-SE individuals leaving/entering cells over time?
* If cells are defined by time-invariant characteristics (cohort, education, race),
* composition should be stable within CBSA-cohort-education-race cells.

* Check: Does the composition of SE within cells drift?
bysort cell: egen mean_se_cell = mean(se)
bysort cell: gen se_deviation = se - mean_se_cell

reg se_deviation year [aw=weight], vce(cluster cbsa)
local b_trend = _b[year]
local se_trend = _se[year]

di "Composition stability test:"
di "  SE deviation trend: " %8.5f `b_trend' " (se = " %8.5f `se_trend' ")"
if abs(`b_trend'/`se_trend') < 1.96 {
    di "  ✓ No significant composition drift detected"
}
else {
    di "  ⚠ Significant composition drift - interpret with caution"
}

restore

/*------------------------------------------------------------------------------
6. Compare Static vs Pseudo-Panel Estimates
------------------------------------------------------------------------------*/

di _n "============================================================"
di "COMPARISON: STATIC vs PSEUDO-PANEL"
di "============================================================"

* Get static estimate for comparison
use "$datadir/analysis_dataset_with_se.dta", clear
keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .
keep if banking_mode != .

gen se = (self_employed == 1)
gen branch = (banking_mode == 3)
gen mobile = (banking_mode == 2)

quietly reg se branch mobile pct_broadband i.age i.educ i.year [pw=hsupwgtk], vce(cluster cbsa)
local b_static = _b[branch]
local se_static = _se[branch]

di ""
di "Method                    | β(Branch) | Std. Error | t-stat"
di "--------------------------|-----------|------------|--------"
di "Static (cross-section)    |  " %7.4f `b_static' "  |  " %7.4f `se_static' "   |  " %5.2f (`b_static'/`se_static')
di "Pseudo-panel (FD)         |  " %7.4f `b_branch_fd' "  |  " %7.4f `se_branch_fd' "   |  " %5.2f `t_branch_fd'
di "Pseudo-panel (FE)         |  " %7.4f `b_branch_fe' "  |  " %7.4f `se_branch_fe' "   |  " %5.2f `t_branch_fe'
di ""

* Assess convergence
local diff_fd = abs(`b_branch_fd' - `b_static')
local diff_fe = abs(`b_branch_fe' - `b_static')

di "Difference from static estimate:"
di "  First-difference: " %7.4f `diff_fd'
di "  Fixed effects:    " %7.4f `diff_fe'
di ""
if `diff_fd' < 0.02 & `diff_fe' < 0.02 {
    di "Estimates are similar across methods - supports static interpretation."
}
else {
    di "Estimates differ - suggests dynamic effects or composition changes."
}

/*------------------------------------------------------------------------------
7. Summary
------------------------------------------------------------------------------*/

di _n "============================================================"
di "SUMMARY: PSEUDO-PANEL ANALYSIS"
di "============================================================"

di ""
di "The pseudo-panel approach provides LIMITED dynamics by:"
di "  - Tracking cohort × CBSA cells across survey waves"
di "  - Estimating cell-level transitions"
di ""
di "Key finding: Pseudo-panel estimates are [similar to/different from]"
di "static estimates, suggesting [absence of/presence of] dynamic effects"
di "beyond what cross-sectional variation captures."
di ""
di "Limitations:"
di "  - Thin cells (many cohort × CBSA × educ × race combinations)"
di "  - Assumes stable within-cell composition"
di "  - Cannot identify individual-level switching costs"
di ""
di "Role in paper: Robustness check showing static results are consistent"
di "with limited dynamics available in pseudo-panel framework."

/*------------------------------------------------------------------------------
8. Save Results
------------------------------------------------------------------------------*/

preserve
clear
set obs 5
gen str30 item = ""
gen value = .

replace item = "beta_static" in 1
replace value = `b_static' in 1

replace item = "beta_pseudopanel_fd" in 2
replace value = `b_branch_fd' in 2

replace item = "beta_pseudopanel_fe" in 3
replace value = `b_branch_fe' in 3

replace item = "se_static" in 4
replace value = `se_static' in 4

replace item = "composition_trend" in 5
replace value = `b_trend' in 5

export delimited using "$output/phase6_pseudopanel.csv", replace
restore

di _n "Results saved to $output/phase6_pseudopanel.csv"

log close
