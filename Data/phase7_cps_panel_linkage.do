/*==============================================================================
Phase 7: CPS Rotation Panel Linkage

Exploiting CPS Rotation Structure for Limited Panel Data
Following Goodstein & Kutzbach (2022)

CPS households are interviewed for 4 months, rotate out for 8 months, then
return for 4 months. The FDIC supplement is administered in June.

Key insight: Some FDIC respondents appear in CPS basic monthly data in
neighboring months. We can match using household identifiers:
- HRHHID (Household ID)
- HRHHID2 (Household ID part 2)
- PULINENO (Person line number)

This gives LIMITED PANEL DATA for a subset of respondents:
- Observe employment status in months surrounding June FDIC survey
- Banking mode observed only in June (from FDIC supplement)

Applications:
(a) Estimate short-run employment transitions conditional on banking mode
(b) Test whether SE persistence differs by banking mode
(c) Construct Heckman-style initial conditions correction
==============================================================================*/

clear all
set more off
set matsize 11000
set seed 20260211

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"

capture log close
log using "$output/phase7_cps_panel_linkage.log", replace text

di _n "============================================================"
di "CPS ROTATION PANEL LINKAGE"
di "Limited Dynamics via Month-to-Month Matching"
di "============================================================"

/*------------------------------------------------------------------------------
1. Understanding CPS Rotation Structure

   Month-in-Sample (MIS) Pattern:
   MIS 1-4: First 4 months of interviews
   MIS 5-8: Second 4 months (after 8-month break)

   June FDIC supplement captures respondents in MIS 1-8.
   We can link to:
   - Previous months (May, April, March) for MIS 2-4 and 6-8
   - Subsequent months (July, August, September) for MIS 1-3 and 5-7

   The 4-8-4 rotation means ~75% of June respondents can be matched
   to at least one adjacent month.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "CPS ROTATION STRUCTURE"
di "============================================================"

di ""
di "Month-in-Sample (MIS) determines linkage possibilities:"
di ""
di "  MIS | Can link to previous | Can link to subsequent"
di "  ----|---------------------|----------------------"
di "   1  |        No           |     Yes (Jul-Sep)"
di "   2  |    Yes (May)        |     Yes (Jul-Sep)"
di "   3  |   Yes (Apr-May)     |     Yes (Jul-Aug)"
di "   4  |  Yes (Mar-May)      |     Yes (Jul only)"
di "   5  |        No           |     Yes (Jul-Sep)"
di "   6  |    Yes (May)        |     Yes (Jul-Sep)"
di "   7  |   Yes (Apr-May)     |     Yes (Jul-Aug)"
di "   8  |  Yes (Mar-May)      |        No"
di ""

/*------------------------------------------------------------------------------
2. Load FDIC/CPS June Data and Identify Rotation Group
------------------------------------------------------------------------------*/

di _n "============================================================"
di "LOADING FDIC/CPS DATA"
di "============================================================"

use "$datadir/analysis_dataset_with_se.dta", clear

keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013
keep if cbsa > 0 & cbsa != .
keep if banking_mode != .

* Create key variables
gen se = (self_employed == 1)
gen branch = (banking_mode == 3)
gen mobile = (banking_mode == 2)

* Check for rotation group variable (month-in-sample)
capture confirm variable hrmis
if _rc != 0 {
    di "Note: Month-in-sample (HRMIS) not in dataset."
    di "Simulating rotation structure for demonstration."

    * Simulate MIS distribution (uniform 1-8 in actual CPS)
    set seed 20260211
    gen hrmis = ceil(runiform() * 8)
}

tab hrmis, missing

* Linkage potential by MIS
gen can_link_prev = inlist(hrmis, 2, 3, 4, 6, 7, 8)
gen can_link_next = inlist(hrmis, 1, 2, 3, 5, 6, 7)
gen can_link_any = (can_link_prev == 1 | can_link_next == 1)

di _n "Linkage potential:"
tab can_link_any

quietly summarize can_link_any
local pct_linkable = r(mean) * 100
di "Approximately " %4.1f `pct_linkable' "% of June respondents can be linked"

/*------------------------------------------------------------------------------
3. Simulate Adjacent Month Data

   In practice, you would:
   1. Download CPS basic monthly files for May, July, etc.
   2. Match on HRHHID + HRHHID2 + PULINENO
   3. Extract employment status from matched months

   Here we simulate the matched panel for demonstration.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "SIMULATING ADJACENT MONTH MATCHES"
di "============================================================"

di "Note: In production, replace this with actual CPS monthly linkage."
di "Required files: CPS basic monthly for May, July, August, September."
di ""

* Simulate previous month employment (for those linkable)
* In reality: merge with May CPS using household IDs

* Employment transition probabilities (based on CPS literature)
* P(SE_t | SE_{t-1}) ≈ 0.92 (high persistence)
* P(SE_t | Wage_{t-1}) ≈ 0.02 (low entry)
* P(SE_t | NotWork_{t-1}) ≈ 0.05 (moderate entry)

set seed 20260211
gen se_prev = .
gen se_next = .

* For those who can link to previous month
* Simulate based on current state (with persistence)
replace se_prev = (runiform() < 0.92) if se == 1 & can_link_prev == 1
replace se_prev = (runiform() < 0.02) if se == 0 & wage_worker == 1 & can_link_prev == 1
replace se_prev = (runiform() < 0.05) if se == 0 & wage_worker == 0 & can_link_prev == 1

* For those who can link to next month
replace se_next = (runiform() < 0.92) if se == 1 & can_link_next == 1
replace se_next = (runiform() < 0.02) if se == 0 & wage_worker == 1 & can_link_next == 1
replace se_next = (runiform() < 0.05) if se == 0 & wage_worker == 0 & can_link_next == 1

* Transition indicators
gen entry_se = (se == 1 & se_prev == 0) if se_prev != .
gen exit_se = (se == 0 & se_prev == 1) if se_prev != .
gen persist_se = (se == 1 & se_prev == 1) if se_prev != .

di "Simulated transitions (previous month):"
tab se se_prev if se_prev != ., row

/*------------------------------------------------------------------------------
4. Estimate Employment Transitions by Banking Mode

   Key question: Does SE persistence differ by banking mode?
   If branch users have higher SE persistence, this supports the
   "relationship lending enables sustained entrepreneurship" hypothesis.
------------------------------------------------------------------------------*/

di _n "============================================================"
di "EMPLOYMENT TRANSITIONS BY BANKING MODE"
di "============================================================"

* Transition rates by banking mode
di _n "SE Persistence Rate by Banking Mode:"
di "(Among those SE in previous month, % still SE in June)"

forvalues b = 1/3 {
    quietly summarize persist_se if banking_mode == `b' & se_prev == 1 [aw=hsupwgtk]
    local persist_`b' = r(mean) * 100
    local n_`b' = r(N)
    local bname = cond(`b'==1, "Unbanked", cond(`b'==2, "Mobile", "Branch"))
    di "  `bname': " %5.1f `persist_`b'' "% (N = " `n_`b'' ")"
}

di _n "SE Entry Rate by Banking Mode:"
di "(Among those NOT SE in previous month, % SE in June)"

forvalues b = 1/3 {
    quietly summarize entry_se if banking_mode == `b' & se_prev == 0 [aw=hsupwgtk]
    local entry_`b' = r(mean) * 100
    local n_`b' = r(N)
    local bname = cond(`b'==1, "Unbanked", cond(`b'==2, "Mobile", "Branch"))
    di "  `bname': " %5.1f `entry_`b'' "% (N = " `n_`b'' ")"
}

/*------------------------------------------------------------------------------
5. Regression Analysis of Transitions

   Model: P(SE_t | SE_{t-1}, Banking Mode, X)

   This tests whether banking mode affects:
   (a) SE persistence (among those previously SE)
   (b) SE entry (among those not previously SE)
------------------------------------------------------------------------------*/

di _n "============================================================"
di "REGRESSION ANALYSIS OF TRANSITIONS"
di "============================================================"

* Sample with valid transitions
keep if se_prev != .

di "Sample for transition analysis: " _N " observations"

* Persistence regression (among previously SE)
di _n "Model 1: SE Persistence (among SE_{t-1} = 1)"
logit se branch mobile pct_broadband age i.educ i.year [pw=hsupwgtk] ///
    if se_prev == 1, vce(cluster cbsa)

local b_branch_persist = _b[branch]
local se_branch_persist = _se[branch]
local or_branch_persist = exp(`b_branch_persist')

di "  Branch coefficient: " %7.4f `b_branch_persist' " (se = " %6.4f `se_branch_persist' ")"
di "  Odds ratio: " %6.3f `or_branch_persist'

margins, dydx(branch) atmeans
local me_branch_persist = r(table)[1,1]
di "  Marginal effect on persistence: " %7.4f `me_branch_persist'

* Entry regression (among not previously SE)
di _n "Model 2: SE Entry (among SE_{t-1} = 0)"
logit se branch mobile pct_broadband age i.educ i.year [pw=hsupwgtk] ///
    if se_prev == 0, vce(cluster cbsa)

local b_branch_entry = _b[branch]
local se_branch_entry = _se[branch]
local or_branch_entry = exp(`b_branch_entry')

di "  Branch coefficient: " %7.4f `b_branch_entry' " (se = " %6.4f `se_branch_entry' ")"
di "  Odds ratio: " %6.3f `or_branch_entry'

margins, dydx(branch) atmeans
local me_branch_entry = r(table)[1,1]
di "  Marginal effect on entry: " %7.4f `me_branch_entry'

/*------------------------------------------------------------------------------
6. Decomposition: State Dependence vs Heterogeneity

   Following Heckman (1981): Observed persistence = True state dependence +
                                                     Spurious (heterogeneity)

   With lagged employment, we can partially separate these:
   - If branch effect persists after controlling for SE_{t-1},
     it reflects ongoing advantage (not just selection)
   - If branch effect shrinks with SE_{t-1} control,
     cross-sectional effect was partly spurious
------------------------------------------------------------------------------*/

di _n "============================================================"
di "STATE DEPENDENCE vs HETEROGENEITY"
di "============================================================"

* Model without lagged SE
di _n "Model A: Cross-sectional (no lagged SE)"
quietly logit se branch mobile pct_broadband age i.educ i.year [pw=hsupwgtk], ///
    vce(cluster cbsa)
local b_branch_xsec = _b[branch]
local se_branch_xsec = _se[branch]

* Model with lagged SE
di "Model B: With lagged SE control"
logit se branch mobile se_prev pct_broadband age i.educ i.year [pw=hsupwgtk], ///
    vce(cluster cbsa)
local b_branch_panel = _b[branch]
local se_branch_panel = _se[branch]
local b_se_prev = _b[se_prev]
local se_se_prev = _se[se_prev]

di _n "Comparison:"
di "  Branch effect (cross-sectional): " %7.4f `b_branch_xsec' " (se = " %6.4f `se_branch_xsec' ")"
di "  Branch effect (with SE_{t-1}):   " %7.4f `b_branch_panel' " (se = " %6.4f `se_branch_panel' ")"
di "  Lagged SE effect:                " %7.4f `b_se_prev' " (se = " %6.4f `se_se_prev' ")"
di ""

local pct_change = (`b_branch_panel' - `b_branch_xsec') / `b_branch_xsec' * 100
di "Change in branch effect with lagged SE control: " %5.1f `pct_change' "%"

if `pct_change' < -20 {
    di "Substantial reduction suggests cross-sectional effect partly reflects"
    di "unobserved heterogeneity rather than pure state dependence."
}
else {
    di "Modest change suggests cross-sectional effect largely reflects"
    di "genuine banking mode advantage, not just selection."
}

/*------------------------------------------------------------------------------
7. Initial Conditions Correction (Heckman 1981)

   The initial conditions problem: SE_{t-1} is endogenous if unobserved
   heterogeneity affects both banking choice and SE.

   Heckman's solution: Model the initial condition as a function of
   pre-sample information (here: demographics, CBSA characteristics).
------------------------------------------------------------------------------*/

di _n "============================================================"
di "INITIAL CONDITIONS CORRECTION"
di "============================================================"

* First stage: model initial SE status
di "First stage: P(SE_{t-1} | pre-sample X)"
probit se_prev pct_broadband age i.educ i.year female married [pw=hsupwgtk], ///
    vce(cluster cbsa)
predict se_prev_hat, pr

* Second stage: include predicted initial condition
di _n "Second stage: P(SE_t | Banking, SE_{t-1}, X, SE_prev_hat)"
logit se branch mobile se_prev se_prev_hat pct_broadband age i.educ i.year ///
    [pw=hsupwgtk], vce(cluster cbsa)

local b_branch_ic = _b[branch]
local se_branch_ic = _se[branch]
local b_se_prev_ic = _b[se_prev]
local b_se_prev_hat = _b[se_prev_hat]

di _n "Initial Conditions Corrected Results:"
di "  Branch effect: " %7.4f `b_branch_ic' " (se = " %6.4f `se_branch_ic' ")"
di "  True state dependence (SE_{t-1}): " %7.4f `b_se_prev_ic'
di "  Correlated heterogeneity (SE_prev_hat): " %7.4f `b_se_prev_hat'

/*------------------------------------------------------------------------------
8. Summary: What Panel Linkage Reveals
------------------------------------------------------------------------------*/

di _n "============================================================"
di "SUMMARY: CPS PANEL LINKAGE FINDINGS"
di "============================================================"

di ""
di "The CPS rotation structure provides limited panel data that reveals:"
di ""
di "1. SE PERSISTENCE by banking mode:"
di "   - Branch users: " %5.1f `persist_1' "% persistence rate"
di "   - Mobile users: " %5.1f `persist_2' "% persistence rate"
di "   - Unbanked:     " %5.1f `persist_3' "% persistence rate"
di ""
di "2. SE ENTRY by banking mode:"
di "   - Branch users: " %5.1f `entry_1' "% entry rate"
di "   - Mobile users: " %5.1f `entry_2' "% entry rate"
di "   - Unbanked:     " %5.1f `entry_3' "% entry rate"
di ""
di "3. STATE DEPENDENCE vs HETEROGENEITY:"
di "   - Cross-sectional branch effect: " %7.4f `b_branch_xsec'
di "   - Panel branch effect (with SE_{t-1}): " %7.4f `b_branch_panel'
di "   - Change: " %5.1f `pct_change' "%"
di ""
di "4. IMPLICATIONS FOR STATIC MODEL:"
di "   The panel evidence [supports/qualifies] the static interpretation."
di "   [If branch effect persists with lagged SE control, static effect"
di "    captures genuine ongoing advantage, not just selection.]"

/*------------------------------------------------------------------------------
9. Save Results
------------------------------------------------------------------------------*/

preserve
clear
set obs 10
gen str40 parameter = ""
gen value = .

replace parameter = "persist_branch" in 1
replace value = `persist_3' in 1

replace parameter = "persist_mobile" in 2
replace value = `persist_2' in 2

replace parameter = "entry_branch" in 3
replace value = `entry_3' in 3

replace parameter = "entry_mobile" in 4
replace value = `entry_2' in 4

replace parameter = "b_branch_xsec" in 5
replace value = `b_branch_xsec' in 5

replace parameter = "b_branch_panel" in 6
replace value = `b_branch_panel' in 6

replace parameter = "b_se_prev" in 7
replace value = `b_se_prev' in 7

replace parameter = "me_branch_persist" in 8
replace value = `me_branch_persist' in 8

replace parameter = "me_branch_entry" in 9
replace value = `me_branch_entry' in 9

drop if parameter == ""
export delimited using "$output/phase7_cps_panel.csv", replace
restore

di _n "Results saved to $output/phase7_cps_panel.csv"

di _n "============================================================"
di "NOTE ON IMPLEMENTATION"
di "============================================================"
di ""
di "This code demonstrates the methodology with simulated transitions."
di "For actual implementation:"
di ""
di "1. Download CPS basic monthly files (May, July, Aug, Sep)"
di "2. Match on: HRHHID + HRHHID2 + PULINENO + (year, state)"
di "3. Extract PEMLR (employment status) from matched months"
di "4. Recode to SE indicator using class-of-worker variable"
di ""
di "See Goodstein & Kutzbach (2022) for detailed linking protocol."

log close
