/*******************************************************************************
* Mobile Banking, Bank Branch Closures, and Self-Employment
* Phase 0: Reduced-Form Analysis
*
* This do-file implements:
*   1. Descriptive statistics and trends
*   2. OLS baseline regressions
*   3. CBSA fixed effects specifications
*   4. IV regressions using broadband as instrument
*   5. Heterogeneity analysis
*******************************************************************************/

clear all
set more off
set matsize 10000

* Install required packages if not already installed
capture which reghdfe
if _rc ssc install reghdfe, replace
capture which ftools
if _rc ssc install ftools, replace
capture which estout
if _rc ssc install estout, replace
capture which ivreghdfe
if _rc ssc install ivreghdfe, replace
capture which ivreg2
if _rc ssc install ivreg2, replace
capture which ranktest
if _rc ssc install ranktest, replace

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global output "$datadir/output"
capture mkdir "$output"

log using "$output/phase0_analysis.log", replace

* Load the analysis dataset with self-employment
use "$datadir/analysis_dataset_with_se.dta", clear

********************************************************************************
* PART 1: SAMPLE RESTRICTIONS AND VARIABLE CONSTRUCTION
********************************************************************************

* Restrict to working-age adults (18-64) who are in the labor force
keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1  // In labor force

* Keep years with mobile banking data (2013+)
keep if year >= 2013

* Keep observations in identifiable CBSAs
keep if cbsa > 0 & cbsa != .

di "Sample size after restrictions: " _N

********************************************************************************
* Create additional variables for analysis
********************************************************************************

* Age categories
gen age_cat = .
replace age_cat = 1 if age >= 18 & age <= 29
replace age_cat = 2 if age >= 30 & age <= 44
replace age_cat = 3 if age >= 45 & age <= 64
label define age_cat_lbl 1 "18-29" 2 "30-44" 3 "45-64"
label values age_cat age_cat_lbl

* Create interaction: Mobile banking × employed
gen mobile_X_employed = mobile_user * employed if mobile_user != . & employed != .

* Standardize broadband for interpretability
egen pct_broadband_std = std(pct_broadband)

* Create CBSA-year panel identifier
egen cbsa_year_id = group(cbsa year)

********************************************************************************
* PART 2: DESCRIPTIVE STATISTICS
********************************************************************************

di "=============================================="
di "PART 2: DESCRIPTIVE STATISTICS"
di "=============================================="

* Table 1: Sample characteristics by year
estpost tabstat self_employed mobile_user banked employed age ///
    college_degree inc_above75k metro ///
    [aw=weight], by(year) stat(mean sd n) columns(statistics)
esttab using "$output/table1_sample_by_year.csv", replace ///
    cells("mean(fmt(3)) sd(fmt(3)) count(fmt(0))") ///
    title("Sample Characteristics by Survey Year")

* Table 2: Mobile banking adoption trends
preserve
collapse (mean) mobile_user mobile_primary offsite_only branch_user ///
         (rawsum) n=weight [aw=weight], by(year)
list year mobile_user mobile_primary offsite_only branch_user n
export delimited "$output/table2_mobile_trends.csv", replace
restore

* Table 3: Self-employment by banking mode
table banking_mode year [aw=weight] if banking_mode != ., ///
    stat(mean self_employed) stat(count self_employed) nformat(%9.3f)

* Table 4: Self-employment by mobile banking status
table mobile_user year [aw=weight] if mobile_user != ., ///
    stat(mean self_employed) stat(count self_employed) nformat(%9.3f)

********************************************************************************
* Figure 1: Mobile banking and self-employment trends
********************************************************************************

preserve
collapse (mean) self_employed mobile_user [aw=weight], by(year)

twoway (connected self_employed year, yaxis(1) lcolor(navy) mcolor(navy)) ///
       (connected mobile_user year, yaxis(2) lcolor(maroon) mcolor(maroon)), ///
       ytitle("Self-Employment Rate", axis(1)) ///
       ytitle("Mobile Banking Adoption", axis(2)) ///
       xtitle("Year") ///
       legend(label(1 "Self-Employment") label(2 "Mobile Banking")) ///
       title("Mobile Banking Adoption and Self-Employment Over Time") ///
       scheme(s2color)
graph export "$output/fig1_trends.png", replace width(1200)
restore

********************************************************************************
* PART 3: BASELINE OLS REGRESSIONS
********************************************************************************

di "=============================================="
di "PART 3: BASELINE OLS REGRESSIONS"
di "=============================================="

* Equation (1) from research plan:
* SE_ijt = α + β₁·MobileBanking_ijt + X_ijt'γ + φ_j + λ_t + ε_ijt

* Column 1: No controls
reg self_employed mobile_user [pw=weight], cluster(cbsa)
eststo m1

* Column 2: Individual controls
reg self_employed mobile_user ///
    age c.age#c.age ///
    i.praceeth3 i.peducgrp i.hhincome ///
    [pw=weight], cluster(cbsa)
eststo m2

* Column 3: Add CBSA and year FE
reghdfe self_employed mobile_user ///
    age c.age#c.age ///
    i.praceeth3 i.peducgrp i.hhincome ///
    [pw=weight], absorb(cbsa year) cluster(cbsa)
eststo m3

* Column 4: Add CBSA-level controls
reghdfe self_employed mobile_user ///
    age c.age#c.age ///
    i.praceeth3 i.peducgrp i.hhincome ///
    pct_broadband unemployment_rate ///
    [pw=weight], absorb(cbsa year) cluster(cbsa)
eststo m4

* Export baseline results
esttab m1 m2 m3 m4 using "$output/table3_baseline_ols.csv", replace ///
    keep(mobile_user age pct_broadband unemployment_rate) ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2 r2_a, labels("Observations" "R-squared" "Adj. R-squared")) ///
    title("Baseline OLS: Self-Employment on Mobile Banking") ///
    mtitles("(1)" "(2)" "(3)" "(4)") ///
    addnotes("Standard errors clustered at CBSA level" ///
             "All regressions weighted by survey weights")

********************************************************************************
* PART 4: BANKING MODE AS OUTCOME (For Structural Model)
********************************************************************************

di "=============================================="
di "PART 4: BANKING MODE ANALYSIS"
di "=============================================="

* Multinomial logit: Banking mode choice
* Banking mode: 1=Unbanked, 2=Mobile/Online Only, 3=Branch User

mlogit banking_mode ///
    age c.age#c.age ///
    i.praceeth3 i.peducgrp i.hhincome ///
    metro pct_broadband ///
    [pw=weight] if banking_mode != ., baseoutcome(3) cluster(cbsa)
eststo mlogit1

* Marginal effects
margins, dydx(pct_broadband) predict(outcome(1)) predict(outcome(2)) predict(outcome(3))

********************************************************************************
* PART 5: JOINT DISTRIBUTION (Banking × Employment)
********************************************************************************

di "=============================================="
di "PART 5: JOINT DISTRIBUTION ANALYSIS"
di "=============================================="

* Create joint outcome variable (for structural model)
* 9 categories: 3 banking modes × 3 employment statuses
gen joint_choice = .
* Banking mode 1 (Unbanked) × Employment
replace joint_choice = 1 if banking_mode == 1 & wage_worker == 1
replace joint_choice = 2 if banking_mode == 1 & self_employed == 1
replace joint_choice = 3 if banking_mode == 1 & employed == 0
* Banking mode 2 (Mobile/Online) × Employment
replace joint_choice = 4 if banking_mode == 2 & wage_worker == 1
replace joint_choice = 5 if banking_mode == 2 & self_employed == 1
replace joint_choice = 6 if banking_mode == 2 & employed == 0
* Banking mode 3 (Branch) × Employment
replace joint_choice = 7 if banking_mode == 3 & wage_worker == 1
replace joint_choice = 8 if banking_mode == 3 & self_employed == 1
replace joint_choice = 9 if banking_mode == 3 & employed == 0

label define joint_lbl 1 "Unbanked-Wage" 2 "Unbanked-SE" 3 "Unbanked-NotWork" ///
                       4 "Mobile-Wage" 5 "Mobile-SE" 6 "Mobile-NotWork" ///
                       7 "Branch-Wage" 8 "Branch-SE" 9 "Branch-NotWork"
label values joint_choice joint_lbl

* Distribution of joint choices
tab joint_choice [aw=weight], m

* Joint distribution by year
table joint_choice year [aw=weight] if joint_choice != ., stat(percent) nformat(%9.2f)

********************************************************************************
* PART 6: IV REGRESSION - Broadband as Instrument for Mobile Banking
********************************************************************************

di "=============================================="
di "PART 6: IV REGRESSION"
di "=============================================="

* First stage: Mobile banking on broadband penetration
reghdfe mobile_user pct_broadband ///
    age c.age#c.age ///
    i.praceeth3 i.peducgrp i.hhincome ///
    [pw=weight], absorb(cbsa year) cluster(cbsa)
eststo fs1
predict mobile_hat, xb

* Check first-stage F-statistic
test pct_broadband
local fs_F = r(F)
di "First-stage F-statistic: " `fs_F'

* Second stage: Self-employment on instrumented mobile banking
ivreghdfe self_employed ///
    (mobile_user = pct_broadband) ///
    age c.age#c.age ///
    i.praceeth3 i.peducgrp i.hhincome ///
    [pw=weight], absorb(cbsa year) cluster(cbsa) first
eststo iv1

* Export IV results
esttab fs1 iv1 using "$output/table4_iv_results.csv", replace ///
    keep(mobile_user pct_broadband age) ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2 widstat, labels("Observations" "R-squared" "First-stage F")) ///
    title("IV Regression: Broadband as Instrument for Mobile Banking") ///
    mtitles("First Stage" "Second Stage (IV)")

********************************************************************************
* PART 7: HETEROGENEITY ANALYSIS
********************************************************************************

di "=============================================="
di "PART 7: HETEROGENEITY ANALYSIS"
di "=============================================="

* By race/ethnicity
eststo clear
foreach race in 1 2 6 {
    local lbl: label (praceeth3) `race'
    reghdfe self_employed mobile_user ///
        age c.age#c.age i.peducgrp i.hhincome ///
        [pw=weight] if praceeth3 == `race', absorb(cbsa year) cluster(cbsa)
    eststo het_race`race'
}

esttab het_race1 het_race2 het_race6 using "$output/table5_het_race.csv", replace ///
    keep(mobile_user) se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2) title("Heterogeneity by Race/Ethnicity") ///
    mtitles("Black" "Hispanic" "White")

* By education
eststo clear
foreach edu in 1 2 3 4 {
    local lbl: label (peducgrp) `edu'
    reghdfe self_employed mobile_user ///
        age c.age#c.age i.praceeth3 i.hhincome ///
        [pw=weight] if peducgrp == `edu', absorb(cbsa year) cluster(cbsa)
    eststo het_edu`edu'
}

esttab het_edu1 het_edu2 het_edu3 het_edu4 using "$output/table5_het_education.csv", replace ///
    keep(mobile_user) se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2) title("Heterogeneity by Education") ///
    mtitles("No HS" "HS Diploma" "Some College" "College+")

* By income
eststo clear
foreach inc in 1 2 3 4 5 {
    reghdfe self_employed mobile_user ///
        age c.age#c.age i.praceeth3 i.peducgrp ///
        [pw=weight] if hhincome == `inc', absorb(cbsa year) cluster(cbsa)
    eststo het_inc`inc'
}

esttab het_inc* using "$output/table5_het_income.csv", replace ///
    keep(mobile_user) se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2) title("Heterogeneity by Income") ///
    mtitles("<15K" "15-30K" "30-50K" "50-75K" "75K+")

* By metro status
eststo clear
reghdfe self_employed mobile_user ///
    age c.age#c.age i.praceeth3 i.peducgrp i.hhincome ///
    [pw=weight] if metro == 1, absorb(cbsa year) cluster(cbsa)
eststo het_metro

reghdfe self_employed mobile_user ///
    age c.age#c.age i.praceeth3 i.peducgrp i.hhincome ///
    [pw=weight] if rural == 1, absorb(cbsa year) cluster(cbsa)
eststo het_rural

esttab het_metro het_rural using "$output/table5_het_geography.csv", replace ///
    keep(mobile_user) se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2) title("Heterogeneity by Geography") ///
    mtitles("Metropolitan" "Rural")

********************************************************************************
* PART 8: CBSA-LEVEL ANALYSIS (For Aggregated Patterns)
********************************************************************************

di "=============================================="
di "PART 8: CBSA-LEVEL ANALYSIS"
di "=============================================="

* Collapse to CBSA-year level
preserve
collapse (mean) self_employed mobile_user banked ///
         pct_broadband unemployment_rate median_hh_income ///
         (rawsum) pop=weight [aw=weight], by(cbsa year)

* CBSA-level regression
xtset cbsa year
xtreg self_employed mobile_user pct_broadband unemployment_rate i.year, fe cluster(cbsa)
eststo cbsa1

* Export CBSA-level results
esttab cbsa1 using "$output/table6_cbsa_level.csv", replace ///
    keep(mobile_user pct_broadband unemployment_rate) ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2_w r2_b r2_o) title("CBSA-Level Panel Regression")

restore

********************************************************************************
* PART 9: CONDITIONAL CHOICE PROBABILITIES (For Phase 2)
********************************************************************************

di "=============================================="
di "PART 9: CCPs FOR STRUCTURAL MODEL"
di "=============================================="

* Estimate CCPs by CBSA × year × demographic cell
* Cells defined by: age_cat × education × race

preserve
collapse (mean) prob_se=self_employed prob_mobile=mobile_user ///
         (count) n=self_employed [aw=weight], ///
         by(cbsa year age_cat peducgrp praceeth3)

* Save CCPs for Phase 2
save "$output/ccps_cbsa_year_demo.dta", replace
export delimited "$output/ccps_cbsa_year_demo.csv", replace

* Summary of CCP variation
summarize prob_se prob_mobile n
restore

********************************************************************************
* SUMMARY STATISTICS FOR PAPER
********************************************************************************

di "=============================================="
di "SUMMARY FOR PAPER"
di "=============================================="

* Key statistics
summarize self_employed mobile_user [aw=weight]

* Mobile banking adoption growth
table year [aw=weight], stat(mean mobile_user) nformat(%9.3f)

* Self-employment by banking mode
table banking_mode [aw=weight] if banking_mode != ., stat(mean self_employed) nformat(%9.3f)

di "=============================================="
di "PHASE 0 ANALYSIS COMPLETE"
di "=============================================="
di "Output files saved to: $output"
di "=============================================="

log close

********************************************************************************
* Export final analysis dataset for Phase 1
********************************************************************************

save "$datadir/analysis_sample_phase0.dta", replace
