/*******************************************************************************
* Continue merging: ACS data and final dataset
*******************************************************************************/

clear all
set more off

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global acs "$datadir/ACS"

log using "$datadir/finish_merge.log", replace

********************************************************************************
* Process ACS Data
********************************************************************************

* Import ACS broadband data
import delimited "$acs/broadband_cbsa_2023.csv", clear varnames(1)

* List variable names
describe

* The last column is CBSA code - rename it using wildcard
rename metropolitan* cbsa

* Rename other variables
rename name cbsa_name
rename s2801_c01_001e total_households
rename s2801_c02_012e pct_broadband
rename s2801_c02_014e pct_cellular

* Clean
destring cbsa total_households pct_broadband pct_cellular, replace force
keep cbsa cbsa_name total_households pct_broadband pct_cellular
save "$acs/acs_broadband.dta", replace

* Import ACS demographics data
import delimited "$acs/demographics_cbsa_2023.csv", clear varnames(1)

rename metropolitan* cbsa
destring cbsa, replace force

* Rename demographic variables
capture rename dp05_0001e population
capture rename dp05_0037pe pct_white_alone
capture rename dp05_0038pe pct_black_alone
capture rename dp05_0039pe pct_asian_alone
capture rename dp05_0044pe pct_two_more_races
capture rename dp05_0071pe pct_hispanic

destring population pct_*, replace force
keep cbsa name population pct_*
save "$acs/acs_demographics.dta", replace

* Import ACS income/unemployment data
import delimited "$acs/income_unemployment_cbsa_2023.csv", clear varnames(1)

rename metropolitan* cbsa
destring cbsa, replace force

capture rename s1901_c01_012e median_hh_income
capture rename s2301_c04_001e unemployment_rate

destring median_hh_income unemployment_rate, replace force
keep cbsa name median_hh_income unemployment_rate
save "$acs/acs_income_unemp.dta", replace

* Import ACS industry data
import delimited "$acs/industry_cbsa_2023.csv", clear varnames(1)

rename metropolitan* cbsa
destring cbsa, replace force

keep cbsa name s2403_*
destring s2403_*, replace force
save "$acs/acs_industry.dta", replace

* Merge all ACS data
use "$acs/acs_broadband.dta", clear
merge 1:1 cbsa using "$acs/acs_demographics.dta", nogen
merge 1:1 cbsa using "$acs/acs_income_unemp.dta", nogen
merge 1:1 cbsa using "$acs/acs_industry.dta", nogen
save "$acs/acs_cbsa_controls.dta", replace

di "ACS data processed: " _N " CBSAs"

********************************************************************************
* Merge with FDIC data
********************************************************************************

use "$datadir/fdic_individual_temp.dta", clear

* Merge ACS controls
merge m:1 cbsa using "$acs/acs_cbsa_controls.dta", keep(master match) nogen

* Label key variables
label var mobile_primary "Mobile banking is primary access method"
label var mobile_user "Used mobile banking in past 12 months"
label var offsite_only "Uses only off-site channels (no branch visits)"
label var branch_user "Uses bank teller at branch"
label var banked "Has bank account"
label var unbanked "Does not have bank account"
label var employed "Employed"
label var unemployed "Unemployed"
label var nilf "Not in labor force"
label var black "Black/African American"
label var hispanic "Hispanic/Latino"
label var asian "Asian"
label var white "White"
label var metro "Lives in metropolitan area"
label var rural "Lives in rural area"
label var college_degree "Has college degree"
label var homeowner "Homeowner"
capture label var pct_broadband "CBSA % households with broadband (ACS)"
capture label var median_hh_income "CBSA median household income (ACS)"
capture label var unemployment_rate "CBSA unemployment rate (ACS)"

********************************************************************************
* Summary Statistics
********************************************************************************

* Sample size by year
di "=== Sample by Year ==="
tab year, m

* Mobile banking by year
di "=== Mobile Banking Adoption by Year ==="
table year, stat(mean mobile_user mobile_primary) stat(count mobile_user)

* Banking status by year
di "=== Banking Status by Year ==="
table year, stat(mean banked unbanked)

********************************************************************************
* Save Final Dataset
********************************************************************************

* Order variables
order cbsa year weight ///
      mobile_primary mobile_user offsite_only branch_user banking_mode ///
      banked unbanked employed unemployed nilf ///
      age black hispanic asian white other_race ///
      no_hs hs_diploma some_college college_degree ///
      inc_below15k inc_15_30k inc_30_50k inc_50_75k inc_above75k ///
      metro rural homeowner ///
      pct_broadband pct_cellular median_hh_income unemployment_rate

compress
save "$datadir/analysis_dataset.dta", replace
export delimited "$datadir/analysis_dataset.csv", replace

* Final summary
describe
summarize mobile_user mobile_primary banked employed age

di "=============================================="
di "Analysis dataset created successfully!"
di "N = " _N " observations"
di "=============================================="

log close
