/*******************************************************************************
* Mobile Banking, Bank Branch Closures, and Self-Employment in the United States
*
* Construct Merged Analysis Dataset
* Author: Research Project
* Date: February 2025
*
* This do-file merges:
*   1. FDIC/CPS Survey data (individual-level)
*   2. FDIC SOD branch data (CBSA-year level)
*   3. ACS controls (CBSA-level)
*******************************************************************************/

clear all
set more off
set maxvar 10000

* Set paths
global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global fdic "$datadir/FDIC_Survey/hhmultiyear"
global sod "$datadir/SOD"
global acs "$datadir/ACS"
global jmp "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/JMP"

log using "$datadir/construct_analysis_dataset.log", replace

********************************************************************************
* PART 1: Load and Process FDIC/CPS Multi-Year Survey Data
********************************************************************************

* Load the FDIC multi-year analysis dataset
import delimited "$fdic/hh_multiyear_analys.csv", clear

* Key variables from FDIC survey:
* - hryear4: Survey year
* - msa13: CBSA code (2013 delineations)
* - gestfips: State FIPS code
* - gtcbsast: Metropolitan status (1=principal city, 2=balance, 3=non-metro)
* - gtmetsta: Metro/non-metro (1=metro, 2=non-metro, 3=not identified)
* - hbankstatv5/v6: Bank account status (1=unbanked, 2=underbanked, 3=fully banked)
* - hbnkaccm: Most common way to access account (5=mobile banking)
* - hbnkaccm5: Used mobile banking to access account
* - hmobonltelstat: Used off-site channels (mobile/online/phone)
* - hbranchstat: Used bank teller (branch visits)
* - pempstat: Employment status (1=employed, 2=unemployed, 3=NILF)
* - praceeth3: Race/ethnicity
* - peducgrp: Education (1-4)
* - hhincome: Family income category (1-5)
* - prtage: Age
* - hhsupwgt: Household survey weight

* Rename key variables for clarity
rename msa13 cbsa
rename hryear4 year
rename gestfips statefips
rename hhsupwgt weight

* Create CBSA-year identifier
gen cbsa_year = cbsa * 10000 + year

* Label survey year
label define yearlbl 2009 "2009" 2011 "2011" 2013 "2013" 2015 "2015" ///
                     2017 "2017" 2019 "2019" 2021 "2021" 2023 "2023"
label values year yearlbl

********************************************************************************
* Create Key Analysis Variables
********************************************************************************

* Mobile banking indicator (primary access method = mobile)
* hbnkaccm = 5 means mobile banking is most common access method
gen mobile_primary = (hbnkaccm == 5) if hbnkaccm != .

* Mobile banking user (any use)
* hbnkaccm5 = 1 means used mobile banking
gen mobile_user = (hbnkaccm5 == 1) if hbnkaccm5 > 0 & hbnkaccm5 != .

* Off-site banking channels only (mobile/online/phone, no branch visits)
gen offsite_only = (hmobonltelstat == 1) if hmobonltelstat > 0 & hmobonltelstat != .
replace offsite_only = (hmobonltelstatv2 == 1) if offsite_only == . & hmobonltelstatv2 > 0 & hmobonltelstatv2 != .

* Branch user indicator
gen branch_user = (hbranchstat == 1 | hbranchstat == 2) if hbranchstat > 0 & hbranchstat != .
replace branch_user = (hbranchstatv2 == 1 | hbranchstatv2 == 2) if branch_user == . & hbranchstatv2 > 0 & hbranchstatv2 != .

* Banked status
gen banked = (hbankstatv6 == 2 | hbankstatv6 == 3) if hbankstatv6 != .
replace banked = (hbankstatv5 == 2 | hbankstatv5 == 3) if banked == . & hbankstatv5 != .
replace banked = (hbankstatv4 == 2 | hbankstatv4 == 3) if banked == . & hbankstatv4 != .
replace banked = (hbankstatv3 == 2 | hbankstatv3 == 3) if banked == . & hbankstatv3 != .
replace banked = (hbankstatv2 == 2 | hbankstatv2 == 3) if banked == . & hbankstatv2 != .
replace banked = (hbankstat == 2 | hbankstat == 3) if banked == . & hbankstat != .

* Unbanked indicator
gen unbanked = (banked == 0) if banked != .

* Employment status
gen employed = (pempstat == 1) if pempstat != .
gen unemployed = (pempstat == 2) if pempstat != .
gen nilf = (pempstat == 3) if pempstat != .

* Race/ethnicity (using praceeth3 - consistent definition across waves)
gen black = (praceeth3 == 1) if praceeth3 != .
gen hispanic = (praceeth3 == 2) if praceeth3 != .
gen asian = (praceeth3 == 3) if praceeth3 != .
gen white = (praceeth3 == 6) if praceeth3 != .
gen other_race = (praceeth3 == 4 | praceeth3 == 5 | praceeth3 == 7) if praceeth3 != .

* Education categories
gen no_hs = (peducgrp == 1) if peducgrp != .
gen hs_diploma = (peducgrp == 2) if peducgrp != .
gen some_college = (peducgrp == 3) if peducgrp != .
gen college_degree = (peducgrp == 4) if peducgrp != .

* Income categories
gen inc_below15k = (hhincome == 1) if hhincome != .
gen inc_15_30k = (hhincome == 2) if hhincome != .
gen inc_30_50k = (hhincome == 3) if hhincome != .
gen inc_50_75k = (hhincome == 4) if hhincome != .
gen inc_above75k = (hhincome == 5) if hhincome != .

* Age
rename prtage age

* Metro status
gen metro = (gtmetsta == 1) if gtmetsta != . & gtmetsta != 3
gen rural = (gtmetsta == 2) if gtmetsta != . & gtmetsta != 3

* Homeownership
gen homeowner = (hhtenure == 1) if hhtenure != .

* Banking mode categories (for structural model)
gen banking_mode = .
replace banking_mode = 1 if unbanked == 1  // Unbanked
replace banking_mode = 2 if banked == 1 & offsite_only == 1  // Mobile/online only
replace banking_mode = 3 if banked == 1 & branch_user == 1  // Branch user
label define bank_mode 1 "Unbanked" 2 "Mobile/Online Only" 3 "Branch User"
label values banking_mode bank_mode

* Save individual-level data
compress
save "$datadir/fdic_individual_temp.dta", replace

********************************************************************************
* PART 2: Process SOD Branch Data
********************************************************************************

* Load existing SOD data (CBSA-level branch counts)
use "$jmp/Revision/New Data/SOD_zip_bank_branch.dta", clear

* Explore the structure
describe
summarize

* The SOD data should have: year, geography (zip/cbsa), number of branches
* We need to aggregate to CBSA-year level if not already

* Check if we need to aggregate from ZIP to CBSA
* If the data is at ZIP level, we need a ZIP-to-CBSA crosswalk

* For now, save what we have
save "$sod/sod_branch_data.dta", replace

********************************************************************************
* PART 2B: Alternative - Create CBSA-level branch data from raw SOD
* (Run this section if you need to download/process raw SOD data)
********************************************************************************

* Note: Raw SOD data should be downloaded from https://www.fdic.gov/resources/sod/
* Download "Branch Office Deposits" files for each year (2009-2023)
* The files contain: CERT, BRNUM, UNINUMBR, DESSION, YEAR, STCNTY, SIMS_LATITUDE,
*                    SIMS_LONGITUDE, STNAME, ZIPBR, CBSA, CBSA_NAME, etc.

********************************************************************************
* PART 3: Import and Process ACS Data
********************************************************************************

* Import ACS broadband data
import delimited "$acs/broadband_cbsa_2023.csv", clear varnames(1)

* List variable names to see what Stata imported
describe, short

* Rename variables - use v5 as CBSA code (5th column)
* Column order: NAME, S2801_C01_001E, S2801_C02_012E, S2801_C02_014E, cbsa_code
capture rename v5 cbsa_code
capture rename metropolitanstatisticalareamicropoli cbsa_code
capture rename metropolitanstatisticareamicropoli cbsa_code

* Rename other variables
capture rename name cbsa_name
capture rename s2801_c01_001e total_households
capture rename s2801_c02_012e pct_broadband
capture rename s2801_c02_014e pct_cellular

* Clean CBSA code
destring cbsa_code, replace force
destring total_households pct_broadband pct_cellular, replace force

* Keep relevant variables
keep cbsa_code cbsa_name total_households pct_broadband pct_cellular
rename cbsa_code cbsa

* Save broadband data
save "$acs/acs_broadband.dta", replace

* Import ACS demographics data
import delimited "$acs/demographics_cbsa_2023.csv", clear varnames(1)

capture rename v5 cbsa_code
capture rename metropolitanstatisticalareamicropoli cbsa_code
capture rename metropolitanstatisticareamicropoli cbsa_code
destring cbsa_code, replace force

* Rename demographic variables
capture rename dp05_0001e population
capture rename dp05_0037pe pct_white_alone
capture rename dp05_0038pe pct_black_alone
capture rename dp05_0039pe pct_asian_alone
capture rename dp05_0044pe pct_two_more_races
capture rename dp05_0071pe pct_hispanic

* Clean and save
destring population pct_*, replace force
keep cbsa_code name population pct_*
rename cbsa_code cbsa

save "$acs/acs_demographics.dta", replace

* Import ACS income/unemployment data
import delimited "$acs/income_unemployment_cbsa_2023.csv", clear varnames(1)

capture rename v4 cbsa_code
capture rename metropolitanstatisticalareamicropoli cbsa_code
capture rename metropolitanstatisticareamicropoli cbsa_code
destring cbsa_code, replace force

capture rename s1901_c01_012e median_hh_income
capture rename s2301_c04_001e unemployment_rate

destring median_hh_income unemployment_rate, replace force
keep cbsa_code name median_hh_income unemployment_rate
rename cbsa_code cbsa

save "$acs/acs_income_unemp.dta", replace

* Import ACS industry data
import delimited "$acs/industry_cbsa_2023.csv", clear varnames(1)

capture rename v7 cbsa_code
capture rename metropolitanstatisticalareamicropoli cbsa_code
capture rename metropolitanstatisticareamicropoli cbsa_code
destring cbsa_code, replace force

* Keep employment by industry variables
keep cbsa_code name s2403_*
rename cbsa_code cbsa

destring s2403_*, replace force

save "$acs/acs_industry.dta", replace

* Merge all ACS data
use "$acs/acs_broadband.dta", clear
merge 1:1 cbsa using "$acs/acs_demographics.dta", nogen
merge 1:1 cbsa using "$acs/acs_income_unemp.dta", nogen
merge 1:1 cbsa using "$acs/acs_industry.dta", nogen

* Note: ACS data is 2023 only; will be merged to all survey years as controls
* For panel variation, would need to download ACS for each year

save "$acs/acs_cbsa_controls.dta", replace

********************************************************************************
* PART 4: Merge All Datasets
********************************************************************************

* Start with individual-level FDIC data
use "$datadir/fdic_individual_temp.dta", clear

* Merge ACS CBSA-level controls
* Note: Using 2023 ACS as cross-sectional controls
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
label var pct_broadband "CBSA % households with broadband (ACS)"
label var median_hh_income "CBSA median household income (ACS)"
label var unemployment_rate "CBSA unemployment rate (ACS)"

********************************************************************************
* PART 5: Create Summary Statistics
********************************************************************************

* Overall sample size by year
tab year, m

* Mobile banking adoption over time
table year, stat(mean mobile_user mobile_primary) stat(count mobile_user)

* Banking status by year
table year, stat(mean banked unbanked)

* Mobile banking by demographics
table praceeth3, stat(mean mobile_user mobile_primary) stat(count mobile_user)
table peducgrp, stat(mean mobile_user mobile_primary) stat(count mobile_user)
table hhincome, stat(mean mobile_user mobile_primary) stat(count mobile_user)

* Mobile banking by metro status
table gtmetsta if gtmetsta < 3, stat(mean mobile_user mobile_primary) stat(count mobile_user)

********************************************************************************
* PART 6: Save Final Analysis Dataset
********************************************************************************

* Order key variables first
order cbsa year weight ///
      mobile_primary mobile_user offsite_only branch_user banking_mode ///
      banked unbanked employed unemployed nilf ///
      age black hispanic asian white other_race ///
      no_hs hs_diploma some_college college_degree ///
      inc_below15k inc_15_30k inc_30_50k inc_50_75k inc_above75k ///
      metro rural homeowner ///
      pct_broadband pct_cellular median_hh_income unemployment_rate

* Compress and save
compress
save "$datadir/analysis_dataset.dta", replace

* Also export to CSV for use in R/Python
export delimited "$datadir/analysis_dataset.csv", replace

* Summary
describe
summarize mobile_user mobile_primary banked employed age

di "=============================================="
di "Analysis dataset created successfully!"
di "N = " _N " observations"
di "Years: 2009-2023 (biennial)"
di "=============================================="

log close

********************************************************************************
* NOTES ON SELF-EMPLOYMENT:
*
* The FDIC analysis file includes pempstat (employed/unemployed/NILF) but
* NOT the class of worker variable needed to identify self-employment.
*
* To add self-employment:
* 1. The raw CPS files in hhmultiyear/src/raw/ contain PEIO1COW (class of worker)
* 2. Values: 1-2 = Government, 3-4 = Private, 5 = Self-employed (not inc),
*            6 = Self-employed (inc), 7 = Without pay
* 3. Need to parse the fixed-width raw files using the layout files
*    and merge by household ID (HRHHID, HRHHID2, PULINENO)
*
* See separate do-file: extract_self_employment.do
********************************************************************************
