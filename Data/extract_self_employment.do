/*******************************************************************************
* Extract Self-Employment from Raw CPS Data
*
* This do-file parses the raw CPS fixed-width files to extract the class of
* worker variable (PEIO1COW) which identifies self-employment status.
*
* PEIO1COW values:
* 1 = Government - Federal
* 2 = Government - State
* 3 = Government - Local
* 4 = Private, for profit
* 5 = Private, nonprofit
* 6 = Self-employed, incorporated
* 7 = Self-employed, unincorporated
* 8 = Without pay
*******************************************************************************/

clear all
set more off

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global fdic "$datadir/FDIC_Survey/hhmultiyear"
global rawdir "$fdic/src/raw"

log using "$datadir/extract_self_employment.log", replace

********************************************************************************
* Function to parse raw CPS file
* Column positions from layout files:
* - HRHHID: 1-15 (household ID)
* - HRHHID2: 71-75 (household ID part 2)
* - PULINENO: 147-148 (person line number)
* - PEIO1COW: varies by year (class of worker, primary job)
* - GTCBSA: varies by year (CBSA code)
* - HRYEAR4: 18-21 (survey year)
********************************************************************************

* Process 2023 data
infix str hrhhid_str 1-15 hryear4 18-21 str hrhhid2_str 71-75 ///
      pulineno 147-148 ///
      gtcbsa 96-100 ///
      peio1cow 432-433 ///
      using "$rawdir/jun23pub.dat", clear

* Convert string IDs to numeric to match FDIC data
destring hrhhid_str, gen(hrhhid) force
destring hrhhid2_str, gen(hrhhid2) force
drop hrhhid_str hrhhid2_str

* Create unique person ID
gen long person_id = hrhhid * 1000000 + hrhhid2 * 100 + pulineno

* Create self-employment indicators
gen self_employed = (peio1cow == 6 | peio1cow == 7) if peio1cow > 0 & peio1cow < 9
gen self_emp_inc = (peio1cow == 6) if peio1cow > 0 & peio1cow < 9
gen self_emp_uninc = (peio1cow == 7) if peio1cow > 0 & peio1cow < 9
gen wage_worker = (peio1cow >= 1 & peio1cow <= 5) if peio1cow > 0 & peio1cow < 9

* Keep relevant variables
keep hrhhid hrhhid2 pulineno hryear4 gtcbsa peio1cow ///
     self_employed self_emp_inc self_emp_uninc wage_worker person_id

gen year = 2023
save "$datadir/self_emp_2023.dta", replace

* Process 2021 data
infix str hrhhid_str 1-15 hryear4 18-21 str hrhhid2_str 71-75 ///
      pulineno 147-148 ///
      gtcbsa 96-100 ///
      peio1cow 432-433 ///
      using "$rawdir/jun21pub.dat", clear

destring hrhhid_str, gen(hrhhid) force
destring hrhhid2_str, gen(hrhhid2) force
drop hrhhid_str hrhhid2_str
gen long person_id = hrhhid * 1000000 + hrhhid2 * 100 + pulineno
gen self_employed = (peio1cow == 6 | peio1cow == 7) if peio1cow > 0 & peio1cow < 9
gen self_emp_inc = (peio1cow == 6) if peio1cow > 0 & peio1cow < 9
gen self_emp_uninc = (peio1cow == 7) if peio1cow > 0 & peio1cow < 9
gen wage_worker = (peio1cow >= 1 & peio1cow <= 5) if peio1cow > 0 & peio1cow < 9
keep hrhhid hrhhid2 pulineno hryear4 gtcbsa peio1cow ///
     self_employed self_emp_inc self_emp_uninc wage_worker person_id
gen year = 2021
save "$datadir/self_emp_2021.dta", replace

* Process 2019 data
infix str hrhhid_str 1-15 hryear4 18-21 str hrhhid2_str 71-75 ///
      pulineno 147-148 ///
      gtcbsa 96-100 ///
      peio1cow 432-433 ///
      using "$rawdir/jun19pub.dat", clear

destring hrhhid_str, gen(hrhhid) force
destring hrhhid2_str, gen(hrhhid2) force
drop hrhhid_str hrhhid2_str
gen long person_id = hrhhid * 1000000 + hrhhid2 * 100 + pulineno
gen self_employed = (peio1cow == 6 | peio1cow == 7) if peio1cow > 0 & peio1cow < 9
gen self_emp_inc = (peio1cow == 6) if peio1cow > 0 & peio1cow < 9
gen self_emp_uninc = (peio1cow == 7) if peio1cow > 0 & peio1cow < 9
gen wage_worker = (peio1cow >= 1 & peio1cow <= 5) if peio1cow > 0 & peio1cow < 9
keep hrhhid hrhhid2 pulineno hryear4 gtcbsa peio1cow ///
     self_employed self_emp_inc self_emp_uninc wage_worker person_id
gen year = 2019
save "$datadir/self_emp_2019.dta", replace

* Process 2017 data
* Note: Column positions may differ - check 2017 layout file
infix str hrhhid_str 1-15 hryear4 18-21 str hrhhid2_str 71-75 ///
      pulineno 147-148 ///
      gtcbsa 96-100 ///
      peio1cow 432-433 ///
      using "$rawdir/jun17pubi.dat", clear

destring hrhhid_str, gen(hrhhid) force
destring hrhhid2_str, gen(hrhhid2) force
drop hrhhid_str hrhhid2_str
gen long person_id = hrhhid * 1000000 + hrhhid2 * 100 + pulineno
gen self_employed = (peio1cow == 6 | peio1cow == 7) if peio1cow > 0 & peio1cow < 9
gen self_emp_inc = (peio1cow == 6) if peio1cow > 0 & peio1cow < 9
gen self_emp_uninc = (peio1cow == 7) if peio1cow > 0 & peio1cow < 9
gen wage_worker = (peio1cow >= 1 & peio1cow <= 5) if peio1cow > 0 & peio1cow < 9
keep hrhhid hrhhid2 pulineno hryear4 gtcbsa peio1cow ///
     self_employed self_emp_inc self_emp_uninc wage_worker person_id
gen year = 2017
save "$datadir/self_emp_2017.dta", replace

* Process 2015 data
infix str hrhhid_str 1-15 hryear4 18-21 str hrhhid2_str 71-75 ///
      pulineno 147-148 ///
      gtcbsa 96-100 ///
      peio1cow 432-433 ///
      using "$rawdir/jun15pubi.dat", clear

destring hrhhid_str, gen(hrhhid) force
destring hrhhid2_str, gen(hrhhid2) force
drop hrhhid_str hrhhid2_str
gen long person_id = hrhhid * 1000000 + hrhhid2 * 100 + pulineno
gen self_employed = (peio1cow == 6 | peio1cow == 7) if peio1cow > 0 & peio1cow < 9
gen self_emp_inc = (peio1cow == 6) if peio1cow > 0 & peio1cow < 9
gen self_emp_uninc = (peio1cow == 7) if peio1cow > 0 & peio1cow < 9
gen wage_worker = (peio1cow >= 1 & peio1cow <= 5) if peio1cow > 0 & peio1cow < 9
keep hrhhid hrhhid2 pulineno hryear4 gtcbsa peio1cow ///
     self_employed self_emp_inc self_emp_uninc wage_worker person_id
gen year = 2015
save "$datadir/self_emp_2015.dta", replace

* Process 2013 data
infix str hrhhid_str 1-15 hryear4 18-21 str hrhhid2_str 71-75 ///
      pulineno 147-148 ///
      gtcbsa 96-100 ///
      peio1cow 432-433 ///
      using "$rawdir/jun13pubi.dat", clear

destring hrhhid_str, gen(hrhhid) force
destring hrhhid2_str, gen(hrhhid2) force
drop hrhhid_str hrhhid2_str
gen long person_id = hrhhid * 1000000 + hrhhid2 * 100 + pulineno
gen self_employed = (peio1cow == 6 | peio1cow == 7) if peio1cow > 0 & peio1cow < 9
gen self_emp_inc = (peio1cow == 6) if peio1cow > 0 & peio1cow < 9
gen self_emp_uninc = (peio1cow == 7) if peio1cow > 0 & peio1cow < 9
gen wage_worker = (peio1cow >= 1 & peio1cow <= 5) if peio1cow > 0 & peio1cow < 9
keep hrhhid hrhhid2 pulineno hryear4 gtcbsa peio1cow ///
     self_employed self_emp_inc self_emp_uninc wage_worker person_id
gen year = 2013
save "$datadir/self_emp_2013.dta", replace

********************************************************************************
* Append all years
********************************************************************************

use "$datadir/self_emp_2023.dta", clear
append using "$datadir/self_emp_2021.dta"
append using "$datadir/self_emp_2019.dta"
append using "$datadir/self_emp_2017.dta"
append using "$datadir/self_emp_2015.dta"
append using "$datadir/self_emp_2013.dta"

* Label variables
label var self_employed "Self-employed (any)"
label var self_emp_inc "Self-employed, incorporated"
label var self_emp_uninc "Self-employed, unincorporated"
label var wage_worker "Wage and salary worker"
label var peio1cow "Class of worker (primary job)"

label define cowlbl 1 "Federal govt" 2 "State govt" 3 "Local govt" ///
                    4 "Private for-profit" 5 "Private nonprofit" ///
                    6 "Self-emp incorporated" 7 "Self-emp unincorporated" ///
                    8 "Without pay"
label values peio1cow cowlbl

* Summary
tab year self_employed, row
tab peio1cow year, col

* Save combined file
compress
save "$datadir/self_employment_all_years.dta", replace

********************************************************************************
* Merge with FDIC analysis data
********************************************************************************

* The FDIC analysis file has: hrhhid, hrhhid2, pulineno (from variable list)
* We can merge on these identifiers plus year

* First, check if merge keys exist in FDIC data
use "$datadir/fdic_individual_temp.dta", clear
describe hrhhid hrhhid2 pulineno

* If merge keys exist, merge self-employment
merge 1:1 hrhhid hrhhid2 pulineno year using "$datadir/self_employment_all_years.dta", ///
      keep(master match) nogen keepusing(self_employed self_emp_inc self_emp_uninc wage_worker peio1cow)

* Create employment categories for structural model
* d ∈ {wage employment, self-employment, not working}
gen emp_status = .
replace emp_status = 1 if wage_worker == 1        // Wage employment
replace emp_status = 2 if self_employed == 1      // Self-employment
replace emp_status = 3 if employed == 0           // Not working (unemployed or NILF)
label define emp_status_lbl 1 "Wage worker" 2 "Self-employed" 3 "Not working"
label values emp_status emp_status_lbl

label var emp_status "Employment status (3 categories)"

* Save updated dataset
save "$datadir/analysis_dataset_with_se.dta", replace

********************************************************************************
* Summary statistics on self-employment
********************************************************************************

* Self-employment rates by year
tab year self_employed [aw=weight], row

* Self-employment by demographics
tab praceeth3 self_employed [aw=weight], row
tab peducgrp self_employed [aw=weight], row
tab hhincome self_employed [aw=weight], row

* Self-employment by banking status
tab banking_mode self_employed [aw=weight] if banking_mode != ., row

* Mobile banking and self-employment
tab mobile_user self_employed [aw=weight] if mobile_user != ., row

di "=============================================="
di "Self-employment data merged successfully!"
di "=============================================="

log close

* Clean up temporary files
erase "$datadir/self_emp_2013.dta"
erase "$datadir/self_emp_2015.dta"
erase "$datadir/self_emp_2017.dta"
erase "$datadir/self_emp_2019.dta"
erase "$datadir/self_emp_2021.dta"
erase "$datadir/self_emp_2023.dta"
