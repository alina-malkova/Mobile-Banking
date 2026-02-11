/*******************************************************************************
* Fix the merged dataset: Add ACS controls to self-employment data
*******************************************************************************/

clear all
set more off

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
global acs "$datadir/ACS"

* Load the dataset with self-employment
use "$datadir/analysis_dataset_with_se.dta", clear

* Check if ACS variables exist
capture confirm variable pct_broadband
if _rc {
    di "ACS variables missing - merging now..."

    * Merge ACS controls
    merge m:1 cbsa using "$acs/acs_cbsa_controls.dta", keep(master match) nogen

    * Label new variables
    capture label var pct_broadband "CBSA % households with broadband (ACS)"
    capture label var median_hh_income "CBSA median household income (ACS)"
    capture label var unemployment_rate "CBSA unemployment rate (ACS)"
    capture label var pct_cellular "CBSA % with cellular data (ACS)"
    capture label var population "CBSA population (ACS)"

    * Save updated dataset
    compress
    save "$datadir/analysis_dataset_with_se.dta", replace

    di "ACS variables merged successfully!"
}
else {
    di "ACS variables already present."
}

* Verify
describe pct_broadband median_hh_income unemployment_rate
summarize pct_broadband median_hh_income unemployment_rate
