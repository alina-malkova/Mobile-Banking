/*******************************************************************************
* Extract Demographics (Sex, Marital Status, Children) from Raw CPS Data
*
* Parses PESEX, PEMARITL, PRNMCHLD from raw CPS fixed-width files and merges
* into analysis_dataset_with_se.dta.
*
* Column positions (consistent 2013-2023, verified from layout files):
*   PESEX:    129-130  (1=Male, 2=Female)
*   PEMARITL: 125-126  (1=Married spouse present, 2=Married spouse absent,
*                        3=Widowed, 4=Divorced, 5=Separated, 6=Never married)
*   PRNMCHLD: 635-636  (Number of own children <18 in household)
*
* Follows same pattern as extract_self_employment.do
*******************************************************************************/

clear all
set more off

global datadir "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/_Research/Mobile_Money_Banking/Mobile banking USA/Data"
global fdic "$datadir/FDIC_Survey/hhmultiyear"
global rawdir "$fdic/src/raw"

capture log close
log using "$datadir/output/extract_demographics.log", replace

di "============================================================"
di "EXTRACTING DEMOGRAPHICS FROM RAW CPS FILES"
di "============================================================"

********************************************************************************
* Parse each year
********************************************************************************

foreach yr in 2023 2021 2019 2017 2015 2013 {

    * Determine file name
    local yy = substr("`yr'", 3, 2)

    if `yr' == 2023 {
        local fname "jun23pub.dat"
    }
    else if `yr' == 2021 {
        local fname "jun21pub.dat"
    }
    else if `yr' == 2019 {
        local fname "jun19pub.dat"
    }
    else if `yr' == 2017 {
        local fname "jun17pubi.dat"
    }
    else if `yr' == 2015 {
        local fname "jun15pubi.dat"
    }
    else if `yr' == 2013 {
        local fname "jun13pubi.dat"
    }

    di _n "Processing `yr' (`fname')..."

    infix str hrhhid_str 1-15 hryear4 18-21 str hrhhid2_str 71-75 ///
          pulineno 147-148 ///
          pesex 129-130 ///
          pemaritl 125-126 ///
          prnmchld 635-636 ///
          using "$rawdir/`fname'", clear

    * Convert string IDs to numeric (match FDIC data format)
    destring hrhhid_str, gen(hrhhid) force
    destring hrhhid2_str, gen(hrhhid2) force
    drop hrhhid_str hrhhid2_str

    * Create demographic indicators
    gen female = (pesex == 2) if pesex >= 1 & pesex <= 2
    gen married = (pemaritl == 1 | pemaritl == 2) if pemaritl >= 1 & pemaritl <= 6
    gen has_children = (prnmchld > 0 & prnmchld < 99) if prnmchld >= 0

    gen year = `yr'

    keep hrhhid hrhhid2 pulineno year pesex pemaritl prnmchld female married has_children

    di "  N = " _N " persons"
    di "  Female: " %5.1f 100*r(mean) "%" _col(30) "(after: sum female [aw=1])"
    quietly summarize female
    di "  Female rate: " %5.1f 100*r(mean) "%"
    quietly summarize married
    di "  Married rate: " %5.1f 100*r(mean) "%"
    quietly summarize has_children
    di "  Has children rate: " %5.1f 100*r(mean) "%"

    save "$datadir/demographics_`yr'.dta", replace
}

********************************************************************************
* Append all years
********************************************************************************

di _n "============================================================"
di "APPENDING ALL YEARS"
di "============================================================"

use "$datadir/demographics_2023.dta", clear
foreach yr in 2021 2019 2017 2015 2013 {
    append using "$datadir/demographics_`yr'.dta"
}

* Label variables
label var pesex "Sex (1=Male, 2=Female)"
label var pemaritl "Marital status"
label var prnmchld "Number of own children <18"
label var female "Female indicator"
label var married "Married indicator"
label var has_children "Has own children <18"

label define sexlbl 1 "Male" 2 "Female"
label values pesex sexlbl

label define marlbl 1 "Married, spouse present" 2 "Married, spouse absent" ///
                    3 "Widowed" 4 "Divorced" 5 "Separated" 6 "Never married"
label values pemaritl marlbl

di "Total observations: " _N

tab year female, row
tab year married, row

compress
save "$datadir/demographics_all_years.dta", replace

********************************************************************************
* Merge into analysis dataset
********************************************************************************

di _n "============================================================"
di "MERGING INTO ANALYSIS DATASET"
di "============================================================"

use "$datadir/analysis_dataset_with_se.dta", clear
di "Analysis dataset N = " _N

* Check if demographics already exist
capture confirm var female
if _rc == 0 {
    di "WARNING: 'female' already exists in dataset. Dropping old version."
    drop female
}
capture confirm var married
if _rc == 0 {
    di "WARNING: 'married' already exists in dataset. Dropping old version."
    drop married
}
capture confirm var has_children
if _rc == 0 {
    di "WARNING: 'has_children' already exists in dataset. Dropping old version."
    drop has_children
}
capture confirm var pesex
if _rc == 0 {
    drop pesex
}
capture confirm var pemaritl
if _rc == 0 {
    drop pemaritl
}
capture confirm var prnmchld
if _rc == 0 {
    drop prnmchld
}

* Merge
merge 1:1 hrhhid hrhhid2 pulineno year ///
    using "$datadir/demographics_all_years.dta", ///
    keep(master match) nogen ///
    keepusing(pesex pemaritl prnmchld female married has_children)

* Check merge quality
di _n "Post-merge summary:"
di "  Total N = " _N
quietly count if female != .
di "  Non-missing female: " r(N) " (" %5.1f 100*r(N)/_N "%)"
quietly count if married != .
di "  Non-missing married: " r(N) " (" %5.1f 100*r(N)/_N "%)"
quietly count if has_children != .
di "  Non-missing has_children: " r(N) " (" %5.1f 100*r(N)/_N "%)"

* Summary statistics
di _n "Demographic Summary (analysis sample, year >= 2013):"
summarize female married has_children if year >= 2013 & age >= 18 & age <= 64

tab pesex if year >= 2013 & age >= 18 & age <= 64, m
tab pemaritl if year >= 2013 & age >= 18 & age <= 64, m
summarize prnmchld if year >= 2013 & age >= 18 & age <= 64

* Self-employment by gender
di _n "Self-employment by gender:"
tab female self_employed [aw=weight] if year >= 2013 & age >= 18 & age <= 64, row

* Self-employment by marital status
di _n "Self-employment by marital status:"
tab married self_employed [aw=weight] if year >= 2013 & age >= 18 & age <= 64, row

* Save updated dataset
save "$datadir/analysis_dataset_with_se.dta", replace

di _n "============================================================"
di "DEMOGRAPHICS MERGED SUCCESSFULLY"
di "============================================================"
di "Variables added: female, married, has_children, pesex, pemaritl, prnmchld"

* Clean up temporary files
foreach yr in 2023 2021 2019 2017 2015 2013 {
    erase "$datadir/demographics_`yr'.dta"
}
erase "$datadir/demographics_all_years.dta"

log close
