* Quick regression for methods paper OLS table
* Full sample, Branch User + Unbanked (mobile reference)

clear all
set more off

use "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/_Research/Mobile_Money_Banking/Mobile banking USA/Data/analysis_dataset_with_se.dta", clear

* Sample restrictions (same as phase0)
keep if age >= 18 & age <= 64
keep if employed == 1 | unemployed == 1
keep if year >= 2013

* Create SE indicator
gen se = (self_employed == 1) if self_employed != .

* Create banking mode dummies (mobile=reference)
gen branch = (banking_mode == 3) if banking_mode != .
cap drop unbanked
gen unbanked = (banking_mode == 1) if banking_mode != .

gen age2 = age * age

di "=== FULL SAMPLE REGRESSIONS FOR METHODS PAPER TABLE ==="
di "N observations with banking_mode:"
count if banking_mode != . & se != .

* Col 1: No controls
reg se branch unbanked [pw=hsupwgtk] if banking_mode != ., cluster(cbsa)
eststo c1

* Col 2: Demographics
reg se branch unbanked ///
    age age2 ///
    i.praceeth3 i.peducgrp i.hhincome ///
    female married has_children ///
    metro ///
    [pw=hsupwgtk] if banking_mode != ., cluster(cbsa)
eststo c2

* Col 3: Add CBSA and year FE
reghdfe se branch unbanked ///
    age age2 ///
    i.praceeth3 i.peducgrp i.hhincome ///
    female married has_children ///
    metro ///
    [pw=hsupwgtk] if banking_mode != ., absorb(cbsa year) cluster(cbsa)
eststo c3

* Col 4: Add CBSA-level controls
reghdfe se branch unbanked ///
    age age2 ///
    i.praceeth3 i.peducgrp i.hhincome ///
    female married has_children ///
    metro ///
    pct_broadband unemployment_rate ///
    [pw=hsupwgtk] if banking_mode != ., absorb(cbsa year) cluster(cbsa)
eststo c4

esttab c1 c2 c3 c4, ///
    keep(branch unbanked age) ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N r2 r2_a) ///
    title("Methods Paper OLS: SE on Banking Mode (Mobile Reference)")
