################################################################################
# Mobile Banking, Bank Branch Closures, and Self-Employment in the United States
#
# Construct Merged Analysis Dataset (R version)
# Date: February 2025
################################################################################

# Required packages
library(tidyverse)
library(haven)      # For reading Stata files
library(jsonlite)   # For reading JSON

# Set paths
datadir <- "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/Mobile banking USA/Data"
fdic_dir <- file.path(datadir, "FDIC_Survey/hhmultiyear")
acs_dir <- file.path(datadir, "ACS")
sod_dir <- file.path(datadir, "SOD")
jmp_dir <- "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/JMP"

################################################################################
# PART 1: Load FDIC/CPS Multi-Year Survey Data
################################################################################

cat("Loading FDIC/CPS survey data...\n")

fdic_data <- read_csv(file.path(fdic_dir, "hh_multiyear_analys.csv"),
                      show_col_types = FALSE)

cat(paste("Loaded", nrow(fdic_data), "observations\n"))

# Rename key variables
fdic_data <- fdic_data %>%
  rename(
    cbsa = msa13,
    year = hryear4,
    statefips = gestfips,
    weight = hhsupwgt,
    age = prtage
  )

################################################################################
# PART 2: Create Key Analysis Variables
################################################################################

cat("Creating analysis variables...\n")

fdic_data <- fdic_data %>%
  mutate(
    # Mobile banking indicators
    mobile_primary = ifelse(hbnkaccm == 5 & !is.na(hbnkaccm), 1, 0),
    mobile_user = ifelse(hbnkaccm5 == 1 & hbnkaccm5 > 0 & !is.na(hbnkaccm5), 1, 0),

    # Off-site banking (mobile/online/phone only)
    offsite_only = case_when(
      hmobonltelstat == 1 & !is.na(hmobonltelstat) & hmobonltelstat > 0 ~ 1,
      hmobonltelstatv2 == 1 & !is.na(hmobonltelstatv2) & hmobonltelstatv2 > 0 ~ 1,
      TRUE ~ 0
    ),

    # Branch user
    branch_user = case_when(
      hbranchstat %in% c(1, 2) & !is.na(hbranchstat) & hbranchstat > 0 ~ 1,
      hbranchstatv2 %in% c(1, 2) & !is.na(hbranchstatv2) & hbranchstatv2 > 0 ~ 1,
      TRUE ~ 0
    ),

    # Banked status (consolidate across versions)
    banked = case_when(
      hbankstatv6 %in% c(2, 3) & !is.na(hbankstatv6) ~ 1,
      hbankstatv5 %in% c(2, 3) & !is.na(hbankstatv5) ~ 1,
      hbankstatv4 %in% c(2, 3) & !is.na(hbankstatv4) ~ 1,
      hbankstatv3 %in% c(2, 3) & !is.na(hbankstatv3) ~ 1,
      hbankstatv2 %in% c(2, 3) & !is.na(hbankstatv2) ~ 1,
      hbankstat %in% c(2, 3) & !is.na(hbankstat) ~ 1,
      TRUE ~ 0
    ),
    unbanked = 1 - banked,

    # Employment status
    employed = ifelse(pempstat == 1 & !is.na(pempstat), 1, 0),
    unemployed = ifelse(pempstat == 2 & !is.na(pempstat), 1, 0),
    nilf = ifelse(pempstat == 3 & !is.na(pempstat), 1, 0),

    # Race/ethnicity (praceeth3)
    black = ifelse(praceeth3 == 1 & !is.na(praceeth3), 1, 0),
    hispanic = ifelse(praceeth3 == 2 & !is.na(praceeth3), 1, 0),
    asian = ifelse(praceeth3 == 3 & !is.na(praceeth3), 1, 0),
    white = ifelse(praceeth3 == 6 & !is.na(praceeth3), 1, 0),
    other_race = ifelse(praceeth3 %in% c(4, 5, 7) & !is.na(praceeth3), 1, 0),

    # Education
    no_hs = ifelse(peducgrp == 1 & !is.na(peducgrp), 1, 0),
    hs_diploma = ifelse(peducgrp == 2 & !is.na(peducgrp), 1, 0),
    some_college = ifelse(peducgrp == 3 & !is.na(peducgrp), 1, 0),
    college_degree = ifelse(peducgrp == 4 & !is.na(peducgrp), 1, 0),

    # Income
    inc_below15k = ifelse(hhincome == 1 & !is.na(hhincome), 1, 0),
    inc_15_30k = ifelse(hhincome == 2 & !is.na(hhincome), 1, 0),
    inc_30_50k = ifelse(hhincome == 3 & !is.na(hhincome), 1, 0),
    inc_50_75k = ifelse(hhincome == 4 & !is.na(hhincome), 1, 0),
    inc_above75k = ifelse(hhincome == 5 & !is.na(hhincome), 1, 0),

    # Metro status
    metro = ifelse(gtmetsta == 1 & gtmetsta != 3 & !is.na(gtmetsta), 1, 0),
    rural = ifelse(gtmetsta == 2 & gtmetsta != 3 & !is.na(gtmetsta), 1, 0),

    # Homeownership
    homeowner = ifelse(hhtenure == 1 & !is.na(hhtenure), 1, 0),

    # Banking mode categories for structural model
    banking_mode = case_when(
      unbanked == 1 ~ 1,                          # Unbanked
      banked == 1 & offsite_only == 1 ~ 2,        # Mobile/Online Only
      banked == 1 & branch_user == 1 ~ 3,         # Branch User
      TRUE ~ NA_real_
    )
  )

################################################################################
# PART 3: Load and Process ACS Data
################################################################################

cat("Loading ACS data...\n")

# Broadband data
acs_broadband <- read_csv(file.path(acs_dir, "broadband_cbsa_2023.csv"),
                          show_col_types = FALSE) %>%
  rename(
    cbsa = `metropolitan statistical area/micropolitan statistical area`,
    total_households = S2801_C01_001E,
    pct_broadband = S2801_C02_012E,
    pct_cellular = S2801_C02_014E
  ) %>%
  mutate(cbsa = as.numeric(cbsa)) %>%
  select(cbsa, NAME, total_households, pct_broadband, pct_cellular)

# Demographics data
acs_demographics <- read_csv(file.path(acs_dir, "demographics_cbsa_2023.csv"),
                             show_col_types = FALSE) %>%
  rename(
    cbsa = `metropolitan statistical area/micropolitan statistical area`,
    population = DP05_0001E,
    pct_white_alone = DP05_0037PE,
    pct_black_alone = DP05_0038PE,
    pct_asian_alone = DP05_0039PE,
    pct_two_more_races = DP05_0044PE,
    pct_hispanic = DP05_0071PE
  ) %>%
  mutate(cbsa = as.numeric(cbsa)) %>%
  select(cbsa, population, starts_with("pct_"))

# Income/unemployment data
acs_income <- read_csv(file.path(acs_dir, "income_unemployment_cbsa_2023.csv"),
                       show_col_types = FALSE) %>%
  rename(
    cbsa = `metropolitan statistical area/micropolitan statistical area`,
    median_hh_income = S1901_C01_012E,
    unemployment_rate = S2301_C04_001E
  ) %>%
  mutate(cbsa = as.numeric(cbsa)) %>%
  select(cbsa, median_hh_income, unemployment_rate)

# Merge ACS datasets
acs_controls <- acs_broadband %>%
  left_join(acs_demographics, by = "cbsa") %>%
  left_join(acs_income, by = "cbsa")

cat(paste("ACS data:", nrow(acs_controls), "CBSAs\n"))

################################################################################
# PART 4: Merge All Datasets
################################################################################

cat("Merging datasets...\n")

# Merge FDIC with ACS
analysis_data <- fdic_data %>%
  left_join(acs_controls, by = "cbsa")

cat(paste("Merged dataset:", nrow(analysis_data), "observations\n"))

################################################################################
# PART 5: Summary Statistics
################################################################################

cat("\n=== Summary Statistics ===\n\n")

# Sample by year
cat("Sample size by year:\n")
print(table(analysis_data$year))

# Mobile banking adoption over time
cat("\nMobile banking adoption by year:\n")
mobile_by_year <- analysis_data %>%
  filter(!is.na(mobile_user)) %>%
  group_by(year) %>%
  summarise(
    n = n(),
    mobile_user_pct = weighted.mean(mobile_user, weight, na.rm = TRUE) * 100,
    mobile_primary_pct = weighted.mean(mobile_primary, weight, na.rm = TRUE) * 100
  )
print(mobile_by_year)

# Banking status by year
cat("\nBanking status by year:\n")
banking_by_year <- analysis_data %>%
  group_by(year) %>%
  summarise(
    banked_pct = weighted.mean(banked, weight, na.rm = TRUE) * 100,
    unbanked_pct = weighted.mean(unbanked, weight, na.rm = TRUE) * 100
  )
print(banking_by_year)

# Mobile banking by demographics
cat("\nMobile banking by race/ethnicity:\n")
mobile_by_race <- analysis_data %>%
  filter(!is.na(mobile_user)) %>%
  mutate(race = case_when(
    black == 1 ~ "Black",
    hispanic == 1 ~ "Hispanic",
    asian == 1 ~ "Asian",
    white == 1 ~ "White",
    TRUE ~ "Other"
  )) %>%
  group_by(race) %>%
  summarise(
    n = n(),
    mobile_user_pct = weighted.mean(mobile_user, weight, na.rm = TRUE) * 100
  )
print(mobile_by_race)

################################################################################
# PART 6: Save Analysis Dataset
################################################################################

cat("\nSaving analysis dataset...\n")

# Select key variables for analysis
analysis_final <- analysis_data %>%
  select(
    # Identifiers
    hrhhid, hrhhid2, pulineno, cbsa, year, weight,
    # Mobile banking
    mobile_primary, mobile_user, offsite_only, branch_user, banking_mode,
    # Banking status
    banked, unbanked,
    # Employment
    employed, unemployed, nilf, pempstat,
    # Demographics
    age, black, hispanic, asian, white, other_race,
    no_hs, hs_diploma, some_college, college_degree,
    inc_below15k, inc_15_30k, inc_30_50k, inc_50_75k, inc_above75k,
    metro, rural, homeowner,
    # ACS controls
    pct_broadband, pct_cellular, median_hh_income, unemployment_rate,
    population, pct_white_alone, pct_black_alone, pct_asian_alone, pct_hispanic,
    # Original survey variables for reference
    hbankstatv6, hbnkaccm, gtmetsta, praceeth3, peducgrp, hhincome, hhtenure
  )

# Save as RDS
saveRDS(analysis_final, file.path(datadir, "analysis_dataset.rds"))

# Save as CSV
write_csv(analysis_final, file.path(datadir, "analysis_dataset.csv"))

cat(paste("\nAnalysis dataset saved!\n"))
cat(paste("Total observations:", nrow(analysis_final), "\n"))
cat(paste("Years: 2009-2023 (biennial)\n"))
cat(paste("CBSAs:", length(unique(analysis_final$cbsa[!is.na(analysis_final$cbsa)])), "\n"))

################################################################################
# NOTE: Self-employment extraction requires parsing raw CPS fixed-width files
# See extract_self_employment.R for details
################################################################################
