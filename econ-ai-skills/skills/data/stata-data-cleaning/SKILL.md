---
name: stata-data-cleaning
description: Clean and transform messy data in Stata with reproducible workflows
workflow_stage: data
compatibility:
  - claude-code
  - cursor
  - codex
  - gemini-cli
author: Awesome Econ AI Community
version: 1.0.0
tags:
  - Stata
  - data-cleaning
  - data-wrangling
  - reproducibility
---

# Stata Data Cleaning

## Purpose

This skill helps economists clean, transform, and prepare datasets for analysis in Stata. It emphasizes reproducibility, proper documentation, and handling common data quality issues found in economic research.

## When to Use

- Cleaning raw survey or administrative data
- Merging multiple data sources
- Handling missing values, duplicates, and outliers
- Creating analysis-ready panel datasets
- Documenting data transformations for replication

## Instructions

### Step 1: Understand the Data

Before generating code, ask the user:
1. What is the data source? (survey, administrative, API, etc.)
2. What is the unit of observation?
3. What are the key variables needed for analysis?
4. Are there known data quality issues to address?

### Step 2: Generate Cleaning Pipeline

Create a Stata do-file that:

1. **Has a clear header** with project info and date
2. **Sets up the environment** (clear all, set memory, log)
3. **Loads and inspects raw data**
4. **Documents each transformation** with comments
5. **Creates a codebook** for the final dataset

### Step 3: Follow Best Practices

- Use `assert` statements to verify data integrity
- Create labeled variables with `label variable`
- Use value labels for categorical variables
- Generate a log file for reproducibility
- Save intermediate files when appropriate

## Example Output

```stata
/*==============================================================================
    Project:    Economic Analysis Data Cleaning
    Author:     [Your Name]
    Date:       [Date]
    Purpose:    Clean raw survey data for regression analysis
    Input:      raw_survey_data.dta
    Output:     cleaned_analysis_data.dta
==============================================================================*/

* ============================================
* 1. SETUP
* ============================================

clear all
set more off
cap log close
log using "logs/data_cleaning_`c(current_date)'.log", replace

* Set working directory
cd "/path/to/project"

* Define globals for paths
global raw_data "data/raw"
global clean_data "data/clean"
global output "output"

* ============================================
* 2. LOAD AND INSPECT RAW DATA
* ============================================

use "${raw_data}/raw_survey_data.dta", clear

* Basic inspection
describe
summarize
codebook, compact

* Check for duplicates
duplicates report id_var
duplicates list id_var if _dup > 0

* ============================================
* 3. VARIABLE CLEANING
* ============================================

* --- Rename variables for clarity ---
rename q1 age
rename q2 income_reported
rename q3 education_level

* --- Clean numeric variables ---
* Replace missing value codes with .
mvdecode age income_reported, mv(-99 -88 -77)

* Cap outliers at 99th percentile
qui sum income_reported, detail
replace income_reported = r(p99) if income_reported > r(p99) & !mi(income_reported)

* --- Clean string variables ---
* Standardize state names
replace state = upper(trim(state))
replace state = "NEW YORK" if inlist(state, "NY", "N.Y.", "N Y")

* --- Create categorical variables ---
gen education_cat = .
replace education_cat = 1 if education_level < 12
replace education_cat = 2 if education_level == 12
replace education_cat = 3 if education_level > 12 & education_level <= 16
replace education_cat = 4 if education_level > 16 & !mi(education_level)

label define edu_lbl 1 "Less than HS" 2 "High School" 3 "College" 4 "Graduate"
label values education_cat edu_lbl

* ============================================
* 4. HANDLE MISSING DATA
* ============================================

* Create missing indicator variables
gen mi_income = mi(income_reported)

* Document missingness
tab mi_income

* Count complete cases
egen complete_case = rownonmiss(age income_reported education_cat)
tab complete_case

* ============================================
* 5. CREATE DERIVED VARIABLES
* ============================================

* Age groups
gen age_group = .
replace age_group = 1 if age >= 18 & age < 30
replace age_group = 2 if age >= 30 & age < 50
replace age_group = 3 if age >= 50 & age < 65
replace age_group = 4 if age >= 65 & !mi(age)

label define age_lbl 1 "18-29" 2 "30-49" 3 "50-64" 4 "65+"
label values age_group age_lbl

* Log income
gen log_income = ln(income_reported + 1)

* ============================================
* 6. DATA VALIDATION
* ============================================

* Assert expected ranges
assert age >= 18 & age <= 120 if !mi(age)
assert income_reported >= 0 if !mi(income_reported)

* Check variable types
assert !mi(id_var)
isid id_var  // Verify unique identifier

* ============================================
* 7. LABEL VARIABLES
* ============================================

label variable age "Age in years"
label variable income_reported "Annual income (USD)"
label variable education_cat "Education category"
label variable log_income "Log of annual income"
label variable mi_income "Missing income indicator"

* ============================================
* 8. FINAL CHECKS AND SAVE
* ============================================

* Keep relevant variables
keep id_var age age_group income_reported log_income ///
     education_cat mi_income state year

* Order variables logically
order id_var year state age age_group income_reported ///
      log_income education_cat mi_income

* Compress to minimize file size
compress

* Save cleaned data
save "${clean_data}/cleaned_analysis_data.dta", replace

* Create codebook
codebook, compact

* Close log
log close

* ============================================
* END OF FILE
* ============================================
```

## Requirements

### Software
- Stata 15+ (some commands require newer versions)

### Recommended User-Written Commands
```stata
ssc install unique     // For unique value checking
ssc install mdesc      // For missing data patterns
ssc install labutil    // For label manipulation
```

## Best Practices

1. **Always start with `clear all`** to ensure clean environment
2. **Use log files** to document all transformations
3. **Comment extensively** - explain WHY, not just WHAT
4. **Use `assert` statements** to catch data errors early
5. **Create a data dictionary** alongside your cleaned data
6. **Version your do-files** and datasets

## Common Pitfalls

- ❌ Not checking for duplicates before merging
- ❌ Forgetting to handle missing value codes (-99, -88, etc.)
- ❌ Not labeling variables and values
- ❌ Overwriting raw data files
- ❌ Not documenting data transformations

## References

- [Stata Data Management Manual](https://www.stata.com/manuals/d.pdf)
- [Gentzkow & Shapiro (2014) Code and Data for the Social Sciences](https://web.stanford.edu/~gentzkow/research/CodeAndData.pdf)
- [DIME Analytics Data Handbook](https://worldbank.github.io/dime-data-handbook/)

## Changelog

### v1.0.0
- Initial release with comprehensive cleaning template
