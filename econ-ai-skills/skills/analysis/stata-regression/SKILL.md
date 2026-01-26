---
name: stata-regression
description: Run regression analyses in Stata with publication-ready output tables.
workflow_stage: analysis
compatibility:
  - claude-code
  - cursor
  - codex
  - gemini-cli
author: Awesome Econ AI Community
version: 1.0.0
tags:
  - stata
  - regression
  - esttab
  - econometrics
---

# Stata Regression

## Purpose

This skill produces reproducible regression analysis workflows in Stata, including model diagnostics and publication-ready tables using `esttab` or `outreg2`.

## When to Use

- Estimating linear or nonlinear regression models in Stata
- Producing tables for academic papers and reports
- Running robustness checks and alternative specifications

## Instructions

Follow these steps to complete the task:

### Step 1: Understand the Context

Before generating any code, ask the user:

- What is the dependent variable and key regressors?
- What controls and fixed effects are required?
- How should standard errors be clustered?
- What output format is needed (LaTeX, Word, or CSV)?

### Step 2: Generate the Output

Based on the context, generate Stata code that:

1. **Loads and checks the data** - Handle missing values and verify variable types
2. **Runs the requested specification** - Use `regress`, `reghdfe`, or `xtreg` as appropriate
3. **Adds robust or clustered standard errors** - Match the study design
4. **Exports tables** - Use `esttab` or `outreg2` with clear labels

### Step 3: Verify and Explain

After generating output:

- Explain what each model estimates
- Highlight assumptions and diagnostics
- Suggest robustness checks or alternative models

## Example Prompts

- "Run OLS with firm and year fixed effects, clustering by firm"
- "Estimate a logit model and export results to LaTeX"
- "Create a regression table with three specifications"

## Example Output

```stata
* ============================================
* Regression Analysis with Stata
* ============================================

* Load data
use "data.dta", clear

* Summary stats
summarize y x1 x2 x3

* Main regression with clustered SEs
regress y x1 x2 x3, vce(cluster firm_id)
eststo model1

* Alternative specification with fixed effects
reghdfe y x1 x2 x3, absorb(firm_id year) vce(cluster firm_id)
eststo model2

* Export table
esttab model1 model2 using "results/regression_table.tex", replace se label
```

## Requirements

### Software

- Stata 17+

### Packages

- `estout` (for `esttab`)
- `reghdfe` (optional, for high-dimensional fixed effects)

Install with:

```stata
ssc install estout
ssc install reghdfe
```

## Best Practices

1. **Match standard errors to the design** (cluster where treatment varies)
2. **Report all model variants** used in the analysis
3. **Document variable definitions** and transformations

## Common Pitfalls

- Not clustering standard errors at the correct level
- Omitting fixed effects when required by the design
- Exporting tables without clear labels and notes

## References

- [Stata Regression Reference Manual](https://www.stata.com/manuals/rregress.pdf)
- [reghdfe documentation](https://github.com/sergiocorreia/reghdfe)
- [estout documentation](https://repec.sowi.unibe.ch/stata/estout/)

## Changelog

### v1.0.0

- Initial release
