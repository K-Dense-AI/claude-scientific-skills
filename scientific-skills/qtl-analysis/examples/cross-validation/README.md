# Cross-Validation for Genomic Prediction

## What This Example Does

This example demonstrates **why cross-validation strategy matters** in genomic prediction. It shows how different CV approaches can give misleadingly high accuracy due to information leakage.

## The Information Leakage Problem

When close relatives appear in both training and test sets, prediction is artificially easy — the model just "memorizes" family patterns. This leads to:
- Overoptimistic accuracy estimates
- Poor generalization to new breeding lines

## CV Strategies Compared

| Strategy | Description | Use When |
|----------|-------------|----------|
| **Standard K-fold** | Random splits | Baseline (often overly optimistic) |
| **Stratified K-fold** | Balanced subpopulations | Population structure present |
| **GroupKFold (family)** | Keep families together | Known family/pedigree structure |
| **Forward validation** | Train on past, test on future | Temporal/breeding program context |
| **GBLUP** | Kernel on GRM | Using relationship matrix |

## Running the Example

```bash
cd scientific-skills/qtl-analysis/examples/cross-validation
python run_cv.py
```

## Input → Process → Output

### Input
- 400 individuals, 800 markers
- 3 subpopulations, 60 families
- 6 time points (years)
- Heritability: 0.6

### Process
- Run 5 strategies: K-fold, Stratified, GroupKFold, Forward, GBLUP
- Each with 5-fold CV

### Output
- `output/cv_comparison.png` — Bar chart comparing strategies
- `output/cv_results.csv` — Per-fold results
- `output/cv_summary.csv` — Mean ± SD per strategy

## Key Insight

GroupKFold and Forward validation give more realistic accuracy estimates. If Standard K-fold shows r=0.8 but GroupKFold shows r=0.4, your model is mostly exploiting family relatedness.
