# Cross-Validation for Genomic Prediction (Genomic Selection)

Runnable, minimal examples of cross-validation (CV) designs for genomic prediction with attention to:
- population structure
- family structure / relatedness
- forward (chronological) validation

## Files

- `genomic_cv_python.py`  Python: RR-BLUP (ridge on markers) + GBLUP (kernel ridge on GRM) + multiple CV splitters
- `genomic_cv_r.R`        R: rrBLUP `mixed.solve()` + caret-based grouping helpers + multiple CV splitters
- `genomic_cv_research_summary.md`  Notes + citations + design guidance
- `GENOMIC_CV_GUIDE.md`   Quick reference

## Run

Python:
```bash
pip install numpy pandas scikit-learn
python genomic_cv_python.py
```

R:
```r
install.packages(c("rrBLUP", "caret"))
source("genomic_cv_r.R")
```

## What The Examples Cover

- K-fold CV: `KFold`
- Stratified CV by subpopulation: `StratifiedKFold`
- Leave-one-out CV: `LeaveOneOut` (included; not run by default)
- Forward validation: train on time <= t, test on time == t+1
- Relatedness-aware CV:
  - `GroupKFold` by known family IDs (best if you have pedigree/family labels)
  - kinship-threshold connected-components splitter (fallback if you only have a GRM)

## Key Takeaway

Random CV can substantially overestimate genomic prediction accuracy when close relatives or population clusters are present in both train and test sets. Use group-wise CV (family/subpopulation) or forward validation when that matches the deployment scenario.
