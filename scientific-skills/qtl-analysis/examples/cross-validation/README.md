# Guide: Cross-Validation Methods for Genomic Prediction

This is a quick reference. The runnable code is in:
- `genomic_cv_python.py`
- `genomic_cv_r.R`

## 1) CV types and when to use them

- K-fold CV (random): baseline; can be optimistic if close relatives are split across folds.
- Stratified k-fold CV (by subpopulation/cluster): preserves subpopulation proportions per fold; still optimistic if relatedness crosses folds.
- Leave-one-out CV (LOOCV): useful for small n; expensive and still leaks if relatives remain in train.
- Forward validation (chronological): train on earlier years/generations, test on later; best match for real breeding roll-forward.
- Relatedness-aware CV:
  - Leave-one-family-out / GroupKFold: best if you have family/pedigree IDs.
  - Kinship-threshold grouping (connected components on GRM): practical fallback when you only have a GRM.

## 2) Python: minimal runnable snippets

### Simulate data
```python
from genomic_cv_python import simulate_genomic_data

data = simulate_genomic_data(
    n_individuals=300,
    n_markers=800,
    n_qtl=40,
    n_subpopulations=3,
    n_families=60,
    n_years=6,
    heritability=0.6,
    seed=42,
)
X = data["genotypes"]
y = data["phenotypes"]
subp = data["subpopulations"]
fam = data["family_ids"]
time = data["time"]
G = data["kinship"]
```

### K-fold and stratified k-fold (scikit-learn)
```python
from genomic_cv_python import kfold_cv_rrblup, stratified_kfold_cv_rrblup

print(kfold_cv_rrblup(X, y, n_splits=5, random_state=42))
print(stratified_kfold_cv_rrblup(X, y, strata=subp, n_splits=5, random_state=42))
```

### Relatedness-aware CV
If you have family IDs:
```python
from genomic_cv_python import group_kfold_cv_rrblup

print(group_kfold_cv_rrblup(X, y, groups=fam, n_splits=5))
```

If you only have a GRM/kinship matrix:
```python
from genomic_cv_python import relatedness_aware_cv_rrblup

print(relatedness_aware_cv_rrblup(X, y, kinship=G, threshold=0.125, n_splits=5, random_state=42))
```

### Forward (chronological) validation
```python
from genomic_cv_python import forward_validation_rrblup

print(forward_validation_rrblup(X, y, time=time, min_train_timepoints=2))
```

### GBLUP in a CV loop
```python
from genomic_cv_python import gblup_kfold_cv

print(gblup_kfold_cv(y, kinship=G, n_splits=5, random_state=42, lambda_param=1.0))
```

## 3) R: minimal runnable snippets

```r
source("genomic_cv_r.R")

data <- simulate_genomic_data(n_individuals=300, n_markers=800, n_qtl=40,
                              n_subpopulations=3, n_families=60, n_years=6,
                              heritability=0.6, seed=42)
X <- data$genotypes
y <- data$phenotypes
subp <- data$subpopulations
fam <- data$family_ids
time <- data$time
K <- data$kinship

standard_kfold_cv_rrblup(X, y, n_folds=5, seed=42)
stratified_kfold_cv_rrblup(X, y, subp, n_folds=5, seed=42)
family_group_kfold_cv_rrblup(X, y, fam, n_folds=5, seed=42)
relatedness_aware_cv_rrblup(X, y, K, threshold=0.125, n_folds=5, seed=42)
gblup_kfold_cv(y, K, n_folds=5, seed=42)
forward_validation_rrblup(X, y, time, min_train_timepoints=2)
```

## 4) Key considerations (genomics-specific)

- Population stratification:
  - stratify folds by inferred subpopulation (PCA/clusters) or known ancestry labels
  - consider reporting across-subpopulation holdout (leave-one-subpopulation-out) if portability matters
- Relatedness and leakage:
  - do not allow close relatives in both train and test when estimating deployment accuracy
  - prefer group-wise splits by family/pedigree (leave-one-family-out)
  - if only GRM available, group by kinship threshold / clustering on GRM
- Forward validation:
  - match the real breeding timeline (train on past, predict future)
  - watch out for confounding by year/environment effects; consider modeling them explicitly

## References (selected)

- Werner CR, Gaynor RC, Gorjanc G, et al. How Population Structure Impacts Genomic Selection Accuracy in Cross-Validation. Front Plant Sci. 2020. doi:10.3389/fpls.2020.592977
- Guo Z, Tucker DM, Basten CJ, et al. The impact of population structure on genomic prediction in stratified populations. Theor Appl Genet. 2014. doi:10.1007/s00122-013-2255-x
- Runcie D, Cheng H. Pitfalls and Remedies for Cross Validation with Multi-trait Genomic Prediction Methods. G3 (Bethesda). 2019. doi:10.1534/g3.119.400598
- Cheng H, Garrick DJ, Fernando RL. Efficient strategies for leave-one-out cross validation for GBLUP. J Anim Sci Biotechnol. 2017. doi:10.1186/s40104-017-0164-6
