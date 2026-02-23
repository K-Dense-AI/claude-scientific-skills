# Cross-validation methods for genomic prediction (genomic selection): research notes

This note focuses on CV designs that are appropriate for genomic prediction when individuals are structured by ancestry/subpopulation and/or relatedness (families, half-sib/full-sib structure), plus forward (chronological) validation.

## 1) Types of cross-validation (CV) used in genomic prediction

### K-fold CV (random)
- What it does: random split into k folds; train on k-1, test on 1.
- Risk in genomics: optimistic accuracy if close relatives or strong population clusters are split across train/test.

### Stratified k-fold CV (by population structure)
- What it does: preserves the distribution of a stratification label (e.g., subpopulation/cluster) in each fold.
- Use when: you want each fold to contain representation from each subpopulation.
- Still not enough when: relatedness is the dominant driver of predictive performance (you can still leak family information across folds).

### Leave-one-out CV (LOOCV)
- What it does: each individual is its own test set once.
- Notes: expensive; and unless you also remove close relatives, LOOCV does not solve relatedness leakage.
- Genomic-specific note: efficient LOOCV computations exist for BLUP/GBLUP-type models.

### Forward validation (chronological)
- What it does: train on earlier years/generations, test on later years/generations.
- Why it matters: this often matches how genomic selection is deployed (predict future candidates).
- Pitfalls: year/environment confounding, genetic trend, and low relatedness between historical training and future candidates.

### Relatedness-aware CV (accounting for family structure)
Common variants:
- Leave-one-family-out / GroupKFold: keep all members of a family (or half-sib/full-sib group) together.
- Kinship-threshold grouping: cluster individuals using a GRM/kinship matrix; keep clusters intact across folds.

## 2) Key considerations (genomic data)

### Population stratification
- If train/test have different ancestry composition, estimated accuracy depends strongly on whether the model can transfer across groups.
- Recommended: report (i) within-subpopulation accuracy and (ii) across-subpopulation/held-out-subpopulation accuracy when portability matters.

### Family structure / relatedness
- If close relatives are present across train/test, models can “predict” primarily via shared haplotypes rather than portable marker effects.
- Recommended: split by family/pedigree IDs when available; otherwise use GRM-based grouping.

### Avoiding leakage
- Standardization/feature selection must be done inside each training fold.
- In breeding data, be careful with precomputed BLUPs: if BLUPs were estimated using all data (including validation environments/years), that can leak information.

## 3) Python and R implementations

Runnable scripts are provided:
- `genomic_cv_python.py`: uses scikit-learn splitters plus custom kinship-aware splitting; RR-BLUP via ridge on markers; GBLUP via kernel ridge on GRM.
- `genomic_cv_r.R`: uses rrBLUP `mixed.solve()` for RR-BLUP/GBLUP; caret grouping helpers; forward validation.

## 4) Selected references (with DOIs)

Population structure / CV scenarios
- Werner CR, Gaynor RC, Gorjanc G, et al. How Population Structure Impacts Genomic Selection Accuracy in Cross-Validation: Implications for Practical Breeding. Front Plant Sci. 2020. doi:10.3389/fpls.2020.592977
- Guo Z, Tucker DM, Basten CJ, et al. The impact of population structure on genomic prediction in stratified populations. Theor Appl Genet. 2014. doi:10.1007/s00122-013-2255-x

Pitfalls and leakage
- Runcie D, Cheng H. Pitfalls and Remedies for Cross Validation with Multi-trait Genomic Prediction Methods. G3 (Bethesda). 2019. doi:10.1534/g3.119.400598
- Perez-Cabal MA, Vazquez AI, Gianola D, et al. Accuracy of genome-enabled prediction in a dairy cattle population using different cross-validation layouts. Front Genet. 2012. doi:10.3389/fgene.2012.00027

Efficient LOOCV for BLUP/GBLUP
- Cheng H, Garrick DJ, Fernando RL. Efficient strategies for leave-one-out cross validation for genomic best linear unbiased prediction. J Anim Sci Biotechnol. 2017. doi:10.1186/s40104-017-0164-6

Forward/temporal validation in breeding
- Bernal-Vasquez AM, Gordillo A, Schmidt M, Piepho HP. Genomic prediction in early selection stages using multi-year data in a hybrid rye breeding program. BMC Genet. 2017. doi:10.1186/s12863-017-0512-8
- Lenz PRN, Beaulieu J, Mansfield SD, et al. Genomic prediction for hastening and improving efficiency of forward selection in conifer polycross mating designs. Heredity. 2020. doi:10.1038/s41437-019-0290-3
- Schrauf MF, et al. Comparing Genomic Prediction Models by Means of Cross Validation. Front Plant Sci. 2021. doi:10.3389/fpls.2021.734512

