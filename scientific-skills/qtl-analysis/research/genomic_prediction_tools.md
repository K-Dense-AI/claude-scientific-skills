<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Open-Source Genomic Prediction Tools: Replacing QTLmax

Comprehensive guide to open-source alternatives for genomic prediction and BLUP (Best Linear Unbiased Prediction) analysis.

---

## Table of Contents

1. [Overview](#overview)
2. [rrBLUP (R)](#1-rrblup-r)
3. [BGLR - Bayesian Generalized Linear Regression (R)](#2-bglr---bayesian-generalized-linear-regression-r)
4. [glmnet - Elastic Net (R & Python)](#3-glmnet---elastic-net-r--python)
5. [sommer - Mixed Models (R)](#4-sommer---mixed-models-r)
6. [Python scikit-learn Alternatives](#5-python-scikit-learn-alternatives)
7. [Tool Comparison](#tool-comparison)
8. [SNP Data Formats](#snp-data-formats)

---

## Overview

### What is Genomic Prediction?

Genomic prediction uses genome-wide markers (SNPs) to predict breeding values or phenotypes. The two main approaches are:

- **RR-BLUP (Ridge Regression BLUP)**: Estimates marker effects using ridge regression
- **GBLUP (Genomic BLUP)**: Uses a genomic relationship matrix to estimate breeding values

### Why Replace QTLmax?

QTLmax is proprietary software for QTL mapping and genomic prediction. Open-source alternatives provide:
- **Cost**: Free and open-source
- **Flexibility**: Customizable code and methods
- **Community**: Active development and support
- **Integration**: Easy integration with R/Python workflows

---

## 1. rrBLUP (R)

**Package**: `rrBLUP`  
**Author**: Jeffrey Endelman  
**Citation**: Endelman (2011), The Plant Genome  
**Best for**: Ridge regression, GBLUP, kinship-based prediction

### Installation

```r
# Install from CRAN
install.packages("rrBLUP")

# Or install from GitHub for development version
# install.packages("devtools")
# devtools::install_github("jendelman/rrBLUP")
```

### Core Functions

| Function | Purpose |
|----------|---------|
| `mixed.solve()` | General mixed model solver (REML/ML) |
| `kinship.BLUP()` | Genomic prediction using kinship matrix |
| `A.mat()` | Calculate additive relationship matrix |
| `GWAS()` | Genome-wide association analysis |

### Minimal Working Example

```r
library(rrBLUP)

# Load example data (wheat dataset)
data(wheat)

# wheat.Y: Phenotypes (599 lines × 4 traits)
# wheat.X: SNP markers (599 lines × 1279 markers)

# Example 1: RR-BLUP (Ridge Regression) - estimate marker effects
# Model: y = Xβ + ε, where β ~ N(0, σ²_u I)

y <- wheat.Y[, 1]  # First trait
X <- wheat.X       # SNP markers (coded as -1, 0, 1)

# Fit ridge regression
rr_result <- mixed.solve(y = y, Z = X, method = "REML")

# Results
print(rr_result$Vu)      # Marker variance component
print(rr_result$Ve)      # Residual variance
print(head(rr_result$u))  # Marker effects (BLUPs)

# Predict breeding values
gebv <- X %*% rr_result$u
print(cor(y, gebv))  # Training accuracy
```

### GBLUP Example (Using Kinship Matrix)

```r
library(rrBLUP)

# Calculate additive relationship matrix
A <- A.mat(wheat.X, min.MAF = 0.05)

# GBLUP: Use relationship matrix instead of markers
# Model: y = μ + u + ε, where u ~ N(0, A σ²_g)

gblup_result <- kinship.BLUP(
  y = wheat.Y[, 1],
  K = A,
  method = "REML"
)

# Results
print(gblup_result$g)      # Genomic breeding values (GEBVs)
print(gblup_result$Vg)   # Genetic variance
print(gblup_result$Ve)   # Residual variance
```

### Cross-Validation Example

```r
library(rrBLUP)

# 5-fold cross-validation
set.seed(123)
n <- nrow(wheat.X)
folds <- sample(rep(1:5, length.out = n))

accuracy <- numeric(5)

for (i in 1:5) {
  # Training and validation sets
  train_idx <- which(folds != i)
  valid_idx <- which(folds == i)
  
  y_train <- wheat.Y[train_idx, 1]
  y_valid <- wheat.Y[valid_idx, 1]
  X_train <- wheat.X[train_idx, ]
  X_valid <- wheat.X[valid_idx, ]
  
  # Fit model on training data
  fit <- mixed.solve(y = y_train, Z = X_train)
  
  # Predict validation set
  pred <- X_valid %*% fit$u
  
  # Accuracy
  accuracy[i] <- cor(y_valid, pred)
}

cat("Cross-validation accuracy:", mean(accuracy), "\n")
cat("Standard deviation:", sd(accuracy), "\n")
```

### Input/Output Formats

**Input**:
- Phenotypes: Numeric vector or matrix (n × t traits)
- Genotypes: Numeric matrix (n × m markers), typically coded as -1, 0, 1
- Covariates: Optional design matrix for fixed effects

**Output**:
- `u`: Marker effects (for RR-BLUP) or breeding values (for GBLUP)
- `Vu`/`Vg`: Variance component for markers/genetic
- `Ve`: Residual variance
- `beta`: Fixed effects

---

## 2. BGLR - Bayesian Generalized Linear Regression (R)

**Package**: `BGLR`  
**Authors**: Gustavo de los Campos, Paulino Perez Rodriguez  
**Citation**: Perez & de los Campos (2014), Genetics  
**Best for**: Bayesian methods, multiple shrinkage priors, complex models

### Installation

```r
# Install from CRAN
install.packages("BGLR")

# Or from GitHub for latest version
# devtools::install_github("gdlc/BGLR-R")
```

### Key Features

| Model | Prior | Use Case |
|-------|-------|----------|
| BRR | Gaussian | Ridge regression (equivalent to RR-BLUP) |
| BayesA | Scaled-t | Heavy-tailed marker effects |
| BayesB | Mixture | Variable selection (some markers have zero effect) |
| BayesC | Mixture | Similar to BayesB but common variance |
| BL | Double exponential | Bayesian LASSO |

### Minimal Working Example

```r
library(BGLR)

# Load example data (mice dataset)
data(mice)

# mice.X: SNP markers (1814 lines × 10,346 markers)
# mice.pheno: Phenotypes with multiple traits

# Prepare data
y <- mice.pheno$Obesity.BMI  # Continuous trait
X <- mice.X                    # SNP markers

# Center and scale markers (important!)
X <- scale(X, center = TRUE, scale = TRUE)

# Example 1: Bayesian Ridge Regression (BRR)
# Equivalent to RR-BLUP but with Bayesian inference

fm_brr <- BGLR(
  y = y,
  ETA = list(
    markers = list(X = X, model = 'BRR')
  ),
  nIter = 12000,      # Total iterations
  burnIn = 2000,      # Burn-in period
  thin = 5,           # Thinning
  verbose = FALSE
)

# Results
print(fm_brr$varE)           # Residual variance
print(fm_brr$ETA$markers$varB)  # Marker variance
print(head(fm_brr$ETA$markers$b))  # Marker effects
print(head(fm_brr$yHat))     # Predicted values

# Prediction accuracy
cat("Training correlation:", cor(y, fm_brr$yHat), "\n")
```

### BayesA Example (Heavy-tailed priors)

```r
library(BGLR)

# BayesA: Each marker has its own variance (scaled-t prior)
fm_bayesA <- BGLR(
  y = y,
  ETA = list(
    markers = list(X = X, model = 'BayesA')
  ),
  nIter = 12000,
  burnIn = 2000,
  thin = 5,
  verbose = FALSE
)

print(fm_bayesA$ETA$markers$varB[1:5])  # Marker-specific variances
```

### BayesB Example (Variable Selection)

```r
library(BGLR)

# BayesB: Some markers have zero effect (spike-slab prior)
# pi = probability of marker having zero effect

fm_bayesB <- BGLR(
  y = y,
  ETA = list(
    markers = list(X = X, model = 'BayesB', probIn = 0.05)
  ),
  nIter = 12000,
  burnIn = 2000,
  thin = 5,
  verbose = FALSE
)

# Check which markers are selected (non-zero effect)
marker_effects <- fm_bayesB$ETA$markers$b
selected <- which(marker_effects != 0)
cat("Number of selected markers:", length(selected), "\n")
```

### Multi-Environment Model

```r
library(BGLR)

# Model with both markers and pedigree
# Load wheat data from BGLR
data(wheat)

y <- wheat.Y[, 1]  # Trait 1
X <- wheat.X       # Markers

# Create a simple pedigree (for demonstration)
# In practice, use actual pedigree
n <- length(y)
pedigree <- data.frame(
  id = 1:n,
  sire = sample(c(0, 1:(n/2)), n, replace = TRUE),
  dam = sample(c(0, 1:(n/2)), n, replace = TRUE)
)

# Fit model with both markers and pedigree
fm_multi <- BGLR(
  y = y,
  ETA = list(
    markers = list(X = X, model = 'BRR'),
    pedigree = list(pedigree = pedigree, model = 'BRR')
  ),
  nIter = 12000,
  burnIn = 2000,
  verbose = FALSE
)

# Variance components
print(fm_multi$ETA$markers$varB)    # Marker variance
print(fm_multi$ETA$pedigree$varU)   # Pedigree variance
```

### Input/Output Formats

**Input**:
- `y`: Phenotype vector (numeric, can have missing values as NA)
- `ETA`: List of random effects (markers, pedigree, etc.)
- `X`: Marker matrix (n × m, numeric)
- `pedigree`: Data frame with id, sire, dam columns

**Output**:
- `yHat`: Predicted values
- `ETA[[i]]$b`: Marker effects for component i
- `varE`: Residual variance
- `ETA[[i]]$varB`: Marker variance for component i

---

## 3. glmnet - Elastic Net (R & Python)

**Package**: `glmnet` (R), `glmnet-python` (Python)  
**Authors**: Friedman, Hastie, Tibshirani, et al.  
**Citation**: Friedman et al. (2010), JSS  
**Best for**: LASSO, Ridge, Elastic Net regularization

### Installation

```r
# R
install.packages("glmnet")
```

```bash
# Python
pip install glmnet-py
# or
pip install scikit-learn  # For sklearn implementation
```

### Minimal Working Example (R)

```r
library(glmnet)

# Load example data
data(QuickStartExample)

# QuickStartExample contains:
# x: 100 x 20 matrix of predictors
# y: 100 vector of responses

x <- QuickStartExample$x
y <- QuickStartExample$y

# Example 1: Ridge Regression (alpha = 0)
# Equivalent to RR-BLUP
ridge_fit <- glmnet(x, y, alpha = 0)

# Cross-validation to find optimal lambda
cv_ridge <- cv.glmnet(x, y, alpha = 0, nfolds = 10)

# Plot CV results
plot(cv_ridge)

# Best lambda
best_lambda <- cv_ridge$lambda.min
cat("Optimal lambda:", best_lambda, "\n")

# Coefficients at best lambda
coef_ridge <- coef(cv_ridge, s = "lambda.min")
print(head(coef_ridge))

# Predictions
pred_ridge <- predict(cv_ridge, newx = x, s = "lambda.min")
cat("Training correlation:", cor(y, pred_ridge), "\n")
```

### Elastic Net Example (R)

```r
library(glmnet)

# Elastic Net: alpha = 0.5 (balance between ridge and lasso)
# alpha = 0: Ridge, alpha = 1: LASSO

enet_fit <- glmnet(x, y, alpha = 0.5)
cv_enet <- cv.glmnet(x, y, alpha = 0.5, nfolds = 10)

# Coefficients
coef_enet <- coef(cv_enet, s = "lambda.min")

# Count non-zero coefficients (selected markers)
non_zero <- sum(coef_enet != 0) - 1  # Exclude intercept
cat("Selected markers:", non_zero, "\n")
```

### Genomic Prediction with glmnet (R)

```r
library(glmnet)

# Simulate SNP data for demonstration
set.seed(123)
n <- 200  # Individuals
m <- 1000 # Markers

# Simulate SNP matrix (0, 1, 2 coding)
X <- matrix(rbinom(n * m, 2, 0.3), nrow = n, ncol = m)

# Simulate QTL (10 major markers)
qtl_idx <- sample(1:m, 10)
qtl_effects <- rnorm(10, 0, 1)
y <- X[, qtl_idx] %*% qtl_effects + rnorm(n, 0, 2)

# Split into training and validation
set.seed(456)
train_idx <- sample(1:n, round(0.8 * n))
valid_idx <- setdiff(1:n, train_idx)

X_train <- X[train_idx, ]
y_train <- y[train_idx]
X_valid <- X[valid_idx, ]
y_valid <- y[valid_idx]

# Fit Elastic Net
cv_fit <- cv.glmnet(X_train, y_train, alpha = 0.5, nfolds = 5)

# Predict validation set
pred_valid <- predict(cv_fit, newx = X_valid, s = "lambda.min")

# Accuracy
accuracy <- cor(y_valid, pred_valid)
cat("Prediction accuracy:", accuracy, "\n")

# Selected markers
coefs <- coef(cv_fit, s = "lambda.min")
selected_markers <- which(coefs[-1] != 0)  # Exclude intercept
cat("True QTL found:", sum(qtl_idx %in% selected_markers), "/", length(qtl_idx), "\n")
```

### Python Implementation

```python
import numpy as np
from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

# Simulate SNP data
np.random.seed(42)
n = 200  # Individuals
m = 1000  # Markers

# Simulate SNP matrix (0, 1, 2 coding)
X = np.random.binomial(2, 0.3, size=(n, m))

# Simulate QTL
qtl_idx = np.random.choice(m, 10, replace=False)
qtl_effects = np.random.normal(0, 1, 10)
y = X[:, qtl_idx] @ qtl_effects + np.random.normal(0, 2, n)

# Split data
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Scale features (important for regularization)
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Fit Elastic Net with cross-validation
enet_cv = ElasticNetCV(
    l1_ratio=0.5,  # alpha in R (0=ridge, 1=lasso)
    cv=5,
    max_iter=10000,
    random_state=42
)

enet_cv.fit(X_train_scaled, y_train)

# Predictions
y_pred = enet_cv.predict(X_test_scaled)

# Accuracy
accuracy = np.corrcoef(y_test, y_pred)[0, 1]
print(f"Prediction accuracy: {accuracy:.4f}")
print(f"Optimal alpha: {enet_cv.alpha_:.6f}")
print(f"Selected features: {np.sum(enet_cv.coef_ != 0)}")

# Check if true QTL were selected
selected = np.where(enet_cv.coef_ != 0)[0]
true_positives = len(set(selected) & set(qtl_idx))
print(f"True QTL found: {true_positives}/{len(qtl_idx)}")
```

### Input/Output Formats

**Input**:
- `x`: Matrix of predictors (n × p, numeric)
- `y`: Response vector (numeric)
- `alpha`: Elastic net mixing parameter (0=ridge, 1=lasso)
- `lambda`: Regularization parameter (optional, auto-selected by CV)

**Output**:
- `lambda.min`: Optimal lambda from cross-validation
- `lambda.1se`: Lambda within 1 SE of minimum
- `coef`: Coefficient matrix
- `predict`: Prediction function

---

## 4. sommer - Mixed Models (R)

**Package**: `sommer`  
**Author**: Giovanny Covarrubias-Pazaran  
**Citation**: Covarrubias-Pazaran (2016), PLoS ONE  
**Best for**: Complex mixed models, multivariate analysis, dominance/epistasis

### Installation

```r
# Install from CRAN
install.packages("sommer")

# Or from GitHub
# devtools::install_github("covaruber/sommer")
```

### Key Features

- Multi-trait models
- Dominance and epistatic effects
- Heterogeneous variance structures
- Spatial models
- Large dataset support (>10,000 observations)

### Minimal Working Example

```r
library(sommer)

# Load example data
data(CPdata)

# CPdata contains:
# CPgeno: SNP markers (genotypes)
# CPpheno: Phenotypes
# CPped: Pedigree

# Prepare data
Y <- CPpheno
Z <- CPgeno

# Calculate additive relationship matrix
A <- A.mat(Z)

# Example 1: Simple GBLUP
# Model: y = μ + u + ε, where Var(u) = A σ²_a

ans <- mmer(
  fixed = color ~ 1,           # Fixed effects
  random = ~ vsr(id, Gu = A),  # Random effects with covariance A
  rcov = ~ units,              # Residual covariance
  data = Y
)

# Results
summary(ans)

# Extract breeding values
gebv <- ans$U$`u:id`$color
cat("Number of GEBVs:", length(gebv), "\n")
```

### RR-BLUP with sommer

```r
library(sommer)

data(CPdata)

# RR-BLUP: Marker-based model
# Model: y = μ + Xb + ε, where b ~ N(0, I σ²_m)

# Use markers directly
ans_rr <- mmer(
  fixed = color ~ 1,
  random = ~ vsr(id, Gu = diag(nrow(CPgeno))),  # IID markers
  rcov = ~ units,
  data = CPpheno
)

summary(ans_rr)
```

### Multi-Trait Model

```r
library(sommer)

data(CPdata)

# Multi-trait GBLUP
A <- A.mat(CPgeno)

# Model multiple traits simultaneously
ans_mt <- mmer(
  fixed = cbind(color, Yield) ~ 1,
  random = ~ vsr(id, Gu = A),
  rcov = ~ units,
  data = CPpheno
)

# Results show genetic correlations
summary(ans_mt)

# Genetic covariance matrix
print(ans_mt$sigma$`u:id`)
```

### Cross-Validation with sommer

```r
library(sommer)

data(CPdata)
A <- A.mat(CPgeno)

# 5-fold CV
set.seed(123)
n <- nrow(CPpheno)
folds <- sample(rep(1:5, length.out = n))

accuracy <- numeric(5)

for (i in 1:5) {
  # Mask validation set
  Y_train <- CPpheno
  Y_train$color[folds == i] <- NA
  
  # Fit model
  ans_cv <- mmer(
    fixed = color ~ 1,
    random = ~ vsr(id, Gu = A),
    rcov = ~ units,
    data = Y_train
  )
  
  # Predictions
  pred <- ans_cv$U$`u:id`$color
  valid_idx <- which(folds == i)
  
  # Accuracy
  accuracy[i] <- cor(CPpheno$color[valid_idx], pred[valid_idx])
}

cat("CV accuracy:", mean(accuracy), "±", sd(accuracy), "\n")
```

### Input/Output Formats

**Input**:
- `fixed`: Formula for fixed effects
- `random`: Formula for random effects (can include covariance structures)
- `rcov`: Formula for residual covariance
- `data`: Data frame with phenotypes
- `Gu`: Relationship matrix (optional)

**Output**:
- `U`: BLUPs for random effects
- `sigma`: Variance components
- `Beta`: BLUEs for fixed effects
- `fitted`: Fitted values

---

## 5. Python scikit-learn Alternatives

For Python users, scikit-learn provides excellent alternatives for genomic prediction.

### Ridge Regression (Python)

```python
import numpy as np
from sklearn.linear_model import Ridge, RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

# Simulate data
np.random.seed(42)
n, m = 200, 1000
X = np.random.binomial(2, 0.3, (n, m))
qtl_idx = np.random.choice(m, 10, replace=False)
qtl_effects = np.random.normal(0, 1, 10)
y = X[:, qtl_idx] @ qtl_effects + np.random.normal(0, 2, n)

# Split data
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Scale features
scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)
X_test_s = scaler.transform(X_test)

# Ridge regression with CV
ridge_cv = RidgeCV(alphas=np.logspace(-3, 3, 100), cv=5)
ridge_cv.fit(X_train_s, y_train)

# Predictions
y_pred = ridge_cv.predict(X_test_s)
accuracy = np.corrcoef(y_test, y_pred)[0, 1]

print(f"Best alpha: {ridge_cv.alpha_:.6f}")
print(f"Prediction accuracy: {accuracy:.4f}")
```

### Random Forest for Genomic Prediction

```python
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score

# Using same simulated data as above

# Random Forest
rf = RandomForestRegressor(
    n_estimators=500,
    max_features='sqrt',
    min_samples_leaf=5,
    random_state=42,
    n_jobs=-1
)

# Cross-validation
scores = cross_val_score(rf, X, y, cv=5, scoring='r2')
print(f"RF CV R²: {scores.mean():.4f} ± {scores.std():.4f}")

# Feature importance
rf.fit(X, y)
importance = rf.feature_importances_
top_markers = np.argsort(importance)[-10:]
print(f"Top 10 markers: {top_markers}")
```

### Support Vector Regression (SVR)

```python
from sklearn.svm import SVR
from sklearn.model_selection import GridSearchCV

# SVR with RBF kernel
param_grid = {
    'C': [0.1, 1, 10, 100],
    'gamma': ['scale', 'auto', 0.001, 0.01, 0.1],
    'epsilon': [0.01, 0.1, 0.5]
}

svr = SVR(kernel='rbf')
grid_search = GridSearchCV(svr, param_grid, cv=5, n_jobs=-1)
grid_search.fit(X_train_s, y_train)

# Best model
best_svr = grid_search.best_estimator_
y_pred = best_svr.predict(X_test_s)
accuracy = np.corrcoef(y_test, y_pred)[0, 1]

print(f"Best parameters: {grid_search.best_params_}")
print(f"SVR accuracy: {accuracy:.4f}")
```

---

## Tool Comparison

| Feature | rrBLUP | BGLR | glmnet | sommer | Python sklearn |
|---------|--------|------|--------|--------|----------------|
| **Language** | R | R | R/Python | R | Python |
| **Method** | REML/ML | Bayesian | Penalized ML | REML/ML | Various |
| **RR-BLUP** | ✅ | ✅ | ✅ (α=0) | ✅ | ✅ |
| **LASSO** | ❌ | ✅ (BayesL) | ✅ (α=1) | ❌ | ✅ |
| **Elastic Net** | ❌ | ❌ | ✅ | ❌ | ✅ |
| **BayesA/B/C** | ❌ | ✅ | ❌ | ❌ | ❌ |
| **Multi-trait** | ❌ | ✅ | ❌ | ✅ | ❌ |
| **Dominance** | ❌ | ✅ | ❌ | ✅ | ❌ |
| **Speed** | Fast | Slow | Fast | Medium | Fast |
| **Large data** | Medium | Poor | Good | Good | Good |
| **Ease of use** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |

### When to Use Each Tool

| Scenario | Recommended Tool |
|----------|-----------------|
| Quick GBLUP/RR-BLUP | **rrBLUP** |
| Bayesian methods, complex priors | **BGLR** |
| Variable selection (LASSO) | **glmnet** |
| Multi-trait, dominance, epistasis | **sommer** |
| Python workflow, ML methods | **sklearn** |
| Large datasets (>10K samples) | **rrBLUP** or **glmnet** |

---

## SNP Data Formats

### Common Input Formats

#### 1. Numeric Matrix (Recommended)

```
# Markers coded as -1, 0, 1 or 0, 1, 2
# Rows = individuals, Columns = markers

       M1  M2  M3  M4  ...
Ind1    0   1   0   2  ...
Ind2    1   0   1   1  ...
Ind3    0   0   2   0  ...
```

**R code to convert**:
```r
# From PLINK .raw format
library(rrBLUP)
geno <- read.table("data.raw", header = TRUE, row.names = 1)
geno <- as.matrix(geno[, -(1:5)])  # Remove first 5 columns (FID, IID, etc.)

# Recode if needed (PLINK uses 0/1/2, rrBLUP prefers -1/0/1)
geno <- geno - 1
```

#### 2. PLINK Format

```
# .bed (binary), .bim (marker info), .fam (sample info)
# Use BGLR's read_ped() or rrBLUP's read.table()

# In BGLR:
library(BGLR)
ped <- read_ped("data")
X <- ped$X
```

#### 3. VCF Format

```r
# Use vcfR package
library(vcfR)
vcf <- read.vcfR("data.vcf.gz")
gt <- extract.gt(vcf)

# Convert to numeric
geno <- apply(gt, 2, function(x) {
  sapply(strsplit(x, "[/|]"), function(y) sum(as.numeric(y)))
})
```

#### 4. HapMap Format

```r
# Read HapMap file
hapmap <- read.table("data.hmp.txt", header = TRUE, sep = "\t")

# Extract genotype columns (skip metadata columns)
geno <- hapmap[, -(1:11)]  # Skip first 11 metadata columns
geno <- t(as.matrix(geno))  # Transpose: samples as rows
```

### Data Preprocessing

```r
# 1. Quality Control - Remove low MAF markers
maf <- colMeans(geno) / 2
# For 0/1/2 coding: MAF = mean/2
# For -1/0/1 coding: MAF = (mean + 1)/2

keep <- maf > 0.05 & maf < 0.95
geno_filtered <- geno[, keep]

# 2. Impute missing values (if any)
# Mean imputation
geno_imputed <- apply(geno_filtered, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})

# 3. Center and scale (for some methods)
geno_scaled <- scale(geno_imputed, center = TRUE, scale = TRUE)

# 4. Calculate relationship matrix (for GBLUP)
library(rrBLUP)
A <- A.mat(geno_filtered, min.MAF = 0.05)
```

### Output Formats

**Marker Effects**:
```r
# Save marker effects
effects <- data.frame(
  Marker = colnames(geno),
  Effect = rr_result$u
)
write.csv(effects, "marker_effects.csv", row.names = FALSE)
```

**Genomic Breeding Values**:
```r
# Save GEBVs
gebv <- data.frame(
  ID = rownames(geno),
  GEBV = gblup_result$g
)
write.csv(gebv, "genomic_breeding_values.csv", row.names = FALSE)
```

**Variance Components**:
```r
# Save variance components
var_comp <- data.frame(
  Component = c("Genetic", "Residual", "Heritability"),
  Value = c(
    gblup_result$Vg,
    gblup_result$Ve,
    gblup_result$Vg / (gblup_result$Vg + gblup_result$Ve)
  )
)
write.csv(var_comp, "variance_components.csv", row.names = FALSE)
```

---

## Summary

### Quick Start Recommendations

1. **For beginners**: Start with **rrBLUP** - simplest interface, well-documented
2. **For variable selection**: Use **glmnet** - LASSO/Elastic Net for marker selection
3. **For complex models**: Use **BGLR** - Bayesian methods with various priors
4. **For multi-trait**: Use **sommer** - handles multiple traits and complex covariances
5. **For Python users**: Use **scikit-learn** - Ridge, ElasticNet, Random Forest

### Complete Workflow Example

```r
# Complete genomic prediction workflow

# 1. Load packages
library(rrBLUP)
library(glmnet)

# 2. Load data
data(wheat)

# 3. Prepare data
y <- wheat.Y[, 1]
X <- wheat.X

# 4. Quality control
maf <- colMeans(X + 1) / 2  # Convert -1/0/1 to 0/1/2
keep <- maf > 0.05 & maf < 0.95
X <- X[, keep]

# 5. Cross-validation
set.seed(123)
n <- length(y)
folds <- sample(rep(1:5, length.out = n))

results <- data.frame(
  rrBLUP = numeric(5),
  glmnet = numeric(5)
)

for (i in 1:5) {
  train <- folds != i
  valid <- folds == i
  
  # rrBLUP
  fit_rr <- mixed.solve(y[train], Z = X[train, ])
  pred_rr <- X[valid, ] %*% fit_rr$u
  results$rrBLUP[i] <- cor(y[valid], pred_rr)
  
  # glmnet
  fit_glm <- cv.glmnet(X[train, ], y[train], alpha = 0.5)
  pred_glm <- predict(fit_glm, newx = X[valid, ], s = "lambda.min")
  results$glmnet[i] <- cor(y[valid], pred_glm)
}

# 6. Compare results
print(colMeans(results))
print(apply(results, 2, sd))
```

---

## References

1. **rrBLUP**: Endelman, J.B. (2011). Ridge regression and other kernels for genomic selection. *The Plant Genome*, 4(3), 250-255.

2. **BGLR**: Pérez, P., & de los Campos, G. (2014). Genome-wide regression and prediction with the BGLR statistical package. *Genetics*, 198(2), 483-495.

3. **glmnet**: Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. *Journal of Statistical Software*, 33(1), 1.

4. **sommer**: Covarrubias-Pazaran, G. (2016). Genome-assisted prediction of quantitative traits using the R package sommer. *PLoS ONE*, 11(6), e0156744.

5. **Review**: de Vlaming, R., & Groenen, P.J. (2014). The current and future use of ridge regression for prediction in quantitative genetics. *BioMed Research International*.

---

*Document created: 2024*  
*Last updated: Compatible with R 4.0+ and Python 3.8+*