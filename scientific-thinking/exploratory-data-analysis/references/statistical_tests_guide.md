# Statistical Tests Guide

Interpretation guidelines for common EDA statistical tests.

## Normality Tests

### Shapiro-Wilk

**Use**: Small to medium samples (n < 5000)

**H0**: Data is normal | **H1**: Data is not normal

**Interpretation**: p > 0.05 → likely normal | p ≤ 0.05 → not normal

**Note**: Very sensitive to sample size; small deviations may be significant in large samples

### Anderson-Darling

**Use**: More powerful than Shapiro-Wilk, emphasizes tails

**Interpretation**: Test statistic > critical value → reject normality

### Kolmogorov-Smirnov

**Use**: Large samples or testing against non-normal distributions

**Interpretation**: p > 0.05 → matches reference | p ≤ 0.05 → differs from reference

## Distribution Characteristics

### Skewness

**Measures asymmetry**:
- ≈ 0: Symmetric
- \> 0: Right-skewed (tail right)
- < 0: Left-skewed (tail left)

**Magnitude**: |s| < 0.5 (symmetric) | 0.5-1 (moderate) | ≥ 1 (high)

**Action**: High skew → consider transformation (log, sqrt, Box-Cox); use median over mean

### Kurtosis

**Measures tailedness** (excess kurtosis, normal = 0):
- ≈ 0: Normal tails
- \> 0: Heavy tails, more outliers
- < 0: Light tails, fewer outliers

**Magnitude**: |k| < 0.5 (normal) | 0.5-1 (moderate) | ≥ 1 (very different)

**Action**: High kurtosis → investigate outliers carefully

## Correlation

### Pearson

**Measures**: Linear relationship (-1 to +1)

**Strength**: |r| < 0.3 (weak) | 0.3-0.5 (moderate) | 0.5-0.7 (strong) | ≥ 0.7 (very strong)

**Assumptions**: Linear, continuous, normal, no outliers, homoscedastic

**Use**: Expected linear relationship, assumptions met

### Spearman

**Measures**: Monotonic relationship (-1 to +1), rank-based

**Advantages**: Robust to outliers, no linearity assumption, works with ordinal, no normality required

**Use**: Outliers present, non-linear monotonic relationship, ordinal data, non-normal

## Outlier Detection

### IQR Method

**Bounds**: Q1 - 1.5×IQR to Q3 + 1.5×IQR

**Characteristics**: Simple, robust, works with skewed data

**Typical Rates**: < 5% (normal) | 5-10% (moderate) | > 10% (high, investigate)

### Z-Score Method

**Definition**: |z| > 3 where z = (x - μ) / σ

**Use**: Normal data, n > 30

**Avoid**: Small samples, skewed data, many outliers (contaminates mean/SD)

## Hypothesis Testing

**Significance Levels**: α = 0.05 (standard) | 0.01 (conservative) | 0.10 (liberal)

**p-value Interpretation**: ≤ 0.001 (***) | ≤ 0.01 (**) | ≤ 0.05 (*) | ≤ 0.10 (weak) | > 0.10 (none)

**Key Considerations**:
- Statistical ≠ practical significance
- Multiple testing → use correction (Bonferroni, FDR)
- Large samples detect trivial effects
- Always report effect sizes with p-values

## Transformations

**Right-skewed**: Log, sqrt, Box-Cox

**Left-skewed**: Square, cube, exponential

**Heavy tails**: Robust scaling, winsorization, log

**Non-constant variance**: Log, Box-Cox

**Common Methods**:
- **Log**: log(x+1) for positive skew, multiplicative relationships
- **Sqrt**: Count data, moderate skew
- **Box-Cox**: Auto-finds optimal (requires positive values)
- **Standardization**: (x-μ)/σ for scaling to unit variance
- **Min-Max**: (x-min)/(max-min) for [0,1] scaling

## Practical Guidelines

**Sample Size**: n < 30 (non-parametric, cautious) | 30-100 (parametric OK) | ≥ 100 (robust) | ≥ 1000 (may detect trivial effects)

**Missing Data**: < 5% (simple methods) | 5-10% (imputation) | > 10% (investigate patterns, advanced methods)

**Reporting**: Include test statistic, p-value, CI, effect size, n, assumption checks
