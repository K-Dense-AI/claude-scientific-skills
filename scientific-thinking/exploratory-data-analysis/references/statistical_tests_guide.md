# Statistical Tests Guide for EDA

This guide provides interpretation guidelines for statistical tests commonly used in exploratory data analysis.

## Normality Tests

### Shapiro-Wilk Test

**Purpose**: Test if a sample comes from a normally distributed population

**When to use**: Best for small to medium sample sizes (n < 5000)

**Interpretation**:
- **Null Hypothesis (H0)**: The data follows a normal distribution
- **Alternative Hypothesis (H1)**: The data does not follow a normal distribution
- **p-value > 0.05**: Fail to reject H0 → Data is likely normally distributed
- **p-value ≤ 0.05**: Reject H0 → Data is not normally distributed

**Notes**:
- Very sensitive to sample size
- Small deviations from normality may be detected as significant in large samples
- Consider practical significance alongside statistical significance

### Anderson-Darling Test

**Purpose**: Test if a sample comes from a specific distribution (typically normal)

**When to use**: More powerful than Shapiro-Wilk for detecting departures from normality

**Interpretation**:
- Compares test statistic against critical values at different significance levels
- If test statistic > critical value at given significance level, reject normality
- More weight given to tails of distribution than other tests

### Kolmogorov-Smirnov Test

**Purpose**: Test if a sample comes from a reference distribution

**When to use**: When you have a large sample or want to test against distributions other than normal

**Interpretation**:
- **p-value > 0.05**: Sample distribution matches reference distribution
- **p-value ≤ 0.05**: Sample distribution differs from reference distribution

## Distribution Characteristics

### Skewness

**Purpose**: Measure asymmetry of the distribution

**Interpretation**:
- **Skewness ≈ 0**: Symmetric distribution
- **Skewness > 0**: Right-skewed (tail extends to right, most values on left)
- **Skewness < 0**: Left-skewed (tail extends to left, most values on right)

**Magnitude interpretation**:
- **|Skewness| < 0.5**: Approximately symmetric
- **0.5 ≤ |Skewness| < 1**: Moderately skewed
- **|Skewness| ≥ 1**: Highly skewed

**Implications**:
- Highly skewed data may require transformation (log, sqrt, Box-Cox)
- Mean is pulled toward tail; median more robust for skewed data
- Many statistical tests assume symmetry/normality

### Kurtosis

**Purpose**: Measure tailedness and peak of distribution

**Interpretation** (Excess Kurtosis, where normal distribution = 0):
- **Kurtosis ≈ 0**: Normal tail behavior (mesokurtic)
- **Kurtosis > 0**: Heavy tails, sharp peak (leptokurtic)
  - More outliers than normal distribution
  - Higher probability of extreme values
- **Kurtosis < 0**: Light tails, flat peak (platykurtic)
  - Fewer outliers than normal distribution
  - More uniform distribution

**Magnitude interpretation**:
- **|Kurtosis| < 0.5**: Normal-like tails
- **0.5 ≤ |Kurtosis| < 1**: Moderately different tails
- **|Kurtosis| ≥ 1**: Very different tail behavior from normal

**Implications**:
- High kurtosis → Be cautious with outliers
- Low kurtosis → Distribution lacks distinct peak

## Correlation Tests

### Pearson Correlation

**Purpose**: Measure linear relationship between two continuous variables

**Range**: -1 to +1

**Interpretation**:
- **r = +1**: Perfect positive linear relationship
- **r = 0**: No linear relationship
- **r = -1**: Perfect negative linear relationship

**Strength guidelines**:
- **|r| < 0.3**: Weak correlation
- **0.3 ≤ |r| < 0.5**: Moderate correlation
- **0.5 ≤ |r| < 0.7**: Strong correlation
- **|r| ≥ 0.7**: Very strong correlation

**Assumptions**:
- Linear relationship between variables
- Both variables continuous and normally distributed
- No significant outliers
- Homoscedasticity (constant variance)

**When to use**: When relationship is expected to be linear and data meets assumptions

### Spearman Correlation

**Purpose**: Measure monotonic relationship between two variables (rank-based)

**Range**: -1 to +1

**Interpretation**: Same as Pearson, but measures monotonic (not just linear) relationships

**Advantages over Pearson**:
- Robust to outliers (uses ranks)
- Doesn't assume linear relationship
- Works with ordinal data
- Doesn't require normality assumption

**When to use**:
- Data has outliers
- Relationship is monotonic but not linear
- Data is ordinal
- Distribution is non-normal

## Outlier Detection Methods

### IQR Method (Interquartile Range)

**Definition**:
- Lower bound: Q1 - 1.5 × IQR
- Upper bound: Q3 + 1.5 × IQR
- Values outside these bounds are outliers

**Characteristics**:
- Simple and interpretable
- Robust to extreme values
- Works well for skewed distributions
- Conservative approach (Tukey's fences)

**Interpretation**:
- **< 5% outliers**: Typical for most datasets
- **5-10% outliers**: Moderate, investigate causes
- **> 10% outliers**: High rate, may indicate data quality issues or interesting phenomena

### Z-Score Method

**Definition**: Outliers are data points with |z-score| > 3

**Formula**: z = (x - μ) / σ

**Characteristics**:
- Assumes normal distribution
- Sensitive to extreme values
- Standard threshold is |z| > 3 (99.7% of data within ±3σ)

**When to use**:
- Data is approximately normally distributed
- Large sample sizes (n > 30)

**When NOT to use**:
- Small samples
- Heavily skewed data
- Data with many outliers (contaminates mean and SD)

## Hypothesis Testing Guidelines

### Significance Levels

- **α = 0.05**: Standard significance level (5% chance of Type I error)
- **α = 0.01**: More conservative (1% chance of Type I error)
- **α = 0.10**: More liberal (10% chance of Type I error)

### p-value Interpretation

- **p ≤ 0.001**: Very strong evidence against H0 (***)
- **0.001 < p ≤ 0.01**: Strong evidence against H0 (**)
- **0.01 < p ≤ 0.05**: Moderate evidence against H0 (*)
- **0.05 < p ≤ 0.10**: Weak evidence against H0
- **p > 0.10**: Little to no evidence against H0

### Important Considerations

1. **Statistical vs Practical Significance**: A small p-value doesn't always mean the effect is important
2. **Multiple Testing**: When performing many tests, use correction methods (Bonferroni, FDR)
3. **Sample Size**: Large samples can detect trivial effects as significant
4. **Effect Size**: Always report and interpret effect sizes alongside p-values

## Data Transformation Strategies

### When to Transform

- **Right-skewed data**: Log, square root, or Box-Cox transformation
- **Left-skewed data**: Square, cube, or exponential transformation
- **Heavy tails/outliers**: Robust scaling, winsorization, or log transformation
- **Non-constant variance**: Log or Box-Cox transformation

### Common Transformations

1. **Log transformation**: log(x) or log(x + 1)
   - Best for: Positive skewed data, multiplicative relationships
   - Cannot use with zero or negative values

2. **Square root transformation**: √x
   - Best for: Count data, moderate positive skew
   - Less aggressive than log

3. **Box-Cox transformation**: (x^λ - 1) / λ
   - Best for: Automatically finds optimal transformation
   - Requires positive values

4. **Standardization**: (x - μ) / σ
   - Best for: Scaling features to same range
   - Centers data at 0 with unit variance

5. **Min-Max scaling**: (x - min) / (max - min)
   - Best for: Scaling to [0, 1] range
   - Preserves zero values

## Practical Guidelines

### Sample Size Considerations

- **n < 30**: Use non-parametric tests, be cautious with assumptions
- **30 ≤ n < 100**: Moderate sample, parametric tests usually acceptable
- **n ≥ 100**: Large sample, parametric tests robust to violations
- **n ≥ 1000**: Very large sample, may detect trivial effects as significant

### Dealing with Missing Data

- **< 5% missing**: Usually not a problem, simple methods OK
- **5-10% missing**: Use appropriate imputation methods
- **> 10% missing**: Investigate patterns, consider advanced imputation or modeling missingness

### Reporting Results

Always include:
1. Test statistic value
2. p-value
3. Confidence interval (when applicable)
4. Effect size
5. Sample size
6. Assumptions checked and violations noted
