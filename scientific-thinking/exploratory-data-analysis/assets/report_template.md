# EDA Report: [Dataset Name]

**Date**: [Date] | **Analyst**: [Name]

## Executive Summary

[Concise summary of key findings and recommendations]

**Key Findings**:
- [Finding 1]
- [Finding 2]
- [Finding 3]

**Recommendations**:
- [Recommendation 1]
- [Recommendation 2]

---

## 1. Dataset Overview

**Source**: [Source name] | **Format**: [CSV/Excel/JSON/etc.] | **Period**: [Date range]

**Structure**: [Rows] observations × [Columns] variables | **Memory**: [Size] MB

**Variable Types**:
- Numeric ([Count]): [List names]
- Categorical ([Count]): [List names]
- Datetime ([Count]): [List names]
- Boolean ([Count]): [List names]

---

## 2. Data Quality

**Completeness**: [Percentage]% | **Duplicates**: [Count] ([%]%)

**Missing Data**:
| Column | Missing % | Assessment |
|--------|-----------|------------|
| [Column 1] | [%] | [High/Medium/Low] |
| [Column 2] | [%] | [High/Medium/Low] |

![Missing Data](path/to/missing_data.png)

**Quality Issues**:
- [Issue 1]
- [Issue 2]

---

## 3. Univariate Analysis

### Numeric: [Variable Name]

**Stats**: Mean: [Value] | Median: [Value] | Std: [Value] | Range: [[Min]-[Max]]

**Distribution**: Skewness: [Value] | Kurtosis: [Value] | Normality: [Yes/No]

**Outliers**: IQR: [Count] ([%]%) | Z-score: [Count] ([%]%)

![Distribution](path/to/distribution.png)

**Insights**: [Key observations]

### Categorical: [Variable Name]

**Stats**: [Count] unique values | Most common: [Value] ([%]%) | Balance: [Balanced/Imbalanced]

| Category | Count | % |
|----------|-------|---|
| [Cat 1] | [Count] | [%] |
| [Cat 2] | [Count] | [%] |

![Distribution](path/to/categorical.png)

**Insights**: [Key observations]

### Temporal: [Variable Name]

**Range**: [Start] to [End] ([Duration]) | **Trend**: [Increasing/Decreasing/Stable] | **Seasonality**: [Yes/No]

![Time Series](path/to/timeseries.png)

**Insights**: [Key observations]

---

## 4. Bivariate Analysis

**Correlation Summary**: [Count] strong positive | [Count] strong negative | [Count] weak/none

![Correlation Heatmap](path/to/correlation_heatmap.png)

**Notable Correlations**:
| Var 1 | Var 2 | Pearson | Spearman | Strength |
|-------|-------|---------|----------|----------|
| [Var 1] | [Var 2] | [Value] | [Value] | [Strong/Moderate/Weak] |

**Insights**: [Multicollinearity issues, feature engineering opportunities]

### Key Relationship: [Var 1] vs [Var 2]

**Type**: [Linear/Non-linear/None] | **r**: [Value] | **p-value**: [Value]

![Scatter Plot](path/to/scatter.png)

**Insights**: [Description and implications]

---

## 5. Multivariate Analysis

![Scatter Matrix](path/to/scatter_matrix.png)

**Patterns**: [Key observations]

**Clustering** (if performed): [Method] | [Count] clusters identified

---

## 6. Outliers

**Overall Rate**: [%]%

| Variable | Outlier % | Method | Action |
|----------|-----------|--------|--------|
| [Var 1] | [%] | [IQR/Z-score] | [Keep/Investigate/Remove] |
| [Var 2] | [%] | [IQR/Z-score] | [Keep/Investigate/Remove] |

![Box Plots](path/to/boxplots.png)

**Investigation**: [Description of significant outliers, causes, validity]

---

## 7. Key Insights

**Data Quality**:
- [Insight with implication]
- [Insight with implication]

**Statistical Patterns**:
- [Insight with implication]
- [Insight with implication]

**Domain/Research Insights**:
- [Insight with implication]
- [Insight with implication]

**Unexpected Findings**:
- [Finding and significance]

---

## 8. Recommendations

**Data Quality Actions**:
- [ ] [Action - priority]
- [ ] [Action - priority]

**Next Steps**:
- [Step with rationale]
- [Step with rationale]

**Feature Engineering**:
- [Opportunity]
- [Opportunity]

**Modeling Considerations**:
- [Consideration]
- [Consideration]

---

## 9. Limitations

**Data**: [Key limitations]

**Analysis**: [Key limitations]

**Assumptions**: [Key assumptions made]

---

## Appendices

### A: Technical Details

**Environment**: Python with pandas, numpy, scipy, matplotlib, seaborn

**Scripts**: [Repository/location]

### B: Variable Dictionary

| Variable | Type | Description | Unit | Range | Missing % |
|----------|------|-------------|------|-------|-----------|
| [Var 1] | [Type] | [Description] | [Unit] | [Range] | [%] |

### C: Statistical Tests

**Normality**:
| Variable | Test | Statistic | p-value | Result |
|----------|------|-----------|---------|--------|
| [Var 1] | Shapiro-Wilk | [Value] | [Value] | [Normal/Non-normal] |

**Correlations**:
| Var 1 | Var 2 | r | p-value | Significant |
|-------|-------|---|---------|-------------|
| [Var 1] | [Var 2] | [Value] | [Value] | [Yes/No] |

### D: Visualizations

1. [Description](path/to/viz1.png)
2. [Description](path/to/viz2.png)
