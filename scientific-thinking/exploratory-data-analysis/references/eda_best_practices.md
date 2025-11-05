# EDA Best Practices

Methodologies for conducting thorough exploratory data analysis.

## 6-Step EDA Framework

### 1. Initial Understanding

**Questions**:
- What does each column represent?
- What is the unit of observation and time period?
- What is the data collection methodology?
- Are there known quality issues?

**Actions**: Load data, inspect structure, review types, document context

### 2. Quality Assessment

**Check**: Missing data patterns, duplicates, outliers, consistency, accuracy

**Red Flags**:
- Missing data >20%
- Unexpected duplicates
- Constant columns
- Impossible values (negative ages, future dates)
- Suspicious patterns (too many round numbers)

### 3. Univariate Analysis

**Numeric**: Central tendency, dispersion, shape (skewness, kurtosis), distribution plots, outliers

**Categorical**: Frequency distributions, unique counts, balance, bar charts

**Temporal**: Time range, gaps, trends, seasonality, time series plots

### 4. Bivariate Analysis

**Numeric vs Numeric**: Scatter plots, correlations (Pearson, Spearman), detect non-linearity

**Numeric vs Categorical**: Group statistics, box plots by category, t-test/ANOVA

**Categorical vs Categorical**: Cross-tabs, stacked bars, chi-square, Cramér's V

### 5. Multivariate Analysis

**Techniques**: Correlation matrices, pair plots, parallel coordinates, PCA, clustering

**Questions**: Groups of correlated features? Reduce dimensionality? Natural clusters? Conditional patterns?

### 6. Insight Generation

**Look for**: Unexpected patterns, strong correlations, quality issues, feature engineering opportunities, domain implications

## Visualization Guidelines

**Chart Selection**:
- Distribution: Histogram, KDE, box/violin plots
- Relationships: Scatter, line, heatmap
- Composition: Stacked bar
- Comparison: Bar, grouped bar

**Best Practices**: Label axes with units, descriptive titles, purposeful color, appropriate scales, avoid clutter

## Statistical Analysis Guidelines

**Check Assumptions**: Normality, homoscedasticity, independence, linearity

**Method Selection**: Parametric when assumptions met, non-parametric otherwise, report effect sizes

**Context Matters**: Statistical ≠ practical significance, domain knowledge trumps statistics, correlation ≠ causation

## Documentation Guidelines

**Notes**: Document assumptions, decisions, issues, findings

**Reproducibility**: Use scripts, version control, document sources, set random seeds

**Reporting**: Clear summaries, supporting visualizations, highlighted insights, actionable recommendations

## Common Pitfalls

1. **Confirmation Bias**: Seek disconfirming evidence, use blind analysis
2. **Ignoring Quality**: Address issues first, document limitations
3. **Over-automation**: Manually inspect subsets, verify results
4. **Neglecting Outliers**: Investigate before removing - may be informative
5. **Multiple Testing**: Use correction (Bonferroni, FDR) or note exploratory nature
6. **Association ≠ Causation**: Use careful language, acknowledge alternatives
7. **Cherry-picking**: Report complete analysis, including negative results
8. **Ignoring Sample Size**: Report effect sizes, CIs, and sample sizes

## Domain-Specific Considerations

**Time Series**: Check stationarity, identify trends/seasonality, autocorrelation, temporal splits

**High-Dimensional**: Dimensionality reduction, feature importance, regularization, domain-guided selection

**Imbalanced**: Report distributions, appropriate metrics, resampling, stratified CV

**Small Samples**: Non-parametric methods, conservative conclusions, CIs, Bayesian approaches

**Big Data**: Intelligent sampling, efficient structures, parallel computing, scalability

## Iterative Process

EDA is iterative: Explore → Questions → Focused Analysis → Insights → New Questions → Deeper Investigation → Synthesis

**Done When**: Understand structure/quality, characterized variables, identified relationships, documented limitations, answered questions, have actionable insights

**Deliverables**: Data understanding, quality issue list, relationship insights, hypotheses, feature ideas, recommendations

## Communication

**Technical Audiences**: Methodological details, statistical tests, assumptions, reproducible code

**Non-Technical Audiences**: Focus on insights, clear visualizations, avoid jargon, concrete recommendations

**Report Structure**: Executive summary → Data overview → Analysis → Insights → Recommendations → Appendix

## Checklists

**Before**: Understand context, define objectives, identify audience, set up environment

**During**: Inspect structure, assess quality, analyze distributions, explore relationships, document continuously

**After**: Verify findings, check alternatives, document limitations, prepare visualizations, ensure reproducibility
