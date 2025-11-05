---
name: exploratory-data-analysis
description: "Analyze datasets to discover patterns, anomalies, and relationships. Use when exploring data files, generating statistical summaries, checking data quality, or creating visualizations. Supports CSV, Excel, JSON, Parquet, and more."
---

# Exploratory Data Analysis

Discover patterns, anomalies, and relationships in tabular data through statistical analysis and visualization.

**Supported formats**: CSV, Excel (.xlsx, .xls), JSON, Parquet, TSV, Feather, HDF5, Pickle

## Standard Workflow

1. Run statistical analysis:
```bash
python scripts/eda_analyzer.py <data_file> -o <output_dir>
```

2. Generate visualizations:
```bash
python scripts/visualizer.py <data_file> -o <output_dir>
```

3. Read analysis results from `<output_dir>/eda_analysis.json`

4. Create report using `assets/report_template.md` structure

5. Present findings with key insights and visualizations

## Analysis Capabilities

### Statistical Analysis

Run `scripts/eda_analyzer.py` to generate comprehensive analysis:

```bash
python scripts/eda_analyzer.py sales_data.csv -o ./output
```

Produces `output/eda_analysis.json` containing:
- Dataset shape, types, memory usage
- Missing data patterns and percentages
- Summary statistics (numeric and categorical)
- Outlier detection (IQR and Z-score methods)
- Distribution analysis with normality tests
- Correlation matrices (Pearson and Spearman)
- Data quality metrics (completeness, duplicates)
- Automated insights

### Visualizations

Run `scripts/visualizer.py` to generate plots:

```bash
python scripts/visualizer.py sales_data.csv -o ./output
```

Creates high-resolution (300 DPI) PNG files in `output/eda_visualizations/`:
- Missing data heatmaps and bar charts
- Distribution plots (histograms with KDE)
- Box plots and violin plots for outliers
- Correlation heatmaps
- Scatter matrices for numeric relationships
- Categorical bar charts
- Time series plots (if datetime columns detected)

### Automated Insights

Access generated insights from the `"insights"` key in the analysis JSON:
- Dataset size considerations
- Missing data warnings (when exceeding thresholds)
- Strong correlations for feature engineering
- High outlier rate flags
- Skewness requiring transformations
- Duplicate detection
- Categorical imbalance warnings

## Reference Materials

### Statistical Interpretation

See `references/statistical_tests_guide.md` for detailed guidance on:
- Normality tests (Shapiro-Wilk, Anderson-Darling, Kolmogorov-Smirnov)
- Distribution characteristics (skewness, kurtosis)
- Correlation methods (Pearson, Spearman)
- Outlier detection (IQR, Z-score)
- Hypothesis testing and data transformations

Use when interpreting statistical results or explaining findings.

### Methodology

See `references/eda_best_practices.md` for comprehensive guidance on:
- 6-step EDA process framework
- Univariate, bivariate, multivariate analysis approaches
- Visualization and statistical analysis guidelines
- Common pitfalls and domain-specific considerations
- Communication strategies for different audiences

Use when planning analysis or handling specific scenarios.

## Report Template

Use `assets/report_template.md` to structure findings. Template includes:
- Executive summary
- Dataset overview
- Data quality assessment
- Univariate, bivariate, and multivariate analysis
- Outlier analysis
- Key insights and recommendations
- Limitations and appendices

Fill sections with analysis JSON results and embed visualizations using markdown image syntax.

## Example: Complete Analysis

User request: "Explore this sales_data.csv file"

```bash
# 1. Run analysis
python scripts/eda_analyzer.py sales_data.csv -o ./output

# 2. Generate visualizations
python scripts/visualizer.py sales_data.csv -o ./output
```

```python
# 3. Read results
import json
with open('./output/eda_analysis.json') as f:
    results = json.load(f)

# 4. Build report from assets/report_template.md
# - Fill sections with results
# - Embed images: ![Missing Data](./output/eda_visualizations/missing_data.png)
# - Include insights from results['insights']
# - Add recommendations
```

## Special Cases

### Dataset Size Strategy

**If < 100 rows**: Note sample size limitations, use non-parametric methods

**If 100-1M rows**: Standard workflow applies

**If > 1M rows**: Sample first for quick exploration, note sample size in report, recommend distributed computing for full analysis

### Data Characteristics

**High-dimensional (>50 columns)**: Focus on key variables first, use correlation analysis to identify groups, consider PCA or feature selection. See `references/eda_best_practices.md` for guidance.

**Time series**: Datetime columns auto-detected, temporal visualizations generated automatically. Consider trends, seasonality, patterns.

**Imbalanced**: Categorical analysis flags imbalances automatically. Report distributions prominently, recommend stratified sampling if needed.

## Output Guidelines

**Format findings as markdown**:
- Use headers, tables, and lists for structure
- Embed visualizations: `![Description](path/to/image.png)`
- Include code blocks for suggested transformations
- Highlight key insights

**Make reports actionable**:
- Provide clear recommendations
- Flag data quality issues requiring attention
- Suggest next steps (modeling, feature engineering, further analysis)
- Tailor communication to user's technical level

## Error Handling

**Unsupported formats**: Request conversion to supported format (CSV, Excel, JSON, Parquet)

**Files too large**: Recommend sampling or chunked processing

**Corrupted data**: Report specific errors, suggest cleaning steps, attempt partial analysis

**Empty columns**: Flag in data quality section, recommend removal or investigation

## Resources

**Scripts** (handle all formats automatically):
- `scripts/eda_analyzer.py` - Statistical analysis engine
- `scripts/visualizer.py` - Visualization generator

**References** (load as needed):
- `references/statistical_tests_guide.md` - Test interpretation and methodology
- `references/eda_best_practices.md` - EDA process and best practices

**Template**:
- `assets/report_template.md` - Professional report structure

## Key Points

- Run both scripts for complete analysis
- Structure reports using the template
- Provide actionable insights, not just statistics
- Use reference guides for detailed interpretations
- Document data quality issues and limitations
- Make clear recommendations for next steps
