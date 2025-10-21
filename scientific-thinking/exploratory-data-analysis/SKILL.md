---
name: exploratory-data-analysis
description: "EDA toolkit. Analyze CSV/Excel/JSON/Parquet files, statistical summaries, distributions, correlations, outliers, missing data, visualizations, markdown reports, for data profiling and insights."
---

# Exploratory Data Analysis

## Overview

EDA is a process for discovering patterns, anomalies, and relationships in data. Analyze CSV/Excel/JSON/Parquet files to generate statistical summaries, distributions, correlations, outliers, and visualizations. All outputs are markdown-formatted for integration into workflows.

## When to Use This Skill

This skill should be used when:
- User provides a data file and requests analysis or exploration
- User asks to "explore this dataset", "analyze this data", or "what's in this file?"
- User needs statistical summaries, distributions, or correlations
- User requests data visualizations or insights
- User wants to understand data quality issues or patterns
- User mentions EDA, exploratory analysis, or data profiling

**Supported file formats**: CSV, Excel (.xlsx, .xls), JSON, Parquet, TSV, Feather, HDF5, Pickle

## Quick Start Workflow

1. **Receive data file** from user
2. **Run comprehensive analysis** using `scripts/eda_analyzer.py`
3. **Generate visualizations** using `scripts/visualizer.py`
4. **Create markdown report** using insights and the `assets/report_template.md` template
5. **Present findings** to user with key insights highlighted

## Core Capabilities

### 1. Comprehensive Data Analysis

Execute full statistical analysis using the `eda_analyzer.py` script:

```bash
python scripts/eda_analyzer.py <data_file_path> -o <output_directory>
```

**What it provides**:
- Auto-detection and loading of file formats
- Basic dataset information (shape, types, memory usage)
- Missing data analysis (patterns, percentages)
- Summary statistics for numeric and categorical variables
- Outlier detection using IQR and Z-score methods
- Distribution analysis with normality tests (Shapiro-Wilk, Anderson-Darling)
- Correlation analysis (Pearson and Spearman)
- Data quality assessment (completeness, duplicates, issues)
- Automated insight generation

**Output**: JSON file containing all analysis results at `<output_directory>/eda_analysis.json`

### 2. Comprehensive Visualizations

Generate complete visualization suite using the `visualizer.py` script:

```bash
python scripts/visualizer.py <data_file_path> -o <output_directory>
```

**Generated visualizations**:
- **Missing data patterns**: Heatmap and bar chart showing missing data
- **Distribution plots**: Histograms with KDE overlays for all numeric variables
- **Box plots with violin plots**: Outlier detection visualizations
- **Correlation heatmap**: Both Pearson and Spearman correlation matrices
- **Scatter matrix**: Pairwise relationships between numeric variables
- **Categorical analysis**: Bar charts for top categories
- **Time series plots**: Temporal trends with trend lines (if datetime columns exist)

**Output**: High-quality PNG files saved to `<output_directory>/eda_visualizations/`

All visualizations are production-ready with:
- 300 DPI resolution
- Clear titles and labels
- Statistical annotations
- Professional styling using seaborn

### 3. Automated Insight Generation

The analyzer automatically generates actionable insights including:

- **Data scale insights**: Dataset size considerations for processing
- **Missing data alerts**: Warnings when missing data exceeds thresholds
- **Correlation discoveries**: Strong relationships identified for feature engineering
- **Outlier warnings**: Variables with high outlier rates flagged
- **Distribution assessments**: Skewness issues requiring transformations
- **Duplicate alerts**: Duplicate row detection
- **Imbalance warnings**: Categorical variable imbalance detection

Access insights from the analysis results JSON under the `"insights"` key.

### 4. Statistical Interpretation

For detailed interpretation of statistical tests and measures, reference:

**`references/statistical_tests_guide.md`** - Comprehensive guide covering:
- Normality tests (Shapiro-Wilk, Anderson-Darling, Kolmogorov-Smirnov)
- Distribution characteristics (skewness, kurtosis)
- Correlation tests (Pearson, Spearman)
- Outlier detection methods (IQR, Z-score)
- Hypothesis testing guidelines
- Data transformation strategies

Load this reference when needing to interpret specific statistical tests or explain results to users.

### 5. Best Practices Guidance

For methodological guidance, reference:

**`references/eda_best_practices.md`** - Detailed best practices including:
- EDA process framework (6-step methodology)
- Univariate, bivariate, and multivariate analysis approaches
- Visualization guidelines
- Statistical analysis guidelines
- Common pitfalls to avoid
- Domain-specific considerations
- Communication tips for technical and non-technical audiences

Load this reference when planning analysis approach or needing guidance on specific EDA scenarios.

## Creating Analysis Reports

Use the provided template to structure comprehensive EDA reports:

**`assets/report_template.md`** - Professional report template with sections for:
- Executive summary
- Dataset overview
- Data quality assessment
- Univariate, bivariate, and multivariate analysis
- Outlier analysis
- Key insights and findings
- Recommendations
- Limitations and appendices

**To use the template**:
1. Copy the template content
2. Fill in sections with analysis results from JSON output
3. Embed visualization images using markdown syntax
4. Populate insights and recommendations
5. Save as markdown for user consumption

## Typical Workflow Example

When user provides a data file:

```
User: "Can you explore this sales_data.csv file and tell me what you find?"

1. Run analysis:
   python scripts/eda_analyzer.py sales_data.csv -o ./analysis_output

2. Generate visualizations:
   python scripts/visualizer.py sales_data.csv -o ./analysis_output

3. Read analysis results:
   Read ./analysis_output/eda_analysis.json

4. Create markdown report using template:
   - Copy assets/report_template.md structure
   - Fill in sections with analysis results
   - Reference visualizations from ./analysis_output/eda_visualizations/
   - Include automated insights from JSON

5. Present to user:
   - Show key insights prominently
   - Highlight data quality issues
   - Provide visualizations inline
   - Make actionable recommendations
   - Save complete report as .md file
```

## Advanced Analysis Scenarios

### Large Datasets (>1M rows)
- Run analysis on sampled data first for quick exploration
- Note sample size in report
- Recommend distributed computing for full analysis

### High-Dimensional Data (>50 columns)
- Focus on most important variables first
- Consider PCA or feature selection
- Generate correlation analysis to identify variable groups
- Reference `eda_best_practices.md` section on high-dimensional data

### Time Series Data
- Ensure datetime columns are properly detected
- Time series visualizations will be automatically generated
- Consider temporal patterns, trends, and seasonality
- Reference `eda_best_practices.md` section on time series

### Imbalanced Data
- Categorical analysis will flag imbalances
- Report class distributions prominently
- Recommend stratified sampling if needed

### Small Sample Sizes (<100 rows)
- Non-parametric methods automatically used where appropriate
- Be conservative in statistical conclusions
- Note sample size limitations in report

## Output Best Practices

**Always output as markdown**:
- Structure findings using markdown headers, tables, and lists
- Embed visualizations using `![Description](path/to/image.png)` syntax
- Use tables for statistical summaries
- Include code blocks for any suggested transformations
- Highlight key insights with bold or bullet points

**Ensure reports are actionable**:
- Provide clear recommendations based on findings
- Flag data quality issues that need attention
- Suggest next steps for modeling or further analysis
- Identify feature engineering opportunities

**Make insights accessible**:
- Explain statistical concepts in plain language
- Use reference guides to provide detailed interpretations
- Include both technical details and executive summary
- Tailor communication to user's technical level

## Handling Edge Cases

**Unsupported file formats**:
- Request user to convert to supported format
- Suggest using pandas-compatible formats

**Files too large to load**:
- Recommend sampling approach
- Suggest chunked processing
- Consider alternative tools for big data

**Corrupted or malformed data**:
- Report specific errors encountered
- Suggest data cleaning steps
- Try to salvage partial analysis if possible

**All missing data in columns**:
- Flag completely empty columns
- Recommend removal or investigation
- Document in data quality section

## Resources Summary

### scripts/
- **`eda_analyzer.py`**: Main analysis engine - comprehensive statistical analysis
- **`visualizer.py`**: Visualization generator - creates all chart types

Both scripts are fully executable and handle multiple file formats automatically.

### references/
- **`statistical_tests_guide.md`**: Statistical test interpretation and methodology
- **`eda_best_practices.md`**: Comprehensive EDA methodology and best practices

Load these references as needed to inform analysis approach and interpretation.

### assets/
- **`report_template.md`**: Professional markdown report template

Use this template structure for creating consistent, comprehensive EDA reports.

## Key Reminders

1. **Always generate markdown output** for textual results
2. **Run both scripts** (analyzer and visualizer) for complete analysis
3. **Use the template** to structure comprehensive reports
4. **Include visualizations** by referencing generated PNG files
5. **Provide actionable insights** - don't just present statistics
6. **Interpret findings** using reference guides
7. **Document limitations** and data quality issues
8. **Make recommendations** for next steps

This skill transforms raw data into actionable insights through systematic exploration, advanced statistics, rich visualizations, and clear communication.
