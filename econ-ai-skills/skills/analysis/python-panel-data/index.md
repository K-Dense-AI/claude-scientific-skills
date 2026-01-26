---
name: python-panel-data
description: Panel data analysis with Python using linearmodels and pandas.
workflow_stage: analysis
compatibility:
  - claude-code
  - cursor
  - codex
  - gemini-cli
author: Awesome Econ AI Community
version: 1.0.0
tags:
  - python
  - pandas
  - linearmodels
  - panel-data
---

# Python Panel Data

## Purpose

This skill helps economists run panel data models in Python using `pandas`, `statsmodels`, and `linearmodels`, with correct fixed effects, clustering, and diagnostics.

## When to Use

- Estimating fixed effects or random effects models
- Running difference-in-differences on panel data
- Creating regression tables and plots in Python

## Instructions

Follow these steps to complete the task:

### Step 1: Understand the Context

Before generating any code, ask the user:

- What is the unit of observation and panel identifiers?
- Which outcomes and regressors are required?
- What fixed effects or time effects are needed?
- How should standard errors be clustered?

### Step 2: Generate the Output

Based on the context, generate Python code that:

1. **Loads and cleans the data** with `pandas`
2. **Sets a MultiIndex** for panel structure
3. **Fits the model** using `linearmodels.PanelOLS` or `RandomEffects`
4. **Outputs results** in a readable table and optional LaTeX

### Step 3: Verify and Explain

After generating output:

- Interpret key coefficients
- Note assumptions (strict exogeneity, parallel trends, etc.)
- Suggest robustness checks (alternative clustering, placebo tests)

## Example Prompts

- "Run a two-way fixed effects model with firm and year effects"
- "Estimate a DiD using state and year fixed effects"
- "Export panel regression results to LaTeX"

## Example Output

```python
# ============================================
# Panel Data Analysis in Python
# ============================================
import pandas as pd
from linearmodels.panel import PanelOLS

# Load data
df = pd.read_csv("panel_data.csv")

# Set panel index
df = df.set_index(["firm_id", "year"])

# Create treatment indicator
df["treat_post"] = df["treated"] * df["post"]

# Two-way fixed effects model
model = PanelOLS.from_formula(
    "outcome ~ 1 + treat_post + EntityEffects + TimeEffects",
    data=df
)
results = model.fit(cov_type="clustered", cluster_entity=True)

print(results.summary)
```

## Requirements

### Software

- Python 3.10+

### Packages

- `pandas`
- `linearmodels`
- `statsmodels`

Install with:

```bash
pip install pandas linearmodels statsmodels
```

## Best Practices

1. **Always verify panel identifiers** and balanced vs unbalanced panels
2. **Cluster standard errors** at the appropriate level
3. **Check for missing data** before estimation

## Common Pitfalls

- Failing to set a proper panel index
- Using pooled OLS when fixed effects are required
- Misinterpreting coefficients without accounting for fixed effects

## References

- [linearmodels documentation](https://bashtage.github.io/linearmodels/)
- [statsmodels documentation](https://www.statsmodels.org/)
- [Wooldridge (2010) Econometric Analysis of Cross Section and Panel Data](https://mitpress.mit.edu/9780262232586/)

## Changelog

### v1.0.0

- Initial release
