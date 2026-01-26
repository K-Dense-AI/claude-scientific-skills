---
name: econ-visualization
description: Create publication-quality charts and graphs for economics papers.
workflow_stage: communication
compatibility:
  - claude-code
  - cursor
  - codex
  - gemini-cli
author: Awesome Econ AI Community
version: 1.0.0
tags:
  - visualization
  - ggplot2
  - charts
  - publication
---

# Econ Visualization

## Purpose

This skill creates publication-quality figures for economics papers, using clean styling, consistent scales, and export-ready formats.

## When to Use

- Building figures for empirical results and descriptive analysis
- Standardizing chart style across a paper or presentation
- Exporting figures to PDF or PNG at journal quality

## Instructions

Follow these steps to complete the task:

### Step 1: Understand the Context

Before generating any code, ask the user:

- What is the dataset and key variables?
- What chart type is needed (line, bar, scatter, event study)?
- What output format and size are required?

### Step 2: Generate the Output

Based on the context, generate code that:

1. **Uses a consistent theme** for academic styling
2. **Labels axes and legends clearly**
3. **Exports figures** at high resolution
4. **Includes reproducible steps** for data preparation

### Step 3: Verify and Explain

After generating output:

- Explain how to regenerate or update the plot
- Suggest alternatives (log scales, faceting, smoothing)
- Note any data transformations used

## Example Prompts

- "Create an event study plot with confidence intervals"
- "Plot GDP per capita over time for three countries"
- "Build a scatter plot with fitted regression line"

## Example Output

```r
# ============================================
# Publication-Quality Figure in R
# ============================================
library(tidyverse)

df <- read_csv("data.csv")

ggplot(df, aes(x = year, y = gdp_per_capita, color = country)) +
  geom_line(size = 1) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "GDP per Capita Over Time",
    x = "Year",
    y = "GDP per Capita (USD)",
    color = "Country"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("figures/gdp_per_capita.pdf", width = 7, height = 4, dpi = 300)
```

## Requirements

### Software

- R 4.0+ or Python 3.10+

### Packages

- For R: `ggplot2`, `scales`, `dplyr`
- For Python: `matplotlib`, `seaborn` (optional alternative)

## Best Practices

1. **Use vector formats** (PDF, SVG) for publication
2. **Keep labels concise** and readable
3. **Document data filters** used in the figure

## Common Pitfalls

- Overcrowded plots without clear labeling
- Inconsistent scales across figures
- Exporting low-resolution images

## References

- [ggplot2 documentation](https://ggplot2.tidyverse.org/)
- [Tufte (2001) The Visual Display of Quantitative Information](https://www.edwardtufte.com/tufte/books_vdqi)

## Changelog

### v1.0.0

- Initial release
