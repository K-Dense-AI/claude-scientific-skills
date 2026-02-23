# QTL Visualization Quick Reference Card

## Python Quick Start

### Installation
```bash
pip install matplotlib seaborn plotly pandas numpy scipy qmplot
```

### Essential Imports
```python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
```

---

## Quick Code Snippets

### 1. Box Plot (Python)
```python
# Basic box plot
sns.boxplot(data=df, x='Genotype', y='Phenotype')
plt.show()

# Grouped box plot
sns.boxplot(data=df, x='Genotype', y='Phenotype', hue='Environment')
plt.show()
```

### 2. Violin Plot (Python)
```python
sns.violinplot(data=df, x='Genotype', y='Phenotype', inner='box')
plt.show()
```

### 3. Density Plot (Python)
```python
sns.kdeplot(data=df, x='Phenotype', hue='Genotype', fill=True)
plt.show()
```

### 4. Heatmap (Python)
```python
# Correlation heatmap
corr = df.corr()
sns.heatmap(corr, annot=True, cmap='coolwarm', center=0)
plt.show()

# Clustermap
sns.clustermap(df.T, cmap='viridis')
plt.show()
```

### 5. Scatter Matrix (Python)
```python
sns.pairplot(df, hue='Genotype')
plt.show()
```

### 6. Manhattan Plot (Python)
```python
# Using qmplot
from qmplot import manhattanplot
manhattanplot(data=gwas_df, chrom="CHR", pos="BP", pvalue="P")
plt.show()

# Manual implementation
def manhattan_plot(df, ax):
    colors = ['blue', 'orange']
    x_pos = []
    offset = 0
    for i, chrom in enumerate(sorted(df['CHR'].unique())):
        chrom_data = df[df['CHR'] == chrom]
        x = chrom_data['BP'] / 1e6 + offset
        ax.scatter(x, -np.log10(chrom_data['P']), 
                  c=colors[i % 2], s=10)
        offset += chrom_data['BP'].max() / 1e6 + 10
    ax.axhline(y=-np.log10(5e-8), color='red', linestyle='--')
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(p)')
```

### 7. QQ Plot (Python)
```python
from scipy import stats
import numpy as np

# For GWAS p-values
sorted_p = np.sort(pvalues)
expected = np.arange(1, len(sorted_p) + 1) / len(sorted_p)
plt.scatter(-np.log10(expected), -np.log10(sorted_p), s=10)
plt.plot([0, max(-np.log10(expected))], [0, max(-np.log10(expected))], 'r--')
plt.xlabel('Expected -log10(p)')
plt.ylabel('Observed -log10(p)')
plt.show()

# For normality check
stats.probplot(data, dist="norm", plot=plt)
plt.show()
```

### 8. LOD Score Plot (Python)
```python
# Genome-wide scan
for chrom in qtl_df['Chromosome'].unique():
    chrom_data = qtl_df[qtl_df['Chromosome'] == chrom]
    plt.plot(chrom_data['Position'], chrom_data['LOD'], 
             label=f'Chr{chrom}')
plt.axhline(y=3, color='red', linestyle='--', label='Threshold')
plt.xlabel('Position')
plt.ylabel('LOD Score')
plt.legend()
plt.show()
```

---

## R Quick Start

### Installation
```r
install.packages(c("ggplot2", "qqman", "corrplot", "GGally"))
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("qtl2")
```

### Essential Imports
```r
library(ggplot2)
library(qqman)
library(corrplot)
library(qtl2)
```

---

## Quick Code Snippets (R)

### 1. Box Plot (R)
```r
# Basic
ggplot(df, aes(x=Genotype, y=Phenotype, fill=Genotype)) +
  geom_boxplot() +
  theme_minimal()

# Grouped
ggplot(df, aes(x=Genotype, y=Phenotype, fill=Environment)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Set2")
```

### 2. Violin Plot (R)
```r
ggplot(df, aes(x=Genotype, y=Phenotype, fill=Genotype)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill="white")
```

### 3. Density Plot (R)
```r
ggplot(df, aes(x=Phenotype, fill=Genotype)) +
  geom_density(alpha=0.5)
```

### 4. Heatmap (R)
```r
# Correlation
corrplot(cor(df), method="color", addCoef.col="black")

# ggplot2
ggplot(melted_df, aes(Var1, Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", high="red", mid="white")
```

### 5. Scatter Matrix (R)
```r
pairs(df)
# or with GGally
ggpairs(df, aes(color=Genotype))
```

### 6. Manhattan Plot (R)
```r
# Using qqman
manhattan(gwas_df, chr="CHR", bp="BP", snp="SNP", p="P")

# With highlighting
manhattan(gwas_df, highlight=significant_snps)

# QQ plot
qq(gwas_df$P)
```

### 7. LOD Plot (R)
```r
# Using qtl2
plot(scan1_output, map)
abline(h=3, col="red", lty=2)

# Single chromosome
plot(scan1_output, map, chr=3)
```

### 8. Effect Plot (R)
```r
ggplot(effects_df, aes(x=Effect, y=QTL)) +
  geom_point() +
  geom_errorbarh(aes(xmin=Lower, xmax=Upper), height=0.2) +
  geom_vline(xintercept=0, linetype="dashed")
```

---

## Common Customizations

### Python
```python
# Figure size
plt.figure(figsize=(10, 6))

# Save figure
plt.savefig('plot.png', dpi=300, bbox_inches='tight')

# Color palettes
sns.set_palette("viridis")  # or "Set1", "Set2", "husl"

# Themes
sns.set_style("whitegrid")  # or "darkgrid", "white", "dark"

# Labels
plt.xlabel('X Label', fontsize=12)
plt.ylabel('Y Label', fontsize=12)
plt.title('Plot Title', fontsize=14, fontweight='bold')

# Legend
plt.legend(title='Legend Title', bbox_to_anchor=(1.05, 1), loc='upper left')
```

### R
```r
# Theme customization
theme_minimal()  # or theme_bw(), theme_classic()

# Labels
labs(x="X Label", y="Y Label", title="Plot Title")

# Colors
scale_fill_brewer(palette="Set2")  # or "Set1", "Pastel1", "Dark2"
scale_color_viridis_c()  # continuous

# Save plot
ggsave("plot.pdf", width=10, height=6, dpi=300)

# Facets
facet_wrap(~Genotype)  # or facet_grid(Environment ~ Genotype)
```

---

## Data Format Requirements

### GWAS/Manhattan Plot Data
```
Required columns:
- SNP: SNP identifier
- CHR: Chromosome number
- BP: Base pair position
- P: P-value
```

### QTL Scan Data
```
Required:
- Chromosome positions (cM or bp)
- LOD scores or -log10(p)
- Marker names
```

### Phenotype Data
```
Typical structure:
- Sample/Individual ID
- Genotype (AA, AB, BB)
- Environment/Treatment
- Phenotype measurements
```

---

## Troubleshooting

### Python
- **ImportError**: Install missing package with `pip install package_name`
- **Empty plot**: Call `plt.show()` at the end
- **Overlapping labels**: Use `plt.tight_layout()` or rotate labels
- **Memory error**: Reduce data size or use `plt.close()` between plots

### R
- **Package not found**: Install with `install.packages()` or `BiocManager::install()`
- **Plot not showing**: Ensure graphics device is open
- **Legend cut off**: Adjust `ggsave()` dimensions or use `theme(legend.position="bottom")`
- **Factor levels**: Use `factor()` to control order in plots

---

## Recommended Color Palettes

### Colorblind-Friendly
- **Python**: `sns.color_palette("colorblind")`, `sns.color_palette("viridis")`
- **R**: `scale_color_brewer(palette="Dark2")`, `viridis::scale_color_viridis()`

### Categorical (≤8 categories)
- **Python**: `"Set1"`, `"Set2"`, `"Dark2"`
- **R**: `"Set1"`, `"Set2"`, `"Dark2"`, `"Paired"`

### Continuous/Diverging
- **Python**: `"viridis"`, `"coolwarm"`, `"RdBu_r"`
- **R**: `"RdBu"`, `"BrBG"`, `"PRGn"`

---

## Publication-Ready Checklist

- [ ] Figure size appropriate (typically 3.5-7 inches width)
- [ ] Resolution ≥300 DPI
- [ ] Font size ≥8pt, axis labels ≥10pt
- [ ] Colorblind-friendly palette
- [ ] Clear axis labels with units
- [ ] Informative title or caption
- [ ] Legend if multiple groups
- [ ] Statistical significance indicators
- [ ] Consistent styling across figures
- [ ] Vector format (PDF/SVG) for publications
