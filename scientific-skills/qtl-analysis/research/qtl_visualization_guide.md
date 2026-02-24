<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Comprehensive QTL and Genomic Data Visualization Guide

## Table of Contents
1. [Python Packages](#python-packages)
2. [R Packages](#r-packages)
3. [Basic Phenotype Visualizations](#basic-phenotype-visualizations)
4. [Advanced Visualizations](#advanced-visualizations)
5. [QTL-Specific Visualizations](#qtl-specific-visualizations)

---

## Python Packages

### 1. Core Packages Installation

```bash
# Essential Python packages for QTL visualization
pip install matplotlib seaborn plotly pandas numpy scipy

# Specialized genomic visualization
pip install qmplot gwaslab bioinfokit

# Interactive visualization
pip install plotly dash
```

### 2. Package Overview

| Package | Purpose | Key Functions |
|---------|---------|---------------|
| matplotlib | Basic plotting | plt.plot(), plt.scatter(), plt.boxplot() |
| seaborn | Statistical visualization | sns.boxplot(), sns.violinplot(), sns.heatmap() |
| plotly | Interactive plots | px.scatter(), px.line(), go.Figure() |
| pandas | Data manipulation | pd.DataFrame(), pd.read_csv(), .corr() |
| qmplot | GWAS plots | manhattan(), qqplot() |

---

## R Packages

### 1. Core Packages Installation

```r
# From CRAN
install.packages(c("ggplot2", "lattice", "corrplot", "RColorBrewer", 
                   "qqman", "dplyr", "tidyr", "gridExtra"))

# From Bioconductor (for genomic data)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("qtl2", "ggmanh", "ComplexHeatmap"))

# Development version of qtl2
install.packages("remotes")
remotes::install_github("rqtl/qtl2")
```

### 2. Package Overview

| Package | Purpose | Key Functions |
|---------|---------|---------------|
| ggplot2 | Primary visualization | ggplot(), geom_point(), geom_boxplot() |
| lattice | Trellis graphics | xyplot(), bwplot(), levelplot() |
| corrplot | Correlation matrices | corrplot(), corrplot.mixed() |
| RColorBrewer | Color palettes | brewer.pal(), display.brewer.all() |
| qqman | GWAS plots | manhattan(), qq() |
| qtl2 | QTL analysis | plot_scan1(), plot_coef(), plot_peaks() |

---

## Basic Phenotype Visualizations

### Python Examples

#### 1. Box Plots - Phenotype Distributions by Group

```python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Create sample QTL phenotype data
np.random.seed(42)
n_samples = 200

phenotype_data = pd.DataFrame({
    'Genotype': np.random.choice(['AA', 'AB', 'BB'], n_samples),
    'Environment': np.random.choice(['Control', 'Stress', 'Optimal'], n_samples),
    'Plant_Height': np.concatenate([
        np.random.normal(100, 15, 70),  # AA
        np.random.normal(120, 18, 65),    # AB
        np.random.normal(95, 12, 65)      # BB
    ]),
    'Yield': np.concatenate([
        np.random.normal(50, 10, 70),     # AA
        np.random.normal(65, 12, 65),     # AB
        np.random.normal(45, 8, 65)       # BB
    ]),
    'Flowering_Time': np.concatenate([
        np.random.normal(45, 5, 70),      # AA
        np.random.normal(42, 4, 65),      # AB
        np.random.normal(48, 6, 65)       # BB
    ])
})

# Basic box plot with matplotlib
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Box plot by genotype
axes[0].boxplot([phenotype_data[phenotype_data['Genotype'] == g]['Plant_Height'].values 
                 for g in ['AA', 'AB', 'BB']], 
                labels=['AA', 'AB', 'BB'])
axes[0].set_xlabel('Genotype')
axes[0].set_ylabel('Plant Height (cm)')
axes[0].set_title('Plant Height by Genotype')
axes[0].grid(True, alpha=0.3)

# Enhanced box plot with seaborn
sns.boxplot(data=phenotype_data, x='Genotype', y='Yield', 
            hue='Environment', palette='Set2', ax=axes[1])
axes[1].set_title('Yield Distribution by Genotype and Environment')
axes[1].legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig('boxplots_phenotype.png', dpi=300, bbox_inches='tight')
plt.show()

# Statistical summary
print("\nDescriptive Statistics by Genotype:")
print(phenotype_data.groupby('Genotype')[['Plant_Height', 'Yield', 'Flowering_Time']].describe())
```

#### 2. Density Plots - Phenotype Distributions

```python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats

# Create sample data
np.random.seed(42)
phenotype_data = pd.DataFrame({
    'Genotype': np.random.choice(['AA', 'AB', 'BB'], 300),
    'Trait_Value': np.concatenate([
        np.random.normal(100, 15, 100),  # AA
        np.random.normal(120, 20, 100),  # AB
        np.random.normal(90, 12, 100)    # BB
    ])
})

# Create figure with multiple density plots
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 1. Basic density plot with matplotlib
for genotype in ['AA', 'AB', 'BB']:
    data = phenotype_data[phenotype_data['Genotype'] == genotype]['Trait_Value']
    axes[0, 0].hist(data, bins=20, alpha=0.5, density=True, label=f'{genotype} (n={len(data)})')
    
    # Add KDE
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(data)
    x_range = np.linspace(data.min(), data.max(), 100)
    axes[0, 0].plot(x_range, kde(x_range), linewidth=2, label=f'{genotype} KDE')

axes[0, 0].set_xlabel('Trait Value')
axes[0, 0].set_ylabel('Density')
axes[0, 0].set_title('Density Plot with Histogram (Matplotlib)')
axes[0, 0].legend()
axes[0, 0].grid(True, alpha=0.3)

# 2. Seaborn KDE plot
sns.kdeplot(data=phenotype_data, x='Trait_Value', hue='Genotype', 
            fill=True, alpha=0.5, palette='viridis', ax=axes[0, 1])
axes[0, 1].set_title('KDE Plot by Genotype (Seaborn)')
axes[0, 1].grid(True, alpha=0.3)

# 3. Stacked density with rug plot
sns.kdeplot(data=phenotype_data[phenotype_data['Genotype'] == 'AA'], 
            x='Trait_Value', fill=True, alpha=0.3, color='blue', ax=axes[1, 0])
sns.kdeplot(data=phenotype_data[phenotype_data['Genotype'] == 'AB'], 
            x='Trait_Value', fill=True, alpha=0.3, color='green', ax=axes[1, 0])
sns.kdeplot(data=phenotype_data[phenotype_data['Genotype'] == 'BB'], 
            x='Trait_Value', fill=True, alpha=0.3, color='red', ax=axes[1, 0])
sns.rugplot(data=phenotype_data, x='Trait_Value', hue='Genotype', 
            palette=['blue', 'green', 'red'], ax=axes[1, 0])
axes[1, 0].set_title('Stacked Density with Rug Plot')
axes[1, 0].grid(True, alpha=0.3)

# 4. Distribution with normal fit
from scipy.stats import norm
data = phenotype_data['Trait_Value']
axes[1, 1].hist(data, bins=25, density=True, alpha=0.7, color='skyblue', edgecolor='black')
mu, std = norm.fit(data)
xmin, xmax = axes[1, 1].get_xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
axes[1, 1].plot(x, p, 'r-', linewidth=2, label=f'Normal fit: μ={mu:.1f}, σ={std:.1f}')
axes[1, 1].set_xlabel('Trait Value')
axes[1, 1].set_ylabel('Density')
axes[1, 1].set_title('Distribution with Normal Fit')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('density_plots.png', dpi=300, bbox_inches='tight')
plt.show()

# Normality test
from scipy.stats import shapiro, normaltest
stat, p_value = shapiro(phenotype_data['Trait_Value'].sample(100))
print(f"\nShapiro-Wilk Test for Normality:")
print(f"Statistic: {stat:.4f}, p-value: {p_value:.4f}")
print(f"Data is {'normally' if p_value > 0.05 else 'not normally'} distributed (α=0.05)")
```

#### 3. Violin Plots - Distributions with Density

```python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Create multi-environment phenotype data
np.random.seed(42)
n = 150
phenotype_data = pd.DataFrame({
    'Genotype': np.random.choice(['AA', 'AB', 'BB'], n),
    'Environment': np.random.choice(['Control', 'Drought', 'Heat'], n),
    'Biomass': np.concatenate([
        np.random.normal(50, 8, 50),   # AA
        np.random.normal(65, 10, 50),  # AB
        np.random.normal(45, 7, 50)    # BB
    ]),
    'Root_Length': np.concatenate([
        np.random.normal(25, 4, 50),   # AA
        np.random.normal(30, 5, 50),   # AB
        np.random.normal(22, 3, 50)    # BB
    ])
})

# Create comprehensive violin plot figure
fig = plt.figure(figsize=(16, 12))

# 1. Basic violin plot
ax1 = plt.subplot(2, 3, 1)
sns.violinplot(data=phenotype_data, x='Genotype', y='Biomass', palette='Set2', ax=ax1)
ax1.set_title('Basic Violin Plot')
ax1.grid(True, alpha=0.3)

# 2. Violin with swarm overlay
ax2 = plt.subplot(2, 3, 2)
sns.violinplot(data=phenotype_data, x='Genotype', y='Biomass', 
               inner=None, palette='pastel', ax=ax2)
sns.swarmplot(data=phenotype_data, x='Genotype', y='Biomass', 
              color='black', alpha=0.5, size=3, ax=ax2)
ax2.set_title('Violin with Data Points')
ax2.grid(True, alpha=0.3)

# 3. Split violin by environment
ax3 = plt.subplot(2, 3, 3)
sns.violinplot(data=phenotype_data, x='Genotype', y='Biomass', 
               hue='Environment', split=True, palette='muted', ax=ax3)
ax3.set_title('Split Violin by Environment')
ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax3.grid(True, alpha=0.3)

# 4. Box inside violin
ax4 = plt.subplot(2, 3, 4)
sns.violinplot(data=phenotype_data, x='Genotype', y='Root_Length', 
               inner='box', palette='coolwarm', ax=ax4)
ax4.set_title('Violin with Box Plot Inside')
ax4.grid(True, alpha=0.3)

# 5. Multiple traits comparison
ax5 = plt.subplot(2, 3, 5)
phenotype_melted = pd.melt(phenotype_data, id_vars=['Genotype'], 
                             value_vars=['Biomass', 'Root_Length'],
                             var_name='Trait', value_name='Value')
sns.violinplot(data=phenotype_melted, x='Genotype', y='Value', 
               hue='Trait', palette='Set1', ax=ax5)
ax5.set_title('Multiple Traits Comparison')
ax5.legend(title='Trait')
ax5.grid(True, alpha=0.3)

# 6. Customized violin with statistics
ax6 = plt.subplot(2, 3, 6)
parts = ax6.violinplot([phenotype_data[phenotype_data['Genotype'] == g]['Biomass'].values 
                        for g in ['AA', 'AB', 'BB']], 
                       positions=[1, 2, 3], showmeans=True, showmedians=True)

# Customize colors
colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']
for pc, color in zip(parts['bodies'], colors):
    pc.set_facecolor(color)
    pc.set_alpha(0.7)

ax6.set_xticks([1, 2, 3])
ax6.set_xticklabels(['AA', 'AB', 'BB'])
ax6.set_xlabel('Genotype')
ax6.set_ylabel('Biomass')
ax6.set_title('Customized Violin (Matplotlib)')
ax6.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('violin_plots.png', dpi=300, bbox_inches='tight')
plt.show()

# Statistical summary
print("\nViolin Plot Statistics:")
print(phenotype_data.groupby('Genotype')[['Biomass', 'Root_Length']].agg(['mean', 'median', 'std']))
```

---

## Advanced Visualizations

### Python Examples

#### 4. Heatmaps - Genotype Data and Phenotype Correlations

```python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Create sample correlation data
np.random.seed(42)
n_markers = 20
n_samples = 50

# Genotype matrix (0, 1, 2 for AA, AB, BB)
genotype_matrix = np.random.choice([0, 1, 2], size=(n_samples, n_markers))
marker_names = [f'Marker_{i+1}' for i in range(n_markers)]
sample_names = [f'Sample_{i+1}' for i in range(n_samples)]

genotype_df = pd.DataFrame(genotype_matrix, 
                          index=sample_names, 
                          columns=marker_names)

# Phenotype data
phenotype_df = pd.DataFrame({
    'Height': np.random.normal(100, 15, n_samples),
    'Weight': np.random.normal(50, 10, n_samples),
    'Yield': np.random.normal(75, 12, n_samples),
    'Quality': np.random.normal(85, 8, n_samples),
    'Disease_Resistance': np.random.normal(60, 20, n_samples)
}, index=sample_names)

# Calculate phenotype correlation matrix
pheno_corr = phenotype_df.corr()

# Create comprehensive heatmap figure
fig = plt.figure(figsize=(18, 14))

# 1. Genotype heatmap
ax1 = plt.subplot(2, 3, 1)
sns.heatmap(genotype_df.iloc[:20, :15], cmap='RdYlBu_r', 
            cbar_kws={'label': 'Genotype (0=AA, 1=AB, 2=BB)'}, 
            xticklabels=True, yticklabels=True, ax=ax1)
ax1.set_title('Genotype Matrix (First 20 samples, 15 markers)')
ax1.set_xlabel('Markers')
ax1.set_ylabel('Samples')

# 2. Phenotype correlation heatmap
ax2 = plt.subplot(2, 3, 2)
mask = np.triu(np.ones_like(pheno_corr, dtype=bool))
sns.heatmap(pheno_corr, mask=mask, annot=True, fmt='.2f', 
            cmap='coolwarm', center=0, square=True,
            cbar_kws={'shrink': 0.8}, ax=ax2)
ax2.set_title('Phenotype Correlation Matrix')

# 3. Clustermap for phenotypes
ax3 = plt.subplot(2, 3, 3)
sns.clustermap(phenotype_df.T, cmap='viridis', 
               col_cluster=True, row_cluster=True,
               figsize=(6, 5), dendrogram_ratio=0.2)
plt.title('Phenotype Clustermap')

# 4. Diverging heatmap for effect sizes
effect_sizes = np.random.randn(10, 8)
ax4 = plt.subplot(2, 3, 4)
sns.heatmap(effect_sizes, cmap='RdBu_r', center=0, 
            annot=True, fmt='.2f', square=True,
            cbar_kws={'label': 'Effect Size'}, ax=ax4)
ax4.set_title('QTL Effect Sizes Heatmap')
ax4.set_xlabel('Environments')
ax4.set_ylabel('Traits')

# 5. Annotated heatmap with custom labels
ax5 = plt.subplot(2, 3, 5)
annotation_data = np.random.rand(8, 6)
row_labels = ['Trait_' + str(i) for i in range(1, 9)]
col_labels = ['Env_' + str(i) for i in range(1, 7)]
sns.heatmap(annotation_data, annot=True, fmt='.2f', 
            cmap='YlOrRd', xticklabels=col_labels,
            yticklabels=row_labels, ax=ax5)
ax5.set_title('Multi-Environment Trial Results')

# 6. Correlation with dendrogram
ax6 = plt.subplot(2, 3, 6)
linkage = sns.clustermap(pheno_corr, cmap='coolwarm', 
                         figsize=(6, 5), dendrogram_ratio=0.2)
plt.title('Correlation with Hierarchical Clustering')

plt.tight_layout()
plt.savefig('heatmaps.png', dpi=300, bbox_inches='tight')
plt.show()

# Print correlation insights
print("\nStrongest Phenotype Correlations:")
corr_pairs = []
for i in range(len(pheno_corr.columns)):
    for j in range(i+1, len(pheno_corr.columns)):
        corr_pairs.append((pheno_corr.columns[i], pheno_corr.columns[j], 
                          pheno_corr.iloc[i, j]))
corr_pairs.sort(key=lambda x: abs(x[2]), reverse=True)
for trait1, trait2, corr in corr_pairs[:3]:
    print(f"  {trait1} - {trait2}: {corr:.3f}")
```

#### 5. Scatter Plot Matrices - Multiple Trait Correlations

```python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Create multi-trait phenotype data
np.random.seed(42)
n = 200

# Generate correlated traits
mean = [100, 50, 75, 85, 60]
cov = [[225, 80, 120, 60, -40],
       [80, 100, 50, 30, -20],
       [120, 50, 144, 70, -30],
       [60, 30, 70, 64, -15],
       [-40, -20, -30, -15, 400]]

traits = np.random.multivariate_normal(mean, cov, n)
phenotype_df = pd.DataFrame(traits, 
                           columns=['Height', 'Weight', 'Yield', 'Quality', 'Disease_Res'])
phenotype_df['Genotype'] = np.random.choice(['AA', 'AB', 'BB'], n)
phenotype_df['Environment'] = np.random.choice(['Control', 'Stress'], n)

# Create scatter plot matrices
fig = plt.figure(figsize=(18, 14))

# 1. Basic pairplot with seaborn
ax1 = plt.subplot(2, 2, 1)
g = sns.pairplot(phenotype_df, hue='Genotype', palette='Set1', 
                 diag_kind='kde', height=2.5)
g.fig.suptitle('Pairplot by Genotype', y=1.02)

# 2. Custom pairplot with regression
ax2 = plt.subplot(2, 2, 2)
g2 = sns.pairplot(phenotype_df, vars=['Height', 'Yield', 'Quality'], 
                  hue='Environment', palette='coolwarm',
                  kind='reg', diag_kind='hist', height=2.5)
g2.fig.suptitle('Pairplot with Regression Lines', y=1.02)

# 3. Correlation matrix scatter
ax3 = plt.subplot(2, 2, 3)
from pandas.plotting import scatter_matrix
scatter_matrix(phenotype_df[['Height', 'Weight', 'Yield', 'Quality']], 
               alpha=0.5, figsize=(6, 6), diagonal='kde', ax=ax3)
ax3.set_title('Scatter Matrix (Pandas)')

# 4. Joint plots for key relationships
fig2, axes = plt.subplots(2, 2, figsize=(12, 12))

# Height vs Yield
sns.scatterplot(data=phenotype_df, x='Height', y='Yield', 
                hue='Genotype', palette='Set1', alpha=0.6, ax=axes[0, 0])
sns.regplot(data=phenotype_df, x='Height', y='Yield', 
            scatter=False, color='red', ax=axes[0, 0])
axes[0, 0].set_title('Height vs Yield by Genotype')

# Weight vs Quality
sns.scatterplot(data=phenotype_df, x='Weight', y='Quality', 
                hue='Environment', palette='coolwarm', alpha=0.6, ax=axes[0, 1])
sns.regplot(data=phenotype_df, x='Weight', y='Quality', 
            scatter=False, color='green', ax=axes[0, 1])
axes[0, 1].set_title('Weight vs Quality by Environment')

# Yield vs Disease Resistance
sns.scatterplot(data=phenotype_df, x='Yield', y='Disease_Res', 
                hue='Genotype', style='Environment', palette='Set2', ax=axes[1, 0])
axes[1, 0].set_title('Yield vs Disease Resistance')

# 3D scatter (simulated with color)
scatter = axes[1, 1].scatter(phenotype_df['Height'], phenotype_df['Yield'], 
                               c=phenotype_df['Quality'], cmap='viridis', 
                               s=50, alpha=0.6)
axes[1, 1].set_xlabel('Height')
axes[1, 1].set_ylabel('Yield')
axes[1, 1].set_title('Height vs Yield (colored by Quality)')
plt.colorbar(scatter, ax=axes[1, 1], label='Quality')

plt.tight_layout()
plt.savefig('scatter_matrices.png', dpi=300, bbox_inches='tight')
plt.show()

# Calculate and display correlations
print("\nTrait Correlations:")
correlations = phenotype_df[['Height', 'Weight', 'Yield', 'Quality', 'Disease_Res']].corr()
print(correlations.round(3))
```

---

## QTL-Specific Visualizations

### Python Examples

#### 6. Manhattan Plots - GWAS/QTL Results

```python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Try to import qmplot, fall back to manual implementation
try:
    from qmplot import manhattanplot
    qmplot_available = True
except ImportError:
    qmplot_available = False
    print("qmplot not available, using manual implementation")

# Create simulated GWAS results
np.random.seed(42)

# Simulate 5 chromosomes with varying marker densities
chromosomes = []
positions = []
p_values = []
snp_ids = []

snp_counter = 1
for chrom in range(1, 6):
    n_markers = np.random.randint(500, 1000)
    chromosomes.extend([chrom] * n_markers)
    positions.extend(np.sort(np.random.randint(0, 100000000, n_markers)))
    
    # Generate p-values with some significant hits
    if chrom == 3:  # Simulate QTL on chromosome 3
        p = np.random.uniform(0.0001, 0.05, n_markers)
        p[400:420] = np.random.uniform(1e-10, 1e-6, 20)  # Peak
    elif chrom == 5:  # Another QTL
        p = np.random.uniform(0.001, 0.1, n_markers)
        p[600:630] = np.random.uniform(1e-8, 1e-5, 30)
    else:
        p = np.random.uniform(0.1, 1.0, n_markers)
    
    p_values.extend(p)
    snp_ids.extend([f'rs{snp_counter + i}' for i in range(n_markers)])
    snp_counter += n_markers

gwas_results = pd.DataFrame({
    'SNP': snp_ids,
    'CHR': chromosomes,
    'BP': positions,
    'P': p_values
})

# Create Manhattan plot
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

if qmplot_available:
    # Using qmplot
    ax1 = axes[0, 0]
    manhattanplot(data=gwas_results, 
                  chrom="CHR", 
                  pos="BP", 
                  pvalue="P",
                  ax=ax1,
                  title="Manhattan Plot (qmplot)")
    ax1.axhline(y=-np.log10(5e-8), color='red', linestyle='--', label='Genome-wide significance')
    ax1.legend()
else:
    # Manual Manhattan plot implementation
    def plot_manhattan(df, ax, title="Manhattan Plot"):
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
        x_pos = []
        y_vals = []
        colors_list = []
        
        offset = 0
        for i, chrom in enumerate(sorted(df['CHR'].unique())):
            chrom_data = df[df['CHR'] == chrom]
            x_pos.extend(chrom_data['BP'] / 1e6 + offset)
            y_vals.extend(-np.log10(chrom_data['P']))
            colors_list.extend([colors[i % len(colors)]] * len(chrom_data))
            offset += chrom_data['BP'].max() / 1e6 + 10
        
        ax.scatter(x_pos, y_vals, c=colors_list, s=10, alpha=0.6)
        ax.axhline(y=-np.log10(5e-8), color='red', linestyle='--', 
                   label='Genome-wide significance (5e-8)')
        ax.axhline(y=-np.log10(1e-5), color='orange', linestyle='--', 
                   label='Suggestive (1e-5)')
        ax.set_xlabel('Chromosome Position (Mb)')
        ax.set_ylabel('-log10(p-value)')
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plot_manhattan(gwas_results, axes[0, 0], "Manhattan Plot (Manual)")

# Enhanced Manhattan with annotations
ax2 = axes[0, 1]
# Highlight significant SNPs
significant = gwas_results[gwas_results['P'] < 1e-5]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

offset = 0
chrom_offsets = {}
for i, chrom in enumerate(sorted(gwas_results['CHR'].unique())):
    chrom_data = gwas_results[gwas_results['CHR'] == chrom]
    x = chrom_data['BP'] / 1e6 + offset
    y = -np.log10(chrom_data['P'])
    
    ax2.scatter(x, y, c=colors[i % len(colors)], s=10, alpha=0.6)
    
    # Highlight significant
    sig_chrom = significant[significant['CHR'] == chrom]
    if len(sig_chrom) > 0:
        sig_x = sig_chrom['BP'] / 1e6 + offset
        sig_y = -np.log10(sig_chrom['P'])
        ax2.scatter(sig_x, sig_y, c='red', s=50, alpha=0.8, edgecolors='black')
    
    chrom_offsets[chrom] = offset
    offset += chrom_data['BP'].max() / 1e6 + 10

ax2.axhline(y=-np.log10(5e-8), color='red', linestyle='--', linewidth=2)
ax2.axhline(y=-np.log10(1e-5), color='orange', linestyle='--', linewidth=2)
ax2.set_xlabel('Genomic Position')
ax2.set_ylabel('-log10(p-value)')
ax2.set_title('Manhattan Plot with Significant Hits Highlighted')
ax2.grid(True, alpha=0.3)

# Chromosome-specific plot
ax3 = axes[1, 0]
chrom3_data = gwas_results[gwas_results['CHR'] == 3]
ax3.scatter(chrom3_data['BP'] / 1e6, -np.log10(chrom3_data['P']), 
          c='blue', s=20, alpha=0.6)
ax3.axhline(y=-np.log10(5e-8), color='red', linestyle='--')
ax3.set_xlabel('Position on Chromosome 3 (Mb)')
ax3.set_ylabel('-log10(p-value)')
ax3.set_title('Chromosome 3 Detail')
ax3.grid(True, alpha=0.3)

# Summary statistics
ax4 = axes[1, 1]
ax4.axis('off')
summary_text = f"""
GWAS Summary Statistics:

Total SNPs: {len(gwas_results):,}
Chromosomes: {gwas_results['CHR'].nunique()}
Genome-wide significant (p < 5e-8): {len(gwas_results[gwas_results['P'] < 5e-8]):,}
Suggestive (p < 1e-5): {len(gwas_results[gwas_results['P'] < 1e-5]):,}

Top Hits:
"""
top_hits = gwas_results.nsmallest(5, 'P')[['SNP', 'CHR', 'BP', 'P']]
for _, row in top_hits.iterrows():
    summary_text += f"  {row['SNP']} (Chr{row['CHR']}:{row['BP']:,}): p={row['P']:.2e}\n"

ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, 
         fontsize=10, verticalalignment='top', fontfamily='monospace')

plt.tight_layout()
plt.savefig('manhattan_plots.png', dpi=300, bbox_inches='tight')
plt.show()

print("\nManhattan plot created successfully!")
print(f"Total markers: {len(gwas_results)}")
print(f"Significant hits (p < 1e-5): {len(significant)}")
```

#### 7. QQ Plots - Assessing Normality

```python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Create sample data for QQ plots
np.random.seed(42)

# 1. Normal distribution
normal_data = np.random.normal(0, 1, 1000)

# 2. Phenotype data (may not be normal)
phenotype_data = np.concatenate([
    np.random.normal(100, 15, 800),
    np.random.normal(130, 20, 200)  # Some outliers
])

# 3. GWAS p-values (should be uniform under null)
gwas_pvalues = np.random.uniform(0, 1, 10000)
# Add some significant hits
gwas_pvalues[:50] = np.random.uniform(0, 1e-6, 50)

# Create QQ plots
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# 1. QQ plot for normal data
ax1 = axes[0, 0]
stats.probplot(normal_data, dist="norm", plot=ax1)
ax1.set_title('QQ Plot: Normal Distribution\n(Should follow line)')
ax1.grid(True, alpha=0.3)

# 2. QQ plot for phenotype
ax2 = axes[0, 1]
stats.probplot(phenotype_data, dist="norm", plot=ax2)
ax2.set_title('QQ Plot: Phenotype Data\n(Deviations indicate non-normality)')
ax2.grid(True, alpha=0.3)

# 3. GWAS QQ plot
ax3 = axes[0, 2]
# Sort p-values
sorted_p = np.sort(gwas_pvalues)
expected = np.arange(1, len(sorted_p) + 1) / len(sorted_p)

# -log10 transformation
observed_log = -np.log10(sorted_p)
expected_log = -np.log10(expected)

ax3.scatter(expected_log, observed_log, s=10, alpha=0.5)
ax3.plot([0, expected_log.max()], [0, expected_log.max()], 
         'r--', label='Expected (null)')
ax3.set_xlabel('Expected -log10(p)')
ax3.set_ylabel('Observed -log10(p)')
ax3.set_title('QQ Plot: GWAS P-values')
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4. Enhanced GWAS QQ with confidence intervals
ax4 = axes[1, 0]
n = len(gwas_pvalues)
sorted_p = np.sort(gwas_pvalues)
expected = np.arange(1, n + 1) / n

observed_log = -np.log10(sorted_p)
expected_log = -np.log10(expected)

# 95% confidence interval
# Under null, p-values are uniform, so -log10(p) follows exponential
ci_upper = -np.log10(stats.beta.ppf(0.975, np.arange(1, n+1), n - np.arange(1, n+1) + 1))
ci_lower = -np.log10(stats.beta.ppf(0.025, np.arange(1, n+1), n - np.arange(1, n+1) + 1))

ax4.fill_between(expected_log, ci_lower, ci_upper, alpha=0.2, color='gray', label='95% CI')
ax4.scatter(expected_log, observed_log, s=10, alpha=0.5, c='blue')
ax4.plot([0, expected_log.max()], [0, expected_log.max()], 'r--', label='Expected')
ax4.set_xlabel('Expected -log10(p)')
ax4.set_ylabel('Observed -log10(p)')
ax4.set_title('GWAS QQ Plot with 95% CI')
ax4.legend()
ax4.grid(True, alpha=0.3)

# 5. Lambda GC calculation
ax5 = axes[1, 1]
# Lambda GC is median(observed chi-square) / median(expected chi-square)
# For p-values, convert to chi-square: qchisq(1-p, 1)
from scipy.stats import chi2
chi2_obs = chi2.ppf(1 - sorted_p, 1)
chi2_exp = chi2.ppf(1 - expected, 1)
lambda_gc = np.median(chi2_obs) / np.median(chi2_exp)

ax5.scatter(expected_log, observed_log, s=10, alpha=0.5)
ax5.plot([0, expected_log.max()], [0, expected_log.max()], 'r--')
ax5.set_xlabel('Expected -log10(p)')
ax5.set_ylabel('Observed -log10(p)')
ax5.set_title(f'QQ Plot: λGC = {lambda_gc:.3f}')
ax5.text(0.5, 0.95, f'Lambda GC: {lambda_gc:.3f}\n(1.0 = no inflation)', 
         transform=ax5.transAxes, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
ax5.grid(True, alpha=0.3)

# 6. Multiple traits QQ comparison
ax6 = axes[1, 2]
traits = ['Trait1', 'Trait2', 'Trait3']
colors = ['blue', 'green', 'red']

for trait, color in zip(traits, colors):
    if trait == 'Trait1':
        pvals = np.random.uniform(0, 1, 1000)
    elif trait == 'Trait2':
        pvals = np.concatenate([np.random.uniform(0, 1, 950), 
                               np.random.uniform(0, 1e-4, 50)])
    else:
        pvals = np.concatenate([np.random.uniform(0, 1, 900), 
                               np.random.uniform(0, 1e-6, 100)])
    
    sorted_p = np.sort(pvals)
    expected = np.arange(1, len(sorted_p) + 1) / len(sorted_p)
    ax6.scatter(-np.log10(expected), -np.log10(sorted_p), 
               s=10, alpha=0.5, label=trait, c=color)

ax6.plot([0, 4], [0, 4], 'k--', label='Expected')
ax6.set_xlabel('Expected -log10(p)')
ax6.set_ylabel('Observed -log10(p)')
ax6.set_title('QQ Plot: Multiple Traits Comparison')
ax6.legend()
ax6.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('qq_plots.png', dpi=300, bbox_inches='tight')
plt.show()

print("\nQQ Plots Summary:")
print(f"Lambda GC (genomic control): {lambda_gc:.3f}")
print(f"  - λGC ≈ 1.0: No population stratification")
print(f"  - λGC > 1.05: Possible population stratification")
print(f"  - λGC < 1.0: Conservative test")
```

#### 8. LOD Score Plots

```python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Create simulated QTL scan results
np.random.seed(42)

# Simulate LOD scores across chromosomes
chromosomes = []
positions = []
lod_scores = []

for chrom in range(1, 6):
    # Each chromosome has different length and marker density
    chrom_length = np.random.randint(80, 150)  # cM
    n_markers = np.random.randint(50, 100)
    
    pos = np.linspace(0, chrom_length, n_markers)
    
    # Generate LOD scores with some peaks
    if chrom == 3:
        # Strong QTL on chromosome 3
        lod = np.random.exponential(2, n_markers)
        peak_pos = 45
        lod += 15 * np.exp(-((pos - peak_pos) / 10) ** 2)
    elif chrom == 5:
        # Moderate QTL on chromosome 5
        lod = np.random.exponential(1.5, n_markers)
        peak_pos = 60
        lod += 8 * np.exp(-((pos - peak_pos) / 8) ** 2)
    else:
        lod = np.random.exponential(1, n_markers)
    
    chromosomes.extend([chrom] * n_markers)
    positions.extend(pos)
    lod_scores.extend(lod)

qtl_results = pd.DataFrame({
    'Chromosome': chromosomes,
    'Position': positions,
    'LOD': lod_scores
})

# Create LOD plots
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# 1. Genome-wide LOD scan
ax1 = axes[0, 0]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
offset = 0
xticks = []
xtick_labels = []

for i, chrom in enumerate(sorted(qtl_results['Chromosome'].unique())):
    chrom_data = qtl_results[qtl_results['Chromosome'] == chrom]
    x = chrom_data['Position'] + offset
    ax1.plot(x, chrom_data['LOD'], color=colors[i % len(colors)], linewidth=1.5)
    
    xticks.append(offset + chrom_data['Position'].mean())
    xtick_labels.append(str(chrom))
    offset += chrom_data['Position'].max() + 10

ax1.axhline(y=3, color='red', linestyle='--', label='Significance threshold (LOD=3)')
ax1.axhline(y=4.5, color='orange', linestyle='--', label='Highly significant (LOD=4.5)')
ax1.set_xticks(xticks)
ax1.set_xticklabels(xtick_labels)
ax1.set_xlabel('Chromosome')
ax1.set_ylabel('LOD Score')
ax1.set_title('Genome-wide QTL Scan')
ax1.legend()
ax1.grid(True, alpha=0.3)

# 2. Individual chromosome plots
ax2 = axes[0, 1]
chrom3_data = qtl_results[qtl_results['Chromosome'] == 3]
ax2.plot(chrom3_data['Position'], chrom3_data['LOD'], 'b-', linewidth=2)
ax2.fill_between(chrom3_data['Position'], 0, chrom3_data['LOD'], alpha=0.3)
ax2.axhline(y=3, color='red', linestyle='--')
ax2.set_xlabel('Position (cM)')
ax2.set_ylabel('LOD Score')
ax2.set_title('Chromosome 3 QTL Profile')
ax2.grid(True, alpha=0.3)

# Find and annotate peak
peak_idx = chrom3_data['LOD'].idxmax()
peak_pos = chrom3_data.loc[peak_idx, 'Position']
peak_lod = chrom3_data.loc[peak_idx, 'LOD']
ax2.annotate(f'Peak: {peak_lod:.1f}\n@{peak_pos:.1f}cM', 
             xy=(peak_pos, peak_lod), xytext=(peak_pos + 10, peak_lod - 2),
             arrowprops=dict(arrowstyle='->', color='red'),
             bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

# 3. Multiple traits comparison
ax3 = axes[1, 0]
chrom = 3
chrom_data = qtl_results[qtl_results['Chromosome'] == chrom]

# Simulate multiple traits
traits = ['Height', 'Yield', 'Quality']
colors_traits = ['blue', 'green', 'red']
for trait, color in zip(traits, colors_traits):
    if trait == 'Height':
        lod = chrom_data['LOD'].values
    elif trait == 'Yield':
        lod = chrom_data['LOD'].values * 0.7 + np.random.normal(0, 0.5, len(chrom_data))
    else:
        lod = chrom_data['LOD'].values * 0.4 + np.random.normal(0, 0.3, len(chrom_data))
    
    ax3.plot(chrom_data['Position'], lod, color=color, linewidth=2, 
             label=trait, alpha=0.7)

ax3.axhline(y=3, color='black', linestyle='--', alpha=0.5)
ax3.set_xlabel('Position (cM)')
ax3.set_ylabel('LOD Score')
ax3.set_title(f'Chromosome {chrom}: Multiple Traits')
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4. 2D scan (epistasis)
ax4 = axes[1, 1]
# Simulate 2D LOD scan
pos_range = np.linspace(0, 100, 50)
lod_2d = np.random.exponential(1, (50, 50))

# Add interaction peak
for i, p1 in enumerate(pos_range):
    for j, p2 in enumerate(pos_range):
        lod_2d[i, j] += 5 * np.exp(-((p1 - 40)**2 + (p2 - 60)**2) / 200)

im = ax4.contourf(pos_range, pos_range, lod_2d, levels=20, cmap='viridis')
plt.colorbar(im, ax=ax4, label='LOD Score')
ax4.set_xlabel('Chromosome 3 Position (cM)')
ax4.set_ylabel('Chromosome 5 Position (cM)')
ax4.set_title('2D QTL Scan (Epistasis)')
ax4.plot(40, 60, 'r*', markersize=15, label='Interaction peak')
ax4.legend()

plt.tight_layout()
plt.savefig('lod_score_plots.png', dpi=300, bbox_inches='tight')
plt.show()

# Summary
print("\nQTL Scan Summary:")
print(f"Chromosomes scanned: {qtl_results['Chromosome'].nunique()}")
print(f"Total markers: {len(qtl_results)}")
print(f"Significant peaks (LOD > 3): {len(qtl_results[qtl_results['LOD'] > 3])}")
print(f"\nTop QTL peaks:")
top_peaks = qtl_results.nlargest(5, 'LOD')[['Chromosome', 'Position', 'LOD']]
for _, row in top_peaks.iterrows():
    print(f"  Chr{row['Chromosome']} @ {row['Position']:.1f}cM: LOD = {row['LOD']:.2f}")
```

#### 9. Effect Size Plots

```python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Create effect size data
np.random.seed(42)

# Simulate QTL effects for multiple traits
n_qtl = 8
traits = ['Height', 'Yield', 'Quality', 'Disease_Res']

qtl_effects = []
for i in range(1, n_qtl + 1):
    for trait in traits:
        effect = np.random.normal(0, 2)
        se = np.random.uniform(0.5, 1.5)
        qtl_effects.append({
            'QTL': f'QTL_{i}',
            'Trait': trait,
            'Effect': effect,
            'SE': se,
            'Chromosome': np.random.randint(1, 6),
            'Position': np.random.uniform(0, 100)
        })

effects_df = pd.DataFrame(qtl_effects)

# Create effect size plots
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# 1. Forest plot of effects
ax1 = axes[0, 0]
qtl_names = effects_df['QTL'].unique()
y_pos = np.arange(len(qtl_names))

for i, trait in enumerate(traits):
    trait_data = effects_df[effects_df['Trait'] == trait]
    effects = [trait_data[trait_data['QTL'] == q]['Effect'].values[0] 
               for q in qtl_names]
    ses = [trait_data[trait_data['QTL'] == q]['SE'].values[0] 
           for q in qtl_names]
    
    x_pos = y_pos + i * 0.2 - 0.3
    ax1.errorbar(effects, x_pos, xerr=ses, fmt='o', label=trait, capsize=3)

ax1.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
ax1.set_yticks(y_pos)
ax1.set_yticklabels(qtl_names)
ax1.set_xlabel('Effect Size')
ax1.set_title('QTL Effect Sizes by Trait (Forest Plot)')
ax1.legend()
ax1.grid(True, alpha=0.3)

# 2. Allele effect plot
ax2 = axes[0, 1]
genotypes = ['AA', 'AB', 'BB']
n_genotypes = len(genotypes)

# Simulate phenotype means by genotype
phenotype_means = {
    'Height': [100, 115, 95],
    'Yield': [50, 65, 45],
    'Quality': [80, 85, 78]
}

x = np.arange(n_genotypes)
width = 0.25
multiplier = 0

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
for i, (trait, means) in enumerate(phenotype_means.items()):
    offset = width * multiplier
    bars = ax2.bar(x + offset, means, width, label=trait, color=colors[i])
    ax2.bar_label(bars, padding=3, fmt='%.0f')
    multiplier += 1

ax2.set_xlabel('Genotype')
ax2.set_ylabel('Phenotype Value')
ax2.set_title('Allele Effects by Genotype')
ax2.set_xticks(x + width)
ax2.set_xticklabels(genotypes)
ax2.legend(loc='upper left')
ax2.grid(True, alpha=0.3, axis='y')

# 3. Effect size heatmap
ax3 = axes[1, 0]
effect_matrix = effects_df.pivot(index='QTL', columns='Trait', values='Effect')
sns.heatmap(effect_matrix, annot=True, fmt='.2f', cmap='RdBu_r', 
            center=0, ax=ax3, cbar_kws={'label': 'Effect Size'})
ax3.set_title('QTL Effect Size Heatmap')

# 4. Effect vs significance
ax4 = axes[1, 1]
effects_df['-log10p'] = -np.log10(np.random.uniform(0.001, 0.1, len(effects_df)))
effects_df['Significant'] = effects_df['-log10p'] > 2

colors = ['blue' if sig else 'gray' for sig in effects_df['Significant']]
ax4.scatter(effects_df['Effect'], effects_df['-log10p'], 
           c=colors, alpha=0.6, s=50)
ax4.axhline(y=2, color='red', linestyle='--', label='p=0.01')
ax4.set_xlabel('Effect Size')
ax4.set_ylabel('-log10(p-value)')
ax4.set_title('Effect Size vs Significance')
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('effect_size_plots.png', dpi=300, bbox_inches='tight')
plt.show()

print("\nEffect Size Summary:")
print(effects_df.groupby('Trait')[['Effect', 'SE']].mean())
```

#### 10. Multi-Environment Trial Visualization

```python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Create multi-environment trial data
np.random.seed(42)

genotypes = ['G1', 'G2', 'G3', 'G4', 'G5']
environments = ['Env1', 'Env2', 'Env3', 'Env4']

# Simulate GxE data
gxe_data = []
for geno in genotypes:
    for env in environments:
        # Different genotypes perform differently in different environments
        base_yield = {'G1': 50, 'G2': 55, 'G3': 45, 'G4': 60, 'G5': 52}[geno]
        env_effect = {'Env1': 0, 'Env2': 5, 'Env3': -3, 'Env4': 8}[env]
        
        # Add GxE interaction
        if geno == 'G4' and env == 'Env4':
            interaction = 10
        elif geno == 'G3' and env == 'Env3':
            interaction = -8
        else:
            interaction = np.random.normal(0, 3)
        
        yield_value = base_yield + env_effect + interaction + np.random.normal(0, 4)
        
        gxe_data.append({
            'Genotype': geno,
            'Environment': env,
            'Yield': yield_value,
            'Stability': np.random.uniform(0.5, 1.5)
        })

gxe_df = pd.DataFrame(gxe_data)

# Create multi-environment visualizations
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# 1. Interaction plot (line plot)
ax1 = axes[0, 0]
for geno in genotypes:
    geno_data = gxe_df[gxe_df['Genotype'] == geno]
    ax1.plot(geno_data['Environment'], geno_data['Yield'], 
             marker='o', label=geno, linewidth=2)

ax1.set_xlabel('Environment')
ax1.set_ylabel('Yield')
ax1.set_title('Genotype × Environment Interaction')
ax1.legend(title='Genotype', bbox_to_anchor=(1.05, 1), loc='upper left')
ax1.grid(True, alpha=0.3)

# 2. Heatmap of GxE
ax2 = axes[0, 1]
gxe_pivot = gxe_df.pivot(index='Genotype', columns='Environment', values='Yield')
sns.heatmap(gxe_pivot, annot=True, fmt='.1f', cmap='YlOrRd', 
            ax=ax2, cbar_kws={'label': 'Yield'})
ax2.set_title('GxE Heatmap')

# 3. Grouped bar chart
ax3 = axes[0, 2]
gxe_pivot.plot(kind='bar', ax=ax3, width=0.8)
ax3.set_xlabel('Genotype')
ax3.set_ylabel('Yield')
ax3.set_title('Yield by Genotype and Environment')
ax3.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')
ax3.tick_params(axis='x', rotation=0)
ax3.grid(True, alpha=0.3, axis='y')

# 4. AMMI-style biplot (simplified)
ax4 = axes[1, 0]
# Calculate means
geno_means = gxe_df.groupby('Genotype')['Yield'].mean()
env_means = gxe_df.groupby('Environment')['Yield'].mean()
overall_mean = gxe_df['Yield'].mean()

# Plot genotype means vs environment means
for geno in genotypes:
    for env in environments:
        value = gxe_df[(gxe_df['Genotype'] == geno) & 
                      (gxe_df['Environment'] == env)]['Yield'].values[0]
        ax4.scatter(env_means[env], value, s=100, alpha=0.6)

# Add regression line
from scipy import stats
slope, intercept, r_value, p_value, std_err = stats.linregress(
    [env_means[e] for e in environments for _ in genotypes],
    gxe_df['Yield']
)
x_line = np.linspace(env_means.min(), env_means.max(), 100)
ax4.plot(x_line, slope * x_line + intercept, 'r--', 
         label=f'Regression (R²={r_value**2:.3f})')

ax4.set_xlabel('Environment Mean')
ax4.set_ylabel('Genotype Performance')
ax4.set_title('AMMI-style Plot')
ax4.legend()
ax4.grid(True, alpha=0.3)

# 5. Stability analysis
ax5 = axes[1, 1]
geno_stats = gxe_df.groupby('Genotype').agg({
    'Yield': ['mean', 'std']
}).reset_index()
geno_stats.columns = ['Genotype', 'Mean_Yield', 'Std_Yield']

scatter = ax5.scatter(geno_stats['Mean_Yield'], geno_stats['Std_Yield'], 
                     s=200, alpha=0.6, c=range(len(genotypes)), cmap='viridis')
for i, geno in enumerate(geno_stats['Genotype']):
    ax5.annotate(geno, (geno_stats['Mean_Yield'].iloc[i], 
                       geno_stats['Std_Yield'].iloc[i]),
                xytext=(5, 5), textcoords='offset points')

ax5.set_xlabel('Mean Yield')
ax5.set_ylabel('Standard Deviation (Stability)')
ax5.set_title('Stability Analysis\n(Lower right = high yield, stable)')
ax5.grid(True, alpha=0.3)

# 6. Environment correlation
ax6 = axes[1, 2]
env_corr = gxe_pivot.T.corr()
sns.heatmap(env_corr, annot=True, fmt='.2f', cmap='coolwarm', 
            center=0, ax=ax6, square=True)
ax6.set_title('Genotype Correlation Across Environments')

plt.tight_layout()
plt.savefig('multienvironment_plots.png', dpi=300, bbox_inches='tight')
plt.show()

print("\nMulti-Environment Summary:")
print("\nMean Yield by Genotype:")
print(gxe_df.groupby('Genotype')['Yield'].mean().sort_values(ascending=False))
print("\nMean Yield by Environment:")
print(gxe_df.groupby('Environment')['Yield'].mean().sort_values(ascending=False))
```

---

## R Examples

### 1. Box Plots with ggplot2

```r
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Create sample data
set.seed(42)
n <- 200
phenotype_data <- data.frame(
  Genotype = sample(c("AA", "AB", "BB"), n, replace = TRUE),
  Environment = sample(c("Control", "Stress", "Optimal"), n, replace = TRUE),
  Plant_Height = c(rnorm(70, 100, 15), rnorm(65, 120, 18), rnorm(65, 95, 12)),
  Yield = c(rnorm(70, 50, 10), rnorm(65, 65, 12), rnorm(65, 45, 8))
)

# Basic box plot
p1 <- ggplot(phenotype_data, aes(x = Genotype, y = Plant_Height, fill = Genotype)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Plant Height by Genotype",
       x = "Genotype", y = "Plant Height (cm)")

# Grouped box plot
p2 <- ggplot(phenotype_data, aes(x = Genotype, y = Yield, fill = Environment)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Yield by Genotype and Environment",
       x = "Genotype", y = "Yield") +
  scale_fill_brewer(palette = "Set2")

# Notched box plot with jittered points
p3 <- ggplot(phenotype_data, aes(x = Genotype, y = Plant_Height, fill = Genotype)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  theme_minimal() +
  labs(title = "Plant Height with Notches and Individual Points")

# Faceted box plot
p4 <- ggplot(phenotype_data, aes(x = Environment, y = Yield, fill = Genotype)) +
  geom_boxplot() +
  facet_wrap(~Genotype) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Yield Distribution by Environment (Faceted)")

# Display plots
print(p1)
print(p2)
print(p3)
print(p4)
```

### 2. Density and Violin Plots with ggplot2

```r
library(ggplot2)
library(ggridges)

# Create sample data
set.seed(42)
phenotype_data <- data.frame(
  Genotype = rep(c("AA", "AB", "BB"), each = 100),
  Trait_Value = c(rnorm(100, 100, 15), rnorm(100, 120, 20), rnorm(100, 90, 12))
)

# Density plot
p1 <- ggplot(phenotype_data, aes(x = Trait_Value, fill = Genotype)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Trait Distribution by Genotype",
       x = "Trait Value", y = "Density") +
  scale_fill_brewer(palette = "Set1")

# Overlapping density with rug
p2 <- ggplot(phenotype_data, aes(x = Trait_Value, color = Genotype)) +
  geom_density(size = 1.2) +
  geom_rug(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Density with Rug Plot")

# Violin plot
p3 <- ggplot(phenotype_data, aes(x = Genotype, y = Trait_Value, fill = Genotype)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  theme_minimal() +
  labs(title = "Violin Plot with Box Plot Overlay")

# Ridgeline plot
p4 <- ggplot(phenotype_data, aes(x = Trait_Value, y = Genotype, fill = Genotype)) +
  geom_density_ridges(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Ridgeline Plot") +
  scale_fill_brewer(palette = "Set2")

print(p1)
print(p2)
print(p3)
print(p4)
```

### 3. Heatmaps with ggplot2 and corrplot

```r
library(ggplot2)
library(corrplot)
library(reshape2)

# Create correlation matrix
set.seed(42)
n <- 100
phenotype_df <- data.frame(
  Height = rnorm(n, 100, 15),
  Weight = rnorm(n, 50, 10),
  Yield = rnorm(n, 75, 12),
  Quality = rnorm(n, 85, 8),
  Disease_Res = rnorm(n, 60, 20)
)

# Add correlations
phenotype_df$Weight <- phenotype_df$Weight + 0.4 * phenotype_df$Height + rnorm(n, 0, 5)
phenotype_df$Yield <- phenotype_df$Yield + 0.5 * phenotype_df$Height + rnorm(n, 0, 8)

# Calculate correlation matrix
corr_matrix <- cor(phenotype_df)

# corrplot visualization
corrplot(corr_matrix, method = "color", type = "upper",
         addCoef.col = "black", tl.col = "black", tl.srt = 45,
         title = "Phenotype Correlation Matrix")

# Mixed correlation plot
corrplot.mixed(corr_matrix, lower = "number", upper = "circle",
               tl.col = "black", title = "Mixed Correlation Plot")

# ggplot2 heatmap
corr_melted <- melt(corr_matrix)
ggplot(corr_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Heatmap", x = "", y = "")
```

### 4. Scatter Plot Matrix with GGally

```r
library(GGally)
library(ggplot2)

# Create data
set.seed(42)
n <- 150
df <- data.frame(
  Height = rnorm(n, 100, 15),
  Weight = rnorm(n, 50, 10),
  Yield = rnorm(n, 75, 12),
  Quality = rnorm(n, 85, 8),
  Genotype = sample(c("AA", "AB", "BB"), n, replace = TRUE)
)
df$Weight <- df$Weight + 0.4 * df$Height + rnorm(n, 0, 5)
df$Yield <- df$Yield + 0.5 * df$Height + rnorm(n, 0, 8)

# Basic pairs plot
pairs(df[, 1:4], main = "Scatter Plot Matrix")

# ggplot2 pairs plot with GGally
ggpairs(df, columns = 1:4, aes(color = Genotype),
        title = "Scatter Plot Matrix by Genotype",
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.5)),
        diag = list(continuous = wrap("barDiag", alpha = 0.7)))

# Customized pairs plot
ggpairs(df, columns = 1:4,
        upper = list(continuous = "density"),
        lower = list(continuous = wrap("points", size = 0.5)),
        diag = list(continuous = "bar")) +
  theme_minimal() +
  ggtitle("Customized Scatter Plot Matrix")
```

### 5. Manhattan Plots with qqman

```r
library(qqman)

# Create sample GWAS data
set.seed(42)
n_snps <- 10000
gwas_results <- data.frame(
  SNP = paste0("rs", 1:n_snps),
  CHR = sample(1:5, n_snps, replace = TRUE),
  BP = sample(1:100000000, n_snps),
  P = runif(n_snps, 0, 1)
)

# Add significant hits
gwas_results$P[sample(1:1000, 20)] <- runif(20, 1e-10, 1e-6)
gwas_results$P[sample(2001:3000, 30)] <- runif(30, 1e-8, 1e-5)

# Basic Manhattan plot
manhattan(gwas_results, chr = "CHR", bp = "BP", snp = "SNP", p = "P",
          main = "Manhattan Plot", col = c("blue4", "orange3"),
          suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8))

# Highlight specific SNPs
snpsOfInterest <- gwas_results$SNP[gwas_results$P < 1e-6]
manhattan(gwas_results, highlight = snpsOfInterest,
          main = "Manhattan Plot with Highlighted SNPs")

# QQ plot
qq(gwas_results$P, main = "QQ Plot of GWAS P-values")

# Enhanced Manhattan with custom colors
manhattan(gwas_results, chr = "CHR", bp = "BP", snp = "SNP", p = "P",
          col = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"),
          chrlabs = c("1", "2", "3", "4", "5"),
          main = "Customized Manhattan Plot")
```

### 6. QTL Plots with qtl2

```r
library(qtl2)

# Create sample QTL data
# In practice, you would load real data with read_cross()
# For demonstration, we'll create a simple scan object

# Simulate a genome scan
set.seed(42)
chromosomes <- c(1:5)
map <- vector("list", 5)
for(i in 1:5) {
  map[[i]] <- seq(0, 100, length.out = 50)
}
names(map) <- as.character(1:5)

# Simulate LOD scores
lod_scores <- vector("list", 5)
for(i in 1:5) {
  lod <- rgamma(50, shape = 2, scale = 1)
  if(i == 3) {
    # Add QTL peak on chromosome 3
    lod <- lod + 15 * dnorm(seq(0, 100, length.out = 50), 45, 10) * 10
  }
  lod_scores[[i]] <- matrix(lod, ncol = 1)
}

# Create scan1 object
scan_result <- list(lod = lod_scores, map = map)
class(scan_result) <- c("scan1", "list")

# Plot LOD scores
plot(scan_result, map, lodcolumn = 1, col = "blue",
     main = "QTL Scan LOD Scores")
abline(h = 3, col = "red", lty = 2, lwd = 2)
legend("topright", legend = "Significance threshold (LOD=3)", 
       col = "red", lty = 2, lwd = 2)

# Plot single chromosome
plot(scan_result, map, chr = 3, col = "darkgreen",
     main = "Chromosome 3 QTL Profile")
abline(h = 3, col = "red", lty = 2)

# Find peaks
peaks <- find_peaks(scan_result, map, threshold = 3)
print(peaks)
```

### 7. Effect Size Plots with ggplot2

```r
library(ggplot2)
library(dplyr)

# Create effect size data
set.seed(42)
n_qtl <- 8
traits <- c("Height", "Yield", "Quality", "Disease_Res")

effects_df <- expand.grid(
  QTL = paste0("QTL_", 1:n_qtl),
  Trait = traits
) %>%
  mutate(
    Effect = rnorm(n(), 0, 2),
    SE = runif(n(), 0.5, 1.5),
    Lower = Effect - 1.96 * SE,
    Upper = Effect + 1.96 * SE
  )

# Forest plot
ggplot(effects_df, aes(x = Effect, y = QTL, color = Trait)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), 
                 height = 0.2, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(title = "QTL Effect Sizes (Forest Plot)",
       x = "Effect Size", y = "QTL") +
  facet_wrap(~Trait, ncol = 1, scales = "free_y")

# Allele effect plot
allele_effects <- data.frame(
  Genotype = rep(c("AA", "AB", "BB"), 3),
  Trait = rep(c("Height", "Yield", "Quality"), each = 3),
  Mean = c(100, 115, 95, 50, 65, 45, 80, 85, 78),
  SE = runif(9, 2, 5)
)

ggplot(allele_effects, aes(x = Genotype, y = Mean, fill = Trait)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(width = 0.9), width = 0.25) +
  facet_wrap(~Trait, scales = "free_y") +
  theme_minimal() +
  labs(title = "Allele Effects by Genotype",
       y = "Phenotype Value") +
  scale_fill_brewer(palette = "Set2")
```

### 8. Multi-Environment Visualization with ggplot2

```r
library(ggplot2)
library(dplyr)
library(tidyr)

# Create GxE data
set.seed(42)
genotypes <- paste0("G", 1:5)
environments <- paste0("Env", 1:4)

gxe_data <- expand.grid(
  Genotype = genotypes,
  Environment = environments
) %>%
  mutate(
    Yield = rnorm(20, 50, 10) + 
           rep(c(0, 5, -3, 8), 5) +  # Environment effects
           rep(c(0, 5, -5, 10, 2), each = 4)  # Genotype effects
  )

# Interaction plot
ggplot(gxe_data, aes(x = Environment, y = Yield, color = Genotype, group = Genotype)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Genotype × Environment Interaction",
       y = "Yield") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Heatmap
ggplot(gxe_data, aes(x = Environment, y = Genotype, fill = Yield)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.1f", Yield)), color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "GxE Heatmap")

# Stability plot
geno_stats <- gxe_data %>%
  group_by(Genotype) %>%
  summarise(
    Mean_Yield = mean(Yield),
    SD_Yield = sd(Yield)
  )

ggplot(geno_stats, aes(x = Mean_Yield, y = SD_Yield, label = Genotype)) +
  geom_point(size = 5, color = "blue") +
  geom_text(vjust = -1) +
  theme_minimal() +
  labs(title = "Stability Analysis",
       subtitle = "Lower right = high yield, stable",
       x = "Mean Yield", y = "Standard Deviation") +
  geom_hline(yintercept = mean(geno_stats$SD_Yield), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = mean(geno_stats$Mean_Yield), linetype = "dashed", alpha = 0.5)
```

---

## Summary Table: Visualization Types and Packages

| Visualization | Python Package | R Package | Use Case |
|---------------|--------------|-----------|----------|
| Box plots | matplotlib, seaborn | ggplot2, lattice | Phenotype distributions |
| Density plots | matplotlib, seaborn | ggplot2, lattice | Distribution shape |
| Violin plots | seaborn | ggplot2, ggridges | Distribution + density |
| Heatmaps | seaborn, matplotlib | corrplot, ComplexHeatmap | Correlations, genotypes |
| Scatter matrices | seaborn, pandas | GGally, pairs | Multi-trait correlations |
| Manhattan plots | qmplot, manual | qqman, ggmanh | GWAS results |
| QQ plots | scipy, matplotlib | qqman, ggplot2 | Normality assessment |
| LOD plots | matplotlib | qtl2 | QTL scans |
| Effect plots | matplotlib, seaborn | ggplot2 | QTL effects |
| GxE plots | seaborn, matplotlib | ggplot2 | Multi-environment |

---

## Best Practices

1. **Color Selection**: Use colorblind-friendly palettes (viridis, colorbrewer)
2. **Resolution**: Save figures at 300 DPI for publications
3. **File Formats**: Use PNG for web, PDF/SVG for publications
4. **Annotations**: Always label axes, add titles, and include legends
5. **Statistical Tests**: Include significance indicators where appropriate
6. **Data Integrity**: Verify data before plotting (check for NAs, outliers)
7. **Reproducibility**: Set random seeds for simulated data

---

## Installation Verification

### Python
```python
# Verify all packages are installed
import importlib
packages = ['matplotlib', 'seaborn', 'plotly', 'pandas', 'numpy', 'scipy']
for pkg in packages:
    try:
        importlib.import_module(pkg)
        print(f"✓ {pkg} installed")
    except ImportError:
        print(f"✗ {pkg} NOT installed")
```

### R
```r
# Verify all packages are installed
packages <- c("ggplot2", "lattice", "corrplot", "RColorBrewer", 
              "qqman", "GGally", "ggridges", "qtl2")
for(pkg in packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "installed\n")
  } else {
    cat("✗", pkg, "NOT installed\n")
  }
}
```
