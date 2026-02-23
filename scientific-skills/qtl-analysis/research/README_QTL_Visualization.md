# QTL and Genomic Data Visualization Guide

Complete guide for visualizing QTL (Quantitative Trait Loci) and genomic data using Python and R.

## 📁 Repository Contents

### Documentation
- **`qtl_visualization_guide.md`** - Comprehensive guide with detailed code examples
- **`qtl_visualization_quick_ref.md`** - Quick reference card for common visualizations
- **`README_QTL_Visualization.md`** - This file

### Code
- **`generate_sample_qtl_data.py`** - Python script to generate sample datasets

### Sample Data (`sample_qtl_data/`)
- `phenotype_data.csv` - 200 samples with multiple traits
- `gwas_results.csv` - 9,006 SNPs across 5 chromosomes
- `qtl_scan.csv` - 381 markers with LOD scores
- `genotype_matrix.csv` - 50 samples × 20 markers
- `gxe_data.csv` - Genotype × Environment data
- `effect_sizes.csv` - QTL effect sizes

---

## 🚀 Quick Start

### Python

```bash
# Install required packages
pip install matplotlib seaborn plotly pandas numpy scipy

# Optional: specialized packages
pip install qmplot gwaslab

# Run sample data generator
python generate_sample_qtl_data.py
```

### R

```r
# Install CRAN packages
install.packages(c("ggplot2", "qqman", "corrplot", "GGally", "ggridges"))

# Install Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("qtl2", "ggmanh"))
```

---

## 📊 Visualization Types Covered

### Basic Phenotype Visualizations
1. **Box Plots** - Distribution by genotype/group
2. **Density Plots** - Distribution shape
3. **Violin Plots** - Distribution with density

### Advanced Visualizations
4. **Heatmaps** - Correlations, genotype data
5. **Scatter Plot Matrices** - Multi-trait relationships
6. **Clustermaps** - Hierarchical clustering

### QTL-Specific Visualizations
7. **Manhattan Plots** - GWAS results
8. **QQ Plots** - Normality assessment
9. **LOD Score Plots** - QTL scans
10. **Effect Size Plots** - QTL effects
11. **Allele Effect Plots** - Genotype effects
12. **Multi-Environment Plots** - G×E interactions

---

## 📖 Usage Examples

### Python Example

```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
df = pd.read_csv('sample_qtl_data/phenotype_data.csv')

# Box plot
sns.boxplot(data=df, x='Genotype', y='Plant_Height', hue='Environment')
plt.title('Plant Height by Genotype and Environment')
plt.savefig('boxplot.png', dpi=300, bbox_inches='tight')
plt.show()

# Correlation heatmap
corr = df[['Plant_Height', 'Yield', 'Quality', 'Disease_Resistance']].corr()
sns.heatmap(corr, annot=True, cmap='coolwarm', center=0)
plt.title('Phenotype Correlations')
plt.savefig('correlation_heatmap.png', dpi=300)
plt.show()
```

### R Example

```r
library(ggplot2)
library(qqman)

# Load data
df <- read.csv('sample_qtl_data/phenotype_data.csv')
gwas <- read.csv('sample_qtl_data/gwas_results.csv')

# Box plot
ggplot(df, aes(x=Genotype, y=Plant_Height, fill=Environment)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle('Plant Height by Genotype and Environment')

# Manhattan plot
manhattan(gwas, chr="CHR", bp="BP", snp="SNP", p="P",
          main="GWAS Manhattan Plot")
```

---

## 🎨 Color Palette Recommendations

### Colorblind-Friendly Options
- **Python**: `sns.color_palette("colorblind")`, `"viridis"`
- **R**: `scale_color_brewer(palette="Dark2")`, `viridis::scale_color_viridis()`

### Categorical Data (≤8 categories)
- `"Set1"`, `"Set2"`, `"Dark2"`, `"Paired"`

### Continuous Data
- `"viridis"`, `"plasma"`, `"coolwarm"`, `"RdBu_r"`

---

## 📏 Publication Guidelines

### Figure Specifications
- **Width**: 3.5-7 inches (single column: 3.5", double column: 7")
- **Resolution**: 300 DPI minimum
- **Font size**: ≥8pt for text, ≥10pt for axis labels
- **Format**: PDF or SVG for publications, PNG for web

### Best Practices
- Use colorblind-friendly palettes
- Include clear axis labels with units
- Add informative titles or captions
- Show statistical significance where applicable
- Maintain consistent styling across figures

---

## 🔧 Troubleshooting

### Python
- **Import errors**: `pip install missing_package`
- **Empty plots**: Add `plt.show()` at the end
- **Overlapping labels**: Use `plt.tight_layout()`

### R
- **Package not found**: Use `install.packages()` or `BiocManager::install()`
- **Plot not displaying**: Check graphics device
- **Legend cut off**: Adjust figure dimensions or legend position

---

## 📚 Additional Resources

### Python Documentation
- [Matplotlib](https://matplotlib.org/)
- [Seaborn](https://seaborn.pydata.org/)
- [Plotly](https://plotly.com/python/)
- [Pandas](https://pandas.pydata.org/)

### R Documentation
- [ggplot2](https://ggplot2.tidyverse.org/)
- [qqman](https://cran.r-project.org/web/packages/qqman/)
- [qtl2](https://kbroman.org/qtl2/)
- [corrplot](https://cran.r-project.org/web/packages/corrplot/)

### QTL Analysis Resources
- [R/qtl2 User Guide](https://kbroman.org/qtl2/assets/vignettes/user_guide.html)
- [GWAS Tutorial](https://github.com/MareesAT/GWA_tutorial/)

---

## 📝 Citation

If you use this guide in your research, please cite the relevant packages:

**Python:**
- Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. Computing in Science & Engineering.
- Waskom, M. L. (2021). seaborn: statistical data visualization. Journal of Open Source Software.

**R:**
- Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer.
- Turner, S. D. (2018). qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots. Journal of Open Source Software.
- Broman, K. W. et al. (2019). R/qtl2: Software for mapping quantitative trait loci with high-dimensional data and multiparent populations. Genetics.

---

## 📧 Support

For questions or issues:
1. Check the troubleshooting section
2. Review the package documentation
3. Search Stack Overflow or Biostars
4. Open an issue on the relevant package's GitHub repository

---

**Last Updated**: February 2025
