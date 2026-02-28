<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Haplotype Analysis with Dendrogram Example

## What is This?

**Haplotype analysis** is a way of looking at groups of DNA letters (SNPs) that tend to be inherited together. Think of it like this: when you inherit DNA from your parents, you don't get individual random letters—you get chunks that stay together.

### Key Terms Explained

- **Haplotype**: A group of SNPs (DNA variants) that are physically close together on a chromosome and tend to be inherited as a group. Like a "sentence" of DNA rather than individual "words."

- **Linkage Disequilibrium (LD)**: A measure of how often two SNPs are inherited together. High LD means they're usually inherited together (linked); low LD means they assort independently.

- **Dendrogram**: A tree-like diagram showing hierarchical relationships. In genetics, it shows which SNPs are most similar to each other.

- **r² (r-squared)**: A statistical measure ranging from 0 to 1 that quantifies the correlation between two SNPs. r² = 1 means perfectly correlated; r² = 0 means independent.

- **Hierarchical Clustering**: A method of grouping similar items together. SNPs with high LD are clustered together.

## Why Would Anyone Do This?

1. **Simplify Complex Data**: Instead of analyzing 50,000 individual SNPs, you can analyze 500 haplotype blocks—much easier!

2. **Find Disease Genes**: Haplotypes can be more informative than individual SNPs for finding disease-causing variants.

3. **Population History**: Haplotype blocks reveal historical recombination events and population migrations.

4. **Better Imputation**: Knowing haplotype structure helps predict missing genotypes more accurately.

## How It Works (The Process)

1. **Generate DNA Data**: We create synthetic DNA data for 100 people across 50 SNPs arranged in blocks.

2. **Calculate LD**: For every pair of SNPs, we calculate r² (how correlated they are).

3. **Cluster SNPs**: Using hierarchical clustering, we group SNPs that have high LD (similar patterns).

4. **Visualize**: We create four plots showing the LD patterns, the dendrogram, and the haplotype blocks.

## Input → Process → Output

### Input
| File | Description | What It Contains |
|------|-------------|------------------|
| `haplotype_genotypes.csv` | DNA data | 100 people × 50 SNPs, with values 0, 1, or 2 (representing 0, 1, or 2 copies of the variant allele) |

**Example Data Preview:**
```
SNP,rs0_0,rs0_1,rs0_2,...
Person_1,2,1,0,...
Person_2,0,2,1,...
...
```
- **0**: No copies of the variant (homozygous reference)
- **1**: One copy (heterozygous)
- **2**: Two copies (homozygous variant)

### Process
1. **Data Generation**: We create 5 haplotype blocks, each with different patterns of SNPs that tend to be inherited together.

2. **LD Calculation**: For each pair of SNPs, we calculate r² using the formula: r² = (correlation coefficient)²

3. **Hierarchical Clustering**: We use Ward's method to cluster SNPs based on their LD distances (distance = 1 - r²).

4. **Visualization**: Four-panel figure showing different aspects of the analysis.

### Output
| File | Description | What It Shows |
|------|-------------|---------------|
| `haplotype_analysis.png` | Four-panel visualization | LD heatmap, dendrogram, haplotype blocks, and LD decay curve |
| `haplotype_blocks.csv` | Block assignments | Which haplotype block each SNP belongs to |
| `ld_matrix.csv` | Full LD matrix | r² values for all SNP pairs |

**Visualization Explanation:**

**Panel 1: LD Heatmap**
- Shows r² between every pair of SNPs (darker = higher LD)
- Yellow dashed lines show block boundaries
- You can see the 5 diagonal blocks—each is a haplotype block

**Panel 2: Dendrogram**
- Tree showing which SNPs are most similar
- The y-axis shows "distance" (1 - r²)
- SNPs that join at the bottom are most similar
- The colored branches show different haplotype blocks

**Panel 3: Haplotype Blocks**
- Each color represents a different haplotype block
- Shows which SNPs cluster together
- 5 distinct colors = 5 haplotype blocks

**Panel 4: LD Decay**
- Shows how LD decreases with physical distance
- Nearby SNPs have high LD; distant SNPs have low LD
- The red line shows the trend

## Running the Example

```bash
cd examples/haplotype-analysis
python run_haplotype.py
```

## Expected Runtime
- Data generation: < 1 second
- LD calculation: ~5 seconds (calculating 2,500 pairs)
- Clustering: ~1 second
- Visualization: ~2 seconds

## Acceptance Criteria
- [x] **5 haplotype blocks identified**: The clustering correctly finds the 5 blocks we created.
- [x] **Dendrogram generated**: Shows the hierarchical relationships between SNPs.
- [x] **LD heatmap shows block boundaries**: The yellow lines match the expected block structure.
- [x] **Block assignments saved**: Each SNP is assigned to its correct block.

## Tools Used
- **scipy.cluster.hierarchy**: For hierarchical clustering and creating dendrograms
- **scipy.spatial.distance**: For converting LD to distances
- **numpy**: For matrix operations and LD calculations
- **matplotlib**: For creating the four-panel visualization

## Real-World Application
In a real study, you might use this to:
- Analyze 500,000 SNPs and reduce them to ~5,000 haplotype blocks
- Find which blocks contain disease-associated variants
- Study recombination hotspots (where haplotype blocks break)
- Compare haplotype structure between populations

## Troubleshooting
**LD values are all 0 or 1**: This can happen with small sample sizes. With 100 samples, we get realistic intermediate values.

**Dendrogram is flat**: This means SNPs aren't forming distinct clusters. Check that you have enough samples and SNPs.

**Blocks don't match expectations**: Try adjusting the clustering threshold or using a different distance metric.
