# Qmapper: Chromosome Ideogram Visualization Example

## What is This?

**Qmapper** (like QTLmax's tool) creates a visual map of chromosomes showing where SNPs are located physically. Think of it like a subway map, but for DNA—showing where the "stations" (SNPs) are along the "lines" (chromosomes).

### Key Terms Explained

- **Ideogram**: A simplified diagram of a chromosome showing its structure. Like a map that shows the main features but not every detail.

- **Chromosome**: A long, tightly coiled DNA molecule. Humans have 23 pairs; plants like rice might have 12.

- **SNP (Single Nucleotide Polymorphism)**: A single position in DNA where people differ. Like finding "color" vs "colour"—same meaning, different spelling.

- **Centromere**: The "waist" of a chromosome where it's pinched. This is important because it affects how chromosomes separate during cell division.

- **Physical Position**: The actual location on the chromosome in base pairs (or megabases, Mb). Unlike genetic position (cM), this is the real physical distance.

- **Manhattan Plot**: A type of graph used in GWAS that shows -log10(p-values) across the genome. Looks like the Manhattan skyline.

## Why Would Anyone Do This?

1. **See the Big Picture**: When you have 10,000 SNPs, it's hard to understand where they all are. An ideogram makes it visual.

2. **Find QTL Locations**: When you find a significant SNP in GWAS, you want to know exactly where it is on the chromosome.

3. **Compare Populations**: Different populations might have different SNP patterns in certain regions.

4. **Publication Figures**: Chromosome ideograms look professional in papers and presentations.

5. **Breeding Programs**: Plant and animal breeders need to know where genes are physically located.

## How It Works (The Process)

1. **Generate Chromosome Data**: Create 5 synthetic chromosomes with realistic lengths (like rice chromosomes: 30-43 Mb each).

2. **Place SNPs**: Generate 100+ SNPs and place them randomly along the chromosomes, with some clustering in QTL regions.

3. **Draw Chromosomes**: Create stylized chromosome shapes with:
   - Rectangular body
   - Centromere "waist" in the middle
   - SNP markers as colored ticks

4. **Color by Significance**: 
   - Red = Genome-wide significant (p < 5×10⁻⁸)
   - Orange = Suggestive (p < 10⁻⁵)
   - Blue = In QTL region
   - Gray = Not significant

5. **Create Combined View**: Show ideogram above Manhattan plot for easy comparison.

## Input → Process → Output

### Input
| File | Description | What It Contains |
|------|-------------|------------------|
| `snp_positions.csv` | SNP data | 100 SNPs with positions, p-values, and significance |

**Example Data Preview:**
```csv
SNP,Chromosome,Position_Mb,P_value,In_QTL
rs0,Chr1,12.34,0.0012,True
rs1,Chr1,15.67,0.0456,False
rs2,Chr2,8.90,0.00001,True
...
```

**What Each Column Means:**
- **SNP**: Unique identifier for the variant
- **Chromosome**: Which chromosome it's on (Chr1-Chr5)
- **Position_Mb**: Position in megabases (millions of base pairs)
- **P_value**: Statistical significance from GWAS (smaller = more significant)
- **In_QTL**: Whether it's in a known QTL region

### Process
1. **Data Generation**: We create 5 chromosomes with lengths based on rice (30-43 Mb), then place SNPs with some clustering.

2. **Ideogram Creation**: 
   - Draw chromosome rectangles
   - Add centromere constrictions
   - Place SNP markers at their physical positions
   - Color by significance

3. **QTL Highlighting**: Semi-transparent red boxes show QTL regions.

4. **Manhattan Alignment**: Create a Manhattan plot that aligns with the ideogram above it.

### Output
| File | Description | What It Shows |
|------|-------------|---------------|
| `ideogram.png` | Chromosome ideograms | 5 chromosomes with SNP positions |
| `qmapper_combined.png` | Combined view | Ideogram + Manhattan alignment |
| `snp_positions.csv` | SNP data | Positions and significance for all SNPs |

**Visualization Explanation:**

**Ideogram.png - Individual Chromosomes**
- Each row is one chromosome
- Gray rectangle = chromosome body
- Dark gray band = centromere
- Vertical lines = SNPs (colored by significance)
- Dashed red box = QTL region
- Bottom label shows chromosome name and length

**Qmapper_combined.png - Aligned View**
- **Top**: Ideogram showing all 5 chromosomes side by side
- **Bottom**: Manhattan plot showing GWAS results
- The x-axis aligns so you can see which SNPs on the ideogram correspond to peaks in the Manhattan plot
- Red/blue horizontal lines = significance thresholds

## Running the Example

```bash
cd examples/qmapper-ideogram
python run_qmapper.py
```

## Expected Runtime
- Data generation: < 1 second
- Ideogram drawing: ~2 seconds
- Combined plot: ~3 seconds
- Total: ~5 seconds

## Acceptance Criteria
- [x] **5 chromosomes drawn**: Each has correct length and centromere.
- [x] **SNPs positioned correctly**: Plotted at their physical positions.
- [x] **QTL regions highlighted**: Semi-transparent boxes visible.
- [x] **Combined view aligned**: Manhattan plot matches ideogram.

## Tools Used
- **matplotlib.patches**: For drawing chromosome rectangles
- **numpy/pandas**: For data generation and handling
- **matplotlib**: For visualization

## Understanding the Output

**Color Coding:**
- **Red SNPs**: p < 5×10⁻⁸ (genome-wide significant)
  - These are the "real" associations
  - Likely to be true genetic effects
  
- **Orange SNPs**: p < 10⁻⁵ (suggestive)
  - Might be real, might be noise
  - Worth investigating further

- **Blue SNPs**: In QTL region
  - Located in a known QTL
  - May not be individually significant

- **Gray SNPs**: Not significant
  - Most SNPs fall here
  - Expected under null hypothesis

**Why the Centromere Matters:**
- Centromeres are constricted regions where the chromosome attaches to spindle fibers
- Genes near centromeres often have reduced recombination
- This affects genetic mapping!

**Physical vs. Genetic Distance:**
- Physical distance: Actual base pairs (Mb)
- Genetic distance: Recombination frequency (cM)
- 1 cM ≈ 1 Mb in most regions, but can vary 10-fold!

## Real-World Application

In a real study, you might use this to:
- Map 10,000 SNPs across a genome
- Identify QTL regions visually
- Compare with known gene locations
- Create publication-quality figures
- Present results to breeders/farmers

## Comparison to QTLmax

QTLmax's Qmapper tool does this interactively:
- Click SNPs to see details
- Zoom in on regions
- Filter by significance
- Export high-resolution images

Our version generates static images but uses the same principles!

## Troubleshooting

**SNPs not showing**: Check that positions are within chromosome lengths.

**Ideogram too small**: Adjust figsize in plt.subplots().

**Colors not right**: Verify p-values are -log10 transformed correctly.

**Misalignment**: Ensure chromosome order matches in both plots.
