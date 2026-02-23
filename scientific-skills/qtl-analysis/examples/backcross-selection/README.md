# Backcross Selection Example

## What is This?

**Backcross selection** is a plant/animal breeding technique where you cross an offspring with one of its parents (or a parent variety) to introduce specific traits from that parent while keeping most of the other parent's desirable characteristics. Think of it like selectively breeding to "import" a specific trait.

### Key Terms Explained

- **Backcross**: Mating an offspring back to one of its parents (or the parent line). If you cross A × B, then cross the result × B, that's a backcross.

- **Donor Parent**: The parent that has the trait you want to introduce (e.g., disease resistance).

- **Recurrent Parent**: The parent that has most of the desirable traits you want to keep (e.g., high yield). You cross back to this parent repeatedly.

- **Genetic Similarity**: How much DNA two individuals share. Calculated from SNP data. 100% = identical; 0% = completely different.

- **Similarity Matrix**: A table showing how similar each individual is to every other individual.

- **BC Generation**: Backcross generation (BC1, BC2, etc.). Each generation is one more backcross to the recurrent parent.

## Why Would Anyone Do This?

1. **Introduce Specific Traits**: Want disease resistance from wild rice but keep the yield of cultivated rice? Backcross it!

2. **Clean Up Background**: Each backcross replaces ~50% of the donor genome with recurrent parent DNA. After 7-8 generations, you have ~99% recurrent parent.

3. **Marker-Assisted Selection**: Use DNA markers to track the target gene so you know which offspring to select.

4. **Shorter Timeline**: Traditional breeding takes decades. Backcrossing + markers can do it in years.

## How It Works (The Process)

Imagine breeding for disease resistance:

**Generation 0 (Parents):**
- Parent A: High yielding, susceptible to disease (recurrent parent - 90% of genome desired)
- Parent B: Disease resistant, low yield (donor parent - 10% of genome desired, just the resistance gene)

**Generation 1 (F1):**
- Cross A × B → 50% A, 50% B

**Generation 2 (BC1):**
- Cross best F1 × A → 75% A, 25% B
- Select offspring with resistance gene + most A-like background

**Generation 3-7 (BC2-BC6):**
- Keep crossing best offspring × A
- Each generation: %A increases, %B decreases
- After BC7: ~99% A genome, but keeping resistance gene

## Input → Process → Output

### Input
| File | Description | What It Contains |
|------|-------------|------------------|
| `parents.csv` | Parent data | Donor and recurrent parent genotypes |
| `offspring_genotypes.csv` | Offspring data | Genotypes of backcross progeny |

**Example Data Preview:**
```csv
Sample,SNP1,SNP2,SNP3,...,Generation,Parent_Type
Parent_A,0,2,1,...,0,Recurrent
Parent_B,2,0,0,...,0,Donor
BC1_001,1,1,0,...,1,BC1
BC2_015,0,2,1,...,2,BC2
...
```

**What Each Column Means:**
- **Sample**: Individual identifier
- **SNP1, SNP2, etc.**: Genotypes (0, 1, or 2)
- **Generation**: Which backcross generation (0=parents, 1=BC1, 2=BC2, etc.)
- **Parent_Type**: Which type of parent or offspring

### Process
1. **Generate Pedigree**: 
   - Start with two parents (100 individuals each)
   - Create F1 generation (50 individuals, 50% from each parent)
   - Create BC1 by crossing F1 × recurrent parent (50 individuals)
   - Continue to BC2, BC3, etc.

2. **Calculate Similarity**: 
   - For each offspring, calculate % similarity to recurrent parent
   - Formula: % matching alleles / total SNPs
   - Expected: 50% in F1, 75% in BC1, 87.5% in BC2, etc.

3. **Select Best**: 
   - Choose offspring with highest similarity to recurrent parent
   - Must retain the target trait (disease resistance)

4. **Track Progress**: 
   - Plot similarity over generations
   - Show donor genome proportion decreasing

### Output
| File | Description | What It Shows |
|------|-------------|---------------|
| `backcross_selection.png` | Four-panel visualization | Pedigree, similarity matrix, progress plot, selection results |
| `similarity_matrix.csv` | Genetic similarities | % similarity of each offspring to recurrent parent |
| `selection_candidates.csv` | Best candidates | Top offspring to select for next generation |

**Visualization Explanation:**

**Panel 1: Pedigree Diagram**
- Shows the crossing scheme
- P1 (Recurrent) × P2 (Donor) → F1
- F1 × P1 → BC1
- BC1 × P1 → BC2, etc.
- Arrows show gene flow

**Panel 2: Similarity Heatmap**
- Rows = offspring
- Columns = SNPs
- Color shows which allele (recurrent=green, donor=red, heterozygous=yellow)
- Over generations, you see more green

**Panel 3: Similarity Progress**
- X-axis = Generation (F1, BC1, BC2, ...)
- Y-axis = % Similarity to recurrent parent
- Each point = one offspring
- Trend line shows expected increase
- Target line shows when to stop (usually 93-99%)

**Panel 4: Selection Plot**
- Histogram of similarity scores
- Red line = selection threshold
- Selected offspring highlighted
- Shows we're choosing the most recurrent-like individuals

## Running the Example

```bash
cd examples/backcross-selection
python run_backcross.py
```

## Expected Runtime
- Data generation: < 1 second
- Similarity calculation: ~2 seconds
- Selection: < 1 second
- Visualization: ~3 seconds
- Total: ~6 seconds

## Acceptance Criteria
- [x] **7 generations created**: From F1 to BC6 (or beyond)
- [x] **Similarity increases**: BC1 ~75%, BC2 ~87.5%, etc.
- [x] **Best candidates selected**: Top 10% chosen for next cross
- [x] **Visualization shows progress**: Clear trend toward recurrent parent

## Tools Used
- **numpy**: For genetic calculations
- **pandas**: For data handling
- **matplotlib**: For visualization
- **scipy**: For statistics

## The Math Behind It

**Expected Recurrent Parent Genome:**
- F1: 50%
- BC1: 75% (50% + 50%)
- BC2: 87.5% (75% + 12.5%)
- BC3: 93.75%
- BC4: 96.875%
- BC5: 98.4375%
- BC6: 99.21875%

Formula: % Recurrent = 1 - (1/2)^(n+1)
where n = backcross generation number

**Selection Strategy:**
- Target: Recover 93-99% recurrent parent genome
- Usually takes 4-7 backcross generations
- Each generation adds ~12.5% of the recurrent parent

## Real-World Application

Common uses in agriculture:
- **Rice**: Introduce submergence tolerance from wild rice
- **Wheat**: Add rust resistance from wild relatives
- **Corn**: Incorporate pest resistance from landraces
- **Tomato**: Transfer disease resistance from wild species

**What You Need:**
1. Two parent varieties
2. DNA markers (SNPs) across the genome
3. A trait you can phenotype (resistance, color, etc.)
4. Time for several generations (years)

## Comparison to QTLmax

QTLmax has a dedicated "Backcross" tab that:
- Calculates similarity automatically from SNP data
- Provides selection recommendations
- Tracks progress across generations
- Shows which offspring to select

Our version does the same calculations with open-source tools!

## Troubleshooting

**Similarity not increasing**: Check that you're crossing with the recurrent parent, not the donor!

**Too much variation**: Use more SNPs (at least 100) for accurate similarity estimates.

**Offspring all look like one parent**: Verify your SNP data is coded correctly (0/1/2).

**Can't distinguish generations**: Make sure SNPs are polymorphic between parents (different allele frequencies).
