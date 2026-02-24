# Rare Variant Tests

## What This Example Does
Compares two gene-level rare variant association strategies:
- Burden test (collapses variants into one weighted score)
- SKAT-like test (variance-component style set test)

## Input -> Process -> Output

### Input
- Simulated rare variants (MAF 0.1% to 4%)
- Gene-set grouping for 25 genes
- Continuous phenotype with a subset of causal genes

### Process
1. Build weighted variant sets by gene
2. Run burden test and SKAT-like test per gene
3. Apply Bonferroni threshold at gene level
4. Compare signal patterns across methods

### Output
- `output/rare_variant_method_comparison.png`
- `output/rare_variant_gene_results.csv`
- `output/rare_variant_summary.csv`

## Run
```bash
python run_rare_variants.py
```

## Interpretation
- Burden works well when effects point in the same direction.
- SKAT-like tests handle mixed effect directions better.
- Compare both to avoid missing biologically plausible genes.
