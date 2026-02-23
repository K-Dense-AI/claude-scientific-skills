# G×E (Genotype × Environment) Genomic Prediction

## What This Example Does

This example demonstrates **multi-environment genomic prediction** — predicting breeding values when the same genotypes are grown across different environments (locations, years, conditions).

## Why It Matters

Genotype × Environment interaction (G×E) is ubiquitous in plant and animal breeding:
- A top-performing line in one environment may perform poorly in another
- Multi-environment models can leverage data across environments to improve predictions
- Reaction norms show how each genotype responds to environmental variation

## Methods Compared

| Method | Description | Best When |
|--------|-------------|-----------|
| **Single-Env** | Train separate model per environment | G×E is strong, environments very different |
| **Multi-Env** | Joint model across all environments | Some correlation across environments, shared genetic effects |

## Reaction Norms

A reaction norm shows how a genotype's phenotype changes across environments:
- Flat line = stable across environments (robust)
- Steep slope = environment-sensitive (variable)

## Running the Example

```bash
cd scientific-skills/qtl-analysis/examples/gxe-prediction
python run_gxe.py
```

## Input → Process → Output

### Input
- 200 individuals, 300 markers, 3 environments
- 25 QTL with both main effects and G×E effects
- Heritability: 0.5

### Process
- Single-environment model: Ridge regression per environment
- Multi-environment model: GBLUP kernel across environments
- 5-fold CV per environment

### Output
- `output/reaction_norms.png` — G×E visualization
- `output/gxe_accuracy.png` — Single vs Multi-Env comparison
- `output/gxe_results.csv` — Per-fold results
- `output/gxe_summary.csv` — Mean accuracy per environment

## Key Insight

Multi-environment models typically outperform single-environment when:
- G×E is moderate (not too strong)
- Training data is limited per environment
- Genetic correlation across environments is >0.3
