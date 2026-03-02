---
phase: 02-advanced-gwas
plan: 03
subsystem: qtl-analysis
tags: [gwas, genomic-control, lambda, rare-variants, burden, skat]
requirements-completed:
  - QTLGW-05
  - QTLGW-06
duration: 6min
completed: 2026-02-24
---

# Phase 2 Plan 3: Genomic Control + Rare Variant Tests Summary

Implemented genomic inflation correction and rare variant set testing examples with runnable scripts and committed outputs.

## Accomplishments
- Added `genomic-control` example that computes lambda and applies genomic control adjustment.
- Added `rare-variant-tests` example with burden and SKAT-like per-gene tests.
- Generated all expected outputs and wrote Input -> Process -> Output READMEs.

## Key Results
- Genomic control: lambda reduced from 1.350 to 1.000 after adjustment.
- Rare variant tests: per-gene significance table and method comparison plot produced.

## Files Created
- `scientific-skills/qtl-analysis/examples/genomic-control/run_genomic_control.py`
- `scientific-skills/qtl-analysis/examples/genomic-control/README.md`
- `scientific-skills/qtl-analysis/examples/genomic-control/output/qq_before_after_gc.png`
- `scientific-skills/qtl-analysis/examples/genomic-control/output/genomic_control_results.csv`
- `scientific-skills/qtl-analysis/examples/genomic-control/output/lambda_summary.csv`
- `scientific-skills/qtl-analysis/examples/rare-variant-tests/run_rare_variants.py`
- `scientific-skills/qtl-analysis/examples/rare-variant-tests/README.md`
- `scientific-skills/qtl-analysis/examples/rare-variant-tests/output/rare_variant_method_comparison.png`
- `scientific-skills/qtl-analysis/examples/rare-variant-tests/output/rare_variant_gene_results.csv`
- `scientific-skills/qtl-analysis/examples/rare-variant-tests/output/rare_variant_summary.csv`

## Requirements Satisfied
- QTLGW-05
- QTLGW-06
