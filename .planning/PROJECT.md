# Meta-QTLmax: Complete Open-Source Breeding Genomics Skills

## What This Is

Two complementary scientific skills that together replace commercial breeding genomics platforms (QTLmax, Golden Helix SVS, Agronomix/Genovix) with open-source alternatives. The first skill (`qtl-analysis`) handles all genomic analysis — GWAS, QTL mapping, genomic prediction, population structure. The second skill (`breeding-trial-management`) handles the breeding program logistics — trial design, field books, germplasm databases, crossing plans, selection indices. Both are CLI-first, example-driven, and follow the existing `claude-scientific-skills` repo conventions.

## Core Value

Researchers and breeders worldwide can run the full breeding genomics pipeline — from genotype QC through genomic prediction to crossing decisions — using only free, open-source tools, with examples that actually work out of the box.

## Requirements

### Validated

- ✓ QTL analysis SKILL.md with 19 working examples — existing
- ✓ Unified CLI (qtl_cli.py) with 6 subcommands — existing
- ✓ QTLmax feature parity research (research/qtlmax_features_analysis.md) — existing
- ✓ Cross-validation guide with Python + R code — existing (research/ only)
- ✓ Tool comparison reference doc — existing
- ✓ Apache-2.0 licensing, Clayton Young attribution — existing

### Active

**QTL Analysis Expansion (~15 new examples/enhancements):**

- [ ] Bayesian GP models (BayesA/B/Cpi) via BGLR
- [ ] Multi-trait GWAS (joint models) via sommer/GAPIT
- [ ] G×E / MET models (genotype-by-environment) via sommer/BGData
- [ ] Elastic Net CV for SNP selection via glmnet/scikit-learn
- [ ] Pedigree kinship (NRM) via AGHmatrix/nadiv
- [ ] Genomic kinship (0-2 NRM scale) via AGHmatrix
- [ ] SNP annotation via snpEff/VEP
- [ ] Genomic control (lambda adjustment) as standalone example
- [ ] Genetic similarity matrix (explicit, beyond kinship)
- [ ] Permutation/FDR threshold correction (BH, Bonferroni, permutation)
- [ ] Sample QC (heterozygosity outliers, sex check, call rate)
- [ ] Rare variant tests (burden/SKAT)
- [ ] Cross-validation promoted to full runnable example
- [ ] Report generation (PDF/HTML from analysis results)
- [ ] Enhanced covariate handling (beyond PCs — sex, batch, environment)

**New Breeding Trial Management Skill:**

- [ ] SKILL.md with overview, installation, CLI reference
- [ ] Trial design — RCBD, alpha lattice, augmented design
- [ ] Field book generation — plot layouts, QR/barcode labels
- [ ] Germplasm database integration — Breedbase/BMS API client
- [ ] Selection indices — multi-trait economic weighting
- [ ] Crossing plans / mate allocation — optimal contributions
- [ ] Phenotype data import — field device → DB pipeline
- [ ] Pedigree management — track lineage, inbreeding coefficients
- [ ] Breeding value ranking — GEBV-based selection lists
- [ ] Unified CLI for breeding operations
- [ ] Preflight system checker (like qtl-analysis check_system.py)

**Shared Infrastructure:**

- [ ] Both skills use synthetic data by default + optional real public datasets
- [ ] Consistent attribution (Clayton Young, Apache-2.0, borealBytes)
- [ ] Each example follows Input→Process→Output with committed outputs
- [ ] READMEs at accessible level (high school / early college)

### Out of Scope

- Interactive web dashboards (Shiny/Streamlit) — skill is CLI-first
- Genome browser hosting (JBrowse/IGV) — infrastructure, not analysis
- Mobile field apps — out of scope for agent skills
- Proprietary tool wrappers (ASReml) — open-source only
- Re-implementing sklearn/matplotlib — reference existing repo skills instead
- Cloud compute integration — skills run locally

## Context

- This lives in the `claude-scientific-skills` repo (K-Dense-AI/claude-scientific-skills)
- The repo already has 144 skills across all scientific domains
- Skills follow the Agent Skills standard (agentskills.io)
- Each skill has SKILL.md, examples/, references/, scripts/
- The QTL skill already has 19 examples, a CLI, and extensive research docs
- Target users: plant/animal breeders, genomics researchers, ag-biotech companies
- Competing with: QTLmax ($500-2000/yr), Golden Helix SVS, Agronomix/Genovix
- OSS backends: PLINK 2, GEMMA, tensorQTL, R/qtl2, GAPIT, BGLR, rrBLUP, sommer, AGHmatrix, agricolae, Breedbase

## Constraints

- **License**: Apache-2.0 for all skill code; backends have their own licenses (GPL-3 for PLINK/GEMMA/R packages)
- **Attribution**: All files must have Clayton Young / borealBytes / Superior Byte Works LLC headers
- **Data**: Synthetic data must be self-contained (no network calls). Real datasets are optional extras.
- **Dependencies**: Python-native where possible, R via Rscript wrappers, no Java required for core features
- **Repo conventions**: Must match existing skill directory structure (SKILL.md, examples/, references/, scripts/)
- **Branch**: Work on `feat/qtl-analysis-skill` branch

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Two skills, not one | Trial management is logistics; QTL is analysis. Different users, different tools. | — Pending |
| Synthetic + real data | Self-contained examples always work; real data for credibility | — Pending |
| Python-first, R wrappers | Most breeding tools are R; skill orchestrates via Rscript | — Pending |
| No ASReml | Commercial license incompatible with open-source goal | ✓ Good |
| Reference existing skills | Don't duplicate scikit-learn, matplotlib, etc. already in repo | ✓ Good |

---
*Last updated: 2026-02-23 after initialization*
