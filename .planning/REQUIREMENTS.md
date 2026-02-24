# Requirements: Meta-QTLmax Breeding Genomics Skills

**Defined:** 2026-02-23
**Core Value:** Researchers and breeders can run the full breeding genomics pipeline using only free, open-source tools, with examples that work out of the box.

## v1 Requirements

Requirements for initial release. Each maps to roadmap phases.

### QTL Expansion — Bayesian & Advanced GP

- [x] **QTLGP-01**: Bayesian GP example (BayesA/B/Cpi) using BGLR with synthetic data, prediction accuracy plots
- [x] **QTLGP-02**: Elastic Net CV example for SNP selection with lambda optimization plot
- [x] **QTLGP-03**: Cross-validation promoted from research/ to full runnable example with committed outputs
- [x] **QTLGP-04**: Multi-environment GP (G×E) example with reaction norms via sommer/BGData

### QTL Expansion — Advanced GWAS

- [x] **QTLGW-01**: Multi-trait GWAS example (joint model across 2+ traits) via sommer or GAPIT
- [x] **QTLGW-02**: G×E GWAS / MET model example with environment-specific effects
- [x] **QTLGW-03**: Enhanced covariate handling example (PCs + sex + batch + environment)
- [x] **QTLGW-04**: Permutation and FDR threshold correction example (BH, Bonferroni, permutation-based)
- [x] **QTLGW-05**: Genomic control (lambda adjustment) as standalone example with before/after QQ plots
- [x] **QTLGW-06**: Rare variant association tests (burden/SKAT) example

### QTL Expansion — Kinship & Relatedness

- [x] **QTLKN-01**: Pedigree kinship (NRM) example from pedigree records via AGHmatrix/nadiv
- [x] **QTLKN-02**: Genomic kinship (0-2 NRM scale) example via AGHmatrix
- [x] **QTLKN-03**: Genetic similarity matrix example (IBS/distance-based) via PLINK --distance

### QTL Expansion — QC & Annotation

- [x] **QTLQC-01**: Sample QC example (heterozygosity outliers, sex check, call rate) via PLINK 2
- [x] **QTLQC-02**: SNP annotation example via snpEff (variant → gene → effect → impact)

### QTL Expansion — Reporting

- [ ] **QTLRP-01**: Report generation — PDF/HTML summary of QC + GWAS + GP results from a single command

### QTL Expansion — Infrastructure

- [ ] **QTLIN-01**: Update SKILL.md to document all new examples and expanded CLI
- [ ] **QTLIN-02**: Update qtl_cli.py with new subcommands for added analyses
- [ ] **QTLIN-03**: Add optional real public dataset examples (Ames panel or NAM as downloads)

### Breeding Trial Management — Core Skill

- [ ] **BREED-01**: SKILL.md with overview, installation, CLI reference, tool selection guide
- [ ] **BREED-02**: Preflight system checker (scripts/check_system.py) for breeding dependencies
- [ ] **BREED-03**: Unified CLI (scripts/breeding_cli.py) with subcommands

### Breeding Trial Management — Trial Design

- [ ] **TRIAL-01**: RCBD (Randomized Complete Block Design) example with plot layout visualization
- [ ] **TRIAL-02**: Alpha-lattice design example for large trials
- [ ] **TRIAL-03**: Augmented design example (unreplicated entries + replicated checks)
- [ ] **TRIAL-04**: Field book generation example — plot maps, QR/barcode labels

### Breeding Trial Management — Germplasm & Pedigree

- [ ] **GERM-01**: Germplasm database integration — Breedbase API client example (read/write accessions)
- [ ] **GERM-02**: Pedigree management example — lineage tracking, inbreeding coefficient calculation
- [ ] **GERM-03**: BMS (Breeding Management System) API client example

### Breeding Trial Management — Selection & Crossing

- [ ] **SEL-01**: Selection index example — multi-trait economic weighting (Smith-Hazel, base index)
- [ ] **SEL-02**: Breeding value ranking example — GEBV-based selection lists with confidence intervals
- [ ] **SEL-03**: Crossing plan / mate allocation example — optimal contributions, avoid inbreeding
- [ ] **SEL-04**: Phenotype data import pipeline — CSV/field device → standardized format

### Breeding Trial Management — Infrastructure

- [ ] **BRDIN-01**: References directory with API docs, data format specs, design theory
- [ ] **BRDIN-02**: All examples follow Input→Process→Output with committed outputs and accessible READMEs

### Shared — Attribution & Quality

- [ ] **SHARED-01**: All files have Apache-2.0 headers with Clayton Young / borealBytes / Superior Byte Works LLC
- [ ] **SHARED-02**: All examples generate synthetic data by default (no network calls required)
- [ ] **SHARED-03**: All READMEs written at high school / early college accessibility level

## v2 Requirements

Deferred to future release. Tracked but not in current roadmap.

### Advanced Features

- **ADV-01**: Multi-trait GP with deep learning (multi-task neural nets)
- **ADV-02**: Genome-wide epistasis detection
- **ADV-03**: Haplotype-based GWAS (beyond single-SNP)
- **ADV-04**: Copy number variation (CNV) analysis integration
- **ADV-05**: Interactive Streamlit/Shiny dashboards for trial results
- **ADV-06**: Real-time field data sync (IoT → Breedbase)
- **ADV-07**: Simulation engine (AlphaSimR-based breeding program simulation)

## Out of Scope

| Feature | Reason |
|---------|--------|
| Genome browser hosting (JBrowse/IGV) | Infrastructure, not analysis — separate concern |
| Mobile field apps | Out of scope for agent skills |
| ASReml integration | Commercial license, open-source only |
| Re-implementing sklearn/matplotlib | Already exist as separate skills in this repo |
| Cloud compute / HPC integration | Skills run locally by design |
| Wet lab protocols | Laboratory work, not computational |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| QTLGP-01 | Phase 1 | Complete |
| QTLGP-02 | Phase 1 | Complete |
| QTLGP-03 | Phase 1 | Complete |
| QTLGP-04 | Phase 1 | Complete |
| QTLGW-01 | Phase 2 | Complete |
| QTLGW-02 | Phase 2 | Complete |
| QTLGW-03 | Phase 2 | Complete |
| QTLGW-04 | Phase 2 | Complete |
| QTLGW-05 | Phase 2 | Complete |
| QTLGW-06 | Phase 2 | Complete |
| QTLKN-01 | Phase 3 | Complete |
| QTLKN-02 | Phase 3 | Complete |
| QTLKN-03 | Phase 3 | Complete |
| QTLQC-01 | Phase 4 | Complete |
| QTLQC-02 | Phase 4 | Complete |
| QTLRP-01 | Phase 5 | Pending |
| QTLIN-01 | Phase 6 | Pending |
| QTLIN-02 | Phase 6 | Pending |
| QTLIN-03 | Phase 6 | Pending |
| BREED-01 | Phase 7 | Pending |
| BREED-02 | Phase 7 | Pending |
| BREED-03 | Phase 7 | Pending |
| TRIAL-01 | Phase 8 | Pending |
| TRIAL-02 | Phase 8 | Pending |
| TRIAL-03 | Phase 8 | Pending |
| TRIAL-04 | Phase 8 | Pending |
| GERM-01 | Phase 9 | Pending |
| GERM-02 | Phase 9 | Pending |
| GERM-03 | Phase 9 | Pending |
| SEL-01 | Phase 10 | Pending |
| SEL-02 | Phase 10 | Pending |
| SEL-03 | Phase 10 | Pending |
| SEL-04 | Phase 10 | Pending |
| BRDIN-01 | Phase 11 | Pending |
| BRDIN-02 | Phase 11 | Pending |
| SHARED-01 | Phase 12 | Pending |
| SHARED-02 | Phase 12 | Pending |
| SHARED-03 | Phase 12 | Pending |

**Coverage:**
- v1 requirements: 38 total
- Mapped to phases: 38
- Unmapped: 0 ✓

---
*Requirements defined: 2026-02-23*
*Last updated: 2026-02-23 after initial definition*
