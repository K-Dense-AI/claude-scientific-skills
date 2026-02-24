# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-23)

**Core value:** Full breeding genomics pipeline using only free, open-source tools
**Current focus:** Phase 2 — Advanced GWAS

## Current Milestone

**v1.0 — Full QTLmax parity + breeding trial management**

| Phase | Status | Description |
|-------|--------|-------------|
| 1 | ✅ Complete | Bayesian & Advanced GP |
| 2 | ✅ Complete | Advanced GWAS |
| 3 | ○ Pending | Kinship & Relatedness |
| 4 | ○ Pending | QC & Annotation |
| 5 | ○ Pending | Report Generation |
| 6 | ○ Pending | QTL Infrastructure |
| 7 | ○ Pending | Breeding Core |
| 8 | ○ Pending | Trial Design |
| 9 | ○ Pending | Germplasm & Pedigree |
| 10 | ○ Pending | Selection & Crossing |
| 11 | ○ Pending | Breeding Infrastructure |
| 12 | ○ Pending | Attribution & Quality |

Progress: █████░░░░░░░░ 26% (10/38 requirements)

## Completed in Phase 1

- ✅ **QTLGP-01**: Bayesian GP example (bayesian-gp/) — BayesA/B/Cpi/GBLUP comparison
- ✅ **QTLGP-02**: Elastic Net CV example (elastic-net-cv/) — SNP selection with lambda optimization
- ✅ **QTLGP-03**: Cross-validation promotion (cross-validation/) — 5 CV strategies
- ✅ **QTLGP-04**: G×E prediction example (gxe-prediction/) — multi-environment GP

## Phase 1 Verification

**Status:** PASSED
**Date:** 2026-02-23
**Requirements:** 4/4 satisfied

All Phase 1 deliverables verified:
- All scripts run successfully with committed outputs
- READMEs follow Input→Process→Output format
- All requirements checked off in REQUIREMENTS.md

## Accumulated Context

### Key Decisions
- Two skills: qtl-analysis (expand) + breeding-trial-management (new)
- Synthetic data by default, real datasets optional
- Python-first, R wrappers for BGLR/sommer/agricolae
- Apache-2.0, Clayton Young attribution
- No ASReml (commercial)

### Current Work
- Phase 3: Kinship & Relatedness planning (QTLKN-01..03)

### Roadmap Evolution
- Phase 02.1 inserted after Phase 2: Progressive PR #55 comment updates after each successful image or major output (URGENT)

### Known Issues
- qtl_cli.py line 197: `np` not defined (missing numpy import in manhattan plot function)

## Session Continuity

Last session: 2026-02-24T00:00:00Z
Stopped at: Phase 2 complete, ready to plan Phase 3
Resume file: 02-01-SUMMARY.md, 02-02-SUMMARY.md, 02-03-SUMMARY.md

---
*Last updated: 2026-02-24 after Phase 2 completion*
