# POD Systematic Review Pipeline - Project Summary

## Overview

This repository contains a **complete, reproducible Python pipeline** for conducting systematic reviews and meta-analyses of **postoperative delirium (POD) in elderly patients undergoing total hip/knee arthroplasty**.

**Repository:** `pod-systematic-review/`
**Branch:** `claude/systematic-review-delirium-7x5sg`
**Version:** 0.1.0
**Python:** 3.11+
**License:** MIT

## Key Features

### ✅ Literature Retrieval
- **Automated PubMed retrieval** via NCBI E-utilities API with configurable queries
- **Import support** for Embase, CINAHL, Web of Science, Cochrane (RIS/BibTeX/CSV formats)
- Handles API rate limiting and retries
- Exports to multiple formats (JSONL, RIS, CSV) for compatibility with EndNote/Covidence

### ✅ Deduplication
- **Exact matching** on DOI and PMID
- **Fuzzy matching** on title (90% threshold) and first author (85% threshold)
- Year validation for fuzzy matches
- Detailed deduplication report with match type statistics

### ✅ Screening
- **Streamlit web application** for citation screening
- **Dual independent review** workflow
- **Conflict adjudication** by third reviewer
- Title/Abstract and Full-Text screening phases
- 12 exclusion reason codes aligned with protocol
- SQLite database for persistence and audit trail
- Real-time progress statistics

### ✅ Data Extraction
- **Structured CSV templates** with validation
- **JSON schema** for programmatic validation
- Captures study characteristics, POD outcomes, and risk factors
- Support for multiple risk factors per study
- Validation functions for completeness checks

### ✅ Risk of Bias Assessment
- **Newcastle-Ottawa Scale (NOS)** for observational studies (0-9 stars)
- **RoB2 tool** for RCTs (5 domains)
- Automated quality rating (Good/Fair/Poor)
- Summary table generation

### ✅ Meta-Analysis
- **Random-effects model** (DerSimonian-Laird method)
- Log-transformation of OR/RR/HR
- Standard error calculation from 95% CI
- **Heterogeneity assessment:** Cochran's Q, I², τ²
- **Forest plots** with study weights and pooled diamond
- **Funnel plots** for publication bias (if k≥10)
- **Egger's test** for funnel plot asymmetry
- **Sensitivity analysis:** Leave-one-out
- **Subgroup analysis:** By surgery type, region, diagnostic method, etc.

### ✅ PRISMA Reporting
- **Automated PRISMA 2020 flow diagram** generation
- Counts table for identification, screening, exclusion, inclusion stages
- Publication-ready figures (PNG, 300 DPI)

### ✅ Reproducibility
- Version-controlled configuration (config.yaml)
- Full audit trail (timestamps, reviewer IDs)
- Makefile for standardized workflow
- Unit tests with pytest
- Type hints and documentation throughout
- Example data for testing without API access

## Repository Structure

```
pod-systematic-review/
├── README.md                # Comprehensive user guide
├── INSTALL.md              # Installation instructions
├── PROJECT_SUMMARY.md      # This file
├── LICENSE                 # MIT License
├── Makefile                # Workflow automation
├── pyproject.toml          # Python package configuration
├── config.yaml             # Review protocol configuration
├── .env.example            # Environment variables template
├── .gitignore              # Git ignore rules
│
├── src/pod_review/         # Source code
│   ├── __init__.py
│   ├── cli.py              # Command-line interface (Click)
│   ├── config.py           # Configuration loader (Pydantic)
│   ├── models.py           # Data models (Pydantic)
│   ├── retrieval.py        # PubMed API (BioPython Entrez)
│   ├── importers.py        # RIS/BibTeX/CSV import
│   ├── exporters.py        # Export to RIS/CSV
│   ├── deduplication.py    # Fuzzy matching (thefuzz)
│   ├── extraction.py       # Data extraction templates
│   ├── risk_of_bias.py     # NOS and RoB2 tools
│   ├── prisma.py           # PRISMA flow diagrams
│   ├── screening/
│   │   ├── __init__.py
│   │   ├── database.py     # SQLite screening database
│   │   └── app.py          # Streamlit screening UI
│   └── meta/
│       ├── __init__.py
│       ├── analysis.py     # Meta-analysis (NumPy, SciPy)
│       └── plots.py        # Forest/funnel plots (Matplotlib)
│
├── tests/                  # Unit tests
│   ├── __init__.py
│   ├── test_deduplication.py
│   └── test_models.py
│
├── examples/               # Example data
│   └── example_pubmed.ris  # Sample RIS file
│
├── data/                   # Data files (gitignored)
│   ├── raw/               # Retrieved records
│   ├── imported/          # Imported records
│   └── processed/         # Deduplicated records
│
└── outputs/               # Analysis outputs (gitignored)
    ├── prisma/            # PRISMA diagrams
    ├── extraction/        # Extraction templates
    ├── rob/               # Risk of bias summaries
    └── meta/              # Meta-analysis results
```

## Core Modules

### 1. Configuration (`config.py`, `config.yaml`)
- **Pydantic Settings** for environment variables (.env)
- **YAML configuration** for review protocol
- PICOS criteria, search strategies, screening settings, extraction schema

### 2. Data Models (`models.py`)
- **Citation:** Full bibliographic metadata
- **ScreeningRecord:** Reviewer decisions with timestamps
- **ExtractionRecord:** Study characteristics + risk factors
- **RiskFactor:** Effect sizes with CI → SE conversion
- **NOSAssessment:** Quality rating calculation
- **MetaAnalysisResult:** Pooled effects + heterogeneity

### 3. Retrieval & Import (`retrieval.py`, `importers.py`, `exporters.py`)
- **PubMedRetriever:** NCBI Entrez API wrapper with rate limiting
- **RISImporter, BibTeXImporter, CSVImporter:** Parse multiple formats
- **Exporters:** Generate RIS/CSV/BibTeX for external tools

### 4. Deduplication (`deduplication.py`)
- **Exact matching:** DOI, PMID
- **Fuzzy matching:** Levenshtein similarity on title/authors
- **Deduplication report:** Match statistics by type

### 5. Screening (`screening/database.py`, `screening/app.py`)
- **SQLite database:** Citations, screening records, progress tracking
- **Streamlit app:** Login, screening, adjudication, statistics
- **Conflict detection:** Automatic flagging of reviewer disagreements

### 6. Extraction (`extraction.py`)
- **Template generation:** CSV with example rows
- **CSV loader:** Parse extraction data with validation
- **JSON schema:** Validation rules
- **Validation functions:** Check completeness, data integrity

### 7. Risk of Bias (`risk_of_bias.py`)
- **NOS template:** 8 items → 0-9 stars → Good/Fair/Poor
- **RoB2 template:** 5 domains → Low/Some concerns/High
- **Summary tables:** Publication-ready CSV

### 8. Meta-Analysis (`meta/analysis.py`, `meta/plots.py`)
- **MetaAnalyzer:** Random-effects pooling (DerSimonian-Laird)
- **Heterogeneity:** Q, I², τ²
- **Publication bias:** Egger's test, funnel plot
- **Sensitivity:** Leave-one-out
- **Subgroup:** Stratified analysis
- **Plots:** Forest plot, funnel plot, sensitivity plot (Matplotlib)

### 9. PRISMA (`prisma.py`)
- **Flow diagram:** PRISMA 2020 compliant (Matplotlib)
- **Counts table:** CSV export

### 10. CLI (`cli.py`)
- **Commands:** retrieve, import, dedup, screen, prisma, extract, meta, rob-summary
- **Click framework:** User-friendly command-line interface

## Clinical Protocol

**Research Question:** What are the risk factors for postoperative delirium in elderly patients (≥65 years) undergoing total hip or knee arthroplasty?

**Inclusion Criteria:**
- Age ≥65 years (or mixed with elderly ≥50% OR subgroup analysis)
- Primary or revision THA/TKA
- POD as outcome with diagnostic method specified (DSM-5, CAM, etc.)
- At least one risk factor examined with statistical analysis
- Cohort, case-control, or secondary RCT analysis
- English or Chinese language
- Peer-reviewed journal

**Exclusion Criteria:**
- Wrong population/intervention/outcome
- Case reports, reviews, opinions
- Grey literature (theses, preprints)
- Insufficient data

**Risk Factors of Interest:**
- Demographics: age, sex, education
- Comorbidities: cardiovascular, respiratory, metabolic, neurological
- Pre-op cognition: MMSE, MoCA, cognitive impairment
- Medications: psychotropics, sedatives, anticholinergics, polypharmacy
- Surgery: type, duration, blood loss, transfusion
- Anesthesia: general vs regional, depth, duration
- Postoperative: ICU, pain management, complications

**Databases:**
- **Automated:** PubMed
- **Manual import:** Embase, CINAHL, Web of Science, Cochrane Library

**Screening:** Dual independent review + adjudication

**Quality Assessment:** NOS (observational), RoB2 (RCTs)

**Meta-Analysis:** Random-effects, publication bias if k≥10, sensitivity & subgroup analyses

## Statistical Methods

### Random-Effects Meta-Analysis (DerSimonian-Laird)

1. **Log-transform** effect sizes: `log_OR = log(OR)`
2. **Calculate SE** from 95% CI: `SE = (log(upper) - log(lower)) / (2 × 1.96)`
3. **Fixed-effect weights:** `w_FE = 1 / variance`
4. **Q statistic:** `Q = Σ w_FE × (effect - pooled_FE)²`
5. **Between-study variance:** `τ² = max(0, (Q - (k-1)) / C)` where `C = Σw - Σw²/Σw`
6. **I² statistic:** `I² = max(0, ((Q - (k-1)) / Q) × 100%)`
7. **Random-effects weights:** `w_RE = 1 / (variance + τ²)`
8. **Pooled effect:** `pooled = Σ(w_RE × effect) / Σw_RE`
9. **Pooled SE:** `SE_pooled = sqrt(1 / Σw_RE)`
10. **95% CI:** `[pooled - 1.96×SE, pooled + 1.96×SE]`

### Publication Bias (if k ≥ 10)

- **Egger's test:** Linear regression of effect on precision (1/SE)
- **Funnel plot:** Effect size vs. SE (inverted y-axis)
- Asymmetry suggests potential publication bias

### Sensitivity Analysis

- **Leave-one-out:** Re-run meta-analysis k times, excluding one study each
- **Exclude high RoB:** Remove studies with high risk of bias
- **Subgroup analysis:** Stratify by surgery type, region, diagnostic method

## Example Workflow

```bash
# 1. Setup
cp .env.example .env
nano .env  # Add NCBI_EMAIL and NCBI_API_KEY
pip install -e ".[dev]"

# 2. Retrieve from PubMed
pod-review retrieve --max-results 10000

# 3. Import from other databases (manual export first)
pod-review import-records embase_export.ris
pod-review import-records cinahl_export.ris
pod-review import-records wos_export.bib

# 4. Deduplicate
pod-review dedup

# 5. Load for screening
pod-review load-screening --input-file data/processed/deduplicated.jsonl

# 6. Screen citations (Streamlit app)
streamlit run src/pod_review/screening/app.py
# - Login as R1, screen title/abstracts
# - Login as R2, screen title/abstracts
# - Login as R3 (adjudicator), resolve conflicts
# - Repeat for full-text phase

# 7. Generate PRISMA
pod-review prisma

# 8. Data extraction
pod-review create-templates
# - Open outputs/extraction/extraction_template.csv
# - Fill in study data and risk factors
# - Save as extraction_completed.csv

# 9. Risk of bias
# - Open outputs/extraction/nos_template.csv
# - Score each study
# - Save as nos_completed.csv
pod-review rob-summary --nos-file outputs/extraction/nos_completed.csv

# 10. Meta-analysis
pod-review meta-analyze \
  outputs/extraction/extraction_completed.csv \
  "Age (per year)" \
  --effect-type OR

# Outputs:
# - outputs/meta/Age_forest.png
# - outputs/meta/Age_funnel.png
# - outputs/meta/Age_sensitivity.png
# - Console: pooled OR, CI, p-value, I², τ²

# 11. Repeat for other risk factors
pod-review meta-analyze extraction_completed.csv "Cognitive impairment" --effect-type OR
pod-review meta-analyze extraction_completed.csv "Surgery duration" --effect-type OR
# etc.
```

## Testing

```bash
# Run all tests
pytest tests/ -v

# With coverage
pytest tests/ --cov=pod_review --cov-report=term-missing

# Test specific module
pytest tests/test_deduplication.py -v
```

**Test Coverage:**
- Deduplication: exact matching (DOI/PMID), fuzzy matching, edge cases
- Models: SE calculation, log transformation, NOS scoring, quality rating

## Dependencies

**Core:**
- biopython (PubMed API)
- requests (HTTP)
- pandas, numpy, scipy (data analysis)
- matplotlib, seaborn (plotting)
- streamlit (web UI)
- pydantic, pydantic-settings (configuration)
- pyyaml (config files)
- click (CLI)
- thefuzz, python-levenshtein (fuzzy matching)
- rispy, bibtexparser (import)

**Development:**
- pytest, pytest-cov (testing)
- mypy (type checking)
- black, ruff (code formatting)

## Future Enhancements

Potential extensions (not yet implemented):

1. **Additional meta-analysis methods:** Paule-Mandel, REML, Hartung-Knapp adjustment
2. **Meta-regression:** Explore sources of heterogeneity
3. **Network meta-analysis:** Compare multiple interventions
4. **Dose-response meta-analysis:** For continuous exposures
5. **Individual patient data (IPD) meta-analysis**
6. **Machine learning screening:** Prioritize citations for screening
7. **Full-text PDF extraction:** Automated data extraction from PDFs
8. **Interactive dashboards:** Plotly/Dash for exploration
9. **API endpoints:** REST API for integration
10. **Docker container:** For easier deployment

## Support

**Documentation:**
- `README.md` - Comprehensive user guide
- `INSTALL.md` - Installation and quick start
- Inline docstrings in all modules
- Type hints for all functions

**Examples:**
- `examples/example_pubmed.ris` - Sample RIS file
- Unit tests demonstrate usage patterns

**Logs:**
- Console logging throughout pipeline
- Database audit trail for screening decisions

## Citation

If you use this pipeline in your research, please cite:

```
[Your systematic review paper once published]

Software:
POD Systematic Review Pipeline, v0.1.0
URL: https://github.com/Circle8211/claude-scientific-skills
Branch: claude/systematic-review-delirium-7x5sg
```

## License

MIT License - see LICENSE file

---

**Created:** 2026-01-07
**Version:** 0.1.0
**Python:** 3.11+
**Authors:** Research Team
**Contact:** [your-email]
