# POD Systematic Review Pipeline

A comprehensive, reproducible systematic review and meta-analysis pipeline for investigating **risk factors for postoperative delirium (POD) in elderly patients (≥65 years) undergoing total hip/knee arthroplasty**.

## Overview

This Python-based pipeline operationalizes a complete systematic review workflow following PRISMA 2020 guidelines, from literature retrieval through meta-analysis and reporting.

### Key Features

- **Automated retrieval** from PubMed via NCBI E-utilities API
- **Flexible import** from Embase, CINAHL, Web of Science, Cochrane (RIS/BibTeX/CSV formats)
- **Intelligent deduplication** using DOI/PMID exact matching + fuzzy title/author matching
- **Web-based screening app** (Streamlit) with dual independent review and adjudication
- **Structured data extraction** with validation
- **Risk of bias assessment** (NOS for observational studies, RoB2 for RCTs)
- **Random-effects meta-analysis** (DerSimonian-Laird) with forest plots, funnel plots, heterogeneity statistics
- **Publication bias testing** (Egger's test) when k≥10
- **Sensitivity & subgroup analyses**
- **PRISMA flow diagram** auto-generation
- **Full audit trail** with SQLite databases

## Protocol Summary

**Population:** Elderly patients ≥65 years undergoing THA/TKA (primary or revision)

**Outcome:** POD incidence with reported diagnostic method (DSM-5, CAM, etc.) and statistical associations with risk factors

**Exposures:** Demographics, comorbidities, pre-op cognition, medications, surgery/anesthesia factors, postoperative factors

**Study Designs:** Prospective/retrospective cohort, case-control, secondary RCT analysis

**Languages:** English + Chinese

**Databases:** PubMed (automated) + Embase/CINAHL/WoS/Cochrane (manual import)

## Installation

### Prerequisites

- Python 3.11 or higher
- pip package manager
- Git

### Setup

1. **Clone the repository**

```bash
git clone <repository-url>
cd pod-systematic-review
```

2. **Install dependencies**

```bash
pip install -e ".[dev]"
```

Or using the Makefile:

```bash
make install
```

3. **Configure environment variables**

```bash
cp .env.example .env
```

Edit `.env` and add your NCBI API credentials:

```bash
NCBI_EMAIL=your_email@example.com
NCBI_API_KEY=your_ncbi_api_key  # Get from https://www.ncbi.nlm.nih.gov/account/
```

**Note:** An NCBI API key is **highly recommended** (increases rate limit from 3 to 10 requests/second). Without it, PubMed retrieval will be slower.

## Quick Start

### End-to-End Workflow

```bash
# 1. Retrieve from PubMed
make retrieve

# 2. (Optional) Import from other databases
pod-review import-records data/my_embase_export.ris --format ris

# 3. Deduplicate
make dedup

# 4. Load into screening database
make load-screening

# 5. Launch screening app
make screen
# Opens at http://localhost:8501

# 6. Generate PRISMA diagram
make prisma

# 7. Create extraction templates
make extract

# 8. After completing extraction, run meta-analysis
pod-review meta-analyze outputs/extraction/extraction_completed.csv "Age" --effect-type OR

# 9. Run tests
make test
```

## Detailed Usage

### 1. Literature Retrieval

#### PubMed (Automated)

The configured PubMed query is defined in `config.yaml` and includes:
- Keywords: postoperative delirium, arthroplasty, elderly, risk factors
- Field tags: [Title/Abstract], [MeSH Terms]
- Language filters: English OR Chinese

```bash
pod-review retrieve --max-results 10000 --output-dir data/raw
```

**Outputs:**
- `data/raw/pubmed_records.jsonl` - Machine-readable format
- `data/raw/pubmed_records.ris` - For EndNote/Covidence
- `data/raw/pubmed_records.csv` - For Excel/review

#### Other Databases (Manual + Import)

For databases without open APIs (Embase, CINAHL, Web of Science, Cochrane):

1. Perform manual search using query strings in `config.yaml`
2. Export results to RIS, BibTeX, or CSV
3. Import into pipeline:

```bash
# RIS format
pod-review import-records my_embase_export.ris --format ris --output-dir data/imported

# BibTeX format
pod-review import-records my_wos_export.bib --format bibtex --output-dir data/imported

# CSV format (auto-detected from extension)
pod-review import-records my_cinahl_export.csv --output-dir data/imported
```

### 2. Deduplication

Combines all imported sources and removes duplicates using:
- **Exact matching** on DOI and PMID
- **Fuzzy matching** on title (90% similarity) and first author (85% similarity)
- **Year validation** for fuzzy matches

```bash
pod-review dedup --input-dir data --output-dir data/processed
```

**Outputs:**
- `data/processed/deduplicated.jsonl`
- `data/processed/deduplicated.ris`
- `data/processed/deduplicated.csv`
- Console report with duplicate statistics

### 3. Screening

#### Load Citations

```bash
pod-review load-screening --input-file data/processed/deduplicated.jsonl
```

Creates SQLite database at `data/screening.db`

#### Launch Screening App

```bash
streamlit run src/pod_review/screening/app.py
# Or: make screen
```

**Streamlit Interface Features:**

1. **Reviewer Login** - Select reviewer ID
2. **Citation Screening**
   - Title/Abstract phase
   - Full-text phase
   - View title, authors, journal, abstract
   - Make decision: Include/Exclude/Unclear
   - Select exclusion reason (E1-E12)
   - Add notes
3. **Adjudication** - Resolve conflicts between reviewers
4. **Statistics** - View screening progress

**Workflow:**
- Two independent reviewers screen each citation
- Conflicts automatically flagged
- Third reviewer adjudicates conflicts
- Full audit trail with timestamps

### 4. PRISMA Flow Diagram

```bash
pod-review prisma --output-dir outputs/prisma
```

**Outputs:**
- `outputs/prisma/prisma_flow.png` - PRISMA 2020 flow diagram
- `outputs/prisma/prisma_counts.csv` - Count table for manuscript

### 5. Data Extraction

#### Create Templates

```bash
pod-review create-templates --output-dir outputs/extraction
```

**Templates created:**
- `extraction_template.csv` - Study characteristics and risk factors
- `nos_template.csv` - Newcastle-Ottawa Scale for observational studies
- `rob2_template.csv` - RoB2 for RCTs

#### Extraction Workflow

1. Open `extraction_template.csv` in Excel/LibreOffice
2. For each included study, fill in:
   - Study characteristics (design, country, sample size, age, sex, surgery type)
   - POD outcome (definition, tool, timing, incidence)
   - Risk factors (name, category, effect type, OR/RR/HR, CI, p-value, adjustment, covariates)
3. For multiple risk factors per study, add additional rows with same `citation_id`
4. Save completed extraction file

#### Risk of Bias Assessment

1. Open `nos_template.csv` or `rob2_template.csv`
2. Score each study according to tool domains
3. NOS: 0-9 stars (Selection, Comparability, Outcome)
4. RoB2: Low/Some concerns/High for each domain

Generate summary tables:

```bash
pod-review rob-summary --nos-file outputs/extraction/nos_completed.csv --output-dir outputs/rob
```

### 6. Meta-Analysis

```bash
pod-review meta-analyze \
  outputs/extraction/extraction_completed.csv \
  "Age (per year)" \
  --effect-type OR \
  --output-dir outputs/meta
```

**Outputs:**
- `outputs/meta/Age_forest.png` - Forest plot
- `outputs/meta/Age_funnel.png` - Funnel plot (if k≥10)
- `outputs/meta/Age_sensitivity.png` - Leave-one-out sensitivity analysis
- Console output with pooled effect, CI, heterogeneity (I², τ², Q)

**Meta-Analysis Methods:**
- Random-effects model (DerSimonian-Laird)
- Log-transformed effect sizes
- 95% confidence intervals
- Heterogeneity: Cochran's Q, I², τ²
- Publication bias: Egger's test (if k≥10), funnel plot

**Sensitivity Analysis:**

```bash
# Leave-one-out analysis is automatic
# For subgroup analysis (requires implementation):
# pod-review meta-analyze <file> <factor> --subgroup surgery_type
```

### 7. Additional Commands

**Generate Rob summary:**

```bash
pod-review rob-summary \
  --nos-file outputs/extraction/nos_completed.csv \
  --rob2-file outputs/extraction/rob2_completed.csv \
  --output-dir outputs/rob
```

## Configuration

### config.yaml

The `config.yaml` file defines all review parameters:

- **Review metadata** - Title, version, reviewers
- **PICOS criteria** - Population, outcomes, exposures, study designs
- **Inclusion/exclusion criteria**
- **Search strategies** - Database-specific queries
- **Screening settings** - Phases, decisions, exclusion reasons
- **Extraction schema** - Fields to extract
- **Risk of bias tools** - NOS and RoB2 parameters
- **Meta-analysis settings** - Methods, publication bias thresholds, subgroup variables

Edit this file to customize the review protocol.

## Directory Structure

```
pod-systematic-review/
├── config.yaml              # Main configuration
├── .env                     # Environment variables (not in git)
├── pyproject.toml           # Python package configuration
├── Makefile                 # Convenience commands
├── README.md                # This file
├── LICENSE                  # MIT License
│
├── src/pod_review/          # Source code
│   ├── __init__.py
│   ├── config.py            # Configuration loader
│   ├── models.py            # Data models (Pydantic)
│   ├── retrieval.py         # PubMed retrieval
│   ├── importers.py         # RIS/BibTeX/CSV import
│   ├── exporters.py         # Export to RIS/CSV
│   ├── deduplication.py     # Duplicate detection
│   ├── prisma.py            # PRISMA diagram
│   ├── extraction.py        # Data extraction
│   ├── risk_of_bias.py      # NOS and RoB2
│   ├── cli.py               # Command-line interface
│   ├── screening/           # Screening module
│   │   ├── __init__.py
│   │   ├── database.py      # SQLite database
│   │   └── app.py           # Streamlit app
│   └── meta/                # Meta-analysis module
│       ├── __init__.py
│       ├── analysis.py      # Statistical methods
│       └── plots.py         # Forest/funnel plots
│
├── data/                    # Data files (gitignored)
│   ├── raw/                 # Retrieved records
│   ├── imported/            # Imported records
│   ├── processed/           # Deduplicated records
│   ├── screening.db         # Screening database
│   └── extraction.db        # Extraction database
│
├── outputs/                 # Analysis outputs (gitignored)
│   ├── prisma/              # PRISMA diagrams
│   ├── extraction/          # Templates and extracted data
│   ├── rob/                 # Risk of bias summaries
│   └── meta/                # Meta-analysis results
│
├── tests/                   # Unit tests
│   ├── __init__.py
│   ├── test_deduplication.py
│   └── test_models.py
│
└── examples/                # Example data
    └── example_pubmed.ris
```

## Data Models

### Citation
- **Fields:** ID, title, authors, year, journal, abstract, DOI, PMID, keywords, MeSH terms
- **Sources:** PubMed API, RIS/BibTeX/CSV import

### ScreeningRecord
- **Fields:** Citation ID, reviewer ID, phase, decision, exclusion reason, notes, timestamp
- **Storage:** SQLite database

### ExtractionRecord
- **Fields:** Study characteristics, POD outcome, risk factors (list)
- **Format:** CSV (wide format, repeatable risk factors)

### RiskFactor
- **Fields:** Name, category, effect type (OR/RR/HR), effect size, CI, p-value, adjustment, covariates
- **Methods:** Calculate SE from CI, log-transform effect size

### NOSAssessment
- **Domains:** Selection (4), Comparability (2), Outcome (3)
- **Scoring:** 0-9 stars → Good (7-9), Fair (4-6), Poor (0-3)

### MetaAnalysisResult
- **Outputs:** Pooled effect, CI, p-value, heterogeneity (Q, I², τ²), Egger's test

## Statistical Methods

### Random-Effects Meta-Analysis

**Method:** DerSimonian-Laird (1986)

**Steps:**
1. Log-transform effect sizes (OR/RR/HR)
2. Calculate within-study variance from 95% CI
3. Estimate between-study variance (τ²)
4. Pool effects using random-effects weights: w_i = 1/(σ²_i + τ²)
5. Calculate pooled effect and 95% CI
6. Assess heterogeneity: Q statistic, I², τ²

**Formulas:**

- SE from 95% CI: `SE = (log(upper) - log(lower)) / (2 × 1.96)`
- τ² (DL method): `τ² = max(0, (Q - (k-1)) / C)` where `C = Σw_i - Σw²_i/Σw_i`
- I² = `max(0, ((Q - (k-1)) / Q) × 100%)`

### Publication Bias

**Egger's Test** (if k ≥ 10):
- Linear regression of effect size on precision (1/SE)
- Tests for funnel plot asymmetry
- P < 0.10 suggests potential bias

**Funnel Plot:**
- X-axis: Effect size
- Y-axis: Standard error (inverted)
- Symmetry expected if no publication bias

### Sensitivity Analysis

**Leave-One-Out:**
- Re-run meta-analysis k times, excluding one study each time
- Assess influence of individual studies on pooled estimate

**Subgroup Analysis:**
- Stratify by: surgery type (THA/TKA), region, POD diagnostic method, study design
- Compare pooled effects across subgroups

## Testing

```bash
# Run all tests
pytest tests/ -v

# With coverage
pytest tests/ -v --cov=pod_review --cov-report=term-missing

# Or use Makefile
make test
```

## Reproducibility

### Version Control
- All code is version-controlled with Git
- Configuration files (config.yaml) stored in repository
- Data files (.gitignored) can be regenerated from scripts

### Documentation
- Inline code documentation (docstrings)
- Type hints for all functions
- This README with step-by-step instructions

### Audit Trail
- All screening decisions logged with timestamps in SQLite
- Reviewer IDs tracked for all assessments
- Import dates recorded for all citations

### Makefile Commands
- Standardized workflow commands
- Easy to reproduce entire pipeline:

```bash
make install
make retrieve
make dedup
make screen
# (Complete screening in Streamlit)
make prisma
make extract
# (Complete extraction in CSV)
make meta
```

## Troubleshooting

### PubMed Retrieval Fails

**Problem:** "HTTP Error 429: Too Many Requests"

**Solution:**
- Add NCBI API key to `.env` file
- Reduces rate limit delays
- Get key at: https://www.ncbi.nlm.nih.gov/account/

### Import Fails

**Problem:** "Cannot determine format for file"

**Solution:**
- Specify format explicitly: `--format ris` or `--format bibtex` or `--format csv`
- Check file is valid RIS/BibTeX/CSV format

### Screening App Won't Start

**Problem:** Streamlit not found

**Solution:**
```bash
pip install streamlit
# Or reinstall package
pip install -e ".[dev]"
```

### Meta-Analysis Fails

**Problem:** "Insufficient studies for [risk factor]"

**Solution:**
- Check risk factor name matches exactly (case-sensitive)
- Ensure at least 2 studies report this factor
- Verify effect_type filter if used

### Missing Dependencies

**Problem:** Import errors

**Solution:**
```bash
pip install -e ".[dev]"
# Or install missing package directly
pip install <package-name>
```

## Contributing

This is a research project. For questions or issues:

1. Check this README and documentation
2. Review `config.yaml` for configuration options
3. Check logs in console output
4. For bugs, provide minimal reproducible example

## Citation

If you use this pipeline, please cite:

```
[Your systematic review paper citation once published]

Software: POD Systematic Review Pipeline
Version: 0.1.0
URL: [repository-url]
```

## License

MIT License - see LICENSE file

## Acknowledgments

- **PRISMA Guidelines:** Page MJ, et al. The PRISMA 2020 statement. BMJ 2021;372:n71
- **Newcastle-Ottawa Scale:** Wells GA, et al. The Newcastle-Ottawa Scale (NOS) for assessing the quality of nonrandomised studies in meta-analyses
- **RoB2:** Sterne JAC, et al. RoB 2: a revised tool for assessing risk of bias in randomised trials. BMJ 2019;366:l4898
- **DerSimonian-Laird:** DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials 1986;7:177-188

## Contact

For questions about this pipeline or the systematic review protocol, contact: [your-email]

---

**Last Updated:** 2026-01-07
**Version:** 0.1.0
