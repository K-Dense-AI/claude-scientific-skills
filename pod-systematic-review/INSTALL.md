# Installation & Quick Start Guide

## System Requirements

- Python 3.11 or higher
- pip (Python package installer)
- 2GB RAM minimum
- Internet connection for PubMed API access

## Installation Steps

### 1. Install Python Dependencies

```bash
# From project root directory
pip install -e ".[dev]"
```

This installs:
- Core dependencies (biopython, requests, pandas, numpy, scipy, matplotlib, streamlit, etc.)
- Development tools (pytest, mypy, black, ruff)

### 2. Configure API Credentials

```bash
# Copy example environment file
cp .env.example .env

# Edit .env with your credentials
nano .env  # or use your preferred editor
```

Add your NCBI credentials:
```
NCBI_EMAIL=your_email@example.com
NCBI_API_KEY=your_ncbi_api_key_here
```

**Get NCBI API Key:**
1. Go to https://www.ncbi.nlm.nih.gov/account/
2. Register or sign in
3. Navigate to Settings → API Key Management
4. Create new API key
5. Copy key to .env file

**Note:** API key is optional but **strongly recommended** (increases rate limit 3x).

### 3. Verify Installation

```bash
# Check CLI is installed
pod-review --help

# Run tests
pytest tests/ -v
```

## First Run Example

### Retrieve Sample Data from PubMed

```bash
# Retrieve up to 100 records (for testing)
pod-review retrieve --max-results 100 --output-dir data/raw
```

### Import Example Data

```bash
# Import the provided example RIS file
pod-review import-records examples/example_pubmed.ris --output-dir data/imported
```

### Deduplicate

```bash
pod-review dedup --input-dir data --output-dir data/processed
```

### Load for Screening

```bash
pod-review load-screening --input-file data/processed/deduplicated.jsonl
```

### Launch Screening App

```bash
streamlit run src/pod_review/screening/app.py
```

Open browser to http://localhost:8501

## Common Issues

### Issue: "ModuleNotFoundError: No module named 'pod_review'"

**Solution:** Install the package in development mode:
```bash
pip install -e .
```

### Issue: "BioPython import error"

**Solution:** Install BioPython:
```bash
pip install biopython
```

### Issue: PubMed retrieval slow

**Solution:** Add NCBI API key to .env file (see step 2 above)

### Issue: Streamlit not starting

**Solution:**
```bash
pip install streamlit
streamlit run src/pod_review/screening/app.py
```

### Issue: Permission denied when creating directories

**Solution:** Create directories manually:
```bash
mkdir -p data/{raw,processed,imported}
mkdir -p outputs/{prisma,extraction,meta,rob}
```

## Alternative: Using Makefile

If `make` is available on your system:

```bash
# Install
make install

# Create .env
make .env  # Then edit .env manually

# Retrieve
make retrieve

# Deduplicate
make dedup

# Screen
make screen

# Run tests
make test
```

## Next Steps

1. Review `config.yaml` and customize for your review
2. Run full PubMed retrieval: `make retrieve`
3. Import records from other databases (Embase, CINAHL, etc.)
4. Complete screening workflow
5. Perform data extraction
6. Run meta-analyses

See README.md for detailed usage instructions.
