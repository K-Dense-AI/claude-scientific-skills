# Econ AI Skills Sync

This directory contains the infrastructure for automatically syncing economics AI skills from the [awesome-econ-ai-stuff](https://github.com/meleantonio/awesome-econ-ai-stuff) repository.

## Available Skills

The following economics AI skills can be synced from the source repository:

| Category | Skill | Description | Status |
|----------|-------|-------------|--------|
| **Ideation** | research-ideation | Generate research questions from economic phenomena | Available |
| **Literature** | lit-review-assistant | Search, summarize, and synthesize papers | Available |
| **Theory** | latex-econ-model | Write economic models in LaTeX | Available |
| **Analysis** | r-econometrics | IV, DiD, RDD analysis in R using `fixest` | Available |
| **Analysis** | python-panel-data | Panel data analysis with Python | Available |
| **Analysis** | stata-regression | Regression analysis in Stata | Available |
| **Data** | api-data-fetcher | Fetch data from FRED, World Bank APIs | Available |
| **Data** | stata-data-cleaning | Clean and transform data in Stata | Available |
| **Writing** | academic-paper-writer | Draft papers with proper structure | Available |
| **Writing** | latex-tables | Generate publication-ready LaTeX tables | Available |
| **Communication** | beamer-presentation | Create Beamer slides | Available |
| **Communication** | econ-visualization | Publication-quality economic charts | Available |

## Directory Structure

```
econ-ai-skills/
├── README.md           # This file
├── config.yml          # Configuration for skill sync
├── manifest.json       # Manifest of installed skills (auto-generated)
├── scripts/
│   └── sync_econ_skills.py  # Main sync script
└── skills/             # Downloaded skills directory
    ├── analysis/
    ├── communication/
    ├── data/
    ├── ideation/
    ├── literature/
    ├── theory/
    └── writing/
```

## Usage

### Manual Sync

Run the sync script manually:

```bash
cd econ-ai-skills

# Dry run (preview changes)
python scripts/sync_econ_skills.py --dry-run

# Full sync
python scripts/sync_econ_skills.py

# Force update all skills
python scripts/sync_econ_skills.py --force

# Generate reports
python scripts/sync_econ_skills.py --output-report report.md --output-manifest manifest.json
```

### Automated Sync (GitHub Actions)

The workflow runs automatically:
- **Daily** at 6:00 AM UTC
- **On push** to `econ-ai-skills/**` or the workflow file
- **Manually** via the GitHub Actions UI (workflow_dispatch)

To trigger manually:
1. Go to Actions > "Sync Econ AI Skills"
2. Click "Run workflow"
3. Optionally enable "dry run" or "force update"

### Configuration

Edit `config.yml` to customize the sync:

```yaml
# Enable/disable specific categories
categories:
  ideation: true
  literature: true
  theory: true
  analysis: true
  data: true
  writing: true
  communication: true

# Include only specific skills (leave empty for all)
include_skills: []

# Exclude specific skills
exclude_skills: []
```

## How It Works

1. **Check Source**: The workflow checks the source repository for updates
2. **Download Skills**: The sync script downloads SKILL.md files and associated content
3. **Compare Hashes**: File hashes are compared to detect actual changes
4. **Install/Update**: New or updated skills are installed locally
5. **Create PR**: Changes are submitted as a pull request for review

## Requirements

- Python 3.9+
- `pyyaml` package
- GitHub token (for API access in CI)

## Source Repository

Skills are sourced from:
- **Repository**: [meleantonio/awesome-econ-ai-stuff](https://github.com/meleantonio/awesome-econ-ai-stuff)
- **Website**: [https://meleantonio.github.io/awesome-econ-ai-stuff/](https://meleantonio.github.io/awesome-econ-ai-stuff/)
- **License**: CC0 1.0 (Public Domain)

## Compatible AI Tools

These skills follow the open SKILL.md standard and are compatible with:
- Claude Code
- Cursor
- Gemini CLI
- GitHub Copilot (partial)
- Windsurf
- Aider (partial)
