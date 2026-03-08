# Local Environment Portability Analysis

> Generated: 2026-03-07
> Total skills scanned: 175 (scientific-skills/) + 1 (skills/team-orchestrator-test)

---

## 1. Hard-coded Local Path Dependencies

Skills with `C:\Users\Jahyun` or machine-specific absolute paths in SKILL.md or scripts.

| Skill | File | Issue | Severity |
|-------|------|-------|----------|
| **weekly-briefing** | SKILL.md | `C:/Users/Jahyun/Desktop/scan_onedrive.ps1`, OneDrive paths (hoseo, korea univ) | HIGH |
| **pptx-reviewer** | SKILL.md | `C:/Users/Jahyun/anaconda3/python.exe`, `C:/Users/Jahyun/tmp_ppt.txt`, `C:/Users/Jahyun/.claude/commands/academic-term-rules.md` | HIGH |
| **conda-env-manager** | SKILL.md | `C:\Users\Jahyun\anaconda3\envs\biosteam\python.exe` | MEDIUM |
| **journal-presentation-maker** | SKILL.md | `C:\Users\Jahyun\OneDrive - hoseo\` paths, team meeting folder paths | MEDIUM |
| **manuscript-writer** | SKILL.md | `C:\Users\Jahyun\OneDrive - hoseo\...\add_comments.py` | LOW |
| **agent-guardian** | SKILL.md | `C:\Users\Jahyun\claude-scientific-skills\scientific-skills\` | LOW |
| **primer-design** | REFERENCE.md, scripts/*.py | `C:/Users/Jahyun/results/primers`, OneDrive paths in Python source, `PycharmProjects` path | HIGH |
| **team-dashboard** | SKILL.md | `jahyunlee082@gmail.com`, Asana user GID | LOW |

### In scripts/ (Python files):

| Skill | File | Issue |
|-------|------|-------|
| **primer-design** | `subst_primer_mode.py` | `Path.home() / "OneDrive - korea univ" / "storage"` |
| **primer-design** | `udh_variant_mutagenesis.py` | `Path.home() / "OneDrive - korea univ"`, `C:/Users/.../Uronate dehydrogenase (UDH)/` |
| **primer-design** | `psxr_cofactor_mutagenesis.py` | `Path.home() / "OneDrive - korea univ"`, `PycharmProjects` path |
| **asana-extended-api** | `asana_api.py`, `mcp_server.py` | `Path.home() / ".claude" / "secrets.json"` (acceptable -- uses `Path.home()`) |

---

## 2. Local File / Config Dependencies (non-API-key)

| Skill | Dependency | Type |
|-------|-----------|------|
| **experiment-hub** | `~/.claude/commands/research-commons.md` | Config file |
| **manuscript-writer** | `~/.claude/commands/research-commons.md` | Config file |
| **research-assistant** | `~/.claude/commands/research-commons.md` | Config file |
| **weekly-briefing** | `~/.claude/projects/C--Users-Jahyun-Desktop/memory/weekly-briefing.md` | Memory file |
| **pptx-reviewer** | `~/.claude/commands/academic-term-rules.md` | Config file |
| **primer-design** | `~/.claude/skills/primer-design/` (MCP server cwd) | Install path |

---

## 3. External API Key Dependencies (EXCLUDED from portability issues)

These are expected external dependencies -- each user must provide their own keys.

| Skill | Key Required |
|-------|-------------|
| alpha-vantage | `ALPHAVANTAGE_API_KEY` |
| adaptyv | `ADAPTYV_API_KEY` |
| datacommons-client | `DC_API_KEY` |
| fda-database | `FDA_API_KEY` |
| fred-economic-data | `FRED_API_KEY` |
| generate-image | `OPENROUTER_API_KEY` |
| infographics | `OPENROUTER_API_KEY` |
| perplexity-search | `OPENROUTER_API_KEY` |
| research-lookup | `OPENROUTER_API_KEY` |
| paper-2-web | `OPENAI_API_KEY` |
| pubmed-database | `NCBI_API_KEY` (optional) |
| pymatgen | `MP_API_KEY` |
| pyzotero | `ZOTERO_API_KEY` |
| rowan | `ROWAN_API_KEY` |
| asana-extended-api | `ASANA_PAT` |
| team-dashboard | `ASANA_PAT` |

---

## 4. Portable Skills (GitHub clone only, no local config needed)

The vast majority of skills (~155/175) are **fully portable** -- they contain only:
- Documentation (SKILL.md with usage instructions)
- Reference materials
- pip-installable Python package guidance
- Public API endpoints (no auth or optional auth)

Examples of fully portable skills:
`alphafold-database`, `anndata`, `astropy`, `biopython`, `biorxiv-database`, `bioservices`, `brenda-database`, `chembl-database`, `cirq`, `clinical-reports`, `clinicaltrials-database`, `cobrapy`, `code-excellence`, `dask`, `datamol`, `deepchem`, `deeptools`, `document-skills/*`, `edgartools`, `ena-database`, `ensembl-database`, `esm`, `etetoolkit`, `exploratory-data-analysis`, `flowio`, `fluidsim`, `gene-database`, `geniml`, `geopandas`, `gget`, `gwas-database`, `histolab`, `hmdb-database`, `kegg-database`, `lamindb`, `latex-posters`, `markdown-mermaid-writing`, `markitdown`, `matchms`, `matplotlib`, `medchem`, `metabolomics-workbench-database`, `modal`, `molfeat`, `networkx`, `neurokit2`, `openalex-database`, `opentargets-database`, `pathml`, `pdb-database`, `peer-review`, `pennylane`, `plotly`, `polars`, `pptx-posters`, `pubchem-database`, `pydeseq2`, `pydicom`, `pyhealth`, `pylabrobot`, `pymc`, `pymoo`, `pyopenms`, `pysam`, `pytdc`, `pytorch-lightning`, `qiskit`, `qutip`, `rdkit`, `reactome-database`, `scanpy`, `scientific-brainstorming`, `scientific-critical-thinking`, `scientific-schematics`, `scientific-slides`, `scientific-visualization`, `scientific-writing`, `scikit-bio`, `scikit-learn`, `scikit-survival`, `scvi-tools`, `seaborn`, `shap`, `simpy`, `solid-principles`, `stable-baselines3`, `statistical-analysis`, `statsmodels`, `string-database`, `sympy`, `tiledbvcf`, `timesfm-forecasting`, `torch_geometric`, `torchdrug`, `transformers`, `umap-learn`, `uniprot-database`, `usfiscaldata`, `uspto-database`, `vaex`, `venue-templates`, `zarr-python`, `zinc-database`

---

## 5. Skills Requiring Additional Setup

### Tier A: Hard-coded user paths (must edit before use)

| Skill | What needs changing | Effort |
|-------|-------------------|--------|
| **weekly-briefing** | OneDrive paths, PowerShell script path, memory file path | Manual edit |
| **pptx-reviewer** | anaconda3 path, tmp file path, academic-term-rules.md | Manual edit |
| **journal-presentation-maker** | OneDrive paths, team meeting folder | Manual edit |
| **primer-design** | OneDrive paths in 3 Python files, output dirs, MCP cwd | Manual edit + code change |
| **conda-env-manager** | anaconda3 path example | Low (documentation only) |
| **manuscript-writer** | OneDrive path reference | Low (one line) |

### Tier B: Require local config files (not in repo)

| Skill | Config needed | Location |
|-------|--------------|----------|
| **experiment-hub** | `research-commons.md` | `~/.claude/commands/` |
| **manuscript-writer** | `research-commons.md` | `~/.claude/commands/` |
| **research-assistant** | `research-commons.md` | `~/.claude/commands/` |
| **pptx-reviewer** | `academic-term-rules.md` | `~/.claude/commands/` |

### Tier C: Require external API keys (standard setup)

See Section 3 above -- users set environment variables per their own accounts.

---

## 6. Structural Solutions

### 6.1 Immediate (automatable, no user decision needed)

1. **Create `setup/user-config.template.json`** -- Template for user-specific paths:
   ```json
   {
     "anaconda_path": "",
     "onedrive_hoseo_path": "",
     "onedrive_korea_univ_path": "",
     "lab_analyses_path": "",
     "tmp_dir": ""
   }
   ```

2. **Move `research-commons.md` and `academic-term-rules.md` into the repo** under `shared-configs/` so they are version-controlled and portable. Skills reference them via relative paths.

3. **Replace hard-coded paths with `Path.home()` patterns** in primer-design Python scripts (already partially done -- some use `Path.home()` but with institution-specific subfolder names).

### 6.2 Requires user decision

1. **OneDrive path strategy**: OneDrive folder names contain institution names ("hoseo", "korea univ"). Options:
   - a) Environment variable: `ONEDRIVE_HOSEO_PATH`, `ONEDRIVE_KU_PATH`
   - b) Config file: `~/.claude/user-paths.json`
   - c) Accept as-is (personal workflow, not meant for sharing)

2. **weekly-briefing skill**: Highly personal (specific OneDrive folders, PowerShell scripts, memory file paths). Options:
   - a) Mark as "personal/non-portable" explicitly
   - b) Refactor to read paths from config

3. **skill-rules.json**: Contains `"Jahyun"` references. Decide if this file should be user-agnostic or remain personal.

4. **team-dashboard**: Contains Asana user GID and email. Personal by nature -- mark as requiring per-user customization.

### 6.3 install.sh / setup script approach

An `install.sh` (already exists at repo root) could be extended to:
- Copy template configs to `~/.claude/commands/`
- Prompt user for OneDrive / anaconda paths
- Generate `~/.claude/user-paths.json`
- Validate API keys are set

---

## 7. Summary

| Category | Count | Action |
|----------|-------|--------|
| Fully portable (clone and use) | ~155 | None needed |
| Need API keys only | ~16 | Standard env var setup |
| Need local config files | 4 | Move configs into repo or document |
| Hard-coded user paths | 6 | Refactor or mark as personal |
| Inherently personal (not shareable) | 3 | weekly-briefing, team-dashboard, journal-presentation-maker |

**Bottom line**: 88% of skills are immediately portable via GitHub clone. The remaining 12% are split between standard API key setup (9%) and truly local-dependent skills (3%) that need path refactoring or explicit "personal" marking.
