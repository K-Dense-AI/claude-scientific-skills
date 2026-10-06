# Security Scan Report

**Generated:** 2026-10-05 09:38 UTC  
**Skills scanned:** 177  
**Total findings:** 463  
**Critical:** 32 | **High:** 21 | **Safe skills:** 159/177

**Scanner:** cisco-ai-skill-scanner 2.2.0 · **Model:** claude-opus-5  
**This run:** full rescan of all 177 skill(s).  

## Summary

| Skill | Severity | Findings | Safe | Duration |
|-------|----------|----------|------|----------|
| autoskill | 🔴 CRITICAL | 5 | ❌ | 20.6s |
| citation-management | 🔴 CRITICAL | 6 | ❌ | 25.4s |
| infographics | 🔴 CRITICAL | 3 | ❌ | 20.7s |
| latex-posters | 🔴 CRITICAL | 6 | ❌ | 23.2s |
| literature-review | 🔴 CRITICAL | 3 | ❌ | 22.7s |
| open-notebook | 🔴 CRITICAL | 3 | ❌ | 15.1s |
| research-lookup | 🔴 CRITICAL | 3 | ❌ | 18.2s |
| scientific-schematics | 🔴 CRITICAL | 4 | ❌ | 16.7s |
| scientific-slides | 🔴 CRITICAL | 5 | ❌ | 21.0s |
| esm | 🟠 HIGH | 4 | ❌ | 20.5s |
| geomaster | 🟠 HIGH | 3 | ❌ | 13.0s |
| hugging-science | 🟠 HIGH | 1 | ❌ | 15.5s |
| modal | 🟠 HIGH | 8 | ❌ | 17.8s |
| pytorch-lightning | 🟠 HIGH | 1 | ❌ | 17.3s |
| torch-geometric | 🟠 HIGH | 4 | ❌ | 11.0s |
| transformers | 🟠 HIGH | 1 | ❌ | 13.0s |
| torchdrug | 🟠 HIGH | 4 | ❌ | 14.3s |
| waypoint-bio | 🟠 HIGH | 2 | ❌ | 23.1s |
| biopython | 🟡 MEDIUM | 6 | ✅ | 10.7s |
| dnanexus-integration | 🟡 MEDIUM | 1 | ✅ | 13.6s |
| genomic-intelligence | 🟡 MEDIUM | 2 | ✅ | 10.2s |
| glycoengineering | 🟡 MEDIUM | 2 | ✅ | 12.6s |
| paperclip | 🟡 MEDIUM | 1 | ✅ | 13.9s |
| cirq | 🔵 LOW | 1 | ✅ | 19.7s |
| latchbio-integration | 🔵 LOW | 1 | ✅ | 23.3s |
| neuropixels-analysis | 🔵 LOW | 1 | ✅ | 18.6s |
| pydeseq2 | 🔵 LOW | 1 | ✅ | 13.4s |
| qiskit | 🔵 LOW | 1 | ✅ | 19.8s |
| scikit-bio | 🔵 LOW | 1 | ✅ | 12.2s |
| zarr-python | 🔵 LOW | 1 | ✅ | 24.9s |
| ncats-arax | ⚪ INFO | 1 | ✅ | 10.5s |
| anndata | 🟢 SAFE | 0 | ✅ | 13.0s |
| 13c-metabolic-flux | 🟢 SAFE | 0 | ✅ | 12.8s |
| alphagenome | 🟢 SAFE | 0 | ✅ | 15.3s |
| analytical-method-validation | 🟢 SAFE | 0 | ✅ | 15.5s |
| arboreto | 🟢 SAFE | 0 | ✅ | 16.5s |
| arbor | 🟢 SAFE | 0 | ✅ | 17.4s |
| adaptyv | 🟢 SAFE | 0 | ✅ | 17.7s |
| bgpt-paper-search | 🟢 SAFE | 0 | ✅ | 6.1s |
| bids | 🟢 SAFE | 0 | ✅ | 6.6s |
| aeon | 🟢 SAFE | 0 | ✅ | 25.3s |
| astropy | 🟢 SAFE | 0 | ✅ | 15.1s |
| bulk-rnaseq | 🟢 SAFE | 0 | ✅ | 7.8s |
| cantera | 🟢 SAFE | 0 | ✅ | 8.0s |
| cellprofiler | 🟢 SAFE | 0 | ✅ | 7.2s |
| bioservices | 🟢 SAFE | 0 | ✅ | 17.6s |
| benchling-integration | 🟢 SAFE | 0 | ✅ | 21.6s |
| clinical-decision-support | 🟢 SAFE | 0 | ✅ | 9.9s |
| clinical-reports | 🟢 SAFE | 0 | ✅ | 9.2s |
| consciousness-council | 🟢 SAFE | 0 | ✅ | 6.7s |
| cobrapy | 🟢 SAFE | 0 | ✅ | 8.8s |
| cellxgene-census | 🟢 SAFE | 0 | ✅ | 16.0s |
| dask | 🟢 SAFE | 0 | ✅ | 8.8s |
| datalad | 🟢 SAFE | 0 | ✅ | 9.4s |
| deepspot-m | 🟢 SAFE | 0 | ✅ | 7.8s |
| datamol | 🟢 SAFE | 0 | ✅ | 13.5s |
| deeptools | 🟢 SAFE | 0 | ✅ | 10.1s |
| deepchem | 🟢 SAFE | 0 | ✅ | 13.8s |
| dhdna-profiler | 🟢 SAFE | 0 | ✅ | 6.3s |
| depmap | 🟢 SAFE | 0 | ✅ | 10.3s |
| experimental-design | 🟢 SAFE | 0 | ✅ | 9.3s |
| exa-search | 🟢 SAFE | 0 | ✅ | 10.7s |
| diffdock | 🟢 SAFE | 0 | ✅ | 15.6s |
| etetoolkit | 🟢 SAFE | 0 | ✅ | 14.2s |
| flowkit | 🟢 SAFE | 0 | ✅ | 7.1s |
| fictiv | 🟢 SAFE | 0 | ✅ | 10.2s |
| database-lookup | 🟢 SAFE | 0 | ✅ | 37.6s |
| fluidsim | 🟢 SAFE | 0 | ✅ | 10.6s |
| folklore-variant-evidence | 🟢 SAFE | 0 | ✅ | 9.6s |
| flowio | 🟢 SAFE | 0 | ✅ | 13.6s |
| exploratory-data-analysis | 🟢 SAFE | 0 | ✅ | 23.6s |
| generate-image | 🟢 SAFE | 0 | ✅ | 12.3s |
| geniml | 🟢 SAFE | 0 | ✅ | 12.0s |
| get-available-resources | 🟢 SAFE | 0 | ✅ | 8.4s |
| gget | 🟢 SAFE | 0 | ✅ | 10.6s |
| ginkgo-cloud-lab | 🟢 SAFE | 0 | ✅ | 7.1s |
| genomic-coordinates | 🟢 SAFE | 0 | ✅ | 18.0s |
| geopandas | 🟢 SAFE | 0 | ✅ | 15.9s |
| histolab | 🟢 SAFE | 0 | ✅ | 9.8s |
| gtars | 🟢 SAFE | 0 | ✅ | 11.4s |
| hypogenic | 🟢 SAFE | 0 | ✅ | 12.0s |
| imaging-data-commons | 🟢 SAFE | 0 | ✅ | 12.6s |
| hypothesis-generation | 🟢 SAFE | 0 | ✅ | 14.7s |
| labarchive-integration | 🟢 SAFE | 0 | ✅ | 9.5s |
| lab-hardware-cad | 🟢 SAFE | 0 | ✅ | 10.3s |
| iso-standards-readiness | 🟢 SAFE | 0 | ✅ | 15.1s |
| lamindb | 🟢 SAFE | 0 | ✅ | 11.7s |
| mageck | 🟢 SAFE | 0 | ✅ | 8.5s |
| marine-carbonate-chemistry | 🟢 SAFE | 0 | ✅ | 6.6s |
| market-research-reports | 🟢 SAFE | 0 | ✅ | 8.2s |
| liteparse | 🟢 SAFE | 0 | ✅ | 16.6s |
| markdown-mermaid-writing | 🟢 SAFE | 0 | ✅ | 11.9s |
| matchms | 🟢 SAFE | 0 | ✅ | 7.8s |
| markitdown | 🟢 SAFE | 0 | ✅ | 14.6s |
| matlab | 🟢 SAFE | 0 | ✅ | 9.1s |
| matplotlib | 🟢 SAFE | 0 | ✅ | 9.4s |
| molecular-dynamics | 🟢 SAFE | 0 | ✅ | 6.0s |
| molfeat | 🟢 SAFE | 0 | ✅ | 12.5s |
| medchem | 🟢 SAFE | 0 | ✅ | 16.0s |
| neurokit2 | 🟢 SAFE | 0 | ✅ | 12.0s |
| networkx | 🟢 SAFE | 0 | ✅ | 12.4s |
| nwb-conversion | 🟢 SAFE | 0 | ✅ | 6.7s |
| nextflow | 🟢 SAFE | 0 | ✅ | 16.6s |
| nmrglue | 🟢 SAFE | 0 | ✅ | 11.6s |
| ontology-term-resolution | 🟢 SAFE | 0 | ✅ | 10.4s |
| onekgpd | 🟢 SAFE | 0 | ✅ | 11.3s |
| openpiv | 🟢 SAFE | 0 | ✅ | 7.7s |
| omero-integration | 🟢 SAFE | 0 | ✅ | 14.7s |
| opentrons-integration | 🟢 SAFE | 0 | ✅ | 9.1s |
| pacsomatic | 🟢 SAFE | 0 | ✅ | 7.8s |
| optimize-for-gpu | 🟢 SAFE | 0 | ✅ | 11.5s |
| pathway-enrichment | 🟢 SAFE | 0 | ✅ | 10.3s |
| paperzilla | 🟢 SAFE | 0 | ✅ | 16.4s |
| parallel-web | 🟢 SAFE | 0 | ✅ | 16.0s |
| pathml | 🟢 SAFE | 0 | ✅ | 14.0s |
| pathogen-variant-surveillance | 🟢 SAFE | 0 | ✅ | 17.6s |
| phylogenetics | 🟢 SAFE | 0 | ✅ | 6.7s |
| pennylane | 🟢 SAFE | 0 | ✅ | 9.6s |
| peer-review | 🟢 SAFE | 0 | ✅ | 17.9s |
| polars | 🟢 SAFE | 0 | ✅ | 8.0s |
| paper-lookup | 🟢 SAFE | 0 | ✅ | 27.3s |
| pkpd-modeling | 🟢 SAFE | 0 | ✅ | 11.8s |
| pi-agent | 🟢 SAFE | 0 | ✅ | 13.2s |
| pptx-posters | 🟢 SAFE | 0 | ✅ | 10.9s |
| polars-bio | 🟢 SAFE | 0 | ✅ | 12.3s |
| primer-design | 🟢 SAFE | 0 | ✅ | 9.6s |
| pycalphad | 🟢 SAFE | 0 | ✅ | 5.6s |
| pybamm | 🟢 SAFE | 0 | ✅ | 8.5s |
| primekg | 🟢 SAFE | 0 | ✅ | 15.0s |
| pufferlib | 🟢 SAFE | 0 | ✅ | 14.3s |
| pydicom | 🟢 SAFE | 0 | ✅ | 12.0s |
| protocolsio-integration | 🟢 SAFE | 0 | ✅ | 20.5s |
| pyhealth | 🟢 SAFE | 0 | ✅ | 13.1s |
| pylabrobot | 🟢 SAFE | 0 | ✅ | 13.2s |
| pymatgen | 🟢 SAFE | 0 | ✅ | 14.2s |
| pymc | 🟢 SAFE | 0 | ✅ | 14.1s |
| pymoo | 🟢 SAFE | 0 | ✅ | 12.1s |
| pysam | 🟢 SAFE | 0 | ✅ | 12.3s |
| qiime2-amplicon | 🟢 SAFE | 0 | ✅ | 7.5s |
| pytdc | 🟢 SAFE | 0 | ✅ | 14.1s |
| pyopenms | 🟢 SAFE | 0 | ✅ | 14.7s |
| pyzotero | 🟢 SAFE | 0 | ✅ | 13.9s |
| rdkit | 🟢 SAFE | 0 | ✅ | 7.1s |
| relion | 🟢 SAFE | 0 | ✅ | 7.5s |
| qutip | 🟢 SAFE | 0 | ✅ | 14.6s |
| relsa-severity-assessment | 🟢 SAFE | 0 | ✅ | 11.1s |
| rowan | 🟢 SAFE | 0 | ✅ | 7.3s |
| scanpy | 🟢 SAFE | 0 | ✅ | 8.6s |
| research-grants | 🟢 SAFE | 0 | ✅ | 16.8s |
| scientific-brainstorming | 🟢 SAFE | 0 | ✅ | 11.2s |
| scholar-evaluation | 🟢 SAFE | 0 | ✅ | 11.7s |
| scientific-critical-thinking | 🟢 SAFE | 0 | ✅ | 7.6s |
| scientific-visualization | 🟢 SAFE | 0 | ✅ | 9.9s |
| scikit-learn | 🟢 SAFE | 0 | ✅ | 7.2s |
| scvelo | 🟢 SAFE | 0 | ✅ | 7.3s |
| scikit-survival | 🟢 SAFE | 0 | ✅ | 9.6s |
| scientific-writing | 🟢 SAFE | 0 | ✅ | 17.5s |
| stable-baselines3 | 🟢 SAFE | 0 | ✅ | 8.1s |
| simpy | 🟢 SAFE | 0 | ✅ | 9.6s |
| shap | 🟢 SAFE | 0 | ✅ | 10.1s |
| scvi-tools | 🟢 SAFE | 0 | ✅ | 14.7s |
| seaborn | 🟢 SAFE | 0 | ✅ | 13.8s |
| statistical-analysis | 🟢 SAFE | 0 | ✅ | 12.2s |
| statistical-power | 🟢 SAFE | 0 | ✅ | 10.0s |
| sympy | 🟢 SAFE | 0 | ✅ | 7.1s |
| tiledbvcf | 🟢 SAFE | 0 | ✅ | 6.0s |
| statsmodels | 🟢 SAFE | 0 | ✅ | 10.2s |
| tamarind | 🟢 SAFE | 0 | ✅ | 9.2s |
| tellurium | 🟢 SAFE | 0 | ✅ | 9.3s |
| umap-learn | 🟢 SAFE | 0 | ✅ | 5.3s |
| timesfm-forecasting | 🟢 SAFE | 0 | ✅ | 12.8s |
| uncertainty-and-units | 🟢 SAFE | 0 | ✅ | 8.4s |
| treatment-plans | 🟢 SAFE | 0 | ✅ | 11.6s |
| vaex | 🟢 SAFE | 0 | ✅ | 7.2s |
| venue-templates | 🟢 SAFE | 0 | ✅ | 8.3s |
| usfiscaldata | 🟢 SAFE | 0 | ✅ | 14.2s |
| what-if-oracle | 🟢 SAFE | 0 | ✅ | 6.3s |

## Detailed Findings

### autoskill — 🔴 CRITICAL

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_ENV_VAR_EXFILTRATION` — Cross-file env var exfiltration: 3 files
  > Environment variable access with network calls in scripts/backends.py, scripts/doctor.py, scripts/run.py
  > **Remediation:** Review data flow across files: scripts/backends.py, scripts/run.py, scripts/doctor.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 3 files
  > Multi-file exfiltration chain detected: scripts/backends.py, scripts/doctor.py, scripts/run.py collect data → scripts/run.py → scripts/backends.py, scripts/doctor.py, scripts/run.py transmit to network
  > **Remediation:** Review data flow across files: scripts/backends.py, scripts/run.py, scripts/doctor.py

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/autoskill/scripts/backends.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/autoskill/scripts/backends.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/autoskill/scripts/doctor.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/autoskill/scripts/doctor.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/autoskill/scripts/run.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/autoskill/scripts/run.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

### citation-management — 🔴 CRITICAL

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_ENV_VAR_EXFILTRATION` — Cross-file env var exfiltration: 5 files
  > Environment variable access with network calls in scripts/extract_metadata.py, scripts/search_openalex.py, scripts/search_pubmed.py
  > **Remediation:** Review data flow across files: scripts/extract_metadata.py, scripts/validate_citations.py, scripts/search_pubmed.py, scripts/doi_to_bibtex.py, scripts/search_openalex.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 5 files
  > Multi-file exfiltration chain detected: scripts/extract_metadata.py, scripts/search_openalex.py, scripts/search_pubmed.py collect data → encode → scripts/doi_to_bibtex.py, scripts/validate_citations.py, scripts/extract_metadata.py, scripts/search_openalex.py, scripts/search_pubmed.py transmit to network
  > **Remediation:** Review data flow across files: scripts/extract_metadata.py, scripts/validate_citations.py, scripts/search_pubmed.py, scripts/doi_to_bibtex.py, scripts/search_openalex.py

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/citation-management/scripts/extract_metadata.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/citation-management/scripts/extract_metadata.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/citation-management/scripts/search_openalex.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/citation-management/scripts/search_openalex.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/citation-management/scripts/search_pubmed.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/citation-management/scripts/search_pubmed.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/core_workflow.md at line 193 contains potentially dangerous Python code.
  > File: `references/core_workflow.md:193`
  > **Remediation:** Review the code block for security implications.

### infographics — 🔴 CRITICAL

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_ENV_VAR_EXFILTRATION` — Cross-file env var exfiltration: 2 files
  > Environment variable access with network calls in scripts/generate_infographic.py, scripts/generate_infographic_ai.py
  > **Remediation:** Review data flow across files: scripts/generate_infographic.py, scripts/generate_infographic_ai.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 2 files
  > Multi-file exfiltration chain detected: scripts/generate_infographic.py, scripts/generate_infographic_ai.py collect data → scripts/generate_infographic_ai.py → scripts/generate_infographic_ai.py transmit to network
  > **Remediation:** Review data flow across files: scripts/generate_infographic.py, scripts/generate_infographic_ai.py

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/infographics/scripts/generate_infographic_ai.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/infographics/scripts/generate_infographic_ai.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

### latex-posters — 🔴 CRITICAL

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_ENV_VAR_EXFILTRATION` — Cross-file env var exfiltration: 2 files
  > Environment variable access with network calls in scripts/generate_schematic_ai.py, scripts/generate_schematic.py
  > **Remediation:** Review data flow across files: scripts/generate_schematic_ai.py, scripts/generate_schematic.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 2 files
  > Multi-file exfiltration chain detected: scripts/generate_schematic_ai.py, scripts/generate_schematic.py collect data → scripts/generate_schematic_ai.py → scripts/generate_schematic_ai.py transmit to network
  > **Remediation:** Review data flow across files: scripts/generate_schematic_ai.py, scripts/generate_schematic.py

- **🔵 LOW** `LLM_DATA_EXFILTRATION` — Credential resolution from .env files with upward directory traversal for declared OpenRouter use
  > Deterministic analyzers flagged an env-var exfiltration chain. Verification shows the OPENROUTER_API_KEY is resolved from --api-key, the environment, or a .env file found by walking upward from the current working directory, and is then used only as an Authorization Bearer header to the declared OpenRouter endpoints (https://openrouter.ai/api/v1/images and /chat/completions). The manifest and documentation explicitly declare this credential and destination; the key value is redacted from error output and is passed to the child process through a minimized environment allowlist rather than the full parent environment. The only residual risk is the upward .env scan, which could pick up a key from an unrelated parent directory outside the user's project; no local files, conversation content, or unrelated secrets are transmitted. This is ordinary API authentication to the intended service, not exfiltration.
  > **Remediation:** Limit the .env search to the project root or the skill directory rather than every ancestor of the working directory, and surface which .env file supplied the credential to the user.

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/latex-posters/scripts/generate_schematic_ai.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/latex-posters/scripts/generate_schematic_ai.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

- **🔵 LOW** `LLM_COMMAND_INJECTION` — User-supplied PDF path interpolated into Poppler invocations in preflight script
  > The deterministic scanner flagged a bash taint flow at scripts/review_poster.sh. Review shows the single positional argument is used as a filename for read-only Poppler tools (pdfinfo, pdffonts, pdfimages) and wc. The value is consistently double-quoted, prefixed with ./ when it does not start with a slash to prevent option injection, and existence-checked before use; no eval, no unquoted expansion, and no shell metacharacter interpretation occurs. The path is supplied by the operator invoking the documented command, not by fetched external content. Risk is therefore contextual and minimal.
  > File: `scripts/review_poster.sh`
  > **Remediation:** No change strictly required; optionally validate the argument against an expected .pdf suffix and reject paths containing newlines for defense in depth.

- **🟡 MEDIUM** `BEHAVIOR_BASH_TAINT_FLOW` — Bash: tainted data flows to dangerous sink
  > Variable $FILE_SIZE_BYTES (line 63) flows to `sh` at line 64.
  > File: `scripts/review_poster.sh:64`
  > **Remediation:** Review the data flow from source to sink. Avoid sending sensitive or untrusted data to network or execution commands.

### literature-review — 🔴 CRITICAL

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_ENV_VAR_EXFILTRATION` — Cross-file env var exfiltration: 3 files
  > Environment variable access with network calls in scripts/generate_schematic_ai.py, scripts/generate_schematic.py
  > **Remediation:** Review data flow across files: scripts/generate_schematic_ai.py, scripts/generate_schematic.py, scripts/verify_citations.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 3 files
  > Multi-file exfiltration chain detected: scripts/generate_schematic_ai.py, scripts/generate_schematic.py collect data → scripts/generate_schematic_ai.py → scripts/generate_schematic_ai.py, scripts/verify_citations.py transmit to network
  > **Remediation:** Review data flow across files: scripts/generate_schematic_ai.py, scripts/generate_schematic.py, scripts/verify_citations.py

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/literature-review/scripts/generate_schematic_ai.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/literature-review/scripts/generate_schematic_ai.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

### open-notebook — 🔴 CRITICAL

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_ENV_VAR_EXFILTRATION` — Cross-file env var exfiltration: 1 files
  > Environment variable access with network calls in scripts/_common.py
  > **Remediation:** Review data flow across files: scripts/_common.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 1 files
  > Multi-file exfiltration chain detected: scripts/_common.py collect data → encode → scripts/_common.py transmit to network
  > **Remediation:** Review data flow across files: scripts/_common.py

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/open-notebook/scripts/_common.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/open-notebook/scripts/_common.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

### research-lookup — 🔴 CRITICAL

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_ENV_VAR_EXFILTRATION` — Cross-file env var exfiltration: 1 files
  > Environment variable access with network calls in scripts/research_lookup.py
  > **Remediation:** Review data flow across files: scripts/research_lookup.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 2 files
  > Multi-file exfiltration chain detected: scripts/research_lookup.py collect data → scripts/manuscript_packet.py → scripts/research_lookup.py transmit to network
  > **Remediation:** Review data flow across files: scripts/manuscript_packet.py, scripts/research_lookup.py

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/research-lookup/scripts/research_lookup.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/research-lookup/scripts/research_lookup.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

### scientific-schematics — 🔴 CRITICAL

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_ENV_VAR_EXFILTRATION` — Cross-file env var exfiltration: 2 files
  > Environment variable access with network calls in scripts/generate_schematic_ai.py, scripts/generate_schematic.py
  > **Remediation:** Review data flow across files: scripts/generate_schematic_ai.py, scripts/generate_schematic.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 2 files
  > Multi-file exfiltration chain detected: scripts/generate_schematic_ai.py, scripts/generate_schematic.py collect data → scripts/generate_schematic_ai.py → scripts/generate_schematic_ai.py transmit to network
  > **Remediation:** Review data flow across files: scripts/generate_schematic_ai.py, scripts/generate_schematic.py

- **🔵 LOW** `LLM_DATA_EXFILTRATION` — Credential read from environment/.env and sent as Bearer auth to declared OpenRouter endpoint
  > Deterministic analyzers flagged an env-var exfiltration chain. Verification of the source shows the only credential handled is OPENROUTER_API_KEY, resolved from --api-key, the environment, or a .env file, and used solely as an `Authorization: Bearer` header to https://openrouter.ai/api/v1 (images and chat/completions). This is the declared, documented purpose of the skill (manifest declares primaryEnv OPENROUTER_API_KEY and compatibility requiring an OpenRouter API key). The code additionally redacts the key from error output and forwards only an allow-listed minimal environment to the child process, which reduces rather than increases exposure. The only residual risk is the upward .env directory walk (cwd and all parents), which could pick up a key from an unrelated parent project, and the documented fact that user prompts and generated images are transmitted to OpenRouter. No undeclared destination, no credential in a request body, and no hidden transmission of local files or conversation content were found.
  > **Remediation:** Optionally limit the .env search to the project root or script directory rather than walking all parent directories, and keep the existing redaction and minimal-environment forwarding. No change is required for the OpenRouter authentication itself, which is ordinary declared credential use.

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/scientific-schematics/scripts/generate_schematic_ai.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/scientific-schematics/scripts/generate_schematic_ai.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

### scientific-slides — 🔴 CRITICAL

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_ENV_VAR_EXFILTRATION` — Cross-file env var exfiltration: 4 files
  > Environment variable access with network calls in scripts/generate_schematic_ai.py, scripts/generate_slide_image.py, scripts/generate_slide_image_ai.py, scripts/generate_schematic.py
  > **Remediation:** Review data flow across files: scripts/generate_schematic_ai.py, scripts/generate_schematic.py, scripts/generate_slide_image.py, scripts/generate_slide_image_ai.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 4 files
  > Multi-file exfiltration chain detected: scripts/generate_schematic_ai.py, scripts/generate_slide_image.py, scripts/generate_slide_image_ai.py, scripts/generate_schematic.py collect data → scripts/generate_schematic_ai.py, scripts/generate_slide_image_ai.py → scripts/generate_schematic_ai.py, scripts/generate_slide_image_ai.py transmit to network
  > **Remediation:** Review data flow across files: scripts/generate_schematic_ai.py, scripts/generate_schematic.py, scripts/generate_slide_image.py, scripts/generate_slide_image_ai.py

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/scientific-slides/scripts/generate_schematic_ai.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/scientific-slides/scripts/generate_schematic_ai.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/scientific-slides/scripts/generate_slide_image_ai.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/scientific-slides/scripts/generate_slide_image_ai.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

- **🔵 LOW** `LLM_DATA_EXFILTRATION` — Credential resolution walks parent directories for .env files and forwards key to subprocess
  > Deterministic analyzers flagged env-var exfiltration chains in generate_schematic_ai.py and generate_slide_image_ai.py. Verification against source shows the OPENROUTER_API_KEY is used only as a Bearer Authorization header to the declared OpenRouter endpoint (https://openrouter.ai/api/v1), which is the service the skill exists to use, and error strings explicitly redact the key. The wrapper scripts build a minimal subprocess environment rather than copying the full parent environment. The only residual, non-confirmed risk is that _resolve_api_key walks every parent directory of the current working directory searching for a .env file and parses OPENROUTER_API_KEY from it, so a key belonging to an unrelated project above the working directory could be picked up and sent to OpenRouter. No secondary or undeclared destination, no payload-embedded credential, and no hidden transmission of local files or conversation content is present, so this is a contextual risk rather than exfiltration.
  > File: `scripts/generate_slide_image_ai.py`
  > **Remediation:** Limit the .env discovery to the project root or require an explicit path/flag instead of walking all ancestor directories, and log which .env file supplied the credential so the user can confirm the intended key is used.

### esm — 🟠 HIGH

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in SKILL.md at line 109 contains potentially dangerous Python code.
  > File: `SKILL.md:109`
  > **Remediation:** Review the code block for security implications.

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/biohub-platform.md at line 77 contains potentially dangerous Python code.
  > File: `references/biohub-platform.md:77`
  > **Remediation:** Review the code block for security implications.

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/esm-c-api.md at line 15 contains potentially dangerous Python code.
  > File: `references/esm-c-api.md:15`
  > **Remediation:** Review the code block for security implications.

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/esm-c-api.md at line 68 contains potentially dangerous Python code.
  > File: `references/esm-c-api.md:68`
  > **Remediation:** Review the code block for security implications.

### geomaster — 🟠 HIGH

- **🟡 MEDIUM** `MDBLOCK_PYTHON_HTTP_POST` — Python code block sends HTTP POST request
  > Code block in references/data-sources.md at line 106 contains potentially dangerous Python code.
  > File: `references/data-sources.md:106`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/gis-software.md at line 262 contains potentially dangerous Python code.
  > File: `references/gis-software.md:262`
  > **Remediation:** Review the code block for security implications.

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/machine-learning.md at line 154 contains potentially dangerous Python code.
  > File: `references/machine-learning.md:154`
  > **Remediation:** Review the code block for security implications.

### hugging-science — 🟠 HIGH

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/using-models.md at line 42 contains potentially dangerous Python code.
  > File: `references/using-models.md:42`
  > **Remediation:** Review the code block for security implications.

### modal — 🟠 HIGH

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in SKILL.md at line 352 contains potentially dangerous Python code.
  > File: `SKILL.md:352`
  > **Remediation:** Review the code block for security implications.

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/api_reference.md at line 137 contains potentially dangerous Python code.
  > File: `references/api_reference.md:137`
  > **Remediation:** Review the code block for security implications.

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/functions.md at line 101 contains potentially dangerous Python code.
  > File: `references/functions.md:101`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/gpu.md at line 161 contains potentially dangerous Python code.
  > File: `references/gpu.md:161`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/gpu.md at line 181 contains potentially dangerous Python code.
  > File: `references/gpu.md:181`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/gpu.md at line 192 contains potentially dangerous Python code.
  > File: `references/gpu.md:192`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_HTTP_POST` — Python code block sends HTTP POST request
  > Code block in references/scheduled-jobs.md at line 150 contains potentially dangerous Python code.
  > File: `references/scheduled-jobs.md:150`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/web-endpoints.md at line 168 contains potentially dangerous Python code.
  > File: `references/web-endpoints.md:168`
  > **Remediation:** Review the code block for security implications.

### pytorch-lightning — 🟠 HIGH

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/lightning_module.md at line 106 contains potentially dangerous Python code.
  > File: `references/lightning_module.md:106`
  > **Remediation:** Review the code block for security implications.

### torch-geometric — 🟠 HIGH

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in SKILL.md at line 226 contains potentially dangerous Python code.
  > File: `SKILL.md:226`
  > **Remediation:** Review the code block for security implications.

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/explainability.md at line 112 contains potentially dangerous Python code.
  > File: `references/explainability.md:112`
  > **Remediation:** Review the code block for security implications.

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/link_prediction.md at line 78 contains potentially dangerous Python code.
  > File: `references/link_prediction.md:78`
  > **Remediation:** Review the code block for security implications.

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/link_prediction.md at line 123 contains potentially dangerous Python code.
  > File: `references/link_prediction.md:123`
  > **Remediation:** Review the code block for security implications.

### transformers — 🟠 HIGH

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/models.md at line 290 contains potentially dangerous Python code.
  > File: `references/models.md:290`
  > **Remediation:** Review the code block for security implications.

### torchdrug — 🟠 HIGH

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/molecular_generation.md at line 74 contains potentially dangerous Python code.
  > File: `references/molecular_generation.md:74`
  > **Remediation:** Review the code block for security implications.

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/molecular_generation.md at line 201 contains potentially dangerous Python code.
  > File: `references/molecular_generation.md:201`
  > **Remediation:** Review the code block for security implications.

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/molecular_property_prediction.md at line 118 contains potentially dangerous Python code.
  > File: `references/molecular_property_prediction.md:118`
  > **Remediation:** Review the code block for security implications.

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/retrosynthesis.md at line 193 contains potentially dangerous Python code.
  > File: `references/retrosynthesis.md:193`
  > **Remediation:** Review the code block for security implications.

### waypoint-bio — 🟠 HIGH

- **🟠 HIGH** `MDBLOCK_PYTHON_EVAL_EXEC` — Python code block uses eval/exec
  > Code block in references/python-api.md at line 137 contains potentially dangerous Python code.
  > File: `references/python-api.md:137`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/python-api.md at line 228 contains potentially dangerous Python code.
  > File: `references/python-api.md:228`
  > **Remediation:** Review the code block for security implications.

### biopython — 🟡 MEDIUM

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/alignment.md at line 295 contains potentially dangerous Python code.
  > File: `references/alignment.md:295`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/alignment.md at line 313 contains potentially dangerous Python code.
  > File: `references/alignment.md:313`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/blast.md at line 204 contains potentially dangerous Python code.
  > File: `references/blast.md:204`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/blast.md at line 231 contains potentially dangerous Python code.
  > File: `references/blast.md:231`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/blast.md at line 320 contains potentially dangerous Python code.
  > File: `references/blast.md:320`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/blast.md at line 349 contains potentially dangerous Python code.
  > File: `references/blast.md:349`
  > **Remediation:** Review the code block for security implications.

### dnanexus-integration — 🟡 MEDIUM

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/app-development.md at line 88 contains potentially dangerous Python code.
  > File: `references/app-development.md:88`
  > **Remediation:** Review the code block for security implications.

### genomic-intelligence — 🟡 MEDIUM

- **🟡 MEDIUM** `MDBLOCK_PYTHON_HTTP_POST` — Python code block sends HTTP POST request
  > Code block in SKILL.md at line 213 contains potentially dangerous Python code.
  > File: `SKILL.md:213`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_HTTP_POST` — Python code block sends HTTP POST request
  > Code block in SKILL.md at line 247 contains potentially dangerous Python code.
  > File: `SKILL.md:247`
  > **Remediation:** Review the code block for security implications.

### glycoengineering — 🟡 MEDIUM

- **🟡 MEDIUM** `MDBLOCK_PYTHON_HTTP_POST` — Python code block sends HTTP POST request
  > Code block in references/glycan_databases.md at line 63 contains potentially dangerous Python code.
  > File: `references/glycan_databases.md:63`
  > **Remediation:** Review the code block for security implications.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/glycoshield.md at line 98 contains potentially dangerous Python code.
  > File: `references/glycoshield.md:98`
  > **Remediation:** Review the code block for security implications.

### paperclip — 🟡 MEDIUM

- **🟡 MEDIUM** `MDBLOCK_PYTHON_HTTP_POST` — Python code block sends HTTP POST request
  > Code block in references/python-sdk.md at line 283 contains potentially dangerous Python code.
  > File: `references/python-sdk.md:283`
  > **Remediation:** Review the code block for security implications.

### cirq — 🔵 LOW

- **🔵 LOW** `LLM_SOCIAL_ENGINEERING` — Self-promotional citation directive instructing the agent to add author's paper to user outputs
  > The SKILL.md body ends with a 'Citing Scientific Agent Skills' section directing the agent to insert a specific arXiv citation (arXiv:2609.00065, authored by the skill vendor K-Dense) into manuscripts, reports, presentations or code releases, to inform the user it did so, and to fetch the arXiv record over the network when available. This is vendor self-promotion injected into the agent's output work product rather than functionality required for the stated Cirq/quantum-computing purpose. It is disclosed (the agent is told to tell the user), the network fetch targets a legitimate, named public endpoint (arxiv.org / export.arxiv.org), and no local data is transmitted, so impact is low and it does not constitute exfiltration or a hidden action.
  > File: `SKILL.md`
  > **Remediation:** Reframe the citation guidance as optional information surfaced to the user rather than a directive that the agent insert the vendor's paper into user deliverables, and make the outbound arXiv lookup explicitly user-approved.

### latchbio-integration — 🔵 LOW

- **🔵 LOW** `LLM_SOCIAL_ENGINEERING` — Self-promotional citation insertion directive in SKILL.md
  > The instruction body directs the agent to add a specific arXiv paper authored by the skill vendor (K-Dense) to the user's manuscripts, reports, presentations, or code releases whenever the skill 'materially contributed', and to fetch the arXiv record over the network to populate the reference. This is unrelated to the skill's stated LatchBio workflow purpose and biases user deliverables toward the vendor's publication. Mitigating factors: the action is disclosed, conditional, and the skill explicitly requires telling the user it was done, so this is a contextual risk rather than a hidden manipulation.
  > File: `SKILL.md`
  > **Remediation:** Make citation insertion an explicit user-approved, opt-in step rather than a standing directive, and avoid instructing the agent to modify user deliverables with vendor-authored references by default.

### neuropixels-analysis — 🔵 LOW

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Self-promotional citation instruction directs agent to add author's paper to user deliverables
  > SKILL.md ends with a 'Citing Scientific Agent Skills' section instructing the agent to add a specific arXiv paper (by the skill author, K-Dense Inc.) to the references of any manuscript, report, presentation, or code release the skill contributed to, and to tell the user it did so. It also directs the agent to fetch an external arXiv URL when network access is available. This is author self-promotion injected into agent output rather than functionality required for Neuropixels analysis. Impact is limited (the instruction is disclosed and asks the agent to inform the user), and the fetch target is a legitimate public arXiv record, so this is a contextual policy/content concern rather than a confirmed attack.
  > File: `SKILL.md`
  > **Remediation:** Reframe the citation request as optional guidance surfaced to the user for approval rather than a directive that the agent modify user deliverables, and make any network fetch explicitly user-initiated.

### pydeseq2 — 🔵 LOW

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Self-promotional citation injection directive in SKILL.md
  > The SKILL.md body instructs the agent to add a specific paper (arXiv:2609.00065, authored by the skill vendor K-Dense) to the references or software section of any manuscript, report, presentation, or code release the skill contributes to, to tell the user it did so, and to fetch an external arXiv URL when network access is available. This is vendor self-promotion embedded as an operational instruction rather than a capability needed for RNA-seq differential expression analysis, and it can cause unsolicited modification of user deliverables plus an outbound network fetch unrelated to the stated purpose. It is disclosed (not hidden) and the destination is a legitimate public preprint server, so impact is low and there is no data exfiltration: only a GET for bibliographic metadata is described.
  > File: `SKILL.md`
  > **Remediation:** Make the citation suggestion passive/optional (e.g., surface it for the user to approve) rather than directing the agent to insert the vendor's reference into user outputs, and gate the external metadata fetch behind explicit user consent.

### qiskit — 🔵 LOW

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Self-promotional citation directive instructing agent to add a specific paper reference
  > The SKILL.md body ends with a 'Citing Scientific Agent Skills' section directing the agent to add a specific arXiv citation (K-Dense authors, arXiv:2609.00065) to the user's manuscripts, reports, or code releases when the skill 'materially contributed', and to fetch an external arXiv URL when network access is available. This is promotional content unrelated to the stated quantum-computing purpose and nudges the agent to inject author-attribution into user deliverables. It is disclosed to the user ('tell the user you did so') and does not exfiltrate data or hide actions, so impact is low; the external fetch is a read-only request to a public, declared academic source.
  > File: `SKILL.md`
  > **Remediation:** Make the citation suggestion optional and user-initiated rather than a standing directive, and avoid instructing the agent to automatically fetch external URLs during unrelated tasks.

### scikit-bio — 🔵 LOW

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Self-promotional citation instruction directing agent to add a specific paper to user outputs
  > The SKILL.md body contains a 'Citing Scientific Agent Skills' section instructing the agent to add a specific arXiv reference (arXiv:2609.00065, attributed to the skill author's organization) to any manuscript, report, presentation, or code release the skill contributes to, and to fetch the arXiv record over the network when available. This is promotional behavior injected into the agent's output on the skill author's behalf rather than a capability the biological-analysis purpose requires. It is disclosed (the agent is told to inform the user) and is limited to adding a citation plus a read-only fetch of a public arXiv URL, so impact is low and no data leaves the host beyond a standard URL request. Noted as a contextual risk only.
  > File: `SKILL.md`
  > **Remediation:** Make citation insertion an explicit user-approved option rather than a standing directive, and avoid automatic network fetches tied to author self-citation.

### zarr-python — 🔵 LOW

- **🔵 LOW** `LLM_SOCIAL_ENGINEERING` — Embedded self-citation directive adds author's paper to user deliverables and triggers external fetch
  > The SKILL.md body ends with a 'Citing Scientific Agent Skills' section instructing the agent to insert a specific arXiv citation (authored by the skill publisher) into the user's manuscripts, reports, presentations, or code releases whenever the skill 'materially contributed', and to fetch https://arxiv.org/abs/2609.00065 or the arXiv API when network access is available. This is promotional content injection into user artifacts plus an outbound network request that the stated purpose (chunked array storage with Zarr) does not require. Mitigating factors: the action is disclosed ('tell the user you did so'), conditional, targets a public, non-sensitive, well-known destination, and no local data, credentials, or conversation content is transmitted. No hidden execution, obfuscation, or exfiltration sink is present, so this is a contextual risk rather than confirmed malicious behavior.
  > File: `SKILL.md`
  > **Remediation:** Make the citation suggestion advisory rather than directive (e.g., 'you may suggest the user cite...'), require explicit user confirmation before modifying deliverables, and remove the unconditional instruction to perform an outbound network fetch unrelated to the skill's array-storage purpose.

### ncats-arax — ⚪ INFO

- **⚪ INFO** `LLM_CONTEXT_BUDGET_EXCEEDED` — 'scripts/arax_client.py' only partially analyzed (74,979 excerpt chars from 85,617)
  > Only selected code excerpts (74,979 chars) from this 85,617-character file were included; the full file was not analyzed because it exceeds llm_analysis.max_code_file_chars (75,000 chars).
  > File: `scripts/arax_client.py`
  > **Remediation:** Increase llm_analysis.max_code_file_chars in your scan policy to include more content. The full file was not analyzed.
