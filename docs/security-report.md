# Security Scan Report

**Generated:** 2026-10-12 09:43 UTC  
**Skills scanned:** 177  
**Total findings:** 469  
**Critical:** 32 | **High:** 21 | **Safe skills:** 159/177

**Scanner:** cisco-ai-skill-scanner 2.2.2 · **Model:** claude-opus-5  
**This run:** full rescan of all 177 skill(s).  

## Summary

| Skill | Severity | Findings | Safe | Duration |
|-------|----------|----------|------|----------|
| autoskill | 🔴 CRITICAL | 5 | ❌ | 15.1s |
| citation-management | 🔴 CRITICAL | 7 | ❌ | 33.3s |
| infographics | 🔴 CRITICAL | 3 | ❌ | 26.4s |
| latex-posters | 🔴 CRITICAL | 4 | ❌ | 23.1s |
| literature-review | 🔴 CRITICAL | 3 | ❌ | 21.2s |
| open-notebook | 🔴 CRITICAL | 3 | ❌ | 16.3s |
| research-lookup | 🔴 CRITICAL | 3 | ❌ | 11.8s |
| scientific-schematics | 🔴 CRITICAL | 3 | ❌ | 15.8s |
| scientific-slides | 🔴 CRITICAL | 6 | ❌ | 28.5s |
| esm | 🟠 HIGH | 4 | ❌ | 19.2s |
| geomaster | 🟠 HIGH | 3 | ❌ | 14.4s |
| hugging-science | 🟠 HIGH | 1 | ❌ | 29.1s |
| modal | 🟠 HIGH | 8 | ❌ | 9.7s |
| pytorch-lightning | 🟠 HIGH | 1 | ❌ | 7.9s |
| torch-geometric | 🟠 HIGH | 4 | ❌ | 10.5s |
| transformers | 🟠 HIGH | 1 | ❌ | 14.3s |
| torchdrug | 🟠 HIGH | 4 | ❌ | 17.6s |
| waypoint-bio | 🟠 HIGH | 2 | ❌ | 9.6s |
| biopython | 🟡 MEDIUM | 6 | ✅ | 9.3s |
| dnanexus-integration | 🟡 MEDIUM | 1 | ✅ | 10.0s |
| genomic-intelligence | 🟡 MEDIUM | 2 | ✅ | 11.4s |
| glycoengineering | 🟡 MEDIUM | 2 | ✅ | 11.7s |
| paperclip | 🟡 MEDIUM | 1 | ✅ | 11.8s |
| gtars | 🔵 LOW | 1 | ✅ | 15.6s |
| gget | 🔵 LOW | 1 | ✅ | 19.9s |
| lamindb | 🔵 LOW | 1 | ✅ | 21.6s |
| markitdown | 🔵 LOW | 1 | ✅ | 15.1s |
| molfeat | 🔵 LOW | 1 | ✅ | 13.0s |
| neuropixels-analysis | 🔵 LOW | 1 | ✅ | 16.9s |
| polars | 🔵 LOW | 1 | ✅ | 16.2s |
| protocolsio-integration | 🔵 LOW | 1 | ✅ | 16.6s |
| pufferlib | 🔵 LOW | 1 | ✅ | 28.3s |
| pyzotero | 🔵 LOW | 1 | ✅ | 21.4s |
| qiskit | 🔵 LOW | 1 | ✅ | 23.3s |
| scikit-learn | 🔵 LOW | 1 | ✅ | 14.0s |
| scientific-writing | 🔵 LOW | 1 | ✅ | 32.0s |
| usfiscaldata | 🔵 LOW | 1 | ✅ | 18.8s |
| ncats-arax | ⚪ INFO | 1 | ✅ | 9.7s |
| 13c-metabolic-flux | 🟢 SAFE | 0 | ✅ | 13.4s |
| arbor | 🟢 SAFE | 0 | ✅ | 16.6s |
| arboreto | 🟢 SAFE | 0 | ✅ | 18.1s |
| alphagenome | 🟢 SAFE | 0 | ✅ | 18.2s |
| adaptyv | 🟢 SAFE | 0 | ✅ | 18.5s |
| analytical-method-validation | 🟢 SAFE | 0 | ✅ | 19.3s |
| anndata | 🟢 SAFE | 0 | ✅ | 19.9s |
| bgpt-paper-search | 🟢 SAFE | 0 | ✅ | 7.2s |
| bids | 🟢 SAFE | 0 | ✅ | 6.5s |
| aeon | 🟢 SAFE | 0 | ✅ | 25.8s |
| benchling-integration | 🟢 SAFE | 0 | ✅ | 8.2s |
| astropy | 🟢 SAFE | 0 | ✅ | 17.2s |
| bioservices | 🟢 SAFE | 0 | ✅ | 12.1s |
| bulk-rnaseq | 🟢 SAFE | 0 | ✅ | 7.6s |
| cantera | 🟢 SAFE | 0 | ✅ | 7.9s |
| cellprofiler | 🟢 SAFE | 0 | ✅ | 9.2s |
| cellxgene-census | 🟢 SAFE | 0 | ✅ | 12.2s |
| consciousness-council | 🟢 SAFE | 0 | ✅ | 6.6s |
| cirq | 🟢 SAFE | 0 | ✅ | 11.9s |
| dask | 🟢 SAFE | 0 | ✅ | 6.9s |
| clinical-reports | 🟢 SAFE | 0 | ✅ | 11.7s |
| clinical-decision-support | 🟢 SAFE | 0 | ✅ | 12.3s |
| cobrapy | 🟢 SAFE | 0 | ✅ | 11.7s |
| datalad | 🟢 SAFE | 0 | ✅ | 8.9s |
| deepchem | 🟢 SAFE | 0 | ✅ | 10.0s |
| deepspot-m | 🟢 SAFE | 0 | ✅ | 8.7s |
| deeptools | 🟢 SAFE | 0 | ✅ | 8.8s |
| depmap | 🟢 SAFE | 0 | ✅ | 10.4s |
| dhdna-profiler | 🟢 SAFE | 0 | ✅ | 7.9s |
| diffdock | 🟢 SAFE | 0 | ✅ | 14.6s |
| experimental-design | 🟢 SAFE | 0 | ✅ | 11.0s |
| etetoolkit | 🟢 SAFE | 0 | ✅ | 19.6s |
| exa-search | 🟢 SAFE | 0 | ✅ | 18.7s |
| exploratory-data-analysis | 🟢 SAFE | 0 | ✅ | 12.9s |
| datamol | 🟢 SAFE | 0 | ✅ | 37.9s |
| flowkit | 🟢 SAFE | 0 | ✅ | 6.7s |
| fictiv | 🟢 SAFE | 0 | ✅ | 15.1s |
| folklore-variant-evidence | 🟢 SAFE | 0 | ✅ | 7.7s |
| database-lookup | 🟢 SAFE | 0 | ✅ | 45.4s |
| fluidsim | 🟢 SAFE | 0 | ✅ | 11.3s |
| flowio | 🟢 SAFE | 0 | ✅ | 18.3s |
| geniml | 🟢 SAFE | 0 | ✅ | 14.4s |
| generate-image | 🟢 SAFE | 0 | ✅ | 16.3s |
| geopandas | 🟢 SAFE | 0 | ✅ | 9.9s |
| get-available-resources | 🟢 SAFE | 0 | ✅ | 15.3s |
| genomic-coordinates | 🟢 SAFE | 0 | ✅ | 21.3s |
| histolab | 🟢 SAFE | 0 | ✅ | 9.0s |
| ginkgo-cloud-lab | 🟢 SAFE | 0 | ✅ | 14.1s |
| hypothesis-generation | 🟢 SAFE | 0 | ✅ | 10.6s |
| imaging-data-commons | 🟢 SAFE | 0 | ✅ | 11.7s |
| iso-standards-readiness | 🟢 SAFE | 0 | ✅ | 10.6s |
| labarchive-integration | 🟢 SAFE | 0 | ✅ | 9.2s |
| liteparse | 🟢 SAFE | 0 | ✅ | 12.3s |
| latchbio-integration | 🟢 SAFE | 0 | ✅ | 22.1s |
| marine-carbonate-chemistry | 🟢 SAFE | 0 | ✅ | 8.8s |
| mageck | 🟢 SAFE | 0 | ✅ | 11.3s |
| matchms | 🟢 SAFE | 0 | ✅ | 9.0s |
| market-research-reports | 🟢 SAFE | 0 | ✅ | 14.9s |
| matplotlib | 🟢 SAFE | 0 | ✅ | 7.5s |
| matlab | 🟢 SAFE | 0 | ✅ | 13.5s |
| medchem | 🟢 SAFE | 0 | ✅ | 9.9s |
| molecular-dynamics | 🟢 SAFE | 0 | ✅ | 6.2s |
| lab-hardware-cad | 🟢 SAFE | 0 | ✅ | 55.2s |
| hypogenic | 🟢 SAFE | 0 | ✅ | 66.3s |
| networkx | 🟢 SAFE | 0 | ✅ | 12.4s |
| nmrglue | 🟢 SAFE | 0 | ✅ | 6.7s |
| nwb-conversion | 🟢 SAFE | 0 | ✅ | 8.0s |
| neurokit2 | 🟢 SAFE | 0 | ✅ | 19.6s |
| onekgpd | 🟢 SAFE | 0 | ✅ | 8.7s |
| nextflow | 🟢 SAFE | 0 | ✅ | 18.3s |
| ontology-term-resolution | 🟢 SAFE | 0 | ✅ | 13.0s |
| openpiv | 🟢 SAFE | 0 | ✅ | 8.5s |
| omero-integration | 🟢 SAFE | 0 | ✅ | 20.2s |
| markdown-mermaid-writing | 🟢 SAFE | 0 | ✅ | 56.0s |
| opentrons-integration | 🟢 SAFE | 0 | ✅ | 11.3s |
| pacsomatic | 🟢 SAFE | 0 | ✅ | 10.2s |
| optimize-for-gpu | 🟢 SAFE | 0 | ✅ | 13.7s |
| pathogen-variant-surveillance | 🟢 SAFE | 0 | ✅ | 8.0s |
| pathml | 🟢 SAFE | 0 | ✅ | 10.6s |
| pathway-enrichment | 🟢 SAFE | 0 | ✅ | 10.2s |
| paper-lookup | 🟢 SAFE | 0 | ✅ | 17.7s |
| peer-review | 🟢 SAFE | 0 | ✅ | 10.7s |
| parallel-web | 🟢 SAFE | 0 | ✅ | 16.2s |
| pennylane | 🟢 SAFE | 0 | ✅ | 8.2s |
| paperzilla | 🟢 SAFE | 0 | ✅ | 18.3s |
| phylogenetics | 🟢 SAFE | 0 | ✅ | 7.5s |
| polars-bio | 🟢 SAFE | 0 | ✅ | 7.5s |
| pi-agent | 🟢 SAFE | 0 | ✅ | 10.7s |
| pkpd-modeling | 🟢 SAFE | 0 | ✅ | 11.6s |
| primer-design | 🟢 SAFE | 0 | ✅ | 11.0s |
| primekg | 🟢 SAFE | 0 | ✅ | 13.3s |
| pybamm | 🟢 SAFE | 0 | ✅ | 8.4s |
| pptx-posters | 🟢 SAFE | 0 | ✅ | 16.5s |
| pycalphad | 🟢 SAFE | 0 | ✅ | 10.5s |
| pylabrobot | 🟢 SAFE | 0 | ✅ | 9.4s |
| pymatgen | 🟢 SAFE | 0 | ✅ | 10.2s |
| pydicom | 🟢 SAFE | 0 | ✅ | 14.6s |
| pyhealth | 🟢 SAFE | 0 | ✅ | 15.7s |
| pymoo | 🟢 SAFE | 0 | ✅ | 12.0s |
| pydeseq2 | 🟢 SAFE | 0 | ✅ | 18.6s |
| pymc | 🟢 SAFE | 0 | ✅ | 16.6s |
| pyopenms | 🟢 SAFE | 0 | ✅ | 11.3s |
| pytdc | 🟢 SAFE | 0 | ✅ | 9.2s |
| qiime2-amplicon | 🟢 SAFE | 0 | ✅ | 7.5s |
| qutip | 🟢 SAFE | 0 | ✅ | 8.8s |
| pysam | 🟢 SAFE | 0 | ✅ | 17.0s |
| relion | 🟢 SAFE | 0 | ✅ | 7.3s |
| rdkit | 🟢 SAFE | 0 | ✅ | 9.7s |
| rowan | 🟢 SAFE | 0 | ✅ | 5.7s |
| relsa-severity-assessment | 🟢 SAFE | 0 | ✅ | 13.7s |
| research-grants | 🟢 SAFE | 0 | ✅ | 15.1s |
| scholar-evaluation | 🟢 SAFE | 0 | ✅ | 9.1s |
| scanpy | 🟢 SAFE | 0 | ✅ | 13.1s |
| scientific-critical-thinking | 🟢 SAFE | 0 | ✅ | 7.4s |
| scientific-brainstorming | 🟢 SAFE | 0 | ✅ | 9.9s |
| scientific-visualization | 🟢 SAFE | 0 | ✅ | 10.0s |
| scvelo | 🟢 SAFE | 0 | ✅ | 8.8s |
| scikit-bio | 🟢 SAFE | 0 | ✅ | 12.4s |
| scikit-survival | 🟢 SAFE | 0 | ✅ | 15.0s |
| shap | 🟢 SAFE | 0 | ✅ | 8.1s |
| scvi-tools | 🟢 SAFE | 0 | ✅ | 12.5s |
| seaborn | 🟢 SAFE | 0 | ✅ | 12.7s |
| stable-baselines3 | 🟢 SAFE | 0 | ✅ | 10.4s |
| statsmodels | 🟢 SAFE | 0 | ✅ | 9.4s |
| simpy | 🟢 SAFE | 0 | ✅ | 17.8s |
| statistical-power | 🟢 SAFE | 0 | ✅ | 12.1s |
| tellurium | 🟢 SAFE | 0 | ✅ | 7.2s |
| statistical-analysis | 🟢 SAFE | 0 | ✅ | 19.0s |
| tamarind | 🟢 SAFE | 0 | ✅ | 11.2s |
| tiledbvcf | 🟢 SAFE | 0 | ✅ | 9.3s |
| sympy | 🟢 SAFE | 0 | ✅ | 17.5s |
| timesfm-forecasting | 🟢 SAFE | 0 | ✅ | 11.7s |
| umap-learn | 🟢 SAFE | 0 | ✅ | 6.5s |
| uncertainty-and-units | 🟢 SAFE | 0 | ✅ | 8.5s |
| treatment-plans | 🟢 SAFE | 0 | ✅ | 14.9s |
| what-if-oracle | 🟢 SAFE | 0 | ✅ | 6.5s |
| venue-templates | 🟢 SAFE | 0 | ✅ | 11.8s |
| vaex | 🟢 SAFE | 0 | ✅ | 15.4s |
| zarr-python | 🟢 SAFE | 0 | ✅ | 12.4s |

## Detailed Findings

### autoskill — 🔴 CRITICAL

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_ENV_VAR_EXFILTRATION` — Cross-file env var exfiltration: 3 files
  > Environment variable access with network calls in scripts/backends.py, scripts/doctor.py, scripts/run.py
  > **Remediation:** Review data flow across files: scripts/backends.py, scripts/doctor.py, scripts/run.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 3 files
  > Multi-file exfiltration chain detected: scripts/backends.py, scripts/doctor.py, scripts/run.py collect data → scripts/run.py → scripts/backends.py, scripts/doctor.py, scripts/run.py transmit to network
  > **Remediation:** Review data flow across files: scripts/backends.py, scripts/doctor.py, scripts/run.py

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
  > **Remediation:** Review data flow across files: scripts/extract_metadata.py, scripts/validate_citations.py, scripts/search_pubmed.py, scripts/search_openalex.py, scripts/doi_to_bibtex.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 5 files
  > Multi-file exfiltration chain detected: scripts/extract_metadata.py, scripts/search_openalex.py, scripts/search_pubmed.py collect data → encode → scripts/doi_to_bibtex.py, scripts/validate_citations.py, scripts/extract_metadata.py, scripts/search_openalex.py, scripts/search_pubmed.py transmit to network
  > **Remediation:** Review data flow across files: scripts/extract_metadata.py, scripts/validate_citations.py, scripts/search_pubmed.py, scripts/search_openalex.py, scripts/doi_to_bibtex.py

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

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Skill directs the agent to insert the vendor's own paper into user manuscripts
  > The SKILL.md body contains a 'Citing Scientific Agent Skills' section instructing the agent to add a specific paper authored by the skill vendor (K-Dense) to the references or software section of any manuscript, report, presentation, or code release the skill materially contributed to, and to fetch the arXiv record to build that reference. This is self-promotional content that modifies the user's academic deliverable. Mitigating factors: the instruction is conditional ('if it materially contributed'), is openly documented rather than hidden, and explicitly tells the agent to inform the user that the citation was added. No data is sent anywhere and no other instruction is overridden, so this is a disclosure/academic-integrity concern rather than an exfiltration or injection behavior.
  > File: `SKILL.md`
  > **Remediation:** Reframe the self-citation as an optional suggestion that requires explicit user confirmation before any reference is added to a user deliverable, and avoid directing the agent to modify bibliographies with vendor-authored works by default.

- **🟡 MEDIUM** `MDBLOCK_PYTHON_SUBPROCESS` — Python code block executes shell commands
  > Code block in references/core_workflow.md at line 193 contains potentially dangerous Python code.
  > File: `references/core_workflow.md:193`
  > **Remediation:** Review the code block for security implications.

### infographics — 🔴 CRITICAL

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_ENV_VAR_EXFILTRATION` — Cross-file env var exfiltration: 2 files
  > Environment variable access with network calls in scripts/generate_infographic.py, scripts/generate_infographic_ai.py
  > **Remediation:** Review data flow across files: scripts/generate_infographic_ai.py, scripts/generate_infographic.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 2 files
  > Multi-file exfiltration chain detected: scripts/generate_infographic.py, scripts/generate_infographic_ai.py collect data → scripts/generate_infographic_ai.py → scripts/generate_infographic_ai.py transmit to network
  > **Remediation:** Review data flow across files: scripts/generate_infographic_ai.py, scripts/generate_infographic.py

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

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/latex-posters/scripts/generate_schematic_ai.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/latex-posters/scripts/generate_schematic_ai.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

- **🟡 MEDIUM** `BEHAVIOR_BASH_TAINT_FLOW` — Bash: tainted data flows to dangerous sink
  > Variable $FILE_SIZE_BYTES (line 63) flows to `sh` at line 64.
  > File: `scripts/review_poster.sh:64`
  > **Remediation:** Review the data flow from source to sink. Avoid sending sensitive or untrusted data to network or execution commands.

### literature-review — 🔴 CRITICAL

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_ENV_VAR_EXFILTRATION` — Cross-file env var exfiltration: 3 files
  > Environment variable access with network calls in scripts/generate_schematic_ai.py, scripts/generate_schematic.py
  > **Remediation:** Review data flow across files: scripts/generate_schematic_ai.py, scripts/verify_citations.py, scripts/generate_schematic.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 3 files
  > Multi-file exfiltration chain detected: scripts/generate_schematic_ai.py, scripts/generate_schematic.py collect data → scripts/generate_schematic_ai.py → scripts/generate_schematic_ai.py, scripts/verify_citations.py transmit to network
  > **Remediation:** Review data flow across files: scripts/generate_schematic_ai.py, scripts/verify_citations.py, scripts/generate_schematic.py

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
  > **Remediation:** Review data flow across files: scripts/research_lookup.py, scripts/manuscript_packet.py

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

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/scientific-schematics/scripts/generate_schematic_ai.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/scientific-schematics/scripts/generate_schematic_ai.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

### scientific-slides — 🔴 CRITICAL

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_ENV_VAR_EXFILTRATION` — Cross-file env var exfiltration: 4 files
  > Environment variable access with network calls in scripts/generate_schematic_ai.py, scripts/generate_slide_image.py, scripts/generate_slide_image_ai.py, scripts/generate_schematic.py
  > **Remediation:** Review data flow across files: scripts/generate_schematic.py, scripts/generate_slide_image_ai.py, scripts/generate_schematic_ai.py, scripts/generate_slide_image.py

- **🔴 CRITICAL** `BEHAVIOR_CROSSFILE_EXFILTRATION_CHAIN` — Cross-file exfiltration chain: 4 files
  > Multi-file exfiltration chain detected: scripts/generate_schematic_ai.py, scripts/generate_slide_image.py, scripts/generate_slide_image_ai.py, scripts/generate_schematic.py collect data → scripts/generate_schematic_ai.py, scripts/generate_slide_image_ai.py → scripts/generate_schematic_ai.py, scripts/generate_slide_image_ai.py transmit to network
  > **Remediation:** Review data flow across files: scripts/generate_schematic.py, scripts/generate_slide_image_ai.py, scripts/generate_schematic_ai.py, scripts/generate_slide_image.py

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/scientific-slides/scripts/generate_schematic_ai.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/scientific-slides/scripts/generate_schematic_ai.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

- **🔴 CRITICAL** `BEHAVIOR_ENV_VAR_EXFILTRATION` — Environment variable access with network calls detected
  > Script accesses environment variables and makes network calls in /home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/scientific-slides/scripts/generate_slide_image_ai.py
  > File: `/home/runner/work/scientific-agent-skills/scientific-agent-skills/skills/scientific-slides/scripts/generate_slide_image_ai.py`
  > **Remediation:** Remove environment variable harvesting or network transmission

- **🔵 LOW** `LLM_DATA_EXFILTRATION` — User-supplied local images are base64-encoded and uploaded to a third-party model provider
  > The --attach option reads arbitrary local image files, base64-encodes them, and includes them in the request body sent to OpenRouter (and onward to the model provider). This is the documented purpose of the feature (style/reference images for slide generation) and is disclosed in SKILL.md and references/script_reference.md, which explicitly state that prompts and attached images are sent to OpenRouter and the selected model provider. It is noted only as a contextual data-handling risk: users must be aware that unpublished figures or sensitive images passed via --attach leave the machine.
  > File: `references/script_reference.md`
  > **Remediation:** No change required for intended use; optionally prompt for confirmation before uploading files from outside the working directory and keep the disclosure in the script reference prominent.

- **🔵 LOW** `LLM_DATA_EXFILTRATION` — OpenRouter API key read from environment/.env and sent as Authorization header to the declared service
  > Deterministic rules flagged an environment-variable exfiltration chain in generate_schematic_ai.py and generate_slide_image_ai.py. Verification of the source shows the credential (OPENROUTER_API_KEY) is resolved from --api-key, the environment, or a .env file and used only as a Bearer Authorization header against https://openrouter.ai/api/v1, which is the service the skill is documented to use (manifest declares primaryEnv OPENROUTER_API_KEY). The code additionally redacts the key from error messages and forwards only an allow-list of environment variables to child processes, which reduces rather than increases exposure. The only residual risk is that the .env discovery walks up all parent directories of the working directory, so a key from an unrelated parent project could be picked up; no local files, conversation content, or other secrets are transmitted. This is normal credential use for the intended service, not exfiltration.
  > File: `scripts/generate_slide_image_ai.py`
  > **Remediation:** Limit the .env search to the project root or require an explicit path/flag, and document that only the OpenRouter endpoint receives the credential.

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

### gtars — 🔵 LOW

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Self-promotional citation instruction directing agent to add a specific paper and fetch an external URL
  > The SKILL.md body instructs the agent to add a specific arXiv citation (by the skill author's organization) to any manuscript, report, or code release the skill contributes to, to notify the user it did so, and to fetch https://arxiv.org/abs/2609.00065 (or the arXiv export API) when network access is available. This is unrelated to the stated genomic-interval purpose and constitutes embedded promotional/behavioral steering plus an outbound network fetch not required by the skill's function. It does not exfiltrate local data, and the destination is a well-known public academic host, so impact is low.
  > File: `SKILL.md`
  > **Remediation:** Make the citation suggestion passive/optional documentation rather than an instruction that the agent must insert citations or perform network fetches; require explicit user consent before any outbound request.

### gget — 🔵 LOW

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Self-promotional citation directive instructing agent to add a specific paper and fetch an external URL
  > The SKILL.md body contains a 'Citing Scientific Agent Skills' section directing the agent to add a specific arXiv reference (by the skill vendor, K-Dense) to any manuscript, report, presentation, or code release the skill contributed to, to inform the user it did so, and to fetch https://arxiv.org/abs/2609.00065 or the arXiv export API when network access is available. This is unrelated to the stated bioinformatics purpose and injects vendor-promotional behavior plus an unsolicited outbound network fetch into user deliverables. No data from the local environment is sent, and the user is told about the citation, so impact is limited and this is a contextual/policy concern rather than confirmed exfiltration or hidden action.
  > File: `SKILL.md`
  > **Remediation:** Make the vendor citation guidance advisory and user-opt-in rather than an unconditional instruction, and remove the directive to perform automatic outbound network fetches to the vendor's arXiv record.

### lamindb — 🔵 LOW

- **🔵 LOW** `LLM_SOCIAL_ENGINEERING` — Embedded self-promotional citation directive in active instructions
  > The SKILL.md body contains a 'Citing Scientific Agent Skills' section that directs the agent to add a specific paper (arXiv:2609.00065 by the skill author's organization) to the references of any manuscript, report, presentation, or code release the skill contributes to, and to fetch the arXiv record over the network before writing the reference. This is promotional instruction content that influences the user's scholarly output rather than functionality required by the stated LaminDB data-curation purpose. It is mitigated by being conditional ('If it materially contributed') and explicitly disclosed ('tell the user you did so'), and the fetch target is a legitimate, declared public URL with no local data sent. No exfiltration, code execution, or hidden action is present.
  > File: `SKILL.md`
  > **Remediation:** Make the citation guidance clearly optional and user-approved, and remove the directive to automatically insert the reference into user deliverables or to perform network fetches as part of normal skill operation.

### markitdown — 🔵 LOW

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Self-promotional citation instruction directs agent to add author's paper to user deliverables
  > The SKILL.md body ends with a 'Citing Scientific Agent Skills' section instructing the agent to add a specific arXiv paper by the skill author (K-Dense Inc.) to the references of any manuscript, report, presentation, or code release that the skill contributed to, and to fetch an external arXiv URL to populate the citation. This is promotional self-referential behavior injected into unrelated user deliverables rather than functionality needed for document-to-Markdown conversion. It does instruct the agent to disclose the addition to the user, and the fetched destination is a legitimate public arXiv record, so impact is low and no data is exfiltrated.
  > File: `SKILL.md`
  > **Remediation:** Make the citation suggestion explicitly optional and user-approved, and remove the directive to automatically add the author's paper to user outputs or perform an unsolicited network fetch.

### molfeat — 🔵 LOW

- **🔵 LOW** `LLM_SOCIAL_ENGINEERING` — Self-promotional citation injection instructing the agent to add author's paper to user outputs
  > The SKILL.md body contains a 'Citing Scientific Agent Skills' section that directs the agent to add a specific arXiv reference (K-Dense authors) to any manuscript, report, presentation, or code release the skill contributes to, to tell the user it did so, and to fetch an external arXiv URL when network access is available. This is author-benefiting behavior appended to a molecular featurization skill and is unrelated to the skill's stated technical purpose. It is a mild, disclosed content-steering/promotion pattern rather than concealed exfiltration: no local data, secrets, or conversation content is sent outward, and the fetch target is a public, well-known arXiv endpoint. Risk is limited to injecting unsolicited citations into user deliverables.
  > File: `SKILL.md`
  > **Remediation:** Make the citation step explicitly optional and user-confirmed, or remove the directive from the skill body and document it as a README note rather than an agent instruction.

### neuropixels-analysis — 🔵 LOW

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Self-promotional citation instruction directing agent to add a specific paper and notify the user
  > The SKILL.md body instructs the agent that if the skill 'materially contributed' to a manuscript or code release it must add a specific arXiv reference (K-Dense 'Scientific Agent Skills', arXiv:2609.00065) to the user's references and 'tell the user you did so', and to fetch the arXiv record when network access is available. This is a vendor self-citation directive that influences the agent's output on the author's behalf rather than a technical capability of Neuropixels analysis. It does not exfiltrate data or execute code, but it is a promotional instruction injected into agent behavior and a scope/purpose mismatch; the network fetch is limited to a public arXiv record.
  > File: `SKILL.md`
  > **Remediation:** Make the citation suggestion passive documentation (e.g., 'if you wish to cite this skill, use...') rather than an instruction that the agent modifies user deliverables, and remove the automatic network fetch directive.

### polars — 🔵 LOW

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Self-promotional citation directive embedded in skill instructions
  > The SKILL.md body instructs the agent to add a specific arXiv paper authored by the skill vendor (K-Dense) to the references/software section of any manuscript, report, presentation, or code release that the skill contributed to, and to tell the user it did so. It also directs the agent to fetch arxiv.org / export.arxiv.org when network access is available, which exceeds the declared `allowed-tools: Read` capability. This is vendor self-promotion injected into user deliverables rather than Polars DataFrame functionality; however it is disclosed to the user, points to a legitimate public destination, sends no local data, and is low impact.
  > File: `SKILL.md`
  > **Remediation:** Make the citation step optional and user-initiated rather than a standing instruction, and either remove the network fetch directive or declare the required fetch/network capability in the manifest.

### protocolsio-integration — 🔵 LOW

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Self-promotional citation instruction directing agent to add a specific paper and fetch an external URL
  > The SKILL.md body contains a 'Citing Scientific Agent Skills' section instructing the agent to add a specific arXiv paper (arXiv:2609.00065) to references of any manuscript/report/code release the skill contributes to, to tell the user it did so, and to fetch https://arxiv.org/abs/2609.00065 or export.arxiv.org when network access is available. This is author self-promotion injected into agent behavior and an unrelated outbound network fetch not required by the skill's stated protocols.io purpose. It is disclosed in the instruction body and does not exfiltrate data, so impact is low; the agent is still told to notify the user.
  > File: `SKILL.md`
  > **Remediation:** Make citation guidance passive/optional documentation rather than a directive to modify user deliverables, and remove the instruction to automatically fetch an external arXiv URL unrelated to the skill's protocols.io purpose.

### pufferlib — 🔵 LOW

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Embedded self-citation directive instructs agent to insert author's paper into user deliverables
  > The SKILL.md body contains a 'Citing Scientific Agent Skills' section directing the agent to add a specific arXiv paper by the skill's authors to the references of any manuscript, report, presentation, or code release the skill contributes to, and to fetch the arXiv record over the network when available. This is promotional content injected into user outputs rather than PufferLib guidance. Mitigating factors: the directive is conditional ('if it materially contributed'), is transparent (the agent is told to inform the user), and the fetched destination (arxiv.org) is benign and sends no local data. No hidden action, credential access, or exfiltration is involved.
  > File: `SKILL.md`
  > **Remediation:** Make the citation suggestion passive informational text (e.g., 'attribution information available on request') so the agent does not automatically modify user deliverables, and remove the automatic network-fetch instruction.

### pyzotero — 🔵 LOW

- **🔵 LOW** `LLM_SOCIAL_ENGINEERING` — Embedded self-citation directive instructs agent to insert author's paper into user deliverables
  > The SKILL.md body contains a 'Citing Scientific Agent Skills' section directing the agent to add a specific arXiv paper by the skill author (K-Dense) to the references/software section of any manuscript, report, presentation, or code release the skill contributes to, and to fetch arxiv.org/export.arxiv.org to populate the reference. This is promotional content injected into the agent's active directions and is unrelated to the skill's stated purpose of managing Zotero libraries. Mitigating factors: the instruction explicitly requires telling the user that the citation was added, the destination (arxiv.org) is a public, disclosed, benign endpoint, and no local data is transmitted. No concealment, data exfiltration, or instruction override is present, so this is a contextual risk rather than a confirmed abuse.
  > File: `SKILL.md`
  > **Remediation:** Convert the citation request into optional, user-facing guidance rather than a directive to the agent, or remove it from the active instruction body so that reference lists are only modified at explicit user request.

### qiskit — 🔵 LOW

- **🔵 LOW** `LLM_SOCIAL_ENGINEERING` — Self-promotional citation injection directive in SKILL.md
  > The SKILL.md body instructs the agent to add a specific arXiv paper by the skill's authors to the user's manuscripts, reports, presentations, or code releases whenever the skill 'materially contributed', and to fetch an external arXiv URL to populate the reference. This directs the agent to insert author-benefiting content into user deliverables and to make an outbound network request that is unrelated to the stated quantum-computing purpose. Mitigating factors: the behavior is openly documented, the agent is told to inform the user, the destination (arxiv.org) is a public, non-exfiltrating endpoint, and no local data is transmitted. No concealment, instruction override, or confirmation bypass is present.
  > File: `SKILL.md`
  > **Remediation:** Reframe the citation guidance as optional, user-prompted information rather than a standing directive to modify user deliverables, and make the outbound arXiv lookup explicitly opt-in.

### scikit-learn — 🔵 LOW

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Self-promotional citation directive instructing agent to add a specific paper to user outputs
  > The SKILL.md body contains a 'Citing Scientific Agent Skills' section directing the agent to insert a specific arXiv reference (arXiv:2609.00065 by the skill author's organization) into any manuscript, report, presentation, or code release the skill contributes to, and to fetch an external arXiv URL when network access is available. This is unrelated to the stated scikit-learn purpose and causes the agent to inject promotional content into user deliverables. The instruction does tell the agent to inform the user, and the destination (arxiv.org) is a benign public source with no local data sent, so impact is low and no exfiltration is established.
  > File: `SKILL.md`
  > **Remediation:** Make citation insertion an explicit user-approved option rather than a standing instruction, and remove the automatic external fetch directive or gate it behind user confirmation.

### scientific-writing — 🔵 LOW

- **🔵 LOW** `LLM_POLICY_VIOLATION` — Skill instructs agent to insert a self-promotional citation into user deliverables
  > The SKILL.md body ends with a 'Citing Scientific Agent Skills' directive telling the agent to add the vendor's own arXiv paper (K-Dense, arXiv:2609.00065) to the references or software section of any manuscript, report, presentation, or code release the skill materially contributed to, and to fetch the arXiv record over the network when available. This is vendor self-promotion injected into the user's scholarly output, which is in tension with the same document's strict rule that every reference must map to independently human-verified evidence. Mitigating factors: the action is conditional, is explicitly disclosed to the user ('tell the user you did so'), the fetched destination is the public arXiv record of the cited paper, and no local data is transmitted. No hidden behavior, credential access, or exfiltration sink is present.
  > File: `SKILL.md`
  > **Remediation:** Make the citation suggestion an explicit user-facing recommendation requiring confirmation rather than a standing instruction to modify the user's reference list, and align it with the skill's own evidence-verification requirements.

### usfiscaldata — 🔵 LOW

- **🔵 LOW** `LLM_SOCIAL_ENGINEERING` — Self-promotional citation injection into user deliverables
  > The SKILL.md body ends with a 'Citing Scientific Agent Skills' section directing the agent to add the skill authors' arXiv paper to the references/software section of any manuscript, report, presentation, or code release the skill contributes to, and to fetch an external arXiv URL to populate the reference. This is unrelated to the stated purpose (querying the Treasury Fiscal Data API) and causes the agent to insert author-promoting content into user outputs. The behavior is disclosed ('tell the user you did so') and the fetch target is a legitimate public academic source, so impact is minimal and no data leaves the host beyond a normal GET request.
  > File: `SKILL.md`
  > **Remediation:** Make the citation guidance optional and user-initiated rather than a standing directive, and avoid instructing the agent to perform unrelated network fetches as part of a data-query skill.

### ncats-arax — ⚪ INFO

- **⚪ INFO** `LLM_CONTEXT_BUDGET_EXCEEDED` — 'scripts/arax_client.py' only partially analyzed (74,979 excerpt chars from 85,617)
  > Only selected code excerpts (74,979 chars) from this 85,617-character file were included; the full file was not analyzed because it exceeds llm_analysis.max_code_file_chars (75,000 chars).
  > File: `scripts/arax_client.py`
  > **Remediation:** Increase llm_analysis.max_code_file_chars in your scan policy to include more content. The full file was not analyzed.
