---
name: multi-agent-research
description: Orchestrate multi-agent scientific research workflows using Claude Code subagents. Use when a research task has distinct phases (data collection, analysis, synthesis) that can run in parallel. For single-phase research, use literature-review or hypothesis-generation instead.
license: MIT
metadata:
    skill-author: K-Dense Inc.
---

# Multi-Agent Research: Parallel Scientific Workflows

## Overview

Orchestrate complex research tasks by decomposing them into independent phases that run as parallel Claude Code subagents. Each phase specializes in one aspect (literature search, data analysis, synthesis) and results are combined into a unified output.

## When to Use This Skill

- Research tasks with 3+ distinct phases
- Tasks where phases can run independently (parallelizable)
- Complex analyses requiring different expertise per phase
- Reproducible research pipelines that need consistent structure

## When NOT to Use This Skill

- Simple literature lookups → use `literature-review`
- Single hypothesis exploration → use `hypothesis-generation`
- Data analysis only → use `exploratory-data-analysis`
- Already using denario for AI-agent research → use `denario`

## Quick Start

```python
# Conceptual multi-agent research pipeline
# Each phase maps to a Task tool invocation in Claude Code

phases = {
    "literature": "Search PubMed for papers on {topic} from last 5 years",
    "data": "Analyze dataset {dataset} for {markers}",
    "synthesis": "Combine literature and data findings into hypothesis"
}

# Phase 1 & 2 run in parallel (independent)
# Phase 3 waits for both to complete (dependent)
```

## Core Patterns

### Pattern 1: Parallel Literature + Data Collection

```
User: Research the role of BRCA1 in triple-negative breast cancer

Claude decomposes into:
├── Agent A: Search PubMed for BRCA1 + TNBC papers (last 3 years)
├── Agent B: Query ClinicalTrials.gov for BRCA1 TNBC trials
├── Agent C: Fetch BRCA1 protein interactions from STRING database
└── Synthesis: Combine all findings into structured report
```

### Pattern 2: Multi-Database Cross-Reference

```
User: Find drug candidates for target EGFR

Claude decomposes into:
├── Agent A: ChEMBL → bioactivity data for EGFR
├── Agent B: DrugBank → approved drugs targeting EGFR
├── Agent C: PubChem → compound properties
├── Agent D: ClinicalTrials → active EGFR trials
└── Synthesis: Cross-reference and rank candidates
```

### Pattern 3: Hypothesis Generation Pipeline

```
Phase 1 (parallel):
├── Literature scan → key findings
├── Database query → relevant data points
└── Existing hypothesis review → gaps identified

Phase 2 (sequential):
├── Gap analysis → where evidence is weak
└── Novel hypothesis generation → testable predictions

Phase 3 (parallel):
├── Feasibility check → can we test this?
├── Literature support → any corroborating evidence?
└── Experimental design → suggested methodology
```

## Implementation in Claude Code

When implementing multi-agent patterns, use the Task tool:

1. **Identify independent phases** -- these run in parallel
2. **Identify dependent phases** -- these wait for inputs
3. **Define clear inputs/outputs** for each phase
4. **Use the Bash subagent** for data-heavy tasks
5. **Use the general-purpose subagent** for research tasks
6. **Combine results** in the main agent thread

### Concrete Task Tool Invocation

Each subagent is launched via the Task tool with a focused prompt. Structure your prompts to return structured data:

```
Task prompt for Agent A:
"Search PubMed for papers on BRCA1 in triple-negative breast cancer
published in the last 3 years. Return a JSON list with fields:
title, authors, year, pmid, key_findings (1-2 sentences)."

Task prompt for Agent B:
"Query ClinicalTrials.gov for active BRCA1 TNBC trials.
Return a JSON list with fields: nct_id, title, phase, status,
enrollment, primary_endpoint."
```

### Defining Phase Dependencies

Map dependencies explicitly before launching agents:

```
Phase graph:
  literature_search ──┐
                      ├──> synthesis
  database_query ─────┘
                      │
  feasibility_check <─┘

Rules:
- Phases with no incoming edges run in parallel (literature_search, database_query)
- Phases with incoming edges wait for all predecessors (synthesis waits for both)
- Chain phases when output of one feeds into the next (feasibility_check after synthesis)
```

### Combining Results in the Main Thread

After all subagents return, the main agent thread merges their outputs:

```
Main agent receives:
- Agent A result: list of 15 papers with key findings
- Agent B result: list of 8 clinical trials
- Agent C result: protein interaction network with 23 nodes

Synthesis step:
1. Cross-reference paper findings with trial endpoints
2. Identify proteins mentioned in both literature and interaction network
3. Flag gaps where trials exist but literature support is weak
4. Produce structured report with sections: Evidence Summary,
   Active Trials, Interaction Network, Knowledge Gaps
```

## Available MCP Servers for Research Agents

| Server | Best For | Example Query |
|--------|----------|---------------|
| biothings-mcp | Gene/variant/drug annotation | Gene function, variant effects |
| gget-mcp | Quick gene/protein lookups | Gene ID, sequence retrieval |
| synergy-age-mcp | Drug combination data | Synergistic drug pairs |
| opengenes-mcp | Aging/longevity genes | Lifespan-associated genes |
| pharmacology-mcp | Drug pharmacology | Drug targets, mechanisms |
| biocontext-kb | Biological context | Pathway and disease context |

## Best Practices

- **Keep subagent prompts focused** -- each agent should have one clear objective. A prompt that says "search literature AND analyze data AND generate hypotheses" defeats the purpose of decomposition.
- **Request structured output** -- ask subagents to return JSON or markdown tables so the synthesis step can parse results programmatically.
- **Include a control agent** -- for large pipelines (5+ agents), dedicate one agent to quality control: checking for contradictions, missing data, or low-confidence findings.
- **Set scope boundaries** -- specify date ranges, organism filters, and database limits to prevent agents from returning overly broad results.
- **Always include a wildtype/control** -- when querying databases for variants, include a reference query so results can be compared against a baseline.

## Common Patterns

### Pattern: Systematic Review Pipeline

```
Phase 1 (parallel):
├── Agent A: PubMed search with MeSH terms → structured paper list
├── Agent B: bioRxiv/medRxiv preprint search → preprint list
└── Agent C: OpenAlex citation network → highly-cited papers

Phase 2 (sequential):
├── Deduplication → merge results, remove duplicates by DOI/PMID
├── Quality screening → filter by study design, sample size
└── Data extraction → pull key metrics from each paper

Phase 3 (parallel):
├── Agent D: Meta-analysis of extracted effect sizes
└── Agent E: Bias assessment using study quality scores

Phase 4:
└── Final report with forest plots and evidence grading
```

### Pattern: Target Validation

```
Phase 1 (parallel):
├── Agent A: Gene expression across tissues (via gget/biothings)
├── Agent B: Known disease associations (via biocontext-kb)
├── Agent C: Existing drugs for this target (via pharmacology-mcp)
└── Agent D: Genetic evidence from GWAS studies

Phase 2:
└── Score target by aggregating druggability, genetic evidence,
    expression specificity, and existing therapeutic coverage
```

## Troubleshooting

| Problem | Cause | Fix |
|---------|-------|-----|
| Agent timeout | Phase scope too broad | Break into smaller sub-tasks with tighter queries |
| Context overflow | Subagent returns too much raw data | Instruct agents to return summaries, not full records |
| Conflicting results between agents | Different data sources disagree | Add a reconciliation step in synthesis; flag conflicts explicitly |
| Missing data from one agent | Database down or no matching records | Specify fallback data sources in the prompt (e.g., "if PubMed returns nothing, try Semantic Scholar") |
| Synthesis is shallow | Agents return unstructured text | Require JSON or tabular output from each agent |
| Duplicate work across agents | Overlapping prompts | Review phase boundaries; each agent should query distinct sources or answer distinct questions |

## See Also

- `literature-review` -- Single-agent literature search and summarization
- `hypothesis-generation` -- Structured hypothesis development from evidence
- `exploratory-data-analysis` -- Data analysis without multi-agent orchestration
- `scientific-brainstorming` -- Creative ideation for research directions
- `denario` -- AI-agent research platform (alternative orchestration approach)
