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

1. **Identify independent phases** — these run in parallel
2. **Identify dependent phases** — these wait for inputs
3. **Define clear inputs/outputs** for each phase
4. **Use the Bash subagent** for data-heavy tasks
5. **Use the general-purpose subagent** for research tasks
6. **Combine results** in the main agent thread

## Available MCP Servers for Research Agents

| Server | Best For |
|--------|----------|
| biothings-mcp | Gene/variant/drug annotation |
| gget-mcp | Quick gene/protein lookups |
| synergy-age-mcp | Drug combination data |
| opengenes-mcp | Aging/longevity genes |
| pharmacology-mcp | Drug pharmacology |
| biocontext-kb | Biological context |

## Troubleshooting

- **Agent timeout**: Break large phases into smaller sub-tasks
- **Context overflow**: Use Bash subagents for data processing, return only summaries
- **Conflicting results**: Include a reconciliation step in synthesis phase
- **Missing data**: Always include fallback data sources per phase
