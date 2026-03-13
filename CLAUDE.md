# CLAUDE.md — Claude Scientific Skills

## Project Overview

Open-source collection of 142 scientific skills for AI agents (Cursor, Claude Code, Codex). Each skill is a self-contained folder under `scientific-skills/` with a `SKILL.md` file following the open Agent Skills standard. Maintained by K-Dense Inc.

## Development Workflow

### 1. Adding or Updating a Skill

Every skill lives in `scientific-skills/<skill-name>/` and MUST contain a `SKILL.md` with YAML frontmatter:

```yaml
---
name: skill-name
description: One-line description. Be specific about when to use vs alternatives.
license: MIT
metadata:
    skill-author: K-Dense Inc.
---
```

**Rules:**
- Keep `description` actionable — explain WHEN to use this skill and how it differs from similar ones
- Include `references/` directory for API docs, guides
- Include `scripts/` directory for example code when applicable
- Never add dependencies that require compilation without documenting prerequisites

### 2. Registering Skills

After adding a skill, register it in `.claude-plugin/marketplace.json` under `plugins[0].skills`:
```
"./scientific-skills/<skill-name>"
```

Then bump the version in `marketplace.json` → `metadata.version`. This triggers the release workflow.

### 3. Version Bumping & Releases

- Version is in `.claude-plugin/marketplace.json` → `metadata.version` (semver)
- Bumping the version on `main` triggers `.github/workflows/release.yml`
- The workflow creates a GitHub release with auto-generated changelog
- **Always bump version when adding/updating skills**

### 4. Commit Style

Follow existing patterns:
```
Add <skill-name> skill for <one-line purpose>
Update <skill-name> skill to v<X.Y.Z>
Bump version number to <X.Y.Z> in marketplace.json
Enhance <skill-name> documentation
```

Do NOT use conventional commits (no `feat:`, `fix:` prefixes). Keep it simple and descriptive.

### 5. Before Creating a PR

- Verify the new SKILL.md has valid YAML frontmatter (name, description, license, metadata)
- Verify the skill is registered in `marketplace.json`
- Verify no broken links in documentation
- Keep PRs focused: one skill per PR when possible

## Architecture

```
.
├── .claude-plugin/marketplace.json   # Skill registry + version (triggers releases)
├── .github/workflows/release.yml     # Auto-release on version bump
├── scientific-skills/                 # 142 skills, each with SKILL.md
│   ├── <skill-name>/
│   │   ├── SKILL.md                  # Required: frontmatter + documentation
│   │   ├── references/               # Optional: API docs, guides
│   │   └── scripts/                  # Optional: example code
├── docs/                             # Project documentation
├── README.md                         # Main documentation
└── LICENSE.md                        # MIT
```

## Common Mistakes — Do NOT Repeat These

- **Do not add skills without registering them in marketplace.json** — they won't be discoverable
- **Do not forget to bump the version** — no release will be created
- **Do not use `interface` in TypeScript** — use `type` instead (if any TS code is written)
- **Do not create overly long SKILL.md files** — keep them focused and scannable. Prefer code examples over long explanations
- **Do not duplicate skills** — check existing skills first. If rdkit already covers molecular fingerprints, don't create a separate fingerprint skill
- **Do not add Python package skills without specifying the install command** — always include `uv pip install <package>` or equivalent

## Skill Categories Quick Reference

| Category | Example Skills | Count |
|----------|---------------|-------|
| Python Packages | rdkit, scanpy, biopython, pytorch-lightning | 55+ |
| Scientific Databases | pubmed, chembl, uniprot, clinvar, cosmic | 28+ |
| Integrations | benchling, dnanexus, latchbio, omero | 15+ |
| Analysis & Communication | literature-review, scientific-writing, visualization | 30+ |
| Research & Clinical | hypothesis-generation, clinical-decision-support | 10+ |

## Key Files

- **Version/registry**: `.claude-plugin/marketplace.json`
- **Release automation**: `.github/workflows/release.yml`
- **All skills**: `scientific-skills/*/SKILL.md`

## Python Environment

- Require Python 3.9+, recommend 3.12+
- Use `uv` for package management (not pip directly)
- Install command: `uv pip install <package>`
- Virtual environments: `uv venv && source .venv/bin/activate`

## Planning & Implementation Strategy

Use plan mode for any task that touches multiple files or requires architectural decisions. The pattern:

1. **Plan first** — Start with `/plan` or ask Claude to plan before coding. A good plan should identify all files to modify, the order of operations, and edge cases.
2. **Validate the plan** — Before implementing, check: does the skill already exist? Is there overlap with existing skills? Will `marketplace.json` need updating?
3. **Implement in order** — SKILL.md first, then `references/` and `scripts/`, then register in `marketplace.json`, then bump version.
4. **If implementation drifts from plan** — Stop and re-plan instead of patching. Return to plan mode immediately.

## Quality Checklist for SKILL.md Files

Good SKILL.md files in this repo follow a consistent pattern. Use this as a template:

```markdown
---
name: <tool-name>
description: <What it does>. <When to use it vs alternatives>. <Key differentiator>.
license: <SPDX identifier>
metadata:
    skill-author: K-Dense Inc.
---

# <Tool Name>: <Subtitle>

## Overview
One paragraph: what it is, when to apply this skill.

## When to Use This Skill
Bullet list of specific scenarios.

## Quick Start
Minimal working code example (under 20 lines).

## Core Capabilities
Organized by use case, each with code examples.

## Common Patterns
Real-world patterns users will actually need.

## Troubleshooting
Known issues and their solutions.
```

**Size guidelines:**
- Target: 200-500 lines (sweet spot for usefulness without bloat)
- Under 150 lines: probably too thin — add more code examples
- Over 800 lines: too long — split into sections, move reference material to `references/`
- Skills like `latex-posters` (1600 lines) are the anti-pattern — be more concise

## Skill Description Writing Guide

The `description` field in frontmatter is the most important line. It determines when Claude auto-invokes the skill. Write it as:

```
<What it does>. <When to use — specific trigger>. <How it differs from similar skills>.
```

**Good examples from this repo:**
- rdkit: "Cheminformatics toolkit for fine-grained molecular control. SMILES/SDF parsing, descriptors... For standard workflows with simpler interface, use datamol."
- scanpy: "Standard single-cell RNA-seq analysis pipeline... For deep learning models use scvi-tools; for data format questions use anndata."

**Bad pattern:** Generic descriptions like "A tool for data analysis" — never say what makes it unique.

## Cross-Skill References

Many skills in this repo overlap. When writing a new skill, always reference alternatives:
- **rdkit** ↔ **datamol** (datamol wraps rdkit with simpler API)
- **scanpy** ↔ **scvi-tools** ↔ **anndata** (different layers of single-cell analysis)
- **biopython** ↔ **gget** ↔ **bioservices** (different levels of database access)
- **pubmed** ↔ **biorxiv** ↔ **openalex** (different literature databases)
- **chembl** ↔ **pubchem** ↔ **drugbank** (different chemical databases)

Add a "When NOT to use this skill" or "See also" section when there's overlap.

## MCP Servers Available in This Environment

This project has access to the following MCP servers. Reference them in skills when relevant:

| Server | Purpose |
|--------|---------|
| `biothings-mcp` | Gene/variant/drug annotation (MyGene, MyVariant, MyChem) |
| `gget-mcp` | Quick gene/protein lookups via gget |
| `synergy-age-mcp` | Synergistic drug combinations database |
| `opengenes-mcp` | Aging/longevity gene database |
| `pharmacology-mcp` | Drug pharmacology data |
| `biocontext-kb` | Biological context knowledge base |

## Parallel Workflow with Git Worktrees

This project uses git worktrees extensively (see `~/.claude-worktrees/`). Each worktree can run an independent Claude session. Use this for:

- Working on multiple skills in parallel (one skill per worktree)
- Testing changes without affecting the main branch
- Running independent Claude sessions that don't interfere

Create a new worktree: `git worktree add ~/.claude-worktrees/claude-scientific-skills/<name> -b <branch-name>`

## Lessons Learned — Update This Section After Corrections

<!-- Add entries here after every correction using the format below -->
<!-- Each entry captures a mistake so it's never repeated -->

- **Path encoding**: Claude Code encodes CWD paths non-trivially for `~/.claude/projects/` directories. Don't assume simple `/` → `-` replacement. Hidden directories (starting with `.`) get the dot stripped, creating double-dashes. Use fuzzy matching on the last two path components instead.
- **SKILL.md size**: Keep under 500 lines. Long skills (800+) cause context bloat when auto-invoked. Move API reference docs to `references/` subdirectory.
- **License field**: Use SPDX identifiers (`MIT`, `BSD-3-Clause`, `Apache-2.0`), not `Unknown`. Check upstream package license before writing the skill.
