---
name: copilot-skills-integration
description: Bridge between Claude Code skills and GitHub Copilot Agent Skills. Use when porting scientific skills to VS Code/Copilot, or when importing Copilot skills into Claude Code. Covers the open Agent Skills standard shared by both platforms.
license: MIT
metadata:
    skill-author: K-Dense Inc.
---

# Copilot Skills Integration: Cross-Agent Skill Compatibility

## Overview

The Agent Skills specification is an open standard shared by Claude Code, GitHub Copilot (VS Code, CLI), and other coding agents. Skills written for one agent work across all compatible agents with minimal or no changes. This skill covers how to write, port, and discover cross-compatible agent skills.

## When to Use This Skill

- Porting a Claude Code skill to work in GitHub Copilot (VS Code)
- Importing a Copilot Agent Skill into Claude Code
- Writing a new skill that works across both platforms
- Publishing skills to the microsoft/skills or anthropics/skills repositories
- Setting up skill discovery and progressive loading

## Quick Start

### Skill Structure (works in both Claude Code and Copilot)

```
my-skill/
├── SKILL.md          # Required: YAML frontmatter + instructions
├── references/       # Optional: API docs, guides
└── scripts/          # Optional: code examples, templates
```

### Minimal SKILL.md

```markdown
---
name: my-skill
description: What it does. When to use it. How it differs from alternatives.
license: MIT
metadata:
    skill-author: Your Name
---

# My Skill

## Instructions
When the user asks about [topic], follow these steps:
1. ...
2. ...
```

## Core Concepts

### Progressive Loading (3 stages)

Both Claude Code and Copilot use progressive loading to keep context efficient:

1. **Discovery** — Agent reads only the `name` and `description` from YAML frontmatter
2. **Instructions** — When relevant, agent loads the full SKILL.md body
3. **Resources** — Agent accesses files in `references/` and `scripts/` only when referenced

This means you can install hundreds of skills without context bloat.

### YAML Frontmatter Fields

| Field | Required | Claude Code | Copilot | Notes |
|-------|----------|-------------|---------|-------|
| `name` | Yes | Yes | Yes | Used as skill identifier |
| `description` | Yes | Yes | Yes | Triggers auto-invocation |
| `license` | Yes | Yes | Optional | SPDX identifier |
| `metadata.skill-author` | Yes | Yes | Optional | Attribution |
| `allowed-tools` | No | Yes | No | Claude Code only: restrict tool access |
| `disable-model-invocation` | No | Yes | No | Claude Code only: user-invoked only |
| `context` | No | Yes | No | Claude Code only: `fork` for isolated subagent |

### Key Differences Between Platforms

| Feature | Claude Code | GitHub Copilot |
|---------|-------------|----------------|
| Install location | `~/.claude/skills/` or `.claude/skills/` | `.github/skills/` or `.vscode/skills/` |
| Slash commands | `.claude/commands/` (`.md` files) | Not yet supported for skills |
| Auto-invocation | Based on `description` match | Based on `description` match |
| Hooks | `PreToolUse`, `PostToolUse`, etc. | Not supported |
| MCP servers | Full support via `settings.json` | Via VS Code extensions |
| Subagents | `context: fork` | Not yet supported |

## Common Patterns

### Porting a Claude Code Skill to Copilot

1. Copy the skill folder to `.github/skills/` in the target repo
2. Remove Claude-specific frontmatter fields (`allowed-tools`, `context`, `disable-model-invocation`)
3. Replace Claude-specific tool references (e.g., `use the Bash tool`) with generic instructions
4. Test in VS Code with GitHub Copilot enabled

```bash
# Copy skill
cp -r ~/.claude/skills/my-skill /path/to/repo/.github/skills/

# Remove Claude-specific fields from frontmatter
# (allowed-tools, context, disable-model-invocation)
```

### Porting a Copilot Skill to Claude Code

1. Copy from `.github/skills/` to `~/.claude/skills/` or `.claude/skills/`
2. Add `license` and `metadata.skill-author` to frontmatter if missing
3. Optionally add `allowed-tools` for security
4. Works immediately — the SKILL.md format is the same

```bash
# Copy skill
cp -r /path/to/repo/.github/skills/my-skill ~/.claude/skills/
```

### Publishing to microsoft/skills Repository

The microsoft/skills repo hosts 132+ skills with 1-click install for VS Code:

```bash
# Fork and clone
gh repo fork microsoft/skills --clone

# Add your skill
mkdir -p skills/my-skill
# Copy SKILL.md and supporting files

# Create PR
gh pr create --title "Add my-skill for [purpose]"
```

### Publishing to anthropics/skills Repository

```bash
gh repo fork anthropics/skills --clone
mkdir -p skills/my-skill
# Copy SKILL.md
gh pr create --title "Add my-skill"
```

## Skill Discovery Resources

| Repository | Skills | Platform |
|-----------|--------|----------|
| [microsoft/skills](https://github.com/microsoft/skills) | 132+ | Copilot, Claude Code |
| [anthropics/skills](https://github.com/anthropics/skills) | Growing | Claude Code, Copilot |
| [K-Dense-AI/claude-scientific-skills](https://github.com/K-Dense-AI/claude-scientific-skills) | 142+ | Claude Code, Copilot |
| [awesome-copilot](https://github.com/github/awesome-copilot) | Community | Copilot |

## Microsoft 365 Copilot Cowork (Enterprise)

Microsoft's Copilot Cowork (powered by Claude) provides agentic automation across M365 apps. This is a separate product from Agent Skills — it runs as a managed service, not as local skills.

Key differences:
- **Agent Skills** = local SKILL.md files loaded by coding agents (Claude Code, Copilot)
- **Copilot Cowork** = cloud service that automates tasks across Outlook, Teams, OneDrive, etc.
- You cannot install Agent Skills into Copilot Cowork (different runtime)

## Troubleshooting

- **Skill not auto-invoked**: Check that `description` clearly states when to use it. Both platforms match on description keywords.
- **Skill loads too much context**: Move reference material to `references/` subdirectory. Only the SKILL.md body loads initially.
- **Cross-platform field errors**: Remove platform-specific frontmatter fields when porting (see table above).
- **Copilot ignores scripts/**: Copilot accesses supporting files only when referenced with relative paths in SKILL.md body, e.g., `[see template](./scripts/template.py)`.
