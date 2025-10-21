---
name: context-initialization
description: "Auto-invoked skill that creates/updates workspace claude.md to instruct Claude to always search for existing methodology skills before attempting any scientific analysis, writing, visualization, or research methodology task."
---

# Scientific Thinking Context Initialization

## Overview

This skill automatically creates or updates a `claude.md` file in the workspace root that instructs Claude to search for existing methodology skills before attempting to perform scientific analysis, writing, or research tasks. This ensures Claude uses established frameworks, best practices, and structured approaches from the repository's methodology skills rather than inventing ad-hoc methods.

## When to Use This Skill

This skill is automatically invoked when:
- Claude Code initializes in this workspace
- User begins any data analysis, visualization, or scientific writing task
- User mentions analysis methods, statistics, hypothesis generation, or peer review
- Any research methodology task is started

**No manual invocation required** - this skill runs automatically.

## What This Skill Does

Creates or updates `claude.md` in the workspace root with instructions for Claude to:

1. **Search first**: Look for relevant methodology skills in `scientific-thinking/` before attempting analysis or writing tasks
2. **Use established frameworks**: Apply documented analysis workflows, writing structures, and visualization templates
3. **Follow best practices**: Use scientific standards, reporting guidelines, and community conventions
4. **Adapt examples**: Leverage templates, scripts, and style files from methodology folders

## Implementation

When invoked, this skill creates/updates the workspace `claude.md` file with a section instructing Claude to search for and use existing methodology skills for any research tasks.

The reference template is available in `references/claude.md`.

## Integration

Works alongside other context-initialization skills:
- `scientific-databases/context-initialization` - for database access
- `scientific-packages/context-initialization` - for Python package usage
- `scientific-integrations/context-initialization` - for lab platform integration

Together, these ensure Claude always leverages existing expertise before attempting scientific tasks.
