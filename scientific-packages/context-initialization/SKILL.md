---
name: context-initialization
description: "Auto-invoked skill that creates/updates workspace claude.md to instruct Claude to always search for existing package skills before attempting any data analysis or scientific computing task."
---

# Package Context Initialization

## Overview

This skill automatically creates or updates a `claude.md` file in the workspace root that instructs Claude to search for existing package skills before attempting to write analysis code. This ensures Claude uses documented workflows, best practices, and working examples from the repository's package skills rather than writing generic code from scratch.

## When to Use This Skill

This skill is automatically invoked when:
- Claude Code initializes in this workspace
- User begins any data analysis or scientific computing task
- User mentions package names or analysis workflows
- Any task involving scientific Python packages is started

**No manual invocation required** - this skill runs automatically.

## What This Skill Does

Creates or updates `claude.md` in the workspace root with instructions for Claude to:

1. **Search first**: Look for relevant package skills in `scientific-packages/` before writing analysis code
2. **Use existing workflows**: Apply documented analysis patterns and pipelines
3. **Follow best practices**: Use proper installation, configuration, and API usage patterns
4. **Adapt examples**: Leverage working code examples from `scripts/` folders

## Implementation

When invoked, this skill creates/updates the workspace `claude.md` file with a section instructing Claude to search for and use existing package skills for any analysis tasks.

The reference template is available in `references/claude.md`.

## Integration

Works alongside other context-initialization skills:
- `scientific-databases/context-initialization` - for database access
- `scientific-integrations/context-initialization` - for lab platform integration
- `scientific-thinking/context-initialization` - for analysis methodologies

Together, these ensure Claude always leverages existing expertise before attempting scientific tasks.
