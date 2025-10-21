---
name: context-initialization
description: "Auto-invoked skill that creates/updates workspace claude.md to instruct Claude to always search for existing integration skills before attempting any laboratory platform or cloud service integration task."
---

# Integration Context Initialization

## Overview

This skill automatically creates or updates a `claude.md` file in the workspace root that instructs Claude to search for existing integration skills before attempting to write platform integration code. This ensures Claude uses documented authentication patterns, API access methods, and platform-specific best practices from the repository's integration skills rather than inventing integration code from scratch.

## When to Use This Skill

This skill is automatically invoked when:
- Claude Code initializes in this workspace
- User begins any laboratory platform or cloud service integration task
- User mentions platform names, LIMS, ELN, automation, or workflows
- Any task involving scientific platforms is started

**No manual invocation required** - this skill runs automatically.

## What This Skill Does

Creates or updates `claude.md` in the workspace root with instructions for Claude to:

1. **Search first**: Look for relevant integration skills in `scientific-integrations/` before writing platform code
2. **Use existing patterns**: Apply documented authentication and API access patterns
3. **Follow best practices**: Use platform-specific conventions and workflows
4. **Adapt examples**: Leverage working integration examples from `scripts/` folders

## Implementation

When invoked, this skill creates/updates the workspace `claude.md` file with a section instructing Claude to search for and use existing integration skills for any platform integration tasks.

The reference template is available in `references/claude.md`.

## Integration

Works alongside other context-initialization skills:
- `scientific-databases/context-initialization` - for database access
- `scientific-packages/context-initialization` - for Python package usage
- `scientific-thinking/context-initialization` - for analysis methodologies

Together, these ensure Claude always leverages existing expertise before attempting scientific tasks.
