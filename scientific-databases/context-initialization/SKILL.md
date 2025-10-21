---
name: context-initialization
description: "Always Auto-invoked skill that creates/updates workspace claude.md to instruct Claude to always search for existing database skills before attempting any database access task."
---

# Database Context Initialization

## Overview

This skill automatically creates or updates a `claude.md` file in the workspace root that instructs Claude to search for existing database skills before attempting to write database access code. This ensures Claude uses documented API patterns, authentication methods, and best practices from the repository's database skills rather than inventing solutions from scratch.

## When to Use This Skill

This skill is automatically invoked when:
- Claude Code initializes in this workspace
- User begins any database-related task
- User mentions database access, queries, or API interactions
- Any scientific data retrieval task is started

**No manual invocation required** - this skill runs automatically.

## What This Skill Does

Creates or updates `claude.md` in the workspace root with instructions for Claude to:

1. **Search first**: Look for relevant database skills in `scientific-databases/` before writing code
2. **Use existing patterns**: Apply documented API access patterns and examples
3. **Follow best practices**: Use rate limits, authentication, and error handling from skills
4. **Adapt examples**: Leverage working code examples from `scripts/` folders

## Implementation

When invoked, this skill creates/updates the workspace `claude.md` file with a section instructing Claude to search for and use existing database skills for any database access tasks.

The reference template is available in `references/claude.md`.

## Integration

Works alongside other context-initialization skills:
- `scientific-packages/context-initialization` - for Python package usage
- `scientific-integrations/context-initialization` - for lab platform integration
- `scientific-thinking/context-initialization` - for analysis methodologies

Together, these ensure Claude always leverages existing expertise before attempting scientific tasks.
