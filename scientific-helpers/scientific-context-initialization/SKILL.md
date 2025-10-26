---
name: context-initialization
description: "Always Auto-invoked skill that creates/updates workspace AGENTS.md to instruct the agent to always search for existing skills before attempting any scientific task."
---

# Context Initialization

## Overview

This skill automatically creates or updates an `AGENTS.md` file in the workspace root that instructs the agent to use existing skills before attempting to solve scientific tasks. This ensures the agent uses documented patterns, authentication methods, working examples, and best practices from the repository's skills rather than inventing solutions from scratch.

## When to Use This Skill

This skill is automatically invoked when:
- The agent initializes in this workspace
- User begins any scientific task (database access, package usage, platform integration, or methodology)
- User mentions specific databases, packages, platforms, or research methods
- Any scientific data retrieval, analysis, or research task is started

**No manual invocation required** - this skill runs automatically.

## What This Skill Does

Creates or updates `AGENTS.md` in the workspace root with instructions for the agent to:

1. **Search first**: Look for relevant skills across all skill categories before writing code
2. **Use existing patterns**: Apply documented API access patterns, workflows, and examples
3. **Follow best practices**: Use rate limits, authentication, configurations, and established methodologies
4. **Adapt examples**: Leverage working code examples from `scripts/` folders and reference documentation

**Important**: If `AGENTS.md` already exists in the workspace, this skill will update it intelligently rather than overwriting it. This preserves any custom instructions or modifications while ensuring the essential skill-search directives are present.

## Skill Categories

This unified context initialization covers four major skill categories:

### Database Access Tasks
- Search `scientific-databases/` for 24+ database skills
- Use documented API endpoints and authentication patterns
- Apply working code examples and best practices
- Follow rate limits and error handling patterns

### Scientific Package Usage
- Search `scientific-packages/` for 40+ Python package skills
- Use installation instructions and API usage examples
- Apply best practices and common patterns
- Leverage working scripts and reference documentation

### Laboratory Platform Integration
- Search `scientific-integrations/` for 6+ platform integration skills
- Use authentication and setup instructions
- Apply API access patterns and platform-specific best practices
- Leverage working integration examples

### Scientific Analysis & Research Methods
- Search `scientific-thinking/` for methodology skills
- Use established data analysis frameworks (EDA, statistical analysis)
- Apply research methodologies (hypothesis generation, brainstorming, critical thinking)
- Leverage communication skills (scientific writing, visualization, peer review)
- Use document processing skills (DOCX, PDF, PPTX, XLSX)

## Implementation

When invoked, this skill manages the workspace `AGENTS.md` file as follows:

- **If `AGENTS.md` does not exist**: Creates a new file using the complete template from `references/AGENTS.md`
- **If `AGENTS.md` already exists**: Updates the file to ensure the essential skill-search directives are present, while preserving any existing custom content or modifications

The file includes sections instructing the agent to search for and use existing skills across all scientific task categories.

The complete reference template is available in `references/AGENTS.md`.

## Benefits

By centralizing context initialization, this skill ensures:
- **Consistency**: The agent always uses the same approach across all skill types
- **Efficiency**: One initialization covers all scientific tasks
- **Maintainability**: Updates to the initialization strategy occur in one place
- **Completeness**: The agent is reminded to search across all available skill categories

