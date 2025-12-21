---
name: mermaid
description: Guide for creating mermaid diagrams. This skill should be used when users want to create a mermaid diagram (or update an existing diagram).
allowed_tools:
  - Bash(mmdc *)
---

# Mermaid Diagram Skill

This skill helps you create, validate, and render Mermaid diagrams from natural language descriptions, code analysis, or by editing existing diagrams.

## Core Workflow

Follow these steps for every diagram creation/modification:

### 1. Understand the Request
- Identify the diagram type needed (flowchart, sequence, class, etc.)
- Extract key elements: entities, relationships, flow, hierarchy
- Determine if creating new, editing existing, or generating from code

### 2. Generate Mermaid Code
- Create syntactically correct Mermaid diagram code
- Use appropriate diagram type syntax
- Apply consistent naming and styling
- Refer to `examples/` directory for type-specific syntax (load on-demand)

### 3. Save to File
- Write code to `.mmd` file
- Generate meaningful filename based on diagram purpose (e.g., `user-authentication-flow.mmd`, `database-schema.mmd`)
- Use kebab-case for filenames

### 4. Validate with CLI
- **ALWAYS** validate using: `mmdc -i <filename>.mmd`
- Check for syntax errors and warnings
- If validation fails, proceed to step 5

### 5. Auto-Correct Errors
- Analyze error messages from `mmdc`
- Common issues:
  - Invalid syntax or keywords
  - Missing quotes around labels with spaces
  - Incorrect arrow syntax
  - Malformed node definitions
- Automatically fix the code
- Re-save and re-validate
- Repeat until validation succeeds

### 6. Render on Request
- Default format: **SVG**
- Command: `mmdc -i <filename>.mmd -o <filename>.svg`
- Alternative formats: PNG (`-o <filename>.png`), PDF (`-o <filename>.pdf`)
- Only render when user explicitly requests a preview/image

## Supported Diagram Types

Load examples from `examples/` directory as needed:

**Basic Diagrams:**
- Flowchart (`examples/flowchart.md`)
- Sequence Diagram (`examples/sequence.md`)
- Class Diagram (`examples/class.md`)
- State Diagram (`examples/state.md`)
- Entity Relationship (`examples/er.md`)

**Planning & Management:**
- Gantt Chart (`examples/gantt.md`)
- User Journey (`examples/journey.md`)
- Timeline (`examples/timeline.md`)
- Kanban (`examples/kanban.md`)

**Data Visualization:**
- Pie Chart (`examples/pie.md`)
- XY Chart (`examples/xy-chart.md`)
- Quadrant Chart (`examples/quadrant.md`)
- Sankey (`examples/sankey.md`)
- Radar (`examples/radar.md`)
- Treemap (`examples/treemap.md`)

**Technical Diagrams:**
- Git Graph (`examples/git.md`)
- C4 Diagram (`examples/c4.md`)
- Requirement Diagram (`examples/requirement.md`)
- Architecture (`examples/architecture.md`)
- Block Diagram (`examples/block.md`)
- Packet (`examples/packet.md`)

**Organizational:**
- Mindmap (`examples/mindmap.md`)
- ZenUML (`examples/zenuml.md`)

## Code Analysis → Diagram

When analyzing code to create diagrams:

**Class Diagrams:**
- Extract classes, methods, properties, inheritance, interfaces
- Show relationships: inheritance, composition, aggregation

**Sequence Diagrams:**
- Track function calls, async operations, API interactions
- Show actors, lifelines, activation boxes

**Flowcharts:**
- Map control flow, conditionals, loops
- Show function entry/exit points

**State Diagrams:**
- Identify states from enums, state machines, status fields
- Map transitions and events

## Templates

Common patterns available in `templates/common-patterns.md` (load on-demand):
- Standard flowchart structures
- API sequence patterns
- Database ER patterns
- Microservice architecture layouts
- State machine templates

## Best Practices

**Styling:**
- Use meaningful node IDs
- Add clear, concise labels
- Apply subgraphs for grouping related elements
- Use classDefs for visual consistency

**Readability:**
- Keep diagrams focused (split large diagrams into smaller ones)
- Use top-to-bottom or left-to-right orientation consistently
- Add comments in code for complex sections

**Validation:**
- NEVER skip validation step
- Always fix errors before presenting to user
- Test rendered output when in doubt

## Error Handling

If `mmdc` validation fails:
1. Read error message carefully
2. Identify line number and issue
3. Apply fix (common fixes in examples)
4. Re-validate
5. Inform user only if repeated attempts fail

## Output Format

Present to user:
```
Created: <filename>.mmd

<Show the mermaid code in a code block>

✓ Validated successfully with mermaid-cli
```

If rendered:
```
Created: <filename>.mmd
Rendered: <filename>.svg

[Show file paths]
```

## Notes

- Assume `mmdc` (mermaid-cli) is installed and available
- Default output format: SVG
- Always validate before presenting
- Fix errors autonomously
- Only load example files when needed for specific diagram type
