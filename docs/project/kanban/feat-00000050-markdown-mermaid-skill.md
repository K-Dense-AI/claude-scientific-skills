# Kanban: feat-00000050 — Markdown and Mermaid Writing Skill

_Feature board tracking PR #50 — `feat/markdown-mermaid-writing-skill` → `main`_

---

## 📋 Board overview

```mermaid
---
config:
  kanban:
    ticketBaseUrl: "https://github.com/K-Dense-AI/claude-scientific-skills/issues/"
---
kanban
    accTitle: PR 50 Markdown Mermaid Writing Skill Board
    accDescr: Kanban board tracking all tasks for the markdown-mermaid-writing skill PR from backlog through done.

    column0[Done]
        task0@{ ticket: "50", title: "Create issue-00000050 feature request doc" }
        task1@{ ticket: "50", title: "Create PR-00000050 record doc" }
        task2@{ ticket: "50", title: "Create skill directory structure" }
        task3@{ ticket: "50", title: "Port markdown style guide (733 lines)" }
        task4@{ ticket: "50", title: "Port mermaid style guide (458 lines)" }
        task5@{ ticket: "50", title: "Port all 24 diagram type guides" }
        task6@{ ticket: "50", title: "Port all 9 document templates" }
        task7@{ ticket: "50", title: "Create SKILL.md with source format philosophy" }
        task8@{ ticket: "50", title: "Add example research report" }
        task9@{ ticket: "50", title: "Fix attribution to SuperiorByteWorks-LLC/agent-project" }
        task10@{ ticket: "50", title: "Fix Mermaid \\n → <br/> in all node labels" }
        task11@{ ticket: "50", title: "Renumber docs to #50 (upstream next number)" }
        task12@{ ticket: "50", title: "Create kanban board" }

    column1[In Review]
        task13@{ ticket: "50", title: "Push branch and open PR to K-Dense-AI/claude-scientific-skills" }

    column2[Backlog]
        task14@{ ticket: "50", title: "Human review and approval" }
        task15@{ ticket: "50", title: "Merge to main" }
```

---

## 🎯 Scope

| Field | Value |
| -------------- | -------------------------------------------------------------------- |
| **PR** | [#50](../pr/pr-00000050-markdown-mermaid-skill.md) |
| **Issue** | [#50](../issues/issue-00000050-markdown-mermaid-skill.md) |
| **Branch** | `feat/markdown-mermaid-writing-skill` |
| **Target repo** | `K-Dense-AI/claude-scientific-skills` |
| **Fork** | `borealBytes/claude-scientific-skills` |
| **Started** | 2026-02-19 |
| **Status** | 🟡 In Review — PR #50 open, awaiting K-Dense team review |

---

## ✅ Completed tasks

| Task | Evidence |
| -------------------------------------------------- | -------------------------------------------------- |
| Issue-00000050 feature request doc | `docs/project/issues/issue-00000050-markdown-mermaid-skill.md` |
| PR-00000050 record doc | `docs/project/pr/pr-00000050-markdown-mermaid-skill.md` |
| Skill directory structure | `scientific-skills/markdown-mermaid-writing/` |
| Markdown style guide (733 lines) | `references/markdown_style_guide.md` |
| Mermaid style guide (458 lines) | `references/mermaid_style_guide.md` |
| 24 diagram type guides | `references/diagrams/*.md` |
| 9 document templates | `templates/*.md` |
| SKILL.md with source format philosophy | `scientific-skills/markdown-mermaid-writing/SKILL.md` |
| Example CRISPR research report | `assets/examples/example-research-report.md` |
| Attribution corrected to SuperiorByteWorks-LLC | All 42 refs updated |
| Mermaid `\n` → `<br/>` in node labels | 12 occurrences fixed across 4 files |
| Docs renumbered to #50 | Matches upstream next available number |

---

## 🔄 In progress

_Nothing currently in progress._

---

## 📋 Backlog

| Task | Blocked by | Notes |
| ----------------------- | ---------- | --------------------------------------------- |
| Human review / approval | PR push | K-Dense team reviews and merges |
| Merge to main | Approval | Standard merge process |

---

## 📊 Progress

```mermaid
pie
    accTitle: PR 50 Task Completion
    accDescr: Pie chart showing 13 tasks done out of 15 total tasks for PR 50

    "Done" : 13
    "In Review" : 1
    "Backlog" : 2
```

---

## 🔗 References

- [PR-00000050](../pr/pr-00000050-markdown-mermaid-skill.md)
- [Issue-00000050](../issues/issue-00000050-markdown-mermaid-skill.md)
- [Upstream PR](https://github.com/K-Dense-AI/claude-scientific-skills/pull/50) (once pushed)
- [Skill source](https://github.com/SuperiorByteWorks-LLC/agent-project)

---

_Last updated: 2026-02-19_
