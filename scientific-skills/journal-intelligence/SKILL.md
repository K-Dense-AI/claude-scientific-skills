---
name: journal-intelligence
description: "Fetch live submission requirements for any scientific journal before writing. Use before scientific-writing or any manuscript task: retrieves author guidelines, word limits, abstract format, figure limits, citation style, AI authorship policy (exact disclosure language), and LaTeX template via WebSearch/WebFetch. Stores results in journal_profile.yaml as the single source of truth for all formatting decisions. Covers Nature family, Cell Press, Elsevier, Wiley, Springer, PLOS, PNAS, ACS, RSC, IEEE, ACM, and any other journal."
---

# Journal Intelligence

Fetches live, current submission requirements for any scientific journal before writing begins. Journal guidelines, word limits, AI authorship policies, and LaTeX templates change without notice — this skill always retrieves the authoritative current version rather than relying on cached or hardcoded information.

Run this before every manuscript session. The output (`journal_profile.yaml`) becomes the single source of truth for all formatting decisions in the `scientific-writing` skill and any other writing tools.

---

## When to Use

- Before starting any manuscript, letter, review, or methods paper
- When you need to know a journal's current word limits, figure limits, or abstract format
- To get the exact AI authorship disclosure language a journal requires
- To download the current LaTeX template before writing
- To check desk-rejection criteria before investing time in a submission

**Requirements:** This skill requires **WebSearch** and **WebFetch** to be enabled in the agent's tool set. It cannot fetch live guidelines without these tools. If they are unavailable, see the fallback behavior below.

---

## What It Does

### Step 1 — Locate Author Guidelines

Searches and fetches the journal's official submission/author guidelines page.

Example URL patterns to try (publisher structures vary — these are starting points, not authoritative):
- Nature family: `https://www.nature.com/[journal]/submission-guidelines`
- Elsevier: `https://www.elsevier.com/journals/[journal]/[issn]/guide-for-authors`
- Cell Press: `https://www.cell.com/[journal]/authors`
- PLOS: `https://journals.plos.org/[journal]/s/submission-guidelines`
- Wiley, Springer, ACS, RSC, IEEE, ACM — searched via WebSearch

If the pattern-based URL fails, fall back to a direct WebSearch for `"[journal name] author guidelines submission"`.

### Step 2 — Extract and Store Requirements

Extracts all format requirements into `journal_profile.yaml`:

```yaml
journal: ""
article_type: ""
guidelines_url: ""
guidelines_fetched: ""        # date fetched

# Length limits
main_text_limit: ""           # words; what is excluded?
abstract_limit: ""            # words; structured or unstructured?
abstract_sections: ""         # required section labels if any
methods_limit: ""
figure_limit: ""
figure_legend_limit: ""
reference_limit: ""
supplementary: ""

# Format
citation_style: ""            # Vancouver numbered / Author-date / etc.
reference_format_example: ""
section_order: ""
methods_placement: ""         # inline / online / after references
special_elements: ""          # Key Resources Table, Significance Statement, eTOC, etc.

# AI policy
ai_policy_summary: ""
ai_disclosure_required: ""
ai_tools_declaration_format: ""
ai_tools_declaration_location: ""

# Template
latex_template_url: ""
latex_template_notes: ""

# Editorial
acceptance_rate: ""
desk_rejection_criteria: ""
editorial_priorities: ""
```

### Step 3 — Fetch AI Authorship Policy

Every major publisher now has an explicit AI tool use policy. Fetched separately and stored in `ai_policy.md` with the exact disclosure language required and its placement in the manuscript.

### Step 4 — Download LaTeX Template

Searches Overleaf gallery, journal page, and publisher template hub. Downloads or links to the current versioned template. Notes document class, required packages, and submission file format.

### Step 5 — Editorial Priorities

Searches for editor blog posts, "what we publish" guides, and desk-rejection criteria. Records what the journal explicitly says it will not send for review.

---

## Output

```
writing_outputs/[project]/journal/
├── journal_profile.yaml          ← source of truth for all formatting
├── journal_profile_summary.md    ← human-readable summary
├── ai_policy.md                  ← exact disclosure language + placement
├── guidelines_fetched.md         ← raw extracted text from guidelines page
└── latex_template/
    ├── [template files]
    └── template_notes.md
```

---

## Usage Examples

```
# Before writing — triggers full journal intelligence pass
Fetch the current submission requirements for Nature Biotechnology, Article format.

# Before a methods paper
What does Genome Biology currently require for a Software article?
Get the word limits, AI policy, and LaTeX template.

# Quick format check before submission
What are the current requirements for Nucleic Acids Research, Database Issue?

# Used automatically by scientific-writing
Write an Article for Nature Methods on our sequencing pipeline.
[journal-intelligence runs first, then scientific-writing uses the profile]
```

---

## Fallback When Guidelines Cannot Be Fetched

If WebSearch and WebFetch are unavailable, or if both fail (network unavailable or journal behind hard paywall):

1. Uses conservative defaults: ≤4,000 words, ≤250-word abstract, Vancouver numbered references
2. Flags every format-dependent decision with `[VERIFY AGAINST AUTHOR GUIDELINES]`
3. Reports explicitly: "Guidelines could not be fetched. Verify at [journal URL] before submission."

---

## Why Hardcoded Requirements Fail

| What changed | When |
|-------------|------|
| Nature Biotechnology article length | 2023 |
| AI authorship disclosure requirements (all major publishers) | 2023–2025 (ongoing) |
| LaTeX template versions | Continuously |
| PLOS data availability mandate | 2024 |
| Cell Press STAR Methods format | 2022 |

A requirement accurate six months ago may now cause desk rejection. Always fetch live.

---

## Recommended workflow with scientific-writing

Run `journal-intelligence` first, then pass the output to `scientific-writing`:

```
Step 1: "Fetch the current submission requirements for [journal], [article type]."
         → produces journal_profile.yaml and journal_profile_summary.md

Step 2: "Write the paper using the requirements in journal_profile.yaml."
         → scientific-writing reads the profile and applies word limits,
           citation format, required sections, and template class from it
```

The two skills are independent — `scientific-writing` does not automatically invoke `journal-intelligence`. Run them sequentially and reference the profile file explicitly in Step 2.
