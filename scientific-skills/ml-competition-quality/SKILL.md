---
name: ml-competition-quality
description: "Review and enforce code quality in tabular ML competition Python code. Use when: reviewing any src/*.py or scripts/*.py for dead code, unclear contracts, single-responsibility violations, implicit types, or silent data mutations; investigating production bugs via the 16-pattern pitfall catalogue; running the 6-point pre-commit quality gate; deciding whether to split or simplify a function. NOT for model training, feature engineering, or pipeline setup."
argument-hint: "Describe your task: e.g. 'review this trainer for quality issues', 'check for dead code', 'diagnose silent scoring bug', 'why is my metric wrapper wrong'"
license: MIT
metadata:
    skill-author: eak
---

# ML Competition — Code Quality & Common Pitfalls

## Overview

This skill covers two defensive layers applied before finalizing any Python change:

1. **Coding rules** — six non-negotiable standards that prevent the most common classes of bugs in fast-iteration competition codebases: no dead code, clear contracts, single responsibility, explicit types and names, predictable data handling, structured logging
2. **Common pitfalls** — 16 production bugs with ❌ BAD / ✅ GOOD patterns; read before finalizing any component

---

## Coding Standards — Non-Negotiable

Apply to every `src/*.py` and `scripts/*.py` file. These are hard quality gates, not suggestions.

| Rule | What it means |
|------|---------------|
| **No dead code** | No unused imports, variables, parameters, or uncalled private helpers (`_name`) |
| **Clear contracts** | Every public function has explicit input/output docs; optional params state their default behavior |
| **Single responsibility** | Functions do one job — if a function loads data AND computes metrics AND logs, split it |
| **Explicit types and names** | Type hints on all signatures; use `train_df`, `class_weights` not `tmp`, `d`, `x1` |
| **Predictable data handling** | Validate required columns up front and `raise` with a specific message; never silently mutate |
| **Structured logging** | Use `logger.*`; no `print`, no commented-out debug blocks, no stale TODOs |

### Pre-commit quality gate — run before every Python change

1. No unused imports, variables, or parameters
2. No uncalled private helper functions (`_name`)
3. Function signatures match actual behavior (not aspirational behavior)
4. Error messages are specific and actionable
5. Logs are concise and useful for iteration debugging
6. Another engineer can understand what a function does in under 30 seconds — if not, rewrite or split

---

## Reference Files

| File | What it covers |
|------|----------------|
| [coding-rules.md](./references/coding-rules.md) | No dead code, clear contracts, single responsibility, explicit types, structured logging — with ❌ / ✅ examples |
| [common-pitfalls.md](./references/common-pitfalls.md) | 16 production bugs with ❌ / ✅ patterns — read before finalizing any component |

---

## See Also

| Skill | When to use it instead |
|-------|------------------------|
| `ml-competition` | Full pipeline overview, task type decision guide, first-principles checklist |
| `ml-competition-setup` | Project structure, RunConfig, process management |
| `ml-competition-features` | Feature engineering, validation strategy |
| `ml-competition-training` | Base model training, competition metrics, output format |
| `ml-competition-tuning` | Optuna hyperparameter tuning |
| `ml-competition-advanced` | Pseudo-labeling, ensemble, post-processing, experiment tracking |
