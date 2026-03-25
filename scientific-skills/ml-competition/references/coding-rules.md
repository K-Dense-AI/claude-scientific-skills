# Coding Rules

## Overview

These are Python coding standards applied as a hard quality gate before finalizing any change in `src/*.py` or `scripts/*.py`. They exist because competition codebases are iterated quickly under deadline pressure — dead code accumulates, function signatures drift from their behavior, and debug artifacts get committed. Each rule below directly prevents a class of bugs observed in production pipelines.

**Six non-negotiable rules:** no dead code; clear contracts; single responsibility; explicit types and names; predictable data handling; structured logging. Each failure mode is illustrated with a ❌ BAD and ✅ GOOD example.

**How to use this file:** Before committing any Python change, run through the 6-point verification checklist at the end. If another engineer cannot understand what a function does in under 30 seconds, the function must be rewritten or split.

---

# Role and Objective

Apply these Python coding standards as a hard quality gate before finalizing any change in `src/*.py` or `scripts/*.py`.

Begin with a concise checklist (3-7 bullets) of what you will do; keep items conceptual, not implementation-level.

# Instructions

## When to Apply

Use these standards when:
- Writing or editing `src/*.py` or `scripts/*.py`
- Reviewing teammate or agent code for maintainability
- Refactoring quickly under iteration pressure
- Investigating issues such as "why is this code confusing?" or "why is this helper unused?"

## Non-Negotiable Rules

### 1. No Dead Code
- No unused imports
- No unused local variables
- No unused function parameters unless required by an interface
- No private helper (`_name`) that is never called

If a helper is intentionally staged for later, do **not** leave it in production files. Remove it and re-add it when needed.

### 2. Clear Contracts
- Every public function must have a clear input/output contract
- Optional parameters must have explicit behavior in the docstring
- If a parameter is currently ignored, remove it now or wire it correctly now

### 3. Single Responsibility
- Functions should do one job
- If a function both loads data and computes metrics and logs summaries, split it
- Prefer small composable helpers over monolithic functions

### 4. Explicit Types and Names
- Use type hints on function signatures
- Use descriptive names such as `train_df` and `class_weights`, not vague names such as `tmp`, `d`, or `x1`
- Avoid one-letter variables except loop indices in tiny scopes

### 5. Predictable Data Handling
- Validate required columns up front and fail fast
- Keep dtype casting deterministic
- Do not silently mutate unrelated columns

### 6. Logging and Errors
- Use structured, useful log messages
- Raise actionable errors for missing columns, files, or types; avoid generic failures
- No debug leftovers such as `print`, commented code blocks, or stale TODOs

## Anti-Patterns to Avoid

```python
# ❌ BAD: private helper is never used

def _cast_categorical_columns(df: pd.DataFrame) -> pd.DataFrame:
    ...

# helper is defined but never called
```

```python
# ❌ BAD: parameter is accepted but ignored

def load_train(fold: int = 0):
    df = pd.read_csv("train.csv")
    return df  # fold has no effect
```

```python
# ✅ GOOD: remove or implement parameter behavior

def load_train() -> tuple[pd.DataFrame, pd.Series]:
    ...
```

# Planning and Verification

If editing code, state assumptions explicitly, follow existing repository style and patterns, and create or run the minimal relevant tests or checks where possible. If tests or checks cannot be run, say so plainly and identify what should be validated later.

Before finishing a change, verify all of the following:
1. No unused imports, variables, or parameters
2. No uncalled private helper functions
3. Function signatures match actual behavior
4. Docstrings describe real behavior, not aspirational behavior
5. Error messages are specific and actionable
6. Logs are concise and useful for iteration debugging

After each substantive edit or validation step, provide a brief validation note stating what was checked and whether the result passed or needs correction.

# Practical Standard

If another engineer cannot understand what a function does in under 30 seconds, rewrite the function or split it.

# Verbosity

Default to concise implementation and review guidance. For code, prefer high verbosity in readability: clear names, explicit typing, comments where helpful, and straightforward control flow.

# Persistence

Keep applying these standards until the Python change is fully resolved. Attempt a conservative first pass autonomously when the required context is present; ask for clarification only if a missing detail blocks correct application of these standards.
Do not leave behind known dead code, confusing contracts, or avoidable technical debt.

# Stop Conditions

Finish only when the code meets all applicable rules above for maintainability, clarity, and production safety.
---

## See Also

| File | Why |
|------|-----|
| [common-pitfalls.md](./common-pitfalls.md) | Production bugs caused by contract violations and dead code |
| [process-management.md](./process-management.md) | Safe launch patterns for training scripts that follow these coding standards |
