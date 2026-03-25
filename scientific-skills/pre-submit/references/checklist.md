# ML Quality Checklist

## Overview

This checklist is the first of three mandatory validation checks. It covers four categories — data leakage, metric correctness, submission format, and robustness — each categorized as CRITICAL or Important.

**CRITICAL items must all pass before reporting any OOF score or generating a submission.** A single CRITICAL failure invalidates your OOF score, even if training completed successfully.

**How to use:** Work through each section top-to-bottom on every new script or after any change to feature engineering, CV splits, or metric code. Mark items as your code is verified. The most common CRITICAL failure in practice is target encoding computed on the full training set (leaks the target into features) — check this first.

---

Fix every **CRITICAL** item before reporting results.

---

## CRITICAL — Data Leakage

- [ ] Target-based encodings (mean / target encoding) computed **inside** each CV fold — never on full train.
- [ ] Temporal features (lag, rolling stats) use only past data — no future leakage.
- [ ] No test rows appear in any training fold.
- [ ] StandardScaler / other transformers fit on **train fold only**, applied to val/test.
- [ ] `train_test_split` NOT used instead of proper k-fold CV.

## CRITICAL — Metric Correctness

- [ ] OOF metric computed on the **full OOF array**, not fold-by-fold averages.
- [ ] Metric function matches the competition definition exactly (`average='macro'` vs `'binary'` etc.).
- [ ] For probability metrics (AUC, log-loss): predictions are **probabilities**, not class labels.
- [ ] Metric direction (maximize / minimize) respected when comparing scores.

## CRITICAL — Submission Format

- [ ] Submission CSV column names match `sample_submission.csv` exactly.
- [ ] Submission row count matches `sample_submission.csv` exactly.
- [ ] No NaN or Inf in prediction column.
- [ ] File saved to the path reported in `submission_file`.

## CRITICAL — Open-Ended Tasks

- [ ] Deliverable runs end-to-end without crashes.
- [ ] All requirements from `README.md` are implemented.
- [ ] Dependencies declared in `pyproject.toml`.
- [ ] Reproducible from scratch: `uv sync && uv run ...`.
- [ ] Submission artifact exists at the reported path.

## Important — Robustness

- [ ] No hard-coded file paths — use `pathlib` and config variables.
- [ ] Random seeds set (`random_state=42`, `np.random.seed(42)`).
- [ ] OOF score printed as `OOF {metric}: {score:.6f}` — evaluator parses this.
- [ ] Script runs end-to-end without manual intervention.

## Quality Score (open-ended tasks)

| Score | Meaning |
| --- | --- |
| 90–100 | Exceeds requirements; polished, documented, tested |
| 70–89 | Meets all stated requirements; no major gaps |
| 50–69 | Meets most requirements; some gaps or rough edges |
| 30–49 | Partial implementation; core functionality works |
| 0–29 | Incomplete; significant requirements unmet |

---

## Failure Patterns and Fixes

| Symptom | Most likely CRITICAL failure | Fix |
|---------|------------------------------|-----|
| OOF AUC >> LB AUC (gap ≥ 0.05) | Target encoding fit on full train | Compute inside CV fold — see [validation-strategy.md](../../ml-competition/references/validation-strategy.md) |
| OOF improves but LB flat | Feature uses future/test info | Check temporal and test-row leakage items |
| Submission rejected / score 0.0 | Column name, row count, or NaN | Re-run Workflow 2 submission validation |
| OOF flat but LB improves | OOF metric formula wrong | Re-check metric implementation against leaderboard formula |
| OOF AUC suspiciously high (> 0.98) on tabular | Target directly or proxy in features | Check ID columns, datetime ordering, calculated targets |

---

## See Also

| File | Why |
|------|-----|
| [../SKILL.md](../SKILL.md) | Workflows 1-3 that use this checklist — run all three in order |
| [../../ml-competition/references/validation-strategy.md](../../ml-competition/references/validation-strategy.md) | Correct CV split and OOF accumulation patterns |
| [../../ml-competition/references/output-format.md](../../ml-competition/references/output-format.md) | Metric → prediction type — governs the "probability vs label" submission check |
| [../../ml-competition/references/common-pitfalls.md](../../ml-competition/references/common-pitfalls.md) | Extended list of production bugs that correspond to checklist failures |
