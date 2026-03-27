---
name: ml-competition-pre-submit
description: Pre-submission quality gate for ML competition pipelines. Runs three checks before any result is reported, code review (data leakage, CV contamination, metric errors), submission CSV format validation, and adversarial train/test distribution shift detection. Invoke before reporting any OOF score or submitting predictions.
license: MIT
metadata:
    skill-author: eak
---
# Validation Skill

## Overview

This skill is a mandatory pre-submission quality gate for tabular ML competition pipelines. It catches the three most costly bugs before they reach the leaderboard: target leakage inflating OOF scores, submission CSV format rejections, and OOF–LB correlation collapse caused by distribution shift.

**Three checks, always in this order — skip none:**

1. **Code review** — data leakage, CV contamination, metric implementation errors (`references/checklist.md`)
2. **Submission format validation** — column names, row count, NaN/Inf, value ranges (Workflow 2)
3. **Adversarial distribution shift** — train vs test AUC to detect shift and identify offending features (`scripts/adversarial_validation.py`)

**Critical principle:** no OOF score is meaningful until all CRITICAL items in the checklist pass. A leaky 0.95 AUC that passes no leakage checks is worthless — it will not exceed 0.70 on the leaderboard.

**When to use:** Before reporting any OOF score or submitting predictions — every single iteration, without exception.

---

## When to Use

- Before reporting any OOF score or submitting predictions — every single iteration.
- When OOF looks suspiciously high relative to the task difficulty.
- When LB score lags OOF by more than 0.01 (CV-LB gap).
- After writing any new feature batch — verify no shift introduced.
- At competition start — run adversarial AUC baseline to calibrate CV strategy.

## Critical Rules

### ✅ DO

- **Always run all three checks** — skip one and you will submit a broken or leaky result.
- **Fix CRITICAL items before reporting** — no OOF score is meaningful until leakage is ruled out.
- **Run adversarial validation once at the start** — know your AUC baseline before adding any features.
- **Check submission format programmatically** — never eyeball column names; trailing spaces and case differences cause silent rejections.
- **Print OOF as `OOF {metric}: {value:.6f}`** — the evaluator agent parses this exact format.

### ❌ DON'T

- **Don't report OOF before running the leakage checklist** — a leaky 0.95 AUC is worthless.
- **Don't skip submission validation** — wrong column order, NaN values, and mismatched row counts are silent failures on most platforms.
- **Don't wait for LB feedback to detect shift** — adversarial AUC catches it before you waste submissions.
- **Don't drop features just because adversarial AUC is high** — first confirm they hurt CV; sometimes shift is benign.

## Anti-Patterns (NEVER)

```python
# ❌ BAD: Computing OOF metric as average of fold scores
fold_scores = []
for tr, val in cv.split(X, y):
    model.fit(X[tr], y[tr])
    fold_scores.append(roc_auc_score(y[val], model.predict_proba(X[val])[:, 1]))
oof_score = np.mean(fold_scores)  # WRONG — biased by unequal fold sizes

# ✅ GOOD: Compute on the full OOF array
oof = np.zeros(len(y))
for tr, val in cv.split(X, y):
    model.fit(X[tr], y[tr])
    oof[val] = model.predict_proba(X[val])[:, 1]
oof_score = roc_auc_score(y, oof)  # correct
```

```python
# ❌ BAD: Submit without checking column names
submission.to_csv("submission.csv", index=False)

# ✅ GOOD: Validate against sample before saving
sample = pd.read_csv("data/sample_submission.csv")
assert list(submission.columns) == list(sample.columns), \
    f"Column mismatch: {submission.columns.tolist()}"
assert len(submission) == len(sample), \
    f"Row count mismatch: {len(submission)} vs {len(sample)}"
assert submission.isnull().sum().sum() == 0, "NaN in submission"
submission.to_csv("artifacts/submission.csv", index=False)
```

## Workflows

The three workflows below map to the three mandatory checks. Run them in order.

---

## Workflow 1 — Code Review Checklist

Run through `references/checklist.md` on every new script. Non-negotiable items:

### Data Leakage

- [ ] Target encodings computed **inside** each CV fold — never on full train.
- [ ] Temporal features use only past data — no future leakage.
- [ ] No test rows in any training fold.
- [ ] StandardScaler / transformers fit on **train fold only**, applied to val/test.
- [ ] `train_test_split` NOT used instead of proper k-fold CV.

### Metric Correctness

- [ ] OOF metric computed on the **full OOF array**, not fold-by-fold averages.
- [ ] Metric function matches competition definition exactly (`average='macro'` vs `'binary'` etc.).
- [ ] For probability metrics (AUC, log-loss): predictions are **probabilities**, not class labels.
- [ ] Metric direction (maximize/minimize) consistent throughout.

### Robustness

- [ ] OOF score printed as `OOF {metric}: {score:.6f}` — evaluator agent parses this.
- [ ] No hard-coded paths — use `pathlib` and config variables.
- [ ] Random seeds set (`random_state=42`, `np.random.seed(42)`).

---

## Workflow 2 — Submission File Validation

```python
import pandas as pd
import numpy as np

def validate_submission(submission_path="artifacts/submission.csv",
                        sample_path="data/sample_submission.csv") -> str:
    submission = pd.read_csv(submission_path)
    sample     = pd.read_csv(sample_path)

    errors = []
    if list(submission.columns) != list(sample.columns):
        errors.append(f"Column mismatch — got {submission.columns.tolist()}, expected {sample.columns.tolist()}")
    if len(submission) != len(sample):
        errors.append(f"Row count mismatch — got {len(submission)}, expected {len(sample)}")
    if submission.isnull().any().any():
        errors.append(f"NaN in: {submission.columns[submission.isnull().any()].tolist()}")

    # Check prediction column value ranges
    pred_col = sample.columns[-1]
    vals = submission[pred_col]
    if vals.dtype in [float, "float64", "float32"]:
        if (vals < 0).any() or (vals > 1).any():
            if (vals < 0).any() or (vals > 1.001).any():
                errors.append(f"Values out of [0,1]: min={vals.min():.4f}, max={vals.max():.4f}")

    return "VALID" if not errors else "INVALID: " + "; ".join(errors)

print(validate_submission())
```

---

## Workflow 3 — Adversarial Validation

```bash
uv run python scripts/adversarial_validation.py \
    --train data/train.csv --test data/test.csv --target <target_col>
# Prints: AUC, verdict, top-20 leaking features
# Saves:  artifacts/adversarial_weights.npy

# Optional: specify a different classifier
uv run python scripts/adversarial_validation.py \
    --train data/train.csv --test data/test.csv --target <target_col> --clf rf
```

Or inline from Python:

```python
from scripts.adversarial_validation import run_adversarial_validation
result = run_adversarial_validation(
    train_path="data/train.csv",
    test_path="data/test.csv",
    target_col="target",
)
# {"auc": 0.58, "verdict": "⚠️ Mild shift ...", "top_features": {...}}

# Custom classifier (any sklearn-compatible object):
from sklearn.ensemble import RandomForestClassifier
result = run_adversarial_validation(
    train_path="data/train.csv",
    test_path="data/test.csv",
    target_col="target",
    clf=RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=-1),
)
```

### AUC interpretation

| AUC        | Verdict           | Action                                         |
| ---------- | ----------------- | ---------------------------------------------- |
| 0.50–0.55 | ✅ No shift       | Proceed normally                               |
| 0.55–0.65 | ⚠️ Mild shift   | Check top features; monitor LB-OOF gap         |
| 0.65–0.80 | ❌ Moderate shift | Drop or transform top leaking features         |
| 0.80–1.00 | 🚨 Severe shift   | Likely ID/time leak — investigate immediately |

### Using adversarial sample weights

```python
weights = np.load("artifacts/adversarial_weights.npy")
model.fit(X_train, y_train, sample_weight=weights)
```

---

## Common Pitfalls and Solutions

### The "Leaky OOF" Problem

OOF AUC is 0.95 but leaderboard gives 0.70. Classic target leakage.

**Check:** `references/checklist.md` data leakage section. Most common culprit: target encoding fit on full train rather than per fold.

### The "False Negative" Adversarial AUC

Adversarial AUC is 0.52 but LB-OOF gap is still large. Shift exists but in non-numeric columns not included in the adversarial classifier.

**Fix:** Include string frequency encodings or hash-encoded string features in the adversarial dataset.

### The "Silent Submission Failure"

Platform returns 0.0 or error despite local validation passing. Common cause: integer IDs in prediction column instead of floats for probability competitions.

**Fix:** Always `submission[pred_col] = submission[pred_col].astype(float)` before saving.

---

## Reference Files

| File                                                          | What it covers                                                                                          |
| ------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| [checklist.md](./references/checklist.md)                        | CRITICAL and Important quality gates — data leakage, metric correctness, submission format, robustness |
| [adversarial_validation.py](./scripts/adversarial_validation.py) | Distribution shift detector: AUC, top leaking features, sample weights for retraining                   |

---

## See Also

| Skill | Why |
| ----- | --- |
| `ml-competition-features` | CV split strategy — GroupKFold / TimeSeriesSplit, OOF array accumulation |
| `ml-competition-advanced` | OOF vs LB divergence diagnosis — run after this checklist passes |
| `ml-competition-training` | Metric → prediction type table — governs what Workflow 2 checks |
| `ml-competition-quality` | 16 production bugs — most are variants of the leakage and metric bugs caught here |
