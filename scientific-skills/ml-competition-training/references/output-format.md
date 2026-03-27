# Output Format & Submission Prediction Type

## Overview

The #1 source of silent score destruction in ML competitions is submitting the wrong prediction type for the competition metric. A model with 0.91 AUC submitting "Yes"/"No" strings instead of probabilities scores ~0.5 — indistinguishable from random — with no error message. This is preventable with a single check before writing any submission code.

This file covers: the core rule (metric determines prediction type, not `train.csv` target dtype); a complete metric-to-prediction-type table; sklearn call patterns for every metric; per-task submission format examples; OOF collection patterns; and the scout checklist to complete before writing submission code.

**The most dangerous confusion:** `train.csv` target column may contain strings ("Yes"/"No", "cat"/"dog") even when the submission requires floats (probabilities). The training target type and the submission prediction type are independent. Always derive the submission format from `sample_submission.csv` and the leaderboard metric name — never from `train.csv`.

**When to use:** Immediately after reading the competition README, before writing any `predict_proba` or `model.predict` call. Also consult when an apparently good model is scoring near 0.5 AUC on the leaderboard despite strong OOF — prediction type mismatch is the first diagnosis to rule out.

---

## The Core Rule

> **The metric determines the prediction type. The target column values in `train.csv` do NOT.**

| Metric | Prediction type | Submission column values |
|---|---|---|
| AUC-ROC, ROC AUC | **Probability** (float 0–1) | `0.0`–`1.0` |
| Average Precision, PR-AUC | **Probability** (float 0–1) | `0.0`–`1.0` |
| Log Loss (binary or multi) | **Probability** (float 0–1 per class) | `0.0`–`1.0` |
| Brier Score | **Probability** (float 0–1) | `0.0`–`1.0` |
| Accuracy, F1, Precision, Recall | **Class label** (as-is) | `0`, `1`, `Yes`, `No`, … |
| Cohen's Kappa (simple) | **Class label** (as-is) | integer or string class |
| Quadratic Weighted Kappa | **Integer class label or ordinal** | `0`, `1`, `2`, … |
| RMSE, MAE, RMSLE, MAPE | **Continuous value** | float |
| RMSLE | **Continuous value ≥ 0** | non-negative float |
| MAP@K, NDCG@K | **Ranked list** (space-separated IDs) | `id1 id2 id3` |
| Mean Columnwise Log Loss | **Probability per class column** | one float col per class |
| Multi-class log loss | **Probability per class column** | one float col per class |
| Mean Average Precision (COCO) | **Bounding box + score** | task-specific |

---

## Metric → sklearn Call

```python
from sklearn.metrics import (
    roc_auc_score,          # AUC-ROC   — expects proba
    average_precision_score, # AP/PR-AUC — expects proba
    log_loss,                # Log loss   — expects proba
    accuracy_score,          # Accuracy   — expects label
    f1_score,                # F1         — expects label
    mean_squared_error,      # RMSE       — expects value
    mean_absolute_error,     # MAE        — expects value
    cohen_kappa_score,       # Kappa      — expects label
)
import numpy as np

# ✅ AUC-ROC — binary, positive-class probability
roc_auc_score(y_true, y_pred_proba)            # y_pred_proba: shape (n,)

# ✅ AUC-ROC — multiclass
roc_auc_score(y_true, y_pred_proba_matrix,
              multi_class="ovr", average="macro")  # shape (n, n_classes)

# ✅ Log loss — binary
log_loss(y_true, y_pred_proba)

# ✅ Log loss — multiclass
log_loss(y_true, y_pred_proba_matrix)          # shape (n, n_classes)

# ✅ RMSE
np.sqrt(mean_squared_error(y_true, y_pred))

# ✅ RMSLE — apply log1p transform before RMSE
np.sqrt(mean_squared_error(np.log1p(y_true), np.log1p(y_pred.clip(0))))

# ✅ Quadratic Weighted Kappa — round predictions to integer classes first
from sklearn.metrics import cohen_kappa_score
cohen_kappa_score(y_true, y_pred_rounded, weights="quadratic")
```

---

## Identifying the Metric from Competition Materials

When reading a README or competition description, look for these phrases:

| Phrase | Metric | Prediction type |
|---|---|---|
| "Area Under the ROC Curve", "AUC", "ROC AUC", "AUROC" | AUC-ROC | probability |
| "Average Precision", "PR-AUC", "mAP" (non-detection) | AP | probability |
| "Log Loss", "Cross-Entropy", "Logloss" | Log loss | probability |
| "Brier Score" | Brier | probability |
| "Accuracy", "Correct Classification Rate" | Accuracy | class label |
| "F1 Score", "F-Measure", "F-beta" | F1 | class label |
| "Root Mean Square Error", "RMSE" | RMSE | value |
| "Mean Absolute Error", "MAE" | MAE | value |
| "Root Mean Square Log Error", "RMSLE" | RMSLE | value ≥ 0 |
| "Mean Absolute Percentage Error", "MAPE" | MAPE | value |
| "Quadratic Weighted Kappa", "QWK", "κ" | QWK | ordinal label |
| "Mean Columnwise Log Loss", "MCLL" | MCLL | proba matrix |
| "Mean Average Precision @ K", "MAP@K" | MAP@K | ranked list |
| "NDCG", "Normalized Discounted Cumulative Gain" | NDCG | score or rank |

---

## Submission Format by Task Type

### Binary classification with probability metric (AUC, log-loss, brier)

```csv
id,target
1,0.23
2,0.87
3,0.04
```
- `target` column: **float in [0, 1]** — never "Yes"/"No", never 0/1 integers only
- Use `model.predict_proba(X)[:, 1]` — NOT `model.predict(X)`

### Binary classification with label metric (accuracy, F1)

```csv
id,target
1,1
2,0
3,1
```
- `target` column: original class labels — check `sample_submission.csv` for exact dtype

### Multiclass classification with per-class probability (log-loss)

```csv
id,class_0,class_1,class_2
1,0.1,0.7,0.2
2,0.5,0.3,0.2
```
- One column per class, values sum to 1.0 per row
- Column order must match `sample_submission.csv` column order exactly

### Multiclass with single-label (accuracy, QWK)

```csv
id,target
1,2
2,0
3,1
```

### Regression

```csv
id,target
1,142.7
2,9.3
```
- RMSLE: clip predictions to ≥ 0; do NOT log-transform before saving

---

## OOF Collection Pattern

> **Always collect full OOF arrays, never fold-by-fold averages.**

```python
import numpy as np

# Binary probability metric (AUC, log-loss)
oof = np.zeros(len(y))
for tr, val in folds:
    model.fit(X[tr], y[tr])
    oof[val] = model.predict_proba(X[val])[:, 1]
# score:
from sklearn.metrics import roc_auc_score
score = roc_auc_score(y, oof)

# Multiclass probability metric (log-loss)
oof = np.zeros((len(y), n_classes))
for tr, val in folds:
    model.fit(X[tr], y[tr])
    oof[val] = model.predict_proba(X[val])   # shape (n_val, n_classes)
score = log_loss(y, oof)

# Regression (RMSE)
oof = np.zeros(len(y))
for tr, val in folds:
    model.fit(X[tr], y[tr])
    oof[val] = model.predict(X[val])
score = np.sqrt(mean_squared_error(y, oof))
```

---

## Scout Checklist — Submission Format Section

When writing `DATA_BRIEFING.md`, always answer ALL of the following:

1. **What is the exact metric name?** (Quote directly from the README.)
2. **Is it a probability metric?** (AUC, log-loss, AP, brier → YES → float predictions required.)
3. **What are the exact column names in `sample_submission.csv`?** (List them.)
4. **What are the value types in those columns?** (float / int / string — check `sample_submission.csv` dtype, not `train.csv` target dtype.)
5. **What do 2–3 example rows look like?** (Copy verbatim from `sample_submission.csv`.)
6. **Quote the README evaluation section** confirming the format.

> ⚠️ **The target in `train.csv` may be strings ("Yes"/"No", "cat"/"dog") even when the submission requires probabilities (floats).** Train target strings are labels for *training* — the submission column is for *predictions*. Always derive the submission format from `sample_submission.csv` and the metric name, never from the train target dtype.

---

## See Also

| File | Why |
|------|-----|
| [competition-metrics.md](./competition-metrics.md) | Metric implementations that match the leaderboard formula identified here |
| [submission-postprocessing.md](./submission-postprocessing.md) | Calibration and clipping applied after the correct prediction type is confirmed |
