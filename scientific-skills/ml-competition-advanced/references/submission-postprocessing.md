# Submission Post-Processing

## Overview

Post-processing transforms raw model outputs into submission-ready predictions. It can recover 0.001–0.003 on the leaderboard without any retraining — but only when the right technique is applied to the right situation. Applied blindly, calibration and clipping can hurt more than they help.

This file covers: probability calibration (Platt scaling, isotonic regression); OOF-optimised prediction clipping; domain constraint enforcement; and the YAML-configurable calibration step. Each technique includes a decision table indicating when it helps vs. when it is harmful.

**The most dangerous misuse:** applying isotonic calibration to a 4+ model ensemble. The ensemble already softens extreme predictions through averaging — additional calibration usually over-corrects and produces a score drop (−0.006 OOF was observed in a 4-model binary blend). Always verify on OOF before enabling calibration. The YAML option `calibrate.enabled` defaults to `false` for this reason.

**When to use:** As the final step before generating the submission file, after ensemble meta-learning. Also consult when a single model's predictions cluster near 0 or 1 (overconfidence signal).

---

## Why Post-Processing Matters

Raw model outputs are often miscalibrated (overconfident probabilities) and may violate domain constraints (monotone orderings, forbidden ranges). Post-processing can add 0.001–0.003 to LB without retraining.

**Order of operations:**
1. Calibrate probabilities (if needed)
2. Apply domain constraints (if any)
3. Clip to submission range
4. Optimise per-target clip bounds on OOF

---

## Probability Calibration

Raw model probabilities often deviate from true probabilities, especially for tree models trained on imbalanced data.

### Platt Scaling (Logistic Calibration)
```python
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression

def platt_calibrate(oof_preds: np.ndarray, y_true: np.ndarray,
                    test_preds: np.ndarray) -> np.ndarray:
    """Fit sigmoid calibration on OOF; apply to test."""
    lr = LogisticRegression(C=1.0)
    lr.fit(oof_preds.reshape(-1, 1), y_true)
    return lr.predict_proba(test_preds.reshape(-1, 1))[:, 1]
```

### Isotonic Regression (non-parametric)
```python
from sklearn.isotonic import IsotonicRegression

def isotonic_calibrate(oof_preds, y_true, test_preds):
    ir = IsotonicRegression(out_of_bounds="clip")
    ir.fit(oof_preds, y_true)
    return ir.predict(test_preds)
```

### When to calibrate
| Signal | Action |
|--------|--------|
| OOF scores well but predictions cluster near 0 or 1 | Calibrate with Platt or isotonic |
| Log-loss metric is prominent in formula | Always try on OOF first — only commit if gain ≥ 0.001 |
| AUC-only competition | Calibration won't help (ranking = invariant to monotone transforms) |
| Ensemble of 4+ models | **Usually harmful** — ensemble already softens extremes. E.g. −0.006 OOF was observed from isotonic calibration on a 4-model binary blend. Verify on OOF first; keep disabled by default. |
| Single model submission | Calibration more likely to help — verify on OOF first |

---

## Prediction Clipping

### OOF-optimised per-target clip

Find the upper-bound clip that maximises OOF competition score:

```python
from scipy.optimize import minimize_scalar

def optimize_clip(oof_preds: np.ndarray, y_true: np.ndarray,
                  metric_fn, min_bound=0.5) -> float:
    """Find clip upper bound that maximises metric_fn on OOF."""
    def neg_score(hi):
        clipped = np.clip(oof_preds, 1e-7, hi)
        return -metric_fn(y_true, clipped)

    result = minimize_scalar(neg_score, bounds=(min_bound, 1.0), method="bounded")
    return result.x

# Per-target clip
clip_map = {}
for t in TARGETS:
    clip_map[t] = optimize_clip(oof[t], y_dict[t], competition_score)
    test_preds[t] = np.clip(test_preds[t], 1e-7, clip_map[t])
```

**Note:** Optimise clipping on OOF — do not use test predictions to decide clips (that's future leakage).

### Safe clipping (no OOF optimisation)
```python
# For probabilistic outputs: avoid exactly 0.0 or 1.0 (log-loss goes infinite)
preds = np.clip(preds, 1e-7, 1 - 1e-7)
```

---

## Domain Constraints

### Monotone hierarchy (ordered probabilities)
When logic requires P(A) ≤ P(B) ≤ P(C):
```python
def enforce_monotone(df, cols: list[str]) -> pd.DataFrame:
    """
    Enforce non-decreasing order across cols for each row.
    E.g. cols = ["P_7d", "P_30d", "P_90d"]
    """
    for col_lo, col_hi in zip(cols[:-1], cols[1:]):
        mask = df[col_lo] > df[col_hi]
        df.loc[mask, col_lo] = df.loc[mask, col_hi]
    return df
```

### Zeroing out impossible predictions
If certain rows have known zero probability (no exposure, no eligibility):
```python
zero_mask = df["eligibility_flag"] == 0
preds[zero_mask] = 0.0
```

### Regression floor/ceiling
```python
# Force predictions into physically meaningful range
preds = np.clip(preds, 0.0, None)   # non-negative (counts, prices)
preds = np.clip(preds, 0.0, 1.0)    # bounded range
```

---

## Submission File Construction

```python
def make_submission(test_preds: dict, test_df: pd.DataFrame, name: str,
                    clip_map: dict | None = None) -> str:
    """
    Build and save submission CSV.
    test_preds: {target: np.array} — one array per target
    """
    os.makedirs(cfg.sub_dir, exist_ok=True)
    sub = test_df[["id"]].copy()   # always start with the ID column

    for t in TARGETS:
        raw = test_preds[t].copy()
        hi  = clip_map[t] if clip_map else (1 - 1e-7)
        sub[t] = np.clip(raw, 1e-7, hi)

    path = os.path.join(cfg.sub_dir, f"{name}.csv")
    sub.to_csv(path, index=False)
    return path
```

**Verify submission format before uploading:**
```python
sample = pd.read_csv("data/SampleSubmission.csv")
assert list(sub.columns) == list(sample.columns), "Column mismatch!"
assert len(sub) == len(sample), "Row count mismatch!"
assert sub["id"].equals(sample["id"]), "ID order mismatch!"
```

---

## Safe Experiment Protocol

Post-processing changes the submission without changing the model. Still follow the same discipline:

1. **One change at a time** — calibration alone, clip alone, constraint alone
2. **Always compute on OOF first** — never tune post-processing on test predictions
3. **Log the OOF gain** before submitting — if OOF gain < 0.0005 it's within noise
4. **Keep the raw (unprocessed) submission** — you may want to revert

```
# Safe order:
# 1. Submit base + no post-processing → establish baseline
# 2. Try calibration → compute OOF gain → submit if gain > noise
# 3. Try optimal clips → compute OOF gain → submit if gain > noise
# 4. Apply domain constraints → submit if logically required
```

---

## See Also

| File | Why |
|------|-----|
| [output-format.md](./output-format.md) | Submission column types and format derived from `sample_submission.csv` |
| [experiment-tracking.md](./experiment-tracking.md) | OOF gain must justify any post-processing before spending a submission slot |
| [ensemble-meta.md](./ensemble-meta.md) | Post-processing is applied to the ensemble output, not to base model outputs |
