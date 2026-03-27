---
name: ml-competition-setup
description: "Set up a tabular ML competition project from scratch or review existing setup. Use when: creating the package layout and directory structure; implementing the RunConfig dataclass singleton; writing the YAML-driven orchestrator with resume-by-existence logic; managing training processes (pre-flight checks, PID tracking, launch/wait/kill patterns, diagnosing fast exits). NOT for feature engineering, model training, or tuning."
argument-hint: "Describe your task: e.g. 'scaffold competition project', 'implement RunConfig', 'debug why train.py reruns completed steps', 'fix process management'"
license: MIT
metadata:
    skill-author: eak
---

# ML Competition — Setup & Process Management

## Overview

This skill covers two foundational concerns that every other part of the pipeline depends on:

1. **Project structure & config** — canonical package layout, `RunConfig` singleton, `config.yaml`, YAML orchestrator with resume-by-existence logic, OOF/submission/tuning naming conventions
2. **Process management** — pre-flight checks before every training launch, PID tracking, launch/wait/kill patterns, diagnosing deceptively fast exits

**Two governing principles:**
- **Single source of truth for config**: all user-tunable knobs live in `config.yaml` and are loaded into `RunConfig` at startup; no magic constants in trainer files
- **Resume by existence**: the orchestrator skips any step whose output file already exists — interrupt and rerun at any point without duplicating work or corrupting results

---

## Separation of Concerns — Non-Negotiable

Every layer owns exactly one thing. Never blur these boundaries.

| Layer | Owns | Key files |
|-------|------|-----------|
| **Config** | All user-tunable knobs | `base/config.py` (RunConfig dataclass + singleton) |
| **Features** | Raw → engineered DataFrames, versioned pkl cache | `base/features.py`, `cache/features_vN.pkl` |
| **Matrices** | NumPy/DataFrame arrays per model family | `base/features.py::build_model_matrices()` |
| **Metrics** | Competition metric, framework-specific wrappers | `base/metrics.py` |
| **Trainers** | Stateless fold-loop engines, injectable params | `base/lgb_trainer.py`, `xgb_trainer.py`, `nn_trainer.py` |
| **Entrypoints** | Load → call trainer → score → save OOF | `train/cat.py`, `lgb.py`, `xgb.py`, `nn.py` |
| **Tuning** | Optuna objectives, run_study, save JSON | `tune/tune_*.py`, `tune.py` |
| **Pseudo** | Retrain on train + pseudo-labeled test | `train/pseudo.py` |
| **Meta** | Weighted blend / stacking / gating on OOF | `train/meta.py`, `train/meta_gating.py` |
| **Orchestrator** | YAML-driven sequential runner + resume logic | `train.py` |

---

## Process Management — Non-Negotiable

Training scripts can run for hours. The most expensive mistake is launching a duplicate process.

### Before every training launch — pre-flight check

```bash
# 1. Is training already running?
RUNNING_PIDS=$(pgrep -f "python scripts/train.py" 2>/dev/null)
[ -n "$RUNNING_PIDS" ] && echo "⚠️  Already running — PIDs: $RUNNING_PIDS" && exit 0

# 2. Are artifacts already fresh? (< 5 min old)
[ -f "artifacts/oof.npy" ] && \
  ARTIFACT_AGE=$(( $(date +%s) - $(stat -c %Y artifacts/oof.npy) )) && \
  [ $ARTIFACT_AGE -lt 300 ] && echo "✅ Artifacts are fresh — skip retraining"
```

### Correct launch pattern

```bash
# TRAIN_PID=$! MUST be on its own line immediately after & — not on the same line
nohup uv run python scripts/train.py > train.log 2>&1 &
TRAIN_PID=$!
echo "Training started — PID: $TRAIN_PID"

# Wait loop
while kill -0 $TRAIN_PID 2>/dev/null; do sleep 60; tail -5 train.log; done
echo "✅ Done"; tail -50 train.log
```

### Hard rules

- **Never launch without the pre-flight check** — two processes writing the same `oof.npy` corrupt results silently
- **`TRAIN_PID=$!` must be on its own line** — assigning on the same compound line as `&` captures an empty string
- **Never pipe `train.py` to `head`** — `head` closes the pipe and kills the process; always redirect to a log file
- **Fast exit ≠ failure** — check artifact timestamps and CPU load before concluding anything went wrong
- **Kill before relaunch** — if code changed and an old process is running, kill it first, then relaunch

See [process-management.md](./references/process-management.md) for all four workflow patterns (pre-flight, launch/wait, kill/relaunch, diagnose fast exit).

---

## Reference Files

| File | What it covers |
|------|----------------|
| [project-structure.md](./references/project-structure.md) | Package layout, `RunConfig` singleton, YAML orchestrator, resume-by-existence logic |
| [process-management.md](./references/process-management.md) | Pre-flight checks, PID tracking, launch/wait/kill patterns, fast-exit diagnosis |

---

## See Also

| Skill | When to use it instead |
|-------|------------------------|
| `ml-competition` | Full pipeline overview, task type decision guide, first-principles checklist |
| `ml-competition-features` | Feature engineering, validation strategy |
| `ml-competition-training` | Model training, competition metrics, output format |
| `ml-competition-tuning` | Optuna hyperparameter tuning |
| `ml-competition-advanced` | Pseudo-labeling, ensemble, post-processing, experiment tracking |
| `ml-competition-quality` | Coding rules, common pitfalls |
