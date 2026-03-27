# Process Management

## Overview

Training scripts routinely run for 30 minutes to several hours. Without disciplined process management, the most common mistake is launching a second training process before the first has finished — two processes writing to the same `oof.npy` file, silently corrupting results with no error message.

This file provides the exact shell patterns for: pre-flight checks (is training already running? are artifacts already fresh?); launching and monitoring a training process; killing and relaunching after a code change; and diagnosing a fast or silent exit. Each workflow is a copy-paste shell block — use them verbatim.

**The single most common anti-pattern:** assigning `TRAIN_PID=$!` on the same line as the `&` command. Bash evaluates `$!` before the background job is registered, so `TRAIN_PID` is empty. Every subsequent `kill -0 $TRAIN_PID` returns immediately, the wait loop never waits, and the next agent step runs while training is still in progress — overwriting outputs.

**When to use:** Before every `uv run python scripts/train.py` call. This is not optional — it is a hard prerequisite, like checking for syntax errors before deploying.

---

Training scripts can take minutes to hours. The single most costly mistake is
launching a second process before the first has finished. This skill gives you
the exact pre-flight checks and wait patterns to avoid that.

## When to Use

- Before calling `uv run python scripts/train.py` — every time, no exceptions.
- When confused about whether training already ran (log looks short, artifacts have old timestamps).
- When you need to rerun training with a code change and must kill the old process first.

## Critical Rules

### ✅ DO

- **Always run the pre-flight check** before launching — one command prevents an entire class of bugs.
- **Capture PID with `TRAIN_PID=$!`** immediately after the `&` — on the very next line.
- **Wait with `kill -0 $TRAIN_PID`** — the only reliable synchronous wait pattern.
- **Kill before relaunch** — if an old process is running and you need to change code, kill it first.
- **Check artifact freshness by timestamp** — a fast script exit usually means the artifacts are already up-to-date, not that something failed.

### ❌ DON'T

- **Don't launch `scripts/train.py` without checking if it's already running** — you will get two processes writing to the same `oof.npy`.
- **Don't use background task IDs** (`TaskOutput`, `TaskStop`, `run_in_background: true`) — they are unreliable for long jobs and break the PID tracking pattern.
- **Don't infer "something went wrong" from a short `train.log`** — check artifact timestamps and CPU load first; the script may have exited because outputs are already up-to-date, or the kernel is buffering stdout because all CPUs are saturated.
- **Don't run `scripts/train.py` interactively (`| head -N`)** — this pipes stdout to head and kills the process when head closes; use `tee` or redirect to a log file.
- **Don't spawn a new training run to "diagnose" a fast exit** — diagnose first, run second.

## Anti-Patterns (NEVER)

```bash
# ❌ BAD: TRAIN_PID=$! on the same line as the & command
nohup uv run python scripts/train.py > train.log 2>&1 & TRAIN_PID=$!
# $! is NOT captured — bash evaluates $! before & has set it for this line.
# TRAIN_PID is empty. Every subsequent kill -0 $TRAIN_PID call silently exits
# immediately, so the wait loop never waits. The model then polls manually.

# ❌ BAD: Inline echo — same problem
nohup uv run python scripts/train.py > train.log 2>&1 & TRAIN_PID=$! && echo "PID: $TRAIN_PID"
# Still broken — $! is empty when assigned on the same compound line.

# ✅ CORRECT: TRAIN_PID=$! MUST be on its own line, immediately after &
nohup uv run python scripts/train.py > train.log 2>&1 &
TRAIN_PID=$!
echo "Training PID: $TRAIN_PID"
# Now $! is the PID of the background job just started. Do NOT put anything
# between the & and the TRAIN_PID=$! assignment.
```

```bash
# ❌ BAD: Interactive pipe — kills the process when head closes
uv run python scripts/train.py 2>&1 | head -80

# ✅ GOOD: Redirect to log, then tail separately
nohup uv run python scripts/train.py > train.log 2>&1 &
TRAIN_PID=$!
sleep 5 && tail -f train.log &   # follow log in background
```

```bash
# ❌ BAD: Assume fast exit = failure; relaunch
# (script exits fast, agent panics and runs it again — now two processes)
uv run python scripts/train.py > train.log 2>&1 &
# 3 seconds later...
uv run python scripts/train.py > train.log 2>&1 &  # ← DUPLICATE

# ✅ GOOD: Check artifacts before concluding anything
ls -la artifacts/oof.npy 2>/dev/null && echo "artifact exists" || echo "missing"
```

---

## Workflow 1 — Pre-flight check (run before every training launch)

```bash
# Step 1: Is training already running?
RUNNING_PIDS=$(pgrep -f "python scripts/train.py" 2>/dev/null)
if [ -n "$RUNNING_PIDS" ]; then
    echo "⚠️  Training already running — PID(s): $RUNNING_PIDS"
    echo "Wait for it to finish or kill it before relaunching."
    ps -p $RUNNING_PIDS -o pid,stat,etime,cmd --no-headers
    exit 0   # stop here — do NOT launch a second process
fi

# Step 2: Check CPU load — if already saturated, launching more work will
# stall stdout buffering and make new logs appear empty for a long time
echo "=== CPU / load ==="
uptime                          # load average (1m / 5m / 15m)
nproc                           # total logical CPUs
top -bn1 | head -5              # snapshot: %us, %sy, load
# Rule of thumb: if load average (1m) > nproc × 1.5, the machine is overloaded.
# Wait for existing work to finish before launching a new training run.

# Step 3: Are artifacts already fresh?
if [ -f "artifacts/oof.npy" ]; then
    ARTIFACT_AGE=$(( $(date +%s) - $(stat -c %Y artifacts/oof.npy) ))
    echo "oof.npy last modified ${ARTIFACT_AGE}s ago"
    if [ $ARTIFACT_AGE -lt 300 ]; then   # < 5 minutes old
        echo "✅ Artifacts are fresh — skip retraining unless code changed."
        # Only proceed if you made a deliberate code change this iteration
    fi
fi

# Step 4: Check that the script exists and is importable
uv run python -c "import importlib.util, sys; spec=importlib.util.spec_from_file_location('train','scripts/train.py'); assert spec, 'train.py not found'" && echo "Script OK"
```

---

## Workflow 2 — Launch and wait (correct pattern)

```bash
# Launch
nohup uv run python scripts/train.py > train.log 2>&1 &
TRAIN_PID=$!
echo "Training started — PID: $TRAIN_PID"

# Verify it actually started (give it 3 seconds)
sleep 3
if ! kill -0 $TRAIN_PID 2>/dev/null; then
    echo "❌ Process exited immediately — check train.log for errors:"
    tail -30 train.log
    # DO NOT relaunch yet — read the error first
fi

# Monitor (check every 60 seconds)
while kill -0 $TRAIN_PID 2>/dev/null; do
    echo "Still running (PID $TRAIN_PID)..."
    tail -5 train.log
    sleep 60
done

echo "✅ Training finished"
tail -50 train.log
```

---

## Workflow 3 — Kill and relaunch after a code change

Only do this when you have deliberately modified `scripts/train.py`, `src/models.py`,
or `src/features.py` and need to rerun with the new code.

```bash
# Kill any running training process
RUNNING_PIDS=$(pgrep -f "python scripts/train.py" 2>/dev/null)
if [ -n "$RUNNING_PIDS" ]; then
    echo "Killing PID(s): $RUNNING_PIDS"
    kill $RUNNING_PIDS
    sleep 2
    # Confirm dead
    pgrep -f "python scripts/train.py" && echo "still running!" || echo "killed"
fi

# Clean stale log so you don't confuse old output with new
mv train.log train.log.bak 2>/dev/null || true

# Now launch fresh (see Workflow 2)
nohup uv run python scripts/train.py > train.log 2>&1 &
TRAIN_PID=$!
echo "Relaunched — PID: $TRAIN_PID"
```

---

## Workflow 4 — Diagnose a fast / silent exit

If training finished in under 10 seconds and you did not expect that:

```bash
# 1. Check CPU load — saturated CPUs delay stdout flushing; the process may
#    still be running even if the log looks empty or has only 1 line
uptime && nproc
# If load (1m) > nproc: wait, then tail the log again before concluding anything

# 2. Check whether the process is actually still alive
pgrep -a -f "python scripts/train.py" || echo "no training process found"

# 3. Check artifact timestamps — was this a legitimate cached result?
ls -la artifacts/oof.npy artifacts/oof_classes.npy 2>/dev/null

# 4. Check the last OOF line — did a previous run already succeed?
grep -E "FINAL OOF|OOF " train.log | tail -5

# 5. Read the full log for tracebacks
cat train.log

# 6. Only if steps 1-5 give no answer: run with visible output to a new log
#    (NOT interactive pipe — that kills the process when the pipe closes)
nohup uv run python scripts/train.py > train_debug.log 2>&1 &
DBG_PID=$!
while kill -0 $DBG_PID 2>/dev/null; do sleep 5; done
cat train_debug.log
```

Only **after** understanding the cause should you modify any code and relaunch.

---

## See Also

| File | Why |
|------|-----|
| [coding-rules.md](./coding-rules.md) | Code quality gate to apply before launching any training script |
| [experiment-tracking.md](./experiment-tracking.md) | Log OOF scores and artifact timestamps produced by the training process |
