# Install Escalation Policy (All Skills)

When a skill requires system-level installation that cannot be completed by the agent in the current environment, the skill response must:

1. State the exact reason it cannot be auto-resolved.
2. Provide a copy-paste one-liner command for the user.
3. Provide a quick verification command.

## Required Response Pattern

- Reason: what is missing and why agent cannot finish it (system package, driver, kernel module, privileged install, reboot).
- One-liner: command user can paste in terminal.
- Verify: command proving the fix worked.

## Common One-Liners

### CUDA-enabled PyTorch (conda)

```bash
conda create -y -n cs-gpu python=3.11 && conda activate cs-gpu && conda install -y pytorch pytorch-cuda=12.4 -c pytorch -c nvidia
```

### NVIDIA driver/tools (Ubuntu/Debian, requires sudo)

```bash
sudo apt-get update && sudo apt-get install -y nvidia-driver-550 nvidia-utils-550
```

### Verification

```bash
nvidia-smi && python - <<'PY'
import torch
print(torch.__version__, torch.version.cuda, torch.cuda.is_available())
print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'no-gpu')
PY
```
