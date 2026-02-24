<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

---
name: interactive-dashboards
description: >
  Build minimal scientific dashboards using Streamlit or Shiny patterns, plus optional GPU/HPC checks.
  Use when you need interactive views instead of static figures.
license: Apache-2.0
metadata:
  author: Clayton Young (borealBytes / Superior Byte Works, LLC)
  version: "0.1.0"
  category: visualization
---

# Interactive Dashboards

## Overview

This skill adds lightweight dashboard scaffolds for scientific workflows.

## Included Examples

- `examples/streamlit-demo/`: interactive phenotype explorer
- `examples/shiny-demo/`: Python Shiny-style interactive trait dashboard
- `scripts/verify_gpu_hpc.py`: local GPU/HPC readiness check for NVIDIA systems

## CUDA Setup (One-Liners)

If `nvidia-smi` works but the checker reports CPU-only torch, run one of these:

```bash
conda create -y -n cs-gpu python=3.11 && conda activate cs-gpu && conda install -y pytorch pytorch-cuda=12.4 -c pytorch -c nvidia
```

Or install into the active conda environment:

```bash
conda install -y pytorch pytorch-cuda=12.4 -c pytorch -c nvidia
```

If your system is missing `nvidia-smi`, this is a system-level driver/toolkit issue the agent cannot fix inside a user env. Use (Ubuntu/Debian):

```bash
sudo apt-get update && sudo apt-get install -y nvidia-driver-550 nvidia-utils-550
```

Then reboot and rerun:

```bash
python scripts/verify_gpu_hpc.py --json
```
