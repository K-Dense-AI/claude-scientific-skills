#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import argparse
import json
import subprocess
import time
import shutil
import matplotlib.pyplot as plt


def probe_conda_env(env_name: str):
    code = (
        "import json,time;"
        "import torch;"
        "r={'cuda_available':bool(torch.cuda.is_available()),"
        "'gpu_name':'none','total_vram_gb':0.0,'small_matmul_seconds':None};"
        "\n"
        "if r['cuda_available']:\n"
        "  d=torch.device('cuda:0');p=torch.cuda.get_device_properties(d);"
        "  r['gpu_name']=p.name;r['total_vram_gb']=round(p.total_memory/(1024**3),2);"
        "  a=torch.randn((2048,2048),device=d);b=torch.randn((2048,2048),device=d);"
        "  torch.cuda.synchronize();t=time.time();_=a@b;torch.cuda.synchronize();"
        "  r['small_matmul_seconds']=round(time.time()-t,4);"
        "print(json.dumps(r))"
    )
    conda_bin = shutil.which("conda")
    if conda_bin is None:
        fallback = Path.home() / "miniconda3" / "bin" / "conda"
        if fallback.exists():
            conda_bin = str(fallback)
    if conda_bin is None:
        return None

    try:
        res = subprocess.run(
            [conda_bin, "run", "-n", env_name, "python", "-c", code],
            capture_output=True,
            text=True,
            check=True,
        )
        return json.loads(res.stdout.strip().splitlines()[-1])
    except Exception:
        return None


def run_check():
    result = {
        "cuda_available": False,
        "gpu_name": "none",
        "total_vram_gb": 0.0,
        "small_matmul_seconds": None,
        "status": "cpu_fallback",
        "detected_system_gpu": None,
        "detected_system_vram_gb": None,
        "config_hint": None,
        "active_backend": "current_env",
    }

    try:
        smi = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=name,memory.total",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        first = smi.stdout.strip().splitlines()[0]
        name, mem_mb = [x.strip() for x in first.split(",")]
        result["detected_system_gpu"] = name
        result["detected_system_vram_gb"] = round(float(mem_mb) / 1024.0, 2)
    except Exception:
        pass

    try:
        import torch

        result["cuda_available"] = bool(torch.cuda.is_available())
        if result["cuda_available"]:
            dev = torch.device("cuda:0")
            props = torch.cuda.get_device_properties(dev)
            result["gpu_name"] = props.name
            result["total_vram_gb"] = round(props.total_memory / (1024**3), 2)

            a = torch.randn((4096, 4096), device=dev)
            b = torch.randn((4096, 4096), device=dev)
            torch.cuda.synchronize()
            t0 = time.time()
            _ = a @ b
            torch.cuda.synchronize()
            result["small_matmul_seconds"] = round(time.time() - t0, 4)
            result["status"] = "gpu_ready"
        else:
            if result["detected_system_gpu"] and (torch.version.cuda is None):
                result["status"] = "misconfigured_torch_cpu_build"
                result["config_hint"] = (
                    "GPU detected by nvidia-smi, but this Python env has CPU-only PyTorch. "
                    "Install CUDA-enabled torch in this environment."
                )
    except Exception as exc:
        result["status"] = f"check_failed: {exc.__class__.__name__}"

    result["recommended_for_12gb"] = bool(
        result["cuda_available"] and result["total_vram_gb"] >= 11.5
    )

    if (not result["cuda_available"]) and result.get("detected_system_gpu"):
        for env_name in ["cs-gpu", "gpu", "pytorch-gpu"]:
            alt = probe_conda_env(env_name)
            if alt and alt.get("cuda_available"):
                result["cuda_available"] = True
                result["gpu_name"] = alt.get("gpu_name", result["gpu_name"])
                result["total_vram_gb"] = float(
                    alt.get("total_vram_gb", result["total_vram_gb"])
                )
                result["small_matmul_seconds"] = alt.get("small_matmul_seconds")
                result["status"] = "gpu_ready_via_conda_env"
                result["active_backend"] = f"conda:{env_name}"
                result["config_hint"] = (
                    f"Detected working CUDA backend in conda env '{env_name}'. "
                    f"Use `conda run -n {env_name} python ...` or activate that env."
                )
                result["recommended_for_12gb"] = bool(result["total_vram_gb"] >= 11.5)
                break

    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", action="store_true", help="print JSON output")
    args = parser.parse_args()

    res = run_check()
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    (out / "gpu_hpc_check.json").write_text(
        json.dumps(res, indent=2) + "\n", encoding="utf-8"
    )

    vals = [
        1.0 if res["cuda_available"] else 0.0,
        min(float(res["total_vram_gb"]), 12.0) / 12.0,
        1.0 if res["recommended_for_12gb"] else 0.0,
    ]
    labels = ["CUDA", "VRAM12", "READY"]
    colors = ["#1f77b4", "#2ca02c", "#9467bd"]

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(9.2, 4.8), gridspec_kw={"width_ratios": [1.65, 1]}
    )

    ax1.bar(labels, [1.0, 1.0, 1.0], color="#e6e6e6", edgecolor="#9a9a9a")
    ax1.bar(labels, vals, color=colors, edgecolor="#3a3a3a")
    ax1.set_ylim(0, 1.1)
    ax1.set_ylabel("Normalized score")
    ax1.set_title("GPU/HPC Readiness Metrics")
    ax1.tick_params(axis="x", labelrotation=0)
    ax1.grid(axis="y", alpha=0.2)
    for idx, val in enumerate(vals):
        ax1.text(idx, val + 0.03, f"{val:.2f}", ha="center", fontsize=9)

    ax2.axis("off")
    lines = [
        f"Status: {res['status']}",
        f"Backend: {res.get('active_backend', 'current_env')}",
        f"GPU (torch): {res['gpu_name']}",
        f"VRAM (GB): {res['total_vram_gb']}",
    ]
    if res.get("detected_system_gpu"):
        lines.append(
            f"GPU (system): {res['detected_system_gpu']} ({res['detected_system_vram_gb']} GB)"
        )
    if res["small_matmul_seconds"] is not None:
        lines.append(f"4K matmul (s): {res['small_matmul_seconds']}")
    else:
        lines.append("4K matmul (s): n/a")

    if res.get("config_hint"):
        lines.append("Hint: GPU exists but active env torch is CPU-only")
        lines.append("Fix: use CUDA-enabled conda env or reinstall torch")
    elif res["cuda_available"]:
        lines.append("Next: run torch/cuDF benchmark profile")
    else:
        lines.append("Next: install CUDA-enabled PyTorch and rerun")

    ax2.text(
        0.02,
        0.95,
        "\n".join(lines),
        va="top",
        fontsize=9.5,
        bbox={
            "facecolor": "#f8f8f8",
            "edgecolor": "#cccccc",
            "boxstyle": "round,pad=0.4",
        },
    )

    fig.suptitle("GPU/HPC Readiness Check", fontsize=14)
    fig.tight_layout()
    fig.savefig(out / "gpu_hpc_check.png", dpi=160)
    plt.close(fig)

    if args.json:
        print(json.dumps(res))
    else:
        print(f"Status: {res['status']}")
        print(f"CUDA available: {res['cuda_available']}")
        print(f"GPU: {res['gpu_name']}")
        print(f"Total VRAM (GB): {res['total_vram_gb']}")
        print(f"4096x4096 matmul seconds: {res['small_matmul_seconds']}")
        print(f"12GB-ready recommendation: {res['recommended_for_12gb']}")


if __name__ == "__main__":
    main()
