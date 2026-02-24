#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import argparse
import json
import time
import matplotlib.pyplot as plt


def run_check():
    result = {
        "cuda_available": False,
        "gpu_name": "none",
        "total_vram_gb": 0.0,
        "small_matmul_seconds": None,
        "status": "cpu_fallback",
    }

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
    except Exception as exc:
        result["status"] = f"check_failed: {exc.__class__.__name__}"

    result["recommended_for_12gb"] = bool(
        result["cuda_available"] and result["total_vram_gb"] >= 11.5
    )
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
    labels = ["CUDA", "VRAM (12GB)", "Ready"]
    colors = ["#1f77b4", "#2ca02c", "#9467bd"]

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(9.2, 4.8), gridspec_kw={"width_ratios": [1.65, 1]}
    )

    ax1.bar(labels, [1.0, 1.0, 1.0], color="#e6e6e6", edgecolor="#9a9a9a")
    ax1.bar(labels, vals, color=colors, edgecolor="#3a3a3a")
    ax1.set_ylim(0, 1.1)
    ax1.set_ylabel("Normalized score")
    ax1.set_title("GPU/HPC Readiness Metrics")
    ax1.grid(axis="y", alpha=0.2)
    for idx, val in enumerate(vals):
        ax1.text(idx, val + 0.03, f"{val:.2f}", ha="center", fontsize=9)

    ax2.axis("off")
    lines = [
        f"Status: {res['status']}",
        f"GPU: {res['gpu_name']}",
        f"VRAM (GB): {res['total_vram_gb']}",
    ]
    if res["small_matmul_seconds"] is not None:
        lines.append(f"4K matmul (s): {res['small_matmul_seconds']}")
    else:
        lines.append("4K matmul (s): n/a")

    if res["cuda_available"]:
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
