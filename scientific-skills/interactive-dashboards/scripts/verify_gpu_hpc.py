#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import argparse
import json
import time


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
