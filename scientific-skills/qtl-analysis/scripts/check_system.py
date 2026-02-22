#!/usr/bin/env python3
"""
QTL Analysis System Requirements Check

Mandatory preflight checker. Run before first use to verify system can run QTL analyses.
Checks: Python version, required packages (tensorQTL, pyQTL), CLI tools (PLINK, GEMMA, R, R/qtl2)
"""

import sys
import subprocess
import shutil
from pathlib import Path


def check_python_version():
    """Check Python >= 3.9"""
    version = sys.version_info
    if version.major < 3 or (version.major == 3 and version.minor < 9):
        return False, f"Python {version.major}.{version.minor} (need >= 3.9)"
    return True, f"Python {version.major}.{version.minor}.{version.micro}"


def check_python_package(package_name):
    """Check if a Python package is installed"""
    try:
        __import__(package_name)
        return True, f"{package_name} installed"
    except ImportError:
        return False, f"{package_name} NOT installed"


def check_cli_tool(command, version_arg="--version"):
    """Check if a CLI tool is available"""
    try:
        result = subprocess.run(
            [command, version_arg], capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0 or result.stdout or result.stderr:
            # Some tools return non-zero but still print version
            version = (result.stdout or result.stderr).strip().split("\n")[0]
            return True, version[:60]
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass
    return False, f"{command} NOT installed"


def check_r_package(package_name):
    """Check if an R package is installed"""
    try:
        result = subprocess.run(
            ["Rscript", "-e", f"library({package_name})"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode == 0:
            return True, f"R/{package_name} installed"
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass
    return False, f"R/{package_name} NOT installed"


def check_gpu():
    """Check if GPU is available for tensorQTL"""
    try:
        import torch

        if torch.cuda.is_available():
            gpu_name = torch.cuda.get_device_name(0)
            vram = torch.cuda.get_device_properties(0).total_memory / 1e9
            return True, f"GPU: {gpu_name} ({vram:.1f} GB VRAM)"
        return False, "GPU: Not available (CPU-only mode)"
    except ImportError:
        return False, "GPU check failed (PyTorch not installed)"


def check_disk_space():
    """Check available disk space"""
    stat = shutil.disk_usage(Path.home())
    free_gb = stat.free / (1024**3)
    if free_gb < 5:
        return False, f"Disk: {free_gb:.1f} GB free (need >= 5 GB)"
    return True, f"Disk: {free_gb:.1f} GB free"


def main():
    print("=" * 60)
    print("QTL Analysis System Requirements Check")
    print("=" * 60)
    print()

    checks = []

    # Python version
    ok, msg = check_python_version()
    checks.append(("Python", ok, msg))

    # Required Python packages
    packages = ["pandas", "numpy", "matplotlib", "scipy", "torch"]
    for pkg in packages:
        ok, msg = check_python_package(pkg)
        checks.append((f"Py/{pkg}", ok, msg))

    # Optional packages (tensorQTL needs GPU, but can be installed)
    optional_packages = ["tensorqtl", "qtl"]
    for pkg in optional_packages:
        ok, msg = check_python_package(pkg)
        checks.append((f"Py/{pkg}", ok, msg, True))  # True = optional

    # CLI tools
    cli_tools = [
        ("plink", "--version"),
        ("plink2", "--version"),
        ("gemma", "-v"),
        ("Rscript", "--version"),
    ]
    for tool, arg in cli_tools:
        ok, msg = check_cli_tool(tool, arg)
        checks.append((tool, ok, msg))

    # R packages
    ok, msg = check_r_package("qtl2")
    checks.append(("R/qtl2", ok, msg))

    # GPU
    ok, msg = check_gpu()
    checks.append(("GPU", ok, msg))

    # Disk
    ok, msg = check_disk_space()
    checks.append(("Disk", ok, msg))

    # Print results
    print()
    required_failures = 0
    optional_failures = 0

    for item in checks:
        if len(item) == 4:  # Optional
            name, ok, msg, optional = item
            if ok:
                print(f"[{name:15s}] ✅ {msg}")
            else:
                print(f"[{name:15s}] ⚠️  {msg} (optional)")
                optional_failures += 1
        else:
            name, ok, msg = item
            if ok:
                print(f"[{name:15s}] ✅ {msg}")
            else:
                print(f"[{name:15s}] ❌ {msg}")
                required_failures += 1

    print()
    print("=" * 60)

    # Verdict
    if required_failures == 0:
        print("✅ System is ready for QTL analysis")
        print()
        print("Capabilities:")
        if any(c[0] == "tensorqtl" and c[1] for c in checks if len(c) == 4):
            print("  - eQTL mapping (tensorQTL, GPU-accelerated)")
        else:
            print("  - eQTL mapping (install tensorQTL for GPU speed)")
        print("  - GWAS (GEMMA LMM, PLINK)")
        if any(c[0] == "R/qtl2" and c[1] for c in checks):
            print("  - Classical QTL (R/qtl2)")
        else:
            print("  - Classical QTL (install R/qtl2)")
        print("  - Population structure (PCA, kinship)")
        return 0
    else:
        print(f"❌ System NOT ready ({required_failures} required components missing)")
        print()
        print("Install missing components:")
        print("  bash scripts/install_deps.sh")
        return 1


if __name__ == "__main__":
    sys.exit(main())
