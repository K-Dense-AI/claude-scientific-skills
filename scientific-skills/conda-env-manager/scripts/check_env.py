#!/usr/bin/env python3
"""
Conda Environment Health Checker

Diagnoses conda environment issues:
  - Lists all environments and verifies they exist on disk
  - Detects conda/pip duplicate packages
  - Checks for known version conflicts
  - Reports Python version per environment
"""

import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


def run_conda(args: List[str], timeout: int = 30) -> Optional[str]:
    """Run a conda command and return stdout."""
    try:
        result = subprocess.run(
            ["conda"] + args,
            capture_output=True, text=True, timeout=timeout
        )
        return result.stdout if result.returncode == 0 else None
    except Exception:
        return None


def list_environments() -> List[Dict[str, Any]]:
    """List all conda environments with existence check."""
    output = run_conda(["env", "list", "--json"])
    if not output:
        return []

    data = json.loads(output)
    envs = []
    for env_path in data.get("envs", []):
        p = Path(env_path)
        python_exe = p / "python.exe" if sys.platform == "win32" else p / "bin" / "python"
        exists = python_exe.exists()

        # Get Python version
        py_version = None
        if exists:
            try:
                result = subprocess.run(
                    [str(python_exe), "--version"],
                    capture_output=True, text=True, timeout=5
                )
                if result.returncode == 0:
                    py_version = result.stdout.strip().replace("Python ", "")
            except Exception:
                pass

        # Determine env name
        name = p.name
        if name == "anaconda3" or name == "miniconda3" or name == "miniforge3":
            name = "base"

        envs.append({
            "name": name,
            "path": str(p),
            "exists": exists,
            "python_version": py_version,
        })

    return envs


def get_packages(env_name: str) -> Dict[str, List[Dict]]:
    """Get packages from both conda and pip for an environment."""
    # Conda packages
    conda_output = run_conda(["list", "-n", env_name, "--json"], timeout=60)
    conda_pkgs = json.loads(conda_output) if conda_output else []

    conda_map = {}
    pip_map = {}

    for pkg in conda_pkgs:
        name = pkg["name"].lower()
        channel = pkg.get("channel", "")
        entry = {
            "name": pkg["name"],
            "version": pkg["version"],
            "channel": channel,
        }
        if channel == "pypi":
            pip_map[name] = entry
        else:
            conda_map[name] = entry

    return {"conda": conda_map, "pip": pip_map}


def find_duplicates(packages: Dict[str, Dict]) -> List[Dict]:
    """Find packages installed by both conda and pip."""
    conda = packages["conda"]
    pip_pkgs = packages["pip"]
    dupes = []

    for name, pip_info in pip_pkgs.items():
        if name in conda:
            conda_info = conda[name]
            dupes.append({
                "package": name,
                "conda_version": conda_info["version"],
                "pip_version": pip_info["version"],
                "conflict": conda_info["version"] != pip_info["version"],
            })

    return dupes


def check_known_conflicts(
    packages: Dict[str, Dict],
    conflicts_file: Optional[str] = None
) -> List[Dict]:
    """Check for known version conflicts."""
    # Load known conflicts
    if conflicts_file and Path(conflicts_file).exists():
        with open(conflicts_file) as f:
            known = yaml.safe_load(f) or []
    else:
        # Built-in known conflicts
        known = [
            {
                "name": "numba-numpy",
                "packages": ["numba", "numpy"],
                "check": "numba_numpy",
                "description": "Numba requires specific NumPy version range",
            },
            {
                "name": "torch-numpy",
                "packages": ["torch", "numpy"],
                "check": "torch_numpy",
                "description": "PyTorch may conflict with latest NumPy",
            },
        ]

    all_pkgs = {**packages["conda"], **packages["pip"]}
    issues = []

    for conflict in known:
        pkgs = conflict["packages"]
        present = {p: all_pkgs.get(p.lower()) for p in pkgs if p.lower() in all_pkgs}

        if len(present) < 2:
            continue  # Not all packages present

        check = conflict.get("check", "")

        if check == "numba_numpy":
            numba_v = present.get("numba", {}).get("version", "")
            numpy_v = present.get("numpy", {}).get("version", "")
            if numpy_v and numba_v:
                np_parts = numpy_v.split(".")
                if len(np_parts) >= 2:
                    np_minor = int(np_parts[1])
                    # numba 0.63.x needs numpy <= 2.3
                    if numba_v.startswith("0.63") and int(np_parts[0]) >= 2 and np_minor > 3:
                        issues.append({
                            "conflict": conflict["name"],
                            "description": conflict["description"],
                            "detail": f"numba {numba_v} requires numpy<=2.3, got {numpy_v}",
                            "fix": "pip install 'numpy<=2.3'",
                        })

        elif check == "torch_numpy":
            numpy_v = present.get("numpy", {}).get("version", "")
            if numpy_v:
                np_parts = numpy_v.split(".")
                if len(np_parts) >= 2 and int(np_parts[0]) >= 2 and int(np_parts[1]) >= 4:
                    issues.append({
                        "conflict": conflict["name"],
                        "description": conflict["description"],
                        "detail": f"numpy {numpy_v} may be too new for PyTorch",
                        "fix": "pip install 'numpy<2.4'",
                    })

    return issues


def diagnose_env(env_name: str, conflicts_file: Optional[str] = None) -> Dict:
    """Full diagnosis of a conda environment."""
    packages = get_packages(env_name)
    duplicates = find_duplicates(packages)
    conflicts = check_known_conflicts(packages, conflicts_file)

    total = len(packages["conda"]) + len(packages["pip"])

    return {
        "env_name": env_name,
        "total_packages": total,
        "conda_packages": len(packages["conda"]),
        "pip_packages": len(packages["pip"]),
        "duplicates": duplicates,
        "conflicts": conflicts,
        "healthy": len(duplicates) == 0 and len(conflicts) == 0,
    }


def main():
    if len(sys.argv) < 2 or sys.argv[1] == "--list":
        # List all environments
        envs = list_environments()
        print(json.dumps({"environments": envs}, indent=2, ensure_ascii=False))
    elif sys.argv[1] == "--check":
        if len(sys.argv) < 3:
            print("Usage: check_env.py --check <env_name> [conflicts_file]")
            sys.exit(1)
        env_name = sys.argv[2]
        conflicts_file = sys.argv[3] if len(sys.argv) > 3 else None
        report = diagnose_env(env_name, conflicts_file)
        print(json.dumps(report, indent=2, ensure_ascii=False))
    else:
        print("Usage:")
        print("  check_env.py --list              List all environments")
        print("  check_env.py --check <env_name>  Diagnose an environment")
        sys.exit(1)


if __name__ == "__main__":
    main()
