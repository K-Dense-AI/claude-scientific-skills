#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

import importlib
import shutil


def check_python(pkg: str) -> bool:
    try:
        importlib.import_module(pkg)
        return True
    except Exception:
        return False


def check_cmd(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def main():
    py_pkgs = ["pandas", "numpy"]
    cli_cmds = ["Rscript"]
    optional_r = ["agricolae", "lme4"]

    print("Breeding Trial Management System Check")
    print("=" * 40)
    for pkg in py_pkgs:
        print(f"python:{pkg}: {'OK' if check_python(pkg) else 'MISSING'}")
    for cmd in cli_cmds:
        print(f"command:{cmd}: {'OK' if check_cmd(cmd) else 'MISSING'}")

    print("Optional R packages:")
    for pkg in optional_r:
        print(f"r:{pkg}: OPTIONAL")

    print("Breedbase connectivity: OPTIONAL (configured in later phases)")


if __name__ == "__main__":
    main()
