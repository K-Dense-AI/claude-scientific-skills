#!/usr/bin/env python3
"""
Biomni Environment Setup and Validation Script

This script helps users set up and validate their Biomni environment,
including checking dependencies, API keys, and data availability.
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple


def check_python_version() -> Tuple[bool, str]:
    """Check if Python version is compatible."""
    version = sys.version_info
    if version.major == 3 and version.minor >= 8:
        return True, f"Python {version.major}.{version.minor}.{version.micro} ✓"
    else:
        return False, f"Python {version.major}.{version.minor} - requires Python 3.8+"


def check_conda_env() -> Tuple[bool, str]:
    """Check if running in biomni conda environment."""
    conda_env = os.environ.get('CONDA_DEFAULT_ENV', None)
    if conda_env == 'biomni_e1':
        return True, f"Conda environment: {conda_env} ✓"
    else:
        return False, f"Not in biomni_e1 environment (current: {conda_env})"


def check_package_installed(package: str) -> bool:
    """Check if a Python package is installed."""
    try:
        __import__(package)
        return True
    except ImportError:
        return False


def check_dependencies() -> Tuple[bool, List[str]]:
    """Check for required and optional dependencies."""
    required = ['biomni']
    optional = ['weasyprint', 'markdown2pdf']

    missing_required = [pkg for pkg in required if not check_package_installed(pkg)]
    missing_optional = [pkg for pkg in optional if not check_package_installed(pkg)]

    messages = []
    success = len(missing_required) == 0

    if missing_required:
        messages.append(f"Missing required packages: {', '.join(missing_required)}")
        messages.append("Install with: pip install biomni --upgrade")
    else:
        messages.append("Required packages: ✓")

    if missing_optional:
        messages.append(f"Missing optional packages: {', '.join(missing_optional)}")
        messages.append("For PDF reports, install: pip install weasyprint")

    return success, messages


def check_api_keys() -> Tuple[bool, Dict[str, bool]]:
    """Check which API keys are configured."""
    api_keys = {
        'ANTHROPIC_API_KEY': os.environ.get('ANTHROPIC_API_KEY'),
        'OPENAI_API_KEY': os.environ.get('OPENAI_API_KEY'),
        'GEMINI_API_KEY': os.environ.get('GEMINI_API_KEY'),
        'GROQ_API_KEY': os.environ.get('GROQ_API_KEY'),
    }

    configured = {key: bool(value) for key, value in api_keys.items()}
    has_any = any(configured.values())

    return has_any, configured


def check_data_directory(data_path: str = './data') -> Tuple[bool, str]:
    """Check if Biomni data directory exists and has content."""
    path = Path(data_path)

    if not path.exists():
        return False, f"Data directory not found at {data_path}"

    # Check if directory has files (data has been downloaded)
    files = list(path.glob('*'))
    if len(files) == 0:
        return False, f"Data directory exists but is empty. Run agent once to download."

    # Rough size check (should be ~11GB)
    total_size = sum(f.stat().st_size for f in path.rglob('*') if f.is_file())
    size_gb = total_size / (1024**3)

    if size_gb < 1:
        return False, f"Data directory exists but seems incomplete ({size_gb:.1f} GB)"

    return True, f"Data directory: {data_path} ({size_gb:.1f} GB) ✓"


def check_disk_space(required_gb: float = 20) -> Tuple[bool, str]:
    """Check if sufficient disk space is available."""
    try:
        import shutil
        stat = shutil.disk_usage('.')
        free_gb = stat.free / (1024**3)

        if free_gb >= required_gb:
            return True, f"Disk space: {free_gb:.1f} GB available ✓"
        else:
            return False, f"Low disk space: {free_gb:.1f} GB (need {required_gb} GB)"
    except Exception as e:
        return False, f"Could not check disk space: {e}"


def test_biomni_import() -> Tuple[bool, str]:
    """Test if Biomni can be imported and initialized."""
    try:
        from biomni.agent import A1
        from biomni.config import default_config
        return True, "Biomni import successful ✓"
    except ImportError as e:
        return False, f"Cannot import Biomni: {e}"
    except Exception as e:
        return False, f"Biomni import error: {e}"


def suggest_fixes(results: Dict[str, Tuple[bool, any]]) -> List[str]:
    """Generate suggestions for fixing issues."""
    suggestions = []

    if not results['python'][0]:
        suggestions.append("➜ Upgrade Python to 3.8 or higher")

    if not results['conda'][0]:
        suggestions.append("➜ Activate biomni environment: conda activate biomni_e1")

    if not results['dependencies'][0]:
        suggestions.append("➜ Install Biomni: pip install biomni --upgrade")

    if not results['api_keys'][0]:
        suggestions.append("➜ Set API key: export ANTHROPIC_API_KEY='your-key'")
        suggestions.append("   Or create .env file with API keys")

    if not results['data'][0]:
        suggestions.append("➜ Data will auto-download on first agent.go() call")

    if not results['disk_space'][0]:
        suggestions.append("➜ Free up disk space (need ~20GB total)")

    return suggestions


def main():
    """Run all environment checks and display results."""
    print("=" * 60)
    print("Biomni Environment Validation")
    print("=" * 60)
    print()

    # Run all checks
    results = {}

    print("Checking Python version...")
    results['python'] = check_python_version()
    print(f"  {results['python'][1]}")
    print()

    print("Checking conda environment...")
    results['conda'] = check_conda_env()
    print(f"  {results['conda'][1]}")
    print()

    print("Checking dependencies...")
    results['dependencies'] = check_dependencies()
    for msg in results['dependencies'][1]:
        print(f"  {msg}")
    print()

    print("Checking API keys...")
    results['api_keys'] = check_api_keys()
    has_keys, key_status = results['api_keys']
    for key, configured in key_status.items():
        status = "✓" if configured else "✗"
        print(f"  {key}: {status}")
    print()

    print("Checking Biomni data directory...")
    results['data'] = check_data_directory()
    print(f"  {results['data'][1]}")
    print()

    print("Checking disk space...")
    results['disk_space'] = check_disk_space()
    print(f"  {results['disk_space'][1]}")
    print()

    print("Testing Biomni import...")
    results['biomni_import'] = test_biomni_import()
    print(f"  {results['biomni_import'][1]}")
    print()

    # Summary
    print("=" * 60)
    all_passed = all(result[0] for result in results.values())

    if all_passed:
        print("✓ All checks passed! Environment is ready.")
        print()
        print("Quick start:")
        print("  from biomni.agent import A1")
        print("  agent = A1(path='./data', llm='claude-sonnet-4-20250514')")
        print("  agent.go('Your biomedical task')")
    else:
        print("⚠ Some checks failed. See suggestions below:")
        print()
        suggestions = suggest_fixes(results)
        for suggestion in suggestions:
            print(suggestion)

    print("=" * 60)

    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
