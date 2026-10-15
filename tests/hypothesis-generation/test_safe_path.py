"""Subprocess regressions for standalone scripts under Python safe-path mode."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "hypothesis-generation"
CLIS = sorted(path.name for path in (SKILL_ROOT / "scripts").glob("*.py")
              if not path.name.startswith("_"))


@pytest.fixture
def installed(tmp_path):
    skill = tmp_path / "installed skill with spaces"
    shutil.copytree(SKILL_ROOT, skill,
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    cwd = tmp_path / "unrelated working directory"
    cwd.mkdir()
    (cwd / "_common.py").write_text('raise RuntimeError("wrong helper")\n')
    return skill, cwd


def launch(installed, script, arguments, mode):
    skill, cwd = installed
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONSAFEPATH", None)
    flags = ["-B"]
    if mode == "flag":
        flags.append("-P")
    elif mode == "environment":
        env["PYTHONSAFEPATH"] = "1"
        # An explicitly supplied conflicting helper must not win either.
        env["PYTHONPATH"] = str(cwd)
    return subprocess.run(
        [sys.executable, *flags, str(skill / "scripts" / script), *arguments],
        cwd=cwd, env=env, capture_output=True, text=True, timeout=15,
    )


@pytest.mark.parametrize("mode", ["flag", "environment"])
@pytest.mark.parametrize("script", CLIS)
def test_every_cli_answers_help_under_safe_path(installed, script, mode):
    result = launch(installed, script, ["--help"], mode)
    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout.lower()
    assert "Traceback" not in result.stderr
    assert not list(installed[0].rglob("*.pyc"))


@pytest.mark.parametrize("mode", ["flag", "environment"])
def test_safe_path_preserves_result_and_exit_status(installed, mode):
    arguments = [str(installed[0] / "assets" / "hypothesis_record_template.json")]
    baseline = launch(installed, 'validate_hypothesis_schema.py', arguments, "normal")
    result = launch(installed, 'validate_hypothesis_schema.py', arguments, mode)
    assert result.returncode == baseline.returncode == 0
    assert "Traceback" not in result.stderr
    assert json.loads(result.stdout) == json.loads(baseline.stdout)
    assert json.loads(result.stdout)["valid"] is True
