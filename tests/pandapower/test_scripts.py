"""Tests for the pandapower power-flow helper.

The helper is a thin, deterministic layer over `pandapower.runpp`: it reads a
case file, solves it, and reports limits. These tests therefore check the three
things that can silently go wrong in that layer, with real solves rather than
mocks:

* the demo case keeps converging to the voltages recorded in `SKILL.md`, so a
  transcription error in the case file or in the bus-by-bus report is caught;
* an overloaded case exits with the violation code instead of quietly reporting
  a healthy grid;
* an N-1 outage that islands a bus is reported as de-energized rather than as
  an undervoltage (pandapower returns NaN for an isolated bus, and a naive
  comparison turns that into a bogus overvoltage).

Skipped when pandapower is not installed; the repo's `--isolated` runner
installs it from `tests/skill-requirements.toml`.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pytest

import skill_contract

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "pandapower"
SCRIPT = SKILL_ROOT / "scripts" / "run_power_flow.py"
DEMO_CASE = SKILL_ROOT / "assets" / "case-6bus-ring.json"
FIXTURES = Path(__file__).resolve().parent / "fixtures"
OVERLOAD_CASE = FIXTURES / "case-overload.json"

pytest.importorskip("pandapower", reason="the helper wraps pandapower.runpp")

CliHelpTests = skill_contract.cli.help_test_case(SKILL_ROOT)


def run_script(case: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), str(case), *arguments],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )


class DemoCaseTests(unittest.TestCase):
    """The documented example must keep producing the documented numbers."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.tmp = tempfile.TemporaryDirectory()
        cls.report_path = Path(cls.tmp.name) / "report.json"
        cls.result = run_script(DEMO_CASE, "--json-out", str(cls.report_path))
        cls.report = json.loads(cls.report_path.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmp.cleanup()

    def test_converges_without_violations(self) -> None:
        self.assertEqual(self.result.returncode, 0, self.result.stdout + self.result.stderr)
        power_flow = self.report["power_flow"]
        self.assertTrue(power_flow["converged"])
        self.assertEqual(power_flow["violations"], [])
        self.assertEqual(power_flow["iterations"], 3)

    def test_bus_and_branch_counts(self) -> None:
        power_flow = self.report["power_flow"]
        self.assertEqual(len(power_flow["buses"]), 5)
        self.assertEqual(
            [branch["element"] for branch in power_flow["branches"]],
            ["line L12", "line L23", "line L34", "line L14", "trafo T45"],
        )

    def test_voltage_profile_matches_skill_documentation(self) -> None:
        voltages = {
            bus["bus"]: bus["vm_pu"] for bus in self.report["power_flow"]["buses"]
        }
        self.assertAlmostEqual(voltages["B1"], 1.02, places=5)
        self.assertAlmostEqual(voltages["B2"], 1.01, places=5)
        self.assertAlmostEqual(voltages["B3"], 1.0076, places=3)
        self.assertAlmostEqual(voltages["B4"], 1.0113, places=3)
        # B5 is the only bus behind the transformer, so its voltage reflects
        # both the transformer impedance and the 10 kV load.
        self.assertAlmostEqual(voltages["B5"], 0.9974, places=3)

    def test_summary_matches_skill_documentation(self) -> None:
        summary = self.report["power_flow"]["summary"]
        self.assertAlmostEqual(summary["vm_min_pu"], 0.99739, places=4)
        self.assertAlmostEqual(summary["max_loading_percent"], 31.799, places=2)
        self.assertAlmostEqual(summary["p_loss_mw"], 0.1769, places=3)
        self.assertEqual(summary["de_energized_buses"], 0)

    def test_versions_are_recorded(self) -> None:
        self.assertIn("pandapower", self.report["versions"])
        self.assertRegex(self.report["versions"]["pandapower"], r"^\d+\.\d+")


class LimitViolationTests(unittest.TestCase):
    def test_overload_exits_with_violation_code(self) -> None:
        result = run_script(OVERLOAD_CASE)
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("overload: line L12", result.stdout)
        self.assertIn("overload: line L23", result.stdout)

    def test_loading_limit_is_configurable(self) -> None:
        # The demo's worst branch sits at 31.8%, so a 25% limit must fail it.
        result = run_script(DEMO_CASE, "--loading-limit", "25")
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("overload: trafo T45", result.stdout)

    def test_voltage_limit_is_configurable(self) -> None:
        # B5 sits at 0.9974 pu, inside the default band but below a 1.00 floor.
        result = run_script(DEMO_CASE, "--vm-min", "1.00")
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("undervoltage: bus B5", result.stdout)


class N1ScreeningTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tmp = tempfile.TemporaryDirectory()
        report_path = Path(cls.tmp.name) / "n1.json"
        cls.result = run_script(DEMO_CASE, "--n1", "--json-out", str(report_path))
        cls.report = json.loads(report_path.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmp.cleanup()

    def test_every_branch_is_screened(self) -> None:
        screenings = self.report["n1"]["screenings"]
        self.assertEqual(
            [item["outage"] for item in screenings],
            ["line L12", "line L23", "line L34", "line L14", "trafo T45"],
        )
        self.assertEqual(self.report["n1"]["non_converged"], 0)

    def test_islanded_bus_is_reported_as_de_energized(self) -> None:
        trafo = next(
            item
            for item in self.report["n1"]["screenings"]
            if item["outage"] == "trafo T45"
        )
        self.assertEqual(trafo["violations"], 1)
        self.assertEqual(trafo["worst"], ["bus B5 de-energized"])
        self.assertEqual(trafo["de_energized_buses"], 1)

    def test_islanding_makes_the_run_fail(self) -> None:
        self.assertEqual(self.result.returncode, 2, self.result.stdout + self.result.stderr)
        kinds = [v["kind"] for v in self.report["power_flow"]["violations"]]
        self.assertEqual(kinds, ["n1_failure"])
        self.assertIn("[FAIL] N-1 trafo T45", self.result.stdout)


class CaseValidationTests(unittest.TestCase):
    def test_unknown_bus_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp) / "case.json"
            case.write_text(
                json.dumps(
                    {
                        "buses": [{"id": "B1", "vn_kv": 110}],
                        "ext_grids": [{"bus": "B1"}],
                        "loads": [{"bus": "B9", "p_mw": 1}],
                    }
                ),
                encoding="utf-8",
            )
            result = run_script(case)
        self.assertEqual(result.returncode, 3, result.stdout + result.stderr)
        self.assertIn("unknown bus 'B9'", result.stdout)

    def test_missing_required_field_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp) / "case.json"
            case.write_text(
                json.dumps({"buses": [{"id": "B1"}]}), encoding="utf-8"
            )
            result = run_script(case)
        self.assertEqual(result.returncode, 3, result.stdout + result.stderr)
        self.assertIn("missing required field 'vn_kv'", result.stdout)

    def test_unreadable_case_is_rejected(self) -> None:
        result = run_script(FIXTURES / "does-not-exist.json")
        self.assertEqual(result.returncode, 3, result.stdout + result.stderr)
        self.assertIn("cannot read case", result.stdout)


if __name__ == "__main__":
    unittest.main()
