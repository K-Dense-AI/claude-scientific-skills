#!/usr/bin/env python3
"""Audit falsifiers, discriminating tests, nulls, and controls without scoring."""

from __future__ import annotations

import argparse
import sys
import re
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _common import (  # noqa: E402
    ValidationError,
    error_exit,
    issue,
    read_json,
    require_enum,
    require_exact_keys,
    require_identifier,
    require_list,
    require_object,
    require_text,
    require_text_list,
    require_unique,
    write_json_report,
)
from validate_hypothesis_schema import load_hypothesis_record, validate_record  # noqa: E402

CONTROL_TYPES = {
    "negative_exposure",
    "negative_outcome",
    "procedural_negative",
    "positive_control",
    "vehicle_or_sham",
    "other",
}
NEGATIVE_TYPES = {"negative_exposure", "negative_outcome", "procedural_negative"}


def _parse_falsifier(raw: Any, context: str) -> dict[str, Any]:
    value = require_object(raw, context)
    require_exact_keys(
        value,
        required={
            "prediction_id",
            "conditions",
            "observable",
            "incompatible_result",
            "assumption_failure_checks",
        },
        context=context,
    )
    return {
        "prediction_id": require_identifier(
            value["prediction_id"], f"{context}.prediction_id"
        ),
        "conditions": require_text(
            value["conditions"], f"{context}.conditions", minimum=10
        ),
        "observable": require_text(
            value["observable"], f"{context}.observable", minimum=5
        ),
        "incompatible_result": require_text(
            value["incompatible_result"],
            f"{context}.incompatible_result",
            minimum=10,
        ),
        "assumption_failure_checks": require_text_list(
            value["assumption_failure_checks"],
            f"{context}.assumption_failure_checks",
            minimum=1,
        ),
    }


def _parse_discriminating_tests(raw: Any, context: str) -> list[dict[str, str]]:
    values = require_list(raw, context, minimum=1, maximum=100)
    parsed: list[dict[str, str]] = []
    identifiers: list[str] = []
    fields = {
        "test_id",
        "rival_hypothesis_id",
        "focal_expected",
        "rival_expected",
        "indeterminate_result",
    }
    for index, raw_value in enumerate(values):
        item_context = f"{context}[{index}]"
        value = require_object(raw_value, item_context)
        require_exact_keys(value, required=fields, context=item_context)
        test_id = require_identifier(value["test_id"], f"{item_context}.test_id")
        identifiers.append(test_id)
        parsed.append(
            {
                "test_id": test_id,
                "rival_hypothesis_id": require_identifier(
                    value["rival_hypothesis_id"],
                    f"{item_context}.rival_hypothesis_id",
                ),
                "focal_expected": require_text(
                    value["focal_expected"],
                    f"{item_context}.focal_expected",
                    minimum=10,
                ),
                "rival_expected": require_text(
                    value["rival_expected"],
                    f"{item_context}.rival_expected",
                    minimum=10,
                ),
                "indeterminate_result": require_text(
                    value["indeterminate_result"],
                    f"{item_context}.indeterminate_result",
                    minimum=10,
                ),
            }
        )
    require_unique(identifiers, context)
    return parsed


def _parse_nulls(raw: Any, context: str) -> list[dict[str, str]]:
    values = require_list(raw, context, minimum=1, maximum=100)
    parsed: list[dict[str, str]] = []
    identifiers: list[str] = []
    fields = {"null_id", "statement", "analysis_id", "interpretation_limit"}
    for index, raw_value in enumerate(values):
        item_context = f"{context}[{index}]"
        value = require_object(raw_value, item_context)
        require_exact_keys(value, required=fields, context=item_context)
        null_id = require_identifier(value["null_id"], f"{item_context}.null_id")
        identifiers.append(null_id)
        parsed.append(
            {
                "null_id": null_id,
                "statement": require_text(
                    value["statement"], f"{item_context}.statement", minimum=10
                ),
                "analysis_id": require_identifier(
                    value["analysis_id"], f"{item_context}.analysis_id"
                ),
                "interpretation_limit": require_text(
                    value["interpretation_limit"],
                    f"{item_context}.interpretation_limit",
                    minimum=10,
                ),
            }
        )
    require_unique(identifiers, context)
    return parsed


def _parse_controls(raw: Any, context: str) -> list[dict[str, str]]:
    values = require_list(raw, context, minimum=1, maximum=100)
    parsed: list[dict[str, str]] = []
    identifiers: list[str] = []
    fields = {
        "control_id",
        "control_type",
        "rationale",
        "expected_result",
        "failure_implication",
    }
    for index, raw_value in enumerate(values):
        item_context = f"{context}[{index}]"
        value = require_object(raw_value, item_context)
        require_exact_keys(value, required=fields, context=item_context)
        control_id = require_identifier(
            value["control_id"], f"{item_context}.control_id"
        )
        identifiers.append(control_id)
        parsed.append(
            {
                "control_id": control_id,
                "control_type": require_enum(
                    value["control_type"],
                    CONTROL_TYPES,
                    f"{item_context}.control_type",
                ),
                "rationale": require_text(
                    value["rationale"], f"{item_context}.rationale", minimum=10
                ),
                "expected_result": require_text(
                    value["expected_result"],
                    f"{item_context}.expected_result",
                    minimum=5,
                ),
                "failure_implication": require_text(
                    value["failure_implication"],
                    f"{item_context}.failure_implication",
                    minimum=10,
                ),
            }
        )
    require_unique(identifiers, context)
    return parsed


def _parse_outcome_interpretation(raw: Any, context: str) -> dict[str, str]:
    value = require_object(raw, context)
    fields = {
        "consistent_with_candidate",
        "challenges_candidate",
        "supports_neither_or_mixed",
    }
    require_exact_keys(value, required=fields, context=context)
    return {
        field: require_text(value[field], f"{context}.{field}", minimum=10)
        for field in fields
    }


def load_checklist(payload: Any) -> dict[str, Any]:
    root = require_object(payload, "checklist")
    require_exact_keys(
        root,
        required={"schema_version", "checklist_id", "record_id", "hypotheses"},
        context="checklist",
    )
    raw_hypotheses = require_list(
        root["hypotheses"], "checklist.hypotheses", minimum=1, maximum=50
    )
    parsed_hypotheses: list[dict[str, Any]] = []
    hypothesis_ids: list[str] = []
    fields = {
        "hypothesis_id",
        "candidate_status",
        "assumptions",
        "boundary_conditions",
        "falsifier",
        "discriminating_tests",
        "nulls",
        "controls",
        "outcome_interpretation",
        "human_review_status",
    }
    for index, raw_hypothesis in enumerate(raw_hypotheses):
        context = f"checklist.hypotheses[{index}]"
        hypothesis = require_object(raw_hypothesis, context)
        require_exact_keys(hypothesis, required=fields, context=context)
        hypothesis_id = require_identifier(
            hypothesis["hypothesis_id"], f"{context}.hypothesis_id"
        )
        hypothesis_ids.append(hypothesis_id)
        parsed_hypotheses.append(
            {
                "hypothesis_id": hypothesis_id,
                "candidate_status": require_enum(
                    hypothesis["candidate_status"],
                    {"candidate"},
                    f"{context}.candidate_status",
                ),
                "assumptions": require_text_list(
                    hypothesis["assumptions"],
                    f"{context}.assumptions",
                    minimum=1,
                ),
                "boundary_conditions": require_text_list(
                    hypothesis["boundary_conditions"],
                    f"{context}.boundary_conditions",
                    minimum=1,
                ),
                "falsifier": _parse_falsifier(
                    hypothesis["falsifier"], f"{context}.falsifier"
                ),
                "discriminating_tests": _parse_discriminating_tests(
                    hypothesis["discriminating_tests"],
                    f"{context}.discriminating_tests",
                ),
                "nulls": _parse_nulls(
                    hypothesis["nulls"], f"{context}.nulls"
                ),
                "controls": _parse_controls(
                    hypothesis["controls"], f"{context}.controls"
                ),
                "outcome_interpretation": _parse_outcome_interpretation(
                    hypothesis["outcome_interpretation"],
                    f"{context}.outcome_interpretation",
                ),
                "human_review_status": require_enum(
                    hypothesis["human_review_status"],
                    {"pending", "complete", "specialist_required"},
                    f"{context}.human_review_status",
                ),
            }
        )
    require_unique(hypothesis_ids, "checklist.hypotheses")
    return {
        "schema_version": require_enum(
            root["schema_version"], {"2.0"}, "checklist.schema_version"
        ),
        "checklist_id": require_identifier(
            root["checklist_id"], "checklist.checklist_id"
        ),
        "record_id": require_identifier(root["record_id"], "checklist.record_id"),
        "hypotheses": parsed_hypotheses,
    }


def _normalize(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.casefold()).strip()


def audit(
    checklist: dict[str, Any], record: dict[str, Any] | None = None
) -> dict[str, Any]:
    errors: list[dict[str, str]] = []
    warnings: list[dict[str, str]] = []
    review_counts: Counter[str] = Counter()
    test_ids: list[str] = []
    null_ids: list[str] = []

    record_ids: dict[str, set[str]] | None = None
    if record is not None:
        record_report = validate_record(record)
        if not record_report["valid"]:
            raise ValidationError(
                "optional hypothesis record must pass schema validation first"
            )
        record_ids = {
            "hypotheses": {
                item["hypothesis_id"] for item in record["hypotheses"]
            },
            "predictions": {
                item["prediction_id"] for item in record["predictions"]
            },
            "analyses": {
                item["analysis_id"] for item in record["analysis_plan"]["analyses"]
            },
            "controls": {
                item["control_id"] for item in record["negative_controls"]
            },
            "nulls": {item["null_id"] for item in record["null_hypotheses"]},
        }

    for hypothesis in checklist["hypotheses"]:
        hypothesis_id = hypothesis["hypothesis_id"]
        review_counts[hypothesis["human_review_status"]] += 1
        if hypothesis["human_review_status"] != "complete":
            warnings.append(issue("HUMAN_REVIEW_INCOMPLETE", hypothesis_id))

        falsifier = hypothesis["falsifier"]
        if record_ids is not None:
            if hypothesis_id not in record_ids["hypotheses"]:
                errors.append(issue("UNKNOWN_HYPOTHESIS_ID", hypothesis_id))
            if falsifier["prediction_id"] not in record_ids["predictions"]:
                errors.append(
                    issue("UNKNOWN_PREDICTION_ID", falsifier["prediction_id"])
                )

        for test in hypothesis["discriminating_tests"]:
            test_ids.append(test["test_id"])
            if test["rival_hypothesis_id"] == hypothesis_id:
                errors.append(issue("FOCAL_HYPOTHESIS_LISTED_AS_RIVAL", test["test_id"]))
            if _normalize(test["focal_expected"]) == _normalize(
                test["rival_expected"]
            ):
                errors.append(
                    issue("FOCAL_AND_RIVAL_EXPECTATIONS_IDENTICAL", test["test_id"])
                )
            if (
                record_ids is not None
                and test["rival_hypothesis_id"] not in record_ids["hypotheses"]
            ):
                errors.append(
                    issue("UNKNOWN_RIVAL_HYPOTHESIS_ID", test["test_id"])
                )

        for null in hypothesis["nulls"]:
            null_ids.append(null["null_id"])
            if record_ids is not None:
                if null["null_id"] not in record_ids["nulls"]:
                    errors.append(issue("UNKNOWN_NULL_ID", null["null_id"]))
                if null["analysis_id"] not in record_ids["analyses"]:
                    errors.append(
                        issue("UNKNOWN_ANALYSIS_ID", null["analysis_id"])
                    )

        negative_controls = [
            control
            for control in hypothesis["controls"]
            if control["control_type"] in NEGATIVE_TYPES
        ]
        if not negative_controls:
            errors.append(issue("NEGATIVE_CONTROL_REQUIRED", hypothesis_id))
        if record_ids is not None:
            for control in negative_controls:
                if control["control_id"] not in record_ids["controls"]:
                    errors.append(
                        issue("UNKNOWN_NEGATIVE_CONTROL_ID", control["control_id"])
                    )

    require_unique(test_ids, "checklist.discriminating_tests")
    require_unique(null_ids, "checklist.nulls")

    return {
        "schema_version": "2.0",
        "checklist_id": checklist["checklist_id"],
        "record_id": checklist["record_id"],
        "valid": not errors,
        "status": (
            "INVALID_CHECKLIST"
            if errors
            else "VALID_PENDING_HUMAN_REVIEW"
            if any(item["human_review_status"] != "complete" for item in checklist["hypotheses"])
            else "VALID_HUMAN_REVIEW_DECLARED_COMPLETE"
        ),
        "errors": errors,
        "warnings": warnings,
        "hypothesis_ids": sorted(
            item["hypothesis_id"] for item in checklist["hypotheses"]
        ),
        "discriminating_test_ids": sorted(test_ids),
        "null_ids": sorted(null_ids),
        "human_review_status_counts": dict(sorted(review_counts.items())),
        "record_cross_check_performed": record is not None,
        "notice": (
            "This audit checks declared falsifiers, rival contrasts, nulls, "
            "controls, and links. It does not prove falsifiability, validate "
            "control assumptions, interpret results, or select a hypothesis."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Audit a bounded local falsification/control JSON checklist without "
            "scientific scoring or candidate selection."
        )
    )
    parser.add_argument("checklist", help="Local falsification/control checklist JSON")
    parser.add_argument(
        "--record", help="Optional local hypothesis record JSON for cross-checks"
    )
    parser.add_argument("-o", "--output", help="Optional local JSON report")
    parser.add_argument(
        "--force", action="store_true", help="Replace an existing output file"
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        record = (
            load_hypothesis_record(read_json(args.record)) if args.record else None
        )
        report = audit(load_checklist(read_json(args.checklist)), record)
        write_json_report(report, args.output, force=args.force)
        return 0 if report["valid"] else 1
    except ValidationError as exc:
        return error_exit(exc)


if __name__ == "__main__":
    raise SystemExit(main())
