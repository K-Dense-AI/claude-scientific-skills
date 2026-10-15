#!/usr/bin/env python3
"""Audit a local operationalization and measurement checklist without scoring."""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _common import (  # noqa: E402
    ValidationError,
    error_exit,
    issue,
    read_json,
    require_bool,
    require_enum,
    require_exact_keys,
    require_identifier,
    require_identifier_list,
    require_list,
    require_object,
    require_text,
    require_unique,
    write_json_report,
)

ROLES = {
    "intervention",
    "exposure",
    "outcome",
    "mediator",
    "confounder",
    "selection",
    "effect_modifier",
    "negative_control",
    "positive_control",
    "other",
}
PLANNING_STATUSES = {"planned", "complete", "not_applicable", "unresolved"}
REVIEW_STATUSES = {"pending", "complete", "specialist_required"}
BOOLEAN_FIELDS = (
    "construct_defined",
    "operational_definition_recorded",
    "population_scope_recorded",
    "unit_or_categories_recorded",
    "timing_recorded",
    "instrument_or_method_recorded",
    "validity_applicability_reviewed",
    "reliability_or_repeatability_plan_recorded",
    "missingness_plan_recorded",
    "limitations_recorded",
)
STATUS_FIELDS = (
    "calibration_or_quality_control_status",
    "measurement_invariance_or_comparability_status",
    "masking_status",
    "threshold_or_cutpoint_status",
)


def load_checklist(payload: Any) -> dict[str, Any]:
    root = require_object(payload, "checklist")
    require_exact_keys(
        root,
        required={
            "schema_version",
            "checklist_id",
            "record_id",
            "human_reviewer",
            "items",
        },
        context="checklist",
    )
    raw_items = require_list(root["items"], "checklist.items", minimum=1, maximum=200)
    parsed_items: list[dict[str, Any]] = []
    identifiers: list[str] = []
    item_fields = {
        "measurement_id",
        "applicability",
        "variable_role",
        "validity_evidence_source_ids",
        "human_review_status",
        "note",
        *BOOLEAN_FIELDS,
        *STATUS_FIELDS,
    }
    for index, raw_item in enumerate(raw_items):
        context = f"checklist.items[{index}]"
        item = require_object(raw_item, context)
        require_exact_keys(item, required=item_fields, context=context)
        measurement_id = require_identifier(
            item["measurement_id"], f"{context}.measurement_id"
        )
        identifiers.append(measurement_id)
        parsed: dict[str, Any] = {
            "measurement_id": measurement_id,
            "applicability": require_enum(
                item["applicability"],
                {"applicable", "not_applicable"},
                f"{context}.applicability",
            ),
            "variable_role": require_enum(
                item["variable_role"], ROLES, f"{context}.variable_role"
            ),
            "validity_evidence_source_ids": require_identifier_list(
                item["validity_evidence_source_ids"],
                f"{context}.validity_evidence_source_ids",
                maximum=100,
            ),
            "human_review_status": require_enum(
                item["human_review_status"],
                REVIEW_STATUSES,
                f"{context}.human_review_status",
            ),
            "note": require_text(
                item["note"], f"{context}.note", allow_empty=True, maximum=2_000
            ),
        }
        for field in BOOLEAN_FIELDS:
            parsed[field] = require_bool(item[field], f"{context}.{field}")
        for field in STATUS_FIELDS:
            parsed[field] = require_enum(
                item[field], PLANNING_STATUSES, f"{context}.{field}"
            )
        parsed_items.append(parsed)
    require_unique(identifiers, "checklist.items")
    return {
        "schema_version": require_enum(
            root["schema_version"], {"2.0"}, "checklist.schema_version"
        ),
        "checklist_id": require_identifier(
            root["checklist_id"], "checklist.checklist_id"
        ),
        "record_id": require_identifier(root["record_id"], "checklist.record_id"),
        "human_reviewer": require_text(
            root["human_reviewer"], "checklist.human_reviewer", minimum=3
        ),
        "items": parsed_items,
    }


def audit(checklist: dict[str, Any]) -> dict[str, Any]:
    errors: list[dict[str, str]] = []
    warnings: list[dict[str, str]] = []
    gaps_by_measurement: dict[str, list[str]] = {}
    role_counts: Counter[str] = Counter()
    review_counts: Counter[str] = Counter()

    for item in checklist["items"]:
        measurement_id = item["measurement_id"]
        role_counts[item["variable_role"]] += 1
        review_counts[item["human_review_status"]] += 1
        gaps: list[str] = []

        if item["applicability"] == "not_applicable":
            if not item["note"]:
                errors.append(issue("NOT_APPLICABLE_RATIONALE_REQUIRED", measurement_id))
            continue

        for field in BOOLEAN_FIELDS:
            if not item[field]:
                gaps.append(field)
        if not item["validity_evidence_source_ids"]:
            errors.append(issue("VALIDITY_SOURCE_REQUIRED", measurement_id))
        for field in STATUS_FIELDS:
            status = item[field]
            if status in {"planned", "unresolved"}:
                gaps.append(field)
            if status == "not_applicable" and not item["note"]:
                errors.append(
                    issue("NOT_APPLICABLE_STATUS_NEEDS_RATIONALE", f"{measurement_id}:{field}")
                )
        if item["human_review_status"] != "complete":
            gaps.append("human_review_status")
        if not item["note"]:
            warnings.append(issue("MEASUREMENT_NOTE_EMPTY", measurement_id))
        if gaps:
            gaps_by_measurement[measurement_id] = sorted(set(gaps))

    if not role_counts["outcome"]:
        warnings.append(issue("NO_OUTCOME_MEASUREMENT_DECLARED", "items"))
    if not (role_counts["intervention"] or role_counts["exposure"]):
        warnings.append(issue("NO_INTERVENTION_OR_EXPOSURE_DECLARED", "items"))

    return {
        "schema_version": "2.0",
        "checklist_id": checklist["checklist_id"],
        "record_id": checklist["record_id"],
        "valid": not errors,
        "status": (
            "INVALID_CHECKLIST"
            if errors
            else "VALID_WITH_MEASUREMENT_GAPS"
            if gaps_by_measurement
            else "VALID_HUMAN_REVIEW_COMPLETE"
        ),
        "errors": errors,
        "warnings": warnings,
        "measurement_count": len(checklist["items"]),
        "role_counts": dict(sorted(role_counts.items())),
        "human_review_status_counts": dict(sorted(review_counts.items())),
        "gap_fields_by_measurement_id": gaps_by_measurement,
        "notice": (
            "This checklist records declared measurement work. It does not test "
            "an instrument, establish construct validity, certify comparability, "
            "or score measurement quality. Qualified human review is required."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Audit a bounded local JSON operationalization checklist and report "
            "measurement IDs and unresolved fields without scientific scoring."
        )
    )
    parser.add_argument("checklist", help="Local operationalization checklist JSON")
    parser.add_argument("-o", "--output", help="Optional local JSON report")
    parser.add_argument(
        "--force", action="store_true", help="Replace an existing output file"
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        report = audit(load_checklist(read_json(args.checklist)))
        write_json_report(report, args.output, force=args.force)
        return 0 if report["valid"] else 1
    except ValidationError as exc:
        return error_exit(exc)


if __name__ == "__main__":
    raise SystemExit(main())
