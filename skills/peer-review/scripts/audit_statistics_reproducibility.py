#!/usr/bin/env python3
"""Audit a structured statistics and reproducibility checklist locally."""

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
    require_text_list,
    require_unique,
    write_json_report,
)

CORE_ITEM_IDS = (
    "question.estimand_alignment",
    "design.unit_and_independence",
    "design.sample_size_precision",
    "design.allocation_randomization",
    "design.blinding",
    "data.inclusion_exclusion",
    "data.missing_data",
    "data.outliers_transformations",
    "analysis.prespecification",
    "analysis.method_design_alignment",
    "analysis.assumptions_diagnostics",
    "analysis.multiplicity",
    "analysis.clustering_repeated_measures",
    "results.effect_sizes_uncertainty",
    "results.denominators_flow",
    "results.complete_outcomes_harms",
    "reproducibility.data_materials_access",
    "reproducibility.code_environment_parameters",
    "reproducibility.provenance_versions",
    "ethics.approval_consent_governance",
    "interpretation.claim_evidence_causality",
    "integrity.deviations_selective_reporting",
)
CATEGORIES = {
    "question",
    "design",
    "data",
    "analysis",
    "results",
    "reproducibility",
    "ethics",
    "interpretation",
    "integrity",
}
APPLICABILITY = {"applicable", "not_applicable"}
STATUSES = {
    "verified_present",
    "partly_documented",
    "missing",
    "not_assessed",
    "not_applicable",
}
SPECIALIST_TRIGGER_IDS = {
    "question.estimand_alignment",
    "data.missing_data",
    "analysis.method_design_alignment",
    "analysis.assumptions_diagnostics",
    "analysis.multiplicity",
    "analysis.clustering_repeated_measures",
    "results.effect_sizes_uncertainty",
    "interpretation.claim_evidence_causality",
}


def load_checklist(payload: Any) -> dict[str, Any]:
    root = require_object(payload, "checklist")
    require_exact_keys(
        root,
        required={
            "schema_version",
            "checklist_id",
            "study_design",
            "specialist_review",
            "items",
        },
        context="checklist",
    )
    schema_version = require_enum(
        root["schema_version"], {"2.0"}, "checklist.schema_version"
    )
    checklist_id = require_identifier(root["checklist_id"], "checklist.checklist_id")
    study_design = require_identifier(root["study_design"], "checklist.study_design")

    specialist = require_object(
        root["specialist_review"], "checklist.specialist_review"
    )
    require_exact_keys(
        specialist,
        required={"needed", "areas", "requested"},
        context="checklist.specialist_review",
    )
    specialist_review = {
        "needed": require_enum(
            specialist["needed"],
            {"yes", "no", "undetermined"},
            "checklist.specialist_review.needed",
        ),
        "areas": require_identifier_list(
            specialist["areas"], "checklist.specialist_review.areas"
        ),
        "requested": require_bool(
            specialist["requested"], "checklist.specialist_review.requested"
        ),
    }

    raw_items = require_list(root["items"], "checklist.items", minimum=1, maximum=200)
    items: list[dict[str, Any]] = []
    item_ids: list[str] = []
    for index, raw_item in enumerate(raw_items):
        context = f"checklist.items[{index}]"
        item = require_object(raw_item, context)
        require_exact_keys(
            item,
            required={
                "id",
                "category",
                "applicability",
                "status",
                "evidence_locations",
                "note",
                "requested_action",
            },
            context=context,
        )
        item_id = require_identifier(item["id"], f"{context}.id")
        item_ids.append(item_id)
        items.append(
            {
                "id": item_id,
                "category": require_enum(
                    item["category"], CATEGORIES, f"{context}.category"
                ),
                "applicability": require_enum(
                    item["applicability"], APPLICABILITY, f"{context}.applicability"
                ),
                "status": require_enum(
                    item["status"], STATUSES, f"{context}.status"
                ),
                "evidence_locations": require_text_list(
                    item["evidence_locations"],
                    f"{context}.evidence_locations",
                    maximum=50,
                ),
                "note": require_text(
                    item["note"],
                    f"{context}.note",
                    allow_empty=True,
                    maximum=2_000,
                ),
                "requested_action": require_text(
                    item["requested_action"],
                    f"{context}.requested_action",
                    allow_empty=True,
                    maximum=2_000,
                ),
            }
        )
    require_unique(item_ids, "checklist.items")
    return {
        "schema_version": schema_version,
        "checklist_id": checklist_id,
        "study_design": study_design,
        "specialist_review": specialist_review,
        "items": items,
    }


def audit(checklist: dict[str, Any]) -> dict[str, Any]:
    errors: list[dict[str, str]] = []
    warnings: list[dict[str, str]] = []
    items_by_id = {item["id"]: item for item in checklist["items"]}
    missing_core = sorted(set(CORE_ITEM_IDS) - set(items_by_id))
    if missing_core:
        errors.extend(issue("CORE_ITEM_MISSING", item_id) for item_id in missing_core)

    status_counts = Counter()
    gaps_by_category: dict[str, list[str]] = {}
    action_missing: list[str] = []
    specialist_triggers: list[str] = []

    for item in checklist["items"]:
        item_id = item["id"]
        status = item["status"]
        status_counts[status] += 1
        if item["applicability"] == "not_applicable":
            if status != "not_applicable":
                errors.append(issue("APPLICABILITY_STATUS_MISMATCH", item_id))
            if not item["note"]:
                errors.append(issue("NOT_APPLICABLE_RATIONALE_REQUIRED", item_id))
            continue
        if status == "not_applicable":
            errors.append(issue("APPLICABILITY_STATUS_MISMATCH", item_id))
            continue
        if status == "verified_present" and not item["evidence_locations"]:
            errors.append(issue("EVIDENCE_LOCATION_REQUIRED", item_id))
        if status == "partly_documented":
            if not item["evidence_locations"]:
                errors.append(issue("PARTIAL_ITEM_NEEDS_EVIDENCE_LOCATION", item_id))
            if not item["requested_action"]:
                errors.append(issue("PARTIAL_ITEM_NEEDS_REQUESTED_ACTION", item_id))
        if status in {"partly_documented", "missing", "not_assessed"}:
            gaps_by_category.setdefault(item["category"], []).append(item_id)
            if not item["requested_action"]:
                action_missing.append(item_id)
                warnings.append(issue("REQUESTED_ACTION_MISSING", item_id))
            if item_id in SPECIALIST_TRIGGER_IDS:
                specialist_triggers.append(item_id)

    specialist = checklist["specialist_review"]
    specialist_recommended = bool(specialist_triggers)
    if specialist["needed"] == "yes" and not specialist["requested"]:
        warnings.append(
            issue("SPECIALIST_REVIEW_NOT_REQUESTED", "specialist_review.requested")
        )
    if specialist["needed"] == "undetermined":
        warnings.append(
            issue("SPECIALIST_REVIEW_UNDETERMINED", "specialist_review.needed")
        )
    if specialist_recommended and specialist["needed"] == "no":
        warnings.append(
            issue("SPECIALIST_REVIEW_DECLARATION_RECHECK", "specialist_review.needed")
        )

    return {
        "schema_version": "2.0",
        "checklist_id": checklist["checklist_id"],
        "valid": not errors,
        "status": (
            "INVALID_CHECKLIST"
            if errors
            else "VALID_WITH_REVIEW_GAPS"
            if gaps_by_category
            else "VALID_NO_RECORDED_GAPS"
        ),
        "errors": errors,
        "warnings": warnings,
        "study_design": checklist["study_design"],
        "item_count": len(checklist["items"]),
        "status_counts": dict(sorted(status_counts.items())),
        "gap_item_ids_by_category": {
            category: sorted(item_ids)
            for category, item_ids in sorted(gaps_by_category.items())
        },
        "item_ids_missing_requested_action": sorted(action_missing),
        "specialist_review": {
            "declared_needed": specialist["needed"],
            "declared_areas": specialist["areas"],
            "declared_requested": specialist["requested"],
            "trigger_item_ids": sorted(specialist_triggers),
            "recheck_recommended": specialist_recommended,
        },
        "notice": (
            "This is a structured completeness and consistency audit, not a "
            "statistical reanalysis, reproducibility claim, quality score, or "
            "publication recommendation. A qualified specialist must evaluate "
            "methods outside the reviewer's competence."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Audit a local statistics/reproducibility checklist and emit only "
            "item identifiers and counts."
        )
    )
    parser.add_argument("checklist", help="Local checklist JSON")
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
