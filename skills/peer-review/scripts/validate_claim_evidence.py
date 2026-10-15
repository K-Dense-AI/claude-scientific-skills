#!/usr/bin/env python3
"""Validate a bounded claim-evidence alignment matrix without echoing prose."""

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
    read_csv_records,
    require_enum,
    require_identifier,
    require_text,
    require_unique,
    split_identifiers,
    write_json_report,
)

FIELDS = (
    "claim_id",
    "location",
    "claim_type",
    "claim_summary",
    "evidence_ids",
    "support_level",
    "alignment_issue",
    "limitation",
    "requested_action",
)
CLAIM_TYPES = {
    "primary_outcome",
    "secondary_outcome",
    "causal",
    "mechanistic",
    "diagnostic",
    "prediction",
    "safety",
    "generalization",
    "methods",
    "other",
}
SUPPORT_LEVELS = {
    "supported",
    "partly_supported",
    "unsupported",
    "not_assessed",
}
ALIGNMENT_ISSUES = {
    "none",
    "direction",
    "magnitude",
    "population",
    "outcome",
    "timepoint",
    "causal_language",
    "scope",
    "uncertainty",
    "selective_reporting",
    "other",
}


def load_matrix(raw_path: str) -> list[dict[str, Any]]:
    rows = read_csv_records(
        raw_path,
        required_fields=FIELDS,
        allowed_fields=FIELDS,
    )
    parsed: list[dict[str, Any]] = []
    claim_ids: list[str] = []
    for line_number, row in enumerate(rows, start=2):
        context = f"matrix row {line_number}"
        claim_id = require_identifier(row["claim_id"], f"{context}.claim_id")
        claim_ids.append(claim_id)
        parsed.append(
            {
                "claim_id": claim_id,
                "location": require_text(
                    row["location"], f"{context}.location", maximum=500
                ),
                "claim_type": require_enum(
                    row["claim_type"], CLAIM_TYPES, f"{context}.claim_type"
                ),
                "claim_summary": require_text(
                    row["claim_summary"],
                    f"{context}.claim_summary",
                    minimum=10,
                    maximum=2_000,
                ),
                "evidence_ids": split_identifiers(
                    row["evidence_ids"],
                    f"{context}.evidence_ids",
                    allow_empty=True,
                ),
                "support_level": require_enum(
                    row["support_level"],
                    SUPPORT_LEVELS,
                    f"{context}.support_level",
                ),
                "alignment_issue": require_enum(
                    row["alignment_issue"],
                    ALIGNMENT_ISSUES,
                    f"{context}.alignment_issue",
                ),
                "limitation": require_text(
                    row["limitation"],
                    f"{context}.limitation",
                    allow_empty=True,
                    maximum=2_000,
                ),
                "requested_action": require_text(
                    row["requested_action"],
                    f"{context}.requested_action",
                    allow_empty=True,
                    maximum=2_000,
                ),
            }
        )
    require_unique(claim_ids, "claim-evidence matrix")
    return parsed


def validate_matrix(rows: list[dict[str, Any]]) -> dict[str, Any]:
    errors: list[dict[str, str]] = []
    warnings: list[dict[str, str]] = []
    support_counts = Counter()
    issue_counts = Counter()
    gap_ids: list[str] = []
    action_missing_ids: list[str] = []

    for row in rows:
        claim_id = row["claim_id"]
        support = row["support_level"]
        alignment_issue = row["alignment_issue"]
        evidence_ids = row["evidence_ids"]
        support_counts[support] += 1
        issue_counts[alignment_issue] += 1

        if support == "supported":
            if not evidence_ids:
                errors.append(issue("SUPPORTED_CLAIM_HAS_NO_EVIDENCE", claim_id))
            if alignment_issue != "none":
                errors.append(issue("SUPPORTED_CLAIM_HAS_ALIGNMENT_ISSUE", claim_id))
        elif support == "partly_supported":
            gap_ids.append(claim_id)
            if not evidence_ids:
                errors.append(issue("PARTIAL_CLAIM_HAS_NO_EVIDENCE", claim_id))
            if alignment_issue == "none":
                errors.append(issue("PARTIAL_CLAIM_NEEDS_ISSUE_CODE", claim_id))
        elif support in {"unsupported", "not_assessed"}:
            gap_ids.append(claim_id)
            if alignment_issue == "none" and support == "unsupported":
                errors.append(issue("UNSUPPORTED_CLAIM_NEEDS_ISSUE_CODE", claim_id))

        if support != "supported" and not row["requested_action"]:
            action_missing_ids.append(claim_id)
            warnings.append(issue("REQUESTED_ACTION_MISSING", claim_id))
        if alignment_issue != "none" and not row["limitation"]:
            warnings.append(issue("LIMITATION_CONTEXT_MISSING", claim_id))
        if row["claim_type"] in {"causal", "mechanistic"} and support == "supported":
            warnings.append(
                issue("CAUSAL_OR_MECHANISTIC_SUPPORT_REQUIRES_EXPERT_REVIEW", claim_id)
            )

    return {
        "schema_version": "2.0",
        "valid": not errors,
        "status": (
            "INVALID_MATRIX"
            if errors
            else "VALID_WITH_ALIGNMENT_GAPS"
            if gap_ids
            else "VALID_NO_RECORDED_GAPS"
        ),
        "errors": errors,
        "warnings": warnings,
        "claim_count": len(rows),
        "support_counts": dict(sorted(support_counts.items())),
        "alignment_issue_counts": dict(sorted(issue_counts.items())),
        "claim_ids_requiring_resolution": sorted(gap_ids),
        "claim_ids_missing_requested_action": sorted(action_missing_ids),
        "notice": (
            "This report checks matrix structure and declared claim-evidence "
            "alignment only. It does not verify evidence truth, reproduce analyses, "
            "or determine manuscript merit. Reports contain identifiers, not claim "
            "or manuscript text."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate a local claim-evidence CSV and emit only identifiers and "
            "counts; raw claim text is never echoed."
        )
    )
    parser.add_argument("matrix", help="Local claim-evidence matrix CSV")
    parser.add_argument("-o", "--output", help="Optional local JSON report")
    parser.add_argument(
        "--force", action="store_true", help="Replace an existing output file"
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        report = validate_matrix(load_matrix(args.matrix))
        write_json_report(report, args.output, force=args.force)
        return 0 if report["valid"] else 1
    except ValidationError as exc:
        return error_exit(exc)


if __name__ == "__main__":
    raise SystemExit(main())
