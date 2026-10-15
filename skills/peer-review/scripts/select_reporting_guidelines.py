#!/usr/bin/env python3
"""Select bundled reporting guidance and audit checklist coverage locally."""

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
    read_csv_records,
    read_json,
    require_enum,
    require_exact_keys,
    require_identifier,
    require_identifier_list,
    require_object,
    require_text,
    require_url,
    write_json_report,
)

ASSET_CATALOG = Path(__file__).resolve().parents[1] / "assets" / "reporting_guidelines.json"
COVERAGE_FIELDS = (
    "guideline_id",
    "item_id",
    "status",
    "location",
    "rationale",
)
COVERAGE_STATUSES = {
    "reported",
    "partly_reported",
    "not_reported",
    "not_applicable",
    "not_assessed",
}
CATALOG_STATUSES = {"current", "legacy_current_qualified"}
CATEGORIES = {"reporting_guideline", "domain_metadata_standard"}
REPORT_KINDS = {"results", "protocol", "abstract", "data_release"}
ITEM_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,63}$")


def load_catalog(path: Path = ASSET_CATALOG) -> dict[str, Any]:
    payload = require_object(read_json(path), "catalog")
    require_exact_keys(
        payload,
        required={"schema_version", "reviewed_on", "notice", "guidelines"},
        context="catalog",
    )
    require_enum(payload["schema_version"], {"2.0"}, "catalog.schema_version")
    require_text(payload["reviewed_on"], "catalog.reviewed_on", maximum=10)
    require_text(payload["notice"], "catalog.notice", maximum=2_000)
    guidelines = payload["guidelines"]
    if not isinstance(guidelines, list) or not guidelines:
        raise ValidationError("catalog.guidelines must be a non-empty array")

    parsed: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, raw_entry in enumerate(guidelines):
        context = f"catalog.guidelines[{index}]"
        entry = require_object(raw_entry, context)
        require_exact_keys(
            entry,
            required={
                "id",
                "name",
                "version",
                "status",
                "category",
                "study_types",
                "report_kinds",
                "required_features",
                "domains",
                "main_item_count",
                "url",
                "notes",
            },
            context=context,
        )
        guideline_id = require_identifier(entry["id"], f"{context}.id")
        if guideline_id in seen:
            raise ValidationError(f"catalog contains duplicate ID: {guideline_id}")
        seen.add(guideline_id)
        item_count = entry["main_item_count"]
        if item_count is not None:
            if type(item_count) is not int or not 1 <= item_count <= 200:
                raise ValidationError(
                    f"{context}.main_item_count must be null or an integer 1-200"
                )
        parsed.append(
            {
                "id": guideline_id,
                "name": require_text(entry["name"], f"{context}.name", maximum=200),
                "version": require_text(
                    entry["version"], f"{context}.version", maximum=100
                ),
                "status": require_enum(
                    entry["status"], CATALOG_STATUSES, f"{context}.status"
                ),
                "category": require_enum(
                    entry["category"], CATEGORIES, f"{context}.category"
                ),
                "study_types": require_identifier_list(
                    entry["study_types"], f"{context}.study_types", minimum=1
                ),
                "report_kinds": [
                    require_enum(item, REPORT_KINDS, f"{context}.report_kinds")
                    for item in entry["report_kinds"]
                ],
                "required_features": require_identifier_list(
                    entry["required_features"], f"{context}.required_features"
                ),
                "domains": require_identifier_list(
                    entry["domains"], f"{context}.domains", minimum=1
                ),
                "main_item_count": item_count,
                "url": require_url(entry["url"], f"{context}.url"),
                "notes": require_text(
                    entry["notes"], f"{context}.notes", maximum=1_000
                ),
            }
        )
    return {
        "schema_version": payload["schema_version"],
        "reviewed_on": payload["reviewed_on"],
        "notice": payload["notice"],
        "guidelines": parsed,
    }


def load_profile(payload: Any) -> dict[str, Any]:
    profile = require_object(payload, "profile")
    require_exact_keys(
        profile,
        required={
            "schema_version",
            "profile_id",
            "study_types",
            "report_kind",
            "features",
            "domains",
        },
        context="profile",
    )
    return {
        "schema_version": require_enum(
            profile["schema_version"], {"2.0"}, "profile.schema_version"
        ),
        "profile_id": require_identifier(profile["profile_id"], "profile.profile_id"),
        "study_types": require_identifier_list(
            profile["study_types"], "profile.study_types", minimum=1
        ),
        "report_kind": require_enum(
            profile["report_kind"], REPORT_KINDS, "profile.report_kind"
        ),
        "features": require_identifier_list(profile["features"], "profile.features"),
        "domains": require_identifier_list(
            profile["domains"], "profile.domains", minimum=1
        ),
    }


def select_guidelines(
    profile: dict[str, Any], catalog: dict[str, Any]
) -> list[dict[str, Any]]:
    study_types = set(profile["study_types"])
    features = set(profile["features"])
    domains = set(profile["domains"])
    selected: list[dict[str, Any]] = []
    for entry in catalog["guidelines"]:
        if not study_types.intersection(entry["study_types"]):
            continue
        if profile["report_kind"] not in entry["report_kinds"]:
            continue
        if not set(entry["required_features"]).issubset(features):
            continue
        if "all" not in entry["domains"] and not domains.intersection(entry["domains"]):
            continue
        selected.append(entry)
    return sorted(
        selected,
        key=lambda item: (
            item["category"] != "reporting_guideline",
            item["id"],
        ),
    )


def load_coverage(raw_path: str | Path) -> list[dict[str, str]]:
    rows = read_csv_records(
        raw_path,
        required_fields=COVERAGE_FIELDS,
        allowed_fields=COVERAGE_FIELDS,
    )
    seen: set[tuple[str, str]] = set()
    parsed: list[dict[str, str]] = []
    for index, row in enumerate(rows, start=2):
        guideline_id = require_identifier(
            row["guideline_id"], f"coverage row {index}.guideline_id"
        )
        item_id = require_text(
            row["item_id"], f"coverage row {index}.item_id", maximum=64
        )
        if not ITEM_ID_RE.fullmatch(item_id):
            raise ValidationError(
                f"coverage row {index}.item_id has an invalid identifier format"
            )
        key = (guideline_id, item_id)
        if key in seen:
            raise ValidationError(
                f"coverage contains duplicate guideline/item pair: "
                f"{guideline_id}/{item_id}"
            )
        seen.add(key)
        status = require_enum(
            row["status"], COVERAGE_STATUSES, f"coverage row {index}.status"
        )
        location = require_text(
            row["location"],
            f"coverage row {index}.location",
            allow_empty=True,
            maximum=500,
        )
        rationale = require_text(
            row["rationale"],
            f"coverage row {index}.rationale",
            allow_empty=True,
            maximum=1_000,
        )
        parsed.append(
            {
                "guideline_id": guideline_id,
                "item_id": item_id,
                "status": status,
                "location": location,
                "rationale": rationale,
            }
        )
    return parsed


def assess(
    profile: dict[str, Any],
    catalog: dict[str, Any],
    coverage_rows: list[dict[str, str]] | None = None,
) -> dict[str, Any]:
    selected = select_guidelines(profile, catalog)
    selected_by_id = {entry["id"]: entry for entry in selected}
    errors: list[dict[str, str]] = []
    warnings: list[dict[str, str]] = []
    if not selected:
        warnings.append(issue("NO_BUNDLED_GUIDELINE_MATCH", "profile.study_types"))

    rows_by_guideline: dict[str, list[dict[str, str]]] = {}
    if coverage_rows is not None:
        for row in coverage_rows:
            rows_by_guideline.setdefault(row["guideline_id"], []).append(row)
            subject = f"{row['guideline_id']}:{row['item_id']}"
            if row["status"] in {"reported", "partly_reported"} and not row["location"]:
                errors.append(issue("COVERAGE_LOCATION_REQUIRED", subject))
            if row["status"] == "not_applicable" and not row["rationale"]:
                errors.append(issue("NOT_APPLICABLE_RATIONALE_REQUIRED", subject))
            if row["guideline_id"] not in selected_by_id:
                warnings.append(issue("COVERAGE_FOR_UNSELECTED_GUIDELINE", subject))

    coverage: list[dict[str, Any]] = []
    for entry in selected:
        rows = rows_by_guideline.get(entry["id"], [])
        statuses = Counter(row["status"] for row in rows)
        observed_ids = {row["item_id"] for row in rows}
        expected_count = entry["main_item_count"]
        unrecorded: list[str] = []
        unexpected: list[str] = []
        record_complete: bool | None = None
        if expected_count is not None:
            expected_ids = {str(number) for number in range(1, expected_count + 1)}
            unrecorded = sorted(
                expected_ids - observed_ids, key=lambda value: int(value)
            )
            unexpected = sorted(observed_ids - expected_ids)
            record_complete = not unrecorded and not unexpected
            if unexpected:
                warnings.append(
                    issue("NON_MAIN_ITEM_IDS_RECORDED", entry["id"])
                )
        elif rows:
            record_complete = None
            warnings.append(
                issue("OFFICIAL_ITEM_SET_MUST_BE_VERIFIED", entry["id"])
            )

        gap_ids = sorted(
            row["item_id"]
            for row in rows
            if row["status"]
            in {"partly_reported", "not_reported", "not_assessed"}
        )
        coverage.append(
            {
                "guideline_id": entry["id"],
                "main_item_count": expected_count,
                "recorded_item_count": len(rows),
                "coverage_record_complete": record_complete,
                "status_counts": dict(sorted(statuses.items())),
                "unrecorded_main_item_ids": unrecorded,
                "unexpected_item_ids": unexpected,
                "reporting_gap_item_ids": gap_ids,
            }
        )

    return {
        "schema_version": "2.0",
        "profile_id": profile["profile_id"],
        "valid": not errors,
        "status": "VALID" if not errors else "INVALID_COVERAGE_RECORD",
        "errors": errors,
        "warnings": warnings,
        "selected_guidelines": [
            {
                "id": entry["id"],
                "name": entry["name"],
                "version": entry["version"],
                "status": entry["status"],
                "category": entry["category"],
                "url": entry["url"],
                "notes": entry["notes"],
            }
            for entry in selected
        ],
        "coverage": coverage if coverage_rows is not None else None,
        "catalog_reviewed_on": catalog["reviewed_on"],
        "notice": (
            "Guideline selection and checklist coverage concern reporting "
            "completeness only. They are not scores and do not establish study "
            "quality, validity, conduct, or manuscript merit. Check the official "
            "guideline, applicable extensions, and target venue policy."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Select current bundled reporting guidance and optionally audit a "
            "local checklist coverage CSV. No network calls are made."
        )
    )
    parser.add_argument("profile", help="Local study profile JSON")
    parser.add_argument(
        "--coverage",
        help=(
            "Optional coverage CSV using aggregate main item IDs (for example, "
            "1 through 30 for CONSORT 2025)"
        ),
    )
    parser.add_argument("-o", "--output", help="Optional local JSON report")
    parser.add_argument(
        "--force", action="store_true", help="Replace an existing output file"
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        profile = load_profile(read_json(args.profile))
        catalog = load_catalog()
        coverage = load_coverage(args.coverage) if args.coverage else None
        report = assess(profile, catalog, coverage)
        write_json_report(report, args.output, force=args.force)
        return 0 if report["valid"] else 1
    except ValidationError as exc:
        return error_exit(exc)


if __name__ == "__main__":
    raise SystemExit(main())
