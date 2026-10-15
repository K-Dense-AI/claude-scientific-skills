#!/usr/bin/env python3
"""Audit a local evidence ledger and dated search boundary without networking."""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from datetime import date
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
    require_https_url,
    require_identifier,
    require_iso_date,
    require_object,
    require_partial_date,
    require_text,
    require_text_list,
    require_unique,
    split_identifiers,
    write_json_report,
)
from validate_hypothesis_schema import load_hypothesis_record, validate_record  # noqa: E402

FIELDS = (
    "source_id",
    "claim_ids",
    "title",
    "authors_or_organization",
    "publication_date",
    "source_type",
    "identifier",
    "url",
    "accessed_on",
    "relation",
    "study_design_or_document_type",
    "limitations",
    "notes",
)
SOURCE_TYPES = {
    "official_guidance",
    "regulation_policy",
    "reporting_guideline",
    "primary_research",
    "primary_method",
    "systematic_review",
    "consensus_report",
    "book",
    "dataset_or_registry",
    "preprint",
    "other",
}
RELATIONS = {"supportive", "challenging", "contextual", "method", "safety", "mixed"}
NOVELTY_STATES = {
    "not_assessed",
    "requires_specialist_review",
    "supported_by_documented_comprehensive_search",
}


def load_ledger(raw_path: str) -> list[dict[str, Any]]:
    rows = read_csv_records(raw_path, fields=FIELDS)
    parsed: list[dict[str, Any]] = []
    source_ids: list[str] = []
    for line_number, row in enumerate(rows, start=2):
        context = f"ledger row {line_number}"
        source_id = require_identifier(row["source_id"], f"{context}.source_id")
        source_ids.append(source_id)
        parsed.append(
            {
                "source_id": source_id,
                "claim_ids": split_identifiers(
                    row["claim_ids"], f"{context}.claim_ids"
                ),
                "title": require_text(
                    row["title"], f"{context}.title", minimum=5, maximum=2_000
                ),
                "authors_or_organization": require_text(
                    row["authors_or_organization"],
                    f"{context}.authors_or_organization",
                    minimum=2,
                    maximum=1_000,
                ),
                "publication_date": require_partial_date(
                    row["publication_date"], f"{context}.publication_date"
                ),
                "source_type": require_enum(
                    row["source_type"], SOURCE_TYPES, f"{context}.source_type"
                ),
                "identifier": require_text(
                    row["identifier"],
                    f"{context}.identifier",
                    minimum=3,
                    maximum=1_000,
                ),
                "url": require_https_url(row["url"], f"{context}.url"),
                "accessed_on": require_iso_date(
                    row["accessed_on"], f"{context}.accessed_on"
                ),
                "relation": require_enum(
                    row["relation"], RELATIONS, f"{context}.relation"
                ),
                "study_design_or_document_type": require_text(
                    row["study_design_or_document_type"],
                    f"{context}.study_design_or_document_type",
                    minimum=3,
                    maximum=1_000,
                ),
                "limitations": require_text(
                    row["limitations"],
                    f"{context}.limitations",
                    minimum=3,
                    maximum=2_000,
                ),
                "notes": require_text(
                    row["notes"],
                    f"{context}.notes",
                    allow_empty=True,
                    maximum=2_000,
                ),
            }
        )
    require_unique(source_ids, "evidence ledger")
    return parsed


def load_search_boundary(payload: Any) -> dict[str, Any]:
    root = require_object(payload, "search_boundary")
    fields = {
        "schema_version",
        "search_boundary_id",
        "searched_on",
        "searched_by",
        "purpose",
        "databases_or_indexes",
        "queries",
        "date_limits",
        "language_limits",
        "inclusion_scope",
        "exclusion_scope",
        "known_limitations",
        "last_result_screened_or_stop_rule",
        "novelty_status",
    }
    require_exact_keys(root, required=fields, context="search_boundary")
    return {
        "schema_version": require_enum(
            root["schema_version"], {"2.0"}, "search_boundary.schema_version"
        ),
        "search_boundary_id": require_identifier(
            root["search_boundary_id"], "search_boundary.search_boundary_id"
        ),
        "searched_on": require_iso_date(
            root["searched_on"], "search_boundary.searched_on"
        ),
        "searched_by": require_text(
            root["searched_by"], "search_boundary.searched_by", minimum=3
        ),
        "purpose": require_text(
            root["purpose"], "search_boundary.purpose", minimum=10
        ),
        "databases_or_indexes": require_text_list(
            root["databases_or_indexes"],
            "search_boundary.databases_or_indexes",
            minimum=1,
            maximum=100,
        ),
        "queries": require_text_list(
            root["queries"], "search_boundary.queries", minimum=1, maximum=200
        ),
        "date_limits": require_text(
            root["date_limits"], "search_boundary.date_limits", minimum=3
        ),
        "language_limits": require_text(
            root["language_limits"], "search_boundary.language_limits", minimum=3
        ),
        "inclusion_scope": require_text(
            root["inclusion_scope"], "search_boundary.inclusion_scope", minimum=10
        ),
        "exclusion_scope": require_text(
            root["exclusion_scope"], "search_boundary.exclusion_scope", minimum=10
        ),
        "known_limitations": require_text_list(
            root["known_limitations"],
            "search_boundary.known_limitations",
            minimum=1,
            maximum=100,
        ),
        "last_result_screened_or_stop_rule": require_text(
            root["last_result_screened_or_stop_rule"],
            "search_boundary.last_result_screened_or_stop_rule",
            minimum=10,
        ),
        "novelty_status": require_enum(
            root["novelty_status"],
            NOVELTY_STATES,
            "search_boundary.novelty_status",
        ),
    }


def _record_claim_ids(record: dict[str, Any]) -> set[str]:
    identifiers = {"OBS1", record["project_id"]}
    identifiers.update(item["hypothesis_id"] for item in record["hypotheses"])
    identifiers.update(item["estimand_id"] for item in record["causal_estimands"])
    identifiers.update(item["prediction_id"] for item in record["predictions"])
    identifiers.update(
        item["alternative_id"] for item in record["alternative_explanations"]
    )
    identifiers.update(item["null_id"] for item in record["null_hypotheses"])
    identifiers.update(item["control_id"] for item in record["negative_controls"])
    identifiers.update(
        item["measurement_id"] for item in record["operationalizations"]
    )
    identifiers.update(
        item["analysis_id"] for item in record["analysis_plan"]["analyses"]
    )
    return identifiers


def audit(
    ledger: list[dict[str, Any]],
    boundary: dict[str, Any],
    record: dict[str, Any] | None = None,
) -> dict[str, Any]:
    errors: list[dict[str, str]] = []
    warnings: list[dict[str, str]] = []
    source_type_counts: Counter[str] = Counter()
    relation_counts: Counter[str] = Counter()
    ledger_source_ids = {row["source_id"] for row in ledger}
    linked_claim_ids: set[str] = set()

    searched_on = date.fromisoformat(boundary["searched_on"])
    for row in ledger:
        source_type_counts[row["source_type"]] += 1
        relation_counts[row["relation"]] += 1
        linked_claim_ids.update(row["claim_ids"])
        if date.fromisoformat(row["accessed_on"]) > searched_on:
            warnings.append(
                issue("SOURCE_ACCESSED_AFTER_SEARCH_BOUNDARY", row["source_id"])
            )

    if not relation_counts["challenging"]:
        warnings.append(issue("NO_CHALLENGING_SOURCE_DECLARED", "ledger"))
    if boundary["novelty_status"] == "supported_by_documented_comprehensive_search":
        warnings.append(
            issue("NOVELTY_STATUS_REQUIRES_QUALIFIED_HUMAN_REVIEW", "search_boundary")
        )

    record_cross_check = record is not None
    if record is not None:
        record_report = validate_record(record)
        if not record_report["valid"]:
            raise ValidationError(
                "optional hypothesis record must pass schema validation first"
            )
        if record["evidence"]["search_boundary_id"] != boundary["search_boundary_id"]:
            errors.append(
                issue("SEARCH_BOUNDARY_ID_MISMATCH", boundary["search_boundary_id"])
            )
        record_source_ids = set(record["evidence"]["source_ids"])
        for source_id in sorted(record_source_ids - ledger_source_ids):
            errors.append(issue("RECORD_SOURCE_MISSING_FROM_LEDGER", source_id))
        for source_id in sorted(ledger_source_ids - record_source_ids):
            warnings.append(issue("LEDGER_SOURCE_NOT_DECLARED_IN_RECORD", source_id))
        allowed_claim_ids = _record_claim_ids(record)
        for claim_id in sorted(linked_claim_ids - allowed_claim_ids):
            errors.append(issue("UNKNOWN_CLAIM_ID", claim_id))

    return {
        "schema_version": "2.0",
        "search_boundary_id": boundary["search_boundary_id"],
        "searched_on": boundary["searched_on"],
        "novelty_status": boundary["novelty_status"],
        "valid": not errors,
        "status": "INVALID_LEDGER" if errors else "VALID_FOR_HUMAN_SOURCE_REVIEW",
        "errors": errors,
        "warnings": warnings,
        "source_count": len(ledger),
        "source_type_counts": dict(sorted(source_type_counts.items())),
        "relation_counts": dict(sorted(relation_counts.items())),
        "source_ids": sorted(ledger_source_ids),
        "linked_claim_ids": sorted(linked_claim_ids),
        "record_cross_check_performed": record_cross_check,
        "notice": (
            "This local audit does not visit URLs, verify source existence or "
            "content, appraise evidence, establish novelty, or determine whether "
            "a source supports a claim. Verify every source and link manually."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Audit a bounded local evidence CSV and search-boundary JSON without "
            "network access or scientific scoring."
        )
    )
    parser.add_argument("ledger", help="Local evidence ledger CSV")
    parser.add_argument("search_boundary", help="Local search-boundary JSON")
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
        report = audit(
            load_ledger(args.ledger),
            load_search_boundary(read_json(args.search_boundary)),
            record,
        )
        write_json_report(report, args.output, force=args.force)
        return 0 if report["valid"] else 1
    except ValidationError as exc:
        return error_exit(exc)


if __name__ == "__main__":
    raise SystemExit(main())
