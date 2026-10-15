#!/usr/bin/env python3
"""Validate a structured hypothesis record without judging scientific merit."""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path, PurePosixPath
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
    require_iso_date,
    require_list,
    require_object,
    require_text,
    require_text_list,
    require_unique,
    write_json_report,
)

QUESTION_TYPES = {"descriptive", "associational", "predictive", "causal", "mechanistic"}
HYPOTHESIS_STATUS = {"candidate"}
MEASUREMENT_ROLES = {
    "intervention",
    "exposure",
    "outcome",
    "mediator",
    "confounder",
    "selection",
    "negative_control",
    "positive_control",
    "effect_modifier",
    "other",
}
RISK_TYPES = {
    "confounding",
    "selection_bias",
    "collider_bias",
    "reverse_causation",
    "measurement_bias",
    "other",
}
CONTROL_TYPES = {
    "negative_exposure",
    "negative_outcome",
    "procedural_negative",
    "positive_control",
    "vehicle_or_sham",
    "other",
}
PLAN_MODES = {"confirmatory", "exploratory"}
GATE_STATES = {"not_applicable", "undetermined", "requires_review", "approved", "blocked"}
FEASIBILITY_STATES = {
    "undetermined",
    "feasible_for_planning",
    "requires_pilot",
    "infeasible",
    "blocked",
}


def _parse_observation(raw: Any) -> dict[str, Any]:
    value = require_object(raw, "record.observation")
    require_exact_keys(
        value,
        required={"statement", "provenance", "source_ids", "uncertainties"},
        context="record.observation",
    )
    return {
        "statement": require_text(
            value["statement"], "record.observation.statement", minimum=10
        ),
        "provenance": require_text(
            value["provenance"], "record.observation.provenance", minimum=10
        ),
        "source_ids": require_identifier_list(
            value["source_ids"], "record.observation.source_ids", minimum=1
        ),
        "uncertainties": require_text_list(
            value["uncertainties"],
            "record.observation.uncertainties",
            minimum=1,
        ),
    }


def _parse_question(raw: Any) -> dict[str, str]:
    value = require_object(raw, "record.research_question")
    fields = {
        "statement",
        "framework",
        "question_type",
        "population_or_system",
        "intervention_or_exposure",
        "comparator",
        "outcome",
        "timeframe",
    }
    require_exact_keys(value, required=fields, context="record.research_question")
    parsed = {
        field: require_text(
            value[field],
            f"record.research_question.{field}",
            minimum=2 if field == "framework" else 5,
        )
        for field in fields
        if field != "question_type"
    }
    parsed["question_type"] = require_enum(
        value["question_type"],
        QUESTION_TYPES,
        "record.research_question.question_type",
    )
    return parsed


def _parse_hypotheses(raw: Any) -> list[dict[str, Any]]:
    entries = require_list(raw, "record.hypotheses", minimum=1, maximum=50)
    parsed: list[dict[str, Any]] = []
    identifiers: list[str] = []
    fields = {
        "hypothesis_id",
        "statement",
        "mechanism",
        "status",
        "source_ids",
        "assumptions",
        "boundary_conditions",
        "uncertainties",
        "prediction_ids",
        "rival_hypothesis_ids",
    }
    for index, raw_entry in enumerate(entries):
        context = f"record.hypotheses[{index}]"
        entry = require_object(raw_entry, context)
        require_exact_keys(entry, required=fields, context=context)
        identifier = require_identifier(entry["hypothesis_id"], f"{context}.hypothesis_id")
        identifiers.append(identifier)
        parsed.append(
            {
                "hypothesis_id": identifier,
                "statement": require_text(
                    entry["statement"], f"{context}.statement", minimum=10
                ),
                "mechanism": require_text(
                    entry["mechanism"], f"{context}.mechanism", minimum=10
                ),
                "status": require_enum(
                    entry["status"], HYPOTHESIS_STATUS, f"{context}.status"
                ),
                "source_ids": require_identifier_list(
                    entry["source_ids"], f"{context}.source_ids", minimum=1
                ),
                "assumptions": require_text_list(
                    entry["assumptions"], f"{context}.assumptions", minimum=1
                ),
                "boundary_conditions": require_text_list(
                    entry["boundary_conditions"],
                    f"{context}.boundary_conditions",
                    minimum=1,
                ),
                "uncertainties": require_text_list(
                    entry["uncertainties"], f"{context}.uncertainties", minimum=1
                ),
                "prediction_ids": require_identifier_list(
                    entry["prediction_ids"], f"{context}.prediction_ids", minimum=1
                ),
                "rival_hypothesis_ids": require_identifier_list(
                    entry["rival_hypothesis_ids"],
                    f"{context}.rival_hypothesis_ids",
                ),
            }
        )
    require_unique(identifiers, "record.hypotheses")
    return parsed


def _parse_estimands(raw: Any) -> list[dict[str, Any]]:
    entries = require_list(raw, "record.causal_estimands", maximum=50)
    parsed: list[dict[str, Any]] = []
    identifiers: list[str] = []
    fields = {
        "estimand_id",
        "linked_hypothesis_ids",
        "population",
        "intervention_or_exposure",
        "comparator",
        "outcome",
        "time_horizon",
        "population_summary",
        "intercurrent_event_strategy",
        "identification_assumptions",
    }
    for index, raw_entry in enumerate(entries):
        context = f"record.causal_estimands[{index}]"
        entry = require_object(raw_entry, context)
        require_exact_keys(entry, required=fields, context=context)
        identifier = require_identifier(entry["estimand_id"], f"{context}.estimand_id")
        identifiers.append(identifier)
        parsed.append(
            {
                "estimand_id": identifier,
                "linked_hypothesis_ids": require_identifier_list(
                    entry["linked_hypothesis_ids"],
                    f"{context}.linked_hypothesis_ids",
                    minimum=1,
                ),
                "population": require_text(entry["population"], f"{context}.population"),
                "intervention_or_exposure": require_text(
                    entry["intervention_or_exposure"],
                    f"{context}.intervention_or_exposure",
                ),
                "comparator": require_text(
                    entry["comparator"], f"{context}.comparator"
                ),
                "outcome": require_text(entry["outcome"], f"{context}.outcome"),
                "time_horizon": require_text(
                    entry["time_horizon"], f"{context}.time_horizon"
                ),
                "population_summary": require_text(
                    entry["population_summary"], f"{context}.population_summary"
                ),
                "intercurrent_event_strategy": require_text(
                    entry["intercurrent_event_strategy"],
                    f"{context}.intercurrent_event_strategy",
                ),
                "identification_assumptions": require_text_list(
                    entry["identification_assumptions"],
                    f"{context}.identification_assumptions",
                    minimum=1,
                ),
            }
        )
    require_unique(identifiers, "record.causal_estimands")
    return parsed


def _parse_predictions(raw: Any) -> list[dict[str, Any]]:
    entries = require_list(raw, "record.predictions", minimum=1, maximum=200)
    parsed: list[dict[str, Any]] = []
    identifiers: list[str] = []
    fields = {
        "prediction_id",
        "hypothesis_id",
        "statement",
        "observable",
        "conditions",
        "expected_pattern",
        "falsifier",
        "rival_hypothesis_ids",
        "measurement_ids",
        "analysis_ids",
    }
    for index, raw_entry in enumerate(entries):
        context = f"record.predictions[{index}]"
        entry = require_object(raw_entry, context)
        require_exact_keys(entry, required=fields, context=context)
        identifier = require_identifier(entry["prediction_id"], f"{context}.prediction_id")
        identifiers.append(identifier)
        parsed.append(
            {
                "prediction_id": identifier,
                "hypothesis_id": require_identifier(
                    entry["hypothesis_id"], f"{context}.hypothesis_id"
                ),
                "statement": require_text(
                    entry["statement"], f"{context}.statement", minimum=10
                ),
                "observable": require_text(
                    entry["observable"], f"{context}.observable", minimum=5
                ),
                "conditions": require_text(
                    entry["conditions"], f"{context}.conditions", minimum=5
                ),
                "expected_pattern": require_text(
                    entry["expected_pattern"], f"{context}.expected_pattern", minimum=5
                ),
                "falsifier": require_text(
                    entry["falsifier"], f"{context}.falsifier", minimum=10
                ),
                "rival_hypothesis_ids": require_identifier_list(
                    entry["rival_hypothesis_ids"],
                    f"{context}.rival_hypothesis_ids",
                    minimum=1,
                ),
                "measurement_ids": require_identifier_list(
                    entry["measurement_ids"],
                    f"{context}.measurement_ids",
                    minimum=1,
                ),
                "analysis_ids": require_identifier_list(
                    entry["analysis_ids"], f"{context}.analysis_ids", minimum=1
                ),
            }
        )
    require_unique(identifiers, "record.predictions")
    return parsed


def _parse_alternatives(raw: Any) -> list[dict[str, Any]]:
    entries = require_list(raw, "record.alternative_explanations", minimum=1, maximum=100)
    parsed: list[dict[str, Any]] = []
    identifiers: list[str] = []
    fields = {
        "alternative_id",
        "statement",
        "linked_hypothesis_ids",
        "risk_types",
        "discriminating_prediction_ids",
    }
    for index, raw_entry in enumerate(entries):
        context = f"record.alternative_explanations[{index}]"
        entry = require_object(raw_entry, context)
        require_exact_keys(entry, required=fields, context=context)
        identifier = require_identifier(entry["alternative_id"], f"{context}.alternative_id")
        identifiers.append(identifier)
        risk_types = require_text_list(
            entry["risk_types"], f"{context}.risk_types", minimum=1, maximum=10
        )
        for risk_index, risk_type in enumerate(risk_types):
            require_enum(risk_type, RISK_TYPES, f"{context}.risk_types[{risk_index}]")
        require_unique(risk_types, f"{context}.risk_types")
        parsed.append(
            {
                "alternative_id": identifier,
                "statement": require_text(
                    entry["statement"], f"{context}.statement", minimum=10
                ),
                "linked_hypothesis_ids": require_identifier_list(
                    entry["linked_hypothesis_ids"],
                    f"{context}.linked_hypothesis_ids",
                    minimum=1,
                ),
                "risk_types": risk_types,
                "discriminating_prediction_ids": require_identifier_list(
                    entry["discriminating_prediction_ids"],
                    f"{context}.discriminating_prediction_ids",
                    minimum=1,
                ),
            }
        )
    require_unique(identifiers, "record.alternative_explanations")
    return parsed


def _parse_nulls(raw: Any) -> list[dict[str, Any]]:
    entries = require_list(raw, "record.null_hypotheses", minimum=1, maximum=100)
    parsed: list[dict[str, Any]] = []
    identifiers: list[str] = []
    fields = {
        "null_id",
        "statement",
        "linked_prediction_ids",
        "rejection_or_compatibility_rule",
    }
    for index, raw_entry in enumerate(entries):
        context = f"record.null_hypotheses[{index}]"
        entry = require_object(raw_entry, context)
        require_exact_keys(entry, required=fields, context=context)
        identifier = require_identifier(entry["null_id"], f"{context}.null_id")
        identifiers.append(identifier)
        parsed.append(
            {
                "null_id": identifier,
                "statement": require_text(
                    entry["statement"], f"{context}.statement", minimum=10
                ),
                "linked_prediction_ids": require_identifier_list(
                    entry["linked_prediction_ids"],
                    f"{context}.linked_prediction_ids",
                    minimum=1,
                ),
                "rejection_or_compatibility_rule": require_text(
                    entry["rejection_or_compatibility_rule"],
                    f"{context}.rejection_or_compatibility_rule",
                    minimum=10,
                ),
            }
        )
    require_unique(identifiers, "record.null_hypotheses")
    return parsed


def _parse_controls(raw: Any) -> list[dict[str, Any]]:
    entries = require_list(raw, "record.negative_controls", minimum=1, maximum=100)
    parsed: list[dict[str, Any]] = []
    identifiers: list[str] = []
    fields = {
        "control_id",
        "control_type",
        "rationale",
        "expected_result",
        "failure_implication",
        "linked_prediction_ids",
    }
    for index, raw_entry in enumerate(entries):
        context = f"record.negative_controls[{index}]"
        entry = require_object(raw_entry, context)
        require_exact_keys(entry, required=fields, context=context)
        identifier = require_identifier(entry["control_id"], f"{context}.control_id")
        identifiers.append(identifier)
        parsed.append(
            {
                "control_id": identifier,
                "control_type": require_enum(
                    entry["control_type"], CONTROL_TYPES, f"{context}.control_type"
                ),
                "rationale": require_text(
                    entry["rationale"], f"{context}.rationale", minimum=10
                ),
                "expected_result": require_text(
                    entry["expected_result"], f"{context}.expected_result", minimum=5
                ),
                "failure_implication": require_text(
                    entry["failure_implication"],
                    f"{context}.failure_implication",
                    minimum=10,
                ),
                "linked_prediction_ids": require_identifier_list(
                    entry["linked_prediction_ids"],
                    f"{context}.linked_prediction_ids",
                    minimum=1,
                ),
            }
        )
    require_unique(identifiers, "record.negative_controls")
    return parsed


def _parse_operationalizations(raw: Any) -> list[dict[str, Any]]:
    entries = require_list(raw, "record.operationalizations", minimum=1, maximum=200)
    parsed: list[dict[str, Any]] = []
    identifiers: list[str] = []
    fields = {
        "measurement_id",
        "construct",
        "variable",
        "role",
        "operational_definition",
        "instrument_or_method",
        "unit",
        "timing",
        "population_or_system",
        "validity_evidence_source_ids",
        "reliability_plan",
        "missingness_plan",
        "blinding_or_masking",
        "threshold_rationale",
    }
    for index, raw_entry in enumerate(entries):
        context = f"record.operationalizations[{index}]"
        entry = require_object(raw_entry, context)
        require_exact_keys(entry, required=fields, context=context)
        identifier = require_identifier(entry["measurement_id"], f"{context}.measurement_id")
        identifiers.append(identifier)
        text_fields = fields - {
            "measurement_id",
            "role",
            "validity_evidence_source_ids",
        }
        parsed_entry: dict[str, Any] = {
            "measurement_id": identifier,
            "role": require_enum(entry["role"], MEASUREMENT_ROLES, f"{context}.role"),
            "validity_evidence_source_ids": require_identifier_list(
                entry["validity_evidence_source_ids"],
                f"{context}.validity_evidence_source_ids",
                minimum=1,
            ),
        }
        for field in text_fields:
            parsed_entry[field] = require_text(
                entry[field],
                f"{context}.{field}",
                minimum=2 if field in {"unit", "variable"} else 5,
            )
        parsed.append(parsed_entry)
    require_unique(identifiers, "record.operationalizations")
    return parsed


def _parse_analysis_plan(raw: Any) -> dict[str, Any]:
    plan = require_object(raw, "record.analysis_plan")
    require_exact_keys(
        plan,
        required={"analyses", "harking_control", "deviation_reporting"},
        context="record.analysis_plan",
    )
    entries = require_list(
        plan["analyses"], "record.analysis_plan.analyses", minimum=1, maximum=200
    )
    parsed: list[dict[str, Any]] = []
    identifiers: list[str] = []
    fields = {
        "analysis_id",
        "prediction_ids",
        "estimand_ids",
        "analysis_population",
        "method",
        "effect_or_summary_measure",
        "uncertainty_method",
        "missing_data_plan",
        "multiplicity_plan",
        "sensitivity_analyses",
        "decision_rule",
        "exploratory_or_confirmatory",
    }
    for index, raw_entry in enumerate(entries):
        context = f"record.analysis_plan.analyses[{index}]"
        entry = require_object(raw_entry, context)
        require_exact_keys(entry, required=fields, context=context)
        identifier = require_identifier(entry["analysis_id"], f"{context}.analysis_id")
        identifiers.append(identifier)
        parsed.append(
            {
                "analysis_id": identifier,
                "prediction_ids": require_identifier_list(
                    entry["prediction_ids"], f"{context}.prediction_ids", minimum=1
                ),
                "estimand_ids": require_identifier_list(
                    entry["estimand_ids"], f"{context}.estimand_ids"
                ),
                "analysis_population": require_text(
                    entry["analysis_population"],
                    f"{context}.analysis_population",
                    minimum=5,
                ),
                "method": require_text(entry["method"], f"{context}.method", minimum=5),
                "effect_or_summary_measure": require_text(
                    entry["effect_or_summary_measure"],
                    f"{context}.effect_or_summary_measure",
                    minimum=3,
                ),
                "uncertainty_method": require_text(
                    entry["uncertainty_method"],
                    f"{context}.uncertainty_method",
                    minimum=5,
                ),
                "missing_data_plan": require_text(
                    entry["missing_data_plan"],
                    f"{context}.missing_data_plan",
                    minimum=5,
                ),
                "multiplicity_plan": require_text(
                    entry["multiplicity_plan"],
                    f"{context}.multiplicity_plan",
                    minimum=5,
                ),
                "sensitivity_analyses": require_text_list(
                    entry["sensitivity_analyses"],
                    f"{context}.sensitivity_analyses",
                    minimum=1,
                ),
                "decision_rule": require_text(
                    entry["decision_rule"], f"{context}.decision_rule", minimum=10
                ),
                "exploratory_or_confirmatory": require_enum(
                    entry["exploratory_or_confirmatory"],
                    PLAN_MODES,
                    f"{context}.exploratory_or_confirmatory",
                ),
            }
        )
    require_unique(identifiers, "record.analysis_plan.analyses")
    return {
        "analyses": parsed,
        "harking_control": require_text(
            plan["harking_control"],
            "record.analysis_plan.harking_control",
            minimum=10,
        ),
        "deviation_reporting": require_text(
            plan["deviation_reporting"],
            "record.analysis_plan.deviation_reporting",
            minimum=10,
        ),
    }


def _parse_evidence(raw: Any) -> dict[str, Any]:
    value = require_object(raw, "record.evidence")
    require_exact_keys(
        value,
        required={
            "ledger_path",
            "source_ids",
            "search_boundary_id",
            "evidence_limitations",
        },
        context="record.evidence",
    )
    ledger_path = require_text(
        value["ledger_path"], "record.evidence.ledger_path", maximum=500
    )
    pure_path = PurePosixPath(ledger_path)
    if pure_path.is_absolute() or ".." in pure_path.parts or "://" in ledger_path:
        raise ValidationError("record.evidence.ledger_path must be a safe relative path")
    if pure_path.suffix.lower() != ".csv":
        raise ValidationError("record.evidence.ledger_path must end in .csv")
    return {
        "ledger_path": ledger_path,
        "source_ids": require_identifier_list(
            value["source_ids"], "record.evidence.source_ids", minimum=1
        ),
        "search_boundary_id": require_identifier(
            value["search_boundary_id"], "record.evidence.search_boundary_id"
        ),
        "evidence_limitations": require_text_list(
            value["evidence_limitations"],
            "record.evidence.evidence_limitations",
            minimum=1,
        ),
    }


def _parse_risk_register(raw: Any) -> dict[str, list[str]]:
    value = require_object(raw, "record.risk_register")
    fields = {
        "confounding",
        "selection_bias",
        "collider_bias",
        "reverse_causation",
        "measurement_bias",
        "other",
    }
    require_exact_keys(value, required=fields, context="record.risk_register")
    return {
        field: require_text_list(
            value[field], f"record.risk_register.{field}", maximum=50
        )
        for field in fields
    }


def _parse_ethics(raw: Any) -> dict[str, Any]:
    value = require_object(raw, "record.ethics_and_feasibility")
    gate_fields = {
        "human_subjects_gate",
        "animal_research_gate",
        "biosafety_gate",
        "dual_use_gate",
        "regulatory_gate",
        "data_governance_gate",
    }
    require_exact_keys(
        value,
        required=gate_fields
        | {"feasibility_status", "required_reviews", "unresolved_blocks"},
        context="record.ethics_and_feasibility",
    )
    parsed: dict[str, Any] = {
        field: require_enum(
            value[field], GATE_STATES, f"record.ethics_and_feasibility.{field}"
        )
        for field in gate_fields
    }
    parsed["feasibility_status"] = require_enum(
        value["feasibility_status"],
        FEASIBILITY_STATES,
        "record.ethics_and_feasibility.feasibility_status",
    )
    parsed["required_reviews"] = require_text_list(
        value["required_reviews"],
        "record.ethics_and_feasibility.required_reviews",
        maximum=50,
    )
    parsed["unresolved_blocks"] = require_text_list(
        value["unresolved_blocks"],
        "record.ethics_and_feasibility.unresolved_blocks",
        maximum=50,
    )
    return parsed


def _parse_ai_use(raw: Any) -> dict[str, Any]:
    value = require_object(raw, "record.ai_use")
    fields = {
        "used",
        "sensitive_or_unpublished_data_sent_externally",
        "local_policy_checked",
        "citation_verification_required",
        "human_accountable",
        "diversity_mitigation",
    }
    require_exact_keys(value, required=fields, context="record.ai_use")
    return {
        "used": require_bool(value["used"], "record.ai_use.used"),
        "sensitive_or_unpublished_data_sent_externally": require_bool(
            value["sensitive_or_unpublished_data_sent_externally"],
            "record.ai_use.sensitive_or_unpublished_data_sent_externally",
        ),
        "local_policy_checked": require_bool(
            value["local_policy_checked"], "record.ai_use.local_policy_checked"
        ),
        "citation_verification_required": require_bool(
            value["citation_verification_required"],
            "record.ai_use.citation_verification_required",
        ),
        "human_accountable": require_bool(
            value["human_accountable"], "record.ai_use.human_accountable"
        ),
        "diversity_mitigation": require_text(
            value["diversity_mitigation"],
            "record.ai_use.diversity_mitigation",
            minimum=10,
        ),
    }


def load_hypothesis_record(payload: Any) -> dict[str, Any]:
    """Parse the exact v2 record schema and normalize bounded values."""
    root = require_object(payload, "record")
    fields = {
        "schema_version",
        "project_id",
        "status",
        "updated_on",
        "human_owner",
        "observation",
        "research_question",
        "hypotheses",
        "causal_estimands",
        "predictions",
        "alternative_explanations",
        "null_hypotheses",
        "negative_controls",
        "operationalizations",
        "analysis_plan",
        "evidence",
        "risk_register",
        "ethics_and_feasibility",
        "ai_use",
    }
    require_exact_keys(root, required=fields, context="record")
    return {
        "schema_version": require_enum(
            root["schema_version"], {"2.0"}, "record.schema_version"
        ),
        "project_id": require_identifier(root["project_id"], "record.project_id"),
        "status": require_enum(
            root["status"], {"draft", "preregistered", "archived"}, "record.status"
        ),
        "updated_on": require_iso_date(root["updated_on"], "record.updated_on"),
        "human_owner": require_text(
            root["human_owner"], "record.human_owner", minimum=3
        ),
        "observation": _parse_observation(root["observation"]),
        "research_question": _parse_question(root["research_question"]),
        "hypotheses": _parse_hypotheses(root["hypotheses"]),
        "causal_estimands": _parse_estimands(root["causal_estimands"]),
        "predictions": _parse_predictions(root["predictions"]),
        "alternative_explanations": _parse_alternatives(
            root["alternative_explanations"]
        ),
        "null_hypotheses": _parse_nulls(root["null_hypotheses"]),
        "negative_controls": _parse_controls(root["negative_controls"]),
        "operationalizations": _parse_operationalizations(root["operationalizations"]),
        "analysis_plan": _parse_analysis_plan(root["analysis_plan"]),
        "evidence": _parse_evidence(root["evidence"]),
        "risk_register": _parse_risk_register(root["risk_register"]),
        "ethics_and_feasibility": _parse_ethics(root["ethics_and_feasibility"]),
        "ai_use": _parse_ai_use(root["ai_use"]),
    }


def _missing_references(
    values: list[str], valid_values: set[str], code: str, field: str
) -> list[dict[str, str]]:
    return [
        issue(code, f"{field}:{value}")
        for value in sorted(set(values) - valid_values)
    ]


def validate_record(record: dict[str, Any]) -> dict[str, Any]:
    """Check cross-links and safety declarations without scoring candidates."""
    errors: list[dict[str, str]] = []
    warnings: list[dict[str, str]] = []

    hypothesis_ids = {item["hypothesis_id"] for item in record["hypotheses"]}
    estimand_ids = {item["estimand_id"] for item in record["causal_estimands"]}
    prediction_ids = {item["prediction_id"] for item in record["predictions"]}
    measurement_ids = {
        item["measurement_id"] for item in record["operationalizations"]
    }
    analysis_ids = {
        item["analysis_id"] for item in record["analysis_plan"]["analyses"]
    }
    control_ids = {item["control_id"] for item in record["negative_controls"]}
    source_ids = set(record["evidence"]["source_ids"])

    used_sources = set(record["observation"]["source_ids"])
    for hypothesis in record["hypotheses"]:
        used_sources.update(hypothesis["source_ids"])
        errors.extend(
            _missing_references(
                hypothesis["prediction_ids"],
                prediction_ids,
                "UNKNOWN_PREDICTION_ID",
                hypothesis["hypothesis_id"],
            )
        )
        errors.extend(
            _missing_references(
                hypothesis["rival_hypothesis_ids"],
                hypothesis_ids,
                "UNKNOWN_RIVAL_HYPOTHESIS_ID",
                hypothesis["hypothesis_id"],
            )
        )
        if hypothesis["hypothesis_id"] in hypothesis["rival_hypothesis_ids"]:
            errors.append(
                issue("HYPOTHESIS_CANNOT_RIVAL_ITSELF", hypothesis["hypothesis_id"])
            )
        if len(hypothesis_ids) > 1 and not hypothesis["rival_hypothesis_ids"]:
            warnings.append(
                issue("RIVAL_HYPOTHESIS_LINK_MISSING", hypothesis["hypothesis_id"])
            )
    for operationalization in record["operationalizations"]:
        used_sources.update(operationalization["validity_evidence_source_ids"])
    errors.extend(
        _missing_references(
            sorted(used_sources),
            source_ids,
            "SOURCE_NOT_DECLARED_IN_EVIDENCE",
            "record.evidence.source_ids",
        )
    )

    if len(hypothesis_ids) == 1:
        warnings.append(issue("SINGLE_CANDIDATE_REQUIRES_RIVAL_REVIEW", "hypotheses"))

    for estimand in record["causal_estimands"]:
        errors.extend(
            _missing_references(
                estimand["linked_hypothesis_ids"],
                hypothesis_ids,
                "UNKNOWN_HYPOTHESIS_ID",
                estimand["estimand_id"],
            )
        )
    if (
        record["research_question"]["question_type"] == "causal"
        and not record["causal_estimands"]
    ):
        errors.append(issue("CAUSAL_QUESTION_REQUIRES_ESTIMAND", "causal_estimands"))

    predictions_by_hypothesis: Counter[str] = Counter()
    for prediction in record["predictions"]:
        prediction_id = prediction["prediction_id"]
        predictions_by_hypothesis[prediction["hypothesis_id"]] += 1
        errors.extend(
            _missing_references(
                [prediction["hypothesis_id"]],
                hypothesis_ids,
                "UNKNOWN_HYPOTHESIS_ID",
                prediction_id,
            )
        )
        errors.extend(
            _missing_references(
                prediction["rival_hypothesis_ids"],
                hypothesis_ids,
                "UNKNOWN_RIVAL_HYPOTHESIS_ID",
                prediction_id,
            )
        )
        errors.extend(
            _missing_references(
                prediction["measurement_ids"],
                measurement_ids,
                "UNKNOWN_MEASUREMENT_ID",
                prediction_id,
            )
        )
        errors.extend(
            _missing_references(
                prediction["analysis_ids"],
                analysis_ids,
                "UNKNOWN_ANALYSIS_ID",
                prediction_id,
            )
        )
        if prediction["hypothesis_id"] in prediction["rival_hypothesis_ids"]:
            errors.append(issue("PREDICTION_RIVAL_IS_FOCAL", prediction_id))
    for hypothesis_id in sorted(hypothesis_ids):
        if predictions_by_hypothesis[hypothesis_id] == 0:
            errors.append(issue("HYPOTHESIS_HAS_NO_PREDICTION", hypothesis_id))

    for alternative in record["alternative_explanations"]:
        errors.extend(
            _missing_references(
                alternative["linked_hypothesis_ids"],
                hypothesis_ids,
                "UNKNOWN_HYPOTHESIS_ID",
                alternative["alternative_id"],
            )
        )
        errors.extend(
            _missing_references(
                alternative["discriminating_prediction_ids"],
                prediction_ids,
                "UNKNOWN_PREDICTION_ID",
                alternative["alternative_id"],
            )
        )

    for null in record["null_hypotheses"]:
        errors.extend(
            _missing_references(
                null["linked_prediction_ids"],
                prediction_ids,
                "UNKNOWN_PREDICTION_ID",
                null["null_id"],
            )
        )

    for control in record["negative_controls"]:
        errors.extend(
            _missing_references(
                control["linked_prediction_ids"],
                prediction_ids,
                "UNKNOWN_PREDICTION_ID",
                control["control_id"],
            )
        )
    if not any(
        control["control_type"]
        in {"negative_exposure", "negative_outcome", "procedural_negative"}
        for control in record["negative_controls"]
    ):
        warnings.append(issue("NEGATIVE_CONTROL_TYPE_REVIEW_REQUIRED", "controls"))

    for analysis in record["analysis_plan"]["analyses"]:
        errors.extend(
            _missing_references(
                analysis["prediction_ids"],
                prediction_ids,
                "UNKNOWN_PREDICTION_ID",
                analysis["analysis_id"],
            )
        )
        errors.extend(
            _missing_references(
                analysis["estimand_ids"],
                estimand_ids,
                "UNKNOWN_ESTIMAND_ID",
                analysis["analysis_id"],
            )
        )

    for category, entries in record["risk_register"].items():
        if not entries:
            warnings.append(issue("RISK_CATEGORY_EMPTY_REQUIRES_RATIONALE", category))

    ethics = record["ethics_and_feasibility"]
    gate_values = {
        key: value for key, value in ethics.items() if key.endswith("_gate")
    }
    unresolved_gate_names = sorted(
        key
        for key, value in gate_values.items()
        if value in {"undetermined", "requires_review", "blocked"}
    )
    if unresolved_gate_names and not ethics["unresolved_blocks"]:
        errors.append(
            issue("UNRESOLVED_GATE_REQUIRES_BLOCK_RECORD", "ethics_and_feasibility")
        )
    if (
        any(value == "requires_review" for value in gate_values.values())
        and not ethics["required_reviews"]
    ):
        errors.append(
            issue("REQUIRED_REVIEW_LIST_MISSING", "ethics_and_feasibility")
        )
    if ethics["feasibility_status"] in {"undetermined", "infeasible", "blocked"}:
        unresolved_gate_names.append("feasibility_status")

    ai_use = record["ai_use"]
    if ai_use["sensitive_or_unpublished_data_sent_externally"]:
        errors.append(
            issue(
                "EXTERNAL_SENSITIVE_DATA_NOT_SUPPORTED",
                "ai_use.sensitive_or_unpublished_data_sent_externally",
            )
        )
    if not ai_use["local_policy_checked"]:
        errors.append(issue("LOCAL_AI_POLICY_NOT_CHECKED", "ai_use.local_policy_checked"))
    if not ai_use["citation_verification_required"]:
        errors.append(
            issue(
                "CITATION_VERIFICATION_MUST_BE_REQUIRED",
                "ai_use.citation_verification_required",
            )
        )
    if not ai_use["human_accountable"]:
        errors.append(
            issue("HUMAN_ACCOUNTABILITY_REQUIRED", "ai_use.human_accountable")
        )

    return {
        "schema_version": "2.0",
        "project_id": record["project_id"],
        "valid": not errors,
        "status": (
            "INVALID_RECORD"
            if errors
            else "VALID_BLOCKED_BY_GATES"
            if unresolved_gate_names
            else "VALID_FOR_HUMAN_REVIEW"
        ),
        "errors": errors,
        "warnings": warnings,
        "counts": {
            "hypotheses": len(hypothesis_ids),
            "causal_estimands": len(estimand_ids),
            "predictions": len(prediction_ids),
            "alternative_explanations": len(record["alternative_explanations"]),
            "null_hypotheses": len(record["null_hypotheses"]),
            "negative_controls": len(control_ids),
            "operationalizations": len(measurement_ids),
            "analyses": len(analysis_ids),
            "declared_sources": len(source_ids),
        },
        "identifiers": {
            "hypothesis_ids": sorted(hypothesis_ids),
            "estimand_ids": sorted(estimand_ids),
            "prediction_ids": sorted(prediction_ids),
            "measurement_ids": sorted(measurement_ids),
            "analysis_ids": sorted(analysis_ids),
            "control_ids": sorted(control_ids),
        },
        "unresolved_gate_fields": sorted(set(unresolved_gate_names)),
        "notice": (
            "This report validates schema, cross-references, and declared safety "
            "controls only. It does not verify evidence, scientific validity, "
            "novelty, ethics approval, causal identification, or hypothesis merit, "
            "and it does not rank or select candidates."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate a bounded local hypothesis JSON record and emit identifiers, "
            "counts, and rule codes without scientific scoring."
        )
    )
    parser.add_argument("record", help="Local hypothesis record JSON")
    parser.add_argument("-o", "--output", help="Optional local JSON report")
    parser.add_argument(
        "--force", action="store_true", help="Replace an existing output file"
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        record = load_hypothesis_record(read_json(args.record))
        report = validate_record(record)
        write_json_report(report, args.output, force=args.force)
        return 0 if report["valid"] else 1
    except ValidationError as exc:
        return error_exit(exc)


if __name__ == "__main__":
    raise SystemExit(main())
