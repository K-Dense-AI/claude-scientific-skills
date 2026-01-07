"""Tests for data models."""

import math
import pytest
from pod_review.models import RiskFactor, EffectType, NOSAssessment


def test_risk_factor_se_calculation() -> None:
    """Test standard error calculation from confidence interval."""
    rf = RiskFactor(
        name="Age",
        effect_type=EffectType.OR,
        effect_size=1.5,
        ci_lower=1.2,
        ci_upper=1.8,
    )

    se = rf.get_se_from_ci()

    assert se is not None
    assert se > 0

    # Verify calculation
    expected_se = (math.log(1.8) - math.log(1.2)) / (2 * 1.96)
    assert abs(se - expected_se) < 0.001


def test_risk_factor_log_effect() -> None:
    """Test log effect size calculation."""
    rf = RiskFactor(
        name="Age",
        effect_type=EffectType.OR,
        effect_size=2.0,
        ci_lower=1.5,
        ci_upper=2.5,
    )

    log_effect = rf.get_log_effect()

    assert abs(log_effect - math.log(2.0)) < 0.001


def test_nos_total_stars() -> None:
    """Test NOS total stars calculation."""
    assessment = NOSAssessment(
        citation_id="test",
        assessor_id="R1",
        representativeness=1,
        selection_nonexposed=1,
        ascertainment_exposure=1,
        outcome_not_present=1,
        comparability=2,
        assessment_outcome=1,
        followup_duration=1,
        followup_completeness=1,
    )

    assert assessment.total_stars == 9
    assert assessment.quality_rating == "Good"


def test_nos_quality_rating() -> None:
    """Test NOS quality rating."""
    # Good quality (7-9 stars)
    good = NOSAssessment(
        citation_id="test",
        assessor_id="R1",
        representativeness=1,
        selection_nonexposed=1,
        ascertainment_exposure=1,
        outcome_not_present=1,
        comparability=2,
        assessment_outcome=1,
        followup_duration=0,
        followup_completeness=0,
    )
    assert good.quality_rating == "Good"

    # Fair quality (4-6 stars)
    fair = NOSAssessment(
        citation_id="test",
        assessor_id="R1",
        representativeness=1,
        selection_nonexposed=1,
        ascertainment_exposure=1,
        outcome_not_present=0,
        comparability=1,
        assessment_outcome=1,
        followup_duration=0,
        followup_completeness=0,
    )
    assert fair.quality_rating == "Fair"

    # Poor quality (0-3 stars)
    poor = NOSAssessment(
        citation_id="test",
        assessor_id="R1",
        representativeness=1,
        selection_nonexposed=1,
        ascertainment_exposure=0,
        outcome_not_present=0,
        comparability=0,
        assessment_outcome=0,
        followup_duration=0,
        followup_completeness=0,
    )
    assert poor.quality_rating == "Poor"
