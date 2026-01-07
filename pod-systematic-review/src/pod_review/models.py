"""Data models for systematic review."""

from datetime import datetime
from enum import Enum
from typing import Any

from pydantic import BaseModel, Field


class EffectType(str, Enum):
    """Type of effect size measure."""

    OR = "OR"  # Odds ratio
    RR = "RR"  # Risk ratio
    HR = "HR"  # Hazard ratio


class StudyDesign(str, Enum):
    """Study design types."""

    PROSPECTIVE_COHORT = "prospective_cohort"
    RETROSPECTIVE_COHORT = "retrospective_cohort"
    CASE_CONTROL = "case_control"
    RCT_SECONDARY = "rct_secondary"


class ScreeningDecision(str, Enum):
    """Screening decision options."""

    INCLUDE = "include"
    EXCLUDE = "exclude"
    UNCLEAR = "unclear"


class ScreeningPhase(str, Enum):
    """Screening phase."""

    TITLE_ABSTRACT = "title_abstract"
    FULLTEXT = "fulltext"


class Citation(BaseModel):
    """Citation/study record."""

    id: str = Field(description="Unique identifier (PMID, DOI, or generated)")
    title: str
    authors: list[str] = Field(default_factory=list)
    year: int | None = None
    journal: str = ""
    volume: str = ""
    issue: str = ""
    pages: str = ""
    abstract: str = ""
    doi: str = ""
    pmid: str = ""
    pmcid: str = ""
    keywords: list[str] = Field(default_factory=list)
    mesh_terms: list[str] = Field(default_factory=list)
    source_database: str = ""
    import_date: datetime = Field(default_factory=datetime.now)
    url: str = ""
    language: str = ""
    metadata: dict[str, Any] = Field(default_factory=dict)

    def get_identifier(self) -> str:
        """Get best available identifier."""
        return self.doi or self.pmid or self.id

    def get_author_string(self) -> str:
        """Get formatted author string."""
        if not self.authors:
            return ""
        if len(self.authors) == 1:
            return self.authors[0]
        elif len(self.authors) <= 3:
            return ", ".join(self.authors)
        else:
            return f"{self.authors[0]} et al."

    def to_citation_text(self) -> str:
        """Generate citation text."""
        parts = []
        if self.authors:
            parts.append(self.get_author_string())
        if self.year:
            parts.append(f"({self.year})")
        if self.title:
            parts.append(self.title.rstrip(".") + ".")
        if self.journal:
            parts.append(f"{self.journal}.")
        if self.volume:
            vol_str = self.volume
            if self.issue:
                vol_str += f"({self.issue})"
            if self.pages:
                vol_str += f":{self.pages}"
            parts.append(vol_str + ".")
        if self.doi:
            parts.append(f"doi:{self.doi}")

        return " ".join(parts)


class ScreeningRecord(BaseModel):
    """Record of screening decision."""

    citation_id: str
    reviewer_id: str
    phase: ScreeningPhase
    decision: ScreeningDecision
    exclusion_reason: str = ""
    notes: str = ""
    timestamp: datetime = Field(default_factory=datetime.now)


class RiskFactor(BaseModel):
    """Risk factor and effect size."""

    name: str
    category: str = ""  # demographic, comorbidity, medication, surgery, anesthesia, postop
    comparison: str = ""
    effect_type: EffectType
    effect_size: float
    ci_lower: float | None = None
    ci_upper: float | None = None
    p_value: float | None = None
    adjustment: str = "crude"  # crude or adjusted
    covariates: list[str] = Field(default_factory=list)
    subgroup: str = ""

    def get_se_from_ci(self) -> float | None:
        """
        Calculate standard error from confidence interval.

        Assumes 95% CI and log-transformed effect size for OR/RR/HR.
        SE = (log(upper) - log(lower)) / (2 * 1.96)
        """
        import math

        if self.ci_lower is None or self.ci_upper is None:
            return None

        try:
            log_upper = math.log(self.ci_upper)
            log_lower = math.log(self.ci_lower)
            se = (log_upper - log_lower) / (2 * 1.96)
            return se
        except (ValueError, ZeroDivisionError):
            return None

    def get_log_effect(self) -> float:
        """Get log-transformed effect size."""
        import math

        return math.log(self.effect_size)


class ExtractionRecord(BaseModel):
    """Data extraction record for a study."""

    citation_id: str
    extractor_id: str

    # Study characteristics
    study_design: StudyDesign
    country: str = ""
    setting: str = ""
    recruitment_period: str = ""
    sample_size: int | None = None
    mean_age: float | None = None
    sd_age: float | None = None
    age_range: str = ""
    female_n: int | None = None
    female_percent: float | None = None
    surgery_type: str = ""  # THA, TKA, mixed
    primary_revision: str = ""  # primary, revision, mixed

    # POD outcome
    pod_definition: str = ""
    pod_assessment_tool: str = ""
    pod_assessment_timing: str = ""
    pod_incidence_n: int | None = None
    pod_incidence_percent: float | None = None

    # Risk factors
    risk_factors: list[RiskFactor] = Field(default_factory=list)

    # Metadata
    extraction_date: datetime = Field(default_factory=datetime.now)
    notes: str = ""


class NOSAssessment(BaseModel):
    """Newcastle-Ottawa Scale assessment."""

    citation_id: str
    assessor_id: str

    # Selection (max 4 stars)
    representativeness: int = Field(ge=0, le=1)  # 0 or 1 star
    selection_nonexposed: int = Field(ge=0, le=1)
    ascertainment_exposure: int = Field(ge=0, le=1)
    outcome_not_present: int = Field(ge=0, le=1)

    # Comparability (max 2 stars)
    comparability: int = Field(ge=0, le=2)

    # Outcome (max 3 stars)
    assessment_outcome: int = Field(ge=0, le=1)
    followup_duration: int = Field(ge=0, le=1)
    followup_completeness: int = Field(ge=0, le=1)

    notes: str = ""
    assessment_date: datetime = Field(default_factory=datetime.now)

    @property
    def total_stars(self) -> int:
        """Calculate total NOS stars (0-9)."""
        return (
            self.representativeness
            + self.selection_nonexposed
            + self.ascertainment_exposure
            + self.outcome_not_present
            + self.comparability
            + self.assessment_outcome
            + self.followup_duration
            + self.followup_completeness
        )

    @property
    def quality_rating(self) -> str:
        """Get quality rating based on stars."""
        stars = self.total_stars
        if stars >= 7:
            return "Good"
        elif stars >= 4:
            return "Fair"
        else:
            return "Poor"


class RoB2Assessment(BaseModel):
    """Risk of Bias 2 assessment for RCTs."""

    citation_id: str
    assessor_id: str

    # Domains (Low/Some concerns/High)
    randomization: str = ""
    deviations: str = ""
    missing_data: str = ""
    measurement: str = ""
    selection_reported: str = ""

    overall: str = ""  # Overall judgment
    notes: str = ""
    assessment_date: datetime = Field(default_factory=datetime.now)


class MetaAnalysisResult(BaseModel):
    """Results from meta-analysis."""

    risk_factor: str
    n_studies: int
    pooled_effect: float  # Pooled log OR/RR/HR
    pooled_ci_lower: float
    pooled_ci_upper: float
    pooled_p_value: float | None = None

    # Heterogeneity
    q_statistic: float
    q_p_value: float
    i_squared: float
    tau_squared: float

    # Publication bias (if applicable)
    egger_p_value: float | None = None

    method: str = "random_effects"
    model_details: dict[str, Any] = Field(default_factory=dict)
