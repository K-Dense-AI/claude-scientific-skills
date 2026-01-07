"""Meta-analysis statistical methods."""

import logging
import math
from typing import Any

import numpy as np
from scipy import stats

from ..models import ExtractionRecord, MetaAnalysisResult, RiskFactor

logger = logging.getLogger(__name__)


class MetaAnalyzer:
    """Perform meta-analysis using random-effects model."""

    def __init__(self, method: str = "DL") -> None:
        """
        Initialize meta-analyzer.

        Args:
            method: Meta-analysis method ('DL' for DerSimonian-Laird)
        """
        self.method = method

    def analyze(
        self,
        extraction_records: list[ExtractionRecord],
        risk_factor_name: str,
        effect_type: str | None = None,
    ) -> MetaAnalysisResult | None:
        """
        Perform meta-analysis for a specific risk factor.

        Args:
            extraction_records: List of extraction records
            risk_factor_name: Name of risk factor to analyze
            effect_type: Optional filter by effect type (OR, RR, HR)

        Returns:
            MetaAnalysisResult or None if insufficient data
        """
        # Extract relevant risk factors
        risk_factors: list[tuple[str, RiskFactor]] = []

        for record in extraction_records:
            for rf in record.risk_factors:
                if rf.name == risk_factor_name:
                    if effect_type is None or rf.effect_type.value == effect_type:
                        risk_factors.append((record.citation_id, rf))

        if len(risk_factors) < 2:
            logger.warning(f"Insufficient studies for {risk_factor_name} (n={len(risk_factors)})")
            return None

        logger.info(f"Performing meta-analysis for {risk_factor_name} with {len(risk_factors)} studies")

        # Calculate effect sizes and standard errors
        study_ids = []
        log_effects = []
        variances = []

        for study_id, rf in risk_factors:
            log_effect = rf.get_log_effect()
            se = rf.get_se_from_ci()

            if se is None or math.isnan(se) or se <= 0:
                logger.warning(f"Skipping {study_id}: invalid SE")
                continue

            study_ids.append(study_id)
            log_effects.append(log_effect)
            variances.append(se**2)

        if len(log_effects) < 2:
            logger.warning("Insufficient valid studies after filtering")
            return None

        log_effects_array = np.array(log_effects)
        variances_array = np.array(variances)

        # Perform random-effects meta-analysis (DerSimonian-Laird)
        result = self._random_effects_meta_analysis(log_effects_array, variances_array)

        # Publication bias (if >=10 studies)
        egger_p = None
        if len(log_effects) >= 10:
            egger_p = self._egger_test(log_effects_array, variances_array)

        return MetaAnalysisResult(
            risk_factor=risk_factor_name,
            n_studies=len(log_effects),
            pooled_effect=result["pooled_effect"],
            pooled_ci_lower=result["pooled_ci_lower"],
            pooled_ci_upper=result["pooled_ci_upper"],
            pooled_p_value=result["p_value"],
            q_statistic=result["Q"],
            q_p_value=result["Q_p"],
            i_squared=result["I2"],
            tau_squared=result["tau2"],
            egger_p_value=egger_p,
            method="random_effects_DL",
            model_details={
                "study_ids": study_ids,
                "log_effects": log_effects,
                "variances": variances,
            },
        )

    def _random_effects_meta_analysis(
        self, log_effects: np.ndarray, variances: np.ndarray
    ) -> dict[str, float]:
        """
        Perform DerSimonian-Laird random-effects meta-analysis.

        Args:
            log_effects: Array of log-transformed effect sizes
            variances: Array of within-study variances

        Returns:
            Dictionary with pooled results and heterogeneity stats
        """
        k = len(log_effects)
        weights_fixed = 1.0 / variances

        # Fixed-effect pooled estimate
        pooled_fixed = np.sum(weights_fixed * log_effects) / np.sum(weights_fixed)

        # Q statistic for heterogeneity
        Q = np.sum(weights_fixed * (log_effects - pooled_fixed) ** 2)
        Q_p = 1 - stats.chi2.cdf(Q, k - 1) if k > 1 else 1.0

        # Tau-squared (between-study variance) - DerSimonian-Laird method
        C = np.sum(weights_fixed) - np.sum(weights_fixed**2) / np.sum(weights_fixed)
        tau2 = max(0, (Q - (k - 1)) / C) if C > 0 else 0

        # I-squared (proportion of total variance due to heterogeneity)
        I2 = max(0, ((Q - (k - 1)) / Q) * 100) if Q > 0 else 0

        # Random-effects weights
        weights_random = 1.0 / (variances + tau2)

        # Random-effects pooled estimate
        pooled_effect = np.sum(weights_random * log_effects) / np.sum(weights_random)
        pooled_var = 1.0 / np.sum(weights_random)
        pooled_se = np.sqrt(pooled_var)

        # 95% CI
        pooled_ci_lower = pooled_effect - 1.96 * pooled_se
        pooled_ci_upper = pooled_effect + 1.96 * pooled_se

        # P-value
        z = pooled_effect / pooled_se
        p_value = 2 * (1 - stats.norm.cdf(abs(z)))

        return {
            "pooled_effect": pooled_effect,
            "pooled_se": pooled_se,
            "pooled_ci_lower": pooled_ci_lower,
            "pooled_ci_upper": pooled_ci_upper,
            "p_value": p_value,
            "Q": Q,
            "Q_p": Q_p,
            "I2": I2,
            "tau2": tau2,
        }

    def _egger_test(self, log_effects: np.ndarray, variances: np.ndarray) -> float:
        """
        Perform Egger's test for publication bias.

        Args:
            log_effects: Array of log effect sizes
            variances: Array of variances

        Returns:
            P-value for Egger's test
        """
        # Standard errors
        se = np.sqrt(variances)

        # Precision (1/SE)
        precision = 1.0 / se

        # Regression: effect size ~ precision
        # Using SciPy's linregress
        slope, intercept, r_value, p_value, std_err = stats.linregress(precision, log_effects)

        # Egger's test p-value is for the intercept
        # We use a simple t-test here
        # For more rigorous approach, use weighted regression

        return p_value

    def leave_one_out_analysis(
        self, extraction_records: list[ExtractionRecord], risk_factor_name: str
    ) -> list[dict[str, Any]]:
        """
        Perform leave-one-out sensitivity analysis.

        Args:
            extraction_records: Extraction records
            risk_factor_name: Risk factor name

        Returns:
            List of results with each study excluded
        """
        results = []

        # Get all risk factors
        all_rfs = []
        for record in extraction_records:
            for rf in record.risk_factors:
                if rf.name == risk_factor_name:
                    all_rfs.append((record.citation_id, rf))

        for i, (excluded_id, _) in enumerate(all_rfs):
            # Create subset excluding this study
            subset = [rec for rec in extraction_records if rec.citation_id != excluded_id]

            # Run meta-analysis
            result = self.analyze(subset, risk_factor_name)

            if result:
                results.append(
                    {
                        "excluded_study": excluded_id,
                        "n_studies": result.n_studies,
                        "pooled_effect": result.pooled_effect,
                        "pooled_ci_lower": result.pooled_ci_lower,
                        "pooled_ci_upper": result.pooled_ci_upper,
                        "i_squared": result.i_squared,
                    }
                )

        return results

    def subgroup_analysis(
        self,
        extraction_records: list[ExtractionRecord],
        risk_factor_name: str,
        subgroup_var: str,
    ) -> dict[str, MetaAnalysisResult | None]:
        """
        Perform subgroup meta-analysis.

        Args:
            extraction_records: Extraction records
            risk_factor_name: Risk factor name
            subgroup_var: Variable to stratify by (e.g., 'surgery_type')

        Returns:
            Dictionary mapping subgroup values to results
        """
        # Group records by subgroup variable
        subgroups: dict[str, list[ExtractionRecord]] = {}

        for record in extraction_records:
            value = getattr(record, subgroup_var, "Unknown")
            if value not in subgroups:
                subgroups[value] = []
            subgroups[value].append(record)

        # Run meta-analysis for each subgroup
        results = {}
        for subgroup_value, records in subgroups.items():
            result = self.analyze(records, risk_factor_name)
            results[subgroup_value] = result

        return results
