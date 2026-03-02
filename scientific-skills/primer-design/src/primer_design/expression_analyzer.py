#!/usr/bin/env python3
"""
Expression Analyzer
====================
CDS 서열의 E. coli 발현 최적화 분석.

분석 항목:
  - GC 함량, 분자량, CAI (Codon Adaptation Index)
  - E. coli K12 rare codon 빈도 및 클러스터 검출
  - Signal peptide 예측 (rule-based)
  - N-terminal MAP (Methionine Aminopeptidase) 제거 예측
  - 발현 숙주 균주 추천 (BL21 / Rosetta / Rosetta 2)
"""

from __future__ import annotations

import math
import re

# ── Genetic code (standard) ─────────────────────────────────────────────────

_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _translate(dna: str) -> str:
    """DNA -> amino acid string (standard code). Stop = '*'."""
    aa = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3].upper()
        aa.append(_CODON_TABLE.get(codon, "?"))
    return "".join(aa)


# ── E. coli K12 Rare Codons ─────────────────────────────────────────────────

ECOLI_RARE_CODONS = {
    "AGG": {"aa": "Arg", "freq_per_1000": 1.2,  "prare": True,  "prare2": True},
    "AGA": {"aa": "Arg", "freq_per_1000": 2.1,  "prare": True,  "prare2": True},
    "CGA": {"aa": "Arg", "freq_per_1000": 3.1,  "prare": False, "prare2": True},
    "CUA": {"aa": "Leu", "freq_per_1000": 3.9,  "prare": True,  "prare2": True},
    "AUA": {"aa": "Ile", "freq_per_1000": 4.1,  "prare": True,  "prare2": True},
    "CCC": {"aa": "Pro", "freq_per_1000": 4.3,  "prare": True,  "prare2": True},
    "GGA": {"aa": "Gly", "freq_per_1000": 8.0,  "prare": True,  "prare2": True},
    "CGG": {"aa": "Arg", "freq_per_1000": 5.4,  "prare": False, "prare2": True},
}

# ── E. coli K12 Codon Usage (relative adaptiveness w_i for CAI) ─────────────
# Source: Kazusa codon usage database — E. coli K12
# w_i = freq / max_freq for that amino acid (pre-computed)

ECOLI_CODON_W = {
    # Phe
    "TTT": 0.58, "TTC": 1.00,
    # Leu
    "TTA": 0.13, "TTG": 0.13, "CTT": 0.10, "CTC": 0.10, "CTA": 0.04, "CTG": 1.00,
    # Ile
    "ATT": 0.49, "ATC": 1.00, "ATA": 0.07,
    # Met
    "ATG": 1.00,
    # Val
    "GTT": 1.00, "GTC": 0.42, "GTA": 0.50, "GTG": 0.37,
    # Ser
    "TCT": 1.00, "TCC": 0.74, "TCA": 0.12, "TCG": 0.15, "AGT": 0.16, "AGC": 0.76,
    # Pro
    "CCT": 0.16, "CCC": 0.12, "CCA": 0.19, "CCG": 1.00,
    # Thr
    "ACT": 0.50, "ACC": 1.00, "ACA": 0.14, "ACG": 0.36,
    # Ala
    "GCT": 1.00, "GCC": 0.63, "GCA": 0.59, "GCG": 0.50,
    # Tyr
    "TAT": 0.59, "TAC": 1.00,
    # Stop
    "TAA": 1.00, "TAG": 0.08, "TGA": 0.30,
    # His
    "CAT": 0.57, "CAC": 1.00,
    # Gln
    "CAA": 0.34, "CAG": 1.00,
    # Asn
    "AAT": 0.49, "AAC": 1.00,
    # Lys
    "AAA": 0.77, "AAG": 1.00,
    # Asp
    "GAT": 0.63, "GAC": 1.00,
    # Glu
    "GAA": 1.00, "GAG": 0.29,
    # Cys
    "TGT": 0.46, "TGC": 1.00,
    # Trp
    "TGG": 1.00,
    # Arg
    "CGT": 1.00, "CGC": 0.60, "CGA": 0.07, "CGG": 0.10, "AGA": 0.04, "AGG": 0.02,
    # Gly
    "GGT": 1.00, "GGC": 0.77, "GGA": 0.11, "GGG": 0.15,
}

# ── Amino acid molecular weights (average) ──────────────────────────────────

_AA_MW = {
    "A": 89.1, "R": 174.2, "N": 132.1, "D": 133.1, "C": 121.2,
    "E": 147.1, "Q": 146.2, "G": 75.0, "H": 155.2, "I": 131.2,
    "L": 131.2, "K": 146.2, "M": 149.2, "F": 165.2, "P": 115.1,
    "S": 105.1, "T": 119.1, "W": 204.2, "Y": 181.2, "V": 117.1,
}
_WATER_MW = 18.02

# ── MAP (Methionine Aminopeptidase) removal rules ───────────────────────────
# E. coli MAP removes N-terminal Met when the second amino acid is small
# (radius of gyration <= 1.29 A)

MAP_REMOVABLE_AA = {"A", "C", "G", "P", "S", "T", "V"}


# ── ExpressionAnalyzer ──────────────────────────────────────────────────────

class ExpressionAnalyzer:
    """CDS 서열의 E. coli 발현 최적화 분석기."""

    def analyze(
        self,
        cds_seq: str,
        organism_source: str = "",
    ) -> dict:
        """CDS 서열 종합 분석.

        Parameters
        ----------
        cds_seq : str
            CDS DNA 서열 (ATG ~ stop codon 포함)
        organism_source : str
            원본 생물 (예: "Agrobacterium tumefaciens"). 정보성 메모 용도.

        Returns
        -------
        dict
            분석 결과 (basic_info, rare_codons, clusters, cai, signal_peptide 등)
        """
        warnings: list[str] = []
        recommendations: list[str] = []

        # 1. Clean input
        seq = cds_seq.upper().replace(" ", "").replace("\n", "").replace("\r", "")
        invalid = set(seq) - {"A", "T", "G", "C"}
        if invalid:
            warnings.append(
                f"Non-standard characters in CDS: {', '.join(sorted(invalid))}"
            )
            seq = re.sub(r"[^ATGC]", "", seq)

        # 2. Start codon check
        has_start = seq[:3] == "ATG"
        if not has_start:
            warnings.append(
                f"CDS does not start with ATG (found: {seq[:3]})"
            )

        # 3. Translate to protein
        protein = _translate(seq)
        # Strip terminal stop
        if protein.endswith("*"):
            protein_no_stop = protein[:-1]
        else:
            protein_no_stop = protein
            warnings.append("CDS does not end with a stop codon")

        # 4. Basic info
        cds_length_bp = len(seq)
        protein_length_aa = len(protein_no_stop)
        gc_count = seq.count("G") + seq.count("C")
        gc_content = gc_count / cds_length_bp if cds_length_bp > 0 else 0.0

        # Molecular weight (sum of aa MW - (n-1)*water)
        mw_sum = sum(_AA_MW.get(aa, 0.0) for aa in protein_no_stop)
        if protein_length_aa > 1:
            mw_sum -= (protein_length_aa - 1) * _WATER_MW
        molecular_weight_kda = mw_sum / 1000.0

        # 5. Rare codon analysis
        codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
        total_codons = len(codons)

        rare_codon_counts: dict[str, int] = {}
        rare_codon_positions: dict[str, list[int]] = {}
        total_rare = 0
        prare_count = 0
        prare2_count = 0

        for idx, codon in enumerate(codons):
            if codon in ECOLI_RARE_CODONS:
                info = ECOLI_RARE_CODONS[codon]
                total_rare += 1
                rare_codon_counts[codon] = rare_codon_counts.get(codon, 0) + 1
                if codon not in rare_codon_positions:
                    rare_codon_positions[codon] = []
                rare_codon_positions[codon].append(idx + 1)  # 1-indexed
                if info["prare"]:
                    prare_count += 1
                if info["prare2"]:
                    prare2_count += 1

        rare_codon_pct = (total_rare / total_codons * 100) if total_codons > 0 else 0.0

        rare_codon_detail = []
        for codon in sorted(rare_codon_counts.keys()):
            info = ECOLI_RARE_CODONS[codon]
            rare_codon_detail.append({
                "codon": codon,
                "amino_acid": info["aa"],
                "count": rare_codon_counts[codon],
                "freq_per_1000": info["freq_per_1000"],
                "positions": rare_codon_positions[codon],
                "prare": info["prare"],
                "prare2": info["prare2"],
            })

        # 6. Rare codon cluster detection (sliding window of 5 codons)
        cluster_windows: list[tuple[int, int]] = []
        for i in range(len(codons) - 4):
            window = codons[i:i + 5]
            rare_in_window = sum(1 for c in window if c in ECOLI_RARE_CODONS)
            if rare_in_window >= 2:
                cluster_windows.append((i + 1, i + 5))  # 1-indexed

        # Merge overlapping clusters
        merged_clusters: list[dict] = []
        if cluster_windows:
            current_start, current_end = cluster_windows[0]
            for start, end in cluster_windows[1:]:
                if start <= current_end:
                    current_end = max(current_end, end)
                else:
                    # Collect rare codons in this merged cluster
                    cluster_rare = []
                    for pos_idx in range(current_start - 1, current_end):
                        if pos_idx < len(codons) and codons[pos_idx] in ECOLI_RARE_CODONS:
                            cluster_rare.append({
                                "codon": codons[pos_idx],
                                "position": pos_idx + 1,
                            })
                    merged_clusters.append({
                        "start": current_start,
                        "end": current_end,
                        "rare_codons": cluster_rare,
                    })
                    current_start, current_end = start, end
            # Final cluster
            cluster_rare = []
            for pos_idx in range(current_start - 1, current_end):
                if pos_idx < len(codons) and codons[pos_idx] in ECOLI_RARE_CODONS:
                    cluster_rare.append({
                        "codon": codons[pos_idx],
                        "position": pos_idx + 1,
                    })
            merged_clusters.append({
                "start": current_start,
                "end": current_end,
                "rare_codons": cluster_rare,
            })

        # 7. CAI calculation
        cai = self.compute_cai(seq)

        # 8. Signal peptide check
        signal_peptide = self.check_signal_peptide(protein_no_stop)

        # 9. MAP removal prediction
        map_removal = False
        map_note = ""
        if protein_length_aa >= 2:
            second_aa = protein_no_stop[1]
            if second_aa in MAP_REMOVABLE_AA:
                map_removal = True
                map_note = (
                    f"N-terminal Met likely removed by MAP "
                    f"(second residue: {second_aa}). "
                    f"Mature protein starts at {second_aa} (position 2)."
                )
            else:
                map_note = (
                    f"N-terminal Met likely retained "
                    f"(second residue: {second_aa}, not small enough for MAP)."
                )

        # Adjusted MW for MAP removal
        mature_mw_kda = molecular_weight_kda
        if map_removal:
            mature_mw_kda = (mw_sum - _AA_MW.get("M", 0.0) + _WATER_MW) / 1000.0

        # 10. Recommendations and warnings
        if rare_codon_pct > 8:
            recommendations.append(
                "High rare codon content (>8%): strongly recommend codon optimization "
                "or use Rosetta 2(DE3)"
            )
        elif rare_codon_pct > 3:
            recommendations.append(
                "Moderate rare codon content (3-8%): Rosetta(DE3) recommended"
            )

        if merged_clusters:
            recommendations.append(
                f"{len(merged_clusters)} rare codon cluster(s) detected: "
                f"may cause ribosomal stalling. Consider codon optimization "
                f"in clustered regions."
            )

        if cai < 0.2:
            recommendations.append(
                "Very low CAI (<0.2): expression likely very poor. "
                "Codon optimization strongly recommended."
            )
        elif cai < 0.4:
            recommendations.append(
                "Low CAI (<0.4): expression may be poor. "
                "Consider codon optimization."
            )

        if signal_peptide["has_signal_peptide"]:
            recommendations.append(
                "Potential signal peptide detected: may cause secretion or "
                "membrane targeting. Consider removing for cytoplasmic expression."
            )

        if gc_content > 0.65:
            warnings.append(
                f"High GC content ({gc_content:.1%}): may form secondary structures "
                f"affecting translation"
            )
        elif gc_content < 0.35:
            warnings.append(
                f"Low GC content ({gc_content:.1%}): mRNA may be unstable in E. coli"
            )

        # Strain recommendation
        strain_rec = self.recommend_strain({
            "rare_codon_pct": rare_codon_pct,
            "rare_codon_clusters": merged_clusters,
            "cai": cai,
            "signal_peptide": signal_peptide,
        })

        return {
            "organism_source": organism_source,
            "basic_info": {
                "cds_length_bp": cds_length_bp,
                "protein_length_aa": protein_length_aa,
                "gc_content": round(gc_content, 4),
                "molecular_weight_kda": round(molecular_weight_kda, 2),
                "mature_mw_kda": round(mature_mw_kda, 2),
                "has_start_codon": has_start,
                "protein_sequence": protein_no_stop,
            },
            "rare_codons": {
                "total_rare": total_rare,
                "total_codons": total_codons,
                "rare_codon_pct": round(rare_codon_pct, 2),
                "prare_count": prare_count,
                "prare2_count": prare2_count,
                "detail": rare_codon_detail,
            },
            "rare_codon_clusters": merged_clusters,
            "cai": round(cai, 4),
            "signal_peptide": signal_peptide,
            "map_removal": {
                "will_be_removed": map_removal,
                "note": map_note,
            },
            "strain_recommendation": strain_rec,
            "recommendations": recommendations,
            "warnings": warnings,
        }

    def recommend_strain(self, analysis: dict) -> dict:
        """발현 숙주 균주 추천.

        Parameters
        ----------
        analysis : dict
            analyze() 결과 또는 rare_codon_pct, cai, rare_codon_clusters,
            signal_peptide 키를 포함하는 dict.

        Returns
        -------
        dict
            primary_strain, alternative_strain, codon_optimization_needed,
            optimization_priority, rationale
        """
        rare_pct = analysis.get("rare_codon_pct", 0.0)
        clusters = analysis.get("rare_codon_clusters", [])
        cai = analysis.get("cai", 1.0)
        signal = analysis.get("signal_peptide", {})

        rationale: list[str] = []
        primary = "BL21(DE3)"
        alternative: str | None = None
        codon_opt = False
        priority = "low"

        # Rare codon percentage thresholds
        if rare_pct > 8:
            primary = "Rosetta 2(DE3)"
            alternative = "Rosetta(DE3)"
            codon_opt = True
            priority = "high"
            rationale.append(
                f"Rare codon content ({rare_pct:.1f}%) exceeds 8%: "
                f"Rosetta 2(DE3) supplies all rare tRNAs"
            )
        elif rare_pct > 3:
            primary = "Rosetta(DE3)"
            alternative = "BL21(DE3)"
            priority = "medium"
            rationale.append(
                f"Moderate rare codon content ({rare_pct:.1f}%): "
                f"Rosetta(DE3) recommended for supplemental rare tRNAs"
            )
        else:
            rationale.append(
                f"Low rare codon content ({rare_pct:.1f}%): "
                f"BL21(DE3) sufficient"
            )

        # Rare codon clusters override
        if clusters:
            if primary == "BL21(DE3)":
                primary = "Rosetta(DE3)"
                alternative = "BL21(DE3)"
            if priority == "low":
                priority = "medium"
            codon_opt = True
            rationale.append(
                f"{len(clusters)} rare codon cluster(s) detected: "
                f"Rosetta strain mandatory or codon optimization required"
            )

        # CAI threshold
        if cai < 0.4:
            codon_opt = True
            if priority != "high":
                priority = "high"
            rationale.append(
                f"Low CAI ({cai:.3f}): codon optimization strongly recommended"
            )
        elif cai < 0.6:
            if not codon_opt:
                codon_opt = False
            rationale.append(
                f"Moderate CAI ({cai:.3f}): acceptable but optimization may improve yield"
            )

        # Signal peptide
        if signal.get("has_signal_peptide", False):
            rationale.append(
                "Signal peptide detected: consider removal for cytoplasmic expression"
            )

        return {
            "primary_strain": primary,
            "alternative_strain": alternative,
            "codon_optimization_needed": codon_opt,
            "optimization_priority": priority,
            "rationale": rationale,
        }

    def check_signal_peptide(self, protein_seq: str) -> dict:
        """Rule-based signal peptide detection.

        Parameters
        ----------
        protein_seq : str
            Amino acid sequence (single-letter, no stop)

        Returns
        -------
        dict
            has_signal_peptide, confidence, hydrophobic_pct,
            n_region_positive, cleavage_site_estimate, notes
        """
        notes: list[str] = []
        hydrophobic_aa = {"L", "I", "V", "F", "A", "W", "M"}
        positive_aa = {"K", "R"}
        small_aa = {"A", "G", "S", "T"}

        if len(protein_seq) < 15:
            return {
                "has_signal_peptide": False,
                "confidence": "none",
                "hydrophobic_pct": 0.0,
                "n_region_positive": 0,
                "cleavage_site_estimate": None,
                "notes": ["Protein too short for signal peptide analysis"],
            }

        # 1. Overall hydrophobic content in first 25 aa
        n_check = min(25, len(protein_seq))
        first_25 = protein_seq[:n_check]
        hydrophobic_count = sum(1 for aa in first_25 if aa in hydrophobic_aa)
        hydrophobic_pct = hydrophobic_count / n_check

        # 2. N-region: positively charged residues (K, R) in first 5 aa
        n_region = protein_seq[:min(5, len(protein_seq))]
        n_region_positive = sum(1 for aa in n_region if aa in positive_aa)

        # 3. H-region: hydrophobic core (positions 5-15)
        h_region_start = 5
        h_region_end = min(15, len(protein_seq))
        h_region = protein_seq[h_region_start:h_region_end]
        h_hydrophobic = sum(1 for aa in h_region if aa in hydrophobic_aa) if h_region else 0
        h_hydrophobic_pct = h_hydrophobic / len(h_region) if h_region else 0.0

        # 4. C-region: small amino acids near cleavage site (positions 15-25)
        cleavage_site: int | None = None
        c_region_start = min(15, len(protein_seq))
        c_region_end = min(25, len(protein_seq))
        c_region = protein_seq[c_region_start:c_region_end]

        # Look for small amino acids (A, G, S, T) that could be cleavage site
        for i in range(len(c_region) - 1, -1, -1):
            if c_region[i] in small_aa:
                cleavage_site = c_region_start + i + 1  # 1-indexed, after the small aa
                break

        # Scoring
        score = 0

        if hydrophobic_pct > 0.40:
            score += 2
            notes.append(
                f"High hydrophobic content in first {n_check} aa: "
                f"{hydrophobic_pct:.0%}"
            )
        elif hydrophobic_pct > 0.30:
            score += 1
            notes.append(
                f"Moderate hydrophobic content in first {n_check} aa: "
                f"{hydrophobic_pct:.0%}"
            )

        if n_region_positive >= 2:
            score += 2
            notes.append(
                f"Positively charged N-region: {n_region_positive} K/R in first 5 aa"
            )
        elif n_region_positive >= 1:
            score += 1
            notes.append(
                f"Weakly positive N-region: {n_region_positive} K/R in first 5 aa"
            )

        if h_hydrophobic_pct > 0.60:
            score += 2
            notes.append(
                f"Strong hydrophobic H-region (pos 5-15): {h_hydrophobic_pct:.0%}"
            )
        elif h_hydrophobic_pct > 0.40:
            score += 1
            notes.append(
                f"Moderate hydrophobic H-region (pos 5-15): {h_hydrophobic_pct:.0%}"
            )

        if cleavage_site is not None:
            score += 1
            notes.append(
                f"Potential cleavage site near position {cleavage_site} "
                f"(small amino acid)"
            )

        # Confidence assignment
        # Require both high overall hydrophobicity AND h-region for medium+
        has_n_and_h = (hydrophobic_pct > 0.40) and (h_hydrophobic_pct > 0.50)
        if score >= 5 and has_n_and_h:
            confidence = "high"
            has_signal = True
        elif score >= 4 and has_n_and_h:
            confidence = "medium"
            has_signal = True
        elif score >= 3:
            confidence = "low"
            has_signal = False
            notes.append("Some signal peptide features but below threshold")
        else:
            confidence = "none"
            has_signal = False

        return {
            "has_signal_peptide": has_signal,
            "confidence": confidence,
            "hydrophobic_pct": round(hydrophobic_pct, 4),
            "n_region_positive": n_region_positive,
            "cleavage_site_estimate": cleavage_site if has_signal else None,
            "notes": notes,
        }

    def compute_cai(self, cds_seq: str) -> float:
        """Codon Adaptation Index (CAI) for E. coli K12.

        CAI = exp( (1/n) * sum(ln(w_i)) )

        Parameters
        ----------
        cds_seq : str
            CDS DNA 서열

        Returns
        -------
        float
            CAI 값 (0 ~ 1). 1에 가까울수록 E. coli에 최적화됨.
        """
        seq = cds_seq.upper().replace(" ", "").replace("\n", "").replace("\r", "")
        codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]

        log_sum = 0.0
        n = 0
        for codon in codons:
            w = ECOLI_CODON_W.get(codon)
            if w is None:
                continue
            # Skip stop codons for CAI calculation
            if _CODON_TABLE.get(codon) == "*":
                continue
            if w <= 0:
                continue
            log_sum += math.log(w)
            n += 1

        if n == 0:
            return 0.0

        return math.exp(log_sum / n)


# ── Tests ────────────────────────────────────────────────────────────────────

def _run_tests():
    """ExpressionAnalyzer 테스트."""
    sep = "=" * 70
    passed = 0
    failed = 0

    def check(label: str, actual, expected, tolerance: float | None = None):
        nonlocal passed, failed
        if tolerance is not None:
            ok = abs(actual - expected) <= tolerance
        else:
            ok = actual == expected
        if ok:
            passed += 1
            print(f"  OK: {label}")
        else:
            failed += 1
            print(f"  FAIL: {label}")
            print(f"    expected: {expected!r}")
            print(f"    actual:   {actual!r}")

    analyzer = ExpressionAnalyzer()

    # ── Test 1: E. coli-optimized GFP analysis ───────────────────────────
    # CAI should be > 0.7, few rare codons
    ecoli_gfp = (
        "ATGAGTAAAGGCGAAGAACTGTTTACCGGCGTGGTGCCGATCCTGGTGGAACTGGATGGCGATGTGAAC"
        "GGCCACAAATTTAGCGTGCGCGGCGAAGGCGAAGGCGATGCGACCAACGGCAAACTGACCCTGAAATTTATC"
        "TGCACCACCGGCAAACTGCCGGTGCCGTGGCCGACCCTGGTGACCACCCTGACCTATGGCGTGCAGTGCTTC"
        "GCGCGCTATCCGGATCACATGAAACGCCACGATTTCTTTAAATCCGCGATGCCGGAAGGCTATGTGCAGGAA"
        "CGCACCATCTTTTTCAAAGACGATGGCACCTACAAAACCCGCGCGGAAGTGAAATTTGAAGGCGATACCCTGG"
        "TGAACCGCATCGAACTGAAAGGCATCGATTTTAAAGAAGATGGCAACATCCTGGGCCACAAACTGGAATACAA"
        "CTTTAACTCCCACAACGTGTACATCACCGCGGACAAACAAAAAAACGGCATCAAAGCGAACTTCAAAATCCGCC"
        "ACAACGTGGAAGATGGCAGCGTGCAGCTGGCGGATCATTATCAGCAGAACACCCCGATCGGCGATGGCCCGGTG"
        "CTGCTGCCGGATAACCACTACCTGTCGACCCAGAGCAAACTGTCGAAAGATCCGAACGAAAAACGCGATCAC"
        "ATGGTGCTGCTGGAATTTGTGACCGCGGCGGGCATCACCCTGGGCATGGATGAACTGTATAAATAA"
    )

    print(sep)
    print("Test 1: E. coli-optimized GFP -- high CAI, few rare codons")
    print(sep)
    result_gfp = analyzer.analyze(ecoli_gfp, organism_source="synthetic (E. coli optimized)")
    print(f"  CDS length: {result_gfp['basic_info']['cds_length_bp']} bp")
    print(f"  Protein length: {result_gfp['basic_info']['protein_length_aa']} aa")
    print(f"  GC content: {result_gfp['basic_info']['gc_content']:.1%}")
    print(f"  MW: {result_gfp['basic_info']['molecular_weight_kda']:.2f} kDa")
    print(f"  CAI: {result_gfp['cai']:.4f}")
    print(f"  Rare codons: {result_gfp['rare_codons']['total_rare']}"
          f" / {result_gfp['rare_codons']['total_codons']}"
          f" ({result_gfp['rare_codons']['rare_codon_pct']:.1f}%)")
    print(f"  Strain: {result_gfp['strain_recommendation']['primary_strain']}")
    check("CAI > 0.7", result_gfp["cai"] > 0.7, True)
    check("rare codon pct < 5%", result_gfp["rare_codons"]["rare_codon_pct"] < 5.0, True)
    check("has start codon", result_gfp["basic_info"]["has_start_codon"], True)
    print()

    # ── Test 2: Plant-origin CDS analysis ─────────────────────────────────
    # Should have many rare codons for E. coli
    plant_cds = (
        "ATGAGAAGATCTTCTTCAAGAAGAATCGCTTCTTCAAGATCCTCAGGATCCTCATCCTCAAGAGGATCTTCA"
        "TCGGATTCTCAAGAGGAAGAATCTCATCCTCAAGATCCGGACGATCTTCAAGATCCTCAGGATCCTCATCCTC"
        "AAGAGGATCTTCATCGGATTCTCAAGATAA"
    )

    print(sep)
    print("Test 2: Plant-origin CDS -- many rare codons expected")
    print(sep)
    result_plant = analyzer.analyze(plant_cds, organism_source="Arabidopsis thaliana")
    print(f"  CDS length: {result_plant['basic_info']['cds_length_bp']} bp")
    print(f"  Protein length: {result_plant['basic_info']['protein_length_aa']} aa")
    print(f"  CAI: {result_plant['cai']:.4f}")
    print(f"  Rare codons: {result_plant['rare_codons']['total_rare']}"
          f" / {result_plant['rare_codons']['total_codons']}"
          f" ({result_plant['rare_codons']['rare_codon_pct']:.1f}%)")
    if result_plant["rare_codons"]["detail"]:
        print("  Rare codon breakdown:")
        for d in result_plant["rare_codons"]["detail"]:
            print(f"    {d['codon']} ({d['amino_acid']}): {d['count']}x "
                  f"at positions {d['positions'][:5]}{'...' if len(d['positions']) > 5 else ''}")
    print(f"  Clusters: {len(result_plant['rare_codon_clusters'])}")
    print(f"  Strain: {result_plant['strain_recommendation']['primary_strain']}")
    check("has rare codons", result_plant["rare_codons"]["total_rare"] > 0, True)
    print()

    # ── Test 3: Strain recommendation ─────────────────────────────────────
    print(sep)
    print("Test 3: Strain recommendation logic")
    print(sep)

    # Low rare codons -> BL21
    rec_low = analyzer.recommend_strain({
        "rare_codon_pct": 1.5,
        "rare_codon_clusters": [],
        "cai": 0.75,
        "signal_peptide": {"has_signal_peptide": False},
    })
    check("low rare -> BL21(DE3)", rec_low["primary_strain"], "BL21(DE3)")
    check("low rare -> no codon opt", rec_low["codon_optimization_needed"], False)

    # Medium rare codons -> Rosetta
    rec_mid = analyzer.recommend_strain({
        "rare_codon_pct": 5.0,
        "rare_codon_clusters": [],
        "cai": 0.55,
        "signal_peptide": {"has_signal_peptide": False},
    })
    check("mid rare -> Rosetta(DE3)", rec_mid["primary_strain"], "Rosetta(DE3)")

    # High rare codons -> Rosetta 2
    rec_high = analyzer.recommend_strain({
        "rare_codon_pct": 12.0,
        "rare_codon_clusters": [],
        "cai": 0.30,
        "signal_peptide": {"has_signal_peptide": False},
    })
    check("high rare -> Rosetta 2(DE3)", rec_high["primary_strain"], "Rosetta 2(DE3)")
    check("high rare -> codon opt needed", rec_high["codon_optimization_needed"], True)
    check("high rare -> high priority", rec_high["optimization_priority"], "high")

    # Clusters force Rosetta
    rec_cluster = analyzer.recommend_strain({
        "rare_codon_pct": 2.0,
        "rare_codon_clusters": [{"start": 10, "end": 14, "rare_codons": []}],
        "cai": 0.65,
        "signal_peptide": {"has_signal_peptide": False},
    })
    check("cluster -> Rosetta", rec_cluster["primary_strain"], "Rosetta(DE3)")
    check("cluster -> codon opt", rec_cluster["codon_optimization_needed"], True)
    print()

    # ── Test 4: CAI calculation accuracy ──────────────────────────────────
    print(sep)
    print("Test 4: CAI calculation")
    print(sep)

    # All optimal codons for E. coli should give CAI ~1.0
    optimal_cds = "ATGAAAGCGTTTCTGTAA"  # M K A F L *
    # ATG(1.0) AAA(0.77) GCG(0.50) TTT(0.58) CTG(1.0) -- stop excluded
    # geometric mean = exp((ln(1)+ln(0.77)+ln(0.50)+ln(0.58)+ln(1))/5)
    expected_cai = math.exp(
        (math.log(1.0) + math.log(0.77) + math.log(0.50) + math.log(0.58) + math.log(1.0)) / 5
    )
    actual_cai = analyzer.compute_cai(optimal_cds)
    print(f"  Test CDS: {optimal_cds}")
    print(f"  Expected CAI: {expected_cai:.4f}")
    print(f"  Actual CAI:   {actual_cai:.4f}")
    check("CAI accuracy", actual_cai, expected_cai, tolerance=0.001)

    # All worst codons -> low CAI
    worst_cds = "ATGAGGCUAAUACCCGGATAA"  # using RNA-like U? no, use DNA
    # Let's do: AGG(0.02) CTA(0.04) ATA(0.07) CCC(0.12) GGA(0.11)
    worst_dna = "ATGAGGCTAATACCCGGATAA"
    worst_cai = analyzer.compute_cai(worst_dna)
    print(f"  Worst CDS CAI: {worst_cai:.4f}")
    check("worst codons CAI < 0.15", worst_cai < 0.15, True)

    # GFP CAI should be high
    gfp_cai = analyzer.compute_cai(ecoli_gfp)
    print(f"  GFP CAI: {gfp_cai:.4f}")
    check("GFP CAI > 0.7", gfp_cai > 0.7, True)
    print()

    # ── Test 5: Signal peptide detection ──────────────────────────────────
    print(sep)
    print("Test 5: Signal peptide detection")
    print(sep)

    # Known signal peptide: MalE-like (MKYLLPTAAAGLLLLAAQPAMA...)
    test_signal_protein = "MKYLLPTAAAGLLLLAAQPAMA" + "DIVLTQSPASLAVSLGQRATIS"
    sp_result = analyzer.check_signal_peptide(test_signal_protein)
    print(f"  Test protein: {test_signal_protein[:30]}...")
    print(f"  Has signal peptide: {sp_result['has_signal_peptide']}")
    print(f"  Confidence: {sp_result['confidence']}")
    print(f"  Hydrophobic%: {sp_result['hydrophobic_pct']:.0%}")
    print(f"  N-region K/R: {sp_result['n_region_positive']}")
    print(f"  Cleavage site: {sp_result['cleavage_site_estimate']}")
    for note in sp_result["notes"]:
        print(f"    - {note}")
    check("signal peptide detected", sp_result["has_signal_peptide"], True)
    check("confidence >= medium",
          sp_result["confidence"] in ("high", "medium"), True)

    # Non-signal peptide: cytoplasmic protein
    cytoplasmic = "MSKGEELFTGVVPILVELDGDVNGHKFSVR"
    sp_cyto = analyzer.check_signal_peptide(cytoplasmic)
    print(f"\n  Cytoplasmic protein: {cytoplasmic[:30]}...")
    print(f"  Has signal peptide: {sp_cyto['has_signal_peptide']}")
    print(f"  Confidence: {sp_cyto['confidence']}")
    # Cytoplasmic protein should not have strong signal peptide
    check("no signal peptide for GFP-like",
          sp_cyto["confidence"] in ("none", "low"), True)
    print()

    # ── Test 6: Rare codon cluster detection ──────────────────────────────
    print(sep)
    print("Test 6: Rare codon cluster detection")
    print(sep)

    # Construct CDS with a known cluster: 3 rare codons in 5-codon window
    # AGG (Arg, rare) + AGA (Arg, rare) + ATG (Met, ok) + CCC (Pro, rare) + GCT (Ala, ok)
    cluster_cds = "ATG" + "AGG" "AGA" "ATG" "CCC" "GCT" + "AAA" * 10 + "TAA"
    result_cluster = analyzer.analyze(cluster_cds)
    print(f"  CDS: {cluster_cds[:60]}...")
    print(f"  Rare codons: {result_cluster['rare_codons']['total_rare']}")
    print(f"  Clusters found: {len(result_cluster['rare_codon_clusters'])}")
    for cl in result_cluster["rare_codon_clusters"]:
        print(f"    Cluster: codons {cl['start']}-{cl['end']}, "
              f"rare: {[rc['codon'] for rc in cl['rare_codons']]}")
    check("cluster detected", len(result_cluster["rare_codon_clusters"]) >= 1, True)

    # Verify cluster contains the expected rare codons
    if result_cluster["rare_codon_clusters"]:
        first_cluster = result_cluster["rare_codon_clusters"][0]
        cluster_codons = {rc["codon"] for rc in first_cluster["rare_codons"]}
        check("cluster contains AGG", "AGG" in cluster_codons, True)
        check("cluster contains AGA", "AGA" in cluster_codons, True)
        check("cluster contains CCC", "CCC" in cluster_codons, True)
    print()

    # ── Test 7: MAP removal prediction ────────────────────────────────────
    print(sep)
    print("Test 7: MAP removal prediction")
    print(sep)

    # Second residue = A (small) -> MAP removes Met
    map_cds = "ATGGCGAAATAA"  # MAK*
    result_map = analyzer.analyze(map_cds)
    check("MAP removes Met (second=A)",
          result_map["map_removal"]["will_be_removed"], True)

    # Second residue = D (not small) -> Met retained
    nomap_cds = "ATGGATAAATAA"  # MDK*
    result_nomap = analyzer.analyze(nomap_cds)
    check("MAP retains Met (second=D)",
          result_nomap["map_removal"]["will_be_removed"], False)
    print()

    # ── Summary ───────────────────────────────────────────────────────────
    print(sep)
    total = passed + failed
    print(f"Results: {passed}/{total} passed" +
          (f", {failed} FAILED" if failed else " -- ALL OK"))
    print(sep)
    return failed == 0


if __name__ == "__main__":
    import sys
    ok = _run_tests()
    sys.exit(0 if ok else 1)
