#!/usr/bin/env python3
"""
Restriction Cloning Primer Designer
=====================================
iPCRDesignerBase를 상속하여 RE(restriction enzyme) 클로닝용 프라이머 설계.

설계 원리 (RE cloning):
  Forward: 5'-[protection 4-6bp]-[RE5 recognition]-[spacer_5prime]-[insert 5' annealing ~20bp]-3'
  Reverse: 5'-[protection 4-6bp]-[RE3 recognition]-[spacer_3prime]-[RC(stop codon)?]-[insert 3' RC annealing ~20bp]-3'

  Tm 계산: annealing 부분만 (RE tail은 template과 결합하지 않으므로 Tm에 기여하지 않음).

특수 사례:
  NdeI/NcoI: recognition site에 ATG 포함 → include_start_codon 시 RE site 자체가 ATG 제공
  Compatible overhang: BamHI+BglII 등 → 비방향성 삽입 경고
  Blunt-end: EcoRV → 방향성 없음 경고

Data sources:
  NEB (New England Biolabs) enzyme buffer compatibility
  NEB minimum protection bases guidelines
"""

from __future__ import annotations

from Bio.Seq import Seq

from .subst_primer_mode import iPCRDesignerBase
from .vector_registry import (
    RESTRICTION_ENZYMES,
    check_reading_frame,
    format_frame_report,
    get_vector,
)


# ── NEB Buffer Compatibility ──────────────────────────────────────────────────

RE_BUFFER_INFO: dict[str, dict] = {
    "BamHI-HF":   {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": True},
    "BamHI":      {"buffer": "NEBuffer 3.1", "activity_cutsmart": 75, "hf": False},
    "XhoI":       {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": False},
    "NheI-HF":    {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": True},
    "NotI-HF":    {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": True},
    "NdeI":       {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": False},
    "EcoRI-HF":   {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": True},
    "HindIII-HF": {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": True},
    "SalI-HF":    {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": True},
    "KpnI-HF":    {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": True},
    "PstI-HF":    {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": True},
    "NcoI-HF":    {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": True},
    "NcoI":       {"buffer": "NEBuffer 3.1", "activity_cutsmart": 50, "hf": False},
    "BglII":      {"buffer": "NEBuffer 3.1", "activity_cutsmart": 50, "hf": False},
    "EcoRV-HF":   {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": True},
    "EcoRV":      {"buffer": "NEBuffer 3.1", "activity_cutsmart": 50, "hf": False},
    "MfeI-HF":    {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": True},
    "SacI-HF":    {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": True},
    "AvrII":      {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": False},
    "FseI":       {"buffer": "NEBuffer 4", "activity_cutsmart": 10, "hf": False},
    "AscI":       {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": False},
    "PacI":       {"buffer": "CutSmart", "activity_cutsmart": 100, "hf": False},
}

# ── Minimum Protection Bases (NEB guidelines) ─────────────────────────────────

RE_MIN_PROTECTION: dict[int, int] = {
    6: 4,   # 6-cutter: 4 bp protection
    8: 6,   # 8-cutter: 6 bp protection
}

# ── Compatible Overhang Pairs ──────────────────────────────────────────────────

COMPATIBLE_OVERHANGS: list[tuple[str, str]] = [
    ("BamHI", "BglII"),   # both produce 5'-GATC overhang
    ("EcoRI", "MfeI"),    # both produce 5'-AATT overhang
    ("NheI", "AvrII"),    # both produce 5'-CTAG overhang (SpeI, XbaI too)
    ("SalI", "XhoI"),     # both produce 5'-TCGA overhang
]


def _get_base_enzyme(name: str) -> str:
    """HF variant 등에서 base enzyme 이름 추출. 'BamHI-HF' -> 'BamHI'."""
    return name.replace("-HF", "")


def _find_all_occurrences(seq: str, pattern: str) -> list[int]:
    """seq 내에서 pattern의 모든 시작 위치 반환."""
    positions = []
    start = 0
    while True:
        idx = seq.upper().find(pattern.upper(), start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


class RestrictionCloningDesigner(iPCRDesignerBase):
    """RE 클로닝용 프라이머 설계 (NEB buffer 호환성 + reading frame 검증 포함)."""

    def design(
        self,
        insert_seq: str,
        re_5prime: str,
        re_3prime: str,
        vector_name: str | None = None,
        target_tm: float = 62.0,
        min_ann_len: int = 18,
        max_ann_len: int = 30,
        include_start_codon: bool = True,
        include_stop_codon: bool = False,
        stop_codon: str = "TAA",
        protection_bases_5: int | str | None = None,
        protection_bases_3: int | str | None = None,
        spacer_5prime: str = "",
        spacer_3prime: str = "",
        auto_frame_check: bool = True,
    ) -> dict:
        """RE 클로닝 프라이머 설계 메인 메서드.

        Parameters
        ----------
        insert_seq : str
            삽입할 CDS 서열 (ATG부터, stop codon 포함/미포함 모두 가능)
        re_5prime, re_3prime : str
            5'/3' restriction enzyme 이름
        vector_name : str or None
            벡터 이름 (reading frame 검증용)
        target_tm : float
            annealing 목표 Tm (degC)
        min_ann_len, max_ann_len : int
            annealing 최소/최대 길이
        include_start_codon : bool
            프라이머에 start codon(ATG) 포함 여부
        include_stop_codon : bool
            reverse primer에 stop codon 포함 여부
        stop_codon : str
            stop codon 서열 (default: TAA)
        protection_bases_5 : int, str, or None
            5' protection bases (None=자동, int=길이, str=리터럴 서열)
        protection_bases_3 : int, str, or None
            3' protection bases (None=자동, int=길이, str=리터럴 서열)
        spacer_5prime, spacer_3prime : str
            RE site 뒤에 추가할 spacer 서열
        auto_frame_check : bool
            vector_name 제공 시 자동 frame check 수행 여부

        Returns
        -------
        dict : f_full, r_full, f_ann, r_ann, f_tail, r_tail, QC 결과 등
        """
        insert_seq = insert_seq.upper().replace(" ", "")
        warnings = []

        # ── 1. Input validation ──────────────────────────────────────────
        if not insert_seq:
            raise ValueError("insert_seq is empty")
        if not all(c in "ATGC" for c in insert_seq):
            raise ValueError("insert_seq contains non-ATGC characters")

        if re_5prime not in RESTRICTION_ENZYMES:
            raise ValueError(
                f"Unknown 5' RE: {re_5prime!r}. "
                f"Available: {', '.join(sorted(RESTRICTION_ENZYMES.keys()))}"
            )
        if re_3prime not in RESTRICTION_ENZYMES:
            raise ValueError(
                f"Unknown 3' RE: {re_3prime!r}. "
                f"Available: {', '.join(sorted(RESTRICTION_ENZYMES.keys()))}"
            )

        re5_info = RESTRICTION_ENZYMES[re_5prime]
        re3_info = RESTRICTION_ENZYMES[re_3prime]
        re5_site = re5_info["recognition"]
        re3_site = re3_info["recognition"]

        # ── 2. Internal RE site scan ─────────────────────────────────────
        internal_re_sites_5 = _find_all_occurrences(insert_seq, re5_site)
        internal_re_sites_3 = _find_all_occurrences(insert_seq, re3_site)

        if internal_re_sites_5:
            warnings.append(
                f"INSERT contains {re_5prime} site ({re5_site}) at position(s): "
                f"{internal_re_sites_5} - may be cut during digestion!"
            )
        if internal_re_sites_3:
            warnings.append(
                f"INSERT contains {re_3prime} site ({re3_site}) at position(s): "
                f"{internal_re_sites_3} - may be cut during digestion!"
            )

        # ── 3. Protection bases determination ────────────────────────────
        protection_5 = self._resolve_protection(protection_bases_5, re5_site)
        protection_3 = self._resolve_protection(protection_bases_3, re3_site)

        # ── 4. NdeI/NcoI special case: ATG in recognition site ───────────
        re5_has_atg = "ATG" in re5_site
        if re5_has_atg and include_start_codon:
            warnings.append(
                f"{re_5prime} recognition site ({re5_site}) contains ATG - "
                f"RE site itself provides start codon"
            )

        # ── 5. Compatible overhang warning ───────────────────────────────
        base_5 = _get_base_enzyme(re_5prime)
        base_3 = _get_base_enzyme(re_3prime)
        for pair in COMPATIBLE_OVERHANGS:
            if (base_5 in pair and base_3 in pair) and base_5 != base_3:
                warnings.append(
                    f"{re_5prime} and {re_3prime} produce compatible sticky ends - "
                    f"insert can ligate in BOTH orientations (non-directional)!"
                )
                break

        # Same enzyme on both sides
        if base_5 == base_3:
            warnings.append(
                f"Same RE ({re_5prime}/{re_3prime}) on both sides - "
                f"insert can ligate in BOTH orientations (non-directional)!"
            )

        # Blunt-end warning
        if re5_info["overhang"] == "blunt" or re3_info["overhang"] == "blunt":
            blunt_re = re_5prime if re5_info["overhang"] == "blunt" else re_3prime
            warnings.append(
                f"{blunt_re} produces blunt ends - "
                f"no directionality from this side"
            )

        # ── 6. Forward primer assembly ───────────────────────────────────
        f_tail = protection_5 + re5_site + spacer_5prime

        # Forward annealing: RE tail does NOT contribute to Tm
        # Fallback: extend range then relax Tm for AT-rich regions
        f_ann, f_tm, f_gc, f_ann_len = self._design_annealing(
            template=insert_seq, anchor=0, direction="+",
            target_tm=target_tm, min_len=min_ann_len, max_len=max_ann_len,
            tail_seq="",
        )
        if f_ann is None:
            # Fallback 1: extend max to 35 bp
            f_ann, f_tm, f_gc, f_ann_len = self._design_annealing(
                template=insert_seq, anchor=0, direction="+",
                target_tm=target_tm, min_len=min_ann_len, max_len=35,
                tail_seq="",
            )
        if f_ann is None:
            # Fallback 2: relax Tm by 4°C with extended range
            f_ann, f_tm, f_gc, f_ann_len = self._design_annealing(
                template=insert_seq, anchor=0, direction="+",
                target_tm=target_tm - 4.0, min_len=min_ann_len, max_len=36,
                tail_seq="",
            )
            if f_ann is not None:
                warnings.append(
                    f"Forward primer Tm ({f_tm:.1f}°C) is below target "
                    f"({target_tm}°C) due to AT-rich region"
                )
        if f_ann is None:
            raise RuntimeError(
                "Forward primer annealing design failed "
                f"(target Tm={target_tm}C, range={min_ann_len}-{max_ann_len} bp)"
            )

        # ── 7. Reverse primer assembly ───────────────────────────────────
        stop_rc = ""
        if include_stop_codon:
            stop_rc = str(Seq(stop_codon.upper()).reverse_complement())

        r_tail = protection_3 + re3_site + spacer_3prime + stop_rc

        # Reverse annealing: RE tail does NOT contribute to Tm
        # Fallback: extend range then relax Tm for AT-rich regions
        r_ann, r_tm, r_gc, r_ann_len = self._design_annealing(
            template=insert_seq, anchor=len(insert_seq), direction="-",
            target_tm=target_tm, min_len=min_ann_len, max_len=max_ann_len,
            tail_seq="",
        )
        if r_ann is None:
            # Fallback 1: extend max to 35 bp
            r_ann, r_tm, r_gc, r_ann_len = self._design_annealing(
                template=insert_seq, anchor=len(insert_seq), direction="-",
                target_tm=target_tm, min_len=min_ann_len, max_len=35,
                tail_seq="",
            )
        if r_ann is None:
            # Fallback 2: relax Tm by 4°C with extended range
            r_ann, r_tm, r_gc, r_ann_len = self._design_annealing(
                template=insert_seq, anchor=len(insert_seq), direction="-",
                target_tm=target_tm - 4.0, min_len=min_ann_len, max_len=36,
                tail_seq="",
            )
            if r_ann is not None:
                warnings.append(
                    f"Reverse primer Tm ({r_tm:.1f}°C) is below target "
                    f"({target_tm}°C) due to AT-rich region"
                )
        if r_ann is None:
            raise RuntimeError(
                "Reverse primer annealing design failed "
                f"(target Tm={target_tm}C, range={min_ann_len}-{max_ann_len} bp)"
            )

        # ── 8. 2-pass hairpin avoidance ──────────────────────────────────
        anneal_temp = min(f_tm, r_tm) + 1.0

        f_qc = self.check_primer(f_tail + f_ann, anneal_temp)
        r_qc = self.check_primer(r_tail + r_ann, anneal_temp)

        if f_qc["verdict"] == "FAIL":
            alt = self._design_annealing(
                template=insert_seq, anchor=0, direction="+",
                target_tm=target_tm, min_len=min_ann_len, max_len=max_ann_len,
                tail_seq="", anneal_temp=anneal_temp,
            )
            if alt[0] is not None and alt[0] != f_ann:
                f_ann, f_tm, f_gc, f_ann_len = alt
                warnings.append(
                    f"F annealing adjusted for hairpin avoidance ({len(f_ann)} bp)"
                )
            else:
                warnings.append(
                    "F hairpin: annealing length adjustment cannot resolve"
                )

        if r_qc["verdict"] == "FAIL":
            alt = self._design_annealing(
                template=insert_seq, anchor=len(insert_seq), direction="-",
                target_tm=target_tm, min_len=min_ann_len, max_len=max_ann_len,
                tail_seq="", anneal_temp=anneal_temp,
            )
            if alt[0] is not None and alt[0] != r_ann:
                r_ann, r_tm, r_gc, r_ann_len = alt
                warnings.append(
                    f"R annealing adjusted for hairpin avoidance ({len(r_ann)} bp)"
                )
            else:
                warnings.append(
                    "R hairpin: annealing length adjustment cannot resolve"
                )

        # ── 9. Assemble final primers ────────────────────────────────────
        f_full = f_tail + f_ann
        r_full = r_tail + r_ann

        # ── 10. Primer-level warnings ────────────────────────────────────
        for label, primer in [("F", f_full), ("R", r_full)]:
            hp = self.homopolymer_run(primer)
            if hp:
                warnings.append(f"{label} primer homopolymer: {hp}")
            if not self.gc_clamp_ok(primer):
                warnings.append(f"{label} primer no 3' G/C clamp")
            if len(primer) > 60:
                warnings.append(f"{label} primer length ({len(primer)} nt) > 60 nt")

        # ── 11. Final QC ─────────────────────────────────────────────────
        anneal_temp = min(f_tm, r_tm) + 1.0
        f_qc = self.check_primer(f_full, anneal_temp)
        r_qc = self.check_primer(r_full, anneal_temp)
        het = self.check_heterodimer(f_full, r_full)

        # ── 12. Reading frame check ──────────────────────────────────────
        frame_check = None
        frame_report = None

        if vector_name is not None and auto_frame_check:
            try:
                frame_check = check_reading_frame(
                    vector_name=vector_name,
                    re_5prime=re_5prime,
                    re_3prime=re_3prime,
                    insert_has_atg=include_start_codon,
                    insert_has_stop=include_stop_codon,
                    insert_cds_bp=len(insert_seq),
                )
                frame_report = format_frame_report(frame_check)
            except ValueError as e:
                warnings.append(f"Frame check failed: {e}")

        # ── 13. Build result ─────────────────────────────────────────────
        result = {
            "f_full": f_full,
            "r_full": r_full,
            "f_ann": f_ann,
            "r_ann": r_ann,
            "f_tail": f_tail,
            "r_tail": r_tail,
            "f_tm": round(f_tm, 1),
            "r_tm": round(r_tm, 1),
            "f_gc": round(f_gc, 1),
            "r_gc": round(r_gc, 1),
            "f_len": len(f_full),
            "r_len": len(r_full),
            "f_ann_len": f_ann_len,
            "r_ann_len": r_ann_len,
            "re_5prime": re_5prime,
            "re_3prime": re_3prime,
            "re_5prime_site": re5_site,
            "re_3prime_site": re3_site,
            "protection_5": protection_5,
            "protection_3": protection_3,
            "spacer_5prime": spacer_5prime,
            "spacer_3prime": spacer_3prime,
            "include_start_codon": include_start_codon,
            "include_stop_codon": include_stop_codon,
            "insert_len": len(insert_seq),
            "anneal_temp": round(anneal_temp, 1),
            "f_qc": f_qc,
            "r_qc": r_qc,
            "het": het,
            "internal_re_sites_5": internal_re_sites_5,
            "internal_re_sites_3": internal_re_sites_3,
            "warnings": warnings,
        }

        if frame_check is not None:
            result["frame_check"] = frame_check
            result["frame_report"] = frame_report

        return result

    def recommend_re_pair(
        self,
        insert_seq: str,
        vector_name: str,
        prefer_hf: bool = True,
    ) -> list[dict]:
        """Insert와 벡터에 최적인 RE 쌍을 추천.

        Parameters
        ----------
        insert_seq : str
            삽입할 CDS 서열
        vector_name : str
            벡터 이름
        prefer_hf : bool
            HF variant 선호 여부

        Returns
        -------
        list[dict] : score 높은 순으로 정렬된 RE 쌍 추천 목록
        """
        insert_seq = insert_seq.upper().replace(" ", "")
        vec = get_vector(vector_name)
        re_sites = vec["re_sites"]

        # 벡터 MCS에 있는 RE 목록 (위치 순)
        re_list = sorted(re_sites.items(), key=lambda x: x[1])

        recommendations = []

        for i, (re5_name, re5_pos) in enumerate(re_list):
            for re3_name, re3_pos in re_list[i + 1:]:
                # 5' RE must be upstream of 3' RE
                if re5_pos >= re3_pos:
                    continue

                # RE info 가져오기
                if re5_name not in RESTRICTION_ENZYMES:
                    continue
                if re3_name not in RESTRICTION_ENZYMES:
                    continue

                re5_rec = RESTRICTION_ENZYMES[re5_name]["recognition"]
                re3_rec = RESTRICTION_ENZYMES[re3_name]["recognition"]

                # Insert 내부에 RE site이 있으면 skip
                if _find_all_occurrences(insert_seq, re5_rec):
                    continue
                if _find_all_occurrences(insert_seq, re3_rec):
                    continue

                # Scoring
                score = 0
                reasons = []

                # Double digest compatibility (same buffer >= 75%)
                buf5 = self._get_buffer_info(re5_name)
                buf3 = self._get_buffer_info(re3_name)
                double_digest_ok = False
                buffer_name = None

                if buf5 and buf3:
                    if buf5["buffer"] == buf3["buffer"]:
                        double_digest_ok = True
                        buffer_name = buf5["buffer"]
                        score += 3
                        reasons.append(f"Same buffer ({buffer_name})")
                    elif (buf5["activity_cutsmart"] >= 75
                          and buf3["activity_cutsmart"] >= 75):
                        double_digest_ok = True
                        buffer_name = "CutSmart"
                        score += 3
                        reasons.append("Both >= 75% in CutSmart")
                    else:
                        reasons.append("Sequential digest recommended")

                # CutSmart buffer
                if buf5 and buf3:
                    if (buf5["activity_cutsmart"] == 100
                            and buf3["activity_cutsmart"] == 100):
                        score += 2
                        reasons.append("Both 100% in CutSmart")

                # Reading frame compatibility
                frame_ok = False
                try:
                    fc = check_reading_frame(
                        vector_name=vector_name,
                        re_5prime=re5_name,
                        re_3prime=re3_name,
                        insert_has_atg=True,
                        insert_has_stop=False,
                    )
                    if fc["in_frame_5prime"] and fc["in_frame_3prime"]:
                        frame_ok = True
                        score += 2
                        reasons.append("Reading frame compatible")
                    else:
                        reasons.append("Frame mismatch")
                except ValueError:
                    reasons.append("Frame check N/A")

                # HF version available
                hf_available = False
                if prefer_hf:
                    hf5 = f"{re5_name}-HF" in RE_BUFFER_INFO
                    hf3 = f"{re3_name}-HF" in RE_BUFFER_INFO
                    if hf5 or hf3:
                        hf_available = True
                        score += 1
                        hf_names = []
                        if hf5:
                            hf_names.append(f"{re5_name}-HF")
                        if hf3:
                            hf_names.append(f"{re3_name}-HF")
                        reasons.append(f"HF available: {', '.join(hf_names)}")

                recommendations.append({
                    "re_5prime": re5_name,
                    "re_3prime": re3_name,
                    "buffer": buffer_name,
                    "double_digest_ok": double_digest_ok,
                    "frame_ok": frame_ok,
                    "score": score,
                    "reason": "; ".join(reasons),
                })

        # Sort by score descending
        recommendations.sort(key=lambda x: -x["score"])
        return recommendations

    # ── Internal helpers ─────────────────────────────────────────────────

    @staticmethod
    def _resolve_protection(
        protection: int | str | None,
        recognition: str,
    ) -> str:
        """Protection bases 결정.

        None -> RE_MIN_PROTECTION 기반 자동 생성 (GCGC... 패턴)
        int  -> 해당 길이만큼 자동 생성
        str  -> 리터럴 서열 사용
        """
        if isinstance(protection, str):
            return protection.upper()

        rec_len = len(recognition)
        if protection is None:
            # 가장 가까운 cutoff 이하 key 선택
            n_bases = RE_MIN_PROTECTION.get(rec_len)
            if n_bases is None:
                # fallback: recognition length 이하의 가장 큰 key
                applicable = [k for k in RE_MIN_PROTECTION if k <= rec_len]
                if applicable:
                    n_bases = RE_MIN_PROTECTION[max(applicable)]
                else:
                    n_bases = 4  # default
        else:
            n_bases = protection

        # GC-rich protection pattern (NEB recommendation)
        pattern = "GCGC"
        return (pattern * ((n_bases // len(pattern)) + 1))[:n_bases]

    @staticmethod
    def _get_buffer_info(enzyme_name: str) -> dict | None:
        """RE_BUFFER_INFO에서 효소 정보 조회 (HF variant 우선)."""
        hf_name = f"{enzyme_name}-HF"
        if hf_name in RE_BUFFER_INFO:
            return RE_BUFFER_INFO[hf_name]
        if enzyme_name in RE_BUFFER_INFO:
            return RE_BUFFER_INFO[enzyme_name]
        return None


# ── 테스트 ─────────────────────────────────────────────────────────────────────

def _run_tests():
    """RestrictionCloningDesigner 테스트."""
    sep = "=" * 70
    line = "-" * 70
    designer = RestrictionCloningDesigner()
    passed = 0
    failed = 0

    def check(label, actual, expected):
        nonlocal passed, failed
        if actual == expected:
            passed += 1
            print(f"  OK: {label}")
        else:
            failed += 1
            print(f"  FAIL: {label}")
            print(f"    expected: {expected!r}")
            print(f"    actual:   {actual!r}")

    # Test insert CDS (no stop codon)
    test_insert = (
        "ATGAAAGCTGCCATTGTTCTGAGCGAATCCGGTGTGCACGAATGTCCGAAAGAATCCAAACTG"
        "GTTCCGGCACTGAACGGTCTGGAAATCGAAGATGAAGGCGTTATCCCGGAATTCTTCAAGGGC"
        "AAACTGGATCGCGGCTTCAAAGCGATCCTGAACAATGCGAAAGATTGGAGCCGCGTGGAAAAC"
        "TACCAGCCGGATCTGATCAACGAACTGAAAGCGAAAGCGACCGATGTGGAATAA"
    )

    # ── Test 1: BamHI + XhoI → pET-28a(+) ────────────────────────────────
    print(f"\n{sep}")
    print("  Test 1: BamHI / XhoI -> pET-28a(+) -- frame check PASS")
    print(sep)

    r1 = designer.design(
        insert_seq=test_insert,
        re_5prime="BamHI",
        re_3prime="XhoI",
        vector_name="pET-28a(+)",
        target_tm=62.0,
        include_start_codon=True,
        include_stop_codon=False,
    )

    print(f"  F: 5'-{r1['f_full']}-3'  ({r1['f_len']} nt, Tm={r1['f_tm']}C)")
    print(f"     tail: {r1['f_tail']}  ann: {r1['f_ann']} ({r1['f_ann_len']} bp)")
    print(f"     QC: {r1['f_qc']['verdict']}")
    print(f"  R: 5'-{r1['r_full']}-3'  ({r1['r_len']} nt, Tm={r1['r_tm']}C)")
    print(f"     tail: {r1['r_tail']}  ann: {r1['r_ann']} ({r1['r_ann_len']} bp)")
    print(f"     QC: {r1['r_qc']['verdict']}")
    print(f"  Anneal temp: {r1['anneal_temp']}C")
    print(f"  Het dimer: dG={r1['het'].get('dg','N/A')} kcal/mol")

    if "frame_report" in r1:
        print(f"\n  Frame Report:")
        for fline in r1["frame_report"].split("\n"):
            print(f"    {fline}")

    if r1["warnings"]:
        print(f"\n  Warnings:")
        for w in r1["warnings"]:
            print(f"    - {w}")

    check("RE 5' = BamHI", r1["re_5prime"], "BamHI")
    check("RE 3' = XhoI", r1["re_3prime"], "XhoI")
    check("5' site = GGATCC", r1["re_5prime_site"], "GGATCC")
    check("3' site = CTCGAG", r1["re_3prime_site"], "CTCGAG")
    check("F primer starts with protection + RE",
          r1["f_full"].startswith(r1["protection_5"] + "GGATCC"), True)
    check("R primer starts with protection + RE",
          r1["r_full"].startswith(r1["protection_3"] + "CTCGAG"), True)
    check("frame_check exists", "frame_check" in r1, True)
    if "frame_check" in r1:
        check("5' in-frame", r1["frame_check"]["in_frame_5prime"], True)
        check("3' in-frame", r1["frame_check"]["in_frame_3prime"], True)
    print()

    # ── Test 2: NheI + NotI → pET-21a(+) with stop codon ────────────────
    print(f"{sep}")
    print("  Test 2: NheI / NotI -> pET-21a(+) with stop codon")
    print(sep)

    r2 = designer.design(
        insert_seq=test_insert,
        re_5prime="NheI",
        re_3prime="NotI",
        vector_name="pET-21a(+)",
        target_tm=62.0,
        include_start_codon=True,
        include_stop_codon=True,
        stop_codon="TAA",
    )

    print(f"  F: 5'-{r2['f_full']}-3'  ({r2['f_len']} nt)")
    print(f"  R: 5'-{r2['r_full']}-3'  ({r2['r_len']} nt)")
    print(f"  F tail: {r2['f_tail']}")
    print(f"  R tail: {r2['r_tail']}")

    if "frame_report" in r2:
        print(f"\n  Frame Report:")
        for fline in r2["frame_report"].split("\n"):
            print(f"    {fline}")

    if r2["warnings"]:
        print(f"\n  Warnings:")
        for w in r2["warnings"]:
            print(f"    - {w}")

    check("include_stop_codon", r2["include_stop_codon"], True)
    # R tail should contain RC(TAA) = TTA
    check("R tail contains RC(stop)", "TTA" in r2["r_tail"], True)
    check("frame_check exists", "frame_check" in r2, True)
    if "frame_check" in r2:
        check("stop warning present",
              any("C-terminal tag will NOT" in w
                  for w in r2["frame_check"]["warnings"]), True)
    print()

    # ── Test 3: Internal RE site detection ───────────────────────────────
    print(f"{sep}")
    print("  Test 3: Internal RE site detection")
    print(sep)

    # Insert that contains BamHI site (GGATCC)
    insert_with_bamhi = (
        "ATGAAAGCTGCCATTGTTCTGAGCGAATCCGGATCCGCACGAATGTCCGAAAGAATCCAAAC"
        "TGGTTCCGGCACTGAACGGTCTGGAAATCGAAGATGAAGGCGTTATCCCGGAATTCTTCAAG"
    )

    r3 = designer.design(
        insert_seq=insert_with_bamhi,
        re_5prime="BamHI",
        re_3prime="XhoI",
        target_tm=62.0,
    )

    print(f"  Internal BamHI sites: {r3['internal_re_sites_5']}")
    if r3["warnings"]:
        print(f"  Warnings:")
        for w in r3["warnings"]:
            print(f"    - {w}")

    check("internal RE site detected", len(r3["internal_re_sites_5"]) > 0, True)
    check("warning present",
          any("may be cut" in w for w in r3["warnings"]), True)
    print()

    # ── Test 4: Custom protection bases ──────────────────────────────────
    print(f"{sep}")
    print("  Test 4: Custom protection bases")
    print(sep)

    # Literal string protection
    r4a = designer.design(
        insert_seq=test_insert,
        re_5prime="BamHI",
        re_3prime="XhoI",
        protection_bases_5="AATTCC",
        protection_bases_3=6,
    )

    print(f"  F tail (literal 'AATTCC'): {r4a['f_tail']}")
    print(f"  R tail (int 6):            {r4a['r_tail']}")

    check("5' protection = AATTCC", r4a["protection_5"], "AATTCC")
    check("3' protection length = 6", len(r4a["protection_3"]), 6)
    check("F tail starts with AATTCC",
          r4a["f_tail"].startswith("AATTCC"), True)

    # Integer protection
    r4b = designer.design(
        insert_seq=test_insert,
        re_5prime="NotI",
        re_3prime="XhoI",
        protection_bases_5=None,  # auto: 8-cutter -> 6 bp
        protection_bases_3=None,  # auto: 6-cutter -> 4 bp
    )

    print(f"\n  NotI (8-cutter) auto protection: {r4b['protection_5']} ({len(r4b['protection_5'])} bp)")
    print(f"  XhoI (6-cutter) auto protection: {r4b['protection_3']} ({len(r4b['protection_3'])} bp)")

    check("NotI auto protection = 6 bp", len(r4b["protection_5"]), 6)
    check("XhoI auto protection = 4 bp", len(r4b["protection_3"]), 4)
    print()

    # ── Test 5: recommend_re_pair() ──────────────────────────────────────
    print(f"{sep}")
    print("  Test 5: recommend_re_pair()")
    print(sep)

    recs = designer.recommend_re_pair(
        insert_seq=test_insert,
        vector_name="pET-28a(+)",
        prefer_hf=True,
    )

    print(f"  Found {len(recs)} RE pair recommendations for pET-28a(+)")
    for i, rec in enumerate(recs[:5]):
        print(f"    {i+1}. {rec['re_5prime']}/{rec['re_3prime']} "
              f"(score={rec['score']}, "
              f"DD={'OK' if rec['double_digest_ok'] else 'No'}, "
              f"frame={'OK' if rec['frame_ok'] else 'No'})")
        print(f"       {rec['reason']}")

    check("recommendations not empty", len(recs) > 0, True)
    check("top recommendation has highest score",
          all(recs[0]["score"] >= r["score"] for r in recs), True)
    print()

    # ── Test 6: NdeI with ATG (special case) ─────────────────────────────
    print(f"{sep}")
    print("  Test 6: NdeI (ATG in recognition site)")
    print(sep)

    r6 = designer.design(
        insert_seq=test_insert,
        re_5prime="NdeI",
        re_3prime="XhoI",
        vector_name="pET-28a(+)",
        include_start_codon=True,
    )

    print(f"  F: 5'-{r6['f_full']}-3'")
    print(f"  F tail: {r6['f_tail']}")
    if r6["warnings"]:
        for w in r6["warnings"]:
            print(f"    - {w}")

    check("NdeI ATG warning",
          any("contains ATG" in w for w in r6["warnings"]), True)
    print()

    # ── Test 7: EcoRV blunt-end warning ──────────────────────────────────
    print(f"{sep}")
    print("  Test 7: EcoRV blunt-end warning")
    print(sep)

    r7 = designer.design(
        insert_seq=test_insert,
        re_5prime="EcoRV",
        re_3prime="XhoI",
    )

    if r7["warnings"]:
        for w in r7["warnings"]:
            print(f"    - {w}")

    check("blunt-end warning",
          any("blunt end" in w.lower() for w in r7["warnings"]), True)
    print()

    # ── Test 8: Compatible overhang warning (BamHI + BglII) ──────────────
    print(f"{sep}")
    print("  Test 8: Compatible overhang warning (BamHI + BglII)")
    print(sep)

    r8 = designer.design(
        insert_seq=test_insert,
        re_5prime="BamHI",
        re_3prime="BglII",
    )

    if r8["warnings"]:
        for w in r8["warnings"]:
            print(f"    - {w}")

    check("compatible overhang warning",
          any("BOTH orientations" in w for w in r8["warnings"]), True)
    print()

    # ── Summary ──────────────────────────────────────────────────────────
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
