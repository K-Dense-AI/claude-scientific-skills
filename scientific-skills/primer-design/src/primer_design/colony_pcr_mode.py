#!/usr/bin/env python3
"""
Colony PCR Primer Designer
===========================
Universal primer 기반 colony PCR 프라이머 추천 및
insert 내부 프라이머 설계.

Colony PCR는 Taq polymerase 기반이므로:
  - Annealing temp = min(F_tm, R_tm) - 5.0
  - Q5/Phusion 보정 없음

Universal primer DB는 Macrogen Standard Primer 목록 기준.
"""

from __future__ import annotations

from dataclasses import dataclass

from Bio.Seq import Seq

from .subst_primer_mode import iPCRDesignerBase
from .vector_registry import get_vector, EXPRESSION_VECTORS
from .snapgene_parser import parse_snapgene


# ── UniversalPrimer dataclass ──────────────────────────────────────────────

@dataclass
class UniversalPrimer:
    """Macrogen 등 표준 프라이머 정보."""
    name: str
    sequence: str       # 5'->3'
    direction: str      # "forward" or "reverse"
    catalog: str        # e.g., "Macrogen #51"


# ── Universal Primer DB (Macrogen Standard Primers) ────────────────────────

UNIVERSAL_PRIMERS = {
    "T7promoter":     UniversalPrimer("T7promoter",     "TAATACGACTCACTATAGGG",    "forward",  "Macrogen #51"),
    "T7terminator":   UniversalPrimer("T7terminator",   "GCTAGTTATTGCTCAGCGG",     "reverse",  "Macrogen #50"),
    "T7":             UniversalPrimer("T7",             "AATACGACTCACTATAG",       "forward",  "Macrogen #49"),
    "pET_RP":         UniversalPrimer("pET-RP",         "CTAGTTATTGCTCAGCGG",      "reverse",  "Macrogen #14"),
    "pET_24a":        UniversalPrimer("pET-24a",        "GGGTTATGCTAGTTATTGCTCAG", "reverse",  "Macrogen #13"),
    "pET_upstream":   UniversalPrimer("pET-upstream",   "ATGCGTCCGGCGTAGAGG",      "forward",  "Macrogen #90"),
    "pMalE":          UniversalPrimer("pMalE",          "TCAGACTGTCGATGAAGC",      "forward",  "Macrogen #15"),
    "DuetDown1":      UniversalPrimer("DuetDown1",      "GATTATGCGGCCGTGTACAA",    "reverse",  "Macrogen #88"),
    "DuetUP2":        UniversalPrimer("DuetUP2",        "TTGTACACGGCCGCATAATC",    "forward",  "Macrogen #89"),
    "STag":           UniversalPrimer("STag 18mer",     "GAACGCCAGCACATGGAC",      "forward",  "Macrogen #34"),
    "BGH_R":          UniversalPrimer("BGH-R",          "TAGAAGGCACAGTCGAGG",      "reverse",  "Macrogen #64"),
}


# ── Vector → Primer Pair Mapping ──────────────────────────────────────────

VECTOR_PRIMER_MAP: dict[str, tuple[str, str]] = {
    "pET-21a(+)":       ("T7promoter",    "T7terminator"),
    "pET-28a(+)":       ("T7promoter",    "T7terminator"),
    "pETDuet-1:MCS1":   ("pET_upstream",  "T7terminator"),
    "pETDuet-1:MCS2":   ("DuetUP2",       "T7terminator"),
    "pACYCDuet-1:MCS1": ("pET_upstream",  "T7terminator"),
    "pACYCDuet-1:MCS2": ("DuetUP2",       "DuetDown1"),
    "pMAL-c6T":         ("pMalE",         "pET_RP"),
}


# ── Vector Flanking Defaults ──────────────────────────────────────────────
# (upstream_bp from primer to insert start, downstream_bp from insert end to primer)
# 값은 각 벡터의 universal primer binding site ~ MCS 사이 대략적인 거리.

VECTOR_FLANKING_DEFAULTS: dict[str, tuple[int, int]] = {
    "pET-21a(+)":       (200, 150),
    "pET-28a(+)":       (200, 150),
    "pETDuet-1:MCS1":   (150, 150),
    "pETDuet-1:MCS2":   (150, 150),
    "pACYCDuet-1:MCS1": (150, 150),
    "pACYCDuet-1:MCS2": (150, 150),
    "pMAL-c6T":         (200, 100),
}

# MCS gap: 빈 벡터(insert 없음)에서 universal primer 사이 MCS 영역 대략적 길이
_MCS_GAP_DEFAULTS: dict[str, int] = {
    "pET-21a(+)":       100,
    "pET-28a(+)":       160,
    "pETDuet-1:MCS1":   90,
    "pETDuet-1:MCS2":   130,
    "pACYCDuet-1:MCS1": 90,
    "pACYCDuet-1:MCS2": 130,
    "pMAL-c6T":         65,
}


# ── ColonyPCRDesigner ─────────────────────────────────────────────────────

class ColonyPCRDesigner(iPCRDesignerBase):
    """Colony PCR 프라이머 추천 (Taq-based).

    Universal primer 기반으로 colony PCR 조건을 추천하고,
    필요시 insert 내부 프라이머도 설계.
    """

    def _resolve_canonical_name(self, vector_name: str) -> str:
        """벡터 이름을 canonical name으로 변환.

        get_vector()는 data dict를 반환하므로,
        EXPRESSION_VECTORS에서 동일 dict 객체를 찾아 canonical key를 반환.
        """
        vec_data = get_vector(vector_name)
        for canonical, data in EXPRESSION_VECTORS.items():
            if data is vec_data:
                return canonical
        return vector_name

    def suggest(
        self,
        vector_name: str,
        insert_length_bp: int,
        vector_flanking_bp: int | tuple[int, int] | None = None,
    ) -> dict:
        """Universal primer 기반 colony PCR 조건 추천.

        Parameters
        ----------
        vector_name : str
            벡터 이름 (fuzzy matching 지원, e.g., "pET28a")
        insert_length_bp : int
            Insert 길이 (bp)
        vector_flanking_bp : int or tuple[int, int] or None
            Primer ~ insert 사이 flanking 거리.
            int: 양쪽 동일, tuple: (upstream, downstream), None: 기본값 사용.

        Returns
        -------
        dict
            f_name, f_seq, f_tm, r_name, r_seq, r_tm,
            expected_band_with_insert, expected_band_empty,
            anneal_temp, flanking_upstream, flanking_downstream, notes
        """
        # 1. Resolve vector
        get_vector(vector_name)  # validate
        canonical = self._resolve_canonical_name(vector_name)
        notes: list[str] = []

        # 2. Primer pair lookup
        if canonical not in VECTOR_PRIMER_MAP:
            available = ", ".join(VECTOR_PRIMER_MAP.keys())
            raise ValueError(
                f"No colony PCR primers mapped for {canonical!r}. "
                f"Available: {available}. Use register_primer_pair() to add."
            )

        f_key, r_key = VECTOR_PRIMER_MAP[canonical]
        f_primer = UNIVERSAL_PRIMERS[f_key]
        r_primer = UNIVERSAL_PRIMERS[r_key]

        # 3. Flanking distances
        if vector_flanking_bp is None:
            flanking = VECTOR_FLANKING_DEFAULTS.get(canonical, (200, 150))
            notes.append(f"Flanking distances: default for {canonical}")
        elif isinstance(vector_flanking_bp, int):
            flanking = (vector_flanking_bp, vector_flanking_bp)
        else:
            flanking = vector_flanking_bp

        flanking_up, flanking_dn = flanking

        # 4. Expected band sizes
        expected_with_insert = flanking_up + insert_length_bp + flanking_dn
        mcs_gap = _MCS_GAP_DEFAULTS.get(canonical, 100)
        expected_empty = flanking_up + mcs_gap + flanking_dn

        # 5. Tm calculation
        f_tm = self.calc_tm(f_primer.sequence)
        r_tm = self.calc_tm(r_primer.sequence)

        # 6. Annealing temp (Taq-based: min Tm - 5)
        anneal_temp = min(f_tm, r_tm) - 5.0

        # 7. Notes
        if abs(f_tm - r_tm) > 5.0:
            notes.append(
                f"Tm difference {abs(f_tm - r_tm):.1f}C between F/R; "
                f"consider touchdown PCR"
            )
        if insert_length_bp > 3000:
            notes.append(
                f"Insert {insert_length_bp} bp is large for Taq colony PCR; "
                f"consider extension time >= {insert_length_bp // 1000 + 1} min"
            )

        return {
            "vector_name": canonical,
            "f_name": f_primer.name,
            "f_seq": f_primer.sequence,
            "f_tm": round(f_tm, 1),
            "r_name": r_primer.name,
            "r_seq": r_primer.sequence,
            "r_tm": round(r_tm, 1),
            "expected_band_with_insert": expected_with_insert,
            "expected_band_empty": expected_empty,
            "anneal_temp": round(anneal_temp, 1),
            "flanking_upstream": flanking_up,
            "flanking_downstream": flanking_dn,
            "notes": notes,
        }

    def suggest_with_internal(
        self,
        vector_name: str,
        insert_seq: str,
        target_tm: float = 58.0,
        target_product_range: tuple[int, int] = (300, 1200),
    ) -> dict:
        """Universal primer + insert 내부 프라이머 함께 추천.

        Insert가 길 때 universal primer만으로는 밴드가 너무 크거나
        구분이 어려울 수 있으므로, insert 내부에 annealing하는
        프라이머 쌍도 함께 설계.

        Parameters
        ----------
        vector_name : str
            벡터 이름
        insert_seq : str
            Insert 전체 서열 (DNA)
        target_tm : float
            내부 프라이머 목표 Tm (Taq 기반, default 58C)
        target_product_range : tuple[int, int]
            내부 프라이머 쌍의 목표 product size 범위 (bp)

        Returns
        -------
        dict
            universal: suggest() 결과,
            internal_f_seq, internal_f_tm,
            internal_r_seq, internal_r_tm,
            internal_product_size, internal_notes
        """
        insert_len = len(insert_seq)

        # 1. Universal primer 추천
        universal = self.suggest(vector_name, insert_len)

        # 2. Internal forward primer: ~100bp from start
        internal_notes: list[str] = []
        f_start = min(100, insert_len // 4)
        int_f_seq, int_f_tm = self._pick_internal_primer(
            insert_seq, f_start, "+", target_tm,
        )

        # 3. Internal reverse primer: ~100bp from end
        r_anchor = max(insert_len - 100, insert_len * 3 // 4)
        int_r_seq, int_r_tm = self._pick_internal_primer(
            insert_seq, r_anchor, "-", target_tm,
        )

        # 4. Product size check
        if int_f_seq is not None and int_r_seq is not None:
            # Forward primer starts at f_start, reverse primer ends at r_anchor
            internal_product = r_anchor - f_start
            min_prod, max_prod = target_product_range

            if internal_product < min_prod:
                internal_notes.append(
                    f"Internal product {internal_product} bp < target min {min_prod} bp; "
                    f"insert may be too short for internal primers"
                )
            elif internal_product > max_prod:
                internal_notes.append(
                    f"Internal product {internal_product} bp > target max {max_prod} bp; "
                    f"consider adjusting primer positions"
                )
        else:
            internal_product = None
            internal_notes.append(
                "Internal primer design failed; insert may be too short"
            )

        return {
            "universal": universal,
            "internal_f_seq": int_f_seq,
            "internal_f_tm": round(int_f_tm, 1) if int_f_tm is not None else None,
            "internal_r_seq": int_r_seq,
            "internal_r_tm": round(int_r_tm, 1) if int_r_tm is not None else None,
            "internal_product_size": internal_product,
            "internal_notes": internal_notes,
        }

    def suggest_from_snapgene(
        self,
        snapgene_path: str,
        vector_name: str,
        cds_feature_name: str | None = None,
    ) -> dict:
        """SnapGene .dna 파일에서 CDS feature를 찾아 suggest() 호출.

        Parameters
        ----------
        snapgene_path : str
            SnapGene .dna 파일 경로
        vector_name : str
            벡터 이름
        cds_feature_name : str or None
            CDS feature 이름. None이면 첫 번째 CDS를 사용.

        Returns
        -------
        dict
            suggest() 결과 + snapgene_info
        """
        sequence, is_circular, features = parse_snapgene(snapgene_path)

        # CDS feature 찾기
        cds_features = [f for f in features if f["type"] == "CDS"]
        if not cds_features:
            # CDS가 없으면 gene type도 시도
            cds_features = [f for f in features if f["type"] in ("CDS", "gene")]

        if not cds_features:
            raise ValueError(
                f"No CDS features found in {snapgene_path}. "
                f"Available features: {[f['name'] for f in features]}"
            )

        if cds_feature_name is not None:
            matched = [f for f in cds_features if f["name"] == cds_feature_name]
            if not matched:
                available_names = [f["name"] for f in cds_features]
                raise ValueError(
                    f"CDS feature {cds_feature_name!r} not found. "
                    f"Available CDS features: {available_names}"
                )
            cds = matched[0]
        else:
            cds = cds_features[0]

        insert_length = cds["end"] - cds["start"]

        result = self.suggest(vector_name, insert_length)
        result["snapgene_info"] = {
            "file": snapgene_path,
            "cds_name": cds["name"],
            "cds_start": cds["start"],
            "cds_end": cds["end"],
            "insert_length_bp": insert_length,
            "is_circular": is_circular,
            "total_sequence_length": len(sequence) if sequence else 0,
        }
        return result

    @classmethod
    def register_primer_pair(
        cls, vector_name: str, forward: str, reverse: str,
    ) -> None:
        """VECTOR_PRIMER_MAP에 프라이머 쌍 등록/덮어쓰기.

        Parameters
        ----------
        vector_name : str
            벡터 canonical name (EXPRESSION_VECTORS의 key)
        forward : str
            Forward primer key (UNIVERSAL_PRIMERS의 key)
        reverse : str
            Reverse primer key (UNIVERSAL_PRIMERS의 key)
        """
        if forward not in UNIVERSAL_PRIMERS:
            raise ValueError(
                f"Forward primer {forward!r} not in UNIVERSAL_PRIMERS. "
                f"Available: {', '.join(UNIVERSAL_PRIMERS.keys())}"
            )
        if reverse not in UNIVERSAL_PRIMERS:
            raise ValueError(
                f"Reverse primer {reverse!r} not in UNIVERSAL_PRIMERS. "
                f"Available: {', '.join(UNIVERSAL_PRIMERS.keys())}"
            )
        VECTOR_PRIMER_MAP[vector_name] = (forward, reverse)

    # ── Internal helper ────────────────────────────────────────────────────

    def _pick_internal_primer(
        self,
        seq: str,
        anchor: int,
        direction: str,
        target_tm: float,
        min_len: int = 18,
        max_len: int = 28,
    ) -> tuple[str | None, float | None]:
        """Insert 서열 내부에서 target Tm에 가장 가까운 프라이머를 선택.

        Parameters
        ----------
        seq : str
            Insert 서열
        anchor : int
            프라이머 시작 위치 (0-indexed)
        direction : str
            "+" (forward) or "-" (reverse)
        target_tm : float
            목표 Tm
        min_len, max_len : int
            프라이머 길이 범위

        Returns
        -------
        tuple : (primer_seq, tm) or (None, None)
        """
        best_seq = None
        best_tm = None
        best_diff = float("inf")

        for length in range(min_len, max_len + 1):
            if direction == "+":
                if anchor + length > len(seq):
                    break
                candidate = seq[anchor:anchor + length]
            else:
                start = anchor - length
                if start < 0:
                    continue
                candidate = str(Seq(seq[start:anchor]).reverse_complement())

            tm = self.calc_tm(candidate)
            diff = abs(tm - target_tm)

            if diff < best_diff:
                best_diff = diff
                best_seq = candidate
                best_tm = tm

            # Tm을 초과하면 더 길게 할 필요 없음
            if tm >= target_tm:
                break

        return best_seq, best_tm


# ── Tests ──────────────────────────────────────────────────────────────────

def _run_tests():
    """ColonyPCRDesigner 테스트."""
    sep = "=" * 70
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

    designer = ColonyPCRDesigner()

    # ── Test 1: pET-28a(+) -> T7promoter / T7terminator ──────────────
    print(f"\n{sep}")
    print("  Test 1: pET-28a(+) primer mapping")
    print(sep)
    r1 = designer.suggest("pET-28a(+)", insert_length_bp=900)
    print(f"  Vector:   {r1['vector_name']}")
    print(f"  Forward:  {r1['f_name']} ({r1['f_seq']}) Tm={r1['f_tm']}C")
    print(f"  Reverse:  {r1['r_name']} ({r1['r_seq']}) Tm={r1['r_tm']}C")
    print(f"  Band (insert): {r1['expected_band_with_insert']} bp")
    print(f"  Band (empty):  {r1['expected_band_empty']} bp")
    print(f"  Anneal:   {r1['anneal_temp']}C")
    check("forward = T7promoter", r1["f_name"], "T7promoter")
    check("reverse = T7terminator", r1["r_name"], "T7terminator")
    check("band with insert = 200+900+150", r1["expected_band_with_insert"], 1250)
    print()

    # ── Test 2: pETDuet-1:MCS1 -> pET_upstream / T7terminator ────────
    print(f"{sep}")
    print("  Test 2: pETDuet-1:MCS1 primer mapping")
    print(sep)
    r2 = designer.suggest("pETDuet-1:MCS1", insert_length_bp=1200)
    print(f"  Forward:  {r2['f_name']}")
    print(f"  Reverse:  {r2['r_name']}")
    check("forward = pET-upstream", r2["f_name"], "pET-upstream")
    check("reverse = T7terminator", r2["r_name"], "T7terminator")
    print()

    # ── Test 3: Band size calculation ─────────────────────────────────
    print(f"{sep}")
    print("  Test 3: Band size calculation (known insert)")
    print(sep)
    r3 = designer.suggest("pET-21a(+)", insert_length_bp=750,
                           vector_flanking_bp=(180, 120))
    expected_band = 180 + 750 + 120
    print(f"  Flanking: {r3['flanking_upstream']} / {r3['flanking_downstream']}")
    print(f"  Expected band: {r3['expected_band_with_insert']} bp")
    check("band = 180+750+120", r3["expected_band_with_insert"], expected_band)
    check("flanking_upstream", r3["flanking_upstream"], 180)
    check("flanking_downstream", r3["flanking_downstream"], 120)
    print()

    # ── Test 4: Fuzzy vector name ─────────────────────────────────────
    print(f"{sep}")
    print("  Test 4: Fuzzy vector name ('pET28a' -> 'pET-28a(+)')")
    print(sep)
    r4 = designer.suggest("pET28a", insert_length_bp=600)
    print(f"  Input: 'pET28a' -> canonical: '{r4['vector_name']}'")
    check("canonical name", r4["vector_name"], "pET-28a(+)")
    check("forward = T7promoter", r4["f_name"], "T7promoter")
    print()

    # ── Test 5: suggest_with_internal() ───────────────────────────────
    print(f"{sep}")
    print("  Test 5: suggest_with_internal() basic test")
    print(sep)
    # 800bp 가상 insert 서열
    insert_seq = (
        "ATGCGTAACCTGGCGATCAAGCTGTTCGACGGTACCGATATCCTGCAGAAATTTGCGCCG"
        "GATCTGAACGAATGGCTGCACATCGGTCCTGCGATTGGCACCGATTTCAATCGCCTGATG"
        "CAGTTCGATGCATCGACCGGCTATCTGAACTCCGTCAAGGCGATGGACAAACTGCGCGGC"
        "GATACCGTGGAAATCGCGCAGCAGCTGGGCGATGAAGTGATCATCGATGCGTCCGGCAAA"
        "ATCGCGTTCAAAGGCACCGATACCGTGATGCTGAGCTATCCGGGCACGCCGGTCGATCCG"
        "GCGCTGACCGGCTGGCGCCTGTTTGAAACCGACAATGTCGATCTGGCAGTCAAAGCGCTG"
        "GGCCTGGATCACATCACCGGCGACTTTGCGGACTACCGCGATCTGACCAAACTGGATCTG"
        "GCGTCGATCTTCAACATCGCGAAAGAAGCGGGCATCCACGACGATACGCTGAGCGCAGTG"
        "ACCGATTTCTCCGGTCTGAACGATGCGACCAACGGCAACATCGCGCTGGCCCAGTTCATC"
        "GACACCAAAGACGGCACCGCGATCCTGACCAACGGCTCGCAGGGCAACAAGCTGACCGAT"
        "GCGATCATCAACGGCAAGACCATTCCGCTGAACGATCTGAACCTGACCGAAGCGGCGCAG"
        "GCGTTTGCGGAAGTGCTGAAAGGCGTCGATCCGGAAACCATCACCCTGGAAGCGATCCGC"
        "GACGGCAAGATCAACCTGCAGGATCGCGTGATCGCGATGGGCGATACCCTGCAGGGCAAC"
        "GCGTCGATCACCGCGTAA"
    )
    r5 = designer.suggest_with_internal("pET-28a(+)", insert_seq)
    print(f"  Universal F: {r5['universal']['f_name']}")
    print(f"  Universal R: {r5['universal']['r_name']}")
    print(f"  Internal F:  {r5['internal_f_seq']}  Tm={r5['internal_f_tm']}C")
    print(f"  Internal R:  {r5['internal_r_seq']}  Tm={r5['internal_r_tm']}C")
    print(f"  Internal product: {r5['internal_product_size']} bp")
    if r5['internal_notes']:
        for n in r5['internal_notes']:
            print(f"  Note: {n}")
    check("internal_f_seq is not None", r5["internal_f_seq"] is not None, True)
    check("internal_r_seq is not None", r5["internal_r_seq"] is not None, True)
    check("internal_product_size > 0",
          r5["internal_product_size"] is not None and r5["internal_product_size"] > 0,
          True)
    print()

    # ── Test 6: register_primer_pair() ────────────────────────────────
    print(f"{sep}")
    print("  Test 6: register_primer_pair()")
    print(sep)
    # 원래 pET-21a(+) 매핑 확인
    orig_f, orig_r = VECTOR_PRIMER_MAP["pET-21a(+)"]
    check("original forward", orig_f, "T7promoter")

    # 새 매핑 등록
    ColonyPCRDesigner.register_primer_pair("pET-21a(+)", "pET_upstream", "pET_RP")
    new_f, new_r = VECTOR_PRIMER_MAP["pET-21a(+)"]
    check("updated forward", new_f, "pET_upstream")
    check("updated reverse", new_r, "pET_RP")

    r6 = designer.suggest("pET-21a(+)", insert_length_bp=500)
    check("suggest uses updated mapping", r6["f_name"], "pET-upstream")

    # 원래대로 복원
    ColonyPCRDesigner.register_primer_pair("pET-21a(+)", "T7promoter", "T7terminator")
    restored_f, _ = VECTOR_PRIMER_MAP["pET-21a(+)"]
    check("restored forward", restored_f, "T7promoter")

    # 잘못된 프라이머 키 테스트
    try:
        ColonyPCRDesigner.register_primer_pair("pET-21a(+)", "INVALID", "T7terminator")
        failed += 1
        print("  FAIL: should raise ValueError for invalid primer key")
    except ValueError:
        passed += 1
        print("  OK: raises ValueError for invalid primer key")
    print()

    # ── Summary ──────────────────────────────────────────────────────
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
