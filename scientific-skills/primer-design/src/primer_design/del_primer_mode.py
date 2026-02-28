#!/usr/bin/env python3
"""
iPCR Deletion Primer Designer
===============================
iPCRDesignerBase를 상속하여 결실(deletion) 프라이머 설계.

설계 원리 (back-to-back overlapping):
  원본: ...AAATTT [DEL_REGION] GGGCCC...
  F: 5'-[left_overlap (k1 bp)]-[annealing_downstream]-3'
  R: 5'-[RC(right_overlap) (k2 bp)]-[annealing_upstream]-3'
  overlap = left_before_del + right_after_del (결실 영역 건너뜀)

  Tm 계산: annealing 부분만 (tail은 결실 영역 건너편이라 초기 cycle에서
  template과 연속 결합 불가).

Frameshift 판정:
  결실 길이 % 3 != 0 → FRAMESHIFT 경고
  1 bp 결실 시 코돈 내 위치까지 판단 (cds_start 제공 시)
"""

from Bio.Seq import Seq

from .subst_primer_mode import iPCRDesignerBase


class iPCRDelDesigner(iPCRDesignerBase):
    """결실 프라이머 설계 (back-to-back overlapping + Tm 기반 품질 판정)."""

    def design(self, seq, del_start, del_end,
               target_tm=61.0, overlap_len=18, min_len=18, max_len=35,
               cds_start=None):
        """Deletion primer design.

        Parameters
        ----------
        seq : str           template 전체 서열
        del_start : int     결실 시작 위치 (0-indexed, inclusive)
        del_end : int       결실 끝 위치 (0-indexed, exclusive)
        target_tm : float   annealing 목표 Tm (degC)
        overlap_len : int   overlap 총 길이 (bp)
        min_len : int       annealing 최소 길이
        max_len : int       annealing 최대 길이
        cds_start : int     reading frame 기준 CDS 시작 위치 (선택)

        Returns
        -------
        dict : f_full, r_full, del_seq, frameshift_warning, f_qc, r_qc, ...
        """
        seq = seq.upper()
        warnings = []

        # 1. 입력 검증
        if del_start < 0 or del_end > len(seq) or del_start >= del_end:
            raise ValueError(
                f"Invalid deletion range: [{del_start}, {del_end}) "
                f"for sequence of length {len(seq)}")

        del_seq = seq[del_start:del_end]
        del_len = del_end - del_start

        # 2. Frameshift 판단
        frameshift_warning = False
        if del_len % 3 == 0:
            warnings.append(
                f"In-frame deletion ({del_len} bp = {del_len // 3} codons removed)")
        else:
            frameshift_warning = True
            if del_len == 1 and cds_start is not None:
                codon_pos = (del_start - cds_start) % 3
                if codon_pos == 2:
                    warnings.append(
                        "FRAMESHIFT: 1 bp deletion at codon 3rd position "
                        "- silent change possible but reading frame shifts")
                else:
                    pos_label = {0: "1st", 1: "2nd"}[codon_pos]
                    warnings.append(
                        f"FRAMESHIFT: 1 bp deletion at codon {pos_label} position "
                        f"- reading frame disrupted")
            elif del_len == 1:
                warnings.append(
                    "FRAMESHIFT: 1 bp deletion "
                    "- provide cds_start to check codon position")
            else:
                warnings.append(
                    f"FRAMESHIFT: {del_len} bp deletion (not multiple of 3) "
                    f"- reading frame disrupted")

        # 3. Overlap 분할
        k1 = overlap_len // 2
        k2 = overlap_len - k1

        if del_start - k1 < 0:
            k1 = del_start
            k2 = overlap_len - k1
            warnings.append(
                f"Overlap adjusted: left side limited to {k1} bp (near sequence start)")
        if del_end + k2 > len(seq):
            k2 = len(seq) - del_end
            warnings.append(
                f"Overlap adjusted: right side limited to {k2} bp (near sequence end)")

        left_overlap = seq[del_start - k1: del_start]
        right_overlap = seq[del_end: del_end + k2]
        overlap_region = left_overlap + right_overlap

        # 4. F primer: left_overlap(5' tail) + annealing downstream
        f_ann, f_tm, f_gc, f_alen = self._design_annealing(
            seq, del_end, "+", target_tm, min_len, max_len)
        if f_ann is None:
            raise RuntimeError("Forward primer annealing design failed")
        f_full = left_overlap + f_ann

        # 5. R primer: RC(right_overlap)(5' tail) + annealing upstream
        r_ann, r_tm, r_gc, r_alen = self._design_annealing(
            seq, del_start, "-", target_tm, min_len, max_len)
        if r_ann is None:
            raise RuntimeError("Reverse primer annealing design failed")
        right_overlap_rc = str(Seq(right_overlap).reverse_complement())
        r_full = right_overlap_rc + r_ann

        # 6. Overlap 검증
        f_5prime = f_full[:k1]
        r_5prime = r_full[:k2]
        reconstructed = f_5prime + str(Seq(r_5prime).reverse_complement())
        overlap_verified = (reconstructed == overlap_region)
        if not overlap_verified:
            warnings.append("Overlap verification FAILED - check primer design")

        # 7. Anneal temp & quality check
        anneal_temp = min(f_tm, r_tm) + 1.0
        f_qc = self.check_primer(f_ann, anneal_temp)
        r_qc = self.check_primer(r_ann, anneal_temp)
        het = self.check_heterodimer(f_full, r_full)

        # 8. Primer-level warnings
        for label, primer, ann in [("F", f_full, f_ann), ("R", r_full, r_ann)]:
            hp = self.homopolymer_run(ann)
            if hp:
                warnings.append(f"{label} primer: homopolymer run ({hp})")
            if not self.gc_clamp_ok(primer):
                warnings.append(f"{label} primer: no G/C clamp at 3' end")

        return {
            "f_full": f_full, "r_full": r_full,
            "f_ann": f_ann, "r_ann": r_ann,
            "f_tm": round(f_tm, 1), "r_tm": round(r_tm, 1),
            "f_gc": round(f_gc, 1), "r_gc": round(r_gc, 1),
            "f_len": len(f_full), "r_len": len(r_full),
            "f_ann_len": f_alen, "r_ann_len": r_alen,
            "del_seq": del_seq, "del_len": del_len,
            "frameshift_warning": frameshift_warning,
            "overlap_region": overlap_region,
            "overlap_verified": overlap_verified,
            "anneal_temp": round(anneal_temp, 1),
            "k1": k1, "k2": k2,
            "f_qc": f_qc, "r_qc": r_qc, "het": het,
            "warnings": warnings,
        }

    @staticmethod
    def print_result(result):
        """결과 출력 포맷터."""
        sep = "=" * 70
        line = "-" * 70
        k1, k2 = result['k1'], result['k2']

        print(sep)
        print("  iPCR Deletion Primer Designer")
        print(sep)
        print(f"\n  Deletion : {result['del_len']} bp")
        print(f"  Deleted  : 5'-{result['del_seq']}-3'")
        if result['frameshift_warning']:
            print(f"  *** FRAMESHIFT WARNING ***")

        if result['warnings']:
            print(f"\n  Warnings:")
            for w in result['warnings']:
                print(f"    - {w}")

        print(f"\n  Overlap: {result['overlap_region']} "
              f"({len(result['overlap_region'])} bp = left {k1}bp + right {k2}bp)")

        # F primer
        print(f"\n{line}")
        f_qc = result['f_qc']
        print(f"  Forward  5'-{result['f_full']}-3'  "
              f"({result['f_len']} nt, Tm={result['f_tm']}C)  "
              f"QC: {f_qc['verdict']}")
        print(f"    tail: {result['f_full'][:k1]} ({k1} nt)  "
              f"ann: {result['f_ann']} ({result['f_ann_len']} nt)")

        # R primer
        r_qc = result['r_qc']
        print(f"  Reverse  5'-{result['r_full']}-3'  "
              f"({result['r_len']} nt, Tm={result['r_tm']}C)  "
              f"QC: {r_qc['verdict']}")
        print(f"    tail: {result['r_full'][:k2]} ({k2} nt)  "
              f"ann: {result['r_ann']} ({result['r_ann_len']} nt)")

        print(f"\n  Anneal temp: {result['anneal_temp']}C")
        print(f"  Overlap verified: "
              f"{'YES' if result['overlap_verified'] else 'NO'}")
        het = result['het']
        print(f"  Heterodimer: dG={het.get('dg','N/A')} kcal/mol")
        print()


# ── 테스트 ──────────────────────────────────────────────────────────────

def _run_tests():
    """iPCRDelDesigner 테스트."""
    sep = "=" * 70
    designer = iPCRDelDesigner()

    test_seq = (
        "ATGAAAGCTGCCATTGTTCTGAGCGAATCCGGTGTGCACGAATGTCCGAAAGAATCCAAACTG"
        "GTTCCGGCACTGAACGGTCTGGAAATCGAAGATGAAGGCGTTATCCCGGAATTCTTCAAGGGC"
    )

    print(f"\n{sep}\n  DELETION PRIMER TESTS\n{sep}")
    print(f"  Test sequence: {len(test_seq)} bp")

    # Test 1: 1 bp deletion (frameshift)
    print(f"\n  Test 1: 1 bp deletion at pos 30 (frameshift)")
    r1 = designer.design(test_seq, del_start=30, del_end=31, cds_start=0)
    designer.print_result(r1)
    assert r1['frameshift_warning'] and r1['overlap_verified']
    print("  -> PASSED")

    # Test 2: 3 bp deletion (in-frame)
    print(f"\n  Test 2: 3 bp deletion at pos 30-33 (in-frame)")
    r2 = designer.design(test_seq, del_start=30, del_end=33, cds_start=0)
    assert not r2['frameshift_warning'] and r2['overlap_verified']
    print("  -> PASSED")

    # Test 3: 6 bp deletion (in-frame, 2 codons)
    print(f"\n  Test 3: 6 bp deletion at pos 30-36 (in-frame)")
    r3 = designer.design(test_seq, del_start=30, del_end=36, cds_start=0)
    assert not r3['frameshift_warning'] and r3['overlap_verified']
    print("  -> PASSED")

    # Test 4: 4 bp deletion (frameshift)
    print(f"\n  Test 4: 4 bp deletion at pos 30-34 (frameshift)")
    r4 = designer.design(test_seq, del_start=30, del_end=34, cds_start=0)
    assert r4['frameshift_warning'] and r4['overlap_verified']
    print("  -> PASSED")

    # Test 5: 1 bp deletion at codon 3rd position
    print(f"\n  Test 5: 1 bp deletion at codon 3rd pos (pos 32)")
    r5 = designer.design(test_seq, del_start=32, del_end=33, cds_start=0)
    assert r5['frameshift_warning']
    assert any("3rd position" in w for w in r5['warnings'])
    print("  -> PASSED")

    print(f"\n{sep}\n  ALL TESTS PASSED\n{sep}")


if __name__ == "__main__":
    _run_tests()
