#!/usr/bin/env python3
"""
Vector Registry & Reading Frame Checker
========================================
자주 사용하는 발현 벡터의 MCS 서열을 저장하고,
RE 클로닝 전략의 reading frame 호환성을 자동 검증.

Ligation 모델:
  Sticky-end ligation 후 양쪽 RE recognition site이 완전히 복원됨.
  따라서 reading frame은 cut position이 아닌 site position 기준으로 계산.

  5' junction: mcs[0 : site5 + len(RE5)] + INSERT_CDS
  3' junction: INSERT_CDS + mcs[site3 : stop]

Data sources:
  pET series — Novagen TB055 (pET System Manual, 11th Edition)
  pMAL-c6T  — NEB #E8202 (manual)
  Duet series — Novagen TB340 (Duet Vectors manual)
"""

from __future__ import annotations

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
    """DNA → amino acid string (standard code). Stop = '*'."""
    aa = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3].upper()
        aa.append(_CODON_TABLE.get(codon, "?"))
    return "".join(aa)


# ── Restriction Enzymes ─────────────────────────────────────────────────────

RESTRICTION_ENZYMES: dict[str, dict] = {
    "NheI":    {"recognition": "GCTAGC", "cut_top": 1, "cut_bottom": 5,
                "overhang": "5prime", "overhang_len": 4},
    "NcoI":    {"recognition": "CCATGG", "cut_top": 1, "cut_bottom": 5,
                "overhang": "5prime", "overhang_len": 4},
    "NdeI":    {"recognition": "CATATG", "cut_top": 2, "cut_bottom": 4,
                "overhang": "5prime", "overhang_len": 2},
    "BamHI":   {"recognition": "GGATCC", "cut_top": 1, "cut_bottom": 5,
                "overhang": "5prime", "overhang_len": 4},
    "BglII":   {"recognition": "AGATCT", "cut_top": 1, "cut_bottom": 5,
                "overhang": "5prime", "overhang_len": 4},
    "EcoRI":   {"recognition": "GAATTC", "cut_top": 1, "cut_bottom": 5,
                "overhang": "5prime", "overhang_len": 4},
    "EcoRV":   {"recognition": "GATATC", "cut_top": 3, "cut_bottom": 3,
                "overhang": "blunt", "overhang_len": 0},
    "MfeI":    {"recognition": "CAATTG", "cut_top": 1, "cut_bottom": 5,
                "overhang": "5prime", "overhang_len": 4},
    "SacI":    {"recognition": "GAGCTC", "cut_top": 5, "cut_bottom": 1,
                "overhang": "3prime", "overhang_len": 4},
    "SalI":    {"recognition": "GTCGAC", "cut_top": 1, "cut_bottom": 5,
                "overhang": "5prime", "overhang_len": 4},
    "PstI":    {"recognition": "CTGCAG", "cut_top": 5, "cut_bottom": 1,
                "overhang": "3prime", "overhang_len": 4},
    "HindIII": {"recognition": "AAGCTT", "cut_top": 1, "cut_bottom": 5,
                "overhang": "5prime", "overhang_len": 4},
    "NotI":    {"recognition": "GCGGCCGC", "cut_top": 2, "cut_bottom": 6,
                "overhang": "5prime", "overhang_len": 4},
    "XhoI":    {"recognition": "CTCGAG", "cut_top": 1, "cut_bottom": 5,
                "overhang": "5prime", "overhang_len": 4},
    "AvrII":   {"recognition": "CCTAGG", "cut_top": 1, "cut_bottom": 5,
                "overhang": "5prime", "overhang_len": 4},
    "KpnI":    {"recognition": "GGTACC", "cut_top": 5, "cut_bottom": 1,
                "overhang": "3prime", "overhang_len": 4},
    "FseI":    {"recognition": "GGCCGGCC", "cut_top": 6, "cut_bottom": 2,
                "overhang": "3prime", "overhang_len": 4},
    "AscI":    {"recognition": "GGCGCGCC", "cut_top": 2, "cut_bottom": 6,
                "overhang": "5prime", "overhang_len": 4},
    "PacI":    {"recognition": "TTAATTAA", "cut_top": 5, "cut_bottom": 3,
                "overhang": "3prime", "overhang_len": 2},
}


# ── Vector MCS Data ─────────────────────────────────────────────────────────
# 각 벡터의 MCS는 ATG(start) ~ Stop codon (inclusive)
# reading frame 0 = ATG의 A가 position 0
# re_sites: recognition seq 첫 뉴클레오타이드의 0-indexed position

def _verify_re_sites(mcs: str, re_sites: dict[str, int], label: str) -> None:
    """MCS 내 모든 RE site 위치 검증 (import 시 assertion)."""
    for enzyme, expected_pos in re_sites.items():
        rec = RESTRICTION_ENZYMES[enzyme]["recognition"]
        actual = mcs.find(rec)
        assert actual == expected_pos, (
            f"{label} {enzyme}: expected pos {expected_pos}, found {actual} "
            f"(recognition={rec})"
        )


# ── pET-21a(+) ──────────────────────────────────────────────────────────────
# T7-tag (N-term) + C-terminal His6
# Source: Novagen TB055
#
# ATG GCT AGC ATG ACT GGT GGA CAG CAA ATG GGT CGC GGA TCC GAA TTC
#  M   A   S   M   T   G   G   Q   Q   M   G   R   G   S   E   F
#  0   3   6   9  12  15  18  21  24  27  30  33  36  39  42  45
#                NheI(3)                         BamHI(36) EcoRI(42)
#
# GAG CTC GTC GAC AAG CTT GCG GCC GCA CTC GAG CAC CAC CAC CAC CAC CAC TGA
#  E   L   V   D   K   L   A   A   A   L   E   H   H   H   H   H   H   *
# 48  51  54  57  60  63  66  69  72  75  78  81  84  87  90  93  96  99
# SacI(48) SalI(54) HindIII(60) NotI(66)    XhoI(75)  C-His6(81-96)   Stop

_PET21A_MCS = (
    "ATG"   # 0   M  — start
    "GCT"   # 3   A  — NheI: GCTAGC starts at 3
    "AGC"   # 6   S
    "ATG"   # 9   M
    "ACT"   # 12  T
    "GGT"   # 15  G
    "GGA"   # 18  G
    "CAG"   # 21  Q
    "CAA"   # 24  Q
    "ATG"   # 27  M
    "GGT"   # 30  G
    "CGC"   # 33  R
    "GGA"   # 36  G  — BamHI: GGATCC starts at 36
    "TCC"   # 39  S
    "GAA"   # 42  E  — EcoRI: GAATTC starts at 42
    "TTC"   # 45  F
    "GAG"   # 48  E  — SacI: GAGCTC starts at 48
    "CTC"   # 51  L
    "GTC"   # 54  V  — SalI: GTCGAC starts at 54
    "GAC"   # 57  D
    "AAG"   # 60  K  — HindIII: AAGCTT starts at 60
    "CTT"   # 63  L
    "GCG"   # 66  A  — NotI: GCGGCCGC starts at 66
    "GCC"   # 69  A
    "GCA"   # 72  A
    "CTC"   # 75  L  — XhoI: CTCGAG starts at 75
    "GAG"   # 78  E
    "CAC"   # 81  H  — C-His6 start
    "CAC"   # 84  H
    "CAC"   # 87  H
    "CAC"   # 90  H
    "CAC"   # 93  H
    "CAC"   # 96  H
    "TGA"   # 99  *  — stop
)

_PET21A_RE_SITES = {
    "NheI": 3, "BamHI": 36, "EcoRI": 42, "SacI": 48,
    "SalI": 54, "HindIII": 60, "NotI": 66, "XhoI": 75,
}
_verify_re_sites(_PET21A_MCS, _PET21A_RE_SITES, "pET-21a(+)")


# ── pET-28a(+) ──────────────────────────────────────────────────────────────
# N-terminal His6 + Thrombin + T7-tag + MCS + C-terminal His6
# Source: Novagen TB055
#
# ATG GGC AGC AGC CAT CAT CAT CAT CAT CAC AGC AGC GGC CTG GTG CCG CGC GGC AGC
#  M   G   S   S   H   H   H   H   H   H   S   S   G   L   V   P   R   G   S
#  0   3   6   9  12  15  18  21  24  27  30  33  36  39  42  45  48  51  54
#  start         N-His6(12-27)                       Thrombin(LVPR|GS,39-54)
#
# CAT ATG GCT AGC ATG ACT GGT GGA CAG CAA ATG GGT CGC GGA TCC GAA TTC
#  H   M   A   S   M   T   G   G   Q   Q   M   G   R   G   S   E   F
# 57  60  63  66  69  72  75  78  81  84  87  90  93  96  99 102 105
# NdeI(57) NheI(63)  T7-tag(69-92)                BamHI(96) EcoRI(102)
#
# GAG CTC GTC GAC AAG CTT GCG GCC GCA CTC GAG CAC CAC CAC CAC CAC CAC TGA
#  E   L   V   D   K   L   A   A   A   L   E   H   H   H   H   H   H   *
# 108 111 114 117 120 123 126 129 132 135 138 141 144 147 150 153 156 159
# SacI(108)SalI(114)HindIII(120)NotI(126)  XhoI(135) C-His6(141-156) Stop

_PET28A_MCS = (
    "ATG"   # 0   M  — start
    "GGC"   # 3   G
    "AGC"   # 6   S
    "AGC"   # 9   S
    "CAT"   # 12  H  — N-His6 start
    "CAT"   # 15  H
    "CAT"   # 18  H
    "CAT"   # 21  H
    "CAT"   # 24  H
    "CAC"   # 27  H  — N-His6 end
    "AGC"   # 30  S
    "AGC"   # 33  S
    "GGC"   # 36  G
    "CTG"   # 39  L  — Thrombin start (LVPR↓GS)
    "GTG"   # 42  V
    "CCG"   # 45  P
    "CGC"   # 48  R
    "GGC"   # 51  G
    "AGC"   # 54  S  — Thrombin end
    "CAT"   # 57  H  — NdeI: CATATG starts at 57
    "ATG"   # 60  M
    "GCT"   # 63  A  — NheI: GCTAGC starts at 63
    "AGC"   # 66  S
    "ATG"   # 69  M  — T7-tag start
    "ACT"   # 72  T
    "GGT"   # 75  G
    "GGA"   # 78  G
    "CAG"   # 81  Q
    "CAA"   # 84  Q
    "ATG"   # 87  M
    "GGT"   # 90  G
    "CGC"   # 93  R  — T7-tag end
    "GGA"   # 96  G  — BamHI: GGATCC starts at 96
    "TCC"   # 99  S
    "GAA"   # 102 E  — EcoRI: GAATTC starts at 102
    "TTC"   # 105 F
    "GAG"   # 108 E  — SacI: GAGCTC starts at 108
    "CTC"   # 111 L
    "GTC"   # 114 V  — SalI: GTCGAC starts at 114
    "GAC"   # 117 D
    "AAG"   # 120 K  — HindIII: AAGCTT starts at 120
    "CTT"   # 123 L
    "GCG"   # 126 A  — NotI: GCGGCCGC starts at 126
    "GCC"   # 129 A
    "GCA"   # 132 A
    "CTC"   # 135 L  — XhoI: CTCGAG starts at 135
    "GAG"   # 138 E
    "CAC"   # 141 H  — C-His6 start
    "CAC"   # 144 H
    "CAC"   # 147 H
    "CAC"   # 150 H
    "CAC"   # 153 H
    "CAC"   # 156 H
    "TGA"   # 159 *  — stop
)

_PET28A_RE_SITES = {
    "NdeI": 57, "NheI": 63, "BamHI": 96, "EcoRI": 102, "SacI": 108,
    "SalI": 114, "HindIII": 120, "NotI": 126, "XhoI": 135,
}
_verify_re_sites(_PET28A_MCS, _PET28A_RE_SITES, "pET-28a(+)")


# ── pMAL-c6T ────────────────────────────────────────────────────────────────
# NEB #E8202: His6-MBP-TEV fusion, cytoplasmic expression
# Source: NEB pMAL-c6T map, reconstructed from documentation
#
# 구조: tac promoter → malE (His6-MBP) → TEV site → polylinker → stop
# 아래 MCS는 TEV cleavage site의 G 잔기(frame 0)부터 stop codon까지.
# 실제 start codon은 malE ATG (MCS 밖, 훨씬 upstream).
#
# GGC ATG CTG ATG GGC GGC CGC GAT ATC GTC GAC GGA TCC GAA TTC CCT GCA GGT AAT AAG CTT TAA
#  G   M   L   M   G   G   R   D   I   V   D   G   S   E   F   P   A   G   N   K   L   *
#  0   3   6   9  12  15  18  21  24  27  30  33  36  39  42  45  48  51  54  57  60  63
#                     NotI(13)  EcoRV(21) SalI(27) BamHI(33) EcoRI(39) PstI(46) HindIII(57)
#
# NOTE: Position 0 (GGC) = TEV cleavage 후 남는 G 잔기, malE와 동일 reading frame.
#       Stop codon (TAA) 이후 추가 stop이 다른 frame에도 존재 (모든 frame 종결).
#       PstI (pos 46)는 codon 경계를 넘으므로 주의.

_PMALC6T_MCS = (
    "GGC"   # 0   G  — TEV cleavage product (frame 0 reference)
    "ATG"   # 3   M
    "CTG"   # 6   L
    "ATG"   # 9   M
    "GGC"   # 12  G
    "GGC"   # 15  G  — NotI: GCGGCCGC starts at 13
    "CGC"   # 18  R
    "GAT"   # 21  D  — EcoRV: GATATC starts at 21
    "ATC"   # 24  I
    "GTC"   # 27  V  — SalI: GTCGAC starts at 27
    "GAC"   # 30  D
    "GGA"   # 33  G  — BamHI: GGATCC starts at 33
    "TCC"   # 36  S
    "GAA"   # 39  E  — EcoRI: GAATTC starts at 39
    "TTC"   # 42  F
    "CCT"   # 45  P  — PstI: CTGCAG starts at 46 (cross-codon)
    "GCA"   # 48  A
    "GGT"   # 51  G
    "AAT"   # 54  N
    "AAG"   # 57  K  — HindIII: AAGCTT starts at 57
    "CTT"   # 60  L
    "TAA"   # 63  *  — stop
)

_PMALC6T_RE_SITES = {
    "NotI": 13, "EcoRV": 21, "SalI": 27, "BamHI": 33,
    "EcoRI": 39, "PstI": 46, "HindIII": 57,
}
_verify_re_sites(_PMALC6T_MCS, _PMALC6T_RE_SITES, "pMAL-c6T")


# ── pETDuet-1 MCS1 ──────────────────────────────────────────────────────────
# Novagen TB340: N-terminal His6 + MCS1 (co-expression vector)
# Source: EcoliWiki (ecoliwiki.org/colipedia/index.php/pETDuet-1)
#
# ATG GGC AGC AGC CAT CAC CAT CAT CAC CAC AGC CAG GAT CCG AAT TCG
#  M   G   S   S   H   H   H   H   H   H   S   Q   D   P   N   S
#  0   3   6   9  12  15  18  21  24  27  30  33  36  39  42  45
#                  N-His6(12-29)           BamHI(35)     EcoRI(41)
#
# AGC TCG GCG CGC CTG CAG GTC GAC AAG CTT GCG GCC GCA TAA
#  S   S   A   R   L   Q   V   D   K   L   A   A   A   *
# 48  51  54  57  60  63  66  69  72  75  78  81  84  87
# SacI(47) AscI(53)PstI(60)SalI(66)HindIII(72)NotI(78)Stop
#
# WARNING: BamHI/EcoRI/SacI는 reading frame +2 offset (His6와 out-of-frame).
#          NcoI는 start codon과 겹침 (MCS 외부, position -2).
#          실험 전 반드시 TB340 원본 또는 SnapGene 파일로 확인할 것.

_PETDUET1_MCS1 = (
    "ATG"   # 0   M  — start (within NcoI site CCATGG)
    "GGC"   # 3   G
    "AGC"   # 6   S
    "AGC"   # 9   S
    "CAT"   # 12  H  — N-His6 start
    "CAC"   # 15  H
    "CAT"   # 18  H
    "CAT"   # 21  H
    "CAC"   # 24  H
    "CAC"   # 27  H  — N-His6 end
    "AGC"   # 30  S
    "CAG"   # 33  Q
    "GAT"   # 36  D  — BamHI: GGATCC at 35 (cross-codon, frame +2)
    "CCG"   # 39  P
    "AAT"   # 42  N  — EcoRI: GAATTC at 41 (cross-codon)
    "TCG"   # 45  S
    "AGC"   # 48  S  — SacI: GAGCTC at 47 (cross-codon)
    "TCG"   # 51  S
    "GCG"   # 54  A  — AscI: GGCGCGCC at 53 (cross-codon)
    "CGC"   # 57  R
    "CTG"   # 60  L  — PstI: CTGCAG at 60 (codon-aligned)
    "CAG"   # 63  Q
    "GTC"   # 66  V  — SalI: GTCGAC at 66 (codon-aligned)
    "GAC"   # 69  D
    "AAG"   # 72  K  — HindIII: AAGCTT at 72 (codon-aligned)
    "CTT"   # 75  L
    "GCG"   # 78  A  — NotI: GCGGCCGC at 78 (codon-aligned)
    "GCC"   # 81  A
    "GCA"   # 84  A
    "TAA"   # 87  *  — stop
)

_PETDUET1_MCS1_RE_SITES = {
    "BamHI": 35, "EcoRI": 41, "SacI": 47, "AscI": 53,
    "PstI": 60, "SalI": 66, "HindIII": 72, "NotI": 78,
}
_verify_re_sites(_PETDUET1_MCS1, _PETDUET1_MCS1_RE_SITES, "pETDuet-1:MCS1")


# ── Duet MCS2 (shared: pETDuet-1, pACYCDuet-1) ───────────────────────────
# Novagen TB340: NdeI(start) + cloning sites + S-tag + stop
# Source: NovoPro GenBank (V11042) / SnapGene / Novagen TB340
#
# ATG GCA GAT CTC AAT TGG ATA TCG GCC GGC CAC GCG ATC GCT GAC GTC
#  M   A   D   L   N   W   I   S   A   G   H   A   I   A   D   V
#  0   3   6   9  12  15  18  21  24  27  30  33  36  39  42  45
#       BglII(5)  MfeI(11) EcoRV(17) FseI(23)
#
# GGT ACC CTC GAG TCT GGT AAA GAA ACC GCT GCT GCG AAA TTT GAA CGC
#  G   T   L   E   S   G   K   E   T   A   A   A   K   F   E   R
# 48  51  54  57  60  63  66  69  72  75  78  81  84  87  90  93
# KpnI(48) XhoI(54)       S-tag(66-108): KETAAAKFERQHMDS
#
# CAG CAC ATG GAC TCG TCT ACT AGC GCA GCT TAA
#  Q   H   M   D   S   S   T   S   A   A   *
# 96  99 102 105 108 111 114 117 120 123 126
#                   ↑S-tag end              stop
#
# NOTE: NdeI (CATATG) spans position -2 to 3, overlapping start codon.
#       BglII/MfeI/EcoRV/FseI are cross-codon (frame +2).
#       Only KpnI/XhoI are codon-aligned.
#       SgfI(GCGATCGC, pos 33) and AatII(GACGTC, pos 42) also present
#       but not included in RESTRICTION_ENZYMES dict.

_DUET_MCS2 = (
    "ATG"   # 0   M  -- start (within NdeI site CATATG)
    "GCA"   # 3   A
    "GAT"   # 6   D  -- BglII: AGATCT at 5 (cross-codon)
    "CTC"   # 9   L
    "AAT"   # 12  N  -- MfeI: CAATTG at 11 (cross-codon)
    "TGG"   # 15  W
    "ATA"   # 18  I  -- EcoRV: GATATC at 17 (cross-codon)
    "TCG"   # 21  S
    "GCC"   # 24  A  -- FseI: GGCCGGCC at 23 (cross-codon)
    "GGC"   # 27  G
    "CAC"   # 30  H
    "GCG"   # 33  A
    "ATC"   # 36  I
    "GCT"   # 39  A
    "GAC"   # 42  D
    "GTC"   # 45  V
    "GGT"   # 48  G  -- KpnI: GGTACC at 48 (codon-aligned)
    "ACC"   # 51  T
    "CTC"   # 54  L  -- XhoI: CTCGAG at 54 (codon-aligned)
    "GAG"   # 57  E
    "TCT"   # 60  S
    "GGT"   # 63  G
    "AAA"   # 66  K  -- S-tag start (KETAAAKFERQHMDS)
    "GAA"   # 69  E
    "ACC"   # 72  T
    "GCT"   # 75  A
    "GCT"   # 78  A
    "GCG"   # 81  A
    "AAA"   # 84  K
    "TTT"   # 87  F
    "GAA"   # 90  E
    "CGC"   # 93  R
    "CAG"   # 96  Q
    "CAC"   # 99  H
    "ATG"   # 102 M
    "GAC"   # 105 D
    "TCG"   # 108 S  -- S-tag end
    "TCT"   # 111 S
    "ACT"   # 114 T
    "AGC"   # 117 S
    "GCA"   # 120 A
    "GCT"   # 123 A
    "TAA"   # 126 *  -- stop
)

_DUET_MCS2_RE_SITES = {
    "BglII": 5, "MfeI": 11, "EcoRV": 17, "FseI": 23,
    "KpnI": 48, "XhoI": 54,
}
_verify_re_sites(_DUET_MCS2, _DUET_MCS2_RE_SITES, "Duet:MCS2")


# ── Vector Registry ─────────────────────────────────────────────────────────

EXPRESSION_VECTORS: dict[str, dict] = {
    "pET-21a(+)": {
        "mcs_seq": _PET21A_MCS,
        "re_sites": _PET21A_RE_SITES,
        "tags": {
            "T7-tag":  {"type": "N-terminal", "region": (0, 33)},
            "C-His6":  {"type": "C-terminal", "region": (81, 96)},
        },
        "stop_pos": 99,
        "aliases": ["pET21a", "pET-21a", "pet21a", "pet-21a(+)"],
    },
    "pET-28a(+)": {
        "mcs_seq": _PET28A_MCS,
        "re_sites": _PET28A_RE_SITES,
        "tags": {
            "N-His6":    {"type": "N-terminal", "region": (12, 27)},
            "Thrombin":  {"type": "N-terminal", "region": (39, 54)},
            "T7-tag":    {"type": "N-terminal", "region": (69, 93)},
            "C-His6":    {"type": "C-terminal", "region": (141, 156)},
        },
        "stop_pos": 159,
        "aliases": ["pET28a", "pET-28a", "pet28a", "pet-28a(+)"],
    },
    "pMAL-c6T": {
        "mcs_seq": _PMALC6T_MCS,
        "re_sites": _PMALC6T_RE_SITES,
        "tags": {},  # MBP, His6 are upstream of polylinker (in malE)
        "stop_pos": 63,
        "aliases": ["pMALc6T", "pmal-c6t", "pmalc6t", "pMAL c6T"],
        "notes": "MCS position 0 = TEV cleavage product (G), not ATG. "
                 "malE ATG is far upstream.",
    },
    "pETDuet-1:MCS1": {
        "mcs_seq": _PETDUET1_MCS1,
        "re_sites": _PETDUET1_MCS1_RE_SITES,
        "tags": {
            "N-His6": {"type": "N-terminal", "region": (12, 27)},
        },
        "stop_pos": 87,
        "aliases": ["pETDuet1", "pETDuet-1", "petduet1", "petduet-1",
                     "pETDuet1:MCS1"],
        "notes": "BamHI/EcoRI/SacI are cross-codon (frame +2). "
                 "Only PstI/SalI/HindIII/NotI are codon-aligned.",
    },
    "pETDuet-1:MCS2": {
        "mcs_seq": _DUET_MCS2,
        "re_sites": _DUET_MCS2_RE_SITES,
        "tags": {
            "S-tag": {"type": "C-terminal", "region": (66, 108)},
        },
        "stop_pos": 126,
        "aliases": ["pETDuet1:MCS2"],
        "notes": "NdeI at position -2 (overlaps start ATG). "
                 "BglII/MfeI/EcoRV/FseI are cross-codon (frame +2). "
                 "Only KpnI/XhoI are codon-aligned.",
    },
    "pACYCDuet-1:MCS1": {
        "mcs_seq": _PETDUET1_MCS1,       # identical to pETDuet-1 MCS1
        "re_sites": _PETDUET1_MCS1_RE_SITES,
        "tags": {
            "N-His6": {"type": "N-terminal", "region": (12, 27)},
        },
        "stop_pos": 87,
        "aliases": ["pACYCDuet1", "pACYCDuet-1", "pacycduet1",
                     "pACYCDuet1:MCS1"],
        "notes": "MCS1 identical to pETDuet-1:MCS1. "
                 "BamHI/EcoRI/SacI are cross-codon (frame +2).",
    },
    "pACYCDuet-1:MCS2": {
        "mcs_seq": _DUET_MCS2,            # shared Duet MCS2 design
        "re_sites": _DUET_MCS2_RE_SITES,
        "tags": {
            "S-tag": {"type": "C-terminal", "region": (66, 108)},
        },
        "stop_pos": 126,
        "aliases": ["pACYCDuet1:MCS2"],
        "notes": "MCS2 identical to pETDuet-1:MCS2. "
                 "NdeI at position -2 (overlaps start ATG).",
    },
}


# ── Public API ──────────────────────────────────────────────────────────────

def get_vector(name: str) -> dict:
    """벡터 조회 (fuzzy matching: "pET21a" → "pET-21a(+)").

    Returns
    -------
    dict with keys: mcs_seq, re_sites, tags, stop_pos, aliases
    """
    if name in EXPRESSION_VECTORS:
        return EXPRESSION_VECTORS[name]

    name_lower = name.lower().replace(" ", "").replace("-", "").replace("(", "").replace(")", "").replace("+", "")
    for canonical, data in EXPRESSION_VECTORS.items():
        canon_norm = canonical.lower().replace(" ", "").replace("-", "").replace("(", "").replace(")", "").replace("+", "")
        if name_lower == canon_norm:
            return data
        for alias in data["aliases"]:
            alias_norm = alias.lower().replace(" ", "").replace("-", "").replace("(", "").replace(")", "").replace("+", "")
            if name_lower == alias_norm:
                return data

    available = ", ".join(EXPRESSION_VECTORS.keys())
    raise ValueError(f"Unknown vector: {name!r}. Available: {available}")


def check_reading_frame(
    vector_name: str,
    re_5prime: str,
    re_3prime: str,
    insert_has_atg: bool = True,
    insert_has_stop: bool = False,
    insert_cds_bp: int | None = None,
) -> dict:
    """RE 클로닝 전략의 reading frame 호환성 검증.

    Ligation 후 양쪽 RE site이 완전히 복원된다는 모델 사용:
      5' junction: mcs[0 : site5 + len(RE5)] → INSERT start
      3' junction: INSERT end → mcs[site3 : stop]

    Parameters
    ----------
    vector_name : str
        벡터 이름 (fuzzy matching 지원)
    re_5prime, re_3prime : str
        5'/3' restriction enzyme 이름
    insert_has_atg : bool
        Insert에 자체 start codon 포함 여부
    insert_has_stop : bool
        Insert에 자체 stop codon 포함 여부
    insert_cds_bp : int or None
        Insert CDS 길이 (bp, stop codon 제외). None이면 길이 검증 skip.

    Returns
    -------
    dict — in_frame_5prime, in_frame_3prime, linker aa, topology, warnings 등
    """
    vec = get_vector(vector_name)
    mcs = vec["mcs_seq"]
    re_sites = vec["re_sites"]
    stop_pos = vec["stop_pos"]
    warnings = []

    # RE site 존재 확인
    for label, enzyme in [("5'", re_5prime), ("3'", re_3prime)]:
        if enzyme not in re_sites:
            available = ", ".join(sorted(re_sites.keys()))
            raise ValueError(
                f"{enzyme} not in {vector_name} MCS. Available: {available}"
            )

    # RE site positions
    site5 = re_sites[re_5prime]
    site3 = re_sites[re_3prime]
    re5_len = len(RESTRICTION_ENZYMES[re_5prime]["recognition"])
    re3_len = len(RESTRICTION_ENZYMES[re_3prime]["recognition"])

    # 5' RE가 3' RE보다 upstream에 있어야 함
    if site5 >= site3:
        raise ValueError(
            f"5' RE ({re_5prime}, pos {site5}) must be upstream of "
            f"3' RE ({re_3prime}, pos {site3})"
        )

    # ── 5' Frame Analysis ───────────────────────────────────────────────
    # Insert는 5' RE site 직후에 시작
    insert_start = site5 + re5_len
    frame_at_insert_start = insert_start % 3
    in_frame_5prime = (frame_at_insert_start == 0)

    # 5' linker: vector ATG ~ 5' RE site 끝 (insert 시작 직전까지)
    linker_5prime_nt = mcs[:insert_start]
    linker_5prime_aa = _translate(linker_5prime_nt)

    # ── 3' Frame Analysis ───────────────────────────────────────────────
    # 3' 쪽: insert CDS 끝 → 3' RE site 시작 → stop codon
    # Ligation 후 3' RE site 전체가 복원되므로, mcs[site3:stop] 전체가 linker
    nt_3prime_to_stop = stop_pos - site3
    in_frame_3prime = (nt_3prime_to_stop % 3 == 0)

    # 3' linker: 3' RE site 시작 ~ stop codon 직전
    linker_3prime_nt = mcs[site3:stop_pos]
    linker_3prime_aa = _translate(linker_3prime_nt)

    # ── Insert CDS 길이 검증 ────────────────────────────────────────────
    if insert_cds_bp is not None:
        if insert_cds_bp % 3 != 0:
            warnings.append(
                f"Insert CDS length ({insert_cds_bp} bp) is not a multiple of 3"
            )
        if not insert_has_stop:
            # C-tag in-frame이려면: insert_cds_bp + nt_3prime_to_stop ≡ 0 (mod 3)
            total_insert_to_stop = insert_cds_bp + nt_3prime_to_stop
            if total_insert_to_stop % 3 != 0:
                warnings.append(
                    f"Insert ({insert_cds_bp} bp) + 3' linker ({nt_3prime_to_stop} bp) "
                    f"= {total_insert_to_stop} bp, not a multiple of 3 "
                    f"-> C-tag out of frame"
                )

    # ── Context-dependent warnings ──────────────────────────────────────
    if insert_has_stop:
        warnings.append(
            "Insert has stop codon -> C-terminal tag will NOT be fused"
        )

    if insert_has_atg and in_frame_5prime:
        warnings.append(
            "Dual ATG: vector reading frame continues into insert which has "
            "its own ATG. N-terminal fusion tag will be translated"
        )
    elif insert_has_atg and not in_frame_5prime:
        warnings.append(
            "5' junction out of frame -> vector N-terminal tag not fused. "
            "Insert's own ATG will serve as translation start"
        )

    if not in_frame_3prime and not insert_has_stop:
        warnings.append(
            "3' junction out of frame -> C-terminal tag will be mistranslated"
        )

    # ── Topology 문자열 ─────────────────────────────────────────────────
    parts = []

    # N-terminal tags (insert_start보다 앞에 있는 것들)
    for tag_name, tag_info in vec["tags"].items():
        if tag_info["type"] == "N-terminal":
            tag_start = tag_info["region"][0]
            if tag_start < insert_start:
                parts.append(f"[{tag_name}]")

    parts.append(f"({re_5prime})")

    # Insert
    insert_label = "ATG-INSERT" if insert_has_atg else "INSERT"
    if insert_has_stop:
        insert_label += "-STOP"
    parts.append(f"[{insert_label}]")

    parts.append(f"({re_3prime})")

    # C-terminal tags
    if not insert_has_stop:
        for tag_name, tag_info in vec["tags"].items():
            if tag_info["type"] == "C-terminal":
                tag_start = tag_info["region"][0]
                if tag_start >= site3:
                    suffix = "" if in_frame_3prime else ":OUT-OF-FRAME"
                    parts.append(f"[{tag_name}{suffix}]")
        parts.append("[STOP]")

    topology = " ".join(parts)

    return {
        "vector_name": vector_name,
        "re_5prime": re_5prime,
        "re_3prime": re_3prime,
        "in_frame_5prime": in_frame_5prime,
        "in_frame_3prime": in_frame_3prime,
        "frame_at_insert_start": frame_at_insert_start,
        "insert_start_pos": insert_start,
        "linker_5prime_nt": linker_5prime_nt,
        "linker_5prime_aa": linker_5prime_aa,
        "linker_3prime_nt": linker_3prime_nt,
        "linker_3prime_aa": linker_3prime_aa,
        "nt_3prime_to_stop": nt_3prime_to_stop,
        "insert_has_stop": insert_has_stop,
        "topology": topology,
        "warnings": warnings,
    }


def format_frame_report(result: dict) -> str:
    """check_reading_frame 결과를 보기 좋게 포맷팅."""
    lines = [
        f"Vector:   {result['vector_name']}",
        f"Strategy: {result['re_5prime']} / {result['re_3prime']}",
        "",
    ]

    # 5' side
    frame5 = result["frame_at_insert_start"]
    if result["in_frame_5prime"]:
        status5 = "IN-FRAME (insert starts at frame 0)"
    else:
        status5 = f"OUT-OF-FRAME (insert starts at +{frame5} nt offset)"
    lines.append(f"5' junction: {status5}")
    if result["linker_5prime_aa"]:
        lines.append(f"  5' leader:  {result['linker_5prime_aa']} ({len(result['linker_5prime_nt'])} nt)")

    # 3' side
    nt3 = result["nt_3prime_to_stop"]
    if result["in_frame_3prime"]:
        status3 = f"IN-FRAME ({nt3} nt to stop, {nt3 // 3} codons)"
    else:
        status3 = f"OUT-OF-FRAME ({nt3} nt to stop, {nt3 % 3} nt remainder)"
    lines.append(f"3' junction: {status3}")
    if result["linker_3prime_aa"] and not result["insert_has_stop"]:
        lines.append(f"  3' linker:  {result['linker_3prime_aa']}")

    lines.append("")
    lines.append(f"Topology: {result['topology']}")

    if result["warnings"]:
        lines.append("")
        for w in result["warnings"]:
            lines.append(f"  ! {w}")

    return "\n".join(lines)


# ── Tests ───────────────────────────────────────────────────────────────────

def _run_tests():
    """vector_registry 테스트."""
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

    # ── Test 1: pET-21a(+) EcoRI/NotI ───────────────────────────────────
    # 채연 리뷰 케이스: C-His6 fusion 확인
    # EcoRI at 42 (6bp) → insert starts at 48 → 48%3=0 ✓
    # NotI at 66, stop at 99 → 33 nt → 33%3=0 ✓
    # 3' linker: GCG GCC GCA CTC GAG CAC×6 = AAALEHHHHHH
    print(sep)
    print("Test 1: pET-21a(+) EcoRI/NotI -- C-His6 fusion")
    print(sep)
    r1 = check_reading_frame(
        "pET-21a(+)", "EcoRI", "NotI",
        insert_has_atg=True, insert_has_stop=False,
    )
    print(format_frame_report(r1))
    print()
    check("5' in-frame", r1["in_frame_5prime"], True)
    check("3' in-frame", r1["in_frame_3prime"], True)
    check("3' linker aa contains AAALEHHHHHH", r1["linker_3prime_aa"], "AAALEHHHHHH")
    check("insert starts at pos 48", r1["insert_start_pos"], 48)
    print()

    # ── Test 2: pET-28a(+) BamHI/XhoI ──────────────────────────────────
    # N-His6 + C-His6 모두 in-frame 확인
    # BamHI at 96 (6bp) → insert starts at 102 → 102%3=0 ✓
    # XhoI at 135, stop at 159 → 24 nt → 24%3=0 ✓
    print(sep)
    print("Test 2: pET-28a(+) BamHI/XhoI -- dual His6")
    print(sep)
    r2 = check_reading_frame(
        "pET-28a(+)", "BamHI", "XhoI",
        insert_has_atg=False, insert_has_stop=False,
    )
    print(format_frame_report(r2))
    print()
    check("5' in-frame", r2["in_frame_5prime"], True)
    check("3' in-frame", r2["in_frame_3prime"], True)
    print()

    # ── Test 3: pET-21a(+) EcoRI/NotI + stop ───────────────────────────
    # Insert에 stop codon → C-His6 안 붙음
    print(sep)
    print("Test 3: pET-21a(+) EcoRI/NotI + stop -- no C-tag")
    print(sep)
    r3 = check_reading_frame(
        "pET-21a(+)", "EcoRI", "NotI",
        insert_has_atg=True, insert_has_stop=True,
    )
    print(format_frame_report(r3))
    print()
    check("stop warning present",
          any("C-terminal tag will NOT" in w for w in r3["warnings"]), True)
    print()

    # ── Test 4: Frame mismatch ──────────────────────────────────────────
    # insert_cds_bp=901 → 3의 배수 아님
    print(sep)
    print("Test 4: insert_cds_bp=901 -- frame mismatch")
    print(sep)
    r4 = check_reading_frame(
        "pET-21a(+)", "EcoRI", "NotI",
        insert_has_atg=True, insert_has_stop=False,
        insert_cds_bp=901,
    )
    print(format_frame_report(r4))
    print()
    check("non-3n warning",
          any("not a multiple of 3" in w for w in r4["warnings"]), True)
    print()

    # ── Test 5: Fuzzy name matching ─────────────────────────────────────
    print(sep)
    print("Test 5: Fuzzy name matching")
    print(sep)
    for name in ["pET21a", "pET-21a", "pet21a", "pET-28a", "pet28a",
                  "pMALc6T", "pmal-c6t", "pETDuet1", "petduet-1",
                  "pACYCDuet1", "pacycduet1", "pACYCDuet1:MCS2"]:
        try:
            v = get_vector(name)
            print(f"  OK: '{name}' -> MCS {len(v['mcs_seq'])} bp")
            passed += 1
        except ValueError as e:
            print(f"  FAIL: '{name}' -> {e}")
            failed += 1
    print()

    # ── Test 6: pET-28a(+) NdeI/XhoI ───────────────────────────────────
    # NdeI → thrombin 이후 insert, XhoI → C-His6
    # NdeI at 57 (6bp) → insert starts at 63 → 63%3=0 ✓
    # XhoI at 135, stop at 159 → 24%3=0 ✓
    print(sep)
    print("Test 6: pET-28a(+) NdeI/XhoI -- N-His6 + thrombin cleavable")
    print(sep)
    r6 = check_reading_frame(
        "pET-28a(+)", "NdeI", "XhoI",
        insert_has_atg=True, insert_has_stop=False,
    )
    print(format_frame_report(r6))
    print()
    check("5' in-frame (NdeI)", r6["in_frame_5prime"], True)
    check("3' in-frame (XhoI)", r6["in_frame_3prime"], True)
    print()

    # ── Test 7: pMAL-c6T NotI/BamHI ─────────────────────────────────────
    # NotI at 13 (8bp) -> insert starts at 21 -> 21%3=0 ✓
    # BamHI at 33, stop at 63 -> 30 nt -> 30%3=0 ✓
    # 3' linker: mcs[33:63] = GGATCC GAATTC CCTGCA GGTAAT AAGCTT TAA
    #            -> GSEFPAGNKL*
    print(sep)
    print("Test 7: pMAL-c6T NotI/BamHI -- MBP fusion")
    print(sep)
    r7 = check_reading_frame(
        "pMAL-c6T", "NotI", "BamHI",
        insert_has_atg=False, insert_has_stop=False,
    )
    print(format_frame_report(r7))
    print()
    check("5' in-frame", r7["in_frame_5prime"], True)
    check("3' in-frame", r7["in_frame_3prime"], True)
    check("insert starts at pos 21", r7["insert_start_pos"], 21)
    print()

    # ── Test 8: pMAL-c6T EcoRI/HindIII ───────────────────────────────
    # EcoRI at 39 (6bp) -> insert starts at 45 -> 45%3=0 ✓
    # HindIII at 57, stop at 63 -> 6 nt -> 6%3=0 ✓
    print(sep)
    print("Test 8: pMAL-c6T EcoRI/HindIII -- compact linker")
    print(sep)
    r8 = check_reading_frame(
        "pMAL-c6T", "EcoRI", "HindIII",
        insert_has_atg=False, insert_has_stop=False,
    )
    print(format_frame_report(r8))
    print()
    check("5' in-frame", r8["in_frame_5prime"], True)
    check("3' in-frame", r8["in_frame_3prime"], True)
    check("3' linker = KL", r8["linker_3prime_aa"], "KL")
    print()

    # ── Test 9: pETDuet-1:MCS1 PstI/NotI ─────────────────────────────
    # Codon-aligned pair: PstI at 60 (6bp) -> insert starts at 66 -> 66%3=0 ✓
    # NotI at 78, stop at 87 -> 9 nt -> 9%3=0 ✓
    print(sep)
    print("Test 9: pETDuet-1:MCS1 PstI/NotI -- codon-aligned")
    print(sep)
    r9 = check_reading_frame(
        "pETDuet-1:MCS1", "PstI", "NotI",
        insert_has_atg=False, insert_has_stop=False,
    )
    print(format_frame_report(r9))
    print()
    check("5' in-frame", r9["in_frame_5prime"], True)
    check("3' in-frame", r9["in_frame_3prime"], True)
    check("3' linker = AAA", r9["linker_3prime_aa"], "AAA")
    print()

    # ── Test 10: pETDuet-1:MCS1 BamHI/NotI ────────────────────────────
    # BamHI at 35 (cross-codon, frame +2): 35+6=41 -> 41%3=2 -> OUT OF FRAME
    # NotI at 78, stop at 87 -> 9 nt -> 9%3=0 ✓
    print(sep)
    print("Test 10: pETDuet-1:MCS1 BamHI/NotI -- 5' out-of-frame")
    print(sep)
    r10 = check_reading_frame(
        "pETDuet-1:MCS1", "BamHI", "NotI",
        insert_has_atg=True, insert_has_stop=False,
    )
    print(format_frame_report(r10))
    print()
    check("5' OUT-of-frame", r10["in_frame_5prime"], False)
    check("frame offset = 2", r10["frame_at_insert_start"], 2)
    check("3' in-frame", r10["in_frame_3prime"], True)
    print()

    # ── Test 11: Duet MCS2 KpnI/XhoI ────────────────────────────────────
    # Codon-aligned pair: KpnI at 48 (6bp) -> insert starts at 54 -> 54%3=0
    # XhoI at 54, stop at 126 -> 72 nt -> 72%3=0
    # S-tag should be in-frame as C-terminal fusion
    print(sep)
    print("Test 11: pETDuet-1:MCS2 KpnI/XhoI -- S-tag fusion")
    print(sep)
    r11 = check_reading_frame(
        "pETDuet-1:MCS2", "KpnI", "XhoI",
        insert_has_atg=False, insert_has_stop=False,
    )
    print(format_frame_report(r11))
    print()
    check("5' in-frame", r11["in_frame_5prime"], True)
    check("3' in-frame", r11["in_frame_3prime"], True)
    check("S-tag in topology",
          "S-tag" in r11["topology"], True)
    print()

    # ── Test 12: Duet MCS2 BglII/XhoI ─────────────────────────────────
    # BglII at 5 (cross-codon): 5+6=11 -> 11%3=2 -> OUT OF FRAME
    # XhoI at 54, stop at 126 -> 72%3=0 -> in-frame
    print(sep)
    print("Test 12: pACYCDuet-1:MCS2 BglII/XhoI -- 5' out-of-frame")
    print(sep)
    r12 = check_reading_frame(
        "pACYCDuet-1:MCS2", "BglII", "XhoI",
        insert_has_atg=True, insert_has_stop=False,
    )
    print(format_frame_report(r12))
    print()
    check("5' OUT-of-frame", r12["in_frame_5prime"], False)
    check("frame offset = 2", r12["frame_at_insert_start"], 2)
    check("3' in-frame", r12["in_frame_3prime"], True)
    print()

    # ── Summary ─────────────────────────────────────────────────────────
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
