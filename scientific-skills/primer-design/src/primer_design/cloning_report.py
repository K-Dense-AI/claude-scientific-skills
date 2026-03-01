"""RE Cloning primer design report & vector construct map (matplotlib).

design_re_cloning_primers() 결과를 시각적 리포트(PNG)로 생성.
- 프라이머 구조 다이어그램 (protection / RE / spacer / annealing 색상 구분)
- QC 요약 테이블
- 클로닝 construct 맵
- Circular vector construct map (원형 플라스미드 맵)
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Arc, Wedge
import numpy as np

# ── Color scheme ─────────────────────────────────────────────────────────────

_C = {
    "protection": "#B0BEC5",
    "re_site":    "#EF5350",
    "spacer":     "#FFE082",
    "annealing":  "#42A5F5",
    "stop_rc":    "#FF7043",
    "cds":        "#66BB6A",
    "pass":       "#4CAF50",
    "warning":    "#FF9800",
    "fail":       "#F44336",
    "header_bg":  "#1565C0",
    "header_tx":  "#FFFFFF",
    "bg":         "#FAFAFA",
    "border":     "#E0E0E0",
}

_FONT = "DejaVu Sans Mono"
_FONT_SANS = "DejaVu Sans"


# ── Public API ───────────────────────────────────────────────────────────────

def generate_cloning_report(
    design_result: dict,
    output_path: Path | str,
    gene_name: str = "Insert",
) -> Path:
    """Generate a visual primer design report as PNG.

    Parameters
    ----------
    design_result : dict
        Output from ``RestrictionCloningDesigner.design()``.
    output_path : Path or str
        Output PNG file path.
    gene_name : str
        Gene name for labeling.

    Returns
    -------
    Path : Written file path.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(16, 11), facecolor="white")

    gs = fig.add_gridspec(
        4, 1,
        height_ratios=[0.7, 2.8, 1.8, 2.0],
        hspace=0.25,
        left=0.04, right=0.96, top=0.96, bottom=0.03,
    )

    # Row 0: Header
    ax0 = fig.add_subplot(gs[0])
    _draw_header(ax0, design_result, gene_name)

    # Row 1: Primer architecture
    ax1 = fig.add_subplot(gs[1])
    _draw_primer_architecture(ax1, design_result)

    # Row 2: QC table
    ax2 = fig.add_subplot(gs[2])
    _draw_qc_table(ax2, design_result)

    # Row 3: Construct map
    ax3 = fig.add_subplot(gs[3])
    _draw_construct_map(ax3, design_result, gene_name)

    fig.savefig(str(output_path), dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return output_path


# ── Drawing helpers ──────────────────────────────────────────────────────────

def _draw_header(ax, dr: dict, gene_name: str):
    """Draw the title/header bar."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    # Background
    ax.add_patch(FancyBboxPatch(
        (0, 0), 1, 1,
        boxstyle="round,pad=0.02",
        facecolor=_C["header_bg"], edgecolor="none",
    ))

    ax.text(
        0.5, 0.62, "RE Cloning Primer Design Report",
        ha="center", va="center",
        fontsize=16, fontweight="bold", color=_C["header_tx"],
        fontfamily=_FONT_SANS,
    )

    insert_len = dr.get("insert_len", 0)
    vector = dr.get("frame_check", {}).get("vector_name", "N/A")
    if vector == "N/A" and "frame_report" in dr:
        # try to extract from frame_report text
        pass

    info = (
        f"Gene: {gene_name}    "
        f"Insert: {insert_len} bp    "
        f"RE: {dr['re_5prime']} / {dr['re_3prime']}    "
        f"Date: {date.today().isoformat()}"
    )
    ax.text(
        0.5, 0.20, info,
        ha="center", va="center",
        fontsize=10, color=_C["header_tx"],
        fontfamily=_FONT,
    )


def _draw_primer_architecture(ax, dr: dict):
    """Draw colored primer structure diagrams for F and R primers."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_axis_off()

    ax.add_patch(FancyBboxPatch(
        (0, 0), 10, 10,
        boxstyle="round,pad=0.1",
        facecolor=_C["bg"], edgecolor=_C["border"], linewidth=1,
    ))

    # ── Forward primer ────────────────────────────────────────────────────
    _draw_single_primer(
        ax, y_center=7.5,
        label="Forward Primer",
        full_seq=dr["f_full"],
        tm=dr["f_tm"],
        gc=dr["f_gc"],
        total_len=dr["f_len"],
        sections=_build_sections(dr, "forward"),
    )

    # ── Reverse primer ────────────────────────────────────────────────────
    _draw_single_primer(
        ax, y_center=3.0,
        label="Reverse Primer",
        full_seq=dr["r_full"],
        tm=dr["r_tm"],
        gc=dr["r_gc"],
        total_len=dr["r_len"],
        sections=_build_sections(dr, "reverse"),
    )


def _build_sections(dr: dict, direction: str) -> list[dict]:
    """Build a list of primer sections for visualization."""
    sections = []

    if direction == "forward":
        prot_len = len(dr["protection_5"])
        re_site = dr["re_5prime_site"]
        re_name = dr["re_5prime"]
        spacer = dr.get("spacer_5prime", "")
        ann = dr["f_ann"]
        ann_len = dr["f_ann_len"]

        if prot_len > 0:
            sections.append({
                "label": f"prot ({prot_len}bp)",
                "seq": dr["protection_5"],
                "color": _C["protection"],
                "width": prot_len,
            })
        sections.append({
            "label": re_name,
            "seq": re_site,
            "color": _C["re_site"],
            "width": len(re_site),
        })
        if spacer:
            sections.append({
                "label": "spacer",
                "seq": spacer,
                "color": _C["spacer"],
                "width": len(spacer),
            })
        sections.append({
            "label": f"annealing ({ann_len}bp)",
            "seq": ann,
            "color": _C["annealing"],
            "width": ann_len,
        })

    else:  # reverse
        prot_len = len(dr["protection_3"])
        re_site = dr["re_3prime_site"]
        re_name = dr["re_3prime"]
        spacer = dr.get("spacer_3prime", "")
        ann = dr["r_ann"]
        ann_len = dr["r_ann_len"]

        # r_tail = protection + RE + spacer + [stop_rc]
        if prot_len > 0:
            sections.append({
                "label": f"prot ({prot_len}bp)",
                "seq": dr["protection_3"],
                "color": _C["protection"],
                "width": prot_len,
            })
        sections.append({
            "label": re_name,
            "seq": re_site,
            "color": _C["re_site"],
            "width": len(re_site),
        })
        if spacer:
            sections.append({
                "label": "spacer",
                "seq": spacer,
                "color": _C["spacer"],
                "width": len(spacer),
            })
        if dr.get("include_stop_codon"):
            sections.append({
                "label": "stop(RC)",
                "seq": "TTA",
                "color": _C["stop_rc"],
                "width": 3,
            })
        sections.append({
            "label": f"annealing ({ann_len}bp)",
            "seq": ann,
            "color": _C["annealing"],
            "width": ann_len,
        })

    return sections


def _draw_single_primer(ax, y_center, label, full_seq, tm, gc, total_len, sections):
    """Draw a single primer with colored section boxes."""
    # Title line
    ax.text(
        0.4, y_center + 1.6,
        f"{label}  ({total_len} nt,  Tm = {tm}°C,  GC = {gc}%)",
        fontsize=11, fontweight="bold", fontfamily=_FONT_SANS,
        va="center",
    )

    # Calculate proportional widths
    total_bp = sum(s["width"] for s in sections)
    bar_x_start = 0.6
    bar_width_total = 8.8
    bar_height = 0.8

    x = bar_x_start

    # 5' label
    ax.text(
        x - 0.35, y_center, "5'",
        fontsize=10, fontweight="bold", ha="right", va="center",
        fontfamily=_FONT,
    )

    for sec in sections:
        w = (sec["width"] / total_bp) * bar_width_total
        rect = FancyBboxPatch(
            (x, y_center - bar_height / 2), w, bar_height,
            boxstyle="round,pad=0.02",
            facecolor=sec["color"],
            edgecolor="#555555",
            linewidth=0.8,
        )
        ax.add_patch(rect)

        # Section label (inside box)
        ax.text(
            x + w / 2, y_center + 0.05,
            sec["label"],
            ha="center", va="center",
            fontsize=8, fontweight="bold",
            fontfamily=_FONT_SANS,
            color="#333333",
        )

        # Sequence snippet (below box, truncated)
        seq_display = sec["seq"]
        if len(seq_display) > 14:
            seq_display = seq_display[:6] + "..." + seq_display[-5:]
        ax.text(
            x + w / 2, y_center - 0.65,
            seq_display,
            ha="center", va="center",
            fontsize=6.5,
            fontfamily=_FONT,
            color="#666666",
        )

        x += w

    # 3' label
    ax.text(
        x + 0.15, y_center, "3'",
        fontsize=10, fontweight="bold", ha="left", va="center",
        fontfamily=_FONT,
    )

    # Full sequence line (below sections)
    full_display = full_seq
    if len(full_display) > 70:
        full_display = full_display[:35] + "..." + full_display[-32:]
    ax.text(
        0.4, y_center - 1.2,
        f"5'-{full_display}-3'",
        fontsize=7, fontfamily=_FONT, color="#888888",
        va="center",
    )


def _draw_qc_table(ax, dr: dict):
    """Draw the QC summary table."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_axis_off()

    ax.add_patch(FancyBboxPatch(
        (0, 0), 10, 10,
        boxstyle="round,pad=0.1",
        facecolor=_C["bg"], edgecolor=_C["border"], linewidth=1,
    ))

    ax.text(
        0.4, 9.2, "Quality Control",
        fontsize=12, fontweight="bold", fontfamily=_FONT_SANS,
    )

    # Table data
    f_qc = dr.get("f_qc", {})
    r_qc = dr.get("r_qc", {})
    het = dr.get("het", {})

    col_labels = ["Primer", "Tm (°C)", "GC%", "Length", "Anneal", "Hairpin Tm", "Verdict"]
    cell_data = [
        [
            "Forward",
            f"{dr['f_tm']}",
            f"{dr['f_gc']}",
            f"{dr['f_len']} nt",
            f"{dr['f_ann_len']} bp",
            f"{f_qc.get('hairpin_tm', 'N/A')}°C",
            f_qc.get("verdict", "N/A"),
        ],
        [
            "Reverse",
            f"{dr['r_tm']}",
            f"{dr['r_gc']}",
            f"{dr['r_len']} nt",
            f"{dr['r_ann_len']} bp",
            f"{r_qc.get('hairpin_tm', 'N/A')}°C",
            r_qc.get("verdict", "N/A"),
        ],
    ]

    # Draw table manually for reliable rendering
    col_widths = [1.1, 0.8, 0.7, 0.8, 0.8, 1.0, 0.9]
    total_w = sum(col_widths)
    scale = 8.5 / total_w
    col_widths = [w * scale for w in col_widths]

    x_start = 0.5
    y_top = 8.2
    row_h = 1.2
    header_h = 1.2

    # Header row
    x = x_start
    for j, label in enumerate(col_labels):
        w = col_widths[j]
        ax.add_patch(FancyBboxPatch(
            (x, y_top - header_h), w, header_h,
            boxstyle="square,pad=0",
            facecolor="#37474F", edgecolor="#263238", linewidth=0.8,
        ))
        ax.text(
            x + w / 2, y_top - header_h / 2, label,
            ha="center", va="center",
            fontsize=9, fontweight="bold", color="white",
            fontfamily=_FONT_SANS,
        )
        x += w

    # Data rows
    for i, row in enumerate(cell_data):
        y_row = y_top - header_h - i * row_h
        row_bg = "#FFFFFF" if i % 2 == 0 else "#F5F5F5"
        x = x_start
        for j, val in enumerate(row):
            w = col_widths[j]
            # Verdict cell coloring
            if j == len(row) - 1:
                if val == "PASS":
                    cell_bg = "#C8E6C9"
                    cell_color = "#2E7D32"
                elif val == "WARNING":
                    cell_bg = "#FFF3E0"
                    cell_color = "#E65100"
                elif val == "FAIL":
                    cell_bg = "#FFCDD2"
                    cell_color = "#C62828"
                else:
                    cell_bg = row_bg
                    cell_color = "#333333"
            else:
                cell_bg = row_bg
                cell_color = "#333333"

            ax.add_patch(FancyBboxPatch(
                (x, y_row - row_h), w, row_h,
                boxstyle="square,pad=0",
                facecolor=cell_bg, edgecolor="#CCCCCC", linewidth=0.5,
            ))
            fw = "bold" if j == 0 or j == len(row) - 1 else "normal"
            ax.text(
                x + w / 2, y_row - row_h / 2, val,
                ha="center", va="center",
                fontsize=9, fontweight=fw, color=cell_color,
                fontfamily=_FONT_SANS,
            )
            x += w

    # Additional info below table
    het_dg = het.get("dg", "N/A")
    het_tm = het.get("tm", "N/A")
    anneal = dr.get("anneal_temp", "N/A")

    info_y = y_top - header_h - len(cell_data) * row_h - 0.8
    info_text = (
        f"Annealing temp: {anneal}°C    |    "
        f"Heterodimer: dG = {het_dg} kcal/mol, Tm = {het_tm}°C"
    )
    ax.text(
        0.5, info_y, info_text,
        fontsize=9.5, fontfamily=_FONT, color="#555555",
    )

    # Warnings
    warnings = dr.get("warnings", [])
    if warnings:
        warn_text = "Warnings: " + "; ".join(warnings[:3])
        if len(warnings) > 3:
            warn_text += f" (+{len(warnings) - 3} more)"
        ax.text(
            0.5, info_y - 1.0,
            warn_text,
            fontsize=8, fontfamily=_FONT_SANS, color=_C["warning"],
            style="italic",
        )


def _draw_construct_map(ax, dr: dict, gene_name: str):
    """Draw linear construct map and reading frame info."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_axis_off()

    ax.add_patch(FancyBboxPatch(
        (0, 0), 10, 10,
        boxstyle="round,pad=0.1",
        facecolor=_C["bg"], edgecolor=_C["border"], linewidth=1,
    ))

    ax.text(
        0.4, 9.2, "Cloning Construct (PCR Product)",
        fontsize=12, fontweight="bold", fontfamily=_FONT_SANS,
    )

    # Construct diagram
    y_bar = 6.5
    bar_height = 1.0
    x_start = 1.0
    x_end = 9.0
    bar_w = x_end - x_start

    insert_len = dr.get("insert_len", 0)
    f_tail_len = len(dr["f_tail"])
    r_tail_len = len(dr["r_tail"])
    total_bp = f_tail_len + insert_len + r_tail_len

    # Proportional widths
    w_ftail = (f_tail_len / total_bp) * bar_w if total_bp else 0.5
    w_insert = (insert_len / total_bp) * bar_w if total_bp else 6.0
    w_rtail = (r_tail_len / total_bp) * bar_w if total_bp else 0.5

    # Minimum widths for visibility
    min_w = 0.6
    if w_ftail < min_w:
        w_ftail = min_w
    if w_rtail < min_w:
        w_rtail = min_w
    w_insert = bar_w - w_ftail - w_rtail

    # 5' RE site box
    ax.add_patch(FancyBboxPatch(
        (x_start, y_bar - bar_height / 2), w_ftail, bar_height,
        boxstyle="round,pad=0.02",
        facecolor=_C["re_site"], edgecolor="#333", linewidth=1,
    ))
    ax.text(
        x_start + w_ftail / 2, y_bar,
        dr["re_5prime"],
        ha="center", va="center",
        fontsize=9, fontweight="bold", color="white",
        fontfamily=_FONT_SANS,
    )

    # Insert CDS box
    cds_x = x_start + w_ftail
    ax.add_patch(FancyBboxPatch(
        (cds_x, y_bar - bar_height / 2), w_insert, bar_height,
        boxstyle="round,pad=0.02",
        facecolor=_C["cds"], edgecolor="#333", linewidth=1,
    ))
    ax.text(
        cds_x + w_insert / 2, y_bar,
        f"{gene_name} CDS ({insert_len} bp)",
        ha="center", va="center",
        fontsize=10, fontweight="bold", color="white",
        fontfamily=_FONT_SANS,
    )

    # 3' RE site box
    re3_x = cds_x + w_insert
    ax.add_patch(FancyBboxPatch(
        (re3_x, y_bar - bar_height / 2), w_rtail, bar_height,
        boxstyle="round,pad=0.02",
        facecolor=_C["re_site"], edgecolor="#333", linewidth=1,
    ))
    ax.text(
        re3_x + w_rtail / 2, y_bar,
        dr["re_3prime"],
        ha="center", va="center",
        fontsize=9, fontweight="bold", color="white",
        fontfamily=_FONT_SANS,
    )

    # 5'→3' direction arrow
    ax.annotate(
        "", xy=(x_end + 0.15, y_bar), xytext=(x_start - 0.15, y_bar),
        arrowprops=dict(
            arrowstyle="->", color="#888", lw=1.5,
            connectionstyle="arc3,rad=0",
        ),
    )
    ax.text(x_start - 0.35, y_bar, "5'", fontsize=9, ha="right",
            va="center", fontweight="bold", fontfamily=_FONT)
    ax.text(x_end + 0.35, y_bar, "3'", fontsize=9, ha="left",
            va="center", fontweight="bold", fontfamily=_FONT)

    # Total PCR product size
    total_pcr = f_tail_len + insert_len + r_tail_len
    ax.text(
        5.0, y_bar - 1.0,
        f"PCR product: {total_pcr} bp",
        ha="center", va="center",
        fontsize=9, fontfamily=_FONT, color="#555",
    )

    # ── Reading frame info ────────────────────────────────────────────────
    fc = dr.get("frame_check")
    if fc:
        y_frame = 3.0
        in5 = fc.get("in_frame_5prime", False)
        in3 = fc.get("in_frame_3prime", False)

        mark5 = "\u2713" if in5 else "\u2717"
        mark3 = "\u2713" if in3 else "\u2717"
        c5 = _C["pass"] if in5 else _C["fail"]
        c3 = _C["pass"] if in3 else _C["fail"]

        ax.text(
            0.4, y_frame + 0.6,
            "Reading Frame:",
            fontsize=10, fontweight="bold", fontfamily=_FONT_SANS,
        )
        ax.text(
            3.0, y_frame + 0.6,
            f"{mark5} 5' junction",
            fontsize=10, color=c5, fontfamily=_FONT_SANS,
        )
        ax.text(
            5.5, y_frame + 0.6,
            f"{mark3} 3' junction",
            fontsize=10, color=c3, fontfamily=_FONT_SANS,
        )

        # Topology
        topology = fc.get("topology", "")
        if topology:
            # Show only first line of topology
            topo_line = topology.split("\n")[0] if "\n" in topology else topology
            if len(topo_line) > 90:
                topo_line = topo_line[:87] + "..."
            ax.text(
                0.4, y_frame - 0.3,
                topo_line,
                fontsize=7.5, fontfamily=_FONT, color="#666",
            )

        # Linker sequences
        linker5 = fc.get("linker_5prime_aa", "")
        linker3 = fc.get("linker_3prime_aa", "")
        if linker5 or linker3:
            linker_text = f"5' linker: {linker5 or 'N/A'}    3' linker: {linker3 or 'N/A'}"
            ax.text(
                0.4, y_frame - 1.0,
                linker_text,
                fontsize=8, fontfamily=_FONT, color="#777",
            )
    else:
        ax.text(
            0.4, 3.0,
            "Reading frame: (no vector specified)",
            fontsize=10, fontfamily=_FONT_SANS, color="#999",
            style="italic",
        )

    # ── Legend ─────────────────────────────────────────────────────────────
    legend_y = 0.7
    legend_items = [
        (_C["protection"], "Protection"),
        (_C["re_site"], "RE site"),
        (_C["spacer"], "Spacer"),
        (_C["annealing"], "Annealing"),
        (_C["cds"], "CDS"),
    ]
    if dr.get("include_stop_codon"):
        legend_items.insert(4, (_C["stop_rc"], "Stop(RC)"))

    x_legend = 0.5
    for color, text in legend_items:
        ax.add_patch(FancyBboxPatch(
            (x_legend, legend_y - 0.15), 0.3, 0.3,
            boxstyle="round,pad=0.02",
            facecolor=color, edgecolor="#999", linewidth=0.5,
        ))
        ax.text(
            x_legend + 0.4, legend_y,
            text,
            fontsize=7.5, fontfamily=_FONT_SANS, va="center", color="#555",
        )
        x_legend += 1.5


# ── Backbone feature data for supported vectors ─────────────────────────────
# Angular fractions (0 = top/12 o'clock, clockwise) + approximate sizes
# Based on GenBank annotations for each vector series.

_BACKBONE_FEATURES: dict[str, dict] = {
    "pET-28a(+)": {
        "total_bp": 5369,
        "resistance": "KanR",
        "features": [
            # (name, frac_start, frac_span, color, label_side)
            ("T7 promoter",    0.065, 0.010, "#E74C3C", "out"),
            ("lac operator",   0.076, 0.008, "#E67E22", "out"),
            ("RBS",            0.085, 0.004, "#F39C12", "out"),
            ("T7 terminator",  0.005, 0.015, "#C0392B", "out"),
            ("lacI",           0.14,  0.20,  "#95A5A6", "out"),
            ("KanR",           0.74,  0.15,  "#E74C3C", "out"),
            ("f1 ori",         0.91,  0.08,  "#3498DB", "out"),
            ("pBR322 ori",     0.61,  0.12,  "#2980B9", "out"),
        ],
    },
    "pET-21a(+)": {
        "total_bp": 5443,
        "resistance": "AmpR",
        "features": [
            ("T7 promoter",    0.065, 0.010, "#E74C3C", "out"),
            ("lac operator",   0.076, 0.008, "#E67E22", "out"),
            ("RBS",            0.085, 0.004, "#F39C12", "out"),
            ("T7 terminator",  0.005, 0.015, "#C0392B", "out"),
            ("lacI",           0.14,  0.20,  "#95A5A6", "out"),
            ("AmpR",           0.72,  0.16,  "#E74C3C", "out"),
            ("f1 ori",         0.91,  0.08,  "#3498DB", "out"),
            ("pBR322 ori",     0.61,  0.10,  "#2980B9", "out"),
        ],
    },
    "pMAL-c6T": {
        "total_bp": 6721,
        "resistance": "AmpR",
        "features": [
            ("Ptac promoter",  0.065, 0.010, "#E74C3C", "out"),
            ("malE (MBP)",     0.08,  0.17,  "#8E44AD", "out"),
            ("TEV site",       0.25,  0.005, "#F39C12", "out"),
            ("rrnB terminator",0.005, 0.015, "#C0392B", "out"),
            ("AmpR",           0.72,  0.14,  "#E74C3C", "out"),
            ("pBR322 ori",     0.58,  0.10,  "#2980B9", "out"),
            ("lacIq",          0.38,  0.15,  "#95A5A6", "out"),
        ],
    },
}

# Duet vectors reuse pET backbone layout
for _duet in ("pETDuet-1:MCS1", "pETDuet-1:MCS2",
              "pACYCDuet-1:MCS1", "pACYCDuet-1:MCS2"):
    _res = "CmR" if "pACYC" in _duet else "AmpR"
    _BACKBONE_FEATURES[_duet] = {
        "total_bp": 5420 if "pET" in _duet else 4008,
        "resistance": _res,
        "features": [
            ("T7 promoter",  0.065, 0.010, "#E74C3C", "out"),
            ("T7 terminator",0.005, 0.015, "#C0392B", "out"),
            ("lacI",         0.14,  0.20,  "#95A5A6", "out"),
            (_res,           0.72,  0.15,  "#E74C3C", "out"),
            ("ori",          0.61,  0.10,  "#2980B9", "out"),
        ],
    }


# ── Color palette for construct map ──────────────────────────────────────────

_CM = {
    "backbone":    "#B0BEC5",
    "insert":      "#2ECC71",
    "tag":         "#9B59B6",
    "re_site":     "#E74C3C",
    "promoter":    "#E74C3C",
    "resistance":  "#E74C3C",
    "origin":      "#3498DB",
    "regulatory":  "#95A5A6",
    "signal_pep":  "#F39C12",
}


# ── Public API: Circular Vector Construct Map ────────────────────────────────

def generate_vector_construct_map(
    vector_name: str,
    gene_name: str = "Insert",
    insert_len: int = 0,
    re_5prime: str = "",
    re_3prime: str = "",
    include_stop: bool = False,
    frame_check: dict | None = None,
    signal_peptide: dict | None = None,
    output_path: Path | str | None = None,
) -> Path:
    """Circular vector construct map (원형 플라스미드 맵) PNG 생성.

    Parameters
    ----------
    vector_name : str
        벡터 이름 (예: "pET-28a(+)")
    gene_name : str
        삽입 유전자 이름
    insert_len : int
        Insert 길이 (bp)
    re_5prime, re_3prime : str
        5'/3' restriction enzyme 이름
    include_stop : bool
        Insert에 stop codon 포함 여부
    frame_check : dict | None
        check_reading_frame() 결과
    signal_peptide : dict | None
        signal peptide 분석 결과 (has_signal_peptide, cleavage_site_estimate 등)
    output_path : Path | str | None
        출력 파일 경로 (.png)

    Returns
    -------
    Path : 생성된 PNG 파일 경로
    """
    if output_path is None:
        output_path = Path.cwd() / f"{gene_name}_{vector_name.replace('(', '').replace(')', '').replace('+', '')}_construct.png"
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10), facecolor="white")
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_aspect("equal")
    ax.set_axis_off()

    # ── Get vector data ─────────────────────────────────────────────────
    # Try to resolve canonical name for backbone lookup
    backbone = None
    for canon_name, bb_data in _BACKBONE_FEATURES.items():
        if canon_name.lower().replace(" ", "").replace("-", "").replace("(", "").replace(")", "").replace("+", "") == \
           vector_name.lower().replace(" ", "").replace("-", "").replace("(", "").replace(")", "").replace("+", ""):
            backbone = bb_data
            break
    if backbone is None:
        # Fallback: generic E. coli expression vector
        backbone = {
            "total_bp": 5000,
            "resistance": "AmpR",
            "features": [
                ("Promoter",    0.065, 0.010, "#E74C3C", "out"),
                ("Terminator",  0.005, 0.015, "#C0392B", "out"),
                ("AmpR",        0.72,  0.15,  "#E74C3C", "out"),
                ("ori",         0.61,  0.10,  "#2980B9", "out"),
            ],
        }

    total_bp = backbone["total_bp"] + insert_len
    R = 1.0  # backbone circle radius
    lw_backbone = 10
    lw_feature = 14
    lw_insert = 18

    # ── Draw backbone circle ─────────────────────────────────────────────
    circle = plt.Circle((0, 0), R, fill=False, edgecolor=_CM["backbone"],
                         linewidth=lw_backbone, zorder=1)
    ax.add_patch(circle)

    # ── Helper: fraction → angle (0=top, clockwise) → matplotlib angle ──
    def frac_to_angle(frac):
        """Fraction (0=top, clockwise) to matplotlib angle (degrees, CCW from right)."""
        return 90 - frac * 360

    def draw_arc(frac_start, frac_span, radius, color, linewidth, zorder=2):
        """Draw a colored arc on the circle."""
        angle_start = frac_to_angle(frac_start + frac_span)
        angle_end = frac_to_angle(frac_start)
        theta = np.linspace(np.radians(angle_start), np.radians(angle_end), 100)
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        ax.plot(x, y, color=color, linewidth=linewidth, solid_capstyle="round", zorder=zorder)

    def label_at_frac(frac, text, radius_offset=0.18, fontsize=8, color="#333",
                      fontweight="normal", ha_override=None):
        """Place a label at a fraction position, pointing outward."""
        angle_rad = np.radians(frac_to_angle(frac))
        r_label = R + radius_offset
        x = r_label * np.cos(angle_rad)
        y = r_label * np.sin(angle_rad)

        # Determine horizontal alignment based on position
        if ha_override:
            ha = ha_override
        elif abs(x) < 0.05:
            ha = "center"
        elif x > 0:
            ha = "left"
        else:
            ha = "right"

        ax.text(x, y, text, fontsize=fontsize, fontweight=fontweight,
                fontfamily=_FONT_SANS, color=color, ha=ha, va="center",
                zorder=10)

    def tick_at_frac(frac, radius, length=0.04, color="#555", lw=1.5):
        """Draw a small radial tick mark."""
        angle_rad = np.radians(frac_to_angle(frac))
        x1 = (radius - length / 2) * np.cos(angle_rad)
        y1 = (radius - length / 2) * np.sin(angle_rad)
        x2 = (radius + length / 2) * np.cos(angle_rad)
        y2 = (radius + length / 2) * np.sin(angle_rad)
        ax.plot([x1, x2], [y1, y2], color=color, linewidth=lw, zorder=5)

    # ── Insert region ────────────────────────────────────────────────────
    # Place insert at the top of the circle (fraction ~0.92 to ~0.08)
    insert_frac = insert_len / total_bp if total_bp > 0 else 0.1
    insert_frac = max(insert_frac, 0.06)  # minimum visibility
    insert_start_frac = 1.0 - insert_frac / 2  # centered at top

    # Draw insert arc (thick, green)
    draw_arc(insert_start_frac, insert_frac, R, _CM["insert"], lw_insert, zorder=3)

    # Insert label
    insert_mid = (insert_start_frac + insert_frac / 2) % 1.0
    label_at_frac(insert_mid, f"{gene_name}\n({insert_len} bp)",
                  radius_offset=0.22, fontsize=11, color="#1B5E20",
                  fontweight="bold")

    # ── RE site markers ──────────────────────────────────────────────────
    if re_5prime:
        tick_at_frac(insert_start_frac, R, length=0.08, color=_CM["re_site"], lw=2.5)
        label_at_frac(insert_start_frac, re_5prime,
                      radius_offset=0.14, fontsize=8, color=_CM["re_site"],
                      fontweight="bold")

    if re_3prime:
        re3_frac = (insert_start_frac + insert_frac) % 1.0
        tick_at_frac(re3_frac, R, length=0.08, color=_CM["re_site"], lw=2.5)
        label_at_frac(re3_frac, re_3prime,
                      radius_offset=0.14, fontsize=8, color=_CM["re_site"],
                      fontweight="bold")

    # ── N-terminal / C-terminal tags (from frame_check) ──────────────────
    if frame_check:
        topology = frame_check.get("topology", "")
        # Parse tag info from topology string
        # N-terminal tags are before (RE5), C-terminal tags after (RE3)

        # Draw small tag arcs adjacent to the insert
        tag_frac_size = 0.015  # small arc for tags

        # N-terminal tags (before insert)
        n_tags = []
        if "[N-His6]" in topology:
            n_tags.append("N-His6")
        if "[Thrombin]" in topology:
            n_tags.append("Thrombin")
        if "[T7-tag]" in topology:
            n_tags.append("T7-tag")

        for i, tag in enumerate(n_tags):
            tag_start = insert_start_frac - (i + 1) * tag_frac_size
            draw_arc(tag_start, tag_frac_size, R, _CM["tag"], lw_feature, zorder=3)
            if i == 0:  # only label the first/nearest tag
                label_at_frac(tag_start + tag_frac_size / 2,
                              " + ".join(n_tags),
                              radius_offset=0.15, fontsize=7.5, color="#6A1B9A")

        # C-terminal tags (after insert)
        c_tags = []
        if "[C-His6]" in topology and ":OUT-OF-FRAME" not in topology:
            c_tags.append("C-His6")
        if "[S-tag]" in topology and ":OUT-OF-FRAME" not in topology:
            c_tags.append("S-tag")

        for i, tag in enumerate(c_tags):
            tag_start = (insert_start_frac + insert_frac + i * tag_frac_size) % 1.0
            draw_arc(tag_start, tag_frac_size, R, _CM["tag"], lw_feature, zorder=3)
            if i == 0:
                label_at_frac((tag_start + tag_frac_size / 2) % 1.0,
                              " + ".join(c_tags),
                              radius_offset=0.15, fontsize=7.5, color="#6A1B9A")

    # ── Signal peptide region (within insert) ─────────────────────────────
    if signal_peptide and signal_peptide.get("has_signal_peptide"):
        cleavage = signal_peptide.get("cleavage_site_estimate", 20)
        if cleavage and insert_len > 0:
            sp_bp = cleavage * 3  # approximate bp for signal peptide
            sp_frac = (sp_bp / total_bp) if total_bp > 0 else 0.02
            sp_frac = min(sp_frac, insert_frac * 0.4)  # max 40% of insert arc
            # Draw signal peptide as a distinct color within the insert arc
            draw_arc(insert_start_frac, sp_frac, R, _CM["signal_pep"], lw_insert - 2, zorder=4)
            # Add "scissors" annotation for cleavage site
            sp_end_frac = (insert_start_frac + sp_frac) % 1.0
            tick_at_frac(sp_end_frac, R, length=0.10, color="#D35400", lw=2)
            sp_angle = np.radians(frac_to_angle(sp_end_frac))
            sx = (R + 0.08) * np.cos(sp_angle)
            sy = (R + 0.08) * np.sin(sp_angle)
            ax.text(sx, sy, "\u2702", fontsize=14, ha="center", va="center",
                    color="#D35400", zorder=10)
            label_at_frac((insert_start_frac + sp_frac / 2) % 1.0,
                          f"SP ({cleavage} aa)",
                          radius_offset=-0.18, fontsize=7.5, color="#D35400",
                          fontweight="bold")

    # ── Backbone features ────────────────────────────────────────────────
    for feat_name, f_start, f_span, f_color, _ in backbone["features"]:
        draw_arc(f_start, f_span, R, f_color, lw_feature, zorder=2)
        label_at_frac(f_start + f_span / 2, feat_name,
                      radius_offset=0.18, fontsize=7.5, color="#444")

    # ── Center text ──────────────────────────────────────────────────────
    ax.text(0, 0.12, vector_name, fontsize=16, fontweight="bold",
            ha="center", va="center", fontfamily=_FONT_SANS, color="#1A237E")
    ax.text(0, -0.05, f"+ {gene_name}", fontsize=12,
            ha="center", va="center", fontfamily=_FONT_SANS, color="#2E7D32")
    ax.text(0, -0.22, f"{total_bp:,} bp", fontsize=11,
            ha="center", va="center", fontfamily=_FONT, color="#555")

    # ── Frame status indicator ───────────────────────────────────────────
    if frame_check:
        in5 = frame_check.get("in_frame_5prime", False)
        in3 = frame_check.get("in_frame_3prime", False)
        mark5 = "\u2713" if in5 else "\u2717"
        mark3 = "\u2713" if in3 else "\u2717"
        c5 = _CM["insert"] if in5 else _CM["re_site"]
        c3 = _CM["insert"] if in3 else _CM["re_site"]

        ax.text(0, -0.42, f"5' frame: {mark5}  |  3' frame: {mark3}",
                fontsize=9, ha="center", va="center", fontfamily=_FONT_SANS,
                color="#555")

    # ── Direction arrow (clockwise) ──────────────────────────────────────
    arrow_frac = 0.5  # bottom of circle
    arrow_angle = np.radians(frac_to_angle(arrow_frac))
    ax.annotate("",
                xy=((R + 0.02) * np.cos(arrow_angle + 0.05),
                    (R + 0.02) * np.sin(arrow_angle + 0.05)),
                xytext=((R + 0.02) * np.cos(arrow_angle - 0.05),
                        (R + 0.02) * np.sin(arrow_angle - 0.05)),
                arrowprops=dict(arrowstyle="->", color="#888", lw=1.5),
                zorder=5)

    # ── Date ─────────────────────────────────────────────────────────────
    ax.text(0, -0.58, f"Generated: {date.today().isoformat()}",
            fontsize=7, ha="center", va="center", fontfamily=_FONT,
            color="#999")

    # ── Legend ────────────────────────────────────────────────────────────
    legend_items = [
        (_CM["insert"], "Insert CDS"),
        (_CM["tag"], "Fusion tag"),
        (_CM["re_site"], "Promoter / Resistance"),
        (_CM["origin"], "Origin of replication"),
        (_CM["regulatory"], "Regulatory"),
    ]
    if signal_peptide and signal_peptide.get("has_signal_peptide"):
        legend_items.insert(1, (_CM["signal_pep"], "Signal peptide"))

    legend_y = -1.35
    legend_x_start = -1.2
    for i, (lc, lt) in enumerate(legend_items):
        x = legend_x_start + i * 0.55
        ax.add_patch(FancyBboxPatch(
            (x, legend_y - 0.03), 0.08, 0.06,
            boxstyle="round,pad=0.01",
            facecolor=lc, edgecolor="#999", linewidth=0.5,
        ))
        ax.text(x + 0.10, legend_y, lt, fontsize=6.5,
                fontfamily=_FONT_SANS, va="center", color="#555")

    fig.savefig(str(output_path), dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return output_path
