"""RE Cloning primer design report image generator (matplotlib).

design_re_cloning_primers() 결과를 시각적 리포트(PNG)로 생성.
- 프라이머 구조 다이어그램 (protection / RE / spacer / annealing 색상 구분)
- QC 요약 테이블
- 클로닝 construct 맵
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

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
