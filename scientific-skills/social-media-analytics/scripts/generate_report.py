#!/usr/bin/env python3
"""
Social Media Analytics Report Generator

Produces a comprehensive analytics report in Markdown or DOCX format from
analysis JSON output.

Usage:
    python generate_report.py --analysis analysis.json --format md --output report.md
    python generate_report.py --analysis analysis.json --format docx --brand-dna brand.json --output report.docx
    python generate_report.py --analysis analysis.json --benchmark benchmark.json --format md --output report.md

Arguments:
    --analysis      Path to analysis JSON (from analyze_performance.py)
    --benchmark     Optional path to benchmark JSON (from competitor_benchmark.py)
    --format        Output format: md or docx (default: md)
    --brand-dna     Optional brand configuration JSON (name, logo, colors)
    --output        Output file path (default: social_media_report.md)
    --title         Report title (default: "Social Media Analytics Report")
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Markdown report generation
# ---------------------------------------------------------------------------

def generate_markdown_report(analysis: dict, benchmark: Optional[dict], title: str) -> str:
    """Generate a full Markdown analytics report."""
    meta = analysis.get("analysis_metadata", {})
    lines = []

    # --- Title & metadata ---
    lines.append(f"# {title}")
    lines.append("")
    lines.append(f"**Account**: {meta.get('account_name', 'N/A')}  ")
    lines.append(f"**Period**: {meta.get('period_start', '?')} to {meta.get('period_end', '?')} ({meta.get('period', '?')})  ")
    lines.append(f"**Generated**: {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}  ")
    lines.append("")

    # --- Executive Summary ---
    lines.append("## Executive Summary")
    lines.append("")
    lines.append(_build_executive_summary(analysis, benchmark))
    lines.append("")

    # --- Platform Breakdown ---
    for platform, pdata in analysis.get("platforms", {}).items():
        lines.append(f"## {platform.title()} Performance")
        lines.append("")

        # Account KPIs table
        kpis = pdata.get("account_kpis", {})
        lines.append("### Key Performance Indicators")
        lines.append("")
        lines.append("| Metric | Value |")
        lines.append("|---|---|")
        kpi_labels = {
            "engagement_rate_by_followers": "Engagement Rate (by followers)",
            "engagement_rate_by_reach": "Engagement Rate (by reach)",
            "follower_growth_rate": "Follower Growth Rate",
            "comment_rate": "Comment Rate",
            "share_rate": "Share Rate",
            "save_rate": "Save Rate",
            "amplification_rate": "Amplification Rate",
            "virality_rate": "Virality Rate",
            "profile_visit_rate": "Profile Visit Rate",
            "website_click_rate": "Website Click Rate",
            "impression_reach_ratio": "Impression/Reach Ratio",
        }
        for key, label in kpi_labels.items():
            val = kpis.get(key)
            if val is not None:
                suffix = "%" if key != "impression_reach_ratio" else "x"
                lines.append(f"| {label} | {val}{suffix} |")
            else:
                lines.append(f"| {label} | -- |")
        lines.append("")

        # Trends
        trends = pdata.get("trends", {})
        if trends and trends.get("status") != "insufficient_data":
            lines.append("### Trend Analysis")
            lines.append("")
            lines.append("| Metric | Direction | Change |")
            lines.append("|---|---|---|")
            for metric, tdata in trends.items():
                if isinstance(tdata, dict) and "direction" in tdata:
                    direction_icon = {
                        "improving": "Improving",
                        "declining": "Declining",
                        "stable": "Stable",
                    }.get(tdata["direction"], tdata["direction"])
                    lines.append(f"| {metric.replace('_', ' ').title()} | {direction_icon} | {tdata.get('change_pct', '--')}% |")
            lines.append("")

        # Content patterns
        patterns = pdata.get("content_patterns", {})
        by_type = patterns.get("by_content_type", {})
        if by_type:
            lines.append("### Content Performance by Format")
            lines.append("")
            lines.append("| Format | Posts | Avg ER | Best ER |")
            lines.append("|---|---|---|---|")
            for ctype, stats in by_type.items():
                avg_er = f"{stats['avg_engagement_rate']}%" if stats.get("avg_engagement_rate") is not None else "--"
                max_er = f"{stats['max_engagement_rate']}%" if stats.get("max_engagement_rate") is not None else "--"
                lines.append(f"| {ctype.title()} | {stats.get('count', 0)} | {avg_er} | {max_er} |")
            lines.append("")

        # Top performers
        top_posts = [p for p in pdata.get("posts", []) if p.get("performance_flag") == "top_performer"]
        if top_posts:
            lines.append("### Top Performing Content")
            lines.append("")
            for post in top_posts[:5]:
                pid = post.get("post_id") or "?"
                ctype = post.get("content_type") or "?"
                er = post["kpis"].get("engagement_rate_by_reach")
                er_str = f"{er}%" if er is not None else "--"
                lines.append(f"- **{pid}** ({ctype}): ER = {er_str}")
            lines.append("")

        # Bottom performers
        bottom_posts = [p for p in pdata.get("posts", []) if p.get("performance_flag") == "bottom_performer"]
        if bottom_posts:
            lines.append("### Underperforming Content")
            lines.append("")
            for post in bottom_posts[:5]:
                pid = post.get("post_id") or "?"
                ctype = post.get("content_type") or "?"
                er = post["kpis"].get("engagement_rate_by_reach")
                er_str = f"{er}%" if er is not None else "--"
                lines.append(f"- **{pid}** ({ctype}): ER = {er_str}")
            lines.append("")

    # --- Benchmark Section ---
    if benchmark:
        lines.append("## Competitor & Industry Benchmark")
        lines.append("")

        industry_bench = benchmark.get("industry_benchmark", {})
        if industry_bench:
            industry_name = benchmark.get("benchmark_metadata", {}).get("industry", "N/A")
            lines.append(f"### Industry Benchmark: {industry_name.title()}")
            lines.append("")
            lines.append("| Platform | Metric | Brand | Benchmark | Gap | Assessment |")
            lines.append("|---|---|---|---|---|---|")
            for platform, bench_data in industry_bench.items():
                for metric, gap_info in bench_data.get("gaps", {}).items():
                    if isinstance(gap_info, dict) and "brand_value" in gap_info:
                        lines.append(
                            f"| {platform.title()} | {metric.replace('_', ' ').title()} "
                            f"| {gap_info['brand_value']}% "
                            f"| {gap_info['benchmark_value']}% "
                            f"| {gap_info['absolute_gap']:+.2f}% "
                            f"| {gap_info['assessment'].replace('_', ' ').title()} |"
                        )
            lines.append("")

        comp_data = benchmark.get("competitor_comparison", {})
        if comp_data:
            lines.append("### Competitor Comparison")
            lines.append("")
            for platform, comps in comp_data.items():
                lines.append(f"**{platform.title()}**:")
                for comp in comps:
                    name = comp.get("competitor", "?")
                    gap = comp.get("engagement_rate_gap", {})
                    if isinstance(gap, dict) and "assessment" in gap:
                        lines.append(f"- vs {name}: {gap['assessment'].replace('_', ' ')}")
                lines.append("")

    # --- Recommendations ---
    lines.append("## Recommendations")
    lines.append("")
    lines.append(_build_recommendations(analysis, benchmark))

    return "\n".join(lines)


def _build_executive_summary(analysis: dict, benchmark: Optional[dict]) -> str:
    """Build the executive summary paragraph."""
    parts = []
    platforms = list(analysis.get("platforms", {}).keys())
    parts.append(f"This report covers performance across {', '.join(p.title() for p in platforms)}.")

    for platform, pdata in analysis.get("platforms", {}).items():
        er = pdata.get("account_kpis", {}).get("engagement_rate_by_followers")
        fgr = pdata.get("account_kpis", {}).get("follower_growth_rate")
        post_count = pdata.get("post_count", 0)
        if er is not None:
            parts.append(
                f"{platform.title()} achieved a {er}% engagement rate across {post_count} posts."
            )
        if fgr is not None:
            parts.append(f"Follower growth rate on {platform.title()} was {fgr}%.")

    # Highlight best platform
    best_platform = None
    best_er = -1
    for platform, pdata in analysis.get("platforms", {}).items():
        er = pdata.get("account_kpis", {}).get("engagement_rate_by_followers")
        if er is not None and er > best_er:
            best_er = er
            best_platform = platform
    if best_platform:
        parts.append(f"{best_platform.title()} was the highest-engaging platform this period.")

    return " ".join(parts)


def _build_recommendations(analysis: dict, benchmark: Optional[dict]) -> str:
    """Generate prioritized recommendations based on analysis."""
    recs = []

    for platform, pdata in analysis.get("platforms", {}).items():
        patterns = pdata.get("content_patterns", {})
        best_format = patterns.get("best_format")
        worst_format = patterns.get("worst_format")
        if best_format and worst_format and best_format != worst_format:
            recs.append(
                f"- **Content Mix ({platform.title()})**: Increase production of "
                f"**{best_format}** content which outperforms other formats. "
                f"Reduce or redesign **{worst_format}** content."
            )

        trends = pdata.get("trends", {})
        for metric, tdata in trends.items():
            if isinstance(tdata, dict) and tdata.get("direction") == "declining":
                recs.append(
                    f"- **Declining {metric.replace('_', ' ').title()} ({platform.title()})**: "
                    f"This metric dropped {abs(tdata.get('change_pct', 0))}% period-over-period. "
                    f"Review recent content changes and audience feedback."
                )

    if benchmark:
        for platform, bench_data in benchmark.get("industry_benchmark", {}).items():
            opp_score = bench_data.get("opportunity_score", 0)
            if opp_score > 30:
                recs.append(
                    f"- **Close Industry Gap ({platform.title()})**: Opportunity score of "
                    f"{opp_score}/100 indicates significant room to improve vs. industry averages. "
                    f"Focus on the metrics with the largest gaps."
                )

    if not recs:
        recs.append("- Performance is on track. Continue current strategy and monitor for changes.")

    return "\n".join(recs)


# ---------------------------------------------------------------------------
# DOCX generation (delegates to document-skills/docx)
# ---------------------------------------------------------------------------

def generate_docx_instructions(analysis: dict, benchmark: Optional[dict],
                                title: str, brand_dna: Optional[dict], output: str) -> str:
    """Return instructions for generating a DOCX report via the docx skill."""
    return (
        f"To generate the DOCX version of this report:\n\n"
        f"1. First generate the Markdown report using this script with --format md.\n"
        f"2. Use the document-skills/docx skill to convert the Markdown to DOCX:\n"
        f"   - Pass the Markdown content as input.\n"
        f"   - If --brand-dna is provided, apply brand colors and logo from: {brand_dna}\n"
        f"   - Output to: {output}\n\n"
        f"3. For charts and visualizations:\n"
        f"   - Use the scientific-visualization skill to generate matplotlib charts.\n"
        f"   - Recommended charts:\n"
        f"     a. Engagement rate by platform (grouped bar chart)\n"
        f"     b. Content format performance comparison (horizontal bar)\n"
        f"     c. Follower growth trend line\n"
        f"     d. Top 10 posts by engagement (bar chart)\n"
        f"     e. Industry benchmark radar chart\n"
        f"   - Embed generated chart images into the DOCX.\n"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate social media analytics report.")
    parser.add_argument("--analysis", type=str, required=True, help="Analysis JSON path")
    parser.add_argument("--benchmark", type=str, default=None, help="Benchmark JSON path")
    parser.add_argument("--format", type=str, default="md", choices=["md", "docx"],
                        help="Output format (default: md)")
    parser.add_argument("--brand-dna", type=str, default=None,
                        help="Brand config JSON (name, logo path, colors)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output file path")
    parser.add_argument("--title", type=str, default="Social Media Analytics Report",
                        help="Report title")
    args = parser.parse_args()

    if args.output is None:
        args.output = f"social_media_report.{args.format}"

    with open(args.analysis, "r", encoding="utf-8") as f:
        analysis = json.load(f)

    benchmark = None
    if args.benchmark:
        with open(args.benchmark, "r", encoding="utf-8") as f:
            benchmark = json.load(f)

    brand_dna = None
    if args.brand_dna:
        with open(args.brand_dna, "r", encoding="utf-8") as f:
            brand_dna = json.load(f)

    if args.format == "md":
        report = generate_markdown_report(analysis, benchmark, args.title)
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(report)
        print(f"Markdown report written to {args.output}")

    elif args.format == "docx":
        # Generate Markdown first, then provide DOCX instructions
        md_output = args.output.replace(".docx", ".md")
        report = generate_markdown_report(analysis, benchmark, args.title)
        with open(md_output, "w", encoding="utf-8") as f:
            f.write(report)
        print(f"Markdown draft written to {md_output}")

        instructions = generate_docx_instructions(analysis, benchmark, args.title, brand_dna, args.output)
        print(f"\n{instructions}")

    print("\nReport generation complete.")
    print("Review the output and adjust recommendations as needed.")


if __name__ == "__main__":
    main()
