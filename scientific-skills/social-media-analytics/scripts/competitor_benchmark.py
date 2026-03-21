#!/usr/bin/env python3
"""
Competitor Benchmarking Tool

Compares brand metrics against competitor profiles and industry benchmarks.
Produces gap analysis, opportunity scores, and actionable insights.

Usage:
    python competitor_benchmark.py \\
        --brand-metrics metrics.json \\
        --competitor-profiles @competitor1,@competitor2 \\
        --output benchmark.json

    python competitor_benchmark.py \\
        --brand-metrics metrics.json \\
        --competitor-data competitor_metrics.json \\
        --industry saas \\
        --output benchmark.json

Arguments:
    --brand-metrics         Path to brand metrics JSON (from collect_metrics.py)
    --competitor-profiles   Comma-separated competitor handles (for collection guidance)
    --competitor-data       Path to competitor metrics JSON (pre-collected)
    --industry              Industry for benchmark comparison (default: none)
    --output                Output JSON path (default: benchmark.json)
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Industry benchmarks (mirrors references/benchmarks_by_industry.md)
# ---------------------------------------------------------------------------

INDUSTRY_BENCHMARKS = {
    "ecommerce": {
        "facebook":  {"engagement_rate": 0.8, "follower_growth_rate": 0.5, "ctr": 1.2, "posts_per_week": 6},
        "instagram": {"engagement_rate": 1.8, "follower_growth_rate": 1.2, "ctr": 0.8, "posts_per_week": 6},
        "linkedin":  {"engagement_rate": 1.5, "follower_growth_rate": 1.0, "ctr": 1.0, "posts_per_week": 3.5},
        "tiktok":    {"engagement_rate": 5.2, "follower_growth_rate": 4.5, "ctr": 0.6, "posts_per_week": 7.5},
    },
    "saas": {
        "facebook":  {"engagement_rate": 0.5, "follower_growth_rate": 0.3, "ctr": 1.5, "posts_per_week": 4},
        "instagram": {"engagement_rate": 1.2, "follower_growth_rate": 0.8, "ctr": 0.5, "posts_per_week": 4},
        "linkedin":  {"engagement_rate": 2.8, "follower_growth_rate": 1.5, "ctr": 2.2, "posts_per_week": 5},
        "tiktok":    {"engagement_rate": 3.8, "follower_growth_rate": 3.0, "ctr": 0.4, "posts_per_week": 4},
    },
    "healthcare": {
        "facebook":  {"engagement_rate": 0.7, "follower_growth_rate": 0.4, "ctr": 0.9, "posts_per_week": 4},
        "instagram": {"engagement_rate": 2.5, "follower_growth_rate": 1.5, "ctr": 0.6, "posts_per_week": 5},
        "linkedin":  {"engagement_rate": 1.8, "follower_growth_rate": 0.8, "ctr": 1.2, "posts_per_week": 2.5},
        "tiktok":    {"engagement_rate": 6.0, "follower_growth_rate": 5.0, "ctr": 0.3, "posts_per_week": 5},
    },
    "finance": {
        "facebook":  {"engagement_rate": 0.4, "follower_growth_rate": 0.3, "ctr": 1.8, "posts_per_week": 4},
        "instagram": {"engagement_rate": 1.0, "follower_growth_rate": 0.6, "ctr": 0.4, "posts_per_week": 3.5},
        "linkedin":  {"engagement_rate": 2.2, "follower_growth_rate": 1.2, "ctr": 2.5, "posts_per_week": 5},
        "tiktok":    {"engagement_rate": 4.0, "follower_growth_rate": 3.5, "ctr": 0.3, "posts_per_week": 3},
    },
    "food": {
        "facebook":  {"engagement_rate": 1.0, "follower_growth_rate": 0.6, "ctr": 0.8, "posts_per_week": 6},
        "instagram": {"engagement_rate": 3.2, "follower_growth_rate": 2.0, "ctr": 0.5, "posts_per_week": 7.5},
        "linkedin":  {"engagement_rate": 1.0, "follower_growth_rate": 0.5, "ctr": 0.6, "posts_per_week": 1.5},
        "tiktok":    {"engagement_rate": 7.5, "follower_growth_rate": 6.0, "ctr": 0.4, "posts_per_week": 7.5},
    },
    "fashion": {
        "facebook":  {"engagement_rate": 0.6, "follower_growth_rate": 0.4, "ctr": 0.7, "posts_per_week": 5},
        "instagram": {"engagement_rate": 2.8, "follower_growth_rate": 1.8, "ctr": 0.6, "posts_per_week": 10.5},
        "linkedin":  {"engagement_rate": 0.8, "follower_growth_rate": 0.4, "ctr": 0.5, "posts_per_week": 1.5},
        "tiktok":    {"engagement_rate": 6.5, "follower_growth_rate": 5.5, "ctr": 0.5, "posts_per_week": 10.5},
    },
    "education": {
        "facebook":  {"engagement_rate": 0.6, "follower_growth_rate": 0.4, "ctr": 1.0, "posts_per_week": 4},
        "instagram": {"engagement_rate": 2.0, "follower_growth_rate": 1.0, "ctr": 0.5, "posts_per_week": 4},
        "linkedin":  {"engagement_rate": 2.5, "follower_growth_rate": 1.3, "ctr": 1.8, "posts_per_week": 4},
        "tiktok":    {"engagement_rate": 5.5, "follower_growth_rate": 4.0, "ctr": 0.3, "posts_per_week": 5},
    },
    "realestate": {
        "facebook":  {"engagement_rate": 0.7, "follower_growth_rate": 0.5, "ctr": 1.4, "posts_per_week": 6},
        "instagram": {"engagement_rate": 1.8, "follower_growth_rate": 1.0, "ctr": 0.7, "posts_per_week": 5.5},
        "linkedin":  {"engagement_rate": 1.5, "follower_growth_rate": 0.8, "ctr": 1.5, "posts_per_week": 2.5},
        "tiktok":    {"engagement_rate": 5.0, "follower_growth_rate": 3.5, "ctr": 0.4, "posts_per_week": 5},
    },
    "travel": {
        "facebook":  {"engagement_rate": 0.8, "follower_growth_rate": 0.5, "ctr": 0.9, "posts_per_week": 5},
        "instagram": {"engagement_rate": 3.0, "follower_growth_rate": 1.5, "ctr": 0.5, "posts_per_week": 7.5},
        "linkedin":  {"engagement_rate": 1.2, "follower_growth_rate": 0.6, "ctr": 0.8, "posts_per_week": 2},
        "tiktok":    {"engagement_rate": 6.8, "follower_growth_rate": 5.0, "ctr": 0.3, "posts_per_week": 7.5},
    },
    "fitness": {
        "facebook":  {"engagement_rate": 0.7, "follower_growth_rate": 0.4, "ctr": 0.8, "posts_per_week": 5},
        "instagram": {"engagement_rate": 2.8, "follower_growth_rate": 1.5, "ctr": 0.5, "posts_per_week": 10.5},
        "linkedin":  {"engagement_rate": 1.0, "follower_growth_rate": 0.5, "ctr": 0.7, "posts_per_week": 1.5},
        "tiktok":    {"engagement_rate": 7.0, "follower_growth_rate": 5.5, "ctr": 0.3, "posts_per_week": 10.5},
    },
    "b2b": {
        "facebook":  {"engagement_rate": 0.3, "follower_growth_rate": 0.2, "ctr": 1.0, "posts_per_week": 2.5},
        "instagram": {"engagement_rate": 0.8, "follower_growth_rate": 0.5, "ctr": 0.3, "posts_per_week": 3},
        "linkedin":  {"engagement_rate": 2.5, "follower_growth_rate": 1.5, "ctr": 2.8, "posts_per_week": 5.5},
        "tiktok":    {"engagement_rate": 3.0, "follower_growth_rate": 2.0, "ctr": 0.2, "posts_per_week": 2},
    },
    "nonprofit": {
        "facebook":  {"engagement_rate": 1.0, "follower_growth_rate": 0.5, "ctr": 1.5, "posts_per_week": 4},
        "instagram": {"engagement_rate": 2.5, "follower_growth_rate": 1.2, "ctr": 0.6, "posts_per_week": 4},
        "linkedin":  {"engagement_rate": 2.0, "follower_growth_rate": 1.0, "ctr": 1.5, "posts_per_week": 3},
        "tiktok":    {"engagement_rate": 5.5, "follower_growth_rate": 4.0, "ctr": 0.3, "posts_per_week": 3.5},
    },
}


# ---------------------------------------------------------------------------
# Gap analysis
# ---------------------------------------------------------------------------

def compute_gap(brand_value: Optional[float], benchmark_value: Optional[float]) -> dict:
    """Compute the gap between brand metric and a benchmark."""
    if brand_value is None or benchmark_value is None:
        return {"status": "insufficient_data"}

    absolute_gap = round(brand_value - benchmark_value, 4)
    relative_gap_pct = round((absolute_gap / benchmark_value) * 100, 2) if benchmark_value != 0 else None

    if relative_gap_pct is not None:
        if relative_gap_pct >= 20:
            assessment = "significantly_above"
        elif relative_gap_pct >= 0:
            assessment = "above"
        elif relative_gap_pct >= -20:
            assessment = "slightly_below"
        else:
            assessment = "significantly_below"
    else:
        assessment = "unknown"

    return {
        "brand_value": brand_value,
        "benchmark_value": benchmark_value,
        "absolute_gap": absolute_gap,
        "relative_gap_pct": relative_gap_pct,
        "assessment": assessment,
    }


def compute_opportunity_score(gaps: dict) -> float:
    """
    Compute an opportunity score (0-100) based on how many metrics are
    below benchmark and by how much.
    """
    below_gaps = []
    for metric, gap_data in gaps.items():
        if isinstance(gap_data, dict) and gap_data.get("relative_gap_pct") is not None:
            if gap_data["relative_gap_pct"] < 0:
                below_gaps.append(abs(gap_data["relative_gap_pct"]))

    if not below_gaps:
        return 0.0

    # Score: weighted by magnitude of underperformance, capped at 100
    raw_score = sum(min(g, 100) for g in below_gaps) / len(below_gaps)
    return round(min(raw_score, 100), 1)


# ---------------------------------------------------------------------------
# Competitor data collection guidance
# ---------------------------------------------------------------------------

def competitor_collection_guidance(handles: list[str]) -> dict:
    """Generate instructions for collecting competitor data."""
    guidance = {}
    for handle in handles:
        guidance[handle] = {
            "instagram": (
                f"1. Navigate to https://www.instagram.com/{handle.lstrip('@')}/\n"
                f"2. Note follower count, following count, post count.\n"
                f"3. Review the last 12-20 posts and record: likes, comments, post type.\n"
                f"4. Calculate average engagement rate: (avg likes + avg comments) / followers * 100.\n"
                f"5. Note posting frequency (posts per week over last 30 days)."
            ),
            "linkedin": (
                f"1. Navigate to the company page for {handle}.\n"
                f"2. Note follower count.\n"
                f"3. Review the last 10-15 posts: reactions, comments, reposts.\n"
                f"4. Calculate average engagement per post.\n"
                f"5. Note posting frequency."
            ),
            "tiktok": (
                f"1. Navigate to https://www.tiktok.com/@{handle.lstrip('@')}\n"
                f"2. Note follower count, total likes.\n"
                f"3. Review the last 10-20 videos: views, likes, comments, shares.\n"
                f"4. Calculate average engagement rate.\n"
                f"5. Note posting frequency."
            ),
            "facebook": (
                f"1. Navigate to the Facebook page for {handle}.\n"
                f"2. Note page likes/followers.\n"
                f"3. Review the last 10-15 posts: reactions, comments, shares.\n"
                f"4. Calculate average engagement rate.\n"
                f"5. Note posting frequency."
            ),
        }
    return guidance


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Benchmark brand against competitors.")
    parser.add_argument("--brand-metrics", type=str, required=True,
                        help="Path to brand metrics JSON")
    parser.add_argument("--competitor-profiles", type=str, default=None,
                        help="Comma-separated competitor handles")
    parser.add_argument("--competitor-data", type=str, default=None,
                        help="Path to pre-collected competitor metrics JSON")
    parser.add_argument("--industry", type=str, default=None,
                        choices=list(INDUSTRY_BENCHMARKS.keys()),
                        help="Industry for benchmark comparison")
    parser.add_argument("--output", type=str, default="benchmark.json",
                        help="Output JSON path")
    args = parser.parse_args()

    with open(args.brand_metrics, "r", encoding="utf-8") as f:
        brand_data = json.load(f)

    competitor_data = None
    if args.competitor_data:
        with open(args.competitor_data, "r", encoding="utf-8") as f:
            competitor_data = json.load(f)

    competitor_handles = []
    if args.competitor_profiles:
        competitor_handles = [h.strip() for h in args.competitor_profiles.split(",")]

    result = {
        "benchmark_metadata": {
            "analyzed_at": datetime.utcnow().isoformat() + "Z",
            "brand": brand_data.get("collection_metadata", {}).get("account_name"),
            "industry": args.industry,
            "competitors": competitor_handles,
        },
        "industry_benchmark": {},
        "competitor_comparison": {},
        "collection_guidance": {},
    }

    # --- Industry benchmark comparison ---
    if args.industry and args.industry in INDUSTRY_BENCHMARKS:
        benchmarks = INDUSTRY_BENCHMARKS[args.industry]
        for platform, pdata in brand_data.get("platforms", {}).items():
            if platform not in benchmarks:
                continue
            bench = benchmarks[platform]
            acct = pdata.get("account_metrics", {})
            followers = acct.get("followers")
            engagement_total = acct.get("engagement_total")
            impressions = acct.get("impressions")

            brand_er = (engagement_total / followers * 100) if followers and engagement_total else None
            brand_ctr = None  # Would need clicks/impressions from account level

            gaps = {
                "engagement_rate": compute_gap(brand_er, bench["engagement_rate"]),
            }
            opp_score = compute_opportunity_score(gaps)

            result["industry_benchmark"][platform] = {
                "industry": args.industry,
                "gaps": gaps,
                "opportunity_score": opp_score,
                "industry_posting_frequency": bench["posts_per_week"],
            }

    # --- Competitor comparison ---
    if competitor_data:
        for platform, pdata in brand_data.get("platforms", {}).items():
            brand_acct = pdata.get("account_metrics", {})
            brand_followers = brand_acct.get("followers")
            brand_engagement = brand_acct.get("engagement_total")
            brand_er = (brand_engagement / brand_followers * 100) if brand_followers and brand_engagement else None

            comp_comparisons = []
            for comp_name, comp_platforms in competitor_data.items():
                comp_pdata = comp_platforms.get(platform, {})
                comp_er = comp_pdata.get("engagement_rate")
                comp_followers = comp_pdata.get("followers")
                comp_growth = comp_pdata.get("follower_growth_rate")

                comp_comparisons.append({
                    "competitor": comp_name,
                    "engagement_rate_gap": compute_gap(brand_er, comp_er),
                    "follower_comparison": {
                        "brand_followers": brand_followers,
                        "competitor_followers": comp_followers,
                    },
                })

            result["competitor_comparison"][platform] = comp_comparisons

    # --- Collection guidance for competitors ---
    if competitor_handles and not competitor_data:
        result["collection_guidance"] = competitor_collection_guidance(competitor_handles)
        print("\nCompetitor data not yet collected. See 'collection_guidance' in output for instructions.")

    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, default=str)

    print(f"\nBenchmark written to {args.output}")
    if args.industry:
        print(f"Industry: {args.industry}")
        for platform, bench in result.get("industry_benchmark", {}).items():
            print(f"  {platform}: opportunity score = {bench['opportunity_score']}")


if __name__ == "__main__":
    main()
