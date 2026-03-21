#!/usr/bin/env python3
"""
Social Media Performance Analyzer

Reads collected metrics JSON and computes KPIs, ranks content, identifies
trends, and flags top/bottom performers.

Usage:
    python analyze_performance.py --metrics metrics.json --output analysis.json
    python analyze_performance.py --metrics metrics.json --output analysis.json --top-pct 10 --bottom-pct 10

Arguments:
    --metrics       Path to metrics JSON (output of collect_metrics.py)
    --output        Output JSON file path (default: analysis.json)
    --top-pct       Percentile threshold for top performers (default: 10)
    --bottom-pct    Percentile threshold for bottom performers (default: 10)
"""

import argparse
import json
import math
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Optional


# ---------------------------------------------------------------------------
# KPI computation helpers
# ---------------------------------------------------------------------------

def safe_div(numerator: Optional[float], denominator: Optional[float]) -> Optional[float]:
    """Divide two numbers, returning None if either is None or denominator is zero."""
    if numerator is None or denominator is None or denominator == 0:
        return None
    return numerator / denominator


def pct(value: Optional[float]) -> Optional[float]:
    """Multiply by 100 to convert ratio to percentage."""
    return round(value * 100, 4) if value is not None else None


def compute_account_kpis(acct: dict) -> dict:
    """Compute account-level KPIs from raw account metrics."""
    m = acct.get("account_metrics", {})
    followers = m.get("followers")
    reach = m.get("reach")
    impressions = m.get("impressions")
    engagement_total = m.get("engagement_total")
    profile_visits = m.get("profile_visits")
    website_clicks = m.get("website_clicks")
    followers_gained = m.get("followers_gained")
    followers_lost = m.get("followers_lost")

    kpis = {}

    # Engagement Rate (by followers)
    kpis["engagement_rate_by_followers"] = pct(safe_div(engagement_total, followers))
    # Engagement Rate (by reach)
    kpis["engagement_rate_by_reach"] = pct(safe_div(engagement_total, reach))

    # Follower Growth Rate
    net_new = None
    if followers_gained is not None and followers_lost is not None:
        net_new = followers_gained - followers_lost
    starting_followers = followers - net_new if followers and net_new is not None else None
    kpis["follower_growth_rate"] = pct(safe_div(net_new, starting_followers))

    # Profile Visit Rate
    kpis["profile_visit_rate"] = pct(safe_div(profile_visits, impressions))

    # Website Click Rate
    kpis["website_click_rate"] = pct(safe_div(website_clicks, profile_visits))

    # Impression-to-Reach Ratio
    kpis["impression_reach_ratio"] = round(safe_div(impressions, reach), 2) if safe_div(impressions, reach) is not None else None

    # Comment Rate
    kpis["comment_rate"] = pct(safe_div(m.get("comments_total"), followers))
    # Share Rate
    kpis["share_rate"] = pct(safe_div(m.get("shares_total"), followers))
    # Save Rate (by reach)
    kpis["save_rate"] = pct(safe_div(m.get("saves_total"), reach))
    # Amplification Rate
    kpis["amplification_rate"] = pct(safe_div(m.get("shares_total"), followers))
    # Virality Rate
    kpis["virality_rate"] = pct(safe_div(m.get("shares_total"), impressions))

    return kpis


def compute_post_kpis(post: dict, followers: Optional[float]) -> dict:
    """Compute per-post KPIs."""
    reach = post.get("reach")
    impressions = post.get("impressions")
    likes = post.get("likes")
    comments = post.get("comments")
    shares = post.get("shares")
    saves = post.get("saves")
    clicks = post.get("clicks")

    total_engagement = sum(v for v in [likes, comments, shares, saves] if v is not None)

    kpis = {}
    kpis["total_engagement"] = total_engagement
    kpis["engagement_rate_by_reach"] = pct(safe_div(total_engagement, reach))
    kpis["engagement_rate_by_followers"] = pct(safe_div(total_engagement, followers))
    kpis["ctr"] = pct(safe_div(clicks, impressions))
    kpis["save_rate"] = pct(safe_div(saves, reach))
    kpis["share_rate"] = pct(safe_div(shares, reach))
    kpis["comment_rate"] = pct(safe_div(comments, reach))
    kpis["virality_rate"] = pct(safe_div(shares, impressions))

    # Video-specific
    kpis["video_completion_rate"] = post.get("video_completion_rate")
    kpis["video_avg_watch_seconds"] = post.get("video_avg_watch_seconds")

    return kpis


# ---------------------------------------------------------------------------
# Ranking and flagging
# ---------------------------------------------------------------------------

def percentile_value(sorted_values: list[float], pct_val: float) -> float:
    """Compute the percentile value from a sorted list."""
    if not sorted_values:
        return 0
    idx = (pct_val / 100) * (len(sorted_values) - 1)
    lower = int(math.floor(idx))
    upper = int(math.ceil(idx))
    if lower == upper:
        return sorted_values[lower]
    frac = idx - lower
    return sorted_values[lower] * (1 - frac) + sorted_values[upper] * frac


def rank_posts(posts_with_kpis: list[dict], top_pct: float, bottom_pct: float) -> list[dict]:
    """Rank posts by engagement rate and flag top/bottom performers."""
    valid = [p for p in posts_with_kpis if p["kpis"].get("engagement_rate_by_reach") is not None]
    if not valid:
        return posts_with_kpis

    sorted_er = sorted([p["kpis"]["engagement_rate_by_reach"] for p in valid])
    top_threshold = percentile_value(sorted_er, 100 - top_pct)
    bottom_threshold = percentile_value(sorted_er, bottom_pct)

    for post in posts_with_kpis:
        er = post["kpis"].get("engagement_rate_by_reach")
        if er is None:
            post["performance_flag"] = "insufficient_data"
        elif er >= top_threshold:
            post["performance_flag"] = "top_performer"
        elif er <= bottom_threshold:
            post["performance_flag"] = "bottom_performer"
        else:
            post["performance_flag"] = "average"

    posts_with_kpis.sort(
        key=lambda p: p["kpis"].get("engagement_rate_by_reach") or 0,
        reverse=True
    )
    for i, post in enumerate(posts_with_kpis):
        post["rank"] = i + 1

    return posts_with_kpis


# ---------------------------------------------------------------------------
# Trend analysis
# ---------------------------------------------------------------------------

def detect_trends(posts_with_kpis: list[dict]) -> dict:
    """Analyze whether key metrics are improving or declining over time."""
    dated = [
        p for p in posts_with_kpis
        if p.get("published_at") and p["kpis"].get("engagement_rate_by_reach") is not None
    ]
    if len(dated) < 4:
        return {"status": "insufficient_data", "detail": "Need at least 4 dated posts for trend analysis."}

    dated.sort(key=lambda p: p["published_at"])
    mid = len(dated) // 2
    first_half = dated[:mid]
    second_half = dated[mid:]

    def avg_metric(posts, metric_key):
        vals = [p["kpis"][metric_key] for p in posts if p["kpis"].get(metric_key) is not None]
        return sum(vals) / len(vals) if vals else None

    trends = {}
    for metric in ["engagement_rate_by_reach", "ctr", "save_rate", "share_rate", "comment_rate"]:
        first_avg = avg_metric(first_half, metric)
        second_avg = avg_metric(second_half, metric)
        if first_avg is not None and second_avg is not None and first_avg > 0:
            change_pct = ((second_avg - first_avg) / first_avg) * 100
            if change_pct > 10:
                direction = "improving"
            elif change_pct < -10:
                direction = "declining"
            else:
                direction = "stable"
            trends[metric] = {
                "first_half_avg": round(first_avg, 4),
                "second_half_avg": round(second_avg, 4),
                "change_pct": round(change_pct, 2),
                "direction": direction,
            }
        else:
            trends[metric] = {"status": "insufficient_data"}

    return trends


# ---------------------------------------------------------------------------
# Content pattern analysis
# ---------------------------------------------------------------------------

def analyze_content_patterns(posts_with_kpis: list[dict]) -> dict:
    """Identify which content attributes correlate with high performance."""
    by_type = {}
    for post in posts_with_kpis:
        ctype = post.get("content_type") or "unknown"
        er = post["kpis"].get("engagement_rate_by_reach")
        if er is None:
            continue
        by_type.setdefault(ctype, []).append(er)

    type_summary = {}
    for ctype, ers in by_type.items():
        type_summary[ctype] = {
            "count": len(ers),
            "avg_engagement_rate": round(sum(ers) / len(ers), 4) if ers else None,
            "max_engagement_rate": round(max(ers), 4) if ers else None,
            "min_engagement_rate": round(min(ers), 4) if ers else None,
        }

    # Sort by average ER descending
    ranked_types = sorted(type_summary.items(), key=lambda x: x[1]["avg_engagement_rate"] or 0, reverse=True)

    return {
        "by_content_type": dict(ranked_types),
        "best_format": ranked_types[0][0] if ranked_types else None,
        "worst_format": ranked_types[-1][0] if ranked_types else None,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Analyze social media performance.")
    parser.add_argument("--metrics", type=str, required=True, help="Input metrics JSON")
    parser.add_argument("--output", type=str, default="analysis.json", help="Output analysis JSON")
    parser.add_argument("--top-pct", type=float, default=10, help="Top performer percentile")
    parser.add_argument("--bottom-pct", type=float, default=10, help="Bottom performer percentile")
    args = parser.parse_args()

    with open(args.metrics, "r", encoding="utf-8") as f:
        data = json.load(f)

    analysis = {
        "analysis_metadata": {
            "analyzed_at": datetime.utcnow().isoformat() + "Z",
            "source_file": args.metrics,
            "period": data.get("collection_metadata", {}).get("period"),
            "period_start": data.get("collection_metadata", {}).get("period_start"),
            "period_end": data.get("collection_metadata", {}).get("period_end"),
            "account_name": data.get("collection_metadata", {}).get("account_name"),
        },
        "platforms": {},
    }

    for platform, pdata in data.get("platforms", {}).items():
        followers = (pdata.get("account_metrics") or {}).get("followers")

        # Account-level KPIs
        account_kpis = compute_account_kpis(pdata)

        # Post-level KPIs
        posts_analyzed = []
        for post in pdata.get("posts", []):
            post_kpis = compute_post_kpis(post, followers)
            posts_analyzed.append({
                "post_id": post.get("post_id"),
                "published_at": post.get("published_at"),
                "content_type": post.get("content_type"),
                "caption_length": post.get("caption_length"),
                "hashtag_count": post.get("hashtag_count"),
                "kpis": post_kpis,
            })

        # Rank and flag
        posts_analyzed = rank_posts(posts_analyzed, args.top_pct, args.bottom_pct)

        # Trends
        trends = detect_trends(posts_analyzed)

        # Content patterns
        patterns = analyze_content_patterns(posts_analyzed)

        # Summary counts
        top_count = sum(1 for p in posts_analyzed if p.get("performance_flag") == "top_performer")
        bottom_count = sum(1 for p in posts_analyzed if p.get("performance_flag") == "bottom_performer")

        analysis["platforms"][platform] = {
            "account_kpis": account_kpis,
            "post_count": len(posts_analyzed),
            "top_performers_count": top_count,
            "bottom_performers_count": bottom_count,
            "posts": posts_analyzed,
            "trends": trends,
            "content_patterns": patterns,
        }

    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(analysis, f, indent=2, default=str)

    print(f"Analysis written to {args.output}")
    print(f"\nPlatforms analyzed: {list(analysis['platforms'].keys())}")
    for platform, panalysis in analysis["platforms"].items():
        print(f"\n  {platform.upper()}:")
        print(f"    Posts analyzed: {panalysis['post_count']}")
        print(f"    Top performers: {panalysis['top_performers_count']}")
        print(f"    Bottom performers: {panalysis['bottom_performers_count']}")
        er = panalysis["account_kpis"].get("engagement_rate_by_followers")
        if er is not None:
            print(f"    Engagement Rate (by followers): {er}%")
        best = panalysis["content_patterns"].get("best_format")
        if best:
            print(f"    Best content format: {best}")


if __name__ == "__main__":
    main()
