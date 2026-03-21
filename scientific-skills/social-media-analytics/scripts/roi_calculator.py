#!/usr/bin/env python3
"""
Social Media ROI Calculator

Calculates return on investment for social media activities by combining
engagement metrics with cost data and attributable revenue.

Usage:
    python roi_calculator.py \\
        --metrics metrics.json \\
        --costs costs.json \\
        --revenue-attributable 12000 \\
        --output roi.json

    python roi_calculator.py \\
        --metrics metrics.json \\
        --costs costs.json \\
        --revenue-attributable 0 \\
        --hourly-rate 50 \\
        --output roi.json

Arguments:
    --metrics               Path to metrics JSON (from collect_metrics.py)
    --costs                 Path to costs JSON (see schema below)
    --revenue-attributable  Total revenue attributable to social media ($)
    --hourly-rate           Hourly rate for time-invested calculation (default: $50)
    --output                Output JSON path (default: roi.json)

Costs JSON schema:
{
    "time_invested_hours": {
        "content_creation": 40,
        "community_management": 20,
        "strategy_planning": 10,
        "analytics_reporting": 5
    },
    "ad_spend": {
        "facebook": 500,
        "instagram": 800,
        "linkedin": 300,
        "tiktok": 200
    },
    "tools_cost": {
        "scheduling_tool": 50,
        "design_tool": 30,
        "analytics_tool": 100
    },
    "other_costs": {
        "influencer_fees": 1000,
        "content_production": 500
    }
}
"""

import argparse
import json
import sys
from datetime import datetime
from typing import Optional


# ---------------------------------------------------------------------------
# Cost aggregation
# ---------------------------------------------------------------------------

def load_costs(costs_path: str, hourly_rate: float) -> dict:
    """Load and aggregate cost data."""
    with open(costs_path, "r", encoding="utf-8") as f:
        costs_raw = json.load(f)

    # Time invested
    time_hours = costs_raw.get("time_invested_hours", {})
    total_hours = sum(v for v in time_hours.values() if isinstance(v, (int, float)))
    time_cost = total_hours * hourly_rate

    # Ad spend
    ad_spend = costs_raw.get("ad_spend", {})
    total_ad_spend = sum(v for v in ad_spend.values() if isinstance(v, (int, float)))

    # Tools
    tools = costs_raw.get("tools_cost", {})
    total_tools = sum(v for v in tools.values() if isinstance(v, (int, float)))

    # Other
    other = costs_raw.get("other_costs", {})
    total_other = sum(v for v in other.values() if isinstance(v, (int, float)))

    total_cost = time_cost + total_ad_spend + total_tools + total_other

    return {
        "breakdown": {
            "time_invested": {
                "hours": total_hours,
                "hourly_rate": hourly_rate,
                "cost": round(time_cost, 2),
                "detail": time_hours,
            },
            "ad_spend": {
                "total": round(total_ad_spend, 2),
                "by_platform": ad_spend,
            },
            "tools": {
                "total": round(total_tools, 2),
                "detail": tools,
            },
            "other": {
                "total": round(total_other, 2),
                "detail": other,
            },
        },
        "total_cost": round(total_cost, 2),
    }


# ---------------------------------------------------------------------------
# ROI calculations
# ---------------------------------------------------------------------------

def safe_div(a: Optional[float], b: Optional[float]) -> Optional[float]:
    if a is None or b is None or b == 0:
        return None
    return a / b


def calculate_overall_roi(total_cost: float, revenue: float) -> dict:
    """Calculate overall content ROI."""
    net_return = revenue - total_cost
    roi_pct = safe_div(net_return, total_cost)
    roi_pct = round(roi_pct * 100, 2) if roi_pct is not None else None

    return {
        "revenue_attributable": round(revenue, 2),
        "total_cost": round(total_cost, 2),
        "net_return": round(net_return, 2),
        "roi_percentage": roi_pct,
        "roi_ratio": round(safe_div(revenue, total_cost), 2) if safe_div(revenue, total_cost) is not None else None,
    }


def calculate_engagement_costs(metrics: dict, costs: dict) -> dict:
    """Calculate cost-per-engagement and cost-per-click across platforms."""
    total_cost = costs["total_cost"]
    ad_spend_by_platform = costs["breakdown"]["ad_spend"]["by_platform"]

    results = {}
    total_engagements = 0
    total_clicks = 0

    for platform, pdata in metrics.get("platforms", {}).items():
        acct = pdata.get("account_metrics", {})
        engagement = acct.get("engagement_total")
        clicks = acct.get("website_clicks")

        platform_ad_spend = ad_spend_by_platform.get(platform, 0)

        platform_result = {
            "engagements": engagement,
            "clicks": clicks,
            "ad_spend": platform_ad_spend,
        }

        if engagement and platform_ad_spend > 0:
            platform_result["cpe_paid"] = round(platform_ad_spend / engagement, 4)
        else:
            platform_result["cpe_paid"] = None

        if clicks and platform_ad_spend > 0:
            platform_result["cpc_paid"] = round(platform_ad_spend / clicks, 4)
        else:
            platform_result["cpc_paid"] = None

        if engagement:
            total_engagements += engagement
        if clicks:
            total_clicks += clicks

        results[platform] = platform_result

    # Overall CPE and CPC (all costs, not just ad spend)
    results["overall"] = {
        "total_engagements": total_engagements,
        "total_clicks": total_clicks,
        "total_cost": total_cost,
        "cpe_all_costs": round(safe_div(total_cost, total_engagements), 4) if safe_div(total_cost, total_engagements) else None,
        "cpc_all_costs": round(safe_div(total_cost, total_clicks), 4) if safe_div(total_cost, total_clicks) else None,
    }

    return results


def calculate_roas(costs: dict, revenue: float) -> dict:
    """Calculate Return on Ad Spend."""
    total_ad_spend = costs["breakdown"]["ad_spend"]["total"]
    if total_ad_spend == 0:
        return {"status": "no_ad_spend", "roas": None}

    roas = round(revenue / total_ad_spend, 2)

    # ROAS by platform (proportional attribution)
    by_platform = {}
    for platform, spend in costs["breakdown"]["ad_spend"]["by_platform"].items():
        if spend > 0 and total_ad_spend > 0:
            platform_share = spend / total_ad_spend
            platform_revenue = revenue * platform_share  # Proportional attribution
            by_platform[platform] = {
                "ad_spend": spend,
                "attributed_revenue": round(platform_revenue, 2),
                "roas": round(safe_div(platform_revenue, spend), 2) if safe_div(platform_revenue, spend) else None,
            }

    return {
        "total_ad_spend": total_ad_spend,
        "total_revenue": revenue,
        "roas": roas,
        "roas_interpretation": _interpret_roas(roas),
        "by_platform": by_platform,
        "note": "Platform-level ROAS uses proportional attribution. For accurate per-platform ROAS, use platform-specific conversion tracking.",
    }


def _interpret_roas(roas: float) -> str:
    if roas >= 5.0:
        return "Excellent. Strong return on ad spend."
    elif roas >= 3.0:
        return "Good. Ad spend is generating healthy returns."
    elif roas >= 2.0:
        return "Average. Returns cover costs but there is room to optimize."
    elif roas >= 1.0:
        return "Below average. Revenue barely covers ad spend. Review targeting and creative."
    else:
        return "Poor. Ad spend exceeds attributed revenue. Consider pausing and restructuring campaigns."


def calculate_content_type_roi(metrics: dict, total_cost: float, revenue: float) -> dict:
    """Estimate ROI by content type based on engagement share."""
    type_engagements = {}
    total_engagement = 0

    for platform, pdata in metrics.get("platforms", {}).items():
        for post in pdata.get("posts", []):
            ctype = post.get("content_type") or "unknown"
            likes = post.get("likes") or 0
            comments = post.get("comments") or 0
            shares = post.get("shares") or 0
            saves = post.get("saves") or 0
            eng = likes + comments + shares + saves
            type_engagements.setdefault(ctype, 0)
            type_engagements[ctype] += eng
            total_engagement += eng

    if total_engagement == 0:
        return {"status": "no_engagement_data"}

    results = {}
    for ctype, eng in sorted(type_engagements.items(), key=lambda x: x[1], reverse=True):
        share = eng / total_engagement
        attributed_cost = total_cost * share
        attributed_revenue = revenue * share
        net = attributed_revenue - attributed_cost
        roi_pct = round((net / attributed_cost) * 100, 2) if attributed_cost > 0 else None

        results[ctype] = {
            "engagements": eng,
            "engagement_share": round(share * 100, 2),
            "attributed_cost": round(attributed_cost, 2),
            "attributed_revenue": round(attributed_revenue, 2),
            "net_return": round(net, 2),
            "roi_percentage": roi_pct,
        }

    results["note"] = (
        "Content type ROI uses engagement-share attribution: each content type is assigned "
        "a share of costs and revenue proportional to its share of total engagements. "
        "This is a heuristic; use conversion tracking for precise attribution."
    )

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Calculate social media ROI.")
    parser.add_argument("--metrics", type=str, required=True, help="Metrics JSON path")
    parser.add_argument("--costs", type=str, required=True, help="Costs JSON path")
    parser.add_argument("--revenue-attributable", type=float, required=True,
                        help="Revenue attributable to social media ($)")
    parser.add_argument("--hourly-rate", type=float, default=50.0,
                        help="Hourly rate for time cost (default: $50)")
    parser.add_argument("--output", type=str, default="roi.json", help="Output JSON path")
    args = parser.parse_args()

    with open(args.metrics, "r", encoding="utf-8") as f:
        metrics = json.load(f)

    costs = load_costs(args.costs, args.hourly_rate)
    revenue = args.revenue_attributable

    result = {
        "roi_metadata": {
            "calculated_at": datetime.utcnow().isoformat() + "Z",
            "account_name": metrics.get("collection_metadata", {}).get("account_name"),
            "period": metrics.get("collection_metadata", {}).get("period"),
            "hourly_rate": args.hourly_rate,
        },
        "cost_breakdown": costs,
        "overall_roi": calculate_overall_roi(costs["total_cost"], revenue),
        "engagement_costs": calculate_engagement_costs(metrics, costs),
        "roas": calculate_roas(costs, revenue),
        "roi_by_content_type": calculate_content_type_roi(metrics, costs["total_cost"], revenue),
    }

    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, default=str)

    print(f"ROI analysis written to {args.output}")
    print(f"\n--- Summary ---")
    print(f"Total Cost:     ${costs['total_cost']:,.2f}")
    print(f"Revenue:        ${revenue:,.2f}")
    overall = result["overall_roi"]
    print(f"Net Return:     ${overall['net_return']:,.2f}")
    print(f"ROI:            {overall['roi_percentage']}%")
    roas = result["roas"].get("roas")
    if roas is not None:
        print(f"ROAS:           {roas}x")
    cpe = result["engagement_costs"].get("overall", {}).get("cpe_all_costs")
    if cpe is not None:
        print(f"Cost/Engagement: ${cpe}")
    cpc = result["engagement_costs"].get("overall", {}).get("cpc_all_costs")
    if cpc is not None:
        print(f"Cost/Click:      ${cpc}")


if __name__ == "__main__":
    main()
