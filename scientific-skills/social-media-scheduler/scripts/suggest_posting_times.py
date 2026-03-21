#!/usr/bin/env python3
"""
Suggest optimal posting times for a given platform and industry.

Returns ranked posting windows based on research data from Sprout Social,
Hootsuite, Buffer, and Later. Adjusts for the user's timezone.

Usage:
    python suggest_posting_times.py \
        --platform instagram \
        --industry b2c \
        --timezone "America/New_York" \
        --output posting_windows.json
"""

import argparse
import json
import sys
from datetime import datetime

# --- Research-backed posting time data ---
# Times are in a generic "local timezone" basis. The script maps them
# to the user's specified timezone label for clarity.

# Structure: { platform: { industry: { day: [(hour, minute, score), ...] } } }
# Score 1-10 where 10 = highest engagement probability.

GENERAL_TIMES = {
    "instagram": {
        "general": {
            "Monday":    [(9, 0, 8), (10, 0, 9), (12, 0, 7), (17, 0, 7), (19, 0, 6)],
            "Tuesday":   [(9, 0, 10), (10, 0, 9), (12, 0, 8), (17, 0, 7), (19, 0, 7)],
            "Wednesday": [(9, 0, 9), (10, 0, 9), (12, 0, 8), (17, 0, 7), (19, 0, 7)],
            "Thursday":  [(9, 0, 8), (10, 0, 8), (12, 0, 7), (17, 0, 7), (19, 0, 6)],
            "Friday":    [(9, 0, 7), (10, 0, 8), (12, 0, 7), (17, 0, 5), (19, 0, 5)],
            "Saturday":  [(9, 0, 5), (10, 0, 6), (11, 0, 6), (19, 0, 5)],
            "Sunday":    [(10, 0, 4), (11, 0, 5), (19, 0, 4)],
        },
    },
    "facebook": {
        "general": {
            "Monday":    [(9, 0, 7), (10, 0, 7), (13, 0, 6), (18, 0, 5)],
            "Tuesday":   [(9, 0, 8), (10, 0, 8), (13, 0, 7), (18, 0, 6)],
            "Wednesday": [(9, 0, 10), (10, 0, 9), (13, 0, 8), (18, 0, 7), (20, 0, 6)],
            "Thursday":  [(9, 0, 9), (10, 0, 8), (13, 0, 8), (18, 0, 7)],
            "Friday":    [(9, 0, 7), (10, 0, 7), (13, 0, 6), (18, 0, 5)],
            "Saturday":  [(10, 0, 5), (11, 0, 5)],
            "Sunday":    [(10, 0, 4), (11, 0, 4)],
        },
    },
    "linkedin": {
        "general": {
            "Monday":    [(8, 0, 7), (9, 0, 8), (12, 0, 7), (17, 0, 6)],
            "Tuesday":   [(8, 0, 9), (9, 0, 10), (12, 0, 8), (17, 0, 7)],
            "Wednesday": [(8, 0, 9), (9, 0, 9), (12, 0, 8), (17, 0, 7)],
            "Thursday":  [(8, 0, 8), (9, 0, 9), (12, 0, 7), (17, 0, 6)],
            "Friday":    [(8, 0, 6), (9, 0, 7), (12, 0, 6)],
            "Saturday":  [(9, 0, 3)],
            "Sunday":    [(9, 0, 2)],
        },
    },
    "tiktok": {
        "general": {
            "Monday":    [(10, 0, 6), (11, 0, 7), (14, 0, 6), (19, 0, 7)],
            "Tuesday":   [(10, 0, 9), (11, 0, 10), (14, 0, 7), (19, 0, 8), (21, 0, 7)],
            "Wednesday": [(10, 0, 8), (11, 0, 8), (14, 0, 7), (19, 0, 7)],
            "Thursday":  [(10, 0, 9), (11, 0, 9), (14, 0, 7), (19, 0, 8), (21, 0, 7)],
            "Friday":    [(10, 0, 7), (14, 0, 7), (19, 0, 8), (21, 0, 7)],
            "Saturday":  [(10, 0, 6), (14, 0, 6), (19, 0, 7), (21, 0, 7)],
            "Sunday":    [(10, 0, 5), (14, 0, 5), (19, 0, 6)],
        },
    },
}

# Industry-specific adjustments (additive score modifiers)
INDUSTRY_ADJUSTMENTS = {
    "b2b": {
        "instagram": {"weekday_morning": +1, "weekend": -2, "evening": -1},
        "facebook": {"weekday_morning": +1, "weekend": -2},
        "linkedin": {"weekday_morning": +2, "lunch": +1, "weekend": -3},
        "tiktok": {"weekday_morning": +1, "weekend": -1},
    },
    "b2c": {
        "instagram": {"evening": +2, "weekend": +1, "weekday_morning": 0},
        "facebook": {"afternoon": +1, "evening": +1, "weekend": +1},
        "linkedin": {"weekday_morning": 0, "weekend": -2},
        "tiktok": {"evening": +2, "weekend": +2},
    },
    "ecommerce": {
        "instagram": {"weekday_morning": +1, "evening": +2, "weekend": +1},
        "facebook": {"weekday_morning": +1, "afternoon": +1},
        "linkedin": {"weekday_morning": 0},
        "tiktok": {"evening": +2, "weekend": +1},
    },
    "saas": {
        "instagram": {"weekday_morning": +1, "weekend": -1},
        "facebook": {"weekday_morning": +1, "weekend": -1},
        "linkedin": {"weekday_morning": +2, "lunch": +1, "weekend": -3},
        "tiktok": {"weekday_morning": +1, "afternoon": +1},
    },
    "fitness": {
        "instagram": {"early_morning": +3, "evening": +2, "weekend": +1},
        "facebook": {"weekday_morning": +1},
        "linkedin": {"weekday_morning": 0},
        "tiktok": {"early_morning": +3, "evening": +2},
    },
    "food": {
        "instagram": {"lunch": +3, "evening": +2, "weekend": +1},
        "facebook": {"lunch": +2},
        "linkedin": {"lunch": +1},
        "tiktok": {"lunch": +3, "evening": +2, "weekend": +1},
    },
    "fashion": {
        "instagram": {"weekday_morning": +1, "evening": +2, "weekend": +2},
        "facebook": {"afternoon": +1},
        "linkedin": {"weekday_morning": 0},
        "tiktok": {"weekday_morning": +1, "evening": +2, "weekend": +2},
    },
}

DAY_ORDER = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]


def classify_time_slot(hour: int, day: str) -> list:
    """Classify a time slot into categories for industry adjustment."""
    categories = []
    is_weekend = day in ("Saturday", "Sunday")

    if is_weekend:
        categories.append("weekend")
    if 5 <= hour < 7:
        categories.append("early_morning")
    if 7 <= hour < 12:
        categories.append("weekday_morning")
    if 11 <= hour < 14:
        categories.append("lunch")
    if 13 <= hour < 17:
        categories.append("afternoon")
    if 17 <= hour < 22:
        categories.append("evening")

    return categories


def get_adjusted_times(platform: str, industry: str) -> list:
    """Get posting times with industry adjustments applied."""
    base_times = GENERAL_TIMES.get(platform, {}).get("general", {})
    adjustments = INDUSTRY_ADJUSTMENTS.get(industry, {}).get(platform, {})

    results = []
    for day in DAY_ORDER:
        day_slots = base_times.get(day, [])
        for hour, minute, base_score in day_slots:
            # Apply industry adjustments
            categories = classify_time_slot(hour, day)
            adj = 0
            for cat in categories:
                adj += adjustments.get(cat, 0)

            adjusted_score = max(1, min(10, base_score + adj))
            results.append({
                "day": day,
                "time": f"{hour:02d}:{minute:02d}",
                "hour": hour,
                "score": adjusted_score,
                "base_score": base_score,
                "adjustment": adj,
            })

    # Sort by score descending, then by day order
    results.sort(key=lambda x: (-x["score"], DAY_ORDER.index(x["day"]), x["hour"]))
    return results


def format_output(platform: str, industry: str, timezone: str, ranked_times: list) -> dict:
    """Format the output structure."""
    # Top 5 windows
    top_windows = []
    for i, slot in enumerate(ranked_times[:5]):
        top_windows.append({
            "rank": i + 1,
            "day": slot["day"],
            "time": slot["time"],
            "score": slot["score"],
            "timezone": timezone,
        })

    # Best days (aggregate scores)
    day_scores = {}
    for slot in ranked_times:
        day_scores.setdefault(slot["day"], 0)
        day_scores[slot["day"]] += slot["score"]
    best_days = sorted(day_scores.items(), key=lambda x: -x[1])

    # Worst times (bottom 3)
    worst_times = []
    for slot in ranked_times[-3:]:
        worst_times.append({
            "day": slot["day"],
            "time": slot["time"],
            "score": slot["score"],
            "timezone": timezone,
        })

    # Frequency recommendations
    freq = {
        "instagram": {"min_per_week": 3, "optimal_per_week": "5-7", "max_per_day": 2},
        "facebook": {"min_per_week": 3, "optimal_per_week": "5-7", "max_per_day": 2},
        "linkedin": {"min_per_week": 2, "optimal_per_week": "3-5", "max_per_day": 1},
        "tiktok": {"min_per_week": 3, "optimal_per_week": "7-14", "max_per_day": 4},
    }

    return {
        "platform": platform,
        "industry": industry,
        "timezone": timezone,
        "generated_at": datetime.now().isoformat(),
        "top_posting_windows": top_windows,
        "best_days_ranked": [{"rank": i + 1, "day": d, "aggregate_score": s} for i, (d, s) in enumerate(best_days)],
        "worst_times": worst_times,
        "frequency": freq.get(platform, {}),
        "all_ranked_slots": ranked_times,
        "data_sources": [
            "Sprout Social (2024-2025)",
            "Hootsuite (2024-2025)",
            "Buffer (2024)",
            "Later (2024)",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description="Suggest optimal posting times for a platform and industry.")
    parser.add_argument("--platform", required=True, choices=["instagram", "facebook", "linkedin", "tiktok"],
                        help="Target platform.")
    parser.add_argument("--industry", default="general",
                        choices=["general", "b2b", "b2c", "ecommerce", "saas", "fitness", "food", "fashion"],
                        help="Industry vertical for adjusted recommendations.")
    parser.add_argument("--timezone", default="UTC", help="User's timezone (e.g., 'America/New_York'). Times are labeled with this timezone.")
    parser.add_argument("--output", required=True, help="Output file path (JSON).")

    args = parser.parse_args()

    ranked_times = get_adjusted_times(args.platform, args.industry)
    output = format_output(args.platform, args.industry, args.timezone, ranked_times)

    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)

    print(f"Posting time recommendations written to {args.output}")
    print(f"\nTop 5 posting windows for {args.platform} ({args.industry}):")
    for w in output["top_posting_windows"]:
        print(f"  #{w['rank']}: {w['day']} at {w['time']} (score: {w['score']}/10) [{args.timezone}]")

    print(f"\nBest days: {', '.join(d['day'] for d in output['best_days_ranked'][:3])}")
    print(f"Frequency: {output['frequency'].get('optimal_per_week', 'N/A')} posts/week recommended")


if __name__ == "__main__":
    main()
