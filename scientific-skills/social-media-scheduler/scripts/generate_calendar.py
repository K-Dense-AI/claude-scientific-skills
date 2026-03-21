#!/usr/bin/env python3
"""
Generate a content calendar from a campaign plan.

Creates a day-by-day calendar with content slots distributed across platforms,
enforcing content mix rules, pillar rotation, and format variety.

Usage:
    python generate_calendar.py \
        --campaign-plan campaign.json \
        --start-date 2026-04-01 \
        --duration-days 30 \
        --platforms instagram,linkedin,tiktok \
        --posts-per-week '{"instagram": 5, "linkedin": 3, "tiktok": 4}' \
        --output calendar.json
"""

import argparse
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

# --- Default posting times by platform (local timezone) ---
DEFAULT_POSTING_TIMES = {
    "instagram": ["09:00", "10:00", "12:00", "17:00", "19:00"],
    "facebook": ["09:00", "10:00", "13:00", "18:00", "20:00"],
    "linkedin": ["08:00", "09:00", "12:00", "17:00"],
    "tiktok": ["10:00", "11:00", "14:00", "19:00", "21:00"],
}

# --- Default format pools by platform ---
DEFAULT_FORMATS = {
    "instagram": ["reel", "carousel", "static", "story"],
    "facebook": ["video", "image", "link", "carousel"],
    "linkedin": ["text", "carousel", "video", "document"],
    "tiktok": ["short-video", "duet", "stitch"],
}

# --- Default content pillars and their target distribution ---
DEFAULT_PILLARS = [
    {"name": "education", "weight": 0.30},
    {"name": "social-proof", "weight": 0.20},
    {"name": "behind-the-scenes", "weight": 0.20},
    {"name": "promotion", "weight": 0.15},
    {"name": "entertainment", "weight": 0.10},
    {"name": "community", "weight": 0.05},
]

# --- Default CTAs ---
DEFAULT_CTAS = ["save", "share", "comment", "link-in-bio", "DM", "shop-now"]

# --- Days ranked by engagement per platform ---
BEST_DAYS = {
    "instagram": ["Tuesday", "Wednesday", "Monday", "Thursday", "Friday", "Saturday", "Sunday"],
    "facebook": ["Wednesday", "Thursday", "Tuesday", "Friday", "Monday", "Saturday", "Sunday"],
    "linkedin": ["Tuesday", "Wednesday", "Thursday", "Monday", "Friday", "Saturday", "Sunday"],
    "tiktok": ["Tuesday", "Thursday", "Wednesday", "Friday", "Saturday", "Monday", "Sunday"],
}

DAY_NAMES = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]


def load_campaign_plan(path: str) -> dict:
    """Load campaign plan from JSON file."""
    with open(path, "r") as f:
        return json.load(f)


def build_phase_map(phases: list, start_date: datetime, duration_days: int) -> list:
    """Map each day in the range to a campaign phase."""
    phase_map = []
    if not phases:
        # No phases defined: all days are 'general'
        for i in range(duration_days):
            phase_map.append("general")
        return phase_map

    # Calculate total phase days
    total_phase_days = sum(p.get("days", 0) for p in phases)
    if total_phase_days < duration_days:
        # Extend the last phase to fill remaining days
        phases[-1]["days"] = phases[-1].get("days", 0) + (duration_days - total_phase_days)

    day_idx = 0
    for phase in phases:
        for _ in range(phase.get("days", 0)):
            if day_idx < duration_days:
                phase_map.append(phase["name"])
                day_idx += 1
    # Fill any remaining days
    while len(phase_map) < duration_days:
        phase_map.append(phases[-1]["name"] if phases else "general")
    return phase_map


def distribute_posts(
    start_date: datetime,
    duration_days: int,
    platform: str,
    posts_per_week: int,
    phase_map: list,
) -> list:
    """Distribute posting slots across the date range for a single platform."""
    slots = []
    best_days = BEST_DAYS.get(platform, DAY_NAMES)

    # Process week by week
    week_start = 0
    while week_start < duration_days:
        week_end = min(week_start + 7, duration_days)
        week_days = []

        for d in range(week_start, week_end):
            current_date = start_date + timedelta(days=d)
            day_name = DAY_NAMES[current_date.weekday()]
            week_days.append((d, current_date, day_name))

        # Sort days in this week by platform preference
        week_days.sort(key=lambda x: best_days.index(x[2]) if x[2] in best_days else 99)

        # Pick the top N days for posting
        posts_this_week = min(posts_per_week, len(week_days))
        selected = week_days[:posts_this_week]
        selected.sort(key=lambda x: x[0])  # Re-sort by date order

        for day_idx, current_date, day_name in selected:
            slots.append({
                "day_index": day_idx,
                "date": current_date.strftime("%Y-%m-%d"),
                "day_of_week": day_name,
                "phase": phase_map[day_idx] if day_idx < len(phase_map) else "general",
            })

        week_start += 7

    return slots


def assign_content(
    slots: list,
    platform: str,
    campaign_plan: dict,
) -> list:
    """Assign content pillars, formats, CTAs, and posting times to slots."""
    formats = campaign_plan.get("formats", {}).get(platform, DEFAULT_FORMATS.get(platform, ["post"]))
    pillars = campaign_plan.get("pillars", [p["name"] for p in DEFAULT_PILLARS])
    posting_times = DEFAULT_POSTING_TIMES.get(platform, ["10:00"])
    ctas = DEFAULT_CTAS.copy()

    entries = []
    last_pillar = None
    format_idx = 0
    time_idx = 0
    cta_idx = 0

    for slot in slots:
        # Rotate pillar (avoid consecutive repeats)
        pillar_pool = [p for p in pillars if p != last_pillar] if last_pillar else pillars
        if not pillar_pool:
            pillar_pool = pillars
        pillar = pillar_pool[len(entries) % len(pillar_pool)]
        last_pillar = pillar

        # Rotate format
        fmt = formats[format_idx % len(formats)]
        format_idx += 1

        # Rotate posting time
        posting_time = posting_times[time_idx % len(posting_times)]
        time_idx += 1

        # Rotate CTA
        cta = ctas[cta_idx % len(ctas)]
        cta_idx += 1

        # Determine visual type based on format
        visual_map = {
            "reel": "talking-head",
            "carousel": "graphic",
            "static": "photo",
            "story": "graphic",
            "video": "animation",
            "image": "photo",
            "link": "graphic",
            "text": "graphic",
            "document": "graphic",
            "short-video": "talking-head",
            "duet": "UGC",
            "stitch": "UGC",
        }
        visual_type = visual_map.get(fmt, "graphic")

        entry = {
            "date": slot["date"],
            "day_of_week": slot["day_of_week"],
            "platform": platform,
            "format": fmt,
            "content_pillar": pillar,
            "hook_placeholder": f"[{pillar} hook for {platform} {fmt}]",
            "visual_type": visual_type,
            "cta": cta,
            "posting_time": posting_time,
            "status": "draft",
            "phase": slot["phase"],
        }
        entries.append(entry)

    return entries


def calendar_to_markdown(entries: list) -> str:
    """Convert calendar entries to a Markdown table."""
    lines = []
    lines.append("| Date | Day | Platform | Format | Pillar | Hook | Visual | CTA | Time | Phase | Status |")
    lines.append("|------|-----|----------|--------|--------|------|--------|-----|------|-------|--------|")

    for e in sorted(entries, key=lambda x: (x["date"], x["platform"])):
        lines.append(
            f"| {e['date']} | {e['day_of_week'][:3]} | {e['platform']} | {e['format']} "
            f"| {e['content_pillar']} | {e['hook_placeholder']} | {e['visual_type']} "
            f"| {e['cta']} | {e['posting_time']} | {e['phase']} | {e['status']} |"
        )

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Generate a content calendar from a campaign plan.")
    parser.add_argument("--campaign-plan", required=True, help="Path to campaign plan JSON file.")
    parser.add_argument("--start-date", required=True, help="Start date in YYYY-MM-DD format.")
    parser.add_argument("--duration-days", type=int, required=True, help="Number of days for the calendar.")
    parser.add_argument("--platforms", required=True, help="Comma-separated list of platforms.")
    parser.add_argument(
        "--posts-per-week",
        required=True,
        help='JSON object mapping platform to posts per week, e.g. \'{"instagram": 5, "linkedin": 3}\'.',
    )
    parser.add_argument("--output", required=True, help="Output file path (JSON). Markdown also generated with .md extension.")

    args = parser.parse_args()

    # Parse inputs
    campaign_plan = load_campaign_plan(args.campaign_plan)
    start_date = datetime.strptime(args.start_date, "%Y-%m-%d")
    platforms = [p.strip().lower() for p in args.platforms.split(",")]
    posts_per_week = json.loads(args.posts_per_week)

    # Build phase map
    phases = campaign_plan.get("phases", [])
    phase_map = build_phase_map(phases, start_date, args.duration_days)

    # Generate calendar entries for each platform
    all_entries = []
    for platform in platforms:
        ppw = posts_per_week.get(platform, 3)
        slots = distribute_posts(start_date, args.duration_days, platform, ppw, phase_map)
        entries = assign_content(slots, platform, campaign_plan)
        all_entries.extend(entries)

    # Sort by date then platform
    all_entries.sort(key=lambda x: (x["date"], x["platform"]))

    # Build output
    output = {
        "campaign_name": campaign_plan.get("campaign_name", "Untitled Campaign"),
        "brand": campaign_plan.get("brand", ""),
        "start_date": args.start_date,
        "duration_days": args.duration_days,
        "platforms": platforms,
        "total_posts": len(all_entries),
        "entries": all_entries,
    }

    # Write JSON
    output_path = Path(args.output)
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"JSON calendar written to {output_path} ({len(all_entries)} posts)")

    # Write Markdown
    md_path = output_path.with_suffix(".md")
    md_content = f"# Content Calendar: {output['campaign_name']}\n\n"
    md_content += f"**Brand**: {output['brand']}  \n"
    md_content += f"**Period**: {args.start_date} to {(start_date + timedelta(days=args.duration_days - 1)).strftime('%Y-%m-%d')}  \n"
    md_content += f"**Total Posts**: {len(all_entries)}  \n\n"
    md_content += calendar_to_markdown(all_entries)
    md_content += "\n"

    with open(md_path, "w") as f:
        f.write(md_content)
    print(f"Markdown calendar written to {md_path}")


if __name__ == "__main__":
    main()
