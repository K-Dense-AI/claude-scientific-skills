#!/usr/bin/env python3
"""
Social Media Campaign Orchestrator

Master orchestration script that generates a complete campaign plan
by coordinating all sub-skills through a structured 7-phase pipeline.

Outputs a campaign plan JSON and creates the full output directory structure
for content production, publishing, and analytics.

Usage:
    python orchestrate_campaign.py \
        --brand-url "https://example.com" \
        --campaign-type brand_awareness \
        --platforms instagram,linkedin \
        --duration-days 30 \
        --goals "followers:+2000,engagement_rate:4%" \
        --output-dir ./campaign_output
"""

import argparse
import json
import os
import sys
from datetime import datetime, timedelta
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CAMPAIGN_TYPES = [
    "product_launch",
    "brand_awareness",
    "engagement_growth",
    "lead_generation",
    "event_promotion",
    "seasonal",
]

SUPPORTED_PLATFORMS = [
    "instagram",
    "linkedin",
    "tiktok",
    "facebook",
    "twitter",
]

CAMPAIGN_STATES = ["init", "plan", "create", "review", "schedule", "monitor", "complete"]

CONTENT_MIX = {
    "value": 0.60,
    "curated": 0.30,
    "promotional": 0.10,
}

DEFAULT_POSTING_CADENCE = {
    "instagram": {"feed": 5, "stories": 7, "reels": 4},  # per week
    "linkedin": {"posts": 4, "articles": 1},
    "tiktok": {"videos": 7},
    "facebook": {"posts": 4, "stories": 5},
    "twitter": {"tweets": 14, "threads": 2},
}

CAMPAIGN_DEFAULTS = {
    "product_launch": {
        "phases": [
            {"name": "tease", "duration_pct": 0.25, "content_focus": "curiosity"},
            {"name": "announce", "duration_pct": 0.35, "content_focus": "excitement"},
            {"name": "sustain", "duration_pct": 0.40, "content_focus": "value"},
        ],
        "recommended_platforms": ["instagram", "tiktok", "linkedin"],
        "content_volume_multiplier": 1.5,
    },
    "brand_awareness": {
        "phases": [
            {"name": "foundation", "duration_pct": 0.25, "content_focus": "brand_story"},
            {"name": "education", "duration_pct": 0.30, "content_focus": "value"},
            {"name": "entertainment", "duration_pct": 0.25, "content_focus": "engagement"},
            {"name": "community", "duration_pct": 0.20, "content_focus": "ugc"},
        ],
        "recommended_platforms": ["instagram", "tiktok", "facebook"],
        "content_volume_multiplier": 1.0,
    },
    "engagement_growth": {
        "phases": [
            {"name": "activate", "duration_pct": 0.30, "content_focus": "interactive"},
            {"name": "build", "duration_pct": 0.40, "content_focus": "community"},
            {"name": "celebrate", "duration_pct": 0.30, "content_focus": "recognition"},
        ],
        "recommended_platforms": ["instagram", "tiktok", "facebook"],
        "content_volume_multiplier": 1.2,
    },
    "lead_generation": {
        "phases": [
            {"name": "attract", "duration_pct": 0.35, "content_focus": "tofu"},
            {"name": "nurture", "duration_pct": 0.40, "content_focus": "mofu"},
            {"name": "convert", "duration_pct": 0.25, "content_focus": "bofu"},
        ],
        "recommended_platforms": ["linkedin", "facebook", "instagram"],
        "content_volume_multiplier": 1.0,
    },
    "event_promotion": {
        "phases": [
            {"name": "before", "duration_pct": 0.50, "content_focus": "anticipation"},
            {"name": "during", "duration_pct": 0.15, "content_focus": "live"},
            {"name": "after", "duration_pct": 0.35, "content_focus": "recap"},
        ],
        "recommended_platforms": ["instagram", "linkedin", "twitter"],
        "content_volume_multiplier": 1.3,
    },
    "seasonal": {
        "phases": [
            {"name": "tease", "duration_pct": 0.25, "content_focus": "anticipation"},
            {"name": "build", "duration_pct": 0.30, "content_focus": "engagement"},
            {"name": "execute", "duration_pct": 0.30, "content_focus": "conversion"},
            {"name": "wrap_up", "duration_pct": 0.15, "content_focus": "gratitude"},
        ],
        "recommended_platforms": ["instagram", "tiktok", "facebook"],
        "content_volume_multiplier": 1.4,
    },
}


# ---------------------------------------------------------------------------
# Core Functions
# ---------------------------------------------------------------------------

def parse_goals(goals_str: str) -> dict:
    """Parse goal string like 'followers:+2000,engagement_rate:4%' into a dict."""
    goals = {}
    if not goals_str:
        return goals
    for pair in goals_str.split(","):
        pair = pair.strip()
        if ":" not in pair:
            continue
        key, value = pair.split(":", 1)
        goals[key.strip()] = value.strip()
    return goals


def calculate_content_slots(
    platforms: list[str],
    duration_days: int,
    multiplier: float,
) -> dict:
    """Calculate total content slots per platform for the campaign duration."""
    weeks = duration_days / 7
    slots = {}
    for platform in platforms:
        cadence = DEFAULT_POSTING_CADENCE.get(platform, {})
        platform_slots = {}
        for fmt, per_week in cadence.items():
            count = int(per_week * weeks * multiplier)
            platform_slots[fmt] = count
        slots[platform] = platform_slots
    return slots


def apply_content_mix(total_posts: int) -> dict:
    """Distribute posts according to the 60/30/10 content mix."""
    return {
        "value": int(total_posts * CONTENT_MIX["value"]),
        "curated": int(total_posts * CONTENT_MIX["curated"]),
        "promotional": total_posts
        - int(total_posts * CONTENT_MIX["value"])
        - int(total_posts * CONTENT_MIX["curated"]),
    }


def generate_phase_timeline(
    campaign_type: str,
    start_date: str,
    duration_days: int,
) -> list[dict]:
    """Generate phase timeline with start/end dates."""
    config = CAMPAIGN_DEFAULTS[campaign_type]
    start = datetime.strptime(start_date, "%Y-%m-%d")
    phases = []
    current_start = start

    for phase_def in config["phases"]:
        phase_days = max(1, int(duration_days * phase_def["duration_pct"]))
        phase_end = current_start + timedelta(days=phase_days - 1)
        phases.append({
            "name": phase_def["name"],
            "content_focus": phase_def["content_focus"],
            "start_date": current_start.strftime("%Y-%m-%d"),
            "end_date": phase_end.strftime("%Y-%m-%d"),
            "duration_days": phase_days,
        })
        current_start = phase_end + timedelta(days=1)

    return phases


def generate_skill_invocations(campaign_type: str, platforms: list[str]) -> list[dict]:
    """Generate the ordered list of skill invocations for the pipeline."""
    return [
        {
            "phase": 1,
            "phase_name": "Brand Discovery",
            "skill": "brand-identity-miner",
            "description": "Extract brand identity from web presence",
            "output": "brand_dna.json",
            "required_input": "brand_url",
        },
        {
            "phase": 2,
            "phase_name": "Strategy Planning",
            "skill": "orchestrator (internal)",
            "description": f"Generate {campaign_type} campaign plan for {', '.join(platforms)}",
            "output": "campaign_plan.json",
            "required_input": "brand_dna.json",
        },
        {
            "phase": 3,
            "phase_name": "Content Planning",
            "skill": "viral-content-creator + social-media-scheduler",
            "description": "Generate content calendar and briefs",
            "output": "content_calendar.json + content_briefs/*.json",
            "required_input": "campaign_plan.json + brand_dna.json",
        },
        {
            "phase": 4,
            "phase_name": "Content Creation",
            "skill": "social-media-copywriting + social-media-visual-generator",
            "description": "Produce copy and visual assets",
            "output": "captions/*.json + visuals/*.png + videos/*.mp4",
            "required_input": "content_briefs/*.json + brand_dna.json",
        },
        {
            "phase": 5,
            "phase_name": "Review & Approval",
            "skill": "orchestrator (user interaction)",
            "description": "Present content for approval with brand consistency check",
            "output": "approved content batch",
            "required_input": "all content assets",
        },
        {
            "phase": 6,
            "phase_name": "Publishing",
            "skill": "social-media-scheduler",
            "description": f"Publish to {', '.join(platforms)} via Composio integrations",
            "output": "publish_log.json",
            "required_input": "approved content + content_calendar.json",
        },
        {
            "phase": 7,
            "phase_name": "Analytics & Optimization",
            "skill": "social-media-analytics",
            "description": "Track performance, generate report, feed back insights",
            "output": "analytics_report.json + optimization_recommendations.json",
            "required_input": "publish_log.json + platform API access",
        },
    ]


def create_output_structure(output_dir: str) -> dict:
    """Create the campaign output directory structure."""
    directories = [
        output_dir,
        os.path.join(output_dir, "content_briefs"),
        os.path.join(output_dir, "captions"),
        os.path.join(output_dir, "visuals"),
        os.path.join(output_dir, "videos"),
        os.path.join(output_dir, "analytics"),
    ]
    created = []
    for d in directories:
        os.makedirs(d, exist_ok=True)
        created.append(d)
    return {"directories_created": created}


def generate_campaign_plan(
    brand_url: str,
    campaign_type: str,
    platforms: list[str],
    duration_days: int,
    goals: dict,
    output_dir: str,
) -> dict:
    """Generate the complete campaign plan JSON."""
    config = CAMPAIGN_DEFAULTS[campaign_type]
    start_date = datetime.now().strftime("%Y-%m-%d")
    end_date = (datetime.now() + timedelta(days=duration_days)).strftime("%Y-%m-%d")

    content_slots = calculate_content_slots(
        platforms, duration_days, config["content_volume_multiplier"]
    )
    total_posts = sum(
        sum(formats.values()) for formats in content_slots.values()
    )
    content_mix = apply_content_mix(total_posts)
    phase_timeline = generate_phase_timeline(campaign_type, start_date, duration_days)
    skill_invocations = generate_skill_invocations(campaign_type, platforms)

    plan = {
        "campaign_id": f"{campaign_type}_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        "brand_url": brand_url,
        "campaign_type": campaign_type,
        "state": "init",
        "created_at": datetime.now().isoformat(),
        "timeline": {
            "start_date": start_date,
            "end_date": end_date,
            "duration_days": duration_days,
        },
        "platforms": platforms,
        "goals": goals,
        "phases": phase_timeline,
        "content_plan": {
            "total_content_pieces": total_posts,
            "content_mix": content_mix,
            "slots_per_platform": content_slots,
        },
        "skill_invocations": skill_invocations,
        "output_directory": os.path.abspath(output_dir),
        "kpi_targets": {
            "engagement_rate": goals.get("engagement_rate", "3%"),
            "follower_growth": goals.get("followers", "+500"),
            "reach_multiplier": "2x baseline",
            "content_completion_rate": "100%",
        },
        "analytics_checkpoints": [
            {
                "name": "7-day check",
                "day": 7,
                "date": (datetime.now() + timedelta(days=7)).strftime("%Y-%m-%d"),
                "action": "early_performance_review",
            },
            {
                "name": "14-day review",
                "day": 14,
                "date": (datetime.now() + timedelta(days=14)).strftime("%Y-%m-%d"),
                "action": "mid_campaign_optimization",
            },
            {
                "name": "30-day report",
                "day": min(30, duration_days),
                "date": (
                    datetime.now() + timedelta(days=min(30, duration_days))
                ).strftime("%Y-%m-%d"),
                "action": "final_campaign_report",
            },
        ],
    }

    return plan


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Social Media Campaign Orchestrator — generate a complete campaign plan."
    )
    parser.add_argument(
        "--brand-url",
        required=True,
        help="Brand website URL (e.g., https://example.com)",
    )
    parser.add_argument(
        "--campaign-type",
        required=True,
        choices=CAMPAIGN_TYPES,
        help="Type of campaign to run",
    )
    parser.add_argument(
        "--platforms",
        required=True,
        help="Comma-separated list of target platforms (e.g., instagram,linkedin)",
    )
    parser.add_argument(
        "--duration-days",
        type=int,
        default=30,
        help="Campaign duration in days (default: 30)",
    )
    parser.add_argument(
        "--goals",
        default="",
        help='Campaign goals as key:value pairs (e.g., "followers:+2000,engagement_rate:4%%")',
    )
    parser.add_argument(
        "--output-dir",
        default="./campaign_output",
        help="Output directory for campaign files (default: ./campaign_output)",
    )
    args = parser.parse_args()

    # Validate platforms
    platforms = [p.strip().lower() for p in args.platforms.split(",")]
    invalid = [p for p in platforms if p not in SUPPORTED_PLATFORMS]
    if invalid:
        print(f"Error: Unsupported platforms: {', '.join(invalid)}", file=sys.stderr)
        print(f"Supported platforms: {', '.join(SUPPORTED_PLATFORMS)}", file=sys.stderr)
        sys.exit(1)

    # Parse goals
    goals = parse_goals(args.goals)

    # Create output structure
    print(f"Creating output directory structure at: {args.output_dir}")
    structure = create_output_structure(args.output_dir)
    for d in structure["directories_created"]:
        print(f"  Created: {d}")

    # Generate campaign plan
    print(f"\nGenerating {args.campaign_type} campaign plan...")
    print(f"  Brand: {args.brand_url}")
    print(f"  Platforms: {', '.join(platforms)}")
    print(f"  Duration: {args.duration_days} days")
    if goals:
        print(f"  Goals: {goals}")

    plan = generate_campaign_plan(
        brand_url=args.brand_url,
        campaign_type=args.campaign_type,
        platforms=platforms,
        duration_days=args.duration_days,
        goals=goals,
        output_dir=args.output_dir,
    )

    # Write campaign plan
    plan_path = os.path.join(args.output_dir, "campaign_plan.json")
    with open(plan_path, "w") as f:
        json.dump(plan, f, indent=2)
    print(f"\nCampaign plan written to: {plan_path}")

    # Summary
    print("\n" + "=" * 60)
    print("CAMPAIGN PLAN SUMMARY")
    print("=" * 60)
    print(f"Campaign ID:    {plan['campaign_id']}")
    print(f"Type:           {plan['campaign_type']}")
    print(f"Duration:       {plan['timeline']['duration_days']} days")
    print(f"Platforms:      {', '.join(plan['platforms'])}")
    print(f"Total content:  {plan['content_plan']['total_content_pieces']} pieces")
    print(f"  Value:        {plan['content_plan']['content_mix']['value']}")
    print(f"  Curated:      {plan['content_plan']['content_mix']['curated']}")
    print(f"  Promotional:  {plan['content_plan']['content_mix']['promotional']}")
    print(f"\nPhases:")
    for phase in plan["phases"]:
        print(f"  {phase['name']:15s} {phase['start_date']} to {phase['end_date']} ({phase['duration_days']}d)")
    print(f"\nSkill Pipeline:")
    for inv in plan["skill_invocations"]:
        print(f"  Phase {inv['phase']}: {inv['phase_name']:25s} -> {inv['skill']}")
    print(f"\nNext step: Run brand-identity-miner on {args.brand_url}")
    print(f"           to generate brand_dna.json in {args.output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
