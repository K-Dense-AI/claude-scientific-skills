#!/usr/bin/env python3
"""
Social Media Metrics Collector

Collects performance metrics from Facebook, Instagram, LinkedIn, and TikTok.
Supports two collection methods:
  1. Composio MCP integration (automated API access)
  2. Chrome automation guidance (manual extraction with structured output)

Usage:
    python collect_metrics.py --platforms instagram,linkedin --period 30d --output metrics.json
    python collect_metrics.py --platforms all --period 7d --output weekly_metrics.json
    python collect_metrics.py --input manual_data.csv --platforms instagram --output metrics.json

Arguments:
    --platforms     Comma-separated list of platforms: facebook, instagram, linkedin, tiktok, all
    --period        Time period: 7d, 14d, 30d, 60d, 90d (default: 30d)
    --output        Output JSON file path (default: metrics.json)
    --input         Optional CSV with manually collected data to normalize into standard format
    --account-name  Brand/account name for labeling (default: "brand")
"""

import argparse
import json
import csv
import sys
from datetime import datetime, timedelta
from pathlib import Path


# ---------------------------------------------------------------------------
# Schema definitions
# ---------------------------------------------------------------------------

PLATFORMS = ["facebook", "instagram", "linkedin", "tiktok"]

PERIOD_DAYS = {
    "7d": 7,
    "14d": 14,
    "30d": 30,
    "60d": 60,
    "90d": 90,
}

ACCOUNT_LEVEL_METRICS = [
    "followers",
    "followers_gained",
    "followers_lost",
    "reach",
    "impressions",
    "profile_visits",
    "website_clicks",
    "engagement_total",
    "likes_total",
    "comments_total",
    "shares_total",
    "saves_total",
    "video_views_total",
    "story_views_total",
]

POST_LEVEL_METRICS = [
    "post_id",
    "published_at",
    "content_type",       # image, carousel, video, reel, story, text, link, document
    "caption_length",
    "hashtag_count",
    "reach",
    "impressions",
    "likes",
    "comments",
    "shares",
    "saves",
    "clicks",
    "video_views",
    "video_avg_watch_seconds",
    "video_completion_rate",
]


def empty_account_metrics(platform: str, period: str, account_name: str) -> dict:
    """Return a blank account-level metrics template."""
    return {
        "platform": platform,
        "account_name": account_name,
        "period": period,
        "collected_at": datetime.utcnow().isoformat() + "Z",
        "account_metrics": {m: None for m in ACCOUNT_LEVEL_METRICS},
        "posts": [],
    }


def empty_post_metrics() -> dict:
    """Return a blank post-level metrics template."""
    return {m: None for m in POST_LEVEL_METRICS}


# ---------------------------------------------------------------------------
# Composio / API collection helpers
# ---------------------------------------------------------------------------

COMPOSIO_INSTRUCTIONS = {
    "instagram": (
        "To collect Instagram metrics via Composio MCP:\n"
        "1. Ensure the Instagram Business/Creator connector is active in Composio.\n"
        "2. Use the `instagram_get_insights` action with fields: impressions, reach, "
        "follower_count, profile_views, website_clicks.\n"
        "3. Use `instagram_get_media` to list recent posts, then `instagram_get_media_insights` "
        "per post for engagement breakdown.\n"
        "4. Normalize the response into the schema defined by this script."
    ),
    "facebook": (
        "To collect Facebook metrics via Composio MCP:\n"
        "1. Ensure the Facebook Pages connector is active.\n"
        "2. Use `facebook_get_page_insights` with metrics: page_impressions, page_engaged_users, "
        "page_fans, page_views_total.\n"
        "3. Use `facebook_get_posts` and iterate `facebook_get_post_insights` per post.\n"
        "4. Normalize into the standard schema."
    ),
    "linkedin": (
        "To collect LinkedIn metrics via Composio MCP:\n"
        "1. Ensure the LinkedIn Pages connector is active.\n"
        "2. Use `linkedin_get_organization_statistics` for follower and visitor analytics.\n"
        "3. Use `linkedin_get_shares` to list posts and `linkedin_get_share_statistics` per post.\n"
        "4. Normalize into the standard schema."
    ),
    "tiktok": (
        "To collect TikTok metrics via Composio MCP:\n"
        "1. Ensure the TikTok Business connector is active.\n"
        "2. Use `tiktok_get_account_analytics` for profile-level data.\n"
        "3. Use `tiktok_get_videos` then `tiktok_get_video_analytics` per video.\n"
        "4. Normalize into the standard schema."
    ),
}

# ---------------------------------------------------------------------------
# Chrome automation guidance
# ---------------------------------------------------------------------------

CHROME_INSTRUCTIONS = {
    "instagram": (
        "Chrome automation steps for Instagram Insights:\n"
        "1. Navigate to https://www.instagram.com/{username}/\n"
        "2. Click 'Professional dashboard' or 'Insights'.\n"
        "3. Select the desired date range.\n"
        "4. Record: Accounts reached, Impressions, Profile visits, Website clicks, "
        "Followers (gained/lost).\n"
        "5. Go to 'Content You Shared' and for each post record: type, reach, likes, "
        "comments, shares, saves, plays (video).\n"
        "6. Export or transcribe data into CSV format matching POST_LEVEL_METRICS."
    ),
    "facebook": (
        "Chrome automation steps for Facebook Page Insights:\n"
        "1. Navigate to https://www.facebook.com/{page}/insights/\n"
        "2. Select 'Overview' and set date range.\n"
        "3. Record: Page reach, Page visits, Post engagement, New followers.\n"
        "4. Navigate to 'Posts' tab for per-post metrics.\n"
        "5. Export CSV from Facebook Insights (use the 'Export Data' button)."
    ),
    "linkedin": (
        "Chrome automation steps for LinkedIn Page Analytics:\n"
        "1. Navigate to https://www.linkedin.com/company/{company}/admin/analytics/\n"
        "2. Select 'Visitors', 'Followers', and 'Content' tabs.\n"
        "3. Record: Impressions, Unique visitors, Follower growth, Engagement rate.\n"
        "4. For each post: impressions, clicks, reactions, comments, shares.\n"
        "5. Use the 'Export' button on each analytics tab."
    ),
    "tiktok": (
        "Chrome automation steps for TikTok Analytics:\n"
        "1. Navigate to https://www.tiktok.com/analytics (Creator/Business account required).\n"
        "2. Select date range under 'Overview'.\n"
        "3. Record: Video views, Profile views, Followers (gained/lost), Likes.\n"
        "4. Go to 'Content' tab for per-video metrics: views, likes, comments, shares, "
        "average watch time, completion rate.\n"
        "5. TikTok does not offer CSV export; transcribe manually or use screen reader tools."
    ),
}


# ---------------------------------------------------------------------------
# CSV import
# ---------------------------------------------------------------------------

def load_csv_input(csv_path: str, platform: str) -> list[dict]:
    """Load manually collected post data from CSV and normalize to schema."""
    posts = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            post = empty_post_metrics()
            for key in post:
                if key in row and row[key] not in ("", None):
                    try:
                        post[key] = float(row[key]) if key != "post_id" and key != "published_at" and key != "content_type" else row[key]
                    except (ValueError, TypeError):
                        post[key] = row[key]
            posts.append(post)
    return posts


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_output(data: dict, output_path: str) -> None:
    """Write collected metrics to JSON."""
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"Metrics written to {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Collect social media metrics.")
    parser.add_argument("--platforms", type=str, default="all",
                        help="Comma-separated platforms or 'all'")
    parser.add_argument("--period", type=str, default="30d",
                        choices=list(PERIOD_DAYS.keys()),
                        help="Reporting period (default: 30d)")
    parser.add_argument("--output", type=str, default="metrics.json",
                        help="Output JSON path")
    parser.add_argument("--input", type=str, default=None,
                        help="Optional CSV with manually collected data")
    parser.add_argument("--account-name", type=str, default="brand",
                        help="Account/brand name label")
    args = parser.parse_args()

    platforms = PLATFORMS if args.platforms == "all" else [
        p.strip().lower() for p in args.platforms.split(",")
    ]
    for p in platforms:
        if p not in PLATFORMS:
            print(f"Error: Unknown platform '{p}'. Choose from: {', '.join(PLATFORMS)}")
            sys.exit(1)

    results = {
        "collection_metadata": {
            "collected_at": datetime.utcnow().isoformat() + "Z",
            "period": args.period,
            "period_days": PERIOD_DAYS[args.period],
            "period_start": (datetime.utcnow() - timedelta(days=PERIOD_DAYS[args.period])).strftime("%Y-%m-%d"),
            "period_end": datetime.utcnow().strftime("%Y-%m-%d"),
            "account_name": args.account_name,
            "platforms": platforms,
        },
        "platforms": {},
    }

    for platform in platforms:
        metrics = empty_account_metrics(platform, args.period, args.account_name)

        if args.input:
            metrics["posts"] = load_csv_input(args.input, platform)
            metrics["collection_method"] = "csv_import"
        else:
            metrics["collection_method"] = "template"
            # Print instructions for the operator or calling agent
            print(f"\n{'='*60}")
            print(f"COLLECTION INSTRUCTIONS FOR {platform.upper()}")
            print(f"{'='*60}")
            print(f"\n--- Option A: Composio MCP ---")
            print(COMPOSIO_INSTRUCTIONS.get(platform, "No Composio instructions available."))
            print(f"\n--- Option B: Chrome Automation ---")
            print(CHROME_INSTRUCTIONS.get(platform, "No Chrome instructions available."))
            print()

        results["platforms"][platform] = metrics

    write_output(results, args.output)

    print(f"\nCollection complete for {len(platforms)} platform(s).")
    print("If using template mode, populate the null fields in the output JSON with actual data.")
    print("Then run: python analyze_performance.py --metrics", args.output)


if __name__ == "__main__":
    main()
