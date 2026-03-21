#!/usr/bin/env python3
"""
analyze_competitors.py - Scrape and analyze competitor brands for positioning.

Fetches competitor websites, extracts brand signals, and creates competitive
positioning analysis including SWOT, positioning maps, content gaps, and
share of voice estimates.

This script combines automated extraction with structured templates.
For production competitor analysis, supplement with Chrome automation MCP
for social media data and WebSearch for market intelligence.

Usage:
    python analyze_competitors.py --competitors "https://comp1.com,https://comp2.com" --output competitor_analysis.json
    python analyze_competitors.py --input raw_competitor_data.json --output competitor_analysis.json

Requirements:
    pip install requests beautifulsoup4
"""

import argparse
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

try:
    import requests
    from bs4 import BeautifulSoup

    HAS_DEPS = True
except ImportError:
    HAS_DEPS = False


def load_env(env_path: str = ".env") -> None:
    """Load environment variables from .env file if it exists."""
    env_file = Path(env_path)
    if env_file.exists():
        with open(env_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    key, _, value = line.partition("=")
                    os.environ.setdefault(key.strip(), value.strip().strip("\"'"))


EMPTY_COMPETITOR_PROFILE = {
    "name": "",
    "url": "",
    "fetched_at": "",
    "brand_signals": {
        "tagline": "",
        "meta_description": "",
        "primary_colors": [],
        "fonts": [],
        "hero_text": [],
        "value_propositions": [],
    },
    "social_profiles": {
        "instagram": "",
        "facebook": "",
        "linkedin": "",
        "tiktok": "",
        "x_twitter": "",
        "youtube": "",
    },
    "social_metrics": {
        "instagram_followers": None,
        "facebook_followers": None,
        "linkedin_followers": None,
        "tiktok_followers": None,
    },
    "content_analysis": {
        "posting_frequency": "",
        "content_themes": [],
        "formats_used": [],
        "engagement_rate": None,
    },
    "positioning": {
        "price_tier": "",
        "personality_axis": "",
        "target_market": "",
        "key_differentiators": [],
    },
    "swot": {
        "strengths": [],
        "weaknesses": [],
        "opportunities": [],
        "threats": [],
    },
}


def fetch_page(url: str, timeout: int = 30) -> tuple[str, int]:
    """Fetch a web page. Returns (html, status_code)."""
    if not HAS_DEPS:
        return "", 0

    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/120.0.0.0 Safari/537.36"
        ),
    }
    try:
        response = requests.get(url, headers=headers, timeout=timeout)
        return response.text, response.status_code
    except requests.RequestException as e:
        print(f"Error fetching {url}: {e}", file=sys.stderr)
        return "", 0


def extract_competitor_signals(url: str) -> dict[str, Any]:
    """Extract basic brand signals from a competitor website.

    Args:
        url: The competitor website URL.

    Returns:
        Dictionary with extracted brand signals.
    """
    profile = json.loads(json.dumps(EMPTY_COMPETITOR_PROFILE))
    profile["url"] = url
    profile["fetched_at"] = datetime.now(timezone.utc).isoformat()

    # Derive company name from domain
    parsed = urlparse(url)
    domain_name = parsed.netloc.replace("www.", "").split(".")[0]
    profile["name"] = domain_name.capitalize()

    if not HAS_DEPS:
        profile["_error"] = "Dependencies not installed. pip install requests beautifulsoup4"
        return profile

    html, status = fetch_page(url)
    if not html:
        profile["_error"] = f"Failed to fetch {url} (status: {status})"
        return profile

    soup = BeautifulSoup(html, "html.parser")

    # Extract meta information
    for tag in soup.find_all("meta"):
        name = tag.get("name", "").lower()
        prop = tag.get("property", "").lower()
        content = tag.get("content", "")
        if name == "description" and content:
            profile["brand_signals"]["meta_description"] = content
        elif prop == "og:title" and content:
            profile["name"] = content.split("|")[0].split("-")[0].strip()

    # Page title as fallback name
    title_tag = soup.find("title")
    if title_tag and not profile["name"]:
        profile["name"] = title_tag.get_text(strip=True).split("|")[0].split("-")[0].strip()

    # Hero text
    h1 = soup.find("h1")
    if h1:
        profile["brand_signals"]["tagline"] = h1.get_text(strip=True)

    hero_section = soup.find(["section", "div"], class_=re.compile(
        r"hero|banner|jumbotron", re.IGNORECASE
    ))
    if hero_section:
        for el in hero_section.find_all(["h1", "h2", "p"])[:5]:
            text = el.get_text(strip=True)
            if text and len(text) > 5:
                profile["brand_signals"]["hero_text"].append(text)

    # Extract colors from inline CSS
    css_parts = []
    for style_tag in soup.find_all("style"):
        if style_tag.string:
            css_parts.append(style_tag.string)
    css_text = "\n".join(css_parts)

    hex_colors = list(set(re.findall(r"#(?:[0-9a-fA-F]{3}){1,2}\b", css_text)))
    profile["brand_signals"]["primary_colors"] = hex_colors[:10]

    # Extract fonts
    fonts = set()
    for match in re.findall(r"font-family\s*:\s*([^;]+)", css_text, re.IGNORECASE):
        for f in match.split(","):
            f = f.strip().strip("\"'")
            if f.lower() not in ("serif", "sans-serif", "monospace", "inherit", "system-ui"):
                fonts.add(f)

    for link in soup.find_all("link", href=True):
        if "fonts.googleapis.com" in link["href"]:
            family_match = re.search(r"family=([^&]+)", link["href"])
            if family_match:
                for fam in family_match.group(1).split("|"):
                    fonts.add(fam.split(":")[0].replace("+", " "))

    profile["brand_signals"]["fonts"] = list(fonts)

    # Detect social profile links
    social_patterns = {
        "instagram": r"instagram\.com/([^/?\s\"']+)",
        "facebook": r"facebook\.com/([^/?\s\"']+)",
        "linkedin": r"linkedin\.com/(company/[^/?\s\"']+)",
        "tiktok": r"tiktok\.com/@([^/?\s\"']+)",
        "x_twitter": r"(?:twitter|x)\.com/([^/?\s\"']+)",
        "youtube": r"youtube\.com/(?:c/|channel/|@)([^/?\s\"']+)",
    }
    page_text = str(soup)
    for platform, pattern in social_patterns.items():
        match = re.search(pattern, page_text, re.IGNORECASE)
        if match:
            profile["social_profiles"][platform] = match.group(1)

    return profile


def create_positioning_map(
    competitors: list[dict[str, Any]],
) -> dict[str, Any]:
    """Create a competitive positioning map structure.

    In practice, the axis positions should be determined by human analysis
    or LLM evaluation of the extracted brand signals.

    Returns:
        Positioning map structure with competitors placed.
    """
    positioning_map = {
        "axes": {
            "x": {"label": "Personality", "left": "Serious/Corporate", "right": "Playful/Casual"},
            "y": {"label": "Price Tier", "bottom": "Budget/Value", "top": "Premium/Luxury"},
        },
        "competitors": [],
        "_note": (
            "Position values (1-10) should be assigned based on analysis "
            "of brand signals, pricing, and communication style. "
            "Use an LLM to evaluate positions from the extracted data."
        ),
    }

    for comp in competitors:
        positioning_map["competitors"].append({
            "name": comp.get("name", "Unknown"),
            "x_position": None,  # 1=Serious, 10=Playful
            "y_position": None,  # 1=Budget, 10=Premium
            "quadrant": "",  # To be filled: TL, TR, BL, BR
        })

    return positioning_map


def identify_content_gaps(competitors: list[dict[str, Any]]) -> dict[str, Any]:
    """Identify content themes and formats that represent opportunities.

    Args:
        competitors: List of competitor profiles with content analysis.

    Returns:
        Content gap analysis structure.
    """
    all_themes: dict[str, int] = {}
    all_formats: dict[str, int] = {}

    for comp in competitors:
        content = comp.get("content_analysis", {})
        for theme in content.get("content_themes", []):
            all_themes[theme] = all_themes.get(theme, 0) + 1
        for fmt in content.get("formats_used", []):
            all_formats[fmt] = all_formats.get(fmt, 0) + 1

    total_comps = len(competitors) if competitors else 1

    return {
        "theme_coverage": {
            theme: {
                "competitors_using": count,
                "saturation": round(count / total_comps * 100, 1),
                "opportunity_type": (
                    "crowded" if count >= total_comps * 0.7
                    else "moderate" if count >= total_comps * 0.3
                    else "underserved"
                ),
            }
            for theme, count in sorted(all_themes.items(), key=lambda x: x[1], reverse=True)
        },
        "format_coverage": {
            fmt: {
                "competitors_using": count,
                "saturation": round(count / total_comps * 100, 1),
            }
            for fmt, count in sorted(all_formats.items(), key=lambda x: x[1], reverse=True)
        },
        "gap_recommendations": [],  # To be filled by LLM analysis
        "_note": (
            "Underserved themes and unused formats represent potential "
            "content gap opportunities. Use LLM analysis to generate "
            "specific recommendations based on the target brand's strengths."
        ),
    }


def estimate_share_of_voice(competitors: list[dict[str, Any]]) -> dict[str, Any]:
    """Estimate relative share of voice across competitors.

    Uses available social metrics to approximate share of voice.
    Actual SOV measurement requires paid tools or extensive data collection.
    """
    sov = {
        "method": "follower_weighted_estimate",
        "platforms": {},
        "overall": {},
        "_note": (
            "Share of voice estimated from available follower counts. "
            "For accurate SOV, use social listening tools or manual "
            "engagement tracking over a 30-day period."
        ),
    }

    platforms = ["instagram", "facebook", "linkedin", "tiktok"]

    for platform in platforms:
        metric_key = f"{platform}_followers"
        total = 0
        entries = []

        for comp in competitors:
            followers = comp.get("social_metrics", {}).get(metric_key)
            if followers and followers > 0:
                entries.append({"name": comp["name"], "followers": followers})
                total += followers

        if entries and total > 0:
            sov["platforms"][platform] = {
                "total_followers": total,
                "breakdown": [
                    {
                        "name": e["name"],
                        "followers": e["followers"],
                        "share_pct": round(e["followers"] / total * 100, 1),
                    }
                    for e in sorted(entries, key=lambda x: x["followers"], reverse=True)
                ],
            }

    return sov


def analyze_competitors(
    urls: list[str],
    pre_collected: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Run full competitor analysis.

    Args:
        urls: List of competitor website URLs.
        pre_collected: Optional pre-collected competitor data (overrides URL scraping).

    Returns:
        Complete competitor analysis document.
    """
    competitors: list[dict[str, Any]]

    if pre_collected:
        competitors = pre_collected
    else:
        competitors = []
        for url in urls:
            print(f"Analyzing: {url}", file=sys.stderr)
            comp = extract_competitor_signals(url)
            competitors.append(comp)

    result = {
        "analyzed_at": datetime.now(timezone.utc).isoformat(),
        "total_competitors": len(competitors),
        "competitors": competitors,
        "positioning_map": create_positioning_map(competitors),
        "content_gaps": identify_content_gaps(competitors),
        "share_of_voice": estimate_share_of_voice(competitors),
        "summary": {
            "competitors_analyzed": [c.get("name", "Unknown") for c in competitors],
            "common_themes": [],
            "differentiation_opportunities": [],
            "_note": (
                "Use LLM analysis to fill in common_themes and "
                "differentiation_opportunities based on the extracted data."
            ),
        },
    }

    return result


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Analyze competitor brands for competitive positioning.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            '  python analyze_competitors.py --competitors "https://comp1.com,https://comp2.com"\n'
            "  python analyze_competitors.py --input raw_data.json --output analysis.json\n"
            "\n"
            "For richer analysis, supplement with:\n"
            "  - Chrome automation MCP for social media metrics\n"
            "  - WebSearch for market intelligence\n"
            "  - firecrawl MCP for deeper website crawling"
        ),
    )

    parser.add_argument(
        "--competitors",
        default="",
        help="Comma-separated competitor website URLs",
    )
    parser.add_argument(
        "--input",
        default="",
        help="Path to pre-collected competitor data JSON (alternative to --competitors)",
    )
    parser.add_argument(
        "--output",
        default="competitor_analysis.json",
        help="Output JSON file path (default: competitor_analysis.json)",
    )

    args = parser.parse_args()

    load_env()

    pre_collected = None
    urls: list[str] = []

    if args.input:
        input_path = Path(args.input)
        if not input_path.exists():
            print(f"Error: Input file not found: {args.input}", file=sys.stderr)
            sys.exit(1)
        with open(input_path) as f:
            data = json.load(f)
        if isinstance(data, list):
            pre_collected = data
        elif isinstance(data, dict) and "competitors" in data:
            pre_collected = data["competitors"]
        else:
            pre_collected = [data]
        print(f"Loaded {len(pre_collected)} competitor(s) from {args.input}", file=sys.stderr)
    elif args.competitors:
        urls = [u.strip() for u in args.competitors.split(",") if u.strip()]
        if not urls:
            print("Error: No valid URLs provided.", file=sys.stderr)
            sys.exit(1)
        print(f"Analyzing {len(urls)} competitor(s)...", file=sys.stderr)
    else:
        parser.print_help()
        sys.exit(1)

    result = analyze_competitors(urls, pre_collected)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"\nAnalysis saved to: {output_path}", file=sys.stderr)
    print(f"Competitors analyzed: {result['total_competitors']}", file=sys.stderr)


if __name__ == "__main__":
    main()
