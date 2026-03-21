#!/usr/bin/env python3
"""
mine_website.py - Extract brand identity signals from a company website.

Fetches the target website and extracts visual identity (colors, fonts, logo),
brand voice (tone, vocabulary, formality), and strategic messaging (tagline,
mission, values, USPs).

This script serves as both a functional extraction tool and a template
documenting the extraction logic. For richer results, use firecrawl MCP
or Chrome automation MCP to fetch pages before running analysis.

Usage:
    python mine_website.py --url https://example.com --output website_analysis.json
    python mine_website.py --url https://example.com --pages about,services --output analysis.json

Requirements:
    pip install requests beautifulsoup4

Environment Variables:
    FIRECRAWL_API_KEY - Optional. If set, uses firecrawl API for better extraction.
"""

import argparse
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urljoin, urlparse

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


def extract_colors_from_css(css_text: str) -> dict[str, list[str]]:
    """Extract color values from CSS text.

    Looks for hex colors, rgb/rgba values, and CSS custom properties
    that appear to be color definitions.

    Returns:
        Dictionary with 'all_colors' list and 'custom_properties' dict.
    """
    colors: dict[str, Any] = {"all_colors": [], "custom_properties": {}}

    # Extract hex colors
    hex_pattern = r"#(?:[0-9a-fA-F]{3}){1,2}\b"
    hex_colors = re.findall(hex_pattern, css_text)
    colors["all_colors"].extend(hex_colors)

    # Extract CSS custom properties that look like colors
    custom_prop_pattern = r"--([\w-]*(?:color|bg|background|brand|primary|secondary|accent)[\w-]*)\s*:\s*([^;]+)"
    for name, value in re.findall(custom_prop_pattern, css_text, re.IGNORECASE):
        colors["custom_properties"][f"--{name}"] = value.strip()

    # Extract rgb/rgba colors
    rgb_pattern = r"rgba?\(\s*\d+\s*,\s*\d+\s*,\s*\d+(?:\s*,\s*[\d.]+)?\s*\)"
    rgb_colors = re.findall(rgb_pattern, css_text)
    colors["all_colors"].extend(rgb_colors)

    # Deduplicate while preserving order
    seen = set()
    unique = []
    for c in colors["all_colors"]:
        normalized = c.lower()
        if normalized not in seen:
            seen.add(normalized)
            unique.append(c)
    colors["all_colors"] = unique

    return colors


def extract_fonts_from_css(css_text: str) -> dict[str, list[str]]:
    """Extract font family declarations from CSS text.

    Returns:
        Dictionary with 'font_families' list and 'google_fonts' list.
    """
    fonts: dict[str, list[str]] = {"font_families": [], "google_fonts": []}

    # Extract font-family declarations
    font_pattern = r"font-family\s*:\s*([^;]+)"
    for match in re.findall(font_pattern, css_text, re.IGNORECASE):
        families = [f.strip().strip("\"'") for f in match.split(",")]
        for family in families:
            if family and family not in fonts["font_families"]:
                # Skip generic families
                if family.lower() not in (
                    "serif",
                    "sans-serif",
                    "monospace",
                    "cursive",
                    "fantasy",
                    "system-ui",
                    "inherit",
                ):
                    fonts["font_families"].append(family)

    return fonts


def extract_fonts_from_html(soup: BeautifulSoup) -> list[str]:
    """Extract Google Fonts or other external font references from HTML."""
    google_fonts = []

    for link in soup.find_all("link", href=True):
        href = link["href"]
        if "fonts.googleapis.com" in href:
            # Parse font names from Google Fonts URL
            family_match = re.search(r"family=([^&]+)", href)
            if family_match:
                families = family_match.group(1).split("|")
                for fam in families:
                    name = fam.split(":")[0].replace("+", " ")
                    if name not in google_fonts:
                        google_fonts.append(name)
        elif "fonts.adobe.com" in href or "use.typekit.net" in href:
            google_fonts.append(f"[Adobe Fonts: {href}]")

    return google_fonts


def extract_meta_info(soup: BeautifulSoup) -> dict[str, str]:
    """Extract meta tags relevant to brand identity."""
    meta = {}

    # Standard meta tags
    for tag in soup.find_all("meta"):
        name = tag.get("name", "").lower()
        prop = tag.get("property", "").lower()
        content = tag.get("content", "")

        if name == "description" and content:
            meta["meta_description"] = content
        elif name == "keywords" and content:
            meta["meta_keywords"] = content
        elif prop == "og:title" and content:
            meta["og_title"] = content
        elif prop == "og:description" and content:
            meta["og_description"] = content
        elif name == "theme-color" and content:
            meta["theme_color"] = content

    # Title tag
    title = soup.find("title")
    if title:
        meta["page_title"] = title.get_text(strip=True)

    return meta


def extract_text_content(soup: BeautifulSoup) -> dict[str, Any]:
    """Extract text content from key page sections for voice analysis."""
    content: dict[str, Any] = {
        "hero_text": [],
        "headings": [],
        "paragraphs": [],
        "cta_text": [],
    }

    # Extract hero section text (first large heading + nearby text)
    hero_candidates = soup.find_all(["section", "div"], class_=re.compile(
        r"hero|banner|jumbotron|masthead|splash", re.IGNORECASE
    ))
    for hero in hero_candidates[:2]:
        for el in hero.find_all(["h1", "h2", "p", "span"]):
            text = el.get_text(strip=True)
            if text and len(text) > 5:
                content["hero_text"].append(text)

    # If no hero section found, grab the first h1
    if not content["hero_text"]:
        h1 = soup.find("h1")
        if h1:
            content["hero_text"].append(h1.get_text(strip=True))

    # Extract all headings
    for level in ["h1", "h2", "h3"]:
        for heading in soup.find_all(level):
            text = heading.get_text(strip=True)
            if text and len(text) > 3:
                content["headings"].append({"level": level, "text": text})

    # Extract paragraph text (first 20 paragraphs for analysis)
    for p in soup.find_all("p")[:20]:
        text = p.get_text(strip=True)
        if text and len(text) > 20:
            content["paragraphs"].append(text)

    # Extract CTA button text
    for btn in soup.find_all(["a", "button"], class_=re.compile(
        r"btn|button|cta", re.IGNORECASE
    )):
        text = btn.get_text(strip=True)
        if text and len(text) > 1:
            content["cta_text"].append(text)

    return content


def extract_logo_info(soup: BeautifulSoup, base_url: str) -> dict[str, str]:
    """Extract logo information from the page header/navigation."""
    logo_info: dict[str, str] = {}

    # Look for logo in header/nav
    header = soup.find(["header", "nav"])
    if header:
        # Check for img with logo-related attributes
        for img in header.find_all("img"):
            src = img.get("src", "")
            alt = img.get("alt", "")
            classes = " ".join(img.get("class", []))

            if any(
                kw in (src + alt + classes).lower()
                for kw in ["logo", "brand", "site-icon"]
            ):
                logo_info["logo_src"] = urljoin(base_url, src)
                logo_info["logo_alt"] = alt
                break

        # Check for SVG logo
        if not logo_info:
            svg = header.find("svg")
            if svg:
                logo_info["logo_type"] = "inline-svg"
                title = svg.find("title")
                if title:
                    logo_info["logo_alt"] = title.get_text(strip=True)

    return logo_info


def fetch_page(url: str, timeout: int = 30) -> tuple[str, int]:
    """Fetch a web page and return (html_content, status_code).

    If FIRECRAWL_API_KEY is set, logs a note that firecrawl MCP
    should be preferred for production use.
    """
    if os.environ.get("FIRECRAWL_API_KEY"):
        print(
            "NOTE: FIRECRAWL_API_KEY detected. For richer extraction, use "
            "firecrawl MCP: firecrawl_scrape(url=url, formats=['markdown', 'html'])",
            file=sys.stderr,
        )

    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/120.0.0.0 Safari/537.36"
        ),
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "Accept-Language": "en-US,en;q=0.5",
    }

    try:
        response = requests.get(url, headers=headers, timeout=timeout)
        return response.text, response.status_code
    except requests.RequestException as e:
        print(f"Error fetching {url}: {e}", file=sys.stderr)
        return "", 0


def extract_inline_styles(soup: BeautifulSoup) -> str:
    """Gather all inline and embedded CSS from the page."""
    css_parts = []

    # Embedded <style> tags
    for style_tag in soup.find_all("style"):
        if style_tag.string:
            css_parts.append(style_tag.string)

    # Inline style attributes on key elements
    for el in soup.find_all(style=True):
        css_parts.append(el["style"])

    return "\n".join(css_parts)


def analyze_website(url: str, additional_pages: list[str] | None = None) -> dict[str, Any]:
    """Run full website analysis and return structured results.

    Args:
        url: The website URL to analyze.
        additional_pages: Optional list of page paths to also analyze (e.g., ['about', 'services']).

    Returns:
        Dictionary containing all extracted brand signals.
    """
    if not HAS_DEPS:
        return {
            "error": "Missing dependencies. Install with: pip install requests beautifulsoup4",
            "note": (
                "Alternatively, use firecrawl MCP or WebFetch to retrieve "
                "the page HTML, then pass it to the analysis functions."
            ),
        }

    parsed_url = urlparse(url)
    base_url = f"{parsed_url.scheme}://{parsed_url.netloc}"

    result: dict[str, Any] = {
        "url": url,
        "analyzed_at": datetime.now(timezone.utc).isoformat(),
        "pages_analyzed": [url],
        "visual_identity": {
            "colors": {},
            "fonts": {"css_fonts": [], "external_fonts": []},
            "logo": {},
        },
        "meta": {},
        "content": {},
        "additional_pages": {},
    }

    # Fetch and parse main page
    html, status = fetch_page(url)
    if not html:
        result["error"] = f"Failed to fetch {url} (status: {status})"
        return result

    soup = BeautifulSoup(html, "html.parser")

    # Extract CSS (inline + embedded)
    css_text = extract_inline_styles(soup)

    # Visual identity
    result["visual_identity"]["colors"] = extract_colors_from_css(css_text)
    result["visual_identity"]["fonts"]["css_fonts"] = extract_fonts_from_css(css_text)[
        "font_families"
    ]
    result["visual_identity"]["fonts"]["external_fonts"] = extract_fonts_from_html(soup)
    result["visual_identity"]["logo"] = extract_logo_info(soup, base_url)

    # Meta information
    result["meta"] = extract_meta_info(soup)

    # Text content for voice analysis
    result["content"] = extract_text_content(soup)

    # Analyze additional pages
    if additional_pages:
        for page_path in additional_pages:
            page_url = urljoin(base_url + "/", page_path.strip("/"))
            result["pages_analyzed"].append(page_url)

            page_html, page_status = fetch_page(page_url)
            if page_html:
                page_soup = BeautifulSoup(page_html, "html.parser")
                result["additional_pages"][page_path] = {
                    "meta": extract_meta_info(page_soup),
                    "content": extract_text_content(page_soup),
                }

    # Summary statistics
    result["summary"] = {
        "total_colors_found": len(
            result["visual_identity"]["colors"].get("all_colors", [])
        ),
        "total_fonts_found": len(
            result["visual_identity"]["fonts"]["css_fonts"]
        )
        + len(result["visual_identity"]["fonts"]["external_fonts"]),
        "total_headings": len(result["content"].get("headings", [])),
        "total_paragraphs": len(result["content"].get("paragraphs", [])),
        "has_logo": bool(result["visual_identity"]["logo"]),
        "pages_analyzed": len(result["pages_analyzed"]),
    }

    return result


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Extract brand identity signals from a website.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python mine_website.py --url https://example.com --output analysis.json\n"
            "  python mine_website.py --url https://example.com --pages about,services\n"
            "\n"
            "For richer extraction, set FIRECRAWL_API_KEY in your .env file\n"
            "or use firecrawl MCP / WebFetch tools to retrieve page content."
        ),
    )

    parser.add_argument(
        "--url",
        required=True,
        help="Website URL to analyze (e.g., https://example.com)",
    )
    parser.add_argument(
        "--output",
        default="website_analysis.json",
        help="Output JSON file path (default: website_analysis.json)",
    )
    parser.add_argument(
        "--pages",
        default="",
        help="Comma-separated additional page paths to analyze (e.g., about,services,blog)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="HTTP request timeout in seconds (default: 30)",
    )

    args = parser.parse_args()

    load_env()

    additional_pages = [p.strip() for p in args.pages.split(",") if p.strip()] if args.pages else None

    print(f"Analyzing website: {args.url}", file=sys.stderr)
    result = analyze_website(args.url, additional_pages)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"Analysis saved to: {output_path}", file=sys.stderr)

    if "error" in result:
        print(f"WARNING: {result['error']}", file=sys.stderr)
        sys.exit(1)

    summary = result.get("summary", {})
    print(f"Colors found: {summary.get('total_colors_found', 0)}", file=sys.stderr)
    print(f"Fonts found: {summary.get('total_fonts_found', 0)}", file=sys.stderr)
    print(f"Pages analyzed: {summary.get('pages_analyzed', 0)}", file=sys.stderr)


if __name__ == "__main__":
    main()
