---
name: brand-identity-miner
description: Extract comprehensive brand identity from websites and social media profiles. Produces a brand_dna.json document containing visual identity (colors, fonts, logo), brand voice, values, target audience, and competitive positioning. Used as foundation for all social media content creation skills.
allowed-tools: [Bash, Read, Write, Edit, WebSearch, WebFetch, Agent]
license: MIT license
metadata:
    skill-author: K-Dense Inc.
---

# Brand Identity Miner

Extract and synthesize a comprehensive brand identity profile from a company's digital presence. This skill crawls websites, analyzes social media profiles, maps competitor positioning, and produces a unified `brand_dna.json` document that serves as the single source of truth for all downstream content creation.

## Overview

Brand Identity Miner performs deep analysis across four dimensions:

1. **Website Mining** -- Extracts visual identity (colors, typography, imagery style), brand voice (copy tone, vocabulary, formality), and strategic messaging (value propositions, mission, target audience signals) from a company's website.

2. **Social Profile Mining** -- Analyzes existing social media presence across Facebook, Instagram, LinkedIn, and TikTok to identify posting patterns, audience engagement, visual consistency, and content themes.

3. **Competitor Analysis** -- Scrapes 3-5 competitor brands to map relative positioning, identify content gaps, and surface differentiation opportunities.

4. **Brand DNA Synthesis** -- Merges all collected data into a structured `brand_dna.json` file that downstream skills consume for content generation, scheduling, and performance tracking.

## When to Use This Skill

- **New client onboarding** -- Run the full pipeline to build a brand profile before creating any content.
- **Brand refresh** -- Re-mine after a rebrand, website redesign, or strategic pivot to update the DNA file.
- **Competitive audit** -- Run competitor analysis independently when market positioning needs reassessment.
- **Content strategy kickoff** -- Generate brand DNA to feed into content calendar and post generation skills.
- **Brand consistency check** -- Compare current social output against the DNA file to identify drift.

## Quick Start

### Mine a single website

```
/brand-identity-miner Analyze the website https://example.com and extract their brand identity into brand_dna.json
```

### Full pipeline with competitors

```
/brand-identity-miner Mine brand identity for https://acme.com, analyze their Instagram @acme and LinkedIn /company/acme, compare against https://competitor1.com and https://competitor2.com, then generate the complete brand DNA.
```

### Update existing brand DNA with new social data

```
/brand-identity-miner Update the brand DNA at ./brand_dna.json with fresh analysis of @acme on Instagram and TikTok
```

## Detailed Workflows

### 1. Website Mining

The website mining phase extracts three categories of brand signals from the target URL.

#### Visual Identity Extraction

Use WebFetch or firecrawl MCP to retrieve the website HTML and CSS. Extract:

- **Primary colors**: Parse CSS custom properties (`--primary`, `--brand-*`), inspect `background-color` and `color` on header, hero, and CTA elements. Look for the dominant 2-3 colors.
- **Secondary colors**: Background tones, border colors, subtle UI element colors.
- **Accent colors**: CTA buttons, hover states, notification badges.
- **Typography**: Parse `font-family` declarations. Identify heading vs. body fonts. Check Google Fonts links or `@font-face` declarations.
- **Logo**: Describe the logo from the `<img>` in the header/nav. Note color, shape, wordmark vs. icon.
- **Photography style**: Analyze hero images and about page photos. Classify as lifestyle, product-focused, abstract, people-centric, or illustration-heavy.
- **Visual mood**: Summarize the overall aesthetic (minimal, bold, corporate, artisan, tech-forward, warm, etc.).

#### Voice Analysis

Scrape text content from these priority pages (in order):

1. Homepage hero section and subheadings
2. About / Our Story page
3. Services or Products overview
4. Blog post titles and first paragraphs (sample 5-10)
5. FAQ page

Analyze the collected copy for:

- **Tone**: Map on four axes -- Formal/Casual, Serious/Playful, Technical/Simple, Authoritative/Friendly.
- **Personality traits**: 3-5 adjectives (e.g., "innovative", "approachable", "bold").
- **Vocabulary level**: Technical (jargon-heavy), Casual (conversational), or Mixed.
- **Formality**: Formal (third person, no contractions), Semiformal (contractions ok, second person), Casual (slang, first person plural).
- **Sample phrases**: Extract 5-10 characteristic phrases that embody the brand voice.
- **Do/Don't say lists**: Infer from copy patterns what language the brand embraces and avoids.

#### Strategic Messaging

From the same pages, extract:

- **Tagline**: The primary headline or positioning statement.
- **Mission**: Usually on the About page.
- **Values**: Look for explicit values sections, or infer from repeated themes.
- **Unique selling propositions**: What the brand claims makes it different.
- **Target audience signals**: Industry terms, problem statements, imagery choices that indicate who they serve.

#### Implementation

```bash
# Using the helper script
python scripts/mine_website.py --url https://example.com --output website_analysis.json

# Or via firecrawl MCP for richer extraction
# firecrawl_scrape(url="https://example.com", formats=["markdown", "html"])
```

The script outputs a JSON file with raw extracted data. The synthesis step later merges this with other sources.

### 2. Social Profile Mining

Analyze existing social media profiles to understand current brand expression and audience response.

#### Supported Platforms

| Platform  | Key Signals                                      |
|-----------|--------------------------------------------------|
| Instagram | Grid aesthetic, caption voice, hashtag strategy   |
| Facebook  | Page about, post mix, community engagement        |
| LinkedIn  | Company description, thought leadership tone      |
| TikTok    | Content themes, audio choices, trend participation |

#### For Each Platform, Extract

- **Bio / About**: The profile description, links, and positioning statement.
- **Visual consistency**: Color palette used in graphics, filter choices, grid layout patterns.
- **Content themes**: Categorize the last 20-30 posts into themes (educational, promotional, behind-the-scenes, user-generated, entertainment).
- **Posting frequency**: Average posts per week, time-of-day patterns.
- **Hashtag strategy**: Most-used hashtags, branded hashtags, community hashtags.
- **Engagement patterns**: Which post types generate the most interaction.
- **Audience signals**: Comment sentiment, follower demographics if visible.

#### Implementation

```bash
python scripts/mine_social_profiles.py \
  --profiles "instagram:@acme" "linkedin:/company/acme" "tiktok:@acme" \
  --output social_analysis.json
```

Social mining typically requires Chrome automation MCP for platforms that block scraping. The script documents the expected data structure; actual extraction uses browser automation tools.

### 3. Competitor Analysis

Map the competitive landscape to identify positioning opportunities and content gaps.

#### Process

1. **Select 3-5 competitors**: Direct competitors in the same market segment.
2. **Mine each competitor website**: Run the website mining workflow on each.
3. **Analyze social presence**: Note follower counts, posting frequency, content themes.
4. **Build positioning map**: Plot competitors on two key axes:
   - **Price/Value axis**: Premium <---> Budget
   - **Personality axis**: Serious/Corporate <---> Playful/Casual
5. **Content gap analysis**: Identify topics competitors cover that the brand does not, and vice versa.
6. **SWOT per competitor**: Strengths, Weaknesses, Opportunities, Threats specific to social media.

#### Implementation

```bash
python scripts/analyze_competitors.py \
  --competitors "https://competitor1.com,https://competitor2.com,https://competitor3.com" \
  --output competitor_analysis.json
```

See `references/competitor_analysis_guide.md` for the complete SWOT framework and positioning methodology.

### 4. Brand DNA Synthesis

Merge all analysis outputs into the final `brand_dna.json` document.

#### Process

1. Load website analysis, social analysis, and competitor analysis JSONs.
2. Resolve conflicts (e.g., different color palettes detected on website vs. social).
3. Apply the Brand Archetype framework to classify the brand (see `references/brand_analysis_framework.md`).
4. Define 3-5 content pillars based on values, audience needs, and competitive whitespace.
5. Generate the unified `brand_dna.json`.

#### Implementation

```bash
python scripts/generate_brand_dna.py \
  --website website_analysis.json \
  --social social_analysis.json \
  --competitors competitor_analysis.json \
  --output brand_dna.json
```

## Brand DNA JSON Schema

The output `brand_dna.json` follows a strict schema documented in `references/brand_dna_schema.md`. Key sections:

- `visual_identity` -- Colors, fonts, logo, photography style
- `voice` -- Tone, personality, vocabulary, do/don't say lists
- `values` and `mission` -- Core beliefs and purpose
- `target_audience` -- Demographics, psychographics, pain points
- `competitors` -- Mapped competitors with positioning
- `content_pillars` -- 3-5 strategic content themes
- `brand_archetype` -- Jungian archetype classification

## Integration with Other Skills

The `brand_dna.json` file is the primary input for downstream social media skills:

| Skill                     | How It Uses brand_dna.json                            |
|---------------------------|-------------------------------------------------------|
| Content Calendar Planner  | Uses content pillars and audience to schedule themes   |
| Post Generator            | Uses voice, tone, and sample phrases for copy          |
| Visual Template Builder   | Uses colors, fonts, and photography style              |
| Hashtag Strategist        | Uses industry, audience, and competitor hashtags        |
| Analytics Reporter        | Compares performance against competitor benchmarks      |

## File Reference

| File                                         | Purpose                                     |
|----------------------------------------------|---------------------------------------------|
| `SKILL.md`                                   | This file -- skill documentation            |
| `references/brand_analysis_framework.md`     | Archetypes, voice dimensions, methodology   |
| `references/competitor_analysis_guide.md`    | SWOT, positioning maps, gap analysis        |
| `references/brand_dna_schema.md`             | Full JSON schema with field descriptions    |
| `scripts/mine_website.py`                    | Website extraction script/template          |
| `scripts/mine_social_profiles.py`            | Social profile analysis script/template     |
| `scripts/analyze_competitors.py`             | Competitor scraping and positioning          |
| `scripts/generate_brand_dna.py`              | Merge all sources into brand_dna.json       |

## Notes

- Always respect `robots.txt` and rate limits when crawling.
- Social platform scraping may require authentication; prefer public data and official APIs where available.
- The brand DNA should be reviewed and refined with the client before using it for content generation.
- Re-run the mining pipeline quarterly or after significant brand changes.
