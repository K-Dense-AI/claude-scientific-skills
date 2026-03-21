---
name: social-media-orchestrator
description: Master orchestrator for social media promotion campaigns. Coordinates brand-identity-miner, viral-content-creator, social-media-visual-generator, social-media-copywriting, social-media-scheduler, and social-media-analytics into cohesive multi-platform campaigns. Entry point for all social media tasks — from brand discovery to content creation, publishing, and performance optimization.
allowed-tools: [Bash, Read, Write, Edit, WebSearch, WebFetch, Agent]
license: MIT license
metadata:
    skill-author: K-Dense Inc.
---

# Social Media Orchestrator

Master orchestrator that turns a brand URL into a complete, multi-platform social media campaign. This skill coordinates 6 specialized sub-skills through a structured 7-phase pipeline, managing the full lifecycle from brand discovery through content creation, publishing, and performance optimization.

## Architecture Overview

```
                         social-media-orchestrator
                                   |
            +----------+-----------+-----------+----------+----------+
            |          |           |           |          |          |
     brand-identity  viral-content  copywriting  visual-gen  scheduler  analytics
        -miner        -creator
            |                                                    |
        brand_dna.json -----> consumed by all skills -----> report.json
```

The orchestrator operates as a state machine, advancing through phases sequentially while allowing feedback loops between analytics and content planning for continuous optimization.

---

## The 7-Phase Campaign Pipeline

### Phase 1: Brand Discovery

**Objective**: Extract a complete brand identity profile from public digital presence.

**Invoke**: `brand-identity-miner` skill

**Inputs**:
- Website URL (required)
- Social media handles (optional — will be discovered if not provided)
- Industry vertical (optional — will be inferred)

**Process**:
1. Scrape the brand website using firecrawl or Chrome automation to extract:
   - Visual identity (colors, typography, logo usage)
   - Voice and tone patterns from copy
   - Product/service offerings
   - Mission, vision, values
2. Analyze existing social media profiles across platforms
3. Research 3-5 direct competitors for positioning context
4. Identify target audience signals from content and engagement patterns

**Output**: `brand_dna.json` — the central contract consumed by every downstream skill. Structure:
```json
{
  "brand_name": "...",
  "tagline": "...",
  "industry": "...",
  "mission": "...",
  "visual_identity": {
    "primary_colors": ["#hex1", "#hex2"],
    "secondary_colors": ["#hex3"],
    "typography": {"heading": "...", "body": "..."},
    "logo_url": "...",
    "visual_style": "minimalist | bold | playful | corporate | luxury"
  },
  "voice": {
    "tone": "professional | casual | witty | authoritative | empathetic",
    "personality_traits": ["innovative", "approachable"],
    "vocabulary_level": "simple | moderate | technical",
    "do_list": ["Use active voice", "Include data points"],
    "dont_list": ["Avoid jargon", "No clickbait"]
  },
  "audience": {
    "primary_demographic": "...",
    "age_range": "25-40",
    "interests": [],
    "pain_points": [],
    "platforms": ["instagram", "linkedin"]
  },
  "competitors": [
    {"name": "...", "strengths": [], "weaknesses": [], "social_handles": {}}
  ],
  "content_pillars": ["education", "behind-the-scenes", "customer-stories"],
  "unique_selling_points": []
}
```

**MCP Tools Used**: firecrawl-scraper, Chrome automation (Claude in Chrome), tavily-web

---

### Phase 2: Strategy Planning

**Objective**: Define campaign structure, goals, platforms, and timeline.

**Reference**: See `references/campaign_frameworks.md` for detailed frameworks.

**Process**:
1. Select campaign type based on business objectives:
   - **Product Launch** — Tease, Announce, Sustain
   - **Brand Awareness** — Top-of-funnel reach
   - **Engagement Growth** — Community building
   - **Lead Generation** — Conversion-focused
   - **Event Promotion** — Before, During, After
   - **Seasonal/Holiday** — Timely cultural moments
2. Define SMART goals:
   - Specific: "Increase Instagram followers by 2,000"
   - Measurable: Tied to a numeric KPI
   - Achievable: Benchmarked against industry averages
   - Relevant: Aligned with business objectives
   - Time-bound: Within the campaign duration
3. Select platforms using decision matrix in `references/platform_ecosystem.md`
4. Define 3-5 content pillars aligned with brand_dna.json
5. Set budget allocation across platforms and content types
6. Establish posting cadence per platform

**Output**: `campaign_plan.json` containing goals, platforms, timeline, content pillars, budget, and KPI targets.

---

### Phase 3: Content Planning

**Objective**: Generate a complete content calendar with individual content briefs.

**Invoke**:
- `viral-content-creator` skill — for content ideation and brief generation
- `social-media-scheduler` skill — for calendar structure and optimal timing

**Process**:
1. Apply the **60/30/10 content mix** rule:
   - **60% Value Content**: Educational, entertaining, or inspiring posts that serve the audience
   - **30% Curated/Shared Content**: Industry news, user-generated content, partner content
   - **10% Promotional**: Direct sales, offers, product features
2. Generate content briefs for each post including:
   - Topic and angle
   - Platform and format (image, carousel, video, story, reel)
   - Key message and CTA
   - Visual direction
   - Hashtag strategy
   - Optimal posting time
3. Map content to campaign phases and dates
4. Ensure variety across formats and themes

**Output**:
- `content_calendar.json` — date-by-date posting schedule
- `content_briefs/*.json` — individual brief per content piece

---

### Phase 4: Content Creation

**Objective**: Produce all creative assets (copy + visuals) aligned to brand identity.

**Invoke**:
- `social-media-copywriting` skill — captions, CTAs, hashtag sets, alt text
- `social-media-visual-generator` skill — images, carousels, video storyboards

**Process**:
1. For each content brief:
   a. Generate platform-specific copy (caption, headline, CTA)
   b. Generate visual assets matching brand_dna.json guidelines
   c. Create platform-adapted variations (e.g., square for Instagram, vertical for Stories/Reels)
2. Brand consistency check: every output validated against brand_dna.json colors, voice, and style
3. Batch production: create content in weekly batches for efficiency

**MCP Tools**:
- Copy: Claude native generation with brand voice constraints
- Images: `fal-generate` (Flux models), `fal-image-edit` (refinement), `generate-image` skill
- Carousels: Gamma MCP (`mcp__gamma__generate`) for multi-slide presentations
- Video: `fal-audio` for voiceover, storyboard generation for video briefs
- Extended: Gemini API via `GOOGLE_API_KEY` in `.env` for additional image/video generation

**Output**:
- `captions/*.json` — copy assets per post
- `visuals/*.png` — image assets
- `videos/*.mp4` — video assets (or storyboards if video production is manual)

---

### Phase 5: Review & Approval

**Objective**: Quality assurance and user sign-off before publishing.

**Process**:
1. Present content batch to user in a structured review format:
   - Post preview (visual + copy together)
   - Platform destination
   - Scheduled date/time
   - Hashtags and mentions
2. Run automated brand consistency check:
   - Color palette adherence in visuals
   - Voice/tone alignment in copy
   - CTA presence and clarity
   - Hashtag relevance and count
3. Quality checklist per post:
   - [ ] Visual quality meets platform standards (resolution, aspect ratio)
   - [ ] Copy is error-free and within character limits
   - [ ] CTA is clear and actionable
   - [ ] Hashtags are relevant and properly researched
   - [ ] Mentions and tags are correct
   - [ ] Alt text provided for accessibility
4. Revision loop: user can request changes, which feed back into Phase 4
5. Batch approval: approve all, approve with exceptions, or reject batch

**Output**: Approved content batch ready for publishing.

---

### Phase 6: Publishing

**Objective**: Distribute approved content to target platforms on schedule.

**Invoke**: `social-media-scheduler` publish workflow

**Process**:
1. Upload assets to platform publishing pipelines
2. Schedule posts according to content calendar timing
3. Support both scheduled and immediate publishing modes
4. Track publish status per post (queued, published, failed)

**Platform Integrations via Composio MCP**:
| Platform   | Integration                | Capabilities                           |
|------------|---------------------------|----------------------------------------|
| Instagram  | `instagram-automation`     | Feed posts, Stories, Reels, Carousels  |
| Facebook   | `facebook`                 | Page posts, Stories, Groups             |
| LinkedIn   | `linkedin-automation`      | Articles, Posts, Document shares        |
| TikTok     | `tiktok-automation`        | Video posts, captions                   |
| Twitter/X  | `twitter` (if available)   | Tweets, threads, polls                  |

**Output**: Publishing log with post IDs, URLs, timestamps, and status per platform.

---

### Phase 7: Analytics & Optimization Loop

**Objective**: Measure performance, extract insights, and optimize future content.

**Invoke**: `social-media-analytics` skill

**Process**:
1. Collect engagement data at defined intervals:
   - **7-day check**: Early performance signals, quick adjustments
   - **14-day review**: Trend identification, mid-campaign pivots
   - **30-day report**: Full campaign analysis and ROI
2. Metrics tracked per post:
   - Reach and impressions
   - Engagement rate (likes, comments, shares, saves)
   - Click-through rate
   - Follower growth attributed to post
   - Conversion events (if tracked)
3. Identify top and bottom performers:
   - Top 20% by engagement rate — analyze what worked
   - Bottom 20% — diagnose issues (timing, format, topic, visual)
4. Generate optimization recommendations:
   - Content type adjustments
   - Posting time refinements
   - Audience targeting shifts
   - Budget reallocation suggestions
5. **Feedback loop**: insights feed directly back into Phase 3 for the next content cycle

**Output**:
- `analytics_report.json` — metrics per post and aggregate
- `campaign_report.md` or `.docx` — final human-readable report
- `optimization_recommendations.json` — actionable changes for next cycle

Use `scripts/campaign_report.py` to generate the final report.

---

## Quick Start Examples

### Example 1: 30-Day Instagram Campaign for a Coffee Brand

```
User: "Create a 30-day Instagram campaign for https://artisancoffee.com"

Orchestrator execution:
1. Phase 1 → brand-identity-miner scrapes artisancoffee.com
   Output: brand_dna.json (earthy tones, artisanal voice, 25-40 demographic)

2. Phase 2 → Strategy: Brand Awareness campaign
   Platforms: Instagram (primary), TikTok (secondary)
   Goals: +1,500 followers, 4% avg engagement rate
   Content pillars: Bean origins, Brewing tutorials, Customer moments

3. Phase 3 → 30 content briefs generated
   Mix: 18 value posts, 9 curated/UGC, 3 promotional
   Calendar: 1 post/day, Stories 3x/week, Reels 2x/week

4. Phase 4 → Copy + visuals produced for all 30 posts
   Style: Warm photography, minimalist text overlays, earthy palette

5. Phase 5 → User reviews batch of 30 posts, approves 27, revises 3

6. Phase 6 → Scheduled via instagram-automation
   Week 1-4 posts queued with optimal timing (7am, 12pm, 6pm)

7. Phase 7 → Day 7 check shows Reels outperform static 3:1
   Adjustment: Increase Reels from 2x to 4x/week in remaining weeks
```

### Example 2: Competitive Analysis + Viral Content Strategy

```
User: "Analyze my competitor https://rivaltech.io and propose viral content for my brand https://mytech.co"

Orchestrator execution:
1. Phase 1 → brand-identity-miner on both URLs
   Output: brand_dna.json (mytech.co) + competitor_analysis.json

2. Phase 2 → Strategy: Engagement Growth (differentiation-focused)
   Identify gaps in competitor content that mytech.co can own

3. Phase 3 → viral-content-creator generates 10 viral content concepts
   Based on: trending formats + competitor gaps + brand strengths
   Each concept includes hook, format, predicted virality score

4. Phase 4 → Produce top 5 concepts as ready-to-post content

5. Phase 5 → User selects final 3 for publication
```

### Example 3: Weekly LinkedIn Posts for a SaaS Product

```
User: "Generate a week of LinkedIn posts for https://cloudsaas.io"

Orchestrator execution (abbreviated pipeline):
1. Phase 1 → brand-identity-miner extracts B2B SaaS identity
   Voice: Professional, data-driven, thought leadership

2. Phase 3 → 5 weekday posts planned
   Mon: Industry insight, Tue: Product tip, Wed: Customer story,
   Thu: Team/culture, Fri: Engagement question

3. Phase 4 → social-media-copywriting generates LinkedIn-optimized copy
   - Hook line (first 2 lines visible before "see more")
   - Value body (3-5 short paragraphs)
   - CTA (comment prompt or link)
   - 3-5 relevant hashtags
   social-media-visual-generator creates supporting graphics

4. Phase 5 → User reviews and approves

5. Phase 6 → Scheduled via linkedin-automation
```

---

## Sub-Skill Integration Map

| Sub-Skill | Responsibility | Input | Output | Path |
|-----------|---------------|-------|--------|------|
| `brand-identity-miner` | Extract brand identity from web presence | URL, social handles | `brand_dna.json` | `scientific-skills/brand-identity-miner/` |
| `viral-content-creator` | Ideate high-engagement content concepts | brand_dna.json, trends | `content_briefs/*.json` | `scientific-skills/viral-content-creator/` |
| `social-media-copywriting` | Write platform-specific copy | content briefs, brand_dna.json | `captions/*.json` | `scientific-skills/social-media-copywriting/` |
| `social-media-visual-generator` | Produce images, carousels, video storyboards | content briefs, brand_dna.json | `visuals/*.png`, `videos/*.mp4` | `scientific-skills/social-media-visual-generator/` |
| `social-media-scheduler` | Calendar generation + publishing | content calendar, assets | `calendar.json`, publish log | `scientific-skills/social-media-scheduler/` |
| `social-media-analytics` | Performance tracking + optimization | publish log, platform APIs | `report.json`, recommendations | `scientific-skills/social-media-analytics/` |

### Data Flow

```
brand_dna.json ─────────────────────────────────────────────────────┐
      │                                                              │
      ├──> campaign_plan.json                                        │
      │         │                                                    │
      │         ├──> content_calendar.json                           │
      │         │         │                                          │
      │         │         ├──> content_briefs/*.json                 │
      │         │         │         │                                │
      │         │         │         ├──> captions/*.json             │
      │         │         │         ├──> visuals/*.png               │
      │         │         │         └──> videos/*.mp4                │
      │         │         │                   │                      │
      │         │         │         [Review & Approval]              │
      │         │         │                   │                      │
      │         │         │         publish_log.json                 │
      │         │         │                   │                      │
      │         │         │         analytics_report.json            │
      │         │         │                   │                      │
      │         │         └───────────────────┘ (optimization loop)  │
      │         │                                                    │
      └─────────┴────────────────────────────────────────────────────┘
                        (brand consistency throughout)
```

---

## MCP Integration Table

All MCP servers and tools used across the orchestrator and its sub-skills:

| MCP Server | Tools | Used In Phase | Purpose |
|------------|-------|---------------|---------|
| `firecrawl-scraper` | `firecrawl_scrape`, `firecrawl_crawl` | Phase 1 | Website content extraction |
| `Claude in Chrome` | `computer`, `read_page`, `navigate`, `javascript_tool` | Phase 1 | Dynamic page interaction, social profile scraping |
| `tavily-web` | `tavily_search` | Phase 1, 2 | Competitor research, trend discovery |
| `fal-generate` | `fal_generate` | Phase 4 | AI image generation (Flux models) |
| `fal-image-edit` | `fal_image_edit` | Phase 4 | Image refinement and editing |
| `fal-audio` | `fal_audio` | Phase 4 | Voiceover generation for video |
| `Gamma` | `mcp__gamma__generate` | Phase 4 | Carousel and presentation creation |
| `Composio` | `instagram-automation`, `facebook`, `linkedin-automation`, `tiktok-automation` | Phase 6 | Platform publishing |
| `Google Gemini` | Image/video generation via API | Phase 4 | Extended visual generation |
| `filesystem` | `read_file`, `write_file`, `list_directory` | All | File management |
| `memory` | `create_entities`, `search_nodes` | All | Campaign state persistence |

---

## Campaign Types Reference

| Campaign Type | Primary Goal | Key KPIs | Typical Duration | Best Platforms |
|--------------|-------------|----------|-----------------|----------------|
| **Product Launch** | Drive awareness and trials for new product | Reach, website clicks, sign-ups | 4-8 weeks | Instagram, TikTok, LinkedIn |
| **Brand Awareness** | Maximize reach among target audience | Impressions, reach, follower growth | 4-12 weeks | Instagram, TikTok, Facebook |
| **Engagement Growth** | Build active community | Engagement rate, comments, shares, saves | Ongoing (monthly cycles) | Instagram, TikTok, Facebook Groups |
| **Lead Generation** | Capture qualified leads | Click-through rate, conversions, cost per lead | 2-6 weeks | LinkedIn, Facebook, Instagram |
| **Event Promotion** | Drive attendance and participation | RSVPs, ticket sales, event mentions | 2-4 weeks pre-event | All platforms |
| **Seasonal/Holiday** | Capitalize on cultural moments | Sales, engagement spikes, share of voice | 1-3 weeks | Instagram, TikTok, Facebook |

---

## Output Directory Structure

When running a full campaign, the orchestrator creates:

```
campaign_output/
  brand_dna.json
  campaign_plan.json
  content_calendar.json
  content_briefs/
    brief_001.json
    brief_002.json
    ...
  captions/
    caption_001.json
    caption_002.json
    ...
  visuals/
    post_001.png
    post_002.png
    ...
  videos/
    reel_001.mp4
    ...
  publish_log.json
  analytics/
    report_7d.json
    report_14d.json
    report_30d.json
  campaign_report.md
  optimization_recommendations.json
```

---

## Usage

### Full Campaign Pipeline
```bash
# Generate complete campaign from brand URL
python scripts/orchestrate_campaign.py \
  --brand-url "https://example.com" \
  --campaign-type "brand_awareness" \
  --platforms "instagram,linkedin" \
  --duration-days 30 \
  --goals "followers:+2000,engagement_rate:4%" \
  --output-dir ./campaign_output

# Generate campaign brief from existing brand DNA
python scripts/brand_brief_generator.py \
  --brand-dna ./campaign_output/brand_dna.json \
  --campaign-type "product_launch" \
  --output ./campaign_output/campaign_brief.json

# Generate final campaign report
python scripts/campaign_report.py \
  --campaign-plan ./campaign_output/campaign_plan.json \
  --analytics ./campaign_output/analytics/report_30d.json \
  --output ./campaign_output/campaign_report.md \
  --format md
```

### Partial Pipeline (specific phases only)
The orchestrator supports running individual phases when the full pipeline is not needed. Each phase can be invoked independently as long as its required inputs exist.
