---
name: social-media-scheduler
description: Create content calendars, determine optimal posting schedules, manage A/B testing, and handle publishing to Facebook, Instagram, LinkedIn, and TikTok via Composio MCP integrations. Works from campaign plans produced by social-media-orchestrator.
allowed-tools: [Bash, Read, Write, Edit, WebSearch, WebFetch, Agent]
license: MIT license
metadata:
    skill-author: K-Dense Inc.
---

# Social Media Scheduler

Scheduling and publishing engine for multi-platform social media campaigns. Transforms campaign plans into actionable content calendars, identifies optimal posting windows, manages A/B testing workflows, and handles publishing through Composio MCP integrations.

## When to Use

- Creating weekly, monthly, or campaign-specific content calendars
- Determining the best times to post on each platform for a given industry
- Setting up A/B tests to optimize hooks, visuals, CTAs, posting times, or formats
- Publishing or scheduling content to Facebook, Instagram, LinkedIn, and TikTok
- Exporting calendars in JSON, CSV, or Markdown table format

## Quick Start

### Generate a Content Calendar

```bash
python scripts/generate_calendar.py \
  --campaign-plan campaign.json \
  --start-date 2026-04-01 \
  --duration-days 30 \
  --platforms instagram,linkedin,tiktok \
  --posts-per-week '{"instagram": 5, "linkedin": 3, "tiktok": 4}' \
  --output calendar_april.json
```

### Get Optimal Posting Times

```bash
python scripts/suggest_posting_times.py \
  --platform instagram \
  --industry b2c \
  --timezone "America/New_York" \
  --output posting_windows.json
```

### Plan an A/B Test

```bash
python scripts/ab_test_planner.py \
  --variable hook \
  --base-content base_post.json \
  --num-variants 3 \
  --duration-days 7 \
  --output ab_test_plan.json
```

### Publish Content

```bash
python scripts/publish_content.py \
  --content-file post.json \
  --platform instagram \
  --schedule "2026-04-01T10:00:00-04:00" \
  --output publish_log.json
```

## Content Calendar Generation

The calendar generator takes a campaign plan (JSON) and distributes content across platforms and days according to configurable rules.

### Workflow

1. **Load campaign plan** -- reads content pillars, formats, platform targets, and campaign phases from JSON.
2. **Assign posting slots** -- distributes posts across the date range respecting per-platform frequency limits.
3. **Apply content mix rules** -- ensures pillar rotation (no two consecutive posts from the same pillar), format variety (carousel, reel, static, story rotated), and CTA diversity.
4. **Inject optimal times** -- assigns posting times from the research-backed optimal windows for each platform and industry.
5. **Export** -- outputs the calendar as JSON and a Markdown table.

### Campaign Plan JSON Schema

```json
{
  "campaign_name": "Spring Launch",
  "brand": "Acme Co",
  "goal": "Drive sign-ups for new product",
  "pillars": ["education", "social-proof", "behind-the-scenes", "promotion"],
  "formats": {
    "instagram": ["reel", "carousel", "story", "static"],
    "linkedin": ["text", "carousel", "video", "document"],
    "tiktok": ["short-video", "duet", "stitch"],
    "facebook": ["video", "image", "link", "carousel"]
  },
  "phases": [
    {"name": "tease", "days": 5},
    {"name": "launch", "days": 3},
    {"name": "sustain", "days": 18},
    {"name": "recap", "days": 4}
  ]
}
```

### Calendar Output Fields

Each calendar entry contains: `date`, `day_of_week`, `platform`, `format`, `content_pillar`, `hook_placeholder`, `visual_type`, `cta`, `posting_time`, `status` (draft/scheduled/published), and `phase`.

## Optimal Posting Windows

Research-backed posting times are stored in `references/optimal_posting_times.md` and consumed by `scripts/suggest_posting_times.py`. Data is sourced from Sprout Social (2024-2025), Hootsuite (2024-2025), Buffer (2024), and Later (2024) research reports.

Key factors considered:

- **Platform-specific peak engagement hours** ranked by day of week
- **Industry vertical adjustments** for B2B, B2C, e-commerce, SaaS, fitness, food, and fashion
- **Weekend vs weekday patterns** per platform
- **Holiday and event considerations** that shift typical patterns
- **Frequency guardrails** (minimum and maximum posts per day and per week)
- **Timezone normalization** so all recommendations output in the user's local time

## A/B Testing Framework

The A/B testing system follows a structured methodology documented in `references/ab_testing_guide.md`.

### Testable Variables

| Variable | What Changes | What Stays Constant |
|----------|-------------|-------------------|
| Hook | Opening line / first 3 seconds | Visual, CTA, posting time |
| Visual | Image style, thumbnail, colors | Caption, CTA, posting time |
| CTA | Call-to-action text or placement | Hook, visual, posting time |
| Posting time | Hour and day of week | Caption, visual, CTA |
| Format | Reel vs carousel vs static | Topic, messaging, posting time |
| Caption length | Short vs medium vs long | Hook, visual, CTA |
| Hashtag set | Different hashtag groups | Caption, visual, posting time |

### Test Workflow

1. Define the variable to test and the hypothesis.
2. Create variants (minimum 2, recommended 3) using `ab_test_planner.py`.
3. Run variants over a minimum of 7 days with equal distribution.
4. Measure the primary metric (engagement rate, click-through, saves, shares).
5. Determine statistical significance (aim for 95% confidence).
6. Roll out the winner and document findings.

## Publishing via Composio MCP

Publishing is handled through Composio MCP tool integrations. The `publish_content.py` script constructs the appropriate tool calls for each platform.

### Supported Platforms and Composio Tools

| Platform | Composio MCP Integration | Capabilities |
|----------|------------------------|-------------|
| Instagram | `instagram-automation` | Feed posts, reels, stories, carousels |
| Facebook | `facebook` | Page posts, stories, reels, link shares |
| LinkedIn | `linkedin-automation` | Text posts, image posts, document shares, articles |
| TikTok | `tiktok-automation` | Video uploads, captions, hashtags |

### Content File Schema (for publish_content.py)

```json
{
  "platform": "instagram",
  "post_type": "carousel",
  "caption": "Your caption text here with #hashtags",
  "media_paths": ["/path/to/slide1.png", "/path/to/slide2.png"],
  "alt_text": "Descriptive alt text for accessibility",
  "schedule_time": "2026-04-01T10:00:00-04:00",
  "first_comment": "Additional hashtags or engagement prompt",
  "location_tag": "New York, NY"
}
```

### Publishing Modes

- **Immediate** (`--schedule now`): Publishes the content right away through the Composio MCP tool.
- **Scheduled** (`--schedule <ISO-datetime>`): Schedules the post for the specified time. The script stores the scheduled job and confirms the platform accepted it.

### Publish Log

Every publish action is logged with: `post_id`, `platform`, `timestamp_utc`, `schedule_time`, `status` (success/failed/scheduled), `response_data`, and `error_message` (if any).

## Calendar Export Formats

### JSON

Full structured output with all fields. Used as input by other scripts and for programmatic processing.

### CSV

Flat table suitable for import into Google Sheets, Excel, or Notion databases. Columns match the calendar output fields.

### Markdown Table

Human-readable table for embedding in campaign briefs, Notion pages, or Slack messages. Grouped by week with platform emoji indicators.

## Reference Files

- `references/optimal_posting_times.md` -- Research-backed posting times by platform, day, hour, and industry
- `references/content_calendar_templates.md` -- Reusable calendar templates for weekly, monthly, campaign, product launch, event, and evergreen schedules
- `references/ab_testing_guide.md` -- Complete A/B testing methodology, naming conventions, measurement frameworks, and interpretation guide
