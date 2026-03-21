---
name: viral-content-creator
description: Create viral-optimized content strategies and detailed content briefs for all social media platforms. Features 150+ viral hooks, trending format templates, engagement psychology frameworks, and platform-specific content planning. Produces structured briefs that feed into visual-generator and copywriting skills.
allowed-tools: [Bash, Read, Write, Edit, WebSearch, WebFetch, Agent]
license: MIT license
metadata:
    skill-author: K-Dense Inc.
---

# Viral Content Creator

## Overview

This skill is a **content strategy engine** that produces structured content briefs, viral angle suggestions, hashtag strategies, and balanced content calendars. It does NOT produce the actual visual or written content -- instead, it generates detailed briefs that feed into the `visual-generator` and `copywriting` skills for execution.

The engine draws on 150+ proven viral hooks organized by psychological trigger, 20+ trending content formats, platform-specific technical specs, engagement psychology frameworks, and battle-tested content structures to transform any topic into a high-engagement content plan.

## When to Use

- **Content ideation**: Generate viral angles and hooks for any topic or brand
- **Content calendar planning**: Build balanced weekly/monthly content mixes
- **Platform adaptation**: Translate one idea into platform-specific briefs
- **Trend research**: Identify trending formats and adapt them to a brand
- **Campaign planning**: Structure multi-post campaigns with narrative arcs
- **Hashtag strategy**: Build tiered hashtag sets for maximum discoverability

## Quick Start

### Generate a Content Brief

```bash
python scripts/generate_content_brief.py \
  --brand-dna "SaaS productivity tool for remote teams" \
  --topic "async communication best practices" \
  --platform instagram \
  --format carousel \
  --output brief.json
```

### Suggest Viral Angles

```bash
python scripts/suggest_viral_angles.py \
  --brand-dna "Sustainable fashion brand targeting Gen Z" \
  --topic "fast fashion environmental impact" \
  --num-angles 7 \
  --output angles.json
```

### Build a Hashtag Strategy

```bash
python scripts/hashtag_optimizer.py \
  --topic "home workout routines" \
  --platform tiktok \
  --brand-dna "Fitness app for busy professionals" \
  --num-hashtags 25 \
  --output hashtags.json
```

### Generate a Content Calendar

```bash
python scripts/content_mixer.py \
  --brand-dna "B2B data analytics platform" \
  --duration-days 30 \
  --platforms instagram,linkedin,tiktok \
  --posts-per-week 5 \
  --output calendar.json
```

## Content Brief Schema

Every brief produced by this skill follows a standard schema:

```json
{
  "brief_id": "uuid",
  "created_at": "ISO-8601 timestamp",
  "brand_context": {
    "brand_dna": "string",
    "content_pillar": "string",
    "target_audience": "string"
  },
  "content_spec": {
    "platform": "instagram|tiktok|linkedin|twitter|youtube|facebook",
    "format": "post|carousel|reel|story|quiz|thread|short",
    "hook": "string - the opening line or concept",
    "angle": "string - the strategic approach",
    "key_messages": ["array of 2-4 core messages"],
    "visual_direction": "string - guidance for visual-generator",
    "caption_notes": "string - guidance for copywriting skill",
    "cta": "string - call to action",
    "target_emotion": "string - primary emotion to evoke"
  },
  "distribution": {
    "hashtags": {
      "high_volume": ["reach hashtags"],
      "medium": ["discoverability hashtags"],
      "niche": ["targeted hashtags"],
      "branded": ["brand-specific hashtags"]
    },
    "posting_time": "suggested time window",
    "cross_post_adaptations": {}
  },
  "metadata": {
    "content_mix_category": "value|curated|promotional",
    "estimated_engagement": "low|medium|high|viral_potential",
    "trend_alignment": "string or null"
  }
}
```

## The 60/30/10 Content Mix Rule

All content calendars follow this proven ratio:

| Category | Percentage | Purpose | Examples |
|----------|-----------|---------|----------|
| **Value** | 60% | Educate, entertain, inspire | Tutorials, tips, behind-the-scenes, storytelling |
| **Curated / Shared** | 30% | Build community, show awareness | Industry news, UGC reposts, collaborations, polls |
| **Promotional** | 10% | Drive conversions | Product launches, offers, testimonials, case studies |

This ratio keeps the audience engaged without feeling sold to. The content mixer script enforces this distribution automatically.

## Content Pillar Rotation

Every brand should define 3-5 content pillars that rotate through their calendar:

1. **Authority** -- Demonstrate expertise (how-tos, frameworks, data)
2. **Relatability** -- Show the human side (behind-the-scenes, mistakes, day-in-the-life)
3. **Community** -- Foster belonging (polls, Q&A, UGC, challenges)
4. **Aspiration** -- Paint the future (transformations, case studies, vision)
5. **Entertainment** -- Create joy (trends, humor, unexpected content)

The content mixer rotates through these pillars to prevent content fatigue and build a multi-dimensional brand presence.

## Trend Research Workflow

Use the tavily-web MCP or WebSearch tool to research current trends:

1. **Identify trending topics** in the brand's niche
2. **Analyze top-performing content** on each target platform
3. **Map trends to brand pillars** -- only pursue trends that align
4. **Check trend lifecycle stage** -- early = high reward, late = low differentiation
5. **Adapt, don't copy** -- put a unique brand spin on every trend

The `suggest_viral_angles.py` script automates steps 3-5 by combining trending formats from `references/trending_formats.md` with the brand DNA to produce original angle suggestions.

## Platform-Specific Content Adaptation

One idea should become multiple platform-native posts:

| Source Idea | Instagram | TikTok | LinkedIn | Twitter/X |
|-------------|-----------|--------|----------|-----------|
| Tutorial | 10-slide carousel | 60s step-by-step reel | Long-form text post with screenshots | Thread with key steps |
| Hot take | Bold text graphic + caption | Green screen reaction | "Broetry" post | Single punchy tweet + poll |
| Case study | Before/after carousel | Transformation video | Data-driven article post | Key stat + thread |
| Behind-the-scenes | Photo dump | Day-in-the-life | Culture/team spotlight | Casual observation tweet |

See `references/platform_formats.md` for exact specs and `references/trending_formats.md` for format-specific templates.

## Integration with Other Skills

### visual-generator

Pass the `visual_direction` field from any content brief to the visual-generator skill to produce platform-ready graphics, carousel slides, or video storyboards.

### copywriting

Pass the `caption_notes`, `hook`, `key_messages`, and `cta` fields to the copywriting skill to produce platform-native captions, threads, or scripts.

### Workflow

```
viral-content-creator (strategy + brief)
    |
    +---> visual-generator (graphics + video)
    |
    +---> copywriting (captions + scripts)
    |
    +---> [manual review + publish]
```

## Reference Files

| File | Contents |
|------|----------|
| `references/viral_hooks_library.md` | 150+ viral hooks organized by psychology type |
| `references/platform_formats.md` | Technical specs for every platform and format |
| `references/trending_formats.md` | 20+ trending content formats with templates |
| `references/engagement_psychology.md` | Psychology frameworks for viral content |
| `references/content_frameworks.md` | Proven content structures and story frameworks |

## Scripts

| Script | Purpose |
|--------|---------|
| `scripts/generate_content_brief.py` | Produce a structured content brief JSON |
| `scripts/suggest_viral_angles.py` | Generate viral angle suggestions for any topic |
| `scripts/hashtag_optimizer.py` | Build tiered hashtag strategies |
| `scripts/content_mixer.py` | Generate balanced content calendars |
