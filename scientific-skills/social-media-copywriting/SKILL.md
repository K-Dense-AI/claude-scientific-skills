---
name: social-media-copywriting
description: Write platform-specific social media copy optimized for engagement. Handles character limits, tone adaptation, CTA placement, emoji strategy, hashtag integration, and bilingual Romanian/English content. Uses brand_dna.json for voice consistency.
allowed-tools: [Bash, Read, Write, Edit]
license: MIT license
metadata:
    skill-author: K-Dense Inc.
---

# Social Media Copywriting Engine

Platform-optimized copywriting system that generates high-engagement social media content across Instagram, LinkedIn, TikTok, Facebook, and Twitter/X. Supports bilingual Romanian/English output with cultural adaptation.

## When to Use

- Writing captions for any social media platform
- Crafting CTAs optimized for specific goals (engagement, traffic, conversion)
- Creating bilingual RO/EN social content
- Adapting a single message across multiple platforms
- Building hashtag strategies per platform
- Writing ad copy, bios, or comment responses
- Generating carousel or reel captions
- Localizing English content for Romanian audiences (or vice versa)

## Quick Start

### 1. Generate a caption

```bash
python scripts/generate_caption.py \
  --brief '{"topic": "product launch", "key_message": "New organic skincare line", "tone": "excited", "target_audience": "women 25-40"}' \
  --platform instagram \
  --language en \
  --brand-dna brand_dna.json \
  --output caption.json
```

### 2. Get a CTA recommendation

```bash
python scripts/generate_cta.py \
  --goal engage \
  --platform instagram \
  --brand-dna brand_dna.json \
  --output cta.json
```

### 3. Localize content to Romanian

```bash
python scripts/localize_content.py \
  --input caption.json \
  --target-language ro \
  --brand-dna brand_dna.json \
  --output caption_ro.json
```

### 4. Optimize emoji placement

```bash
python scripts/emoji_optimizer.py \
  --input "Check out our new collection" \
  --platform instagram \
  --brand-voice playful \
  --output optimized.json
```

## Brand DNA File Format

Create a `brand_dna.json` file to maintain voice consistency:

```json
{
  "brand_name": "YourBrand",
  "voice": {
    "tone": ["friendly", "expert", "approachable"],
    "personality": "A knowledgeable friend who makes complex topics simple",
    "formality": "casual-professional",
    "humor_level": "moderate",
    "emoji_density": "medium"
  },
  "audience": {
    "primary": "women 25-40",
    "secondary": "health-conscious millennials",
    "pain_points": ["lack of time", "information overload", "budget concerns"],
    "aspirations": ["wellness", "self-improvement", "work-life balance"]
  },
  "language": {
    "primary": "en",
    "secondary": "ro",
    "ro_formality": "tu",
    "banned_words": ["cheap", "buy now", "limited time"],
    "preferred_words": ["discover", "transform", "elevate"]
  },
  "hashtags": {
    "branded": ["#YourBrand", "#YourBrandLife"],
    "always_include": ["#YourBrand"],
    "never_use": ["#ad", "#sponsored"]
  }
}
```

## Platform-Specific Guidelines

### Instagram
- **Caption sweet spot:** 150-300 characters for feed posts, up to 2200 max
- **First line is everything:** Only ~125 characters show before "more"
- **Hashtags:** 5-15 in a comment or after line breaks (30 max)
- **Emoji:** Liberal use, 1-3 per sentence acceptable
- **Tone:** Visual, aspirational, personal
- **Best formats:** Storytelling, listicles, question hooks

### LinkedIn
- **Caption sweet spot:** 150-300 characters for maximum engagement
- **First line hook:** Only ~140 characters before "see more"
- **Hashtags:** 3-5 max, placed at the end
- **Emoji:** Minimal, professional (checkmarks, arrows, pointing hands)
- **Tone:** Professional, insightful, thought-leadership
- **Best formats:** Hot takes, lessons learned, data-driven insights

### TikTok
- **Caption limit:** 4000 characters (but shorter is better)
- **Optimal length:** Under 150 characters
- **Hashtags:** 3-5, mix trending + niche
- **Emoji:** Moderate, trend-dependent
- **Tone:** Casual, authentic, trend-aware
- **Best formats:** Hook-driven, conversational, challenge-oriented

### Facebook
- **Caption sweet spot:** 80-150 characters for link posts, up to 500 for stories
- **Truncation:** Around 477 characters before "See More"
- **Hashtags:** 1-3 max (less important here)
- **Emoji:** Moderate use
- **Tone:** Conversational, community-oriented
- **Best formats:** Questions, personal stories, community polls

### Twitter/X
- **Character limit:** 280 (4000 for premium)
- **Optimal length:** 71-100 characters for highest engagement
- **Hashtags:** 1-2 max
- **Emoji:** Strategic, 1-2 per tweet
- **Tone:** Witty, concise, opinionated
- **Best formats:** Hot takes, threads, quote tweets

## Copywriting Frameworks

### AIDA (Attention-Interest-Desire-Action)
Best for: Product launches, promotional posts

```
[ATTENTION] Bold hook statement
[INTEREST] Expand with benefits/features
[DESIRE] Social proof or emotional appeal
[ACTION] Clear CTA
```

### PAS (Problem-Agitation-Solution)
Best for: Educational content, service promotion

```
[PROBLEM] Name the pain point
[AGITATION] Amplify the frustration
[SOLUTION] Present your offer as the answer
```

### BAB (Before-After-Bridge)
Best for: Testimonials, transformation stories

```
[BEFORE] Describe the current state
[AFTER] Paint the desired outcome
[BRIDGE] Show how you get them there
```

### 4U (Useful-Urgent-Unique-Ultra-specific)
Best for: Twitter/X, short-form content

### SPIN (Situation-Problem-Implication-Need-payoff)
Best for: LinkedIn thought leadership, B2B content

## Caption Structure

Every high-performing caption follows this skeleton:

```
1. HOOK (first line, visible before fold)
2. BODY (value delivery, 2-5 sentences)
3. CTA (single clear action)
4. HASHTAGS (platform-appropriate count)
```

## Hashtag Strategy

| Platform   | Count | Placement              | Type Mix                         |
|------------|-------|------------------------|----------------------------------|
| Instagram  | 5-15  | First comment or below  | 30% branded, 40% niche, 30% broad |
| LinkedIn   | 3-5   | End of caption          | 50% industry, 50% topic          |
| TikTok     | 3-5   | In caption              | 50% trending, 50% niche          |
| Facebook   | 1-3   | In caption              | Branded only                     |
| Twitter/X  | 1-2   | In tweet                | Trending or branded               |

## Emoji Usage Guidelines

| Brand Voice    | Emoji Density | Style                                    |
|----------------|---------------|------------------------------------------|
| Corporate      | Low (0-1)     | Professional only (checkmarks, arrows)   |
| Professional   | Low-Med (1-2) | Subtle emphasis, no faces                |
| Friendly       | Medium (2-4)  | Warm, relatable, occasional faces        |
| Playful        | High (3-6)    | Expressive, faces, objects, creative     |
| Gen-Z/Trendy   | High (4-8)    | Ironic, layered, trend-specific          |

## Bilingual RO/EN Support

This skill supports full bilingual content creation:

- **Parallel content:** Generate the same post in both languages
- **Cultural adaptation:** Romanian posts are culturally adapted, not literally translated
- **Tu/Dumneavoastra:** Formality level set in brand_dna.json
- **Mixed-language posts:** Support for RO caption with EN hashtags (common pattern)

See `references/multilingual_guide.md` for Romanian-specific conventions.

## Reference Files

| File | Contents |
|------|----------|
| `references/platform_copy_specs.md` | Character limits, formatting, best practices per platform |
| `references/emotional_triggers.md` | Copywriting frameworks, power words, persuasion psychology |
| `references/cta_library.md` | 100+ CTAs organized by goal and platform |
| `references/caption_templates.md` | Fill-in-the-blank templates in EN and RO |
| `references/multilingual_guide.md` | Romanian social media conventions and localization guide |

## Scripts

| Script | Purpose |
|--------|---------|
| `scripts/generate_caption.py` | Generate platform-optimized captions |
| `scripts/generate_cta.py` | Select and customize CTAs |
| `scripts/localize_content.py` | Culturally adapt content between EN and RO |
| `scripts/emoji_optimizer.py` | Add strategic emoji placement |
