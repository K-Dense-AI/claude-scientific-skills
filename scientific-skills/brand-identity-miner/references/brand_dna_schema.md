# Brand DNA JSON Schema

Complete specification for the `brand_dna.json` output file. This document serves as the contract between the Brand Identity Miner skill and all downstream content creation skills.

## Schema Version

```
schema_version: "1.0.0"
```

## Full Schema

```json
{
  "schema_version": "1.0.0",
  "generated_at": "ISO-8601 timestamp",
  "generated_by": "brand-identity-miner",

  "company_name": "",
  "tagline": "",
  "website_url": "",
  "industry": "",
  "brand_archetype": "",

  "visual_identity": {
    "primary_colors": ["#hex"],
    "secondary_colors": ["#hex"],
    "accent_colors": ["#hex"],
    "fonts": {
      "heading": "",
      "body": "",
      "accent": ""
    },
    "logo_description": "",
    "photography_style": "",
    "visual_mood": ""
  },

  "voice": {
    "tone": [],
    "personality_traits": [],
    "vocabulary_level": "technical|casual|mixed",
    "formality": "formal|semiformal|casual",
    "voice_dimensions": {
      "formal_casual": 5,
      "serious_playful": 5,
      "technical_simple": 5,
      "authoritative_friendly": 5
    },
    "do_say": [],
    "dont_say": [],
    "sample_phrases": []
  },

  "values": [],
  "mission": "",
  "unique_selling_propositions": [],

  "target_audience": {
    "demographics": {
      "age_range": "",
      "gender": "",
      "location": "",
      "income_level": "",
      "education": "",
      "occupation": ""
    },
    "psychographics": {
      "interests": [],
      "values": [],
      "lifestyle": "",
      "media_consumption": []
    },
    "pain_points": [],
    "aspirations": []
  },

  "competitors": [
    {
      "name": "",
      "url": "",
      "positioning": "",
      "strengths": [],
      "weaknesses": [],
      "social_profiles": {},
      "estimated_sov": ""
    }
  ],

  "differentiators": [],
  "content_pillars": [
    {
      "name": "",
      "description": "",
      "audience_need": "",
      "example_topics": [],
      "formats": [],
      "frequency": ""
    }
  ],

  "social_profiles": {
    "instagram": "",
    "facebook": "",
    "linkedin": "",
    "tiktok": "",
    "x_twitter": "",
    "youtube": ""
  }
}
```

## Field Descriptions

### Top-Level Fields

| Field              | Type      | Required | Description                                                              |
|--------------------|-----------|----------|--------------------------------------------------------------------------|
| `schema_version`   | string    | Yes      | Schema version for forward compatibility. Currently `"1.0.0"`.           |
| `generated_at`     | string    | Yes      | ISO-8601 timestamp of when the brand DNA was generated.                  |
| `generated_by`     | string    | Yes      | Always `"brand-identity-miner"`.                                         |
| `company_name`     | string    | Yes      | Legal or commonly known name of the company.                             |
| `tagline`          | string    | Yes      | Primary brand tagline or slogan. Extracted from website hero section.    |
| `website_url`      | string    | Yes      | Primary website URL.                                                     |
| `industry`         | string    | Yes      | Industry or vertical (e.g., "SaaS", "Health & Wellness", "E-commerce"). |
| `brand_archetype`  | string    | Yes      | One of the 12 Jungian archetypes. See `brand_analysis_framework.md`.     |

### visual_identity

Visual brand elements extracted from website and social profiles.

| Field                | Type     | Required | Description                                                         |
|----------------------|----------|----------|---------------------------------------------------------------------|
| `primary_colors`     | string[] | Yes      | 2-3 hex color codes. Dominant brand colors from logo, headers, CTAs.|
| `secondary_colors`   | string[] | Yes      | 2-4 hex codes. Supporting background and UI colors.                 |
| `accent_colors`      | string[] | No       | 1-2 hex codes. High-contrast emphasis colors.                       |
| `fonts.heading`      | string   | Yes      | Typeface name for headings (e.g., "Montserrat", "Playfair Display").|
| `fonts.body`         | string   | Yes      | Typeface name for body text (e.g., "Inter", "Open Sans").           |
| `fonts.accent`       | string   | No       | Optional accent or display typeface.                                |
| `logo_description`   | string   | Yes      | Text description of the logo (shape, color, wordmark vs. icon).     |
| `photography_style`  | string   | Yes      | One of: lifestyle, product-focused, abstract, people-centric, illustration, data-visualization, stock-minimal, ugc-style. |
| `visual_mood`        | string   | Yes      | 2-4 word description of the overall aesthetic (e.g., "minimal and modern", "warm and artisan"). |

### voice

Brand communication style and language patterns.

| Field                             | Type     | Required | Description                                                      |
|-----------------------------------|----------|----------|------------------------------------------------------------------|
| `tone`                            | string[] | Yes      | 3-5 tone descriptors (e.g., ["confident", "warm", "direct"]).    |
| `personality_traits`              | string[] | Yes      | 3-5 adjectives describing brand personality.                     |
| `vocabulary_level`                | string   | Yes      | One of: `"technical"`, `"casual"`, `"mixed"`.                    |
| `formality`                       | string   | Yes      | One of: `"formal"`, `"semiformal"`, `"casual"`.                  |
| `voice_dimensions.formal_casual`  | integer  | Yes      | 1-10 scale. 1 = very formal, 10 = very casual.                  |
| `voice_dimensions.serious_playful`| integer  | Yes      | 1-10 scale. 1 = very serious, 10 = very playful.                |
| `voice_dimensions.technical_simple`| integer | Yes      | 1-10 scale. 1 = highly technical, 10 = very simple.             |
| `voice_dimensions.authoritative_friendly`| integer | Yes | 1-10 scale. 1 = very authoritative, 10 = very friendly.         |
| `do_say`                          | string[] | Yes      | 5-10 phrases or language patterns the brand should use.          |
| `dont_say`                        | string[] | Yes      | 5-10 phrases or language patterns the brand should avoid.        |
| `sample_phrases`                  | string[] | Yes      | 5-10 representative phrases extracted from existing brand copy.  |

### values

| Field    | Type     | Required | Description                                           |
|----------|----------|----------|-------------------------------------------------------|
| `values` | string[] | Yes      | 3-7 core brand values (e.g., ["innovation", "transparency", "sustainability"]). |

### mission

| Field     | Type   | Required | Description                                |
|-----------|--------|----------|--------------------------------------------|
| `mission` | string | Yes      | 1-2 sentence brand mission statement.      |

### unique_selling_propositions

| Field                        | Type     | Required | Description                                          |
|------------------------------|----------|----------|------------------------------------------------------|
| `unique_selling_propositions`| string[] | Yes      | 3-5 specific claims about what makes the brand unique. |

### target_audience

| Field                              | Type     | Required | Description                                          |
|------------------------------------|----------|----------|------------------------------------------------------|
| `demographics.age_range`           | string   | Yes      | Target age range (e.g., "25-45").                    |
| `demographics.gender`              | string   | No       | Target gender if applicable, or "all".               |
| `demographics.location`            | string   | Yes      | Geographic focus (e.g., "US", "Global", "Urban EU"). |
| `demographics.income_level`        | string   | No       | Income bracket (e.g., "middle", "upper-middle").     |
| `demographics.education`           | string   | No       | Education level (e.g., "college-educated").           |
| `demographics.occupation`          | string   | No       | Professional context (e.g., "tech professionals").   |
| `psychographics.interests`         | string[] | Yes      | 5-10 interest areas.                                 |
| `psychographics.values`            | string[] | Yes      | 3-5 personal values of the target audience.          |
| `psychographics.lifestyle`         | string   | Yes      | Brief lifestyle description.                         |
| `psychographics.media_consumption` | string[] | No       | Platforms and media the audience consumes.            |
| `pain_points`                      | string[] | Yes      | 3-7 problems the audience faces that the brand solves.|
| `aspirations`                      | string[] | Yes      | 3-5 goals or desires the audience has.               |

### competitors

Array of 3-5 competitor profiles.

| Field               | Type     | Required | Description                                            |
|---------------------|----------|----------|--------------------------------------------------------|
| `name`              | string   | Yes      | Competitor company name.                               |
| `url`               | string   | Yes      | Competitor website URL.                                |
| `positioning`       | string   | Yes      | Brief positioning statement (e.g., "premium B2B SaaS").|
| `strengths`         | string[] | Yes      | 3-5 competitive strengths (social media focused).      |
| `weaknesses`        | string[] | Yes      | 3-5 competitive weaknesses (social media focused).     |
| `social_profiles`   | object   | No       | Map of platform to handle/URL.                         |
| `estimated_sov`     | string   | No       | Estimated share of voice percentage.                   |

### differentiators

| Field             | Type     | Required | Description                                                |
|-------------------|----------|----------|------------------------------------------------------------|
| `differentiators` | string[] | Yes      | 3-5 specific ways the brand differs from listed competitors.|

### content_pillars

Array of 3-5 content pillar definitions.

| Field            | Type     | Required | Description                                              |
|------------------|----------|----------|----------------------------------------------------------|
| `name`           | string   | Yes      | Short pillar name (e.g., "Industry Insights").           |
| `description`    | string   | Yes      | 1-2 sentence description of the content theme.           |
| `audience_need`  | string   | Yes      | Which audience need this pillar addresses.               |
| `example_topics` | string[] | Yes      | 3-5 specific topic ideas.                                |
| `formats`        | string[] | Yes      | Recommended content formats (carousel, video, blog, etc.)|
| `frequency`      | string   | Yes      | Recommended posting frequency for this pillar.           |

### social_profiles

| Field       | Type   | Required | Description                   |
|-------------|--------|----------|-------------------------------|
| `instagram` | string | No       | Instagram handle or URL.      |
| `facebook`  | string | No       | Facebook page URL.            |
| `linkedin`  | string | No       | LinkedIn company page URL.    |
| `tiktok`    | string | No       | TikTok handle or URL.         |
| `x_twitter` | string | No       | X (Twitter) handle or URL.    |
| `youtube`   | string | No       | YouTube channel URL.          |

## Validation Rules

1. All hex color codes must be valid 6-digit hex with `#` prefix (e.g., `#1A2B3C`).
2. `brand_archetype` must be one of: Innocent, Explorer, Sage, Hero, Outlaw, Magician, Regular Guy, Lover, Jester, Caregiver, Creator, Ruler.
3. `vocabulary_level` must be one of: technical, casual, mixed.
4. `formality` must be one of: formal, semiformal, casual.
5. Voice dimension scores must be integers from 1 to 10.
6. `content_pillars` array must contain 3-5 entries.
7. `competitors` array should contain 3-5 entries for a complete analysis.
8. `values` array should contain 3-7 entries.
9. All URL fields must be valid URLs starting with `http://` or `https://`.

## Example

```json
{
  "schema_version": "1.0.0",
  "generated_at": "2026-03-18T10:30:00Z",
  "generated_by": "brand-identity-miner",

  "company_name": "Greenleaf Analytics",
  "tagline": "Data-driven sustainability for modern businesses",
  "website_url": "https://greenleafanalytics.com",
  "industry": "SaaS / Sustainability Tech",
  "brand_archetype": "Sage",

  "visual_identity": {
    "primary_colors": ["#2D6A4F", "#1B4332"],
    "secondary_colors": ["#B7E4C7", "#D8F3DC", "#F5F5F5"],
    "accent_colors": ["#FF6B35"],
    "fonts": {
      "heading": "DM Sans",
      "body": "Inter",
      "accent": ""
    },
    "logo_description": "Stylized leaf icon composed of data points, paired with a clean sans-serif wordmark in dark green",
    "photography_style": "data-visualization",
    "visual_mood": "clean, professional, nature-tech hybrid"
  },

  "voice": {
    "tone": ["knowledgeable", "optimistic", "direct", "approachable"],
    "personality_traits": ["analytical", "forward-thinking", "trustworthy", "pragmatic"],
    "vocabulary_level": "mixed",
    "formality": "semiformal",
    "voice_dimensions": {
      "formal_casual": 4,
      "serious_playful": 3,
      "technical_simple": 4,
      "authoritative_friendly": 6
    },
    "do_say": [
      "data shows",
      "measurable impact",
      "let's look at the numbers",
      "sustainable growth",
      "actionable insights"
    ],
    "dont_say": [
      "trust us",
      "revolutionary",
      "game-changer",
      "synergy",
      "at the end of the day"
    ],
    "sample_phrases": [
      "Your sustainability data, finally making sense.",
      "We believe what gets measured gets improved.",
      "Built for teams who take impact seriously."
    ]
  },

  "values": ["transparency", "sustainability", "data-integrity", "pragmatism", "continuous improvement"],
  "mission": "Empower businesses to measure, understand, and reduce their environmental impact through accessible data analytics.",
  "unique_selling_propositions": [
    "Only platform that integrates supply chain and operational carbon data in one dashboard",
    "Automated ESG reporting that saves 40+ hours per quarter",
    "Built by sustainability scientists, not just engineers"
  ],

  "target_audience": {
    "demographics": {
      "age_range": "30-50",
      "gender": "all",
      "location": "North America, Western Europe",
      "income_level": "upper-middle",
      "education": "college-educated, often MBA or STEM background",
      "occupation": "sustainability managers, CSOs, operations directors"
    },
    "psychographics": {
      "interests": ["climate tech", "ESG investing", "operational efficiency", "corporate responsibility"],
      "values": ["environmental stewardship", "evidence-based decisions", "professional growth"],
      "lifestyle": "busy professionals balancing business performance with sustainability goals",
      "media_consumption": ["LinkedIn", "industry newsletters", "sustainability podcasts", "Harvard Business Review"]
    },
    "pain_points": [
      "Scattered sustainability data across spreadsheets and systems",
      "Difficulty translating environmental metrics into business language",
      "Increasing regulatory pressure for ESG reporting",
      "Lack of benchmarks to compare sustainability performance"
    ],
    "aspirations": [
      "Become a recognized leader in corporate sustainability",
      "Reduce environmental footprint while maintaining profitability",
      "Simplify reporting to save time and reduce errors"
    ]
  },

  "competitors": [
    {
      "name": "EcoTrack Pro",
      "url": "https://ecotrackpro.com",
      "positioning": "Enterprise ESG platform for Fortune 500",
      "strengths": ["Strong enterprise features", "Regulatory compliance depth"],
      "weaknesses": ["Complex onboarding", "Weak social media presence", "No SMB offering"],
      "social_profiles": {
        "linkedin": "/company/ecotrackpro",
        "instagram": "@ecotrackpro"
      },
      "estimated_sov": "35%"
    }
  ],

  "differentiators": [
    "Mid-market focus vs. enterprise-only competitors",
    "Scientist-founded credibility",
    "Supply chain integration unique in the market",
    "Self-serve onboarding (no 6-month implementation)"
  ],

  "content_pillars": [
    {
      "name": "Data Decoded",
      "description": "Breaking down complex sustainability metrics into understandable, actionable insights",
      "audience_need": "Understanding what their environmental data actually means",
      "example_topics": ["Carbon accounting basics", "Scope 3 demystified", "ESG metrics that matter"],
      "formats": ["carousel", "infographic", "linkedin-article"],
      "frequency": "2x per week"
    },
    {
      "name": "Regulation Radar",
      "description": "Timely updates on ESG regulations and what they mean for businesses",
      "audience_need": "Staying compliant with evolving requirements",
      "example_topics": ["EU CSRD breakdown", "SEC climate disclosure rules", "ISSB standards"],
      "formats": ["short-video", "text-post", "newsletter"],
      "frequency": "1x per week"
    },
    {
      "name": "Impact Stories",
      "description": "Real customer case studies showing measurable sustainability improvements",
      "audience_need": "Proof that sustainability efforts lead to results",
      "example_topics": ["Customer reduced emissions 30%", "ROI of sustainability reporting"],
      "formats": ["video-testimonial", "carousel", "blog"],
      "frequency": "1x per week"
    }
  ],

  "social_profiles": {
    "instagram": "@greenleafanalytics",
    "facebook": "",
    "linkedin": "/company/greenleaf-analytics",
    "tiktok": "",
    "x_twitter": "@greenleafdata",
    "youtube": ""
  }
}
```
