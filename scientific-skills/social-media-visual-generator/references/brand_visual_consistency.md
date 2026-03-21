# Brand Visual Consistency

Guidelines for maintaining brand identity across all AI-generated social media visuals.

---

## Encoding Brand Colors in AI Image Prompts

AI image generators interpret color instructions through natural language. Direct hex codes are generally ignored by diffusion models; instead, describe colors descriptively and reinforce them through scene context.

### Effective Color Prompting Strategies

**Name the colors explicitly:**
```
"A workspace scene dominated by deep navy blue (#1A2744) surfaces and warm coral (#FF6B6B) accents"
```

**Use material references to anchor colors:**
```
"Matte navy blue ceramic desk accessories, coral-colored notebook, brushed gold pen holder"
```

**Describe lighting to reinforce palette:**
```
"Warm golden hour lighting casting coral tones across navy blue architecture"
```

**Use negative prompts to exclude off-brand colors:**
```
Negative: "neon, fluorescent, pastel pink, bright green, rainbow"
```

### brand_dna.json Color Fields

The `brand_dna.json` file should include:

```json
{
  "colors": {
    "primary": "#1A2744",
    "primary_name": "deep navy blue",
    "secondary": "#FF6B6B",
    "secondary_name": "warm coral",
    "accent": "#F5A623",
    "accent_name": "golden amber",
    "neutral_light": "#F8F9FA",
    "neutral_dark": "#2D3436",
    "background": "#FFFFFF"
  },
  "color_prompt_prefix": "Deep navy blue and warm coral color palette with golden amber accents.",
  "color_negative_prompt": "neon colors, fluorescent, pastel pink, lime green, rainbow"
}
```

---

## Style Keyword Libraries

Select and combine keywords based on brand personality. Store the chosen set in `brand_dna.json` under `style_keywords`.

### Minimalist
`clean lines, generous whitespace, simple composition, muted tones, flat design, negative space, uncluttered, geometric simplicity, restrained palette, zen-like calm`

### Bold
`high contrast, saturated colors, strong typography feel, dynamic composition, graphic punch, heavy shapes, impactful, striking, dramatic lighting, powerful`

### Organic
`natural textures, earth tones, flowing shapes, botanical elements, handcrafted feel, soft edges, warm and inviting, raw materials, linen texture, terracotta`

### Tech
`futuristic, sleek surfaces, gradient meshes, glass morphism, digital interfaces, holographic, circuit patterns, dark mode aesthetic, neon accents on dark, sharp edges`

### Luxury
`rich textures, gold accents, deep jewel tones, marble surfaces, velvet, soft lighting, editorial quality, haute couture, refined, sophisticated shadows`

### Playful
`bright primary colors, rounded shapes, hand-drawn elements, confetti, pop art influence, whimsical, cartoon-like, bouncy, lighthearted, comic style`

### Corporate
`professional, polished, clean backgrounds, structured grid, business context, muted blues and grays, trustworthy, institutional, formal, boardroom aesthetic`

### Warm / Human
`candid feeling, soft natural light, genuine smiles, earth tones, community, togetherness, lifestyle photography style, authentic, relatable, cozy`

### brand_dna.json Style Fields

```json
{
  "style_keywords": ["minimalist", "clean lines", "generous whitespace", "geometric simplicity"],
  "style_prompt_suffix": "Minimalist aesthetic with clean lines, generous whitespace, and geometric simplicity.",
  "style_negative": "cluttered, busy, ornate, maximalist, chaotic"
}
```

---

## Typography in AI-Generated Images

Most AI image generators produce unreliable text. Follow these guidelines:

### In-Generation Text (Avoid When Possible)
- Limit to 1-3 words maximum in the prompt.
- Use common, short words (AI handles these better).
- Specify the text in quotes: `with the text "SALE" in bold sans-serif`.
- Accept that text may need post-processing correction.

### Post-Processing Text (Recommended)
- Generate images without text, then add text using `add_brand_overlay.py` or design tools.
- Store brand fonts in `brand_dna.json`:

```json
{
  "typography": {
    "heading_font": "Montserrat",
    "body_font": "Open Sans",
    "accent_font": "Playfair Display",
    "heading_weight": "Bold",
    "body_weight": "Regular"
  }
}
```

### Text Placement Rules
- Headlines: top third or center of image.
- Body text: never in AI-generated images (add in post-processing).
- CTA text: bottom third, high contrast background.
- Always maintain minimum 4.5:1 contrast ratio (WCAG AA).

---

## Maintaining Mood Board Consistency

When generating multiple images for a campaign or content series:

### Seed Management
- Use a base seed for the first image in a series.
- Increment by 1 for each subsequent image (seed, seed+1, seed+2, ...).
- Record seeds in a campaign manifest for reproducibility.

### Prompt Template System
Build prompts from a consistent template:

```
[brand_color_prefix] + [style_keywords] + [subject_description] + [composition_direction] + [brand_suffix]
```

Example:
```
"Deep navy and coral palette. Minimalist, clean lines. A person working at a standing desk in a modern office. Centered composition, eye-level angle. Professional and approachable mood, no text."
```

### Reference Image Anchoring
When available, use image-to-image or style reference features:
- fal.ai: Use `image_url` parameter for style reference.
- FLUX: Use img2img mode with low denoising (0.3-0.5) to maintain base style.

### Campaign Consistency Checklist
1. Same color palette prompt prefix across all images.
2. Same style keyword set for all images.
3. Same negative prompt to exclude off-brand elements.
4. Same or similar composition direction (centered, rule of thirds, etc.).
5. Same lighting direction and quality description.
6. Sequential seeds within the same batch.

---

## Color Psychology in Social Media

Use these associations when selecting brand-aligned visual treatments:

| Color | Associations | Best For |
|---|---|---|
| Blue | Trust, stability, professionalism | Finance, tech, healthcare, B2B |
| Red | Energy, urgency, passion | Sales, food, entertainment, CTA buttons |
| Green | Growth, health, nature | Wellness, sustainability, finance (growth) |
| Yellow | Optimism, warmth, attention | Retail, food, youth brands |
| Orange | Creativity, enthusiasm, action | Tech, creative agencies, SaaS |
| Purple | Luxury, wisdom, creativity | Beauty, premium brands, education |
| Pink | Playfulness, compassion, femininity | Fashion, beauty, lifestyle |
| Black | Sophistication, power, elegance | Luxury, fashion, premium tech |
| White | Purity, simplicity, cleanliness | Healthcare, minimalist brands, tech |

---

## Brand Template System

Organize reusable visual templates by content type:

```
brand_templates/
  post/
    quote_template.json
    announcement_template.json
    tip_template.json
  carousel/
    listicle_template.json
    storytelling_template.json
    tutorial_template.json
  story/
    poll_template.json
    countdown_template.json
    swipe_template.json
  thumbnail/
    talking_head_template.json
    text_focus_template.json
    collage_template.json
```

Each template JSON contains:
- `prompt_template`: Base prompt with `{subject}` placeholder.
- `dimensions`: Target width and height.
- `text_zones`: Array of text placement regions with font, size, color.
- `brand_elements`: Which brand overlays to apply (logo, border, CTA).

---

## Visual Hierarchy for Social Media

Ensure every visual communicates its message within 1-3 seconds of viewing:

### The 3-Second Rule
1. **Primary element** (takes 60% of visual weight): The main subject or headline.
2. **Secondary element** (takes 25% of visual weight): Supporting context or subheading.
3. **Tertiary element** (takes 15% of visual weight): Brand mark, CTA, or attribution.

### Hierarchy Techniques
- **Scale:** Largest element draws attention first.
- **Contrast:** Highest contrast element is seen first.
- **Color:** Brand accent color on the primary CTA or focal point.
- **Position:** Top-left to bottom-right reading flow (for LTR audiences).
- **Isolation:** Surround the key element with whitespace.

### Mobile-First Composition
- Design at actual pixel size, then zoom to 50% to test mobile legibility.
- No text smaller than 24px at 1080px width.
- Touch targets (CTA areas) at least 80x80px equivalent.
- Single focal point per image -- avoid splitting attention.
