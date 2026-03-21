#!/usr/bin/env python3
"""
Generate a video storyboard for Reels, TikToks, and Shorts.

Produces a structured JSON storyboard with frame-by-frame shot descriptions,
timing, text overlays, transitions, and music direction based on proven templates.

Usage:
    python generate_video_storyboard.py \
        --template 30s --topic "3 morning habits that changed my productivity" \
        --brand-dna brand_dna.json --output output/storyboard.json

    python generate_video_storyboard.py \
        --template 60s --topic "How to build a personal brand from scratch" \
        --brand-dna brand_dna.json --output output/storyboard.json \
        --generate-keyframes --keyframe-dir output/keyframes/
"""

import argparse
import json
import sys
from pathlib import Path

TEMPLATES = {
    "15s": {
        "name": "15-Second Reel: Hook - Content - CTA",
        "total_duration": 15,
        "structure": [
            {
                "section": "HOOK",
                "start": 0,
                "end": 3,
                "duration": 3,
                "shot_type": "Close-up or text card",
                "purpose": "Stop the scroll with a bold statement, question, or surprising visual.",
                "text_overlay_style": "Large, centered, 2-5 words",
                "transition_in": "Hard cut",
                "transition_out": "Cut",
                "prompt_template": "Bold attention-grabbing visual for social media reel. {topic_hook}. {brand_style} Close-up composition, dramatic lighting.",
            },
            {
                "section": "CONTENT",
                "start": 3,
                "end": 7,
                "duration": 4,
                "shot_type": "Medium shot or B-roll",
                "purpose": "Deliver the core message or show the product.",
                "text_overlay_style": "Key phrase, lower third",
                "transition_in": "Cut",
                "transition_out": "Smooth cut",
                "prompt_template": "Social media content visual. {topic_content}. {brand_style} Medium shot, clean composition.",
            },
            {
                "section": "PROOF",
                "start": 7,
                "end": 11,
                "duration": 4,
                "shot_type": "Close-up or demo",
                "purpose": "Supporting visual, quick demo, or result.",
                "text_overlay_style": "Supporting stat or text, centered",
                "transition_in": "Smooth cut",
                "transition_out": "Cut",
                "prompt_template": "Proof or result visual. {topic_proof}. {brand_style} Close-up detail shot.",
            },
            {
                "section": "CTA",
                "start": 11,
                "end": 15,
                "duration": 4,
                "shot_type": "Medium shot or text card",
                "purpose": "Direct instruction: follow, save, link in bio, comment.",
                "text_overlay_style": "CTA text, large. Handle/username.",
                "transition_in": "Cut",
                "transition_out": "Fade or brand end card",
                "prompt_template": "Clean branded end card for social media. {brand_style} Centered composition, space for text overlay.",
            },
        ],
        "music_direction": "Upbeat, trending audio. Beat drop at 0:00. Energy peak at 0:07-0:11. Clean ending on beat.",
    },
    "30s": {
        "name": "30-Second Reel: Hook - Problem - Solution - CTA",
        "total_duration": 30,
        "structure": [
            {
                "section": "HOOK",
                "start": 0,
                "end": 3,
                "duration": 3,
                "shot_type": "Direct to camera or text card",
                "purpose": "Provocative question or bold claim to stop scrolling.",
                "text_overlay_style": "Hook text, large, bold",
                "transition_in": "Hard cut",
                "transition_out": "Cut",
                "prompt_template": "Attention-grabbing hook visual. {topic_hook}. {brand_style} Bold, dramatic composition.",
            },
            {
                "section": "PROBLEM",
                "start": 3,
                "end": 8,
                "duration": 5,
                "shot_type": "B-roll or demonstration",
                "purpose": "Show or explain the pain point. Relatable scenario.",
                "text_overlay_style": "Problem statement, lower third",
                "transition_in": "Cut",
                "transition_out": "Cut",
                "prompt_template": "Visual depicting a common problem or frustration. {topic_problem}. {brand_style} Muted tones, relatable setting.",
            },
            {
                "section": "PIVOT",
                "start": 8,
                "end": 12,
                "duration": 4,
                "shot_type": "Direct to camera or product shot",
                "purpose": "Bridge from problem to solution. 'Here's the fix.'",
                "text_overlay_style": "Pivot phrase, centered",
                "transition_in": "Whip pan or quick zoom",
                "transition_out": "Cut",
                "prompt_template": "Transition moment visual, shift from problem to solution. {topic_pivot}. {brand_style} Energetic, optimistic shift.",
            },
            {
                "section": "SOLUTION",
                "start": 12,
                "end": 20,
                "duration": 8,
                "shot_type": "Demo, walkthrough, or results",
                "purpose": "Show the solution in action. 2-3 quick steps if needed.",
                "text_overlay_style": "Step labels or key benefit text",
                "transition_in": "Cut",
                "transition_out": "Cuts between steps",
                "prompt_template": "Solution demonstration visual. {topic_solution}. {brand_style} Bright, clear, instructional composition.",
            },
            {
                "section": "PROOF",
                "start": 20,
                "end": 25,
                "duration": 5,
                "shot_type": "Results shot or before/after",
                "purpose": "Show the outcome, results, or transformation.",
                "text_overlay_style": "Result stat or before/after label",
                "transition_in": "Side-by-side or morph",
                "transition_out": "Soft cut",
                "prompt_template": "Results or transformation visual. {topic_proof}. {brand_style} Bright, satisfying, clean.",
            },
            {
                "section": "CTA",
                "start": 25,
                "end": 30,
                "duration": 5,
                "shot_type": "Direct to camera or end card",
                "purpose": "Follow for more, save this, link in bio, or engagement question.",
                "text_overlay_style": "CTA + handle. Brand colors.",
                "transition_in": "Cut",
                "transition_out": "Fade to brand card",
                "prompt_template": "Branded CTA end card. {brand_style} Clean, centered, space for text and logo.",
            },
        ],
        "music_direction": "Trending or original audio. Build energy to 0:20, peak at solution. Softer at CTA.",
    },
    "60s": {
        "name": "60-Second TikTok: Hook - Story - Twist - CTA",
        "total_duration": 60,
        "structure": [
            {
                "section": "HOOK",
                "start": 0,
                "end": 3,
                "duration": 3,
                "shot_type": "Close-up or POV",
                "purpose": "Pattern interrupt. Unexpected visual or mid-action start.",
                "text_overlay_style": "Hook text, large",
                "transition_in": "Hard cut, no fade",
                "transition_out": "Smooth cut",
                "prompt_template": "Unexpected, attention-grabbing opening frame. {topic_hook}. {brand_style} Dynamic POV or close-up.",
            },
            {
                "section": "CONTEXT",
                "start": 3,
                "end": 8,
                "duration": 5,
                "shot_type": "Medium shot or B-roll",
                "purpose": "Set the scene. Who, what, where.",
                "text_overlay_style": "Context label if needed",
                "transition_in": "Smooth cut",
                "transition_out": "Cut",
                "prompt_template": "Scene-setting establishing shot. {topic_context}. {brand_style} Medium shot, clear environment.",
            },
            {
                "section": "RISING_ACTION",
                "start": 8,
                "end": 18,
                "duration": 10,
                "shot_type": "Mixed shots, quick cuts",
                "purpose": "Build the story. Show the process, journey, or challenge.",
                "text_overlay_style": "Key phrases at each cut",
                "transition_in": "Quick cuts",
                "transition_out": "Cut",
                "prompt_template": "Action sequence visual showing progress or journey. {topic_rising}. {brand_style} Dynamic, multiple angles.",
            },
            {
                "section": "COMPLICATION",
                "start": 18,
                "end": 25,
                "duration": 7,
                "shot_type": "Close-up or reaction",
                "purpose": "The obstacle, the unexpected, the 'but then...' moment.",
                "text_overlay_style": "Tension text or 'But...' label",
                "transition_in": "Dramatic pause or slow-mo",
                "transition_out": "Cut",
                "prompt_template": "Moment of tension or complication. {topic_complication}. {brand_style} Dramatic lighting, close framing.",
            },
            {
                "section": "RESOLUTION",
                "start": 25,
                "end": 35,
                "duration": 10,
                "shot_type": "Demo or walkthrough",
                "purpose": "The breakthrough, solution, or reveal. Slower pacing.",
                "text_overlay_style": "Resolution text, step labels if tutorial",
                "transition_in": "Reveal transition (zoom out, uncover)",
                "transition_out": "Cut",
                "prompt_template": "Breakthrough or reveal moment. {topic_resolution}. {brand_style} Bright, satisfying, clear.",
            },
            {
                "section": "PAYOFF",
                "start": 35,
                "end": 45,
                "duration": 10,
                "shot_type": "Results or montage",
                "purpose": "Final result, transformation, or outcome montage.",
                "text_overlay_style": "Result stats, 'Final result' label",
                "transition_in": "Montage cuts or split screen",
                "transition_out": "Soft cut",
                "prompt_template": "Final results montage. {topic_payoff}. {brand_style} Vibrant, polished, multiple angles.",
            },
            {
                "section": "REFLECTION",
                "start": 45,
                "end": 52,
                "duration": 7,
                "shot_type": "Direct to camera",
                "purpose": "Key takeaway. Why it matters. Personal connection.",
                "text_overlay_style": "Takeaway text, clean",
                "transition_in": "Soft cut",
                "transition_out": "Cut",
                "prompt_template": "Reflective, personal moment. {topic_reflection}. {brand_style} Warm, intimate close-up.",
            },
            {
                "section": "CTA",
                "start": 52,
                "end": 60,
                "duration": 8,
                "shot_type": "End card or direct to camera",
                "purpose": "Engagement driver. Follow, comment, save, part 2.",
                "text_overlay_style": "CTA + handle + follow prompt",
                "transition_in": "Cut",
                "transition_out": "Brand end card or loop point",
                "prompt_template": "Branded end card with loop potential. {brand_style} Clean, centered, brand colors.",
            },
        ],
        "music_direction": "Story-driven audio. Tension at 0:18-0:25, release at 0:25-0:35. Volume dips under voiceover, rises during montage.",
    },
    "90s": {
        "name": "90-Second Educational Template",
        "total_duration": 90,
        "structure": [
            {
                "section": "HOOK",
                "start": 0,
                "end": 5,
                "duration": 5,
                "shot_type": "Direct to camera or motion graphic",
                "purpose": "Bold claim or common misconception to challenge.",
                "text_overlay_style": "Hook text, bold and large",
                "transition_in": "Hard cut",
                "transition_out": "Smooth cut",
                "prompt_template": "Bold educational hook visual. {topic_hook}. {brand_style} Authoritative, clean.",
            },
            {
                "section": "FOUNDATION",
                "start": 5,
                "end": 12,
                "duration": 7,
                "shot_type": "B-roll or diagram",
                "purpose": "Establish baseline knowledge. What audience thinks they know.",
                "text_overlay_style": "Definition or common belief text",
                "transition_in": "Smooth cut",
                "transition_out": "Cut",
                "prompt_template": "Educational foundation visual or diagram. {topic_foundation}. {brand_style} Clean, informational.",
            },
            {
                "section": "MISCONCEPTION",
                "start": 12,
                "end": 20,
                "duration": 8,
                "shot_type": "Direct to camera or animation",
                "purpose": "Reveal the gap: 'But here's what actually happens...'",
                "text_overlay_style": "'Myth:' or 'What most people think:' label",
                "transition_in": "Visual break effect (glitch, crack)",
                "transition_out": "Cut",
                "prompt_template": "Misconception reveal visual. {topic_misconception}. {brand_style} Visual disruption effect.",
            },
            {
                "section": "EXPLANATION",
                "start": 20,
                "end": 35,
                "duration": 15,
                "shot_type": "Mixed: diagrams, demos, B-roll",
                "purpose": "Core educational content. 2-3 sub-points with visuals.",
                "text_overlay_style": "Numbered points or keyword labels",
                "transition_in": "Cut between sub-points",
                "transition_out": "Cut",
                "prompt_template": "Educational explanation visual with diagrams. {topic_explanation}. {brand_style} Clear, structured, informational.",
            },
            {
                "section": "EXAMPLE",
                "start": 35,
                "end": 50,
                "duration": 15,
                "shot_type": "Walkthrough or demo",
                "purpose": "Concrete example or case study illustrating the explanation.",
                "text_overlay_style": "Example label, key data highlighted",
                "transition_in": "Smooth cuts within example",
                "transition_out": "Cut",
                "prompt_template": "Concrete example or case study visual. {topic_example}. {brand_style} Specific, detailed, real-world context.",
            },
            {
                "section": "APPLICATION",
                "start": 50,
                "end": 60,
                "duration": 10,
                "shot_type": "Direct to camera or infographic",
                "purpose": "Practical implication: 'What this means for you.'",
                "text_overlay_style": "'What this means:' header, action items",
                "transition_in": "Cut",
                "transition_out": "Cut",
                "prompt_template": "Practical application visual. {topic_application}. {brand_style} Action-oriented, clear.",
            },
            {
                "section": "ADVANCED_TIP",
                "start": 60,
                "end": 75,
                "duration": 15,
                "shot_type": "Quick-cut montage or demo",
                "purpose": "Bonus insight for retained viewers. Rewards watching.",
                "text_overlay_style": "'Bonus:' or 'Pro tip:' label",
                "transition_in": "Energetic cut or zoom",
                "transition_out": "Cut",
                "prompt_template": "Advanced pro tip visual. {topic_advanced}. {brand_style} Energetic, expert-level.",
            },
            {
                "section": "SUMMARY_CTA",
                "start": 75,
                "end": 90,
                "duration": 15,
                "shot_type": "Direct to camera + end card",
                "purpose": "Recap key point in one sentence. CTA to follow, save, check link.",
                "text_overlay_style": "Summary sentence, CTA, handle",
                "transition_in": "Soft cut",
                "transition_out": "Fade to brand end card",
                "prompt_template": "Summary and CTA branded end card. {brand_style} Clean, authoritative, brand colors.",
            },
        ],
        "music_direction": "Calm, focused background. Documentary tone. Subtle energy lift at advanced tip. Clean fade at end.",
    },
}


def load_brand_dna(path: str) -> dict:
    """Load brand DNA configuration."""
    brand_path = Path(path)
    if not brand_path.exists():
        return {}
    with open(brand_path, "r") as f:
        return json.load(f)


def get_brand_style_string(brand: dict) -> str:
    """Extract a style string from brand DNA for prompt insertion."""
    parts = []
    color_prefix = brand.get("color_prompt_prefix", "")
    if color_prefix:
        parts.append(color_prefix.strip())
    style_suffix = brand.get("style_prompt_suffix", "")
    if style_suffix:
        parts.append(style_suffix.strip())
    return " ".join(parts) if parts else "Professional, clean aesthetic."


def generate_storyboard(template_key: str, topic: str, brand: dict) -> dict:
    """Generate a complete storyboard from a template and topic."""
    if template_key not in TEMPLATES:
        print(f"Error: Unknown template '{template_key}'.")
        print(f"Available: {', '.join(TEMPLATES.keys())}")
        sys.exit(1)

    template = TEMPLATES[template_key]
    brand_style = get_brand_style_string(brand)

    storyboard = {
        "template": template["name"],
        "topic": topic,
        "total_duration_seconds": template["total_duration"],
        "dimensions": {"width": 1080, "height": 1920},
        "aspect_ratio": "9:16",
        "music_direction": template["music_direction"],
        "brand_style": brand_style,
        "frames": [],
        "safe_zones": {
            "top_avoid": "0-250px (platform UI: username, follow, music info)",
            "bottom_avoid": "1620-1920px (caption area, CTA buttons, navigation)",
            "right_avoid": "960-1080px in bottom half (like, comment, share icons)",
            "safe_content_area": "135-945px horizontal, 250-1620px vertical (810x1370px)",
        },
    }

    for section in template["structure"]:
        # Build the keyframe prompt
        prompt = section["prompt_template"].format(
            brand_style=brand_style,
            topic_hook=f"Topic: {topic}. Hook the viewer immediately.",
            topic_content=f"Core message about: {topic}.",
            topic_proof=f"Evidence or result related to: {topic}.",
            topic_problem=f"The problem or pain point related to: {topic}.",
            topic_pivot=f"Transition to solution for: {topic}.",
            topic_solution=f"The solution or fix for: {topic}.",
            topic_context=f"Setting and context for: {topic}.",
            topic_rising=f"Building narrative around: {topic}.",
            topic_complication=f"Unexpected challenge in: {topic}.",
            topic_resolution=f"Breakthrough moment for: {topic}.",
            topic_payoff=f"Final results of: {topic}.",
            topic_reflection=f"Key lesson from: {topic}.",
            topic_foundation=f"Baseline knowledge about: {topic}.",
            topic_misconception=f"Common misconception about: {topic}.",
            topic_explanation=f"Core explanation of: {topic}.",
            topic_example=f"Real-world example of: {topic}.",
            topic_application=f"Practical application of: {topic}.",
            topic_advanced=f"Advanced insight on: {topic}.",
        )

        frame = {
            "section": section["section"],
            "time_range": f"{section['start']:.0f}s - {section['end']:.0f}s",
            "duration_seconds": section["duration"],
            "shot_type": section["shot_type"],
            "purpose": section["purpose"],
            "text_overlay": {
                "style": section["text_overlay_style"],
                "suggested_text": f"[Write {section['section'].lower()} text for: {topic}]",
            },
            "transition": {
                "in": section["transition_in"],
                "out": section["transition_out"],
            },
            "keyframe_prompt": prompt,
        }
        storyboard["frames"].append(frame)

    return storyboard


def main():
    parser = argparse.ArgumentParser(
        description="Generate a video storyboard for Reels, TikToks, and Shorts."
    )
    parser.add_argument(
        "--template",
        required=True,
        choices=list(TEMPLATES.keys()),
        help="Storyboard duration template.",
    )
    parser.add_argument(
        "--topic",
        required=True,
        help="Video topic or title.",
    )
    parser.add_argument(
        "--brand-dna",
        default="brand_dna.json",
        help="Path to brand_dna.json file.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output path for the storyboard JSON.",
    )
    parser.add_argument(
        "--generate-keyframes",
        action="store_true",
        help="Generate key frame images for each section.",
    )
    parser.add_argument(
        "--keyframe-dir",
        default=None,
        help="Directory for generated key frame images.",
    )

    args = parser.parse_args()

    brand = load_brand_dna(args.brand_dna)

    storyboard = generate_storyboard(args.template, args.topic, brand)

    # Save storyboard
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(storyboard, f, indent=2)

    print(f"Storyboard generated: {args.output}")
    print(f"Template: {storyboard['template']}")
    print(f"Duration: {storyboard['total_duration_seconds']}s")
    print(f"Frames: {len(storyboard['frames'])}")
    print()

    for frame in storyboard["frames"]:
        print(f"  [{frame['time_range']}] {frame['section']}: {frame['purpose'][:60]}...")

    if args.generate_keyframes:
        kf_dir = Path(args.keyframe_dir) if args.keyframe_dir else output_path.parent / "keyframes"
        kf_dir.mkdir(parents=True, exist_ok=True)
        print(f"\nKey frame prompts saved in storyboard JSON.")
        print(f"To generate images, run generate_social_image.py for each frame's keyframe_prompt")
        print(f"with --platform tiktok --format video and output to {kf_dir}/")

        # Save a convenience script
        keyframe_manifest = {
            "keyframe_dir": str(kf_dir),
            "frames": [],
        }
        for i, frame in enumerate(storyboard["frames"]):
            keyframe_manifest["frames"].append({
                "index": i + 1,
                "section": frame["section"],
                "prompt": frame["keyframe_prompt"],
                "output": str(kf_dir / f"keyframe_{i + 1:02d}_{frame['section'].lower()}.png"),
            })
        manifest_path = kf_dir / "keyframe_manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(keyframe_manifest, f, indent=2)
        print(f"Keyframe manifest saved: {manifest_path}")


if __name__ == "__main__":
    main()
