#!/usr/bin/env python3
"""
Culturally adapt social media content between English and Romanian.

This script handles cultural adaptation (not just translation), including
formality level, idiom adaptation, hashtag localization, and cultural
reference substitution.

Usage:
    python localize_content.py \
        --input caption.json \
        --target-language ro \
        --brand-dna brand_dna.json \
        --output caption_ro.json
"""

import argparse
import json
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Localization rules and mappings
# ---------------------------------------------------------------------------

# Common English social media phrases with Romanian cultural equivalents
PHRASE_ADAPTATIONS_EN_TO_RO = {
    "link in bio": "link in bio",  # Kept as-is (universally understood)
    "save this for later": "salveaza pentru mai tarziu",
    "drop a": "lasa un",
    "tag someone who": "tagueaza pe cineva care",
    "double tap if": "da dublu tap daca",
    "share with a friend": "trimite unui prieten",
    "let us know in the comments": "spune-ne in comentarii",
    "follow for more": "urmareste pentru mai multe",
    "swipe up": "gliseaza in sus",
    "click the link below": "acceseaza linkul de mai jos",
    "what do you think?": "tu ce crezi?",
    "agree or disagree?": "esti de acord sau nu?",
    "hot take": "parere controversata",
    "game changer": "o schimbare radicala",
    "level up": "treci la nivelul urmator",
    "stay tuned": "ramai pe fir",
    "no brainer": "alegerea e clara",
    "behind the scenes": "din culise",
    "glow up": "transformare spectaculoasa",
    "pro tip": "sfat de expert",
    "plot twist": "surpriza",
    "real talk": "sincer sa fiu",
}

PHRASE_ADAPTATIONS_RO_TO_EN = {v: k for k, v in PHRASE_ADAPTATIONS_EN_TO_RO.items()}

# False friends to catch and warn about
FALSE_FRIENDS = {
    "actually": {
        "wrong": "actual",
        "correct": "de fapt",
        "note": "'actual' in Romanian means 'current'",
    },
    "eventually": {
        "wrong": "eventual",
        "correct": "in cele din urma",
        "note": "'eventual' in Romanian means 'possibly'",
    },
    "sensible": {
        "wrong": "sensibil",
        "correct": "rezonabil",
        "note": "'sensibil' in Romanian means 'sensitive'",
    },
    "library": {
        "wrong": "librarie",
        "correct": "biblioteca",
        "note": "'librarie' in Romanian means 'bookstore'",
    },
    "magazine": {
        "wrong": "magazin",
        "correct": "revista",
        "note": "'magazin' in Romanian means 'store'",
    },
    "sympathetic": {
        "wrong": "simpatic",
        "correct": "compasional",
        "note": "'simpatic' in Romanian means 'nice/likeable'",
    },
    "resume": {
        "wrong": "rezuma",
        "correct": "CV",
        "note": "'a rezuma' in Romanian means 'to summarize'",
    },
}

# Formality transformations (tu <-> dumneavoastra patterns)
FORMALITY_PATTERNS_TU = {
    "descoperiti": "descopera",
    "aflati": "afla",
    "incercati": "incearca",
    "vizitati": "viziteaza",
    "urmariti": "urmareste",
    "comentati": "comenteaza",
    "distribuiti": "distribuie",
    "salvati": "salveaza",
    "scrieti": "scrie",
    "spuneti": "spune",
    "alegeti": "alege",
    "dumneavoastra": "tu",
    "va invitam": "te invitam",
    "va prezentam": "iti prezentam",
    "va oferim": "iti oferim",
    "va rugam": "te rugam",
}

FORMALITY_PATTERNS_DVS = {v: k for k, v in FORMALITY_PATTERNS_TU.items()}


def load_json_file(path: str) -> dict:
    """Load a JSON file."""
    file_path = Path(path)
    if not file_path.exists():
        print(f"Error: File not found: {path}", file=sys.stderr)
        sys.exit(1)
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)


def load_brand_dna(path: str | None) -> dict:
    """Load brand DNA configuration."""
    if not path:
        return {
            "language": {"ro_formality": "tu"},
            "voice": {"tone": ["friendly"]},
        }
    return load_json_file(path)


def detect_language(text: str) -> str:
    """Simple heuristic to detect if text is Romanian or English."""
    ro_indicators = [
        " si ", " sau ", " este ", " sunt ", " pentru ", " care ",
        " din ", " mai ", " acum ", " foarte ", " bine ", " asta ",
        " noastra ", " despre ", " acest ",
    ]
    en_indicators = [
        " the ", " and ", " is ", " are ", " for ", " that ",
        " with ", " this ", " from ", " your ", " our ",
        " have ", " about ", " will ",
    ]

    text_lower = f" {text.lower()} "
    ro_score = sum(1 for ind in ro_indicators if ind in text_lower)
    en_score = sum(1 for ind in en_indicators if ind in text_lower)

    return "ro" if ro_score > en_score else "en"


def adapt_phrases(text: str, target_language: str) -> tuple[str, list[str]]:
    """Replace common social media phrases with cultural equivalents."""
    changes = []
    if target_language == "ro":
        adaptations = PHRASE_ADAPTATIONS_EN_TO_RO
    else:
        adaptations = PHRASE_ADAPTATIONS_RO_TO_EN

    for source, target in adaptations.items():
        if source.lower() in text.lower():
            pattern = re.compile(re.escape(source), re.IGNORECASE)
            text = pattern.sub(target, text)
            changes.append(f"'{source}' -> '{target}'")

    return text, changes


def check_false_friends(text: str, target_language: str) -> list[dict]:
    """Check for false friend pitfalls in the text."""
    warnings = []
    if target_language == "ro":
        for en_word, info in FALSE_FRIENDS.items():
            if info["wrong"].lower() in text.lower():
                warnings.append({
                    "word": info["wrong"],
                    "suggestion": info["correct"],
                    "note": info["note"],
                })
    return warnings


def adjust_formality(text: str, formality: str) -> tuple[str, list[str]]:
    """Adjust Romanian text formality (tu vs dumneavoastra)."""
    changes = []
    if formality == "tu":
        patterns = FORMALITY_PATTERNS_TU
    else:
        patterns = FORMALITY_PATTERNS_DVS

    for source, target in patterns.items():
        if source.lower() in text.lower():
            pattern = re.compile(re.escape(source), re.IGNORECASE)
            text = pattern.sub(target, text)
            changes.append(f"Formality: '{source}' -> '{target}'")

    return text, changes


def localize_hashtags(hashtags: list[str], target_language: str) -> list[str]:
    """Adapt hashtags for the target language market."""
    ro_equivalents = {
        "#MadeInRomania": "#FacutInRomania",
        "#SmallBusiness": "#AfacereMica",
        "#SkinCare": "#IngrijirePiele",
        "#HealthyLiving": "#ViataSanatoasa",
        "#Wellness": "#Wellness",  # Kept as-is (used in RO)
        "#Fitness": "#Fitness",  # Kept as-is
        "#Marketing": "#Marketing",  # Kept as-is
        "#Startup": "#Startup",  # Kept as-is
    }

    en_equivalents = {v: k for k, v in ro_equivalents.items()}

    localized = []
    equivalents = ro_equivalents if target_language == "ro" else en_equivalents

    for tag in hashtags:
        if tag in equivalents:
            localized.append(equivalents[tag])
        else:
            localized.append(tag)

    return localized


def estimate_length_change(source_lang: str, target_lang: str, text_length: int) -> dict:
    """Estimate character count change for the target language."""
    # Romanian text tends to be 10-20% longer than English
    if source_lang == "en" and target_lang == "ro":
        factor = 1.15
    elif source_lang == "ro" and target_lang == "en":
        factor = 0.87
    else:
        factor = 1.0

    estimated = int(text_length * factor)
    return {
        "source_length": text_length,
        "estimated_target_length": estimated,
        "expansion_factor": factor,
        "note": (
            "Romanian text is typically 10-20% longer than English. "
            "Verify against platform character limits."
            if target_lang == "ro"
            else "English text is typically 10-15% shorter than Romanian."
        ),
    }


def localize_content(input_data: dict, target_language: str, brand_dna: dict) -> dict:
    """
    Main localization pipeline.

    Processes a caption JSON and produces a localized version with
    cultural adaptations, formality adjustments, and warnings.
    """
    formality = brand_dna.get("language", {}).get("ro_formality", "tu")
    all_changes = []
    all_warnings = []

    # Process caption sections
    caption = input_data.get("caption", {})
    localized_caption = {}

    for section_name, section_data in caption.items():
        if isinstance(section_data, dict) and "text" in section_data:
            text = section_data["text"]

            # Adapt phrases
            adapted_text, phrase_changes = adapt_phrases(text, target_language)
            all_changes.extend(phrase_changes)

            # Check false friends
            ff_warnings = check_false_friends(adapted_text, target_language)
            all_warnings.extend(ff_warnings)

            # Adjust formality if target is Romanian
            if target_language == "ro":
                adapted_text, form_changes = adjust_formality(adapted_text, formality)
                all_changes.extend(form_changes)

            localized_section = dict(section_data)
            localized_section["text"] = adapted_text
            localized_section["original_text"] = text
            localized_caption[section_name] = localized_section
        else:
            localized_caption[section_name] = section_data

    # Localize hashtags
    original_hashtags = input_data.get("hashtags", {}).get("suggested", [])
    branded_hashtags = input_data.get("hashtags", {}).get("branded", [])
    localized_tags = localize_hashtags(original_hashtags, target_language)

    # Estimate length changes
    total_text = " ".join(
        s.get("text", "") if isinstance(s, dict) else str(s)
        for s in caption.values()
    )
    source_lang = "en" if target_language == "ro" else "ro"
    length_estimate = estimate_length_change(source_lang, target_language, len(total_text))

    # Build output
    output = {
        "metadata": {
            **input_data.get("metadata", {}),
            "language": target_language,
            "source_language": source_lang,
            "formality": formality if target_language == "ro" else "n/a",
            "localization_type": "cultural_adaptation",
        },
        "caption": localized_caption,
        "hashtags": {
            **input_data.get("hashtags", {}),
            "suggested": localized_tags,
            "branded": branded_hashtags,
        },
        "emoji_placements": input_data.get("emoji_placements", {}),
        "localization_report": {
            "adaptations_made": all_changes,
            "false_friend_warnings": all_warnings,
            "length_estimate": length_estimate,
            "notes": [
                f"Formality set to '{formality}'" if target_language == "ro" else None,
                "Text requires human review for natural flow",
                "Hashtags include both localized and universal tags",
            ],
        },
    }

    # Clean None values from notes
    output["localization_report"]["notes"] = [
        n for n in output["localization_report"]["notes"] if n
    ]

    return output


def main():
    parser = argparse.ArgumentParser(
        description="Culturally adapt social media content between EN and RO",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to caption JSON file to localize",
    )
    parser.add_argument(
        "--target-language",
        required=True,
        choices=["en", "ro"],
        help="Target language for localization",
    )
    parser.add_argument(
        "--brand-dna",
        default=None,
        help="Path to brand_dna.json for voice and formality settings",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output file path (default: stdout)",
    )

    args = parser.parse_args()

    # Check for .env
    env_path = Path(".env")
    if env_path.exists():
        print("Note: .env file detected in working directory", file=sys.stderr)

    # Load inputs
    input_data = load_json_file(args.input)
    brand_dna = load_brand_dna(args.brand_dna)

    # Localize
    localized = localize_content(input_data, args.target_language, brand_dna)

    # Output
    output_json = json.dumps(localized, indent=2, ensure_ascii=False)

    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(output_json)
        print(f"Localized content written to {args.output}", file=sys.stderr)
    else:
        print(output_json)


if __name__ == "__main__":
    main()
