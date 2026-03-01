"""Load agent prompts from markdown files."""

from __future__ import annotations

from pathlib import Path

PROMPTS_DIR = Path(__file__).parent.parent / "prompts"

# Domain agents that share the commons preamble
COMMONS_AGENTS = {
    "literature",
    "bioinformatics",
}


def load_prompt(filename: str) -> str:
    """Load a single prompt file."""
    path = PROMPTS_DIR / filename
    return path.read_text(encoding="utf-8").strip()


def load_prompt_with_commons(filename: str) -> str:
    """Load a prompt with the commons preamble prepended."""
    commons = load_prompt("commons.md")
    specific = load_prompt(filename)
    return f"{commons}\n\n---\n\n{specific}"


def load_agent_prompt(agent_name: str) -> str:
    """Load the appropriate prompt for an agent, with commons if applicable."""
    filename = agent_name.replace("-", "_") + ".md"
    stem = agent_name.replace("-", "_")
    if stem in COMMONS_AGENTS:
        return load_prompt_with_commons(filename)
    return load_prompt(filename)
