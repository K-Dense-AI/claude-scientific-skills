"""Test that all prompt files load correctly."""

from __future__ import annotations

import pytest

from research_agent.agents.prompts import load_prompt, load_prompt_with_commons

PROMPT_FILES = [
    "commons.md",
    "orchestrator.md",
    "bioinformatics.md",
    "chemistry.md",
    "clinical.md",
    "computation.md",
    "literature.md",
    "visualization.md",
    "workflow.md",
]


@pytest.mark.parametrize("filename", PROMPT_FILES)
def test_prompt_loads(filename):
    content = load_prompt(filename)
    assert len(content) > 50, f"Prompt {filename} is too short ({len(content)} chars)"


def test_commons_prepended_to_literature():
    content = load_prompt_with_commons("literature.md")
    # commons content should appear first
    assert "APA 7th Edition" in content
    # literature-specific content should appear after
    assert "논문 검색" in content


def test_commons_prepended_to_bioinformatics():
    content = load_prompt_with_commons("bioinformatics.md")
    assert "APA 7th Edition" in content
    assert "유전체" in content


def test_commons_not_in_computation():
    from research_agent.agents.prompts import load_agent_prompt

    content = load_agent_prompt("computation")
    assert "APA 7th Edition" not in content
    assert "ML/DL" in content


def test_commons_not_in_workflow():
    from research_agent.agents.prompts import load_agent_prompt

    content = load_agent_prompt("workflow")
    assert "APA 7th Edition" not in content
    assert "Git" in content
