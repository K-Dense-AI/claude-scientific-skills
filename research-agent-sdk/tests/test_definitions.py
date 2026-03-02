"""Test domain agent definitions load correctly."""

from __future__ import annotations

import pytest

from research_agent.config import AgentConfig

EXPECTED_AGENTS = [
    "bioinformatics",
    "chemistry",
    "clinical",
    "computation",
    "literature",
    "visualization",
    "workflow",
]


@pytest.fixture
def config():
    return AgentConfig()


def test_all_agents_defined(config):
    from research_agent.agents.definitions import build_agents

    agents = build_agents(config)
    for name in EXPECTED_AGENTS:
        assert name in agents, f"Missing agent: {name}"
        assert agents[name].description, f"Empty description: {name}"
        assert agents[name].prompt, f"Empty prompt: {name}"


def test_agent_count(config):
    from research_agent.agents.definitions import build_agents

    agents = build_agents(config)
    assert len(agents) == len(EXPECTED_AGENTS)


def test_tool_scoping(config):
    from research_agent.agents.definitions import build_agents

    agents = build_agents(config)

    # Every agent should have skill registry tools
    for name in EXPECTED_AGENTS:
        agent_tools = agents[name].tools or []
        assert "skill_search" in agent_tools, f"{name} missing skill_search"
        assert "skill_load" in agent_tools, f"{name} missing skill_load"
        assert "skill_list" in agent_tools, f"{name} missing skill_list"

    # bioinformatics should have WebSearch
    bio_tools = agents["bioinformatics"].tools or []
    assert "WebSearch" in bio_tools

    # literature should have WebSearch and WebFetch
    lit_tools = agents["literature"].tools or []
    assert "WebSearch" in lit_tools
    assert "WebFetch" in lit_tools

    # workflow should have Edit and Bash
    wf_tools = agents["workflow"].tools or []
    assert "Edit" in wf_tools
    assert "Bash" in wf_tools


def test_model_assignments(config):
    assert config.get_model("bioinformatics") == "opus"
    assert config.get_model("literature") == "opus"
    assert config.get_model("chemistry") == "sonnet"
    assert config.get_model("computation") == "sonnet"
    assert config.get_model("visualization") == "sonnet"
    assert config.get_model("workflow") == "sonnet"
    assert config.get_model("clinical") == "sonnet"
