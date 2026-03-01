"""Skill Registry MCP Server — exposes skill_search, skill_load, skill_list as MCP tools."""

from __future__ import annotations

from pathlib import Path

from claude_agent_sdk import create_sdk_mcp_server, tool

from research_agent.domains import DOMAIN_SKILLS
from research_agent.skills.registry import SkillRegistry

# Module-level singleton, initialized lazily
_registry: SkillRegistry | None = None


def _get_registry() -> SkillRegistry:
    """Get or create the singleton SkillRegistry."""
    global _registry  # noqa: PLW0603
    if _registry is None:
        skills_dir = Path(__file__).resolve().parents[3] / "scientific-skills"
        _registry = SkillRegistry(skills_dir=skills_dir, domain_map=DOMAIN_SKILLS)
    return _registry


@tool(
    "skill_search",
    "Search for scientific skills by keyword query. Returns matching skill names, descriptions, and domains.",
    {
        "type": "object",
        "properties": {
            "query": {
                "type": "string",
                "description": "Search keywords (e.g., 'single cell RNA', 'molecular docking', 'time series')",
            },
            "domain": {
                "type": "string",
                "description": "Optional domain filter",
                "enum": ["bioinformatics", "chemistry", "clinical", "computation", "literature", "visualization", "workflow"],
            },
            "limit": {
                "type": "integer",
                "description": "Maximum results (default: 10, max: 30)",
                "default": 10,
            },
        },
        "required": ["query"],
    },
)
async def skill_search(args: dict) -> dict:
    """Search skills and return formatted results."""
    registry = _get_registry()
    query = args["query"]
    domain = args.get("domain")
    limit = min(args.get("limit", 10), 30)

    results = registry.search(query, domain=domain, limit=limit)

    if not results:
        text = f"No skills found for query: '{query}'"
        if domain:
            text += f" in domain: {domain}"
        return {"content": [{"type": "text", "text": text}]}

    lines = [f"Found {len(results)} skill(s) for '{query}':\n"]
    for entry in results:
        desc = entry.description[:150] if entry.description else "(no description)"
        lines.append(f"- **{entry.name}** [{entry.domain}]: {desc}")

    return {"content": [{"type": "text", "text": "\n".join(lines)}]}


@tool(
    "skill_load",
    "Load the full content of a scientific skill by name. Returns the complete SKILL.md with instructions and code examples.",
    {
        "type": "object",
        "properties": {
            "skill_name": {
                "type": "string",
                "description": "Exact skill name (e.g., 'scanpy', 'rdkit', 'pubmed-database')",
            },
        },
        "required": ["skill_name"],
    },
)
async def skill_load(args: dict) -> dict:
    """Load and return full SKILL.md content."""
    registry = _get_registry()
    skill_name = args["skill_name"]
    content = registry.load(skill_name)

    if content is None:
        candidates = registry.search(skill_name, limit=3)
        if candidates:
            suggestions = ", ".join(c.name for c in candidates)
            text = f"Skill '{skill_name}' not found. Did you mean: {suggestions}?"
        else:
            text = f"Skill '{skill_name}' not found. Use skill_search to find available skills."
        return {"content": [{"type": "text", "text": text}]}

    return {"content": [{"type": "text", "text": content}]}


@tool(
    "skill_list",
    "List all available scientific skills, optionally filtered by domain. Returns skill names and short descriptions.",
    {
        "type": "object",
        "properties": {
            "domain": {
                "type": "string",
                "description": "Optional domain filter",
                "enum": ["bioinformatics", "chemistry", "clinical", "computation", "literature", "visualization", "workflow"],
            },
        },
        "required": [],
    },
)
async def skill_list(args: dict) -> dict:
    """List all skills with optional domain filter."""
    registry = _get_registry()
    domain = args.get("domain")
    entries = registry.list_skills(domain=domain)

    if not entries:
        text = "No skills found"
        if domain:
            text += f" in domain: {domain}"
        return {"content": [{"type": "text", "text": text}]}

    header = f"Skills in domain '{domain}'" if domain else "All skills"
    lines = [f"{header} ({len(entries)} total):\n"]
    for entry in entries:
        if entry.description:
            desc = entry.description[:100] + "..." if len(entry.description) > 100 else entry.description
            lines.append(f"- {entry.name}: {desc}")
        else:
            lines.append(f"- {entry.name}")

    return {"content": [{"type": "text", "text": "\n".join(lines)}]}


def create_skill_registry_server():
    """Create the MCP server configuration for the skill registry."""
    return create_sdk_mcp_server(
        name="skill-registry",
        version="0.1.0",
        tools=[skill_search, skill_load, skill_list],
    )


# MCP tool names for adding to allowed_tools
SKILL_REGISTRY_TOOLS = ["skill_search", "skill_load", "skill_list"]
