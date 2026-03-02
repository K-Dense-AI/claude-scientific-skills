"""Central configuration for domain agents, tools, and pipeline mapping."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# Skill registry tool names — automatically added to every domain agent
SKILL_REGISTRY_TOOLS = ["skill_search", "skill_load", "skill_list"]

# Domain agent tool assignments (base tools per domain)
DOMAIN_TOOLS: dict[str, list[str]] = {
    "bioinformatics": ["Read", "Write", "Edit", "Bash", "Glob", "Grep", "WebSearch", "WebFetch"],
    "chemistry": ["Read", "Write", "Edit", "Bash", "Glob", "Grep", "WebSearch"],
    "clinical": ["Read", "Write", "Edit", "Bash", "Glob", "Grep", "WebSearch", "WebFetch"],
    "computation": ["Read", "Write", "Edit", "Bash", "Glob", "Grep"],
    "literature": ["Read", "Write", "Edit", "Bash", "Glob", "Grep", "WebSearch", "WebFetch"],
    "visualization": ["Read", "Write", "Edit", "Bash", "Glob", "Grep"],
    "workflow": ["Read", "Write", "Edit", "Bash", "Glob", "Grep"],
}

# Model assignments per domain
OPUS_DOMAINS = {"bioinformatics", "literature"}
# Everything else defaults to sonnet

# Pipeline phase -> domain agent mapping
PIPELINE_PHASES = [
    {
        "phase": "P1",
        "label": "문헌 검색 및 리뷰",
        "domain": "literature",
        "instruction": (
            "논문 검색, 문헌 리뷰, 사실관계 확인을 수행하라. "
            "skill_load('research-assistant') 또는 skill_load('pubmed-database')로 관련 스킬을 로드하여 사용하라."
        ),
    },
    {
        "phase": "P1.5",
        "label": "코드 사례 조사",
        "domain": "workflow",
        "instruction": (
            "GitHub 구현체를 탐색하고 매핑 문서를 작성하라. "
            "skill_load('reference-surveyor')로 절차를 따르라."
        ),
    },
    {
        "phase": "P2",
        "label": "코드 구현",
        "domain": "workflow",
        "instruction": (
            "매핑 문서 기반으로 코드를 구현하라. "
            "skill_load('code-implementer')를 로드하고, ruff+mypy+pytest 도구 체인을 적용하라."
        ),
    },
    {
        "phase": "P3",
        "label": "코드 검증",
        "domain": "workflow",
        "instruction": (
            "레퍼런스 재현 테스트 및 디버그 루프를 수행하라. "
            "skill_load('code-validator')를 로드하라."
        ),
    },
    {
        "phase": "P4",
        "label": "원고 작성",
        "domain": "literature",
        "instruction": (
            "학술 원고를 분석하고 교정하라. "
            "skill_load('manuscript-writer')와 skill_load('scientific-writing')를 로드하라."
        ),
    },
    {
        "phase": "P5",
        "label": "개선 분석",
        "domain": "workflow",
        "instruction": (
            "사이클 에러·회귀 패턴을 분석하고 개선점을 도출하라. "
            "skill_load('improvement-analyzer')를 로드하라."
        ),
    },
]

PROMPTS_DIR = Path(__file__).parent / "prompts"


@dataclass
class AgentConfig:
    """Runtime configuration."""

    prompts_dir: Path = field(default_factory=lambda: PROMPTS_DIR)
    default_cwd: Path = field(default_factory=Path.cwd)
    pipeline_phases: list[dict] = field(default_factory=lambda: list(PIPELINE_PHASES))

    def get_model(self, domain_name: str) -> str:
        if domain_name in OPUS_DOMAINS:
            return "opus"
        return "sonnet"

    def get_tools(self, domain_name: str) -> list[str]:
        """Return full tool list: domain-specific tools + skill registry tools."""
        base_tools = DOMAIN_TOOLS.get(domain_name, [])
        return base_tools + list(SKILL_REGISTRY_TOOLS)
