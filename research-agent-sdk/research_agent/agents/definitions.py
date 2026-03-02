"""Build AgentDefinition instances for 7 domain agents."""

from __future__ import annotations

from claude_agent_sdk import AgentDefinition

from research_agent.agents.prompts import load_agent_prompt
from research_agent.config import AgentConfig

# Agent descriptions tell the orchestrator WHEN to delegate to each agent
AGENT_DESCRIPTIONS: dict[str, str] = {
    "bioinformatics": (
        "유전체, 단백질, 단일세포, 대사 모델링, 실험실 자동화, 바이오 데이터베이스. "
        "바이오인포매틱스 관련 질문이나 실험 데이터 분석 시 사용."
    ),
    "chemistry": (
        "분자 설계, 약물 발견, 대사체학, 재료과학, 화학 데이터베이스. "
        "화학/약학 관련 질문이나 분자 분석 시 사용."
    ),
    "clinical": (
        "임상 보고서, 임상시험, 치료 계획, FDA, 약물유전체학. "
        "의료/임상 관련 질문 시 사용."
    ),
    "computation": (
        "ML/DL, 통계, 최적화, 양자컴퓨팅, 시뮬레이션, 대용량 데이터 처리. "
        "계산/모델링 관련 질문이나 데이터 분석 시 사용."
    ),
    "literature": (
        "논문 검색, 문헌 리뷰, 학술 작문, 인용 관리, 피어 리뷰, Discussion 구조화. "
        "문헌 기반 질문, P1 단계, P4 단계, Discussion 작업 시 사용."
    ),
    "visualization": (
        "그래프, Figure, 프레젠테이션, 포스터, 다이어그램, 인포그래픽. "
        "시각화/발표 자료 작성 시 사용."
    ),
    "workflow": (
        "코드 구현, 검증, Git 관리, 사이클 분석, 실험 관리, 코드 사례 조사. "
        "P1.5, P2, P3, P5 단계 또는 코드/Git 관련 작업 시 사용."
    ),
}


def build_agents(config: AgentConfig) -> dict[str, AgentDefinition]:
    """Build all domain agent definitions with tools and model assignments."""
    agents: dict[str, AgentDefinition] = {}

    for name, description in AGENT_DESCRIPTIONS.items():
        prompt = load_agent_prompt(name)
        tools = config.get_tools(name)
        model = config.get_model(name)

        agents[name] = AgentDefinition(
            description=description,
            prompt=prompt,
            tools=tools if tools else None,
            model=model,
        )

    return agents
