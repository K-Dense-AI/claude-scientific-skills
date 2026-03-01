# Research Agent SDK

Claude Agent SDK 기반 연구 파이프라인 에이전트 프로젝트.

## 구조
- `research_agent/agents/` — 에이전트 정의 및 프롬프트 로더
- `research_agent/prompts/` — 에이전트별 시스템 프롬프트 (.md)
- `research_agent/utils/` — 훅, 메시지 핸들러
- `research_agent/orchestrator.py` — 파이프라인 오케스트레이터

## 개발
```bash
conda activate research-agent
pip install -e .
pytest tests/
```

## 코드 품질
- ruff: line-length=120, py311
- mypy: warn_return_any
- pytest: tests/ 디렉토리
