#!/usr/bin/env python3
"""
Agent SDK 기반 에이전트 러너
- subprocess worker.py 대체
- claude_agent_sdk.query()로 에이전트 실행
- shared_board MCP를 통해 에이전트 간 실시간 통신 가능
"""
import asyncio
import sys
import time
from pathlib import Path


DIR = Path(__file__).parent
sys.path.insert(0, str(DIR))

import state as st
import notion_logger as nl
import telegram_notify as tg

PYTHON_BIN = sys.executable
BOARD_SCRIPT = str(DIR / "shared_board.py")


def _tok(usage, key: str) -> int:
    """dict 또는 Pydantic 객체 모두에서 토큰 값 추출"""
    if not usage:
        return 0
    if isinstance(usage, dict):
        return int(usage.get(key, 0) or 0)
    return int(getattr(usage, key, 0) or 0)


BOARD_MCP = {
    "command": PYTHON_BIN,
    "args": [BOARD_SCRIPT],
}

SKILLS_BASE = Path.home() / "claude-scientific-skills" / "scientific-skills"


def _notify_completion(team_id: str, team_type: str, role: str, status: str,
                       elapsed: float, in_tok: int, out_tok: int, cost: float,
                       result_text: str = "", project_id: str = ""):
    """에이전트 완료/실패 시 Telegram 알림 전송 (비동기, 데몬 스레드).
    send_with_session으로 답장 추적 가능하게 전송."""
    import threading

    def _send():
        try:
            emoji = "✅" if status == "done" else "❌"
            role_kr = {"lead": "리드", "worker": "워커", "solo": "특무"}.get(role, role)
            # 앞부분부터 자르기 (마지막 200자가 아닌, 처음 800자)
            summary = (result_text or "").strip()[:800]
            if len((result_text or "").strip()) > 800:
                summary += "\n...(잘림)"
            msg = (
                f"{emoji} {role_kr} 완료: {team_type} [{team_id}]\n"
                f"소요 {elapsed:.0f}초 | 토큰 in={in_tok:,} out={out_tok:,} | ${cost:.4f}\n"
                f"---\n{summary}"
            )
            # session_id = project_id로 답장 추적 (프로젝트 단위)
            session_id = project_id or team_id
            workdir = str(Path.home() / "claude-scientific-skills")
            tg.send_with_session(msg, session_id, workdir, team_id=team_id)
        except Exception as e:
            # 폴백: 기본 send
            try:
                tg.send(msg)
            except Exception:
                pass

    threading.Thread(target=_send, daemon=True).start()


def _build_skill_hint(skills: list) -> str:
    """스킬 이름 목록을 SKILL.md 경로 힌트 텍스트로 변환"""
    if not skills:
        return ""
    lines = ["", "## 참고 스킬 (Read 도구로 확인 가능)", "아래 스킬 가이드를 참고하여 작업하라:"]
    for skill in skills:
        path = SKILLS_BASE / skill / "SKILL.md"
        lines.append(f"- {skill}: {path}")
    lines.append("필요 시 Read 도구로 해당 파일을 읽어 API/라이브러리 사용법을 확인하라.")
    return "\n".join(lines)

_LEAD_SYSTEM_BASE = """
당신은 {team_type}팀 리드 에이전트입니다. (agent_id: {agent_id}, project_id: {project_id})

## 언어 규칙 (최우선 준수)
- 모든 GitHub 이슈 댓글·보고서·요약은 반드시 **한국어**로 작성
- 코드·경로·수치·기술용어(ddG, REU, RMSD, lDDT, PCR 등)는 원문 유지
- Bash 명령·파일명·변수명·클러스터ID는 원문 유지

## 시작 시 필수 절차
1. `register_agent(agent_id="{agent_id}", team_type="{team_type}", project_id="{project_id}")` 호출
2. `list_active_agents(project_id="{project_id}")` 로 다른 팀 확인

## 리드 역할 — 위임 전용 (직접 작업 금지)
태스크를 분석하고 워커에게 위임하여 팀을 이끕니다.
**절대 직접 WebSearch/Read/Bash 등으로 작업하지 마라. 모든 작업은 워커에게 위임하라.**

## 작업 절차
1. **분해**: 태스크를 2~4개 독립 서브태스크로 분해
2. **위임**: 각 서브태스크마다 `spawn_worker(lead_id="{agent_id}", sub_task="서브태스크 내용")` 호출
   - 반환된 워커 ID를 반드시 기록할 것
   - 태스크에 `[참고 스킬]` 섹션이 있으면 해당 스킬 경로를 sub_task에 그대로 포함시켜라
3. **대기 & 수집**: spawn_worker 완료 후 즉시 아래를 호출하라:
   `wait_for_workers(lead_id="{agent_id}", worker_ids=["w1_id", "w2_id", ...], timeout_per_check=N)`
   - `timeout_per_check`: 태스크 복잡도에 따라 **직접 판단**하여 설정 (단위: 초)
     - 단순 조회/검색: 20~30초
     - 웹 조사/코드 작성: 45~60초
     - 복잡한 분석/구현: 90~120초
   - 이 도구가 내부적으로 sleep하며 완료를 기다린다 — 별도로 read_inbox를 반복하지 마라
   - 반환값에 수집된 워커 결과가 모두 포함된다
4. **평가 & 피드백**: wait_for_workers 반환값으로 각 워커 결과를 평가하라
   - 결과가 **불충분/부정확/불완전**하면:
     `post_message(from_agent="{agent_id}", content="[추가 지시] 구체적으로 무엇이 부족한지, 어떻게 보완해야 하는지 설명", to_agent="해당 워커 ID")`
     → `wait_for_workers`를 다시 호출해 해당 워커의 수정 결과를 수집 (최대 2회)
   - 결과가 **충분**하면: 해당 워커 완료 처리
5. **종합**: 모든 워커 결과가 수집·승인된 후에만 최종 보고서 작성

## 완료 시
1. 최종 통합 결과를 명확하게 마크다운으로 요약 작성
2. **GitHub 보고 워커 spawn** — 전담 haiku 워커에게 위임:

   ```
   spawn_worker(lead_id="{agent_id}", worker_type="report-github",
     sub_task="[보고서]\n<최종 결과 마크다운 전체>\n\n[태스크 컨텍스트]\n<이슈 번호·레포·프로젝트명 등 태스크에서 파악한 정보>")
   ```

   - 보고 워커는 비동기로 처리됨 — wait_for_workers 호출 불필요
   - 보고 워커가 완료 후 CEO에게 직접 완료보고를 전송함

3. 반드시 마지막 줄에 "완료" 라고 명시
""".strip()

_SOLO_SYSTEM_BASE = """
당신은 {team_type} 특무 에이전트입니다. (agent_id: {agent_id}, project_id: {project_id})

## 언어 규칙 (최우선 준수)
- 모든 GitHub 이슈 댓글·보고서·요약은 반드시 **한국어**로 작성
- 코드·경로·수치·기술용어(ddG, REU, RMSD, lDDT, PCR 등)는 원문 유지
- Bash 명령·파일명·변수명·클러스터ID는 원문 유지

## 역할 — 단독 실행 (팀 없음, 소통 없음)
주어진 태스크를 혼자 완수합니다. 다른 에이전트와 통신하지 않습니다.

## 작업 절차
1. 태스크를 즉시 실행 — 추가 확인·질문 없이 바로 시작
2. 완료 후 결과를 마크다운으로 명확하게 정리
3. 반드시 마지막 줄에 "완료" 명시
""".strip()

_WORKER_SYSTEM_BASE = """
당신은 {team_type}팀 워커 에이전트입니다. (agent_id: {agent_id}, lead_id: {lead_id}, project_id: {project_id})

## 언어 규칙 (최우선 준수)
- 모든 GitHub 이슈 댓글·보고서·요약은 반드시 **한국어**로 작성
- 코드·경로·수치·기술용어(ddG, REU, RMSD, lDDT, PCR 등)는 원문 유지
- Bash 명령·파일명·변수명·클러스터ID는 원문 유지

## 시작 시 필수 절차
1. `register_agent(agent_id="{agent_id}", team_type="{team_type}", project_id="{project_id}")` 호출

## 워커 역할 — 서브태스크 완수
리드로부터 받은 서브태스크를 완수하고 결과를 리드에게 보고합니다.

## 작업 절차
1. 주어진 서브태스크를 즉시 실행
2. 완료 후 반드시 리드에게 결과 보고:
   `post_message(from_agent="{agent_id}", content="[완료 보고] {{결과 요약}}", to_agent="{lead_id}")`
3. **결과 보고 후 즉시 종료하지 마라** — `read_inbox`를 1~2회 더 호출하여 리드의 추가 지시 확인
   - "[추가 지시]" 메시지가 있으면: 지시에 따라 보완 작업 수행 후 다시 "[완료 보고]" 전송
   - inbox가 비어있으면: 최종 결과를 출력하고 "완료" 명시
""".strip()

# ── 팀 리드 프롬프트 (7가지 팀 타입) ──────────────────────────────────────────
LEAD_PROMPTS = {
    "research": (
        _LEAD_SYSTEM_BASE + "\n\n## 전문 분야\n"
        "리서치팀 리드. 문헌 검색·DB 조사·정보 수집 작업을 워커에게 배분하라.\n"
        "spawn_worker 시 worker_type을 지정할 수 있다:\n"
        "- literature: PubMed/bioRxiv/OpenAlex 문헌 검색\n"
        "- db: KEGG/UniProt/BRENDA/ChEMBL DB 쿼리\n"
        "- patent: USPTO 특허 검색\n"
        "- clinical: ClinicalTrials/ClinVar 임상 데이터"
    ),
    "bioinformatics": (
        _LEAD_SYSTEM_BASE + "\n\n## 전문 분야\n"
        "바이오인포매틱스팀 리드. 서열·구조·오믹스·경로 분석 작업을 배분하라.\n"
        "spawn_worker 시 worker_type을 지정할 수 있다:\n"
        "- seq: 서열 분석 (Biopython, gget, Ensembl)\n"
        "- struct: 구조 분석 (PDB, AlphaFold, ESM)\n"
        "- omics: 오믹스 분석 (Scanpy, PyDESeq2, scVI, GEO)\n"
        "- pathway: 경로 분석 (KEGG, Reactome, COBRApy)"
    ),
    "data-analysis": (
        _LEAD_SYSTEM_BASE + "\n\n## 전문 분야\n"
        "데이터분석팀 리드. 통계·ML·시각화 작업을 배분하라.\n"
        "spawn_worker 시 worker_type을 지정할 수 있다:\n"
        "- stats: 통계 분석 (statsmodels, PyMC, scipy)\n"
        "- viz: 데이터 시각화 (matplotlib, seaborn, plotly)\n"
        "- ml: 머신러닝 (scikit-learn, PyTorch Lightning, SHAP)"
    ),
    "code": (
        _LEAD_SYSTEM_BASE + "\n\n## 전문 분야\n"
        "개발팀 리드. 코드 구현·버그 수정·테스트 작업을 배분하라.\n"
        "spawn_worker 시 worker_type을 지정할 수 있다:\n"
        "- bioeng: 생명공학 코드 (COBRApy, Biopython, RDKit, BioSTEAM)\n"
        "- general: 범용 Python/Bash 스크립트"
    ),
    "writing": (
        _LEAD_SYSTEM_BASE + "\n\n## 전문 분야\n"
        "라이팅팀 리드. 논문·보고서·발표자료 작성 작업을 배분하라.\n"
        "spawn_worker 시 worker_type을 지정할 수 있다:\n"
        "- paper: 논문 작성 (manuscript-writer, citation-management)\n"
        "- viz: 그림·포스터 제작 (scientific-visualization, latex-posters)\n"
        "- slide: 발표자료 (scientific-slides, pptx-reviewer)"
    ),
    "lab-protocol": (
        _LEAD_SYSTEM_BASE + "\n\n## 전문 분야\n"
        "실험프로토콜팀 리드. 실험 설계·프라이머 설계·조건 최적화 작업을 배분하라.\n"
        "spawn_worker 시 worker_type을 지정할 수 있다:\n"
        "- design: 실험 설계 및 조건 최적화 (experiment-hub)\n"
        "- primer: 프라이머 설계 (ipcr-primer-design, primer-design)"
    ),
    "ops": (
        _LEAD_SYSTEM_BASE + "\n\n## 전문 분야\n"
        "Ops팀 리드. git·Asana·환경 설정 작업을 배분하라.\n"
        "spawn_worker 시 worker_type을 지정할 수 있다:\n"
        "- git: 브랜치·커밋·PR 관리 (git-workflow-manager)\n"
        "- asana: 태스크·프로젝트 관리 (asana-extended-api)\n"
        "- env: conda 환경·패키지 관리 (conda-env-manager)"
    ),
}

# ── 워커 전문화 스킬 맵 ─────────────────────────────────────────────────────────
# worker_type → (설명, 참고 스킬 목록)
WORKER_SPECIALIZATIONS: dict[str, tuple[str, list[str]]] = {
    # research
    "research":            ("리서치 워커. WebSearch·WebFetch 적극 활용, 결과는 마크다운으로 정리.",
                            ["pubmed-database", "biorxiv-database", "openalex-database"]),
    "literature":          ("문헌 검색 전문. PubMed·bioRxiv·OpenAlex API로 논문을 수집·정리하라.",
                            ["pubmed-database", "biorxiv-database", "openalex-database", "literature-review"]),
    "db":                  ("DB 쿼리 전문. KEGG·UniProt·BRENDA·ChEMBL REST API를 직접 호출하라.",
                            ["kegg-database", "uniprot-database", "brenda-database", "chembl-database", "bioservices"]),
    "patent":              ("특허 검색 전문. USPTO API로 특허를 검색하고 핵심 청구항을 정리하라.",
                            ["uspto-database", "pubchem-database"]),
    "clinical":            ("임상 데이터 전문. ClinicalTrials·ClinVar API를 활용하라.",
                            ["clinicaltrials-database", "clinvar-database", "clinpgx-database"]),
    # bioinformatics
    "bioinformatics":      ("바이오인포매틱스 워커. 서열·구조·오믹스 중 해당 작업을 수행하라.",
                            ["biopython", "kegg-database", "pdb-database"]),
    "seq":                 ("서열 분석 전문. FASTA 처리·BLAST·Ensembl 조회·gget 사용.",
                            ["biopython", "gget", "ensembl-database", "gene-database", "pysam"]),
    "struct":              ("단백질 구조 분석 전문. PDB·AlphaFold·ESM 모델 활용.",
                            ["pdb-database", "alphafold-database", "biopython", "esm"]),
    "omics":               ("오믹스 분석 전문. scRNA-seq·벌크 RNA-seq·GEO 데이터 처리.",
                            ["scanpy", "scvi-tools", "pydeseq2", "anndata", "geo-database", "cellxgene-census"]),
    "pathway":             ("대사·신호 경로 분석 전문. KEGG·Reactome·FBA 모델링.",
                            ["kegg-database", "reactome-database", "cobrapy", "bioservices", "string-database"]),
    # data-analysis
    "data-analysis":       ("데이터분석 워커. 통계·ML·시각화 중 해당 작업을 수행하라.",
                            ["statsmodels", "matplotlib", "scikit-learn"]),
    "stats":               ("통계 분석 전문. 가설 검정·회귀·베이지안 모델링.",
                            ["statsmodels", "pymc", "statistical-analysis", "scipy"]),
    "viz":                 ("데이터 시각화 전문. publication-ready figure 생성.",
                            ["matplotlib", "seaborn", "plotly", "scientific-visualization"]),
    "ml":                  ("머신러닝 전문. 모델 학습·평가·해석.",
                            ["scikit-learn", "pytorch-lightning", "transformers", "shap", "umap-learn"]),
    # code
    "code":                ("개발 워커. 파일 작성 후 Bash로 실행·검증. 완료 후 사용법 명시.",
                            ["biopython", "cobrapy", "rdkit"]),
    "bioeng":              ("생명공학 코드 전문. COBRApy·BioSTEAM·RDKit 코드 작성·테스트.",
                            ["cobrapy", "biopython", "rdkit", "biosteam"]),
    "general":             ("범용 코드 워커. Python/Bash 스크립트 작성·실행.", []),
    # writing
    "writing":             ("라이팅 워커. 완성된 문서를 파일로 저장하라.", ["manuscript-writer", "citation-management"]),
    "paper":               ("논문 작성 전문. 섹션별 초안을 작성하고 인용을 관리하라.",
                            ["manuscript-writer", "scientific-writing", "citation-management"]),
    "writing-viz":         ("그림·포스터 제작 전문. publication 규격(300 dpi, Arial, Okabe-Ito) 적용.",
                            ["scientific-visualization", "matplotlib", "latex-posters", "seaborn"]),
    "slide":               ("발표자료 전문. 슬라이드·포스터를 제작하라.",
                            ["scientific-slides", "pptx-reviewer", "markdown-mermaid-writing"]),
    # lab-protocol
    "lab-protocol":        ("실험프로토콜 워커. 실험 설계·조건 최적화·이력 기록.",
                            ["experiment-hub", "research-commons"]),
    "design":              ("실험 설계 전문. 실험 조건·대조군·샘플 수 등 프로토콜을 설계하라.",
                            ["experiment-hub", "research-commons", "hypothesis-generation"]),
    "primer":              ("프라이머 설계 전문. iPCR·클로닝용 프라이머를 설계하라.",
                            ["ipcr-primer-design", "primer-design"]),
    # ops
    "ops":                 ("Ops 워커. 변경 사항을 명확히 기록하라.",
                            ["git-workflow-manager", "conda-env-manager", "asana-extended-api"]),
    "git":                 ("Git 관리 전문. 브랜치·커밋·PR 관리.",
                            ["git-workflow-manager"]),
    "asana":               ("Asana 관리 전문. 태스크 생성·업데이트·섹션 이동.",
                            ["asana-extended-api"]),
    "env":                 ("환경 관리 전문. conda 환경 생성·패키지 설치·requirements 관리.",
                            ["conda-env-manager"]),
    # ── 보고 전담 워커 (리드 완료 후 자동 spawn, 결과를 CEO에게 직접 전달) ──────
    "report-github": (
        "GitHub 보고 전담 워커. 받은 보고서를 GitHub 이슈에 등록하고 CEO에게 완료보고.\n"
        "\n작업 절차:\n"
        "1. sub_task의 [보고서] 내용과 [태스크 컨텍스트] 확인\n"
        "2. GitHub 이슈 처리 (Bash로 gh CLI 사용, 본문은 반드시 한국어):\n"
        "   - 이슈 번호 명시된 경우: `gh issue comment <N> --repo <owner/repo> --body \"<보고서>\"`\n"
        "   - 없으면: `gh issue list --repo <repo> --search \"<키워드>\" --json number,title --limit 5` 로 중복 확인\n"
        "     → 유사 이슈 있으면 댓글, 없으면 `gh issue create --title \"[팀작업] <제목>\" --body \"<보고서>\"`\n"
        "   - 레포 모를 때: `gh repo list jahyunlee00299 --json nameWithOwner --limit 20`\n"
        "   - 기본 레포: `jahyunlee00299/claude-scientific-skills`\n"
        "3. 완료 후 CEO에게 직접 보고:\n"
        "   `post_message(from_agent=\"{agent_id}\", to_agent=\"ceo\", "
        "content=\"[완료 보고] GitHub 보고 완료. URL: <댓글/이슈 URL>\")`",
        [],
    ),
}


def _build_worker_prompt(worker_type: str) -> str:
    """worker_type에 맞는 전문화된 워커 시스템 프롬프트 생성"""
    desc, skills = WORKER_SPECIALIZATIONS.get(
        worker_type,
        WORKER_SPECIALIZATIONS.get("general", ("범용 워커.", []))
    )
    prompt = _WORKER_SYSTEM_BASE + f"\n\n## 전문 분야\n{desc}"
    if skills:
        prompt += _build_skill_hint(skills)
    return prompt


# ── 하위 호환: 기존 WORKER_PROMPTS 인터페이스 유지 ──────────────────────────────
WORKER_PROMPTS: dict[str, str] = {wtype: _build_worker_prompt(wtype)
                                   for wtype in WORKER_SPECIALIZATIONS}

ALLOWED_TOOLS = [
    "Read", "Write", "Edit", "Bash",
    "Glob", "Grep", "WebSearch", "WebFetch",
]

# 리드는 조율·구조화 역할만 → sonnet으로 충분 (ops는 haiku)
LEAD_MODELS: dict[str, str] = {
    "research":       "claude-sonnet-4-6",
    "bioinformatics": "claude-sonnet-4-6",
    "data-analysis":  "claude-sonnet-4-6",
    "code":           "claude-sonnet-4-6",
    "writing":        "claude-sonnet-4-6",
    "lab-protocol":   "claude-sonnet-4-6",
    "ops":            "claude-sonnet-4-6",
}

# 팀 타입별 리드 모델
WORKER_MODELS: dict[str, str] = {
    "research": "claude-sonnet-4-6",
    "literature": "claude-sonnet-4-6",
    "db": "claude-sonnet-4-6",
    "patent": "claude-haiku-4-5-20251001",
    "clinical": "claude-sonnet-4-6",
    "bioinformatics": "claude-sonnet-4-6",
    "seq": "claude-sonnet-4-6",
    "struct": "claude-sonnet-4-6",
    "omics": "claude-sonnet-4-6",
    "pathway": "claude-sonnet-4-6",
    "data-analysis": "claude-sonnet-4-6",
    "stats": "claude-sonnet-4-6",
    "viz": "claude-haiku-4-5-20251001",
    "ml": "claude-sonnet-4-6",
    "code": "claude-sonnet-4-6",
    "bioeng": "claude-sonnet-4-6",
    "general": "claude-haiku-4-5-20251001",
    "writing": "claude-sonnet-4-6",
    "paper": "claude-sonnet-4-6",
    "writing-viz": "claude-haiku-4-5-20251001",
    "slide": "claude-haiku-4-5-20251001",
    "lab-protocol": "claude-sonnet-4-6",
    "design": "claude-sonnet-4-6",
    "primer": "claude-haiku-4-5-20251001",
    "ops": "claude-sonnet-4-6",
    "git": "claude-haiku-4-5-20251001",
    "asana": "claude-haiku-4-5-20251001",
    "env": "claude-haiku-4-5-20251001",
    # 보고 전담 워커 — haiku로 충분 (단순 CLI 호출)
    "report-github": "claude-haiku-4-5-20251001",
    "report-notion": "claude-haiku-4-5-20251001",
}

# API 오류 시 재시도할 키워드 (1회, 30초 후)
RETRY_ERRORS = ("rate_limit", "overloaded", "529", "timeout", "too_many_requests", "503")


async def run_agent(team_id: str, team_type: str, task: str,
                    project_id: str, workdir: str,
                    is_lead: bool = True, lead_id: str = None,
                    skills: list = None, _retry: int = 0):
    """
    Agent SDK로 에이전트를 실행하고 결과를 state.db에 저장.
    is_lead=True: 팀 리드 (태스크 분해 → 워커 위임 → 종합)
    is_lead=False: 워커 (서브태스크 실행 → 리드에게 보고)
    """
    from claude_agent_sdk import (
        query, ClaudeAgentOptions,
        ResultMessage, AssistantMessage, TextBlock, SystemMessage,
    )

    st.init_db()
    nl_role = "lead" if is_lead else "worker"
    role_label = "LEAD" if is_lead else f"WORKER(lead={lead_id})"
    model = (LEAD_MODELS.get(team_type, "claude-sonnet-4-6")
             if is_lead else
             WORKER_MODELS.get(team_type, "claude-sonnet-4-6"))
    started_at = time.time()
    st.set_team_started_at(team_id, started_at)
    st.update_team_status(team_id, "running")
    st.append_team_output(team_id, f"[AGENT SDK 시작: {team_type.upper()} {role_label}]\n")
    nl.log_event("agent_start", task[:80],
                 agent_id=team_id, project_id=project_id,
                 team_type=team_type, role=nl_role, model=model,
                 status="running", task=task[:300])

    if is_lead:
        system_prompt = LEAD_PROMPTS.get(team_type, LEAD_PROMPTS["research"]).format(
            agent_id=team_id,
            team_type=team_type,
            project_id=project_id,
        )
    else:
        base = WORKER_PROMPTS.get(team_type, WORKER_PROMPTS["research"]).format(
            agent_id=team_id,
            team_type=team_type,
            project_id=project_id,
            lead_id=lead_id or "",
        )
        system_prompt = base + _build_skill_hint(skills or [])

    prompt = (
        "아래 태스크를 지금 즉시 실행하라. "
        "추가 질문이나 확인 없이 바로 작업을 시작하라.\n\n"
        f"[태스크]\n{task}"
    )

    try:
        model = (LEAD_MODELS.get(team_type, "claude-sonnet-4-6")
                 if is_lead else
                 WORKER_MODELS.get(team_type, "claude-sonnet-4-6"))
        # ResultMessage 없이 종료될 때를 대비한 usage 누적 추적
        _acc_in_tok = 0
        _acc_out_tok = 0
        _acc_cost = 0.0
        _last_result_text = ""
        async for message in query(
            prompt=prompt,
            options=ClaudeAgentOptions(
                cwd=workdir,
                system_prompt=system_prompt,
                allowed_tools=ALLOWED_TOOLS,
                mcp_servers={"board": BOARD_MCP},
                permission_mode="bypassPermissions",
                env={"CLAUDECODE": ""},
                model=model,
                max_turns=50,
            ),
        ):
            if isinstance(message, AssistantMessage):
                for block in message.content:
                    if isinstance(block, TextBlock):
                        st.append_team_output(team_id, block.text)
                        _last_result_text = block.text
                # AssistantMessage에서 usage 누적
                if hasattr(message, "usage") and message.usage:
                    _acc_in_tok += _tok(message.usage, "input_tokens")
                    _acc_out_tok += _tok(message.usage, "output_tokens")

            elif isinstance(message, ResultMessage):
                if message.result:
                    st.append_team_output(team_id, f"\n[최종 결과]\n{message.result}")
                elapsed = time.time() - started_at
                in_tok = _tok(message.usage, "input_tokens")
                out_tok = _tok(message.usage, "output_tokens")
                cost = message.total_cost_usd or 0.0
                st.update_team_usage(
                    team_id,
                    input_tokens=in_tok,
                    output_tokens=out_tok,
                    model_used=model,
                    elapsed_seconds=elapsed,
                    total_cost_usd=cost,
                )
                st.update_team_status(team_id, "done")
                role_label = "리드" if is_lead else "워커"
                st.add_notification(
                    team_id,
                    f"✅ [{team_type}] {role_label} 완료 (소요 {elapsed:.0f}초, "
                    f"토큰 in={in_tok:,} out={out_tok:,})"
                )
                st.append_team_output(
                    team_id,
                    f"\n[완료 | 소요 {elapsed:.1f}초 | "
                    f"토큰 in={in_tok:,} out={out_tok:,} | "
                    f"모델={model}]\n"
                )
                nl_role2 = "리드" if is_lead else "워커"
                nl.log_event("agent_done", task[:80],
                             agent_id=team_id, project_id=project_id,
                             team_type=team_type, role=nl_role, model=model, status="done",
                             elapsed=elapsed, tokens_in=in_tok, tokens_out=out_tok, cost=cost,
                             task=task[:300], result=(message.result or "")[:300])
                _notify_completion(team_id, team_type, nl_role, "done",
                                   elapsed, in_tok, out_tok, cost,
                                   (message.result or "")[:800],
                                   project_id=project_id)
                return

        # ResultMessage 없이 스트림 종료 시 — 누적 usage 활용
        elapsed = time.time() - started_at
        st.update_team_usage(
            team_id,
            input_tokens=_acc_in_tok,
            output_tokens=_acc_out_tok,
            model_used=model,
            elapsed_seconds=elapsed,
            total_cost_usd=_acc_cost,
        )
        st.update_team_status(team_id, "done")
        st.add_notification(
            team_id,
            f"✅ [{team_type}] 완료 (소요 {elapsed:.0f}초, "
            f"토큰 in={_acc_in_tok:,} out={_acc_out_tok:,})"
        )
        st.append_team_output(
            team_id,
            f"\n[완료 | 소요 {elapsed:.1f}초 | "
            f"토큰 in={_acc_in_tok:,} out={_acc_out_tok:,} | "
            f"모델={model}]\n"
        )
        nl.log_event("agent_done", task[:80],
                     agent_id=team_id, project_id=project_id,
                     team_type=team_type, role=nl_role, model=model, status="done",
                     elapsed=elapsed, tokens_in=_acc_in_tok, tokens_out=_acc_out_tok,
                     cost=_acc_cost, task=task[:300],
                     result=_last_result_text[:300])
        _notify_completion(team_id, team_type, nl_role, "done",
                           elapsed, _acc_in_tok, _acc_out_tok, _acc_cost,
                           _last_result_text[:800],
                           project_id=project_id)

    except Exception as e:
        err_str = str(e).lower()
        if _retry == 0 and any(kw in err_str for kw in RETRY_ERRORS):
            delay = 30
            st.append_team_output(team_id, f"\n[API 오류 감지, {delay}초 후 재시도 (1/1)] {e}\n")
            await asyncio.sleep(delay)
            return await run_agent(team_id, team_type, task, project_id, workdir,
                                   is_lead=is_lead, lead_id=lead_id, skills=skills, _retry=1)
        st.update_team_status(team_id, "failed", error=str(e))
        st.add_notification(team_id, f"❌ [{team_type}] 실패: {str(e)[:80]}")
        st.append_team_output(team_id, f"\n[ERROR] {e}\n")
        nl.log_event("agent_failed", task[:80],
                     agent_id=team_id, project_id=project_id,
                     team_type=team_type, role=nl_role, model=model, status="failed",
                     task=task[:300], result=str(e)[:300])
        _notify_completion(team_id, team_type, nl_role, "failed",
                           time.time() - started_at, 0, 0, 0.0, str(e)[:800],
                           project_id=project_id)


async def run_solo_agent(team_id: str, team_type: str, task: str,
                         project_id: str, workdir: str, skills: list = None,
                         _retry: int = 0):
    """
    특무 에이전트 — board MCP 없이 단독 실행.
    리드/워커 계층 없음, 에이전트 간 소통 없음.
    단순·단일 태스크에 최적화.
    """
    from claude_agent_sdk import (
        query, ClaudeAgentOptions,
        ResultMessage, AssistantMessage, TextBlock,
    )

    st.init_db()
    solo_model = WORKER_MODELS.get(team_type, "claude-sonnet-4-6")
    started_at = time.time()
    st.set_team_started_at(team_id, started_at)
    st.update_team_status(team_id, "running")
    st.append_team_output(team_id, f"[특무 에이전트 시작: {team_type.upper()}]\n")
    nl.log_event("agent_start", task[:80],
                 agent_id=team_id, project_id=project_id,
                 team_type=team_type, role="solo", model=solo_model,
                 status="running", task=task[:300])

    system_prompt = (
        _SOLO_SYSTEM_BASE.format(agent_id=team_id, team_type=team_type, project_id=project_id)
        + _build_skill_hint(skills or [])
    )
    # worker_type 전문화 설명 추가
    if team_type in WORKER_SPECIALIZATIONS:
        desc, extra_skills = WORKER_SPECIALIZATIONS[team_type]
        system_prompt += f"\n\n## 전문 분야\n{desc}"
        if extra_skills and not skills:
            system_prompt += _build_skill_hint(extra_skills)

    prompt = (
        "아래 태스크를 지금 즉시 실행하라. "
        "추가 질문이나 확인 없이 바로 작업을 시작하라.\n\n"
        f"[태스크]\n{task}"
    )

    model = WORKER_MODELS.get(team_type, "claude-sonnet-4-6")

    try:
        # ResultMessage 없이 종료될 때를 대비한 usage 누적 추적
        _acc_in_tok = 0
        _acc_out_tok = 0
        _acc_cost = 0.0
        _last_result_text = ""
        async for message in query(
            prompt=prompt,
            options=ClaudeAgentOptions(
                cwd=workdir,
                system_prompt=system_prompt,
                allowed_tools=ALLOWED_TOOLS,
                mcp_servers={},          # board MCP 없음 — 소통 불필요
                permission_mode="bypassPermissions",
                env={"CLAUDECODE": ""},
                model=model,
                max_turns=30,
            ),
        ):
            if isinstance(message, AssistantMessage):
                for block in message.content:
                    if isinstance(block, TextBlock):
                        st.append_team_output(team_id, block.text)
                        _last_result_text = block.text
                # AssistantMessage에서 usage 누적
                if hasattr(message, "usage") and message.usage:
                    _acc_in_tok += _tok(message.usage, "input_tokens")
                    _acc_out_tok += _tok(message.usage, "output_tokens")

            elif isinstance(message, ResultMessage):
                if message.result:
                    st.append_team_output(team_id, f"\n[최종 결과]\n{message.result}")
                elapsed = time.time() - started_at
                in_tok = _tok(message.usage, "input_tokens")
                out_tok = _tok(message.usage, "output_tokens")
                cost = message.total_cost_usd or 0.0
                st.update_team_usage(
                    team_id,
                    input_tokens=in_tok,
                    output_tokens=out_tok,
                    model_used=model,
                    elapsed_seconds=elapsed,
                    total_cost_usd=cost,
                )
                st.update_team_status(team_id, "done")
                st.add_notification(
                    team_id,
                    f"✅ [특무/{team_type}] 완료 (소요 {elapsed:.0f}초, "
                    f"토큰 in={in_tok:,} out={out_tok:,})"
                )
                st.append_team_output(
                    team_id,
                    f"\n[완료 | 소요 {elapsed:.1f}초 | "
                    f"토큰 in={in_tok:,} out={out_tok:,} | "
                    f"모델={model}]\n"
                )
                nl.log_event("agent_done", task[:80],
                             agent_id=team_id, project_id=project_id,
                             team_type=team_type, role="solo", model=solo_model, status="done",
                             elapsed=elapsed, tokens_in=in_tok, tokens_out=out_tok, cost=cost,
                             task=task[:300], result=(message.result or "")[:300])
                _notify_completion(team_id, team_type, "solo", "done",
                                   elapsed, in_tok, out_tok, cost,
                                   (message.result or "")[:800],
                                   project_id=project_id)
                return

        # ResultMessage 없이 스트림 종료 시 — 누적 usage 활용
        elapsed = time.time() - started_at
        st.update_team_usage(
            team_id,
            input_tokens=_acc_in_tok,
            output_tokens=_acc_out_tok,
            model_used=model,
            elapsed_seconds=elapsed,
            total_cost_usd=_acc_cost,
        )
        st.update_team_status(team_id, "done")
        st.add_notification(
            team_id,
            f"✅ [특무/{team_type}] 완료 (소요 {elapsed:.0f}초, "
            f"토큰 in={_acc_in_tok:,} out={_acc_out_tok:,})"
        )
        st.append_team_output(
            team_id,
            f"\n[완료 | 소요 {elapsed:.1f}초 | "
            f"토큰 in={_acc_in_tok:,} out={_acc_out_tok:,} | "
            f"모델={model}]\n"
        )
        nl.log_event("agent_done", task[:80],
                     agent_id=team_id, project_id=project_id,
                     team_type=team_type, role="solo", model=solo_model, status="done",
                     elapsed=elapsed, tokens_in=_acc_in_tok, tokens_out=_acc_out_tok,
                     cost=_acc_cost, task=task[:300],
                     result=_last_result_text[:300])
        _notify_completion(team_id, team_type, "solo", "done",
                           elapsed, _acc_in_tok, _acc_out_tok, _acc_cost,
                           _last_result_text[:800],
                           project_id=project_id)

    except Exception as e:
        err_str = str(e).lower()
        if _retry == 0 and any(kw in err_str for kw in RETRY_ERRORS):
            delay = 30
            st.append_team_output(team_id, f"\n[API 오류 감지, {delay}초 후 재시도 (1/1)] {e}\n")
            await asyncio.sleep(delay)
            return await run_solo_agent(team_id, team_type, task, project_id, workdir,
                                        skills=skills, _retry=1)
        st.update_team_status(team_id, "failed", error=str(e))
        st.add_notification(team_id, f"❌ [특무/{team_type}] 실패: {str(e)[:80]}")
        st.append_team_output(team_id, f"\n[ERROR] {e}\n")
        nl.log_event("agent_failed", task[:80],
                     agent_id=team_id, project_id=project_id,
                     team_type=team_type, role="solo", model=solo_model, status="failed",
                     task=task[:300], result=str(e)[:300])
        _notify_completion(team_id, team_type, "solo", "failed",
                           time.time() - started_at, 0, 0, 0.0, str(e)[:800],
                           project_id=project_id)
