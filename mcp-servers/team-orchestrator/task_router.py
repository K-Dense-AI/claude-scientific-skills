"""프로젝트 설명을 분석해 팀 배정 결정 - 키워드 기반 + LLM fallback"""
import json
import os
import re
import subprocess
from pathlib import Path

SKILL_BASE = str(Path.home() / "claude-scientific-skills" / "scientific-skills")

# 팀 유형별 기본 스킬 목록 (키워드 분류 시 상위 3개 주입)
TEAM_DEFAULT_SKILLS: dict[str, list[str]] = {
    "research":        ["pubmed-database", "biorxiv-database", "openalex-database"],
    "bioinformatics":  ["biopython", "kegg-database", "pdb-database"],
    "data-analysis":   ["statsmodels", "matplotlib", "scikit-learn"],
    "code":            ["cobrapy", "biopython", "rdkit"],
    "writing":         ["manuscript-writer", "citation-management", "scientific-writing"],
    "lab-protocol":    ["experiment-hub", "research-commons", "hypothesis-generation"],
    "ops":             ["git-workflow-manager", "conda-env-manager", "asana-extended-api"],
    "planning":        ["hypothesis-generation", "scientific-brainstorming", "research-commons"],
    "sci-review":      ["scientific-critical-thinking", "peer-review", "literature-review"],
    "code-review":     ["solid-principles", "code-validator", "git-workflow-manager"],
}

ROUTER_PROMPT = """당신은 프로젝트 관리 전문가입니다.
주어진 프로젝트 설명을 분석해 아래 팀 유형 중 적합한 것에 서브태스크를 배정하세요.

팀 유형 (10가지):
- research: 문헌 검색, 학술 DB 조회, 정보 수집, 문헌 리뷰
- bioinformatics: 서열 분석, 단백질 구조, 오믹스(scRNA-seq/bulk), 대사경로 모델링
- data-analysis: 통계 분석, 머신러닝, 데이터 시각화, EDA
- code: 코드 구현, 버그 수정, 테스트, 리팩토링, 파이프라인 개발
- writing: 논문 초안, 보고서, 발표자료, 포스터 제작
- lab-protocol: 실험 설계, 조건 최적화, 프라이머 설계, 프로토콜 작성
- ops: git 관리, Asana 업데이트, conda 환경 설정, 배포
- planning: 작업 계획 수립, 범위 명확화, 태스크 분해, 실현 가능성 검토, 리스크 평가
- sci-review: 과학적 검토, 방법론 타당성, 선행연구 비교, 주장 비판, 연구 격차 분석
- code-review: 코드 리뷰, 보안 취약점 스캔, SOLID 원칙 검토, 성능 분석, 테스트 커버리지

각 팀에서 쓸 수 있는 스킬 (1-3개 선택):
research: pubmed-database, biorxiv-database, openalex-database, kegg-database, uniprot-database, brenda-database, gene-database, clinicaltrials-database, chembl-database
bioinformatics: biopython, gget, ensembl-database, pdb-database, alphafold-database, esm, scanpy, pydeseq2, scvi-tools, anndata, geo-database, kegg-database, reactome-database, cobrapy, string-database
data-analysis: statsmodels, pymc, matplotlib, seaborn, plotly, scientific-visualization, scikit-learn, pytorch-lightning, transformers, shap, umap-learn
code: cobrapy, biopython, rdkit, biosteam, scikit-learn, pytorch-lightning
writing: manuscript-writer, citation-management, scientific-writing, scientific-visualization, latex-posters, pptx-reviewer, scientific-slides
lab-protocol: experiment-hub, research-commons, hypothesis-generation, ipcr-primer-design, primer-design
ops: git-workflow-manager, conda-env-manager, asana-extended-api
planning: hypothesis-generation, scientific-brainstorming, research-commons
sci-review: scientific-critical-thinking, peer-review, literature-review, pubmed-database, openalex-database
code-review: solid-principles, code-validator, git-workflow-manager

배정 규칙:
1. 독립적으로 병렬 실행 가능한 것은 각각 다른 팀으로
2. 한 팀의 결과가 다른 팀에 필요하면 dependency 명시 (depends_on에 팀 이름 기입)
3. 같은 유형 작업 여러 개는 하나의 팀으로 합치기
4. 팀당 작업은 구체적이고 실행 가능하게 (2-3문장)
5. skills: 해당 태스크에 가장 적합한 스킬 1-3개 (없으면 빈 배열)

코드블록 없이 순수 JSON만 응답하세요:
{
  "project_summary": "한 줄 요약",
  "teams": [
    {
      "name": "팀 이름",
      "type": "research|bioinformatics|data-analysis|code|writing|lab-protocol|ops",
      "task": "구체적인 작업 내용 (2-3문장)",
      "depends_on": [],
      "skills": ["skill-name-1", "skill-name-2"]
    }
  ]
}"""


RESEARCH_KEYWORDS = [
    "조사", "분석", "논문", "DB", "찾기", "리뷰", "문헌", "수집", "파악", "확인",
    "mutation", "kinetic", "keq", "ref", "산분해",
]
BIOINFORMATICS_KEYWORDS = [
    "sequence", "서열", "alignment", "blast", "clustering", "structural", "structure",
    "단백질", "protein", "genome", "유전체", "transcriptome", "RNA-seq", "scRNA",
    "pathway", "경로", "대사", "metabolic", "FBA", "PCR", "cloning", "primer",
    "eta", "루트", "AlphaFold", "PDB", "omics", "오믹스",
]
DATA_ANALYSIS_KEYWORDS = [
    "통계", "statistical", "regression", "머신러닝", "machine learning", "ML",
    "시각화", "visualization", "plot", "graph", "figure", "EDA", "clustering",
    "Standard curve", "BO", "학습셋", "모델", "예측", "classification",
]
CODE_KEYWORDS = [
    "코드", "구현", "버그", "스크립트", "레포", "추가", "입력", "표시", "모듈",
    "PeakPicker", "Kinetic-modeling", "UDH_Clustering", "biosteam",
    "pH", "브랜치", "Tracker", "파이프라인",
]
LAB_PROTOCOL_KEYWORDS = [
    "실험", "프로토콜", "protocol", "조건", "최적화", "optimization",
    "프라이머", "primer", "iPCR", "cloning", "설계", "샘플",
]
OPS_KEYWORDS = ["git", "Asana", "환경", "설정", "배포", "하드코딩", "스캔", "커밋", "푸시", "conda"]
WRITING_KEYWORDS = ["문서", "보고서", "논문 초안", "발표", "작성", "PFD", "포스터", "슬라이드"]
PLANNING_KEYWORDS = [
    "계획", "플래닝", "planning", "설계", "범위", "scope", "roadmap", "로드맵",
    "전략", "방향", "우선순위", "어떻게", "무엇을", "방법론", "아키텍처",
    "태스크 분해", "분해", "실현 가능성", "feasibility", "리스크", "risk",
]
SCI_REVIEW_KEYWORDS = [
    "검토", "리뷰", "review", "타당성", "과학적", "scientific",
    "통계 검토", "논리", "근거", "evidence", "claim", "주장", "비판",
    "peer review", "선행연구 비교", "문헌 비교", "방법론 검토", "실험 검토",
]
CODE_REVIEW_KEYWORDS = [
    "코드 리뷰", "code review", "security", "보안", "취약점", "SOLID",
    "리팩토링 검토", "성능 검토", "버그 탐지", "품질", "quality",
    "테스트 커버리지", "coverage", "코드 품질", "오와스프", "owasp",
]

CLAUDE_BIN = str(Path.home() / ".local" / "bin" / "claude.exe")


_KEYWORD_MAP = [
    ("research",       RESEARCH_KEYWORDS),
    ("bioinformatics", BIOINFORMATICS_KEYWORDS),
    ("data-analysis",  DATA_ANALYSIS_KEYWORDS),
    ("code",           CODE_KEYWORDS),
    ("lab-protocol",   LAB_PROTOCOL_KEYWORDS),
    ("ops",            OPS_KEYWORDS),
    ("writing",        WRITING_KEYWORDS),
    ("planning",       PLANNING_KEYWORDS),
    ("sci-review",     SCI_REVIEW_KEYWORDS),
    ("code-review",    CODE_REVIEW_KEYWORDS),
]


def _score_task(task_text: str) -> dict:
    txt = task_text.lower()
    scores = {t: 0 for t, _ in _KEYWORD_MAP}
    for team_type, keywords in _KEYWORD_MAP:
        for kw in keywords:
            if kw.lower() in txt:
                scores[team_type] += 1
    return scores


def _classify_task(task_text: str) -> str:
    scores = _score_task(task_text)
    return max(scores, key=scores.get)


def _is_ambiguous(tasks: list) -> bool:
    """키워드 점수가 모두 0이거나 최고점이 동점인 태스크가 있으면 True"""
    for task in tasks:
        scores = _score_task(task)
        max_score = max(scores.values())
        if max_score == 0:
            return True
        top_types = [t for t, s in scores.items() if s == max_score]
        if len(top_types) > 1:
            return True
    return False


def _classify_with_llm(description: str) -> dict | None:
    """LLM(claude.exe)으로 태스크 분류 — 키워드 기반 모호 시 fallback"""
    prompt = f"{ROUTER_PROMPT}\n\n[프로젝트 설명]\n{description}"
    clean_env = {k: v for k, v in os.environ.items() if k != "CLAUDECODE"}
    try:
        result = subprocess.run(
            [CLAUDE_BIN, "-p", prompt, "--dangerously-skip-permissions",
             "--model", "claude-opus-4-6"],
            capture_output=True, timeout=60, env=clean_env,
        )
        result = type('r', (), {
            'returncode': result.returncode,
            'stdout': result.stdout.decode('utf-8', errors='replace'),
        })()
        if result.returncode != 0:
            return None
        output = result.stdout.strip()
        json_match = re.search(r"\{.*\}", output, re.DOTALL)
        if json_match:
            return json.loads(json_match.group())
    except Exception:
        pass
    return None


def analyze_and_split(description: str) -> dict:
    """
    프로젝트 설명 -> 팀 배정 계획 반환
    - 명확한 키워드가 있으면 키워드 기반 분류 (빠름)
    - 모호하면 LLM fallback (정확, depends_on 자동 생성 포함)
    Returns: { project_summary, teams: [{name, type, task, depends_on}] }
    """
    lines = description.strip().split("\n")

    tasks = []
    current_task = []
    for line in lines:
        if line.strip().startswith("[") and "]" in line:
            if current_task:
                tasks.append("\n".join(current_task))
            current_task = [line]
        elif current_task:
            current_task.append(line)
    if current_task:
        tasks.append("\n".join(current_task))

    # [태스크] 패턴이 없으면 description 전체를 하나의 태스크로 취급
    if not tasks:
        tasks = [description.strip()]

    # 모호한 태스크가 있으면 LLM fallback 시도
    if _is_ambiguous(tasks):
        llm_result = _classify_with_llm(description)
        if llm_result and llm_result.get("teams"):
            return llm_result

    # 키워드 기반 분류 — 태스크당 별도 팀 생성 (같은 유형도 병렬 실행)
    import re as _re
    type_names = {
        "research": "Research", "code": "Dev", "ops": "Ops", "writing": "Writer",
        "bioinformatics": "BioInfo", "data-analysis": "DataAna", "lab-protocol": "LabProto",
        "planning": "Planning", "sci-review": "SciReview", "code-review": "CodeReview",
    }
    type_counters: dict = {}
    teams = []
    for task_block in tasks:
        t = _classify_task(task_block)
        type_counters[t] = type_counters.get(t, 0) + 1
        # 태스크 첫 줄에서 이름 추출 (예: "[리서치1: xxx]" → "리서치1")
        first_line = task_block.strip().splitlines()[0] if task_block.strip() else ""
        label_match = _re.match(r"\[([^\]]{1,30})", first_line)
        name = label_match.group(1).strip() if label_match else f"{type_names[t]}{type_counters[t]}"
        teams.append({
            "name": name,
            "type": t,
            "task": task_block.strip(),
            "depends_on": [],
            "skills": TEAM_DEFAULT_SKILLS.get(t, [])[:3],
        })

    return {
        "project_summary": f"Lab sprint: {len(tasks)} tasks distributed to {len(teams)} agents in parallel",
        "teams": teams,
    }
