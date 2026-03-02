---
name: asana-extended-api
description: |
  Asana MCP 도구에 없는 REST API를 Python으로 직접 호출하는 확장 스킬.
  태스크를 여러 프로젝트에 동시 소속(multi-homing), 섹션 생성/관리, 태그 추가/제거,
  포트폴리오 관리, 프로젝트 복제, 배치 작업 등을 지원한다.

  다음 상황에서 반드시 이 스킬을 사용할 것:
  - "태스크를 프로젝트에 추가해줘", "이 태스크를 팀원관리에도 넣어줘" 등 태스크 multi-homing 요청
  - "섹션 만들어줘", "섹션에 태스크 이동해줘" 등 섹션 관련 요청
  - "태그 달아줘", "태그 만들어줘" 등 태그 관련 요청
  - "프로젝트 복제해줘", "포트폴리오에 추가해줘" 등 MCP에서 지원하지 않는 작업
  - "여러 태스크를 한번에 옮겨줘" 등 배치/일괄 작업 요청
  - addProject, removeProject, createSection, addTaskToSection, addTag, removeTag,
    duplicateProject, batch 등의 API 키워드가 언급될 때
---

# Asana Extended API

Asana MCP 도구가 커버하지 않는 REST API를 직접 호출한다.

## 실행 방법 (우선순위)

### 방법 1: Chrome JS (권장, PAT 불필요)

브라우저에서 Asana 로그인 세션을 활용. **핵심 헤더: `X-Allow-Asana-Client: 1`**

```javascript
(async () => {
  const resp = await fetch('https://app.asana.com/api/1.0/tasks/TASK_GID/addProject?opt_silent=true', {
    method: 'POST',
    credentials: 'include',
    headers: {
      'Content-Type': 'application/json',
      'Accept': 'application/json',
      'X-Requested-With': 'XMLHttpRequest',
      'X-Allow-Asana-Client': '1',
    },
    body: JSON.stringify({data: {project: 'PROJECT_GID'}})
  });
  return JSON.stringify({status: resp.status});
})();
```

**주의**: GET 요청은 `credentials: 'include'`만으로 동작하지만,
POST/PUT/DELETE는 반드시 `X-Allow-Asana-Client: 1` 헤더가 있어야 401이 아닌 200을 반환한다.

배치 처리 시 `await new Promise(r => setTimeout(r, 300))`으로 300ms 간격 유지.

### 방법 2: Python + PAT (MCP 서버 또는 CLI)

환경변수 `ASANA_PAT`가 설정된 경우 `scripts/asana_api.py` 헬퍼 또는 MCP 서버 사용.

```
ASANA_PAT 발급: https://app.asana.com/0/my-apps → Personal access token → New
설정: setx ASANA_PAT "2/xxxx/xxxx:xxxx"  (Windows)
```

`scripts/asana_api.py` 헬퍼 모듈 또는 `scripts/mcp_server.py` (stdio MCP 서버)로 호출.

## 지원 API

### 1. Task → Project (multi-homing)

태스크를 추가 프로젝트에 소속시키거나 제거한다. 원본 프로젝트에서 삭제되지 않는다.

```python
# 태스크를 프로젝트에 추가
POST /tasks/{task_gid}/addProject
Body: {"data": {"project": "<project_gid>", "section": "<section_gid>(optional)"}}

# 태스크를 프로젝트에서 제거
POST /tasks/{task_gid}/removeProject
Body: {"data": {"project": "<project_gid>"}}
```

**사용 예:** 김유담의 "AcP Synthesis" 태스크를 Nucleoside APIs 프로젝트에 두면서 팀원관리 프로젝트에도 추가.

### 2. Section 관리

프로젝트 내 섹션을 생성/수정/삭제하고, 태스크를 특정 섹션으로 이동한다.

```python
# 섹션 생성
POST /projects/{project_gid}/sections
Body: {"data": {"name": "섹션 이름"}}

# 섹션에 태스크 추가 (다른 섹션에서 제거됨)
POST /sections/{section_gid}/addTask
Body: {"data": {"task": "<task_gid>"}}

# 섹션 삭제 (빈 섹션만 가능)
DELETE /sections/{section_gid}
```

### 3. Tag 관리

태그를 만들고 태스크에 부착/제거한다.

```python
# 태그 생성
POST /workspaces/{workspace_gid}/tags
Body: {"data": {"name": "태그명"}}

# 태스크에 태그 추가
POST /tasks/{task_gid}/addTag
Body: {"data": {"tag": "<tag_gid>"}}

# 태스크에서 태그 제거
POST /tasks/{task_gid}/removeTag
Body: {"data": {"tag": "<tag_gid>"}}
```

### 4. Portfolio 관리

포트폴리오를 생성하고 프로젝트를 추가/제거한다.

```python
# 포트폴리오에 항목 추가
POST /portfolios/{portfolio_gid}/addItem
Body: {"data": {"item": "<project_gid>"}}

# 포트폴리오에서 항목 제거
POST /portfolios/{portfolio_gid}/removeItem
Body: {"data": {"item": "<project_gid>"}}
```

### 5. Project 복제

프로젝트를 태스크 포함하여 비동기 복제한다.

```python
POST /projects/{project_gid}/duplicate
Body: {"data": {"name": "복제된 프로젝트", "include": ["members","task_notes","task_assignee","task_subtasks","task_attachments","task_dates","task_dependencies","task_followers","task_tags","task_projects"]}}
```

### 6. Batch 작업

여러 API 호출을 하나의 요청으로 묶는다 (최대 약 10개).

```python
POST /batch
Body: {
  "data": {
    "actions": [
      {
        "method": "POST",
        "relative_path": "/tasks/TASK_GID_1/addProject",
        "data": {"project": "PROJECT_GID"}
      },
      {
        "method": "POST",
        "relative_path": "/tasks/TASK_GID_2/addProject",
        "data": {"project": "PROJECT_GID"}
      }
    ]
  }
}
```

## 스크립트 생성 패턴

사용자 요청에 맞는 Python 스크립트를 workspace 폴더에 생성한다.
아래 `scripts/asana_api.py` 모듈을 import하여 사용하면 편리하다.

```python
# 사용 예시
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from asana_api import AsanaAPI

api = AsanaAPI()  # ASANA_PAT 환경변수에서 자동 로드

# 태스크를 프로젝트에 추가
api.add_task_to_project("TASK_GID", "PROJECT_GID")

# 섹션 생성
section = api.create_section("PROJECT_GID", "새 섹션")

# 태스크를 섹션으로 이동
api.add_task_to_section("SECTION_GID", "TASK_GID")

# 일괄 addProject
task_gids = ["GID1", "GID2", "GID3"]
api.batch_add_to_project(task_gids, "PROJECT_GID")
```

## 주의사항

- Asana Free 플랜에서도 모든 REST API가 동작한다 (search_tasks만 Premium 필요).
- Rate limit: 분당 1500 요청. batch 내 각 action은 1 요청으로 계산된다.
- multi-homing: 태스크당 최대 20개 프로젝트에 소속 가능.
- 섹션 삭제는 빈 섹션만 가능하고, 프로젝트의 마지막 섹션은 삭제할 수 없다.
