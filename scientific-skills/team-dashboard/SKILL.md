---
name: team-dashboard
description: "Ribose 팀원별 Asana 현황 조회 및 디스커션 연동 스킬. (1) '김유담 현황' 같이 팀원명 언급 시 Asana 태스크 자동 조회, (2) 프로젝트별 진행도/마감 현황, (3) 관련 프로젝트 참조자료 크로스레퍼런스, (4) 디스커션 중 실시간 현황 참조. '팀원 현황', '진행도', '누가 뭐하고 있어', '마감 임박' 등 트리거."
---

# Team Dashboard (팀원별 Asana 현황 조회)

Ribose team Asana 워크스페이스에서 팀원별·프로젝트별 현황을 조회하고, 디스커션/미팅에서 실시간 참조하는 스킬.

## MCP 도구 사용 방법

**MCP 도구는 메인 모델이 직접 호출한다. 서브에이전트에 위임하지 않는다.**
(서브에이전트는 MCP 도구에 접근 불가)

사용 가능한 MCP 도구 (`asana-extended` 서버):
- `asana_get_tasks(project_id, opt_fields, limit)` — 프로젝트 태스크 목록
- `asana_get_task(task_gid, opt_fields)` — 태스크 상세 (notes 포함)
- `asana_get_project_task_counts(project_id)` — 완료율
- `asana_search_tasks_in_workspace(workspace_gid, text, opt_fields)` — 키워드 검색

## Workspace & Team 정보

```
Workspace: Ribose team (GID: 1202946172929462)
```

### 팀원 매핑 (이름 → GID)

```
이자현 (PI)    → 1202946172861105  jahyunlee082@gmail.com
김유담         → 1206785467601510  fogree@korea.ac.kr
김지원         → 1206415808608485  gogokiru29@naver.com
장서현         → 1209566337407843  tjgus1098@naver.com
김호준         → 1211203914550274  hojun000331@naver.com
윤주           → 1202947125854982  wooknjoo50@gmail.com
강임수         → 1203628275591657  rkddlatn7473@naver.com
박성현         → 1205439025541060  delate2137@naver.com
김서혁         → 1205439025542118  kshpiano116@naver.com
신채연         → 1205439025541058  kiv0818@naver.com
장정화         → 1205432756567100  wjdghk0242@naver.com
표건           → 1205439025541054  gen6773@naver.com
lilik krisna   → 1205439025542122  lilikkrisnamukti@gmail.com
박준하         → 1203628275610421  ahee4747@naver.com
이주현         → 1204954542012371  leeju01@korea.ac.kr
김선희         → 1205439025541056  33shkim@naver.com
```

### 프로젝트 매핑

```
주요 활성 프로젝트:
├─ Enzymatic Nucleoside APIs       → 1203520533771201  (김유담 중심)
├─ D-Gal to DL-Tag               → 1207549530067901  (김지원 중심)
├─ [한미] Aldaric acids           → 1209503367534283  (장서현+김호준)
├─ Cloning and Strains            → 1205504956716782  (다수 팀원)
├─ L-Ribose production            → 1203828166717600  (이자현)
├─ D-Ribose production            → 1202968410293444  (이자현)
├─ Methods & maintenance          → 1202946550171010  (공용)

참조용 프로젝트:
├─ [경북대] Troubleshooting        → 1205432756567083  (장정화)
├─ 고려대_학부인턴                  → 1206415808608498  (인턴)
└─ Python                         → 1206118091459210  (교육)
```

---

## MCP 미연결 시 폴백: Chrome JS

Asana MCP 도구(`asana_get_tasks` 등)가 현재 세션에 없거나 ASANA_PAT가 없을 때,
**브라우저 로그인 세션을 활용한 Chrome JS 스크립트**를 생성해 사용자에게 제공한다.

### 폴백 감지 조건

MCP 도구 목록에 `asana_get_tasks`가 없으면 → Chrome JS 방식으로 전환

### Chrome JS 공통 패턴

```javascript
// 핵심: credentials:'include' + X-Allow-Asana-Client:'1' 헤더
const r = await fetch(URL, {
  credentials: 'include',
  headers: { 'X-Allow-Asana-Client': '1', 'Accept': 'application/json' }
});
const json = await r.json();
```

- Asana가 로그인된 브라우저에서 F12 → Console에 붙여넣어 실행
- PAT 불필요, 브라우저 세션 쿠키로 인증
- GET 요청은 헤더 없이도 동작, POST는 반드시 `X-Allow-Asana-Client: 1` 필요

### 폴백 스크립트 템플릿

#### M1 폴백: 팀원 현황 조회

```javascript
// {팀원명} 현황 조회 — Asana 콘솔(F12)에서 실행
(async () => {
  const PROJECTS = [{gid: '{PROJECT_GID}', name: '{PROJECT_NAME}'}]; // 팀원의 주요 프로젝트
  const ASSIGNEE_GID = '{ASSIGNEE_GID}'; // 팀원 GID
  const FIELDS = 'name,assignee.name,assignee.gid,due_on,completed,completed_at,created_at,gid';
  const today = new Date().toISOString().slice(0,10);
  const week = new Date(Date.now()+7*864e5).toISOString().slice(0,10);

  let all = [];
  for (const p of PROJECTS) {
    let offset = null;
    do {
      let url = `https://app.asana.com/api/1.0/projects/${p.gid}/tasks?opt_fields=${FIELDS}&limit=100`;
      if (offset) url += `&offset=${offset}`;
      const j = await (await fetch(url, {credentials:'include', headers:{'X-Allow-Asana-Client':'1'}})).json();
      (j.data||[]).filter(t => t.assignee?.gid === ASSIGNEE_GID).forEach(t => all.push({...t, _project: p.name}));
      offset = j.next_page?.offset || null;
    } while (offset);
  }

  const overdue = all.filter(t => !t.completed && t.due_on && t.due_on < today);
  const thisWeek = all.filter(t => !t.completed && t.due_on && t.due_on >= today && t.due_on <= week);
  const inProgress = all.filter(t => !t.completed && !t.due_on);
  const recentDone = all.filter(t => t.completed && t.completed_at >= new Date(Date.now()-7*864e5).toISOString());

  console.log('🔴 마감 지연:', overdue.map(t=>({name:t.name, due:t.due_on, project:t._project})));
  console.log('🟡 이번주 마감:', thisWeek.map(t=>({name:t.name, due:t.due_on})));
  console.log('🟢 진행중:', inProgress.map(t=>({name:t.name, project:t._project})));
  console.log('✅ 최근완료:', recentDone.map(t=>({name:t.name, done:t.completed_at?.slice(0,10)})));
  console.log(`완료율: ${all.filter(t=>t.completed).length}/${all.length}`);
})();
```

#### M3 폴백: 키워드 검색

```javascript
// {키워드} 검색 — Asana 콘솔(F12)에서 실행
(async () => {
  const PROJECTS = [
    {gid: '1202946550171010', name: 'Methods & maintenance'},
    {gid: '1205504956716782', name: 'Cloning and Strains'},
    {gid: '1203520533771201', name: 'Enzymatic Nucleoside APIs'},
    // 추가 프로젝트 GID 삽입
  ];
  const KEYWORDS = ['{키워드1}', '{키워드2}'];
  const FIELDS = 'name,assignee.name,completed,completed_at,created_at,gid,notes';

  const results = [];
  for (const p of PROJECTS) {
    let offset = null;
    do {
      let url = `https://app.asana.com/api/1.0/projects/${p.gid}/tasks?opt_fields=${FIELDS}&limit=100`;
      if (offset) url += `&offset=${offset}`;
      const j = await (await fetch(url, {credentials:'include', headers:{'X-Allow-Asana-Client':'1'}})).json();
      (j.data||[]).forEach(t => {
        const hay = (t.name + ' ' + (t.notes||'')).toLowerCase();
        if (KEYWORDS.some(k => hay.includes(k.toLowerCase())))
          results.push({project: p.name, gid: t.gid, name: t.name,
            assignee: t.assignee?.name||'-', completed: t.completed,
            date: (t.completed_at||t.created_at||'').slice(0,10), notes: (t.notes||'').slice(0,200)});
      });
      offset = j.next_page?.offset || null;
    } while (offset);
  }
  console.table(results);
  return JSON.stringify(results, null, 2);
})();
```

#### M4 폴백: 마감 임박 전체 조회

```javascript
// 전체 팀 마감 임박 조회 — Asana 콘솔(F12)에서 실행
(async () => {
  const PROJECTS = [
    {gid:'1203520533771201', name:'Enzymatic Nucleoside APIs'},
    {gid:'1207549530067901', name:'D-Gal to DL-Tag'},
    {gid:'1209503367534283', name:'Aldaric acids'},
    {gid:'1205504956716782', name:'Cloning and Strains'},
    {gid:'1203828166717600', name:'L-Ribose'},
    {gid:'1202968410293444', name:'D-Ribose'},
  ];
  const today = new Date().toISOString().slice(0,10);
  const week = new Date(Date.now()+7*864e5).toISOString().slice(0,10);
  const two = new Date(Date.now()+14*864e5).toISOString().slice(0,10);
  const FIELDS = 'name,assignee.name,due_on,completed,gid';

  const overdue=[], thisWeek=[], nextWeek=[];
  for (const p of PROJECTS) {
    const j = await (await fetch(
      `https://app.asana.com/api/1.0/projects/${p.gid}/tasks?opt_fields=${FIELDS}&limit=100`,
      {credentials:'include', headers:{'X-Allow-Asana-Client':'1'}}
    )).json();
    (j.data||[]).filter(t=>!t.completed && t.due_on).forEach(t => {
      const row = {project:p.name, name:t.name, due:t.due_on, assignee:t.assignee?.name||'-'};
      if (t.due_on < today) overdue.push(row);
      else if (t.due_on <= week) thisWeek.push(row);
      else if (t.due_on <= two) nextWeek.push(row);
    });
  }
  console.log('🔴 마감 지연:'); console.table(overdue);
  console.log('🟡 이번 주:'); console.table(thisWeek);
  console.log('🔵 다음 주:'); console.table(nextWeek);
})();
```

### 폴백 적용 규칙

1. MCP 도구 없음 감지 → "MCP 미연결 상태입니다. Chrome JS로 조회합니다." 안내
2. 요청 내용(팀원명/키워드/기능)에 맞게 위 템플릿의 `{변수}` 부분을 채워서 제공
3. 스크립트 실행 결과를 사용자가 붙여넣으면 → M1/M3/M4 출력 형식으로 파싱·정리
4. 팀원 GID는 위 "팀원 매핑" 테이블에서 조회, 프로젝트 GID는 "프로젝트 매핑" 테이블에서 조회

---

## 핵심 기능

### M1: 팀원 현황 조회

#### 트리거
- "{이름} 현황", "{이름} 뭐하고 있어?"
- "김유담 태스크", "장서현 진행도"
- "팀원 현황 전체"

#### 실행 절차

```
1. 팀원명 → GID 매핑 (위 테이블)
2. Asana API 호출:
   ├─ asana_get_tasks(project={관련 프로젝트}, opt_fields=name,assignee.name,due_on,completed,completed_at,modified_at,permalink_url)
   │  → 해당 팀원의 assignee 태스크 필터링
   ├─ 여러 프로젝트 순회 (팀원의 주요 프로젝트 우선)
   └─ 미완료 태스크 중 due_on 기준 정렬
3. 출력:
   ├─ 진행 중(미완료) 태스크 목록 + 마감일
   ├─ 최근 완료 태스크 (7일 이내)
   ├─ 마감 지연 태스크 (due_on < today & !completed)
   └─ 완료율 (완료/전체)
```

#### 출력 형식

```
═══════════════════════════════════════════════════
        {팀원명} 현황 ({date})
═══════════════════════════════════════════════════
완료율: {completed}/{total} ({rate}%)

🔴 마감 지연:
├─ [{task_name}]({permalink_url}) (마감: {due_on}, {days}일 지연)
└─ [{task_name}]({permalink_url}) (마감: {due_on}, {days}일 지연)

🟡 이번 주 마감:
├─ [{task_name}]({permalink_url}) (마감: {due_on})
└─ [{task_name}]({permalink_url}) (마감: {due_on})

🟢 진행 중:
├─ [{task_name}]({permalink_url}) [{project_name}]
└─ [{task_name}]({permalink_url}) [{project_name}]

✅ 최근 완료 (7일):
├─ [{task_name}]({permalink_url}) (완료: {completed_at})
└─ [{task_name}]({permalink_url}) (완료: {completed_at})
═══════════════════════════════════════════════════
```

### M2: 프로젝트 진행도 조회

#### 트리거
- "{프로젝트명} 진행도", "Nucleoside 현황"
- "프로젝트별 진행도", "전체 진행률"

#### 실행 절차

```
1. 프로젝트명 → GID 매핑
2. asana_get_project_task_counts(project_id, opt_fields=num_tasks,num_incomplete_tasks,num_completed_tasks)
3. asana_get_tasks(project, opt_fields=...) → 팀원별 분류
4. 출력: 프로젝트 진행률 + 팀원별 기여도
```

### M3: 참조자료 크로스레퍼런스

#### 트리거
- 디스커션 중 method, protocol, 기존 결과 언급 시
- "{프로젝트}에서 {키워드} 관련 태스크 찾아줘"
- "인턴이 했던 {실험} 결과"

#### 실행 절차

```
1. 키워드 추출
2. 관련 프로젝트 태스크 검색:
   ├─ Methods & maintenance (공용 프로토콜)
   ├─ 고려대_학부인턴 (기초 실험 데이터)
   ├─ [경북대] (트러블슈팅 기록)
   └─ 해당 프로젝트 내 관련 태스크
3. 태스크 상세 조회 (asana_get_task → notes/html_notes)
4. 관련 정보 요약 반환
```

#### 출력 형식

```
📎 참조자료: "{키워드}"
──────────────────────────────────────
출처: {project_name} > [{task_name}]({permalink_url})
담당: {assignee} | 완료: {date}
내용: {notes 요약}
──────────────────────────────────────
```

### M4: 마감 임박 알림

#### 트리거
- "마감 임박", "이번 주 마감"
- "지연된 태스크", "overdue"

#### 실행 절차

```
1. 전체 활성 프로젝트 순회
2. 미완료 태스크 중 due_on 필터:
   ├─ 🔴 지연 (due_on < today)
   ├─ 🟡 이번 주 (today ≤ due_on ≤ today+7)
   └─ 🔵 다음 주 (today+7 < due_on ≤ today+14)
3. 팀원별 그룹핑 후 출력
```

---

## research-discussion 연동

### 디스커션 중 자동 호출 규칙

디스커션(research-discussion) 진행 중 다음 조건에서 team-dashboard가 자동 호출된다:

| 디스커션 상황 | 자동 호출 | 용도 |
|--------------|----------|------|
| 팀원명 + 실험 결과 논의 | M1 (팀원 현황) | 해당 팀원의 관련 태스크 컨텍스트 제공 |
| "이전에 했던 실험" 언급 | M3 (크로스레퍼런스) | 관련 과거 실험 기록 검색 |
| 액션 아이템 도출 | → Asana 태스크 생성 | disc:action → asana_create_task |
| 미팅 시작 시 | M4 (마감 임박) | 팀 전체 긴급 현황 브리핑 |

### 연동 데이터 흐름

```
research-discussion M1(해석)
    │ 팀원명/프로젝트 감지
    ↓
team-dashboard M1/M3 호출
    │ Asana 태스크 + 참조자료 반환
    ↓
research-discussion에서 컨텍스트로 활용
    │ (해석 시 과거 실험 조건 참조)
    ↓
디스커션 완료 → 액션 아이템 도출
    ↓
asana_create_task(project_id, name, assignee, due_on, notes)
    │ 디스커션 결론 → Asana 태스크 자동 생성
    ↓
disc:action에도 기록
```

---

## M5: 마감기한 추천 및 복잡도 분석

### 트리거
- "마감 언제가 적당해?", "기한 추천"
- "이 실험 얼마나 걸릴까?", "일정 잡아줘"
- 태스크 생성 시 due_on 미지정 → 자동 추천

### 실험 복잡도 분류 기준

```
┌─────────────────────────────────────────────────────────────┐
│ Level 1: 단순 (1-3일)                                        │
│ ├─ HPLC 분석, Bradford assay, SDS-PAGE                      │
│ ├─ 기존 프로토콜 반복 실험                                    │
│ ├─ Stock solution 준비, Buffer 조제                           │
│ └─ 문헌 조사, 데이터 정리                                     │
├─────────────────────────────────────────────────────────────┤
│ Level 2: 보통 (3-7일)                                        │
│ ├─ Protein expression & purification (단일 효소)             │
│ ├─ Enzyme activity assay (조건 1-2개)                        │
│ ├─ PCR + Cloning (기존 벡터/프라이머)                         │
│ ├─ 단일 효소 반응 최적화 (1 변수)                             │
│ └─ Cell stock 제작                                           │
├─────────────────────────────────────────────────────────────┤
│ Level 3: 복잡 (1-2주)                                        │
│ ├─ Multi-enzyme cascade reaction 셋업                        │
│ ├─ Cloning (새 유전자 + 신규 벡터 구축)                       │
│ ├─ 효소 반응 최적화 (다변수, factorial design)                │
│ ├─ Co-expression system 구축                                 │
│ └─ Substrate/Product 정제 및 분석법 개발                      │
├─────────────────────────────────────────────────────────────┤
│ Level 4: 장기 (2-4주+)                                       │
│ ├─ 신규 반응 경로 설계 + 검증                                 │
│ ├─ Enzyme engineering (mutation library screening)            │
│ ├─ Scale-up (mL → L)                                         │
│ ├─ Biomass 전처리 → 반응 → 분석 전체 파이프라인               │
│ └─ 논문 작성 + Figure + SI 완성                               │
└─────────────────────────────────────────────────────────────┘
```

### 추천 로직

```
입력: 태스크 이름/설명 + 팀원 + 프로젝트 컨텍스트
    ↓
Step 1: 복잡도 판별
  ├─ 태스크명에서 키워드 추출 (cloning, expression, purification, reaction, optimization 등)
  ├─ 프로젝트 특성 반영 (Nucleoside > Aldaric > Cloning 순 복잡도 경향)
  └─ Level 1-4 분류
    ↓
Step 2: 팀원 역량/이력 반영
  ├─ 과거 유사 태스크 완료 소요일 분석 (completed_at - created_at)
  ├─ 현재 진행 중 태스크 수 (병렬 작업 부하)
  ├─ 최근 완료 속도 (주당 완료 태스크 수)
  └─ 경험 레벨: 해당 유형 태스크 완료 이력 수
    ↓
Step 3: 일정 조정 요소
  ├─ 선행 태스크 의존성 (예: cloning → expression → reaction)
  ├─ 장비 대기 (HPLC, 인큐베이터 등 공유 장비)
  ├─ 주말/공휴일 제외
  ├─ 동시 진행 태스크 수 (3개 이상이면 +30% 버퍼)
  └─ 실험 실패 확률 반영 (신규 실험 +50% 버퍼, 반복 실험 +20%)
    ↓
Step 4: 추천 출력
```

### 출력 형식

```
═══════════════════════════════════════════════════
        마감기한 추천
═══════════════════════════════════════════════════
태스크: {task_name}
담당: {assignee}

복잡도: Level {N} ({label})
근거:
├─ {keyword 1} → {복잡도 요소}
├─ {keyword 2} → {복잡도 요소}
└─ 유사 과거 실험 평균 소요: {N}일

팀원 부하:
├─ 현재 진행 중: {N}건
├─ 이번 주 마감: {N}건
└─ 주간 완료 속도: {N}건/주

추천 기한: {recommended_due_on}
├─ 기본 소요: {base_days}일
├─ 부하 보정: +{buffer_days}일
├─ 선행 태스크 대기: +{dependency_days}일
└─ 최종: {total_days}일 → {date}

대안:
├─ 낙관적: {optimistic_date} ({base_days}일)
├─ 현실적: {realistic_date} ({total_days}일) ★ 추천
└─ 보수적: {conservative_date} ({total_days * 1.5}일)
═══════════════════════════════════════════════════
```

### 팀원별 실험 속도 기준값 (Asana 이력 기반)

```
이자현: 주 4-5건 완료, Level 3까지 단독 수행 가능
김유담: 주 2-3건 완료, Nucleoside 특화 (Level 3 독립)
김지원: 주 1-2건 완료, D-Gal 및 cloning 중심
장서현: 주 1-2건 완료, Aldaric acids + cloning
김호준: 신규 (기준값 없음 → 보수적 추정)
```

### 자동 마감일 관리

디스커션에서 액션 아이템 생성 시 자동 적용:

```
disc:action 도출
    ↓
복잡도 자동 판별 (M5 Step 1-3)
    ↓
asana_create_task(
  name="{action_item}",
  assignee="{팀원 GID}",
  due_on="{추천 기한}",
  notes="복잡도 Level {N}, 추천 근거: {summary}"
)
    ↓
사용자에게 추천 기한 제시 → 승인/수정
```

---

## 작업 생성 원칙 (Task Creation Principles)

### 기본 규칙

1. **워크스페이스**: 항상 `Ribose team (GID: 1202946172929462)` 사용
   - 다른 워크스페이스(FBE, FBE-FBR 공동연구 등)에 생성하지 않도록 주의
2. **담당자**: 특별히 언급 없으면 `me` (이자현) 기본 할당
3. **마감일**: 명시적 언급 없으면 M5 복잡도 분석으로 자동 추천

### 체크박스 항목 → 서브태스크로 생성

Asana notes에 번호 목록으로 적는 대신, **체크박스가 필요한 항목은 서브태스크(subtask)로 생성**한다.

```python
# ❌ 잘못된 방법: notes에 목록으로 작성
asana_create_task(
  name="작업명",
  notes="1. 항목A\n2. 항목B"
)

# ✅ 올바른 방법: 각 항목을 서브태스크로 생성
parent_task = asana_create_task(name="작업명", workspace=RIBOSE_GID, assignee="me")
asana_create_task(name="항목A", parent=parent_task.gid, assignee="me")
asana_create_task(name="항목B", parent=parent_task.gid, assignee="me")
```

**적용 기준:**
- 사용자가 "체크박스로", "항목별로", "리스트로" 언급 시 → 서브태스크
- 디스커션에서 세부 확인사항이 2개 이상 있는 액션 아이템 → 서브태스크
- 단순 메모성 설명 → notes에 작성

---

## API 호출 최적화 (토큰 절약)

### 원칙
1. **opt_fields 항상 지정**: 필요한 필드만 요청 (name,assignee.name,due_on,completed,completed_at,modified_at,permalink_url) — **permalink_url 필수 포함**
2. **limit 설정**: 프로젝트당 최대 100건, 필요시 pagination
3. **캐시 활용**: 같은 세션 내 동일 프로젝트 재조회 최소화
4. **점진적 조회**: 전체 → 특정 프로젝트 → 특정 태스크 순으로 필요한 만큼만

### 조회 우선순위 (팀원별)
```
김유담 → Enzymatic Nucleoside APIs > Cloning
김지원 → D-Gal to DL-Tag > Cloning > Aldaric
장서현 → Aldaric acids > Cloning
김호준 → Aldaric acids > Cloning
이자현 → 전체 (PI)
```

---

## 주의사항

- Asana Free 플랜: advanced search (search_tasks) 사용 불가 → get_tasks + 클라이언트 필터링 사용
- 팀원별 프로젝트 매핑은 주기적 업데이트 필요
- 민감한 연구 데이터는 태스크 notes에 상세 기록하지 않음
- 디스커션 중 Asana 조회는 핵심 정보만 (너무 많은 API 호출 자제)
