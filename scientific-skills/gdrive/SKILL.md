---
name: gdrive
description: Search and access Google Drive files using the mcp__gdrive__search MCP tool. Use when you need to find documents, spreadsheets, presentations, or other files stored in Google Drive.
---

# Google Drive Search

Search and retrieve files from Google Drive via the connected MCP gdrive server.

## Core Tool

**`mcp__gdrive__search`** — Google Drive 파일 검색 도구 (MCP gdrive 서버 필요)

## Usage

### Basic File Search

```
mcp__gdrive__search(query="<검색어>")
```

반환값: 파일명, 파일 ID, MIME 타입, 수정일, 웹 링크 등

### Common Search Examples

```
# 논문 PDF 찾기
mcp__gdrive__search(query="논문 site:drive.google.com type:pdf")

# 스프레드시트 검색
mcp__gdrive__search(query="실험 데이터 mimeType:application/vnd.google-apps.spreadsheet")

# 특정 폴더 내 문서
mcp__gdrive__search(query="회의록 in:folder")

# 키워드로 문서 찾기
mcp__gdrive__search(query="UDH clustering results")
```

---

## 주요 활용 시나리오

| 시나리오 | 쿼리 예시 |
|----------|-----------|
| 논문/보고서 PDF | `"실험명 report"` |
| Google 스프레드시트 | `"데이터 분석"` |
| Google Docs 문서 | `"회의록 2026"` |
| 프레젠테이션 | `"발표자료 슬라이드"` |
| 특정 프로젝트 파일 | `"프로젝트명"` |

---

## 주의사항

- **MCP gdrive 서버가 연결된 세션에서만 작동** (`mcp__gdrive__search` 도구 활성화 필요)
- 검색 결과는 Google Drive 검색 엔진 기준으로 반환 (파일명, 내용 전체 텍스트 포함)
- 파일 다운로드/수정은 이 스킬 범위 외 — 별도 Google Drive API 연동 필요
- 검색어는 한국어/영어 모두 지원

---

## 연동 스킬

- `pdf` — Drive에서 찾은 PDF 파일 분석 시 활용
- `literature-review` — 논문 파일 검색 후 리뷰 작성
- `scientific-writing` — 문서 검색 후 작성 작업 연계
