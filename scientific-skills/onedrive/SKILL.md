---
name: onedrive
description: Access and search Microsoft OneDrive files using the MCP onedrive server. Use when you need to find or work with files stored in OneDrive/SharePoint.
---

# OneDrive File Access

Search and access files from Microsoft OneDrive via the connected MCP onedrive server.

## Core Tool

**MCP onedrive 서버 도구** — OneDrive/SharePoint 파일 접근 (MCP onedrive 서버 필요)

## Usage

### File Search

OneDrive MCP 서버가 제공하는 도구로 파일 검색 및 접근:

```
# 파일명 또는 키워드로 검색
onedrive_search(query="<검색어>")

# 특정 폴더 내 파일 목록
onedrive_list(path="<폴더 경로>")
```

### Common Search Examples

```
# 논문 폴더 검색
onedrive_search(query="논문 2026")

# 발표자료 찾기
onedrive_search(query="발표 pptx")

# 실험 데이터 파일
onedrive_search(query="실험 데이터 xlsx")

# 특정 OneDrive 계정 경로
onedrive_list(path="OneDrive - 호서대학교/논문")
```

---

## 주요 활용 시나리오

| 시나리오 | 설명 |
|----------|------|
| 논문/보고서 파일 | OneDrive에 저장된 PDF, DOCX 검색 |
| 대학교 OneDrive | 호서대학교/고려대학교 OneDrive 파일 접근 |
| 발표자료 | PPTX 파일 검색 및 내용 확인 |
| 실험 데이터 | Excel/CSV 파일 검색 |
| SharePoint 문서 | 팀 공유 문서 접근 |

---

## 연결된 OneDrive 계정 (현재 환경)

- `OneDrive - 호서대학교` — `c:\Users\Jahyun\OneDrive - 호서대학교\`
- `OneDrive - 고려대학교` — `c:\Users\Jahyun\OneDrive - 고려대학교\저장소`

> 로컬 동기화 경로에 직접 접근하려면 Bash 도구로 파일 읽기 가능.

---

## 주의사항

- **MCP onedrive 서버가 연결된 세션에서만 원격 OneDrive API 작동**
- 로컬 동기화 파일은 Bash/Read 도구로 직접 접근 가능 (MCP 불필요)
- 파일 수정/업로드는 별도 Microsoft Graph API 연동 필요
- SharePoint 문서 라이브러리 접근 시 권한 확인 필요

---

## 연동 스킬

- `pdf` — OneDrive에서 찾은 PDF 파일 분석
- `docx` — Word 문서 읽기/편집
- `xlsx` — Excel 파일 처리
- `pptx` — PowerPoint 파일 분석
