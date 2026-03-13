# Configuration Guide / 설정 가이드

## 파일 구조

```
apps/ai-code-execution/
├── app.py           # Streamlit 앱 (MD 편집기)
├── requirements.txt # streamlit>=1.32.0 만 필요
├── README.ko.md     # 한글 사용 설명서
└── CONFIGURE.md     # 이 파일
```

---

## 로컬 실행

```bash
cd apps/ai-code-execution
pip install streamlit
streamlit run app.py
```

API 키 불필요 — 외부 서비스 없음.

---

## Streamlit Cloud 배포

1. [share.streamlit.io](https://share.streamlit.io) → **New app**
2. 설정:
   - Repository: `jahyunlee00299/claude-scientific-skills`
   - Branch: `feat/team-orchestrator-v2`
   - Main file path: `apps/ai-code-execution/app.py`
3. **Deploy** — Secrets 불필요

---

## 앱이 하는 일

1. NO CODE / WITH CODE 탭에서 `.md` 템플릿 제공
2. 사용자가 내용 편집
3. **📥 Download** 버튼으로 파일 저장
4. 사용자가 claude.ai에 `.md` + 데이터 파일 첨부 후 전송
5. Claude가 코드 작성·실행·결과 제공
