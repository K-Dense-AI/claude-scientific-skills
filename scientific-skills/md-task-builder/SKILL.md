# MD Task Builder

기존 스킬을 **Claude 웹용 `.md` 지시 파일**로 변환합니다.
로컬 라이브러리 없이 claude.ai에 첨부하면 Python으로 실행됩니다.

## 사용법

```
/md-task-builder [스킬명]
예: /md-task-builder primer-design
```

## 사용 가능한 템플릿

| 파일 | 용도 | 커버하는 기능 |
|------|------|-------------|
| [templates/primer-design.md](templates/primer-design.md) | 프라이머 설계 | SDM (iPCR back-to-back), RE 클로닝 |

## 사용 방법

1. `templates/` 에서 필요한 `.md` 파일 열기
2. 플레이스홀더(`[...]`) 를 실제 값으로 교체
3. [claude.ai](https://claude.ai) → 📎 → `.md` 첨부 → 전송
4. Claude가 Python 코드 작성 + 실행 + 결과 출력

## 새 템플릿 추가 원칙

- 알고리즘 가이드 포함 (Claude가 라이브러리 없이 구현 가능하게)
- 사용자 입력 섹션과 알고리즘 섹션 명확히 분리
- QC 기준 명시
- 출력 형식 지정
