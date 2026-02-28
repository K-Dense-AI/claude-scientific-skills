Base directory for this skill: C:\Users\gogok\.claude\skills\primer-design

# Primer Design Suite

iPCR mutagenesis(치환·결실), RE 클로닝, Colony PCR, E. coli 발현 분석, Macrogen 주문서 생성 통합 스킬.

**소스코드**: `./src/primer_design/`
**상세 문서**: `./REFERENCE.md` ← 파라미터·반환값·예시 등 상세 내용은 이 파일 참조

## MCP 서버

```json
{
  "mcpServers": {
    "primer-design": {
      "command": "python",
      "args": ["-m", "src.primer_design.mcp_server"],
      "cwd": "~/.claude/skills/primer-design"
    }
  }
}
```

## 기능 & Tool 매핑

| 기능 | MCP Tool / 클래스 | 용도 |
|------|-------------------|------|
| 아미노산 치환 (SDM) | `iPCRSubstDesigner` | back-to-back iPCR mutagenesis |
| 코돈 결실 | `iPCRDelDesigner` | in-frame/frameshift deletion |
| RE 클로닝 프라이머 | `design_re_cloning_primers` | 발현벡터 삽입용 프라이머 + frame/발현 검증 |
| RE 쌍 추천 | `recommend_re_pair` | 최적 RE 조합 자동 탐색 |
| Colony PCR | `suggest_colony_pcr` | 유니버설 프라이머 + band size 예측 |
| 발현 분석 | `analyze_expression` | CAI, 희귀 코돈, 균주 추천 |
| Reading frame 검증 | `check_reading_frame_tool` | 클로닝 전 frame 호환성 |
| 벡터/RE 목록 | `list_vectors`, `list_restriction_enzymes` | 지원 벡터·효소 조회 |
| Macrogen 주문서 | `generate_macrogen_order` | XLSX 자동 생성 |

## 방법론 선택

| 상황 | 권장 |
|------|------|
| 아미노산 치환 (1–수십 nt) | `iPCRSubstDesigner` |
| 코돈 결실 | `iPCRDelDesigner` |
| 유전자 → 발현벡터 삽입 | `design_re_cloning_primers` |
| RE 조합 모를 때 | `recommend_re_pair` 먼저 |
| 클로닝 후 확인 | `suggest_colony_pcr` |
| 삽입 > 100 nt | Gibson Assembly 고려 |

## 검증 체크리스트 (필수)

- Overlap 상보성: `overlap_verified: True`
- Reading frame 보존 확인
- RE cloning: `expression_check.verdict != "FAIL"`
- QC: PASS 또는 WARNING 허용, FAIL이면 재설계

## 상세 사용법

**파라미터, 반환값, 코드 예시, 설계 원리, 문제 해결** 등은 `REFERENCE.md` 참조:
```
Read ~/.claude/skills/primer-design/REFERENCE.md
```
