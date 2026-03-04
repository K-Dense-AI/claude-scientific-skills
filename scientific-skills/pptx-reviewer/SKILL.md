name: pptx-reviewer
description: PowerPoint 파일에 학술 피드백 댓글(말풍선) 자동 추가 스킬. PPTX 파일을 읽고 슬라이드별 내용/노트를 분석하여 학술적 엄밀성 기준으로 피드백을 PowerPoint 네이티브 댓글(Comment) 형식으로 삽입. 랩미팅 PPT 리뷰, 논문 발표 초고 검토에 사용.

---

# PPTX Reviewer — PowerPoint 학술 댓글 자동 삽입

## 개요

PPTX 파일을 읽고 각 슬라이드의 내용과 발표자 노트를 분석하여, PowerPoint 네이티브 댓글(Comment 말풍선) 형식으로 학술적 피드백을 삽입합니다.

## 트리거

- "PPT에 댓글 달아줘"
- "PPT 피드백 달아줘" / "슬라이드 리뷰해줘"
- "학술적으로 피드백해줘"
- "PPT 오타 찾아줘" / "오타 수정해줘"

## 실행 절차

### Step 1. 파일 읽기 (markitdown — 내용/노트 분석용)

```python
from markitdown import MarkItDown
md = MarkItDown()
result = md.convert(pptx_path)
with open("C:/Users/Jahyun/tmp_ppt.txt", 'w', encoding='utf-8', errors='replace') as f:
    f.write(result.text_content)
```

### Step 1.5. 정밀 텍스트 추출 및 오타 탐지 (python-pptx)

markitdown은 유니코드 특수문자로 인해 중간에 잘릴 수 있으므로, 오타 탐지는 python-pptx로 별도 추출.

```python
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
from pptx import Presentation

prs = Presentation(pptx_path)
for i, slide in enumerate(prs.slides, 1):
    print(f'\n=== Slide {i} ===')
    for shape in slide.shapes:
        if shape.has_text_frame:
            for para in shape.text_frame.paragraphs:
                text = para.text.strip()
                if text:
                    print(repr(text))
```

**오타 수정 시 주의사항 — run 분리 문제**:
- PPTX에서 텍스트는 run 단위로 분리 저장되므로 `shape.text`에서는 보이는 오타가 단일 run에 없을 수 있음
- 오타가 run 경계에 걸린 경우: 개별 run 출력으로 확인 후 해당 run만 수정

```python
# run 분리 확인
for para in shape.text_frame.paragraphs:
    runs = [r.text for r in para.runs]
    if any('오타단어' in r for r in runs):
        print('runs:', runs)

# 오타가 단일 run에 있는 경우 — replace 사용
for run in para.runs:
    if '오타' in run.text:
        run.text = run.text.replace('오타', '수정')

# 오타가 run 분리된 경우 — 해당 run만 직접 수정
for run in para.runs:
    if run.text == 'dependen':   # 'PQQ '와 분리된 경우
        run.text = 'dependent'
```

오타 수정 후 `prs.save(out_path)`로 저장.

### Step 2. 슬라이드별 분석

각 슬라이드(본문 + 노트 모두) 확인 항목:

**내용 완성도**
- "나중에 보완", "추후 수정" 등의 미완성 표시 → 수정 필수 표기
- 슬라이드 본문이 없고 노트만 있는 경우 → 시각화 필요

**학술적 엄밀성 체크리스트**
- [ ] 통계: n수, error bar, p-value 존재 여부
- [ ] 방법론: 실험 조건 가정 근거 명시 여부
- [ ] 문헌 인용: 가설/결론 근거 문헌 여부
- [ ] 효소 kinetics: Km 추정 시 [S] 범위가 Km 주변을 포함하는지 (0.2×Km ~ 10×Km)
- [ ] 분석 신뢰성: detection limit, S/N ratio 고려 여부
- [ ] 데이터 보정: 보정 방법론 가정 명시 여부
- [ ] 결론 타당성: 데이터→결론 논리 검토

### Step 3. 댓글 삽입 — PowerPoint COM API (권장)

**win32com 방식 (최신 PowerPoint 365에서 정상 동작)**:

```python
import win32com.client, shutil, time

shutil.copy2(src_path, dest_path)

ppt = win32com.client.Dispatch("PowerPoint.Application")
ppt.Visible = True  # False 불가 (PowerPoint COM 제한)

try:
    pres = ppt.Presentations.Open(dest_path, ReadOnly=False, Untitled=False, WithWindow=False)
    time.sleep(1)

    # COMMENTS: [(slide_1based, left_pt, top_pt, text), ...]
    for slide_num, left, top, text in COMMENTS:
        slide = pres.Slides(slide_num)
        slide.Comments.Add(
            Left=left, Top=top,
            Author="학술리뷰",
            AuthorInitials="리",
            Text=text
        )
    pres.Save()
    pres.Close()
finally:
    ppt.Quit()
```

**주의사항**:
- `Width`, `Height`는 `Add()`에서 지원하지 않음 — Left/Top/Author/AuthorInitials/Text만 사용
- `ppt.Visible = False` 불가 — PowerPoint가 잠깐 화면에 뜸
- 저장 경로: `{원본파일명}_리뷰.pptx` (원본 절대 덮어쓰지 않음)
- Python 실행: `C:/Users/Jahyun/anaconda3/python.exe`
- 스크립트는 `C:/Users/Jahyun/` 에 저장 후 실행 (한글 경로 /tmp/ 접근 불가)

### Step 4. 출력

- 저장 경로: `{원본파일명}_리뷰.pptx`
- PowerPoint에서 열고 **검토(Review) > 댓글(Comments)** 패널에서 확인

## 댓글 위치 가이드

```
기본 위치: Left=10, Top=10 (포인트 단위)
슬라이드 우측: Left=600, Top=10
슬라이드 하단: Left=10, Top=450
```

## 학술 피드백 템플릿

```
[학술 피드백] {주제}:
- {문제점 1}: {개선안}
- {문제점 2}: {개선안}
```

```
[학술 피드백 + 수정 필수] {문제}:
- {조치 사항}
```

```
[검토 필요] {사항}:
- 결정 후 → A 또는 B 선택
```
