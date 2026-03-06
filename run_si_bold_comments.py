"""SI 파일 bold 서식 분석 및 Word 코멘트 추가 - 한글 경로 안전 버전"""
import sys, io, re, zipfile, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import docx
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from lxml import etree

COMMENT_AUTHOR = "Claude"
COMMENT_DATE = "2026-03-02T00:00:00Z"
W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
W = W_NS

# 하드코딩된 경로 (string literal → UTF-8로 정확히 읽힘)
SI_PATH = "C:/Users/Jahyun/OneDrive - 호서대학교/논문/LRib/LRib Manuscripts/투고 이전/JACS/260302_L-Rib_SI_JACS.docx"


def make_comment_xml(comment_id, text):
    comment = OxmlElement('w:comment')
    comment.set(qn('w:id'), str(comment_id))
    comment.set(qn('w:author'), COMMENT_AUTHOR)
    comment.set(qn('w:date'), COMMENT_DATE)
    p = OxmlElement('w:p')
    p.append(OxmlElement('w:pPr'))
    r = OxmlElement('w:r')
    r.append(OxmlElement('w:rPr'))
    t = OxmlElement('w:t')
    t.text = text
    t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    r.append(t)
    p.append(r)
    comment.append(p)
    return comment


def add_comment_ref(para, comment_id):
    p = para._p
    crs = OxmlElement('w:commentRangeStart')
    crs.set(qn('w:id'), str(comment_id))
    p.insert(0, crs)
    cre = OxmlElement('w:commentRangeEnd')
    cre.set(qn('w:id'), str(comment_id))
    p.append(cre)
    r = OxmlElement('w:r')
    rPr = OxmlElement('w:rPr')
    rStyle = OxmlElement('w:rStyle')
    rStyle.set(qn('w:val'), 'CommentReference')
    rPr.append(rStyle)
    r.append(rPr)
    ref = OxmlElement('w:commentReference')
    ref.set(qn('w:id'), str(comment_id))
    r.append(ref)
    p.append(r)


def get_bold_and_nonbold(para):
    bold, nonbold = [], []
    for run in para.runs:
        t = run.text
        if not t: continue
        b = run._r.find('{%s}rPr/{%s}b' % (W, W))
        is_b = b is not None and b.get(qn('w:val'), 'true').lower() not in ('false', '0')
        (bold if is_b else nonbold).append(t)
    return ''.join(bold), ''.join(nonbold)


# ── 분석 ──────────────────────────────────────────
doc = docx.Document(SI_PATH)
paras = doc.paragraphs
new_comments = []

SKIP = {'BA_Title','BB_Author_Name','BC_Author_Address',
        'FA_Corresponding_Author_Footnote','FA_Author_Info_Subtitle',
        'TD_Acknowledgments','TF_References_Section','TE_Supporting_Information',
        'VC_Scheme_Title','SchemeCaption','EndNote Bibliography'}

print("=== SI 캡션 bold 분석 ===")
for i, para in enumerate(paras):
    s = para.text.strip()
    if not re.match(r'^(Figure|Table|Scheme)\s*[S\d]', s): continue

    bold_t, non_bold_t = get_bold_and_nonbold(para)
    cap_match = re.match(r'^(Figure|Table|Scheme)\s*[\dS]+\.?\s*', s)
    label = cap_match.group(0).strip() if cap_match else ''

    if not bold_t.strip():
        status = "전체_nonbold"
        new_comments.append((i, f"[서식-캡션] 번호({label})에 bold 없음 → ACS: 번호만 bold 필요"))
    elif bold_t.strip() == label:
        status = "OK"
    elif not non_bold_t.strip():
        status = "전체bold_ERROR"
        new_comments.append((i, f"[서식-캡션] 전체 캡션 bold → ACS: 번호({label})만 bold, 나머지 plain"))
    else:
        status = "혼합bold"
        if bold_t.strip() != label:
            new_comments.append((i, f"[서식-캡션] bold 불일치: BOLD='{bold_t[:50]}' (ACS: 번호만 bold)"))

    print(f"  [{i}] {status}: '{s[:60]}'")

print(f"\n=== SI 본문 부분 bold 분석 ===")
for i, para in enumerate(paras):
    style = para.style.name
    if style in SKIP: continue
    s = para.text.strip()
    if not s or re.match(r'^(Figure|Table|Scheme)\s*[S\d]', s): continue

    bold_t, non_bold_t = get_bold_and_nonbold(para)
    if not bold_t.strip() or not non_bold_t.strip(): continue
    if re.fullmatch(r'[\d\s\-()]+', bold_t.strip()): continue  # ACS 연도/번호 정상
    if bold_t.strip() == s and len(s) < 80: continue  # heading

    print(f"  [{i}] Style={style!r} BOLD='{bold_t[:60]}'")
    new_comments.append((i, f"[서식-본문bold] '{bold_t[:80]}' 의도적 강조인지 확인. 트래킹 변경 잔재 가능성"))

print(f"\n추가할 신규 코멘트: {len(new_comments)}개")
if not new_comments:
    print("SI 서식 정상 - 추가 코멘트 없음")
    sys.exit(0)

# ── 기존 코멘트 파악 ──────────────────────────────
with zipfile.ZipFile(SI_PATH, 'r') as z:
    existing_xml = z.read('word/comments.xml') if 'word/comments.xml' in z.namelist() else None

existing_elems = []
max_id = -1
if existing_xml:
    ex_tree = etree.fromstring(existing_xml)
    for c in ex_tree.findall('{%s}comment' % W):
        cid = int(c.get('{%s}id' % W, -1))
        max_id = max(max_id, cid)
        existing_elems.append(c)

start_id = max_id + 1
print(f"기존 코멘트 {len(existing_elems)}개 (max id={max_id}), 신규 시작 id={start_id}")

# paragraph에 참조 삽입
for offset, (para_idx, _) in enumerate(new_comments):
    add_comment_ref(paras[para_idx], start_id + offset)

# 새 comments.xml 빌드
NSMAP = {'w': W_NS, 'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
new_root = etree.Element('{%s}comments' % W_NS, nsmap=NSMAP)
for el in existing_elems:
    new_root.append(el)
for offset, (_, text) in enumerate(new_comments):
    new_root.append(make_comment_xml(start_id + offset, text))

comments_bytes = etree.tostring(new_root, xml_declaration=True, encoding='UTF-8', standalone=True)

# zip 조작
tmp = SI_PATH + '.boldtmp.docx'
doc.save(tmp)

CT_NS = 'http://schemas.openxmlformats.org/package/2006/content-types'
PKG_NS = 'http://schemas.openxmlformats.org/package/2006/relationships'

with zipfile.ZipFile(tmp, 'r') as zin:
    with zipfile.ZipFile(SI_PATH, 'w', zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            if item.filename == 'word/comments.xml':
                continue
            data = zin.read(item.filename)

            if item.filename == 'word/_rels/document.xml.rels':
                tree = etree.fromstring(data)
                ids = {r.get('Id') for r in tree}
                if 'rIdComments' not in ids:
                    rel = etree.SubElement(tree, '{%s}Relationship' % PKG_NS)
                    rel.set('Id', 'rIdComments')
                    rel.set('Type', 'http://schemas.openxmlformats.org/officeDocument/2006/relationships/comments')
                    rel.set('Target', 'comments.xml')
                data = etree.tostring(tree, xml_declaration=True, encoding='UTF-8', standalone=True)

            if item.filename == '[Content_Types].xml':
                tree = etree.fromstring(data)
                parts = {el.get('PartName') for el in tree.findall('{%s}Override' % CT_NS)}
                if '/word/comments.xml' not in parts:
                    el = etree.SubElement(tree, '{%s}Override' % CT_NS)
                    el.set('PartName', '/word/comments.xml')
                    el.set('ContentType', 'application/vnd.openxmlformats-officedocument.wordprocessingml.comments+xml')
                data = etree.tostring(tree, xml_declaration=True, encoding='UTF-8', standalone=True)

            zout.writestr(item, data)
        zout.writestr('word/comments.xml', comments_bytes)

os.remove(tmp)
total = len(existing_elems) + len(new_comments)
print(f"\n완료! 총 코멘트 {total}개 (기존 {len(existing_elems)}개 + 신규 {len(new_comments)}개)")
