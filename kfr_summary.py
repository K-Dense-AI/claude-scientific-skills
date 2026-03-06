from fpdf import FPDF

FONT_PATH = "C:/Windows/Fonts/malgun.ttf"
FONT_BOLD = "C:/Windows/Fonts/malgunbd.ttf"

class PDF(FPDF):
    def header(self):
        self.set_font("Malgun", "B", 11)
        self.set_fill_color(30, 80, 150)
        self.set_text_color(255, 255, 255)
        self.cell(0, 10, "5-Ketofructose Reductase (KFR) — Gluconobacter cerinus", align="C", fill=True, new_x="LMARGIN", new_y="NEXT")
        self.set_text_color(0, 0, 0)
        self.ln(3)

    def footer(self):
        self.set_y(-12)
        self.set_font("Malgun", "", 8)
        self.set_text_color(120, 120, 120)
        self.cell(0, 8, f"Page {self.page_no()}", align="C")

    def section_title(self, title):
        self.set_font("Malgun", "B", 11)
        self.set_fill_color(220, 230, 245)
        self.set_text_color(20, 60, 130)
        self.cell(0, 8, f"  {title}", fill=True, new_x="LMARGIN", new_y="NEXT")
        self.set_text_color(0, 0, 0)
        self.ln(2)

    def body_text(self, text, indent=5):
        self.set_font("Malgun", "", 10)
        self.set_x(self.l_margin + indent)
        self.multi_cell(0, 6, text)
        self.ln(1)

    def bullet(self, text, indent=8):
        self.set_font("Malgun", "", 10)
        self.set_x(self.l_margin + indent)
        self.multi_cell(0, 6, f"  - {text}")

    def table_row(self, col1, col2, fill=False):
        self.set_font("Malgun", "", 10)
        if fill:
            self.set_fill_color(240, 245, 255)
        else:
            self.set_fill_color(255, 255, 255)
        self.set_x(self.l_margin + 5)
        self.cell(70, 7, col1, border=1, fill=True)
        self.cell(105, 7, col2, border=1, fill=True, new_x="LMARGIN", new_y="NEXT")

    def table_header_row(self, col1, col2):
        self.set_font("Malgun", "B", 10)
        self.set_fill_color(180, 200, 235)
        self.set_x(self.l_margin + 5)
        self.cell(70, 7, col1, border=1, fill=True)
        self.cell(105, 7, col2, border=1, fill=True, new_x="LMARGIN", new_y="NEXT")


pdf = PDF()
pdf.add_font("Malgun", "", FONT_PATH)
pdf.add_font("Malgun", "B", FONT_BOLD)
pdf.add_page()
pdf.set_auto_page_break(auto=True, margin=15)

# ── 1. 반응 및 역할
pdf.section_title("1. 반응 및 생물학적 역할")
pdf.body_text("반응식:")
pdf.set_font("Malgun", "", 10)
pdf.set_x(pdf.l_margin + 10)
pdf.set_fill_color(248, 248, 248)
pdf.multi_cell(0, 7, "5-keto-D-fructose  +  NADPH  +  H+   <->   D-fructose  +  NADP+", fill=True)
pdf.ln(2)
pdf.bullet("G. cerinus의 D-fructose 대사 경로에서 막결합 fructose dehydrogenase (fructose → 5-KF)의 역방향 환원을 담당")
pdf.bullet("세포 내 NADPH-linked reductase로 처음 정제 (1966, Avigad et al., J Biol Chem)")
pdf.bullet("5-KF는 G. cerinus가 fructose를 산화하면서 생성하는 독특한 dicarbonyl sugar")
pdf.ln(3)

# ── 2. 핵심 생화학적 특성
pdf.section_title("2. 핵심 생화학적 특성")
pdf.ln(1)
pdf.table_header_row("파라미터", "값")
rows = [
    ("Keq (환원 방향)",          "2 x 1010  (fructose 생성 방향으로 극단적 편향)"),
    ("Km  NADPH",                "1.8 x 10-6 M  (1.8 uM) — 매우 높은 친화도"),
    ("Km  5-KF (기질)",           "4.5 x 10-3 M  (4.5 mM)"),
    ("Km  NADP+",                "1.3 x 10-3 M  (1.3 mM)"),
    ("Km  fructose (역반응)",     "7.0 x 10-2 M  (70 mM) — 매우 낮은 친화도"),
    ("Vmax (환원 방향)",          "~920 U/mg protein"),
    ("분자량",                    "71 kDa  (homodimer, ~30 kDa x 2)"),
    ("보조인자 특이성",            "NADPH 전용  (NADH 불사용)"),
    ("최적 pH",                   "7.4"),
]
for i, (c1, c2) in enumerate(rows):
    pdf.table_row(c1, c2, fill=(i % 2 == 0))
pdf.ln(4)

# ── 3. 구조적 분류
pdf.section_title("3. 구조적 분류 및 유전자 정보")
pdf.bullet("Shikimate dehydrogenase (SDH) family의 novel subclass")
pdf.bullet("유전자: GLF_2050  (Gluconobacter sp. strain CHM43 기준)")
pdf.bullet("촉매 잔기: Lys72, Asp108  — hydride transfer dyad (SDH와 보존됨)")
pdf.bullet("Asn21: shikimate 결합을 입체적으로 차단 → 5-KF 기질 특이성 부여")
pdf.bullet("substrate-binding 잔기 9개 중 4개가 canonical SDH와 상이")
pdf.bullet("결정구조: SDH와 유사한 fold이나 기질 결합 포켓 형태가 다름")
pdf.ln(3)

# ── 4. DHA scavenging 전략에서의 의미
pdf.section_title("4. DHA Scavenging 전략 관점에서의 평가")
pdf.ln(1)

pdf.set_font("Malgun", "B", 10)
pdf.set_x(pdf.l_margin + 5)
pdf.cell(0, 7, "  전략 아이디어", new_x="LMARGIN", new_y="NEXT")
pdf.body_text(
    "Fructose를 NAD(P)H 재생 시약으로 활용하여 DHA를 scavenging하는 개념.\n"
    "이론적 흐름: fructose + NADP+ → 5-KF + NADPH  (KFR 역반응 이용)\n"
    "생성된 NADPH로 DHA를 환원 (glycerol 등으로 전환).",
    indent=8)

pdf.set_font("Malgun", "B", 10)
pdf.set_x(pdf.l_margin + 5)
pdf.cell(0, 7, "  [!] 문제점 (Keq 우려 사항)", new_x="LMARGIN", new_y="NEXT")
pdf.bullet("KFR의 Keq = 2 x 1010 이 환원(5-KF → fructose) 방향")
pdf.bullet("역방향 Keq_rev ~= 5 x 10-11 → 산화 방향은 열역학적으로 극히 불리")
pdf.bullet("fructose의 Km (역반응) = 70 mM로 매우 높아 기질 포화도 확보 어려움")
pdf.bullet("=> KFR 단독으로는 NADPH 재생효소로 기능 불가능")
pdf.ln(2)

pdf.set_font("Malgun", "B", 10)
pdf.set_x(pdf.l_margin + 5)
pdf.cell(0, 7, "  결론 및 대안", new_x="LMARGIN", new_y="NEXT")
pdf.body_text(
    "KFR은 NADPH 소비효소 (5-KF 환원). Fructose를 NAD(P)H 재생에 활용하려면\n"
    "역방향을 촉매하는 별도의 fructose 5-dehydrogenase (NADP+-dependent)\n"
    "또는 coupled oxidation system이 별도로 필요함.\n"
    "→ 후보 2 이상에서 다른 효소 조합 검토 필요.",
    indent=8)
pdf.ln(3)

# ── 5. 참고문헌
pdf.section_title("5. 참고문헌")
refs = [
    "Avigad G, Englard S, Pifko S (1966). 5-Keto-D-fructose IV: A specific NADPH-linked reductase\n"
    "     from Gluconobacter cerinus. J Biol Chem 241(2):373-8. PMID: 4379259",
    "Lee S et al. (2021). The 5-Ketofructose Reductase of Gluconobacter sp. Strain CHM43\n"
    "     Is a Novel Class in the Shikimate Dehydrogenase Family. ACS Chem Biol. PMC8425407",
    "Ameyama M et al. (1981). 5-Keto-d-Fructose: Formation and Utilization in G. cerinus.\n"
    "     J Bacteriol. PMC246855",
]
for i, r in enumerate(refs):
    pdf.set_font("Malgun", "", 9)
    pdf.set_x(pdf.l_margin + 8)
    pdf.multi_cell(0, 5, f"[{i+1}] {r}")
    pdf.ln(1)

output_path = "c:/Users/Jahyun/claude-scientific-skills/KFR_Gluconobacter_cerinus_summary.pdf"
pdf.output(output_path)
print(f"Saved: {output_path}")
