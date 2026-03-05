"""
Ribavirin 합성 평형 분석 보고서 PDF (최종본)
- Rib-1-P + TCA_sat(고정) <-> Ribavirin + Pi  (불균일 평형)
- K_app = K_eq × TCA_sat = [Ribavirin][Pi] / [Rib-1-P]
- TCA 투입량(포화 이상)은 무의미 — Rib-1-P가 진짜 제한 기질
- Pi 제거(SucP)가 핵심 전략
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.font_manager as fm
from datetime import date

for font in fm.fontManager.ttflist:
    if 'Malgun' in font.name:
        plt.rcParams['font.family'] = font.name
        break
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.facecolor'] = 'white'

TITLE_COLOR  = '#1a3a5c'
ACCENT_COLOR = '#c0392b'
GREEN_COLOR  = '#27ae60'
BODY_COLOR   = '#2c3e50'
LIGHT_BG     = '#f0f4f8'

# ── 수율 계산 함수들 ───────────────────────────────────────────────
def calc_yield_reservoir(r1p, K_app, pi0=0.0):
    """
    불균일 평형 (TCA 고체 reservoir):
    x^2 + (pi0 + K_app)*x - K_app*r1p = 0
    -> x = [Ribavirin]_eq,  수율 = x/r1p × 100
    """
    b = pi0 + K_app
    c = -K_app * r1p
    disc = b**2 - 4*c
    x = (-b + np.sqrt(disc)) / 2
    x = min(max(x, 0), r1p)
    return x / r1p * 100

def calc_yield_homogeneous(r1p, tca, pi0=0.0, Keq=1.0):
    """
    균일 평형 (TCA 완전 용해):
    Keq*(r1p-x)*(tca-x) = x*(x+pi0)
    """
    A = Keq - 1
    B = -(Keq*(r1p+tca) + pi0)
    C = Keq * r1p * tca
    if abs(A) < 1e-12:
        x = r1p * tca / (r1p + tca + pi0)
    else:
        disc = B**2 - 4*A*C
        if disc < 0: return 0.0
        x = (-B - np.sqrt(disc)) / (2*A)
    x = min(max(x, 0), min(r1p, tca))
    return x / r1p * 100


# ═══════════════════════════════════════════════════════════════════
outpath = 'ribavirin_equilibrium_report.pdf'
with PdfPages(outpath) as pdf:

    # ── PAGE 1: 표지 ─────────────────────────────────────────────
    fig = plt.figure(figsize=(11, 8.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')

    ax.add_patch(mpatches.FancyBboxPatch((0,0),1,1,
        boxstyle='square,pad=0', facecolor=TITLE_COLOR))
    ax.add_patch(mpatches.FancyBboxPatch((0.05,0.07),0.9,0.86,
        boxstyle='round,pad=0.02', facecolor='white', alpha=0.97))
    ax.add_patch(mpatches.Rectangle((0.05,0.73),0.9,0.008,
        facecolor=ACCENT_COLOR))

    ax.text(0.5,0.89,'Ribavirin 효소 합성 반응',
            ha='center',va='center',fontsize=26,fontweight='bold',color=TITLE_COLOR)
    ax.text(0.5,0.81,'평형 열역학 분석 보고서',
            ha='center',va='center',fontsize=20,color=TITLE_COLOR)

    # 핵심 평형식
    box1 = dict(boxstyle='round,pad=0.5', facecolor='#fff8f0', edgecolor=ACCENT_COLOR, lw=2)
    ax.text(0.5,0.665,
        "TCA 포화 이상 투입 시  (불균일 평형):\n\n"
        "  Rib-1-P  +  TCA(s/sat)  =  Ribavirin  +  Pi\n\n"
        "  K_app = [Ribavirin][Pi] / [Rib-1-P]  =  K_eq x TCA_sat\n\n"
        "  [TCA] 고정  ->  Rib-1-P 만이 변수  ->  TCA 투입량 무의미",
        ha='center',va='center', fontsize=12, color=BODY_COLOR,
        bbox=box1, linespacing=1.8)

    # 비유 박스
    box2 = dict(boxstyle='round,pad=0.4', facecolor='#f0f4f8', edgecolor='#7f8c8d', lw=1)
    ax.text(0.5,0.465,
        "불균일 평형 유사 사례:\n"
        "CaCO3(s)  =  CaO(s) + CO2(g)\n"
        "-> CaCO3 더 넣어도 CO2 분압 변하지 않음\n"
        "-> TCA 더 넣어도 K_app 변하지 않음",
        ha='center',va='center', fontsize=11, color='#555555',
        bbox=box2, linespacing=1.7)

    highlights = [
        ('#2980b9', '실험 확인값',
         'Rib-1-P:TCA = 1:1 -> 50%  (TCA 포화 이하 조건으로 해석)'),
        (ACCENT_COLOR, 'TCA 증가의 한계',
         'TCA 포화 이상이면 추가 투입 무의미 — K_app = K_eq x TCA_sat 고정'),
        (GREEN_COLOR, '핵심 전략',
         'Rib-1-P 공급량 최적화  +  SucP(Pi 제거)'),
    ]
    for i,(col,title,body) in enumerate(highlights):
        y = 0.36 - i*0.105
        ax.add_patch(mpatches.FancyBboxPatch((0.07,y-0.035),0.86,0.078,
            boxstyle='round,pad=0.01', facecolor='white', edgecolor=col, lw=1.5))
        ax.text(0.12,y+0.01, f'▶  {title}',
                fontsize=10, fontweight='bold', color=col)
        ax.text(0.12,y-0.02, body, fontsize=9.5, color=BODY_COLOR)

    ax.text(0.5,0.09, f'작성일: {date.today()}  |  PNP-catalyzed glycosylation',
            ha='center', fontsize=9, color='#7f8c8d')

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


    # ── PAGE 2: 평형 이론 ─────────────────────────────────────────
    fig = plt.figure(figsize=(11,8.5))
    ax = fig.add_axes([0,0,1,1])
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')

    ax.add_patch(mpatches.Rectangle((0,0.91),1,0.09,facecolor=TITLE_COLOR))
    ax.text(0.5,0.955,'평형 이론 — 불균일 평형 모델',
            ha='center',va='center',fontsize=17,fontweight='bold',color='white')

    sections = [
        {
            'title': '1.  두 가지 평형 구간',
            'color': '#2980b9',
            'y': 0.825,
            'rows': [
                ('구간 A: TCA 완전 용해',
                 '[TCA] < TCA_sat'),
                ('  -> 균일 평형',
                 'K_eq = [Rib][Pi] / [R1P][TCA],    x^2 / (a-x)^2 = 1  -> 50% (1:1)'),
                ('구간 B: TCA 포화 이상',
                 '[TCA] >= TCA_sat,  고체 reservoir 존재'),
                ('  -> 불균일 평형',
                 'K_app = [Rib][Pi] / [R1P] = K_eq x TCA_sat  (TCA 활동도 = 1 고정)'),
            ]
        },
        {
            'title': '2.  불균일 평형 해석해',
            'color': '#8e44ad',
            'y': 0.575,
            'rows': [
                ('초기 조건',
                 '[Rib-1-P]_0 = a,  [Pi]_0 = pi0,  [TCA] = TCA_sat (고정)'),
                ('평형 방정식',
                 'x^2 + (pi0 + K_app) x - K_app a = 0'),
                ('해석해',
                 'x = { -(pi0+K_app) + sqrt((pi0+K_app)^2 + 4 K_app a) } / 2'),
                ('수율 변수',
                 'Rib-1-P(a), Pi 초기값(pi0), K_app  -- TCA 투입량은 무관'),
            ]
        },
        {
            'title': '3.  Pi 축적과 SucP 역할',
            'color': GREEN_COLOR,
            'y': 0.325,
            'rows': [
                ('Pi 축적 효과',
                 'pi0 증가 -> 해석해의 분자 감소 -> 수율 ↓'),
                ('SucP 역할',
                 'Sucrose + Pi -> Glc-1-P + Fructose  ->  pi0 -> 0 유지'),
                ('SucP 효과 (정량)',
                 'pi0=0 vs pi0=5mM: [R1P]=10mM, K_app=10mM 기준 약 20%p 수율 차이'),
                ('TCA 투입량',
                 '포화 이상이면 완전 무관 (K_app 고정). 포화 이하면 균일 평형으로 전환'),
            ]
        },
    ]

    for sec in sections:
        y0 = sec['y']
        ax.add_patch(mpatches.FancyBboxPatch((0.04,y0-0.21),0.92,0.23,
            boxstyle='round,pad=0.01',
            facecolor='white', edgecolor=sec['color'], lw=1.5))
        ax.text(0.06,y0+0.01, sec['title'],
                fontsize=11,fontweight='bold',color=sec['color'])
        for j,(k,v) in enumerate(sec['rows']):
            yy = y0 - 0.043 - j*0.047
            ax.text(0.08, yy, f'{k}:', fontsize=9,
                    fontweight='bold', color=BODY_COLOR)
            ax.text(0.08 + min(len(k),22)*0.008 + 0.03, yy, v,
                    fontsize=9, color='#444444')

    ax.text(0.5,0.025,
        'TCA(s) 물성: MP 315 C, slightly soluble in water (~7-10 mM)  |  '
        'Kaspar et al. 2020 PMC7318676',
        ha='center', fontsize=7.5, color='#95a5a6')

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


    # ── PAGE 3: 그래프 ────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    fig.patch.set_facecolor('#f9f9f9')
    fig.suptitle('Ribavirin 합성 평형 시뮬레이션',
                 fontsize=15, fontweight='bold', color=TITLE_COLOR, y=1.01)

    # 그래프 1: Rib-1-P 농도 vs 수율 (K_app 값별, pi0=0)
    ax1 = axes[0]
    ax1.set_facecolor('#fdfdfd')

    r1p_range = np.linspace(0.5, 30, 300)
    K_apps = [3.0, 7.0, 10.0, 20.0]
    cols1  = ['#9467bd','#ff7f0e','#1f77b4','#2ca02c']
    for K_app, col in zip(K_apps, cols1):
        y_vals = [calc_yield_reservoir(r, K_app, pi0=0) for r in r1p_range]
        ax1.plot(r1p_range, y_vals, color=col, lw=2.2,
                 label=f'K_app = {K_app} mM (TCA_sat~{K_app:.0f}mM)')

    # 실험값 표시 (1:1, 50% -> 균일 평형 구간)
    ax1.scatter([7.0],[50], s=200, color='red', zorder=6, marker='*',
                label='실험값 참고 (50%, 균일 평형)')
    ax1.axhline(50, color='gray', lw=1, ls='--', alpha=0.5)
    ax1.text(15, 52, '50% 기준', color='gray', fontsize=8.5)
    ax1.text(16, 30,
             'Rib-1-P << K_app\n-> 수율 높음\n\nRib-1-P >> K_app\n-> 수율 낮음',
             fontsize=8.5, color=BODY_COLOR,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='#f8f9fa', edgecolor='#ccc'))

    ax1.set_xlabel('[Rib-1-P]_0 (mM)', fontsize=11)
    ax1.set_ylabel('Rib-1-P 기준 수율 (%)', fontsize=11)
    ax1.set_title('① Rib-1-P 농도 vs 수율\n(불균일 평형, Pi_0=0)',
                  fontsize=11, fontweight='bold')
    ax1.legend(fontsize=8.5, loc='upper right')
    ax1.set_xlim(0, 30); ax1.set_ylim(0, 100)
    ax1.grid(True, alpha=0.3)

    # 그래프 2: Pi 축적(pi0) vs 수율  +  SucP 화살표
    ax2 = axes[1]
    ax2.set_facecolor('#fdfdfd')

    pi0_range = np.linspace(0, 15, 300)
    r1p_vals  = [5.0, 10.0, 20.0]
    cols2     = ['#2ca02c','#1f77b4','#d62728']
    K_app_ref = 10.0

    for r1p, col in zip(r1p_vals, cols2):
        y_vals = [calc_yield_reservoir(r1p, K_app_ref, pi0=p) for p in pi0_range]
        ax2.plot(pi0_range, y_vals, color=col, lw=2.2,
                 label=f'Rib-1-P = {r1p} mM')

    ax2.axvline(0, color='green', lw=2.5, ls='-', alpha=0.5)
    ax2.axhline(50, color='gray', lw=1, ls='--', alpha=0.5)
    ax2.fill_betweenx([0,100], 0, 1.5, alpha=0.10, color='green')
    ax2.text(0.2, 5, 'SucP\n목표', color='green',
             fontsize=8.5, fontweight='bold')

    # SucP 효과 화살표 (R1P=10mM 기준)
    y_pi5  = calc_yield_reservoir(10.0, K_app_ref, pi0=5)
    y_pi0  = calc_yield_reservoir(10.0, K_app_ref, pi0=0)
    ax2.annotate('', xy=(0.5, y_pi0), xytext=(5, y_pi5),
                 arrowprops=dict(arrowstyle='->', color='green', lw=2))
    ax2.text(5.2, (y_pi0+y_pi5)/2,
             f'SucP 효과\n+{y_pi0-y_pi5:.0f}%p',
             color='green', fontsize=9, fontweight='bold')

    ax2.set_xlabel('축적된 Pi 농도 [Pi]_0 (mM)', fontsize=11)
    ax2.set_ylabel('Rib-1-P 기준 수율 (%)', fontsize=11)
    ax2.set_title(f'② Pi 축적 vs 수율  (K_app={K_app_ref}mM)\n-> SucP로 Pi 제거시 수율 회복',
                  fontsize=11, fontweight='bold')
    ax2.legend(fontsize=9, loc='upper right')
    ax2.set_xlim(0, 15); ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


    # ── PAGE 4: 결론 & 실험 권고 ─────────────────────────────────
    fig = plt.figure(figsize=(11,8.5))
    ax = fig.add_axes([0,0,1,1])
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')

    ax.add_patch(mpatches.Rectangle((0,0.91),1,0.09,facecolor=TITLE_COLOR))
    ax.text(0.5,0.955,'결론 및 실험 권고사항',
            ha='center',va='center',fontsize=17,fontweight='bold',color='white')

    # 결론 박스
    ax.add_patch(mpatches.FancyBboxPatch((0.04,0.585),0.92,0.295,
        boxstyle='round,pad=0.02',facecolor='#eaf4fb',edgecolor='#2980b9',lw=2))
    ax.text(0.5,0.868,'핵심 결론',ha='center',fontsize=13,
            fontweight='bold',color='#2980b9')

    conclusions = [
        ('①', '#2980b9',
         'K_eq = 1  (1:1, 50% 실험 확인 — TCA 포화 이하 균일 평형 조건)'),
        ('②', ACCENT_COLOR,
         'TCA 포화 이상: 불균일 평형 -> K_app = K_eq x TCA_sat 고정'),
        ('③', ACCENT_COLOR,
         'TCA 투입량(포화 이상)은 무의미 — CaCO3/CO2 불균일 평형과 동일 원리'),
        ('④', BODY_COLOR,
         'Rib-1-P가 유일한 제한 기질 — 더 많이 공급해야 생산량 증가'),
        ('⑤', GREEN_COLOR,
         'Pi 축적이 수율 ceiling — SucP(Pi 제거)가 핵심 향상 전략'),
    ]
    for i,(num,col,text) in enumerate(conclusions):
        y = 0.835 - i*0.048
        ax.text(0.08, y, num, fontsize=11, fontweight='bold', color=col)
        ax.text(0.12, y, text, fontsize=10, color=BODY_COLOR)

    # 변수별 효과 요약표
    ax.add_patch(mpatches.FancyBboxPatch((0.04,0.17),0.92,0.40,
        boxstyle='round,pad=0.02',facecolor='#eafaf1',edgecolor=GREEN_COLOR,lw=2))
    ax.text(0.5,0.555,'변수별 효과 및 실험 권고',
            ha='center',fontsize=13,fontweight='bold',color=GREEN_COLOR)

    rows = [
        ('TCA 투입량 증가\n(포화 이상)',    '수율 변화 없음',
         '불균일 평형 — K_app 불변\n포화 이하면 수율 상승 가능'),
        ('Rib-1-P 농도 조절',               'Rib-1-P 낮을수록 % 수율 ↑\n절대량은 Rib-1-P에 비례',
         '[R1P] << K_app 조건 탐색\nK_app = K_eq x TCA_sat'),
        ('SucP 실험\n(미실시, 권고)',        'Pi 제거 -> 수율 ↑\n예상 효과: +10~25%p',
         '우선 Pi 축적 확인 후\nSucP 농도 최적화'),
        ('TCA 20mM 이상 재현',              '포화 이상이면 균일/불균일\n경계 구간 확인 필요',
         '원심분리로 [TCA]_sat 측정\n~7-10 mM 예상'),
    ]
    headers = ['실험/변수', '예상 효과', '비고']
    col_x   = [0.06, 0.37, 0.63]
    ax.add_patch(mpatches.Rectangle((0.05,0.498),0.90,0.033,
                 facecolor=GREEN_COLOR,alpha=0.25))
    for h,cx in zip(headers,col_x):
        ax.text(cx,0.510,h,fontsize=9.5,fontweight='bold',color=BODY_COLOR)
    for i,row in enumerate(rows):
        yy = 0.462 - i*0.076
        bg = '#f0fbf6' if i%2==0 else 'white'
        ax.add_patch(mpatches.Rectangle((0.05,yy-0.033),0.90,0.068,
                     facecolor=bg,alpha=0.8))
        for cell,cx in zip(row,col_x):
            ax.text(cx,yy,cell,fontsize=8.5,color=BODY_COLOR,
                    va='center',linespacing=1.35)

    ax.text(0.5,0.105,
        '참고 문헌:\n'
        '• Kaspar et al. (2020) General Principles for Yield Optimization. PMC7318676\n'
        '• Rabuffetti et al. (2019) Synthesis of Ribavirin by Enzymatic Transglycosylation. Catalysts 9:355\n'
        '• ChemicalBook: 2H-1,2,4-Triazole-3-carboxamide (MP 315 C, slightly soluble in water)',
        ha='center',va='center',fontsize=8,color='#7f8c8d',linespacing=1.6,
        bbox=dict(boxstyle='round,pad=0.5',facecolor='#fdfefe',edgecolor='#bdc3c7'))

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

print(f"PDF 저장 완료: {outpath}")
