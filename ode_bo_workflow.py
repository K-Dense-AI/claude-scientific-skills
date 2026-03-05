"""ODE Kinetic Modeling + Bayesian Optimization Workflow Diagram"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(12, 20))
ax.set_xlim(0, 10)
ax.set_ylim(0, 22)
ax.axis('off')
fig.patch.set_facecolor('white')

# ─── 색상 정의 ──────────────────────────────────────────────
C = {
    'phase1': '#DDEEFF',   # 연파랑 - 데이터 수집
    'phase2': '#C8DCFF',   # 파랑 - ODE 모델
    'phase3': '#FFE4C0',   # 오렌지 - 파라미터 추정
    'phase4': '#FFF3C0',   # 노랑 - 검증
    'phase5': '#D4EFC0',   # 초록 - BO 루프
    'phase6': '#B8E0B8',   # 진초록 - 최적 조건
    'decision': '#F5E6FF', # 보라 - 결정 다이아몬드
    'border1': '#2E75B6',
    'border2': '#2E75B6',
    'border3': '#C55A11',
    'border4': '#C99B00',
    'border5': '#375623',
    'border6': '#1F4E1F',
    'arrow': '#404040',
    'feedback': '#888888',
}

def rounded_box(ax, x, y, w, h, color, border_color, lw=1.8, radius=0.25, alpha=1.0):
    box = FancyBboxPatch(
        (x - w/2, y - h/2), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={radius}",
        facecolor=color, edgecolor=border_color,
        linewidth=lw, alpha=alpha, zorder=3
    )
    ax.add_patch(box)
    return box

def diamond(ax, x, y, w, h, color, border_color, lw=1.8):
    pts = np.array([[x, y+h/2], [x+w/2, y], [x, y-h/2], [x-w/2, y]])
    patch = plt.Polygon(pts, closed=True, facecolor=color,
                         edgecolor=border_color, linewidth=lw, zorder=3)
    ax.add_patch(patch)

def arrow_down(ax, x, y_start, y_end, color='#404040'):
    ax.annotate('', xy=(x, y_end), xytext=(x, y_start),
                arrowprops=dict(arrowstyle='->', color=color,
                                lw=2.0, mutation_scale=18))

def text(ax, x, y, s, fontsize=9, bold=False, color='#1a1a1a', ha='center', va='center',
         wrap_width=None):
    weight = 'bold' if bold else 'normal'
    ax.text(x, y, s, ha=ha, va=va, fontsize=fontsize,
            fontweight=weight, color=color, zorder=5,
            linespacing=1.4)

def phase_label(ax, x, y, label, color):
    ax.text(x, y, label, ha='left', va='center', fontsize=8.5,
            fontweight='bold', color=color, zorder=5,
            style='italic')

# ════════════════════════════════════════════════════════════════
# TITLE
# ════════════════════════════════════════════════════════════════
ax.text(5, 21.5, 'ODE Kinetic Modeling → Bayesian Optimization Workflow',
        ha='center', va='center', fontsize=13, fontweight='bold', color='#1a1a1a')
ax.text(5, 21.0, 'for Biocatalytic Process Optimization',
        ha='center', va='center', fontsize=10, color='#555555')

# ════════════════════════════════════════════════════════════════
# PHASE 1 — Experimental Data Collection
# ════════════════════════════════════════════════════════════════
Y1 = 19.8
phase_label(ax, 0.3, Y1+0.7, 'Phase 1', C['border1'])
rounded_box(ax, 5, Y1, 8.5, 1.6, C['phase1'], C['border1'], lw=2.2)
text(ax, 5, Y1+0.45, 'Experimental Data Collection', fontsize=10.5, bold=True, color=C['border1'])
text(ax, 5, Y1-0.05, 'Initial experiments via DoE or Latin Hypercube Sampling (LHS)', fontsize=8.5)
text(ax, 5, Y1-0.45, 'n = 10–20 initial points  ·  concentration time-courses  ·  reaction conditions', fontsize=8)

arrow_down(ax, 5, Y1-0.8, Y1-1.35)

# ════════════════════════════════════════════════════════════════
# PHASE 2 — ODE Model Construction
# ════════════════════════════════════════════════════════════════
Y2 = 17.5
phase_label(ax, 0.3, Y2+0.8, 'Phase 2', C['border2'])
rounded_box(ax, 5, Y2, 8.5, 1.8, C['phase2'], C['border2'], lw=2.2)
text(ax, 5, Y2+0.55, 'ODE Model Construction', fontsize=10.5, bold=True, color=C['border2'])
text(ax, 5, Y2+0.05, 'Mechanistic rate equations  ·  Ordered Bi-Bi / cascade kinetics', fontsize=8.5)
text(ax, 5, Y2-0.32, 'State variables: substrates, products, cofactors (NADPH/NAD⁺), enzymes', fontsize=8)
text(ax, 5, Y2-0.65, 'Parameters: kcat, Km, Ki  ·  Cofactor decomposition sub-models', fontsize=8)

arrow_down(ax, 5, Y2-0.9, Y2-1.45)

# ════════════════════════════════════════════════════════════════
# PHASE 3 — Parameter Estimation
# ════════════════════════════════════════════════════════════════
Y3 = 15.0
phase_label(ax, 0.3, Y3+0.9, 'Phase 3', C['border3'])
rounded_box(ax, 5, Y3, 8.5, 2.0, C['phase3'], C['border3'], lw=2.2)
text(ax, 5, Y3+0.70, 'Parameter Estimation', fontsize=10.5, bold=True, color=C['border3'])
text(ax, 5, Y3+0.22, 'Numerical ODE integration:  scipy.integrate.odeint  /  Runge-Kutta', fontsize=8.5)
text(ax, 5, Y3-0.18, 'Objective:  min Σ(y_sim − y_obs)²   ·   log-transformed parameters', fontsize=8.5)
text(ax, 5, Y3-0.55, 'Optimizer:  lmfit Nelder-Mead  or  scipy trust-constr', fontsize=8.5)

# 결과 박스 (작게)
rounded_box(ax, 5, Y3-0.90, 5.5, 0.42, '#FFF8F0', C['border3'], lw=1.5, radius=0.15)
text(ax, 5, Y3-0.90, '→  Calibrated kinetic parameters  {kcat, Km, Ki, kd} ± σ',
     fontsize=8, color=C['border3'])

arrow_down(ax, 5, Y3-1.1, Y3-1.65)

# ════════════════════════════════════════════════════════════════
# PHASE 4 — Model Validation  (diamond)
# ════════════════════════════════════════════════════════════════
Y4 = 12.6
phase_label(ax, 0.3, Y4+1.0, 'Phase 4', C['border4'])
text(ax, 5, Y4+0.78, 'Model Validation', fontsize=10.5, bold=True, color=C['border4'])
text(ax, 5, Y4+0.38, 'cross-validation  ·  sensitivity analysis  ·  confidence intervals', fontsize=8)

diamond(ax, 5, Y4-0.1, 3.8, 1.05, C['decision'], '#7B4CAF', lw=2.0)
text(ax, 5, Y4-0.1, 'Model\nAcceptable?', fontsize=9, bold=True, color='#4A235A')

# YES 아래로
ax.annotate('', xy=(5, Y4-0.65), xytext=(5, Y4-1.15),
            arrowprops=dict(arrowstyle='->', color='#375623', lw=2.0, mutation_scale=18))
text(ax, 5.35, Y4-0.88, 'Yes', fontsize=8.5, bold=True, color='#375623')

# NO 오른쪽 → 위로 → Phase 2로 피드백
ax.annotate('', xy=(8.7, Y4-0.1), xytext=(6.9, Y4-0.1),
            arrowprops=dict(arrowstyle='->', color='#C55A11', lw=1.6,
                            mutation_scale=15, connectionstyle='arc3,rad=0'))
text(ax, 7.8, Y4+0.18, 'No', fontsize=8.5, bold=True, color='#C55A11')
# 피드백 라인 (오른쪽으로 올라가서 Phase 2 근처로)
ax.plot([8.7, 8.7], [Y4-0.1, Y2], color='#C55A11', lw=1.6, ls='--', zorder=2)
ax.annotate('', xy=(9.3/2+4.25*2/2, Y2), xytext=(8.7, Y2),
            arrowprops=dict(arrowstyle='->', color='#C55A11', lw=1.6, mutation_scale=15))
text(ax, 9.05, (Y4+Y2)/2, 'Revise\nmodel', fontsize=7.5, color='#C55A11', ha='center')

arrow_down(ax, 5, Y4-1.15, Y4-1.7)

# ════════════════════════════════════════════════════════════════
# PHASE 5 — Bayesian Optimization Loop
# ════════════════════════════════════════════════════════════════
Y5_top = 10.7
phase_label(ax, 0.3, Y5_top+0.25, 'Phase 5', C['border5'])

# 큰 BO 루프 박스
rounded_box(ax, 5, Y5_top - 2.3, 8.5, 5.2, C['phase5'], C['border5'], lw=2.5, radius=0.35)
text(ax, 5, Y5_top-0.05, 'Bayesian Optimization Loop', fontsize=10.5, bold=True, color=C['border5'])

# 모델 수식 박스
rounded_box(ax, 5, Y5_top-0.72, 7.0, 0.52, '#F0FFF0', C['border5'], lw=1.5, radius=0.12)
text(ax, 5, Y5_top-0.72, 'y_obs(x) = f_ODE(x) + ε_GP(x)     '
     '[ODE prior + GP residual]',
     fontsize=8.5, color='#1F4E1F')

# BO 사이클 5개 스텝 (사각형 cycle)
cx, cy = 5.0, Y5_top - 2.45
steps = [
    ('①  ODE Prior\nPrediction',   cx,      cy+1.35, C['phase2'], C['border2']),
    ('②  GP Surrogate\nFit',       cx+2.5,  cy+0.45, C['phase5'], C['border5']),
    ('③  EI Acquisition\nFunction', cx+2.5, cy-0.65, C['phase4'], C['border4']),
    ('④  Suggest Next\nConditions', cx,      cy-1.35, C['phase3'], C['border3']),
    ('⑤  Run\nExperiment',         cx-2.5,  cy-0.45, '#E8E8FF', '#4444AA'),
    ('⑥  Update\nGP model',        cx-2.5,  cy+0.65, C['phase5'], C['border5']),
]

bw, bh = 2.05, 0.72
for label, bx, by, fc, ec in steps:
    rounded_box(ax, bx, by, bw, bh, fc, ec, lw=1.6, radius=0.18)
    text(ax, bx, by, label, fontsize=8, bold=False, color='#1a1a1a')

# 사이클 화살표
cycle_pts = [(cx, cy+1.35+bh/2),       # ①아래
             (cx+2.5, cy+0.45+bh/2),   # ②위
             (cx+2.5, cy-0.65+bh/2),   # ③위
             (cx, cy-1.35+bh/2),       # ④위
             (cx-2.5, cy-0.45+bh/2),   # ⑤위
             (cx-2.5, cy+0.65+bh/2),   # ⑥위
             (cx, cy+1.35-bh/2)]       # ①아래로 다시

arrow_color = '#375623'
for i in range(len(cycle_pts)-1):
    x0, y0 = cycle_pts[i]
    x1, y1 = cycle_pts[i+1]
    # 중간 경로
    if i == 0:   # ①→② (오른쪽 위)
        ax.annotate('', xy=(x1, cy+0.45+bh/2), xytext=(cx+bw/2, cy+1.35),
                    arrowprops=dict(arrowstyle='->', color=arrow_color, lw=1.5,
                                    mutation_scale=14, connectionstyle='arc3,rad=-0.2'))
    elif i == 1: # ②→③
        ax.annotate('', xy=(cx+2.5, cy-0.65+bh/2), xytext=(cx+2.5, cy+0.45-bh/2),
                    arrowprops=dict(arrowstyle='->', color=arrow_color, lw=1.5, mutation_scale=14))
    elif i == 2: # ③→④
        ax.annotate('', xy=(cx+bw/2, cy-1.35), xytext=(cx+2.5, cy-0.65-bh/2),
                    arrowprops=dict(arrowstyle='->', color=arrow_color, lw=1.5,
                                    mutation_scale=14, connectionstyle='arc3,rad=-0.2'))
    elif i == 3: # ④→⑤
        ax.annotate('', xy=(cx-2.5+bw/2, cy-0.45), xytext=(cx-bw/2, cy-1.35),
                    arrowprops=dict(arrowstyle='->', color=arrow_color, lw=1.5,
                                    mutation_scale=14, connectionstyle='arc3,rad=-0.2'))
    elif i == 4: # ⑤→⑥
        ax.annotate('', xy=(cx-2.5, cy+0.65-bh/2), xytext=(cx-2.5, cy-0.45+bh/2),
                    arrowprops=dict(arrowstyle='->', color=arrow_color, lw=1.5, mutation_scale=14))
    elif i == 5: # ⑥→①
        ax.annotate('', xy=(cx-bw/2, cy+1.35), xytext=(cx-2.5+bw/2, cy+0.65),
                    arrowprops=dict(arrowstyle='->', color=arrow_color, lw=1.5,
                                    mutation_scale=14, connectionstyle='arc3,rad=-0.2'))

# 수렴 다이아몬드
Yd = Y5_top - 4.45
diamond(ax, 5, Yd, 3.5, 0.95, C['decision'], '#7B4CAF', lw=1.8)
text(ax, 5, Yd, 'Converged?\n(~25–35 exp)', fontsize=8.5, bold=True, color='#4A235A')

# BO 루프 내부 → 다이아몬드
arrow_down(ax, 5, Y5_top-3.82, Yd+0.48)

# No 루프백
ax.plot([6.75, 8.2, 8.2], [Yd, Yd, cy+1.35],
        color='#7B4CAF', lw=1.6, ls='--', zorder=2)
ax.annotate('', xy=(cx+bw/2, cy+1.35), xytext=(8.2, cy+1.35),
            arrowprops=dict(arrowstyle='->', color='#7B4CAF', lw=1.6, mutation_scale=14))
text(ax, 8.5, (Yd+cy+1.35)/2, 'No\n(next\niter.)', fontsize=7.5, color='#7B4CAF', ha='center')

# Yes 아래로
arrow_down(ax, 5, Yd-0.475, Yd-1.05)
text(ax, 5.4, Yd-0.75, 'Yes', fontsize=8.5, bold=True, color=C['border5'])

# ════════════════════════════════════════════════════════════════
# PHASE 6 — Optimal Conditions
# ════════════════════════════════════════════════════════════════
Y6 = Y5_top - 5.95
phase_label(ax, 0.3, Y6+0.72, 'Phase 6', C['border6'])
rounded_box(ax, 5, Y6, 8.5, 1.6, C['phase6'], C['border6'], lw=2.5)
text(ax, 5, Y6+0.50, 'Optimal Process Conditions', fontsize=10.5, bold=True, color=C['border6'])
text(ax, 5, Y6+0.05, 'Maximized yield / titer / enantiomeric excess (ee%)', fontsize=8.5)
text(ax, 5, Y6-0.32, 'Pareto front (multi-objective)  ·  validation experiment', fontsize=8.5)
text(ax, 5, Y6-0.62, '→ Report: optimal [E], [S], pH, T, reaction time', fontsize=8, color=C['border6'])

# ════════════════════════════════════════════════════════════════
# LEGEND
# ════════════════════════════════════════════════════════════════
lx, ly = 0.35, 3.6
ax.text(lx, ly, 'Phase Legend', fontsize=9, fontweight='bold', color='#333333')
legend_items = [
    ('Phase 1  Data Collection',    C['phase1'],   C['border1']),
    ('Phase 2  ODE Model',          C['phase2'],   C['border2']),
    ('Phase 3  Parameter Est.',     C['phase3'],   C['border3']),
    ('Phase 4  Validation',         C['phase4'],   C['border4']),
    ('Phase 5  BO Loop',            C['phase5'],   C['border5']),
    ('Phase 6  Optimal Conditions', C['phase6'],   C['border6']),
]
for i, (label, fc, ec) in enumerate(legend_items):
    bx = lx
    by = ly - 0.52 - i * 0.48
    rect = FancyBboxPatch((bx, by-0.16), 0.38, 0.32,
                          boxstyle="round,pad=0,rounding_size=0.06",
                          facecolor=fc, edgecolor=ec, linewidth=1.5, zorder=3)
    ax.add_patch(rect)
    ax.text(bx+0.52, by, label, fontsize=8, va='center', color='#1a1a1a')

# 참고 논문
ref_y = 0.85
ax.text(0.3, ref_y, 'References:', fontsize=7.5, fontweight='bold', color='#555555')
refs = [
    'petBOA (Comput. Phys. Commun. 2024)  ·  Biotechnol. Bioeng. 2025 (BO in bioprocess review)',
    'ACS Synth. Biol. 2023 (ODE-based circuit optimization)  ·  Anal. Chem. 2022 (Bayesian enzyme kinetics)',
]
for i, r in enumerate(refs):
    ax.text(0.3, ref_y - 0.32*(i+1), r, fontsize=7, color='#777777')

# ════════════════════════════════════════════════════════════════
plt.tight_layout(pad=0.3)
out = 'c:/Users/Jahyun/claude-scientific-skills/ode_bo_workflow.png'
plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
print(f"Saved: {out}")
plt.close()
