import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

# === 파라미터 ===
C_sat = 8.0          # TCA 포화농도 (mM)
K_eq = 0.8           # 평형상수 (임의값, 시각화 목적)
S0 = 1.0             # 기질 초기농도 (mM, normalized)

# === TCA 투입량 범위 ===
tca_added = np.linspace(0, 30, 500)  # mM

# === [TCA(aq)] 계산: 불균일 평형 적용 ===
tca_aq = np.where(tca_added < C_sat, tca_added, C_sat)

# === Yield 계산 (단순 1:1 평형 가정) ===
# K_eq = yield / ((1-yield) * [TCA])  →  yield = K_eq*[TCA] / (1 + K_eq*[TCA])
yield_pct = K_eq * tca_aq / (1 + K_eq * tca_aq) * 100

# K_eq가 큰 경우 (강한 반응) 비교용
K_high = 3.0
yield_high = K_high * tca_aq / (1 + K_high * tca_aq) * 100

# === 플롯 ===
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
fig.suptitle("Heterogeneous Equilibrium: TCA Input vs. Yield", fontsize=14, fontweight='bold', y=1.02)

# ── 왼쪽: Yield vs TCA added ──
ax = axes[0]

ax.plot(tca_added, yield_pct, 'steelblue', lw=2.5, label=f'K_eq = {K_eq}')
ax.plot(tca_added, yield_high, 'tomato', lw=2.5, linestyle='--', label=f'K_eq = {K_high}')

# 포화농도 경계선
ax.axvline(C_sat, color='gray', linestyle=':', lw=1.8, alpha=0.8)
ax.axhspan(0, 100, xmin=C_sat/30, alpha=0.04, color='orange')

# 영역 레이블
ax.text(C_sat*0.45, 92, 'Homogeneous\n[TCA] ∝ input', ha='center', fontsize=9,
        color='steelblue', style='italic')
ax.text(C_sat*1.7, 92, 'Heterogeneous\n[TCA] = const', ha='center', fontsize=9,
        color='darkorange', style='italic')

# Yield plateau 표시 (K_eq = 0.8 기준)
y_plateau = K_eq * C_sat / (1 + K_eq * C_sat) * 100
ax.axhline(y_plateau, color='steelblue', linestyle=':', lw=1.2, alpha=0.6)
ax.annotate(f'Plateau\n({y_plateau:.1f}%)', xy=(25, y_plateau), fontsize=8.5,
            color='steelblue', va='center',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='steelblue', alpha=0.7))

y_plateau_h = K_high * C_sat / (1 + K_high * C_sat) * 100
ax.axhline(y_plateau_h, color='tomato', linestyle=':', lw=1.2, alpha=0.6)
ax.annotate(f'Plateau\n({y_plateau_h:.1f}%)', xy=(25, y_plateau_h+2), fontsize=8.5,
            color='tomato', va='center',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='tomato', alpha=0.7))

ax.set_xlabel('TCA Input (mM)', fontsize=11)
ax.set_ylabel('Yield (%)', fontsize=11)
ax.set_xlim(0, 30)
ax.set_ylim(0, 100)
ax.legend(fontsize=10, loc='lower right')
ax.set_title('Yield vs. TCA Input', fontsize=12)
ax.annotate(f'Solubility\nlimit\n~{C_sat} mM', xy=(C_sat, 50), xytext=(C_sat+3, 45),
            arrowprops=dict(arrowstyle='->', color='gray'), fontsize=8.5, color='gray')

# ── 오른쪽: [TCA(aq)] vs TCA added ──
ax2 = axes[1]

ax2.plot(tca_added, tca_aq, 'darkorange', lw=2.5)
ax2.fill_between(tca_added, tca_aq, alpha=0.15, color='darkorange')

# 포화 이전: 선형 구간
ax2.plot([0, C_sat], [0, C_sat], 'gray', lw=1.2, linestyle='--', alpha=0.5, label='Linear (ideal)')

# 포화 농도 경계
ax2.axvline(C_sat, color='gray', linestyle=':', lw=1.8, alpha=0.8)
ax2.axhline(C_sat, color='gray', linestyle=':', lw=1.2, alpha=0.5)

# 고체 TCA 존재 표현
for x_pos in np.linspace(C_sat+1, 29, 6):
    ax2.annotate('', xy=(x_pos, C_sat*0.15), xytext=(x_pos, C_sat*0.55),
                arrowprops=dict(arrowstyle='->', color='saddlebrown', lw=1.2))
ax2.text(22, C_sat*0.65, 'TCA(s) dissolves\nto maintain [TCA(aq)]', fontsize=8.5,
         color='saddlebrown', ha='center',
         bbox=dict(boxstyle='round', fc='wheat', ec='saddlebrown', alpha=0.7))

# 영역 음영
ax2.axvspan(0, C_sat, alpha=0.06, color='steelblue', label='Homogeneous')
ax2.axvspan(C_sat, 30, alpha=0.06, color='darkorange', label='Heterogeneous')

ax2.set_xlabel('TCA Input (mM)', fontsize=11)
ax2.set_ylabel('[TCA(aq)] (mM)', fontsize=11)
ax2.set_xlim(0, 30)
ax2.set_ylim(0, 16)
ax2.legend(fontsize=9, loc='upper left')
ax2.set_title('[TCA(aq)] vs. TCA Input', fontsize=12)
ax2.annotate(f'C_sat ≈ {C_sat} mM', xy=(C_sat, C_sat), xytext=(C_sat+2, C_sat+1.5),
            arrowprops=dict(arrowstyle='->', color='gray'), fontsize=8.5, color='gray')

plt.tight_layout()
plt.savefig('tca_yield_heterogeneous.png', dpi=150, bbox_inches='tight')
print("저장 완료: tca_yield_heterogeneous.png")
plt.show()
