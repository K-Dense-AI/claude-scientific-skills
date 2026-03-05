import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# 한글 폰트 설정
for font in fm.fontManager.ttflist:
    if 'Malgun' in font.name or 'NanumGothic' in font.name or 'Gulim' in font.name:
        plt.rcParams['font.family'] = font.name
        break
plt.rcParams['axes.unicode_minus'] = False

def calc_yield(K_N, donor_ratio, tca=1.0):
    """
    equilibrium yield of ribavirin (% of TCA converted)
    K_N = K1/K2,  donor_ratio = [donor]0 / [TCA]0
    """
    a = donor_ratio * tca
    b = tca
    if abs(K_N - 1.0) < 1e-9:
        x = a * b / (a + b)
    else:
        disc = K_N**2 * (a + b)**2 - 4 * K_N * (K_N - 1) * a * b
        if disc < 0:
            return 0.0
        x = (K_N * (a + b) - np.sqrt(disc)) / (2 * (K_N - 1))
    x = min(x, min(a, b))
    return x / b * 100

# ── Figure ──
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor('#f9f9f9')

# ── Plot 1: Donor:TCA 비율 vs 수율 (K_N별) ──────────────────────────
ax1 = axes[0]
ax1.set_facecolor('#fdfdfd')

donor_ratios = np.linspace(0.2, 10, 300)
K_N_values   = [0.5, 1.0, 2.0, 4.0, 10.0]
colors       = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4', '#9467bd']
labels       = ['K_N = 0.5', 'K_N = 1.0  ← 현재 시스템', 'K_N = 2.0', 'K_N = 4.0', 'K_N = 10.0']

for K_N, color, label in zip(K_N_values, colors, labels):
    lw = 3.0 if K_N == 1.0 else 1.8
    ls = '-'
    y = [calc_yield(K_N, r) for r in donor_ratios]
    ax1.plot(donor_ratios, y, color=color, linewidth=lw, linestyle=ls, label=label)

# 확인된 실험 포인트
ax1.scatter([1.0], [50.0], s=200, color='red', zorder=6, marker='*', label='실험 확인값 (1:1, 50%)')

# 참고선
ax1.axhline(50, color='gray', lw=1, ls='--', alpha=0.5)
ax1.axvline(1.0, color='gray', lw=1, ls='--', alpha=0.5)

# SucP 화살표 (K_N=1 → 높은 K_N 이동 효과)
ax1.annotate('', xy=(4.0, calc_yield(4.0, 4.0)),
             xytext=(4.0, calc_yield(1.0, 4.0)),
             arrowprops=dict(arrowstyle='->', color='green', lw=2))
ax1.text(4.15, 72, 'SucP로 Pi 제거\n→ 실질 K_N ↑\n→ 수율 상승', color='green',
         fontsize=9, va='center')

ax1.set_xlabel('Donor : TCA 비율', fontsize=12)
ax1.set_ylabel('이론 수율 (%, TCA 기준)', fontsize=12)
ax1.set_title('Donor:TCA 비율 vs 이론 수율\n(K_N값 별)', fontsize=13, fontweight='bold')
ax1.legend(loc='lower right', fontsize=9)
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 100)
ax1.grid(True, alpha=0.3)

# ── Plot 2: TCA 증가가 왜 효과 없는가 (K_N=1, Donor 고정) ──────────────
ax2 = axes[1]
ax2.set_facecolor('#fdfdfd')

tca_ratios = np.linspace(0.2, 10, 300)   # TCA:Donor 배수

def yields_tca_sweep(K_N, tca_ratio, donor=1.0):
    a = donor
    b = tca_ratio * donor
    if abs(K_N - 1.0) < 1e-9:
        x = a * b / (a + b)
    else:
        disc = K_N**2*(a+b)**2 - 4*K_N*(K_N-1)*a*b
        x = (K_N*(a+b) - np.sqrt(max(disc, 0))) / (2*(K_N-1))
    x = min(x, min(a, b))
    return x / b * 100, x / a * 100   # TCA 기준 수율, Donor 기준 수율

y_tca_kn1, y_don_kn1 = zip(*[yields_tca_sweep(1.0, t) for t in tca_ratios])
y_tca_kn4, y_don_kn4 = zip(*[yields_tca_sweep(4.0, t) for t in tca_ratios])

ax2.plot(tca_ratios, y_tca_kn1, color='#ff7f0e', lw=2.5,
         label='TCA 기준 수율  (K_N=1)')
ax2.plot(tca_ratios, y_don_kn1, color='#ff7f0e', lw=2.5, ls='--',
         label='Donor 기준 수율 (K_N=1)')
ax2.plot(tca_ratios, y_tca_kn4, color='#1f77b4', lw=2.0, alpha=0.7,
         label='TCA 기준 수율  (K_N=4, 참고)')

# 1:1 포인트
ax2.scatter([1.0], [50], s=200, color='red', zorder=6, marker='*')
ax2.axhline(50, color='gray', lw=1, ls='--', alpha=0.5)
ax2.axvline(1.0, color='gray', lw=1, ls='--', alpha=0.5)

# 해석 텍스트
ax2.annotate('TCA ↑ →\nTCA 기준 수율 ↓\n(Donor가 limiting)',
             xy=(4, list(y_tca_kn1)[int(4/10.2*300)]),
             xytext=(5.5, 30),
             arrowprops=dict(arrowstyle='->', color='#d62728'),
             color='#d62728', fontsize=9)

ax2.annotate('Donor 기준 수율은\n오히려 ↑\n(TCA 과잉이므로)',
             xy=(4, list(y_don_kn1)[int(4/10.2*300)]),
             xytext=(5.5, 68),
             arrowprops=dict(arrowstyle='->', color='#1f77b4'),
             color='#1f77b4', fontsize=9)

ax2.set_xlabel('TCA : Donor 배수 (TCA 과잉 배율)', fontsize=12)
ax2.set_ylabel('이론 수율 (%)', fontsize=12)
ax2.set_title('TCA 과잉 시 수율 변화 (K_N=1, Donor 고정)\n→ TCA 아무리 늘려도 TCA 기준 수율 안 오름',
              fontsize=12, fontweight='bold')
ax2.legend(loc='center right', fontsize=9)
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 100)
ax2.grid(True, alpha=0.3)

fig.text(0.5, 0.01,
    '결론: K_N≈1 시스템에서 수율 향상 전략 → Donor 과잉(좌) + SucP로 Pi 제거(좌, 녹색 화살표)',
    ha='center', fontsize=10.5, color='darkgreen', style='italic')

plt.tight_layout(rect=[0, 0.04, 1, 1])
outpath = 'ribavirin_equilibrium.png'
plt.savefig(outpath, dpi=150, bbox_inches='tight')
print(f"저장 완료: {outpath}")
plt.show()
