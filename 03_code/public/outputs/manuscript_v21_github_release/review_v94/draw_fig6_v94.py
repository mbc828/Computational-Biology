# Public source snapshot: replace <PROJECT_ROOT> with the workspace root before running.
from pathlib import Path
import csv
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(r'<PROJECT_ROOT>/outputs/manuscript_v21_github_release')
FIGURE = ROOT / '02_figures/main/Figure_06.png'
DATA = ROOT / '04_data/processed/wetlab/antibacterial_gradient_reporting_v94.csv'

# Figure 6 v93 comprises panels A–D (3976 × 2330 px), a 36-px separator,
# and panel E (3976 × 1100 px). Only panel E is regenerated here.
TOP_HEIGHT = 2330
PANEL_HEIGHT = 1100
SEPARATOR = 36
EXPECTED_SIZE = (3976, 3466)

with DATA.open(encoding='utf-8-sig', newline='') as handle:
    rows = list(csv.DictReader(handle))
assert len(rows) == 4
counts = [np.array([int(row[f'plate_{i}_colonies']) for i in (1, 2, 3)], dtype=float) for row in rows]
assert np.array_equal(counts[0], [726, 711, 672])
assert np.array_equal(counts[1], [668, 654, 641])
assert np.array_equal(counts[2], [476, 481, 597])
assert np.array_equal(counts[3], [0, 0, 0])
assert rows[3]['interpretation'] == 'not_detected'

# One colony corresponds to 1e6 CFU/mL under the recorded plating/dilution rule.
means = np.array([row.mean() * 0.01 for row in counts])  # units: 1e8 CFU/mL
sds = np.array([row.std(ddof=1) * 0.01 for row in counts])
control_mean_cfu = counts[0].mean() * 1e6
lod_cfu = 1e6
survival_upper_pct = lod_cfu / control_mean_cfu * 100
assert survival_upper_pct < 0.15
assert np.isclose(means[0], 7.03)
assert np.isclose(means[1], 6.5433333333)
assert np.isclose(means[2], 5.18)

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
    'axes.linewidth': 1.0,
    'axes.labelsize': 13,
    'xtick.labelsize': 12,
    'ytick.labelsize': 11,
})
labels = ['Control\n(new)', '16', '32', '64']
colors = ['#222222', '#C8C8C8', '#969696']
fig, ax = plt.subplots(figsize=(13.253333, 3.666667), dpi=300)
x = np.arange(4)
ax.bar(x[:3], means[:3], yerr=sds[:3], capsize=5, width=0.62, color=colors,
       edgecolor='black', linewidth=0.9,
       error_kw={'elinewidth': 1.1, 'capthick': 1.1})
offsets = [-0.12, 0.0, 0.12]
for i in range(3):
    vals = counts[i] * 0.01
    ax.scatter(x[i] + np.array(offsets), vals, s=32, marker='o',
               facecolors='white', edgecolors='black', linewidths=0.9, zorder=5)
for i in range(3):
    ax.text(x[i], means[i] + sds[i] + 0.18, f'{means[i]:.2f}',
            ha='center', va='bottom', fontsize=11)
ax.text(x[3], 0.20, 'ND\n<0.01\nSurvival <0.15%', ha='center', va='bottom',
        fontsize=10.5, fontweight='bold', color='#9E2A2B')
ax.set_xticks(x, labels)
ax.set_xlim(-0.55, 3.55)
ax.set_xlabel('Baicalein (μg/mL)', labelpad=8)
ax.set_ylabel('CFU/mL (×10⁸)', labelpad=8)
ax.set_ylim(0, 8.25)
ax.set_yticks(np.arange(0, 8.1, 2))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='x', length=0, pad=5)
ax.tick_params(axis='y', direction='out', length=4)
ax.text(0.985, 0.94,
        'Mean ± SD; n = 3 plates\nND = not detected; LOD = 0.01 ×10⁸ CFU/mL',
        transform=ax.transAxes, ha='right', va='top', fontsize=10.5)
fig.text(0.008, 0.94, 'E', fontsize=25, fontweight='bold', ha='left', va='top')
fig.subplots_adjust(left=0.075, right=0.985, top=0.80, bottom=0.255)
fig.canvas.draw()
panel = Image.frombuffer('RGBA', fig.canvas.get_width_height(), fig.canvas.buffer_rgba()).convert('RGB')
plt.close(fig)
if panel.size != (EXPECTED_SIZE[0], PANEL_HEIGHT):
    panel = panel.resize((EXPECTED_SIZE[0], PANEL_HEIGHT), Image.Resampling.LANCZOS)
assert panel.size == (EXPECTED_SIZE[0], PANEL_HEIGHT), panel.size

source = Image.open(FIGURE).convert('RGB')
assert source.size == EXPECTED_SIZE, source.size
original = source.crop((0, 0, EXPECTED_SIZE[0], TOP_HEIGHT))
combined = Image.new('RGB', EXPECTED_SIZE, 'white')
combined.paste(original, (0, 0))
combined.paste(panel, (0, TOP_HEIGHT + SEPARATOR))
combined.save(FIGURE, format='PNG', optimize=True, dpi=(300, 300))
print(FIGURE)
print('dimensions:', combined.size, 'bytes:', FIGURE.stat().st_size)
print('64 μg/mL: ND; <1e6 CFU/mL; relative survival <0.15%')
