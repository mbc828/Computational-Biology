# Public source snapshot: replace path placeholders before running.
from pathlib import Path
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(r'<PROJECT_ROOT>/outputs/manuscript_v21_github_release')
original_path = ROOT / 'review_v87/assets/fig6_original.png'
panel_path = ROOT / 'review_v87/assets/fig6E_gradient.png'
out_path = ROOT / 'fig6_antibacterial_gradient_v87.png'

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
    'axes.linewidth': 1.0,
    'axes.labelsize': 13,
    'xtick.labelsize': 12,
    'ytick.labelsize': 11,
})
means = np.array([7.03, 6.5433333333, 5.18, 0.0])
sds = np.array([0.2787471973, 0.1350308607, 0.6846166811, 0.0])
replicates = [
    np.array([7.26, 7.11, 6.72]),
    np.array([6.68, 6.54, 6.41]),
    np.array([4.76, 4.81, 5.97]),
    np.array([0.0, 0.0, 0.0]),
]
labels = ['Control\n(new)', '16', '32', '64']
colors = ['#222222', '#C8C8C8', '#969696', '#606060']

fig, ax = plt.subplots(figsize=(13.253333, 3.666667), dpi=300)
x = np.arange(4)
bars = ax.bar(x, means, yerr=sds, capsize=5, width=0.62, color=colors,
              edgecolor='black', linewidth=0.9, error_kw={'elinewidth': 1.1, 'capthick': 1.1})
offsets = [-0.12, 0.0, 0.12]
for i, vals in enumerate(replicates):
    ax.scatter(x[i] + np.array(offsets), vals, s=32, marker='o',
               facecolors='white', edgecolors='black', linewidths=0.9, zorder=5)
for i, (m, sd) in enumerate(zip(means, sds)):
    y = m + sd + 0.18 if m > 0 else 0.24
    ax.text(x[i], y, f'{m:.2f}' if m > 0 else '0', ha='center', va='bottom', fontsize=11)
ax.set_xticks(x, labels)
ax.set_xlabel('Baicalein (μg/mL)', labelpad=8)
ax.set_ylabel('CFU/mL (×10⁸)', labelpad=8)
ax.set_ylim(0, 8.25)
ax.set_yticks(np.arange(0, 8.1, 2))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='x', length=0, pad=5)
ax.tick_params(axis='y', direction='out', length=4)
ax.text(0.985, 0.94, 'Mean ± SD; n = 3 plates', transform=ax.transAxes,
        ha='right', va='top', fontsize=11)
fig.text(0.008, 0.94, 'E', fontsize=25, fontweight='bold', ha='left', va='top')
fig.subplots_adjust(left=0.075, right=0.985, top=0.80, bottom=0.255)
fig.savefig(panel_path, dpi=300, facecolor='white')
plt.close(fig)

original = Image.open(original_path).convert('RGB')
panel = Image.open(panel_path).convert('RGB')
if panel.width != original.width:
    panel = panel.resize((original.width, panel.height), Image.Resampling.LANCZOS)
separator = 36
combined = Image.new('RGB', (original.width, original.height + separator + panel.height), 'white')
combined.paste(original, (0, 0))
combined.paste(panel, (0, original.height + separator))
combined.save(out_path, format='PNG', optimize=True, dpi=(300, 300))
print(out_path)
print('dimensions:', combined.size, 'bytes:', out_path.stat().st_size)
