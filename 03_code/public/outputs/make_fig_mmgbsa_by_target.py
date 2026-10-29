# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""表4 图化：分靶点 MM/GBSA 分组柱状图（英文，SCI 风格）
新增图：MM/GBSA binding energy by target (mean ± sample SD, n=3 trajectories)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False
OUT = r'<PROJECT_ROOT>\outputs'

data = {
    'IKKβ': [('Baicalein', 32.82, 3.51), ('Wogonin', 31.18, 2.59), ('Berberine', 30.57, 2.47),
             ('Oroxylin A', 33.48, 3.96), ('Apigenin', 30.41, 2.80), ('Coptisine', 27.33, 2.61),
             ('Emodin', 24.77, 3.86), ('IKK-16 (anchor)', 35.66, 4.50)],
    'NLRP3': [('Baicalein', 27.61, 2.47), ('Berberine', 27.64, 2.87)],
    'TLR4–MD2': [('Baicalein', 27.02, 3.53), ('Oroxylin A', 27.03, 2.52), ('Wogonin', 27.82, 1.25)],
}
colors = {'Baicalein': '#d62728', 'Wogonin': '#1f77b4', 'Berberine': '#2ca02c',
          'Oroxylin A': '#9467bd', 'Apigenin': '#ff7f0e', 'Coptisine': '#17becf',
          'Emodin': '#8c564b', 'IKK-16 (anchor)': '0.75'}

fig, ax = plt.subplots(figsize=(10.5, 5.6), dpi=300)
xpos = 0
xticks, xlabels = [], []
for ti, (target, rows) in enumerate(data.items()):
    start = xpos
    for name, m, s in rows:
        hatch = '//' if 'IKK-16' in name else ''
        ax.bar(xpos, m, yerr=s, capsize=3, width=0.78, color=colors[name],
               edgecolor='black', linewidth=0.9, hatch=hatch, error_kw=dict(lw=1.0))
        xpos += 1
    xticks.append((start + xpos - 1) / 2)
    xlabels.append(target + f'\n(n={len(rows)})')
    xpos += 0.9  # 靶点间距
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels, fontsize=11)
ax.set_ylabel('|ΔG| MM/GBSA (kcal/mol)', fontsize=11.5)
ax.spines[['top', 'right']].set_visible(False)
ax.tick_params(direction='out', length=4, width=1)
handles = [plt.Rectangle((0, 0), 1, 1, facecolor=colors[k], edgecolor='black',
                         hatch='//' if 'IKK-16' in k else '') for k in colors]
ax.legend(handles, list(colors), fontsize=8.5, ncol=4, frameon=False,
          loc='upper left', bbox_to_anchor=(0.0, 1.10))
fig.tight_layout()
fig.savefig(OUT + r'\fig_new_mmgbsa_by_target.png', dpi=300, bbox_inches='tight',
            pad_inches=0.15, facecolor='white')
print('saved fig_new_mmgbsa_by_target.png')
