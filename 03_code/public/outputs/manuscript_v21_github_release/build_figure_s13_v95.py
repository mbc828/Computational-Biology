# -*- coding: utf-8 -*-
"""Rebuild Supplementary Figure S13 from the audited rank-1 redocking RMSD values."""
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('Figure_S13.png')
out.parent.mkdir(parents=True, exist_ok=True)

records = [
    ('IKKβ (4KIK)', 4.46),
    ('NLRP3 (7PZC)', 5.78),
    ('CASP1 (6PZP)', 6.40),
    ('TLR4–MD2 (3FXI)', 8.47),
    ('IL-17A (5HI4)', 8.56),
]
labels = [x[0] for x in records]
values = np.array([x[1] for x in records])
y = np.arange(len(records))

plt.rcParams.update({'font.family': 'Arial', 'axes.unicode_minus': False})
fig, ax = plt.subplots(figsize=(8.0, 7.37), dpi=300)
fig.subplots_adjust(left=0.30, right=0.955, top=0.82, bottom=0.16)
ax.barh(y, values, color='#c96f74', height=0.58)
ax.set_yticks(y, labels, fontsize=12)
ax.invert_yaxis()
ax.set_xlim(0, 10.0)
ax.set_ylim(4.85, -0.72)
ax.set_xticks(np.arange(0, 11, 2))
ax.tick_params(axis='x', labelsize=11, width=1.2, length=7)
ax.tick_params(axis='y', width=1.2, length=7)
ax.axvline(2.0, color='#20252b', linestyle=(0, (4, 3)), linewidth=1.35)
for yi, value in zip(y, values):
    ax.text(value + 0.13, yi, f'{value:.2f}', va='center', fontsize=11.5, color='#20252b')
ax.set_xlabel('Rank-1 redock RMSD (Å)', fontsize=13, labelpad=12)
ax.set_title('Redocking benchmark', loc='left', fontsize=24, fontweight='bold', pad=30, color='#20252b')
ax.text(9.82, 4.48, 'dashed line: 2.0 Å criterion — none passed',
        ha='right', va='center', fontsize=10.2, color='#bd6269')
for side in ('top', 'right'):
    ax.spines[side].set_visible(False)
for side in ('left', 'bottom'):
    ax.spines[side].set_color('#20252b')
    ax.spines[side].set_linewidth(1.25)

fig.savefig(out, dpi=300, facecolor='white')
plt.close(fig)
print(out)
print(out.stat().st_size)
