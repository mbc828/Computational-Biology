# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""残基 x 成分接触频率热图：合并 18 条 fingerprint.csv（每成分 3 重复取均值）"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = r'<MD_DATA_ROOT>\fes_ikkb_20261022'
OUT = os.path.dirname(os.path.abspath(__file__))
LIGS = ['Baicalein', 'Wogonin', 'Berberine', 'OroxylinA', 'Apigenin', 'Emodin']

data = {}
res_order = []
for lig in LIGS:
    acc = {}
    n = 0
    for rep in (1, 2, 3):
        p = os.path.join(ROOT, f'{lig}_rep{rep}', 'fingerprint.csv')
        if not os.path.exists(p):
            print('missing:', p)
            continue
        n += 1
        with open(p, encoding='utf-8-sig') as f:
            for row in csv.DictReader(f):
                rk = row.get('residue') or list(row.values())[0]
                fv = float(row.get('contact_freq') or row.get('freq') or list(row.values())[1])
                acc[rk] = acc.get(rk, 0.0) + fv
    if n:
        data[lig] = {k: v / n for k, v in acc.items()}
        for k in acc:
            if k not in res_order:
                res_order.append(k)

# 残基按自然排序（按编号）
import re
def reskey(r):
    m = re.search(r'(\d+)', r)
    return int(m.group(1)) if m else 9999
res_order.sort(key=reskey)

M = np.array([[data.get(lig, {}).get(r, 0.0) for r in res_order] for lig in LIGS])
fig, ax = plt.subplots(figsize=(max(8, len(res_order) * 0.42), 4.2))
im = ax.imshow(M, aspect='auto', cmap='YlGnBu', vmin=0, vmax=1)
ax.set_xticks(range(len(res_order)))
ax.set_xticklabels(res_order, rotation=90, fontsize=8)
ax.set_yticks(range(len(LIGS)))
ax.set_yticklabels(LIGS, fontsize=10, fontstyle='italic')
ax.set_xlabel('IKKβ residue')
ax.set_ylabel('Ligand')
ax.set_title('Per-residue contact frequency (mean of 3 replicas, 50–100 ns, cutoff 4.5 Å)',
             fontsize=10)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
cb.set_label('Contact frequency', fontsize=9)
fig.tight_layout()
out = os.path.join(OUT, 'fig_fingerprint_heatmap.png')
fig.savefig(out, dpi=300)
print('saved:', out, M.shape)
