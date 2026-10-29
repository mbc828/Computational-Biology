# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""FES 二维自由能景观图：6 成分 IKKβ 复合物（50-100 ns，3 重复合并）
轴: x=配体重原子RMSD(Backbone拟合, nm), y=配体-口袋接触数(0.45 nm)
F = -kT ln(P/Pmax)，单位 kBT
读取 /mnt/m/shxdpji/fes_ikkb_20261022/<Ligand>_rep<1-3>/{rmsd_lig,contacts}.xvg
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob, os

BASE = r'<MD_DATA_ROOT>\fes_ikkb_20261022'
LIGS = ['Baicalein', 'Wogonin', 'Berberine', 'OroxylinA', 'Apigenin', 'Emodin']
NICE = {'OroxylinA': 'Oroxylin A'}

def read_xvg(p):
    t, v = [], []
    for line in open(p, encoding='utf-8', errors='ignore'):
        if line.startswith(('#', '@')) or not line.strip():
            continue
        a = line.split()
        t.append(float(a[0])); v.append(float(a[1]))
    return np.array(t), np.array(v)

fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.6), dpi=300)
for ax, lig in zip(axes.flat, LIGS):
    rms, cnt = [], []
    for rep in (1, 2, 3):
        d = os.path.join(BASE, f'{lig}_rep{rep}')
        if not os.path.exists(os.path.join(d, 'rmsd_lig.xvg')):
            continue
        _, r = read_xvg(os.path.join(d, 'rmsd_lig.xvg'))
        _, c = read_xvg(os.path.join(d, 'contacts.xvg'))
        n = min(len(r), len(c))
        rms.append(r[:n]); cnt.append(c[:n])
    r = np.concatenate(rms); c = np.concatenate(cnt)
    H, xe, ye = np.histogram2d(r, c, bins=(45, 45))
    P = H / H.sum()
    Pmax = P[P > 0].max()
    with np.errstate(divide='ignore'):
        F = -np.log(P / Pmax)
    F[~np.isfinite(F)] = np.nan
    F = np.ma.masked_invalid(F)
    xc = 0.5 * (xe[:-1] + xe[1:]); yc = 0.5 * (ye[:-1] + ye[1:])
    cf = ax.contourf(xc, yc, F.T, levels=np.arange(0, 8.5, 0.5), cmap='turbo', extend='max')
    ax.set_title(NICE.get(lig, lig), fontsize=11)
    ax.tick_params(direction='out', labelsize=8.5)
    print(lig, 'frames:', len(r), 'RMSD range: %.2f-%.2f' % (r.min(), r.max()))

for ax in axes[1, :]:
    ax.set_xlabel('Ligand RMSD (nm)', fontsize=10)
for ax in axes[:, 0]:
    ax.set_ylabel('Pocket contacts (<0.45 nm)', fontsize=10)
fig.suptitle('Free-energy landscapes of six IKKβ-bound components (50–100 ns, three replicates pooled; F in kBT)',
             fontsize=12)
fig.subplots_adjust(left=0.06, right=0.90, top=0.90, bottom=0.08, wspace=0.24, hspace=0.30)
cax = fig.add_axes([0.92, 0.12, 0.015, 0.74])
cb = fig.colorbar(cf, cax=cax)
cb.set_label('ΔF (kBT)', fontsize=10)
fig.savefig(r'<PROJECT_ROOT>\outputs\fig_new_fes_ikkb_6components.png',
            dpi=300, bbox_inches='tight', pad_inches=0.15, facecolor='white')
print('saved fig_new_fes_ikkb_6components.png')
