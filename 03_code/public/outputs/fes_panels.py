# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""逐成分 FES 面板：fes2d_<lig>.png（等值面）+ fes3d_<lig>.png（3D 曲面）
供四联复合图拼装。数据：<MD_DATA_ROOT>\fes_ikkb_20261022\<Ligand>_rep<1-3>
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
from scipy.ndimage import gaussian_filter
import os

BASE = r'<MD_DATA_ROOT>\fes_ikkb_20261022'
OUT = r'<PROJECT_ROOT>\outputs\fes_panels'
os.makedirs(OUT, exist_ok=True)
LIGS = ['Baicalein', 'Wogonin', 'Berberine', 'OroxylinA', 'Apigenin', 'Emodin']

def read_xvg(p):
    t, v = [], []
    for line in open(p, encoding='utf-8', errors='ignore'):
        if line.startswith(('#', '@')) or not line.strip():
            continue
        a = line.split()
        t.append(float(a[0])); v.append(float(a[1]))
    return np.array(t), np.array(v)

def pooled(lig):
    rms, cnt = [], []
    for rep in (1, 2, 3):
        d = os.path.join(BASE, f'{lig}_rep{rep}')
        if not os.path.exists(os.path.join(d, 'rmsd_lig.xvg')):
            continue
        _, r = read_xvg(os.path.join(d, 'rmsd_lig.xvg'))
        _, c = read_xvg(os.path.join(d, 'contacts.xvg'))
        n = min(len(r), len(c))
        rms.append(r[:n]); cnt.append(c[:n])
    return np.concatenate(rms), np.concatenate(cnt)

plt.rcParams['font.family'] = 'Arial'

for lig in LIGS:
    r, c = pooled(lig)
    H, xe, ye = np.histogram2d(r, c, bins=(42, 42))
    Hs = gaussian_filter(H, sigma=1.2)
    Ps = Hs / Hs.sum()
    Pmax = Ps.max()
    with np.errstate(divide='ignore'):
        F = -np.log(Ps / Pmax)
    F = np.clip(F, 0, 8)
    mask = Ps < Pmax * np.exp(-8.2)   # 8 kBT 以外区域留白
    Fm = np.ma.masked_array(F, mask=mask)
    xc = 0.5 * (xe[:-1] + xe[1:]); yc = 0.5 * (ye[:-1] + ye[1:])
    X, Y = np.meshgrid(xc, yc)
    lv = np.arange(0, 8.5, 0.5)

    # 2D 等值面
    fig, ax = plt.subplots(figsize=(4.0, 3.6), dpi=300)
    cf = ax.contourf(X, Y, Fm.T, levels=lv, cmap='turbo', extend='max')
    ax.contour(X, Y, Fm.T, levels=lv[::2], colors='k', linewidths=0.15, alpha=0.4)
    ax.set_xlabel('Ligand RMSD (nm)', fontsize=9)
    ax.set_ylabel('Pocket contacts', fontsize=9)
    ax.set_title('Gibbs Energy Landscape', fontsize=9.5)
    ax.tick_params(direction='out', labelsize=8)
    cb = fig.colorbar(cf, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label('kBT', fontsize=8); cb.ax.tick_params(labelsize=7)
    fig.savefig(os.path.join(OUT, f'fes2d_{lig}.png'), bbox_inches='tight',
                pad_inches=0.06, facecolor='white')
    plt.close(fig)

    # 3D 曲面（参考图样式：F 向下为谷 → 用 -F 使谷为漏斗）
    fig = plt.figure(figsize=(4.4, 3.8), dpi=300)
    ax3 = fig.add_subplot(111, projection='3d')
    ax3.plot_surface(X, Y, F.T, cmap='turbo', rstride=1, cstride=1,
                     linewidth=0, antialiased=True, vmin=0, vmax=8)
    ax3.view_init(elev=32, azim=-128)
    ax3.set_xlabel('RMSD', fontsize=7.5, labelpad=-4)
    ax3.set_ylabel('Contacts', fontsize=7.5, labelpad=-4)
    ax3.set_zlabel('kBT', fontsize=7.5, labelpad=-4)
    ax3.tick_params(labelsize=6, pad=-3)
    ax3.set_box_aspect((1.15, 1.15, 0.75))
    fig.savefig(os.path.join(OUT, f'fes3d_{lig}.png'), bbox_inches='tight',
                pad_inches=0.04, facecolor='white')
    plt.close(fig)
    print('panels done:', lig)
print('ALL_PANELS_OK', OUT)
