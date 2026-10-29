# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""图A: 16成分×16端点毒理热图+跨平台一致率; 图B: IKKβ五成分MD六维雷达图"""
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

D = '<PROJECT_ROOT>/admet_xval'
OUT = '<PROJECT_ROOT>/outputs'

df = pd.read_csv(D + '/internal_admet_tox.csv')
xv = pd.read_csv(D + '/xval_comparison.csv')

# 成分按药材分组排序
order = ['Baicalein', 'Baicalin', 'Wogonin', 'Wogonoside', 'Oroxylin A', 'Chrysin', 'Apigenin',
         'Berberine', 'Coptisine', 'Palmatine', 'Jatrorrhizine', 'Magnoflorine',
         'Emodin', 'Rhein', 'Aloe-emodin', 'Chrysophanol']
cn = {'Baicalein': '黄芩素', 'Baicalin': '黄芩苷', 'Wogonin': '汉黄芩素', 'Wogonoside': '汉黄芩苷',
      'Oroxylin A': '千层纸素A', 'Chrysin': '白杨素', 'Apigenin': '芹菜素',
      'Berberine': '小檗碱', 'Coptisine': '黄连碱', 'Palmatine': '巴马汀', 'Jatrorrhizine': '药根碱',
      'Magnoflorine': '木兰花碱', 'Emodin': '大黄素', 'Rhein': '大黄酸',
      'Aloe-emodin': '芦荟大黄素', 'Chrysophanol': '大黄酚'}
df = df.set_index('Compound').loc[order].reset_index()

eps = [('Ames', 'Ames 致突变'), ('Micronucleus', '微核'), ('hERG', 'hERG'), ('DILI', 'DILI'),
       ('Hepatotoxicity', '肝毒性'), ('Nephrotoxicity', '肾毒性'), ('Mitochondrial tox', '线粒体毒性'),
       ('Synoviocyte tox', '滑膜细胞毒性'), ('Chondrocyte tox', '软骨细胞毒性'),
       ('Macrophage viability', '巨噬细胞活力'), ('CYP1A2', 'CYP1A2 抑制'), ('CYP2C9', 'CYP2C9 抑制'),
       ('CYP2C19', 'CYP2C19 抑制'), ('CYP2D6', 'CYP2D6 抑制'), ('CYP3A4', 'CYP3A4 抑制'), ('P-gp', 'P-gp 抑制')]

def num(s):
    m = re.match(r'\s*([\d.]+)', str(s))
    return float(m.group(1)) if m else np.nan

M = np.array([[num(v) for v in df[c]] for c, _ in eps])

# 跨平台一致率
def ep_short(s):
    return s.split(' ↔')[0].strip()
xv['ep'] = xv['端点对(内部↔ADMETlab2)'].map(ep_short)
agree = xv.groupby('ep')['方向一致'].apply(lambda s: (s == '✓').mean()).to_dict()
xv_eps = ['Ames', 'hERG', 'DILI', 'Hepatotoxicity', 'CYP1A2', 'CYP2C9', 'CYP2C19', 'CYP2D6', 'CYP3A4', 'P-gp']
xv_cn = {'Ames': 'Ames', 'hERG': 'hERG', 'DILI': 'DILI', 'Hepatotoxicity': '肝毒性', 'CYP1A2': 'CYP1A2',
         'CYP2C9': 'CYP2C9', 'CYP2C19': 'CYP2C19', 'CYP2D6': 'CYP2D6', 'CYP3A4': 'CYP3A4', 'P-gp': 'P-gp'}
print('一致率:', {k: round(agree.get(k, float('nan')), 3) for k in xv_eps})

# ---------- 图A ----------
fig = plt.figure(figsize=(14.5, 6.8))
gs = fig.add_gridspec(1, 2, width_ratios=[3.1, 1.15], wspace=0.12)
ax = fig.add_subplot(gs[0])
im = ax.imshow(M, cmap='YlOrRd', vmin=0, vmax=1, aspect='auto')
ax.set_xticks(range(16))
ax.set_xticklabels([cn[c] for c in order], rotation=45, ha='right', fontsize=9)
ax.set_yticks(range(len(eps)))
ax.set_yticklabels([c for _, c in eps], fontsize=9)
for i in range(len(eps) + 1):
    ax.axhline(i - 0.5, color='white', lw=0.8)
for j in range(17):
    ax.axvline(j - 0.5, color='white', lw=0.8)
# 药材分隔线
for j in (7, 12):
    ax.axvline(j - 0.5, color='black', lw=2)
# 关键信号标注
ax.add_patch(Rectangle((order.index('Berberine') - 0.5, 2 - 0.5), 1, 1, fill=False, ec='blue', lw=2.5))
ax.text(order.index('Berberine'), 2 - 0.72, '*', ha='center', va='bottom', fontsize=14, color='blue', fontweight='bold')
for c in ['Emodin', 'Rhein', 'Aloe-emodin', 'Chrysophanol']:
    ax.add_patch(Rectangle((order.index(c) - 0.5, -0.5), 1, 1, fill=False, ec='purple', lw=2))
ax.text(11.2, -1.3, '紫框：蒽醌 Ames（跨平台判阳性，与文献一致）  蓝框*：小檗碱 hERG（两平台冲突，待膜片钳）',
        fontsize=8.5, color='black', ha='left')
ax.set_title('A  内部规则模型风险概率（0–1，描述性）', fontsize=11, loc='left')
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
cb.set_label('预测概率', fontsize=9)

ax2 = fig.add_subplot(gs[1])
ypos = np.arange(len(xv_eps))[::-1]
vals = [agree.get(k, np.nan) * 100 for k in xv_eps]
colors = ['#2ca02c' if v >= 70 else ('#ff7f0e' if v >= 40 else '#d62728') for v in vals]
ax2.barh(ypos, vals, color=colors, alpha=0.85, height=0.62)
for y, v in zip(ypos, vals):
    ax2.text(v + 2, y, f'{v:.0f}%', va='center', fontsize=9)
ax2.set_yticks(ypos)
ax2.set_yticklabels([xv_cn[k] for k in xv_eps], fontsize=9)
ax2.set_xlim(0, 118)
ax2.set_xlabel('双平台方向一致率（%）', fontsize=9)
ax2.set_title('B  与 ADMETlab 2.0 一致率', fontsize=11, loc='left')
ax2.axvline(70, color='grey', ls='--', lw=0.8)
ax2.text(0.0, -0.55, '注：微核/肾毒性/线粒体/滑膜/软骨/巨噬细胞端点为\n内部独有，未经交叉验证；HIA 一致率 81%（13/16）。',
         fontsize=8.5, transform=ax2.transAxes)
fig.tight_layout()
fig.savefig(OUT + '/fig_tox_xval_heatmap.png', dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.close(fig)
print('saved fig_tox_xval_heatmap.png')

# ---------- 图B 雷达 ----------
data = {  # ΔG, -dRMSF, RMSD, D166occ, I165A ddG, |Vina|
    '黄芩素':   [32.82, 0.0008, 0.244, 1.000, 6.46, 9.540],
    '汉黄芩素': [31.18, 0.0093, 0.250, 0.816, 0.61, 9.380],
    '小檗碱':   [30.57, 0.0117, 0.437, 0.953, 7.36, 9.795],
    '千层纸素A': [33.48, 0.0279, 0.388, 0.643, 5.98, 9.090],
    '芹菜素':   [30.41, -0.0661, 0.299, 0.994, 0.07, 9.622],
}
axes_lab = ['结合能\n|MM/GBSA|', '口袋稳定化\n(−ΔRMSF)', '轨迹稳定\n(1−RMSD/max)',
            'D166 接触\n占有率', 'I165 锚定依赖\n(ΔΔG)', '对接强度\n|Vina|']
arr = np.array(list(data.values()), dtype=float)
norm = np.zeros_like(arr)
norm[:, 0] = arr[:, 0] / 33.48
norm[:, 1] = (arr[:, 1] - (-0.0661)) / (0.0279 + 0.0661)
norm[:, 2] = 1 - arr[:, 2] / 0.437
norm[:, 3] = arr[:, 3]
norm[:, 4] = arr[:, 4] / 7.36
norm[:, 5] = arr[:, 5] / 9.795

ang = np.linspace(0, 2 * np.pi, 6, endpoint=False).tolist()
ang += ang[:1]
colors5 = ['#d62728', '#1f77b4', '#2ca02c', '#9467bd', '#ff7f0e']
fig, axr = plt.subplots(figsize=(7.6, 7.0), subplot_kw=dict(polar=True))
for (name, _), row, c in zip(data.items(), norm, colors5):
    v = row.tolist() + row[:1].tolist()
    axr.plot(ang, v, color=c, lw=2, label=name)
    axr.fill(ang, v, color=c, alpha=0.06)
axr.set_xticks(ang[:-1])
axr.set_xticklabels(axes_lab, fontsize=10)
axr.set_ylim(0, 1.05)
axr.set_yticks([0.25, 0.5, 0.75, 1.0])
axr.set_yticklabels(['0.25', '0.50', '0.75', '1.00'], fontsize=8)
axr.legend(loc='upper right', bbox_to_anchor=(1.32, 1.10), fontsize=10)
axr.set_title('IKKβ 五成分多维计算画像（各轴独立归一化，外=更强）', fontsize=12, pad=24)
fig.tight_layout()
fig.savefig(OUT + '/fig_radar_ikkb_5components.png', dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.close(fig)
print('saved fig_radar_ikkb_5components.png')

# 雷达图源数据存 provenance
src = pd.DataFrame(arr, index=list(data.keys()),
                   columns=['MMGBSA_abs', 'neg_dRMSF_pocket', 'ligand_RMSD', 'D166_occupancy', 'I165A_ddG', 'Vina_abs'])
src.to_csv(OUT + '/manuscript_v21_github_release/tables_data/radar_ikkb_5components_20261018.csv')
print('radar source csv saved')
