# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""六成分 IKKβ 雷达图 正式版（黄连碱 I165A ΔΔG=+3.46 已回填，实线；英文标签，供附录A.7）"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '<PROJECT_ROOT>/outputs'

# axes: |MM/GBSA|, pocket stabilization (-dRMSF), traj stability (1-RMSD/max), D166 occupancy, I165A ddG, |Vina|
data = {
    'Baicalein':  [32.82, 0.0008, 0.244, 1.000, 6.46, 9.540],
    'Wogonin':    [31.18, 0.0093, 0.250, 0.816, 0.61, 9.380],
    'Berberine':  [30.57, 0.0117, 0.437, 0.953, 7.36, 9.795],
    'Oroxylin A': [33.48, 0.0279, 0.388, 0.643, 5.98, 9.090],
    'Apigenin':   [30.41, -0.0661, 0.299, 0.994, 0.07, 9.622],
    'Coptisine':  [27.33, 0.0534, 0.634, 0.578, 3.46, 11.160],  # ΔRMSF为S1口袋−全局差分；RMSD/占有率为S1口径；ΔΔG已回填
}
axes_lab = ['Binding energy\n|MM/GBSA|', 'Pocket stabilization\n(-dRMSF)', 'Trajectory stability\n(1-RMSD/max)',
            'D166 contact\noccupancy', 'I165 anchor\ndependence (ddG)', 'Docking strength\n|Vina|']

arr = np.array(list(data.values()), dtype=float)
norm = np.full_like(arr, np.nan)
norm[:, 0] = arr[:, 0] / 33.48
norm[:, 1] = (arr[:, 1] - (-0.0661)) / (0.0534 + 0.0661)
norm[:, 2] = 1 - arr[:, 2] / 0.634
norm[:, 3] = arr[:, 3]
norm[:, 4] = arr[:, 4] / 7.36
norm[:, 5] = arr[:, 5] / 11.16

ang = np.linspace(0, 2 * np.pi, 6, endpoint=False).tolist()
ang += ang[:1]
colors6 = ['#d62728', '#1f77b4', '#2ca02c', '#9467bd', '#ff7f0e', '#17becf']

fig, axr = plt.subplots(figsize=(8.4, 7.4), subplot_kw=dict(polar=True))
for (name, row), normed, c in zip(data.items(), norm, colors6):
    v = normed.tolist() + [normed[0]]
    lw = 2.2 if name == 'Coptisine' else 1.8
    axr.plot(ang, v, color=c, lw=lw, label=name)
    axr.fill(ang, v, color=c, alpha=0.06)
axr.set_xticks(ang[:-1])
axr.set_xticklabels(axes_lab, fontsize=10)
axr.set_ylim(0, 1.05)
axr.set_yticks([0.25, 0.5, 0.75, 1.0])
axr.set_yticklabels(['0.25', '0.50', '0.75', '1.00'], fontsize=8)
axr.legend(loc='upper right', bbox_to_anchor=(1.38, 1.12), fontsize=9)
axr.set_title('Multi-dimensional computational profile of six IKKβ-binding components\n(each axis independently normalized, outward = stronger; Coptisine structural axes from batch S1)',
              fontsize=10.5, pad=26)
fig.tight_layout()
fig.savefig(OUT + '/figA7_radar_ikkb_6components_final_EN.png', dpi=300, bbox_inches='tight', pad_inches=0.2)
print('saved figA7_radar_ikkb_6components_final_EN.png')

# 同步更新源数据CSV
csv = '<PROJECT_ROOT>/outputs/manuscript_v21_github_release/tables_data/radar_ikkb_6components_20261018.csv'
txt = open(csv, encoding='utf-8').read()
txt = txt.replace('黄连碱,27.3300,-0.3191,0.6339,0.5778,,11.1600', '黄连碱,27.3300,-0.3191,0.6339,0.5778,3.4600,11.1600')
txt = txt.replace('I165A_ddG(pending)', 'I165A_ddG')
open(csv, 'w', encoding='utf-8').write(txt)
print('csv updated')
