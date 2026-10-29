# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""ELISA 新版三指标组合图（2026-10 重测版）
分组：Control / HKSA / HKSA+BAI(6,12,24,48) / HKSA+IKK-16 / HKSA+BAI24+IKK-16 / IKK-16
浓度一律按各表内嵌公式由 OD 计算；单位 pg/mL（IL-1β 按用户前定口径 pg/mL）。
"""
import numpy as np
import openpyxl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

BASE = r'<WETLAB_DATA_ROOT>\elisa\计算pji'
OUT = r'<PROJECT_ROOT>\outputs'

# 每块板的配置：文件、sheet、公式系数、主分组 OD 区、剂量梯度 OD 区、y轴标题
CFG = {
    'IL-1β': dict(file='IL-1b.xlsx', sheet='il-1β 3.29', a=7.716, b=20.272, c=0.5275,
                  main_cols=range(3, 8), main_rows=range(20, 26),
                  dose_cols={6: 8, 12: 9, 48: 10}, dose_rows=range(20, 25)),
    'IL-6': dict(file='IL-6.xlsx', sheet='IL-6 4.6', a=0.5599, b=39.067, c=3.1522,
                 main_cols=range(3, 8), main_rows=range(20, 26),
                 dose_cols={6: 8, 12: 9, 48: 10}, dose_rows=range(20, 25)),
    'TNF-α': dict(file='TNF-a.xlsx', sheet='TNF-α 3.29', a=32.194, b=275.45, c=3.9883,
                  main_cols=range(4, 9), main_rows=range(20, 26),
                  dose_cols={6: 9, 12: 10, 48: 11}, dose_rows=range(20, 25)),
}

GROUPS = ['Control', 'IKK-16', 'HKSA', 'HKSA+IKK-16',
          'HKSA+BAI(6)', 'HKSA+BAI(12)', 'HKSA+BAI(24)', 'HKSA+BAI(48)',
          'HKSA+BAI(24)+IKK-16']
COLORS = ['#4d4d4d', '#9ecae1', '#c00000', '#6baed6',
          '#c7e9c0', '#74c476', '#238b45', '#00441b',
          '#3182bd']


def get(ws, r, c):
    v = ws.cell(row=r, column=c).value
    return float(v) if v not in (None, '', ' ') else np.nan


def extract(key):
    cfg = CFG[key]
    ws = openpyxl.load_workbook(BASE + '\\' + cfg['file'], data_only=True)[cfg['sheet']]
    a, b, c0 = cfg['a'], cfg['b'], cfg['c']
    f = lambda od: a * od * od + b * od + c0
    # 主分组：行序 = Control, HKSA, IKK-16(不诱导), HKSA+IKK-16, HKSA+BAI24, HKSA+BAI24+IKK-16
    rows = list(cfg['main_rows'])
    data = {}
    names = ['Control', 'HKSA', 'IKK-16', 'HKSA+IKK-16', 'HKSA+BAI(24)', 'HKSA+BAI(24)+IKK-16']
    for i, r in enumerate(rows):
        ods = [get(ws, r, c) for c in cfg['main_cols']]
        data[names[i]] = [f(od) for od in ods if not np.isnan(od)]
    for dose, col in cfg['dose_cols'].items():
        ods = [get(ws, r, col) for r in cfg['dose_rows']]
        data[f'HKSA+BAI({dose})'] = [f(od) for od in ods if not np.isnan(od)]
    return data


def welch(x, y):
    t = stats.ttest_ind(x, y, equal_var=True)
    return min(t.pvalue * 5, 1.0)  # 5 个计划内对比的 Bonferroni 校正


def star(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


all_data = {k: extract(k) for k in CFG}

# ---- 打印汇总 ----
for k in CFG:
    d = all_data[k]
    print(f'==== {k} (pg/mL) ====')
    for g in GROUPS:
        v = np.array(d[g])
        print(f'  {g:26s} n={len(v)}  mean={v.mean():8.2f}  sd={v.std(ddof=1):6.2f}')
    for pair in [('Control', 'HKSA'), ('HKSA', 'HKSA+BAI(24)'), ('HKSA', 'HKSA+IKK-16'),
                 ('HKSA+IKK-16', 'HKSA+BAI(24)+IKK-16'), ('IKK-16', 'Control')]:
        p = welch(d[pair[0]], d[pair[1]])
        print(f'  {pair[0]} vs {pair[1]}: p_adj={p:.4f} {star(p)}')
    hk, ct = np.mean(d['HKSA']), np.mean(d['Control'])
    b24 = np.mean(d['HKSA+BAI(24)'])
    print(f'  induction +{(hk-ct)/ct*100:.1f}%  BAI24 inhibition {(hk-b24)/(hk-ct)*100:.1f}%')

# ---- 绘图 1×3 ----
fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2), dpi=300)
plt.subplots_adjust(left=0.055, right=0.99, top=0.86, bottom=0.16, wspace=0.30)

for ax, key, L in zip(axes, ['IL-1β', 'IL-6', 'TNF-α'], 'ABC'):
    d = all_data[key]
    means = [np.mean(d[g]) for g in GROUPS]
    sems = [np.std(d[g], ddof=1) for g in GROUPS]  # mean±SD 与正文统计口径一致
    x = np.arange(len(GROUPS))
    ax.bar(x, means, yerr=sems, capsize=3, width=0.68, color=COLORS,
           edgecolor='black', linewidth=0.8, error_kw=dict(lw=1.1))
    ax.set_xticks(x)
    ax.set_xticklabels([g.replace('HKSA+', '+').replace('HKSA', 'HKSA') for g in GROUPS],
                       fontsize=8, rotation=32, ha='right')
    ax.set_ylabel(f'{key} (pg/mL)', fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(direction='out', length=3.5, width=0.9)
    ymax = max(m + s for m, s in zip(means, sems)) * 1.60
    ax.set_ylim(0, ymax)
    h = ymax * 0.030
    # 关键比较括号
    iC, iH = GROUPS.index('Control'), GROUPS.index('HKSA')
    iB24, iIKK = GROUPS.index('HKSA+BAI(24)'), GROUPS.index('HKSA+IKK-16')
    iBIKK = GROUPS.index('HKSA+BAI(24)+IKK-16')
    ytop = max(means[iC], means[iH]) + sems[iH] + ymax * 0.04
    pairs = [
        (iC, iH, star(welch(d['Control'], d['HKSA'])), ytop),
        (iH, iB24, star(welch(d['HKSA'], d['HKSA+BAI(24)'])), ytop + ymax * 0.085),
        (iH, iIKK, star(welch(d['HKSA'], d['HKSA+IKK-16'])), ytop + ymax * 0.17),
        (iIKK, iBIKK, star(welch(d['HKSA+IKK-16'], d['HKSA+BAI(24)+IKK-16'])), ytop + ymax * 0.255),
    ]
    for i, j, lab, y in pairs:
        ax.plot([i, i, j, j], [y, y + h, y + h, y], lw=1.0, c='k')
        ax.text((i + j) / 2, y + h + ymax * 0.008, lab, ha='center', fontsize=9.5)
    ax.text(-0.10, 1.04, L, transform=ax.transAxes, fontsize=15, fontweight='bold', va='top')

fig.savefig(OUT + r'\fig_elisa_hksa_3cytokines.png', dpi=300,
            bbox_inches='tight', pad_inches=0.12, facecolor='white')
print('saved fig_elisa_hksa_3cytokines.png')
