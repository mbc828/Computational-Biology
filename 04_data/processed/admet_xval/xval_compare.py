# -*- coding: utf-8 -*-
"""xval_compare.py — 内部规则模型 vs ADMETlab 2.0 一致性对照
输出: admet_xval/xval_comparison.csv (逐端点逐成分) + 汇总打印
方向判定: 内部概率≥0.5 判阳性; ADMETlab2 +/++/+++ 判阳性, -/--/--- 判阴性
GI吸收: 按 ADMETlab2 HIA 定义, 吸收<30% 判 'HIA+(差吸收)'
"""
import csv, re
from pathlib import Path

# ---------- 读内部模型 ----------
internal = {}   # internal[compound][endpoint] = raw string
with open('admet_xval/internal_admet_tox.csv', encoding='utf-8-sig') as f:
    rd = csv.DictReader(f)
    for row in rd:
        internal[row['Compound']] = row

# ---------- 读 ADMETlab2 ----------
al2 = {}        # al2[(compound, endpoint)] = label
with open('admet_xval/admetlab2_master.csv', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        al2[(row['Compound'], row['Endpoint'])] = row['Value'].strip()

def fnum(s):
    m = re.search(r'([\d.]+)', str(s))
    return float(m.group(1)) if m else None

def al2_pos(label):
    return 1 if label.startswith('+') else 0

COMPS = ['Baicalein','Berberine','Wogonoside','Wogonin','Baicalin','Oroxylin A',
         'Apigenin','Chrysin','Coptisine','Palmatine','Jatrorrhizine','Magnoflorine',
         'Emodin','Rhein','Aloe-emodin','Chrysophanol']

# (内部列, ADMETlab2端点, 内部阳性规则函数)
PAIRS = [
    ('Ames',          'AMES Toxicity',      lambda v: int(fnum(v) >= 0.5)),
    ('hERG',          'hERG Blockers',      lambda v: int(fnum(v) >= 0.5)),
    ('DILI',          'DILI',               lambda v: int(fnum(v) >= 0.5)),
    ('Hepatotoxicity','H-HT',               lambda v: int(fnum(v) >= 0.5)),
    ('CYP1A2',        'CYP1A2 inhibitor',   lambda v: int(fnum(v) >= 0.5)),
    ('CYP2C9',        'CYP2C9 inhibitor',   lambda v: int(fnum(v) >= 0.5)),
    ('CYP2C19',       'CYP2C19 inhibitor',  lambda v: int(fnum(v) >= 0.5)),
    ('CYP2D6',        'CYP2D6 inhibitor',   lambda v: int(fnum(v) >= 0.5)),
    ('CYP3A4',        'CYP3A4 inhibitor',   lambda v: int(fnum(v) >= 0.5)),
    ('P-gp',          'Pgp-inhibitor',      lambda v: int(fnum(v) >= 0.5)),
    ('GI absorption', 'HIA',                lambda v: int(fnum(v) < 30)),  # HIA+ = 吸收<30%
]

rows, pair_stat = [], {}
for comp in COMPS:
    for icol, aep, rule in PAIRS:
        iv = internal[comp][icol]
        av = al2.get((comp, aep), '?')
        ipos = rule(iv)
        apos = al2_pos(av) if av != '?' else None
        agree = (ipos == apos) if apos is not None else None
        rows.append([comp, f'{icol} ↔ {aep}', iv, '阳性' if ipos else '阴性',
                     av, '阳性' if apos else '阴性', '✓' if agree else '✗'])
        key = f'{icol} ↔ {aep}'
        pair_stat.setdefault(key, [0, 0])
        pair_stat[key][1] += 1
        if agree: pair_stat[key][0] += 1

with open('admet_xval/xval_comparison.csv', 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['Compound', '端点对(内部↔ADMETlab2)', '内部值', '内部判定', 'AL2值', 'AL2判定', '方向一致'])
    w.writerows(rows)

total_a = sum(v[0] for v in pair_stat.values()); total_n = sum(v[1] for v in pair_stat.values())
print('=== 逐端点一致率 ===')
for k, (a, n) in pair_stat.items():
    print(f'{k:32s} {a}/{n}  ({100*a/n:.0f}%)')
print(f'{"总体":34s} {total_a}/{total_n} ({100*total_a/total_n:.1f}%)')

# ---------- 逐成分 ----------
print('\n=== 逐成分一致率 ===')
comp_stat = {}
for r in rows:
    comp_stat.setdefault(r[0], [0, 0])
    comp_stat[r[0]][1] += 1
    if r[6] == '✓': comp_stat[r[0]][0] += 1
for c, (a, n) in comp_stat.items():
    print(f'{c:15s} {a}/{n} ({100*a/n:.0f}%)')

# ---------- PPB 数值相关 (Spearman) ----------
ippb = [fnum(internal[c]['Protein binding']) for c in COMPS]
appb = [float(al2[(c, 'PPB')].replace('%', '')) for c in COMPS]
def rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0]*len(xs)
    for rk, i in enumerate(order): r[i] = rk
    return r
import math
ri, ra = rank(ippb), rank(appb)
n = len(ri)
mi = sum(ri)/n; ma = sum(ra)/n
cov = sum((ri[i]-mi)*(ra[i]-ma) for i in range(n))
si = math.sqrt(sum((x-mi)**2 for x in ri)); sa = math.sqrt(sum((x-ma)**2 for x in ra))
rho = cov/(si*sa)
print(f'\n=== PPB 数值比较 ===')
print('内部PPB%   :', ippb)
print('AL2 PPB%   :', appb)
print(f'Spearman ρ = {rho:.3f}  (系统偏移: AL2平均 {sum(appb)/n:.1f}% vs 内部 {sum(ippb)/n:.1f}%)')

# ---------- ADMETlab2 新增独立发现 (内部模型没有的端点) ----------
print('\n=== ADMETlab2 新增风险信号 (内部模型未覆盖) ===')
for ep in ['Carcinogencity', 'Respiratory Toxicity', 'Rat Oral Acute Toxicity',
           'Skin Sensitization', 'BBB Penetration', 'NR-AhR', 'SR-p53']:
    hits = [(c, al2.get((c, ep), '?')) for c in COMPS]
    pos = [f'{c}({v})' for c, v in hits if v.startswith('+')]
    print(f'{ep:26s} 阳性: {", ".join(pos) if pos else "无"}')
