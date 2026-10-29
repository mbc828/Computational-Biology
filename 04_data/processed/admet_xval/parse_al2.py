# -*- coding: utf-8 -*-
"""parse_al2.py — 解析 16 个 ADMETlab 2.0 结果页 → 主表 CSV
输出 admet_xval/admetlab2_master.csv：Compound, Category, Endpoint, Value, FlagColor"""
import csv, re
from pathlib import Path

RES = Path('admet_xval/results_al2')

def parse(path):
    html = path.read_text(encoding='utf-8', errors='ignore')
    # 类别区块标题
    out = []
    # 按 <tr> 切分，记录最近出现的类别标题
    blocks = re.split(r'(<div class="alert alert-info[^"]*"[^>]*>)', html)
    cur_cat = ''
    for seg in blocks:
        m = re.match(r'<div class="alert alert-info[^"]*"[^>]*>\s*([^<]+?)\s*</div>', seg)
        if m:
            cur_cat = m.group(1).strip(); continue
        # 数值端点：<td class="text-left">NAME</td> <td>VALUE</td> <td> <i class="fas fa-circle text-COLOR">
        for nm, val, color in re.findall(
            r'<td class="text-left"[^>]*>([^<]+)</td>\s*<td[^>]*>\s*([^<]+?)\s*</td>\s*<td[^>]*>\s*<i class="fas fa-circle text-(\w+)"',
            seg):
            out.append((cur_cat, nm.strip(), val.strip(), color))
        # 分类端点：<td class="text-left">NAME</td> <td> <span class="font-weight-bold">LABEL</span>
        # 第三格色点可选（CYP 端点无色点，为 <td></td>）
        for nm, lab, color in re.findall(
            r'<td class="text-left"[^>]*>([^<]+)</td>\s*<td[^>]*>\s*<span class="font-weight-bold">([^<]*)</span>\s*</td>\s*(?:<td[^>]*>\s*<i class="fas fa-circle text-(\w+)")?',
            seg):
            out.append((cur_cat, nm.strip(), lab.strip(), color))
    return out

rows = []
for f in sorted(RES.glob('*.html')):
    comp = f.stem
    got = parse(f)
    for cat, nm, val, color in got:
        rows.append((comp, cat, nm, val, color))
    print(f'{comp}: {len(got)} 端点')

with open('admet_xval/admetlab2_master.csv', 'w', newline='', encoding='utf-8-sig') as fp:
    w = csv.writer(fp)
    w.writerow(['Compound', 'Category', 'Endpoint', 'Value', 'FlagColor'])
    w.writerows(rows)
print('主表:', len(rows), '行 ->', 'admet_xval/admetlab2_master.csv')
eps = sorted({r[2] for r in rows})
print('端点清单(%d):' % len(eps))
print(' | '.join(eps))
