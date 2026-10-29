# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""IKKβ 配体诱导动力学稳定化 + 接触占有指纹 + SAR 能量分解（2026-10-12）
数据全部为既有轨迹/能量批次的后处理，无新 MD。
输出：outputs/动力学稳定化_SAR_20261012/ 下的 F12 三联图、TSV 与摘要 JSON。
"""
import json, re, csv, io
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = ['Microsoft YaHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

ANA = Path(r'<MD_DATA_ROOT>\analysis_all_completed\per_rep')
PBC = Path(r'<MD_DATA_ROOT>\MD_SHXD_PJI\05_analysis\core_3rep_mmpbsa_pbc\per_rep')
FILL = Path(r'<MD_DATA_ROOT>\mmgbsa_ikkb_fill_20261003\per_rep')
ALA = Path(r'<MD_DATA_ROOT>\mmgbsa_ala_scan_20261004\ALA_SCAN_SUMMARY.csv')
OUT = Path(r'<PROJECT_ROOT>\outputs\动力学稳定化_SAR_20261012')
OUT.mkdir(parents=True, exist_ok=True)

SYSTEMS = {
    'apo':       [('5090', 'v2_IKKb_apo', r) for r in (1, 2, 3)],
    'Baicalein': [('5090', 'v2_IKKb_Baicalein', r) for r in (1, 2, 3)],
    'Wogonin':   [('5090', 'v2_IKKb_Wogonin', r) for r in (1, 2, 3)],
    'Berberine': [('5090', 'v2_IKKb_Berberine', r) for r in (1, 2, 3)],
    'OroxylinA': [('ikkb_fill_20261003', 'v2_IKKb_OroxylinA', r) for r in (1, 2, 3)],
    'Apigenin':  [('5090', 'v2_IKKb_Apigenin', 1), ('ikkb_fill_20261003', 'v2_IKKb_Apigenin', 2), ('ikkb_fill_20261003', 'v2_IKKb_Apigenin', 3)],
}
CN = dict(Baicalein='黄芩素', Wogonin='汉黄芩素', Berberine='小檗碱', OroxylinA='千层纸素A', Apigenin='芹菜素')
POCKET = [21, 29, 44, 96, 98, 99, 152, 165, 166]  # L21 V29 K44 M96 Y98 C99 V152 I165 D166
POCKET_LABEL = {21: 'L21', 29: 'V29', 44: 'K44', 96: 'M96', 98: 'Y98', 99: 'C99', 152: 'V152', 165: 'I165', 166: 'D166'}

def read_xvg(p):
    d = np.loadtxt(p, comments=['#', '@'])
    return d[:, 0].astype(int), d[:, 1]

def rep_dir(src, system, rep):
    return ANA / src / system / ('rep%d' % rep)

# ---------- 1. ΔRMSF ----------
rmsf = {}
for sysname, reps in SYSTEMS.items():
    arr = []
    for src, s, r in reps:
        f = rep_dir(src, s, r) / 'rmsf_ca_steady.xvg'
        assert f.is_file(), f
        idx, val = read_xvg(f)
        arr.append(pd.Series(val, index=idx) if False else (idx, val))
    # 对齐残基号
    base = arr[0][0]
    for idx, _ in arr:
        assert np.array_equal(idx, base), '%s 残基号不一致' % sysname
    M = np.vstack([v for _, v in arr])  # 3 × nres
    rmsf[sysname] = dict(res=base, mean=M.mean(0), sd=M.std(0, ddof=1))

res = rmsf['apo']['res']
dRMSF = {}
for s in SYSTEMS:
    if s == 'apo':
        continue
    d = rmsf[s]['mean'] - rmsf['apo']['mean']
    sd = np.sqrt(rmsf[s]['sd'] ** 2 + rmsf['apo']['sd'] ** 2)
    dRMSF[s] = dict(d=d, sd=sd)

pocket_mask = np.isin(res, POCKET)
drmsf_summary = {}
for s, v in dRMSF.items():
    sig = np.abs(v['d']) > 2 * v['sd']
    drmsf_summary[s] = dict(
        pocket_mean=float(v['d'][pocket_mask].mean()),
        global_mean=float(v['d'].mean()),
        n_stabilized=int(((v['d'] < 0) & sig).sum()),
        n_total=int(len(res)))

# ---------- 2. 接触占有指纹 ----------
occ = {}   # sys -> resnr -> (mean_occ, n_reps_ge_0.5)
for sysname, reps in SYSTEMS.items():
    if sysname == 'apo':
        continue
    per_rep = []
    for src, s, r in reps:
        f = rep_dir(src, s, r) / 'contact_occupancy_residue.tsv'
        assert f.is_file(), f
        rows = list(csv.DictReader(io.open(f, encoding='utf-8'), delimiter='\t'))
        per_rep.append({int(x['resnr']): float(x['occupancy']) for x in rows})
    allres = sorted(set().union(*[set(x) for x in per_rep]))
    occ[sysname] = {rn: (float(np.mean([x.get(rn, 0.0) for x in per_rep])),
                         int(sum(1 for x in per_rep if x.get(rn, 0.0) >= 0.5))) for rn in allres}

# 选取指纹展示的残基：任一体系 mean≥0.3 或口袋残基
fingerprint_res = sorted({rn for s in occ for rn, (m, _) in occ[s].items() if m >= 0.30} | set(POCKET))

# ala-scan ΔΔG
ala = {}
for r in csv.DictReader(io.open(ALA, encoding='utf-8')):
    ala[(r['system'].replace('v2_IKKb_', ''), r['mutation'])] = float(r['ddg_mean_kcal_mol'])

# ---------- 3. SAR 能量分量 ----------
def parse_delta_dat(p):
    txt = io.open(p, encoding='utf-8', errors='replace').read()
    m = re.search(r'Delta \(Complex - Receptor - Ligand\):(.*?)\Z', txt, re.S)
    blk = m.group(1)
    out = {}
    for line in blk.splitlines():
        mm = re.match(r'\s*Δ?(VDWAALS|EEL|EGB|ESURF|GGAS|GSOLV|TOTAL)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', line)
        if mm:
            out[mm.group(1)] = float(mm.group(2))
    return out

sar = {}
for sysname in ['Baicalein', 'Wogonin', 'Berberine', 'OroxylinA', 'Apigenin']:
    vals = []
    for r in (1, 2, 3):
        if sysname in ('Baicalein', 'Wogonin', 'Berberine'):
            f = PBC / ('v2_IKKb_' + sysname) / ('rep%d' % r) / 'FINAL_RESULTS_MMPBSA_PBC_51FRAMES.dat'
        else:
            f = FILL / ('v2_IKKb_' + sysname) / ('rep%d' % r) / 'FINAL_RESULTS_MMGBSA.dat'
        assert f.is_file(), f
        vals.append(parse_delta_dat(f))
    sar[sysname] = {k: (float(np.mean([v[k] for v in vals])), float(np.std([v[k] for v in vals], ddof=1)))
                    for k in ['VDWAALS', 'EEL', 'EGB', 'ESURF', 'TOTAL']}

# ---------- 图 F12 ----------
fig = plt.figure(figsize=(16, 5.4), dpi=200)
colors = dict(Baicalein='#c8363b', Wogonin='#e07b00', Berberine='#2e7d32', OroxylinA='#7b5ea7', Apigenin='#b0762a')

# A: ΔRMSF 剖面
axA = fig.add_subplot(1, 3, 1)
for s, v in dRMSF.items():
    axA.plot(res, v['d'], color=colors[s], lw=1.1, label=CN[s], alpha=0.9)
axA.axhline(0, color='black', lw=0.8)
top_labels = [165, 166, 152, 99]
row = 0
for rn in POCKET:
    axA.axvline(rn, color='#999999', lw=0.6, ls=':', zorder=0)
for k, rn in enumerate(POCKET):
    y0, y1 = (-0.30, -0.22) if rn in top_labels else (-0.22, -0.145)
    axA.annotate(POCKET_LABEL[rn], xy=(rn, 0), xytext=(rn, y0 + (0 if k % 2 else 0.0)),
                 fontsize=6.8, rotation=90, va='bottom', ha='center', color='#555555',
                 annotation_clip=False)
axA.set_ylim(-0.31, 0.20)
axA.set_xlabel('残基编号（配体接触残基全部 ≤175）')
axA.set_ylabel('ΔRMSF = holo − apo (nm)')
axA.set_title('A  配体诱导的柔性变化（三重复均值）', fontsize=11)
axA.legend(fontsize=7.5, frameon=False, loc='upper right')

# B: 接触占有热图
axB = fig.add_subplot(1, 3, 2)
sysnames = ['Baicalein', 'Wogonin', 'Berberine', 'OroxylinA', 'Apigenin']
M = np.array([[occ[s].get(rn, (0, 0))[0] for rn in fingerprint_res] for s in sysnames])
im = axB.imshow(M, cmap='YlGnBu', vmin=0, vmax=1, aspect='auto')
axB.set_xticks(range(len(fingerprint_res)))
labels = []
for rn in fingerprint_res:
    lab = POCKET_LABEL.get(rn, str(rn))
    labels.append(lab)
axB.set_xticklabels(labels, rotation=90, fontsize=7.5)
for j, rn in enumerate(fingerprint_res):
    if rn in (165, 166):
        axB.get_xticklabels()[j].set_color('#c8363b')
        axB.get_xticklabels()[j].set_fontweight('bold')
axB.set_yticks(range(5)); axB.set_yticklabels([CN[s] for s in sysnames], fontsize=9)
for i in range(5):
    for j, rn in enumerate(fingerprint_res):
        n_consistent = occ[sysnames[i]].get(rn, (0, 0))[1]
        if M[i, j] >= 0.5:
            axB.text(j, i, '●' if n_consistent == 3 else '○', ha='center', va='center',
                     fontsize=6, color='white' if M[i, j] > 0.7 else '#333333')
axB.set_title('B  配体—残基接触占有率（●=三重复一致）', fontsize=11)
axB.set_xlabel('红标为 ala-scan 验证位点：I165（五配体 ΔΔG +2.0~+2.9）；D166（选择性影响芹菜素）', fontsize=7.5)
fig.colorbar(im, ax=axB, fraction=0.04, label='占有率')

# C: SAR 能量分量
axC = fig.add_subplot(1, 3, 3)
comp_names = ['VDWAALS', 'EEL', 'EGB', 'ESURF']
comp_cn = ['范德华', '静电', '极性溶剂化(不利)', '非极性溶剂化']
x = np.arange(5); w = 0.19
cmap4 = ['#4a6fa5', '#c8363b', '#e07b00', '#2e7d32']
for k, (cn_k, lab) in enumerate(zip(comp_names, comp_cn)):
    means = [sar[s][cn_k][0] for s in sysnames]
    sds = [sar[s][cn_k][1] for s in sysnames]
    axC.bar(x + (k - 1.5) * w, means, w, yerr=sds, capsize=2, color=cmap4[k], label=lab, alpha=0.9)
axC.plot(x, [sar[s]['TOTAL'][0] for s in sysnames], 'kD-', ms=5, lw=1.2, label='ΔG 合计')
axC.set_xticks(x); axC.set_xticklabels([CN[s] for s in sysnames], fontsize=8.5, rotation=18)
axC.axhline(0, color='black', lw=0.8)
axC.set_ylabel('能量分量 (kcal/mol)')
axC.set_title('C  结合能量分解的化学基础（n=3）', fontsize=11)
axC.legend(fontsize=7.5, frameon=False)

fig.tight_layout()
fpng = OUT / 'F12_动力学稳定化_接触指纹_SAR.png'
fig.savefig(fpng, dpi=200, facecolor='white')
print('saved', fpng)

# ---------- 表与摘要 ----------
with io.open(OUT / 'drmsf_summary.tsv', 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['system', 'pocket_mean_dRMSF_nm', 'global_mean_dRMSF_nm', 'n_stabilized_residues', 'n_residues'])
    for s, v in drmsf_summary.items():
        w.writerow([s, round(v['pocket_mean'], 4), round(v['global_mean'], 4), v['n_stabilized'], v['n_total']])

with io.open(OUT / 'occupancy_matrix.tsv', 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['resnr', 'label'] + sysnames)
    for j, rn in enumerate(fingerprint_res):
        w.writerow([rn, POCKET_LABEL.get(rn, '')] + [round(float(M[i, j]), 3) for i in range(5)])

with io.open(OUT / 'sar_components.tsv', 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['system', 'component', 'mean', 'sd'])
    for s in sysnames:
        for k in ['VDWAALS', 'EEL', 'EGB', 'ESURF', 'TOTAL']:
            w.writerow([s, k, round(sar[s][k][0], 2), round(sar[s][k][1], 2)])

summary = dict(drmsf=drmsf_summary,
               d166_occupancy={s: occ[s].get(166, (0, 0))[0] for s in sysnames},
               i165_occupancy={s: occ[s].get(165, (0, 0))[0] for s in sysnames},
               ala_scan_ddg={'%s_%s' % k: v for k, v in ala.items()},
               sar_total={s: sar[s]['TOTAL'] for s in sysnames})
(OUT / 'analysis_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=1), encoding='utf-8')

print('\nΔRMSF 口袋区均值：')
for s, v in drmsf_summary.items():
    print('  %-10s 口袋 %+.4f nm | 全局 %+.4f | 显著稳定化残基 %d/%d' % (s, v['pocket_mean'], v['global_mean'], v['n_stabilized'], v['n_total']))
print('\nI165 占有率：', {CN[s]: round(occ[s].get(165, (0, 0))[0], 2) for s in sysnames})
print('D166 占有率：', {CN[s]: round(occ[s].get(166, (0, 0))[0], 2) for s in sysnames})
print('\nSAR ΔG 分量（均值）：')
for s in sysnames:
    print('  %-10s vdW %+.1f | EEL %+.1f | EGB %+.1f | ESURF %+.1f | TOTAL %+.1f' %
          (s, sar[s]['VDWAALS'][0], sar[s]['EEL'][0], sar[s]['EGB'][0], sar[s]['ESURF'][0], sar[s]['TOTAL'][0]))
