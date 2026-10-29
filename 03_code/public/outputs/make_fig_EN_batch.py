# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""Generate English versions of fig04 (tox heatmap), fig20 (radar), fig24 (integrated schematic).
Outputs are letterboxed onto white canvases matching the exact pixel sizes embedded in the v60 docx,
so zip-level media replacement keeps layout undistorted.
"""
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch, Circle
from PIL import Image

WS = '<PROJECT_ROOT>'
D = WS + '/admet_xval'
OUT = WS + '/outputs'

TARGETS = {  # output name -> (W,H) embedded pixel size
    'fig04_tox_xval_heatmap_EN.png': (3799, 2669),
    'fig20_radar_ikkb_5components_EN.png': (2481, 1952),
    'fig24_schematic_EN.png': (3200, 1889),
}

def letterbox(path, target):
    W, H = target
    im = Image.open(path).convert('RGB')
    s = min(W / im.width, H / im.height)
    nw, nh = int(im.width * s), int(im.height * s)
    im = im.resize((nw, nh), Image.LANCZOS)
    canvas = Image.new('RGB', (W, H), 'white')
    canvas.paste(im, ((W - nw) // 2, (H - nh) // 2))
    canvas.save(path, dpi=(300, 300))
    print('letterboxed', path, canvas.size)

# ================= fig04: tox heatmap + concordance =================
df = pd.read_csv(D + '/internal_admet_tox.csv')
xv = pd.read_csv(D + '/xval_comparison.csv')

order = ['Baicalein', 'Baicalin', 'Wogonin', 'Wogonoside', 'Oroxylin A', 'Chrysin', 'Apigenin',
         'Berberine', 'Coptisine', 'Palmatine', 'Jatrorrhizine', 'Magnoflorine',
         'Emodin', 'Rhein', 'Aloe-emodin', 'Chrysophanol']
df = df.set_index('Compound').loc[order].reset_index()

eps = [('Ames', 'Ames mutagenicity'), ('Micronucleus', 'Micronucleus'), ('hERG', 'hERG'),
       ('DILI', 'DILI'), ('Hepatotoxicity', 'Hepatotoxicity'), ('Nephrotoxicity', 'Nephrotoxicity'),
       ('Mitochondrial tox', 'Mitochondrial tox'), ('Synoviocyte tox', 'Synoviocyte tox'),
       ('Chondrocyte tox', 'Chondrocyte tox'), ('Macrophage viability', 'Macrophage viability'),
       ('CYP1A2', 'CYP1A2 inhibition'), ('CYP2C9', 'CYP2C9 inhibition'),
       ('CYP2C19', 'CYP2C19 inhibition'), ('CYP2D6', 'CYP2D6 inhibition'),
       ('CYP3A4', 'CYP3A4 inhibition'), ('P-gp', 'P-gp inhibition')]

def num(s):
    m = re.match(r'\s*([\d.]+)', str(s))
    return float(m.group(1)) if m else np.nan

M = np.array([[num(v) for v in df[c]] for c, _ in eps])

def ep_short(s):
    return s.split(' ↔')[0].strip()
xv['ep'] = xv['端点对(内部↔ADMETlab2)'].map(ep_short)
agree = xv.groupby('ep')['方向一致'].apply(lambda s: (s == '✓').mean()).to_dict()
xv_eps = ['Ames', 'hERG', 'DILI', 'Hepatotoxicity', 'CYP1A2', 'CYP2C9', 'CYP2C19', 'CYP2D6', 'CYP3A4', 'P-gp']

fig = plt.figure(figsize=(13.6, 9.3))
gs = fig.add_gridspec(1, 2, width_ratios=[2.85, 1.15], wspace=0.42,
                      left=0.105, right=0.985, top=0.90, bottom=0.13)
ax = fig.add_subplot(gs[0])
im = ax.imshow(M, cmap='YlOrRd', vmin=0, vmax=1, aspect='auto')
ax.set_xticks(range(16))
ax.set_xticklabels(order, rotation=42, ha='right', fontsize=8.6)
ax.set_yticks(range(len(eps)))
ax.set_yticklabels([c for _, c in eps], fontsize=8.6)
for i in range(len(eps) + 1):
    ax.axhline(i - 0.5, color='white', lw=0.8)
for j in range(17):
    ax.axvline(j - 0.5, color='white', lw=0.8)
for j in (7, 12):
    ax.axvline(j - 0.5, color='black', lw=2)
# herb group labels above columns
for x0, x1, lab in ((0, 6, 'Scutellaria baicalensis'), (7, 11, 'Coptis chinensis'), (12, 15, 'Rheum palmatum')):
    ax.text((x0 + x1) / 2, -1.55, lab, ha='center', va='bottom', fontsize=8.2, style='italic', color='#444444')
# key-signal annotations
ax.add_patch(Rectangle((order.index('Berberine') - 0.5, 2 - 0.5), 1, 1, fill=False, ec='blue', lw=2.5))
ax.text(order.index('Berberine'), 2 - 0.72, '*', ha='center', va='bottom', fontsize=14, color='blue', fontweight='bold')
for c in ['Emodin', 'Rhein', 'Aloe-emodin', 'Chrysophanol']:
    ax.add_patch(Rectangle((order.index(c) - 0.5, -0.5), 1, 1, fill=False, ec='purple', lw=2))
ax.set_ylim(len(eps) - 0.5, -2.3)
ax.set_title('A  In-house rule-model risk probabilities (0–1, descriptive)', fontsize=10.5, loc='left', pad=14)
cb = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.06)
cb.set_label('Predicted probability', fontsize=8.5)
cb.ax.tick_params(labelsize=8)

ax2 = fig.add_subplot(gs[1])
ypos = np.arange(len(xv_eps))[::-1]
vals = [agree.get(k, np.nan) * 100 for k in xv_eps]
colors = ['#2ca02c' if v >= 70 else ('#ff7f0e' if v >= 40 else '#d62728') for v in vals]
ax2.barh(ypos, vals, color=colors, alpha=0.88, height=0.62)
for y, v in zip(ypos, vals):
    ax2.text(v + 2, y, f'{v:.0f}%', va='center', fontsize=8.2)
ax2.set_yticks(ypos)
ax2.set_yticklabels(xv_eps, fontsize=8.6)
ax2.set_xlim(0, 118)
ax2.set_xlabel('Cross-platform directional\nconcordance (%)', fontsize=8.5)
ax2.set_title('B  Concordance with ADMETlab 2.0', fontsize=10.5, loc='left', pad=14)
ax2.axvline(70, color='grey', ls='--', lw=0.8)
ax2.tick_params(labelsize=8.2)
fig.text(0.105, 0.055,
         'Purple boxes: anthraquinone Ames (cross-platform positive, consistent with literature);  '
         'blue box *: berberine hERG (platform conflict; patch-clamp pending).',
         fontsize=7.8, ha='left')
fig.text(0.105, 0.022,
         'Note: micronucleus, nephrotoxicity, mitochondrial, synoviocyte, chondrocyte and macrophage-viability '
         'endpoints are in-house only (not cross-validated); HIA concordance 81% (13/16).',
         fontsize=7.8, ha='left')
fig.savefig(OUT + '/fig04_tox_xval_heatmap_EN.png', dpi=300)
plt.close(fig)
letterbox(OUT + '/fig04_tox_xval_heatmap_EN.png', TARGETS['fig04_tox_xval_heatmap_EN.png'])

# ================= fig20: 5-component radar =================
data = {  # ΔG, -dRMSF, RMSD, D166occ, I165A ddG, |Vina|
    'Baicalein':  [32.82, 0.0008, 0.244, 1.000, 6.46, 9.540],
    'Wogonin':    [31.18, 0.0093, 0.250, 0.816, 0.61, 9.380],
    'Berberine':  [30.57, 0.0117, 0.437, 0.953, 7.36, 9.795],
    'Oroxylin A': [33.48, 0.0279, 0.388, 0.643, 5.98, 9.090],
    'Apigenin':   [30.41, -0.0661, 0.299, 0.994, 0.07, 9.622],
}
axes_lab = ['Binding energy\n|MM/GBSA|', 'Pocket stabilization\n(−ΔRMSF)', 'Trajectory stability\n(1−RMSD/max)',
            'D166 contact\noccupancy', 'I165 anchor\ndependence (ΔΔG)', 'Docking strength\n|Vina|']
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
fig, axr = plt.subplots(figsize=(8.6, 6.8), subplot_kw=dict(polar=True))
fig.subplots_adjust(left=0.16, right=0.74, top=0.86, bottom=0.08)
for name, row, c in zip(data.keys(), norm, colors5):
    v = row.tolist() + row[:1].tolist()
    axr.plot(ang, v, color=c, lw=2, label=name)
    axr.fill(ang, v, color=c, alpha=0.06)
axr.set_xticks(ang[:-1])
axr.set_xticklabels(axes_lab, fontsize=9)
axr.tick_params(axis='x', pad=13)
axr.set_rlabel_position(30)
axr.set_ylim(0, 1.05)
axr.set_yticks([0.25, 0.5, 0.75, 1.0])
axr.set_yticklabels(['0.25', '0.50', '0.75', '1.00'], fontsize=7.5)
axr.legend(loc='upper left', bbox_to_anchor=(1.04, 1.06), fontsize=9, frameon=False)
axr.set_title('Multi-axis computational profile of five IKKβ-stage compounds\n(axes independently normalized; outer = stronger)',
              fontsize=10, pad=26)
fig.savefig(OUT + '/fig20_radar_ikkb_5components_EN.png', dpi=300)
plt.close(fig)
letterbox(OUT + '/fig20_radar_ikkb_5components_EN.png', TARGETS['fig20_radar_ikkb_5components_EN.png'])

# ================= fig24: integrated schematic =================
fig, ax = plt.subplots(figsize=(13.55, 8.0))
ax.set_xlim(0, 135.5); ax.set_ylim(0, 80); ax.axis('off')
fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)

C_BLUE = '#dce9f7'; E_BLUE = '#2c5f8a'
C_RED = '#f5d5d5'; E_RED = '#b03a3a'
C_GREEN = '#e2f0da'; E_GREEN = '#4e8542'
C_ORG = '#fdeeda'; E_ORG = '#d0892c'
C_PINK = '#fbe4e4'; E_PINK = '#c0504d'
C_YEL = '#fdf6dc'; E_YEL = '#a08c3a'

def box(x, y, w, h, text, fc, ec, fs=8.6, bold=False, dashed=False, tc='black'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.4,rounding_size=1.2',
                                fc=fc, ec=ec, lw=1.6, linestyle='--' if dashed else '-', zorder=3))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs,
            color=tc, fontweight='bold' if bold else 'normal', zorder=4)

def arrow(x1, y1, x2, y2, color='#555555', lw=1.8, ls='-', style='-|>'):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=13,
                                 color=color, lw=lw, linestyle=ls, zorder=2))

def blunt(x1, y1, x2, y2, color, lw=2.4):
    ax.plot([x1, x2], [y1, y2], color=color, lw=lw, zorder=2)
    dx, dy = x2 - x1, y2 - y1
    n = (dx**2 + dy**2) ** 0.5
    px, py = -dy / n * 1.6, dx / n * 1.6
    ax.plot([x2 - px, x2 + px], [y2 - py, y2 + py], color=color, lw=lw, zorder=2)

def tag(x, y, text, fc, fs=7.0):
    ax.text(x, y, text, fontsize=fs, ha='center', va='center', zorder=5,
            bbox=dict(boxstyle='round,pad=0.28', fc=fc, ec='none', alpha=0.95))

ax.set_title('Baicalein in HKSA-stimulated macrophages: integration of acquired data\n'
             'and the compound–structure hypothesis', fontsize=12, fontweight='bold', pad=8)

# --- cell boundary (membrane) ---
ax.add_patch(FancyBboxPatch((16, 6), 92, 52, boxstyle='round,pad=0.5,rounding_size=6',
                            fc='#f4f8fc', ec='#8aa8c8', lw=1.8, zorder=1))
ax.text(20, 55.5, 'Cell membrane', fontsize=8, color='#5a7a9a', style='italic', rotation=0)

# --- HKSA stimulus ---
for (cx, cy, r) in ((6.5, 40, 2.2), (10, 43.5, 2.4), (12.5, 39, 2.0), (8.5, 35.5, 2.3), (12, 32.5, 1.9), (5, 33, 1.8)):
    ax.add_patch(Circle((cx, cy), r, fc='#cfcfcf', ec='#888888', lw=1.2, zorder=3))
ax.text(8.5, 47.5, 'HKSA', fontsize=10, fontweight='bold', ha='center')
ax.text(8.5, 45.2, 'heat-inactivated\nS. aureus', fontsize=7.2, ha='center', va='top')
tag(11.5, 28.3, 'S. aureus–associated inflammation model', C_GREEN, fs=6.6)
arrow(14.5, 38, 21.5, 38, lw=2.2)

# --- cascade boxes ---
box(22, 34.5, 15, 7.5, 'TLR2 / NOD\npattern recognition', C_BLUE, E_BLUE, bold=True)
box(41.5, 34.5, 14.5, 7.5, 'MyD88\nIRAK–TRAF6', C_BLUE, E_BLUE, bold=True)
box(60.5, 34.5, 17, 7.5, 'IKKβ kinase complex\n(computational candidate node)', C_RED, E_RED, bold=True)
box(81.5, 34.5, 17.5, 7.5, 'p-IκBα (S32/36)↑\n→ IκBα degradation', C_BLUE, E_BLUE, bold=True)
arrow(37, 38.2, 41.5, 38.2); arrow(56, 38.2, 60.5, 38.2); arrow(77.5, 38.2, 81.5, 38.2)
tag(79.5, 31.6, 'WB: p-IκBα/IκBα +165% (n=3)', '#d7e8f7')

# --- baicalein (green, top) ---
box(46, 62, 26, 12, 'Baicalein (24 μM, WB + ELISA)\nsource: Sanhuang Xiexin Decoction\n(lead formula)',
    C_GREEN, E_GREEN, bold=True)
tag(51, 59.6, 'WB (n=3)', '#d7e8f7'); tag(61.5, 59.6, 'ELISA (n=5)', '#fce6c9')
blunt(59, 61.6, 66.5, 42.4, E_GREEN)

# --- IKK-16 (orange, top right) ---
box(96, 62, 22, 11, 'IKK-16 (0.5 μM)\nhighly selective IKK inhibitor', C_ORG, E_ORG, bold=True)
tag(107, 59.3, 'no additive effect after pathway blockade (non-additive)', '#fce6c9')
blunt(98, 61.6, 76.5, 42.4, E_ORG)

# --- nucleus ---
ax.add_patch(Circle((79, 16.5), 10.5, fc='#f7e8ee', ec='#b06a8a', lw=1.8, zorder=1))
ax.text(79, 18.8, 'NF-κB\ntranscription', fontsize=9.5, fontweight='bold', ha='center', zorder=4)
ax.text(79, 13.4, 'IL1B · IL6 · TNF', fontsize=8, ha='center', zorder=4)
ax.text(66.5, 26.6, 'Nucleus', fontsize=8, color='#8a5a70', style='italic')

# --- p-IκBα -> nucleus, p65 translocation ---
arrow(90, 34.2, 84.5, 26.8, lw=2.0)
box(45, 24.5, 20, 7, 'p65 nuclear translocation\n(nuclear/cytoplasmic gray ratio)', 'white', '#999999',
    fs=8.2, dashed=True)
tag(55, 21.6, 'WB: nuclear/cytoplasmic p65 ratio decreased (exploratory)', '#e8e8e8')
arrow(65, 27, 71.5, 22.5, color='#999999', ls='--')

# --- cytokine output ---
box(110, 24, 23.5, 8.5, 'Supernatant cytokines\nIL-1β · IL-6 · TNF-α ↑', C_PINK, E_PINK, bold=True)
tag(121.8, 21.2, 'ELISA (n=5)', '#fce6c9')
arrow(89.5, 17.5, 109.6, 27.5, lw=2.0)

# --- pocket hypothesis (yellow dashed) ---
box(4, 3.5, 40, 16.5, '', C_YEL, E_YEL, dashed=True)
ax.text(24, 17.2, 'IKKβ catalytic pocket (structural hypothesis, to be validated)',
        fontsize=8.6, fontweight='bold', ha='center')
ax.text(24, 13.6, 'Baicalein / Wogonin / Berberine / Oroxylin A / Apigenin', fontsize=7.6, ha='center')
ax.text(24, 10.9, 'anchoring residues: I165·V29·L21 / V152·M96', fontsize=7.6, ha='center')
ax.text(24, 8.2, 'MM/GBSA −30.4 to −33.5 kcal/mol (model estimate)', fontsize=7.6, ha='center')
ax.text(24, 5.5, 'Ala-scan candidate I165; first mutant MD:\n1 system consistent / 1 inconsistent',
        fontsize=7.2, ha='center', va='top')
arrow(44, 12, 62, 33.9, color=E_YEL, ls='--', lw=1.8)

# --- legend ---
lx, ly = 111, 4.5
ax.add_patch(FancyBboxPatch((lx - 1.5, ly - 1.2), 25.5, 14.5, boxstyle='round,pad=0.3,rounding_size=1',
                            fc='white', ec='#bbbbbb', lw=1.0, zorder=3))
items = [
    ('arrow', '#555555', 'pathway flow (canonical cascade)'),
    ('blunt', E_GREEN, 'inhibition (measured, this study)'),
    ('dash', '#999999', 'computational hypothesis / exploratory'),
]
for k, (kind, col, txt) in enumerate(items):
    yy = ly + 11.2 - k * 3.1
    if kind == 'arrow':
        arrow(lx, yy, lx + 3.6, yy, color=col, lw=1.8)
    elif kind == 'blunt':
        blunt(lx, yy, lx + 3.6, yy, col, lw=2.2)
    else:
        ax.plot([lx, lx + 3.6], [yy, yy], color=col, lw=1.6, ls='--', zorder=4)
    ax.text(lx + 4.6, yy, txt, fontsize=6.8, va='center', zorder=4)
tag(lx + 5.2, ly + 0.9, 'blue = WB (n=3)', '#d7e8f7', fs=6.6)
tag(lx + 16.4, ly + 0.9, 'orange = ELISA (n=5)', '#fce6c9', fs=6.6)

fig.savefig(OUT + '/fig24_schematic_EN.png', dpi=300)
plt.close(fig)
letterbox(OUT + '/fig24_schematic_EN.png', TARGETS['fig24_schematic_EN.png'])
print('ALL DONE')
