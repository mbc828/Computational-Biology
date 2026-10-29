# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""Prepare two new figures for v43 NP merge:
1) F7 mechanism schematic — crop of img05 panel B (Chinese title bands removed)
2) F12 real docking figure — regenerated in English from docking_master_30.csv
   with audited redock RMSD values (v42 canonical: IKKb 4.46, NLRP3 5.78, MD2 8.47,
   IL-17A 8.56, CASP1 6.40; all FAIL the 2.0 A criterion)
"""
import csv, os
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = r"<PROJECT_ROOT>"
IMGDIR = os.path.join(ROOT, "outputs", "_v7_imgs")
FIGDIR = os.path.join(ROOT, "outputs", "manuscript_v21_github_release", "figures")
os.makedirs(FIGDIR, exist_ok=True)

# ---------- 1) mechanism schematic crop ----------
src = Image.open(os.path.join(IMGDIR, "img05_rId14.png"))
# original 4002x1763; panel B diagram box determined by region inspection
crop = src.crop((2050, 545, 3965, 1545))
crop.save(os.path.join(FIGDIR, "F7_proposed_mechanism_v43.png"), optimize=True)
print("mechanism crop:", crop.size)

# ---------- 2) docking figure ----------
rows = list(csv.DictReader(open(os.path.join(ROOT, "SHXD_MD", "docking_local", "docking_master_30.csv"), encoding="utf-8-sig")))
targets = ["IKKb", "MD2", "NLRP3", "IL17A", "CASP1"]
tlabels = ["IKK\u03b2", "MD2", "NLRP3", "IL-17A", "CASP1"]
comps = ["Emodin", "Berberine", "Apigenin", "Baicalein", "Wogonin", "OroxylinA"]
clabels = ["Emodin", "Berberine", "Apigenin", "Baicalein", "Wogonin", "Oroxylin A"]
M = np.full((len(comps), len(targets)), np.nan)
for r in rows:
    if r["ligand"] in comps and r["target"] in targets:
        M[comps.index(r["ligand"]), targets.index(r["target"])] = float(r["vina_score_kcal_mol"])

# audited rank-1 redock RMSD (v42 canonical text + redock logs)
redock = [("IKK\u03b2 (4KIK)", 4.46), ("NLRP3 (7PZC)", 5.78), ("CASP1 (6PZP)", 6.40),
          ("TLR4\u2013MD2 (3FXI)", 8.47), ("IL-17A (5H4I)", 8.56)]

fig = plt.figure(figsize=(8.6, 3.6), dpi=300)
gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1.0], wspace=0.52, left=0.125, right=0.99, top=0.85, bottom=0.16)

# Panel A: heatmap
axA = fig.add_subplot(gs[0])
vmax, vmin = -6.0, -10.0
im = axA.imshow(M, cmap="RdYlGn_r", vmin=vmin, vmax=vmax, aspect="auto")
axA.set_xticks(range(len(targets))); axA.set_xticklabels(tlabels, fontsize=8)
axA.set_yticks(range(len(comps))); axA.set_yticklabels(clabels, fontsize=8)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        if not np.isnan(v):
            axA.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=6.5,
                     color="black")
axA.set_title("A  Vina score matrix (6\u00d75)", fontsize=9.5, fontweight="bold", loc="left")
cb = fig.colorbar(im, ax=axA, fraction=0.028, pad=0.02)
cb.set_label("Vina score (kcal/mol)", fontsize=6.5); cb.ax.tick_params(labelsize=6)
axA.set_xticks(np.arange(-0.5, len(targets), 1), minor=True)
axA.set_yticks(np.arange(-0.5, len(comps), 1), minor=True)
axA.grid(which="minor", color="white", linewidth=0.8)
axA.tick_params(which="minor", length=0)

# Panel B: redock RMSD bars
axB = fig.add_subplot(gs[1])
names = [x[0] for x in redock]; vals = [x[1] for x in redock]
colors = ["#c44e52"] * len(vals)
ypos = np.arange(len(vals))[::-1]
axB.barh(ypos, vals, color=colors, alpha=0.85, height=0.62)
for y, v in zip(ypos, vals):
    axB.text(v + 0.12, y, "%.2f" % v, va="center", fontsize=7.5)
axB.axvline(2.0, color="black", ls="--", lw=1.0)
axB.set_yticks(ypos); axB.set_yticklabels(names, fontsize=8)
axB.set_xlabel("Rank-1 redock RMSD (\u00c5)", fontsize=8)
axB.set_xlim(0, 10.2)
axB.set_ylim(-0.7, 4.7)
axB.set_title("B  Redocking benchmark", fontsize=9.5, fontweight="bold", loc="left")
axB.text(9.9, -0.52, "dashed line: 2.0 \u00c5 criterion \u2014 none passed",
         fontsize=7, ha="right", color="#c44e52")
axB.tick_params(labelsize=7.5)
for s in ["top", "right"]:
    axB.spines[s].set_visible(False)

fig.savefig(os.path.join(FIGDIR, "F12_docking_real_v43.png"), bbox_inches="tight", pad_inches=0.08)
fig.savefig(os.path.join(FIGDIR, "F12_docking_real_v43.pdf"), bbox_inches="tight", pad_inches=0.08)
print("docking figure saved")
