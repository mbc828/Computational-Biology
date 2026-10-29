# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""make_fig_wb.py — WB 六面板组图（Prism 风格）
A: p-IκBα/IκBα/p-p65/p65/Actin 条带; B: p-IκBα/IκBα 比率; C: p-p65/p65 比率
D: Nuclear/Cytoplasmic p65/Actin 条带; E: 核 p65 占比; F: 质 p65 占比
"""
import json, os
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

BASE = os.path.join(r"<PROJECT_ROOT>", "outputs", "wet_figures_20261014")
PPTX = os.path.join(r"<PROJECT_ROOT>", "outputs", "pptx_strips")
OUT = os.path.join(BASE, "figures")
os.makedirs(OUT, exist_ok=True)
stats = json.load(open(os.path.join(r"<PROJECT_ROOT>", "outputs", "manuscript_v21_github_release", "04_data", "wb_audit", "wb_stats_v92.json"), encoding="utf-8"))

# PPT 裁好的条带（泳道顺序 Control/PJI+SHXD/PJI，rep1 列）
STRIP_MAP = {
    "p-IkBa":   (PPTX, "022", "image1.png"),
    "IkBa":     (PPTX, "022", "image4.png"),
    "p-p65":    (PPTX, "044", "image1.png"),
    "p65":      (PPTX, "044", "image4.png"),
    "Nuc-p65":  (PPTX, "022", "image7.png"),
    "Cyto-p65": (PPTX, "022", "image10.png"),
    "Actin":    (PPTX, "022", "image13.png"),
}

rcParams["font.family"] = "Arial"
rcParams["axes.linewidth"] = 1.2

GROUPS = ["Control", "PJI+SHXD", "PJI"]   # stats.json 键名（内部）
XLAB = ["Control", "HKSA\n+Baicalein", "HKSA"]  # 图中显示标签（WB 干预为黄芩素单体）
BAR_FC = "#dbe9f6"  # 浅蓝
BAR_HATCH = ".."

def load_strip(name):
    """读取 PPT 裁好的条带（已转正），保持原始宽高比。"""
    base, sub, fn = STRIP_MAP[name]
    return Image.open(os.path.join(base, sub, "ppt", "media", fn)).convert("L")

def strip_panel(fig, rect, names, labels, group_labels=True):
    """在 fig 的 rect 区域内从顶向下紧凑堆叠条带（按原始宽高比，不拉伸）。"""
    ax = fig.add_axes(rect)
    ax.axis("off")
    ax_w_in = rect[2]*13.6
    ax_h_in = rect[3]*8.6
    x0, x1 = 0.13, 0.995             # 条带拉长, 保持宽高比
    disp_w_in = (x1-x0)*ax_w_in
    gap_in = 0.025                     # 条带间小间隙
    y_top = 1.0 - 0.005                # 从面板顶部开始
    for nm, lb in zip(names, labels):
        im = load_strip(nm)
        hf = (disp_w_in*im.height/im.width)/ax_h_in
        yc = y_top - hf/2
        ax.imshow(im, extent=[x0, x1, yc-hf/2, yc+hf/2], aspect="auto",
                  cmap="gray", vmin=0, vmax=255, transform=ax.transAxes)
        ax.text(0.125, yc, lb, transform=ax.transAxes, ha="right",
                va="center", fontsize=12, fontweight="normal", fontstyle="normal")
        y_top = yc - hf/2 - gap_in/ax_h_in
    if group_labels:
        # 水平两行标签，整体放在条带堆下方并留足间隙，不压图
        for k, g in enumerate(XLAB):
            x = x0 + (k+0.5)*((x1-x0)/3)
            ax.text(x, y_top - 0.035, g, transform=ax.transAxes, ha="center",
                    va="top", fontsize=11.5, linespacing=1.1)

def bar_panel(fig, rect, metric, ylab, ymax=None):
    ax = fig.add_axes(rect)
    blk = stats["WB"][metric]
    means = [blk["groups"][g]["mean"] for g in GROUPS]
    sds = [blk["groups"][g]["sd"] for g in GROUPS]
    x = np.arange(3)
    ax.bar(x, means, width=0.62, color=BAR_FC, edgecolor="black",
           hatch=BAR_HATCH, linewidth=1.2, zorder=3)
    ax.errorbar(x, means, yerr=sds, fmt="none", ecolor="black",
                elinewidth=1.2, capsize=4, zorder=4)
    ax.set_xticks(x); ax.set_xticklabels(XLAB, fontsize=11.5, rotation=38,
                                         rotation_mode="anchor", ha="right")
    ax.set_ylabel(ylab, fontsize=12.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=1.2, labelsize=11)
    mx = max(m+s for m, s in zip(means, sds))
    ymax = ymax or mx*1.55
    ax.set_ylim(0, ymax)
    # 显著性括号：三组两两 (0-1, 1-2, 0-2)
    pw = blk["pairwise"]
    def star(p):
        return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
    pairs = [(0, 1, pw["Control vs PJI+SHXD"]["p"]),
             (1, 2, pw["PJI+SHXD vs PJI"]["p"]),
             (0, 2, pw["Control vs PJI"]["p"])]
    hs = [mx*1.18, mx*1.30, mx*1.46]
    for (a, b, p), hh in zip(pairs, hs):
        ax.plot([a, a, b, b], [hh-mx*0.03, hh, hh, hh-mx*0.03], lw=1.1, c="black")
        ax.text((a+b)/2, hh+mx*0.005, star(p), ha="center", va="bottom",
                fontsize=11.5, fontweight="bold")

fig = plt.figure(figsize=(13.6, 8.6))
# 上排: A 条带(总) / B / C；下排: D 条带(核质) / E / F
# 柱状图高度 = 同行条带堆叠高度 (行1: 0.32, 行2: 0.19), 顶边逐行对齐
strip_panel(fig, [0.075, 0.58, 0.24, 0.32],
            ["p-IkBa", "IkBa", "p-p65", "p65", "Actin"],
            ["p-IκBα", "IκBα", "p-p65", "p65", "Actin"])
bar_panel(fig, [0.40, 0.58, 0.20, 0.32], "p-IkBa/IkBa", "p-IκBα/IκBα ratio", ymax=1.55)
bar_panel(fig, [0.71, 0.58, 0.20, 0.32], "p-p65/p65", "p-p65/p65 ratio", ymax=1.25)
strip_panel(fig, [0.075, 0.245, 0.24, 0.19],
            ["Nuc-p65", "Cyto-p65", "Actin"],
            ["Nuclear p65", "Cytoplasmic p65", "Actin"])
bar_panel(fig, [0.40, 0.245, 0.20, 0.19], "Nuclear p65 fraction",
          "Nuclear p65 level", ymax=1.25)
bar_panel(fig, [0.71, 0.245, 0.20, 0.19], "Cytoplasmic p65 fraction",
          "Cytoplasmic p65 level", ymax=1.25)

# 面板序号：统一用图坐标，同一行同一高度；B/E 放在条带图与柱状图间隙正中，C/F 同理
TAG_FS = 20
ROW1_TAG_Y = 0.928   # 行1 序号统一高度（axes 顶 0.90 上方）
ROW2_TAG_Y = 0.465   # 行2 序号统一高度（axes 顶 0.435 上方；与行1标签留足间距）
tags = [("A", 0.030, ROW1_TAG_Y, "left"),
        ("B", 0.357, ROW1_TAG_Y, "center"),
        ("C", 0.655, ROW1_TAG_Y, "center"),
        ("D", 0.030, ROW2_TAG_Y, "left"),
        ("E", 0.357, ROW2_TAG_Y, "center"),
        ("F", 0.655, ROW2_TAG_Y, "center")]
for tag, tx, ty, ha in tags:
    fig.text(tx, ty, tag, fontsize=TAG_FS, fontweight="bold", va="top", ha=ha)

out = os.path.join(OUT, "Fig_WB_NFkB_6panel.png")
fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.12)
print("saved:", out)
