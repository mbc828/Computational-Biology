# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""图13 v41版：整合模型示意——实线元素=本研究已获得数据，虚线元素=计算假说/待验证。
   核/质 p65 标探索性；删去"遮蔽成立""↓78%""锚定Ki"等过强表述。"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyBboxPatch, FancyArrowPatch, Circle
import numpy as np

plt.rcParams['font.family'] = ['Microsoft YaHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

fig, ax = plt.subplots(figsize=(16, 9.45), dpi=200)
ax.set_xlim(0, 16); ax.set_ylim(0, 9.45); ax.axis('off')
fig.patch.set_facecolor('white')

C_CYT = '#eef4fb'; C_MEM = '#7a9cc6'; C_NUC = '#f6e8f0'; C_NUCMEM = '#b57a9c'
C_BOX = '#dce9f7'; C_BOXE = '#4a6fa5'; C_IKK = '#f8d7da'; C_IKKE = '#c8363b'
C_SHXD = '#d9ead3'; C_SHXDE = '#2e7d32'; C_IKK16 = '#fde8cf'; C_IKK16E = '#e07b00'

def box(x, y, w, h, text, fc=C_BOX, ec=C_BOXE, fs=11, bold=False, tc='black', lw=1.5, r=0.12, ls='-'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0.06,rounding_size={r}',
                                fc=fc, ec=ec, lw=lw, ls=ls, zorder=5))
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fs,
            color=tc, fontweight='bold' if bold else 'normal', zorder=6, linespacing=1.35)

def arrow(x1, y1, x2, y2, color='#444444', lw=2.2, style='-|>', ms=16, ls='-', z=4):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=ms,
                                 color=color, lw=lw, linestyle=ls, zorder=z))

def blunt(x1, y1, x2, y2, color, lw=2.6):
    ax.plot([x1, x2], [y1, y2], color=color, lw=lw, zorder=7, solid_capstyle='round')
    ax.plot([x2], [y2], marker='_', ms=22, color=color, mew=3.2, zorder=8)

# ---- 细胞体 ----
ax.add_patch(Ellipse((7.1, 4.4), 11.4, 7.6, fc=C_CYT, ec=C_MEM, lw=3, zorder=1))
ax.add_patch(Ellipse((7.1, 4.4), 11.9, 8.1, fc='none', ec=C_MEM, lw=1.2, alpha=0.45, zorder=1))
ax.text(2.0, 7.9, '细胞膜', fontsize=10, color=C_MEM, style='italic')
ax.add_patch(Circle((10.3, 3.3), 1.75, fc=C_NUC, ec=C_NUCMEM, lw=2.5, zorder=2))
ax.text(10.3, 5.35, '细胞核', fontsize=10, color=C_NUCMEM, style='italic', ha='center')

# ---- 胞外：HKSA ----
for dx, dy in [(-0.35, 0.3), (0.35, 0.32), (0, -0.05), (-0.3, -0.42), (0.38, -0.4), (0.02, -0.72)]:
    ax.add_patch(Circle((0.85 + dx, 8.0 + dy), 0.16, fc='#c9c9c9', ec='#777777', lw=1, zorder=3))
ax.text(0.85, 8.85, 'HKSA\n热灭活金葡菌', ha='center', fontsize=10.5, color='#555555', linespacing=1.3)
ax.text(0.85, 6.6, '金葡菌相关炎症模型', ha='center', fontsize=8.5, color='white',
        bbox=dict(boxstyle='round,pad=0.25', fc='#5a8a5a', ec='none'))

# ---- 膜受体与胞质级联（经典通路背景） ----
box(2.35, 6.7, 1.85, 0.95, 'TLR2 / NOD\n模式识别', fs=10.5, bold=True)
arrow(1.35, 7.6, 2.35, 7.2)
box(4.75, 6.7, 2.0, 0.95, 'MyD88\nIRAK·TRAF6', fs=10.5)
arrow(4.2, 7.18, 4.75, 7.18)
box(7.35, 6.55, 2.5, 1.25, 'IKKβ 激酶复合物\n（计算候选节点）', fc=C_IKK, ec=C_IKKE, fs=12, bold=True, lw=2.6)
arrow(6.75, 7.18, 7.35, 7.18)
box(10.45, 6.7, 2.35, 0.95, 'p-IκBα (S32/36)↑\n→ IκBα 降解', fs=10)
arrow(9.85, 7.18, 10.45, 7.18)
ax.text(11.6, 6.35, 'WB：p-IκBα/IκBα 增量↓65%（n=3）', fontsize=8.5, color='white', ha='center',
        bbox=dict(boxstyle='round,pad=0.25', fc='#3a6ea5', ec='none'))

# p65 核转位（探索性：核/质灰度直比，未分区归一化）
box(8.0, 4.35, 2.2, 1.0, 'p65 核转位\n（核/质灰度直比）', fs=10.5, ls='--', ec='#888888', fc='#eeeeee')
arrow(11.6, 6.7, 9.6, 5.35)
arrow(10.1, 4.85, 10.1, 4.35, ls='--', color='#888888')
ax.text(5.9, 3.9, 'WB：核/质 p65 直比下降（探索性）', fontsize=8.5, color='white', ha='center',
        bbox=dict(boxstyle='round,pad=0.25', fc='#888888', ec='none'))

# 核内转录
ax.text(10.3, 3.55, 'NF-κB 转录', ha='center', fontsize=11.5, fontweight='bold', color='#8a2e44', zorder=6)
ax.text(10.3, 3.0, 'IL1B · IL6 · TNF', ha='center', fontsize=10, color='#8a2e44', zorder=6)

# 细胞因子释放
box(12.6, 2.5, 2.9, 1.1, '上清细胞因子\nIL-1β·IL-6·TNF-α ↑', fs=10.5, fc='#f3e2e2', ec='#a94442')
arrow(11.3, 2.8, 12.6, 2.95)
ax.text(14.05, 2.05, 'ELISA：增量抑制 45–65%（n=5）', fontsize=8.5, color='white', ha='center',
        bbox=dict(boxstyle='round,pad=0.25', fc='#b0762a', ec='none'))

# ---- 干预：黄芩素（WB）＋整方（ELISA/抗菌）+ IKK-16（本研究测量） ----
box(4.35, 8.25, 3.4, 1.05, '黄芩素 Baicalein（24 μM，WB·ELISA）\n来源：三黄泻心汤（先导方剂）', fc=C_SHXD, ec=C_SHXDE, fs=9.5, bold=True)
blunt(6.6, 8.25, 7.7, 7.85, C_SHXDE)
ax.text(4.55, 7.95, 'WB（n=3）·ELISA（n=5）＝黄芩素（本研究测量）', fontsize=8, color='white',
        bbox=dict(boxstyle='round,pad=0.22', fc='#2e7d32', ec='none'))

box(11.3, 8.35, 2.6, 0.95, 'IKK-16（0.5 μM）\n高选择性 IKK 抑制剂', fc=C_IKK16, ec=C_IKK16E, fs=10.5, bold=True)
blunt(12.0, 8.35, 9.3, 7.85, C_IKK16E)
ax.text(12.75, 8.0, '预阻断后未检测到黄芩素额外效应（非加和）', fontsize=8, color='white',
        bbox=dict(boxstyle='round,pad=0.22', fc='#e07b00', ec='none'))

# ---- 口袋插图（左下，结构假说，虚线） ----
ax.add_patch(FancyBboxPatch((0.45, 0.45), 5.5, 2.6, boxstyle='round,pad=0.08,rounding_size=0.15',
                            fc='#faf7ee', ec='#b8a25e', lw=1.8, ls='--', zorder=4))
ax.text(3.2, 2.78, 'IKKβ 催化口袋（结构假说，待验证）', ha='center', fontsize=10.5,
        fontweight='bold', color='#8a6d1f', zorder=6)
t = np.linspace(0, np.pi, 100)
px = 1.35 + 1.05*np.cos(t); py = 0.75 + 1.3*np.sin(t)
ax.plot(px, py, color='#8a6d1f', lw=3, zorder=5)
ax.text(1.35, 1.15, 'I165·V29·L21\nV152·M96', ha='center', fontsize=8.5, color='#8a6d1f', zorder=6, linespacing=1.3)
ax.text(3.6, 1.95, '黄芩素·汉黄芩素·小檗碱\n千层纸素A·芹菜素', fontsize=9.5, color='#333333', zorder=6, linespacing=1.4)
ax.text(3.6, 1.0, 'MM/GBSA −30.4~−33.5 kcal/mol（模型估计）\nAla-scan 候选 I165；首批突变体MD\n1 体系一致 / 1 体系不一致', fontsize=8.5,
        color='#555555', zorder=6, linespacing=1.45)
arrow(6.0, 1.85, 7.5, 6.5, color='#b8a25e', lw=1.6, ls='--', ms=12)

# ---- 图例 ----
lx, ly = 12.45, 0.35
ax.add_patch(FancyBboxPatch((lx-0.15, ly), 3.45, 1.95, boxstyle='round,pad=0.06', fc='white', ec='#bbbbbb', lw=1, zorder=8))
ax.plot([lx, lx+0.45], [ly+1.68, ly+1.68], color='#444444', lw=2.2, zorder=9)
ax.text(lx+0.6, ly+1.68, '通路流向（经典级联背景）', fontsize=9, va='center', zorder=9)
ax.plot([lx, lx+0.45], [ly+1.31, ly+1.31], color=C_SHXDE, lw=2.6, zorder=9)
ax.plot([lx+0.45], [ly+1.31], marker='_', ms=16, color=C_SHXDE, mew=3, zorder=9)
ax.text(lx+0.6, ly+1.31, '抑制（本研究测量）', fontsize=9, va='center', zorder=9)
ax.plot([lx, lx+0.45], [ly+0.94, ly+0.94], color='#b8a25e', lw=1.6, ls='--', zorder=9)
ax.text(lx+0.6, ly+0.94, '计算假说 / 待验证环节', fontsize=9, va='center', zorder=9)
ax.add_patch(FancyBboxPatch((lx+0.02, ly+0.50), 0.42, 0.26, boxstyle='round,pad=0.03', fc='white', ec='#888888', lw=1.4, ls='--', zorder=9))
ax.text(lx+0.6, ly+0.63, '探索性（未完成归一化）', fontsize=9, va='center', zorder=9)
ax.text(lx+0.05, ly+0.18, '蓝标=WB(n=3)  橙标=ELISA(n=5)', fontsize=8.5, va='center', color='#555555', zorder=9)

ax.set_title('黄芩素对 HKSA 诱导巨噬细胞炎症的影响：已获得数据与成分结构假说的整合示意',
             fontsize=15, fontweight='bold', pad=14)
fig.tight_layout()
out = r'<PROJECT_ROOT>\outputs\manuscript_v21_github_release\figures\F13_integrated_model_v41.png'
fig.savefig(out, dpi=200, facecolor='white')
print('saved', out)
