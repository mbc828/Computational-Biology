# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""抗菌实验定量：平板菌落计数 + 活/死染色红色通道定量
输出: outputs/antibacterial_counts_20261022.csv
"""
import numpy as np
from PIL import Image
from scipy import ndimage
import os, csv

BASE = r'<WETLAB_DATA_ROOT>\抗菌实验\细菌平板'
LD = r'<WETLAB_DATA_ROOT>\抗菌实验\活死染色'

def count_white_field(path, min_area=8, max_area=2000):
    """暗视野白板照片：亮菌落 on 暗背景。限圆形平皿区域。"""
    im = Image.open(path).convert('L')
    a = np.array(im, dtype=float)
    # 平皿区域：亮度>25 的连通主体取最大圆盘近似
    mask0 = a > 25
    lab, n = ndimage.label(mask0)
    sizes = ndimage.sum(mask0, lab, range(1, n + 1))
    plate = lab == (1 + int(np.argmax(sizes)))
    plate = ndimage.binary_fill_holes(plate)
    # 菌落：平皿内显著亮点
    thr = np.percentile(a[plate], 75) + 1.2 * a[plate].std()
    m = (a > max(thr, 90)) & plate
    m = ndimage.binary_opening(m, iterations=1)
    lab2, n2 = ndimage.label(m)
    sz = ndimage.sum(m, lab2, range(1, n2 + 1))
    keep = (sz >= min_area) & (sz <= max_area)
    return int(keep.sum()), a.shape

def count_phone(path):
    """手机照片：手动裁平皿区域后数亮点。"""
    im = Image.open(path).convert('L')
    a = np.array(im, dtype=float)
    h, w = a.shape
    # 平皿位于画面中部偏下（约 x 25-58%, y 38-78%，按观察）
    crop = a[int(0.38*h):int(0.78*h), int(0.25*w):int(0.58*w)]
    thr = crop.mean() + 2.0 * crop.std()
    m = crop > max(thr, 120)
    m = ndimage.binary_opening(m, iterations=2)
    lab, n = ndimage.label(m)
    sz = ndimage.sum(m, lab, range(1, n + 1))
    keep = (sz >= 20) & (sz <= 5000)
    return int(keep.sum()), a.shape

def red_fraction(path):
    """活死染色红通道：红点面积占比(%)与点数。"""
    im = Image.open(path).convert('RGB')
    a = np.array(im, dtype=float)
    r = a[..., 0]
    m = r > 80
    frac = m.mean() * 100
    lab, n = ndimage.label(m)
    sz = ndimage.sum(m, lab, range(1, n + 1))
    spots = int((sz >= 6).sum())
    return frac, spots

rows = []
print('== 平板菌落计数（CFU/皿）==')
for grp, files, fn in [
    ('Control', ['Control/2_260805_230216_White.jpg', 'Control/2_260805_233820_White.jpg', 'Control/2_260805_234034_White.jpg'], count_white_field),
    ('Baicalein', ['中药黄芩素/1_260805_234744_White.jpg', '中药黄芩素/1_260805_235825_White.jpg', '中药黄芩素/2_260805_235952_White.jpg'], count_white_field),
    ('SHXD formula 125', ['中药整方/中药整方2.jpg', '中药整方/中药整方3.jpg'], count_phone),
]:
    for f in files:
        c, shape = fn(os.path.join(BASE, f))
        rows.append(('plate', grp, os.path.basename(f), c))
        print(f'{grp:20s} {os.path.basename(f):32s} CFU={c}')

print('\n== 活/死染色（红通道=死菌）==')
for grp, files in [
    ('Control', ['control 1_ch01.tif', 'control 2_ch01.tif']),
    ('Baicalein', ['中药1_ch01.tif', '中药2_ch01.tif']),
    ('SHXD formula 125', ['整方1.jpg', '整方2.jpg', '整方3.jpg']),
]:
    for f in files:
        frac, spots = red_fraction(os.path.join(LD, f))
        rows.append(('live_dead_red', grp, f, f'{frac:.2f}% area, {spots} spots'))
        print(f'{grp:20s} {f:22s} red_area={frac:.2f}% spots={spots}')

out = r'<PROJECT_ROOT>\outputs\antibacterial_counts_20261022.csv'
with open(out, 'w', newline='', encoding='utf-8-sig') as fh:
    w = csv.writer(fh)
    w.writerow(['assay', 'group', 'file', 'result'])
    w.writerows(rows)
print('\nsaved', out)
