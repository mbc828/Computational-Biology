# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""拼装 6 成分 x 4 联复合大图：FES2D | FES3D | 复合物结构 | 相互作用指纹"""
import os
from PIL import Image, ImageDraw, ImageFont

OUT = os.path.dirname(os.path.abspath(__file__))
FES = os.path.join(OUT, 'fes_panels')
ROOT = r'<MD_DATA_ROOT>\fes_ikkb_20261022'
SKETCH = FES  # interaction_sketch.py 与 fes 面板同目录输出

LIGS = ['Baicalein', 'Wogonin', 'Berberine', 'OroxylinA', 'Apigenin', 'Emodin']
COLS = ['2D free-energy landscape', '3D energy surface',
        'IKKβ–ligand complex', 'Residue interaction map']

CELL_W, CELL_H = 1100, 900
PAD_L, PAD_T, GAP = 340, 120, 18

def cell_paths(lig):
    return [
        os.path.join(FES, f'fes2d_{lig}.png'),
        os.path.join(FES, f'fes3d_{lig}.png'),
        os.path.join(ROOT, f'render_{lig}.png'),
        os.path.join(SKETCH, f'interact_{lig}.png'),
    ]

def fit(img, w, h):
    iw, ih = img.size
    s = min(w / iw, h / ih)
    return img.resize((max(1, int(iw * s)), max(1, int(ih * s))), Image.LANCZOS)

def font(sz):
    for p in [r'<LOCAL_DRIVE>/Windows\Fonts\arialbd.ttf', r'<LOCAL_DRIVE>/Windows\Fonts\arial.ttf']:
        if os.path.exists(p):
            return ImageFont.truetype(p, sz)
    return ImageFont.load_default()

W = PAD_L + CELL_W * 4 + GAP * 3
H = PAD_T + CELL_H * 6 + GAP * 5
canvas = Image.new('RGB', (W, H), 'white')
dr = ImageDraw.Draw(canvas)
f_col, f_row, f_lab = font(44), font(48), font(64)

missing = []
for j, ct in enumerate(COLS):
    x = PAD_L + j * (CELL_W + GAP) + CELL_W // 2
    bb = dr.textbbox((0, 0), ct, font=f_col)
    dr.text((x - (bb[2] - bb[0]) // 2, 30), ct, fill='black', font=f_col)

for i, lig in enumerate(LIGS):
    y = PAD_T + i * (CELL_H + GAP) + CELL_H // 2
    bb = dr.textbbox((0, 0), lig, font=f_row)
    dr.text((PAD_L - 40 - (bb[2] - bb[0]), y - (bb[3] - bb[1]) // 2), lig,
            fill='black', font=f_row)
    for j, p in enumerate(cell_paths(lig)):
        x0 = PAD_L + j * (CELL_W + GAP)
        y0 = PAD_T + i * (CELL_H + GAP)
        if not os.path.exists(p):
            missing.append(p)
            dr.rectangle([x0, y0, x0 + CELL_W, y0 + CELL_H], outline='#cccccc')
            continue
        img = fit(Image.open(p).convert('RGB'), CELL_W, CELL_H)
        ox, oy = x0 + (CELL_W - img.size[0]) // 2, y0 + (CELL_H - img.size[1]) // 2
        canvas.paste(img, (ox, oy))
        dr.text((x0 + 10, y0 + 6), chr(ord('A') + i * 4 + j), fill='black', font=f_lab)

out = os.path.join(OUT, 'fig_fes_composite_6components.png')
canvas.save(out, dpi=(300, 300))
print('saved:', out, canvas.size)
if missing:
    print('MISSING:')
    for m in missing:
        print(' ', m)
