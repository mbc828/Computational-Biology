# Public source snapshot: replace path placeholders before running.
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(r'<PROJECT_ROOT>')
V82 = ROOT / 'outputs/manuscript_v21_github_release/figure_consolidation_v82'
MEDIA = V82 / 'work/unpacked/word/media'
OUT = ROOT / 'outputs/manuscript_v21_github_release/figure_consolidation_v84/composite_figures'
OUT.mkdir(parents=True, exist_ok=True)
FONT = r'<LOCAL_DRIVE>/Windows\Fonts\arial.ttf'
FONT_B = r'<LOCAL_DRIVE>/Windows\Fonts\arialbd.ttf'


def im(n):
    return Image.open(MEDIA / f'image{n}.png').convert('RGB')


def erase(img, boxes):
    draw = ImageDraw.Draw(img)
    for box in boxes:
        draw.rectangle(box, fill='white')
    return img


def relabel(img, replacements, size=70):
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype(FONT_B, size)
    for x, y, w, h, letter in replacements:
        draw.rectangle((x - 8, y - 5, x + w + 8, y + h + 5), fill='white')
        if letter:
            draw.text((x, y), letter, font=font, fill='#20252b')
    return img


def with_header(img, text, height=150, size=62, x=130):
    panel = Image.new('RGB', (img.width, img.height + height), 'white')
    panel.paste(img, (0, height))
    ImageDraw.Draw(panel).text((x, 25), text, font=ImageFont.truetype(FONT_B, size), fill='#20252b')
    return panel


def canvas(rows, width=3600, margin=40, gap=45):
    placed = []
    y = margin
    for images in rows:
        cell_w = (width - 2 * margin - gap * (len(images) - 1)) // len(images)
        scaled = []
        for image in images:
            h = round(image.height * cell_w / image.width)
            scaled.append(image.resize((cell_w, h), Image.Resampling.LANCZOS))
        row_h = max(image.height for image in scaled)
        placed.append((y, cell_w, scaled, row_h))
        y += row_h + gap
    sheet = Image.new('RGB', (width, y - gap + margin), 'white')
    for y0, cell_w, images, row_h in placed:
        for j, image in enumerate(images):
            x = margin + j * (cell_w + gap)
            sheet.paste(image, (x, y0))
    return sheet


def label_sheet(sheet, labels, size=75):
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.truetype(FONT_B, size)
    for x, y, letter in labels:
        draw.text((x, y), letter, font=font, fill='#20252b')
    return sheet


def save(sheet, name):
    path = OUT / f'{name}.png'
    sheet.save(path, optimize=True, dpi=(300, 300))
    preview = sheet.copy()
    preview.thumbnail((1800, 1800), Image.Resampling.LANCZOS)
    preview.save(OUT / f'{name}_preview.jpg', quality=92, optimize=True)
    print(name, sheet.size, path.stat().st_size, 'preview', preview.size)


# Figure 2: mechanism/workflow on the left; docking and endpoint energy stacked on the right.
a1 = im(1)
dock = im(13)
b13 = with_header(dock.crop((0, 110, dock.width // 2, dock.height)),
                  'Vina score matrix (6×5)', 140, 68)
ikk = im(14).crop((0, 700, 4251, 2100))
erase(ikk, [(450, 1270, 1300, 1400)])
xaxis = im(14).crop((500, 3520, 4251, 3800))
xd = ImageDraw.Draw(xaxis)
xd.line((200, 0, 3500, 0), fill='#20252b', width=3)
for xx in (460, 1080, 1670, 2280, 2890, 3500):
    xd.line((xx, 0, xx, 20), fill='#20252b', width=3)
ikk_with_axis = Image.new('RGB', (4251, ikk.height + xaxis.height), 'white')
ikk_with_axis.paste(ikk, (0, 0))
ikk_with_axis.paste(xaxis, (500, ikk.height))
c14 = with_header(ikk_with_axis, 'IKKβ endpoint energy estimates', 180, 140, 420)

width, margin, gap, left_w = 3600, 40, 45, 2050
right_w = width - 2 * margin - gap - left_w
a = a1.resize((left_w, round(a1.height * left_w / a1.width)), Image.Resampling.LANCZOS)
b = b13.resize((right_w, round(b13.height * right_w / b13.width)), Image.Resampling.LANCZOS)
c = c14.resize((right_w, round(c14.height * right_w / c14.width)), Image.Resampling.LANCZOS)
right_h = b.height + gap + c.height
sheet2 = Image.new('RGB', (width, max(a.height, right_h) + 2 * margin), 'white')
sheet2.paste(a, (margin, margin))
right_x = margin + left_w + gap
sheet2.paste(b, (right_x, margin))
c_y = margin + b.height + gap
sheet2.paste(c, (right_x, c_y))
label_sheet(sheet2, [(55, 55, 'A'), (right_x + 15, 55, 'B'), (right_x + 15, c_y + 15, 'C')], 75)
save(sheet2, 'Figure_2_virtual_screening')

# Figure 3: four baicalein panels on top, three cross-compound fingerprints below.
fes_dir = ROOT / 'outputs/fes_panels'
fes2d = Image.open(fes_dir / 'fes2d_Baicalein.png').convert('RGB')
fes3d = Image.open(fes_dir / 'fes3d_Baicalein.png').convert('RGB')
# Coordinates follow outputs/assemble_composite.py: third column, first row of S8.
complex_pose = im(18).crop((2576, 120, 3676, 1020))
erase(complex_pose, [(0, 0, 105, 105)])
interaction = Image.open(fes_dir / 'interact_Baicalein.png').convert('RGB')
# Remove the original A/B/C panel letters without altering titles or data.
fingerprints = relabel(im(16), [(195, 0, 30, 55, ''),
                                (1205, 0, 30, 55, ''),
                                (2347, 0, 30, 55, '')], 40)
sheet3 = canvas([[fes2d, fes3d, complex_pose, interaction], [fingerprints]])
top_cell_w = (3600 - 80 - 3 * 45) // 4
top_row_h = max(round(source.height * top_cell_w / source.width)
                for source in (fes2d, fes3d, complex_pose, interaction))
bottom_y = 40 + top_row_h + 45
labels3 = [(55 + i * (top_cell_w + 45), 55, letter) for i, letter in enumerate('ABCD')]
bottom_cell_w = (3600 - 80 - 2 * 45) // 3
labels3 += [(55 + i * (bottom_cell_w + 45), bottom_y + 15, letter)
            for i, letter in enumerate('EFG')]
label_sheet(sheet3, labels3, 70)
save(sheet3, 'Figure_3_baicalein_biophysics')

# Figure 5: erase all legacy labels, assemble, then apply one common label style.
wb = im(23)
erase(wb, [
    (105, 20, 205, 120), (1415, 20, 1515, 120), (2640, 20, 2745, 120),
    (105, 1250, 205, 1355), (1415, 1250, 1515, 1355), (2640, 1250, 2745, 1355),
])
el = im(24)
erase(el, [(0, 0, 150, 145), (1615, 0, 1790, 145), (3200, 0, 3460, 155)])
sheet5 = canvas([[wb], [el]])
wb_scale = 3520 / wb.width
wb_y = 40
el_y = 40 + round(wb.height * wb_scale) + 45
wb_positions = [(130, 37, 'A'), (1436, 37, 'B'), (2670, 37, 'C'),
                (130, 1279, 'D'), (1436, 1279, 'E'), (2670, 1279, 'F')]
labels5 = [(40 + round(x * wb_scale), wb_y + round(y * wb_scale), letter)
           for x, y, letter in wb_positions]
el_scale = 3520 / el.width
labels5 += [(40 + round(x * el_scale), el_y + round(8 * el_scale), letter)
            for x, letter in ((25, 'G'), (1640, 'H'), (3270, 'I'))]
label_sheet(sheet5, labels5, 80)
save(sheet5, 'Figure_5_WB_ELISA')
print('ALL_V84_COMPOSITES_DONE')
