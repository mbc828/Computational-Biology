# Public source snapshot: replace path placeholders before running.
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(r'<PROJECT_ROOT>')
V82 = ROOT / 'outputs/manuscript_v21_github_release/figure_consolidation_v82'
MEDIA = V82 / 'work/unpacked/word/media'
OUT = ROOT / 'outputs/manuscript_v21_github_release/figure_consolidation_v83/composite_figures'
OUT.mkdir(parents=True, exist_ok=True)
FONT = r'<LOCAL_DRIVE>/Windows\Fonts\arial.ttf'
FONT_B = r'<LOCAL_DRIVE>/Windows\Fonts\arialbd.ttf'


def im(n):
    return Image.open(MEDIA / f'image{n}.png').convert('RGB')


def relabel(img, replacements, size=70):
    d = ImageDraw.Draw(img)
    f = ImageFont.truetype(FONT_B, size)
    for x, y, w, h, letter in replacements:
        d.rectangle((x - 8, y - 5, x + w + 8, y + h + 5), fill='white')
        d.text((x, y), letter, font=f, fill='#20252b')
    return img


def add_label(img, letter, x=25, y=20, size=75):
    return relabel(img, [(x, y, size, size, letter)], size)


def with_header(img, text, height=150, size=62):
    panel = Image.new('RGB', (img.width, img.height + height), 'white')
    panel.paste(img, (0, height))
    d = ImageDraw.Draw(panel)
    d.text((35, 25), text, font=ImageFont.truetype(FONT_B, size), fill='#20252b')
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
        row_h = max(x.height for x in scaled)
        placed.append((y, cell_w, scaled, row_h))
        y += row_h + gap
    sheet = Image.new('RGB', (width, y - gap + margin), 'white')
    for y0, cell_w, images, row_h in placed:
        for j, image in enumerate(images):
            x = margin + j * (cell_w + gap) + (cell_w - image.width) // 2
            sheet.paste(image, (x, y0))
    return sheet


def save(sheet, name):
    path = OUT / f'{name}.png'
    sheet.save(path, optimize=True, dpi=(300, 300))
    prev = sheet.copy()
    prev.thumbnail((1800, 1800), Image.Resampling.LANCZOS)
    prev.save(OUT / f'{name}_preview.jpg', quality=92, optimize=True)
    print(name, sheet.size, path.stat().st_size, 'preview', prev.size)


# Shared clinical panels A-D.
a9 = im(9).crop((0, 330, 4251, 2500))
go = im(10)
cplot = go.crop((0, 620, 3700, 2600))
cb = go.crop((3700, 950, 4251, 4050))
cb = cb.resize((round(cb.width * cplot.height / cb.height), cplot.height), Image.Resampling.LANCZOS)
cgo = Image.new('RGB', (cplot.width + cb.width, cplot.height), 'white')
cgo.paste(cplot, (0, 0)); cgo.paste(cb, (cplot.width, 0))
c10 = with_header(cgo, 'C  UP: top 10 of 664 significant terms', 140, 66)
dplot = go.crop((0, 2850, 3700, 4600))
dgo = Image.new('RGB', (cgo.width, dplot.height), 'white')
dgo.paste(dplot, (0, 0))
d10 = with_header(dgo, 'D  DOWN: all 7 significant terms', 140, 66)

# Figure 1E-F: actual prespecified marker scores, not deconvolution fractions.
score_path = ROOT / 'outputs/supplementary_20260923/primary/runs/run_20260923T035406_845247Z_522b828cbcc0/results/marker_score/marker_scores.txt'
scores = pd.read_csv(score_path, sep='\t')
plt.rcParams.update({'font.family': 'Arial', 'font.size': 9})
marker_panels = []
for letter, key, title, fdr in [
    ('E', 'neutrophil', 'Neutrophil marker score', r'BH FDR = 2.51 $\times$ 10$^{-16}$'),
    ('F', 'myeloid', 'Myeloid marker score', r'BH FDR = 4.56 $\times$ 10$^{-8}$'),
]:
    fig, ax = plt.subplots(figsize=(5.0, 4.1), dpi=300)
    data = [scores.loc[scores.group == g, key].dropna().to_numpy() for g in ('AF', 'PJI')]
    vp = ax.violinplot(data, positions=[0, 1], widths=0.72, showextrema=False)
    for body, color in zip(vp['bodies'], ('#b9d9e8', '#f2c28b')):
        body.set_facecolor(color); body.set_edgecolor('#56616b'); body.set_alpha(0.75); body.set_linewidth(0.8)
    rng = np.random.default_rng(20261028)
    for x, vals, color in zip((0, 1), data, ('#1677a6', '#c56a16')):
        jitter = rng.uniform(-0.12, 0.12, len(vals))
        ax.scatter(np.full(len(vals), x) + jitter, vals, s=11, color=color, alpha=0.72, edgecolors='none', zorder=3)
        mean = vals.mean(); sd = vals.std(ddof=1)
        ax.errorbar(x, mean, yerr=sd, fmt='_', markersize=18, capsize=5, lw=1.4, color='#20252b', zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels([f'AF\n(n={len(data[0])})', f'PJI\n(n={len(data[1])})'])
    ax.set_ylabel('Gene z-score mean')
    ax.set_title(f'{letter}  {title}', loc='left', fontweight='bold', fontsize=11)
    ax.text(0.03, 0.97, fdr, transform=ax.transAxes, va='top', fontsize=8.5)
    ax.grid(axis='y', color='#e2e5e8', lw=0.6); ax.set_axisbelow(True)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    fig.tight_layout()
    p = OUT / f'marker_{key}.png'
    fig.savefig(p, bbox_inches='tight', pad_inches=0.08, facecolor='white')
    plt.close(fig)
    marker_panels.append(Image.open(p).convert('RGB'))
save(canvas([[a9], [c10, d10], marker_panels]), 'Figure_1_clinical_transcriptome_immune')

# Figure 2: workflow -> docking matrix -> replicate IKKbeta endpoint energies.
a1 = add_label(im(1), 'A', 30, 25, 80)
dock = im(13)
b13 = with_header(dock.crop((0, 110, dock.width // 2, dock.height)),
                  'B  Vina score matrix (6×5)', 140, 68)
ikk = im(14).crop((0, 700, 4251, 2100))
ImageDraw.Draw(ikk).rectangle((450, 1270, 1300, 1400), fill='white')
xaxis = im(14).crop((500, 3520, 4251, 3800))
xd = ImageDraw.Draw(xaxis)
xd.line((200, 0, 3500, 0), fill='#20252b', width=3)
for xx in (460, 1080, 1670, 2280, 2890, 3500):
    xd.line((xx, 0, xx, 20), fill='#20252b', width=3)
ikk_with_axis = Image.new('RGB', (4251, ikk.height + xaxis.height), 'white')
ikk_with_axis.paste(ikk, (0, 0)); ikk_with_axis.paste(xaxis, (500, ikk.height))
c14 = with_header(ikk_with_axis, 'C  IKKβ endpoint energy estimates', 140, 68)
save(canvas([[a1], [b13, c14]]), 'Figure_2_virtual_screening')

# Figure 3: baicalein free-energy landscape and structural/energy fingerprints.
fes_dir = ROOT / 'outputs/fes_panels'
fes2d = add_label(Image.open(fes_dir / 'fes2d_Baicalein.png').convert('RGB'), 'A', 20, 15, 75)
fes3d = add_label(Image.open(fes_dir / 'fes3d_Baicalein.png').convert('RGB'), 'B', 20, 15, 75)
def16 = relabel(im(16), [(195, 0, 30, 55, 'C'),
                         (1205, 0, 30, 55, 'D'),
                         (2347, 0, 30, 55, 'E')], 40)
save(canvas([[fes2d, fes3d], [def16]]), 'Figure_3_baicalein_biophysics')

# Figure 5: complete WB A-F plus ELISA G-I; all measurements remain unchanged.
wb = im(23)
el = im(24)
ImageDraw.Draw(el).rectangle((3200, 0, 3600, 150), fill='white')
el = relabel(el, [(10, 8, 115, 115, 'G'),
                  (1640, 8, 115, 115, 'H'),
                  (3270, 8, 115, 115, 'I')], 92)
save(canvas([[wb], [el]]), 'Figure_5_WB_ELISA')

# Supplementary Figure S12: original prespecified module forest plots.
e11 = with_header(im(11).crop((160, 650, 2140, 2200)), 'A  TLR2 / TLR4 / MyD88', 140, 68)
f11 = with_header(im(11).crop((160, 2400, 2140, 3800)), 'B  Canonical NF-κB', 140, 68)
save(canvas([[e11, f11]]), 'Supplementary_Figure_S12_modules')

# Supplementary Figure S13: network and descriptive centrality.
network = im(6).crop((3975, 0, 7775, 3700))
bd = ImageDraw.Draw(network)
bd.rectangle((0, 3100, 1100, 3600), fill='white')
bd.ellipse((70, 3340, 140, 3410), fill='#dceaf5', outline='#6b8ba4', width=3)
bd.text((175, 3342), 'Scutellaria compound', font=ImageFont.truetype(FONT, 52), fill='#20252b')
network = network.crop((0, 180, network.width, network.height))
ImageDraw.Draw(network).rectangle((0, 0, network.width, 160), fill='white')
network = with_header(network, 'A  Complete compound–target hypothesis network', 150, 62)
compound_src = Image.open(V82 / 'composite_figures/compound_centrality.png').convert('RGB')
compound = with_header(compound_src.crop((0, 120, compound_src.width, compound_src.height)),
                       'B  Compound centrality and connectivity', 130, 60)
target = im(7).crop((1400, 150, 3602, 1658))
ImageDraw.Draw(target).rectangle((0, 450, 25, 700), fill='white')
target = with_header(target, 'C  Top target centrality scores', 130, 60)
save(canvas([[network], [compound, target]]), 'Supplementary_Figure_S13_network')

# Supplementary Figure S14: redocking methodology QC only.
redock = dock.crop((round(dock.width * 0.58), 140, dock.width, dock.height))
redock = with_header(redock, 'Redocking benchmark', 140, 68)
save(canvas([[redock]], width=2400), 'Supplementary_Figure_S14_redocking')
print('ALL_V83_COMPOSITES_DONE')
