# Public source snapshot: replace path placeholders before running.
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(r"<PROJECT_ROOT>\outputs\manuscript_v21_github_release")
DESIGN = ROOT / "review_v89" / "design"
DESIGN.mkdir(parents=True, exist_ok=True)
source = ROOT.parent.parent / "wet_figures_20261014" / "make_fig_wb.py"
# Resolve explicitly because ROOT is the manuscript release directory.
source = ROOT.parents[1] / "wet_figures_20261014" / "make_fig_wb.py"
if not source.exists():
    source = Path(r"<PROJECT_ROOT>\outputs\wet_figures_20261014\make_fig_wb.py")
code = source.read_text(encoding="utf-8")
code = code.replace('OUT = os.path.join(BASE, "figures")', f'OUT = r"{DESIGN}"')
code = code.replace('fontsize=12, fontweight="bold", fontstyle="italic"',
                    'fontsize=12, fontweight="normal", fontstyle="normal"')
code = code.replace('Fig_WB_NFkB_6panel.png', 'Figure5_WB_Arial.png')
assert 'fontstyle="italic"' not in code
exec(compile(code, str(source), "exec"), {"__file__": str(source), "__name__": "__main__"})

wb = Image.open(DESIGN / "Figure5_WB_Arial.png").convert("RGB")
el = Image.open(ROOT / "figure_consolidation_v82" / "work" / "unpacked" / "word" / "media" / "image24.png").convert("RGB")
FONT_B = r"<LOCAL_DRIVE>/Windows\Fonts\arialbd.ttf"

def relabel(img, replacements, size=70):
    d = ImageDraw.Draw(img)
    f = ImageFont.truetype(FONT_B, size)
    for x, y, w, h, letter in replacements:
        d.rectangle((x - 8, y - 5, x + w + 8, y + h + 5), fill="white")
        d.text((x, y), letter, font=f, fill="#20252b")
    return img

ImageDraw.Draw(el).rectangle((3200, 0, 3600, 150), fill="white")
el = relabel(el, [(10, 8, 115, 115, "G"), (1640, 8, 115, 115, "H"), (3270, 8, 115, 115, "I")], 92)

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
    sheet = Image.new("RGB", (width, y - gap + margin), "white")
    for y0, cell_w, images, row_h in placed:
        for j, image in enumerate(images):
            x = margin + j * (cell_w + gap) + (cell_w - image.width) // 2
            sheet.paste(image, (x, y0))
    return sheet

sheet = canvas([[wb], [el]])
target_size = (3600, 3106)
if sheet.size != target_size:
    if sheet.width != target_size[0] or sheet.height > target_size[1]:
        raise RuntimeError(f"Unexpected composite size: {sheet.size}")
    padded = Image.new("RGB", target_size, "white")
    padded.paste(sheet, (0, 0))
    sheet = padded
out = DESIGN / "Figure5_WB_ELISA_Arial.png"
sheet.save(out, optimize=True, dpi=(300, 300))
print("WB", wb.size, "ELISA", el.size, "composite", sheet.size)
print(out)
