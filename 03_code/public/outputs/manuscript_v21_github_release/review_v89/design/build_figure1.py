# Public source snapshot: replace path placeholders before running.
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
import yaml

ROOT = Path(r"<PROJECT_ROOT>\outputs\manuscript_v21_github_release\review_v89")
DESIGN = ROOT / "design" / "fig1_integrated_workflow"
MEDIA = DESIGN / "media"
PAGES = DESIGN / "pages"
MEDIA.mkdir(parents=True, exist_ok=True)
PAGES.mkdir(parents=True, exist_ok=True)
SRC = ROOT / "work" / "unpacked" / "word" / "media" / "image42.png"
OUT = ROOT / "design" / "Figure1_integrated_workflow.png"
W, H = 3600, 1883
SPLIT = 1980

C = {
    "navy": "#244F70", "navy_light": "#EAF1F6",
    "teal": "#2D7B73", "teal_light": "#EAF4F1",
    "gold": "#A87922", "gold_light": "#FBF4E6",
    "plum": "#72558B", "plum_light": "#F2EDF6",
    "red": "#A84B4B", "red_light": "#FAEEEE",
    "ink": "#17232D", "muted": "#5E6B75", "line": "#AAB6BF",
    "white": "#FFFFFF", "pale": "#F7F9FA"
}

def rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

canvas = Image.new("RGB", (W, H), rgb(C["white"]))
draw = ImageDraw.Draw(canvas)
orig = Image.open(SRC).convert("RGB")
assert orig.size == (W, H), orig.size

# Preserve the quantitative panels B and C without resampling; remove only the matrix-dimension label.
b = orig.crop((SPLIT, 0, W, 1100))
ImageDraw.Draw(b).rectangle((2940-SPLIT, 58, 3180-SPLIT, 166), fill=rgb(C["white"]))
c = orig.crop((SPLIT, 1080, W, H))
b.save(MEDIA / "panel_B_original.png")
c.save(MEDIA / "panel_C_original.png")
canvas.paste(b, (SPLIT, 0))
canvas.paste(c, (SPLIT, 1080))

font_regular = r"<LOCAL_DRIVE>/Windows\Fonts\arial.ttf"
font_bold = r"<LOCAL_DRIVE>/Windows\Fonts\arialbd.ttf"
def font(size, bold=False):
    return ImageFont.truetype(font_bold if bold else font_regular, size)

elements = []
def shape(eid, bounds, fill, border=None, shape_name="rect", adjustments=None):
    x, y, w, h = bounds
    draw.rectangle((x, y, x+w, y+h), fill=rgb(fill), outline=rgb(border) if border else None, width=3 if border else 1)
    item = {"elementId": eid, "elementType": "shape", "bounds": bounds,
            "shapeName": shape_name, "fill": {"type": "solid", "color": fill}}
    if border:
        item["border"] = {"style": "solid", "width": 3, "color": border}
    if adjustments is not None:
        item["adjustments"] = adjustments
    elements.append(item)

def text(eid, bounds, value, size=36, color=None, bold=False, align="center"):
    color = color or C["ink"]
    x, y, w, h = bounds
    f = font(size, bold)
    # Pillow and PPTD share explicit line breaks; no implicit layout guessing.
    draw.multiline_text((x+w/2, y+h/2), value, font=f, fill=rgb(color),
                        anchor="mm", align="center" if align == "center" else "left", spacing=max(4, int(size*0.18)))
    elements.append({"elementId": eid, "elementType": "text", "bounds": bounds,
                     "content": {"text": value, "fontSize": size, "fontFamily": "Arial",
                                 "color": color, "bold": bold, "lineHeight": 1.12,
                                 "align": [align, "middle"]}})

def line(eid, bounds, points, viewbox, color=None, width=4, arrow_end=False):
    color = color or C["line"]
    x, y, w, h = bounds
    pts = []
    for px, py in points:
        pts.append((x + px/viewbox[0]*w, y + py/viewbox[1]*h))
    draw.line(pts, fill=rgb(color), width=width)
    if arrow_end:
        x2, y2 = pts[-1]
        x1, y1 = pts[-2]
        import math
        ang = math.atan2(y2-y1, x2-x1)
        al, aw = 22, 12
        p1 = (x2-al*math.cos(ang)+aw*math.sin(ang), y2-al*math.sin(ang)-aw*math.cos(ang))
        p2 = (x2-al*math.cos(ang)-aw*math.sin(ang), y2-al*math.sin(ang)+aw*math.cos(ang))
        draw.polygon([(x2, y2), p1, p2], fill=rgb(color))
    elements.append({"elementId": eid, "elementType": "line", "bounds": bounds,
                     "viewBox": list(viewbox), "points": " ".join(f"{px},{py}" for px, py in points),
                     "curve": "sharp", "arrow": [None, "arrow"] if arrow_end else [None, None],
                     "border": {"style": "solid", "width": width, "color": color}})

# Panel heading.
text("panel_A", [48, 28, 100, 90], "A", 72, C["ink"], True)
text("title", [170, 34, 1740, 64], "From SHXD chemical space to experimental validation", 50, C["ink"], True)
text("subtitle", [170, 101, 1740, 52], "Network overlap guides candidate prioritization; evidence layers remain distinct.", 30, C["muted"])

# Connectors are drawn first so that no line crosses text.
line("arrow_top_left", [520, 448, 280, 112], [(0, 0), (280, 112)], (280, 112), C["navy"], 4, True)
line("arrow_top_right", [1240, 448, 220, 112], [(220, 0), (0, 112)], (220, 112), C["teal"], 4, True)
line("arrow_network_candidate", [1038, 842, 2, 66], [(1, 0), (1, 66)], (2, 66), C["navy"], 4, True)
line("arrow_candidate_structure", [1038, 1108, 2, 62], [(1, 0), (1, 62)], (2, 62), C["gold"], 4, True)
line("arrow_structure_split", [1038, 1418, 2, 64], [(1, 0), (1, 64)], (2, 64), C["plum"], 4, True)
line("split_horizontal", [520, 1480, 940, 2], [(0, 1), (940, 1)], (940, 2), C["plum"], 4, False)
line("split_left", [518, 1480, 2, 32], [(1, 0), (1, 32)], (2, 32), C["plum"], 4, True)
line("split_right", [1458, 1480, 2, 32], [(1, 0), (1, 32)], (2, 32), C["plum"], 4, True)

# Input layers.
shape("input_shxd", [100, 190, 840, 258], C["navy_light"], C["navy"])
shape("input_shxd_rule", [100, 190, 840, 12], C["navy"])
text("input_shxd_title", [140, 218, 760, 58], "SHXD chemical records", 44, C["navy"], True)
text("input_shxd_body", [140, 285, 760, 125], "58 formula records · 98 plasma records\nUPLC-HRMS chemistry and exposure clues\n16 candidate compounds", 34, C["ink"])

shape("input_pji", [1040, 190, 840, 258], C["teal_light"], C["teal"])
shape("input_pji_rule", [1040, 190, 840, 12], C["teal"])
text("input_pji_title", [1080, 218, 760, 58], "PJI disease evidence", 44, C["teal"], True)
text("input_pji_body", [1080, 280, 760, 140], "GSE255786: 53 PJI / 40 AF\nPXD038928 · PXD058084 · literature\nedgeR differential expression and enrichment", 32, C["ink"])

# Network intersection and candidate nomination.
shape("network_band", [290, 560, 1600, 282], C["pale"], C["navy"])
shape("network_rule", [290, 560, 14, 282], C["navy"])
text("network_title", [350, 582, 1480, 58], "Network pharmacology intersection", 46, C["navy"], True)
text("network_body", [350, 650, 1480, 72], "Herb–compound–target network + disease evidence + topology", 34, C["ink"])
text("network_counts", [340, 724, 1500, 46], "51-target set  →  45 detected  →  33 FDR < 0.05  →  28 strict UP", 34, C["navy"], True)
text("network_nodes", [350, 779, 1480, 42], "Direct overlap: IL1B · IL6 · CXCL8 · MMP3", 32, C["muted"])

shape("candidate_band", [190, 910, 1700, 198], C["gold_light"], C["gold"])
text("candidate_title", [250, 926, 1580, 48], "Candidate monomers", 42, C["gold"], True)
text("candidate_list", [240, 980, 1600, 48], "Baicalein · Berberine · Emodin · Wogonin · Apigenin · Oroxylin A", 32, C["ink"])
text("candidate_focus", [250, 1038, 1580, 44], "Baicalein advanced to cell-based functional validation", 34, C["gold"], True)

shape("structure_band", [190, 1170, 1700, 248], C["plum_light"], C["plum"])
shape("structure_rule", [190, 1170, 14, 248], C["plum"])
text("structure_title", [250, 1190, 1580, 52], "Structure-based evaluation", 42, C["plum"], True)
text("structure_body1", [250, 1254, 1580, 52], "Vina  →  triplicate 100 ns MD  →  MM/GBSA", 36, C["ink"])
text("structure_body2", [250, 1315, 1580, 70], "Contact fingerprints · free-energy landscapes\nRepresentative poses · computational mutants", 32, C["muted"])

# Two experimental arms.
shape("host_arm", [100, 1512, 840, 252], C["red_light"], C["red"])
text("host_title", [140, 1532, 760, 52], "Host inflammatory response", 38, C["red"], True)
text("host_body", [140, 1596, 760, 128], "PMA-differentiated THP-1 + HKSA\nBaicalein · IKK-16 perturbation\nWB + ELISA\nIKK/NF-κB-related readouts", 32, C["ink"])

shape("pathogen_arm", [1040, 1512, 840, 252], C["teal_light"], C["teal"])
text("pathogen_title", [1080, 1532, 760, 52], "Pathogen-side readouts", 38, C["teal"], True)
text("pathogen_body", [1080, 1596, 760, 128], "Commercial MRSA reference strain\nBNCC330041 · BeNa Culture Collection\nBaicalein + SHXD formula\nCFU counting + live/dead staining", 32, C["ink"])

text("integrated_readout", [170, 1780, 1710, 58], "Integrated model: SHXD-derived baicalein and IKK/NF-κB-related inflammatory regulation", 30, C["navy"], True)
line("panel_separator", [1954, 35, 2, 1800], [(1, 0), (1, 1800)], (2, 1800), "#D7DEE3", 2, False)

# Original B/C panels are embedded as separate, unaltered evidence images.
elements.extend([
    {"elementId": "panel_B_original", "elementType": "image", "bounds": [SPLIT, 0, W-SPLIT, 1100],
     "src": "media/panel_B_original.png", "fit": {"mode": "fill"}},
    {"elementId": "panel_C_original", "elementType": "image", "bounds": [SPLIT, 1080, W-SPLIT, H-1080],
     "src": "media/panel_C_original.png", "fit": {"mode": "fill"}}
])

main = {
    "version": "v2", "title": "SHXD-PJI Figure 1 integrated workflow",
    "size": [W, H],
    "theme": {"colors": C, "textStyles": {"body": {"fontFamily": "Arial", "fontSize": 34, "color": "$ink"}}},
    "pages": ["pages/figure1.page"]
}
page = {"pageType": "content", "background": {"type": "solid", "color": "#FFFFFF"}, "elements": elements}
(DESIGN / "fig1_integrated_workflow.kanva.pptd").write_text(yaml.safe_dump(main, allow_unicode=True, sort_keys=False), encoding="utf-8")
(PAGES / "figure1.page").write_text(yaml.safe_dump(page, allow_unicode=True, sort_keys=False), encoding="utf-8")
canvas.save(OUT, dpi=(600, 600), optimize=True)
print(OUT)
print(DESIGN / "fig1_integrated_workflow.kanva.pptd")
print("elements", len(elements), "size", canvas.size)
