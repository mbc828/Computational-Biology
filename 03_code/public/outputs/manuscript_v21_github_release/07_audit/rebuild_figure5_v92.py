# Public source snapshot: replace path placeholders before running.
from pathlib import Path
import shutil

from PIL import Image

ROOT = Path(r"<PROJECT_ROOT>\outputs\manuscript_v21_github_release")
FIG = ROOT / "02_figures" / "main" / "Figure_05.png"
WB = ROOT.parent / "wet_figures_20261014" / "figures" / "Fig_WB_NFkB_6panel.png"
AUDIT = ROOT / "07_audit" / "figure5_v92"
AUDIT.mkdir(parents=True, exist_ok=True)
BACKUP = AUDIT / "Figure_05_before_stats_unification.png"
if not BACKUP.exists():
    shutil.copy2(FIG, BACKUP)
shutil.copy2(WB, AUDIT / "Figure_05_WB_panels_stats_unified.png")

old = Image.open(BACKUP).convert("RGB")
wb = Image.open(WB).convert("RGB")
assert old.size == (3600, 3106), old.size
width = 3600
margin = 40
gap = 45
cell_w = width - 2 * margin
wb_h = round(wb.height * cell_w / wb.width)
wb_resized = wb.resize((cell_w, wb_h), Image.Resampling.LANCZOS)
y_elisa = margin + wb_h + gap
elisa_h = old.height - margin - y_elisa
assert elisa_h > 0, (old.size, wb.size, y_elisa)
# Preserve the already-reviewed ELISA row exactly; only the WB panels are replaced.
elisa = old.crop((margin, y_elisa, margin + cell_w, y_elisa + elisa_h))
canvas = Image.new("RGB", old.size, "white")
canvas.paste(wb_resized, (margin, margin))
canvas.paste(elisa, (margin, y_elisa))
canvas.save(FIG, optimize=True, dpi=(300, 300))
print(f"WB source={wb.size}; resized={wb_resized.size}; ELISA top={y_elisa}; ELISA height={elisa_h}")
print(FIG)
