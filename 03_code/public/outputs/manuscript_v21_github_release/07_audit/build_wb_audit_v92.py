# Public source snapshot: replace path placeholders before running.
from pathlib import Path
from hashlib import sha256
from math import sqrt

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.dimensions import ColumnDimension
from PIL import Image
from scipy import stats

ROOT = Path(r"<PROJECT_ROOT>\outputs\manuscript_v21_github_release")
RAW_DIR = ROOT / "06_source_material" / "wetlab_wb_raw_20261030"
SOURCE = RAW_DIR / "10.14灰度值(1).xlsx"
OUT_DIR = ROOT / "04_data" / "wb_audit"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT = OUT_DIR / "WB原始全膜与灰度核对_v1_20261030.xlsx"

src_wb = load_workbook(SOURCE, data_only=True, read_only=True)
src = src_wb["Sheet1"]
# Original source columns: rep1 C:E, rep2 G:I, rep3 K:M.
source_columns = {1: (3, 4, 5), 2: (7, 8, 9), 3: (11, 12, 13)}
source_rows = {
    "p-IκBα": 2,
    "IκBα": 3,
    "p-p65": 5,
    "p65": 6,
    "Nuclear p65": 8,
    "Cytoplasmic p65": 9,
    "β-actin": 10,
}
groups = ["Control", "HKSA+Baicalein", "HKSA"]
proteins = list(source_rows)
raw = {}
for rep, columns in source_columns.items():
    for group_index, group in enumerate(groups):
        for protein, row in source_rows.items():
            value = src.cell(row, columns[group_index]).value
            if not isinstance(value, (int, float)):
                raise ValueError(f"Missing numeric source value: rep{rep} {group} {protein}")
            raw[(rep, group, protein)] = float(value)
src_wb.close()

metrics = [
    ("p-IκBα/IκBα", "ratio", "p-IκBα", "IκBα"),
    ("p-p65/p65", "ratio", "p-p65", "p65"),
    ("Nuclear/Cytoplasmic p65", "ratio", "Nuclear p65", "Cytoplasmic p65"),
    ("Nuclear p65 fraction", "fraction", "Nuclear p65", "Cytoplasmic p65"),
    ("Cytoplasmic p65 fraction", "fraction_cyto", "Nuclear p65", "Cytoplasmic p65"),
]

def value_for(metric, rep, group):
    name, kind, a, b = metric
    av = raw[(rep, group, a)]
    bv = raw[(rep, group, b)]
    if kind == "ratio":
        return av / bv
    if kind == "fraction":
        return av / (av + bv)
    return bv / (av + bv)

# Independent numerical QA, separate from spreadsheet formulas.
print("INDEPENDENT_QA")
all_values = {}
for metric in metrics:
    name = metric[0]
    all_values[name] = {}
    print(f"\n{name}")
    for group in groups:
        vals = [value_for(metric, rep, group) for rep in (1, 2, 3)]
        all_values[name][group] = vals
        print(group, "values=", ", ".join(f"{v:.9f}" for v in vals),
              f"mean={sum(vals)/3:.9f}", f"sample_SD={stats.tstd(vals):.9f}")
    f_stat, anova_p = stats.f_oneway(*(all_values[name][g] for g in groups))
    print(f"ANOVA F={f_stat:.9f}; P={anova_p:.12g}")
    # Reproduce the manuscript's pairwise equal-variance Student t tests with threefold Bonferroni adjustment.
    for a, b in ((groups[0], groups[1]), (groups[0], groups[2]), (groups[1], groups[2])):
        t_stat, raw_p = stats.ttest_ind(all_values[name][a], all_values[name][b], equal_var=True)
        bonf_p = min(1.0, 3 * raw_p)
        print(f"{a} vs {b}: t={t_stat:.9f}; df=4; Bonferroni_P={bonf_p:.12g}")

control = sum(all_values["p-IκBα/IκBα"]["Control"]) / 3
baicalein = sum(all_values["p-IκBα/IκBα"]["HKSA+Baicalein"]) / 3
hksa = sum(all_values["p-IκBα/IκBα"]["HKSA"]) / 3
inhibition = (hksa - baicalein) / (hksa - control)
print(f"\np-IκBα/IκBα induction-normalized inhibition={inhibition:.9f} ({inhibition*100:.3f}%)")
print("Manuscript p-IκBα Control vs Baicalein rounded Bonferroni P: 0.003")

wb = Workbook()
wb.remove(wb.active)
navy = "244F70"
light_blue = "EAF1F6"
light_grey = "F3F5F7"
white = "FFFFFF"
black = "17232D"
thin_grey = Side(style="thin", color="B7C1C9")
medium_navy = Side(style="medium", color=navy)


def setup(ws, widths):
    ws.sheet_view.showGridLines = False
    ws.freeze_panes = "B4"
    for col, width in widths.items():
        ws.column_dimensions[col].width = width
    ws.sheet_properties.pageSetUpPr.fitToPage = True
    ws.page_setup.fitToWidth = 1
    ws.page_setup.fitToHeight = 0


def title(ws, text, end_col):
    ws.merge_cells(start_row=1, start_column=2, end_row=1, end_column=end_col)
    c = ws.cell(1, 2, text)
    c.font = Font(name="Arial", size=16, bold=True, color=navy)
    c.alignment = Alignment(vertical="center")
    ws.row_dimensions[1].height = 28


def header_row(ws, row, start_col, end_col):
    for c in ws.iter_cols(min_col=start_col, max_col=end_col, min_row=row, max_row=row):
        cell = c[0]
        cell.font = Font(name="Arial", size=10, bold=True, color=white)
        cell.fill = PatternFill("solid", fgColor=navy)
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = Border(bottom=medium_navy)
    ws.row_dimensions[row].height = 30


def body_style(ws, min_row, max_row, min_col, max_col):
    for row in ws.iter_rows(min_row=min_row, max_row=max_row, min_col=min_col, max_col=max_col):
        for cell in row:
            cell.font = Font(name="Arial", size=10, color=black)
            cell.alignment = Alignment(vertical="center", wrap_text=True)
            cell.border = Border(bottom=Side(style="hair", color="D9DEE2"))

# Raw grey values in a compact, source-traceable layout.
ws = wb.create_sheet("原始灰度")
setup(ws, {"A": 3, "B": 24, "C": 15, "D": 18, "E": 15, "F": 15, "G": 18, "H": 15, "I": 15, "J": 18, "K": 15, "L": 18, "M": 18})
title(ws, "WB 原始灰度值（来源工作簿逐格转录）", 13)
headers = ["蛋白"] + [f"rep{rep}\n{group}" for rep in (1, 2, 3) for group in groups]
for index, text in enumerate(headers, 2):
    ws.cell(3, index, text)
header_row(ws, 3, 2, 14)
for r, protein in enumerate(proteins, 4):
    ws.cell(r, 2, protein)
    out_col = 3
    for rep in (1, 2, 3):
        for group in groups:
            cell = ws.cell(r, out_col, raw[(rep, group, protein)])
            cell.number_format = "0"
            source_col = source_columns[rep][groups.index(group)]
            cell.comment = None
            out_col += 1
body_style(ws, 4, 10, 2, 14)
ws.cell(12, 2, "来源")
ws.cell(12, 3, "06_source_material/wetlab_wb_raw_20261030/10.14灰度值(1).xlsx")
ws.cell(13, 2, "来源位置")
ws.cell(13, 3, "Sheet1：rep1=C:E；rep2=G:I；rep3=K:M；组序均为 Control／HKSA+Baicalein／HKSA")
ws.cell(14, 2, "内参口径")
ws.cell(14, 3, "原表 actin 行；按用户确认在稿件中写为 β-actin，不改写为 GAPDH")
for row in (12, 13, 14):
    ws.cell(row, 2).font = Font(name="Arial", bold=True, color=navy)
    ws.cell(row, 3).alignment = Alignment(wrap_text=True, vertical="top")
ws.merge_cells("C12:N12")
ws.merge_cells("C13:N13")
ws.merge_cells("C14:N14")

# Cell map provides exact audit trail to source workbook coordinates.
ws = wb.create_sheet("原始单元格映射")
setup(ws, {"A": 3, "B": 24, "C": 12, "D": 20, "E": 16, "F": 16})
title(ws, "原始灰度单元格映射", 6)
for col, value in enumerate(["蛋白", "重复", "分组", "来源单元格", "原始数值"], 2):
    ws.cell(3, col, value)
header_row(ws, 3, 2, 6)
r = 4
for protein, source_row in source_rows.items():
    for rep in (1, 2, 3):
        for gi, group in enumerate(groups):
            source_col = source_columns[rep][gi]
            ws.cell(r, 2, protein)
            ws.cell(r, 3, f"rep{rep}")
            ws.cell(r, 4, group)
            ws.cell(r, 5, f"Sheet1!{get_column_letter(source_col)}{source_row}")
            ws.cell(r, 6, raw[(rep, group, protein)])
            r += 1
body_style(ws, 4, r - 1, 2, 6)
ws.auto_filter.ref = f"B3:F{r-1}"

# Formula-driven replicate calculations.
ws = wb.create_sheet("逐重复计算")
setup(ws, {"A": 3, "B": 30, "C": 12, "D": 18, "E": 20, "F": 18})
title(ws, "逐重复派生指标（全部公式）", 6)
for col, value in enumerate(["指标", "重复", "Control", "HKSA+Baicalein", "HKSA"], 2):
    ws.cell(3, col, value)
header_row(ws, 3, 2, 6)
metric_rows = {}
raw_row_by_protein = {p: i + 4 for i, p in enumerate(proteins)}
r = 4
for metric in metrics:
    name, kind, a, b = metric
    metric_rows[name] = []
    for rep in (1, 2, 3):
        ws.cell(r, 2, name)
        ws.cell(r, 3, f"rep{rep}")
        for gi, group in enumerate(groups):
            raw_col = get_column_letter(3 + (rep - 1) * 3 + gi)
            ar = raw_row_by_protein[a]
            br = raw_row_by_protein[b]
            if kind == "ratio":
                formula = f"='原始灰度'!{raw_col}{ar}/'原始灰度'!{raw_col}{br}"
            elif kind == "fraction":
                formula = f"='原始灰度'!{raw_col}{ar}/('原始灰度'!{raw_col}{ar}+'原始灰度'!{raw_col}{br})"
            else:
                formula = f"='原始灰度'!{raw_col}{br}/('原始灰度'!{raw_col}{ar}+'原始灰度'!{raw_col}{br})"
            cell = ws.cell(r, 4 + gi, formula)
            cell.number_format = "0.000000"
        metric_rows[name].append(r)
        r += 1
body_style(ws, 4, r - 1, 2, 6)

# Group summaries by formulas only.
ws = wb.create_sheet("组均值与SD")
setup(ws, {"A": 3, "B": 30, "C": 14, "D": 18, "E": 20, "F": 18})
title(ws, "组均值与样本标准差（n=3，全部公式）", 6)
for col, value in enumerate(["指标", "统计量", "Control", "HKSA+Baicalein", "HKSA"], 2):
    ws.cell(3, col, value)
header_row(ws, 3, 2, 6)
summary_rows = {}
r = 4
for metric in metrics:
    name = metric[0]
    mean_row, sd_row = r, r + 1
    summary_rows[name] = {"mean": mean_row, "sd": sd_row}
    ws.cell(mean_row, 2, name)
    ws.cell(mean_row, 3, "Mean")
    ws.cell(sd_row, 3, "Sample SD")
    rows = metric_rows[name]
    for gi in range(3):
        col = get_column_letter(4 + gi)
        source_col = get_column_letter(4 + gi)
        refs = ",".join(f"'逐重复计算'!{source_col}{rr}" for rr in rows)
        ws.cell(mean_row, 4 + gi, f"=AVERAGE({refs})").number_format = "0.000000"
        ws.cell(sd_row, 4 + gi, f"=STDEV.S({refs})").number_format = "0.000000"
    r += 2
body_style(ws, 4, r - 1, 2, 6)

# ANOVA table, formula-driven from replicate calculations.
ws = wb.create_sheet("ANOVA")
setup(ws, {"A": 3, "B": 30, "C": 12, "D": 12, "E": 14, "F": 14, "G": 14, "H": 14, "I": 16})
title(ws, "单因素方差分析（公式）", 9)
for col, value in enumerate(["指标", "N", "df between", "df error", "SS between", "SS error", "MS error", "F", "P"], 2):
    ws.cell(3, col, value)
header_row(ws, 3, 2, 10)
anova_rows = {}
r = 4
for metric in metrics:
    name = metric[0]
    anova_rows[name] = r
    ws.cell(r, 2, name)
    ws.cell(r, 3, 9)
    ws.cell(r, 4, 2)
    ws.cell(r, 5, 6)
    mean_row = summary_rows[name]["mean"]
    sd_row = summary_rows[name]["sd"]
    # Three observations per group; SS_error=(n-1)*SD^2 per group.
    ws.cell(r, 6, f"=3*SUMPRODUCT(('组均值与SD'!D{mean_row}:F{mean_row}-AVERAGE('组均值与SD'!D{mean_row}:F{mean_row}))^2)")
    ws.cell(r, 7, f"=2*SUMPRODUCT('组均值与SD'!D{sd_row}:F{sd_row}^2)")
    ws.cell(r, 8, f"=G{r}/E{r}")
    ws.cell(r, 9, f"=(F{r}/C{r})/H{r}")
    ws.cell(r, 10, f"=F.DIST.RT(I{r},D{r},E{r})")
    for c in range(6, 11):
        ws.cell(r, c).number_format = "0.00000000"
    r += 1
body_style(ws, 4, r - 1, 2, 10)

# Pairwise equal-variance Student t tests with Bonferroni adjustment, matching the manuscript pipeline.
ws = wb.create_sheet("两两比较")
setup(ws, {"A": 3, "B": 30, "C": 26, "D": 14, "E": 14, "F": 14, "G": 10, "H": 16, "I": 18})
title(ws, "两两 Student t 检验与 Bonferroni 校正（公式）", 9)
for col, value in enumerate(["指标", "比较", "均值1", "均值2", "t", "df", "原始 P", "Bonferroni P"], 2):
    ws.cell(3, col, value)
header_row(ws, 3, 2, 9)
pair_rows = {}
r = 4
for metric in metrics:
    name = metric[0]
    mean_row = summary_rows[name]["mean"]
    sd_row = summary_rows[name]["sd"]
    for pair_name, c1, c2 in (("Control vs HKSA+Baicalein", "D", "E"), ("Control vs HKSA", "D", "F"), ("HKSA+Baicalein vs HKSA", "E", "F")):
        pair_rows[(name, pair_name)] = r
        ws.cell(r, 2, name)
        ws.cell(r, 3, pair_name)
        ws.cell(r, 4, f"='组均值与SD'!{c1}{mean_row}")
        ws.cell(r, 5, f"='组均值与SD'!{c2}{mean_row}")
        ws.cell(r, 6, f"=(D{r}-E{r})/SQRT('组均值与SD'!{c1}{sd_row}^2/3+'组均值与SD'!{c2}{sd_row}^2/3)")
        ws.cell(r, 7, 4)
        ws.cell(r, 8, f"=T.DIST.2T(ABS(F{r}),G{r})")
        ws.cell(r, 9, f"=MIN(1,3*H{r})")
        for c in (4, 5, 6, 8, 9):
            ws.cell(r, c).number_format = "0.00000000"
        r += 1
body_style(ws, 4, r - 1, 2, 9)

# Comparison against values explicitly printed in manuscript v91.
manuscript = {
    "p-IκBα/IκBα": {"Control": (0.265, 0.026), "HKSA+Baicalein": (0.513, 0.042), "HKSA": (0.965, 0.022)},
    "p-p65/p65": {"Control": (0.124, 0.004), "HKSA+Baicalein": (0.381, 0.012), "HKSA": (0.775, 0.003)},
    "Nuclear/Cytoplasmic p65": {"Control": (0.232, 0.013), "HKSA+Baicalein": (1.005, 0.080), "HKSA": (3.822, 0.139)},
    "Nuclear p65 fraction": {"Control": (0.188, 0.009), "HKSA+Baicalein": (0.501, 0.020), "HKSA": (0.792, 0.006)},
}
ws = wb.create_sheet("现稿核对")
setup(ws, {"A": 3, "B": 30, "C": 20, "D": 14, "E": 14, "F": 16, "G": 16, "H": 14, "I": 14, "J": 14})
title(ws, "与 v91 正文数值逐项核对", 10)
for col, value in enumerate(["指标", "分组", "v91 Mean", "v91 SD", "重算 Mean", "重算 SD", "Mean 差", "SD 差", "三位小数核对"], 2):
    ws.cell(3, col, value)
header_row(ws, 3, 2, 10)
r = 4
for name, by_group in manuscript.items():
    for gi, group in enumerate(groups):
        old_mean, old_sd = by_group[group]
        mean_row = summary_rows[name]["mean"]
        sd_row = summary_rows[name]["sd"]
        source_col = get_column_letter(4 + gi)
        ws.cell(r, 2, name)
        ws.cell(r, 3, group)
        ws.cell(r, 4, old_mean)
        ws.cell(r, 5, old_sd)
        ws.cell(r, 6, f"='组均值与SD'!{source_col}{mean_row}")
        ws.cell(r, 7, f"='组均值与SD'!{source_col}{sd_row}")
        ws.cell(r, 8, f"=F{r}-D{r}")
        ws.cell(r, 9, f"=G{r}-E{r}")
        ws.cell(r, 10, f'=IF(AND(ROUND(F{r},3)=D{r},ROUND(G{r},3)=E{r}),"一致","需核查")')
        for c in range(4, 10):
            ws.cell(r, c).number_format = "0.000000"
        r += 1
body_style(ws, 4, r - 1, 2, 10)

# Executive summary and inhibition calculation.
ws = wb.create_sheet("核对摘要", 0)
setup(ws, {"A": 3, "B": 30, "C": 18, "D": 18, "E": 20, "F": 18, "G": 18, "H": 20})
title(ws, "WB 原始全膜与灰度核对摘要", 8)
summary_notes = [
    ("数据源", SOURCE.name),
    ("内参", "β-actin（按用户确认沿用；原始表为 actin 行）"),
    ("重复", "rep1／rep2／rep3；每组 n=3"),
    ("分组顺序", "Control／HKSA+Baicalein／HKSA"),
    ("原始全膜", "rep1／rep2／rep3 TIFF 已归档并生成缩略核对图；核对图不用于重新定量"),
    ("核对范围", "由原灰度表重算比率、核占比、均值、样本 SD 与 ANOVA/Bonferroni；未重新圈选条带"),
]
for r, (label, text) in enumerate(summary_notes, 3):
    ws.cell(r, 2, label).font = Font(name="Arial", bold=True, color=navy)
    ws.cell(r, 3, text)
    ws.merge_cells(start_row=r, start_column=3, end_row=r, end_column=8)
    ws.cell(r, 3).alignment = Alignment(wrap_text=True, vertical="top")
ws.cell(11, 2, "p-IκBα/IκBα 诱导增量抑制率")
ws.cell(11, 3, "公式")
cm = summary_rows["p-IκBα/IκBα"]["mean"]
ws.cell(11, 4, f"=('组均值与SD'!F{cm}-'组均值与SD'!E{cm})/('组均值与SD'!F{cm}-'组均值与SD'!D{cm})")
ws.cell(11, 4).number_format = "0.0%"
ws.cell(12, 2, "v91 正文")
ws.cell(12, 4, 0.646).number_format = "0.0%"
ws.cell(13, 2, "核对")
ws.cell(13, 4, '=IF(ROUND(D11*100,1)=ROUND(D12*100,1),"一致","需核查")')
ws.cell(15, 2, "精确 P 值提示")
ws.cell(15, 3, "p-IκBα/IκBα：Control vs HKSA+Baicalein 按等方差 Student t 检验并作 3 次 Bonferroni 校正，精确 P=0.002944，三位小数为 0.003，与 v91 正文一致。旧绘图统计文件中的 Welch 未校正 P 不再用于图 5。")
ws.merge_cells("C15:H16")
ws.cell(15, 3).alignment = Alignment(wrap_text=True, vertical="top")
ws.cell(18, 2, "均值与 SD 总结")
ws.cell(18, 2).font = Font(name="Arial", bold=True, color=navy)
ws.cell(18, 3, '=IF(COUNTIF(\'现稿核对\'!J4:J15,"需核查")=0,"正文列出的四项指标均值与 SD 均按三位小数一致","存在需核查项")')
ws.merge_cells("C18:H18")
ws.cell(20, 2, "主要比较 P 值")
ws.cell(20, 2).font = Font(name="Arial", bold=True, color=navy)
for col, value in enumerate(["指标", "比较", "Bonferroni P"], 3):
    ws.cell(21, col, value)
header_row(ws, 21, 3, 5)
r = 22
for pair_name in ("Control vs HKSA+Baicalein", "Control vs HKSA", "HKSA+Baicalein vs HKSA"):
    ws.cell(r, 3, "p-IκBα/IκBα")
    ws.cell(r, 4, pair_name)
    ws.cell(r, 5, f"='两两比较'!I{pair_rows[('p-IκBα/IκBα', pair_name)]}")
    ws.cell(r, 5).number_format = "0.00000000"
    r += 1
body_style(ws, 22, r - 1, 3, 5)

# File inventory and checksums for the archived raw source material.
ws = wb.create_sheet("原文件清单")
setup(ws, {"A": 3, "B": 44, "C": 16, "D": 18, "E": 68})
title(ws, "WB 原始资料文件清单与 SHA-256", 5)
for col, value in enumerate(["文件名", "字节数", "像素尺寸", "SHA-256"], 2):
    ws.cell(3, col, value)
header_row(ws, 3, 2, 5)
r = 4
for path in sorted(RAW_DIR.iterdir(), key=lambda p: p.name.lower()):
    if not path.is_file():
        continue
    dims = ""
    if path.suffix.lower() in {".tif", ".tiff", ".jpg", ".jpeg", ".png"}:
        with Image.open(path) as im:
            dims = f"{im.width}×{im.height}"
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    ws.cell(r, 2, path.name)
    ws.cell(r, 3, path.stat().st_size)
    ws.cell(r, 4, dims)
    ws.cell(r, 5, digest.hexdigest())
    r += 1
body_style(ws, 4, r - 1, 2, 5)
ws.auto_filter.ref = f"B3:E{r-1}"

# Sources and interpretation boundaries.
ws = wb.create_sheet("来源说明")
setup(ws, {"A": 3, "B": 24, "C": 100})
title(ws, "来源与核对口径", 3)
source_notes = [
    ("灰度工作簿", "06_source_material/wetlab_wb_raw_20261030/10.14灰度值(1).xlsx"),
    ("原始膜目录", "06_source_material/wetlab_wb_raw_20261030/"),
    ("正文基线", "01_manuscript/SHXD_PJI_中文全稿_v91_MD参数与ELISA重复伦理核定.docx"),
    ("缩略核对图", "07_audit/wb_raw_20261030/"),
    ("计算定义", "p-IκBα/IκBα 与 p-p65/p65 为同次重复内直比；核/质 p65 为 Nuclear/Cytoplasmic；核占比为 Nuclear/(Nuclear+Cytoplasmic)。"),
    ("统计定义", "均值与样本 SD，n=3；one-way ANOVA；两两比较采用等方差 Student t 检验（每对 df=4），3 次比较 Bonferroni 校正，与现稿统计管线一致。"),
    ("内参说明", "用户确认沿用此前稿件的 β-actin；附件中其他内参记录不覆盖本次选择。"),
    ("外部来源 URL", "无；本工作簿只使用本地原始实验文件与现稿。"),
]
for r, (label, text) in enumerate(source_notes, 3):
    ws.cell(r, 2, label).font = Font(name="Arial", bold=True, color=navy)
    ws.cell(r, 3, text)
    ws.cell(r, 3).alignment = Alignment(wrap_text=True, vertical="top")
    ws.row_dimensions[r].height = 34

# Workbook-level usability and recalculation settings.
for sheet in wb.worksheets:
    sheet.sheet_properties.pageSetUpPr.fitToPage = True
    sheet.sheet_view.zoomScale = 90
wb.calculation.fullCalcOnLoad = True
wb.calculation.forceFullCalc = True
wb.calculation.calcMode = "auto"
wb.save(OUT)
print(f"\nWORKBOOK={OUT}")
