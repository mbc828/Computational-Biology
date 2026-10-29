# Public source snapshot: replace path placeholders before running.
# -*- coding: utf-8 -*-
"""Rebuild the 9-group ELISA figure from the latest copied workbooks.

Concentrations are independently recalculated from OD450 with each source cell's
own quadratic formula and checked against the saved formula cache. No source
workbook is modified.
"""
from pathlib import Path
import hashlib
import json
import math
import re

import numpy as np
import openpyxl
from openpyxl import Workbook
from openpyxl.chart import BarChart, Reference
from openpyxl.chart.label import DataLabelList
from openpyxl.drawing.image import Image as XLImage
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

BASE = Path(r"<PROJECT_ROOT>\outputs\elisa_redraw_20261026")
SOURCE = BASE / "source"
PNG = BASE / "Fig_ELISA_9groups_20261026.png"
PDF = BASE / "Fig_ELISA_9groups_20261026.pdf"
XLSX = BASE / "ELISA三因子_9组重算与制图_20261026.xlsx"
JSON = BASE / "elisa_recalculated_20261026.json"

PLATES = {
    "IL-1β": {
        "file": "IL-1b.xlsx", "sheet": "il-1β 3.29", "first_row": 44,
        "label_col": 2, "value_col": 3, "source_unit": "ng/ml", "plot_unit": "pg/mL",
    },
    "IL-6": {
        "file": "IL-6.xlsx", "sheet": "IL-6 4.6", "first_row": 52,
        "label_col": 1, "value_col": 2, "source_unit": "pg/ml", "plot_unit": "pg/mL",
    },
    "TNF-α": {
        "file": "TNF-a.xlsx", "sheet": "TNF-α 3.29", "first_row": 49,
        "label_col": 1, "value_col": 2, "source_unit": "pg/ml", "plot_unit": "pg/mL",
    },
}
SOURCE_ORDER = [
    "Control", "HKSA", "IKK-16", "HKSA+IKK-16", "HKSA+BAI(24)",
    "HKSA+BAI(24)+IKK-16", "HKSA+BAI(6)", "HKSA+BAI(12)", "HKSA+BAI(48)",
]
GROUPS = [
    "Control", "IKK-16", "HKSA", "HKSA+IKK-16", "HKSA+BAI(6)",
    "HKSA+BAI(12)", "HKSA+BAI(24)", "HKSA+BAI(48)", "HKSA+BAI(24)+IKK-16",
]
COLORS = ["#4D4D4D", "#9ECAE1", "#C00000", "#6BAED6", "#C7E9C0",
          "#74C476", "#238B45", "#00441B", "#3182BD"]
XLABELS = ["Control", "IKK-16", "HKSA", "+IKK-16", "+BAI (6)", "+BAI (12)",
           "+BAI (24)", "+BAI (48)", "+BAI (24)\n+IKK-16"]
PAIRS = [
    ("Control", "HKSA"),
    ("HKSA", "HKSA+BAI(24)"),
    ("HKSA", "HKSA+IKK-16"),
    ("HKSA+IKK-16", "HKSA+BAI(24)+IKK-16"),
    ("IKK-16", "Control"),
]
FORMULA_RE = re.compile(
    r"^=(?P<a>\d+(?:\.\d+)?)\*(?P<cell>[A-Z]{1,3}\d+)\*(?P=cell)"
    r"\+(?P<b>\d+(?:\.\d+)?)\*(?P=cell)\+(?P<c>\d+(?:\.\d+)?)$"
)


def significance(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def extract_plate(analyte, cfg):
    path = SOURCE / cfg["file"]
    wbf = openpyxl.load_workbook(path, data_only=False)
    wbv = openpyxl.load_workbook(path, data_only=True)
    wsf = wbf[cfg["sheet"]]
    wsv = wbv[cfg["sheet"]]
    raw = {}
    source_labels = {}
    source_cells = {}
    coefficients = set()
    max_cache_delta = 0.0
    for offset, group in enumerate(SOURCE_ORDER):
        row = cfg["first_row"] + offset
        source_labels[group] = wsf.cell(row=row, column=cfg["label_col"]).value
        ods, values, cells = [], [], []
        for rep in range(5):
            col = cfg["value_col"] + rep
            cell = wsf.cell(row=row, column=col)
            match = FORMULA_RE.fullmatch(str(cell.value))
            if not match:
                raise ValueError(f"Unexpected formula at {cfg['file']} {cell.coordinate}: {cell.value}")
            a = float(match.group("a"))
            b = float(match.group("b"))
            c = float(match.group("c"))
            coefficients.add((a, b, c))
            od_cell = match.group("cell")
            od = wsf[od_cell].value
            if not isinstance(od, (int, float)):
                raise ValueError(f"Non-numeric OD at {cfg['file']} {od_cell}: {od!r}")
            value = a * float(od) ** 2 + b * float(od) + c
            cached = wsv[cell.coordinate].value
            if not isinstance(cached, (int, float)):
                raise ValueError(f"Missing formula cache at {cfg['file']} {cell.coordinate}")
            max_cache_delta = max(max_cache_delta, abs(value - float(cached)))
            ods.append(float(od))
            values.append(value)
            cells.append(f"{cfg['sheet']}!{cell.coordinate};OD={od_cell}")
        raw[group] = {"od": ods, "concentration": values, "source_cells": cells}
    if len(coefficients) != 1:
        raise ValueError(f"Inconsistent coefficients in {cfg['file']}: {coefficients}")
    if max_cache_delta > 1e-8:
        raise ValueError(f"Formula/cache mismatch in {cfg['file']}: {max_cache_delta}")
    a, b, c = next(iter(coefficients))
    return {
        "analyte": analyte, "source_file": cfg["file"], "source_sheet": cfg["sheet"],
        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "source_unit": cfg["source_unit"], "plot_unit": cfg["plot_unit"],
        "coefficients": {"a": a, "b": b, "c": c},
        "source_labels": source_labels, "groups": raw,
        "max_formula_cache_delta": max_cache_delta,
    }


def analyse(plate):
    groups = plate["groups"]
    result = {"groups": {}, "comparisons": {}}
    for group in GROUPS:
        values = np.asarray(groups[group]["concentration"], dtype=float)
        result["groups"][group] = {
            "n": int(values.size), "mean": float(values.mean()),
            "sd": float(values.std(ddof=1)), "values": values.tolist(),
        }
    arrays = [np.asarray(groups[g]["concentration"], dtype=float) for g in GROUPS]
    f_stat, anova_p = stats.f_oneway(*arrays)
    result["anova"] = {"F": float(f_stat), "df_between": 8, "df_within": 36, "p": float(anova_p)}
    for group_a, group_b in PAIRS:
        a = np.asarray(groups[group_a]["concentration"], dtype=float)
        b = np.asarray(groups[group_b]["concentration"], dtype=float)
        test = stats.ttest_ind(a, b, equal_var=True)
        p_raw = float(test.pvalue)
        p_adj = min(1.0, p_raw * len(PAIRS))
        result["comparisons"][f"{group_a} vs {group_b}"] = {
            "group_a": group_a, "group_b": group_b, "t": float(test.statistic),
            "p_raw": p_raw, "p_bonferroni": p_adj, "significance": significance(p_adj),
        }
    mean = {g: result["groups"][g]["mean"] for g in GROUPS}
    result["induction_percent"] = (mean["HKSA"] - mean["Control"]) / mean["Control"] * 100
    result["bai24_induction_normalized_inhibition_percent"] = (
        (mean["HKSA"] - mean["HKSA+BAI(24)"]) / (mean["HKSA"] - mean["Control"]) * 100
    )
    dose_means = [mean[f"HKSA+BAI({dose})"] for dose in (6, 12, 24, 48)]
    result["dose_means_6_12_24_48"] = dose_means
    result["dose_strictly_decreasing"] = bool(np.all(np.diff(dose_means) < 0))
    return result


def make_figure(results):
    plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["axes.linewidth"] = 0.9
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.25), dpi=300)
    plt.subplots_adjust(left=0.055, right=0.992, top=0.86, bottom=0.19, wspace=0.31)
    for panel_index, (ax, analyte, tag) in enumerate(zip(axes, PLATES, "ABC")):
        result = results[analyte]["analysis"]
        means = np.array([result["groups"][g]["mean"] for g in GROUPS])
        sds = np.array([result["groups"][g]["sd"] for g in GROUPS])
        x = np.arange(len(GROUPS))
        ax.bar(x, means, yerr=sds, capsize=3.2, width=0.68, color=COLORS,
               edgecolor="black", linewidth=0.8,
               error_kw=dict(lw=1.0, capthick=1.0), zorder=2)
        rng = np.random.default_rng(20261026 + panel_index)
        for i, group in enumerate(GROUPS):
            values = np.asarray(result["groups"][group]["values"])
            ax.scatter(np.full(values.size, x[i]) + rng.uniform(-0.11, 0.11, values.size),
                       values, s=12, color="black", alpha=0.78, linewidth=0, zorder=5)
        ax.set_xticks(x)
        ax.set_xticklabels(XLABELS, fontsize=8.4, rotation=32, ha="right", rotation_mode="anchor")
        ax.set_ylabel(f"{analyte} (pg/mL)", fontsize=11.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(direction="out", length=3.5, width=0.9, labelsize=9)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
        ymax = float(np.max(means + sds) * 1.66)
        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.65, len(GROUPS) - 0.35)
        h = ymax * 0.025
        index = {g: GROUPS.index(g) for g in GROUPS}
        ytop = max(means[index["Control"]], means[index["HKSA"]]) + sds[index["HKSA"]] + ymax * 0.035
        shown_pairs = [PAIRS[0], PAIRS[1], PAIRS[2], PAIRS[3]]
        offsets = [0.0, 0.085, 0.170, 0.255]
        for (group_a, group_b), offset in zip(shown_pairs, offsets):
            key = f"{group_a} vs {group_b}"
            p_adj = result["comparisons"][key]["p_bonferroni"]
            i, j = index[group_a], index[group_b]
            y = ytop + ymax * offset
            ax.plot([i, i, j, j], [y, y + h, y + h, y], lw=1.0, color="black", clip_on=False)
            ax.text((i + j) / 2, y + h + ymax * 0.006, significance(p_adj),
                    ha="center", va="bottom", fontsize=9.8)
        ax.text(-0.12, 1.07, tag, transform=ax.transAxes, fontsize=16,
                fontweight="bold", va="top")
    fig.savefig(PNG, dpi=300, bbox_inches="tight", pad_inches=0.12, facecolor="white")
    fig.savefig(PDF, bbox_inches="tight", pad_inches=0.12, facecolor="white")
    plt.close(fig)


def style_table_header(ws, row, start_col, end_col):
    fill = PatternFill("solid", fgColor="1F4E78")
    for col in range(start_col, end_col + 1):
        cell = ws.cell(row=row, column=col)
        cell.font = Font(name="Arial", size=10, bold=True, color="FFFFFF")
        cell.fill = fill
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)


def add_analyte_sheet(wb, analyte, plate):
    ws = wb.create_sheet(analyte)
    ws.sheet_view.showGridLines = False
    ws.merge_cells("A1:J1")
    ws["A1"] = f"{analyte} — raw OD450 and formula-driven concentrations"
    ws["A1"].font = Font(name="Arial", size=14, bold=True, color="1F4E78")
    ws["A1"].alignment = Alignment(horizontal="left")
    coef = plate["coefficients"]
    ws["A2"], ws["B2"] = "a", coef["a"]
    ws["D2"], ws["E2"] = "b", coef["b"]
    ws["G2"], ws["H2"] = "c", coef["c"]
    ws["A3"], ws["B3"] = "Source unit", plate["source_unit"]
    ws["D3"], ws["E3"] = "Figure unit", plate["plot_unit"]
    ws["G3"], ws["H3"] = "Replicates", 5
    ws["A4"], ws["B4"] = "Source file", f"source/{plate['source_file']} | sheet={plate['source_sheet']}"
    ws.merge_cells("B4:J4")
    ws["A5"] = "Concentration formula"
    ws["B5"] = "C = a × OD450² + b × OD450 + c; every concentration below is formula-driven"
    ws.merge_cells("B5:J5")
    for row in (6, 13):
        ws.cell(row=row, column=1, value="Replicate")
        for col, group in enumerate(GROUPS, start=2):
            ws.cell(row=row, column=col, value=group)
        style_table_header(ws, row, 1, 10)
    ws["A6"] = "OD450 / Replicate"
    ws["A13"] = "Concentration (pg/mL) / Replicate"
    for rep in range(5):
        ws.cell(row=7 + rep, column=1, value=rep + 1)
        ws.cell(row=14 + rep, column=1, value=rep + 1)
        for gi, group in enumerate(GROUPS, start=2):
            od = plate["groups"][group]["od"][rep]
            ws.cell(row=7 + rep, column=gi, value=od)
            od_cell = f"{get_column_letter(gi)}{7 + rep}"
            ws.cell(row=14 + rep, column=gi,
                    value=f"=$B$2*{od_cell}^2+$E$2*{od_cell}+$H$2")
    ws["A20"] = "Mean (pg/mL)"
    ws["A21"] = "SD"
    ws["A22"] = "n"
    for col in range(2, 11):
        letter = get_column_letter(col)
        ws.cell(row=20, column=col, value=f"=AVERAGE({letter}14:{letter}18)")
        ws.cell(row=21, column=col, value=f"=STDEV({letter}14:{letter}18)")
        ws.cell(row=22, column=col, value=f"=COUNT({letter}14:{letter}18)")
    for row in (20, 21, 22):
        ws.cell(row=row, column=1).font = Font(name="Arial", bold=True)
        for col in range(2, 11):
            ws.cell(row=row, column=col).number_format = "0.000"
    ws["A25"] = "Planned comparisons (two-sided, equal-variance t tests; Bonferroni family = 5)"
    ws.merge_cells("A25:F25")
    ws["A25"].font = Font(name="Arial", bold=True, color="1F4E78")
    headers = ["Comparison", "Group A", "Group B", "Raw P", "Bonferroni P", "Significance"]
    for col, header in enumerate(headers, start=1):
        ws.cell(row=26, column=col, value=header)
    style_table_header(ws, 26, 1, 6)
    for row, (group_a, group_b) in enumerate(PAIRS, start=27):
        col_a = get_column_letter(2 + GROUPS.index(group_a))
        col_b = get_column_letter(2 + GROUPS.index(group_b))
        ws.cell(row=row, column=1, value=f"{group_a} vs {group_b}")
        ws.cell(row=row, column=2, value=group_a)
        ws.cell(row=row, column=3, value=group_b)
        ws.cell(row=row, column=4, value=f"=TTEST({col_a}14:{col_a}18,{col_b}14:{col_b}18,2,2)")
        ws.cell(row=row, column=5, value=f"=MIN(1,D{row}*5)")
        ws.cell(row=row, column=6,
                value=f'=IF(E{row}<0.001,"***",IF(E{row}<0.01,"**",IF(E{row}<0.05,"*","ns")))')
        ws.cell(row=row, column=4).number_format = "0.000000"
        ws.cell(row=row, column=5).number_format = "0.000000"
    ws["H25"] = "One-way ANOVA across all 9 groups"
    ws["H25"].font = Font(name="Arial", bold=True, color="1F4E78")
    anova_rows = [
        (26, "SS_between", "=SUMPRODUCT(B22:J22,B20:J20,B20:J20)-SUMPRODUCT(B22:J22,B20:J20)^2/SUM(B22:J22)"),
        (27, "SS_within", "=SUMPRODUCT(B22:J22-1,B21:J21,B21:J21)"),
        (28, "df_between", "=COUNT(B20:J20)-1"),
        (29, "df_within", "=SUM(B22:J22)-COUNT(B20:J20)"),
        (30, "MS_between", "=I26/I28"),
        (31, "MS_within", "=I27/I29"),
        (32, "F", "=I30/I31"),
        (33, "P", "=FDIST(I32,I28,I29)"),
    ]
    for row, label, formula in anova_rows:
        ws.cell(row=row, column=8, value=label)
        ws.cell(row=row, column=9, value=formula)
        ws.cell(row=row, column=9).number_format = "0.000000"
    chart = BarChart()
    chart.type = "col"
    chart.style = 10
    chart.title = f"{analyte}: mean concentration by group"
    chart.y_axis.title = "Concentration (pg/mL)"
    chart.x_axis.title = "Experimental group"
    data = Reference(ws, min_col=1, max_col=10, min_row=20, max_row=20)
    categories = Reference(ws, min_col=2, max_col=10, min_row=6, max_row=6)
    chart.add_data(data, titles_from_data=True, from_rows=True)
    chart.set_categories(categories)
    chart.height = 9.5
    chart.width = 25.5
    chart.legend = None
    chart.dLbls = DataLabelList()
    chart.dLbls.showVal = True
    ws.add_chart(chart, "A36")
    ws.freeze_panes = "B7"
    ws.column_dimensions["A"].width = 33
    for col in range(2, 11):
        ws.column_dimensions[get_column_letter(col)].width = 23
    for col in ("H", "I"):
        ws.column_dimensions[col].width = 19
    for row in range(1, 60):
        for cell in ws[row]:
            cell.font = cell.font.copy(name="Arial", size=cell.font.sz or 10)
            cell.alignment = cell.alignment.copy(vertical="center")
    ws.row_dimensions[6].height = 45
    ws.row_dimensions[13].height = 45
    ws.auto_filter.ref = "A26:F31"
    return ws


def make_workbook(results):
    wb = Workbook()
    readme = wb.active
    readme.title = "README"
    readme.sheet_view.showGridLines = False
    readme["A1"] = "ELISA three-cytokine analysis — latest 9-group dataset"
    readme["A1"].font = Font(name="Arial", size=16, bold=True, color="1F4E78")
    readme.merge_cells("A1:F1")
    notes = [
        "Source: IL-1b.xlsx, IL-6.xlsx and TNF-a.xlsx copied from the user-specified folder on 2026-10-26; source workbooks were not modified.",
        "Concentrations were independently reconstructed from each source cell's quadratic formula and checked against its saved cached value (maximum absolute difference is reported below).",
        "All workbook concentration, mean, SD, ANOVA and planned-comparison cells are formula-driven. Python calculations provide an independent cross-check for the figure.",
        "Groups: Control; IKK-16 without HKSA; HKSA; HKSA+IKK-16; HKSA+baicalein 6/12/24/48 μM; HKSA+baicalein 24 μM+IKK-16. BAI = baicalein.",
        "Statistics: mean ± sample SD, n = 5; one-way ANOVA across 9 groups; five planned two-sided equal-variance t tests with Bonferroni adjustment (raw P × 5).",
        "Significance: *P<0.05, **P<0.01, ***P<0.001, ns = not significant. The four mechanistic brackets are shown in the figure; Control vs IKK-16 is retained in the workbook.",
        "Unit note: the IL-1β source sheet still labels its standard curve ng/ml. Following the previously confirmed manuscript convention, the figure uses pg/mL without rescaling numerical values. Verify the kit/reporting unit before submission.",
        "Figure: bars show mean ± SD; black dots show all five replicate values. Full-resolution PNG and vector PDF accompany this workbook.",
    ]
    for row, note in enumerate(notes, start=3):
        readme.cell(row=row, column=1, value=f"{row - 2}. {note}")
        readme.merge_cells(start_row=row, start_column=1, end_row=row, end_column=6)
        readme.cell(row=row, column=1).alignment = Alignment(wrap_text=True, vertical="top")
        readme.cell(row=row, column=1).font = Font(name="Arial", size=10)
        readme.row_dimensions[row].height = 32 if row != 9 else 44
    readme["A13"] = "Source integrity and formula-cache cross-check"
    readme["A13"].font = Font(name="Arial", size=11, bold=True, color="1F4E78")
    for col, header in enumerate(["Analyte", "Source file", "SHA-256", "Max |recalculated − cached|"], start=1):
        readme.cell(row=14, column=col, value=header)
    style_table_header(readme, 14, 1, 4)
    for row, analyte in enumerate(PLATES, start=15):
        plate = results[analyte]["plate"]
        readme.cell(row=row, column=1, value=analyte)
        readme.cell(row=row, column=2, value=plate["source_file"])
        readme.cell(row=row, column=3, value=plate["source_sha256"])
        readme.cell(row=row, column=4, value=plate["max_formula_cache_delta"])
        readme.cell(row=row, column=4).number_format = "0.00E+00"
    readme.column_dimensions["A"].width = 28
    readme.column_dimensions["B"].width = 25
    readme.column_dimensions["C"].width = 70
    readme.column_dimensions["D"].width = 30
    readme.column_dimensions["E"].width = 18
    readme.column_dimensions["F"].width = 18
    figure_ws = wb.create_sheet("Figure", 1)
    figure_ws.sheet_view.showGridLines = False
    figure_ws["A1"] = "Figure 24 replacement candidate — ELISA 9-group panel"
    figure_ws["A1"].font = Font(name="Arial", size=14, bold=True, color="1F4E78")
    figure_ws["A2"] = "BAI = baicalein; numbers in parentheses are μM. '+' groups contain HKSA. Bars = mean ± SD; dots = five replicate values."
    figure_ws["A3"] = "Suggested caption: Baicalein suppresses HKSA-induced cytokine release and shows non-additivity with IKK-16. One-way ANOVA plus five planned Bonferroni-adjusted t tests; ***P<0.001; ns, not significant."
    figure_ws["A2"].alignment = Alignment(wrap_text=True)
    figure_ws["A3"].alignment = Alignment(wrap_text=True)
    image = XLImage(str(PNG))
    ratio = image.height / image.width
    image.width = 1500
    image.height = int(1500 * ratio)
    figure_ws.add_image(image, "A5")
    figure_ws.column_dimensions["A"].width = 210
    figure_ws.row_dimensions[2].height = 32
    figure_ws.row_dimensions[3].height = 45
    for analyte in PLATES:
        add_analyte_sheet(wb, analyte, results[analyte]["plate"])
    wb.calculation.fullCalcOnLoad = True
    wb.calculation.forceFullCalc = True
    wb.calculation.calcMode = "auto"
    wb.save(XLSX)


def main():
    results = {}
    for analyte, cfg in PLATES.items():
        plate = extract_plate(analyte, cfg)
        analysis = analyse(plate)
        results[analyte] = {"plate": plate, "analysis": analysis}
        print(f"\n==== {analyte} (figure unit {plate['plot_unit']}; source unit {plate['source_unit']}) ====")
        print(f"Formula C={plate['coefficients']['a']}*OD^2+{plate['coefficients']['b']}*OD+{plate['coefficients']['c']}")
        print(f"Max formula-cache delta: {plate['max_formula_cache_delta']:.3e}")
        for group in GROUPS:
            item = analysis["groups"][group]
            print(f"  {group:28s} n={item['n']} mean={item['mean']:.6f} SD={item['sd']:.6f}")
        print(f"  ANOVA F(8,36)={analysis['anova']['F']:.6f}; P={analysis['anova']['p']:.6g}")
        for pair_name, comparison in analysis["comparisons"].items():
            print(f"  {pair_name}: raw P={comparison['p_raw']:.6g}; Bonferroni P={comparison['p_bonferroni']:.6g} {comparison['significance']}")
        print(f"  HKSA induction={analysis['induction_percent']:.3f}%")
        print(f"  BAI24 induction-normalized inhibition={analysis['bai24_induction_normalized_inhibition_percent']:.3f}%")
        print(f"  Dose means 6/12/24/48={','.join(f'{x:.6f}' for x in analysis['dose_means_6_12_24_48'])}; strictly decreasing={analysis['dose_strictly_decreasing']}")
    JSON.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    make_figure(results)
    make_workbook(results)
    print(f"\nSAVED PNG: {PNG}")
    print(f"SAVED PDF: {PDF}")
    print(f"SAVED XLSX: {XLSX}")
    print(f"SAVED JSON: {JSON}")


if __name__ == "__main__":
    main()
