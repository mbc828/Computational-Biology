# Public source snapshot: replace path placeholders before running.
#!/usr/bin/env python3
"""Validate six full-length WT(S1)/I165A(v2)–Coptisine trajectories.

Exact 80–100 ns tail reproduction is required only when the old and full
analyses use the same source trajectory.  S1 WT is a distinct batch from old
v2 WT and is therefore a descriptive full-length control, not an exact tail
cross-check or a paired energy replacement.  No old result is modified.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).resolve().parent
FULL = BASE / "full_structure"
OLD = BASE / "structure"
MM = {"WT": BASE / "inputs" / "mmgbsa_wt", "I165A": BASE / "inputs" / "mmgbsa"}
RESULTS = BASE / "full_results"
FIGURES = BASE / "full_figures"
AUDIT = BASE / "audit"
for directory in (RESULTS, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)

SYSTEMS = ("WT", "I165A")
FOLDER = {"WT": "wt", "I165A": "mut"}
COLORS = {"WT": "#286A9F", "I165A": "#C65D3B"}
POCKET = [21, 29, 44, 96, 98, 99, 152, 165, 166]
EXPECTED_TIMES = np.round(np.arange(0.0, 100.0001, 0.4), 6)
EXPECTED_TAIL = np.round(np.arange(80.0, 100.0001, 0.4), 6)
SERIES_NAMES = (
    "rmsd_ligand", "rmsd_backbone", "mindist", "contact_count",
    "sc165_heavy_mindist", "sc165_heavy_contact_count",
)

spec = importlib.util.spec_from_file_location(
    "contact_parser_full", BASE / "inputs" / "scripts" / "parse_contact_occupancy.py"
)
contact_parser = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = contact_parser
spec.loader.exec_module(contact_parser)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(rows: list[dict], path: Path) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def mean_sd(values) -> tuple[float, float]:
    array = np.asarray(list(values), dtype=float)
    if len(array) != 3:
        raise ValueError("Summary requires exactly three computational replicates")
    return float(array.mean()), float(array.std(ddof=1))


def read_xvg(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        values = [float(value) for value in line.split()]
        if len(values) != 2:
            raise ValueError(f"{path}: expected two columns, got {len(values)}")
        rows.append(values)
    array = np.asarray(rows, dtype=float)
    if array.ndim != 2 or not len(array) or not np.isfinite(array).all():
        raise ValueError(f"{path}: empty or non-finite data")
    return array


def validate_series(array: np.ndarray, path: Path, tail: bool = False) -> None:
    expected = EXPECTED_TAIL if tail else EXPECTED_TIMES
    if array.shape != (len(expected), 2):
        raise ValueError(f"{path}: expected {(len(expected), 2)}, got {array.shape}")
    if not np.allclose(array[:, 0], expected, atol=0.002, rtol=0):
        raise ValueError(f"{path}: wrong time grid")


def read_static_ndx(path: Path, only: set[str] | None = None) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = {}
    current_name = None
    keep = False
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            current_name = line[1:-1].strip()
            keep = only is None or current_name in only
            if keep:
                if current_name in groups:
                    raise ValueError(f"{path}: duplicate group {current_name}")
                groups[current_name] = []
            continue
        if keep and current_name is not None:
            groups[current_name].extend(int(token) for token in line.split())
    if only is not None and set(groups) != only:
        raise ValueError(f"{path}: missing groups {sorted(only - set(groups))}")
    return groups


def parse_check_log(path: Path, source: bool, system: str) -> dict:
    text = path.read_text(encoding="utf-8")
    atoms = re.search(r"# Atoms\s+(\d+)", text)
    coords = re.search(r"^Coords\s+(\d+)\s+([\d.]+)", text, re.M)
    first = re.search(r"Reading frame\s+0\s+time\s+([\d.]+)", text)
    last = re.search(r"Last frame\s+(\d+)\s+time\s+([\d.]+)", text)
    result = {
        "atoms": int(atoms.group(1)) if atoms else None,
        "frames": int(coords.group(1)) if coords else None,
        "timestep_ps": float(coords.group(2)) if coords else None,
        "first_time_ps": float(first.group(1)) if first else None,
        "last_frame_index": int(last.group(1)) if last else None,
        "last_time_ps": float(last.group(2)) if last else None,
        "fatal": "Fatal error" in text,
    }
    if source:
        frames, timestep = (2001, 50.0) if system == "WT" else (10001, 10.0)
    else:
        frames, timestep = 251, 400.0
    expected = {
        "frames": frames, "timestep_ps": timestep, "first_time_ps": 0.0, "fatal": False,
    }
    # gmx check 2026.2 prints Last frame for the long source trajectories and
    # short (51-frame) windows, but not for these 251-frame work trajectories.
    # For work trajectories the exact 251-point 0–100 ns grid is independently
    # enforced below for every XVG and for the dynamic contact index.
    if source or last is not None:
        expected.update({"last_frame_index": frames - 1, "last_time_ps": 100000.0})
    for key, value in expected.items():
        if result[key] != value:
            raise ValueError(f"{path}: {key}={result[key]!r}, expected {value!r}")
    if result["atoms"] is None:
        raise ValueError(f"{path}: missing atom count")
    result["last_frame_evidence"] = "gmx check" if last is not None else "per-series 0–100 ns time grids"
    return result


def tail_of(array: np.ndarray) -> np.ndarray:
    tail = array[array[:, 0] >= 79.999]
    if tail.shape != (51, 2) or not np.allclose(tail[:, 0], EXPECTED_TAIL, atol=0.002, rtol=0):
        raise ValueError("Cannot extract the expected 80–100 ns tail")
    return tail


def portable_path(path: str) -> str:
    value = path.replace("\\", "/")
    if len(value) > 2 and value[1] == ":":
        value = f"/mnt/{value[0].lower()}/{value[3:]}"
    return value


def main() -> None:
    validation_rows = []
    structure_rows = []
    pocket_rows = []
    all_contact_rows = []
    comparison_rows = []
    series = {}

    for system in SYSTEMS:
        print(f"[TOPOLOGY] {system}", flush=True)
        mapping = contact_parser.parse_topology_dump(FULL / FOLDER[system] / "topology_dump.txt")
        protein_residues = {
            resind: residue for resind, residue in mapping.residues.items()
            if residue.moltype_name.startswith("Protein_chain_")
        }
        if len(protein_residues) != 655:
            raise ValueError(f"{system}: expected 655 protein residues")
        residue_by_nr = defaultdict(list)
        for resind, residue in protein_residues.items():
            residue_by_nr[residue.input_resnr].append(resind)
        for resnr in POCKET:
            if len(residue_by_nr[resnr]) != 1:
                raise ValueError(f"{system}: pocket residue {resnr} is not unique")

        for rep in (1, 2, 3):
            print(f"[VALIDATE FULL] {system} rep{rep}", flush=True)
            directory = FULL / FOLDER[system] / f"rep{rep}"
            old_directory = OLD / FOLDER[system] / f"rep{rep}"
            manifest = load_json(MM[system] / f"rep{rep}" / "manifest.json")
            raw = load_json(directory / "raw_copy.json")
            source = Path(raw["source"])
            staged = Path(raw["staged_copy"])
            old_source = portable_path(manifest["original_trajectory"]["path"])
            same_source = portable_path(raw["source"]) == old_source
            if raw["system"] != FOLDER[system] or raw["rep"] != rep:
                raise ValueError(f"{system} rep{rep}: raw-copy identity mismatch")
            if not source.exists() or source.stat().st_size != raw["source_size"]:
                raise ValueError(f"{system} rep{rep}: source missing/size changed")
            if not staged.exists() or staged.stat().st_size != raw["staged_size"]:
                raise ValueError(f"{system} rep{rep}: staged copy missing/size changed")
            if raw["source_size"] != raw["staged_size"]:
                raise ValueError(f"{system} rep{rep}: source/staged sizes differ")
            if ": OK" not in (directory / "verify_staged_hash.log").read_text(encoding="utf-8"):
                raise ValueError(f"{system} rep{rep}: staged SHA256 check missing/failed")
            if system == "I165A" and rep in (2, 3):
                backfill = load_json(AUDIT / f"backfill_finalized_rep{rep}.json")
                if backfill["status"] != "PASS" or backfill["full_sha256"] != raw["source_sha256"]:
                    raise ValueError(f"{system} rep{rep}: final backfill/source hash mismatch")
                if not backfill["window_backup_preserved"] or not Path(backfill["window_backup"]).exists():
                    raise ValueError(f"{system} rep{rep}: window backup not preserved")

            source_check = parse_check_log(directory / "check_source_full.log", True, system)
            work_check = parse_check_log(directory / "check_work_full.log", False, system)
            if source_check["atoms"] != work_check["atoms"]:
                raise ValueError(f"{system} rep{rep}: source/work atom mismatch")
            source_tpr = source.parent / "md.tpr"
            source_tpr_hash = sha256(source_tpr)
            tpr_hash = sha256(directory / "md.tpr")
            if tpr_hash != source_tpr_hash:
                raise ValueError(f"{system} rep{rep}: analysis TPR is not the source-run TPR")
            if raw.get("source_tpr_sha256", source_tpr_hash) != source_tpr_hash:
                raise ValueError(f"{system} rep{rep}: raw-copy TPR lineage mismatch")
            if same_source and tpr_hash != manifest["source_tpr_sha256"]:
                raise ValueError(f"{system} rep{rep}: same-source TPR does not match accepted MM/GBSA input")
            if work_check["atoms"] != manifest["atom_counts"]["System"] or work_check["atoms"] != mapping.declared_atom_count:
                raise ValueError(f"{system} rep{rep}: work/manifest/topology atom mismatch")
            ndx_hash = sha256(directory / "analysis.ndx")
            if same_source and ndx_hash != sha256(old_directory / "analysis.ndx"):
                raise ValueError(f"{system} rep{rep}: same-source analysis index differs from old window")
            static = read_static_ndx(directory / "analysis.ndx", {"Protein", "Cop", "Backbone", "C-alpha"})
            for group in ("Protein", "Cop"):
                if len(static[group]) != manifest["atom_counts"][group]:
                    raise ValueError(f"{system} rep{rep}: {group} index size mismatch")
                for atom in static[group]:
                    atom_map = mapping.atoms[atom]
                    is_protein = atom_map.moltype_name.startswith("Protein_chain_")
                    if is_protein != (group == "Protein"):
                        raise ValueError(f"{system} rep{rep}: wrong atom in {group}")
            if len(static["C-alpha"]) != 655 or len(static["Backbone"]) != 1965:
                raise ValueError(f"{system} rep{rep}: C-alpha/backbone index mismatch")

            arrays = {}
            for name in SERIES_NAMES:
                path = directory / f"{name}.xvg"
                array = read_xvg(path)
                validate_series(array, path)
                arrays[name] = array
                series[(system, rep, name)] = array

                old_path = old_directory / f"{name}.xvg"
                old_array = read_xvg(old_path)
                validate_series(old_array, old_path, tail=True)
                new_tail = tail_of(array)
                difference = np.abs(new_tail[:, 1] - old_array[:, 1])
                tolerance = 0.0 if "contact_count" in name else 0.000002
                max_difference = float(difference.max())
                if same_source and max_difference > tolerance:
                    raise ValueError(
                        f"{system} rep{rep}/{name}: same-source tail mismatch max={max_difference}, tolerance={tolerance}"
                    )
                comparison_rows.append({
                    "system": system, "rep": rep, "metric": name,
                    "same_source_trajectory": same_source, "frames_compared": 51,
                    "max_absolute_difference": max_difference, "tolerance": tolerance,
                    "status": "PASS" if same_source else "NOT_APPLICABLE_DISTINCT_SOURCE",
                })

            rmsf = read_xvg(directory / "rmsf_ca.xvg")
            if rmsf.shape != (655, 2) or len(set(rmsf[:, 0])) != 655:
                raise ValueError(f"{system} rep{rep}: invalid full-length C-alpha RMSF")
            rmsf_by_nr = {int(round(row[0])): float(row[1]) for row in rmsf}
            if any(resnr not in rmsf_by_nr for resnr in POCKET):
                raise ValueError(f"{system} rep{rep}: pocket RMSF residue missing")

            dynamic_path = directory / "contact_residues_dynamic.ndx"
            frames = contact_parser.parse_dynamic_ndx(dynamic_path)
            if len(frames) != 251:
                raise ValueError(f"{system} rep{rep}: full dynamic contact index is not 251 frames")
            header_times = [float(value) for value in re.findall(
                r"_f\d+_t([\d.]+)", dynamic_path.read_text(encoding="utf-8")
            )]
            if len(header_times) != 251 or not np.allclose(header_times, EXPECTED_TIMES * 1000, atol=2.0, rtol=0):
                raise ValueError(f"{system} rep{rep}: full dynamic contact time-grid mismatch")
            for frame_index, atom_indices in enumerate(frames):
                for atom in atom_indices:
                    if atom not in mapping.atoms or not mapping.atoms[atom].moltype_name.startswith("Protein_chain_"):
                        raise ValueError(f"{system} rep{rep}: non-protein/unknown contact in frame {frame_index}")

            old_frames = contact_parser.parse_dynamic_ndx(old_directory / "contact_residues_dynamic.ndx")
            dynamic_equal = len(old_frames) == 51 and frames[-51:] == old_frames
            if same_source and not dynamic_equal:
                raise ValueError(f"{system} rep{rep}: same-source dynamic-contact tail differs from old window")
            dynamic_status = "PASS" if same_source else "NOT_APPLICABLE_DISTINCT_SOURCE"
            comparison_rows.append({
                "system": system, "rep": rep, "metric": "dynamic_contact_atom_sets",
                "same_source_trajectory": same_source, "frames_compared": 51,
                "max_absolute_difference": 0.0 if dynamic_equal else "",
                "tolerance": 0.0, "status": dynamic_status,
            })

            occupancy = contact_parser.build_contact_occupancy(mapping, frames)
            contact_parser.write_contact_occupancy(occupancy, directory / "contact_occupancy_residue_full.tsv")
            occupancy_by_nr = {int(row["input_resnr"]): row for row in occupancy}
            for row in occupancy:
                all_contact_rows.append({
                    "system": system, "rep": rep, "resnr": int(row["input_resnr"]),
                    "resname": row["resname"], "hit_frames": int(row["hit_frames"]),
                    "total_frames": 251, "occupancy": float(row["occupancy"]),
                })

            sc_groups = read_static_ndx(directory / "sc165_cop_heavy.ndx")
            if len(sc_groups) != 2:
                raise ValueError(f"{system} rep{rep}: SC165/Cop index must contain two groups")
            sc_atoms, cop_atoms = list(sc_groups.values())
            if not sc_atoms or not cop_atoms:
                raise ValueError(f"{system} rep{rep}: empty SC165/Cop selection")
            for atom in sc_atoms:
                atom_map = mapping.atoms[atom]
                if atom_map.input_resnr != 165 or atom_map.atom_name in {"N", "CA", "C", "O"} or atom_map.atom_name.startswith("H"):
                    raise ValueError(f"{system} rep{rep}: invalid SC165 heavy atom")
            for atom in cop_atoms:
                atom_map = mapping.atoms[atom]
                if atom_map.resname != "Cop" or atom_map.atom_name.startswith("H"):
                    raise ValueError(f"{system} rep{rep}: invalid Cop heavy atom")

            for resnr in POCKET:
                resind = residue_by_nr[resnr][0]
                residue = protein_residues[resind]
                row = occupancy_by_nr.get(resnr)
                pocket_rows.append({
                    "system": system, "rep": rep, "resnr": resnr, "resname": residue.resname,
                    "contact_hit_frames": int(row["hit_frames"]) if row else 0,
                    "total_frames": 251,
                    "contact_occupancy": float(row["occupancy"]) if row else 0.0,
                    "ca_rmsf_nm": rmsf_by_nr[resnr],
                })

            record = {
                "system": system, "rep": rep, "n_frames": 251,
                "source_frames": source_check["frames"], "source_timestep_ps": source_check["timestep_ps"],
            }
            for name in SERIES_NAMES:
                values = arrays[name][:, 1]
                tail_values = tail_of(arrays[name])[:, 1]
                record[f"{name}_full_mean"] = float(values.mean())
                record[f"{name}_full_within_sd"] = float(values.std(ddof=1))
                record[f"{name}_full_max"] = float(values.max())
                record[f"{name}_tail_mean"] = float(tail_values.mean())
            record["sc165_heavy_contact_occupancy_full"] = float(np.mean(arrays["sc165_heavy_mindist"][:, 1] <= 0.45))
            record["sc165_heavy_contact_occupancy_tail"] = float(np.mean(tail_of(arrays["sc165_heavy_mindist"])[:, 1] <= 0.45))
            record["pocket_ca_rmsf_mean_nm"] = float(np.mean([rmsf_by_nr[resnr] for resnr in POCKET]))
            record["res165_ca_rmsf_nm"] = rmsf_by_nr[165]
            structure_rows.append(record)
            validation_rows.append({
                "system": system, "rep": rep, "status": "PASS",
                "analysis_window": "0–100 ns", "work_frames": 251, "work_interval_ps": 400,
                "source": portable_path(raw["source"]), "old_window_source": old_source,
                "same_source_as_old_window": same_source,
                "source_frames": source_check["frames"], "source_interval_ps": source_check["timestep_ps"],
                "atoms": work_check["atoms"], "protein_atoms": len(static["Protein"]),
                "cop_atoms": len(static["Cop"]), "source_size": raw["source_size"],
                "source_sha256": raw["source_sha256"], "work_xtc_sha256": sha256(directory / "work.xtc"),
                "source_tpr_sha256": source_tpr_hash, "tpr_sha256": tpr_hash,
                "analysis_ndx_sha256": ndx_hash, "dynamic_contact_frames": len(frames),
                "old_80_100ns_tail_crosscheck": dynamic_status,
            })
        del mapping

    write_csv(validation_rows, AUDIT / "full_analysis_input_validation.csv")
    (AUDIT / "full_analysis_input_validation.json").write_text(
        json.dumps(validation_rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(structure_rows, RESULTS / "full_structure_metrics_per_replicate.csv")
    write_csv(pocket_rows, RESULTS / "full_pocket_residue_metrics_per_replicate.csv")
    write_csv(all_contact_rows, RESULTS / "full_all_contacting_residues_per_replicate.csv")
    write_csv(comparison_rows, AUDIT / "full_tail_vs_existing_window_crosscheck.csv")

    summary_metrics = [
        "rmsd_ligand_full_mean", "rmsd_backbone_full_mean", "mindist_full_mean",
        "contact_count_full_mean", "sc165_heavy_mindist_full_mean",
        "sc165_heavy_contact_occupancy_full", "pocket_ca_rmsf_mean_nm", "res165_ca_rmsf_nm",
    ]
    structure_summary = []
    for system in SYSTEMS:
        rows = [row for row in structure_rows if row["system"] == system]
        record = {"system": system, "n_computational_replicates": 3}
        for metric in summary_metrics:
            mean, sd = mean_sd(row[metric] for row in rows)
            record[f"{metric}_mean_of_reps"] = mean
            record[f"{metric}_between_rep_sd"] = sd
        structure_summary.append(record)
    write_csv(structure_summary, RESULTS / "full_structure_system_summary.csv")

    contrast_rows = []
    for rep in (1, 2, 3):
        wt = next(row for row in structure_rows if row["system"] == "WT" and row["rep"] == rep)
        mut = next(row for row in structure_rows if row["system"] == "I165A" and row["rep"] == rep)
        record = {
            "rep_label": rep,
            "comparison_basis": "Descriptive only: S1 WT and v2 I165A are distinct batches; labels are not pairing.",
        }
        for metric in summary_metrics:
            record[f"i165a_minus_wt_{metric}"] = mut[metric] - wt[metric]
        contrast_rows.append(record)
    write_csv(contrast_rows, RESULTS / "full_structure_replicate_label_contrasts.csv")

    pocket_summary = []
    for resnr in POCKET:
        record = {"resnr": resnr}
        for system in SYSTEMS:
            rows = [row for row in pocket_rows if row["system"] == system and row["resnr"] == resnr]
            record[f"{system}_resname"] = rows[0]["resname"]
            for metric in ("contact_occupancy", "ca_rmsf_nm"):
                mean, sd = mean_sd(row[metric] for row in rows)
                record[f"{system}_{metric}_mean"] = mean
                record[f"{system}_{metric}_between_rep_sd"] = sd
        record["i165a_minus_wt_contact_occupancy"] = record["I165A_contact_occupancy_mean"] - record["WT_contact_occupancy_mean"]
        record["i165a_minus_wt_ca_rmsf_nm"] = record["I165A_ca_rmsf_nm_mean"] - record["WT_ca_rmsf_nm_mean"]
        pocket_summary.append(record)
    write_csv(pocket_summary, RESULTS / "full_pocket_residue_summary.csv")

    applicable = [row for row in comparison_rows if row["same_source_trajectory"]]
    if len(applicable) != 21 or not all(row["status"] == "PASS" for row in applicable):
        raise ValueError("Expected 21 passing same-source tail comparisons (3 reps × 7 checks)")
    summary = {
        "scope": "WT S1 and I165A v2 Coptisine; distinct source batches; 0–100 ns; 251 frames/replicate; 400 ps interval; three computational replicates per system",
        "validation": {
            "systems": 2, "replicates": 6,
            "all_full_inputs_and_source_tpr_lineage_passed": len(validation_rows) == 6,
        },
        "old_window_crosscheck": {
            "same_source_comparisons": len(applicable),
            "same_source_all_passed": all(row["status"] == "PASS" for row in applicable),
            "distinct_source_not_applicable": sum(row["status"] == "NOT_APPLICABLE_DISTINCT_SOURCE" for row in comparison_rows),
            "rmsf_comparison": "Not applicable: full-length and 80–100 ns RMSF have different averaging windows.",
        },
        "structure": {"systems": structure_summary, "replicate_label_contrasts": contrast_rows, "pocket_residues": pocket_summary},
        "cautions": [
            "Full-length WT uses S1; old 80–100 ns WT energy/structure uses v2. They are not the same trajectories.",
            "Full-length structural analysis does not recalculate, replace, or extend the existing 80–100 ns MM/GBSA energy comparison.",
            "Replicate labels across S1 WT and v2 I165A do not establish paired statistical contrasts.",
            "The 251 time samples are correlated frames, not independent biological replicates.",
            "Between-replicate SD is descriptive (n=3); no inferential P values are calculated.",
            "RMSD is relative to each source-run TPR reference after the specified backbone fit.",
            "Whole-residue contacts do not prove side-chain-specific direct interactions or experimental binding sites.",
            "Residue-165 side-chain selections naturally contain different atoms in Ile and Ala.",
        ],
    }
    (RESULTS / "full_analysis_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 9,
        "axes.titlesize": 9, "axes.linewidth": 0.8, "xtick.labelsize": 8,
        "ytick.labelsize": 8, "legend.fontsize": 8, "figure.dpi": 150,
        "savefig.dpi": 400, "savefig.bbox": "tight", "axes.spines.top": False,
        "axes.spines.right": False,
    })

    def panel_label(ax, label):
        ax.text(-0.16, 1.06, label, transform=ax.transAxes, fontsize=12,
                fontweight="bold", va="top", ha="left")

    fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.5), constrained_layout=True)
    specs = [
        ("rmsd_ligand", "Coptisine RMSD after backbone fit", "RMSD (nm)"),
        ("rmsd_backbone", "Protein backbone RMSD", "RMSD (nm)"),
        ("mindist", "Protein–Coptisine minimum distance", "Distance (nm)"),
        ("contact_count", "Protein–Coptisine contacts < 0.45 nm", "Contact pairs"),
    ]
    for index, (name, title, ylabel) in enumerate(specs):
        ax = axes.flat[index]
        for system in SYSTEMS:
            arrays = [series[(system, rep, name)] for rep in (1, 2, 3)]
            times = arrays[0][:, 0]
            values = np.vstack([array[:, 1] for array in arrays])
            for rep_values in values:
                ax.plot(times, rep_values, color=COLORS[system], alpha=0.17, linewidth=0.55)
            mean = values.mean(axis=0)
            sd = values.std(axis=0, ddof=1)
            ax.plot(times, mean, color=COLORS[system], linewidth=1.35, label=system)
            ax.fill_between(times, mean - sd, mean + sd, color=COLORS[system], alpha=0.12, linewidth=0)
        ax.axvline(80, color="#666666", linestyle="--", linewidth=0.65)
        ax.set(xlim=(0, 100), title=title, xlabel="Time (ns)", ylabel=ylabel)
        panel_label(ax, chr(ord("A") + index))
    axes.flat[0].legend(frameon=False, ncol=2, loc="best")
    fig.text(0.5, -0.025, "Thin lines: replicates; thick line: mean; band: between-replicate SD; dashed line: 80 ns",
             ha="center", va="top", fontsize=7, color="#555555")
    fig.savefig(FIGURES / "fig_full_0_100ns_structure_timeseries.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.8), constrained_layout=True)
    x = np.arange(len(POCKET))
    width = 0.36
    labels = []
    for resnr in POCKET:
        row = next(row for row in pocket_summary if row["resnr"] == resnr)
        labels.append("I/A165" if resnr == 165 else f"{row['WT_resname'].title()}{resnr}")
    for panel_index, (metric, ylabel, title) in enumerate([
        ("contact_occupancy", "Whole-residue contact occupancy (%)", "Pocket contacts over 0–100 ns"),
        ("ca_rmsf_nm", "Cα RMSF (nm)", "Full-length pocket flexibility"),
    ]):
        ax = axes[panel_index]
        for offset, system in ((-width / 2, "WT"), (width / 2, "I165A")):
            rows = [next(row for row in pocket_summary if row["resnr"] == resnr) for resnr in POCKET]
            scale = 100 if metric == "contact_occupancy" else 1
            means = [row[f"{system}_{metric}_mean"] * scale for row in rows]
            sds = [row[f"{system}_{metric}_between_rep_sd"] * scale for row in rows]
            ax.bar(x + offset, means, width=width, yerr=sds, capsize=2,
                   color=COLORS[system], alpha=0.8, label=system, error_kw={"elinewidth": 0.7})
        ax.set_xticks(x, labels, rotation=45, ha="right")
        ax.set(ylabel=ylabel, title=title)
        panel_label(ax, chr(ord("A") + panel_index))
    axes[0].set_ylim(0, 112)
    axes[0].legend(frameon=False, ncol=2)

    ax = axes[2]
    for index, system in enumerate(SYSTEMS):
        rows = [row for row in structure_rows if row["system"] == system]
        values = np.array([row["sc165_heavy_mindist_full_mean"] for row in rows])
        mean, sd = mean_sd(values)
        ax.bar(index, mean, width=0.56, color=COLORS[system], alpha=0.58,
               edgecolor=COLORS[system], yerr=sd, capsize=4, error_kw={"elinewidth": 0.8})
        ax.scatter(index + np.array([-0.08, 0, 0.08]), values, s=25,
                   color=COLORS[system], edgecolor="white", linewidth=0.5, zorder=3)
    ax.set_xticks([0, 1], ["Ile165 (WT)", "Ala165 (I165A)"])
    ax.set(ylabel="Heavy-atom minimum distance (nm)", title="Residue 165 side chain to Coptisine")
    panel_label(ax, "C")
    fig.savefig(FIGURES / "fig_full_pocket_contacts_rmsf_res165.png")
    plt.close(fig)

    print(f"[DONE] {len(validation_rows)} full-length inputs passed validation", flush=True)
    print(f"[DONE] {len(applicable)} same-source old-window tail checks passed; "
          f"{len(comparison_rows) - len(applicable)} distinct-source comparisons marked N/A", flush=True)
    print(f"[DONE] Results: {RESULTS}", flush=True)
    print(f"[DONE] Figures: {FIGURES}", flush=True)


if __name__ == "__main__":
    main()
