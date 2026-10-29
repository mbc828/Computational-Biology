# Public source snapshot: replace path placeholders before running.
#!/usr/bin/env python3
"""Final read-only checks for the backfill, full-length results, and cleanup."""
from __future__ import annotations

import csv
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent
AUDIT = BASE / "audit"
RESULTS = BASE / "full_results"


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def load_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    required = [
        RESULTS / "FULL_STRUCTURE_ANALYSIS_20261025.md",
        RESULTS / "full_structure_metrics_per_replicate.csv",
        RESULTS / "full_structure_system_summary.csv",
        RESULTS / "full_vs_tail_descriptive_summary.csv",
        RESULTS / "full_pocket_residue_summary.csv",
        RESULTS / "full_analysis_summary.json",
        AUDIT / "INTEGRITY_REPORT_20261025.md",
        AUDIT / "full_analysis_input_validation.json",
        AUDIT / "full_tail_vs_existing_window_crosscheck.csv",
        AUDIT / "backfill_finalized.json",
        AUDIT / "staging_cleanup.json",
        BASE / "full_figures" / "fig_full_0_100ns_structure_timeseries.png",
        BASE / "full_figures" / "fig_full_pocket_contacts_rmsf_res165.png",
        BASE.parent / "TRAJECTORY_AUDIT_20261019.md",
    ]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing/empty delivery file: {path}")

    validation = load_json(AUDIT / "full_analysis_input_validation.json")
    assert len(validation) == 6 and all(row["status"] == "PASS" for row in validation)
    cross = load_csv(AUDIT / "full_tail_vs_existing_window_crosscheck.csv")
    applicable = [row for row in cross if row["same_source_trajectory"] == "True"]
    distinct = [row for row in cross if row["same_source_trajectory"] == "False"]
    assert len(applicable) == 21 and all(row["status"] == "PASS" for row in applicable)
    assert all(float(row["max_absolute_difference"]) == 0.0 for row in applicable)
    assert len(distinct) == 21 and all(row["status"] == "NOT_APPLICABLE_DISTINCT_SOURCE" for row in distinct)

    directions = {row["metric"]: row for row in load_csv(RESULTS / "full_vs_tail_descriptive_summary.csv")}
    assert directions["sc165_heavy_mindist"]["full_and_tail_same_descriptive_direction"] == "YES"
    assert int(directions["sc165_heavy_mindist"]["same_label_direction_count_not_paired"]) == 3
    assert directions["sc165_heavy_contact_occupancy"]["full_and_tail_same_descriptive_direction"] == "YES"
    assert int(directions["sc165_heavy_contact_occupancy"]["same_label_direction_count_not_paired"]) == 3
    assert directions["rmsd_ligand"]["full_and_tail_same_descriptive_direction"] == "NO"
    assert int(directions["contact_count"]["same_label_direction_count_not_paired"]) == 2

    for rep in (2, 3):
        backfill = load_json(AUDIT / f"backfill_finalized_rep{rep}.json")
        assert backfill["status"] == "PASS"
        full = Path(backfill["full_path"])
        backup = Path(backfill["window_backup"])
        assert full.is_file() and full.stat().st_size == int(backfill["full_size"])
        assert backup.is_file() and backup.stat().st_size == int(backfill["window_size"])
        assert backfill["window_backup_preserved"] is True

    cleanup = load_json(AUDIT / "staging_cleanup.json")
    assert cleanup["status"] == "PASS" and cleanup["deleted_staging_files"] == 6
    assert cleanup["freed_bytes"] == 42795801688
    assert not (BASE / "full_structure" / "_staging").exists()
    for row in cleanup["files"]:
        assert Path(row["source_retained"]).is_file()
        assert Path(row["work_trajectory_retained"]).is_file()
        assert not Path(row["staged_path"]).exists()

    print("FINAL DELIVERY CHECK: PASS")
    print("6 full-length analyses; 21 same-source checks at zero difference; 21 distinct-source N/A")
    print("rep2/rep3 full trajectories and old-window backups retained on M drive")
    print("6 staging duplicates removed; all 6 sources and 6 work trajectories retained")
    print(f"Verified {len(required)} report/table/figure/audit deliverables")


if __name__ == "__main__":
    main()
