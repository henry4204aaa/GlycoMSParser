"""
Smoke test for Session 2 features — import parsing, score_a_min filter, export.
"""

import csv
import sys
import tempfile
from pathlib import Path

from msp_curated_library_db import CuratedLibraryDB, DuplicateEntryError
from msp_curated_library_ui import _read_csv_rows, _extract_metadata_from_method, ImportPreviewDialog

PASS = 0
FAIL = 0


def check(name: str, condition: bool, detail: str = "") -> None:
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  PASS  {name}")
    else:
        FAIL += 1
        print(f"  FAIL  {name}  {detail}")


def _write_mock_cga_tsv(path: Path) -> None:
    """Create a mock CGA pseudolabel TSV with realistic columns."""
    cols = [
        "MS2scan_no", "protonatedmass", "composition", "theoretical_mass",
        "ppm_error", "ion score", "ion hit count", "ion hits m/z",
        "score_b", "score_b_support_score", "score_b_composition_penalty",
        "score_b_unexpected_penalty", "score_b_motif_summary",
        "peaklist", "peakintensity",
    ]
    rows = [
        {
            "MS2scan_no": "1001", "protonatedmass": "1200.5",
            "composition": "H5N4F1S2", "theoretical_mass": "1200.48",
            "ppm_error": "2.1", "ion score": "0.75", "ion hit count": "8",
            "ion hits m/z": "200.1;300.2;400.3",
            "score_b": "0.82", "score_b_support_score": "0.9",
            "score_b_composition_penalty": "0.05",
            "score_b_unexpected_penalty": "0.03",
            "score_b_motif_summary": "core-fuc(1);Lewis-x(1)",
            "peaklist": "200.1;300.2;400.3;500.4;600.5",
            "peakintensity": "1000;2000;3000;500;800",
        },
        {
            "MS2scan_no": "1002", "protonatedmass": "1400.6",
            "composition": "H6N5F1", "theoretical_mass": "1400.58",
            "ppm_error": "1.5", "ion score": "0.60", "ion hit count": "6",
            "ion hits m/z": "250.1;350.2",
            "score_b": "", "score_b_support_score": "",
            "score_b_composition_penalty": "",
            "score_b_unexpected_penalty": "",
            "score_b_motif_summary": "",
            "peaklist": "250.1;350.2;450.3;550.4",
            "peakintensity": "1500;2500;3500;600",
        },
        {
            "MS2scan_no": "1003", "protonatedmass": "900.3",
            "composition": "H3N2", "theoretical_mass": "900.28",
            "ppm_error": "3.0", "ion score": "0.90", "ion hit count": "10",
            "ion hits m/z": "180.1;220.2;260.3",
            "score_b": "0.91", "score_b_support_score": "0.95",
            "score_b_composition_penalty": "0.02",
            "score_b_unexpected_penalty": "0.01",
            "score_b_motif_summary": "high-mannose(1)",
            "peaklist": "180.1;220.2;260.3;340.4;420.5",
            "peakintensity": "5000;4000;3000;2000;1000",
        },
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_mock_converted_csv(path: Path) -> None:
    """Create a mock converted CSV (no annotation)."""
    cols = ["MS2scan_no", "protonatedmass", "peaklist", "peakintensity"]
    rows = [
        {"MS2scan_no": "2001", "protonatedmass": "800.4",
         "peaklist": "100.1;200.2;300.3", "peakintensity": "500;1000;1500"},
        {"MS2scan_no": "2002", "protonatedmass": "1000.5",
         "peaklist": "150.1;250.2;350.3", "peakintensity": "600;1200;1800"},
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        db_path = tmpdir / "test_lib.db"
        db = CuratedLibraryDB(db_path)

        # --- 1. Create mock CGA TSV and read it ---
        print("\n1. Read mock CGA TSV")
        tsv_path = tmpdir / "mock_cga.tsv"
        _write_mock_cga_tsv(tsv_path)
        rows = _read_csv_rows(tsv_path)
        check("read 3 rows", len(rows) == 3, f"got {len(rows)}")
        check("has composition col", "composition" in rows[0])
        check("has score_b col", "score_b" in rows[0])

        # --- 2. Simulate CGA import (b2 path) ---
        print("\n2. Simulate CGA import (b2 row conversion)")
        meta = {"sample_id": "test_sample", "ion_mode": "negative", "derivatization": "permethylated"}

        # Use ImportPreviewDialog._row_to_entry logic directly
        preview = ImportPreviewDialog.__new__(ImportPreviewDialog)
        preview.meta = meta
        preview.import_mode = "cga"
        preview.session_user = "TST"

        entries = []
        for row in rows:
            entry = preview._row_to_entry(row)
            if entry is not None:
                entries.append(entry)

        check("converted 3 entries", len(entries) == 3, f"got {len(entries)}")
        check("entry 0 composition", entries[0]["composition"] == "H5N4F1S2")
        check("entry 0 score_a", entries[0]["score_a"] == 0.75, f"got {entries[0].get('score_a')}")
        check("entry 0 score_b", entries[0]["score_b"] == 0.82, f"got {entries[0].get('score_b')}")
        check("entry 0 annotation_source", entries[0]["annotation_source"] == "CGA-RouteA+ScoreB")
        # Entry 1 has no score_b
        check("entry 1 no score_b", entries[1]["score_b"] is None)
        check("entry 1 annotation_source", entries[1]["annotation_source"] == "CGA-RouteA")

        # --- 3. Insert into DB ---
        print("\n3. Insert CGA entries into DB")
        ids = db.insert_entries_batch(entries)
        check("inserted 3", len(ids) == 3, f"got {len(ids)}")

        # --- 4. Score A min filter ---
        print("\n4. Score A min filter")
        filtered = db.query_entries({"score_a_min": 0.7})
        check("score_a >= 0.7 returns 2", len(filtered) == 2, f"got {len(filtered)}")

        filtered2 = db.query_entries({"score_a_min": 0.85})
        check("score_a >= 0.85 returns 1", len(filtered2) == 1, f"got {len(filtered2)}")

        # --- 5. Simulate converted CSV import (b1 path) ---
        print("\n5. Simulate converted CSV import (b1)")
        csv_path = tmpdir / "mock_converted.csv"
        _write_mock_converted_csv(csv_path)
        csv_rows = _read_csv_rows(csv_path)

        preview_b1 = ImportPreviewDialog.__new__(ImportPreviewDialog)
        preview_b1.meta = {"sample_id": "converted_test", "ion_mode": "negative", "derivatization": "native"}
        preview_b1.import_mode = "converted_csv"
        preview_b1.session_user = "TST"

        b1_entries = [preview_b1._row_to_entry(r) for r in csv_rows]
        b1_entries = [e for e in b1_entries if e is not None]
        check("b1 converted 2 entries", len(b1_entries) == 2, f"got {len(b1_entries)}")
        check("b1 annotation_source", b1_entries[0]["annotation_source"] == "manual-validated")
        check("b1 empty composition", b1_entries[0]["composition"] == "")

        b1_ids = db.insert_entries_batch(b1_entries)
        check("b1 inserted 2", len(b1_ids) == 2, f"got {len(b1_ids)}")

        # --- 6. Export trainable CSV ---
        print("\n6. Export trainable CSV")
        export_path = tmpdir / "export.csv"
        all_ids = [e["entry_id"] for e in db.query_entries()]
        count = db.export_trainable_csv(all_ids, export_path)
        check("exported 5 rows", count == 5, f"got {count}")
        check("export file exists", export_path.exists())

        # Verify columns
        with open(export_path, encoding="utf-8") as f:
            reader = csv.DictReader(f)
            export_cols = reader.fieldnames
            export_rows = list(reader)
        check("has peaklist col", "peaklist" in export_cols)
        check("has peakintensity col", "peakintensity" in export_cols)
        check("has composition col", "composition" in export_cols)
        check("export row count", len(export_rows) == 5, f"got {len(export_rows)}")

        # --- 7. Method JSON metadata extraction ---
        print("\n7. Method JSON metadata extraction")
        # Test raw metadata JSON format
        meta_result = _extract_metadata_from_method({
            "Glycan Type": "N",
            "Mass Analyzer charge mode": "-",
            "Derivatization Type": "PerMe",
            "Raw filename": "test_sample.raw",
        })
        check("metadata ion_mode", meta_result.get("ion_mode") == "negative")
        check("metadata derivatization", meta_result.get("derivatization") == "PerMe")

        # Test v1 method JSON format
        meta_v1 = _extract_metadata_from_method({
            "json_type": "glycomsp.method",
            "schema_version": "1.0.0",
            "sample": {"sample_name": "zebra_brain"},
            "metadata": {"Mass Analyzer charge mode": "+", "Derivatization Type": "native"},
        })
        check("v1 sample_id", meta_v1.get("sample_id") == "zebra_brain")
        check("v1 ion_mode", meta_v1.get("ion_mode") == "positive")

        # --- 8. Duplicate detection ---
        print("\n8. Duplicate detection on re-import")
        dup_count = 0
        for entry in entries:
            entry_copy = dict(entry)
            entry_copy.pop("entry_id", None)
            try:
                db.insert_entry(entry_copy)
            except DuplicateEntryError:
                dup_count += 1
        check("all 3 duplicates caught", dup_count == 3, f"got {dup_count}")

        # --- 9. Stats after all imports ---
        print("\n9. Stats verification")
        stats = db.get_stats()
        check("total entries 5", stats["total_entries"] == 5, f"got {stats['total_entries']}")
        check("unique samples 2", stats["unique_samples"] == 2, f"got {stats['unique_samples']}")

    # Summary
    print(f"\n{'='*40}")
    print(f"Results: {PASS} passed, {FAIL} failed")
    if FAIL:
        print("IMPORT SMOKE TEST FAILED")
        sys.exit(1)
    else:
        print("IMPORT SMOKE TEST PASSED")


if __name__ == "__main__":
    main()
