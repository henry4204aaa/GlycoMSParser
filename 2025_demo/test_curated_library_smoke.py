"""
Smoke test for msp_curated_library_db.py — exercises core DB operations.
"""

import sys
import tempfile
from pathlib import Path

from msp_curated_library_db import CuratedLibraryDB, DuplicateEntryError

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


def make_entry(n: int, sample: str = "sample_A", comp: str = "H5N4F1S2") -> dict:
    return {
        "sample_id": sample,
        "ms2_scan_no": 1000 + n,
        "precursor_mz": 1200.5 + n,
        "precursor_charge": -1,
        "precursor_adduct": "[M-H]-",
        "ion_mode": "negative",
        "derivatization": "permethylated",
        "composition": comp,
        "theoretical_mass": 1200.0 + n,
        "ppm_error": 2.3,
        "annotation_source": "CGA-RouteA",
        "confidence": "tentative",
        "reviewed_by": "TST",
        "created_by": "TST",
        "last_modified_by": "TST",
        "peaklist": "200.1;300.2;400.3",
        "peakintensity": "1000;2000;3000",
    }


def main() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "test_lib.db"
        db = CuratedLibraryDB(db_path)

        # --- 1. Insert 3 entries ---
        print("\n1. Insert 3 entries")
        ids = []
        for i in range(3):
            eid = db.insert_entry(make_entry(i, comp=f"H{5+i}N4F1S2"))
            ids.append(eid)
        check("inserted 3", len(ids) == 3)
        check("sequential IDs", ids == ["LIB-000001", "LIB-000002", "LIB-000003"],
              f"got {ids}")

        # --- 2. Duplicate detection ---
        print("\n2. Duplicate detection")
        dup = make_entry(0, comp="H5N4F1S2")  # same sample + scan + comp as entry 0
        try:
            db.insert_entry(dup)
            check("duplicate raises", False, "no exception raised")
        except DuplicateEntryError:
            check("duplicate raises", True)

        # --- 3. Query with filter ---
        print("\n3. Query with filter")
        all_entries = db.query_entries()
        check("query all returns 3", len(all_entries) == 3, f"got {len(all_entries)}")

        filtered = db.query_entries({"composition": "H6"})
        check("filter composition H6", len(filtered) == 1, f"got {len(filtered)}")
        if filtered:
            check("correct entry", filtered[0]["entry_id"] == "LIB-000002")

        # --- 4. Update an entry ---
        print("\n4. Update an entry")
        db.update_entry("LIB-000001", {
            "confidence": "confirmed",
            "last_modified_by": "HCT",
            "diagnostic_notes": "Verified by expert",
        })
        updated = db.get_entry("LIB-000001")
        check("confidence updated", updated["confidence"] == "confirmed",
              f"got {updated['confidence']}")
        check("notes updated", updated["diagnostic_notes"] == "Verified by expert")

        # --- 5. Delete an entry ---
        print("\n5. Delete an entry")
        db.delete_entry("LIB-000003")
        remaining = db.query_entries()
        check("delete leaves 2", len(remaining) == 2, f"got {len(remaining)}")
        check("deleted entry gone", db.get_entry("LIB-000003") is None)

        # --- 6. Stats ---
        print("\n6. Stats")
        stats = db.get_stats()
        check("total 2", stats["total_entries"] == 2, f"got {stats['total_entries']}")
        check("unique samples 1", stats["unique_samples"] == 1)
        check("unique comps 2", stats["unique_compositions"] == 2)
        check("confirmed count 1",
              stats["confidence_counts"].get("confirmed") == 1,
              f"got {stats['confidence_counts']}")

        # --- 7. Batch insert ---
        print("\n7. Batch insert (with one duplicate)")
        batch = [
            make_entry(10, sample="sample_B", comp="H3N2"),
            make_entry(0, comp="H5N4F1S2"),  # duplicate
            make_entry(11, sample="sample_B", comp="H4N3"),
        ]
        batch_ids = db.insert_entries_batch(batch)
        check("batch inserts 2 of 3", len(batch_ids) == 2, f"got {len(batch_ids)}")

        # --- 8. Export trainable CSV ---
        print("\n8. Export trainable CSV")
        csv_path = Path(tmpdir) / "export.csv"
        count = db.export_trainable_csv(["LIB-000001", "LIB-000002"], csv_path)
        check("exported 2 rows", count == 2, f"got {count}")
        check("csv file exists", csv_path.exists())

        # --- 9. Unique values ---
        print("\n9. Unique values for dropdowns")
        samples = db.get_unique_values("sample_id")
        check("unique samples list", set(samples) == {"sample_A", "sample_B"},
              f"got {samples}")

    # Summary
    print(f"\n{'='*40}")
    print(f"Results: {PASS} passed, {FAIL} failed")
    if FAIL:
        print("SMOKE TEST FAILED")
        sys.exit(1)
    else:
        print("SMOKE TEST PASSED")


if __name__ == "__main__":
    main()
