"""Smoke tests for the v0.4 delta of the GlycoMSP Curated Spectrum Library.

Covers the four new items introduced by v0.4:

1. Import preview filter helper (``filter_import_rows``) — Item 1.
2. Batch delete + batch set-confidence DB methods — Item 2.
3. Matplotlib graceful fallback flag (``_HAS_MATPLOTLIB``) — Item 3.
4. Ion-list fragment-name resolution with priority-ranked fallback, plus
   ``lookup_fragment_names`` — Item 4.

Designed to be run as a plain script (``python3 test_curated_library_v04.py``)
without pytest — the v0.1..v0.3 smoke tests in the archive use the same
PASS/FAIL counter format so that new checks integrate cleanly.
"""

from __future__ import annotations

import os
import sys
import tempfile
import time
from pathlib import Path

# Ensure the repo root is on sys.path when the script is invoked from
# elsewhere.
_REPO = Path(__file__).resolve().parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from msp_curated_library_db import (  # noqa: E402  (import after sys.path tweak)
    CuratedLibraryDB,
    _load_ion_df,
    filter_import_rows,
    lookup_fragment_names,
)


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


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def make_entry(n: int, comp: str = "H5N4F1S2", confidence: str = "tentative") -> dict:
    """Build a minimal library entry for insertion tests."""
    return {
        "sample_id": f"sample_{n % 2}",
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
        "confidence": confidence,
        "reviewed_by": "TST",
        "created_by": "TST",
        "last_modified_by": "TST",
        "peaklist": "200.1;300.2;400.3",
        "peakintensity": "1000;2000;3000",
    }


def write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            f.write(",".join(row) + "\n")


# ---------------------------------------------------------------------------
# Item 1 — filter_import_rows
# ---------------------------------------------------------------------------


def test_item1_filter() -> None:
    print("\n[Item 1] Import preview filter helper")

    rows = [
        {"MS2scan_no": "1", "composition": "H5N4",   "ion score": "0.90", "score_b": "0.80", "ppm_error": "1.5",  "sample_id": "A"},
        {"MS2scan_no": "2", "composition": "",       "ion score": "0.50", "score_b": "0.40", "ppm_error": "3.0",  "sample_id": "A"},
        {"MS2scan_no": "3", "composition": "H5N4F1", "ion score": "0.20", "score_b": "0.10", "ppm_error": "12.0", "sample_id": "B"},
        {"MS2scan_no": "4", "composition": "  ",     "ion score": "0.75", "score_b": "0.60", "ppm_error": "0.5",  "sample_id": "A"},
    ]

    # Composition-present filter should drop rows 2 and 4 (blank or whitespace).
    keep = filter_import_rows(rows, composition_present=True)
    check(
        "composition_present keeps non-empty composition rows",
        keep == [0, 2],
        f"got {keep}",
    )

    # Score B >= 0.5 keeps rows 1 and 4.
    keep = filter_import_rows(rows, score_b_min=0.5)
    check(
        "score_b_min keeps qualifying rows",
        keep == [0, 3],
        f"got {keep}",
    )

    # Score A >= 0.6 keeps rows 1 and 4.
    keep = filter_import_rows(rows, score_a_min=0.6)
    check("score_a_min uses 'ion score' column", keep == [0, 3], f"got {keep}")

    # ppm max 5 keeps rows 1, 2, 4.
    keep = filter_import_rows(rows, ppm_max=5.0)
    check("ppm_max filters by abs ppm_error", keep == [0, 1, 3], f"got {keep}")

    # Sample dropdown restricts to "A".
    keep = filter_import_rows(rows, sample="A")
    check("sample filter restricts by sample_id", keep == [0, 1, 3], f"got {keep}")

    # Combined filter: composition present + score_b >= 0.5 → row 1 only.
    keep = filter_import_rows(rows, composition_present=True, score_b_min=0.5)
    check(
        "combined filters AND together",
        keep == [0],
        f"got {keep}",
    )

    # Filter referencing a column absent from every row passes through as
    # no-op (does not raise). Build rows without score_b and call with
    # score_b_min set.
    rows_no_sb = [{"MS2scan_no": "1", "composition": "H5N4"}]
    raised = False
    try:
        keep = filter_import_rows(rows_no_sb, score_b_min=0.5)
    except Exception as exc:
        raised = True
        keep = []
    check(
        "missing score_b column passes through without raising",
        not raised and keep == [0],
        f"raised={raised}, keep={keep}",
    )

    # Missing sample column with a sample filter should still pass rows
    # through (per spec, missing column = no-op).
    rows_no_sample = [{"MS2scan_no": "1", "composition": "H5N4"}]
    keep = filter_import_rows(rows_no_sample, sample="A")
    check(
        "missing sample column with sample filter passes through",
        keep == [0],
        f"got {keep}",
    )


# ---------------------------------------------------------------------------
# Item 2 — batch delete and batch set-confidence
# ---------------------------------------------------------------------------


def test_item2_batch_delete_and_confidence() -> None:
    print("\n[Item 2] Batch delete + batch set-confidence")

    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "lib.db"
        db = CuratedLibraryDB(db_path)

        inserted_ids = []
        for i in range(5):
            entry = make_entry(i, comp=f"H{i}N4")
            inserted_ids.append(db.insert_entry(entry))

        # Attach PDFs to first two entries.
        for i in range(2):
            src_pdf = Path(tmpdir) / f"src_{i}.pdf"
            src_pdf.write_text("dummy pdf")
            db.attach_pdf(inserted_ids[i], src_pdf)

        pdf_paths_before = [
            db.db_path.parent / "pdf_references" / f"{eid}.pdf"
            for eid in inserted_ids[:2]
        ]
        check(
            "setup: both PDFs on disk before batch delete",
            all(p.exists() for p in pdf_paths_before),
            f"{[p.exists() for p in pdf_paths_before]}",
        )

        to_delete = inserted_ids[:3]  # includes both PDF-attached ones.
        deleted = db.delete_entries_batch(to_delete)
        check("batch delete removed 3 rows", deleted == 3, f"deleted={deleted}")

        remaining = {e["entry_id"] for e in db.query_entries()}
        check(
            "batch delete left the remaining 2 rows intact",
            remaining == set(inserted_ids[3:]),
            f"remaining={remaining}",
        )
        check(
            "batch delete also removed PDFs from disk",
            not any(p.exists() for p in pdf_paths_before),
            f"{[p.exists() for p in pdf_paths_before]}",
        )

        # --- batch set-confidence ---
        new_ids = []
        for i in range(5):
            entry = make_entry(100 + i, comp=f"H{i}N2", confidence="tentative")
            new_ids.append(db.insert_entry(entry))

        target_ids = new_ids[:3]
        # Capture original last_modified_at so we can confirm it advanced.
        before = {
            e["entry_id"]: e["last_modified_at"]
            for e in db.query_entries()
            if e["entry_id"] in target_ids
        }
        # Sleep just long enough for SQLite's CURRENT_TIMESTAMP (second
        # granularity) to tick over.
        time.sleep(1.1)

        updated = db.update_entries_confidence_batch(
            target_ids, "confirmed", "REV",
        )
        check("batch set-confidence updated 3 rows", updated == 3, f"updated={updated}")

        entries_after = {
            e["entry_id"]: e for e in db.query_entries()
            if e["entry_id"] in target_ids
        }
        all_confirmed = all(
            e["confidence"] == "confirmed" for e in entries_after.values()
        )
        check("all target rows now confidence=confirmed", all_confirmed)
        all_reviewer = all(
            e["last_modified_by"] == "REV" for e in entries_after.values()
        )
        check("all target rows have last_modified_by=REV", all_reviewer)
        timestamps_advanced = all(
            entries_after[eid]["last_modified_at"] > before[eid]
            for eid in target_ids
        )
        check(
            "last_modified_at advanced via trigger",
            timestamps_advanced,
            f"before={before}, after={ {k: v['last_modified_at'] for k,v in entries_after.items()} }",
        )

        # Untouched rows should still have original confidence/reviewer.
        untouched = [e for e in db.query_entries() if e["entry_id"] in new_ids[3:]]
        check(
            "untouched rows keep confidence=tentative",
            all(e["confidence"] == "tentative" for e in untouched),
        )
        check(
            "untouched rows keep last_modified_by=TST",
            all(e["last_modified_by"] == "TST" for e in untouched),
        )


# ---------------------------------------------------------------------------
# Item 3 — matplotlib graceful fallback
# ---------------------------------------------------------------------------


def test_item3_matplotlib_fallback() -> None:
    print("\n[Item 3] Matplotlib graceful fallback")

    # The UI module may pull in tkinter at import time. On headless CI
    # boxes without a display this can raise ``TclError`` — we skip the
    # test cleanly in that case rather than flagging a false failure.
    try:
        import msp_curated_library_ui as ui_mod
    except Exception as exc:
        print(f"  SKIP  could not import msp_curated_library_ui ({exc})")
        return

    check(
        "_HAS_MATPLOTLIB flag is exposed",
        hasattr(ui_mod, "_HAS_MATPLOTLIB"),
    )

    # Simulate a matplotlib-absent environment by forcing the flag False.
    # We then call ``_on_spectrum_preview`` on a minimal dummy object —
    # we only need the early-out path that shows the messagebox to not
    # raise and not import matplotlib.
    original = ui_mod._HAS_MATPLOTLIB
    ui_mod._HAS_MATPLOTLIB = False

    # Intercept the warning popup so the test can run headless. We patch
    # ``messagebox.showwarning`` to record the call.
    import tkinter.messagebox as mb
    called = {"n": 0, "title": None, "body": None}
    orig_showwarning = mb.showwarning

    def _stub(title, msg, **kwargs):
        called["n"] += 1
        called["title"] = title
        called["body"] = msg
        return None

    mb.showwarning = _stub

    class _DummyDB:
        def get_entry(self, *_a, **_kw):
            # Returned value is irrelevant — the early-out fires first.
            return {"entry_id": "LIB-1"}

    class _DummyWindow:
        db = _DummyDB()
        _ion_list_path = ""

    # Pre-test: matplotlib must not be imported by the preview path when
    # the flag is False. Clear any stale module-level import under a
    # different name; we only verify we don't *add* one via this call.
    had_mpl_before = "matplotlib.pyplot" in sys.modules

    try:
        # Bind the unbound method to our dummy receiver.
        ui_mod.CuratedLibraryWindow._on_spectrum_preview(
            _DummyWindow(), "LIB-1",
        )
        raised = False
    except Exception as exc:
        raised = True
        print(f"  DEBUG preview raised: {exc}")
    finally:
        mb.showwarning = orig_showwarning
        ui_mod._HAS_MATPLOTLIB = original

    check("preview call with _HAS_MATPLOTLIB=False does not raise", not raised)
    check(
        "messagebox.showwarning called once with expected title",
        called["n"] == 1 and called["title"] == "Spectrum preview unavailable",
        f"called={called}",
    )
    check(
        "warning body mentions pip install matplotlib",
        called["body"] is not None and "pip install matplotlib" in called["body"],
        f"body={called['body']!r}",
    )

    had_mpl_after = "matplotlib.pyplot" in sys.modules
    check(
        "preview early-out did not newly import matplotlib.pyplot",
        had_mpl_after == had_mpl_before,
    )


# ---------------------------------------------------------------------------
# Item 4 — ion list header variants + fragment lookup with priority
# ---------------------------------------------------------------------------


def test_item4_ion_list_variants() -> None:
    print("\n[Item 4] Ion-list fragment-name header variants")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # (a) only `name`
        p = tmp / "a.csv"
        write_csv(p, ["mass", "name"], [["100.0", "alpha"], ["200.0", "beta"]])
        df = _load_ion_df(p)
        check(
            "(a) only 'name' column - fragment_name populated",
            list(df["fragment_name"]) == ["alpha", "beta"],
            f"got {list(df['fragment_name'])}",
        )

        # (b) only `glycotope_name`
        p = tmp / "b.csv"
        write_csv(
            p, ["mass", "glycotope_name"],
            [["100.0", "LacNAc"], ["200.0", "terminal HexNAc"]],
        )
        df = _load_ion_df(p)
        check(
            "(b) only 'glycotope_name' - fragment_name populated",
            list(df["fragment_name"]) == ["LacNAc", "terminal HexNAc"],
        )

        # (c) only `structural_identifier`
        p = tmp / "c.csv"
        write_csv(
            p, ["mass", "structural_identifier"],
            [["100.0", "core-1"], ["200.0", "core-2"]],
        )
        df = _load_ion_df(p)
        check(
            "(c) only 'structural_identifier' - fragment_name populated",
            list(df["fragment_name"]) == ["core-1", "core-2"],
        )

        # (d) both structural_identifier and glycotope_name → SI wins.
        p = tmp / "d.csv"
        write_csv(
            p,
            ["mass", "structural_identifier", "glycotope_name"],
            [
                ["100.0", "SI-1", "GN-1"],
                ["200.0", "SI-2", "GN-2"],
            ],
        )
        df = _load_ion_df(p)
        check(
            "(d) structural_identifier wins over glycotope_name",
            list(df["fragment_name"]) == ["SI-1", "SI-2"],
        )

        # (e) SI + glycotope_name, SI partly blank → falls through.
        p = tmp / "e.csv"
        write_csv(
            p,
            ["mass", "structural_identifier", "glycotope_name"],
            [
                ["100.0", "",    "GN-1"],
                ["200.0", "SI-2", "GN-2"],
                ["300.0", "",    ""],
            ],
        )
        df = _load_ion_df(p)
        check(
            "(e) SI blank falls back to glycotope_name, SI present wins, both blank returns empty",
            list(df["fragment_name"]) == ["GN-1", "SI-2", ""],
            f"got {list(df['fragment_name'])}",
        )

        # (f) Three or more variants, sparse.
        p = tmp / "f.csv"
        write_csv(
            p,
            [
                "mass", "structural_identifier", "glycotope",
                "fragment_name", "label",
            ],
            [
                ["100.0", "SI-1", "",     "",      "L-1"],
                ["200.0", "",     "GT-2", "",      "L-2"],
                ["300.0", "",     "",     "FN-3",  "L-3"],
                ["400.0", "",     "",     "",      "L-4"],
            ],
        )
        df = _load_ion_df(p)
        check(
            "(f) sparse multiple variants resolved by priority",
            list(df["fragment_name"]) == ["SI-1", "GT-2", "FN-3", "L-4"],
            f"got {list(df['fragment_name'])}",
        )

        # (g) No accepted variants present at all.
        p = tmp / "g.csv"
        write_csv(p, ["mass", "whatever"], [["100.0", "x"], ["200.0", "y"]])
        df = _load_ion_df(p)
        check(
            "(g) no accepted variants - fragment_name column exists and empty",
            "fragment_name" in df.columns
            and list(df["fragment_name"]) == ["", ""],
            f"columns={list(df.columns)}, values={list(df['fragment_name'])}",
        )

        # (h) Mixed-case / whitespace headers.
        p = tmp / "h.csv"
        write_csv(
            p,
            ["mass", "ION NAME", "Structural_Identifier", " glycotope_name "],
            [
                ["100.0", "ion-a", "SI-a", "GN-a"],
                ["200.0", "ion-b", "",     "GN-b"],
                ["300.0", "ion-c", "",     ""],
            ],
        )
        df = _load_ion_df(p)
        check(
            "(h) SI > GN > ion_name priority survives case/whitespace normalization",
            list(df["fragment_name"]) == ["SI-a", "GN-b", "ion-c"],
            f"got {list(df['fragment_name'])}",
        )

        # lookup_fragment_names behaviour — exact-match, near-miss, and miss.
        # Reuse the last df (contains ion-a at 100.0, GN-b at 200.0, ion-c at 300.0).
        names = lookup_fragment_names(df, [100.0, 200.4, 500.0], tolerance=0.5)
        check(
            "lookup_fragment_names exact match returns name",
            names[0] == "SI-a",
            f"got {names!r}",
        )
        check(
            "lookup_fragment_names within tolerance returns name",
            names[1] == "GN-b",
            f"got {names!r}",
        )
        check(
            "lookup_fragment_names out-of-tolerance returns empty string",
            names[2] == "",
            f"got {names!r}",
        )

        # When the ion_df has no fragment_name column at all, the helper
        # returns a list of blanks.
        import pandas as pd  # local import — only needed for this check.
        bare_df = pd.DataFrame({"mass": [100.0, 200.0]})
        names = lookup_fragment_names(bare_df, [100.0, 200.0])
        check(
            "lookup_fragment_names with no fragment_name column returns blanks",
            names == ["", ""],
            f"got {names!r}",
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    print("=" * 64)
    print("Curated Library v0.4 smoke tests")
    print("=" * 64)

    test_item1_filter()
    test_item2_batch_delete_and_confidence()
    test_item3_matplotlib_fallback()
    test_item4_ion_list_variants()

    total = PASS + FAIL
    print("\n" + "=" * 64)
    print(f"Results: {PASS}/{total} PASS, {FAIL} FAIL")
    print("=" * 64)
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
