version = "1.00"
last_update = 20260408
# I decided the spec and needed components from GlycoMSP, then let the Claude Code (opus4.6) do the vibe coding
# QA and user tests are performed manually beside cli tests, and the code review has been done by author to confirm the behaviors as expected

"""
msp_curated_library_db.py — Database manager for GlycoMSP Curated Spectrum Library.

Handles all SQLite operations: schema creation, CRUD, filtering, stats,
PDF attachment, trainable CSV export, and scoring backfill.
No GUI code — no tkinter imports.
"""

from __future__ import annotations

import csv
import os
import shutil
import sqlite3
from pathlib import Path


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

_SCHEMA_SQL = """\
PRAGMA journal_mode = WAL;
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS entries (
    -- Group A: Spectrum Identity
    entry_id            TEXT PRIMARY KEY,
    sample_id           TEXT NOT NULL,
    ms2_scan_no         TEXT NOT NULL,
    precursor_mz        REAL NOT NULL,
    precursor_charge    INTEGER NOT NULL,
    precursor_adduct    TEXT NOT NULL,

    -- Group B: Experimental Conditions
    ion_mode            TEXT NOT NULL CHECK(ion_mode IN ('positive', 'negative')),
    derivatization      TEXT NOT NULL,
    instrument          TEXT,
    collision_energy    TEXT,
    lc_column_type      TEXT,
    retention_time      REAL,

    -- Group C: Composition Assignment
    composition         TEXT NOT NULL,
    theoretical_mass    REAL NOT NULL,
    ppm_error           REAL NOT NULL,
    glytoucan_id        TEXT,

    -- Group D: Scoring
    annotation_source   TEXT NOT NULL,
    score_a             REAL,
    ion_hit_count       INTEGER,
    score_b             REAL,
    score_b_support     REAL,
    score_b_comp_penalty    REAL,
    score_b_unexp_penalty   REAL,
    score_b_motif_summary   TEXT,

    -- Group E: Human Validation
    confidence          TEXT NOT NULL DEFAULT 'tentative'
                        CHECK(confidence IN ('confirmed', 'probable', 'tentative')),
    reviewed_by         TEXT NOT NULL,
    diagnostic_notes    TEXT,
    pdf_reference       TEXT,

    -- Group F: Tracking
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    created_by          TEXT NOT NULL,
    last_modified_at    TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    last_modified_by    TEXT NOT NULL,

    -- Group G: Peak Data
    peaklist            TEXT NOT NULL,
    peakintensity       TEXT NOT NULL,
    ion_hits_mz         TEXT,

    -- Constraints
    UNIQUE(sample_id, ms2_scan_no, composition)
);

CREATE TRIGGER IF NOT EXISTS update_modified_timestamp
AFTER UPDATE ON entries
BEGIN
    UPDATE entries SET last_modified_at = CURRENT_TIMESTAMP
    WHERE entry_id = NEW.entry_id;
END;

CREATE INDEX IF NOT EXISTS idx_composition ON entries(composition);
CREATE INDEX IF NOT EXISTS idx_sample ON entries(sample_id);
CREATE INDEX IF NOT EXISTS idx_confidence ON entries(confidence);
CREATE INDEX IF NOT EXISTS idx_score_b ON entries(score_b);
CREATE INDEX IF NOT EXISTS idx_score_a ON entries(score_a);
"""

# Columns in insertion order (excludes created_at / last_modified_at which are auto-set)
_INSERT_COLUMNS = [
    "entry_id", "sample_id", "ms2_scan_no", "precursor_mz", "precursor_charge",
    "precursor_adduct", "ion_mode", "derivatization", "instrument",
    "collision_energy", "lc_column_type", "retention_time",
    "composition", "theoretical_mass", "ppm_error", "glytoucan_id",
    "annotation_source", "score_a", "ion_hit_count", "score_b",
    "score_b_support", "score_b_comp_penalty", "score_b_unexp_penalty",
    "score_b_motif_summary",
    "confidence", "reviewed_by", "diagnostic_notes", "pdf_reference",
    "created_by", "last_modified_by",
    "peaklist", "peakintensity", "ion_hits_mz",
]

_ALL_COLUMNS = [
    "entry_id", "sample_id", "ms2_scan_no", "precursor_mz", "precursor_charge",
    "precursor_adduct", "ion_mode", "derivatization", "instrument",
    "collision_energy", "lc_column_type", "retention_time",
    "composition", "theoretical_mass", "ppm_error", "glytoucan_id",
    "annotation_source", "score_a", "ion_hit_count", "score_b",
    "score_b_support", "score_b_comp_penalty", "score_b_unexp_penalty",
    "score_b_motif_summary",
    "confidence", "reviewed_by", "diagnostic_notes", "pdf_reference",
    "created_at", "created_by", "last_modified_at", "last_modified_by",
    "peaklist", "peakintensity", "ion_hits_mz",
]


class DuplicateEntryError(Exception):
    """Raised when inserting a row that violates the UNIQUE constraint."""


# ---------------------------------------------------------------------------
# Data format helpers (no GUI, no tkinter)
# ---------------------------------------------------------------------------

def _semicolons_to_list(text: str | None) -> list[float]:
    """Convert a peak string to a list of floats.

    Accepts multiple storage formats:
    - Semicolon-separated: ``"200.1;300.2;400.3"``
    - Python tuple/list string: ``"(200.1, 300.2, 400.3)"`` or ``"[200.1, 300.2]"``
    - ``None`` or empty → ``[]``

    Returns:
        List of floats, or empty list if input is empty/None.
    """
    if not text:
        return []
    s = str(text).strip()
    # Detect Python tuple/list string format: "(100.0, 200.0, ...)" or "[100.0, 200.0, ...]"
    if s.startswith("(") or s.startswith("["):
        import ast
        try:
            parsed = ast.literal_eval(s)
            if isinstance(parsed, (list, tuple)):
                return [float(x) for x in parsed]
        except (ValueError, SyntaxError):
            pass
    # Fall back to semicolon-separated
    return [float(x) for x in s.split(";") if x.strip()]


def _semicolons_to_list_string(text: str | None) -> str:
    """Convert peak string to ``"[200.1, 300.2, 400.3]"`` format.

    Accepts semicolon-separated or Python tuple/list strings.
    Output is always ``ast.literal_eval``-parseable as required by
    ``build_features_from_peaks_log10_plus1()``.
    """
    if not text:
        return "[]"
    s = str(text).strip()
    # If already a Python tuple/list string, parse and re-format as list
    if s.startswith("(") or s.startswith("["):
        import ast
        try:
            parsed = ast.literal_eval(s)
            if isinstance(parsed, (list, tuple)):
                return "[" + ", ".join(str(float(x)) for x in parsed) + "]"
        except (ValueError, SyntaxError):
            pass
    # Fall back to semicolon-separated
    parts = [x.strip() for x in s.split(";") if x.strip()]
    return "[" + ", ".join(parts) + "]"


def _load_ion_df(path: str | Path):
    """Load an ion list from CSV or Excel and return a DataFrame with a ``mass`` column.

    Supports:
    - CSV / TSV files with a ``mass`` (or ``mz``, ``ion_mz``, ``m/z``) column
    - Excel files (.xlsx / .xls) — prefers a sheet named ``ionlist``

    Args:
        path: File path to the ion list.

    Returns:
        pandas DataFrame with a single ``mass`` column of floats.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If no mass-like column is found.
        ImportError: If openpyxl is needed but missing.
    """
    import pandas as pd

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Ion list not found: {p}")

    ext = p.suffix.lower()

    if ext in (".xlsx", ".xls"):
        sheets = pd.read_excel(str(p), sheet_name=None)
        df = None
        for name in sheets:
            if name.strip().lower() in {"ionlist", "ions", "ion_list"}:
                df = sheets[name]
                break
        if df is None:
            for sheet_df in sheets.values():
                cols_l = {c.strip().lower() for c in sheet_df.columns}
                if any(c in cols_l for c in {"mass", "mz", "ion_mz", "m/z"}):
                    df = sheet_df
                    break
        if df is None:
            raise ValueError(f"No sheet with a mass-like column found in: {p}")
    else:
        sep = "\t" if ext == ".tsv" else ","
        df = pd.read_csv(str(p), sep=sep, encoding="utf-8")

    # Normalize mass column name
    col_map = {c.strip().lower(): c for c in df.columns}
    for key in ("mass", "fragmentation_mass", "mz", "ion_mz", "m/z"):
        if key in col_map:
            if col_map[key] != "mass":
                df = df.rename(columns={col_map[key]: "mass"})
            break

    if "mass" not in df.columns:
        raise ValueError(f"No mass-like column found in: {p} (columns: {list(df.columns)})")

    return df[["mass"]].dropna().astype(float)


# ---------------------------------------------------------------------------
# Scoring functions (no GUI, no tkinter)
# ---------------------------------------------------------------------------

def compute_score_a(entry: dict, ion_df, ppm: float = 10.0) -> dict:
    """Compute Score A (ion hit ratio) for a single curated library entry.

    Args:
        entry: Full entry dict from the database.
        ion_df: pandas DataFrame with ``mass`` column (pre-loaded via ``_load_ion_df``).
        ppm: PPM tolerance for ion matching (default 10.0).

    Returns:
        Dict with keys ``score_a``, ``ion_hit_count``, ``ion_hits_mz``.
    """
    import pandas as pd
    from msp_insilicomarker_withGPT import attach_ion_score_on_matched

    row_df = pd.DataFrame([{
        "MS2scan_no": entry["ms2_scan_no"],
        "peaklist": _semicolons_to_list(entry.get("peaklist")),
        "peakintensity": _semicolons_to_list(entry.get("peakintensity")),
    }])

    scored = attach_ion_score_on_matched(
        row_df, ion_df, ppm_value=ppm,
        scan_col="MS2scan_no", ion_mass_col="mass",
    )
    row = scored.iloc[0]

    score_a = row.get("ion score")
    hit_count = row.get("ion hit count", 0)
    hits_mz = row.get("ion hits m/z", "")

    return {
        "score_a": float(score_a) if score_a is not None else None,
        "ion_hit_count": int(hit_count) if hit_count is not None else 0,
        "ion_hits_mz": str(hits_mz) if hits_mz else "",
    }


def compute_score_b(
    entry: dict,
    cfg,
    ppm_tolerance: float = 20.0,
    candidate_flags: dict | None = None,
) -> dict:
    """Compute Score B (motif-based) for a single curated library entry.

    Args:
        entry: Full entry dict from the database.
        cfg: ``ScoreBConfig`` object (pre-loaded via ``load_scoreb_workbook``).
        ppm_tolerance: PPM tolerance for ion matching (default 20.0).
        candidate_flags: Optional dict of flag_id → bool for MotifPolicy.

    Returns:
        Dict with keys ``score_b``, ``score_b_support``, ``score_b_comp_penalty``,
        ``score_b_unexp_penalty``, ``score_b_motif_summary``.
    """
    import pandas as pd
    from msp_CGA_structscore import (
        parse_compact_composition_label,
        build_peak_pairs_from_row,
        extract_observed_ion_hits_from_peaks,
        SpectrumEvidence,
        evaluate_global_motif_evidence,
        evaluate_candidate_score_b,
        _build_score_b_motif_summary,
        _normalize_runtime_mode,
    )

    candidate_composition = parse_compact_composition_label(entry["composition"])

    mzs = _semicolons_to_list(entry.get("peaklist"))
    intensities = _semicolons_to_list(entry.get("peakintensity"))
    row = pd.Series({
        "peaklist": mzs,
        "peakintensity": intensities,
        "MS2scan_no": entry.get("ms2_scan_no"),
    })
    peak_pairs = build_peak_pairs_from_row(row)

    observed_hits = extract_observed_ion_hits_from_peaks(
        peak_pairs, cfg, ppm_tolerance=ppm_tolerance,
    )

    charge_mode = _normalize_runtime_mode(
        "NEG" if entry.get("ion_mode") == "negative" else "POS",
        kind="charge",
    )
    derivatization_mode = _normalize_runtime_mode(
        entry.get("derivatization", "any"),
        kind="derivatization",
    )

    spectrum = SpectrumEvidence(
        observed_hits=observed_hits,
        charge_mode=charge_mode,
        derivatization=derivatization_mode,
    )
    global_evidence = evaluate_global_motif_evidence(cfg, spectrum)
    result = evaluate_candidate_score_b(
        cfg, global_evidence,
        candidate_composition=candidate_composition,
        candidate_flags=candidate_flags,
    )
    motif_summary = _build_score_b_motif_summary(result)

    return {
        "score_b": result.final_score_b,
        "score_b_support": result.motif_support_score,
        "score_b_comp_penalty": result.composition_penalty,
        "score_b_unexp_penalty": result.unexpected_penalty,
        "score_b_motif_summary": motif_summary,
    }


def compute_scores_batch(
    entries: list[dict],
    ion_list_path: str | None = None,
    workbook_path: str | None = None,
    ppm_a: float = 10.0,
    ppm_b: float = 20.0,
    overwrite: bool = False,
    do_score_a: bool = True,
    do_score_b: bool = True,
) -> list[tuple[str, dict]]:
    """Compute scores for a batch of entries, loading resources once.

    Args:
        entries: List of entry dicts from the database.
        ion_list_path: Path to ion list CSV/XLSX (required for Score A).
        workbook_path: Path to Score B workbook (optional).
        ppm_a: PPM for Score A ion matching.
        ppm_b: PPM for Score B ion matching.
        overwrite: If False, skip entries that already have the score.
        do_score_a: Whether to compute Score A.
        do_score_b: Whether to compute Score B.

    Returns:
        List of ``(entry_id, scores_dict)`` tuples for entries that were scored.
    """
    ion_df = None
    if do_score_a and ion_list_path:
        ion_df = _load_ion_df(ion_list_path)

    cfg = None
    if do_score_b and workbook_path:
        from msp_CGA_structscore import load_scoreb_workbook
        cfg = load_scoreb_workbook(workbook_path)

    results: list[tuple[str, dict]] = []

    for entry in entries:
        eid = entry["entry_id"]
        scores: dict = {}

        # Score A
        if do_score_a and ion_df is not None:
            if overwrite or entry.get("score_a") is None:
                if entry.get("peaklist"):
                    try:
                        scores.update(compute_score_a(entry, ion_df, ppm=ppm_a))
                    except Exception as exc:
                        print(f"[score_a] failed for {eid}: {exc}")

        # Score B
        if do_score_b and cfg is not None:
            if overwrite or entry.get("score_b") is None:
                comp = entry.get("composition", "")
                if comp and entry.get("peaklist"):
                    try:
                        scores.update(compute_score_b(
                            entry, cfg, ppm_tolerance=ppm_b,
                        ))
                    except Exception as exc:
                        print(f"[score_b] failed for {eid}: {exc}")

        if scores:
            results.append((eid, scores))

    return results


# ---------------------------------------------------------------------------
# CuratedLibraryDB
# ---------------------------------------------------------------------------

class CuratedLibraryDB:
    """SQLite-backed manager for the GlycoMSP Curated Spectrum Library.

    Args:
        db_path: Path to the SQLite database file. Created if it does not exist.
    """

    def __init__(self, db_path: str | Path) -> None:
        self.db_path = Path(db_path).resolve()
        self.pdf_dir = self.db_path.parent / "pdf_references"
        self.pdf_dir.mkdir(parents=True, exist_ok=True)
        self._init_schema()

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _connect(self) -> sqlite3.Connection:
        """Return a new connection with row_factory set to sqlite3.Row."""
        conn = sqlite3.connect(str(self.db_path))
        conn.row_factory = sqlite3.Row
        return conn

    def _init_schema(self) -> None:
        """Create tables, triggers, and indexes if they don't exist."""
        with self._connect() as conn:
            conn.executescript(_SCHEMA_SQL)

    # ------------------------------------------------------------------
    # Entry ID generation
    # ------------------------------------------------------------------

    def generate_entry_id(self) -> str:
        """Generate the next sequential entry ID (e.g. ``LIB-000001``)."""
        with self._connect() as conn:
            row = conn.execute(
                "SELECT entry_id FROM entries ORDER BY entry_id DESC LIMIT 1"
            ).fetchone()
            if row is None:
                return "LIB-000001"
            last_num = int(row["entry_id"].split("-")[1])
            return f"LIB-{last_num + 1:06d}"

    # ------------------------------------------------------------------
    # CRUD
    # ------------------------------------------------------------------

    def insert_entry(self, data: dict) -> str:
        """Insert a single library entry.

        Args:
            data: Dict of column values. ``entry_id`` is auto-generated if absent.

        Returns:
            The ``entry_id`` of the inserted row.

        Raises:
            DuplicateEntryError: If (sample_id, ms2_scan_no, composition)
                already exists.
        """
        if "entry_id" not in data or not data["entry_id"]:
            data["entry_id"] = self.generate_entry_id()
        cols = [c for c in _INSERT_COLUMNS if c in data]
        placeholders = ", ".join(["?"] * len(cols))
        col_names = ", ".join(cols)
        values = [data[c] for c in cols]
        try:
            with self._connect() as conn:
                conn.execute(
                    f"INSERT INTO entries ({col_names}) VALUES ({placeholders})",
                    values,
                )
        except sqlite3.IntegrityError as exc:
            if "UNIQUE constraint" in str(exc):
                raise DuplicateEntryError(
                    f"Duplicate entry: sample_id={data.get('sample_id')}, "
                    f"ms2_scan_no={data.get('ms2_scan_no')}, "
                    f"composition={data.get('composition')}"
                ) from exc
            raise
        return data["entry_id"]

    def insert_entries_batch(self, entries: list[dict]) -> list[str]:
        """Insert multiple entries, handling duplicates per row."""
        inserted_ids: list[str] = []
        for entry in entries:
            try:
                eid = self.insert_entry(entry)
                inserted_ids.append(eid)
            except DuplicateEntryError:
                continue
        return inserted_ids

    def update_entry(self, entry_id: str, updates: dict) -> None:
        """Update specific fields of an entry."""
        if not updates:
            return
        set_clause = ", ".join(f"{col} = ?" for col in updates)
        values = list(updates.values()) + [entry_id]
        with self._connect() as conn:
            conn.execute(
                f"UPDATE entries SET {set_clause} WHERE entry_id = ?",
                values,
            )

    def delete_entry(self, entry_id: str) -> None:
        """Delete an entry and its associated PDF (if any)."""
        entry = self.get_entry(entry_id)
        if entry and entry.get("pdf_reference"):
            pdf_path = self.db_path.parent / entry["pdf_reference"]
            if pdf_path.exists():
                pdf_path.unlink()
        with self._connect() as conn:
            conn.execute("DELETE FROM entries WHERE entry_id = ?", (entry_id,))

    def get_entry(self, entry_id: str) -> dict | None:
        """Fetch a single entry by ID."""
        with self._connect() as conn:
            row = conn.execute(
                "SELECT * FROM entries WHERE entry_id = ?", (entry_id,)
            ).fetchone()
            if row is None:
                return None
            return dict(row)

    def query_entries(self, filters: dict | None = None) -> list[dict]:
        """Query entries with optional filters.

        Supported filter keys: composition, sample_id, confidence,
        score_a_min, score_b_min, reviewed_by.
        """
        clauses: list[str] = []
        params: list[object] = []
        if filters:
            if "composition" in filters and filters["composition"]:
                clauses.append("composition LIKE ?")
                params.append(f"%{filters['composition']}%")
            if "sample_id" in filters and filters["sample_id"]:
                clauses.append("sample_id = ?")
                params.append(filters["sample_id"])
            if "confidence" in filters and filters["confidence"]:
                clauses.append("confidence = ?")
                params.append(filters["confidence"])
            if "score_a_min" in filters and filters["score_a_min"] is not None:
                clauses.append("score_a >= ?")
                params.append(float(filters["score_a_min"]))
            if "score_b_min" in filters and filters["score_b_min"] is not None:
                clauses.append("score_b >= ?")
                params.append(float(filters["score_b_min"]))
            if "reviewed_by" in filters and filters["reviewed_by"]:
                clauses.append("reviewed_by = ?")
                params.append(filters["reviewed_by"])

        sql = "SELECT * FROM entries"
        if clauses:
            sql += " WHERE " + " AND ".join(clauses)
        sql += " ORDER BY entry_id"

        with self._connect() as conn:
            rows = conn.execute(sql, params).fetchall()
            return [dict(r) for r in rows]

    # ------------------------------------------------------------------
    # Stats
    # ------------------------------------------------------------------

    def get_stats(self) -> dict:
        """Return summary statistics about the library."""
        with self._connect() as conn:
            total = conn.execute("SELECT COUNT(*) FROM entries").fetchone()[0]
            samples = conn.execute(
                "SELECT COUNT(DISTINCT sample_id) FROM entries"
            ).fetchone()[0]
            compositions = conn.execute(
                "SELECT COUNT(DISTINCT composition) FROM entries"
            ).fetchone()[0]
            conf_rows = conn.execute(
                "SELECT confidence, COUNT(*) as cnt FROM entries GROUP BY confidence"
            ).fetchall()
            confidence_counts = {r["confidence"]: r["cnt"] for r in conf_rows}
        return {
            "total_entries": total,
            "unique_samples": samples,
            "unique_compositions": compositions,
            "confidence_counts": confidence_counts,
        }

    # ------------------------------------------------------------------
    # Unique values (for UI dropdowns)
    # ------------------------------------------------------------------

    def get_unique_values(self, column: str) -> list[str]:
        """Return sorted unique non-null values for a column."""
        if column not in _ALL_COLUMNS:
            raise ValueError(f"Invalid column: {column}")
        with self._connect() as conn:
            rows = conn.execute(
                f"SELECT DISTINCT {column} FROM entries "
                f"WHERE {column} IS NOT NULL ORDER BY {column}"
            ).fetchall()
            return [str(r[0]) for r in rows]

    # ------------------------------------------------------------------
    # PDF attachment
    # ------------------------------------------------------------------

    def attach_pdf(self, entry_id: str, source_path: str | Path) -> str:
        """Copy a PDF into pdf_references/ and update the entry."""
        source = Path(source_path)
        dest_name = f"{entry_id}.pdf"
        dest = self.pdf_dir / dest_name
        shutil.copy2(str(source), str(dest))
        rel_path = f"pdf_references/{dest_name}"
        self.update_entry(entry_id, {"pdf_reference": rel_path})
        return rel_path

    # ------------------------------------------------------------------
    # Export (Option A: pre-compute feature matrix)
    # ------------------------------------------------------------------

    def export_trainable_csv(
        self,
        entry_ids: list[str],
        output_path: str | Path,
        ion_list_path: str | Path,
    ) -> int:
        """Export selected entries as a wide-format ML-ready CSV.

        Pre-computes the feature matrix via ``build_features_from_peaks_log10_plus1()``
        so the output is directly usable by the ML training pipeline.

        Output columns: MS2scan_no, Structure, [m/z feature columns...],
        selection_source, protonatedmass.

        Args:
            entry_ids: List of entry IDs to export.
            output_path: Destination CSV file path.
            ion_list_path: Path to the ion list CSV/XLSX.

        Returns:
            Number of rows written.
        """
        import pandas as pd
        from ml_ng_utils import build_features_from_peaks_log10_plus1

        ion_df = _load_ion_df(ion_list_path)
        ion_masses = ion_df["mass"].astype(float).tolist()

        # Fetch entries from DB
        entries = []
        for eid in entry_ids:
            e = self.get_entry(eid)
            if e is not None:
                entries.append(e)

        if not entries:
            return 0

        # Build DataFrame with ML-parseable peak format
        rows = []
        for e in entries:
            rows.append({
                "MS2scan_no": e["ms2_scan_no"],
                "Structure": e.get("composition", ""),
                "peaklist": _semicolons_to_list_string(e.get("peaklist")),
                "peakintensity": _semicolons_to_list_string(e.get("peakintensity")),
                "selection_source": e.get("annotation_source", ""),
                "protonatedmass": e.get("precursor_mz", 0.0),
            })
        df = pd.DataFrame(rows)

        # Extract features
        feat_df = build_features_from_peaks_log10_plus1(df, ion_masses, ppm=20.0)

        # Combine identity columns + features + metadata columns
        out = pd.concat([
            df[["MS2scan_no", "Structure"]],
            feat_df,
            df[["selection_source", "protonatedmass"]],
        ], axis=1)
        out.to_csv(str(output_path), index=False)
        return len(out)
