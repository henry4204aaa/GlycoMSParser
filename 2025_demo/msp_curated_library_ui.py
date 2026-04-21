

"""
msp_curated_library_ui.py — Tkinter UI for GlycoMSP Curated Spectrum Library.

Main window (Toplevel), import dialogs, filter bar, entry table, edit dialog,
spectrum preview, session settings, export, and stats dialog.
All database operations go through CuratedLibraryDB — no direct sqlite3 usage.
"""

from __future__ import annotations

version = "1.04"
last_update = 20260421
# 20260421: 1.04 (20260421) aggregated functional update based on user feedback. Version up to v0.4 (labeled as 1.04)
# I decided the spec and needed components from GlycoMSP, then let the Claude Code (opus4.6) do the vibe coding
# QA and user tests are performed manually beside cli tests, and the code review has been done by author to confirm the behaviors as expected
# Fixed several issues - overflow full peaklists in man add (clipboard_decoder) and reading csv not parsing properly issue
# Can be called from mspfileloaderv14db.py

import configparser
import csv
import json
import os
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog, ttk
from pathlib import Path

from msp_curated_library_db import (
    CuratedLibraryDB, DuplicateEntryError,
    compute_scores_batch, _load_ion_df, lookup_fragment_names,
    filter_import_rows,
)

# Graceful matplotlib import: the module as a whole must remain usable when
# matplotlib is missing; only the spectrum preview dialog depends on it.
try:  # pragma: no cover - import-time check
    import matplotlib  # noqa: F401
    _HAS_MATPLOTLIB = True
except ImportError:  # pragma: no cover - only exercised without matplotlib
    _HAS_MATPLOTLIB = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(val, default=None) -> float | None:
    """Convert a value to float, returning *default* on failure."""
    if val is None:
        return default
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def _read_csv_rows(path: str | Path) -> list[dict]:
    """Read a CSV/TSV and return a list of dicts. Auto-detects delimiter."""
    p = Path(path)
    with open(p, encoding="utf-8", errors="replace") as f:
        sample = f.read(4096)
    first_line = sample.split("\n", 1)[0]
    dialect_sep = "\t" if "\t" in first_line else ","
    with open(p, encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter=dialect_sep)
        return list(reader)


def _load_method_json(path: str | Path) -> dict:
    """Load a method or metadata JSON and return it as a dict.

    Handles both v1 (glycomsp.method) and legacy formats.
    """
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _extract_metadata_from_method(mj: dict) -> dict:
    """Extract sample_id, ion_mode, derivatization from a method/metadata JSON.

    Supports v1 schema, legacy pseudolabel, and raw metadata JSON formats.
    """
    meta: dict = {}

    # v1 method schema
    if mj.get("json_type") == "glycomsp.method":
        sample = mj.get("sample", {})
        meta["sample_id"] = sample.get("sample_name", "")
        params = mj.get("parameters", {})
        # Check for metadata nested inside
        md = mj.get("metadata", {})
        meta["ion_mode"] = _ion_mode_from_charge(
            md.get("Mass Analyzer charge mode", "")
        )
        meta["derivatization"] = md.get("Derivatization Type", "")
        return meta

    # Legacy pseudolabel method
    if mj.get("dataset_type") == "pseudolabel":
        meta["sample_id"] = mj.get("sample_name", "")
        return meta

    # Raw metadata JSON (e.g. .raw.json)
    if "Glycan Type" in mj or "Mass Analyzer charge mode" in mj:
        meta["ion_mode"] = _ion_mode_from_charge(
            mj.get("Mass Analyzer charge mode", "")
        )
        meta["derivatization"] = mj.get("Derivatization Type", "")
        meta["sample_id"] = mj.get("Raw filename", "")
        return meta

    # Legacy MAS method with samples dict
    if "samples" in mj and isinstance(mj["samples"], dict):
        sample_names = list(mj["samples"].keys())
        if sample_names:
            meta["sample_id"] = sample_names[0]
        return meta

    return meta


def _ion_mode_from_charge(charge_str: str) -> str:
    """Convert charge mode string ('+', '-') to 'positive'/'negative'."""
    s = str(charge_str).strip()
    if s == "-" or "neg" in s.lower():
        return "negative"
    return "positive"


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------

class CuratedLibraryWindow(tk.Toplevel):
    """Main window for the Curated Spectrum Library.

    Launched as a Toplevel from the GlycoMSP main application.

    Args:
        master: Parent Tk widget.
    """

    # Treeview column definitions: (col_id, heading, width, anchor)
    _TABLE_COLS = [
        ("entry_id",    "ID",          100, "w"),
        ("sample_id",   "Sample",      140, "w"),
        ("ms2_scan_no", "Scan",         70, "center"),
        ("composition", "Composition", 140, "w"),
        ("score_a",     "ScoreA",       70, "center"),
        ("score_b",     "ScoreB",       70, "center"),
        ("confidence",  "Conf.",        90, "center"),
        ("reviewed_by", "Reviewer",     80, "center"),
        ("last_modified_at", "Modified", 100, "center"),
    ]

    def __init__(self, master: tk.Tk | tk.Toplevel | None = None) -> None:
        super().__init__(master)
        self.title("Curated Spectrum Library")
        self.geometry("1060x640")
        self.minsize(860, 500)
        self.resizable(True, True)

        self.db: CuratedLibraryDB | None = None
        self.session_user: str = ""
        self._config_path: Path | None = None  # .ini alongside .db
        self._ion_list_path: str = ""
        self._scoreb_workbook_path: str = ""

        # Prompt for session user
        self._prompt_session_user()

        # Build UI
        self._build_header()
        self._build_toolbar()
        self._build_settings_bar()
        self._build_filter_bar()
        self._build_table()
        self._build_status_bar()

    # ------------------------------------------------------------------
    # Session user
    # ------------------------------------------------------------------

    def _prompt_session_user(self) -> None:
        """Ask the user for their reviewer name/initials on first launch."""
        name = simpledialog.askstring(
            "Session User",
            "Enter your reviewer name or initials:",
            parent=self,
        )
        self.session_user = (name or "").strip() or "unknown"
        self._update_title()

    def _update_title(self) -> None:
        db_name = self.db.db_path.name if self.db else "No DB"
        self.title(
            f"Curated Spectrum Library \u2014 {db_name}  [User: {self.session_user}]"
        )

    # ------------------------------------------------------------------
    # Config persistence (.ini alongside .db)
    # ------------------------------------------------------------------

    def _load_config(self) -> None:
        """Load session settings from .ini file next to the DB."""
        if self.db is None:
            return
        self._config_path = self.db.db_path.parent / "curated_library.ini"
        cfg = configparser.ConfigParser()
        if self._config_path.exists():
            cfg.read(str(self._config_path), encoding="utf-8")
        sec = "session" if cfg.has_section("session") else None
        if sec:
            stored_user = cfg.get(sec, "user", fallback="")
            if stored_user and self.session_user == "unknown":
                self.session_user = stored_user
            self._ion_list_path = cfg.get(sec, "ion_list", fallback="")
            self._scoreb_workbook_path = cfg.get(sec, "scoreb_workbook", fallback="")
        self._refresh_settings_display()

    def _save_config(self) -> None:
        """Persist current session settings to the .ini file."""
        if self._config_path is None:
            return
        cfg = configparser.ConfigParser()
        cfg["session"] = {
            "user": self.session_user,
            "ion_list": self._ion_list_path,
            "scoreb_workbook": self._scoreb_workbook_path,
        }
        with open(self._config_path, "w", encoding="utf-8") as f:
            cfg.write(f)

    # ------------------------------------------------------------------
    # Header row (v0.4 Item 2: Stats + User indicator live here)
    # ------------------------------------------------------------------

    def _build_header(self) -> None:
        """Build the top-most header row: title, Stats button, User indicator.

        v0.4 moves ``Stats`` out of the toolbar (which now hosts the
        destructive ``Delete Selected`` button) into this header row,
        immediately to the left of the session-user indicator.
        """
        header = ttk.Frame(self)
        header.pack(fill="x", padx=8, pady=(8, 0))

        ttk.Label(
            header,
            text="Curated Spectrum Library",
            font=("TkDefaultFont", 11, "bold"),
        ).pack(side="left", padx=(4, 0))

        # Right side: User (rightmost), then Stats to its left.
        self._user_label = ttk.Label(header, text=f"User: {self.session_user}")
        self._user_label.pack(side="right", padx=(4, 4))
        ttk.Button(header, text="Stats", command=self._on_stats).pack(
            side="right", padx=4,
        )

    # ------------------------------------------------------------------
    # Toolbar
    # ------------------------------------------------------------------

    def _build_toolbar(self) -> None:
        toolbar = ttk.Frame(self)
        toolbar.pack(fill="x", padx=8, pady=(4, 4))

        ttk.Button(
            toolbar, text="Open/Create DB", command=self._on_open_create_db
        ).pack(side="left", padx=4)

        # Import dropdown
        import_btn = ttk.Menubutton(toolbar, text="Import\u25bc")
        import_menu = tk.Menu(import_btn, tearoff=False)
        import_menu.add_command(
            label="From converted CSV (fresh, no annotation)",
            command=self._on_import_converted_csv,
        )
        import_menu.add_command(
            label="From CGA method (has composition + scores)",
            command=self._on_import_cga,
        )
        import_menu.add_command(
            label="From MAS method (manual annotation)",
            command=self._on_import_mas,
        )
        import_menu.add_separator()
        import_menu.add_command(
            label="Single entry (manual)",
            command=self._on_import_single,
        )
        import_btn["menu"] = import_menu
        import_btn.pack(side="left", padx=4)

        self._export_btn = ttk.Button(
            toolbar, text="Export Trainable CSV", command=self._on_export_trainable,
        )
        self._export_btn.pack(side="left", padx=4)

        # v0.4 Item 2: batch delete. Disabled until the user selects rows.
        self._delete_selected_btn = ttk.Button(
            toolbar,
            text="Delete Selected",
            command=self._on_delete_selected,
            state="disabled",
        )
        self._delete_selected_btn.pack(side="left", padx=4)

    # ------------------------------------------------------------------
    # Settings bar (ion list + Score B workbook)
    # ------------------------------------------------------------------

    def _build_settings_bar(self) -> None:
        sbar = ttk.LabelFrame(self, text="Session Settings")
        sbar.pack(fill="x", padx=8, pady=(0, 4))

        row = ttk.Frame(sbar)
        row.pack(fill="x", padx=6, pady=4)

        ttk.Label(row, text="Ion list:").pack(side="left", padx=(0, 4))
        self._ion_list_var = tk.StringVar()
        ttk.Label(row, textvariable=self._ion_list_var, width=36, relief="sunken").pack(
            side="left", padx=(0, 4)
        )
        ttk.Button(row, text="Browse\u2026", command=self._on_browse_ion_list).pack(
            side="left", padx=2
        )
        ttk.Button(row, text="Clear", command=self._on_clear_ion_list).pack(
            side="left", padx=(0, 12)
        )

        ttk.Label(row, text="Score B workbook:").pack(side="left", padx=(0, 4))
        self._scoreb_wb_var = tk.StringVar()
        ttk.Label(row, textvariable=self._scoreb_wb_var, width=36, relief="sunken").pack(
            side="left", padx=(0, 4)
        )
        ttk.Button(row, text="Browse\u2026", command=self._on_browse_scoreb_wb).pack(
            side="left", padx=2
        )
        ttk.Button(row, text="Clear", command=self._on_clear_scoreb_wb).pack(
            side="left", padx=2
        )

    def _refresh_settings_display(self) -> None:
        name_or_none = lambda p: Path(p).name if p else "(none)"
        self._ion_list_var.set(name_or_none(self._ion_list_path))
        self._scoreb_wb_var.set(name_or_none(self._scoreb_workbook_path))
        self._user_label.configure(text=f"User: {self.session_user}")
        self._update_title()

    def _on_browse_ion_list(self) -> None:
        p = filedialog.askopenfilename(
            title="Select Ion List",
            filetypes=[("CSV/Excel", "*.csv *.xlsx"), ("All files", "*.*")],
            parent=self,
        )
        if p:
            self._ion_list_path = p
            self._refresh_settings_display()
            self._save_config()

    def _on_clear_ion_list(self) -> None:
        self._ion_list_path = ""
        self._refresh_settings_display()
        self._save_config()

    def _on_browse_scoreb_wb(self) -> None:
        p = filedialog.askopenfilename(
            title="Select Score B Workbook",
            filetypes=[("Excel", "*.xlsx *.xlsm"), ("All files", "*.*")],
            parent=self,
        )
        if p:
            self._scoreb_workbook_path = p
            self._refresh_settings_display()
            self._save_config()

    def _on_clear_scoreb_wb(self) -> None:
        self._scoreb_workbook_path = ""
        self._refresh_settings_display()
        self._save_config()

    # ------------------------------------------------------------------
    # Filter bar
    # ------------------------------------------------------------------

    def _build_filter_bar(self) -> None:
        fbar = ttk.LabelFrame(self, text="Filters")
        fbar.pack(fill="x", padx=8, pady=4)

        # Row 1
        row1 = ttk.Frame(fbar)
        row1.pack(fill="x", padx=6, pady=(4, 2))

        ttk.Label(row1, text="Composition:").pack(side="left", padx=(0, 4))
        self._filt_composition = ttk.Entry(row1, width=16)
        self._filt_composition.pack(side="left", padx=(0, 12))

        ttk.Label(row1, text="Sample:").pack(side="left", padx=(0, 4))
        self._filt_sample = ttk.Combobox(row1, width=18, state="readonly")
        self._filt_sample.pack(side="left", padx=(0, 12))

        ttk.Label(row1, text="Confidence:").pack(side="left", padx=(0, 4))
        self._filt_confidence = ttk.Combobox(
            row1, width=12, state="readonly",
            values=["", "confirmed", "probable", "tentative"],
        )
        self._filt_confidence.set("")  # BUG FIX: default to empty/all
        self._filt_confidence.pack(side="left", padx=(0, 12))

        # Row 2
        row2 = ttk.Frame(fbar)
        row2.pack(fill="x", padx=6, pady=(2, 4))

        ttk.Label(row2, text="Score A min:").pack(side="left", padx=(0, 4))
        self._filt_score_a = ttk.Entry(row2, width=8)
        self._filt_score_a.pack(side="left", padx=(0, 12))

        ttk.Label(row2, text="Score B min:").pack(side="left", padx=(0, 4))
        self._filt_score_b = ttk.Entry(row2, width=8)
        self._filt_score_b.pack(side="left", padx=(0, 12))

        ttk.Label(row2, text="Reviewer:").pack(side="left", padx=(0, 4))
        self._filt_reviewer = ttk.Combobox(row2, width=12, state="readonly")
        self._filt_reviewer.pack(side="left", padx=(0, 12))

        ttk.Button(row2, text="Apply", command=self._on_apply_filters).pack(
            side="left", padx=4
        )
        ttk.Button(row2, text="Clear", command=self._on_clear_filters).pack(
            side="left", padx=4
        )

    # ------------------------------------------------------------------
    # Entry table
    # ------------------------------------------------------------------

    def _build_table(self) -> None:
        table_frame = ttk.Frame(self)
        table_frame.pack(fill="both", expand=True, padx=8, pady=4)

        col_ids = [c[0] for c in self._TABLE_COLS]
        self.tree = ttk.Treeview(
            table_frame, columns=col_ids, show="headings", selectmode="extended",
        )
        for col_id, heading, width, anchor in self._TABLE_COLS:
            self.tree.heading(col_id, text=heading)
            self.tree.column(col_id, width=width, anchor=anchor, minwidth=50)

        # Scrollbars
        vsb = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(
            table_frame, orient="horizontal", command=self.tree.xview
        )
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        table_frame.columnconfigure(0, weight=1)
        table_frame.rowconfigure(0, weight=1)

        # Bindings
        self.tree.bind("<Double-1>", self._on_double_click)
        self.tree.bind("<Button-2>", self._on_right_click)  # macOS
        self.tree.bind("<Button-3>", self._on_right_click)  # Windows/Linux
        self.tree.bind("<<TreeviewSelect>>", self._on_table_selection_changed)

    # ------------------------------------------------------------------
    # Status bar
    # ------------------------------------------------------------------

    def _build_status_bar(self) -> None:
        self._status_var = tk.StringVar(value="No database loaded")
        status_bar = ttk.Label(
            self, textvariable=self._status_var, relief="sunken", anchor="w",
        )
        status_bar.pack(fill="x", side="bottom", padx=8, pady=(0, 4))

    def _refresh_status(self, filtered_count: int | None = None) -> None:
        if self.db is None:
            self._status_var.set("No database loaded")
            return
        stats = self.db.get_stats()
        total = stats["total_entries"]
        samples = stats["unique_samples"]
        comps = stats["unique_compositions"]
        showing = filtered_count if filtered_count is not None else total
        self._status_var.set(
            f"{total} entries | {samples} samples | {comps} compositions | "
            f"Showing: {showing}/{total}"
        )

    # ------------------------------------------------------------------
    # Table data
    # ------------------------------------------------------------------

    def _populate_table(self, entries: list[dict] | None = None) -> None:
        """Clear and repopulate the Treeview from a list of entry dicts."""
        self.tree.delete(*self.tree.get_children())
        if self.db is None:
            return
        if entries is None:
            entries = self.db.query_entries()
        col_ids = [c[0] for c in self._TABLE_COLS]
        for entry in entries:
            values = []
            for cid in col_ids:
                val = entry.get(cid)
                if val is None:
                    val = "\u2014"  # em dash
                elif cid in ("score_a", "score_b"):
                    try:
                        val = f"{float(val):.2f}"
                    except (ValueError, TypeError):
                        val = "\u2014"
                elif cid == "last_modified_at" and val:
                    val = str(val)[:10]
                values.append(val)
            self.tree.insert("", "end", iid=entry["entry_id"], values=values)
        self._refresh_status(len(entries))

    def _populate_filter_dropdowns(self) -> None:
        """Refresh filter combobox values from the current database."""
        if self.db is None:
            return
        samples = [""] + self.db.get_unique_values("sample_id")
        self._filt_sample["values"] = samples

        reviewers = [""] + self.db.get_unique_values("reviewed_by")
        self._filt_reviewer["values"] = reviewers

    # ------------------------------------------------------------------
    # Button handlers — Open/Create DB
    # ------------------------------------------------------------------

    def _on_open_create_db(self) -> None:
        """Open an existing .db or create a new one."""
        path = filedialog.asksaveasfilename(
            title="Open or create curated library database",
            defaultextension=".db",
            filetypes=[("SQLite Database", "*.db"), ("All files", "*.*")],
            confirmoverwrite=False,
            parent=self,
        )
        if not path:
            return
        self.db = CuratedLibraryDB(path)
        self._load_config()
        self._update_title()
        self._populate_filter_dropdowns()
        self._populate_table()

    # ------------------------------------------------------------------
    # Button handlers — Stats
    # ------------------------------------------------------------------

    def _on_stats(self) -> None:
        if self.db is None:
            messagebox.showinfo("Stats", "No database loaded.", parent=self)
            return
        StatsDialog(self, self.db)

    # ------------------------------------------------------------------
    # Button handlers — Filters
    # ------------------------------------------------------------------

    def _on_apply_filters(self) -> None:
        """Read filter bar values and refresh the table."""
        if self.db is None:
            return
        filters: dict = {}
        comp = self._filt_composition.get().strip()
        if comp:
            filters["composition"] = comp
        sample = self._filt_sample.get().strip()
        if sample:
            filters["sample_id"] = sample
        conf = self._filt_confidence.get().strip()
        if conf:
            filters["confidence"] = conf

        score_a_text = self._filt_score_a.get().strip()
        if score_a_text:
            try:
                filters["score_a_min"] = float(score_a_text)
            except ValueError:
                messagebox.showwarning(
                    "Invalid filter", "Score A min must be a number.", parent=self,
                )
                return

        score_b_text = self._filt_score_b.get().strip()
        if score_b_text:
            try:
                filters["score_b_min"] = float(score_b_text)
            except ValueError:
                messagebox.showwarning(
                    "Invalid filter", "Score B min must be a number.", parent=self,
                )
                return
        reviewer = self._filt_reviewer.get().strip()
        if reviewer:
            filters["reviewed_by"] = reviewer

        entries = self.db.query_entries(filters if filters else None)
        self._populate_table(entries)

    def _on_clear_filters(self) -> None:
        """Reset all filters and show all entries."""
        self._filt_composition.delete(0, "end")
        self._filt_sample.set("")
        self._filt_confidence.set("")
        self._filt_score_a.delete(0, "end")
        self._filt_score_b.delete(0, "end")
        self._filt_reviewer.set("")
        if self.db is not None:
            self._populate_table()

    # ------------------------------------------------------------------
    # Table interactions
    # ------------------------------------------------------------------

    def _selected_entry_ids(self) -> list[str]:
        return list(self.tree.selection())

    def _selected_entry_id(self) -> str | None:
        sel = self.tree.selection()
        return sel[0] if sel else None

    def _on_double_click(self, event: tk.Event) -> None:
        entry_id = self._selected_entry_id()
        if entry_id and self.db:
            EditEntryDialog(
                self, self.db, entry_id, self.session_user,
                self._on_edit_done,
                ion_list_path=self._ion_list_path,
                scoreb_wb_path=self._scoreb_workbook_path,
            )

    def _on_right_click(self, event: tk.Event) -> None:
        row_id = self.tree.identify_row(event.y)
        if row_id:
            # Only change selection if the clicked row is not already
            # selected — preserves multi-select for batch operations.
            if row_id not in self.tree.selection():
                self.tree.selection_set(row_id)
        else:
            return

        sel = list(self.tree.selection())
        multi_sel = len(sel) > 1
        delete_label = f"Delete ({len(sel)})" if multi_sel else "Delete"
        single_state = "disabled" if multi_sel else "normal"

        menu = tk.Menu(self, tearoff=False)
        menu.add_command(
            state=single_state,
            label="Edit\u2026", command=lambda: self._on_context_edit(row_id)
        )
        menu.add_command(
            label=delete_label, command=self._on_delete_selected,
        )
        menu.add_command(
            label="Set Confidence for Selected...",
            command=self._on_set_confidence_selected,
        )
        menu.add_separator()
        menu.add_command(
            state=single_state,
            label="Attach PDF\u2026",
            command=lambda: self._on_context_attach_pdf(row_id),
        )
        menu.add_command(
            state=single_state,
            label="Copy composition",
            command=lambda: self._on_context_copy_comp(row_id),
        )
        menu.add_separator()
        menu.add_command(
            label="Compute Scores for Selected",
            command=self._on_compute_scores,
        )
        menu.add_command(
            state=single_state,
            label="Spectrum Preview\u2026",
            command=lambda: self._on_spectrum_preview(row_id),
        )
        menu.tk_popup(event.x_root, event.y_root)

    def _on_context_edit(self, entry_id: str) -> None:
        if self.db:
            EditEntryDialog(
                self, self.db, entry_id, self.session_user,
                self._on_edit_done,
                ion_list_path=self._ion_list_path,
                scoreb_wb_path=self._scoreb_workbook_path,
            )

    def _on_context_delete(self, entry_id: str) -> None:
        # Legacy single-row delete helper, retained for callers that pass
        # an explicit entry_id. New UI paths go through
        # ``_on_delete_selected`` which is batch-aware.
        if self.db is None:
            return
        if messagebox.askyesno(
            "Confirm Delete",
            f"Delete {entry_id}? This cannot be undone.",
            parent=self,
        ):
            self.db.delete_entry(entry_id)
            self._populate_table()
            self._populate_filter_dropdowns()
            self._on_table_selection_changed(None)
    # 20260421 code review: Claude Code missed definition
    def _on_table_selection_changed(self, event=None) -> None:
        """Enable or disable row-scoped actions based on current Treeview selection.

        Called by:
        - <<TreeviewSelect>> binding (Tkinter passes an event object)
        - _populate_table / filter refresh paths (call with None)
        """
        has_selection = bool(self.tree.selection())
        self._delete_selected_btn.configure(state="normal" if has_selection else "disabled")

    def _on_delete_selected(self) -> None:
        """Batch-delete every currently selected entry (v0.4 Item 2).

        Confirms with the total count, then calls
        ``CuratedLibraryDB.delete_entries_batch`` which removes the rows
        and their attached PDFs under a single SQL transaction.
        """
        if self.db is None:
            return
        sel = list(self.tree.selection())
        if not sel:
            return
        n = len(sel)
        if not messagebox.askyesno(
            "Confirm Delete",
            f"Delete {n} entr{'y' if n == 1 else 'ies'}? This cannot be "
            f"undone. Associated PDFs in pdf_references/ will also be "
            f"removed.",
            parent=self,
        ):
            return
        self.db.delete_entries_batch(sel)
        self._populate_table()
        self._populate_filter_dropdowns()
        self._on_table_selection_changed(None)

    def _on_set_confidence_selected(self) -> None:
        """Batch-update ``confidence`` + ``last_modified_by`` (v0.4 Item 2)."""
        if self.db is None:
            return
        sel = list(self.tree.selection())
        if not sel:
            return

        def _apply(confidence: str, reviewer: str) -> None:
            effective_reviewer = reviewer.strip() or self.session_user
            self.db.update_entries_confidence_batch(
                sel, confidence, effective_reviewer,
            )
            self._populate_table()
            self._populate_filter_dropdowns()

        BatchConfidenceDialog(
            self, n_selected=len(sel),
            session_user=self.session_user,
            on_apply=_apply,
        )

    def _on_context_attach_pdf(self, entry_id: str) -> None:
        if self.db is None:
            return
        path = filedialog.askopenfilename(
            title="Select PDF to attach",
            filetypes=[("PDF files", "*.pdf"), ("All files", "*.*")],
            parent=self,
        )
        if path:
            self.db.attach_pdf(entry_id, path)
            self._populate_table()

    def _on_context_copy_comp(self, entry_id: str) -> None:
        if self.db is None:
            return
        entry = self.db.get_entry(entry_id)
        if entry:
            self.clipboard_clear()
            self.clipboard_append(entry["composition"])
            self.update()

    def _on_edit_done(self) -> None:
        """Callback after edit dialog saves or deletes."""
        self._populate_table()
        self._populate_filter_dropdowns()

    # ==================================================================
    # IMPORT PATHS
    # ==================================================================

    def _require_db(self) -> bool:
        if self.db is None:
            messagebox.showinfo("Import", "Open or create a database first.", parent=self)
            return False
        return True

    # ---- b1: From converted CSV (no annotation) ----

    def _on_import_converted_csv(self) -> None:
        if not self._require_db():
            return
        path = filedialog.askopenfilename(
            title="Select Converted CSV",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            parent=self,
        )
        if not path:
            return
        rows = _read_csv_rows(path)
        if not rows:
            messagebox.showwarning("Import", "CSV is empty.", parent=self)
            return

        # Collect metadata (optional method JSON, then manual fallback)
        meta = self._ask_import_metadata(default_sample=Path(path).stem)

        if meta is None:
            return  # user cancelled

        ImportPreviewDialog(
            self, self.db, rows, meta,
            import_mode="converted_csv",
            session_user=self.session_user,
            on_done=self._on_import_done,
        )

    # ---- b2: From CGA method (has composition + scores) ----

    def _on_import_cga(self) -> None:
        if not self._require_db():
            return
        path = filedialog.askopenfilename(
            title="Select CGA Pseudolabel TSV",
            filetypes=[("TSV files", "*.tsv"), ("CSV files", "*.csv"), ("All files", "*.*")],
            parent=self,
        )
        if not path:
            return
        rows = _read_csv_rows(path)
        if not rows:
            messagebox.showwarning("Import", "File is empty.", parent=self)
            return

        meta = self._ask_import_metadata(default_sample=Path(path).stem)
        if meta is None:
            return

        ImportPreviewDialog(
            self, self.db, rows, meta,
            import_mode="cga",
            session_user=self.session_user,
            on_done=self._on_import_done,
        )

    # ---- b3: From MAS method (manual annotation) ----

    def _on_import_mas(self) -> None:
        if not self._require_db():
            return
        path = filedialog.askopenfilename(
            title="Select MAS-annotated CSV",
            filetypes=[("CSV files", "*.csv"), ("TSV files", "*.tsv"), ("All files", "*.*")],
            parent=self,
        )
        if not path:
            return
        rows = _read_csv_rows(path)
        if not rows:
            messagebox.showwarning("Import", "File is empty.", parent=self)
            return

        meta = self._ask_import_metadata(default_sample=Path(path).stem)
        if meta is None:
            return

        ImportPreviewDialog(
            self, self.db, rows, meta,
            import_mode="mas",
            session_user=self.session_user,
            on_done=self._on_import_done,
        )

    # ---- b5: Single entry (manual) ----

    def _on_import_single(self) -> None:
        if not self._require_db():
            return
        NewEntryDialog(self, self.db, self.session_user, self._on_import_done)

    # ---- Metadata dialog shared by b1/b2/b3 ----

    def _ask_import_metadata(self, default_sample: str = "") -> dict | None:
        """Show a dialog to collect sample_id, ion_mode, derivatization.

        Offers an optional [Select method JSON] to auto-fill.
        Returns dict or None if cancelled.
        """
        dlg = tk.Toplevel(self)
        dlg.title("Import Metadata")
        dlg.geometry("440x240")
        dlg.resizable(False, False)
        dlg.transient(self)
        dlg.grab_set()

        result: dict = {}

        frame = ttk.Frame(dlg, padding=12)
        frame.pack(fill="both", expand=True)

        # Method JSON button
        def _load_method():
            p = filedialog.askopenfilename(
                title="Select Method JSON",
                filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
                parent=dlg,
            )
            if not p:
                return
            try:
                mj = _load_method_json(p)
                meta = _extract_metadata_from_method(mj)
                if meta.get("sample_id"):
                    sample_var.set(meta["sample_id"])
                if meta.get("ion_mode"):
                    ion_mode_var.set(meta["ion_mode"])
                if meta.get("derivatization"):
                    deriv_var.set(meta["derivatization"])
            except Exception as exc:
                messagebox.showwarning("Method JSON", f"Could not parse:\n{exc}", parent=dlg)

        ttk.Button(frame, text="Select method JSON\u2026", command=_load_method).grid(
            row=0, column=0, columnspan=2, sticky="w", pady=(0, 8),
        )

        ttk.Label(frame, text="Sample ID:").grid(row=1, column=0, sticky="w", pady=4)
        sample_var = tk.StringVar(value=default_sample)
        ttk.Entry(frame, textvariable=sample_var, width=30).grid(
            row=1, column=1, sticky="ew", padx=(8, 0), pady=4,
        )

        ttk.Label(frame, text="Ion mode:").grid(row=2, column=0, sticky="w", pady=4)
        ion_mode_var = tk.StringVar(value="positive")
        ttk.Combobox(
            frame, textvariable=ion_mode_var, state="readonly", width=12,
            values=["positive", "negative"],
        ).grid(row=2, column=1, sticky="w", padx=(8, 0), pady=4)

        ttk.Label(frame, text="Derivatization:").grid(row=3, column=0, sticky="w", pady=4)
        deriv_var = tk.StringVar(value="permethylated")
        ttk.Combobox(
            frame, textvariable=deriv_var, width=20,
            values=["permethylated", "native", "2-AB", "2-AA", "other"],
        ).grid(row=3, column=1, sticky="w", padx=(8, 0), pady=4)

        cancelled = {"v": False}

        def _ok():
            sid = sample_var.get().strip()
            if not sid:
                messagebox.showwarning("Metadata", "Sample ID is required.", parent=dlg)
                return
            result["sample_id"] = sid
            result["ion_mode"] = ion_mode_var.get().strip() or "positive"
            result["derivatization"] = deriv_var.get().strip() or "native"
            dlg.destroy()

        def _cancel():
            cancelled["v"] = True
            dlg.destroy()

        btn_row = ttk.Frame(frame)
        btn_row.grid(row=4, column=0, columnspan=2, pady=(12, 0))
        ttk.Button(btn_row, text="OK", command=_ok).pack(side="right", padx=4)
        ttk.Button(btn_row, text="Cancel", command=_cancel).pack(side="right", padx=4)

        frame.columnconfigure(1, weight=1)
        dlg.wait_window()
        if cancelled["v"] or not result:
            return None
        return result

    def _on_import_done(self) -> None:
        self._populate_table()
        self._populate_filter_dropdowns()

    # ==================================================================
    # EXPORT TRAINABLE CSV (b4)
    # ==================================================================

    def _on_export_trainable(self) -> None:
        if self.db is None:
            messagebox.showinfo("Export", "No database loaded.", parent=self)
            return

        # Get visible entry IDs from the treeview
        all_ids = list(self.tree.get_children())
        if not all_ids:
            messagebox.showinfo("Export", "No entries to export.", parent=self)
            return

        # If user has a selection, use that; otherwise use all visible
        sel = list(self.tree.selection())
        export_ids = sel if sel else all_ids
        total_visible = len(all_ids)

        # Prompt for ion list (pre-fill from session settings)
        ion_path = self._ion_list_path
        if not ion_path or not Path(ion_path).exists():
            ion_path = filedialog.askopenfilename(
                title="Select Ion List for Feature Extraction",
                filetypes=[("CSV/Excel", "*.csv *.xlsx"), ("All files", "*.*")],
                parent=self,
            )
        if not ion_path:
            messagebox.showinfo(
                "Export",
                "Ion list is required for trainable CSV export.",
                parent=self,
            )
            return
        # Update session setting if changed
        if ion_path != self._ion_list_path:
            self._ion_list_path = ion_path
            self._refresh_settings_display()
            self._save_config()

        # Confirmation dialog
        comps = set()
        for eid in export_ids:
            entry = self.db.get_entry(eid)
            if entry and entry.get("composition"):
                comps.add(entry["composition"])

        msg = (
            f"Export {len(export_ids)} of {total_visible} visible entries?\n"
            f"Unique compositions: {len(comps)}\n"
            f"Ion list: {Path(ion_path).name}\n\n"
            f"Tip: Select specific rows to narrow the export."
        )
        if not messagebox.askyesno("Export Trainable CSV", msg, parent=self):
            return

        out_path = filedialog.asksaveasfilename(
            title="Save Trainable CSV",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            parent=self,
        )
        if not out_path:
            return

        try:
            count = self.db.export_trainable_csv(export_ids, out_path, ion_path)
            messagebox.showinfo(
                "Export", f"Exported {count} entries to:\n{out_path}", parent=self,
            )
        except Exception as exc:
            messagebox.showerror(
                "Export Error", f"Export failed:\n{exc}", parent=self,
            )

    # ==================================================================
    # SCORING BACKFILL
    # ==================================================================

    def _on_compute_scores(self) -> None:
        """Open the scoring backfill dialog for selected entries."""
        if self.db is None:
            return
        selected = self._selected_entry_ids()
        if not selected:
            messagebox.showinfo("Scoring", "Select entries first.", parent=self)
            return

        entries = []
        for eid in selected:
            entry = self.db.get_entry(eid)
            if entry:
                entries.append(entry)

        if not entries:
            return

        ScoringBackfillDialog(
            self, self.db, entries, self.session_user,
            ion_list_path=self._ion_list_path,
            scoreb_wb_path=self._scoreb_workbook_path,
            on_done=self._on_scoring_done,
        )

    def _on_scoring_done(self, ion_path: str, wb_path: str) -> None:
        """Callback after scoring dialog completes. Updates session settings."""
        if ion_path and ion_path != self._ion_list_path:
            self._ion_list_path = ion_path
        if wb_path and wb_path != self._scoreb_workbook_path:
            self._scoreb_workbook_path = wb_path
        self._refresh_settings_display()
        self._save_config()
        self._populate_table()

    # ==================================================================
    # SPECTRUM PREVIEW
    # ==================================================================

    def _on_spectrum_preview(self, entry_id: str) -> None:
        if self.db is None:
            return
        if not _HAS_MATPLOTLIB:
            messagebox.showwarning(
                "Spectrum preview unavailable",
                "matplotlib is not installed. Install it with "
                "`pip install matplotlib` and restart the app to enable preview.",
                parent=self,
            )
            return
        entry = self.db.get_entry(entry_id)
        if not entry:
            return
        SpectrumPreviewDialog(self, entry, ion_list_path=self._ion_list_path)


# ---------------------------------------------------------------------------
# Import Preview Dialog
# ---------------------------------------------------------------------------

class ImportPreviewDialog(tk.Toplevel):
    """Preview and select rows from an imported file before inserting into DB.

    Args:
        master: Parent window.
        db: Database manager.
        rows: List of dicts from CSV/TSV reader.
        meta: Metadata dict (sample_id, ion_mode, derivatization).
        import_mode: One of 'converted_csv', 'cga', 'mas'.
        session_user: Current session user.
        on_done: Callback after import completes.
    """

    # Column mappings for the preview Treeview per mode
    _PREVIEW_COLS = {
        "converted_csv": [
            ("MS2scan_no", "Scan", 70),
            ("protonatedmass", "Precursor m/z", 120),
            ("peaklist", "Peak count", 80),
        ],
        "cga": [
            ("MS2scan_no", "Scan", 70),
            ("composition", "Composition", 140),
            ("ppm_error", "ppm", 70),
            ("ion score", "Score A", 70),
            ("score_b", "Score B", 70),
        ],
        "mas": [
            ("MS2scan_no", "Scan", 70),
            ("composition", "Composition", 140),
            ("protonatedmass", "Precursor m/z", 120),
        ],
    }

    def __init__(
        self,
        master: CuratedLibraryWindow,
        db: CuratedLibraryDB,
        rows: list[dict],
        meta: dict,
        import_mode: str,
        session_user: str,
        on_done: callable | None = None,
    ) -> None:
        super().__init__(master)
        self.db = db
        self.rows = rows
        self.meta = meta
        self.import_mode = import_mode
        self.session_user = session_user
        self.on_done = on_done

        mode_labels = {
            "converted_csv": "Converted CSV (no annotation)",
            "cga": "CGA Pseudolabel TSV",
            "mas": "MAS Annotation",
        }
        self.title(f"Import Preview \u2014 {mode_labels.get(import_mode, import_mode)}")
        self.geometry("780x560")
        self.resizable(True, True)
        self.transient(master)
        self.grab_set()

        self._has_scoreb = import_mode == "cga" and any(
            r.get("score_b") for r in rows
        )

        # Default: all rows shown, none selected. The selection set is
        # authoritative across filter changes: rows hidden by a later
        # filter change do NOT drop out of the import.
        self._shown_indices: list[int] = list(range(len(rows)))
        self._selected_indices: set[int] = set()

        self._build_ui()

    def _build_ui(self) -> None:
        info = ttk.Label(
            self,
            text=f"Sample: {self.meta.get('sample_id', '?')}  |  "
                 f"Ion mode: {self.meta.get('ion_mode', '?')}",
        )
        info.pack(fill="x", padx=10, pady=(8, 2))

        self._build_filter_bar()

        cols_spec = self._PREVIEW_COLS.get(self.import_mode, self._PREVIEW_COLS["cga"])
        self._cols_spec = cols_spec
        col_ids = [c[0] for c in cols_spec]

        tree_frame = ttk.Frame(self)
        tree_frame.pack(fill="both", expand=True, padx=10, pady=4)

        self._ptree = ttk.Treeview(
            tree_frame, columns=col_ids, show="headings", selectmode="extended",
        )
        for cid, heading, width in cols_spec:
            self._ptree.heading(cid, text=heading)
            self._ptree.column(cid, width=width, anchor="center", minwidth=50)

        vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=self._ptree.yview)
        self._ptree.configure(yscrollcommand=vsb.set)
        self._ptree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        tree_frame.columnconfigure(0, weight=1)
        tree_frame.rowconfigure(0, weight=1)

        self._ptree.bind("<<TreeviewSelect>>", self._on_tree_selection_changed)
        self._render_rows()

        self._summary_var = tk.StringVar()
        ttk.Label(self, textvariable=self._summary_var).pack(
            fill="x", padx=10, pady=(2, 4),
        )

        btn_frame = ttk.Frame(self)
        btn_frame.pack(fill="x", padx=10, pady=(4, 10))

        ttk.Button(btn_frame, text="Select All", command=self._select_all).pack(
            side="left", padx=4,
        )
        ttk.Button(btn_frame, text="Deselect All", command=self._deselect_all).pack(
            side="left", padx=4,
        )
        ttk.Button(
            btn_frame,
            text="Select all with composition",
            command=self._select_with_composition,
        ).pack(side="left", padx=4)

        ttk.Label(btn_frame, text="Score B >=").pack(side="left", padx=(8, 2))
        self._quick_scoreb_var = tk.StringVar()
        ttk.Entry(btn_frame, textvariable=self._quick_scoreb_var, width=6).pack(
            side="left",
        )
        ttk.Button(
            btn_frame, text="Select", command=self._select_with_scoreb,
        ).pack(side="left", padx=(2, 8))

        ttk.Button(btn_frame, text="Cancel", command=self.destroy).pack(
            side="right", padx=4,
        )
        ttk.Button(btn_frame, text="Import Selected", command=self._do_import).pack(
            side="right", padx=4,
        )

        self._update_summary()

    # ------------------------------------------------------------------
    # Filter bar
    # ------------------------------------------------------------------

    def _build_filter_bar(self) -> None:
        """Build the shared import-preview filter bar (v0.4 Item 1)."""
        fbar = ttk.LabelFrame(self, text="Filters")
        fbar.pack(fill="x", padx=10, pady=(4, 2))

        row1 = ttk.Frame(fbar)
        row1.pack(fill="x", padx=6, pady=(4, 2))

        self._fp_comp_present = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            row1, text="Composition present", variable=self._fp_comp_present,
        ).pack(side="left", padx=(0, 12))

        ttk.Label(row1, text="Score A min:").pack(side="left", padx=(0, 2))
        self._fp_score_a = ttk.Entry(row1, width=6)
        self._fp_score_a.pack(side="left", padx=(0, 12))

        ttk.Label(row1, text="Score B min:").pack(side="left", padx=(0, 2))
        self._fp_score_b = ttk.Entry(row1, width=6)
        self._fp_score_b.pack(side="left", padx=(0, 12))

        ttk.Label(row1, text="ppm max:").pack(side="left", padx=(0, 2))
        self._fp_ppm = ttk.Entry(row1, width=6)
        self._fp_ppm.pack(side="left", padx=(0, 12))

        row2 = ttk.Frame(fbar)
        row2.pack(fill="x", padx=6, pady=(2, 4))

        ttk.Label(row2, text="Sample:").pack(side="left", padx=(0, 2))
        self._fp_sample = ttk.Combobox(row2, state="readonly", width=24)
        unique_samples = sorted({
            str(
                r.get("sample_id") or r.get("Sample ID") or r.get("sample") or ""
            ).strip()
            for r in self.rows
        })
        unique_samples = [s for s in unique_samples if s]
        self._fp_sample["values"] = ["(all)", *unique_samples]
        self._fp_sample.set("(all)")
        self._fp_sample.pack(side="left", padx=(0, 12))

        ttk.Button(row2, text="Apply", command=self._on_apply_filters).pack(
            side="left", padx=4,
        )
        ttk.Button(row2, text="Clear", command=self._on_clear_filters).pack(
            side="left", padx=4,
        )

    def _on_apply_filters(self) -> None:
        """Recompute ``_shown_indices`` from the filter widgets.

        Non-numeric entries in numeric boxes are treated as "not set"
        (silent no-op) per Item 1's simplest-path rule.
        """
        def _parse_float(widget) -> float | None:
            text = widget.get().strip()
            if not text:
                return None
            try:
                return float(text)
            except ValueError:
                return None

        sample_val = self._fp_sample.get().strip()
        sample_arg: str | None = (
            sample_val if sample_val and sample_val != "(all)" else None
        )

        self._shown_indices = filter_import_rows(
            self.rows,
            composition_present=bool(self._fp_comp_present.get()),
            score_a_min=_parse_float(self._fp_score_a),
            score_b_min=_parse_float(self._fp_score_b),
            ppm_max=_parse_float(self._fp_ppm),
            sample=sample_arg,
        )
        self._render_rows()
        self._update_summary()

    def _on_clear_filters(self) -> None:
        self._fp_comp_present.set(False)
        self._fp_score_a.delete(0, "end")
        self._fp_score_b.delete(0, "end")
        self._fp_ppm.delete(0, "end")
        self._fp_sample.set("(all)")
        self._shown_indices = list(range(len(self.rows)))
        self._render_rows()
        self._update_summary()

    # ------------------------------------------------------------------
    # Rendering + selection bookkeeping
    # ------------------------------------------------------------------

    def _render_rows(self) -> None:
        """Populate the Treeview with the currently shown rows only.

        Rows hidden by the filter stay in ``_selected_indices`` but do not
        appear in the Treeview. The authoritative selection is re-applied
        to any rows still visible.
        """
        self._ptree.delete(*self._ptree.get_children())
        for i in self._shown_indices:
            row = self.rows[i]
            values = []
            for cid, _, _ in self._cols_spec:
                val = row.get(cid, "")
                if cid == "peaklist" and val:
                    try:
                        s = str(val).strip()
                        if s.startswith("(") or s.startswith("["):
                            import ast
                            val = str(len(ast.literal_eval(s)))
                        else:
                            val = str(len(s.split(";")))
                    except Exception:
                        val = "?"
                values.append(val or "\u2014")
            self._ptree.insert("", "end", iid=str(i), values=values)

        visible_selected = [
            str(i) for i in self._selected_indices if i in self._shown_indices
        ]
        if visible_selected:
            self._ptree.selection_set(visible_selected)

    def _on_tree_selection_changed(self, _event: tk.Event) -> None:
        """Sync ``_selected_indices`` with the Treeview.

        Only rows currently visible contribute. Rows hidden by the filter
        retain whatever selection state they had before the filter change.
        """
        visible_set = set(self._shown_indices)
        tree_sel = {int(iid) for iid in self._ptree.selection()}
        hidden_selected = {
            i for i in self._selected_indices if i not in visible_set
        }
        self._selected_indices = hidden_selected | tree_sel
        self._update_summary()

    def _set_selection(self, new_selection: set[int]) -> None:
        self._selected_indices = set(new_selection)
        visible_selected = [
            str(i) for i in self._selected_indices if i in self._shown_indices
        ]
        self._ptree.selection_set(visible_selected)
        self._update_summary()

    def _update_summary(self) -> None:
        self._summary_var.set(
            f"Rows: {len(self.rows)}  |  Shown: {len(self._shown_indices)}  "
            f"|  Selected: {len(self._selected_indices)}"
        )

    # ------------------------------------------------------------------
    # Quick-select actions (all operate on Shown rows only)
    # ------------------------------------------------------------------

    def _select_all(self) -> None:
        self._set_selection(self._selected_indices | set(self._shown_indices))

    def _deselect_all(self) -> None:
        shown = set(self._shown_indices)
        self._set_selection(self._selected_indices - shown)

    def _select_with_composition(self) -> None:
        to_add = set()
        for i in self._shown_indices:
            comp = str(
                self.rows[i].get("composition")
                or self.rows[i].get("selected_composition")
                or ""
            ).strip()
            if comp:
                to_add.add(i)
        self._set_selection(self._selected_indices | to_add)

    def _select_with_scoreb(self) -> None:
        """Select shown rows with ``score_b >= threshold``.

        Blank / non-numeric threshold is a silent no-op. Missing
        ``score_b`` cell means the row is simply not added.
        """
        text = self._quick_scoreb_var.get().strip()
        if not text:
            return
        try:
            threshold = float(text)
        except ValueError:
            return
        to_add = set()
        for i in self._shown_indices:
            raw = self.rows[i].get("score_b")
            if raw is None or str(raw).strip() == "":
                continue
            try:
                if float(raw) >= threshold:
                    to_add.add(i)
            except (ValueError, TypeError):
                continue
        self._set_selection(self._selected_indices | to_add)

    def _do_import(self) -> None:
        # The authoritative selection set survives filter tweaks:
        # import every index we have, even if currently hidden.
        sel_indices = sorted(self._selected_indices)
        if not sel_indices:
            messagebox.showwarning("Import", "No rows selected.", parent=self)
            return

        entries = []
        for idx in sel_indices:
            row = self.rows[idx]
            entry = self._row_to_entry(row)
            if entry is not None:
                entries.append(entry)

        if not entries:
            messagebox.showwarning("Import", "No valid rows to import.", parent=self)
            return

        # Insert with duplicate handling
        inserted = 0
        duplicates = 0
        for entry in entries:
            try:
                self.db.insert_entry(entry)
                inserted += 1
            except DuplicateEntryError:
                duplicates += 1

        msg = f"Imported {inserted} entries."
        if duplicates:
            msg += f"\n{duplicates} duplicates skipped."
        messagebox.showinfo("Import", msg, parent=self)

        self.destroy()
        if self.on_done:
            self.on_done()

    def _row_to_entry(self, row: dict) -> dict | None:
        """Convert an imported CSV/TSV row to a library entry dict."""
        meta = self.meta

        # Scan number — try several column names
        scan = (
            row.get("MS2scan_no")
            or row.get("ms2_scan_no")
            or row.get("Scan")
            or row.get("scan")
        )
        if scan is None:
            return None
        try:
            scan = int(float(str(scan)))
        except (ValueError, TypeError):
            return None

        # Precursor m/z
        prec_mz = _safe_float(
            row.get("protonatedmass") or row.get("precursor_mz") or row.get("PrecursorMZ"),
            0.0,
        )

        # Peak data
        peaklist = row.get("peaklist", "")
        peakintensity = row.get("peakintensity", "")
        if not peaklist:
            return None

        # Base entry
        entry: dict = {
            "sample_id": meta["sample_id"],
            "ms2_scan_no": scan,
            "precursor_mz": prec_mz,
            "precursor_charge": -1 if meta["ion_mode"] == "negative" else 1,
            "precursor_adduct": "[M-H]-" if meta["ion_mode"] == "negative" else "[M+H]+",
            "ion_mode": meta["ion_mode"],
            "derivatization": meta["derivatization"],
            "peaklist": peaklist,
            "peakintensity": peakintensity,
            "created_by": self.session_user,
            "last_modified_by": self.session_user,
            "reviewed_by": self.session_user,
        }

        if self.import_mode == "converted_csv":
            # No composition — placeholder
            entry["composition"] = ""
            entry["theoretical_mass"] = 0.0
            entry["ppm_error"] = 0.0
            entry["annotation_source"] = "manual-validated"
            entry["confidence"] = "tentative"

        elif self.import_mode == "cga":
            comp = row.get("composition") or row.get("selected_composition") or ""
            entry["composition"] = comp
            entry["theoretical_mass"] = _safe_float(row.get("theoretical_mass"), 0.0)
            entry["ppm_error"] = _safe_float(row.get("ppm_error"), 0.0)
            entry["score_a"] = _safe_float(row.get("ion score") or row.get("score_a"))
            entry["ion_hit_count"] = (
                int(float(row["ion hit count"]))
                if row.get("ion hit count")
                else None
            )
            entry["ion_hits_mz"] = row.get("ion hits m/z") or row.get("ion_hits_mz")

            # Score B columns
            entry["score_b"] = _safe_float(row.get("score_b"))
            entry["score_b_support"] = _safe_float(row.get("score_b_support_score") or row.get("score_b_support"))
            entry["score_b_comp_penalty"] = _safe_float(
                row.get("score_b_composition_penalty") or row.get("score_b_comp_penalty")
            )
            entry["score_b_unexp_penalty"] = _safe_float(
                row.get("score_b_unexpected_penalty") or row.get("score_b_unexp_penalty")
            )
            entry["score_b_motif_summary"] = row.get("score_b_motif_summary")

            has_sb = entry["score_b"] is not None
            entry["annotation_source"] = "CGA-RouteA+ScoreB" if has_sb else "CGA-RouteA"
            entry["confidence"] = "tentative"

        elif self.import_mode == "mas":
            comp = row.get("composition") or row.get("selected_composition") or ""
            entry["composition"] = comp
            entry["theoretical_mass"] = _safe_float(row.get("theoretical_mass"), 0.0)
            entry["ppm_error"] = _safe_float(row.get("ppm_error"), 0.0)
            entry["annotation_source"] = "manual-RouteB"
            entry["confidence"] = "probable"

        return entry


# ---------------------------------------------------------------------------
# New Entry Dialog (b5 — single manual entry)
# ---------------------------------------------------------------------------

class NewEntryDialog(tk.Toplevel):
    """Dialog for manually creating a single library entry.

    Args:
        master: Parent window.
        db: Database manager.
        session_user: Current session user.
        on_done: Callback after insertion.
    """

    def __init__(
        self,
        master: CuratedLibraryWindow,
        db: CuratedLibraryDB,
        session_user: str,
        on_done: callable | None = None,
    ) -> None:
        super().__init__(master)
        self.db = db
        self.session_user = session_user
        self.on_done = on_done

        self.title("New Library Entry (Manual)")
        self.geometry("520x600")
        self.resizable(True, True)
        self.transient(master)
        self.grab_set()

        self._vars: dict[str, tk.StringVar] = {}
        self._build_form()
        self._build_paste_bar()

    def _build_form(self) -> None:
        canvas = tk.Canvas(self, borderwidth=0)
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        frame = ttk.Frame(canvas, padding=12)

        frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        fields = [
            ("sample_id", "Sample ID:", ""),
            ("ms2_scan_no", "MS2 Scan No:", ""),
            ("precursor_mz", "Precursor m/z:", ""),
            ("precursor_charge", "Precursor Charge:", "-1"),
            ("precursor_adduct", "Precursor Adduct:", "[M-H]-"),
            ("ion_mode", "Ion Mode:", "negative"),
            ("derivatization", "Derivatization:", "permethylated"),
            ("composition", "Composition:", ""),
            ("theoretical_mass", "Theoretical Mass:", ""),
            ("ppm_error", "ppm Error:", ""),
            ("confidence", "Confidence:", "tentative"),
            ("annotation_source", "Annotation Source:", "manual-validated"),
            ("peaklist", "Peaklist (;-sep):", ""),
            ("peakintensity", "Peak Intensity (;-sep):", ""),
            ("diagnostic_notes", "Diagnostic Notes:", ""),
        ]

        for row_i, (key, label, default) in enumerate(fields):
            ttk.Label(frame, text=label).grid(row=row_i, column=0, sticky="w", padx=4, pady=3)
            var = tk.StringVar(value=default)
            self._vars[key] = var
            ttk.Entry(frame, textvariable=var, width=40).grid(
                row=row_i, column=1, sticky="ew", padx=4, pady=3,
            )

        frame.columnconfigure(1, weight=1)

        # Buttons — created here but packed later (after paste bar) by _build_paste_bar
        self._btn_frame = ttk.Frame(self)
        ttk.Button(self._btn_frame, text="Insert", command=self._on_insert).pack(side="right", padx=4)
        ttk.Button(self._btn_frame, text="Cancel", command=self.destroy).pack(side="right", padx=4)

    def _build_paste_bar(self) -> None:
        """Add a 'Paste from Clipboard' bar, then pack the button frame below it."""
        paste_frame = ttk.LabelFrame(self, text="Quick Fill from Spectrum Clipboard")
        paste_frame.pack(fill="x", padx=10, pady=(0, 4))

        ttk.Button(
            paste_frame, text="Paste from Clipboard",
            command=self._on_paste_clipboard,
        ).pack(side="left", padx=8, pady=6)

        self._paste_status = ttk.Label(paste_frame, text="Thermo Xcalibur supported", foreground="gray")
        self._paste_status.pack(side="left", padx=4, pady=6)

        # Now pack the button frame (Insert / Cancel) below the paste bar
        self._btn_frame.pack(fill="x", padx=10, pady=(4, 10))

    # Maximum peaks to display in Entry widgets before truncating.
    # Full data is stored in _paste_full_peaks for insertion.
    _PEAK_DISPLAY_LIMIT = 500

    def _on_paste_clipboard(self) -> None:
        """Read clipboard, decode vendor spectrum, and pre-fill form fields."""
        try:
            clip_text = self.clipboard_get()
        except tk.TclError:
            self._paste_status.configure(text="Clipboard is empty", foreground="red")
            return

        from msp_clipboard_decoder import parse_clipboard
        result = parse_clipboard(clip_text)

        if result["format"] == "unsupported":
            self._paste_status.configure(text="Unsupported format", foreground="red")
            messagebox.showwarning(
                "Clipboard Format",
                result.get("message", "Unrecognised clipboard format."),
                parent=self,
            )
            return

        entry = result["entry"]
        peak_count = entry.get("peak_count", 0)

        # Warn about very large peak lists (likely profile-mode, not centroided)
        if peak_count > 2000:
            proceed = messagebox.askyesno(
                "Large Peak List",
                f"The pasted spectrum contains {peak_count:,} peaks.\n\n"
                "This is likely profile-mode (not centroided) data. "
                "Profile spectra are very large and may not be suitable "
                "for direct annotation or scoring.\n\n"
                "Consider using centroided data instead.\n\n"
                "Proceed anyway?",
                parent=self,
            )
            if not proceed:
                self._paste_status.configure(
                    text=f"Cancelled ({peak_count:,} peaks — profile?)", foreground="orange",
                )
                return

        # Store full peak data internally (avoid rendering huge strings in Entry)
        self._paste_full_peaks = {
            "peaklist": entry.get("peaklist", ""),
            "peakintensity": entry.get("peakintensity", ""),
        }

        # Map decoded fields to form variables
        field_map = {
            "sample_id": "sample_id",
            "ms2_scan_no": "ms2_scan_no",
            "precursor_mz": "precursor_mz",
            "precursor_charge": "precursor_charge",
            "precursor_adduct": "precursor_adduct",
            "ion_mode": "ion_mode",
            "diagnostic_notes": "diagnostic_notes",
            "confidence": "confidence",
            "annotation_source": "annotation_source",
        }

        filled = 0
        for src_key, var_key in field_map.items():
            val = str(entry.get(src_key, "")).strip()
            if val and var_key in self._vars:
                self._vars[var_key].set(val)
                filled += 1

        # For peak fields: show truncated preview in Entry, keep full data for insert
        for peak_key in ("peaklist", "peakintensity"):
            full_val = entry.get(peak_key, "")
            if full_val and peak_key in self._vars:
                parts = full_val.split(";")
                if len(parts) > self._PEAK_DISPLAY_LIMIT:
                    display = ";".join(parts[:self._PEAK_DISPLAY_LIMIT]) + f"  [...{len(parts)} total]"
                else:
                    display = full_val
                self._vars[peak_key].set(display)
                filled += 1

        # Show success with peak count
        fmt_name = result["format"].capitalize()
        extra = ""
        if peak_count > 2000:
            extra = " (profile?)"
        self._paste_status.configure(
            text=f"{fmt_name}: {peak_count:,} peaks{extra} ({filled} fields)",
            foreground="green",
        )

    def _on_insert(self) -> None:
        data: dict = {}
        for key, var in self._vars.items():
            val = var.get().strip()
            if val:
                data[key] = val

        # If clipboard paste stored full peak data, use it instead of truncated display
        if hasattr(self, "_paste_full_peaks"):
            for peak_key in ("peaklist", "peakintensity"):
                full = self._paste_full_peaks.get(peak_key, "")
                if full:
                    data[peak_key] = full

        # Required field checks
        for req in ["sample_id", "ms2_scan_no", "composition", "peaklist", "peakintensity"]:
            if req not in data:
                messagebox.showwarning("Missing Field", f"{req} is required.", parent=self)
                return

        # Type conversions — scan number can be int or composite ID (e.g. "AV12_5467-32694")
        try:
            data["ms2_scan_no"] = int(data["ms2_scan_no"])
        except ValueError:
            pass  # keep as string for composite/averaged spectra
        data["precursor_mz"] = _safe_float(data.get("precursor_mz"), 0.0)
        data["precursor_charge"] = int(_safe_float(data.get("precursor_charge"), -1))
        data["theoretical_mass"] = _safe_float(data.get("theoretical_mass"), 0.0)
        data["ppm_error"] = _safe_float(data.get("ppm_error"), 0.0)

        data.setdefault("precursor_adduct", "[M-H]-")
        data.setdefault("ion_mode", "negative")
        data.setdefault("derivatization", "permethylated")
        data.setdefault("annotation_source", "manual-validated")
        data.setdefault("confidence", "tentative")
        data["created_by"] = self.session_user
        data["last_modified_by"] = self.session_user
        data["reviewed_by"] = self.session_user

        try:
            eid = self.db.insert_entry(data)
            messagebox.showinfo("Inserted", f"Created {eid}.", parent=self)
            self.destroy()
            if self.on_done:
                self.on_done()
        except DuplicateEntryError:
            messagebox.showwarning(
                "Duplicate",
                "An entry with this sample + scan + composition already exists.",
                parent=self,
            )


# ---------------------------------------------------------------------------
# Edit Entry Dialog
# ---------------------------------------------------------------------------

class EditEntryDialog(tk.Toplevel):
    """Dialog for editing a single library entry.

    Args:
        master: Parent window.
        db: Database manager instance.
        entry_id: ID of the entry to edit.
        session_user: Current session user name.
        on_done: Callback invoked after save or delete.
        ion_list_path: Current session ion list (for display).
        scoreb_wb_path: Current session Score B workbook (for display).
    """

    def __init__(
        self,
        master: tk.Toplevel,
        db: CuratedLibraryDB,
        entry_id: str,
        session_user: str,
        on_done: callable | None = None,
        ion_list_path: str = "",
        scoreb_wb_path: str = "",
    ) -> None:
        super().__init__(master)
        self.db = db
        self.entry_id = entry_id
        self.session_user = session_user
        self.on_done = on_done

        self.entry = db.get_entry(entry_id)
        if self.entry is None:
            messagebox.showerror(
                "Error", f"Entry {entry_id} not found.", parent=master
            )
            self.destroy()
            return

        self.title(f"Edit Entry: {entry_id}")
        self.geometry("500x520")
        self.resizable(True, True)
        self.transient(master)
        self.grab_set()

        self._ion_list_path = ion_list_path
        self._scoreb_wb_path = scoreb_wb_path
        self._build_form()

    def _build_form(self) -> None:
        e = self.entry
        pad = {"padx": 10, "pady": 4, "sticky": "w"}

        frame = ttk.Frame(self)
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        row = 0

        # Composition (editable)
        ttk.Label(frame, text="Composition:").grid(row=row, column=0, **pad)
        self._comp_var = tk.StringVar(value=e.get("composition", ""))
        ttk.Entry(frame, textvariable=self._comp_var, width=30).grid(
            row=row, column=1, sticky="ew", padx=10, pady=4,
        )
        row += 1

        # Confidence (dropdown)
        ttk.Label(frame, text="Confidence:").grid(row=row, column=0, **pad)
        self._conf_var = tk.StringVar(value=e.get("confidence", "tentative"))
        ttk.Combobox(
            frame, textvariable=self._conf_var, width=14, state="readonly",
            values=["confirmed", "probable", "tentative"],
        ).grid(row=row, column=1, sticky="w", padx=10, pady=4)
        row += 1

        # GlyTouCan ID (editable)
        ttk.Label(frame, text="GlyTouCan ID:").grid(row=row, column=0, **pad)
        self._glyto_var = tk.StringVar(value=e.get("glytoucan_id", "") or "")
        ttk.Entry(frame, textvariable=self._glyto_var, width=20).grid(
            row=row, column=1, sticky="w", padx=10, pady=4,
        )
        row += 1

        # Read-only scores
        ttk.Label(frame, text="Score A:").grid(row=row, column=0, **pad)
        sa = f"{e['score_a']:.4f}" if e.get("score_a") is not None else "\u2014"
        ttk.Label(frame, text=sa).grid(row=row, column=1, sticky="w", padx=10, pady=4)
        row += 1

        ttk.Label(frame, text="Score B:").grid(row=row, column=0, **pad)
        sb = f"{e['score_b']:.4f}" if e.get("score_b") is not None else "\u2014"
        ttk.Label(frame, text=sb).grid(row=row, column=1, sticky="w", padx=10, pady=4)
        row += 1

        # Annotation source (read-only)
        ttk.Label(frame, text="Annotation source:").grid(row=row, column=0, **pad)
        ttk.Label(frame, text=e.get("annotation_source", "\u2014")).grid(
            row=row, column=1, sticky="w", padx=10, pady=4,
        )
        row += 1

        # Scoring reference (read-only)
        ttk.Label(frame, text="Scored with:").grid(row=row, column=0, **pad)
        ref_parts = []
        if self._ion_list_path:
            ref_parts.append(Path(self._ion_list_path).name)
        if self._scoreb_wb_path:
            ref_parts.append(Path(self._scoreb_wb_path).name)
        ref_text = " + ".join(ref_parts) if ref_parts else "(not configured)"
        ttk.Label(frame, text=ref_text).grid(
            row=row, column=1, sticky="w", padx=10, pady=4,
        )
        row += 1

        # PDF reference
        ttk.Label(frame, text="PDF Reference:").grid(row=row, column=0, **pad)
        pdf_frame = ttk.Frame(frame)
        pdf_frame.grid(row=row, column=1, sticky="ew", padx=10, pady=4)
        self._pdf_var = tk.StringVar(value=e.get("pdf_reference", "") or "")
        ttk.Label(pdf_frame, textvariable=self._pdf_var, width=28).pack(side="left")
        ttk.Button(pdf_frame, text="Browse\u2026", command=self._on_browse_pdf).pack(
            side="left", padx=4
        )
        row += 1

        # Diagnostic notes (Text widget)
        ttk.Label(frame, text="Diagnostic Notes:").grid(
            row=row, column=0, padx=10, pady=4, sticky="nw",
        )
        self._notes_text = tk.Text(frame, width=36, height=5, wrap="word")
        self._notes_text.grid(row=row, column=1, sticky="nsew", padx=10, pady=4)
        self._notes_text.insert("1.0", e.get("diagnostic_notes", "") or "")
        row += 1

        frame.columnconfigure(1, weight=1)
        frame.rowconfigure(row - 1, weight=1)

        # Buttons
        btn_frame = ttk.Frame(self)
        btn_frame.pack(fill="x", padx=10, pady=(0, 10))
        ttk.Button(btn_frame, text="Save", command=self._on_save).pack(
            side="right", padx=4
        )
        ttk.Button(btn_frame, text="Cancel", command=self.destroy).pack(
            side="right", padx=4
        )
        ttk.Button(btn_frame, text="Delete", command=self._on_delete).pack(
            side="left", padx=4
        )

    def _on_browse_pdf(self) -> None:
        path = filedialog.askopenfilename(
            title="Select PDF",
            filetypes=[("PDF files", "*.pdf"), ("All files", "*.*")],
            parent=self,
        )
        if path:
            rel = self.db.attach_pdf(self.entry_id, path)
            self._pdf_var.set(rel)

    def _on_save(self) -> None:
        updates: dict = {}
        new_comp = self._comp_var.get().strip()
        if new_comp != self.entry.get("composition", ""):
            updates["composition"] = new_comp
        new_conf = self._conf_var.get().strip()
        if new_conf != self.entry.get("confidence", ""):
            updates["confidence"] = new_conf
        new_glyto = self._glyto_var.get().strip() or None
        if new_glyto != self.entry.get("glytoucan_id"):
            updates["glytoucan_id"] = new_glyto
        new_notes = self._notes_text.get("1.0", "end-1c").strip() or None
        if new_notes != self.entry.get("diagnostic_notes"):
            updates["diagnostic_notes"] = new_notes

        if updates:
            updates["last_modified_by"] = self.session_user
            updates["reviewed_by"] = self.session_user
            self.db.update_entry(self.entry_id, updates)
        self.destroy()
        if self.on_done:
            self.on_done()

    def _on_delete(self) -> None:
        if messagebox.askyesno(
            "Confirm Delete",
            f"Delete {self.entry_id}? This cannot be undone.",
            parent=self,
        ):
            self.db.delete_entry(self.entry_id)
            self.destroy()
            if self.on_done:
                self.on_done()


# ---------------------------------------------------------------------------
# Stats Dialog
# ---------------------------------------------------------------------------

class StatsDialog(tk.Toplevel):
    """Simple dialog showing library statistics."""

    def __init__(self, master: tk.Toplevel, db: CuratedLibraryDB) -> None:
        super().__init__(master)
        self.title("Library Statistics")
        self.geometry("320x240")
        self.resizable(False, False)
        self.transient(master)
        self.grab_set()

        stats = db.get_stats()

        frame = ttk.Frame(self, padding=16)
        frame.pack(fill="both", expand=True)

        row = 0
        for label, val in [
            ("Total entries:", stats["total_entries"]),
            ("Unique samples:", stats["unique_samples"]),
            ("Unique compositions:", stats["unique_compositions"]),
        ]:
            ttk.Label(frame, text=label, font=("TkDefaultFont", 11)).grid(
                row=row, column=0, sticky="w", pady=2,
            )
            ttk.Label(
                frame, text=str(val), font=("TkDefaultFont", 11, "bold")
            ).grid(row=row, column=1, sticky="e", padx=(12, 0), pady=2)
            row += 1

        # Confidence breakdown
        ttk.Separator(frame, orient="horizontal").grid(
            row=row, column=0, columnspan=2, sticky="ew", pady=8,
        )
        row += 1
        ttk.Label(
            frame, text="By confidence:", font=("TkDefaultFont", 10, "bold")
        ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(0, 4))
        row += 1
        for conf_level in ["confirmed", "probable", "tentative"]:
            count = stats["confidence_counts"].get(conf_level, 0)
            ttk.Label(frame, text=f"  {conf_level}:").grid(
                row=row, column=0, sticky="w"
            )
            ttk.Label(frame, text=str(count)).grid(
                row=row, column=1, sticky="e", padx=(12, 0)
            )
            row += 1

        # BUG FIX: close button — pack in self (outside frame) with explicit text
        close_btn = ttk.Button(self, text="Close", command=self.destroy)
        close_btn.pack(side="bottom", pady=(4, 10))


# ---------------------------------------------------------------------------
# Batch Confidence Dialog (v0.4 Item 2)
# ---------------------------------------------------------------------------


class BatchConfidenceDialog(tk.Toplevel):
    """Confirmation dialog for batch-updating ``confidence`` across rows.

    Presented from the right-click ``Set Confidence for Selected...`` menu
    item. Collects a confidence level (``confirmed`` / ``probable`` /
    ``tentative``) and an optional reviewer override; empty reviewer field
    falls back to the session user.

    Args:
        master: Parent window.
        n_selected: Count of rows that will be updated. Shown in the
            Apply-button label and the intro text.
        session_user: Default reviewer when the override entry is blank.
        on_apply: Callback ``(confidence, reviewer)`` invoked on Apply.
    """

    def __init__(
        self,
        master: tk.Toplevel,
        n_selected: int,
        session_user: str,
        on_apply,
    ) -> None:
        super().__init__(master)
        self.title("Set Confidence for Selected")
        self.geometry("360x220")
        self.resizable(False, False)
        self.transient(master)
        self.grab_set()

        self._on_apply = on_apply

        frame = ttk.Frame(self, padding=14)
        frame.pack(fill="both", expand=True)

        ttk.Label(
            frame,
            text=f"Apply to {n_selected} selected entr"
                 f"{'y' if n_selected == 1 else 'ies'}.",
        ).pack(anchor="w", pady=(0, 8))

        ttk.Label(frame, text="Confidence:").pack(anchor="w")
        self._conf_var = tk.StringVar(value="probable")
        for value in ("confirmed", "probable", "tentative"):
            ttk.Radiobutton(
                frame, text=value, value=value, variable=self._conf_var,
            ).pack(anchor="w", padx=16)

        rev_row = ttk.Frame(frame)
        rev_row.pack(fill="x", pady=(8, 0))
        ttk.Label(rev_row, text="Reviewer override:").pack(side="left")
        self._reviewer_var = tk.StringVar(value="")
        ttk.Entry(rev_row, textvariable=self._reviewer_var, width=16).pack(
            side="left", padx=6,
        )
        ttk.Label(
            rev_row,
            text=f"(blank = {session_user})",
            foreground="gray",
        ).pack(side="left")

        btns = ttk.Frame(frame)
        btns.pack(fill="x", pady=(12, 0))
        ttk.Button(btns, text="Cancel", command=self.destroy).pack(
            side="right", padx=4,
        )
        ttk.Button(
            btns, text=f"Apply to {n_selected}", command=self._do_apply,
        ).pack(side="right", padx=4)

    def _do_apply(self) -> None:
        confidence = self._conf_var.get()
        reviewer = self._reviewer_var.get().strip()
        try:
            self._on_apply(confidence, reviewer)
        finally:
            self.destroy()


# ---------------------------------------------------------------------------
# Spectrum Preview Dialog
# ---------------------------------------------------------------------------

class SpectrumPreviewDialog(tk.Toplevel):
    """Read-only spectrum preview as a stick plot.

    Uses matplotlib embedded in Tkinter via FigureCanvasTkAgg. Ion-hit peaks
    are rendered in red and, when the current ion list supplies names,
    labeled with the resolved ``fragment_name`` (or the m/z rounded to two
    decimals as a fallback). A checkbutton toggles the labels off for dense
    profile-mode spectra.

    Args:
        master: Parent window.
        entry: Dict of entry column values.
        ion_list_path: Path to the active ion list (session setting) — used
            to resolve fragment names. When empty or unreadable, the dialog
            silently falls back to m/z-only labels.
    """

    _LABEL_ROTATION_DEG = 30  # v0.4 spec: default rotation, stick-tip anchor.
    _HIT_TOLERANCE = 0.5      # m/z tolerance for hit highlighting + name lookup.

    def __init__(
        self,
        master: tk.Toplevel,
        entry: dict,
        ion_list_path: str = "",
    ) -> None:
        super().__init__(master)
        self.entry = entry
        self._ion_list_path = ion_list_path
        self._ion_df = None  # Lazily loaded; None on failure or empty path.
        self._show_labels_var = tk.BooleanVar(value=True)
        # Populated during _build_peak_data; reused by plot + peak table.
        self._hit_name_by_mz: dict[float, str] = {}
        self._canvas_widget = None  # matplotlib canvas widget (for re-render).

        eid = entry.get("entry_id", "?")
        scan = entry.get("ms2_scan_no", "?")
        self.title(f"Spectrum Preview: {eid} (Scan {scan})")
        self.geometry("750x720")
        self.resizable(True, True)
        self.transient(master)

        # Info bar at the top
        info_frame = ttk.Frame(self)
        info_frame.pack(fill="x", padx=10, pady=(8, 4))

        prec = entry.get("precursor_mz")
        adduct = entry.get("precursor_adduct", "")
        charge = entry.get("precursor_charge", "")
        comp = entry.get("composition", "")
        ppm = entry.get("ppm_error")
        sa = entry.get("score_a")
        sb = entry.get("score_b")

        info_parts = [f"Precursor: {prec} {adduct} (z={charge})"]
        if comp:
            info_parts.append(f"Composition: {comp}")
        if ppm is not None:
            info_parts.append(f"ppm: {ppm}")
        if sa is not None:
            info_parts.append(f"Score A: {sa}")
        if sb is not None:
            info_parts.append(f"Score B: {sb}")
        ttk.Label(info_frame, text="   ".join(info_parts)).pack(anchor="w")

        # Label toggle — placed between info bar and plot so the user can
        # toggle before/after reading the data.
        ctrl_frame = ttk.Frame(self)
        ctrl_frame.pack(fill="x", padx=10, pady=(0, 2))
        ttk.Checkbutton(
            ctrl_frame,
            text="Show fragment labels",
            variable=self._show_labels_var,
            command=self._on_toggle_labels,
        ).pack(side="left")

        # Load ion_df once (lazy failures are silently absorbed — label
        # rendering degrades to m/z-only in that case).
        self._load_ion_df_safe()
        self._prepare_hit_name_map()

        try:
            self._build_matplotlib_plot()
        except ImportError:
            self._build_text_fallback()

        # Peak table below the plot
        self._build_peak_table()

        ttk.Button(self, text="Close", command=self.destroy).pack(
            side="bottom", pady=(0, 8)
        )

    def _load_ion_df_safe(self) -> None:
        """Attempt to load the active ion list; on any failure, leave ``None``."""
        if not self._ion_list_path:
            return
        if not Path(self._ion_list_path).exists():
            return
        try:
            self._ion_df = _load_ion_df(self._ion_list_path)
        except Exception:
            # Any parse error (missing sheet, malformed CSV) just disables
            # label resolution — the preview itself must still open.
            self._ion_df = None

    def _prepare_hit_name_map(self) -> None:
        """Resolve ``hit_mz -> fragment_name`` for every ion-hit peak.

        Uses the stored ``ion_hits_mz`` (from scoring) as the authoritative
        hit list and looks up each against the active ion list. Peaks that
        are not hits never get a label.
        """
        _, _, ion_hits = self._parse_peaks()
        if not ion_hits:
            return
        names = lookup_fragment_names(
            self._ion_df, ion_hits, tolerance=self._HIT_TOLERANCE,
        )
        self._hit_name_by_mz = {
            mz: name for mz, name in zip(ion_hits, names)
        }

    def _parse_peaks(self) -> tuple[list[float], list[float], list[float]]:
        """Parse peaklist, peakintensity, and ion_hits_mz from entry."""
        def _parse_list(raw: str | None) -> list[float]:
            if not raw:
                return []
            s = str(raw).strip()
            # Handle both semicolon-separated and Python list format
            if s.startswith("[") or s.startswith("("):
                try:
                    import ast
                    return [float(x) for x in ast.literal_eval(s)]
                except Exception:
                    pass
            return [float(x) for x in s.replace(",", ";").split(";") if x.strip()]

        mzs = _parse_list(self.entry.get("peaklist"))
        intensities = _parse_list(self.entry.get("peakintensity"))
        ion_hits = _parse_list(self.entry.get("ion_hits_mz"))
        return mzs, intensities, ion_hits

    def _build_matplotlib_plot(self) -> None:
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure

        mzs, intensities, ion_hits = self._parse_peaks()

        fig = Figure(figsize=(7, 4), dpi=96)
        ax = fig.add_subplot(111)
        self._fig = fig
        self._ax = ax

        if mzs and intensities and len(mzs) == len(intensities):
            ion_hit_set = self._build_ion_hit_set(mzs, ion_hits)
            self._ion_hit_set = ion_hit_set

            for mz, intensity in zip(mzs, intensities):
                color = "#d62728" if mz in ion_hit_set else "#1f77b4"
                ax.vlines(mz, 0, intensity, colors=color, linewidth=1.0)

            if self._show_labels_var.get():
                self._draw_fragment_labels(ax, mzs, intensities, ion_hit_set)

            ax.set_xlabel("m/z")
            ax.set_ylabel("Intensity")

            from matplotlib.lines import Line2D
            handles = [Line2D([0], [0], color="#1f77b4", lw=2, label="All peaks")]
            if ion_hit_set:
                handles.append(
                    Line2D([0], [0], color="#d62728", lw=2, label="Ion hits")
                )
            ax.legend(handles=handles, loc="upper right", fontsize=8)
        else:
            ax.text(
                0.5, 0.5, "No peak data available",
                ha="center", va="center", transform=ax.transAxes,
            )

        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.draw()
        self._canvas = canvas
        self._canvas_widget = canvas.get_tk_widget()
        self._canvas_widget.pack(fill="both", expand=True, padx=8, pady=(8, 4))

    def _build_ion_hit_set(
        self, mzs: list[float], ion_hits: list[float],
    ) -> set[float]:
        """Match observed peak m/z to ion-hit m/z within display tolerance."""
        if not ion_hits:
            return set()
        out: set[float] = set()
        for mz in mzs:
            for hit_mz in ion_hits:
                if abs(mz - hit_mz) <= self._HIT_TOLERANCE:
                    out.add(mz)
                    break
        return out

    def _draw_fragment_labels(
        self,
        ax,
        mzs: list[float],
        intensities: list[float],
        ion_hit_set: set[float],
    ) -> None:
        """Draw label text above every ion-hit stick.

        Label rule (v0.4):
            - fragment_name known for this hit → fragment_name
            - no fragment_name, only m/z → ``"{m/z:.2f}"``
            - non-hit peaks → skipped entirely

        Labels rotated 30°, anchored at the stick tip. Minor overlap is
        acceptable for v0.4; no manual placement UI.
        """
        for mz, intensity in zip(mzs, intensities):
            if mz not in ion_hit_set:
                continue
            label = self._label_for_hit(mz)
            if not label:
                continue
            ax.text(
                mz,
                intensity,
                label,
                rotation=self._LABEL_ROTATION_DEG,
                rotation_mode="anchor",
                ha="left",
                va="bottom",
                fontsize=7,
                color="#d62728",
            )

    def _label_for_hit(self, mz: float) -> str:
        """Pick the label text for a single ion-hit m/z."""
        # Match an entry in self._hit_name_by_mz by the nearest ion-hit mz.
        best_key: float | None = None
        best_diff = self._HIT_TOLERANCE
        for hit_mz in self._hit_name_by_mz:
            diff = abs(mz - hit_mz)
            if diff <= best_diff:
                best_diff = diff
                best_key = hit_mz
        name = self._hit_name_by_mz.get(best_key, "") if best_key is not None else ""
        if name:
            return name
        return f"{mz:.2f}"

    def _on_toggle_labels(self) -> None:
        """Redraw the plot when the fragment-label checkbox toggles."""
        if not hasattr(self, "_ax") or not hasattr(self, "_canvas"):
            return
        mzs, intensities, ion_hits = self._parse_peaks()
        if not (mzs and intensities and len(mzs) == len(intensities)):
            return
        ion_hit_set = self._ion_hit_set
        # Remove only the text labels; keep vlines and legend as-is.
        for text in list(self._ax.texts):
            text.remove()
        if self._show_labels_var.get():
            self._draw_fragment_labels(self._ax, mzs, intensities, ion_hit_set)
        self._canvas.draw_idle()

    def _build_peak_table(self) -> None:
        """Build a scrollable peak list panel with ion-hit highlighting and
        a ``Fragment`` column populated from the resolved fragment names.

        Rows where no fragment name is available (BY-style hits or non-hits)
        display a blank Fragment cell.
        """
        mzs, intensities, ion_hits = self._parse_peaks()
        if not mzs:
            return

        ion_hit_indices: set[int] = set()
        if ion_hits:
            for i, mz in enumerate(mzs):
                for hit_mz in ion_hits:
                    if abs(mz - hit_mz) <= self._HIT_TOLERANCE:
                        ion_hit_indices.add(i)
                        break

        table_frame = ttk.LabelFrame(
            self,
            text=f"Peak List ({len(mzs)} peaks, {len(ion_hit_indices)} ion hits)",
        )
        table_frame.pack(fill="both", expand=True, padx=8, pady=(2, 4))

        text = tk.Text(table_frame, wrap="none", height=8, font=("Courier", 9))
        yscroll = ttk.Scrollbar(table_frame, orient="vertical", command=text.yview)
        text.configure(yscrollcommand=yscroll.set)
        yscroll.pack(side="right", fill="y")
        text.pack(fill="both", expand=True)

        text.tag_configure(
            "ion_hit", foreground="#d62728", font=("Courier", 9, "bold"),
        )

        header = (
            f"  {'#':>4}  {'m/z':>14}  {'Intensity':>14}  "
            f"{'Ion Hit':>8}  {'Fragment':<20}\n"
        )
        text.insert("end", header)
        text.insert("end", "  " + "-" * 68 + "\n")

        for i, (mz, intensity) in enumerate(zip(mzs, intensities)):
            is_hit = i in ion_hit_indices
            mark = "  *" if is_hit else ""
            # Fragment column is populated only for hits with a known name.
            fragment = self._label_for_hit(mz) if is_hit else ""
            # Non-named hits return the m/z-as-label fallback, which we do
            # *not* want to repeat in the Fragment cell; strip it.
            if fragment and fragment == f"{mz:.2f}":
                fragment = ""
            line = (
                f"  {i+1:>4}  {mz:>14.4f}  {intensity:>14.1f}"
                f"{mark:>8}  {fragment:<20}\n"
            )
            if is_hit:
                text.insert("end", line, "ion_hit")
            else:
                text.insert("end", line)

        text.configure(state="disabled")

    def _build_text_fallback(self) -> None:
        """Show a text summary if matplotlib is not installed."""
        mzs, intensities, _ = self._parse_peaks()
        text = tk.Text(self, wrap="word", height=20)
        text.pack(fill="both", expand=True, padx=8, pady=(8, 4))
        text.insert("1.0", f"Spectrum data ({len(mzs)} peaks):\n\n")
        for mz, i in zip(mzs, intensities):
            text.insert("end", f"  m/z={mz:.4f}  I={i:.1f}\n")
        text.configure(state="disabled")


# ---------------------------------------------------------------------------
# Scoring Backfill Dialog
# ---------------------------------------------------------------------------

class ScoringBackfillDialog(tk.Toplevel):
    """Pre-check and configuration dialog for scoring backfill.

    Shows entry summary, checkboxes for Score A/B, file pickers,
    overwrite toggle, and a progress display during computation.

    Args:
        master: Parent window.
        db: Database manager.
        entries: List of entry dicts to score.
        session_user: Current session user.
        ion_list_path: Pre-filled ion list path from session settings.
        scoreb_wb_path: Pre-filled Score B workbook path from session settings.
        on_done: Callback(ion_path, wb_path) after scoring completes.
    """

    def __init__(
        self,
        master: CuratedLibraryWindow,
        db: CuratedLibraryDB,
        entries: list[dict],
        session_user: str,
        ion_list_path: str = "",
        scoreb_wb_path: str = "",
        on_done: callable | None = None,
    ) -> None:
        super().__init__(master)
        self.db = db
        self.entries = entries
        self.session_user = session_user
        self.on_done = on_done

        self.title("Scoring Backfill")
        self.geometry("520x420")
        self.resizable(False, False)
        self.transient(master)
        self.grab_set()

        # Pre-compute summary stats
        self._has_peaklist = sum(1 for e in entries if e.get("peaklist"))
        self._has_composition = sum(1 for e in entries if e.get("composition"))
        self._has_score_a = sum(1 for e in entries if e.get("score_a") is not None)
        self._has_score_b = sum(1 for e in entries if e.get("score_b") is not None)

        self._ion_path_var = tk.StringVar(value=ion_list_path)
        self._wb_path_var = tk.StringVar(value=scoreb_wb_path)
        self._do_score_a = tk.BooleanVar(value=True)
        self._do_score_b = tk.BooleanVar(value=bool(scoreb_wb_path))
        self._overwrite = tk.BooleanVar(value=False)

        self._build_ui()

    def _build_ui(self) -> None:
        n = len(self.entries)
        frame = ttk.Frame(self, padding=12)
        frame.pack(fill="both", expand=True)

        # Summary
        summary = ttk.LabelFrame(frame, text="Selection Summary")
        summary.pack(fill="x", pady=(0, 8))
        labels = [
            f"Selected entries: {n}",
            f"With peaklist: {self._has_peaklist}",
            f"With composition: {self._has_composition}",
            f"Already have Score A: {self._has_score_a}",
            f"Already have Score B: {self._has_score_b}",
        ]
        for txt in labels:
            ttk.Label(summary, text=txt).pack(anchor="w", padx=8, pady=1)

        # Options
        opts = ttk.LabelFrame(frame, text="Options")
        opts.pack(fill="x", pady=(0, 8))
        ttk.Checkbutton(opts, text="Compute Score A", variable=self._do_score_a).pack(
            anchor="w", padx=8, pady=2,
        )
        ttk.Checkbutton(opts, text="Compute Score B", variable=self._do_score_b).pack(
            anchor="w", padx=8, pady=2,
        )
        ttk.Checkbutton(opts, text="Overwrite existing scores", variable=self._overwrite).pack(
            anchor="w", padx=8, pady=2,
        )

        # File pickers
        files = ttk.LabelFrame(frame, text="Files")
        files.pack(fill="x", pady=(0, 8))

        r0 = ttk.Frame(files)
        r0.pack(fill="x", padx=8, pady=4)
        ttk.Label(r0, text="Ion list:").pack(side="left")
        ttk.Entry(r0, textvariable=self._ion_path_var, width=36).pack(
            side="left", padx=4, fill="x", expand=True,
        )
        ttk.Button(r0, text="Browse\u2026", command=self._browse_ion).pack(side="left")

        r1 = ttk.Frame(files)
        r1.pack(fill="x", padx=8, pady=4)
        ttk.Label(r1, text="Score B workbook:").pack(side="left")
        ttk.Entry(r1, textvariable=self._wb_path_var, width=30).pack(
            side="left", padx=4, fill="x", expand=True,
        )
        ttk.Button(r1, text="Browse\u2026", command=self._browse_wb).pack(side="left")

        # Progress
        self._progress_var = tk.StringVar(value="Ready")
        ttk.Label(frame, textvariable=self._progress_var).pack(
            anchor="w", pady=(4, 0),
        )

        # Buttons
        btn_frame = ttk.Frame(frame)
        btn_frame.pack(fill="x", pady=(8, 0))
        self._run_btn = ttk.Button(btn_frame, text="Run Scoring", command=self._run)
        self._run_btn.pack(side="right", padx=4)
        ttk.Button(btn_frame, text="Cancel", command=self.destroy).pack(
            side="right", padx=4,
        )

    def _browse_ion(self) -> None:
        p = filedialog.askopenfilename(
            title="Select Ion List",
            filetypes=[("CSV/Excel", "*.csv *.xlsx"), ("All files", "*.*")],
            parent=self,
        )
        if p:
            self._ion_path_var.set(p)

    def _browse_wb(self) -> None:
        p = filedialog.askopenfilename(
            title="Select Score B Workbook",
            filetypes=[("Excel", "*.xlsx *.xlsm"), ("All files", "*.*")],
            parent=self,
        )
        if p:
            self._wb_path_var.set(p)

    def _run(self) -> None:
        ion_path = self._ion_path_var.get().strip()
        wb_path = self._wb_path_var.get().strip()
        do_a = self._do_score_a.get()
        do_b = self._do_score_b.get()
        overwrite = self._overwrite.get()

        if do_a and not ion_path:
            messagebox.showwarning(
                "Missing", "Ion list path is required for Score A.", parent=self,
            )
            return
        if do_b and not wb_path:
            messagebox.showwarning(
                "Missing", "Score B workbook path is required for Score B.", parent=self,
            )
            return

        self._run_btn.configure(state="disabled")
        self._progress_var.set("Computing scores...")
        self.update_idletasks()

        try:
            results = compute_scores_batch(
                self.entries,
                ion_list_path=ion_path if do_a else None,
                workbook_path=wb_path if do_b else None,
                overwrite=overwrite,
                do_score_a=do_a,
                do_score_b=do_b,
            )
        except Exception as exc:
            messagebox.showerror("Scoring Error", f"Scoring failed:\n{exc}", parent=self)
            self._run_btn.configure(state="normal")
            self._progress_var.set("Failed")
            return

        # Write scores to DB
        new_a = 0
        new_b = 0
        for eid, scores in results:
            if scores.get("score_a") is not None:
                new_a += 1
            if scores.get("score_b") is not None:
                new_b += 1
            scores["last_modified_by"] = self.session_user
            # Update annotation_source if Score B was added
            entry = self.db.get_entry(eid)
            if entry and scores.get("score_b") is not None:
                src = entry.get("annotation_source", "")
                if "ScoreB" not in src:
                    scores["annotation_source"] = (src + "+ScoreB") if src else "scored"
            self.db.update_entry(eid, scores)

        skip_a = len(self.entries) - new_a if do_a else 0
        skip_b = len(self.entries) - new_b if do_b else 0
        no_comp = len(self.entries) - self._has_composition

        lines = []
        if do_a:
            lines.append(f"Score A: {new_a} new, {skip_a} skipped")
        if do_b:
            lines.append(f"Score B: {new_b} new, {skip_b} skipped ({no_comp} had no composition)")

        messagebox.showinfo("Scoring Complete", "\n".join(lines), parent=self)

        self.destroy()
        if self.on_done:
            self.on_done(ion_path, wb_path)


# ---------------------------------------------------------------------------
# Standalone launch (for development / testing)
# ---------------------------------------------------------------------------

def main() -> None:
    """Launch the Curated Library window standalone for development."""
    root = tk.Tk()
    root.withdraw()
    win = CuratedLibraryWindow(root)
    win.protocol("WM_DELETE_WINDOW", root.destroy)
    root.mainloop()


if __name__ == "__main__":
    main()
