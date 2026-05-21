import os
version = "1.09.3"
last_update = 20260521
import msprawextractor as mspext
import mzmlreader as mspmzmlext
import threading
from tkinter import ttk
import time #for testing
import random #for random failure
#import datetime
from datetime import datetime
import configparser #load config
import json
import tkinter.simpledialog as simpledialog
from tkinter import filedialog
from tkinter import messagebox
import shutil
import traceback
import mspvalidator_merger as mspval
import pandas as pd
# 20260520 fix B-28 (Patch D, revised per Codex same-day): opt out of pandas 3.x
# future.infer_string=True default. Mirrors the same opt-out in
# mspvalidator_merger.py top (Patch C) for defense in depth — if v14 is launched
# standalone or before mspval is imported, this still applies. Restores legacy
# 'object' dtype for string columns; required for v14's read_csv on trainable /
# prediction / converted CSVs and any downstream list-into-string assignment.
# Guard: tolerate older pandas where the option does not exist (see mspval Patch C).
try:
    pd.options.future.infer_string = False
except (AttributeError, KeyError):
    pass
import platform
import msp_insilicomarker_withGPT as marker
# --- Add near top-level imports ---
import os, joblib
import compnewv4 as compv4  # assumes dev/test calls are guarded by if __name__ == "__main__"
# Ion mining feature (safe to import even if file is absent)
# 20250907 to solve ion suggest missing issue
import importlib, sys, os, traceback, io
from pathlib import Path
import pathcanon
from typing import Dict
import re
import uuid
# ---- GlycoMSP startup diagnostics ---- #
import traceback, importlib.util

# If user double-clicked on Windows, stdin often isn't a real TTY. Stop from closing to allow error print
def _pause_if_no_tty():
    try:
        if os.name == "nt" and not sys.stdin.isatty():
            input("\nPress Enter to close…")
    except Exception:
        pass

def _check_mod(name, alt=None):
    """Return (present: bool, version_or_reason: str)."""
    modname = alt or name
    spec = importlib.util.find_spec(modname)
    if not spec:
        return False, "not found"
    try:
        m = __import__(modname)
        v = getattr(m, "__version__", "")
        return True, (v or "present")
    except Exception as e:
        return False, f"import failed: {e.__class__.__name__}: {e}"

def _quick_diag(note=None, file=None):
    f = file if file is not None else sys.stdout
    print("=== GlycoMSP quick environment check ===", file=f)
    print(f"Python: {sys.version.split()[0]}  ({sys.executable})", file=f)
    print(f"Platform: {platform.platform()}", file=f)
    if note:
        print(note, file=f)

    wanted = [
        ("tkinter", None),      # stdlib (bundled on python.org installers)
        ("numpy",   None),
        ("pandas",  None),
        ("joblib",  None),
        ("sklearn", "sklearn"), # scikit-learn
    ]
    for name, mod in wanted:
        ok, info = _check_mod(name, alt=mod)
        tag = "OK  " if ok else "MISS"
        shown = (mod or name)
        print(f"{tag} {shown:<10} {info}", file=f)

    if os.name == "nt":
        ok, info = _check_mod("pymsfilereader")
        tag = "OK  " if ok else "MISS"
        print(f"{tag} pymsfilereader {info}  (needed only for RAW conversion)", file=f)

    if sys.version_info >= (3, 13):
        print("\n[Notice] Running on Python 3.13.",
              "If wheels for numpy/pandas/scikit-learn are missing, prefer Python 3.11/3.12 for now.",
              file=f)

def _excepthook(exctype, value, tb):
    log_path = None
    try:
        if os.name == "nt":
            log_dir = Path(os.environ.get("LOCALAPPDATA", os.path.expanduser("~"))) / "GlycoMSP"
        else:
            log_dir = Path(os.path.expanduser("~")) / ".glycomsp"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / f"startup_error_{datetime.now():%Y%m%d_%H%M%S_%f}.log"
    except Exception:
        pass

    buf = io.StringIO()
    banner = "\n[StartupError] Unhandled exception during launch/import:\n"
    print(banner, file=sys.stdout)
    buf.write(banner + "\n")
    traceback.print_exception(exctype, value, tb, file=sys.stdout)
    traceback.print_exception(exctype, value, tb, file=buf)
    print(file=sys.stdout)
    buf.write("\n")
    _quick_diag(file=sys.stdout)
    _quick_diag(file=buf)

    if log_path is not None:
        try:
            log_path.write_text(buf.getvalue(), encoding="utf-8")
            print(f"\n[StartupError] Crash details saved to: {log_path}", file=sys.stdout)
        except Exception as e:
            print(f"\n[StartupError] Failed to save crash log: {e}", file=sys.stdout)

    _pause_if_no_tty()
    sys.exit(1)
sys.excepthook = _excepthook

# Allow a no-import diagnostics mode:
if "--diag" in sys.argv or "--diag" in sys.argv:
    _quick_diag(note="(ran with --diag; skipped importing the full GUI)")
    _pause_if_no_tty()
    sys.exit(0)

# ---- End diagnostics prelude ----
# =========================
# =========================
# =========================
# =========================
# ---- JSON type/schema enforcement (GlycoMSP v1) ----

JSON_TYPE_METADATA  = "glycomsp.metadata"
JSON_TYPE_METHOD    = "glycomsp.method"
JSON_TYPE_EXPERIMENT= "glycomsp.experiment"
SCHEMA_V1 = "1.0.0"
import hashlib
from typing import Optional

def file_sha256(path: str) -> Optional[str]:
    try:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1024*1024), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return None

def _now_iso_utc():
    # ISO8601 UTC string with seconds
    return datetime.utcnow().replace(microsecond=0).isoformat() + "Z"

def _ensure_header(d: dict, json_type: str):
    # For newly written JSONs: enforce the required header.
    d.setdefault("json_type", json_type)
    d.setdefault("schema_version", SCHEMA_V1)
    d.setdefault("created_utc", _now_iso_utc())
    d["updated_utc"] = _now_iso_utc()
    return d

def _legacy_classify_metadata(d: dict) -> bool:
    # Conservative: require at least the core metadata keys
    if not isinstance(d, dict):
        return False
    need = ["Experiment Title", "Glycan Type", "Mass Analyzer charge mode"]
    return all(k in d for k in need)

def _legacy_classify_method(d: dict) -> bool:
    # Legacy per-sample method JSON: {"experiment":..., "samples":{ sample: {csv/excel/metadata...}}}
    if not isinstance(d, dict):
        return False
    if "samples" not in d or not isinstance(d.get("samples"), dict) or not d["samples"]:
        return False
    # Inspect one sample entry
    _, sample_blob = next(iter(d["samples"].items()))
    if not isinstance(sample_blob, dict):
        return False
    # Legacy keys seen in v11: "csv","excel","metadata" (or "json" in some contexts)
    return ("csv" in sample_blob) and (("metadata" in sample_blob) or ("json" in sample_blob)) and ("excel" in sample_blob)

def _legacy_classify_method_pseudolabel(d: dict) -> bool:
    if not isinstance(d, dict):
        return False
    if d.get("dataset_type") != "pseudolabel":
        return False
    parents = d.get("parents")
    ionlist = d.get("ionlist")
    if not isinstance(parents, dict) or not isinstance(ionlist, dict):
        return False
    if not parents.get("converted_csv"):
        return False
    if not ionlist.get("path"):
        return False
    if not d.get("experiment_title") or not d.get("sample_name"):
        return False
    return True

def _legacy_classify_experiment(d: dict) -> bool:
    # Current exp.json writer uses {"experiment":..., "generated_on":..., "samples":{...}}
    if not isinstance(d, dict):
        return False
    return ("experiment" in d) and ("samples" in d) and isinstance(d.get("samples"), dict)

def _load_json_raw(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def load_typed_json(path: str, expected_type: str, *, allow_legacy: bool = True, context: str = "") -> dict:
    d = _load_json_raw(path)

    jt = d.get("json_type")
    sv = d.get("schema_version")

    # New-style: must match exactly
    if jt is not None or sv is not None:
        if jt != expected_type:
            raise ValueError(f"{context}Wrong json_type: expected '{expected_type}', got '{jt}' (file={path})")
        if not isinstance(sv, str) or not sv.strip():
            raise ValueError(f"{context}Missing/invalid schema_version (file={path})")
        return d

    # Legacy fallback (optional)
    if not allow_legacy:
        raise ValueError(f"{context}Missing json_type/schema_version (file={path})")

    if expected_type == JSON_TYPE_METADATA and _legacy_classify_metadata(d):
        d["json_type"] = JSON_TYPE_METADATA
        d["schema_version"] = SCHEMA_V1
        return d

    if expected_type == JSON_TYPE_METHOD:
        if _legacy_classify_method(d) or _legacy_classify_method_pseudolabel(d):
            d["json_type"] = JSON_TYPE_METHOD
            d["schema_version"] = SCHEMA_V1
            return d

    if expected_type == JSON_TYPE_EXPERIMENT and _legacy_classify_experiment(d):
        d["json_type"] = JSON_TYPE_EXPERIMENT
        d["schema_version"] = SCHEMA_V1
        return d

    raise ValueError(f"{context}Unrecognized legacy JSON; cannot safely classify (file={path})")

# ---- json schema ends ---- #
# ===========================

# a fail guardsafe import of msp_ion_mining
# For both MAS and CGA(supported now?)? It is called in dataset preparation workflow, do we really need early-assignment here?
def _load_ion_module():
    try:
        return importlib.import_module("msp_ion_mining")
    except Exception as e:
        # 2nd chance: try alongside this file and one level up
        here = os.path.dirname(__file__)
        for cand in (os.path.join(here, "msp_ion_mining.py"),
                     os.path.join(os.path.dirname(here), "msp_ion_mining.py")):
            if os.path.exists(cand):
                try:
                    spec = importlib.util.spec_from_file_location("msp_ion_mining", cand)
                    mod = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(mod)
                    sys.modules["msp_ion_mining"] = mod
                    return mod
                except Exception as e2:
                    print("[IonSuggest] import failed from", cand, "→", repr(e2))
                    traceback.print_exc()
        print("[IonSuggest] import failed:", repr(e))
        traceback.print_exc()
        return None

_ionmod = _load_ion_module()
if _ionmod:
    export_ion_suggestions_csv = _ionmod.export_ion_suggestions_csv
    SuggestParams = _ionmod.SuggestParams
else:
    export_ion_suggestions_csv = None
    class SuggestParams: pass
# ===========================

# 20260520 fix B-30: tkinter filedialog macOS-safety helpers.
# Tk on macOS crashes (SIGABRT via setAllowedFileTypes / NSSavePanel nil-UTI)
# when filetypes patterns or defaultextension use compound double-dot extensions
# (*.exp.json, *.method.json) because no UTI exists for them. Helpers collapse
# compound patterns to bare *.json on macOS and post-process SaveAs returns to
# preserve the .exp.json / .method.json naming convention.
# Per Codex round-2 verdict on packet codex_B30_macos_saveasdialog_review.md:
#  - REPLACE (don't append) compound patterns with safe ones; fallback entries
#    do not rescue because Tk crashes on the first bad UTI conversion.
#  - Dedupe collapsed patterns preserving first label.
#  - Skip multi-pattern strings ("*.xls *.xlsx") via " " not in pattern guard.
#  - Both Save and Open dialogs need patching; _macos_append_compound_suffix is
#    SAVE-ONLY (never call on askopenfilename returns — would mutate user-picked
#    existing-file path).
# Windows/Linux behavior is unchanged (helpers are pass-through).
_IS_MACOS = platform.system() == "Darwin"


def _macos_safe_filetypes(filetypes):
    """Replace compound double-dot patterns (e.g. '*.exp.json', '*.method.json')
    with '*.json' on macOS, then dedupe by pattern preserving the first label.
    Pass-through on Windows/Linux. Returns a new list.

    Skips multi-pattern strings like '*.xls *.xlsx' which are not compound
    suffixes but multiple patterns in one Tk pattern string.
    """
    if not _IS_MACOS:
        return filetypes
    safe = []
    seen_patterns = set()
    for label, pattern in filetypes:
        if (
            isinstance(pattern, str)
            and " " not in pattern
            and pattern.startswith("*.")
            and pattern.count(".") > 1
        ):
            new_pattern = "*.json"
        else:
            new_pattern = pattern
        if new_pattern in seen_patterns:
            continue
        seen_patterns.add(new_pattern)
        safe.append((label, new_pattern))
    return safe


def _macos_safe_defaultextension(ext):
    """Collapse compound double-dot defaultextension to '.json' on macOS.
    SaveAs callers should pair this with _macos_append_compound_suffix on the
    returned path. Pass-through on Windows/Linux or empty input."""
    if not _IS_MACOS or not ext:
        return ext
    if ext.count(".") > 1:
        return ".json"
    return ext


def _macos_append_compound_suffix(path, compound_suffix):
    """SAVE-ONLY. If macOS user-returned SaveAs path doesn't end in the compound
    suffix, append it (handling the 'user picked .json' case by upgrading to
    the compound form). No-op on Windows/Linux or empty path. compound_suffix
    should start with '.' (e.g. '.exp.json').

    Per Codex round-2 verdict on B-30: never call on askopenfilename returns
    — the user selected a specific existing file and the path must not be
    rewritten. Use only on asksaveasfilename() returns.
    """
    if not _IS_MACOS or not path:
        return path
    if path.endswith(compound_suffix):
        return path
    if path.endswith(".json"):
        return path[:-5] + compound_suffix
    return path + compound_suffix
# ===========================

# Changelogs:
# v1.5 (future) allow multiple methods exist under one sample (need 1.2 update first to satisfy requirements)
# v1.3 (future) start cleaning unneeded code blocks, move changelog to wiki and other versionfiles.
# v1.2 (future) fix the tree selection/display logic (would probably bundled with v1.1 update)
# v1.1 [stable version for publish] fix the old macos crash issue in tkinter, confirmed full pipeline executable on both Windows (fixed package version) and MacOS (py3.12+, latest packages) 
# v1.09.5 refactoring v1.09.23 fix conversion early fire issue. If conversion at Thermo COM level fails, need to close app to release. Marked for future fix.
# v1.09.3 deep review and add comments for further refactor work (Claude Code involved)
# v1.09.1 dead code cleanup
# v1.09 (manuscript version): CGA score b implemented
# v1.07: CGA extra score b prepared.
# v1.05: fixed json file relationship definition and support legacy json file load (v11). 
# v1.02: json fix B (now metadata, method, and experiment json should have proper relationship)
# v1.01: json fix A (compatible with legacy, previous mixed json file)
# v1.00: Able to write manuscript although some bug persists.
# v0.9999: fix minor bugs and display issues
# v0.9998: really fix gate issue (confirmed on manual datasets)
# v0.9996: fix gate not working on ones learned with protonated mass
# v0.9995: not sure if gate fixed completely
# v0.99924: allow prediction results having 1. normal 2. mass gated (need theo mass from in silico csv) 3. learned w/ protonated mass 4. 2+3
# v0.9993: allow protonated mass getting learned
# v0.99921: allow trainable PL dataset retain protonated mass (so next version, in ML we can treat it as feature)
# v0.9992: add precursor -composition mass check 
# v0.9991: fix minor bugs
# v0.999: fixed PL unlabel dataset functionality
# v0.9987: hide the old version of converting PL datasets to trainable data in ML analysis tab (intermediate function. can be removed safely)
# v0.9985-86: fix the negative mode pseudolabeling issue (the in silico glycan list was fine, now fixing PL workflow itself)
# v0.9984: fix the issue when user cwd (output directory for file conversion) falls to default when directly open the program by clicking v10.py
# v0.9983: add Pseudolabeling method (generate in silico glycan list) support on negative mode and GlcA (HexA) - testing 
# v0.9981: add ML parameters save/load feature. Add error tracker for configuring first time.
# v0.998: update requirements.txt (the python version and packages needs to be updated). PS. python 3.13 has errors
# v0.997: add integrity check placeholder, fix the path issue
# v0.9969: fix pseudolabeling metadata logics
# v0.9967 try to fix batch conversion issue
# v0.9966: incorporate random sampling manner also in PL trainable dataset generation 
# v0.9965: finished GUI pseudolabeling -> trainable csv 
# v0.9963: fix OG in pseudolabeling
# v0.996: Add ion (feature) mining method 
# v0.993: able to apply pseudolabeling function (functional but may have bugs)
# v0.9925: add split and stratified method when creating training set
# v0.992: Link to pseudolabeling function
# v0.991: Make the software functional to work on MacOS
# v0.99: Add Glypick-like autoannotation back (need to change UI)
# v0.91: ML added. Lacking combining data and include "Non-glycan labels for training"
# v0.9: add ML window
# v0.81: adding derivatization flag
# v0.8: able to export merged trainable datasets
# v0.65: introduce real validation function
# v0.64c: pseudo-link (no real validation) but link the files
# v0.64b: load method file and skip unneeded file conversion for filling missing metadata part
# v0.64a: record current assignment and link information (in one method file or update json and keep sample json in method?)
# v0.63: dataset management + link and validate data (placeholder)
# v0.61: build functional dataset management window (if new version failed, rollback to this)
# v0.6: support saving converted files by selecting user defined folders
# for trace errors
# v0.58: UI refine. (v0.57 works perfectly)
# v0.56: trace errors and also delete temp csv while metadata assingment and rename csv process failed
# v0.56  change the logic to avoid main window freeze during conversion
# v0.55: refine the logic. Go to simulation before real run.
# version: 0.54
# add persistent metadata edit window
# version: 0.53
# add metadata class
# version: 0.52
# add quit destroy function
# version: 0.51
# add logger
# test of adding tkinter windows with the help of GPT4o
# add GUI, add progress bar and making buttons work

# version: 0.4
# added mzml support function (waiting mzml extractor implementation - will be noted as v0.5 then, while MSP will be v0.6)
# version: 0.3
# added excel reader for annotated file. Will add validation here or in annotationreader.py
# version: 0.2
# date: 20240730
# about this file: original trytolistoutheaders.py and GlycoMSP_demo.py modulated to support future standardized processing workflow
# import tkinter part
#mind that we plan to add GUI support in project managing, the tkinter detection may be moved to funtion: GUI in MSPinit.py in future
# re-organized the code structure and readability using copilot and ChatGPT4-o
# no sensitive contents were sent to the server for this part
# ============================

# tkinter-safe init
try: 
    import tkinter as tk
    from tkinter import filedialog, messagebox
except ImportError:
    raise ImportError('Please install tkinter to enable GUI-based file selection')
else:
    has_tkinter = True

# raw/mzml file conversion init
selected_files = {}  # Dictionary to keep track of selected files

#20250907 preventing key error in GUI #need real test to see behavior changes
FILETYPE_TO_KEY = {
    "csv": "csv",
    "excel": "excel",
    "json": "json",        # metadata as json(legacy)
    "method": "json",      # method json stored as "json"
    "metadata": "metadata" # metadata json stored as "metadata"
}
def normalize_ftype(ft: str) -> str:
    return FILETYPE_TO_KEY.get(ft.lower().strip(), ft.lower().strip())

# 20260413 Code review: This is Placeholder #
# --- Integrity checker (optional) ---
try:
    import mspprojintegritykeeper as mspik
except Exception:
    # Soft fallback so GUI still works if the module is missing.
    from types import SimpleNamespace
    def _noimpl(*a, **k): return None
    mspik = SimpleNamespace(
        writehashtojson=_noimpl,
        file_hashcheck=_noimpl,
        relocatemissingfile=_noimpl,
    )
# =================

# App Logger block
class AppLogger:
    def __init__(self):
        self.entries = []
        self.debug_logs = []
        self.gui_writer = None  # Optional live display hook

    def log(self, message):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        entry = f"[{timestamp}] {message}"
        print(entry)
        self.entries.append(entry)
        if self.gui_writer:
            self.gui_writer(entry)

    def debug(self, message):
        self.debug_logs.append(f"[DEBUG] {message}")

    def save(self, filename, include_debug=False):
        with open(filename, "w", encoding="utf-8") as f: #add utf-8 encode to get arrow supported
            for line in self.entries:
                f.write(line + "\n")
            if include_debug and self.debug_logs:
                f.write("-----\ndebug enabled-----\n")
                for line in self.debug_logs:
                    f.write(line + "\n")
                f.write("-----\ndebug log ends-----\n")
        return filename

logger = AppLogger()
# =====================

#20250822 replace composition window with pseudolabel window
# --- Pseudo-Labeling Setup window (replaces with CGASetupWindow) ---
class CGASetupWindow(tk.Toplevel):
    """
    Edit in-silico generation flags and preview essential metadata.
    - Defaults injected via `default_flags` (dict copied before editing)
    - Read-only metadata summary (Glycan Type / Charge / Derivatization)
    - Save/Load flags presets (JSON)
    - Calls `on_submit({"flags": ..., "metadata": {...}})` on Generate
    """
    def __init__(self, master,
            meta_json_path=None, meta_prefill=None, default_flags=None, on_submit=None,
            on_generate=None, on_link_existing=None, on_attach_ionlist=None,
            on_attach_scoreb_workbook=None, on_start=None,
            initial_insilico=None, initial_ionlist=None, initial_scoreb_workbook=None):

        super().__init__(master)
        self.title("Constraint-based Glycan Annotation Setup")
        self.geometry("900x800")
        self.minsize(840, 720)
        self.resizable(True, True)
        self.meta_json_path = meta_json_path
        self.on_submit = on_submit
        self.flags = dict(default_flags or {})  # shallow copy; safe to mutate locally
        self.on_generate = on_generate
        self.on_link_existing = on_link_existing
        self.on_attach_ionlist = on_attach_ionlist
        self.on_attach_scoreb_workbook = on_attach_scoreb_workbook
        self.on_start = on_start
        self.insilico_path_var = tk.StringVar(value=initial_insilico or "")
        self.ionlist_path_var  = tk.StringVar(value=initial_ionlist or "")
        self.scoreb_workbook_path_var = tk.StringVar(value=initial_scoreb_workbook or "")
        # --- Metadata (editable) ---
        meta_frame = ttk.LabelFrame(self, text="Metadata (can override here)")
        meta_frame.pack(fill="x", padx=12, pady=(12, 6))

        ## score B validator helper
        def _validate_scoreb_workbook(path: str):
            """
            Part 1 UI validation only:
            - file exists
            - Excel extension
            - optional deep validation via ScoreBLoader if import is available
            Returns: (ok: bool, message: str)
            """
            path = (path or "").strip()
            if not path:
                return False, "No Score B workbook selected."

            if not os.path.exists(path):
                return False, f"Workbook not found:\n{path}"

            ext = os.path.splitext(path)[1].lower()
            if ext not in {".xlsx", ".xlsm", ".xltx", ".xltm"}:
                return False, "Score B workbook must be an Excel workbook (.xlsx/.xlsm)."

            # 20260518 B-12 fold-in (issue c): split import-try from validation-try (was silently downgrading real ImportError from transitive deps)
            try:
                from msp_CGA_structscore import ScoreBLoader
            except ImportError:
                # Fallback if module import is not ready in this environment yet
                return True, "Workbook file looks acceptable (deep validation unavailable in this build)."
            try:
                # Preferred: deep workbook validation from your Score B module
                ScoreBLoader(path).load()
                return True, "Workbook validated successfully."
            except Exception as e:
                return False, f"Score B workbook validation failed:\n{e}"

        self._validate_scoreb_workbook = _validate_scoreb_workbook

        # hold current + original for a Reset action
        self.meta_vars = {
            "Glycan Type": tk.StringVar(value=""),
            "Mass Analyzer charge mode": tk.StringVar(value=""),
            "Derivatization Type": tk.StringVar(value=""),
        }
        self._meta_original = {"Glycan Type":"", "Mass Analyzer charge mode":"", "Derivatization Type":""}

        # 1) prefill from dict if provided (preferred)
        if meta_prefill:
            self.meta_vars["Glycan Type"].set(meta_prefill.get("Glycan Type", ""))
            self.meta_vars["Mass Analyzer charge mode"].set(meta_prefill.get("Mass Analyzer charge mode", ""))
            self.meta_vars["Derivatization Type"].set(meta_prefill.get("Derivatization Type", ""))
        for k in self._meta_original:
            self._meta_original[k] = self.meta_vars[k].get()

        # render editable controls r= row
        r = 0
        ttk.Label(meta_frame, text="Glycan Type:").grid(row=r, column=0, sticky="w", padx=10, pady=4)
        ttk.Combobox(meta_frame, state="readonly", width=10,
                    values=["N", "O"],
                    textvariable=self.meta_vars["Glycan Type"]).grid(row=r, column=1, sticky="w", padx=8, pady=4)
        r += 1

        ttk.Label(meta_frame, text="Mass Analyzer charge mode:").grid(row=r, column=0, sticky="w", padx=10, pady=4)
        ttk.Combobox(meta_frame, state="readonly", width=10,
                    values=["+", "-"],
                    textvariable=self.meta_vars["Mass Analyzer charge mode"]).grid(row=r, column=1, sticky="w", padx=8, pady=4)
        r += 1

        ttk.Label(meta_frame, text="Derivatization Type:").grid(row=r, column=0, sticky="w", padx=10, pady=4)
        ttk.Combobox(meta_frame, state="readonly", width=20,
                    values=["PerMe", "PerMe(Freeend)", "PerMe(Reduced)", "None"],
                    textvariable=self.meta_vars["Derivatization Type"]).grid(row=r, column=1, sticky="w", padx=8, pady=4)
        # tiny reset link
        def _reset_meta():
            for k, v in self._meta_original.items():
                self.meta_vars[k].set(v)

        ttk.Button(meta_frame, text="Reset from file", command=_reset_meta).grid(row=0, column=2, rowspan=3, padx=8)
        #20250910
        # --- Flag editor (dynamic by glycan type) ---
        flags_frame = ttk.LabelFrame(self, text="In-Silico Generation Flags")
        flags_frame.pack(fill="both", expand=True, padx=12, pady=6)

        # Full spec (we'll subset by glycan type)
        flag_spec = {
            # monitoring / flow
            "debug": ("bool", None),
            "dev": ("bool", None),
            "force_exit": ("bool", None),
            # composition options
            "alphagal_like": ("bool", None),
            "allowldnc": ("bool", None),
            "allowleby": ("bool", None),
            "allow5ac": ("bool", None),
            "allow5gc": ("bool", None),
            "allowkdn": ("bool", None),
            "allowfuc": ("bool", None),
            "allowpsa": ("int", (0, 3)),
            "allowldnf": ("bool", None),
            #20250925 extra
            # >>> NEW: HexA / sulfate / phosphate <<<
            "allowHexA": ("bool", None),
            "HexA_range": ("range", (0, 4)),   # 0–4 is a safe UI cap; adjust if you like
            "SO3": ("bool", None),
            "SO3_range": ("range", (0, 2)),
            "PO3H": ("bool", None),
            "PO3H_range": ("range", (0, 2)),

            # iteration logic
            "arm_count": ("int", (0, 8)),
            "internal_minrep": ("int", (0, 6)),
            "internal_maxrep": ("int", (0, 10)),
            "topology": ("bool", None),
            # NG-only core flags (hide for O)
            "corefuc": ("bool", None),
            "bicorefuc": ("bool", None),
            "highman": ("bool", None),
            "perman": ("bool", None),
            "hybrid": ("bool", None),
            # optional composition check
            "compcheck": ("bool", None),
            "Hex_range": ("range", (0, 24)),
            "HexNAc_range": ("range", (0, 20)),
            "Neu5Ac_range": ("range", (0, 10)),
            "Neu5Gc_range": ("range", (0, 10)),
            "KDN_range": ("range", (0, 10)),
            "Fucose_range": ("range", (0, 10)),

        }

        # Which keys to show per glycan type
        NG_KEYS = {
            "debug","dev","force_exit",
            "alphagal_like","allowldnc","allowleby","allow5ac","allow5gc","allowkdn","allowfuc","allowpsa","allowldnf",
            # NEW:
            "allowHexA","HexA_range","SO3","SO3_range","PO3H","PO3H_range",
            #
            "arm_count","internal_minrep","internal_maxrep","topology",
            "corefuc","bicorefuc","highman","perman","hybrid",
            "compcheck","Hex_range","HexNAc_range","Neu5Ac_range","Neu5Gc_range","KDN_range","Fucose_range"
        }
        OG_KEYS = {
            "debug","dev","force_exit","alphagal_like","allowldnc","allowleby",
            "allow5ac","allow5gc","allowkdn","allowfuc","allowpsa",
            # NEW:
            "allowHexA","HexA_range","SO3","SO3_range","PO3H","PO3H_range",
            #"arm_count",
            "internal_minrep","internal_maxrep","topology",
            "compcheck","Hex_range","HexNAc_range","Neu5Ac_range","Neu5Gc_range","KDN_range","Fucose_range"
        }
        self.flag_vars = {}
        # Define the keys group for the CGA window
        # Keys that belong to the "Advanced options" collapsible group
        ADVANCED_KEYS = {"debug", "dev", "force_exit", "topology", "perman"}
        # Keys that belong to the "glyco-feature restraints" sub-frame
        GLYCO_FEATURE_KEYS = {
            "alphagal_like","allowldnc","allowleby","allow5ac","allow5gc",
            "allowkdn","allowfuc","allowpsa","allowHexA", "allowldnf"
        }
        INTERNAL_EXT_KEYS = {"arm_count", "internal_minrep", "internal_maxrep"}
        NG_CORE_KEYS = {"corefuc", "bicorefuc", "highman", "hybrid"}
        SUBSTITUENT_NEG_KEYS = {"SO3", "SO3_range", "PO3H", "PO3H_range"}
        COMPOSITION_KEYS = {
            "HexA_range", "Hex_range", "HexNAc_range",
            "Neu5Ac_range", "Neu5Gc_range", "KDN_range",
            "Fucose_range", "compcheck"
        }
        # controls visibility of "Advanced options" panel
        self.show_advanced_var = tk.BooleanVar(value=False)

        left = ttk.Frame(flags_frame)
        right = ttk.Frame(flags_frame)

        advanced_frame = ttk.LabelFrame(flags_frame, text="Advanced options")

        def _coerce_int_like(val, default=0):
            return int(val) if isinstance(val, (int, float, str)) and str(val).strip() != "" else int(default)

        def add_bool(parent, key, row):
            var = tk.BooleanVar(value=bool(self.flags.get(key, False)))
            ttk.Checkbutton(parent, text=key, variable=var).grid(row=row, column=0, sticky="w", pady=3)
            self.flag_vars[key] = var
            if key == "hybrid":
                ttk.Label(parent, text="(components auto-filled on submit)").grid(row=row, column=1, sticky="w")

        def add_int(parent, key, row, lo, hi):
            label = key
            if key == "arm_count":
                label = "arm_count (N-glycans)"
            ttk.Label(parent, text=label + ":").grid(row=row, column=0, sticky="w")
            var = tk.IntVar(value=_coerce_int_like(self.flags.get(key, 0), 0))
            ttk.Spinbox(parent, from_=lo, to=hi, textvariable=var, width=6).grid(row=row, column=1, sticky="w", padx=6)
            self.flag_vars[key] = var

        def add_range(parent, key, row, lo, hi):
            ttk.Label(parent, text=key + ":").grid(row=row, column=0, sticky="w")
            default = self.flags.get(key, [0, 0])
            vmin = tk.IntVar(value=_coerce_int_like(default[0] if isinstance(default, (list, tuple)) else 0, 0))
            vmax = tk.IntVar(value=_coerce_int_like(default[1] if isinstance(default, (list, tuple)) else 0, 0))
            wrap = ttk.Frame(parent)
            wrap.grid(row=row, column=1, sticky="w")
            ttk.Spinbox(wrap, from_=lo, to=hi, textvariable=vmin, width=5).pack(side="left")
            ttk.Label(wrap, text=" to ").pack(side="left")
            ttk.Spinbox(wrap, from_=lo, to=hi, textvariable=vmax, width=5).pack(side="left")
            self.flag_vars[key] = (vmin, vmax)

        def _current_keys():
            gtype = self.meta_vars["Glycan Type"].get().strip().upper()
            return NG_KEYS if gtype == "N" else OG_KEYS

        #move OG core panel here to avoid called before assignment exceptions
        # --- O-glycan core types (multi-select) ---
        self.og_core_vars = {i: tk.BooleanVar(value=False) for i in (0, 1, 2, 3, 4)}
        def _selected_coretypes():
            sel = [i for i, v in self.og_core_vars.items() if v.get()]
            return sel if sel else [1, 2, 3, 4]  # sensible default

        def _snapshot_flags():
            snap = {}
            # flags
            for k, w in self.flag_vars.items():
                if isinstance(w, tuple):
                    val = (int(w[0].get()), int(w[1].get()))
                elif isinstance(w, tk.BooleanVar):
                    val = bool(w.get())
                else:
                    val = int(w.get())
                snap[k] = val
                # keep latest GUI value so we can restore even when widgets are rebuilt/hidden
                self.flags[k] = val
            # OG core types (added below)
            snap["_ogcore"] = {i: v.get() for i, v in self.og_core_vars.items()}
            return snap

        def _restore_flags(snap):
            if not snap: return
            for k, w in self.flag_vars.items():
                if k not in snap: continue
                v = snap[k]
                if isinstance(w, tuple):
                    w[0].set(int(v[0] if isinstance(v, (list, tuple)) else 0))
                    w[1].set(int(v[1] if isinstance(v, (list, tuple)) else 0))
                elif isinstance(w, tk.BooleanVar):
                    w.set(bool(v))
                else:
                    w.set(int(v))
            for i, val in snap.get("_ogcore", {}).items():
                if i in self.og_core_vars:
                    self.og_core_vars[i].set(bool(val))

        def _rebuild_flag_panel(*_):
            # preserve state
            snap = _snapshot_flags()
            # clear frames
            for child in left.winfo_children():
                child.destroy()
            for child in right.winfo_children():
                child.destroy()
            for child in advanced_frame.winfo_children():
                child.destroy()
            self.flag_vars.clear()

            # choose keys for current glycan type
            keys = [k for k in flag_spec.keys() if k in _current_keys()]

            # split into glyco-feature / main / advanced
            glyco_feature_keys = [k for k in keys if k in GLYCO_FEATURE_KEYS]
            internal_keys = [k for k in keys if k in INTERNAL_EXT_KEYS]
            ngcore_keys = [k for k in keys if k in NG_CORE_KEYS]
            substituent_neg_keys = [k for k in keys if k in SUBSTITUENT_NEG_KEYS]
            composition_keys = [k for k in keys if k in COMPOSITION_KEYS]
            main_keys = [
                k for k in keys
                if (
                    k not in ADVANCED_KEYS
                    and k not in GLYCO_FEATURE_KEYS
                    and k not in INTERNAL_EXT_KEYS
                    and k not in NG_CORE_KEYS
                    and k not in SUBSTITUENT_NEG_KEYS
                    and k not in COMPOSITION_KEYS
                )
            ]
            adv_keys = [k for k in keys if k in ADVANCED_KEYS]

            # remaining main keys go into left/right columns
            half = (len(main_keys) + 1) // 2
            left_keys, right_keys = main_keys[:half], main_keys[half:]

            # generic renderer
            def render_column(parent, keys_subset, start_row=0):
                r = start_row
                for k in keys_subset:
                    ftype, extra = flag_spec[k]
                    if ftype == "bool":
                        add_bool(parent, k, r)
                    elif ftype == "int":
                        lo, hi = extra
                        add_int(parent, k, r, lo, hi)
                    elif ftype == "range":
                        lo, hi = extra
                        add_range(parent, k, r, lo, hi)
                    r += 1
                return r

            # --- glyco-feature restraints sub-frame (left side) ---
            left_start_row = 0
            # row 0: glyco-feature restraints
            if glyco_feature_keys:
                glyco_frame = ttk.LabelFrame(left, text="glyco-feature restraints")
                glyco_frame.grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 4))
                render_column(glyco_frame, glyco_feature_keys, start_row=0)
                left_start_row = 1 # other left-column flags start under this frame

            # row 1: substituent(neg)
            if substituent_neg_keys:
                subneg_frame = ttk.LabelFrame(left, text="substituent(neg)")
                subneg_frame.grid(row=left_start_row, column=0, columnspan=2, sticky="w", pady=(0, 4))
                render_column(subneg_frame, substituent_neg_keys, start_row=0)
                left_start_row += 1

            # --- internal extension/elongation sub-frame (right side) ---
            right_start_row = 0
            # row 0: internal extension/elongation
            if internal_keys:
                internal_frame = ttk.LabelFrame(right, text="internal extension/elongation")
                internal_frame.grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 4))
                render_column(internal_frame, internal_keys, start_row=0)
                right_start_row = 1  # next content goes to row 1+
            # row 1: NG Core Flags (N-glycan only keys, but safe to render normally)
            if ngcore_keys:
                ngcore_frame = ttk.LabelFrame(right, text="NG Core Flags")
                ngcore_frame.grid(row=right_start_row, column=0, columnspan=2, sticky="w", pady=(0, 4))
                render_column(ngcore_frame, ngcore_keys, start_row=0)
                right_start_row += 1
            # row 2: composition restraints
            if composition_keys:
                comp_frame = ttk.LabelFrame(right, text="composition restraints")
                comp_frame.grid(row=right_start_row, column=0, columnspan=2, sticky="w", pady=(0, 4))
                render_column(comp_frame, composition_keys, start_row=0)
                right_start_row += 1
                # OG core types should appear under the right-side stack (after composition restraints)
                # O-glycan core types (only render when Glycan Type == O)
                if self.meta_vars["Glycan Type"].get().strip().upper() == "O":
                    ogcore_frame = ttk.LabelFrame(right, text="O-glycan core types (select 1–4)")
                    ogcore_frame.grid(
                        row=right_start_row,
                        column=0,
                        columnspan=2,
                        sticky="w",
                        pady=(0, 4)
                    )

                    core_labels = {
                        0: "Tn Antigen",
                        1: "Core 1",
                        2: "Core 2",
                        3: "Core 3",
                        4: "Core 4",
                    }

                    r = 0
                    for i in (0, 1, 2, 3, 4):
                        ttk.Checkbutton(
                            ogcore_frame,
                            text=core_labels[i],
                            variable=self.og_core_vars[i]
                        ).grid(row=r // 2, column=r % 2, sticky="w", padx=6)
                        r += 1

                    right_start_row += 1

            # main columns (remaining flags)
            render_column(left, left_keys, start_row=left_start_row)
            render_column(right, right_keys, start_row=right_start_row)

            # advanced panel (optional)
            if self.show_advanced_var.get() and adv_keys:
                if not advanced_frame.winfo_manager():
                    advanced_frame.pack(fill="x", padx=10, pady=(0, 8))
                render_column(advanced_frame, adv_keys, start_row=0)
            else:
                advanced_frame.pack_forget()

            # restore state
            _restore_flags(snap)

        # After creating the Glycan Type combobox (named via self.meta_vars["Glycan Type"])
        gly_cb = meta_frame.grid_slaves(row=0, column=1)[0]  # the Combobox you just created
        gly_cb.bind("<<ComboboxSelected>>", _rebuild_flag_panel)

        # Toggle to show/hide advanced flags
        advanced_toggle = ttk.Checkbutton(
            flags_frame,
            text="Show advanced options (debug /future options)",
            variable=self.show_advanced_var,
            command=_rebuild_flag_panel,
        )
        advanced_toggle.pack(anchor="w", padx=12, pady=(0, 4))

        # main flag columns sit under the toggle
        left.pack(side="left", fill="both", expand=True, padx=(10, 5), pady=8)
        right.pack(side="left", fill="both", expand=True, padx=(5, 10), pady=8)

        # Initial render
        _rebuild_flag_panel()

        # Small panel to show current links (insilico / ion list)
        links_frame = ttk.LabelFrame(self, text="Files selected for CGA analysis")
        links_frame.pack(fill="x", padx=12, pady=(6, 0))

        # allow the path column to expand
        links_frame.columnconfigure(1, weight=1)

        ttk.Label(links_frame, text="In-silico glycan list (compositions):").grid(row=0, column=0, sticky="w", padx=10, pady=4)
        ttk.Label(
            links_frame,
            textvariable=self.insilico_path_var,
            justify="left",
            wraplength=760  # adjust if you change geometry width
        ).grid(row=0, column=1, sticky="we", padx=8, pady=4)

        ttk.Label(links_frame, text="Fragmentation ion list (features):").grid(row=1, column=0, sticky="w", padx=10, pady=4)
        ttk.Label(
            links_frame,
            textvariable=self.ionlist_path_var,
            justify="left",
            wraplength=760
        ).grid(row=1, column=1, sticky="we", padx=8, pady=4)

        #newly added score b
        ttk.Label(links_frame, text="Score B reference:").grid(row=2, column=0, sticky="w", padx=10, pady=4)
        ttk.Label(
            links_frame,
            textvariable=self.scoreb_workbook_path_var,
            justify="left",
            wraplength=760
        ).grid(row=2, column=1, sticky="we", padx=8, pady=4)

        # --- Footer buttons (top row: preparation) ---
        btns_top = ttk.Frame(self)
        btns_top.pack(fill="x", padx=8, pady=(6, 2))
        ttk.Button(btns_top, text="Load CGA settings", command=self.load_flags).pack(side="left", padx=4)
        ttk.Button(btns_top, text="Save CGA settings", command=self.save_flags).pack(side="left", padx=4)

        ttk.Separator(btns_top, orient="vertical").pack(side="left", fill="y", padx=8)
        #20250910
        def _on_generate():
            payload = {"flags": self.collect_flags(),
                    "metadata": {
                        "Glycan Type": self.meta_vars["Glycan Type"].get(),
                        "Mass Analyzer charge mode": self.meta_vars["Mass Analyzer charge mode"].get(),
                        "Derivatization Type": self.meta_vars["Derivatization Type"].get(),
                        "_meta_json": self.meta_json_path or "",
                        "_overrides_applied": any(self.meta_vars[k].get() != self._meta_original.get(k, "")
                                                    for k in self._meta_original)
                    },
                    "coretype": _selected_coretypes(),   # NEW
                    }
            cb = self.on_submit or self.on_generate
            if cb:
                maybe_path = cb(payload)
                if isinstance(maybe_path, str) and os.path.exists(maybe_path):
                    self.insilico_path_var.set(maybe_path)

        ttk.Button(btns_top, text="Generate In-Silico CSV", command=_on_generate).pack(side="left", padx=4)

        # Link existing in-silico (no auto popups unless user clicks)
        def _link_existing():
            p = filedialog.askopenfilename(title="Select existing in-silico CSV",
                                        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
            if not p:
                return
            self.insilico_path_var.set(p)
            cb = self.on_link_existing
            
            if cb:
                cb(p)

        ttk.Button(btns_top, text="Link in silico glycan list", command=_link_existing).pack(side="left", padx=4)

        # Attach ion list (optional)
        def _attach_ionlist():
            p = filedialog.askopenfilename(title="Attach Ion List (CSV/XLSX)",
                                        filetypes=[("CSV/XLSX", "*.csv;*.xlsx;*.xls"), ("All files", "*.*")])
            if not p:
                return
            self.ionlist_path_var.set(p)
            cb = self.on_attach_ionlist
            if cb:
                cb(p)

        ttk.Button(btns_top, text="Add Ion List", command=_attach_ionlist).pack(side="left", padx=4)

        #score b
        def _attach_scoreb_workbook():
            p = filedialog.askopenfilename(
                title="Select Score B workbook",
                filetypes=[("Excel workbooks", "*.xlsx *.xlsm *.xltx *.xltm"), ("All files", "*.*")]
            )
            if not p:
                return

            self.scoreb_workbook_path_var.set(p)
            cb = self.on_attach_scoreb_workbook
            if cb:
                cb(p)

        def _validate_scoreb_selected():
            ok, msg = self._validate_scoreb_workbook(self.scoreb_workbook_path_var.get())
            if ok:
                messagebox.showinfo("Score B Workbook", msg)
            else:
                messagebox.showerror("Score B Workbook", msg)

        def _clear_scoreb_workbook():
            self.scoreb_workbook_path_var.set("")
            cb = self.on_attach_scoreb_workbook
            if cb:
                cb("")

        ttk.Button(btns_top, text="Add Score B Workbook", command=_attach_scoreb_workbook).pack(side="left", padx=4)
        ttk.Button(btns_top, text="Validate Score B Workbook", command=_validate_scoreb_selected).pack(side="left", padx=4)
        ttk.Button(btns_top, text="Clear Score B Workbook", command=_clear_scoreb_workbook).pack(side="left", padx=4)
        ttk.Separator(btns_top, orient="vertical").pack(side="left", fill="y", padx=8)

        # --- Footer buttons (bottom row: execution) ---
        btns_bottom = ttk.Frame(self)
        btns_bottom.pack(fill="x", padx=8, pady=(2, 10))

        # Start CGA analysis (enabled if converted CSV + insilico present)
        #modified with score b
        def _start():
            scoreb_path = self.scoreb_workbook_path_var.get().strip()

            # Optional validation at launch time
            if scoreb_path:
                ok, msg = self._validate_scoreb_workbook(scoreb_path)
                if not ok:
                    messagebox.showerror("Score B Workbook", msg)
                    return

            payload = {
                "flags": self.collect_flags(),
                "metadata": {
                    "Glycan Type": self.meta_vars["Glycan Type"].get(),
                    "Mass Analyzer charge mode": self.meta_vars["Mass Analyzer charge mode"].get(),
                    "Derivatization Type": self.meta_vars["Derivatization Type"].get(),
                    "_meta_json": self.meta_json_path or "",
                    "_overrides_applied": any(
                        self.meta_vars[k].get() != self._meta_original.get(k, "")
                        for k in self._meta_original
                    ),
                },
                "coretype": _selected_coretypes(),
                "insilico_csv": self.insilico_path_var.get().strip(),
                "ionlist_path": self.ionlist_path_var.get().strip(),
                "score_b_workbook_path": scoreb_path,
            }
            cb = self.on_start
            if cb:
                cb(payload)

        ttk.Button(btns_bottom, text="Start CGA analysis", command=_start).pack(side="left", padx=2)
        ttk.Button(btns_bottom, text="Close", command=self.destroy).pack(side="left", padx=4)


    def collect_flags(self):
        out = {}
        for k, v in self.flag_vars.items():
            if isinstance(v, tuple):
                vmin, vmax = v
                out[k] = [int(vmin.get()), int(vmax.get())]
            elif isinstance(v, tk.BooleanVar):
                out[k] = bool(v.get())
            else:
                out[k] = int(v.get())
        # 20260518 B-13: persist OG core types (previously lost on save/load)
        out["coretype"] = [i for i, v in self.og_core_vars.items() if v.get()]
        return out

    def load_flags(self):
        path = filedialog.askopenfilename(title="Load Flags JSON", filetypes=[("JSON", "*.json")])
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
            for k, v in data.items():
                if k not in self.flag_vars:
                    continue
                w = self.flag_vars[k]
                if isinstance(w, tuple):
                    w[0].set(int(v[0] if isinstance(v, (list, tuple)) else 0))
                    w[1].set(int(v[1] if isinstance(v, (list, tuple)) else 0))
                elif isinstance(w, tk.BooleanVar):
                    w.set(bool(v))
                else:
                    w.set(int(v))
            # 20260518 B-13: restore OG core types (gracefully absent in legacy saves)
            cores = data.get("coretype", [])
            for i, v in self.og_core_vars.items():
                v.set(i in cores)
            messagebox.showinfo("Flags Loaded", f"Loaded: {os.path.basename(path)}")
        except Exception as e:
            messagebox.showerror("Load Failed", str(e))

    def save_flags(self):
        path = filedialog.asksaveasfilename(title="Save Flags JSON",
                                            defaultextension=".json",
                                            filetypes=[("JSON", "*.json")])
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(self.collect_flags(), f, indent=2)
            messagebox.showinfo("Saved", f"Saved: {os.path.basename(path)}")
        except Exception as e:
            messagebox.showerror("Save Failed", str(e))

# Constraint-based Glycan Annotation Setup Window Ends
# ===========================

# Generic Tk widget liveness check (used by MetadataEditorWindow + reserved for future windows)
def widget_alive(w):
    try:
        return (w is not None) and int(w.winfo_exists()) == 1
    except Exception:
        return False

# Metadata Editor (Shows when conversion of raw/mzml is finished, or filling missing metadata in prepare dataset window)
class MetadataEditorWindow:
    def __init__(self, parent, raw_file_list, on_each_metadata_ready_callback, on_finish=None, output_dir=None, skip_conversion=False, input_kind=None):
        self.parent = parent
        self.raw_file_list = raw_file_list
        self.callback = on_each_metadata_ready_callback
        self.on_finish = on_finish
        self.current_index = 0
        self.last_metadata = None
        self.metadata = {}
        self.output_dir = output_dir
        self.skip_conversion = skip_conversion
        self.input_kind = input_kind
        self.entries = {}
        self.files = {}   # ← add this line in 20250916
        fields = [
            "Experiment Title", "Experiment Description", "Author running this analysis",
            "Raw data acquired date", "Glycan Type", "Mass Analyzer charge mode", "Derivatization Type"
        ]

        self.window = tk.Toplevel(parent)
        # 20250917
        self._convert_in_progress = False          # NEW: track background conversion
        self._error_mode = False            # ← NEW: allow close after failures
        self._batch_halted = False          # ← NEW: treat batch as stopped after error

        self.window.protocol("WM_DELETE_WINDOW", self._on_user_close)  # NEW: intercept user close (X)

        self.window.title("Enter Experiment Metadata")
        self.window.geometry("600x500")
        # Row 0: Status label placeholder (spans both columns)
        self.status_label = tk.Label(self.window, text="", fg="blue")
        self.status_label.grid(row=0, column=0, columnspan=2, pady=(10, 5))
        for i, field in enumerate(fields, start=1):
            label = tk.Label(self.window, text=field)
            label.grid(row=i, column=0, padx=10, pady=5, sticky="w")

            if field == "Mass Analyzer charge mode":
                var = tk.StringVar()
                entry = ttk.Combobox(self.window, textvariable=var, values=["+", "-"], state="readonly", width=10)
                entry.current(0)
            elif field == "Derivatization Type":
                var = tk.StringVar()
                entry = ttk.Combobox(
                    self.window, textvariable=var,
                    values=["PerMe(Freeend)", "PerMe(Reduced)", "Others"],  # Add more if needed
                    state="readonly", width=20
                )
                entry.current(0)

            elif field == "Glycan Type":
                var = tk.StringVar()
                entry = ttk.Combobox(self.window, textvariable=var, values=["N", "O", "N+O"], state="readonly", width=10) #add the values when going to support new glycan types
                entry.current(0)
            else:
                entry = tk.Entry(self.window, width=50)

            entry.grid(row=i, column=1, padx=10, pady=5)
            self.entries[field] = entry
        #20260414 some minor fixes 
        if not self.skip_conversion:
            self.load_file(self.raw_file_list[self.current_index])

        # Buttons
        button_frame = tk.Frame(self.window)
        button_row = len(fields) + 1 #so the button adjusted itself in future if we add more contents into metadata
        button_frame.grid(row=button_row, column=0, columnspan=2, pady=10)
        tk.Button(button_frame, text="Import Metadata", command=self.import_metadata).grid(row=0, column=0, padx=5)
        tk.Button(button_frame, text="Use Last Metadata", command=self.use_last_metadata).grid(row=0, column=1, padx=5)
        tk.Button(button_frame, text="Clear", command=self.clear_fields).grid(row=0, column=2, padx=5)
        tk.Button(button_frame, text="Generate Dataset", command=self.generate).grid(row=0, column=3, padx=5)
        tk.Button(button_frame, text="Cancel", command=self._cancel_editor).grid(row=0, column=4, padx=5)  # 20260519 [B-26]: was self.window.destroy (orphaned tmps + bar)

        label = tk.Label(self.window, text="Raw File (optional):")
        label.grid(row=button_row + 1, column=0, padx=10, pady=5, sticky="w")

        raw_file_frame = tk.Frame(self.window)
        raw_file_frame.grid(row=button_row + 1, column=1, padx=10, pady=5, sticky="w")

        self.raw_file_entry = tk.Entry(raw_file_frame, width=40)
        self.raw_file_entry.pack(side=tk.LEFT)

        def select_raw_file():
            path = filedialog.askopenfilename(title="Select RAW File", filetypes=[("RAW files", "*.raw")])
            if path:
                self.raw_file_entry.delete(0, tk.END)
                self.raw_file_entry.insert(0, path)

        tk.Button(raw_file_frame, text="Browse", command=select_raw_file).pack(side=tk.LEFT, padx=5)

    def _cancel_editor(self):
        """20260519 [B-26] Unified cancel handler: cleanup mzml tmps + invoke on_finish + destroy window.
        Why: bare window.destroy from the Cancel button left ms2tmp_/ms3tmp_ orphaned (mzml path)
        and the indeterminate progress bar animating indefinitely (both paths). on_finish ==
        reset_main_status callback wired at MetadataEditorWindow construction, which calls
        progress.stop() + sets status='Idle' + re-enables the convert button.
        """
        try:
            for p in [self.files.get("ms2tmp"), self.files.get("ms3tmp")]:
                if p and os.path.exists(p):
                    try:
                        os.remove(p)
                        # 20260519 [B-26] log success — diagnostic value for Henry + testers (especially when current extraction overwrote a pre-existing orphan from a prior session, then this cleanup removes it)
                        logger.log(f"[Cleanup] Removed leftover temp file (cancel path): {p}")
                    except OSError as e:
                        logger.log(f"[Cleanup] cancel: failed to remove {p}: {e}")
        except Exception:
            pass
        try:
            if callable(getattr(self, "on_finish", None)):
                self.on_finish()
        except Exception:
            pass
        try:
            self.window.destroy()
        except Exception:
            pass

    #20251001 replace
    def on_conversion_complete(self):
        # Status (your existing line)
        self._safe_set_status("Conversion complete. Ready for metadata.", fg="green")
        # Safely reflect current raw path into the (readonly) Entry if it still exists
        try:
            entry = getattr(self, "raw_file_entry", None)
            if entry and entry.winfo_exists():
                entry.config(state="normal")
                entry.delete(0, "end")
                entry.insert(0, self.current_raw_file or "")
                entry.config(state="readonly")
        except tk.TclError:
            # Widget already destroyed; nothing else to do here
            return

        # Re-enable entries and buttons (guard each widget in case the window is closing)
        try:
            for entry_widget in getattr(self, "entries", {}).values():
                if entry_widget and entry_widget.winfo_exists():
                    entry_widget.config(state="normal")
            for child in self.window.winfo_children():
                if isinstance(child, tk.Button) and child.winfo_exists():
                    child.config(state="normal")
        except tk.TclError:
            # Window or widgets may be gone; abort cleanly
            return

        # Pre-fill if last used (as before)
        if getattr(self, "last_metadata", None):
            try:
                self.fill_fields_from_metadata(self.last_metadata)
            except Exception:
                # prefill is best-effort; continue even if it fails
                pass

        # Set title if window still exists
        try:
            if self.window and self.window.winfo_exists():
                self.window.title(f"Metadata for: {os.path.basename(self.current_raw_file)}")
        except tk.TclError:
            return

    def load_file(self, raw_file):
        self._error_mode = False        # ← NEW
        self._batch_halted = False          # ← NEW: treat batch as stopped after error
        self.current_raw_file = raw_file
        self.window.title(f"Converting: {os.path.basename(raw_file)}")
        self.clear_fields()

        #self.status_label.config(text="Converting raw file... please wait.", fg="blue")
        self._safe_set_status("Converting raw file... please wait.", fg="blue")

        # Disable form and buttons temporarily
        for entry in self.entries.values():
            entry.config(state="disabled")
        for child in self.window.winfo_children():
            if isinstance(child, tk.Button):
                child.config(state="disabled")

        # 20250927
        def background_conversion():
            if self.skip_conversion:
                return

            self._convert_in_progress = True
            try:
                # ask extractor to write into the chosen output folder (if any)
                if getattr(self, "input_kind", "raw") == "mzml":
                    res = mspmzmlext.extract_mzML(raw_file,
                                                  outdir=self.output_dir,
                                                  debug=True)
                else:
                    res = mspext.convert_raw_to_csv(
                        raw_file,
                        outdir=self.output_dir,
                        debug=True
                    )

                def _norm_path(p: str | None) -> str | None:
                    """Make absolute; if it's a bare name, assume output_dir (else CWD)."""
                    if not p:
                        return None
                    # If caller returned a bare filename, anchor it
                    if not os.path.dirname(p):
                        base_dir = self.output_dir or os.getcwd()
                        p = os.path.join(base_dir, p)
                    return os.path.abspath(p)

                # Accept either (ms2, ms3) or {"ms2tmp": ..., "ms3tmp": ...}
                if isinstance(res, dict):
                    tmp_ms2 = _norm_path(res.get("ms2tmp"))
                    tmp_ms3 = _norm_path(res.get("ms3tmp"))
                elif isinstance(res, (list, tuple)) and len(res) >= 2:
                    tmp_ms2 = _norm_path(res[0])
                    tmp_ms3 = _norm_path(res[1])
                else:
                    raise RuntimeError("Unexpected return from convert_raw_to_csv")

                # Remember where the temps are for finalize()
                self.files["ms2tmp"] = tmp_ms2
                self.files["ms3tmp"] = tmp_ms3
                logger.log(f"[convert] MS2 temp exists immediately after conversion: {os.path.exists(tmp_ms2) if tmp_ms2 else None}")
                logger.log(f"[convert] MS3 temp exists immediately after conversion: {os.path.exists(tmp_ms3) if tmp_ms3 else None}")

                self._convert_in_progress = False
                self._safe_after(0, self.on_conversion_complete)

            except Exception as e:
                self._convert_in_progress = False
                logger.log(f"[FATAL] Exception during raw file conversion: {e}")
                # 20260518 B-10: marshal to Tk main thread (was direct call from worker)
                self._safe_after(0, lambda: self._enter_error_mode("An error occurred. Please check the log."))

                try:
                    for p in [self.files.get("ms2tmp"), self.files.get("ms3tmp")]:
                        if p and os.path.exists(p):
                            os.remove(p)
                            logger.log(f"[Cleanup] Removed leftover temp file: {p}")
                except Exception as ce:
                    logger.log(f"[WARNING] Cleanup failed: {ce}")

                self._safe_after(
                    0,
                    lambda err=e: messagebox.showerror(
                        "Thermo Library Error",
                        f"Raw file could not be processed.\n\nDetails:\n{err}"
                    )
                )
                self._safe_after(0, self._update_main_status_for_batch_state)
                self._safe_after(500, self._attempt_close_after_error)
            finally:
                # double-guard (in case of unexpected code paths)
                self._convert_in_progress = False
                self._safe_after(0, self._update_main_status_for_batch_state)

        threading.Thread(target=background_conversion, daemon=True).start()
        
    def import_metadata(self):
        path = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
        if path:
            with open(path, "r") as f:
                self.metadata = json.load(f)
            self.fill_fields_from_metadata(self.metadata)

    def use_last_metadata(self):
        if self.last_metadata:
            self.fill_fields_from_metadata(self.last_metadata)
        else:
            messagebox.showinfo("No previous metadata", "You haven't submitted any metadata yet.")

    def fill_fields_from_metadata(self, meta):
        for field, entry in self.entries.items():
            entry.delete(0, tk.END)
            entry.insert(0, meta.get(field, ""))

    def clear_fields(self):
        for entry in self.entries.values():
            entry.delete(0, tk.END)

    def _safe_set_status(self, text, fg=None):
        try:
            if widget_alive(getattr(self, "status_label", None)):
                if fg is None:
                    self.status_label.config(text=text)
                else:
                    self.status_label.config(text=text, fg=fg)
        except Exception as e:
            logger.log(f"[UI] status update skipped: {e}")

    def _safe_after(self, ms, func, *args, **kwargs):
        host = getattr(self, "parent", None) or getattr(self, "window", None)
        if widget_alive(host):
            host.after(ms, func, *args, **kwargs)
    #20250917
    
    def _is_batch_active(self) -> bool:
        """
        True if a conversion thread is running OR more files remain in this batch.
        """
        #to allow main status update when batch error
        if getattr(self, "_batch_halted", False):   # ← NEW
            return False
        more_pending = self.current_index < (len(self.raw_file_list) - 1)
        return bool(self._convert_in_progress or more_pending)

    def _on_user_close(self):
        """
        User clicked the window X. If batch is still running, block close and
        update main status. Otherwise, allow close and set main status Idle.
        """
        #early close if the error happens in batch processing
        if getattr(self, "_error_mode", False):
            # After a failure, always allow the editor to close
            self._safe_after(0, self._update_main_status_for_batch_state)
            try:
                self.window.destroy()
            except Exception:
                pass
            return
                                                                                
        if self._is_batch_active():                                                   
            if messagebox.askyesno(                                                   
                "Conversion Running",                                                 
                "Conversion is still in progress. Force close?\n\n"                   
                "(The background thread will be abandoned. No files will be corrupted.)"):                                                    
                self._convert_in_progress = False                                     
                self._batch_halted = True                                             
                try:
                    self.window.destroy()                                             
                except Exception:                                 
                    pass
                # 20260519 [B-26] stop indeterminate progress bar after force-close; status="Idle (cancelled)" below preserves yellow indicator
                try:
                    progress.stop()
                except Exception:
                    pass
                set_main_status("Idle (cancelled)", fg="orange")
                return                                                                
            else:
                self._safe_set_status("Waiting for conversion to finish...",          
        fg="orange")                                                                  
                return

        # else: safe to close
        set_main_status("Idle", fg="blue")
        try:
            self.window.destroy()
        except Exception:
            pass         
    def _attempt_close_after_error(self):
        """
        Called after showing an error. If batch is NOT active, close the window.
        If batch is active, keep it open (blocked by policy above) but make the state clear.
        """
        #Auto-close after error regardless of batch
        if getattr(self, "_error_mode", False):
            # Failure overrides the batch-close lock; return to main window
            self._safe_after(0, self._update_main_status_for_batch_state)
            if widget_alive(getattr(self, "window", None)):
                try:
                    self.window.destroy()
                except Exception:
                    pass
            return

        if self._is_batch_active():
            set_main_status("Batch conversion running", fg="blue")
            self._safe_set_status("Conversion failed. Batch is still running — editor locked.", fg="red")
        else:
            set_main_status("Idle", fg="blue")
            if widget_alive(getattr(self, "window", None)):
                try:
                    self.window.destroy()
                except Exception:
                    pass        

    def _enter_error_mode(self, msg="Conversion failed."):
        self._error_mode = True
        self._convert_in_progress = False
        self._batch_halted = True                      # ← NEW
        self._safe_set_status(f"{msg} You can close this editor.", fg="red")
        try:
            self.window.title(f"Error: {os.path.basename(self.current_raw_file)}")
        except Exception:
            pass
        # flip main window NOW
        self._safe_after(0, self._update_main_status_for_batch_state)

    def _update_main_status_for_batch_state(self):
        try:
            more_pending = self.current_index < (len(self.raw_file_list) - 1)
        except Exception:
            more_pending = False
        if getattr(self, "_batch_halted", False):           # ← NEW
            set_main_status("Idle (error)", fg="red")
        elif self._convert_in_progress or more_pending:
            set_main_status("Batch conversion running", fg="blue")
        else:
            set_main_status("Idle", fg="blue")
    #20250927
    def generate(self):
        if not self.skip_conversion:
            if self._convert_in_progress:
                messagebox.showwarning("Please wait", "Conversion is still running.")
                return

            if not self.files.get("ms2tmp") or not self.files.get("ms3tmp"):
                messagebox.showwarning(
                    "Conversion not finished",
                    "Temporary conversion files are not ready yet. Please wait until conversion completes."
                )
                logger.log("[generate] blocked: temp paths not ready yet")
                return

        if self.output_dir is None:
            selected_dir = filedialog.askdirectory(
                title="Select folder to save metadata and output files",
                initialdir=os.getcwd()
            )
            if not selected_dir:
                messagebox.showwarning("Cancelled", "You must select an output folder.")
                return
            self.output_dir = selected_dir
            logger.log(f"[Metadata] Output will be saved to: {self.output_dir}")


        # collect metadata
        for field, entry in self.entries.items():
            self.metadata[field] = entry.get()

        extractdate = time.strftime("%Y%m%d_%H%M%S")
        self.metadata["Parameters when GlycoMSP launched"] = ["Autofilled by GUI"]
        self.metadata["Date of file extracted from raw file"] = extractdate

        # resolve raw path
        raw_path = ""
        if hasattr(self, "raw_file_entry"):
            raw_path = (self.raw_file_entry.get() or "").strip()
        if not raw_path:
            raw_path = getattr(self, "current_raw_file", "") or ""
        raw_path = raw_path or "not linked"
        self.metadata["Original raw file path"] = raw_path
        self.metadata["Raw filename"] = (
            os.path.splitext(os.path.basename(raw_path))[0]
            if raw_path != "not linked" else "(not linked)"
        )

        mansavename = simpledialog.askstring("Project Name", "Enter a name for this batch/project:")
        if not mansavename:
            messagebox.showwarning("Missing Name", "You must enter a project/batch name.")
            return

        raw_base = os.path.splitext(os.path.basename(raw_path))[0] if raw_path != "not linked" else "no_raw"
        savename = f"{mansavename}_{extractdate}_{raw_base}"

        # save metadata json
        output_json_path = os.path.join(self.output_dir, savename + ".json")
        # fixed 20260305 P1-C: msprawextractor.savemetadata() does NOT stamp headers, so we must do it here
        # Verified: msprawextractor.savemetadata() and mzmlreader.savemetadata() both do plain json.dump with no key injection.
        # The prior comment claiming "already implemented in mzmlreader and masprawextractor" was incorrect.
        # enforce typed metadata json header (overwrite if wrong / missing)
        self.metadata["json_type"] = "glycomsp.metadata"
        self.metadata["schema_version"] = "1.0.0"
        jsonfile = mspext.savemetadata(self.metadata, output_json_path)
        logger.log(f"Metadata saved to {jsonfile}")
        self.last_metadata = self.metadata.copy()

        # final file targets
        final_ms2 = os.path.join(self.output_dir, f"ms2_{savename}.csv")
        final_ms3 = os.path.join(self.output_dir, f"ms3_{savename}.csv")

        success = True

        if self.skip_conversion:
            logger.log("[debug] filling missing metadata only; no csv rename/check")
        else:
            try:
                Path(final_ms2).parent.mkdir(parents=True, exist_ok=True)

                # prefer absolute paths captured during background_conversion
                temp_ms2 = self.files.get("ms2tmp")
                temp_ms3 = self.files.get("ms3tmp")
                logger.log(f"[generate] project name = {mansavename}")
                logger.log(f"[generate] savename = {savename}")
                logger.log(f"[generate] final_ms2 = {final_ms2}")
                logger.log(f"[generate] final_ms3 = {final_ms3}")
                logger.log(f"[generate] temp_ms2 = {temp_ms2}")
                logger.log(f"[generate] temp_ms3 = {temp_ms3}")
                logger.log(f"[generate] ms2 exists before promote = {os.path.exists(temp_ms2) if temp_ms2 else None}")
                logger.log(f"[generate] ms3 exists before promote = {os.path.exists(temp_ms3) if temp_ms3 else None}")
                # fallback: construct likely locations
                raw_stem = os.path.splitext(os.path.basename(raw_path))[0] if raw_path != "not linked" else "no_raw"
                if not temp_ms2:
                    temp_ms2 = os.path.join(self.output_dir, f"ms2tmp_{raw_stem}.csv") #cand->temp_ms2 avoid cwd guess for user clicking msp
                    #temp_ms2 = cand if os.path.exists(cand) else os.path.abspath(f"ms2tmp_{raw_stem}.csv")
                if not temp_ms3:
                    temp_ms3 = os.path.join(self.output_dir, f"ms3tmp_{raw_stem}.csv") #same as ms2
                    #temp_ms3 = cand if os.path.exists(cand) else os.path.abspath(f"ms3tmp_{raw_stem}.csv")

                # promote temps → finals
                try:
                    promote_temp_to_final(temp_ms2, final_ms2, logger)
                except Exception as e:
                    logger.log(f"[WARNING] MS2 finalize failed: {e}")
                    success = False

                try:
                    promote_temp_to_final(temp_ms3, final_ms3, logger)
                except Exception as e:
                    logger.log(f"[WARNING] MS3 finalize failed: {e}")
                    success = False

                if success:
                    self.files["ms2"] = os.path.abspath(final_ms2)
                    self.files["ms3"] = os.path.abspath(final_ms3)
                    self.files.pop("ms2tmp", None)
                    self.files.pop("ms3tmp", None)

            except Exception as e:
                logger.log(f"[FATAL] Finalize failed: {e}")
                success = False

            if success:
                mspext.finalize_extraction(self.current_raw_file, self.metadata, savename)
                logname = logger.save(os.path.join(self.output_dir, savename + ".log"), include_debug=True)
                logger.log(f"Saved log to: {logname}")

            # advance batch
            self.current_index += 1
            if self.current_index < len(self.raw_file_list):
                self.load_file(self.raw_file_list[self.current_index])
            else:
                messagebox.showinfo("Metadata", "All metadata have been completed.")
                self.window.destroy()
                if self.on_finish:
                    self.on_finish()

        # callback after metadata is created
        if self.callback:
            self.callback(raw_path, self.metadata, savename)
        # If we're only generating metadata (no conversion/batch), close immediately
        if getattr(self, "skip_conversion", False):
            try:
                print("[Should close window]20260126 GPT suggested fix on metadata/method writing close winodow fix. Not sure if it works")
                self.window.destroy()
            except Exception:
                print("[Exception happened when closing window after method data generation]20260126 GPT suggested fix on metadata/method writing close winodow fix. Not sure if it works")
                pass
            return
   
# MetadataEditor Window Ends here
# ==============================

# Decouple of AppLogger from tkinter (intended, for future CLI support cases etc)
def write_to_gui(message):
    text_widget.insert(tk.END, message + "\n")
    text_widget.see(tk.END)

logger.gui_writer = write_to_gui
# AppLogger calls using tkinter ends
# ===============================

#version info reader
#will add same one in extractor but may depreciate after moving metadata editing part to GUI
#20260412 code-review: not active for now. About window carries the work. Considering CHANGE this to CLI or other clickable elements in about window in future)
def get_version_info_for(module_name, manifest_path=r".\project_version_manifest.ini"):
    config = configparser.ConfigParser()
    config.read(manifest_path)
    print(f"current path: {manifest_path}")

    if module_name in config:
        version = config[module_name].get("version", "N/A")
        last_update = config[module_name].get("last_update", "N/A")
        return version, last_update
    return None, None

# Data Loading helpers : csv (globally) and ion list (for CGA score A and ...?)
#newly added 20250824 for reading csv (is that essential?)
def robust_read_csv(path, prefer_tab=False):
    """
    Try several parsing strategies (TSV first if prefer_tab=True).
    Returns a pandas DataFrame or raises the last error.
    """
    import csv

    tries = []

    # Prefer TSV if requested (your converted file)
    if prefer_tab:
        tries += [
            dict(engine="python", sep="\t", encoding="utf-8-sig", low_memory=False),
            dict(engine="c",      sep="\t", encoding="utf-8-sig", low_memory=False),
        ]

    # Standard CSV attempts
    tries += [
        dict(engine="c",      encoding="utf-8-sig", low_memory=False),
        dict(engine="python", sep=None, encoding="utf-8-sig", low_memory=False),  # sniff delimiter
        dict(engine="python", sep=r',(?=(?:[^"]*"[^"]*")*[^"]*$)', encoding="utf-8-sig", low_memory=False),
        dict(engine="python", delimiter=",", quoting=csv.QUOTE_NONE, escapechar="\\",
             encoding="utf-8-sig", on_bad_lines="skip", low_memory=False),
    ]

    last_err = None
    for kw in tries:
        try:
            if kw.get("engine") == "python":
                kw.pop("low_memory", None)
            df = pd.read_csv(path, **kw)
            df.columns = [str(c).strip() for c in df.columns]
            return df
        except Exception as e:
            last_err = e
    raise last_err

# 20260420 promoted to module-level utility
#It strips whitespace and removes the Unicode BOM character (\ufeff) from all column headers
def clean_cols(df):
    if df is None or df.empty: 
        return df
    df = df.copy()
    df.columns = [str(c).strip().replace("\ufeff","") for c in df.columns]
    return df

#fix in 20251001
def _read_ion_df(path):
    """Read ion list from .xlsx/.xls (preferred) or .csv (comma-separated)."""
    if not path:
        return None
    p = str(path).lower()
    if p.endswith((".xlsx", ".xls")):
        try:
            sheets = pd.read_excel(path, sheet_name=None)
            # pick sheet named like "ionlist" first, else the first with a mass-like column
            preferred = None
            for name in sheets:
                if name.strip().lower() in {"ionlist", "ions", "ion_list"}:
                    preferred = sheets[name]
                    break
            if preferred is None:
                for df in sheets.values():
                    cols_l = {c.strip().lower() for c in df.columns}
                    if any(c in cols_l for c in {"mass", "mz", "ion_mz", "m/z"}):
                        preferred = df
                        break
            if preferred is None:
                return None
            df = preferred
        except ImportError as e:
            messagebox.showerror(
                "Missing dependency",
                "Reading .xlsx needs 'openpyxl'.\n\nPlease install:\n\npip install openpyxl"
            )
            raise            
    else:
        try:
            df = pd.read_csv(path, engine="python")
        except Exception:
            raise ValueError(f"Unsupported ion-list file: {path}")

    # normalize a mass column name
    col_map = {c.lower(): c for c in df.columns}
    for key in ("mass", "mz", "ion_mz", "m/z"):
        if key in col_map:
            if key != "mass":
                df = df.rename(columns={col_map[key]: "mass"})
            break
    return df[["mass"]].dropna() if "mass" in df.columns else None
# Data loading helper ends
# =======================

# About Window showing the version info
def open_about_window():
    about_win = tk.Toplevel(root)
    about_win.title("About GlycoMSP")
    about_win.geometry("500x400")

    tk.Label(about_win, text="About GlycoMSP GUI", font=("Arial", 14, "bold")).pack(pady=10)
    tk.Label(about_win, text="Made by Huan-Chuan TSENG", font=("Arial", 10)).pack(pady=5)
    
    # Create a scrollable text box
    about_text = tk.Text(about_win, wrap=tk.WORD, height=20, width=60)
    about_text.pack(padx=10, pady=5, fill="both", expand=True)

    scrollbar = tk.Scrollbar(about_text)
    scrollbar.pack(side="right", fill="y")
    about_text.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=about_text.yview)

    # Load version info from manifest, mind the path issues (on demo folder now)
    config = configparser.ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__),"project_version_manifest.ini"))

    if "version_check" in config:
        ts = config["version_check"].get("checked_on", "Unknown")
        about_text.insert(tk.END, f"Version manifest generated on: {ts}\n\n")

    for section in config.sections():
        if section != "version_check":
            version = config[section].get("version", "N/A")
            updated = config[section].get("last_update", "N/A")
            about_text.insert(tk.END, f"[{section}]\n  Version: {version}\n  Last Update: {updated}\n\n")

    about_text.config(state=tk.DISABLED)

    tk.Label(about_win, text="Component version listed above", font=("Arial", 10)).pack(pady=5)
    tk.Button(about_win, text="Close", command=about_win.destroy).pack(pady=10)
# About Window Ends
# ===================

# Function to validate file path, called when selecting raw/mzml
def validate_file_path(file_path):
    if not os.path.exists(file_path):
        print(f"Error: The file {file_path} does not exist.")
        return False
    if not os.path.isfile(file_path):
        print(f"Error: The path {file_path} is not a file.")
        return False
    return True

# Main GUI-related action
def update_display():
    text_widget.delete(1.0, tk.END)
    for ftype, path in selected_files.items():
        text_widget.insert(tk.END, f"{ftype}: {path}\n")

def clear_files():
    selected_files.clear()
    update_display()

def reset_main_status():
    status_var.set("Idle")
    progress.stop()
    convert_button.config(state="normal")
# Main GUI-related action ends
# ========================

# Spectral file processing for downstream csv-based pre-processing
def select_file(filetype):
    filetypes_dict = {
        "raw": [("Raw file", "*.raw")],
        "mzml": [("mzML file", "*.mzML")],
    }
    filepaths = filedialog.askopenfilenames(filetypes=filetypes_dict.get(filetype, [("All files", "*.*")]))
    if filepaths:
        valid_paths = [fp for fp in filepaths if validate_file_path(fp)]
        selected_files.setdefault(filetype, []).extend(valid_paths)
        update_display()
    elif not filepaths:
        print(f"No valid {filetype} spectral file(s) selected.")

# path resolve
def same_drive(a: str, b: str) -> bool:
    da = os.path.splitdrive(os.path.abspath(a))[0].lower()
    db = os.path.splitdrive(os.path.abspath(b))[0].lower()
    return da == db

# Finalize the ms2 and ms3 csv from temp name to proper exp + timestamp output
def promote_temp_to_final(src: str, dst: str, logger):
    """Promote temp→final.
    - If same drive/volume: os.replace (atomic).
    - Else: copy2 then remove src (Windows cross-drive).
    """
    dst = os.path.abspath(dst)
    src = os.path.abspath(src)
    Path(dst).parent.mkdir(parents=True, exist_ok=True)

    if not os.path.exists(src):
        raise FileNotFoundError(f"Temp file missing: {src}")

    if same_drive(src, dst):
        os.replace(src, dst)
        logger.log(f"[finalize] Promoted (atomic) {src} → {dst}")
    else:
        shutil.copy2(src, dst)
        os.remove(src)
        logger.log(f"[finalize] Promoted (copy+delete) {src} → {dst} (cross-drive)")

# Spectral file processing block Ends
# ========================


# Metadata-related blocks
# Callback to log medatata is ready detection when prep the conversion
# 20260414 code review: add savename display and suggestion from Claude Code
def on_metadata_ready(raw_file, metadata, savename):
    logger.log(f"Confirmed metadata for {raw_file} → {savename}")
    # metadata is a full dict — not useful in a one-liner log. Leave it unused but in the signature. 
    # It's there for future use (e.g., auto-assigning to experiment_projects).     

def launch_metadata_batch():
    # need to rebuild to accept mzML
    # allow mzML conversion without pymsreader. (change the hierarchy to adapt)
    """
    Entry point for the 'Convert Raw to csv' button.

    Current behavior:
      - RAW: existing workflow (metadata batch + mspext raw conversion)
      - mzML: placeholder (inform user and halt), module hook will be inserted later
    """

    has_raw  = bool(selected_files.get("raw"))
    has_mzml = bool(selected_files.get("mzml"))

    # 0) mzML placeholder branch (must come BEFORE pymsreader checks)
    if has_mzml and not has_raw:
        mzml_list = selected_files["mzml"]
        try:
            import pymzml
        except Exception as e:
            messagebox.showerror("pymzML Module Missing", f"Could not import mzmlreader.\n\n{e}")
            return
        status_var.set("Converting mzML file to csv...check the metadata assignment window.")
        progress.start()
        MetadataEditorWindow(
            root,
            mzml_list,
            on_metadata_ready,
            on_finish=reset_main_status,
            input_kind="mzml"
            )
        # 20260519 [B-23-E] removed Preview notification window. mzML is now functional, 
        # the concern of minimal requirements for mzML conversion can be found at user manuals
        return

    # (Optional) prevent ambiguous mixed selection
    if has_mzml and has_raw:
        messagebox.showwarning(
            "Ambiguous Input",
            "Both RAW and mzML files are selected.\n\n"
            "Please convert one format at a time (clear one selection and retry)."
        )
        return

    if not mspext.pymsreader:
        logger.log("[ERROR] Raw file conversion is not available — pymsfilereader missing.")
        messagebox.showerror(
            "Raw Conversion Not Available",
            "Thermo MSFileReader or pymsfilereader is not installed.\n"
            "Please run on a compatible Windows system with the required libraries."
        )
        return

    if "raw" not in selected_files or not selected_files["raw"]:
        messagebox.showwarning("No Raw File selected", "Please select at least one raw file.")
        return

    rawfilelist = selected_files["raw"]
    status_var.set("Converting raw file to csv...check the metadata assignment window.")
    progress.start()
    MetadataEditorWindow(root, rawfilelist, on_metadata_ready, on_finish=reset_main_status)


# Only the essentials are required; Derivatization may be absent in older files.
REQUIRED_META_KEYS = ("Glycan Type", "Mass Analyzer charge mode")

def _is_metadata_dict(d: dict) -> bool:
    return isinstance(d, dict) and all(k in d for k in REQUIRED_META_KEYS)

def _load_json_safely(path):
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None

def _resolve_metadata_for_sample(files: dict, sample_name: str, csv_path: str | None):
    """
    Returns: (meta_path, meta_dict) or (None, None)
    Strategy:
      1) explicit files["metadata"] if it validates
      2) if files["json"] looks like a metadata JSON (.raw.json) → validate & use
      3) if files["json"] is a method JSON → follow samples[sample]["metadata"] relative to method file
      4) else → None, caller decides whether to prompt
    """
    # 1) explicit pointer
    p = files.get("metadata")
    d = _load_json_safely(p)
    if _is_metadata_dict(d):
        return p, d

    # 2) sometimes 'json' already points at the metadata (.raw.json)
    j = files.get("json")
    if isinstance(j, str) and j.lower().endswith(".raw.json"):
        d2 = _load_json_safely(j)
        if _is_metadata_dict(d2):
            return j, d2

    # 3) method file → per-sample 'metadata' entry
    md = _load_json_safely(j)
    if isinstance(md, dict) and "samples" in md:
        entry = md["samples"].get(sample_name)
        if entry is None:
            # tolerant case-insensitive match if the tree cleaned the name
            lowmap = {k.lower(): k for k in md["samples"].keys()}
            k2 = lowmap.get((sample_name or "").lower())
            entry = md["samples"].get(k2) if k2 else None

        if isinstance(entry, dict):
            rel_meta = entry.get("metadata") or entry.get("json")  # tolerate older key
            if rel_meta:
                base = os.path.dirname(j) if j else ""
                cand = rel_meta if os.path.isabs(rel_meta) else os.path.normpath(os.path.join(base, rel_meta))
                d3 = _load_json_safely(cand)
                if _is_metadata_dict(d3):
                    return cand, d3

    # 4) give up (launcher will optionally ask user once)
    return None, None
# Metadata-related support block ends
# ============================

#20260414 code review: it is a fallback (safe net) when marker fails. Around line 5500~? this fallback produces ion_scoring_status = "ok(fallback) in AppLogger"
def _fallback_simple_ion_scoring(matched_df, ion_df, ppm_value):
    """Minimal ion score using validator.findingions + marker.score_counter (anchors-aware)."""
    from mspvalidator_merger import findingions, expandpeaklist

    out = matched_df.copy()
    # ensure peaks are lists (the validator expects python lists, not strings)
    if isinstance(out["peaklist"].iloc[0], str) or isinstance(out["peakintensity"].iloc[0], str):
        out = expandpeaklist(out)

    ionlist_mz = pd.to_numeric(ion_df["mass"], errors="coerce").dropna().to_numpy()

    scores, counts, hits_str = [], [], []
    for _, row in out.iterrows():
        # findingions returns a full listing: [(ion_mz, logI_plus1), ...] length == len(ionlist)
        hit_pairs = findingions(row, ion_df[["mass"]], ppm_value)
        matched_mz = [float(mz) for (mz, logi1) in hit_pairs if float(logi1) > 1.0]
        try:
            # prefer your module’s anchors-aware score if present
            score = float(getattr(marker, "score_counter")(matched_mz, ionlist_mz))
        except Exception:
            # simple fraction fallback
            score = (len(matched_mz) / max(1, len(ionlist_mz))) if len(ionlist_mz) else 0.0
        scores.append(score)
        counts.append(len(matched_mz))
        hits_str.append(";".join(f"{mz:.6f}" for mz in matched_mz))

    out["ion score"] = scores
    out["ion hit count"] = counts
    out["ion hits m/z"] = hits_str
    return out

#20250929 for combining datasets
# --- Combine Trainable Datasets with helper functions useful for other modules ---

# UID related helpers
def _infer_uid_series(df: pd.DataFrame, src_base: str):
    """
    Returns a UID series. Prefer MS2scan_no if present; otherwise index-based.
    UID format: {src_base}#scan{MS2scan_no}  OR  {src_base}#row{n}
    """
    # normalize possible scan column variants
    scan_candidatecols = ["MS2scan_no", "MS2scan", "ScanNo", "scan", "scan_no", "ms2_scan_no"]
    scan_col = next((c for c in scan_candidatecols if c in df.columns), None)

    if scan_col is not None:
        # convert to string for safe concatenation
        return df[scan_col].astype(str).map(lambda s: f"{src_base}#scan{s}")
    else:
        # stable index numbering
        # NOTE: don't use df.index because it may be non-contiguous after prev ops
        return pd.Series([f"{src_base}#row{i}" for i in range(1, len(df) + 1)], index=df.index)

def _ensure_uid_and_origin(df: pd.DataFrame, src_path: str) -> pd.DataFrame:
    """
    - Adds UID if no 'UID'/'uid' column exists.
    - Adds Origin_File and Origin_Basename columns for traceability.
    Returns a new DataFrame (doesn't modify the input df in-place).
    """
    out = df.copy()
    src_abs = os.path.abspath(src_path)
    src_base = os.path.splitext(os.path.basename(src_path))[0]

    # Respect existing UID (case-insensitive)
    uid_col = None
    for c in out.columns:
        if str(c).strip().lower() == "uid":
            uid_col = c
            break
    if uid_col is None:
        out.insert(0, "UID", _infer_uid_series(out, src_base))

    # Always add origin columns (safe to overwrite with same values)
    out.insert(1, "Origin_Basename", src_base)
    out.insert(2, "Origin_File", src_abs)
    return out

def _outer_union_concat(dfs: list[pd.DataFrame]) -> pd.DataFrame:
    """
    Concatenate with outer join to keep superset of columns.
    Missing columns are NA; order: ensure UID/Origin_* stay in front.
    """
    if not dfs:
        return pd.DataFrame()
    big = pd.concat(dfs, axis=0, join="outer", ignore_index=True)

    # Put UID + Origin_* first if present
    front = [c for c in ["UID", "Origin_Basename", "Origin_File"] if c in big.columns]
    rest  = [c for c in big.columns if c not in front]
    return big[front + rest]
# UID related helpers end
# =======================

# Combine Trainable Datasets called in ML analysis window
def combine_trainable_datasets_ui(parent=None):
    """
    UI entry point:
      1) Ask for multiple trainable CSVs
      2) Load each robustly (sniff TSV if needed)
      3) Ensure UID + Origin_*
      4) Concatenate (outer union), prompt for save name
      5) Save combined CSV and a small .txt log
    """
    from tkinter import filedialog, messagebox

    filepaths = filedialog.askopenfilenames(
        title="Select trainable CSV files to combine",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
    )
    if not filepaths:
        return

    dfs = []
    errors = []
    for p in filepaths:
        try:
            # Heuristic: many of your trainables are comma-CSV; but support TSV just in case
            prefer_tab = False
            # if you want to force TSV sniff for certain prefixes:
            # prefer_tab = os.path.basename(p).lower().startswith(("ms2_", "trainable_"))
            df = robust_read_csv(p, prefer_tab=prefer_tab)
            df = _ensure_uid_and_origin(df, p)
            dfs.append(df)
        except Exception as e:
            errors.append(f"{p} -> {e}")

    if not dfs:
        messagebox.showerror("Combine Datasets", "No files could be read.\n\n" + "\n".join(errors))
        return

    combined = _outer_union_concat(dfs)

    # Ask user where to save
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    default_name = f"combined_trainable_{ts}.csv"
    out_csv = filedialog.asksaveasfilename(
        title="Save combined CSV as…",
        defaultextension=".csv",
        initialfile=default_name,
        filetypes=[("CSV files", "*.csv")]
    )
    if not out_csv:
        return

    try:
        combined.to_csv(out_csv, index=False)
    except Exception as e:
        messagebox.showerror("Save Failed", f"Could not save combined CSV:\n{e}")
        return

    # Write a tiny log
    log_path = os.path.splitext(out_csv)[0] + "_combine_log.txt"
    try:
        with open(log_path, "w", encoding="utf-8") as f:
            f.write("GlycoMSP Combine Trainable Datasets Log\n")
            f.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Output CSV: {os.path.abspath(out_csv)}\n\n")
            f.write("Sources:\n")
            for p, df in zip(filepaths, dfs):
                f.write(f"  - {os.path.abspath(p)}  (rows after UID/origin inject: {len(df)})\n")
            f.write("\nRow count (combined): " + str(len(combined)) + "\n")
    except Exception as e:
        # Log failure is non-fatal
        print(f"[combine] Failed to write log: {e}")

    messagebox.showinfo(
        "Combine Datasets",
        "Done!\n\n"
        f"Combined CSV:\n{out_csv}\n\n"
        f"Log:\n{log_path}"
    )
# Combine Trainable Datasets ends
# ==============================

def _utc_now_iso():
    # RFC3339-ish without microseconds
    return datetime.utcnow().replace(microsecond=0).isoformat() + "Z"

# Avoid status symbol contaminate the tree contexts
def clean_sample_name(text):
    for symbol in ["✅", "⚠️", "❌", "⛔"]:
        if text.startswith(symbol):
            text = text[len(symbol):].strip()

    if text.startswith("Sample:"):
        text = text[len("Sample:"):].strip()

    # Remove validation tag suffix
    if " (metadata missing)" in text:
        text = text.replace(" (metadata missing)", "")
    if " (validation failed)" in text:
        text = text.replace(" (validation failed)", "")
    if " (unvalidated)" in text:
        text = text.replace(" (unvalidated)", "")

    return text.strip()
# ====================

# normalized method 


# 20260417 Code review: removed redundant definition and import
# patch in 20260127 to apply new data structure of method json (method v1)
# resolve paths in 3 ways: Empty/None/non-string → returns None; Already absolute returns as-is; Relative path → joins with base and normalizes
def _resolve_path(base, p):
    """Resolve relative paths against method file folder; keep absolute paths (Windows or POSIX) as-is."""
    if not p:
        return None
    if not isinstance(p, str):
        return None
    p = p.strip()
    if not p:
        return None
    # POSIX absolute
    if os.path.isabs(p):
        return p
    # Windows drive absolute, even when running on non-Windows
    if re.match(r"^[A-Za-z]:[\\/]", p):
        return p

    return os.path.normpath(os.path.join(base, p))
# ==

def normalize_method_json(method_obj: dict, method_path: str, base_dir: str = None):
    """
    Normalize ANY accepted method-json family into:
    (exp_name, sample_name, tree_entry_dict, v1_method_dict)

    - tree_entry_dict matches your TreeView keys:
        csv, excel, metadata, json, ionlist_path, insilico_csv, pseudolabel_csv, (optional) trainable_csv
    - v1_method_dict matches the v1 schema you approved (json_type/schema_version + method/sample/inputs/artifacts/parameters)
    """
    if not isinstance(method_obj, dict):
        print(f"[DEBUG][Error] method_obj must be a dict")
        raise ValueError("method_obj must be a dict")

    base = base_dir or os.path.dirname(method_path)
    # ----------------------------
    # A) Detect method family
    # ----------------------------
    jt = method_obj.get("json_type")
    # New v1 schema (preferred)
    is_v1 = (jt == JSON_TYPE_METHOD) and isinstance(method_obj.get("inputs"), dict)
    # Legacy MAS schema (experiment + samples dict)
    is_legacy_mas = (not is_v1) and isinstance(method_obj.get("samples"), dict) and bool(method_obj.get("samples"))
    # Legacy CGA/pseudolabel schema (dataset_type/parents)
    is_legacy_cga = (not is_v1) and (
        method_obj.get("dataset_type") == "pseudolabel" or isinstance(method_obj.get("parents"), dict))

    if not (is_v1 or is_legacy_mas or is_legacy_cga):
        logger.log("[ERROR 1] Unrecognized method JSON structure")
        raise ValueError("Unrecognized method JSON structure (normalize_method_json)")

    # ----------------------------
    # B) Extract exp_name / sample_name + paths into a unified internal record
    # ----------------------------
    if is_v1:
        exp_name = (
            (method_obj.get("sample") or {}).get("experiment_title")
            or (method_obj.get("method") or {}).get("name")
            or "Recovered"
        )
        sample_name = ((method_obj.get("sample") or {}).get("sample_name")
                    or "RecoveredSample")

        inputs = method_obj.get("inputs") or {}
        artifacts = method_obj.get("artifacts") or {}

        csv_path = _resolve_path(base, (inputs.get("converted_csv") or {}).get("path"))
        meta_path = _resolve_path(base, (inputs.get("metadata_json") or {}).get("path"))
        excel_path = _resolve_path(base, (inputs.get("annotation_excel") or {}).get("path"))

        ion_path = _resolve_path(base, (inputs.get("ion_list") or {}).get("path"))
        insilico_path = _resolve_path(base, (inputs.get("insilico_glycan_list") or {}).get("path"))
        scoreb_path = _resolve_path(base, (inputs.get("score_b_workbook") or {}).get("path"))
        #change name from pseudolabels_tsv to CGAresult_tsv
        pl_path = _resolve_path(base, (artifacts.get("CGAresult_tsv") or artifacts.get("pseudolabels_tsv") or {}).get("path"))  
        train_path = _resolve_path(base, (artifacts.get("trainable_csv") or {}).get("path"))

        method_family = (method_obj.get("method") or {}).get("family") or ("CGA" if ion_path or insilico_path else "MAS")

    elif is_legacy_mas:
        exp_name = method_obj.get("experiment") or "Recovered"
        # legacy may contain multiple samples; import them one-by-one upstream
        # here we normalize ONLY the first sample (caller can loop externally if desired)
        sample_name = next(iter(method_obj["samples"].keys()))
        files = method_obj["samples"][sample_name] or {}

        csv_path = _resolve_path(base, files.get("csv"))
        excel_path = _resolve_path(base, files.get("excel"))
        meta_path = _resolve_path(base, files.get("metadata"))

        ion_path = _resolve_path(base, files.get("ionlist_path"))
        insilico_path = _resolve_path(base, files.get("insilico_csv"))
        scoreb_path = _resolve_path(base, files.get("score_b_workbook_path"))
        pl_path = _resolve_path(base, files.get("pseudolabel_csv"))
        train_path = _resolve_path(base, files.get("trainable_csv"))

        method_family = "MAS"  # legacy MAS method files represent MAS by default

    else:  # legacy CGA/pseudolabel
        exp_name = method_obj.get("experiment_title") or "Recovered"
        sample_name = method_obj.get("sample_name") or os.path.splitext(os.path.basename(method_path))[0]
        parents = method_obj.get("parents") or {}
        ionlist = method_obj.get("ionlist") or {}

        csv_path = _resolve_path(base, parents.get("converted_csv"))
        meta_path = _resolve_path(base, parents.get("metadata_json"))
        excel_path = _resolve_path(base, parents.get("annotation_excel"))  # usually absent
        ion_path = _resolve_path(base, ionlist.get("path"))
        insilico_path = _resolve_path(base, parents.get("insilico_csv"))
        scoreb_path = _resolve_path(base, parents.get("score_b_workbook_path")) #files -> parents. Mistake? Found and corrected by Claude Code
        pl_path = _resolve_path(base, parents.get("pseudolabels_tsv") or parents.get("pseudolabel_csv"))
        train_path = _resolve_path(base, parents.get("trainable_csv"))

        method_family = "CGA"

    # normalize sample name the same way GUI expects
    sample_name = clean_sample_name(sample_name)

    # ----------------------------
    # C) Produce TreeView entry (your internal representation)
    # ----------------------------
    tree_entry = {
        "json": method_path,   # method file path belongs here
        "csv": csv_path,
        "excel": excel_path,
        "metadata": meta_path,
        "ionlist_path": ion_path,
        "insilico_csv": insilico_path,
        "score_b_workbook_path": scoreb_path,
        "pseudolabel_csv": pl_path,
        "trainable_csv": train_path,
    }
    # drop empty
    tree_entry = {k: v for k, v in tree_entry.items() if v}

    # ----------------------------
    # D) Produce canonical v1 method dict (in-memory)
    # ----------------------------
    if is_v1:
        v1 = method_obj
        # Ensure headers exist (enforced)
        v1["json_type"] = JSON_TYPE_METHOD
        v1["schema_version"] = v1.get("schema_version") or SCHEMA_V1
        # Ensure timestamps
        v1.setdefault("created_utc", _utc_now_iso())
        v1["updated_utc"] = _utc_now_iso()
        v1.setdefault("uid", str(uuid.uuid4()))
        return exp_name, sample_name, tree_entry, v1

    # Build v1 from legacy shapes
    v1 = {
        "json_type": JSON_TYPE_METHOD,
        "schema_version": SCHEMA_V1,
        "uid": str(uuid.uuid4()),
        "created_utc": _utc_now_iso(),
        "updated_utc": _utc_now_iso(),
        "method": {
            "family": method_family,
            "name": f"{exp_name}:{sample_name}:{method_family}",
            "description": "",
            "tags": []
        },
        "sample": {
            "sample_name": sample_name,
            "experiment_title": exp_name
        },
        "inputs": {
            "converted_csv": {"path": csv_path} if csv_path else None,
            "metadata_json": {"path": meta_path} if meta_path else None,
            "annotation_excel": {"path": excel_path} if excel_path else None,
            "ion_list": {"path": ion_path} if ion_path else None,
            "score_b_workbook": {"path": scoreb_path} if scoreb_path else None,
            "insilico_glycan_list": {"path": insilico_path} if insilico_path else None
        },
        "parameters": {
            "mas": {},
            "cga": {},
            "scoring": {},
            "ml": {}
        },
        "artifacts": {
            "CGAresult_tsv": {"path": pl_path} if pl_path else None, #pseudolabels_tsv -> CGAresult_tsv
            "trainable_csv": {"path": train_path} if train_path else None,
            "unlabeled_csv": None,
            "reports": []
        },
        # 20260521 B-32 (Codex defer + document, NOT a semantic change): the
        # `validation` block here is reserved for future method/file integrity
        # checks (v1.2+), NOT current GUI sample-validation state. Current MAS
        # sample validation is workspace-scoped and persisted in exp.json
        # `samples[].validated` (see export_experiment_json). CGA currently has
        # no equivalent user-facing validation flow; it keeps this scaffold as "unknown".
        "validation": {
            "status": "unknown",
            "checked_utc": None,
            "items": []
        },
        "software": {
            "glycomsp": {"version": "", "commit": ""},
            "extractor": {"name": "", "version": ""}
        },
        "operator": {"name": "", "note": ""}
    }

    # remove nulls in inputs/artifacts for cleanliness
    v1["inputs"] = {k: v for k, v in v1["inputs"].items() if v is not None}
    v1["artifacts"] = {k: v for k, v in v1["artifacts"].items() if v is not None}

    return exp_name, sample_name, tree_entry, v1

# Prepare Dataset Window and related modules
def open_prepare_dataset_window():
    subwin = tk.Toplevel(root)
    subwin.title("Prepare Dataset")
    subwin.geometry("800x550")
    #drag
    drag_data = {
    "item": None,
    "filetype": None,
    "filename": None,
    "from_exp": None,
    "from_sample": None,
}
    
    # --- Nested project data ---
    experiment_projects = {}  # {experiment: {"samples": {sample_name: {csv, excel, json}}}}
    linked_validated_samples = set() #validated samples
    validation_failed_samples = set() #failed sample
    #added 20250412
    experiment_method_paths = {}  # Store .exp.json path per experiment
    sample_method_folder = None  # Global path for saving per-sample method.json files
    experiment_status_labels = {}  # GUI labels for status display, indexed by experiment
    # fixed 20260305 P1-B: uuid import and _utc_now_iso must be at top of outer function, not after inner function definitions that use them


    #20250917 start fixing json issues
    # --- BEGIN: exp.json save/load helpers (generic; PL-ready) ---
    def _abs(p, base_dir):
        if not p: return p
        return os.path.normpath(p if os.path.isabs(p) else os.path.join(base_dir, p))

    # returns experiment title only I think. Role: STATE-MANAGER
    def _current_exp_title():
        """Return selected experiment name; if selection is inside sample/method/file, climb to Experiment node."""
        sel = tree.selection()
        if sel:
            node = sel[0]
            while node:
                txt = tree.item(node, "text") or ""
                if txt.startswith("Experiment:"):
                    return txt.replace("Experiment:", "").strip().split(" (")[0].strip()
                node = tree.parent(node)
        # fallback to first root
        roots = tree.get_children()
        if not roots:
            return None
        txt = tree.item(roots[0], "text") or ""
        return txt.replace("Experiment:", "").strip().split(" (")[0].strip()

    #patch to avoid dropping exp info when load -> save new exp json from legacy format
    #20260415 code review: in future this criteria might change
    # --- Legacy format migration (deprecate after v1.1- publish) ---
    def _infer_method_kinds(files: dict):
        kinds = []
        has_csv = bool(files.get("csv"))
        has_meta = bool(files.get("metadata"))
        if has_csv and has_meta and files.get("excel"):
            kinds.append("MAS")
        if has_csv and has_meta and files.get("ionlist_path") and files.get("insilico_csv"):
            kinds.append("CGA")
        return kinds

    # force method file appears before saving exp json. Produce one method entry per family
    def _ensure_methods_before_saving_exp_v1(exp_title: str):
        """
        Silent upgrade:
        - if sample has legacy links but no method list, generate method v1 files (MAS/CGA) into method folder
        - populate files["_methods"] so saving exp v1 never loses info
        """
        samples = experiment_projects.get(exp_title, {}).get("samples", {})
        if not samples:
            return

        for sample_name, files in samples.items():
            if not isinstance(files, dict):
                continue

            # If methods already present, keep them
            methods = files.get("_methods")
            if isinstance(methods, list) and methods:
                continue

            # If at least one existing method path is linked, promote it into _methods
            if files.get("json"):
                files["_methods"] = [{
                    "method_id": str(uuid.uuid4()),
                    "method_name": os.path.basename(files["json"]),
                    "family": "UNKNOWN",
                    "path": files["json"],
                    "status": "unknown",
                    "notes": "promoted from active method"
                }]
                continue

            # Otherwise: legacy-only links exist → generate method(s)
            kinds = _infer_method_kinds(files)
            if not kinds:
                # nothing we can upgrade
                continue

            files["_methods"] = []

            for kind in kinds:
                try:
                    # Auto-save path in method folder (no dialogs)
                    save_method_v1_for_sample(
                        exp_title,
                        sample_name,
                        out_path=None,
                        kind=kind,
                        auto=True
                    )
                    #fixed 20260312 Test-B B4: prevent duplicate method entries
                    # save_method_v1_for_sample now calls _register_or_update_method_ref internally,
                    # which already appended the correct entry to files["_methods"].
                    # The manual append that was here created a duplicate — removed.
                except Exception as e:
                    logger.log(f"[EXP v1][UPGRADE] Failed to auto-generate {kind} method for {exp_title}/{sample_name}: {e}")

            # Set active method to the first created one (so current UI keeps working)
            if files["_methods"]:
                files["json"] = files["_methods"][0]["path"]


    # --- exp.json save/load helpers ---
    def export_experiment_json(exp_title, out_path):
        """
        Experiment JSON v1:
        - stores only method references + per-method validation status (snapshot + authoritative)
        - does NOT duplicate method contents
        - supports multiple methods per sample via sample["_methods"] (internal cache)
        """
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        samples = experiment_projects.get(exp_title, {}).get("samples", {}) or {}
        _ensure_methods_before_saving_exp_v1(exp_title)
        payload = {
            "experiment": {
                "title": exp_title,
                "description": "",
                "tags": []
            },
            "workspace": {
                "preferred_method_folder": sample_method_folder or "",
                "preferred_exp_path": experiment_method_paths.get(exp_title, "") or "",
                "ui_order": list(samples.keys())
            },
            "samples": [],
            "validation_policy": {
                "hash_alg": "sha256",
                "allow_missing_method": True,
                "allow_hash_mismatch": True
            }
        }
        payload = _ensure_header(payload, JSON_TYPE_EXPERIMENT)

        for sname, files in samples.items():
            files = files or {}

            # Prefer multi-method cache if present; else fall back to active method path (files["json"])
            methods = []
            if isinstance(files.get("_methods"), list) and files["_methods"]:
                methods = files["_methods"]
            else:
                mp = files.get("json")
                if mp:
                    methods = [{
                        "method_name": os.path.basename(mp),
                        "family": _detect_legacy_family(files) or "UNKNOWN",
                        "path": mp,
                        "status": "unknown",
                        "notes": ""
                    }]
            out_methods = []
            for m in methods:
                mpath = m.get("path") or m.get("method_path") or m.get("json")
                if not isinstance(mpath, str) or not mpath.strip():
                    continue

                mpath = mpath.strip()
                exists = os.path.exists(mpath)
                h = file_sha256(mpath) if exists else None

                # if you already store per-method status, keep it; else derive
                status = m.get("status")
                if not status:
                    status = "missing" if not exists else "unknown"

                out_methods.append({
                    "method_id": m.get("method_id") or str(uuid.uuid4()),
                    "method_name": m.get("method_name") or os.path.basename(mpath),
                    "family": m.get("family") or "UNKNOWN",
                    "method_ref": {
                        "path": pathcanon.to_posix_str(mpath),
                        "hash": {"alg": "sha256", "value": h} if h else None,
                        "last_seen_utc": _now_iso_utc()
                    },
                    "validation": {
                        "status": status,
                        "checked_utc": _now_iso_utc(),
                        "notes": m.get("notes", "")
                    }
                })

            # fixed 20260305 P1-A bug: placeholder must be set before appending to payload
            if not out_methods:
                    # Create a non-destructive placeholder so v1 exp.json never drops the sample completely
                out_methods = [{
                    "method_id": str(uuid.uuid4()),
                    "method_name": f"{sname}.INCOMPLETE",
                    "family": "UNKNOWN",
                    "method_ref": {
                        "path": None,
                        "last_seen_utc": _now_iso_utc()
                    },
                    "validation": {
                        "status": "incomplete",
                        "checked_utc": _now_iso_utc(),
                        "notes": "No method JSON available yet (missing required inputs to auto-generate)."
                    }
                }]
            #fixed 20260306 Test-A: persist validation status in exp.json
            payload["samples"].append({
                "sample_name": sname,
                "validated": (exp_title, sname) in linked_validated_samples,
                "methods": out_methods
            })

        # prune None recursively (so hash=None doesn't clutter JSON)
        def _prune_none(x):
            if isinstance(x, dict):
                return {k: _prune_none(v) for k, v in x.items() if v is not None}
            if isinstance(x, list):
                return [_prune_none(v) for v in x]
            return x

        payload = _prune_none(payload)

        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)

        experiment_method_paths[exp_title] = out_path
        if exp_title in experiment_status_labels:  # 20260515 fix B-06: mirror tuple-unpack pattern from load_experiment_method_file (was AttributeError on tuple).
            path_label, _ = experiment_status_labels[exp_title]
            if path_label:
                path_label.set(f"EXP file: {out_path}")
        logger.log(f"[EXP v1] Saved → {out_path}")
        return out_path

    def import_experiment_json(in_path):
        base = os.path.dirname(in_path)

        data = load_typed_json(
            in_path,
            expected_type=JSON_TYPE_EXPERIMENT,
            allow_legacy=True,
            context="[EXP] "
        )

        # v1: {"experiment": {"title": ...}, "samples": [ ... ]}
        # legacy: {"experiment": "...", "samples": { ... }}

        if isinstance(data.get("experiment"), dict):
            exp_title = data["experiment"].get("title") or "Unnamed Experiment"
        else:
            exp_title = data.get("experiment") or "Unnamed Experiment"

        experiment_projects.setdefault(exp_title, {"samples": {}})
        dst = experiment_projects[exp_title]["samples"]

        # ---- v1 loader ----
        if isinstance(data.get("samples"), list):
            for s in data["samples"]:
                if not isinstance(s, dict):
                    continue
                sname = s.get("sample_name")
                if not sname:
                    continue

                # Ensure sample entry exists
                dst.setdefault(sname, {"csv": None, "excel": None, "metadata": None, "json": None})
                files = dst[sname]

                methods = s.get("methods") or []
                norm_methods = []
                active_method_path = None

                for m in methods:
                    if not isinstance(m, dict):
                        continue
                    ref = m.get("method_ref") or {}
                    mpath = ref.get("path") or m.get("path")
                    if not isinstance(mpath, str) or not mpath.strip():
                        continue

                    mpath_abs = _abs(str(pathcanon.to_native_path(mpath)), base)
                    exists = os.path.exists(mpath_abs)

                    # Validate by hash if present
                    status = (m.get("validation") or {}).get("status") or ("missing" if not exists else "unknown")
                    notes = (m.get("validation") or {}).get("notes") or ""
                    expected_hash = ((ref.get("hash") or {}).get("value")) if isinstance(ref.get("hash"), dict) else None

                    if exists and expected_hash:
                        current_hash = file_sha256(mpath_abs)
                        if current_hash and current_hash != expected_hash:
                            status = "warn"
                            notes = (notes + " | hash mismatch").strip(" |")

                    norm_methods.append({
                        "method_id": m.get("method_id"),
                        "method_name": m.get("method_name") or os.path.basename(mpath_abs),
                        "family": m.get("family") or "UNKNOWN",
                        "path": mpath_abs,
                        "status": status,
                        "notes": notes
                    })

                    if active_method_path is None and exists:
                        active_method_path = mpath_abs

                # store list for future UI ("Sample -> Methods")
                files["_methods"] = norm_methods

                #fixed 20260306 Test-A: restore validation status from exp.json into linked_validated_samples set
                if s.get("validated"):
                    linked_validated_samples.add((exp_title, sname))

                # choose an active method (first existing); fallback to first path even if missing
                if active_method_path is None and norm_methods:
                    active_method_path = norm_methods[0]["path"]

                # If we have an active method, load it and populate current sample files
                if active_method_path and os.path.exists(active_method_path):
                    try:
                        #md = load_typed_json(active_method_path, expected_type=JSON_TYPE_METHOD, allow_legacy=True, context="[METHOD] ")
                        #norm = normalize_method_json(md, active_method_path, base_dir=os.path.dirname(active_method_path))
                        # norm is expected to map into your internal keys: csv/excel/metadata/json/ionlist_path/etc.
                        #for k, v in (norm or {}).items():
                        #    files[k] = v
                        #files["json"] = active_method_path  # keep active method linked
                        md = load_typed_json(
                            active_method_path,
                            expected_type=JSON_TYPE_METHOD,
                            allow_legacy=True,
                            context="[METHOD] "
                        )
                        exp2, sample2, tree_entry, v1_obj = normalize_method_json(
                            md,
                            active_method_path,
                            base_dir=os.path.dirname(active_method_path)
                        )

                        # tree_entry is the dict mapping into internal keys: csv/excel/metadata/json/ionlist_path/etc.
                        for k, v in (tree_entry or {}).items():
                            files[k] = v

                        files["json"] = active_method_path  # keep active method linked
                
                    except Exception as e:
                        logger.log(f"[EXP v1] Failed to load method for sample={sname}: {e}")

                #fixed 20260313 LEGACY-1: detect family for v1 samples with no registered methods
                # Safety net: v1-format file has empty methods array but has CGA/MAS artifacts.
                if not files.get("_methods"):
                    _fam_auto = _detect_legacy_family(files)
                    if _fam_auto:
                        logger.log(f"[LEGACY] {sname}: no methods in v1 JSON → auto-detecting {_fam_auto}")
                        _ensure_method_stub(files, family=_fam_auto, sample_name=sname)
                        try:
                            _leg_path = save_method_v1_for_sample(exp_name=exp_title, sample_name=sname, kind=_fam_auto, auto=True)
                            logger.log(f"[LEGACY] Created: {_leg_path}")
                        except Exception as e_leg:
                            logger.log(f"[LEGACY] Conversion failed for {sname}: {e_leg}")

                #fixed 20260313 ML-1: auto-create MAS method for legacy v11 samples
                # Validated samples loaded from v11 JSON may lack a MAS method entry.
                # Pre-register stub and attempt auto-save so the tree shows a MAS node.
                if s.get("validated") and files.get("excel"):
                    _ensure_method_stub(files, family="MAS", sample_name=sname)
                    try:
                        save_method_v1_for_sample(exp_name=exp_title, sample_name=sname, kind="MAS", auto=True)
                    except Exception as e_ml1:
                        logger.log(f"[ML-1] auto-create MAS failed for {sname}: {e_ml1}")

            experiment_method_paths[exp_title] = in_path
            logger.log(f"[EXP v1] Loaded ← {in_path}")
            refresh_tree()
            return exp_title

        # ---- legacy loader (keep your old behavior) ----
        for sname, sample_blob in (data.get("samples") or {}).items():
            if isinstance(sample_blob, dict) and "files" in sample_blob:
                files = sample_blob.get("files") or {}
            else:
                files = sample_blob or {}
            resolved = {k: _abs(str(pathcanon.to_native_path(v)), base) for k, v in files.items() if isinstance(v, str) and v.strip()}
            dst[sname] = resolved

            #fixed 20260313 LEGACY-1: auto-detect family and create method for old-format samples
            # Old deployed format has no _methods array; detect from artifact key patterns.
            _fam_auto = _detect_legacy_family(dst[sname])
            if _fam_auto:
                logger.log(f"[LEGACY] {sname}: old format detected → {_fam_auto}")
                print(f"[LEGACY DEBUG] {sname}: family={_fam_auto}, "
                      f"insilico_csv={bool(dst[sname].get('insilico_csv'))}, "
                      f"ionlist={bool(dst[sname].get('ionlist_path'))}, "
                      f"excel={bool(dst[sname].get('excel'))}")
                _ensure_method_stub(dst[sname], family=_fam_auto, sample_name=sname)
                try:
                    _leg_path = save_method_v1_for_sample(exp_name=exp_title, sample_name=sname, kind=_fam_auto, auto=True)
                    logger.log(f"[LEGACY] Created: {_leg_path}")
                except Exception as e_leg:
                    logger.log(f"[LEGACY] Conversion failed for {sname}: {e_leg}")

        experiment_method_paths[exp_title] = in_path
        logger.log(f"[EXP legacy] Loaded ← {in_path}")
        refresh_tree()
        return exp_title
    # --- END: exp.json save/load helpers ---
    # =====================================

    #20250917 future placeholder: integrity check
    # Reuse an existing toolbar if present; otherwise create one
    try:
        toolbar
    except NameError:
        toolbar = ttk.Frame(subwin)
        toolbar.pack(fill="x", padx=10, pady=(8, 0))

    def _check_integrity_clicked():
        # For now: just a friendly placeholder dialog
        tk.messagebox.showinfo(
            "Project integrity",
            "This feature is in development.\n"
            "Planned: compute & store file hashes, check missing&changed files,\n"
            "and help you to build portable or rebuild migrated project files."
        )
    ttk.Button(toolbar, text="Check integrity",state="disabled", command=_check_integrity_clicked).pack(side="left", padx=4)


    #Main sample node resolver, important
    def _get_selected_context():
        """
        Returns:
        exp_name, sample_name, method_path (may be None if no method)
        Works if user clicks Sample/Method/File nodes.
        """
        sel = tree.selection()
        if not sel:
            return None, None, None

        node = sel[0]

        # climb to sample
        cur = node
        sample_node = None
        while cur:
            t = tree.item(cur, "text")
            if "Sample:" in t:
                sample_node = cur
                break
            cur = tree.parent(cur)
        if not sample_node:
            return None, None, None

        exp_node = tree.parent(sample_node)
        exp_text = tree.item(exp_node, "text")
        exp_name = exp_text.replace("Experiment:", "").strip().split(" (")[0].strip()

        sample_text = tree.item(sample_node, "text")
        sample_name = clean_sample_name(sample_text)

        # climb to method node (optional)
        cur = node
        method_node = None
        while cur:
            t = tree.item(cur, "text")
            if t.startswith("Method:"):
                method_node = cur
                break
            cur = tree.parent(cur)

        files = experiment_projects.get(exp_name, {}).get("samples", {}).get(sample_name, {})
        method_path = None

        if method_node:
            # find "Method JSON:" child basename
            basename = None
            for child in tree.get_children(method_node):
                ct = tree.item(child, "text")
                if ct.startswith("Method JSON:"):
                    basename = ct.replace("Method JSON:", "").strip()
                    break
            if basename:
                # match by basename in _methods
                for m in (files.get("_methods") or []):
                    p = m.get("path")
                    if p and os.path.basename(p) == basename:
                        method_path = p
                        break
                # fallback: active
                if method_path is None and files.get("json") and os.path.basename(files["json"]) == basename:
                    method_path = files["json"]

        if method_path is None:
            method_path = files.get("json")  # active method fallback

        return exp_name, sample_name, method_path
    # Main sample node resolver ends
    # =============================

    # ensure method node exists, if no, add an unknown for it
    def _ensure_method_stub(files: dict, family: str, sample_name: str):
        """Ensure there is a Method node to hang CGA/MAS artifacts under in the tree."""
        if not isinstance(files, dict):
            return

        methods = files.get("_methods")
        if not isinstance(methods, list):
            methods = []
            files["_methods"] = methods

        famU = (family or "").upper()

        for m in methods:
            if (m.get("family") or "").upper() == famU:
                return

        methods.append({
            "method_id": str(uuid.uuid4()),
            "method_name": f"{sample_name}.{famU}.method.v1 (unsaved)",
            "family": famU,
            "path": None,
            "status": "unknown",
            "notes": ""
        })

    # --- Legacy format migration (deprecate after v1.1 publish) ---
    def _detect_legacy_family(files: dict):
        """Auto-detect method family from artifact patterns in a legacy exp.json sample.
        Returns 'CGA', 'MAS', or None.
        CGA: has insilico_csv or ionlist_path (constraint-based annotation artifacts)
        MAS: has excel (manual annotation spreadsheet)
        """
        if not isinstance(files, dict):
            return None
        if files.get("insilico_csv") or files.get("ionlist_path"):
            return "CGA"  # Has CGA-specific artifacts
        if files.get("excel"):
            return "MAS"  # Has manual annotation spreadsheet
        return None

    # main updater for method (called when method file needs to be updated)
    def _register_or_update_method_ref(files: dict, out_path: str, family: str):
        """Attach a saved method JSON path to the right method entry and keep legacy 'json' in sync."""
        famU = (family or "").upper()
        methods = files.get("_methods")
        if not isinstance(methods, list):
            methods = []
            files["_methods"] = methods

        # update existing
        for m in methods:
            if (m.get("family") or "").upper() == famU:
                m["path"] = out_path
                m["method_name"] = os.path.basename(out_path)
                m["status"] = m.get("status") or "unknown"
                break
        else:
            methods.append({
                "method_id": str(uuid.uuid4()),
                "method_name": os.path.basename(out_path),
                "family": famU,
                "path": out_path,
                "status": "unknown",
                "notes": ""
            })

        # legacy compatibility: lots of code still reads files["json"]
        files["json"] = out_path

    # 20260521 fix B-31 (MAS+CGA, Codex Option B + Q7): resolve method.json output path.
    # Reuse existing registered path for this (sample, family) if it still exists on disk;
    # else return a stable name (no timestamp). Replaces the prior timestamped-on-every-save
    # convention that accumulated orphan method.json files in sample_method_folder.
    # Hash mismatch on overwrite is acceptable per validation_policy.allow_hash_mismatch=True
    # (warn-only UX until exp.json is re-saved).
    def _resolve_method_out_path(files: dict, sample_name: str, family: str, base_dir: str) -> str:
        famU = (family or "").upper()
        methods = files.get("_methods") or []
        for m in methods:
            if (m.get("family") or "").upper() == famU:
                mp = m.get("path")
                if isinstance(mp, str) and mp.strip() and os.path.exists(mp):
                    return mp
        return os.path.join(base_dir, f"{sample_name}.{famU}.method.json")

    # Method export writer
    def _export_method_v1(force_family: str):
        exp_name, sample_name, method_path = _get_selected_context() 
        if not exp_name or not sample_name:
            return
        # Build v1 dict from current links (no writing yet)
        try:
            v1 = build_method_v1_from_tree(exp_name, sample_name, force_family=force_family)
        except Exception as e:
            messagebox.showerror("Method v1", str(e))
            return
        # Ask user where to save
        csv_path = experiment_projects[exp_name]["samples"][sample_name].get("csv")
        default_dir = os.path.dirname(csv_path) if csv_path else os.getcwd()
        default_name = f"{sample_name}.{force_family}.method.json"  #method.v1.json -> method.json
        out_path = filedialog.asksaveasfilename(
            title="Save Method v1 JSON",
            initialdir=default_dir,
            initialfile=default_name,
            defaultextension=".json",
            filetypes=[("JSON files", "*.json")]
        )
        if not out_path:
            return
        # Save
        try:
            # 20260521 fix B-31 sister (logger gap close): manual-export observability.
            # filedialog already prompts "overwrite?" if file exists, but a forensic
            # log entry is still useful since users may intentionally pick the same name.
            _path_existed = os.path.exists(out_path)
            with open(out_path, "w", encoding="utf-8") as f:
                json.dump(v1, f, indent=2, ensure_ascii=False)
            files = experiment_projects[exp_name]["samples"][sample_name]
            _register_or_update_method_ref(files, out_path, family=force_family)
            refresh_tree()
            logger.log(f"[Method v1][{force_family}][Manual Export] {'Updated existing' if _path_existed else 'Created new'}: {out_path}")
        except Exception as e:
            messagebox.showerror("Method v1", f"Failed to save:\n{out_path}\n{e}")
            return
        messagebox.showinfo("Method v1", f"Saved:\n{out_path}")
    # Method export writer ends
    # =========================

    # Some definitions/parameters will be called in Widgets within Prepare Dataset panel
    include_mass_feat_var = tk.BooleanVar(value=False)
    #20250905 added for negative label
    # --- Negative sampling (Prepare Dataset) ---
    add_negatives_var   = tk.BooleanVar(value=False)
    neg_ratio_var       = tk.DoubleVar(value=3.0)   # max neg : pos
    #neg_scorecol_var    = tk.StringVar(value="score")
    #neg_scorethr_var    = tk.DoubleVar(value=0.05)  # keep scans with max(score) < thr
    neg_markermin_var   = tk.IntVar(value=1)        # require < this many marker hits
    min_hits_var      = tk.IntVar(value=3)      # keep scans with < min_hits hits
    ppm_tol_var       = tk.StringVar(value="10")# ppm tolerance; text with validation
    n_features_var    = tk.IntVar(value=0)      # total ion features (from ion_df)
    gate_hint_var     = tk.StringVar(value="Ion features not loaded yet")
    # --- Ion mining (ion suggestion in source code) options ---
    ion_suggest_enable_var   = tk.BooleanVar(value=False)
    ion_suggest_ppm_var      = tk.DoubleVar(value=10.0)
    ion_suggest_dafloor_var  = tk.DoubleVar(value=0.03)
    ion_suggest_minsupp_var  = tk.IntVar(value=5)     # min glycan support per ion
    ion_suggest_topk_var     = tk.IntVar(value=60)    # how many to export
    last_suggest_csv_var = tk.StringVar(value="")
    # ===========================

    # --- Treeview UI: panel of sample nodes showing ---
    tree = ttk.Treeview(subwin)
    tree.heading("#0", text="Dataset Explorer", anchor="w")
    tree.pack(expand=True, fill="both", padx=10, pady=10)

    # --- Tree logic block---
    # Check validation when select a tree node. Updates the Link and Merge button states based on what the user clicked in the tree.
    def on_tree_select(event):
        exp_name, sample_name, method_path = _get_selected_context()
        if not exp_name or not sample_name:
            link_button.config(state="disabled")
            merge_button.config(state="disabled")
            return

        files = experiment_projects.get(exp_name, {}).get("samples", {}).get(sample_name, {}) or {}
        has_merge_inputs = bool(files.get("csv") and files.get("excel") and files.get("metadata"))
        is_validated = (exp_name, sample_name) in linked_validated_samples

        if not is_validated:
            link_button.config(state="normal")
            merge_button.config(state="disabled")
        else:
            link_button.config(state="disabled")
            merge_button.config(state="normal" if has_merge_inputs else "disabled")
        # fixed 20260306 Test-A: removed circular binding from handler (was inside else block, never triggered at init)

    # Tree-view selection event listener(handler)
    tree.bind("<<TreeviewSelect>>", on_tree_select)  # fixed 20260306 Test-A: added TreeviewSelect binding to setup

    # Link & Validation on MAS method sample. Considering CGA compatibility
    def try_link_selected_sample():
        sel = tree.selection()
        if not sel:
            return
        node = sel[0] # only select the first one selected (single selection)
        # climb until we find a Sample node (or root)
        while node:
            t = tree.item(node, "text")
            if "Sample:" in t:   # works with ✅/⚠️/❌ prefixes since you already clean those
                break
            node = tree.parent(node)

        if not node:
            return  # not inside a sample

        sample_id = node
        exp_id = tree.parent(sample_id)
        sample_name = clean_sample_name(tree.item(sample_id, "text"))
        exp_name = tree.item(exp_id, "text").replace("Experiment: ", "").split(" (")[0].strip()
        link_and_validate_sample(exp_name, sample_name)

    # --- Tree refresh logic ---
    def refresh_tree():
        # 20260415 code review: attempt to fix elements drifting issue
        # remember current selection context                  
        old_exp, old_sample, _ = _get_selected_context()
        tree.delete(*tree.get_children())

        for exp_title, exp_data in experiment_projects.items():
            exp_node = tree.insert("", "end", text=f"Experiment: {exp_title}", open=True)

            for sample_name, files in exp_data.get("samples", {}).items():
                has_csv = bool(files.get("csv"))
                has_excel = bool(files.get("excel"))
                #has_json = bool(files.get("json")) -> now it's clearly metadata
                has_meta = bool(files.get("metadata"))

                # Decide how to display the sample label
                if (exp_title, sample_name) in linked_validated_samples:
                    sample_display = f"✅ Sample: {sample_name}"
                elif (exp_title, sample_name) in validation_failed_samples:
                    sample_display = f"⛔ Sample: {sample_name} (validation failed)"
                #elif has_csv and has_excel and has_json:
                elif has_csv and has_excel and has_meta:
                    sample_display = f"⚠️ Sample: {sample_name} (unvalidated)"
                #elif not has_json:
                elif not has_meta:
                    sample_display = f"❌ Sample: {sample_name} (metadata missing)"
                else:
                    sample_display = f"Sample: {sample_name}"

                sample_node = tree.insert(exp_node, "end", text=sample_display, open=True)
                #20260202 update 1
                # ---- (A) Show sample-level anchors under Sample ----
                # These are "sample identity" / shared inputs
                for ftype in ["csv", "excel", "metadata"]:
                    if files.get(ftype):
                        if ftype == "metadata":
                            label = "Metadata"
                        else:
                            label = ftype.upper()
                        tree.insert(sample_node, "end", text=f"{label}: {os.path.basename(files[ftype])}")

                # ---- (B) Build method list ----
                # New: allow multiple methods per sample.
                # Internal cache (future): files["_methods"] = [ {method_name, family, path, ...}, ... ]
                methods = []
                if isinstance(files.get("_methods"), list) and files["_methods"]:
                    methods = files["_methods"]
                else:
                    # fallback: current flat structure has only one active method path in files["json"]
                    if files.get("json"):
                        methods = [{
                            "method_name": os.path.basename(files["json"]),
                            "family": "UNKNOWN",
                            "path": files["json"],
                        }]

                # ---- (C) Insert Method nodes ----
                if not methods:
                    # no method yet → show a placeholder so users see what's missing
                    tree.insert(sample_node, "end", text="Method: (none)")
                else:
                    for idx, m in enumerate(methods, start=1):
                        mname = m.get("method_name") or f"Method {idx}"
                        fam = m.get("family") or "UNKNOWN"
                        method_node = tree.insert(sample_node, "end", text=f"Method: {mname} [{fam}]", open=False)

                        # method JSON itself
                        mpath = m.get("path")
                        if mpath:
                            tree.insert(method_node, "end", text=f"Method JSON: {os.path.basename(mpath)}")

                        # fixed 20260305 P2-B: TODO Phase 1.5 — per-method artifact display
                        # Currently shows active artifacts under first method node only.
                        # Full fix requires per-method artifact storage in files["_methods"] entries (Phase 1.5 Step 2).
                        # Each method entry (m) already carries m.get("path"), but the artifact keys
                        # ("insilico_csv", "ionlist_path", "pseudolabel_csv") are stored flat in the
                        # outer "files" dict and are not yet separated per-method.
                        # Keeping idx == 1 guard until per-method artifact storage is implemented.
                        if idx == 1:
                            for ftype in ["insilico_csv", "ionlist_path", "score_b_workbook_path", "pseudolabel_csv", "trainable_csv"]:
                                if files.get(ftype):
                                    if ftype == "insilico_csv":
                                        label = "In-silico CSV"
                                    elif ftype == "ionlist_path":
                                        label = "Ion List"
                                    elif ftype == "pseudolabel_csv":
                                        label = "CGA result tsv"
                                    elif ftype == "score_b_workbook_path":
                                        label = "CGA score B reference"
                                    elif ftype == "trainable_csv":
                                        label = "CGA-based trainable csv"
                                    else:
                                        label = ftype.upper()
                                    tree.insert(method_node, "end", text=f"{label}: {os.path.basename(files[ftype])}")

                        #fixed 20260316: MAS artifact display — load inputs from method JSON for per-method parity
                        # Reads converted_csv / annotation_excel / metadata_json from the saved method file.
                        # Unlike CGA (which reads from the flat files dict due to Phase 1.5 limitation),
                        # MAS reads directly from the method JSON so multi-method display is correct:
                        # each saved method has its own inputs dict, no idx==1 guard needed.
                        if fam == "MAS" and mpath and os.path.exists(mpath):
                            try:
                                with open(mpath, "r", encoding="utf-8") as _mf:
                                    _md = json.load(_mf)
                                _inputs = _md.get("inputs", {})
                                for _key, _label in [
                                    ("converted_csv",    "Converted CSV"),
                                    ("annotation_excel", "Annotation Excel"),
                                    ("metadata_json",    "Metadata"),
                                ]:
                                    _entry = _inputs.get(_key)
                                    if isinstance(_entry, dict) and _entry.get("path"):
                                        tree.insert(method_node, "end",
                                                     text=f"{_label}: {os.path.basename(_entry['path'])}")
                            except Exception:
                                pass  # silently skip if method JSON is missing or malformed
        # 20260415 restore selection after rebuild
        if old_exp and old_sample:                                                   
            for exp_iid in tree.get_children():           
                if tree.item(exp_iid, "text").replace("Experiment:","").strip().split(" (")[0].strip() == old_exp:
                    for sample_iid in tree.get_children(exp_iid):                    
                        if clean_sample_name(tree.item(sample_iid, "text")) == old_sample:
                            tree.selection_set(sample_iid)                           
                            tree.see(sample_iid)          
                            return
    # Tree refresh logic ends
    # =======================


    # --- Add sample (customized empty node) to selected experiment ---
    def add_sample():
        sel = tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select an experiment to add a sample.")
            return

        item_id = sel[0]
        item_text = tree.item(item_id, "text")
        if not item_text.startswith("Experiment:"):
            messagebox.showinfo("Invalid Selection", "You must select an experiment, not a sample or file.")
            return

        exp_title = item_text.replace("Experiment: ", "")
        sample_name = simpledialog.askstring("New Sample", "Enter sample name:")
        if not sample_name:
            return
        # 20260415 code review: Plan to support multiple method on same sample, the logic might need revision
        if sample_name in experiment_projects[exp_title]["samples"]:
            messagebox.showwarning("Duplicate Sample", f"Sample '{sample_name}' already exists.")
            return
        experiment_projects[exp_title]["samples"][sample_name] = {"csv": None, "excel": None, "metadata": None, "json": None} #20260415 added "metadata": None
        refresh_tree()

    # Temp stable solution to remove sample or experiment with no entry inside
    def clean_unassigned_samples():
        exp = "Unassigned"
        if exp not in experiment_projects:
            return
        samples = experiment_projects[exp]["samples"]
        to_delete = []

        for s_name, files in samples.items():
            if not files.get("csv") and not files.get("excel") and not files.get("json"):
                to_delete.append(s_name)

        for s in to_delete:
            del samples[s]

        if not samples:
            del experiment_projects[exp]  # if no samples left, remove experiment too
        refresh_tree()

    # Opens MetadataEditorWindow to create post-metadata when user clicks Link Sample button on a sample with both CSV and EXCEL (MAS workflow)
    def store_metadata_back_to_sample(exp_name, sample_name, rawfile, metadata, savename, outdir=None):
        base_dir = outdir or os.path.dirname(rawfile)
        meta_path = os.path.join(base_dir, savename + ".json")

        # Fallback if sample doesn't exist yet (e.g., metadata was created before file assignment)
        # Ensure experiment and sample containers exist
        if exp_name not in experiment_projects:
            experiment_projects[exp_name] = {"samples": {}, "generated_on": datetime.now().strftime("%Y-%m-%d %H:%M")}
        if sample_name not in experiment_projects[exp_name]["samples"]:
            experiment_projects[exp_name]["samples"][sample_name] = {"csv": None, "excel": None, "metadata": None, "json": None} #20260415 added "metadata": None

        # Now safe to reference sample
        sample = experiment_projects[exp_name]["samples"][sample_name]
        sample.setdefault("metadata", None)

        # Enforce header fields on newly generated metadata
        metadata["json_type"] = "glycomsp.metadata"
        metadata["schema_version"] = "1.0.0"
        # Store as metadata (not method)
        sample["metadata"] = meta_path

        refresh_tree() # 20260415 code review: is refresh needed indeed?

        # Try to rename sample to raw name (from metadata), only if different
        try:
            raw_base = os.path.splitext(os.path.basename(metadata.get("Raw filename", "")))[0]
        except Exception:
            raw_base = sample_name

        if raw_base != sample_name and raw_base not in experiment_projects[exp_name]["samples"]:
            experiment_projects[exp_name]["samples"][raw_base] = sample
            del experiment_projects[exp_name]["samples"][sample_name]
            sample_name = raw_base

        # 20260415 code review: automatically callback and link sample once metadata is generated
        if sample.get("csv") and sample.get("excel") and sample.get("metadata"):   
            logger.log(f"[Validation][auto]: start link and validation using created {sample['metadata']}\n")      
            link_and_validate_sample(exp_name, sample_name)

        refresh_tree()                                     
    
    # Real window caller that calls store_metadata_back_to_sample
    def open_metadata_editor_for_sample(exp_name, sample_name):
        from tkinter import simpledialog

        sample = experiment_projects[exp_name]["samples"][sample_name]

        # Let user select raw file if needed (optional)
        raw_file_path = filedialog.askopenfilename(
            title="Select RAW file for this sample (optional)",
            filetypes=[("Raw files", "*.raw")]
        )
        raw_file_path = raw_file_path or "(not linked)"
        print(f"[debug] set to not linked and need to check if we can add one back in metadata editor window")

        # Ask user for output folder
        outdir = filedialog.askdirectory(title="Select output folder to save metadata")
        if not outdir:
            messagebox.showwarning("Cancelled", "Metadata creation cancelled.")
            return

        editor = MetadataEditorWindow(
            parent=root,
            raw_file_list=[raw_file_path],
            on_each_metadata_ready_callback=lambda rf, md, name: store_metadata_back_to_sample(exp_name, sample_name, rf, md, name, outdir),#20260126
            output_dir=outdir,
            skip_conversion=True
        )
        # 20260415 code review: previous patch in __init__ of MetadataEditorWindow already declares this, no need to prevent too early load file before rendering in tk widget
        # Delay execution of load_file until window is fully initialized
        #if not editor.skip_conversion:
        #    editor.window.after(10, lambda: editor.load_file(raw_file_path))
    #  Post-metadata creating widget block ends
    #  =========================================

    # Negative Sampling Options block
    # using _get_selected_context to return sample dict for downstream work (later blocks) referencing certain entry within it. 
    def current_selected_files():
        exp_name, sample_name, method_path = _get_selected_context()
        if not exp_name or not sample_name:
            return None
        try:
            return experiment_projects[exp_name]["samples"][sample_name]
        except KeyError:
            return None

    # Negative Sampling Excel picker
    def _fetch_ion_count_for_dialog():
        # Try: current sample’s Excel; else let user pick one.
        from tkinter import filedialog, messagebox
        files = current_selected_files()
        excel_path = None
        try:
            excel_path = files.get("excel")
        except Exception:
            excel_path = None

        if not excel_path:
            excel_path = filedialog.askopenfilename(
                title="Select manual annotation Excel (ionlist)",
                filetypes=[("Excel", "*.xlsx *.xls")]
            )
            if not excel_path:
                return  # user cancelled

        try:
            n, _ = mspval.peek_ion_feature_count(excel_path)
            n_features_var.set(int(n))
            _update_gate_hint()
        except Exception as e:
            messagebox.showwarning("Ion list", f"Could not read ionlist:\n{e}")

    # Live update the tk var value in Negative Sampling Options when the ion(feature) counts are known
    def _update_gate_hint():
        try:
            k = int(min_hits_var.get())
        except Exception:
            k = 3
        n = int(n_features_var.get())
        if n > 0:
            pct = 100.0 * k / n
            gate_hint_var.set(f"Gating non-glycan: ion matches < {k} hits (~{pct:.1f}% ; {k} of {n})")
        else:
            gate_hint_var.set(f"Gating non-glycan: ion matches < {k} hits (feature count unknown yet)")

    # Prepare Negative Sampling Options Panel that passes/update the negative sampling accross MAS (merge MAS) and CGA (CGA to trainable) 
    def open_negative_options_dialog():
        dlg = tk.Toplevel(root)        # you already switched to root
        dlg.title("Negative Sampling Options")
        dlg.resizable(False, False)
        dlg.grab_set()

        # Try to auto-load ion feature count from the selected sample
        try:
            files = current_selected_files()
            excel_or_ion_path = None
            if files:
                # Prefer annotated Excel; if your pipeline sometimes stores an external ion file, check those too
                excel_or_ion_path = (
                    files.get("excel") or
                    files.get("ion_sheet_file") or
                    files.get("ionlist_path")
                )
            if excel_or_ion_path:
                n, _ = mspval.peek_ion_feature_count(excel_or_ion_path)
                n_features_var.set(int(n))
        except Exception as e:
            print("[neg opts] auto ionlist load failed:", e)

        # Update the hint line (uses n_features_var + min_hits_var)
        _update_gate_hint()

        tk.Checkbutton(dlg, text="Add Non-glycan entries (easy negatives)",
                    variable=add_negatives_var).grid(row=0, column=0, columnspan=2,
                                                        sticky="w", padx=10, pady=(10,6))

        # ratio row
        tk.Label(dlg, text="Max ratio (neg:pos):").grid(row=1, column=0, sticky="e", padx=10)
        tk.Spinbox(dlg, from_=0.0, to=10.0, increment=0.5, width=6,
                textvariable=neg_ratio_var).grid(row=1, column=1, sticky="w", padx=6, pady=2)

        # ion-hit gate
        tk.Label(dlg, text="Keep scans with < min hits to ion list:").grid(row=2, column=0, sticky="e", padx=10)
        tk.Spinbox(dlg, from_=0, to=50, increment=1, width=6,
                textvariable=min_hits_var).grid(row=2, column=1, sticky="w", padx=6, pady=2)

        # ppm tolerance with float validation
        def _valid_float(s: str) -> bool:
            # allow empty (user mid-typing), or digits/one dot (no lone '.')
            return (s == "") or (re.fullmatch(r"\d+(\.\d+)?", s) is not None)
        vcmd = dlg.register(lambda P: _valid_float(P))

        tk.Label(dlg, text="Ion match tolerance (ppm):").grid(row=3, column=0, sticky="e", padx=10)
        tk.Entry(dlg, textvariable=ppm_tol_var, width=8,
                validate="key", validatecommand=(vcmd, "%P")
                ).grid(row=3, column=1, sticky="w", padx=6, pady=2)

        tk.Button(dlg, text="Load ion list now…",
          command=_fetch_ion_count_for_dialog).grid(row=4, column=0, columnspan=2,
                                                            sticky="w", padx=10, pady=(6,0))

        # react when min-hits changes
        def _on_hits_change(*_): _update_gate_hint()
        min_hits_var.trace_add("write", _on_hits_change)

        # initialize text
        _update_gate_hint()

        # live hint label (single placement at clean row; see freeze v1.10 S6.5 cluster_triage_2026_05_18.md issue m)
        tk.Label(dlg, textvariable=gate_hint_var, fg="#555").grid(row=6, column=0, columnspan=2,
                                                                sticky="w", padx=10, pady=(6,10))
        tk.Button(dlg, text="Close", command=dlg.destroy).grid(row=7, column=0, columnspan=2, pady=(4,10))
    # Negative Sampling Options block end
    # ===================================

    # Ion mining block (suggest features - fragment ions)
    def open_ion_mining_dialog():
        dlg = tk.Toplevel(root)
        dlg.title("Ion Mining")
        dlg.resizable(False, False)
        dlg.grab_set()

        tk.Checkbutton(dlg, text="Generate ion suggestions on merge",
                    variable=ion_suggest_enable_var).grid(row=0, column=0, columnspan=2,
                                                            sticky="w", padx=10, pady=(10,6))

        tk.Label(dlg, text="ppm tolerance:").grid(row=1, column=0, sticky="e", padx=10)
        tk.Spinbox(dlg, from_=1.0, to=50.0, increment=0.5, width=6,
                textvariable=ion_suggest_ppm_var).grid(row=1, column=1, sticky="w", padx=6, pady=2)

        tk.Label(dlg, text="DA floor (low-m/z merge):").grid(row=2, column=0, sticky="e", padx=10)
        tk.Spinbox(dlg, from_=0.00, to=0.10, increment=0.005, width=6,
                format="%.3f", textvariable=ion_suggest_dafloor_var).grid(row=2, column=1, sticky="w", padx=6, pady=2)

        tk.Label(dlg, text="Min glycan support per ion:").grid(row=3, column=0, sticky="e", padx=10)
        tk.Spinbox(dlg, from_=1, to=100, increment=1, width=6,
                textvariable=ion_suggest_minsupp_var).grid(row=3, column=1, sticky="w", padx=6, pady=2)

        tk.Label(dlg, text="Top-K suggestions:").grid(row=4, column=0, sticky="e", padx=10)
        tk.Spinbox(dlg, from_=5, to=500, increment=5, width=6,
                textvariable=ion_suggest_topk_var).grid(row=4, column=1, sticky="w", padx=6, pady=(2,10))

        tk.Button(dlg, text="Close", command=dlg.destroy).grid(row=5, column=0, columnspan=2, pady=(2,10))
    # Ion mining block Ends
    # =====================

    # Validate sample status and link the sample (update status in tree, should be captured and saved to method json)
    def link_and_validate_sample(exp_name, sample_name):

        sample_name = clean_sample_name(sample_name)
        if sample_name not in experiment_projects[exp_name]["samples"]:
            logger.log(f"[TREEVIEW Error] Clean sample name '{sample_name}' not found under '{exp_name}'")
            messagebox.showerror("Invalid Sample", f"Sample not found in experiment: {sample_name}")
            return
        sample = experiment_projects[exp_name]["samples"][sample_name]

        # Check file presence
        if not sample.get("csv") or not sample.get("excel"):
            logger.log(f"[Error 0] Missing either converted spectral data or annotation sheet")
            messagebox.showerror("Missing File", "Sample must have both a CSV and Excel file before linking.")
            return
        
        # Validate CSV format
        if not mspval.validate_csv_structure(sample["csv"]): #when it returns False meaning failed
            messagebox.showerror("CSV Validation Failed", f"File: {os.path.basename(sample['csv'])}\nSee log for details.")
            logger.log(f"[Validation Failed] Spectral data CSV failed: {sample['csv']}\n")
            validation_failed_samples.add((exp_name, sample_name))
            refresh_tree()
            return  # ← this prevents continuing to metadata
        else:
            logger.log(f"[Validation] Spectral data CSV successful: {sample['csv']}\n")
        #Validate Excel format (assume sheet 'MSlist' exists and ion list is valid) 

        try:
            xl = pd.read_excel(sample["excel"], sheet_name=None)
            if not isinstance(xl, dict):
                raise ValueError("Excel read error: file may be corrupted or unreadable.")
            if "MSlist" not in xl:
                raise ValueError("Sheet 'MSlist' not found in Excel.")
        except Exception as e:
            messagebox.showerror("Excel Validation Failed", f"{sample['excel']}\n\n{e}")
            logger.log(f"[Validation Failed] Annotation sheet failed: {sample['excel']}\n{e}")
            validation_failed_samples.add((exp_name, sample_name))
            refresh_tree()
            return  # ← this prevents continuing to metadata

        #try:
        anno = pd.ExcelFile(sample["excel"])
        MSlistdf = pd.read_excel(anno, sheet_name="MSlist")
        ionlistdf = pd.read_excel(anno, sheet_name="ionlist")
        if not mspval.validate_annotation_structure(MSlistdf):
            messagebox.showerror("Excel Validation Failed", f"File: {os.path.basename(sample['excel'])}\nSee log for details.")
            logger.log(f"[Validation Failed] MAS Annotation failed: {sample['excel']}\n")
            validation_failed_samples.add((exp_name, sample_name))
            refresh_tree()
            return  # ← this prevents continuing to metadata
        else:
            logger.log(f"[Validation] MAS Annotation is valid: {sample['excel']}\n")
        
        # 20260416 code review: allow score B workbook ionlist validation
        if "mass" in ionlistdf.columns:                                                      
            ion_df = ionlistdf[["mass"]]                          
        elif "fragmentation_mass" in ionlistdf.columns:                                      
            ion_df = ionlistdf[["fragmentation_mass"]].rename(columns={"fragmentation_mass":
        "mass"})                                                                             
        else:                                                     
            logger.log(f"[Validation] ion list has no 'mass' or 'fragmentation_mass' column")
            return  

        # validate the mass column contains clean float values
        if ion_df["mass"].dtype != "float64":
            logger.log(f"[Validation] ion list has invalid values: {sample['excel']}\n")
            return
        else:
            logger.log(f"[Validation] ion list validation successful: {sample['excel']}\n")

        # Check metadata
        #if not sample.get("json"):
        if not sample.get("metadata"):
            proceed = messagebox.askyesno("Metadata Missing", "No metadata found. Would you like to create it now?")
            if not proceed:
                return
            open_metadata_editor_for_sample(exp_name, sample_name)
            return

        # All checks passed → mark as validated
        linked_validated_samples.add((exp_name, sample_name))

        try:
            save_method_v1_for_sample(exp_name, sample_name, kind="MAS", auto=True) #using v1 method generator
        except Exception as e:
            logger.log(f"[Method v1][ERROR 1] Failed to auto-export Method v1 for {exp_name}/{sample_name}: {e}")
            messagebox.showwarning("Method v1", f"Validated, but failed to export Method v1.\n\n{e}")
            print("try to use old write method file. Notice that it generates .exp.json for unknown reasons")
            write_method_file(exp_name)
        messagebox.showinfo("Validated", f"Sample '{sample_name}' under '{exp_name}' is now validated.")
        logger.log(f"Sample '{sample_name}' under '{exp_name}' is now linked and validated")
        refresh_tree()

    # Merge MAS dataset into trainable (ML applicable) format. Not tracked after dataset production.
    def try_merge_selected_sample():
        exp_name, sample_name, method_path = _get_selected_context()
        if not exp_name or not sample_name:
            return

        files = experiment_projects[exp_name]["samples"][sample_name]
        #fixed 20260309 Test-A Secondary: metadata key with backward compatibility
        # v12 samples store metadata under "metadata"; old v11 legacy stored it under "json"
        if not all([files.get("csv"), files.get("excel"),
                    files.get("metadata") or files.get("json")]):
            messagebox.showerror("Error", "Sample is missing required files.")
            return

        #fixed 20260313 ML-2: ensure MAS method exists before merge
        # Legacy samples may not have a MAS method node; create one now so merge can proceed.
        if not files.get("_methods") or not any(m.get("family") == "MAS" for m in files.get("_methods", [])):
            _ensure_method_stub(files, family="MAS", sample_name=sample_name)
            save_method_v1_for_sample(exp_name=exp_name, sample_name=sample_name, kind="MAS", auto=True)

        # Ask user for output folder
        outdir = filedialog.askdirectory(title="Select output folder to save merged dataset")
        if not outdir:
            return
        # 20260416 code review fixed wrong reference on method json, should be metadata json to get derivatization field string
        metadata_path = files.get('metadata')
        if not metadata_path:
            messagebox.showerror(
                "Merge Failed\n","Could not locate the metadata."
            )
            logger.log(f"[ERROR 1] Metadata lost, merge abort")
            return

        try:
            with open(metadata_path, "r") as f:
                metadata = json.load(f)
            derivatization_type = metadata.get("Derivatization Type", "Others")
            if derivatization_type == "Others":
                value = simpledialog.askfloat("Custom Derivatization Mass","Enter custom derivatization mass (e.g., 50.1234):")
                if value is None:
                    return  # user cancelled
                else:
                    derivatization_type = value
                    # pass it into directassign_files or your protonated mass logic
                    # need to save the info back later
                
            today = datetime.now().strftime("%Y%m%d_%H%M%S")
            outname = f"{sample_name}_merged_{today}.csv"
            outpath = os.path.join(outdir, outname)
            #quick fix on path

            pre_df, iondfindex, ion_df = mspval.directassign_files(files["excel"], files["csv"],derivatization_type, debug = False)
            try:
                n_feats = int(ion_df["mass"].notna().sum()) if ("mass" in ion_df.columns) else int(len(ion_df))
            except Exception:
                n_feats = 0
            n_features_var.set(n_feats)
            _update_gate_hint()

            if add_negatives_var.get():
                try:
                    ppm = float(ppm_tol_var.get() or 10)
                except Exception:
                    ppm = 10.0
                neg_df = mspval.sample_real_negatives(
                    raw_tsv_path = files["csv"],
                    annotated_scans = pre_df["MS2scan_no"],
                    ion_df = ion_df,
                    ppm_tol = ppm,
                    min_hits = int(min_hits_var.get()),
                    max_neg_ratio = float(neg_ratio_var.get())
                )
                added = len(neg_df)
                if added:
                    pre_df = pd.concat([pre_df, neg_df], ignore_index=True)
                    print(f"[Prepare] Added {added} Non-glycan negatives (ppm={ppm}, min_hits<{int(min_hits_var.get())}, cap {neg_ratio_var.get():.1f}×).")
                else:
                    print("[Prepare] No eligible negatives found with current gates.")

            # Calling Ion Mining 
            if ion_suggest_enable_var.get():
                global export_ion_suggestions_csv, SuggestParams, _ionmod
                if export_ion_suggestions_csv is None:
                    _ionmod = _load_ion_module()
                    if _ionmod:
                        export_ion_suggestions_csv = _ionmod.export_ion_suggestions_csv
                        SuggestParams = _ionmod.SuggestParams

                if export_ion_suggestions_csv is None:
                    messagebox.showwarning("Ion suggestions",
                                        "Ion module failed to import. Check console for the exact error.")
                    return
                else:
                    try:
                        # choose output path in the same folder as merged CSV
                        suggest_csv = os.path.join(outdir, f"{sample_name}_ion_suggestions.csv")
                        params = SuggestParams(
                            ppm=float(ion_suggest_ppm_var.get()),
                            da_floor=float(ion_suggest_dafloor_var.get()),
                            min_cluster_count=3,
                            min_support_glycan=int(ion_suggest_minsupp_var.get()),
                            top_k=int(ion_suggest_topk_var.get())
                        )
                        # majority label is Non-glycan in our pipeline
                        export_ion_suggestions_csv(
                            pre_df, ion_df, out_csv=suggest_csv, params=params,
                            label_col="Structure", majority_label="Non-glycan"
                        )
                        print(f"[Prepare] Ion mining report saved: {suggest_csv}")
                        messagebox.showinfo("Ion Mining",
                            f"Suggested ions written to:\n{os.path.basename(suggest_csv)}")
                        last_suggest_csv_var.set(suggest_csv)
                    except Exception as e:
                        messagebox.showwarning("Ion Mining", f"Suggestion failed:\n{e}")

            mspval.createnormailzedionlistcsv(iondfindex, pre_df,ion_df, outpath)
            messagebox.showinfo("Merge Complete", f"Dataset saved:\n{os.path.basename(outpath)}")
            #fixed 20260313 ML-2: update method after merge completes
            # 20260416 add trainable_csv tracking (in future we may need subnode, or separated ML method under same MAS/CGA method)
            files["trainable_csv"] = outpath 

            try:
                save_method_v1_for_sample(exp_name=exp_name, sample_name=sample_name, kind="MAS", auto=True)
                refresh_tree()
            except Exception as e_method:
                logger.log(f"[ML-2] post-merge method update failed: {e_method}")
        except Exception as e:
            messagebox.showerror("Merge Failed", f"Error:\n{str(e)}")
    # MAS to trainable merge workflow block ends
    # =========================================

    # CGA analysis result (previously called Pseudolabel) data to trainable csv functionality block
    # ---------- CGA → TRAINABLE (one-pass) ----------
    import ast, re
    from datetime import datetime

    # 20250912@mark fix MS2Scan_no missing in pseudo -> trainable csv
    SCAN_CANDIDATES = ("MS2scan_no", "MS2Scan_no", "ScanNum", "scan", "Scan", "unique_ID")

    # 20260416 code review: change to add exceptions from _NONMASS list to MASSLIKE list to avoid future contamination
    # in CGA -> Trainable 6) Feature rebuilding, when use "Reuse existing wide features CSV", this function will be triggered                                                                       
    def _infer_ion_masses_from_wide_df(df):                                              
        """Extract ion m/z values from wide-format feature CSV column headers.           
        Ion mass columns are named by their m/z value (e.g., '344.1726').                
        Non-numeric headers (metadata, labels) are skipped automatically.
        """                                                                              
        masses = []                                           
        for c in df.columns:                                                             
            try:                                              
                masses.append(float(c))                                                  
            except (ValueError, TypeError):                   
                continue
        return sorted(set(masses))

    # Keep MS2scan_no consistent and controllable in df, so later dataframe edit operation won't get crazy.
    def _ensure_scan(df, fallback=None):
        """
        Guarantee a canonical string column 'MS2scan_no' exists in df.
        If a candidate exists, copy→cast→drop dups.
        Else, if fallback (Series/array) matches length, insert it.
        Returns df (copy).
        """
        if df is None or df.empty:
            return df
        out = df.copy()
        found = next((c for c in SCAN_CANDIDATES if c in out.columns), None)
        if found:
            out["MS2scan_no"] = out[found].astype(str)
            for c in SCAN_CANDIDATES:
                if c in out.columns and c != "MS2scan_no":
                    out.drop(columns=[c], inplace=True)
            return out
        if fallback is not None and len(fallback) == len(out):
            out.insert(0, "MS2scan_no", pd.Series(fallback, index=out.index, dtype="string"))
            return out
        
        # last resort: fail early with context
        print("[CGA→Train][Error] _ensure_scan failed to assign MS2scan_no; df cols:", list(df.columns)[:20], "len=", len(df))
        logger.log(f"[CGA→Train][Error] _ensure_scan failed to guarantee MS2scan_no in df")
        raise KeyError("MS2scan_no")


    # (Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc) → F, H, N, S, G, K string compact output. MIND the ORDER of input and output
    def _comp_tuple_to_label(x):
        """
        Accepts a tuple/list/str and returns compact label like F1H4N2S3 (omit zeros).
        Order: (Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc) → H,N,S,G,K,F in label
        If you use F,H,N,S,G,K as your canonical, adjust the mapping below.
        """
        if x is None or x == "" or (isinstance(x, float) and str(x) == "nan"):
            return None
        if isinstance(x, str):
            for parser in (json.loads, ast.literal_eval):
                try:
                    x = parser(x)
                    break
                except Exception:
                    pass
            if isinstance(x, str):  # fallback simple split
                parts = [p for p in x.replace("(","").replace(")","").split(",") if p.strip()!=""]
                # safe check of not integer
                for p in parts:                                                                      
                    if float(p) % 1 != 0:                                                            
                        raise ValueError(f"Fractional monosaccharide count in composition: '{p}' (full input: {x})") 
                x = [int(float(p)) for p in parts] if parts else []
        if isinstance(x, (list, tuple)):
            # Assume library order: Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc
            # Label order (O/N both): F, H, N, S(=NeuAc), G(=NeuGc), K(=KDN)
            if len(x) < 6:
                x = list(x) + [0]*(6-len(x))
            hex, hexc, neuac, neugc, kdn, fuc = [int(v) for v in x[:6]]
            parts = []
            if fuc:   parts.append(f"F{fuc}")
            if hex:  parts.append(f"H{hex}")
            if hexc:  parts.append(f"N{hexc}")
            if neuac: parts.append(f"S{neuac}")
            if neugc: parts.append(f"G{neugc}")
            if kdn:   parts.append(f"KDN{kdn}")   # 20260517 fix B-01: 3-char KDN token (was "K{kdn}"); matches _MANUAL_EXT validator + future trainable-CSV consumers
            return "".join(parts) if parts else None
        return None

    # count peak hits only, no hit peak reference (that one is moved to score B iirc, and here count serves as "minimal hit" to assign it as glycan, <hits will be assigned as Non-glycans)
    def _count_hits_to_ionlist(peaklist, peakintensity, ion_masses, ppm: float) -> int:
        """
        Count how many reference ion_masses have at least one observed peak within +/- ppm window.
        Accepts list/ndarray OR stringified lists/tuples for peaks & intensities.
        """
        import numpy as np, json, ast

        def _to_array(x):
            # Already array-like?
            if isinstance(x, (list, tuple)):
                return np.asarray(x, dtype=float)
            # String → try JSON, else Python literal, else comma/semicolon split
            s = str(x).strip()
            for parser in (json.loads, ast.literal_eval):
                try:
                    v = parser(s)
                    return np.asarray(v, dtype=float)
                except Exception:
                    pass
            # fallback: split by common separators
            try:
                parts = [p for p in s.replace(';', ',').split(',') if p.strip() != ""]
                return np.asarray([float(p) for p in parts], dtype=float) if parts else np.asarray([], dtype=float)
            except Exception:
                return np.asarray([], dtype=float)

        peaks = _to_array(peaklist)
        ints  = _to_array(peakintensity)
        if peaks.size == 0 or ints.size == 0 or peaks.size != ints.size:
            return 0

        hits = 0
        for m in ion_masses or []:
            tol = abs(float(m)) * float(ppm) / 1e6
            if np.any(np.abs(peaks - float(m)) <= tol):
                hits += 1
        return hits
    # ========

    def build_trainable_from_CGA(
        sample_name: str,
        cga_path: str,
        ion_file_path: str | None,
        ion_sheet_name: str | None,
        salvage_path: str | None,
        thresholds: dict,
        neg_opts: dict,
        feature_mode: str,          # "extract" or "reuse"
        wide_feat_csv: str | None,  # used when feature_mode == "reuse"
        output_path: str | None,
        logger=None,
        precursor_gate_ppm: float | None = None,   # e.g., 10.0 to enable, None to keep legacy behavior
        converted_csv_path: str | None = None,     # converted TSV with MS2scan_no + protonatedmass
        include_mass_feature: bool = False,               # <— NEW #20251002
        ):
        """
        One-pass builder:
        - Load CGA (pseudolabeled) long TSV
        - Normalize labels (compact FHNSGKDN)
        - Threshold selection (min ion score, max |ppm|, Top-N)
        - Optional salvage override
        - Build/Reuse features + optional negatives
        - Save trainable CSV
        Returns: (outpath, summary_dict)
        """

        import numpy as np
        # apply CLI-friendly log declaration
        log = (logger.log if logger else print)

        log(f"[CGA→Trainable] converting sample={sample_name}")

        # 1) Load CGA TSV (long form)
        pl = robust_read_csv(cga_path, prefer_tab=True) # Always tsv
        pl = clean_cols(pl) # Prevents COM elements from Window MS
        if pl is None or pl.empty:
            log(f"[CGA→Trainable][Error 0] invalid CGA tsv from {sample_name}")
            raise RuntimeError("CGA TSV is empty or unreadable.")
        # unify case/aliases early
        aliases = {c.lower(): c for c in pl.columns}
        def has(col): return col in pl.columns
        def has_lower(col): return col.lower() in aliases

        #20251002 probably include mass to trainable datasets
        # --- [MS1 attach] make sure precursor mass is present when gate or feature needs it
        # detect CGA column (case-sensitive)
        pl_scan_col = next((c for c in SCAN_CANDIDATES if c in pl.columns), None)
        # If failed, use case-insensitve aliases. if the DataFrame has "ms2scan_no" (all lowercase), Pass 1 misses it but Pass 2 finds "ms2scan_no" in aliases and       
        # returns the original column name from the DataFrame
        if not pl_scan_col and any(k.lower() in aliases for k in SCAN_CANDIDATES):
            pl_scan_col = aliases[next(k for k in (s.lower() for s in SCAN_CANDIDATES) if k in aliases)]
        if not pl_scan_col:
            raise RuntimeError("CGA tsv lacks an MS2 scan column (e.g., 'MS2scan_no').")
        pl[pl_scan_col] = pl[pl_scan_col].astype(str)

        #debug lines
        print("[CGA→Train][DEBUG] scan_col in CGA tsv for MS2 scan no:", pl_scan_col)
        print("[CGA→Train][DEBUG] MS1 gate argument: include_mass_feature:", include_mass_feature)
        print("[CGA→Train][DEBUG] column names in CGA tsv (first 20):", pl.columns.tolist()[:20])
        print("[CGA→Train][DEBUG] 'protonatedmass' in CGA tsv before MS1 attach:", "protonatedmass" in pl.columns)
        print("[CGA→Train][DEBUG] converted_csv_path:", converted_csv_path)
        #####

        # When CGA tsv lacks protonatedmass and user checked "Include precursor mass as feature, this block fires"
        need_ms1 = bool(include_mass_feature) or (precursor_gate_ppm is not None)
        if need_ms1 and "protonatedmass" not in pl.columns:
            if not converted_csv_path or not os.path.exists(converted_csv_path):
                log(f"[CGA→Train][ERROR 0]: CGA tsv missing when protonatedmass as feature is selected")
                raise RuntimeError("MS1 mass required (gate/feature), but converted MS2 CSV/TSV is missing.")
            conv = robust_read_csv(converted_csv_path, prefer_tab=True)
            # find scan + mass columns in converted file
            conv_scan = next((c for c in SCAN_CANDIDATES if c in conv.columns), None)
            conv_mz   = next((c for c in ("protonatedmass","ProtonatedMass","precursor_mass","mz","MZ") if c in conv.columns), None)
            if not conv_scan or not conv_mz:
                log(f"[CGA→Train][ERROR 1]: Missing MS2scan_no and protonatedmass-like columns")
                raise RuntimeError("Converted CSV/TSV must include MS2scan_no and protonatedmass-like columns.")
            conv = conv[[conv_scan, conv_mz]].dropna()
            conv[conv_scan] = conv[conv_scan].astype(str)
            pl = pl.merge(conv.rename(columns={conv_mz: "protonatedmass"}), left_on=pl_scan_col, right_on=conv_scan, how="left").drop(columns=[conv_scan])
        #print("[PL→Train][dbg] after MS1 attach: 'protonatedmass' in PL =", "protonatedmass" in pl.columns)
        #if "protonatedmass" in pl.columns:
        #    print("[PL→Train][dbg] protonatedmass head:", pl["protonatedmass"].head(5).tolist())

        if precursor_gate_ppm is not None:
            # Need the library mass per composition; accept any available column name
            lib_mass_col = next((c for c in ("theoretical_mass","Mass","mass","TheoMass", "Theoretical_mass", "theoreticalmass", "TheoreticalMass") if c in pl.columns), None)
            if lib_mass_col is None:
                log(f"[CGA→Train][ERROR 1]: Missing theoretical mass column")
                raise RuntimeError("CGA lack a theoretical mass column to compare against with.")
            #ppm calculation is healthy: denominator is theo mass pl[lib_mass_col]
            with np.errstate(divide="ignore", invalid="ignore"):
                pl["ppm_precursor"] = 1e6 * (pl["protonatedmass"] - pl[lib_mass_col]) / pl[lib_mass_col]

            before = len(pl)
            pl = pl[pl["ppm_precursor"].abs() <= float(precursor_gate_ppm)]
            after = len(pl)
            if logger: logger.log(f"[CGA→Train] precursor gate {precursor_gate_ppm} ppm: kept {after}/{before} rows")

        # 20260416 Code review: refactored to inactivate duplicated scan_col to scan_col. Keep dead block until we make sure it is fine. -> Tested. Looks good, remove in next version
        #scan_col = next((c for c in SCAN_CANDIDATES if has(c) or has_lower(c)), None)
        #if scan_col and scan_col not in pl.columns and has_lower(scan_col):
        #    scan_col = aliases[scan_col.lower()]

        # Standardize 'composition'
        if not has("composition"):
            if has("comp_str"):  # some runs write comp_str
                pl.rename(columns={"comp_str": "composition"}, inplace=True)
            elif has_lower("comp_str"):
                pl.rename(columns={aliases["comp_str"]: "composition"}, inplace=True)
            elif has("comp_tuple") or has_lower("comp_tuple"):
                ct = "comp_tuple" if has("comp_tuple") else aliases["comp_tuple"]
                pl["composition"] = pl[ct].apply(_comp_tuple_to_label)
            elif has("pseudo compositions") or has_lower("pseudo compositions"):
                pc = "pseudo compositions" if has("pseudo compositions") else aliases["pseudo compositions"]
                # if it’s a single best composition string, take it; if it’s a list, take first
                pl["composition"] = pl[pc].apply(
                    lambda v: v if isinstance(v, str) 
                    else (v[0] if isinstance(v, (list, tuple)) and v else None)
                )
        # final sanity
        if "composition" not in pl.columns:
            print("[CGA→Train][ERROR 1] CGA headers:", pl.columns.tolist()[:30])
            raise RuntimeError("No composition/comp_tuple column in CGA file.")     
        log(f"[CGA→Train] composition column: source={'composition' if has('composition') else 'derived'}, non-null={int(pl['composition'].notna().sum())}/{len(pl)}")
        # ============================================

        # 2) Thresholding / selection (Updated 20260330)
        min_score = float(thresholds.get("min_ion_score", 0.0))
        max_abs_ppm = float(thresholds.get("max_abs_ppm", 20.0))
        topn = int(thresholds.get("topn", 1)) # get top 1 composition in same MS2scan_no. Notice that score A is identical for isomers, so w/o score b, first from mspinsilico? wins.
        use_score = ("ion score" in pl.columns)
        if "ppm_error" not in pl.columns:
            pl["ppm_error"] = np.nan #20260416 changed to nan to avoid weird 9e9 in output... do we have that output? I think save issue is visible in ML part #9e9  

        # gates "did Score B select this row?" 
        def _truthy(v):
            if pd.isna(v):
                return False
            if isinstance(v, bool):
                return v
            s = str(v).strip().lower()
            return s in {"true", "1", "yes", "y", "t"}

        # gates "did Score B actually assign a composition?"
        def _nonempty_label(v):
            if pd.isna(v):
                return False
            s = str(v).strip()
            return bool(s) and s.lower() not in {"nan", "none", "null"}

        # First build the eligible pool using the same quality thresholds as before
        eligible = pl.copy()

        # when ion score is presented (should have one after CGA, it is score A), pick those > min_score, default value in GUI is 0.07
        if use_score:
            eligible = eligible[eligible["ion score"] >= min_score]
        # 20260416 nan version, if not working rollback
        if "ppm_error" in pl.columns and pl["ppm_error"].notna().any():
            eligible = eligible[eligible["ppm_error"].abs() <= max_abs_ppm]
        #eligible = eligible[eligible["ppm_error"].abs() <= max_abs_ppm] #old #9e9 version use this line

        # Fix B-25-A: when every row was filtered out by Score A eligibility,
        # the rest of this function (Score B re-ranking, legacy fallback, feature
        # build, MS1 merge, negative-pool, CSV write) all assume a non-empty `sel`
        # with Score B columns present. Defensively guarding each of those is
        # whack-a-mole; abort with an actionable RuntimeError instead. Mirrors
        # the equivalent v14db abort. C-001 lock is preserved: Score A is the
        # eligibility gate; if nothing passes Score A there is nothing for
        # Score B to re-rank.
        if eligible.empty:
            raise RuntimeError(
                f"[CGA→Train] 0 rows passed Score A eligibility "
                f"(min_ion_score={min_score}, max_abs_ppm={max_abs_ppm}).\n"
                f"This usually means the PL was produced by the Score B pipeline "
                f"(which writes ion_score=0) and the UI min_ion_score is > 0, OR "
                f"the CGA parameters used to produce this PL don't match what the "
                f"trainable build expects.\n"
                f"Please re-confirm CGA parameters / the min_ion_score threshold "
                f"and retry."
            )

        # Legacy fallback selection (current behavior), score A based selection is copied to legacy_sel
        order_cols = [pl_scan_col] + (["ion score"] if use_score else ["ppm_error"]) #replaced pl_scan_col
        ascending  = [True] + ([False] if use_score else [True])
        legacy_sel = (
            eligible
            .sort_values(order_cols, ascending=ascending)
            .groupby(pl_scan_col, as_index=False) #replaced pl_scan_col that means get Top N of same MS2scan_no
            .head(topn)
            .copy()
        )

        # New: prefer Score B-selected rows when available
        has_scoreb_selected = "score_b_selected" in eligible.columns

        if has_scoreb_selected:
            # Fix B-25-B: cast mask to bool BEFORE indexing. On older pandas (1.x / early 2.x),
            # `Series.map(callable)` on an empty Series loses dtype and returns object-dtype,
            # then `df[empty_object_series]` is interpreted as `df[[]]` (column-selection)
            # rather than a boolean row-mask — silently dropping all columns. astype(bool)
            # + .loc[] removes both ambiguities.
            _mask = eligible["score_b_selected"].map(_truthy).astype(bool)
            scoreb_sel = eligible.loc[_mask].copy()

            # Fix B-25-C: re-check column presence on scoreb_sel directly, not on a
            # value cached from eligible.columns. The cached check went stale during the
            # older-pandas column-drop quirk fixed above; rechecking here is cheap insurance
            # in case any future code path drops the column between eligible and scoreb_sel.
            if "selected_composition" in scoreb_sel.columns:
                scoreb_sel = scoreb_sel[scoreb_sel["selected_composition"].map(_nonempty_label)].copy()
                scoreb_sel["Structure"] = scoreb_sel["selected_composition"].astype(str).str.strip()
            else:
                # defensive fallback to avoid header inconsistency
                if "composition" in scoreb_sel.columns:
                    scoreb_sel["Structure"] = scoreb_sel["composition"].astype(str).str.strip()

            # one Score B-selected row per scan; use rank if available
            if not scoreb_sel.empty:
                if "score_b_rank" in scoreb_sel.columns:
                    scoreb_sel = (
                        scoreb_sel
                        .sort_values([pl_scan_col, "score_b_rank"], ascending=[True, True]) #replaced pl_scan_col
                        .groupby(pl_scan_col, as_index=False) #pl_scan_col
                        .head(1)
                        .copy()
                    )
                else:
                    scoreb_sel = (
                        scoreb_sel
                        .sort_values([pl_scan_col], ascending=[True]) #pl_scan_col
                        .groupby(pl_scan_col, as_index=False) #pl_scan_col
                        .head(1)
                        .copy()
                    )

                if "selection_source" not in scoreb_sel.columns:
                    scoreb_sel["selection_source"] = "score_b"
                else:
                    scoreb_sel["selection_source"] = scoreb_sel["selection_source"].replace("", "score_b").fillna("score_b")

                scoreb_scan_ids = set(scoreb_sel[pl_scan_col].astype(str))#pl_scan_col

                # fallback only for scans not covered by Score B-selected rows
                fallback_sel = legacy_sel[~legacy_sel[pl_scan_col].astype(str).isin(scoreb_scan_ids)].copy() #pl_scan_col

                # IMPORTANT: legacy fallback must still define Structure
                if "composition" in fallback_sel.columns:
                    fallback_sel["Structure"] = fallback_sel["composition"].astype(str).str.strip()
                    fallback_sel.loc[
                        fallback_sel["Structure"].str.lower().isin({"nan", "none", "null", ""}),
                        "Structure"
                    ] = pd.NA
                else:
                    fallback_sel["Structure"] = pd.NA

                if "selection_source" not in fallback_sel.columns:
                    fallback_sel["selection_source"] = "score_a_fallback"
                else:
                    fallback_sel["selection_source"] = (
                        fallback_sel["selection_source"]
                        .replace("", "score_a_fallback")
                        .fillna("score_a_fallback")
                    )
                log(f"[CGA→Train][DEBUG] scoreb_sel rows: {len(scoreb_sel)} null_Structure: {int(scoreb_sel['Structure'].isna().sum()) if 'Structure' in scoreb_sel.columns else 'NA'}")
                log(f"[CGA→Train][DEBUG] fallback_sel rows: {len(fallback_sel)} null_Structure: {int(fallback_sel['Structure'].isna().sum()) if 'Structure' in fallback_sel.columns else 'NA'}")

                sel = pd.concat([scoreb_sel, fallback_sel], ignore_index=True, sort=False)

                log(f"[CGA→Train][AUTODETECT] ranking method change to: Score B preferred, Score A fallback")
                log(f"[CGA→Train] Score B-selected scans used: {len(scoreb_scan_ids)}")
                log(f"[CGA→Train] Score A fallback rows used: {len(fallback_sel)}")
            else:
                sel = legacy_sel.copy()
                log("[CGA→Train][AUTODETECT] no valid Score B-selected rows found; ranking method is: Legacy Score A based")
        else:
            sel = legacy_sel.copy()
            log("[CGA→Train][AUTODETECT] Score B columns absent; ranking method is: Legacy Score A based")
        #20250912 fix critical root cause: no Structure column if we have pure pl datasets. It's composition!
        # --- normalize label column on the selection table ---
        # we want a guaranteed 'Structure' column to merge later
        if "Structure" not in sel.columns:
            if "composition" in sel.columns:
                sel["Structure"] = sel["composition"].astype(str)
            elif "comp_tuple" in sel.columns:
                sel["Structure"] = sel["comp_tuple"].apply(_comp_tuple_to_label)
            elif "pseudo compositions" in sel.columns:
                def _first_str(v):
                    if isinstance(v, str): 
                        return v
                    if isinstance(v, (list, tuple)) and v:
                        return v[0]
                    return None
                sel["Structure"] = sel["pseudo compositions"].apply(_first_str)
            else:
                # last resort — make it present to avoid KeyError, will become None on merge
                sel["Structure"] = None
        #20250912@mark
        # ---- normalize scan key to string to avoid dtype mismatches ----
        # keep scan key as string to avoid dtype mismatches
        sel[pl_scan_col] = sel[pl_scan_col].astype(str) #pl_scan_col

        #debug track
        print("[CGA→Train][DEBUG]sel cols:", [pl_scan_col, "Structure"],  #pl_scan_col
            "null_Struct:", int(sel["Structure"].isna().sum()))
        # 3) Ion list
        ion_df = _read_ion_df(ion_file_path) if ion_file_path else None
        if ion_df is not None and ion_sheet_name:
            # (optional) if your _read_ion_df can select sheet, pass it there; otherwise ignore
            pass

        # 4) Features
        # If EXTRACT (REBUILD): need long-form peaks and a numeric ion mass list
        # If REUSE : use the provided wide feature CSV and infer numeric masses (for negatives later)
        import numpy as np
        from ml_ng_utils import build_features_from_peaks_log10_plus1

        ion_masses = None
        if feature_mode == "extract":
            if ion_df is None or "mass" not in ion_df.columns:
                log(f"[ERROR 1] Require an ion list with mass column to define features for peak extraction (is file corrupted?)")
                raise RuntimeError("Ion list with a 'mass' column is required to to define features for peak extraction")
            ion_masses = ion_df["mass"].astype(float).tolist()

            # require long-form peaks to build features
            if not {"peaklist","peakintensity"}.issubset(pl.columns):
                log(f"[ERROR 1] CGA TSV requires 'peaklist' and 'peakintensity' for extracting features (is file corrupted?)")
                raise RuntimeError("CGA TSV must contain 'peaklist' and 'peakintensity' to extract and rebuild features.")
        else:  # feature_mode == "reuse"
            if not wide_feat_csv:                                                            
                log(f"[ERROR 1] Require an trainable wide feature CSV file to define features for peak extraction (is file corrupted?)")
                raise RuntimeError("Provide a wide feature CSV when feature_mode='reuse'.")
                                                                        
            # 1) load the wide feature matrix                     
            feat = robust_read_csv(wide_feat_csv)
            if feat is None or feat.empty:                                                   
                log(f"[ERROR 1] Trainable wide feature CSV is unreadable or empty (file corrupted?)")
                raise RuntimeError("Wide feature CSV unreadable or empty.")
                                                                       
            # 2) detect scan column                                                          
            scan_feat = next((c for c in SCAN_CANDIDATES if c in feat.columns), None)
            if not scan_feat:                                                                
                log(f"[ERROR 1] No MS2 scan column founr in trainable wide CSV. (file corrupted?)")
                raise RuntimeError("No scan column found in wide feature CSV.")


            feat[scan_feat] = feat[scan_feat].astype(str)                                    
        
            # 3) infer ion masses from headers (for negatives later)                         
            ion_masses = _infer_ion_masses_from_wide_df(feat)     
            if not ion_masses:                                                               
                log(f"[ERROR 1] Missing features from trainable CSV headers")
                raise RuntimeError("No numeric ion masses inferred from wide feature CSV headers.")
            
            pos_feat = _ensure_scan(feat, fallback=None) 

        #20251030 fix MS2scan_no name and type error
        pos_scans_df = sel[[pl_scan_col,"peaklist","peakintensity"]].rename(columns={pl_scan_col:"MS2scan_no"}) #pl_scan_col
        pos_scans_df["MS2scan_no"] = pos_scans_df["MS2scan_no"].astype(str)

        # Build features from extracted list
        pos_feat = build_features_from_peaks_log10_plus1(pos_scans_df, ion_masses, ppm=float(thresholds.get("ion_ppm", 10.0)))

        # --- NEW: robust scan handling even when features are empty ---
        scan_series = pos_scans_df["MS2scan_no"].astype("string")
        if pos_feat is None or pos_feat.empty:
            # Create a minimal shell so downstream code can proceed gracefully

            pos_feat = pd.DataFrame({"MS2scan_no": scan_series})
            pos_feat["_scan_fallback"] = scan_series
        else:
            # Ensure canonical scan column then record fallback for later merges
            pos_feat = _ensure_scan(pos_feat, fallback=scan_series)
            pos_feat["_scan_fallback"] = pos_feat["MS2scan_no"].astype(str).values
        # --------------------------------------------------------------
        log(f"[CGA→Train][info] features count: {0 if ion_masses is None else len(ion_masses)}")

        #debug print
        print("[CGA→Train][DEBUG] pos_feat pre-merge has:", 
        [c for c in ("MS2scan_no","Structure","_scan_fallback") if c in pos_feat.columns])

        # attach label next; do NOT subset columns yet
        # updated 20260330
        # 20260517 fix B-17: also carry label metadata columns (Glycanannotation2 / GlyToucan ID / WURCS) into pos_feat when present in CGA TSV `sel`. Without this, the finalize block at line 4432 would always produce empty metadata columns even when the upstream CGA TSV had user-supplied or B-17-Edit-1-seeded values (Codex Edit 2 flag).
        _b17_meta_cols = [c for c in ("Glycanannotation2", "GlyToucan ID", "WURCS") if c in sel.columns]
        merge_cols = [pl_scan_col, "Structure"] + (["selection_source"] if "selection_source" in sel.columns else []) + _b17_meta_cols
        pos_feat = pos_feat.merge(sel[merge_cols],
                                left_on="MS2scan_no", right_on=pl_scan_col,
                                how="left").drop(columns=[pl_scan_col]) #pl_scan_col
        # 20260417 A inspection of fixing 100% warn on this [PL→Train][dbg][mass] pos_feat missing MS2scan_no; rebuilding from fallback… /[CGA→Train][WARN] pos_feat missing MS2scan_no, run fallback...
        if "MS2scan_no_x" in pos_feat.columns:                                               
            print("[CGA→Train][DEBUG] MS2scan_no_x found after merge, renaming to MS2scan_no")                                                                         
            pos_feat.rename(columns={"MS2scan_no_x": "MS2scan_no"}, inplace=True)
            pos_feat.drop(columns=["MS2scan_no_y"], errors="ignore", inplace=True)           
        elif "MS2scan_no_y" in pos_feat.columns:                  
            print("[CGA→Train][DEBUG] MS2scan_no_y found after merge (no _x), renaming to MS2scan_no")                                                                         
            pos_feat.rename(columns={"MS2scan_no_y": "MS2scan_no"}, inplace=True)

        ## extra debug to eliminate nan after score b implementation
        # --- DEBUG: inspect NaN Structure rows ---
        nan_mask = pos_feat["Structure"].isna() | (
            pos_feat["Structure"].astype(str).str.strip().str.lower().isin({"nan", "none", "null", ""})
        )

        nan_rows = pos_feat[nan_mask]

        if len(nan_rows) > 0:
            unique_scans = nan_rows["MS2scan_no"].nunique() if "MS2scan_no" in nan_rows.columns else 0
            log(f"[CGA→Train][WARN] NaN Structure rows detected: {len(nan_rows)}")
            log(f"[CGA→Train][WARN] Unique MS2 scans affected: {unique_scans}")
            # show a few examples for inspection
            preview_cols = [c for c in ["MS2scan_no", "Structure"] if c in nan_rows.columns]
            try:
                preview = nan_rows[preview_cols].head(5).to_string(index=False)
                log(f"[CGA→Train][WARN] Sample NaN rows:\n{preview}")
            except Exception:
                pass

        # --- DROP NaN Structure rows (recommended) ---
        pos_feat = pos_feat[~nan_mask].copy()
        # optional deeper trace: check if nan MS2 scans still exist in sel
        if "MS2scan_no" in pos_feat.columns:
            sel_scans = set(sel[pl_scan_col].astype(str)) #pl_scan_col
            nan_scan_ids = set(nan_rows["MS2scan_no"].astype(str))
            missing_in_sel = nan_scan_ids - sel_scans
            log(f"[CGA→Train][ERROR 1] NaN MS2scan_no missing from selection: {len(missing_in_sel)}, is file corrupted?")

        # === attach MS1 (protonatedmass) + delta_ppm into wide features when requested ===
        log(f"[CGA→Train][setting] Extract mode, including MS1 precursor mass as feature: {include_mass_feature}")

        # === Include precursor mass as a feature (rebuild) ===========================
        if include_mass_feature:
            # Build a (scan, protonatedmass) table from CGA, renaming scan to MS2scan_no
            pm_src_scan = pl_scan_col  # detected earlier from CGA
            if "protonatedmass" not in pl.columns:
                log("[CGA→Train][WARN] CGA file missing 'protonatedmass' after MS1 attach; skipping mass feature inclusion")
            else:
                pm = (pl[[pm_src_scan, "protonatedmass"]]
                        .dropna()
                        .drop_duplicates(pm_src_scan)
                        .rename(columns={pm_src_scan: "MS2scan_no"}))
                pm["MS2scan_no"] = pm["MS2scan_no"].astype(str)
                log(f"[CGA→Train][info] MS2 with MS1 precursor feature count: {len(pm)} ")
                print("[CGA→Train][DEBUG] head:", pm.head(3).to_dict("records"))

            # Make sure pos_feat actually has MS2scan_no; if not, rebuild it
            # Why MS2scan_no is missing, rebuild 100% will be triggered, is probably coming from drop(columns=[pl_scan_col], 20260417 fix attempted
            if "MS2scan_no" not in pos_feat.columns:
                log("[CGA→Train][WARN] pos_feat missing MS2scan_no, run fallback...")
                pos_feat = _ensure_scan(pos_feat, fallback=pos_feat.get("_scan_fallback"))
                log(f"[CGA→Train][DEBUG] after ensure_scan fallback, has MS2scan_no now?{'MS2scan_no' in pos_feat.columns}" )

            # Final guard (fail soft if still missing)
            if "MS2scan_no" in pos_feat.columns:
                # Cast to string to avoid dtype merge issues
                pos_feat["MS2scan_no"] = pos_feat["MS2scan_no"].astype(str)
                # Merge
                pos_feat = pos_feat.merge(pm, on="MS2scan_no", how="left")
                print("[CGA→Train][DEBUG] merged mass → pos_feat cols preview",
                    [c for c in pos_feat.columns[:25]])
                print("[CGA→Train][DEBUG] protonatedmass non-null count:",
                    int(pos_feat["protonatedmass"].notna().sum()))
            else:
                print("[CGA→Train][ERROR 1] ABORT mass merge: no MS2scan_no in pos_feat")
        # ============================================================================

        # ensure scan column exists (from the earlier patch) then cast to str
        if "MS2scan_no" not in pos_feat.columns and len(pos_feat) == len(pos_scans_df):
            pos_feat.insert(0, "MS2scan_no", pos_scans_df["MS2scan_no"].values)
        pos_feat["MS2scan_no"] = pos_feat["MS2scan_no"].astype(str)

        if "Structure" not in pos_feat.columns:
            pos_feat = pos_feat.merge(
                sel[[pl_scan_col, "Structure"]],
                left_on="MS2scan_no",
                right_on=pl_scan_col,
                how="left",
                suffixes=("", "_pl")  # avoid _x/_y confusion
            ).drop(columns=[pl_scan_col]) #pl_scan_col
        else:
            # Normalize any legacy duplicates from previous runs
            if "Structure_x" in pos_feat.columns and "Structure_y" in pos_feat.columns:
                pos_feat["Structure"] = pos_feat["Structure_x"].fillna(pos_feat["Structure_y"])
                pos_feat.drop(columns=["Structure_x","Structure_y"], inplace=True, errors="ignore")
            elif "Structure_pl" in pos_feat.columns:
                pos_feat["Structure"] = pos_feat.get("Structure").fillna(pos_feat["Structure_pl"])
                pos_feat.drop(columns=["Structure_pl"], inplace=True, errors="ignore")
        #debug print
        print("[CGA→Train][DEBUG] pos_feat post-merge has:", 
            [c for c in ("MS2scan_no","Structure","Structure_pl","Structure_x","Structure_y","_scan_fallback") 
            if c in pos_feat.columns])


        # 5) Optionally add negatives (easy non-glycan) using the same ion set
        NEG_LABEL = "Non-glycan"  # keep consistent with your ML prep; used by stage-5 negatives + stage-6 summary
        pos_ids_set = set(sel[pl_scan_col].astype(str)) #pl_scan_col
        add_negs = bool(neg_opts.get("enable", False))
        final_df = pos_feat.copy()

        if add_negs:
            # --- params ---
            min_hits  = int(neg_opts.get("min_hits", 3))
            ng_ppm    = float(neg_opts.get("ion_ppm", thresholds.get("ion_ppm", 10.0)))
            max_ratio = float(neg_opts.get("max_ratio", 3.0))
            mode      = str(neg_opts.get("sampling", "random")).lower()   # "random" | "first"
            seed      = int(neg_opts.get("seed", 42))

            # --- sanity: ion masses source ---
            if ion_df is None or ion_df.empty:
                print(f"[DEBUG] ion_df is missing after reference, try ion_masses instead")
                ion_df = pd.DataFrame({"mass": ion_masses})

            # --- collect candidate negatives from long pseudolabel TSV ---
            #   (scans with < min_hits glycan-ion matches and not already positive)
            neg_rows = []
            for r in pl[[pl_scan_col, "peaklist", "peakintensity"]].dropna().itertuples(index=False): #pl_scan_col
                sid = str(getattr(r, pl_scan_col)) #pl_scan_col
                if sid in pos_ids_set:
                    continue
                hits = _count_hits_to_ionlist(
                    getattr(r, "peaklist"), getattr(r, "peakintensity"),
                    ion_df["mass"].tolist(), ppm=ng_ppm
                )
                if hits < min_hits:
                    neg_rows.append((sid, getattr(r, "peaklist"), getattr(r, "peakintensity")))

            # --- optional shuffle FIRST, then cap by ratio ---
            if neg_rows:
                if mode == "random":
                    import random as _rnd
                    _rnd.Random(seed).shuffle(neg_rows)   # in-place
                n_pos  = int((pos_feat["Structure"] != "None").sum()) if "Structure" in pos_feat.columns else len(pos_feat)
                keep   = min(len(neg_rows), int(max_ratio * max(1, n_pos)))
                neg_rows = neg_rows[:keep]

                # build DataFrame only for kept rows
                neg_df = pd.DataFrame(neg_rows, columns=["MS2scan_no", "peaklist", "peakintensity"])
                neg_df["MS2scan_no"] = neg_df["MS2scan_no"].astype(str)

                # DEBUG: prove we have scan ids before featurizing
                log(f"[CGA→Train][DEBUG][Neg sampling] Row counts after capping: {len(neg_df)}")
                print(f"[DEBUG][Neg sampling] first scans: ", neg_df["MS2scan_no"].head(5).tolist())

                # --- build features for negatives (only kept rows) ---
                neg_feat = build_features_from_peaks_log10_plus1(neg_df, ion_masses, ppm=ng_ppm)
                neg_feat = _ensure_scan(neg_feat, fallback=neg_df["MS2scan_no"])
                neg_feat["Structure"] = NEG_LABEL
                neg_feat["_scan_fallback"] = neg_feat["MS2scan_no"].astype(str).values

                print("[CGA→Train][DEBUG][Neg sampling] after features:",
                    [c for c in ("MS2scan_no","Structure","_scan_fallback") if c in neg_feat.columns],
                    "first scans:", neg_feat["MS2scan_no"].head(5).tolist())
                
        # 6) Finalize + save
        if output_path is None or output_path.strip() == "":
            outdir = os.path.dirname(cga_path)
            # suffix reflects whether MS1 is part of the feature set
            mass_suffix = "_withMass" if include_mass_feature else "_noMass"
            outname = f"{sample_name}_trainable_fromPL_{datetime.now().strftime('%Y%m%d_%H%M%S')}{mass_suffix}.csv"
            output_path = os.path.join(outdir, outname)

        log(f"[CGA→Train] include_mass_feature={include_mass_feature} → {os.path.basename(output_path)}")

        # assemble
        # --- assemble ---
        final_df = pos_feat if 'neg_feat' not in locals() else \
                pd.concat([pos_feat, neg_feat], ignore_index=True, sort=False)

        # Debug: check for protonatedmass survival before ensure_scan
        print("[CGA→Train][DEBUG] after concat: has protonatedmass?","protonatedmass" in final_df.columns)
        if "protonatedmass" in final_df.columns:
            print("[CGA→Train][DEBUG] final protonatedmass non-null:",
                int(final_df["protonatedmass"].notna().sum()))
        print("[CGA→Train][DEBUG] after concat: cols=", final_df.columns.tolist()[:20])
        if "protonatedmass" in final_df.columns:
            print("[CGA→Train][DEBUG] protonatedmass head:", final_df["protonatedmass"].head().tolist())

        # --- make sure we can reconstruct the scan column ---
        combined_fallback = final_df["_scan_fallback"].astype(str) if "_scan_fallback" in final_df.columns else None

        final_df = _ensure_scan(final_df, fallback=combined_fallback)

        # Debug after ensure_scan
        print("[CGA→Train][DEBUG] after ensure_scan: cols=", final_df.columns.tolist()[:20])
        if "protonatedmass" in final_df.columns:
            print("[CGA→Train][DEBUG] protonatedmass head:", final_df["protonatedmass"].head().tolist())

        # drop helper only after ensure_scan
        final_df.drop(columns=["_scan_fallback"], errors="ignore", inplace=True)

        log(f"[CGA→Train][DEBUG] end-before-order: cols= {final_df.columns.tolist()[:30]}")

        if not include_mass_feature:
            final_df.drop(columns=["protonatedmass","delta_ppm","ppm_precursor","precursor_gate_comp"],
                        errors="ignore", inplace=True)
            print("[CGA→Train][DEBUG][auto] MS1 feature disabled → dropped MS1 columns in trainable csv product")

        # order columns for saving
        # 20260517 fix B-17: ensure label metadata columns exist (empty-fill if absent post-merge — covers pre-Edit-1 CGA TSVs / reuse-mode trainable CSVs without metadata); reorder so they sit adjacent to Structure (label-metadata block), not at the tail of the feature block. Excludes them from feature_cols so wildcard whitelist doesn't classify them as ML features.
        for _b17_label_col in ("Glycanannotation2", "GlyToucan ID", "WURCS"):
            if _b17_label_col not in final_df.columns:
                final_df[_b17_label_col] = ""
        feature_cols = [c for c in final_df.columns if c not in ("MS2scan_no","Structure","Glycanannotation2","GlyToucan ID","WURCS")]
        if include_mass_feature and "protonatedmass" in final_df.columns:
            # 20260417 code review: fix the guarantee of protonatedmass should be at the beginning of feature rows, next to "MS2scan_no","Structure"
            feature_cols = [c for c in feature_cols if c != "protonatedmass"]
            feature_cols.insert(0, "protonatedmass")
        else:
            # explicitly drop if user disabled
            feature_cols = [c for c in feature_cols if c != "protonatedmass"]

        final_df = final_df[["MS2scan_no","Structure","Glycanannotation2","GlyToucan ID","WURCS"] + feature_cols]

        final_df.to_csv(output_path, index=False)
        log(f"[CGA→Train] saved: {output_path}")
        # summary
        classes = final_df["Structure"].value_counts().to_dict()
        summary = {
            "rows": int(len(final_df)),
            "cols": int(len(final_df.columns)),
            "classes": classes,
            "negatives_added": int((final_df["Structure"] == NEG_LABEL).sum()),
            "feature_mode": feature_mode,
            "ion_masses": len(ion_masses or []),
            "thresholds": thresholds,
            "neg_opts": neg_opts
        }
        return output_path, summary
    # ---------- end CGA (PSEUDOLABEL) → TRAINABLE ----------


    # Added 20250906 ion mining viewer window, not that useful so keep as-is for now
    def open_ion_suggestions_viewer():

        from tkinter import filedialog, messagebox, ttk

        path = last_suggest_csv_var.get().strip()
        if not path or not os.path.exists(path):
            # let user pick if we don't have a saved path yet
            path = filedialog.askopenfilename(
                title="Open ion mining suggestions CSV",
                filetypes=[("CSV files","*.csv"), ("All files","*.*")]
            )
            if not path:
                return

        try:
            df = pd.read_csv(path)
            if df.empty:
                messagebox.showinfo("Ion suggestions", "No suggestions in file.")
                return
            # sort by lift desc, then support_glycan desc
            df = df.sort_values(["lift","support_glycan"], ascending=[False, False]).reset_index(drop=True)
            top = df.head(10)
        except Exception as e:
            messagebox.showerror("Ion suggestions", f"Could not read CSV:\n{e}")
            return

        dlg = tk.Toplevel(root)
        dlg.title(f"Ion suggestions — {os.path.basename(path)}")
        dlg.geometry("620x260")
        dlg.resizable(True, True)
        dlg.grab_set()

        cols = ["mz","support_glycan","support_non","lift","odds_ratio","already_in_list","recommended"]
        tree = ttk.Treeview(dlg, columns=cols, show="headings", height=8)
        for c in cols:
            tree.heading(c, text=c)
            width = 90 if c in ("mz","lift","odds_ratio") else 110
            tree.column(c, width=width, anchor="center")

        # insert rows
        for _, r in top.iterrows():
            mz = f"{float(r.get('mz', 0.0)):.4f}"
            sg = int(r.get("support_glycan", 0))
            sn = int(r.get("support_non", 0))
            lift = f"{float(r.get('lift', 0.0)):.2f}"
            orat = f"{float(r.get('odds_ratio', 0.0)):.2f}"
            ai = str(bool(r.get("already_in_list", False)))
            rec = str(bool(r.get("recommended", True)))
            tree.insert("", "end", values=[mz, sg, sn, lift, orat, ai, rec])

        tree.pack(fill="both", expand=True, padx=8, pady=(8,4))

        def copy_mz():
            sel = tree.selection()
            if not sel: return
            mz_val = tree.item(sel[0], "values")[0]
            dlg.clipboard_clear()
            dlg.clipboard_append(mz_val)
            dlg.update()  # keep it on the clipboard
        def open_csv_folder():
            import subprocess, os
            folder = os.path.dirname(path)
            try:
                if os.name == "nt":
                    os.startfile(folder)
                elif sys.platform == "darwin":
                    subprocess.Popen(["open", folder])
                else:
                    subprocess.Popen(["xdg-open", folder])
            except Exception:
                pass

        btns = tk.Frame(dlg); btns.pack(pady=4)
        tk.Button(btns, text="Copy selected m/z", command=copy_mz).pack(side="left", padx=6)
        tk.Button(btns, text="Open folder", command=open_csv_folder).pack(side="left", padx=6)
        tk.Button(btns, text="Close", command=dlg.destroy).pack(side="left", padx=6)
    # ==================

    # Run CGA analysis block
    # Create jsonl as run log when a CGA analysis is completed (generates scored A + B(optionally) tsv)
    def append_runlog(files: dict, entry: dict):
        """Append one JSON line to a per-sample ops log."""
        try:
            # Prefer insilico/converted csv folder, else cwd
            base = (files.get("insilico_csv") or files.get("csv") or os.getcwd())
            folder = os.path.dirname(base)
            logdir = os.path.join(folder, "logs")
            os.makedirs(logdir, exist_ok=True)
            # Per-sample rolling file is simple; switch to dated if you prefer
            sample = "sample"
            # try to extract a readable name from any path we have
            for k in ("csv", "insilico_csv"):
                p = files.get(k)
                if p:
                    sample = os.path.splitext(os.path.basename(p))[0]
                    break
            logfile = os.path.join(logdir, f"{sample}_ops.jsonl")
            with open(logfile, "a", encoding="utf-8") as f:
                f.write(json.dumps(entry, ensure_ascii=False) + "\n")
        except Exception:
            traceback.print_exc()

    # 20260417 code review: keep the function name for now is fine
    def run_pseudolabeling(sample_name: str, files: dict, meta_overrides: dict, parent=None,
                        ppm_value: float = 20.0, keep_top_n_per_scan: int = 3,
                        ion_ppm: float = 10.0, anchors_required: int = 2):
        import pandas as pd, numpy as np, traceback, os

        """
        label_style: "short" -> A/s/p   ;  "long" -> A/Sul/Phos
        """

        csv_path = files.get("csv")
        ins_path = files.get("insilico_csv")
        ion_path = files.get("ionlist_path")

        if not csv_path or not os.path.exists(csv_path):
            logger.log(f"[WARN] Missing MS2 spectral data csv for CGA analysis")
            messagebox.showwarning("Converted CSV missing", "Link a converted CSV for this sample.")
            return
        if not ins_path or not os.path.exists(ins_path):
            logger.log(f"[WARN] Missing in-silico glycan list csv for CGA analysis")
            messagebox.showwarning("In-silico CSV missing", "Generate or link an in-silico CSV first.")
            return

        # 1) Load inputs
        df  = robust_read_csv(csv_path, prefer_tab=True)   # converted TSV
        lib = robust_read_csv(ins_path)                    # in-silico CSV

        # column heuristics (converted)
        scan_col = next((c for c in ["MS2scan_no","unique_ID","ScanNum","scan","Scan"] if c in df.columns), None)
        mass_col = next((c for c in ["protonatedmass","ProtonatedMass","precursor_mass","mz","MZ"] if c in df.columns), None)
        if not scan_col or not mass_col:
            logger.log(f"[ERROR 1] MS2 spectral data csv missing MS2scan_no and/or protonatedmass for CGA analysis.")
            messagebox.showerror("Columns not found",
                "Could not find scan/mass columns in converted file (need e.g., MS2scan_no + protonatedmass).")
            return
        
        #20250929 supports extra
        # 2) Normalize in-silico library → comp_tuple/comp_str + sorted Mass
        #    Include HexA/SO3/PO3H pass-through (create zeros if missing)
        COMP_COLS_BASE = ("Hex","HexNAc","NeuAc","NeuGc","KDN","Fuc")
        MOD_COLS = ("HexA","SO3","PO3H")

        # Ensure modifiers exist in lib even if the CSV lacked them
        for _c in MOD_COLS:
            if _c not in lib.columns:
                lib[_c] = 0
        try:
            libn = marker.normalize_insilico(
                lib,
                comp_cols=COMP_COLS_BASE,   # keep tuple = 6-core only
                mass_col="Mass",
                add_legacy_repr=True,
            ).copy()
            logger.log(f"[DEBUG] in silico library is loaded and normalized")
        except Exception as e:
            logger.log(f"[ERROR 1] in silico library normalization failed: {e}")
            messagebox.showerror("In-silico library normalization failed", f"Could not normalize the in-silico glycan library.\nSee log for details.\n\n{e}")
            return
        
        # IMPORTANT: keep modifier columns on the normalized table
        # normalize_insilico typically returns a row-per-entry frame preserving order,
        # so we can attach auxiliary columns directly by index alignment:
        for _c in MOD_COLS:
            if _c not in libn.columns and _c in lib.columns:
                libn[_c] = lib[_c].values

        libn = libn.sort_values("Mass").reset_index(drop=True)

        # 3) Precursor matching (top-N by |ppm| per scan) → tidy table
        rows = []
        for r in df[[scan_col, mass_col]].dropna().itertuples(index=False):
            scan, obs = getattr(r, scan_col), float(getattr(r, mass_col))
            hits = marker.find_compositions_for_mass(
                obs_mass=obs,
                insilico_sorted=libn,
                mass_col="Mass",
                ppm=ppm_value,
                mass_transform=None,  # Mass already protonated in your lib
                comp_cols=COMP_COLS_BASE,
            )
            if not hits.empty:
                hits = hits.reindex(hits["ppm_error"].abs().sort_values().index)
                if keep_top_n_per_scan:
                    hits = hits.head(keep_top_n_per_scan)
                hits.insert(0, "MS2scan_no", scan)

                # Keep composition tuple/string + modifier counts
                _keep_cols = [
                    "MS2scan_no","observed_mass","theoretical_mass","ppm_error",
                    "comp_str","comp_tuple",
                    # base composition counts (some versions of normalize may include them)
                    "Hex","HexNAc","NeuAc","NeuGc","KDN","Fuc",
                    # modifiers (we ensured they exist)
                    "HexA","SO3","PO3H",
                ]
                _keep_cols = [c for c in _keep_cols if c in hits.columns]
                rows.append(hits[_keep_cols])            

        #20250929
        matched = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
            columns=[
                "MS2scan_no","observed_mass","theoretical_mass","ppm_error",
                "comp_str","comp_tuple",
                "Hex","HexNAc","NeuAc","NeuGc","KDN","Fuc",
                "HexA","SO3","PO3H",
            ]
        )
        matched = matched.rename(columns={"comp_str": "composition"})

        # 20250920 Try to add ion mining to PL
        # --- helper: build pre_df (positives from PL + optional sampled negatives) ---
        def _pre_df_for_ion_suggest_from_pl(matched, csv_path, ion_df):
            # positives: take peaks from 'matched' (already carries peaklist/peakintensity in PL path)
            need = ["MS2scan_no","peaklist","peakintensity"]
            # error prevention, but why GPT add this? Safer but I didn't ask for it
                # normalize the scan column name if it differs
            if "MS2scan_no" not in matched.columns:
                if "MS2Scan_no" in matched.columns:
                    matched = matched.rename(columns={"MS2Scan_no": "MS2scan_no"})
                elif "scan" in matched.columns:
                    matched = matched.rename(columns={"scan": "MS2scan_no"})
                else:
                    raise KeyError("No MS2scan_no / MS2Scan_no / scan column in pseudolabeled matches.")
                
            pos = matched[need].dropna().copy()
            pos["Structure"] = "Glycan"   # anything != "Non-glycan" counts as glycan for miner
            # negatives: reuse existing gates
            if add_negatives_var.get():
                try:
                    ppm = float(ppm_tol_var.get() or 10.0)
                except Exception:
                    ppm = 10.0
                neg = mspval.sample_real_negatives(
                    raw_tsv_path=csv_path,
                    annotated_scans=pos["MS2scan_no"],
                    ion_df=ion_df,
                    ppm_tol=ppm,
                    min_hits=int(min_hits_var.get()),
                    max_neg_ratio=float(neg_ratio_var.get()),
                )
                if len(neg):
                    pos = pd.concat([pos, neg], ignore_index=True)
            return pos

        #20250909 neg ver
        # --- canonical label (HexA/SO3/PO3H-aware) ---
        try:
            matched["composition"] = matched.apply(
                lambda r: marker.canonical_label_from_row(r, style="short"), axis=1
            )
        except Exception as e:
            logger.log(f"[ERROR 2][CGA] failed to build CGA analysis result; reason={e}")
            print(f"[ERROR 2][CGA] failed to build CGA analysis result; reason={e}")
            matched["composition"] = ""

        # If downstream expects 'composition', mirror it (optional) #changed already above
        #matched["composition"] = matched["Predicted_Label"]

        if matched.empty:
            messagebox.showinfo("CGA Aborted", "No precursor matches within tolerance.")
            return

        # 4) Bring peaklist into matched (needed for ion scoring)

        #temp change to including ion score full format so we commented this part
        peaklist_col = next((c for c in ["peaklist","peaks","ms2_peaks","peak_list","mz_list"] if c in df.columns), None)
        peakint_col  = next((c for c in ["peakintensity","peak_intensity","intensitylist","intensities"] if c in df.columns), None)

        if peaklist_col and peakint_col:
            # merge without creating _x/_y when the scan key matches
            if scan_col == "MS2scan_no":
                matched = matched.merge(df[[scan_col, peaklist_col, peakint_col]], on="MS2scan_no", how="left")
            else:
                matched = matched.merge(df[[scan_col, peaklist_col, peakint_col]],
                                        left_on="MS2scan_no", right_on=scan_col, how="left") \
                                .drop(columns=[scan_col], errors="ignore")
            # standardize names so downstream is stable
            matched = matched.rename(columns={peaklist_col: "peaklist", peakint_col: "peakintensity"})
        else:
            # no peaks → ion scoring will be skipped
            pass

        # --- 4.5) Ion mining from PL (requires peaklist/peakintensity now present) ---
        if ion_suggest_enable_var.get() and "peaklist" in matched.columns and "peakintensity" in matched.columns:
            # ensure ion mining module is available
            global export_ion_suggestions_csv, SuggestParams, _ionmod
            if export_ion_suggestions_csv is None:
                _ionmod = _load_ion_module()
                if _ionmod:
                    export_ion_suggestions_csv = _ionmod.export_ion_suggestions_csv
                    SuggestParams = _ionmod.SuggestParams

            if export_ion_suggestions_csv is not None:
                try:
                    ion_df_for_suggest = _read_ion_df(files.get("ionlist_path"))
                    pre_df = _pre_df_for_ion_suggest_from_pl(matched, files["csv"], ion_df_for_suggest)

                    outdir = os.path.dirname(csv_path)
                    suggest_csv = os.path.join(outdir, f"{sample_name}_ion_suggestions.csv")
                    params = SuggestParams(
                        ppm=float(ion_suggest_ppm_var.get()),
                        da_floor=float(ion_suggest_dafloor_var.get()),
                        min_cluster_count=3,
                        min_support_glycan=int(ion_suggest_minsupp_var.get()),
                        top_k=int(ion_suggest_topk_var.get())
                    )
                    export_ion_suggestions_csv(
                        pre_df, ion_df_for_suggest, out_csv=suggest_csv, params=params,
                        label_col="Structure", majority_label="Non-glycan"
                    )
                    logger.log(f"[CGA] Ion mining suggestions saved: {suggest_csv}")
                    print(f"[CGA] Ion mining suggestions saved: {suggest_csv}")
                    last_suggest_csv_var.set(suggest_csv)
                    messagebox.showinfo("Ion suggestions",
                                        f"Suggested ions written to:\n{os.path.basename(suggest_csv)}")
                except Exception as e:
                    logger.log(f"[CGA][ERROR 2] Ion mining failed, skipped")
                    messagebox.showwarning("Ion suggestions", f"Suggestion failed:\n{e}")
            else:
                logger.log(f"[CGA][ERROR 2] Ion mining failed")
                messagebox.showwarning("Ion mining failed",
                                    "Ion module failed to import. Check console for the exact error.")

        # 5) Ion scoring (optional)
        ion_scoring_status, n_with_scores = "skipped", 0
        ion_df = _read_ion_df(ion_path) if ion_path else None
        print("[ion] path:", ion_path,
            "ion_df_rows:", 0 if (ion_df is None) else len(ion_df),
            "has_peaks:", bool("peaklist" in matched.columns and "peakintensity" in matched.columns))

        if ion_df is not None and not ion_df.empty and \
        "peaklist" in matched.columns and "peakintensity" in matched.columns:
            try:
                # 5a) try the module’s scorer first
                print("[ion] calling attach_ion_score_on_matched ...")
                matched2 = marker.attach_ion_score_on_matched(
                    matched_df=matched,
                    ion_df=ion_df,
                    ppm_value=ion_ppm,
                    scan_col="MS2scan_no",
                    ion_mass_col="mass",
                )
                # unify column names the module might use
                rename_map = {
                    "ion_score": "ion score",
                    "ion_hit_count": "ion hit count",
                    "ion_hits_mz": "ion hits m/z",
                    "ion_hits_intensity": "ion hits intensity",
                    "ion_hits_logI": "ion hits logI",
                    "ion_hits_relI": "ion hits relI",
                    "pseudo_compositions": "pseudo compositions",
                }
                for src, dst in rename_map.items():
                    if src in matched2.columns and dst not in matched2.columns:
                        matched2.rename(columns={src: dst}, inplace=True)

                # accept if it actually added a score column
                if "ion score" in matched2.columns or "ion hit count" in matched2.columns:
                    matched = matched2
                    ion_scoring_status = "ok"
                    n_with_scores = int((matched.get("ion hit count", 0) > 0).sum()) if "ion hit count" in matched.columns else 0
                else:
                    logger.log("[CGA][WARN] ion scorer returned no ion columns; falling back to simple score_counter")
                    print("[CGA][WARN] ion scorer returned no ion columns; falling back to simple score_counter")
                    raise RuntimeError("no_ion_columns")

            except Exception as e:
                import traceback; traceback.print_exc()
                # 5b) fallback — always produce basic ion columns with score_counter
                logger.log("[CGA→Train][EXCEPTION] Using fallback ion scoring (simple fraction of logI>1 matches / reference ion count). ")
                print("[CGA→Train][EXCEPTION] Using fallback ion scoring (simple fraction of logI>1 matches / reference ion count). ")
                try:
                    matched = _fallback_simple_ion_scoring(matched, ion_df, ion_ppm)
                    ion_scoring_status = "ok(fallback)"
                    n_with_scores = int((matched["ion hit count"] > 0).sum())
                except Exception:
                    traceback.print_exc()
                    ion_scoring_status = "failed"

        # 6) Merge back onto the original converted file (one row per composition match)

        scan_right = "MS2scan_no"
        if scan_right not in matched.columns:
            logger.log(f"[CGA][DEBUG] MS2scan_no has _x or _y in upstream drop process")
            print(f"[CGA][DEBUG] MS2scan_no has _x or _y in upstream drop process")
            for alt in ("MS2scan_no_x", "MS2scan_no_y", "ScanNum", "scan", "Scan", "unique_ID"):
                if alt in matched.columns:
                    matched = matched.rename(columns={alt: "MS2scan_no"})
                    break

        if "MS2scan_no" not in matched.columns:
            messagebox.showerror("Merge error", "No scan column found in matched table.")
            return
        # pick a composition column that exists
        comp_col_out = "composition" if "composition" in matched.columns else \
                    ("pseudo compositions" if "pseudo compositions" in matched.columns else None)

        enrich_cols = ["MS2scan_no", "theoretical_mass", "ppm_error"]#, "observed_mass"] <- mind this will still appear internally when doing calculation
        if comp_col_out:
            enrich_cols.insert(1, comp_col_out)
        for extra in ("ion score", "ion hit count", "ion hits m/z", "ion hits intensity", "ion hits logI", "ion hits relI"): #"pseudo compositions"  # if your scorer populates it
            if extra in matched.columns:
                enrich_cols.append(extra)
        #debug use
        print("matched cols:", matched.columns.tolist()[:30])
        print("df cols:", df.columns.tolist()[:30])
        # final join
        out = df.merge(
            matched[enrich_cols],
            left_on=scan_col, right_on="MS2scan_no",
            how="left"
        )
        # just in case it sneaks in from elsewhere:
        out.drop(columns=["observed_mass"], errors="ignore", inplace=True)
        # 20260517 fix B-17: seed label metadata columns with empty defaults so CGA-produced TSVs match the v5 MAS schema. Preserves user-supplied values when upstream `df` already carries them (MAS-converted-reused-for-CGA path).
        # FUTURE-API-HOOK: primary Glycosmos composition API integration point — when the external API caller/collector lands post-freeze, fill GlyToucan ID + WURCS from `selected_composition` / `composition` here, replacing the empty-string defaults below. Single composition string → (glytoucan_id, wurcs) round trip; no structural change to this block needed.
        for _b17_label_col in ("Glycanannotation2", "GlyToucan ID", "WURCS"):
            if _b17_label_col not in out.columns:
                out[_b17_label_col] = ""
       #20250930 fix win11 issue
       # 7) Save TSV next to converted CSV  (ABSOLUTE + explicit encoding)
        outdir  = os.path.dirname(os.path.abspath(csv_path))
        outname = f"{sample_name}_pseudolabels_{datetime.now().strftime('%Y%m%d_%H%M%S')}.tsv" #added HMS to avoid overwriting
        outpath = os.path.abspath(os.path.join(outdir, outname))
        out.to_csv(outpath, index=False, sep="\t", encoding="utf-8")

        # store absolute path so refresh_tree can always find it
        files["pseudolabel_csv"] = outpath
        logger.log(f"CGA analysis saved to: {outpath}")
        print(f"CGA analysis saved to: {outpath}")
        # 7) Save TSV next to converted CSV

        append_runlog(files, {
            "ts": datetime.now().isoformat(timespec="seconds"),
            "action": "pseudolabel_run",
            "sample": sample_name,
            "inputs": {"converted_csv": csv_path, "insilico_csv": ins_path, "ionlist": ion_path or "none"},
            "params": {"ppm_value": ppm_value, "keep_top_n_per_scan": keep_top_n_per_scan, "ion_ppm": ion_ppm},
            "output": {"pseudolabel_tsv": outpath},
            "summary": {
                "n_scans": int(len(df)),
                "n_labeled": int(matched["MS2scan_no"].nunique()),
                "n_rows": int(len(out)),
                "n_with_ion_scores": int(n_with_scores),
                "ion_scoring": ion_scoring_status,
            },
        })

        # ALWAYS hop back to Tk main thread for UI work
        if parent:
            subwin.after(0, lambda p=outpath: messagebox.showinfo(
                "CGA analysis complete", f"Saved and linked:\n{p}"
            ))
        subwin.after(0, refresh_tree)
        return outpath     # 20260501 fix issue (v): expose outpath
    # CGA analysis block ends
    # =============
    # 20260329 new block for score b
    def run_score_b_enrichment(
        pseudolabel_tsv_path: str,
        score_b_workbook_path: str,
        *,
        candidate_flags: dict | None = None,
        meta_overrides: dict | None = None,
        parent=None,
        ppm_tolerance: float = 20.0,
    ):
        import os
        import traceback

        if not pseudolabel_tsv_path or not os.path.exists(pseudolabel_tsv_path):
            messagebox.showwarning("Score B", "CGA TSV not found for Score B enrichment.")
            return pseudolabel_tsv_path

        if not score_b_workbook_path or not os.path.exists(score_b_workbook_path):
            messagebox.showwarning("Score B", "Score B workbook not found.")
            return pseudolabel_tsv_path

        meta_overrides = meta_overrides or {}
        charge_mode = meta_overrides.get("Mass Analyzer charge mode", "any")
        derivatization = meta_overrides.get("Derivatization Type", "any")

        try:
            from msp_CGA_structscore import enrich_cga_tsv_with_score_b
            logger.log("[CGA][Score B] enrichment")
            return enrich_cga_tsv_with_score_b(
                pseudolabel_tsv_path,
                score_b_workbook_path,
                candidate_flags=candidate_flags,
                ppm_tolerance=ppm_tolerance,
                charge_mode=charge_mode,
                derivatization=derivatization,
            )
            
        except Exception as e:
            logger.log("[CGA][Score B][ERROR 2] enrichment failed")
            print("[Score B] enrichment failed")
            traceback.print_exc()
            messagebox.showwarning(
                "Score B",
                f"Score B enrichment failed.\n\nOriginal TSV was kept.\n\nReason:\n{e}"
            )
            return pseudolabel_tsv_path
    # ==========

    # --- Assign file to experiment/sample, It's the low-level "put this file here" operation ---
    def assign_file(filetype, filepath, exp_title="Unassigned", sample_name="Unassigned"):
        if exp_title not in experiment_projects:
            experiment_projects[exp_title] = {"samples": {}}
        if sample_name not in experiment_projects[exp_title]["samples"]:
            experiment_projects[exp_title]["samples"][sample_name] = {"csv": None, "excel": None, "metadata": None,  "json": None} #20260415 added "metadata": None
        experiment_projects[exp_title]["samples"][sample_name][filetype] = filepath
        refresh_tree()
    # =============

    # metadata extractor block
    # extract raw filename from json
    def extract_rawname_from_metadata(json_path):
        try:
            meta = load_typed_json(
                json_path,
                expected_type=JSON_TYPE_METADATA,
                allow_legacy=True,
                context="[Metadata] "
            ) #20260126
            rawbase = os.path.basename(meta.get("Raw filename", ""))
            return os.path.splitext(rawbase)[0]
        except:
            return None
        
    def extract_title_from_metadata(json_path):
        try:
            meta = load_typed_json(
                                        json_path,
                                        expected_type=JSON_TYPE_METADATA,
                                        allow_legacy=True,
                                        context="[Metadata] "
                                    )#20260126
            return meta.get("Experiment Title")
        except:
            return None
    # metadata extractor block ends
    # ============================

    # --- File handlers ---
    def handle_csv_selection(filepaths):
        for path in filepaths:
            sample_id = os.path.splitext(os.path.basename(path))[0]
            assign_file("csv", path, "Unassigned", sample_id)

    def handle_excel_selection(filepaths):
        for path in filepaths:
            sample_id = os.path.splitext(os.path.basename(path))[0]
            assign_file("excel", path, "Unassigned", sample_id)

    def handle_metadata_selection(filepaths):
        for path in filepaths:
            title = extract_title_from_metadata(path) or "Unassigned"
            rawname = extract_rawname_from_metadata(path) or "Unassigned"
            # IMPORTANT: store as metadata, not json(method)
            assign_file("metadata", path, title, "Unassigned")

            # Auto-rename sample (move from 'Unassigned' to rawname)
            if title in experiment_projects:
                samples = experiment_projects[title]["samples"]
                if "Unassigned" in samples:
                    if rawname in samples:
                        messagebox.showwarning("Sample Exists", f"Sample '{rawname}' already exists. Skipping rename.")
                    else:
                        samples[rawname] = samples.pop("Unassigned")
                        refresh_tree()
            # DO NOT auto-write method here.
            # Metadata import should not create method/exp side-effects.
    # End of file handlers
    # ====================

    # 20260417 Code review: add json filter for user selection
    def select_files_generic(filetype_key, allow_multiple=False, on_select_callback=None):
        filetypes_dict = {                                                                   
            "csv": [("CSV files", "*.csv")],                                                 
            "excel": [("Excel files", "*.xls *.xlsx")],                                      
            "json": [("JSON files", "*.json")],                                              
            "method": [("Method JSON", "*.method.json"), ("JSON files", "*.json")],
            "experiment": [("Experiment JSON", "*.exp.json"), ("JSON files", "*.json")],     
            "metadata": [("Metadata JSON", "*.json")],            
            "all": [("All files", "*.*")]                                                    
        }     
        # 20260520 fix B-30: macOS-safe filetypes wrap (closes indirect path for
        # filetype_key="method" / "experiment" without touching call sites).
        filetypes = _macos_safe_filetypes(filetypes_dict.get(filetype_key, filetypes_dict["all"]))
        # 20260417 add crash prevention on certain macos with danger tkinter build
        if allow_multiple:                                                                   
            try:
                filepaths = filedialog.askopenfilenames(filetypes=filetypes)                 
            except Exception as e:                                                           
                logger.log(f"[File Dialog] crashed: {e}")
                filepaths = []                                                               
        else:                                                     
            try:                                                                             
                filepath = filedialog.askopenfilename(filetypes=filetypes)
            except Exception as e:
                logger.log(f"[File Dialog] crashed: {e}")                                    
                filepath = ""
            filepaths = [filepath] if filepath else []  
        # consider data validation from integrity checker 
        if filepaths and on_select_callback:
            on_select_callback(filepaths)
    # =============== 

    # Update experiment-level path display, currently less related to UX
    def update_status_display(exp_name):
        path_label, sample_label = experiment_status_labels.get(exp_name, (None, None))
        if path_label:
            path = experiment_method_paths.get(exp_name)
            path_label.set(f"Experiment method file: {path if path else 'None'}")

        if sample_label:
            sample_label.set(f"Sample method will be saved at: {sample_method_folder if sample_method_folder else 'None'}")

    # TREEVIEW OPERATION
    # mouse dragging event - select
    def on_drag_start(event):
        item_id = tree.identify_row(event.y)
        if not item_id:
              return
        item_text = tree.item(item_id, "text")
    
        if ":" in item_text:  # This is a file node by checking it has a : or not
            sample_id = tree.parent(item_id)
            exp_id = tree.parent(sample_id)
            if not sample_id or not exp_id:
                return  # Avoid broken context
            ft_raw = item_text.split(":")[0].strip()
            drag_data["filetype"] = normalize_ftype(ft_raw)              # <— was .lower()
            drag_data["filename"] = item_text.split(":")[1].strip()
            drag_data["item"] = item_id
            drag_data["from_sample"] = clean_sample_name(tree.item(sample_id, "text"))
            drag_data["from_exp"] = tree.item(exp_id, "text").replace("Experiment: ", "")
            #if debug: (gives raw dict)
        else:
            drag_data["item"] = None

    # mouse dragging event - release
    def on_drag_release(event):
        dest_id = tree.identify_row(event.y)
        if not dest_id or not drag_data["item"]:
            return

        dest_text = tree.item(dest_id, "text")
        # 20260417 code review: replaced with a simpler version. Old one accept symbols existing before Sample: #if not dest_text.startswith("Sample:") and not any(dest_text.startswith(sym + " Sample:") for sym in ["✅", "⚠️", "❌", "⛔"]):
        if "Sample:" not in dest_text:                            
            return  
        to_sample = clean_sample_name(dest_text)
        to_exp_id = tree.parent(dest_id)
        to_exp = tree.item(to_exp_id, "text").replace("Experiment: ", "")
        # Perform move
        move_file(
            from_exp=drag_data["from_exp"],
            from_sample=drag_data["from_sample"],
            ftype=drag_data["filetype"],
            filename=drag_data["filename"],
            to_exp=to_exp,
            to_sample=to_sample
        )
        drag_data["item"] = None
        refresh_tree()


    # mouse right-clicking event
    def on_right_click(event):
        item_id = tree.identify_row(event.y)
        if not item_id:
            return
        selected_text = tree.item(item_id, "text")
        sample_id = tree.parent(item_id)
        exp_id = tree.parent(sample_id)
        sample_name = clean_sample_name(tree.item(sample_id, "text"))#.replace("Sample: ", "").split(" (")[0]
        exp_name = tree.item(exp_id, "text").replace("Experiment: ", "")

        if ":" not in selected_text:
            return  # not a file entry
        
        # Right-click menu option: Remove file
        def remove_file(exp_name, sample_name, filetype):
            entry = experiment_projects[exp_name]["samples"][sample_name]
            entry[filetype] = None
            # Clear validation state if needed
            linked_validated_samples.discard((exp_name, sample_name))
            validation_failed_samples.discard((exp_name, sample_name))
            refresh_tree()

        # get correct keys of filetype and filename
        filetype_raw = selected_text.split(":")[0].strip().lower()
        filetype = normalize_ftype(filetype_raw)                     # <— normalize before use
        filename = selected_text.split(":")[1].strip()

        menu = tk.Menu(subwin, tearoff=0)
        move_menu = tk.Menu(menu, tearoff=0)

        for e_name, e_data in experiment_projects.items():
            for s_name in e_data["samples"]:
                move_menu.add_command(
                    label=f"{e_name} → {s_name}",
                    command=lambda en=e_name, sn=s_name: move_file(from_exp=exp_name, from_sample=clean_sample_name(sample_name), ftype=filetype, filename=filename, to_exp=en, to_sample=clean_sample_name(sn))
                )
        menu.add_cascade(label="Move to...", menu=move_menu)
        
        # Right-click on sample: enable linking
        if "Sample:" in selected_text:
            parent_id = tree.parent(item_id)
            sample_name = selected_text.replace("Sample:", "").split(" (")[0]
            exp_name = tree.item(parent_id, "text").replace("Experiment:", "")
            menu.add_command(
                label="Link and Validate Sample",
                command=lambda: link_and_validate_sample(exp_name, sample_name)
            )

        # Add remove option if it's a valid file
        menu.add_command(label=f"Remove {filetype.upper()}",command=lambda: remove_file(exp_name, clean_sample_name(sample_name), filetype))  
        menu.post(event.x_root, event.y_root)

    # move file shared by drag release and right click menu
    def move_file(from_exp, from_sample, ftype, filename, to_exp, to_sample):
        entry = experiment_projects[from_exp]["samples"][from_sample][ftype]
        if entry and os.path.basename(entry) == filename:
            experiment_projects[from_exp]["samples"][from_sample][ftype] = None
            assign_file(ftype, entry, to_exp, to_sample)
    # TREEVIEW OPERATION block ends
    # ============================

    # --- write to method ---

    def save_method_v1_for_sample(exp_name, sample_name, out_path=None, *, kind="MAS", auto=False):
        """
        Save Method v1 for a sample.

        Backward compatible:
        - old style: save_method_v1_for_sample(exp, sample, out_path, auto=False)
        - new style: save_method_v1_for_sample(exp, sample, kind="MAS", auto=True)  (out_path auto-resolved)

        kind: "MAS" or "CGA"
        auto:
        - if True and out_path is None, write into sample_method_folder (if set) else CSV folder fallback
        - if False and out_path is None, raise (caller should ask via dialog)
        """
        sample_name = clean_sample_name(sample_name)

        if exp_name not in experiment_projects or sample_name not in experiment_projects[exp_name]["samples"]:
            raise ValueError(f"Unknown sample: {exp_name}/{sample_name}")

        sample = experiment_projects[exp_name]["samples"][sample_name]

        #fixed 20260313 Test-2: pre-register stub so MAS is visible in tree even if write fails
        # Mirrors CGA behavior: CGA calls _ensure_method_stub before writing (lines 5951/5984/5997);
        # MAS had no fallback, causing CGA→MAS ordering to hide MAS node in tree.
        _ensure_method_stub(sample, family=kind, sample_name=sample_name)

        #fixed 20260312 Test-B B1.2: split cache by family, prevents cross-contamination
        # Cache is keyed by family so MAS and CGA never share the same slot.
        _cache_key = f"_method_v1_cache_{kind.upper()}"
        v1 = sample.get(_cache_key)
        if not isinstance(v1, dict):
            # Build from current tree links (strict checks inside)
            v1 = build_method_v1_from_tree(exp_name, sample_name, force_family=kind)
            sample[_cache_key] = v1

        # Update required headers/timestamps
        v1["json_type"] = "glycomsp.method"
        v1["schema_version"] = "1.0.0"
        v1["updated_utc"] = _utc_now_iso()

        # Decide output path
        if out_path is None:
            if not auto:
                raise ValueError("out_path is required when auto=False")

            base_dir = sample_method_folder or (
                os.path.dirname(sample["csv"]) if sample.get("csv") else os.getcwd()
            )
            # 20260521 fix B-31 (Codex Option B + Q7): reuse existing registered method
            # path for this (sample, family) if present; else stable name {sample}.{KIND}.method.json.
            # Was: timestamp embedded in filename → new file per save → orphan accumulation
            # (smoke4_5 evidence: 3 method.json files within 10-minute window).
            out_path = _resolve_method_out_path(sample, sample_name, kind, base_dir)

        # Write
        #fixed 20260313 Test-2: ensure parent directory exists (sample_method_folder may not exist yet)
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        # 20260521 fix B-31 sister (logger gap close): observability for overwrite-by-default
        # auto-save; record whether this call updated an existing file or created a new one.
        _path_existed = os.path.exists(out_path)
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(v1, f, indent=2, ensure_ascii=False)

        #fixed 20260312 Test-B B4: register in _methods so MAS survives exp.json save/reload
        # Previously only sample["json"] was set, leaving _methods with a path=None stub.
        # _register_or_update_method_ref also sets files["json"] for legacy compat (redundant
        # sample["json"] = out_path line removed).
        _register_or_update_method_ref(sample, out_path, family=kind)
        refresh_tree()
        logger.log(f"[Method v1][{kind}] {'Updated existing' if _path_existed else 'Created new'}: {out_path}")
        return out_path

    # Legacy writing method file. Only fires when method v1 can't execute anyway.
    def write_method_file(exp_name, auto=False):
        method = {
            "experiment": exp_name,
            "generated_on": datetime.now().strftime("%Y-%m-%d %H:%M"),
            "samples": {}
        }

        for sname, files in experiment_projects[exp_name]["samples"].items():
            #if not all([files.get("csv"), files.get("excel"), files.get("json")]):
            if not all([files.get("csv"), files.get("excel"), files.get("metadata")]): #change json ->metadata specifically
                continue
            try:
                #with open(files["json"], "r") as f:
                with open(files["metadata"], "r", encoding="utf-8") as f:
                    meta = json.load(f)
            except Exception:
                meta = {}

            entry = {
                "csv": os.path.abspath(os.path.normpath(files["csv"])),#files["csv"],  #os.path.basename(files["csv"]),
                "excel": os.path.abspath(os.path.normpath(files["excel"])),#files["excel"], #os.path.basename(files["excel"]),
                #"metadata": os.path.abspath(os.path.normpath(files["json"])),#files["json"], #os.path.basename(files["json"]),
                "metadata": os.path.abspath(os.path.normpath(files["metadata"])),
                "raw_file": os.path.abspath(os.path.normpath(meta.get("Raw filename", "not linked"))),#meta.get("Raw filename", "not linked"),
                "validated": (exp_name, sname) in linked_validated_samples
            }
            if "ion_sheet" in files:
                entry["ion_sheet"] = files["ion_sheet"]
            if "ion_sheet_file" in files:
                entry["ion_sheet_file"] = files["ion_sheet_file"]

            method["samples"][sname] = entry

            # Save per-sample method.json
            try:
                sample_method = {
                    "experiment": exp_name,
                    "samples": {sname: entry},
                    "generated_on": method["generated_on"]
                }
                if sample_method_folder:
                    sample_path = os.path.join(sample_method_folder, f"{sname}.method.json")
                else:
                    #sample_path = os.path.join(os.path.dirname(files["json"]), f"{sname}.method.json")
                    sample_path = os.path.join(os.path.dirname(files["metadata"]), f"{sname}.method.json")
                with open(sample_path, "w") as sf:
                    json.dump(sample_method, sf, indent=4)
                logger.log(f"[Method] Per-sample method saved: {sample_path}")
            except Exception as e:
                logger.log(f"[WARNING] Failed to save per-sample method for {sname}: {e}")

        # Save experiment-level .exp.json
        if auto and exp_name in experiment_method_paths:
            output_path = experiment_method_paths[exp_name]
        else:
            # 20260520 fix B-30: macOS-safe filetypes + defaultextension + post-process
            # to preserve compound .exp.json suffix on macOS where compound patterns crash Tk.
            output_path = filedialog.asksaveasfilename(
                defaultextension=_macos_safe_defaultextension(".exp.json"),
                initialfile=f"{exp_name}.exp.json",
                filetypes=_macos_safe_filetypes([("Experiment Method JSON", "*.exp.json")])
            )
            output_path = _macos_append_compound_suffix(output_path, ".exp.json")
            if not output_path:
                return
            experiment_method_paths[exp_name] = output_path

        try:
            with open(output_path, "w") as f:
                json.dump(method, f, indent=4)
            logger.log(f"[Method] Experiment method saved to: {output_path}")
        except Exception as e:
            messagebox.showerror("Save Failed", f"Could not save method file:\n{e}")

        update_status_display(exp_name)
    # ===


    # build method file from treeview
    def build_method_v1_from_tree(exp_name: str, sample_name: str, *, force_family: str = None) -> dict:
        """
        Build a Method v1 dict from the current experiment_projects tree entry.
        force_family: None | "MAS" | "CGA"
        """
        if exp_name not in experiment_projects:
            raise ValueError(f"Unknown experiment: {exp_name}")
        if sample_name not in experiment_projects[exp_name]["samples"]:
            raise ValueError(f"Unknown sample: {sample_name}")

        s = experiment_projects[exp_name]["samples"][sample_name]

        csv_path = s.get("csv")
        meta_path = s.get("metadata")
        excel_path = s.get("excel")

        ion_path = s.get("ionlist_path")
        insilico_path = s.get("insilico_csv")
        scoreb_path = s.get("score_b_workbook_path")
        pl_path = s.get("pseudolabel_csv")
        train_path = s.get("trainable_csv")
        unlabeled_path = s.get("unlabeled_csv")

        if not csv_path:
            raise ValueError("Missing converted CSV (csv).")

        # Decide family
        if force_family in ("MAS", "CGA"):
            family = force_family
        else:
            # auto: MAS if excel exists; CGA if ion/insilico exists
            family = "MAS" if excel_path else ("CGA" if (ion_path or insilico_path or pl_path) else "MAS")

        # Enforce required inputs per family (for v1 standard)
        if not meta_path:
            raise ValueError("Missing metadata JSON (metadata). Link or generate metadata first.")

        if family == "MAS":
            if not excel_path:
                raise ValueError("MAS method requires annotation Excel (excel).")
        elif family == "CGA":
            if not ion_path:
                raise ValueError("CGA method requires ion list (ionlist_path).")
            if not insilico_path:
                raise ValueError("CGA method requires in-silico glycan list (insilico_csv).")
        else:
            raise ValueError(f"Unknown family: {family}")

        # Build v1 dict (minimal but compliant)
        v1 = {
            "json_type": "glycomsp.method",
            "schema_version": "1.0.0",
            "uid": str(uuid.uuid4()),
            "created_utc": _utc_now_iso(),
            "updated_utc": _utc_now_iso(),

            "method": {
                "family": family,
                "name": f"{exp_name}:{sample_name}:{family}",
                "description": "",
                "tags": []
            },

            "sample": {
                "sample_name": sample_name,
                "experiment_title": exp_name
            },

            "inputs": {
                "converted_csv": {"path": os.path.normpath(csv_path)},
                "metadata_json": {"path": os.path.normpath(meta_path)}
            },

            "parameters": {
                "mas": {},
                "cga": {},
                "scoring": {},
                "ml": {}
            },

            # FUTURE-API-HOOK: when Glycosmos composition API integration lands post-freeze, add a sibling `glytoucan_resolution` block here with shape {"status": "manual"|"api"|"not_attempted", "api_endpoint": "https://doc.glycosmos.org/api/composition", "last_resolved_utc": <iso>}. Deferred for freeze v1.10 per Henry Q1: no schema bump without a reachable consumer.
            "artifacts": {
                "reports": []
            },

            # 20260521 B-32 (Codex defer + document, NOT a semantic change):
            # see normalize_method_json — this `validation` block is reserved
            # for future method/file integrity checks (v1.2+), NOT current GUI
            # sample-validation state (which lives in exp.json
            # `samples[].validated` via export_experiment_json).
            "validation": {
                "status": "unknown",
                "checked_utc": None,
                "items": []
            },

            "software": {
                "glycomsp": {"version": "", "commit": ""},
                "extractor": {"name": "", "version": ""}
            },

            "operator": {"name": "", "note": ""}
        }

        # Conditional inputs
        if family == "MAS":
            v1["inputs"]["annotation_excel"] = {"path": os.path.normpath(excel_path)}
        else:
            v1["inputs"]["ion_list"] = {"path": os.path.normpath(ion_path)}
            v1["inputs"]["insilico_glycan_list"] = {"path": os.path.normpath(insilico_path)}
            if scoreb_path:
                v1["inputs"]["score_b_workbook"] = {"path": os.path.normpath(scoreb_path)}
        #fixed 20260312 Test-B B1.1: gate all CGA artifacts/inputs, prevents cross-contamination
        # ion_list/insilico_glycan_list already gated above (in conditional inputs block).
        # pseudolabels_tsv, trainable_csv, unlabeled_csv are CGA pipeline outputs only;
        # MAS produces no pseudolabeling artifacts, so never include them in a MAS method.
        if family == "CGA":
            if pl_path:
                v1["artifacts"]["CGAresult_tsv"] = {"path": os.path.normpath(pl_path)} #pseudolabels_tsv -> CGAresult_tsv
            if train_path:
                v1["artifacts"]["trainable_csv"] = {"path": os.path.normpath(train_path)}
            if unlabeled_path:
                v1["artifacts"]["unlabeled_csv"] = {"path": os.path.normpath(unlabeled_path)}

        return v1
    # ==

    # load method file(s)
    def load_method_file(paths=None):
        if not paths:
            # 20260520 fix B-30: macOS-safe filetypes wrap (compound .method.json crashes macOS Tk)
            paths = filedialog.askopenfilenames(
                title="Select One or More Sample Method Files",
                filetypes=_macos_safe_filetypes([("Sample Method JSON", "*.method.json"), ("JSON files", "*.json")])
            )
        if not paths:
            return

        loaded = 0
        skipped = 0
        for path in paths:
            try:
                method = load_typed_json(
                    path,
                    expected_type=JSON_TYPE_METHOD,
                    allow_legacy=True,
                    context="[Method Import] "
                )
            except Exception as e:
                logger.log(f"[Method Import][ERROR 1] skipped unsupported or malformed json: {path}: {e}")
                messagebox.showerror("Error", f"Failed to load method file:\n{path}\n{e}")
                continue

            try:
                exp_name, sample_name, tree_entry, v1_obj = normalize_method_json(method, path)
            except Exception as e:
                skipped += 1
                logger.log(f"[Method Import][ERROR 1] skipped probably malformed or broken json: {path}: {e}")
                continue

            if exp_name not in experiment_projects:
                experiment_projects[exp_name] = {"samples": {}}

            # If already exists, skip (or you can decide overwrite policy later)
            if sample_name in experiment_projects[exp_name]["samples"]:
                skipped += 1
                logger.log(f"[Method Import][WARN] skipped {sample_name}: already exists in {exp_name}")
                continue

            experiment_projects[exp_name]["samples"][sample_name] = tree_entry

            # optional: store the normalized v1 object for later "Save method v1"
            #fixed 20260312 Test-B B1.2: family-keyed cache; family derived from the loaded v1 object
            _loaded_family = ((v1_obj or {}).get("method") or {}).get("family") or "UNKNOWN"
            experiment_projects[exp_name]["samples"][sample_name][f"_method_v1_cache_{_loaded_family.upper()}"] = v1_obj

            loaded += 1

        refresh_tree()
        msg = f"Imported {loaded} sample(s) successfully.\nSkipped: {skipped}"
        if skipped > 0:
            msg += "\n(Skipped files were already loaded or invalid)"
        logger.log(f"[Method Import] Imported {loaded} sample(s) successfully. Skipped: {skipped}")    
        messagebox.showinfo("Method Import", msg)
    # == 

    def change_experiment_method_path():
        exp = _current_exp_title() or "Unassigned" # 20260418 code review fix Unassigned bug
        # 20260520 fix B-30: macOS-safe SaveAs wrap (see save_method_v1_for_sample for rationale).
        path = filedialog.asksaveasfilename(
            title="Select path to save experiment method file",
            defaultextension=_macos_safe_defaultextension(".exp.json"),
            initialfile=f"{exp}.exp.json",
            filetypes=_macos_safe_filetypes([("Experiment Method JSON", "*.exp.json")])
        )
        path = _macos_append_compound_suffix(path, ".exp.json")
        if path:
            experiment_method_paths[exp] = path                                          
            update_status_display(exp) 

    def change_sample_method_folder():
        nonlocal sample_method_folder # 20260418 code review: global -> nonlocal
        path = filedialog.askdirectory(title="Select folder to save sample method files")
        if path:
            sample_method_folder = path
            update_status_display(_current_exp_title() or "Unassigned") # 20260418 code review

    #20250917 to avoid hard-fixing the items in exp json so the PL workflow can be saved as well
    def load_experiment_method_file():
        """Repurposed: load an experiment (.exp.json) and rebuild the tree."""
        # 20260520 fix B-30: macOS-safe filetypes wrap (compound .exp.json crashes macOS Tk)
        path = filedialog.askopenfilename(
            title="Open Experiment (.exp.json)",
            filetypes=_macos_safe_filetypes([("Experiment JSON", "*.exp.json"), ("JSON", "*.json")])
        )
        if not path:
            return
        exp_title = import_experiment_json(path)
        if exp_title in experiment_status_labels:  # 20260418 Code review fix: Calling .config(text=...) on a tuple crashes with AttributeError.    
            path_label, _ = experiment_status_labels[exp_title]                              
            if path_label:
                path_label.set(f"EXP file: {path}")  

    def save_current_experiment_method():
        """Repurposed: save the *experiment* (.exp.json), capturing manual + PL keys."""
        exp = _current_exp_title()
        if not exp:
            messagebox.showwarning("No experiment", "Please select or create an Experiment first.")
            return
        # use existing path if known; otherwise ask
        path = experiment_method_paths.get(exp)
        if not path:
            safe = exp.replace(" ", "_")
            stamp = datetime.now().strftime("%Y%m%d_%H%M%S") # introduced in 20260418 code review
            # 20260520 fix B-30: macOS-safe SaveAs wrap (see save_method_v1_for_sample for rationale).
            path = filedialog.asksaveasfilename(
                title="Save Experiment (.exp.json)",
                defaultextension=_macos_safe_defaultextension(".exp.json"),
                initialfile=f"{safe}_{stamp}.exp.json",
                filetypes=_macos_safe_filetypes([("Experiment JSON", "*.exp.json"), ("JSON", "*.json")])
            )
            path = _macos_append_compound_suffix(path, ".exp.json")
            if not path:
                return
        export_experiment_json(exp, path)
        # 20260418 Code review: add auto-mounting method folder                        
        experiment_method_paths[exp] = path
        update_status_display(exp)

    # 20260417 code review: platform-based right-click binding + fix MacOS unclickable issue by adding button 2
    def bind_right_click(widget, callback):                                              
        widget.bind("<Button-3>", callback)           # Windows & Linux right-click
        if platform.system() == "Darwin":                                                
            widget.bind("<Button-2>", callback)       # macOS two-finger/right-click
            widget.bind("<Control-Button-1>", callback)  # macOS Ctrl+click 

    # -- CGA (pseudolabeling)
    def launch_pseudo_labeling():
        #import compnewv4 as compv4  # assumes dev/test calls are guarded by if __name__ == "__main__"
        exp_name, sample_name, method_path = _get_selected_context()
        if not exp_name or not sample_name:
            messagebox.showwarning("No Selection", "Select a sample or method in the tree first.")
            return
        files = experiment_projects[exp_name]["samples"][sample_name]
        csv_path = files.get("csv")
        meta_path = files.get("metadata")

        if not csv_path or not meta_path:
            logger.log(f"[launch][CGA] Missing coupled csv and metadata json file")
            messagebox.showerror("Missing Files", "This sample must have both CSV and Metadata (.json) linked.")
            return
        
        # auto-resolve without prompting
        meta_path, meta_dict = _resolve_metadata_for_sample(files, sample_name, csv_path)

        if meta_path and meta_dict:
            files["metadata"] = meta_path  # cache for next time
        else:
            picked = filedialog.askopenfilename(
                title="Select metadata JSON (must contain Glycan Type / Charge; Derivatization optional)",
                filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
            )
            if picked:
                d = _load_json_safely(picked)
                if _is_metadata_dict(d):
                    meta_path, meta_dict = picked, d
                    files["metadata"] = picked
                else:
                    logger.log(f"[launch][CGA][WARN 1] Metadata invalid or lacking Glycan Type / Mass Analyzer charge mode")
                    messagebox.showwarning("Metadata not valid",
                                        "Selected file lacks Glycan Type / Mass Analyzer charge mode.")
                    meta_path, meta_dict = None, None
            else:
                logger.log(f"[launch][CGA][WARN 0] Metadata not found, derivatization and charges may be missing")
                messagebox.showwarning(
                    "Metadata not found",
                    "Could not locate metadata for this sample. You can proceed, but defaults may be wrong."
                )



        def on_generate(payload):
            
            # 20250910 not sure if GPT is asking here
            flags = dict(payload["flags"])
            #glycan_type = payload["metadata"]["Glycan Type"].strip().upper()
            # Generate NG/OG and link the path; also append run log
            #flags = dict(compv4.NG_flags); flags.update(payload.get("flags", {}))
            glycan_type = (payload["metadata"].get("Glycan Type") or "").strip().upper()
            # NEW: carry derivatization into flags
            deriv_raw = (payload["metadata"].get("Derivatization Type") or "").strip().lower()
            flags["derivatization_type"] = payload["metadata"].get("Derivatization Type", "PerMe")
            flags["reduced"] = ("reduced" in deriv_raw)
            # 20250925 neg mode support + hexA
            # NEW: Ionization mode (+ / −) from metadata into flags for compnewv4
            mode_meta = (payload["metadata"].get("Mass Analyzer charge mode") or "+").strip()
            flags["mode"] = "-" if mode_meta in ("-", "−") else "+"
            # NEW: Decide concrete counts for SO3 / PO3H from toggles + ranges
            # Policy: if toggle is ON, use the upper bound of the range; else 0.
            def _pick_count(toggle_key, range_key):
                if flags.get(toggle_key):
                    r = flags.get(range_key)
                    if isinstance(r, (list, tuple)) and len(r) == 2:
                        return int(r[1])
                return 0

            so3_count  = _pick_count("SO3",  "SO3_range")
            po3h_count = _pick_count("PO3H", "PO3H_range")

            # Stash chosen counts in flags so NG/OGlaunch (compnewv4) can read them
            flags["SO3_count"]  = so3_count
            flags["PO3H_count"] = po3h_count


            outdir = filedialog.askdirectory(title="Select output folder for in-silico CSV")
            if not outdir:
                return None
            outname = f"{sample_name}_insilico_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
            outpath = os.path.join(outdir, outname)

            try:
                if glycan_type == "N":
                    # compnewv4.NGlaunch should:
                    #  - read flags["mode"]
                    #  - pass so3=flags["SO3_count"], po3h=flags["PO3H_count"]
                    #  - honor HexA / HexA_range inside the generator
                    # although we didn't make any changes here, automatically passed from other functions?
                    compv4.NGlaunch(user_flags=flags, filename=outpath, debug=flags.get("debug", False))
                elif glycan_type == "O":
                    # haven't add neg mode and hexA support in OG
                    cores = payload.get("coretype") or [1,2,3,4]
                    compv4.OGlaunch(user_flags=flags, coretype=cores,filename=outpath, debug=flags.get("debug", False))
                else:
                    messagebox.showwarning("Glycan Type Missing",
                                        "Select 'N' or 'O' in the metadata panel.")
                    return None

                files["insilico_csv"] = outpath
                #20260204
                _ensure_method_stub(files, family="CGA", sample_name=sample_name)
                append_runlog(files, {
                    "ts": datetime.now().isoformat(timespec="seconds"),
                    "action": "insilico_generate",
                    "sample": sample_name,
                    "flags": flags,
                    "metadata_used": {
                        "Glycan Type": payload["metadata"].get("Glycan Type", ""),
                        "Mass Analyzer charge mode": payload["metadata"].get("Mass Analyzer charge mode", ""),
                        "Derivatization Type": payload["metadata"].get("Derivatization Type", ""),
                        "overrides_applied": payload["metadata"].get("_overrides_applied", False),
                    },
                    "inputs": {
                        "method_json": files.get("json"),
                        "metadata_json": files.get("metadata"),
                        "converted_csv": files.get("csv"),
                    },
                    "output": {"insilico_csv": outpath},
                    # Helpful to see what was actually used:
                    "chosen_counts": {"SO3": so3_count, "PO3H": po3h_count},
                })
                messagebox.showinfo("In-silico CSV generated", f"Saved and linked:\n{outpath}")
                logger.log(f"In-silico CSV generated, Saved and linked at: {outpath}")
                refresh_tree()
                return outpath  # so the window can show it immediately

            except Exception as e:
                traceback.print_exc()
                logger.log(f"Generation failed {str(e)}")
                messagebox.showerror("Generation failed", str(e))
                return None

        def on_link_existing(path):
            files["insilico_csv"] = path
             #20260204
            _ensure_method_stub(files, family="CGA", sample_name=sample_name)
            append_runlog(files, {
                "ts": datetime.now().isoformat(timespec="seconds"),
                "action": "insilico_link_existing",
                "sample": sample_name,
                "output": {"insilico_csv": path},
            })
            refresh_tree()

        def on_attach_ionlist(path):
            #20260204
            files["ionlist_path"] = path
            #20260204
            _ensure_method_stub(files, family="CGA", sample_name=sample_name)
            append_runlog(files, {
                "ts": datetime.now().isoformat(timespec="seconds"),
                "action": "add_ion_list",
                "sample": sample_name,
                "output": {"ionlist_path": path},
            })            
            refresh_tree()

        def on_attach_scoreb_workbook(path):
            if path:
                files["score_b_workbook_path"] = path
            else:
                files.pop("score_b_workbook_path", None)

            _ensure_method_stub(files, family="CGA", sample_name=sample_name)
            append_runlog(files, {
                "ts": datetime.now().isoformat(timespec="seconds"),
                "action": "set_score_b_workbook",
                "sample": sample_name,
                "output": {"score_b_workbook_path": path or ""},
            })
            refresh_tree()

        def on_start(payload):
            # keep current CGA-related selections synchronized first
            files["insilico_csv"] = payload.get("insilico_csv", "").strip()
            files["ionlist_path"] = payload.get("ionlist_path", "").strip()

            score_b_workbook_path = (payload.get("score_b_workbook_path") or "").strip()
            if score_b_workbook_path:
                files["score_b_workbook_path"] = score_b_workbook_path
            else:
                files.pop("score_b_workbook_path", None)

            out_tsv = run_pseudolabeling(
                sample_name=sample_name,
                files=files,
                meta_overrides=payload.get("metadata", {}),
                parent=root
            )

            if not out_tsv:
                # Since run_pseudolabeling won't return any value, we reuse the stored value from files["pseudolabel_csv"]
                out_tsv = files.get("pseudolabel_csv")

            # Part 2: Score B runtime enrichment
            if out_tsv and score_b_workbook_path:
                out_tsv = run_score_b_enrichment(
                    out_tsv,
                    score_b_workbook_path,
                    candidate_flags=payload.get("flags", {}),
                    meta_overrides=payload.get("metadata", {}),
                    parent=root,
                    ppm_tolerance=20.0,
                )

                append_runlog(files, {
                    "ts": datetime.now().isoformat(timespec="seconds"),
                    "action": "score_b_enrichment",
                    "sample": sample_name,
                    "input": {
                        "pseudolabel_tsv_before_enrichment": out_tsv,
                        "score_b_workbook_path": score_b_workbook_path,
                    },
                    "params": {
                        "ppm_tolerance": 20.0,
                    },
                    "output": {
                        "pseudolabel_tsv_after_enrichment": out_tsv,
                    },
                })

            if out_tsv:
                files["pseudolabel_csv"] = out_tsv

            # Always refresh file links first (so builder sees latest insilico/ion/tsv)
            refresh_tree()

            # If run failed / cancelled, stop here
            if not out_tsv:
                return

            try:
                #fixed 20260312 Test-B B2: CGA auto-save with CSV-folder fallback (mirrors MAS logic)
                # (1) Build v1 dict and keep in memory so Export buttons work before autosave
                v1 = build_method_v1_from_tree(exp_name, sample_name, force_family="CGA")
                files["_method_v1_cache_CGA"] = v1  #fixed 20260312 Test-B B1.2: family-keyed cache

                # (2) Autosave: prefer sample_method_folder, fall back to CSV directory (like MAS)
                base_dir = sample_method_folder or (
                    os.path.dirname(files["csv"]) if files.get("csv") else None
                )
                if base_dir:
                    # 20260521 fix B-31 sister (Codex Option B + Q7, folded with MAS):
                    # reuse existing registered CGA method path if present; else stable
                    # name. Was: timestamped filename → CGA orphan accumulation per
                    # scoring run, same pattern as MAS auto-save.
                    out_path = _resolve_method_out_path(files, sample_name, "CGA", base_dir)
                    #fixed 20260313 Test-2: ensure parent directory exists
                    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
                    # 20260521 fix B-31 sister (logger gap close): same observability
                    # rationale as MAS auto-save site in save_method_v1_for_sample.
                    _path_existed = os.path.exists(out_path)
                    with open(out_path, "w", encoding="utf-8") as f:
                        json.dump(v1, f, indent=2, ensure_ascii=False)
                    _register_or_update_method_ref(files, out_path, family="CGA")
                    logger.log(f"[Method v1][CGA] {'Updated existing' if _path_existed else 'Created new'}: {out_path}")
                else:
                    logger.log("[CGA] No method folder and no CSV path — CGA method JSON not auto-saved.")

            except Exception:
                traceback.print_exc()
                messagebox.showwarning(
                    "Method v1",
                    "CGA finished, but method JSON was not generated.\nSee console for details."
                )

            refresh_tree()

        # Open the setup window with the resolved metadata
        CGASetupWindow(
            root,
            meta_json_path=meta_path,
            meta_prefill=meta_dict or {},
            default_flags=compv4.NG_flags,
            on_generate=on_generate,
            on_link_existing=on_link_existing,
            on_attach_ionlist=on_attach_ionlist,
            on_attach_scoreb_workbook=on_attach_scoreb_workbook,
            on_start=on_start,
            initial_insilico=files.get("insilico_csv"),
            initial_ionlist=files.get("ionlist_path"),
            initial_scoreb_workbook=files.get("score_b_workbook_path"),
        )
    # CGA setup & data passing block ends

    # 20250911 CGA → Trainable window
    def open_pl_to_trainable_modal(root, sample_name, files, logger):
        import tkinter as tk
        from tkinter import filedialog, messagebox

        win = tk.Toplevel(root)
        win.title("CGA → Trainable (one-pass)")
        win.grab_set()

        # --- Inputs
        frm = ttk.Frame(win, padding=10); frm.pack(fill="both", expand=True)

        pseudo_var = tk.StringVar(value=files.get("pseudolabel_csv",""))
        ion_var    = tk.StringVar(value=files.get("ionlist_path",""))
        ion_sheet  = tk.StringVar(value="")
        salvage_var= tk.StringVar(value="")
        wide_var   = tk.StringVar(value=files.get("features_csv",""))

        # NEW — precursor gate controls (default ON)
        apply_precursor_gate_var = tk.BooleanVar(value=True)
        precursor_ppm_var        = tk.StringVar(value="10")   # sensible default; adjust if you prefer
        converted_csv_var        = tk.StringVar(value=files.get("csv","")) 
        # 20260418 converted_csv -> csv. Edge case since run_pseudolabeling reads the converted CSV directly from files["csv"]. External tsv or reuse mode might touch this part
        def browse(var, exts=(("All","*.*"),)):
            p = filedialog.askopenfilename(filetypes=exts)
            if p: var.set(p)

        row=0
        ttk.Label(frm, text="1) CGA analysis report TSV/CSV (long):").grid(row=row, column=0, sticky="w"); 
        ttk.Entry(frm, textvariable=pseudo_var, width=70).grid(row=row, column=1, sticky="we")
        ttk.Button(frm, text="Choose…", command=lambda: browse(pseudo_var,(("TSV/CSV","*.tsv *.csv"),))).grid(row=row, column=2); row+=1

        ttk.Label(frm, text="2) Ion sheet (CSV/XLSX) or Manual Excel:").grid(row=row, column=0, sticky="w")
        ttk.Entry(frm, textvariable=ion_var, width=70).grid(row=row, column=1, sticky="we")
        ttk.Button(frm, text="Choose…", command=lambda: browse(ion_var,(("CSV/XLSX","*.csv *.xlsx"),))).grid(row=row, column=2); row+=1
        ttk.Label(frm, text="   Sheet name (optional):").grid(row=row, column=0, sticky="e")
        ttk.Entry(frm, textvariable=ion_sheet, width=20).grid(row=row, column=1, sticky="w"); row+=1

        ttk.Label(frm, text="3) Salvage composition file (optional):").grid(row=row, column=0, sticky="w")
        ttk.Entry(frm, textvariable=salvage_var, width=70).grid(row=row, column=1, sticky="we")
        ttk.Button(frm, text="Choose…", command=lambda: browse(salvage_var,(("TSV/CSV","*.tsv *.csv"),))).grid(row=row, column=2); row+=1

        # --- Thresholds
        min_score   = tk.DoubleVar(value=0.07)
        max_ppm     = tk.DoubleVar(value=20.0)
        topn        = tk.IntVar(value=1)
        ion_ppm     = tk.DoubleVar(value=10.0)

        thresh_box = ttk.LabelFrame(frm, text="4) Thresholds / selection"); thresh_box.grid(row=row, column=0, columnspan=3, sticky="we", pady=(8,4))
        ttk.Label(thresh_box, text="min ion score").grid(row=0, column=0, sticky="e")
        ttk.Entry(thresh_box, textvariable=min_score, width=6).grid(row=0, column=1, sticky="w")
        ttk.Label(thresh_box, text="max |ppm_error|").grid(row=0, column=2, sticky="e")
        ttk.Entry(thresh_box, textvariable=max_ppm, width=6).grid(row=0, column=3, sticky="w")
        ttk.Label(thresh_box, text="Top N per scan").grid(row=0, column=4, sticky="e")
        ttk.Spinbox(thresh_box, from_=1, to=10, textvariable=topn, width=5).grid(row=0, column=5, sticky="w")
        ttk.Label(thresh_box, text="ion list ppm").grid(row=0, column=6, sticky="e")
        ttk.Entry(thresh_box, textvariable=ion_ppm, width=6).grid(row=0, column=7, sticky="w"); row+=1

        # --- Negatives
        neg_box = ttk.LabelFrame(frm, text="5) Negative sampling")
        neg_box.grid(row=row, column=0, columnspan=3, sticky="we", pady=(4,4))
        neg_enable = tk.BooleanVar(value=True)#(value=False)
        ttk.Checkbutton(neg_box, text="Add Non-glycan entries (easy negatives)", variable=neg_enable).grid(row=0, column=0, columnspan=3, sticky="w")
        max_ratio = tk.DoubleVar(value=3.0)
        min_hits  = tk.IntVar(value=3)
        ttk.Label(neg_box, text="Max ratio (neg:pos)").grid(row=1, column=0, sticky="e")
        ttk.Entry(neg_box, textvariable=max_ratio, width=5).grid(row=1, column=1, sticky="w")
        ttk.Label(neg_box, text="Keep scans with < min hits to ion list").grid(row=1, column=2, sticky="e")
        ttk.Entry(neg_box, textvariable=min_hits, width=5).grid(row=1, column=3, sticky="w"); row+=1
        # --- Negative sampling options extra ---
        neg_sampling_var = tk.BooleanVar(value=True)   # True = random, False = first
        neg_seed_var = tk.IntVar(value=42)

        # use grid (not pack) because neg_box uses grid
        neg_sampling_frame = ttk.Frame(neg_box)
        neg_sampling_frame.grid(row=2, column=0, columnspan=4, sticky="w", pady=(4, 0))

        # lay out the children with grid, too
        chk = ttk.Checkbutton(
            neg_sampling_frame,
            text="Shuffle negatives (random)",
            variable=neg_sampling_var
        )
        chk.grid(row=0, column=0, sticky="w", padx=(0, 8))

        ttk.Label(neg_sampling_frame, text="Seed:").grid(row=0, column=1, sticky="w")
        ttk.Spinbox(
            neg_sampling_frame,
            from_=-2_147_483_648, to=2_147_483_647,
            textvariable=neg_seed_var,
            width=8
        ).grid(row=0, column=2, sticky="w", padx=(4, 0))



        # --- 6) Feature building
        feat_box = ttk.LabelFrame(frm, text="6) Feature building")
        feat_box.grid(row=row, column=0, columnspan=3, sticky="we", pady=(4,8))

        mode = tk.StringVar(value="extract")
        ttk.Radiobutton(feat_box, text="Extract features from long-form peaks and build feature matrix (log10(I) + 1)",
                        variable=mode, value="extract").grid(row=0, column=0, columnspan=3, sticky="w")
        ttk.Radiobutton(feat_box, text="Reuse existing trainable CSV features",
                        variable=mode, value="reuse").grid(row=1, column=0, columnspan=3, sticky="w")
        ttk.Entry(feat_box, textvariable=wide_var, width=70).grid(row=2, column=0, sticky="we")
        ttk.Button(feat_box, text="Choose…",
                command=lambda: browse(wide_var, (("CSV","*.csv"),))).grid(row=2, column=1, sticky="w")

        # make the entry expand
        feat_box.columnconfigure(0, weight=1)
        frm.columnconfigure(0, weight=1)

        # >>> advance the row so the next frame goes *below* section 6
        row += 1

        # --- 7) Include mass as feature
        mass_box = ttk.LabelFrame(frm, text="7) Include mass as feature")
        mass_box.grid(row=row, column=0, columnspan=3, sticky="we", pady=(4,12))

        include_mass_feat_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(mass_box, text="Include precursor mass as feature",
                        variable=include_mass_feat_var).grid(row=0, column=0, sticky="w")
        mass_box.columnconfigure(0, weight=1)

        # --- Output + Run
        out_lbl = ttk.Label(frm, text="Output: (auto-named)"); out_lbl.grid(row=row+1, column=0, sticky="w")
        status  = ttk.Label(frm, text="", foreground="gray"); status.grid(row=row+1, column=1, sticky="w")

        def run_once():
            try:
                thresholds = {"min_ion_score": min_score.get(), "max_abs_ppm": max_ppm.get(),
                            "topn": topn.get(), "ion_ppm": ion_ppm.get()}
                neg_opts = {"enable": neg_enable.get(), "max_ratio": max_ratio.get(),
                            "min_hits": min_hits.get(), "ion_ppm": ion_ppm.get(),     # NEW: sampling + seed for randomness control
                            "sampling": "random" if neg_sampling_var.get() else "first","seed": int(neg_seed_var.get()),}
                # NEW: wire UI ➜ function
                pg_ppm   = float(precursor_ppm_var.get()) if apply_precursor_gate_var.get() else None
                conv_path = converted_csv_var.get().strip() or None
                #debug line
                print("[CGA→Train][DEBUG] include_mass_feature UI:", include_mass_feat_var.get())
                outpath, summary = build_trainable_from_CGA(
                    sample_name=sample_name,
                    cga_path=pseudo_var.get().strip(),
                    ion_file_path=ion_var.get().strip(),
                    ion_sheet_name=ion_sheet.get().strip() or None,
                    salvage_path=salvage_var.get().strip() or None,
                    thresholds=thresholds,
                    neg_opts=neg_opts,
                    feature_mode=mode.get(),
                    wide_feat_csv=wide_var.get().strip() or None,
                    output_path=None,
                    logger=logger,
                    # pass the gate params:
                    precursor_gate_ppm=pg_ppm,
                    converted_csv_path=conv_path,
                    include_mass_feature=bool(include_mass_feat_var.get()),
                    #converted_csv_path=converted_csv_path_entry.get().strip() or None,
                )
                status.config(text=outpath) 
                messagebox.showinfo("Done", f"Saved trainable CSV:\n{outpath}\n\nSummary:\nrows={summary['rows']} cols={summary['cols']}\nclasses={summary['classes']}")
                # 20260417 code review: add trainable csv (CGA route) with tree update
                files["trainable_csv"] = outpath
                print(f"[DEBUG] files['trainable_csv'] set to: {files.get('trainable_csv')}")
                refresh_tree()
            except Exception as e:
                import traceback; traceback.print_exc()
                messagebox.showerror("Failed", str(e))

        ttk.Button(frm, text="Build Trainable CSV", command=run_once).grid(row=row+2, column=1, pady=8)

    def try_pl_to_trainable():
        exp_name, sample_name, method_path = _get_selected_context()
        if not exp_name or not sample_name:
            messagebox.showerror("No Selection", "Please select a sample or method first.")
            return
        
        files = experiment_projects[exp_name]["samples"][sample_name]
        open_pl_to_trainable_modal(root, sample_name, files, logger)



    # --- Button panel ---
    button_frame = tk.Frame(subwin)
    button_frame.pack(pady=5)

    tk.Button(button_frame, text="Add MS2 CSV(s)", command=lambda: select_files_generic("csv", True, handle_csv_selection)).grid(row=0, column=0, padx=5)
    tk.Button(button_frame, text="Add Excels (ion list/man annotation)", command=lambda: select_files_generic("excel", True, handle_excel_selection)).grid(row=0, column=1, padx=5)
    tk.Button(button_frame, text="Add Sample using metadata", command=lambda: select_files_generic("json", True, handle_metadata_selection)).grid(row=0, column=2, padx=5)
    tk.Button(button_frame, text="Add Sample", command=add_sample).grid(row=1, column=0, padx=5)
    tk.Button(button_frame, text="Clean up empty unassigned sample tags", command=clean_unassigned_samples).grid(row=1, column=1, padx=5)
    link_button = tk.Button(button_frame, text="Link Sample", state="disabled", command=lambda: try_link_selected_sample())
    link_button.grid(row=1, column=2, padx=5) #why it was gone?
    merge_button = tk.Button(button_frame, text="Merge Sample", state="disabled", command=lambda: try_merge_selected_sample())
    merge_button.grid(row=1, column=3, padx=5)
    tk.Button(button_frame, text="Negative Sampling Setting",
          command=open_negative_options_dialog).grid(row=2, column=0, padx=5)
    tk.Button(button_frame, text="Ion Mining",
          command=open_ion_mining_dialog).grid(row=2, column=1, padx=5)
    tk.Button(button_frame, text="Browse Ion Mining result",
          command=open_ion_suggestions_viewer).grid(row=2, column=2, padx=5)
    tk.Button(button_frame, text="Load Method", command=load_method_file).grid(row=2, column=3, padx=5)
    tk.Button(button_frame, text="CGA manager", command=lambda:launch_pseudo_labeling()).grid(row=3, column=0, padx=5, pady=5) 
    #GPT said without () it only passes the function, and work only if clicked
    tk.Button(button_frame, text="CGA → Trainable",
          command=try_pl_to_trainable).grid(row=3, column=1, padx=5)


    tk.Button(subwin, text="Close", command=subwin.destroy).pack(pady=10)

    #status label
    status_frame = tk.Frame(subwin, relief=tk.SUNKEN, borderwidth=1)
    status_frame.pack(fill=tk.X, padx=10, pady=(0, 5))
    exp_status_var = tk.StringVar(value="Experiment method file: None")
    sample_status_var = tk.StringVar(value="Sample method will be saved at: None")
    tk.Label(status_frame, textvariable=exp_status_var, anchor="w").pack(fill=tk.X, pady=2)
    tk.Label(status_frame, textvariable=sample_status_var, anchor="w").pack(fill=tk.X, pady=2)
    # Store label references for this experiment
    experiment_status_labels["Unassigned"] = (exp_status_var, sample_status_var)

    btn_frame = tk.Frame(subwin)
    btn_frame.pack(fill=tk.X)

    tk.Button(btn_frame, text="Set Experiment Method File", command=change_experiment_method_path).pack(side=tk.LEFT, padx=10, pady=5)
    tk.Button(btn_frame, text="Set Sample Method Folder", command=change_sample_method_folder).pack(side=tk.LEFT, padx=10, pady=5)
    #bottom place for "global" exp method file
    tk.Button(btn_frame, text="Load .exp.json", command=load_experiment_method_file).pack(side=tk.LEFT, padx=10)
    tk.Button(btn_frame, text="Save .exp.json", command=save_current_experiment_method).pack(side=tk.LEFT)
    tk.Button(btn_frame, text="Export Method v1 (MAS)", command=lambda: _export_method_v1("MAS")).pack(side="left", padx=4)
    tk.Button(btn_frame, text="Export Method v1 (CGA)", command=lambda: _export_method_v1("CGA")).pack(side="left", padx=4)

    # --- Right-click bind ---
    #tree.bind("<Button-3>", on_right_click)
    tree.bind("<ButtonPress-1>", on_drag_start)
    tree.bind("<ButtonRelease-1>", on_drag_release)
    bind_right_click(tree, on_right_click)
    # Prepare dataset Window Setting ENDS
    # ==============================

def open_ml_analysis_window():
    #20250901 add split
    from sklearn.model_selection import train_test_split
    import copy
    import re
    import numpy as np
    #v0.9923~0.9929 utilities
    from ml_ng_utils import (
        collect_ng_candidates,
        cap_non_glycan,
        predict_with_threshold,
        build_features_from_peaks_log10_plus1,  # optional if you need it directly
    )
    from ml_ng_utils_extras import resample_by_strategy, tau_sweep_summary
    # ML Utility and dict definition block
    # Initialize _ml_state default dict (For Random Forest, in future supporting more models, edit BUILTIN_ML)
    _ml_state = {
        "current": {
            "model": {
                "type": "RandomForest",
                "n_estimators": 400,
                "max_depth": None,
                "min_samples_split": 2,
                "min_samples_leaf": 1,
                "random_state": 42,
                "class_weight": "balanced",
            },
            "train": {
                "test_size": 0.20, # 0.25 in code review  for tracking
                "val_size": 0.00, # 0.10 in code review for tracking
                "stratify": True,
                "real_world_test": True, # 20260418 False -> True
                "threshold": {"enabled": False, "tau": 0.65, "margin": 0.05},
            },
            "min_samples_per_class": 5,
            "balance": {"enabled": True, "majority_label": "Non-glycan", "majority_factor": 3},
        },
        # These are set by open_ml_params_window so Train/Test can push live updates:
        "editor_txt": None,
        "preview_txt": None,
    }

    import hashlib, datetime, platform, sys
    from typing import Optional, Dict, Any

    # 20260418 Keep first BUILTIN_ML as init, and remove others, also change reading it to snapshot_train_vars(), which captures real tkinter values
    BUILTIN_ML: Dict[str, Any] = {
        "model": {"type": "RandomForest", "n_estimators": 400, "max_depth": None, "class_weight": "balanced"}, # was 500, others are 400, in code review, set to 401 for tracking
        "split": {"test_size": 0.2, "val_size": 0.0, "random_state": 42, "stratified": True},
        "filters": {"min_samples_per_class": 5, "drop_rare_in_test": True}, # in code review, min_samples_per_class set to 15 for tracking
        "negatives": {"enabled": True, "max_ratio": 3.0},
        "features": {"exclude_cols": ["MS2scan_no", "protonatedmass"], "use_ion_suggestions": False},
        "thresholds": {"tau": 0.65, "margin": 0.05}, # in code review, tau set to 0.60 for tracking
    }
    # ML window state, grouped for better management. Duplicated are removed (lives inside ML window scope) 
    train_csv_path = None                                                                
    predict_input_path = None                                 
    model_file_path = None
    linked_exp_json = None
    effective_ml_params = {}  # what Train uses
    ml_summary_var = tk.StringVar(value="Params: (using built-ins)")

    # tk variables storing train/test + RF tk parameters (shared between settings dialog and the rest)
    n_estimators_var     = tk.IntVar(value=400)     # default trees number, for debug during code review session, change 400 to 333
    class_weight_var     = tk.BooleanVar(value=True)  #  => "balanced" if True else None
    test_split_var       = tk.DoubleVar(value=0.20)   # default test set split ratio % = 0.8 train/0.2 test, no val for glycomics (save entries)
    val_split_var        = tk.DoubleVar(value=0.00)   # validation set ratio %. For pooled datasets, for debug during code review session, set to 0.1
    min_samples_var      = tk.IntVar(value=5)         # tiny class threshold, default is 5, while more plausible count will be 8 or 10 on pooled dataset
    use_balance_var      = tk.BooleanVar(value=True)  # enable balancing pipeline, default to True
    majority_label_var   = tk.StringVar(value="Non-glycan") # Set Non-glycan as major labels (negatives) against other composition-fitted spectra entries
    majority_factor_var  = tk.IntVar(value=3)      # cap = factor * max(minor)
    use_stratify_var     = tk.BooleanVar(value=True)
    # Cap only the training fold (leave Val/Test uncapped)
    real_world_test_var  = tk.BooleanVar(value=True) # Only caps training sets. Default to True. 
    # Beta versions  (include publication data) seems using False flag but the tkinter rendering issue made it unclear -- I saw it was True until we set debug prints
    # Optional Confidence thresholding (post-prediction) 20250902
    enable_thresh_var    = tk.BooleanVar(value=False)   # default ON # 20260418 True to False
    thresh_val_var       = tk.DoubleVar(value=0.65)    # τ in [0,1]
    margin_val_var       = tk.DoubleVar(value=0.05)  # δ for majority-support check

    # ML live update utility block 
    #read current values
    def _on_change(*_):
        # removed _write_train_ui_into_state
        _emit_ml_state_changed_ui_refresh()
        # 20260518 B-09 fold-in (Codex Finding 38): trigger summary refresh on any tk-var write
        _update_ml_summary()
        #_debug_print_train_params()
    
    # for live updating all tk values. Exactly infinite monitoring loop
    for v in (test_split_var, val_split_var, min_samples_var, use_balance_var,
          majority_label_var, majority_factor_var, use_stratify_var,
          n_estimators_var, class_weight_var, real_world_test_var,
          enable_thresh_var, thresh_val_var, margin_val_var):
        v.trace_add("write", _on_change)
    # ML live update utility block ends
    # ======

    # ---- ML state bootstrapper ----
    def _ensure_ml_state():
        """Guarantee _ml_state exists with the keys we expect."""
        nonlocal _ml_state # 20260418 global -> nonlocal fix by Claude Code
        try:
            _ml_state  # noqa: F401
        except NameError:
            _ml_state = {}

        if not isinstance(_ml_state, dict):
            _ml_state = {}

        # 20260418 code review: add debug lines
        if "editor" not in _ml_state:
            print(f"[ML][DEBUG][init/fallback] editor field is empty or not exist. Copy from BUILTIN_ML")
            _ml_state["editor"] = copy.deepcopy(BUILTIN_ML)
        if "effective" not in _ml_state:
            print(f"[ML][DEBUG][init/fallback] effective field is empty or not exist. Copy from editor")
            _ml_state["effective"] = copy.deepcopy(_ml_state["editor"])
        # 20260419 code review: I can't find any source_file related field in either BUILTIN_ML or _ml_state. Future feature?
        #if "source_files" not in _ml_state:
        #    print(f"[ML][DEBUG][init?] source_files do not exist, reset to empty []")
        #    _ml_state["source_files"] = []  # paths you load/merge from
        return _ml_state

    # Handle the refresh/editing value in ML param editor
    def _emit_ml_state_changed_ui_refresh():
        """Refresh the ML Parameters editor/preview panes if that window is open."""
        _ensure_ml_state()
        t1 = _ml_state.get("editor_txt")
        t2 = _ml_state.get("preview_txt")
        if not (t1 or t2):
            return
        try:
            if t1:
                t1.configure(state="normal")
                t1.delete("1.0", "end")
                t1.insert("1.0", json.dumps(_ml_state["editor"], indent=2))
                t1.configure(state="normal")
            if t2:
                t2.configure(state="normal")
                t2.delete("1.0", "end")
                t2.insert("1.0", json.dumps(snapshot_train_vars(), indent=2))  # 20260418, now right window read directly
                #json.dumps(_ml_state.get("effective", _ml_state["editor"]), indent=2))
                t2.configure(state="normal")
        except Exception:
            pass

    # handles ml dict merge, called in train_model and mostly in open_ml_params_window
    def _merge_ml(*layers):
        out = {}
        for layer in layers:
            if not layer: 
                continue
            for k, v in layer.items():
                if isinstance(v, dict) and isinstance(out.get(k), dict):
                    out[k] = {**out[k], **v}
                else:
                    out[k] = v
        return out

    # For collecting training environments (not parameters)
    def collect_versions() -> Dict[str, str]:
        v = {
            "python": platform.python_version(),
            "platform": f"{platform.system()} {platform.release()}",
        }
        try:
            import sklearn
            v["sklearn"] = sklearn.__version__
        except Exception:
            pass
        # add your own GUI/app version constant if you keep one
        try:
            from importlib.metadata import version as _v
            v["joblib"] = _v("joblib")
            v["numpy"] = _v("numpy")
            v["pandas"] = _v("pandas")
        except Exception:
            pass
        # If you track GUI version somewhere:
        try:
            v["gms_gui"] = version  # define elsewhere if available
        except Exception:
            pass
        return v

    # Saving report and enviroments to json file along with model file
    def snapshot_training_run(artifact_dir: str,
                            effective_params: Dict[str, Any],
                            column_order: list[str],
                            classes: list[str],
                            inputs: Dict[str, Any],
                            hashes: Dict[str, Optional[str]],
                            versions: Dict[str, str]) -> str:
        os.makedirs(artifact_dir, exist_ok=True)
        payload = {
            "training_run": {
                "effective_params": effective_params,
                "column_order": column_order,
                "classes": classes,
                "versions": versions,
                "seeds": {"random_state": effective_params.get("split", {}).get("random_state", None)},
                "inputs": inputs,
                "hashes": hashes,
                "created_at": datetime.datetime.now().astimezone().isoformat(),
            }
        }
        out = os.path.join(artifact_dir, "training_run.json")
        with open(out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)
            f.write("\n")
        return out
    # ML Utility block ends
    # ====================

    #20250929
    def _drop_meta_and_get_X(df: pd.DataFrame) -> pd.DataFrame:
        """
        Returns a numeric-only feature matrix X.
        Drops obviously non-feature metadata columns like UID and Origin_*.
        Also attempts to coerce object columns to numeric if they are mostly numeric.
        """

        # Identify meta columns (those SHOULD NOT be features, others will be passed to model.fix(X,Y))
        meta_exact = {"uid", "predicted_label", "label_str", "class_str",
                    "origin_file", "origin_basename"}
        meta_prefix = ("origin_",)
        # drop meta columns
        drop_cols = []
        for c in df.columns:
            lc = str(c).strip().lower()
            if lc in meta_exact or any(lc.startswith(p) for p in meta_prefix):
                drop_cols.append(c)
        df2 = df.drop(columns=drop_cols, errors="ignore").copy()

        # try to coerce mostly numeric object columns (if over 95% is numeric, do the conversion)
        for c in list(df2.columns):
            if df2[c].dtype == object:
                coerced = pd.to_numeric(df2[c], errors="coerce")
                # keep coerced if it's largely numeric (>=95% not NaN)
                if coerced.notna().mean() >= 0.95:
                    df2[c] = coerced

        # keep strictly numeric columns for X, if a string-based feature is going to be fed, tokenize it.
        X = df2.select_dtypes(include=[np.number])
        return X    

    # Label consistency block
    # Always output consistent human-readable string like "F1H5N4S2" or "Non-glycan", regardless of whether the input was a tuple, a clean label, or a messy variant
    # Pure manual: FHNSGKDN order, no neg-mode extras
    _MANUAL_PURE = re.compile(r"^(?:F\d+)?H\d+N\d+(?:S\d+)?(?:G\d+)?(?:KDN\d+)?$", re.IGNORECASE)

    # Manual + PL extensions allowed:
    # optional F..., then H...N..., optional S/G/KDN, optional A..., optional trailing s/p suffixes (lowercase)
    _MANUAL_EXT = re.compile(r"^(?:F\d+)?H\d+N\d+(?:S\d+)?(?:G\d+)?(?:KDN\d+)?(?:A\d+)?(?:[sp]\d+)?$", re.IGNORECASE)

    def _looks_tuplelike_for_structure(x) -> bool:
        # tuple/list/array or string like "(5,4,0,0,0,1)"
        if isinstance(x, (tuple, list)):
            return True
        if isinstance(x, str) and re.fullmatch(r"\(\s*\d+(?:\s*,\s*\d+){5}\s*\)", x):
            return True
        return False

    def normalize_structure_for_training(val, parse_tuple_to_manual):
        """
        - If it's tuple-like -> convert using parse_tuple_to_manual (your pretrain_normalizer.parse_structure_to_manual)
        - If it's already manual (pure or PL-extended, e.g., includes A / trailing s/p) -> KEEP AS-IS
        - Else if it's composition-like but messy -> try parse_tuple_to_manual (it will reorder FHNSGKDN)
        - Else -> passthrough (e.g., 'Non-glycan')
        - Notice that it calls external module. Usage example: lambda x: normalize_structure_for_training(x, normalizer.parse_structure_to_manual)
        """
        s = str(val).strip()
        if _looks_tuplelike_for_structure(val):
            return parse_tuple_to_manual(val)

        if _MANUAL_EXT.match(s) or _MANUAL_PURE.match(s):
            return s  # keep A / s / p (and normal manual) untouched
        # Fallback: if it looks like a composition string (any letter+digits), try parser
        if re.fullmatch(r"(?:[A-Za-z]+?\d+)+", s):
            try:
                return parse_tuple_to_manual(s)
            except Exception:
                return s
        return s
    # Label consistency block Ends
    # ================

    # ML utilities part 2?
    # model loader that handles three file formats: .joblib, .pkl (both via joblib.load), and .skops (via skops.io.load). Returns (model, loader_name) 
    def load_model_any(model_path: str):
        ext = os.path.splitext(model_path)[1].lower()
        print(f"ext {ext}") # debug 20260420
        if ext in (".joblib", ".pkl"):
            model = joblib.load(model_path)
            loader = "joblib"
        elif ext == ".skops":
            print("loading skops model file")
            try:
                from skops.io import load as sk_load, get_untrusted_types
            except Exception as e:
                raise ImportError(
                    "This model is a .skops file but 'skops' is not installed. "
                    "Install it in this environment: pip install skops"
                ) from e
            untrusted = get_untrusted_types(file=model_path)
            model = sk_load(model_path, trusted=None)#)True)
            loader = "skops"
        else:
            raise ValueError(f"Unsupported model file extension: {ext}")
        if not hasattr(model, "predict"):
            raise TypeError(
                f"Loaded object is {type(model).__name__} and has no .predict(). "
                "Did you select the *_labelencoder.joblib by mistake?"
            )
        return model, loader

    # 20260521 fix B-33 (helpers, Codex-approved closure scope): read trained sklearn
    # version from training_run.json sidecar + current sklearn version, for use by
    # the cross-version notice (Branch 2) and version-mismatch guard (Branch 3) in
    # run_prediction. Codex correction: cannot rely on the `import sklearn` at
    # collect_versions (v14:6441) — that import is scoped to collect_versions, not
    # to run_prediction. _current_sklearn_version() isolates a local import.
    def _read_trained_sklearn(model_path):
        """Read trained-with sklearn version from training_run.json sidecar; None if unavailable."""
        tr_path = os.path.join(os.path.dirname(model_path), "training_run.json")
        if not os.path.exists(tr_path):
            return None
        try:
            with open(tr_path, "r", encoding="utf-8") as f:
                return json.load(f).get("training_run", {}).get("versions", {}).get("sklearn")
        except Exception:
            return None

    def _current_sklearn_version():
        """Return current sklearn version string; 'unknown' if sklearn import fails."""
        try:
            import sklearn
            return sklearn.__version__
        except Exception:
            return "unknown"

    # Recovers the exact list of feature column names that the model was trained on, avoid misalignment error
    def get_training_features(model, model_path: str):
        # 1) native sklearn attribute (best)
        feats = getattr(model, "feature_names_in_", None)
        # 2) skops metadata (if exported with metadata)
        if feats is None:
            meta = getattr(model, "__skops_metadata__", None)
            if isinstance(meta, dict):
                feats = meta.get("feature_names") or meta.get("feature_names_in_")
        # 3) sidecar JSON saved at train time
        if feats is None:
            sidecar = (model_path
                    .replace("_rf_model.joblib", "_features.json")
                    .replace(".skops", "_features.json"))
            if os.path.exists(sidecar):
                with open(sidecar, "r", encoding="utf-8") as f:
                    try:
                        feats = json.load(f)
                    except Exception:
                        pass
        return [str(c) for c in feats] if feats is not None else None

    # Class Balance, Split, and Capping
    def balance_and_split(
        df: pd.DataFrame,
        label_col="Structure",
        majority_label="Non-glycan",
        min_count=12,
        majority_factor=3,
        feature_exclude=("MS2scan_no","ID","Source","IUPACname(optional)","Glycanannotation2","GlyToucan ID","WURCS","unique_ID"),  # 20260517 fix B-17: WURCS added — never treated as ML feature; carried as reference metadata only.
        test_size=0.2,
        val_size=0.0,
        stratify=True,
        cap_training_only=False,     # << NEW
        random_state=42
    ):
        from sklearn.model_selection import train_test_split

        # --- 0) drop tiny classes globally so test never contains unseen labels
        counts = df[label_col].value_counts()
        keep = counts[counts >= min_count].index
        dropped_rare = df[~df[label_col].isin(keep)]
        df1 = df[df[label_col].isin(keep)].copy()

        # --- helpers
        def _feat_cols(dff):
            return [c for c in dff.columns if c not in set(feature_exclude) | {label_col}]

        def _cap_majority_df(dff, cap):
            if majority_label not in dff[label_col].unique():
                return dff, None
            maj = dff[dff[label_col] == majority_label]
            if len(maj) <= cap:
                return dff, None
            maj_keep = maj.sample(n=int(cap), random_state=random_state)
            dff2 = pd.concat([maj_keep, dff[dff[label_col] != majority_label]], axis=0)
            return dff2, int(cap)

        # --- robust split math
        ts = float(test_size); vs = float(val_size)
        if ts < 0 or vs < 0: raise ValueError("Test/Validation must be ≥ 0.")
        if ts == 0 and vs == 0: raise ValueError("At least one of test/val must be > 0.")
        holdout = ts + vs
        if holdout >= 0.999: raise ValueError(f"Test+Validation ({holdout:.2f}) must be < 1.0.")

        feat_cols = _feat_cols(df1)
        X_all, y_all = df1[feat_cols], df1[label_col]
        if y_all.nunique() < 2: raise ValueError("After filtering, fewer than 2 classes remain.")

        if not cap_training_only:
            # ----- MODE A: cap BEFORE splitting (affects all folds) -----
            if majority_label in y_all.unique():
                minor_counts = df1[df1[label_col] != majority_label][label_col].value_counts()
                max_minor = int(minor_counts.max()) if not minor_counts.empty else 0
                cap = max(majority_factor * max_minor, 1)
                df1, cap_applied = _cap_majority_df(df1, cap)
                feat_cols = _feat_cols(df1)
                X_all, y_all = df1[feat_cols], df1[label_col]
            else:
                cap_applied = None

            strat = y_all if stratify else None
            try:
                X_train, X_tmp, y_train, y_tmp = train_test_split(
                    X_all, y_all, test_size=holdout, stratify=strat, random_state=random_state
                )
            except ValueError:
                X_train, X_tmp, y_train, y_tmp = train_test_split(
                    X_all, y_all, test_size=holdout, stratify=None, random_state=random_state
                )

            if vs > 0:
                rel_test = min(max(ts / holdout, 1e-6), 1 - 1e-6)
                strat_tmp = y_tmp if stratify else None
                try:
                    X_val, X_test, y_val, y_test = train_test_split(
                        X_tmp, y_tmp, test_size=rel_test, stratify=strat_tmp, random_state=random_state
                    )
                except ValueError:
                    X_val, X_test, y_val, y_test = train_test_split(
                        X_tmp, y_tmp, test_size=rel_test, stratify=None, random_state=random_state
                    )
            else:
                X_val, y_val = X_tmp.iloc[0:0], y_tmp.iloc[0:0]
                X_test, y_test = X_tmp, y_tmp

            info = {
                "mode": "cap_before_split",
                "kept_label_counts": y_all.value_counts().to_dict(),
                "dropped_rare_counts": dropped_rare[label_col].value_counts().to_dict(),
                "majority_cap_applied_to": majority_label if cap_applied else None,
                "majority_cap": int(cap_applied) if cap_applied else None,
                "feature_cols": feat_cols,
            }
            #Future proof, can be hidden
            # after building `info` in either mode:
            info["cap_applied_to"] = info.get("majority_cap_applied_to") or info.get("train_majority_cap_applied_to")
            info["cap_value"]      = info.get("majority_cap") or info.get("train_majority_cap")

            return X_train, y_train, X_val, y_val, X_test, y_test, info

        else:
            # ----- MODE B: cap TRAINING ONLY (real-world test) -----
            strat = y_all if stratify else None
            try:
                X_train, X_tmp, y_train, y_tmp = train_test_split(
                    X_all, y_all, test_size=holdout, stratify=strat, random_state=random_state
                )
            except ValueError:
                X_train, X_tmp, y_train, y_tmp = train_test_split(
                    X_all, y_all, test_size=holdout, stratify=None, random_state=random_state
                )

            if vs > 0:
                rel_test = min(max(ts / holdout, 1e-6), 1 - 1e-6)
                strat_tmp = y_tmp if stratify else None
                try:
                    X_val, X_test, y_val, y_test = train_test_split(
                        X_tmp, y_tmp, test_size=rel_test, stratify=strat_tmp, random_state=random_state
                    )
                except ValueError:
                    X_val, X_test, y_val, y_test = train_test_split(
                        X_tmp, y_tmp, test_size=rel_test, stratify=None, random_state=random_state
                    )
            else:
                X_val, y_val = X_tmp.iloc[0:0], y_tmp.iloc[0:0]
                X_test, y_test = X_tmp, y_tmp

            # compute cap from TRAIN labels only
            if majority_label in y_train.unique():
                minor_counts_train = y_train[y_train != majority_label].value_counts()
                max_minor_train = int(minor_counts_train.max()) if not minor_counts_train.empty else 0
                cap_train = max(majority_factor * max_minor_train, 1)

                train_df = df1.loc[X_train.index].copy()
                train_df_capped, cap_applied = _cap_majority_df(train_df, cap_train)
                # rebuild X_train / y_train from capped rows
                X_train = train_df_capped[feat_cols]
                y_train = train_df_capped[label_col]
            else:
                cap_applied = None

            info = {
                "mode": "cap_training_only",
                "kept_label_counts": y_all.value_counts().to_dict(),
                "dropped_rare_counts": dropped_rare[label_col].value_counts().to_dict(),
                "train_majority_cap_applied_to": majority_label if cap_applied else None,
                "train_majority_cap": int(cap_applied) if cap_applied else None,
                "feature_cols": feat_cols,
            }
            #Future proof, can be hidden
            # after building `info` in either mode:
            info["cap_applied_to"] = info.get("majority_cap_applied_to") or info.get("train_majority_cap_applied_to")
            info["cap_value"]      = info.get("majority_cap") or info.get("train_majority_cap")
            return X_train, y_train, X_val, y_val, X_test, y_test, info
    # ==

    # --- helper: guess method-json-derived folder name for packing ---
    def _guess_method_basename_for_pack(predict_input_path: str, df: pd.DataFrame) -> str:
        """
        Find the correct *.method.json for packaging.
        Priority:
        1) explicit column in CSV (method_json / method_path / method)
        2) unique neighbor *.method.json
        3) neighbor match by sample_name / experiment_title substring
        4) ASK USER to pick one (no silent 'newest' fallback)
        5) CSV stem as last resort
        """
        import os, glob
        folder = os.path.dirname(predict_input_path)

        # 1) explicit hint column
        for col in df.columns:
            if col.lower() in ("method_json", "method_path", "method", "json_path"):
                try:
                    cand = str(df[col].dropna().iloc[0]).strip()
                    if cand and cand.lower().endswith(".method.json") and os.path.exists(cand):
                        print(f"[ML-PACK] method basename from CSV column '{col}':{Path(cand).stem}")
                        return Path(cand).stem
                except Exception:
                    pass

        # 2) neighbors
        candidates = sorted(glob.glob(os.path.join(folder, "*.method.json")))
        if len(candidates) == 1:
            print(f"[ML-PACK] method basename from unique neighbor: {Path(candidates[0]).stem}")
            return Path(candidates[0]).stem

        if len(candidates) > 1:
            # 3) try match by hints
            hints = []
            for hcol in ("sample_name", "experiment_title"):
                if hcol in df.columns:
                    try:
                        hints.append(str(df[hcol].dropna().iloc[0]))
                    except Exception:
                        pass
            for c in candidates:
                name = os.path.basename(c)
                if any(h and (h in name) for h in hints):
                    print(f"[ML-PACK] method basename from hint match '{hints}' → {Path(c).stem}")
                    return Path(c).stem

            # 4) ask user explicitly
            try:
                from tkinter import filedialog
                # 20260520 fix B-30: macOS-safe filetypes wrap (compound .method.json crashes macOS Tk)
                sel = filedialog.askopenfilename(
                    title="Select method.json for packaging",
                    initialdir=folder,
                    filetypes=_macos_safe_filetypes([("Method JSON", "*.method.json")])
                )
                if sel:
                    print(f"[ML-PACK] method basename from user dialog: {Path(sel).stem}")
                    return Path(sel).stem
            except Exception:
                pass  # fall through to #5

        # 5) last resort
        print(f"[ML-PACK] method basename fallback to CSV stem: {Path(predict_input_path).stem}")
        return Path(predict_input_path).stem
    # ===

    # 20260418 code review: possible for CLI implementation in future. Built 20250915. Called by create_unlabeled_from_method
    def _pick_file_cli_or_gui(title="Select a file", patterns=(("CSV", "*.csv"), ("All files", "*.*"))):
        # Try GUI first
        try:
            from tkinter import filedialog
            path = filedialog.askopenfilename(title=title, filetypes=patterns)
            if path:
                return path
        except Exception:
            pass
        # Fallback to CLI prompt
        try:
            print(f"{title}: enter full path (or leave blank to cancel)")
            path = input("> ").strip()
            return path or None
        except Exception:
            return None    


    # read ML parameters written to json file. Field: ml_defaults
    def _read_ml_defaults_from_exp(path):
        try:
            with open(path, "r", encoding="utf-8") as f: data = json.load(f)
            return data.get("ml_defaults")
        except Exception:
            return None
    # write ML parameters to json file. Field: ml_defaults
    def _write_ml_defaults_to_exp(path, params: dict) -> bool:
        try:
            data = {}
            if os.path.exists(path):
                with open(path, "r", encoding="utf-8") as f: data = json.load(f)
            data["ml_defaults"] = params
            with open(path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False); f.write("\n")
            return True
        except Exception as e:
            messagebox.showerror("Save failed", str(e)); return False

    # ========== 20260419 start from here

    # 20260418 code review: reliable helper debug printer to track real values for training
    def _debug_print_train_params():                 
        print(f"[ML-PARAMS] n_estimators={n_estimators_var.get()}, test={test_split_var.get()}, val={val_split_var.get()},", end = "")                                    
        print(f"min_samples={min_samples_var.get()},balance={use_balance_var.get()}, majority={majority_label_var.get()}x{majority_factor_var.get()},",end ="")
        print(f"stratify={use_stratify_var.get()}, real_world={real_world_test_var.get()}, class_weight={class_weight_var.get()}, thresh_on={enable_thresh_var.get()},tau={thresh_val_var.get()}, margin={margin_val_var.get()}")

    # ---------------------------------------------------------

    def select_train_csv():
        nonlocal train_csv_path, linked_exp_json
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if not path:
            return
        train_csv_path = os.path.abspath(path)

        exp_json_path = None
        base_dir = os.path.dirname(path)
        for fname in os.listdir(base_dir):
            if fname.endswith(".exp.json"):
                exp_json_path = os.path.join(base_dir, fname)
                break
        linked_exp_json = exp_json_path

        if exp_json_path and os.path.exists(exp_json_path):
            try:
                with open(exp_json_path, "r") as f:
                    exp_data = json.load(f)
                exp_data["train_csv"] = train_csv_path
                with open(exp_json_path, "w") as f:
                    json.dump(exp_data, f, indent=4)
            except Exception as e:
                print(f"[ML] Failed to update experiment JSON: {e}")

        origin_info.configure(state="normal")
        origin_info.delete(1.0, "end")
        origin_info.insert("end", f"Loaded file: {os.path.basename(path)}\n")
        if exp_json_path:
            origin_info.insert("end", f"Linked .exp.json: {os.path.basename(exp_json_path)}\n")
        origin_info.insert("end", f"Path: {path}")
        origin_info.configure(state="disabled")

     # ---- ML state (analysis window scope) ----
    def _read_ml_from_json(path):
        if not path or not os.path.exists(path):
            return None
        try:
            with open(path, "r", encoding="utf-8") as f:
                d = json.load(f) or {}
            ml = d.get("ml") or {}
            return ml.get("parameters") or ml
        except Exception:
            return None

    # default model + train (used when nothing else provided)
    # Removed BUILTIN_ML(4)
    
    def snapshot_train_vars() -> dict:
        """Live values from the Train/Test dialog & quick RF options."""
        return {
            "model": {
                "type": "RandomForest",
                "n_estimators": int(n_estimators_var.get()),
                "max_depth": None,
                "min_samples_split": 2,
                "min_samples_leaf": 1,
                "random_state": 42,
                "class_weight": ("balanced" if class_weight_var.get() else None),
            },
            "train": {
                "test_size": float(test_split_var.get()),
                "val_size": float(val_split_var.get()),
                "stratify": bool(use_stratify_var.get()),
                "real_world_test": bool(real_world_test_var.get()),
                "threshold": {
                    "enabled": bool(enable_thresh_var.get()),
                    "tau": float(thresh_val_var.get()),
                    "margin": float(margin_val_var.get()),
                },
                # Optional: you can store these for provenance
                "balance": {
                    "enabled": bool(use_balance_var.get()),
                    "majority_label": majority_label_var.get(),
                    "majority_factor": int(majority_factor_var.get()),
                },
            },
        }
    # Show live updates on Train tab Step 4
    def _update_ml_summary(): #effective_ml_params or BUILTIN_ML -> effective_ml_params or snapshot_train_vars()
        try:
            src = effective_ml_params or snapshot_train_vars()
            m = src.get("model", {})
            t = src.get("train", {})
            b = t.get("balance", {})
            # 20260518 B-09: read live slider tk-var (was schema-mismatch top-level lookup) + add separator
            ml_summary_var.set(
            f'Params: RF Trees: {m.get("n_estimators")}, '
            f'Split test={t.get("test_size")}, val={t.get("val_size")}, '
            f'min_class_size={int(min_samples_var.get())}, '
            f'real_world={t.get("real_world_test")}, '
            f'balance={b.get("enabled")}'
            )
        except Exception:
            ml_summary_var.set("Params: (using built-ins)") 

    def open_ml_params_window():
        win = tk.Toplevel(root)
        win.title("ML Parameters")
        win.geometry("1200x600")
        win.transient(root)
        win.grab_set()

        # --- left editor ---
        left = ttk.LabelFrame(win, text="Editor (JSON)")
        left.pack(side="left", fill="both", expand=True, padx=(10,5), pady=10)
        editor_txt = tk.Text(left, wrap="none", height=20, width=45)
        editor_txt.pack(fill="both", expand=True, padx=8, pady=8)

        seed = _ensure_ml_state().get("editor") or snapshot_train_vars() # Replace BUILTIN_ML to snapshot_train_vars()
        editor_txt.insert("1.0", json.dumps(seed, indent=2))

        # --- right preview (effective after merge) ---
        right = ttk.LabelFrame(win, text="Effective (Builtins ← Files ← Editor)")
        right.pack(side="left", fill="both", expand=True, padx=(5,10), pady=10)
        preview_txt = tk.Text(right, wrap="none", height=20, width=45, state="disabled")
        preview_txt.pack(fill="both", expand=True, padx=8, pady=8)
        debug_var = tk.StringVar(value="")
        ttk.Label(right, textvariable=debug_var).pack(anchor="w", padx=8, pady=(0,6))

        st = _ensure_ml_state()
        st["editor_txt"] = editor_txt
        st["preview_txt"] = preview_txt
        _emit_ml_state_changed_ui_refresh()


        # --- local helpers that use the widgets above ---
        def _get_editor_json_or_empty():
            try:
                return json.loads(editor_txt.get("1.0", "end").strip() or "{}")
            except Exception:
                return {}
            
        def _get_ml_context():
            # If you already wrote a context getter elsewhere, you can call it here.
            return {"exp_json": linked_exp_json, "method_json": None}

        def _refresh_effective():
            # ← the key: include snapshot_train_vars() in the merge
            ctx = _get_ml_context()
            exp_defs  = _read_ml_from_json(ctx.get("exp_json"))
            meth_defs = _read_ml_from_json(ctx.get("method_json"))
            editor    = _get_editor_json_or_empty()
            eff = _merge_ml(snapshot_train_vars(),exp_defs,meth_defs,editor)  #removed BUILTIN_ML
            # editor first…snapshot_train_vars() # …then UI snapshot overrides
            #eff = _merge_ml(BUILTIN_ML, exp_defs, meth_defs, snapshot_train_vars(), editor)
            #test, remove if works
            try:
                ts = float(test_split_var.get())
                assert abs((eff.get("train",{}).get("test_size", -1)) - ts) < 1e-9, "merge precedence wrong"
            except Exception:
                pass
            #
            preview_txt.config(state="normal")
            preview_txt.delete("1.0", "end")
            preview_txt.insert("1.0", json.dumps(eff, indent=2))
            preview_txt.config(state="disabled")
            try:
                debug_var.set(
                    f"live: test={test_split_var.get():.2f}, "
                    f"val={val_split_var.get():.2f}, trees={n_estimators_var.get()}"
                )
            except Exception:
                pass

        ## adapted from build_ml_params_panel
        def _get_edit():
            try:
                txt = editor_txt.get("1.0","end").strip() or "{}"
                return json.loads(txt)
            except Exception as e:
                messagebox.showerror("Invalid JSON", f"Editor JSON parse error:\n{e}")
                return None
            
        # copied from build_ml_params_panel
        def _pull_from_exp():
            ctx = _get_ml_context()
            exp_path = ctx.get("exp_json")
            if not exp_path:
                messagebox.showwarning("No .exp.json found", "Select a training CSV that has an .exp.json in the same folder, or use 'Load preset' instead."); return
            defs = _read_ml_defaults_from_exp(exp_path)
            if defs is None:
                messagebox.showinfo("No defaults", "This experiment has no ml_defaults yet.")
                return
            editor_txt.delete("1.0","end"); editor_txt.insert("1.0", json.dumps(defs, indent=2, ensure_ascii=False))
            _refresh_effective()

        def _apply_to_exp():
            ctx = _get_ml_context()
            exp_path = ctx.get("exp_json")
            if not exp_path:
                messagebox.showwarning("No .exp.json found", "Select a training CSV that has an .exp.json in the same folder, or use 'Load preset' instead."); return
            # 20260517 fix B-03: use _get_edit() (validating helper) so editor JSON is parsed before being passed
            # to _write_ml_defaults_to_exp (which expects a dict, not a raw string). Mirrors _save_preset_as.
            d = _get_edit()
            if d is None: return
            if _write_ml_defaults_to_exp(exp_path, d):
                messagebox.showinfo("Saved", f"Updated ml_defaults in:\n{exp_path}")
            _refresh_effective()

        def _load_preset():
            path = filedialog.askopenfilename(title="Load params JSON", filetypes=[("JSON","*.json"), ("All files","*.*")])
            if not path: return
            try:
                with open(path, "r", encoding="utf-8") as f: d = json.load(f)
            except Exception as e:
                messagebox.showerror("Load failed", str(e)); return
            # accept whole-file payloads (exp files) or plain params
            d2 = d.get("ml_defaults", d)
            editor_txt.delete("1.0","end"); editor_txt.insert("1.0", json.dumps(d2, indent=2, ensure_ascii=False))
            _refresh_effective()

        def _save_preset_as():
            path = filedialog.asksaveasfilename(title="Save params JSON", defaultextension=".json",
                                                filetypes=[("JSON","*.json"), ("All files","*.*")])
            if not path: return
            d = _get_edit()
            if d is None: return
            try:
                with open(path, "w", encoding="utf-8") as f:
                    json.dump(d, f, indent=2, ensure_ascii=False); f.write("\n")
                messagebox.showinfo("Saved", f"Saved preset to:\n{path}")
            except Exception as e:
                messagebox.showerror("Save failed", str(e))

        def _link_exp_json():                                     
            nonlocal linked_exp_json                                                         
            # 20260520 fix B-30: macOS-safe filetypes wrap (compound .exp.json crashes macOS Tk)
            path = filedialog.askopenfilename(title="Select experiment JSON",
        filetypes=_macos_safe_filetypes([("Experiment JSON", "*.exp.json"), ("JSON", "*.json")]))
            if path:
                linked_exp_json = path                                                       
                messagebox.showinfo("Linked", f"Experiment linked:\n{os.path.basename(path)}") 
        # absorbed block ends, below are native buttons/function in open_ml_params_window
        # ===================        
        def _load_from_method():
            # prefer current method from context; otherwise ask
            path = _get_ml_context().get("method_json") or filedialog.askopenfilename(
                title="Choose method JSON", filetypes=[("JSON","*.json")]
            )
            if not path: return
            defs = _read_ml_from_json(path)
            if not defs:
                messagebox.showwarning("ML Parameters", "No ML block found in that JSON.")
                return
            editor_txt.delete("1.0", "end")
            editor_txt.insert("1.0", json.dumps(defs, indent=2))
            _refresh_effective()

        # 20260419 Code review refactored
        def _validate_and_use():
            nonlocal effective_ml_params
            try:
                ed = json.loads(editor_txt.get("1.0", "end"))
            except Exception as e:
                messagebox.showerror("Invalid JSON", str(e))
                return
            ctx = _get_ml_context()
            exp_defs  = _read_ml_from_json(ctx.get("exp_json"))
            meth_defs = _read_ml_from_json(ctx.get("method_json"))
            eff = _merge_ml(snapshot_train_vars(), exp_defs, meth_defs, ed) #_merge_ml(BUILTIN_ML, exp_defs, meth_defs, snapshot_train_vars(), editor)
            effective_ml_params = eff
            _ensure_ml_state()["editor"] = ed        # persist what’s in the left pane
            _ensure_ml_state()["effective_params"] = eff # optional: keep a copy
            _update_ml_summary()
            _refresh_effective()
            messagebox.showinfo("Ready", "Effective ML parameters set for training.")
            #win.destroy() # if users prefer auto-close, just uncomment this line

        # 20260419 Code review refactored
        def _save_effective_to_method():
            # pick current method or let user create one
            path = _get_ml_context().get("method_json") or filedialog.asksaveasfilename(
                title="Save or choose method JSON", defaultextension=".json",
                filetypes=[("JSON","*.json")]
            )
            if not path: return
            # ensure effective is up-to-date with what’s in the editor
            # 20260517 fix B-03 fold-B: use _get_edit() instead of _get_editor_json_or_empty() so invalid editor JSON
            # fail-closes (surfaces messagebox.showerror) instead of silently saving snapshot+exp+meth without the editor override. Mirrors B-03 main _apply_to_exp fix.
            editor = _get_edit()
            if editor is None: return
            ctx = _get_ml_context()
            exp_defs  = _read_ml_from_json(ctx.get("exp_json"))
            meth_defs = _read_ml_from_json(ctx.get("method_json"))
            eff = _merge_ml(snapshot_train_vars(), exp_defs, meth_defs, editor)
            # write back to method JSON
            os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
            data = {}
            if os.path.exists(path):
                try:
                    with open(path, "r", encoding="utf-8") as f:
                        data = json.loads(f.read().strip() or "{}")
                except Exception:
                    data = {}
            data.setdefault("ml", {})
            data["ml"]["updated"] = datetime.now().strftime("%Y-%m-%d %H:%M")
            data["ml"]["parameters"] = eff
            with open(path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)

            # reflect as current effective
            nonlocal effective_ml_params
            effective_ml_params = eff
            _update_ml_summary()
            _refresh_effective()
            messagebox.showinfo("Saved", f"ML parameters written to:\n{os.path.basename(path)}")

        def _on_params_close():
            st = _ensure_ml_state()
            st["editor_txt"]  = None
            st["preview_txt"] = None
            win.destroy()

        win.protocol("WM_DELETE_WINDOW", _on_params_close)
        # --- footer buttons ---
        btns_top = ttk.Frame(win)
        btns_top.pack(fill="x", padx=10, pady=(0,2))     
        btns_mid = ttk.Frame(win)
        btns_mid.pack(fill="x", padx=10, pady=(0,6))                                     
        btns_bot = ttk.Frame(win)
        btns_bot.pack(fill="x", padx=10, pady=(0,10))                                        
                                                                                            
        # top row: core actions                                                              
        ttk.Button(btns_mid, text="Validate & Use",                                          
        command=_validate_and_use).pack(side="left", padx=4)                                 
        ttk.Button(btns_top, text="Load from method…",            
        command=_load_from_method).pack(side="left", padx=4)                                 
        ttk.Button(btns_top, text="Save to method…",
        command=_save_effective_to_method).pack(side="left", padx=4)                         
        ttk.Button(btns_top, text="Close", command=_on_params_close).pack(side="right",
        padx=4)                                                                              
                                                                    
        # bottom row: presets + experiment                                                   
        ttk.Button(btns_mid, text="Load preset…", command=_load_preset).pack(side="left",
        padx=4)                                                                              
        ttk.Button(btns_mid, text="Save preset as…",              
        command=_save_preset_as).pack(side="left", padx=4)                                   
        ttk.Button(btns_bot, text="Pull from experiment",
        command=_pull_from_exp).pack(side="left", padx=4)                                    
        ttk.Button(btns_bot, text="Apply to experiment",          
        command=_apply_to_exp).pack(side="left", padx=4) 
        ttk.Button(btns_bot, text="Link exp manually",          
        command=_link_exp_json).pack(side="left", padx=4) 

        # live preview while typing in the left editor
        editor_txt.bind("<KeyRelease>", lambda _=None: _refresh_effective())
        # seed the right pane immediately on window open
        _refresh_effective()
    # =======================

    # ---------- UPDATED: parameters window ----------
    def open_train_settings():
        settings = tk.Toplevel()
        settings.title("Train/Test Parameters")
        settings.geometry("360x720")
        _ensure_ml_state()
        # splits
        tk.Label(settings, text="Test split").pack(pady=(10,0))
        tk.Scale(settings, from_=0.05, to=0.5, resolution=0.05, orient="horizontal",
                 variable=test_split_var).pack(fill="x", padx=10)

        tk.Label(settings, text="Validation split").pack(pady=(10,0))
        tk.Scale(settings, from_=0.00, to=0.3, resolution=0.05, orient="horizontal",
                 variable=val_split_var).pack(fill="x", padx=10)

        # --- Live split summary + note ---------------------------------- added 20250902
        note_lbl = tk.Label(settings, text="Note: Test + Validation must be < 1.0",
                            fg="#666666")
        note_lbl.pack(pady=(4,0))

        split_pct_lbl = tk.Label(settings, text="", font=("TkDefaultFont", 9, "bold"))
        split_pct_lbl.pack(pady=(2,6))


        def _push_vars_into_state(*_):
            # 20260418 _write_train_ui_into_state is removed
            if _settings_initializing:                                                       
                return  
            _emit_ml_state_changed_ui_refresh()  # now it updates directly from tk values
            _debug_print_train_params()

        def _on_settings_close():                                 
            for var, tid in trace_ids:                                                       
                try:                                              
                    var.trace_remove("write", tid)
                except Exception:                                                            
                    pass
            print(f"[DEBUG][WindowClose]", end="")
            _debug_print_train_params()
            settings.destroy()     

        # 2) wire all variables for live updates
        trace_ids = []
        for var in (
            test_split_var, val_split_var, use_stratify_var, real_world_test_var,
            enable_thresh_var, thresh_val_var, margin_val_var,
            min_samples_var, use_balance_var, majority_label_var, majority_factor_var,
            n_estimators_var, class_weight_var,
        ):
            tid = var.trace_add("write", _push_vars_into_state)                              
            trace_ids.append((var, tid))

        # Call _on_settings_close when window closes
        settings.protocol("WM_DELETE_WINDOW", _on_settings_close) 

        # 3) also bind Scales so dragging updates continuously
        def _scale_cb(_val):
            _push_vars_into_state()

        for sc in (test_split_var, val_split_var):  # use your actual Scale objects
            try:
                sc.configure(command=_scale_cb)
            except Exception:
                pass

        def _update_split_labels(*_):
            try:
                if not split_pct_lbl.winfo_exists():
                    return
            except Exception:
                return
            ts = float(test_split_var.get())
            vs = float(val_split_var.get())
            tr = 1.0 - (ts + vs)
            split_pct_lbl.configure(
                text=f"Train {max(tr,0)*100:.0f}%  |  Val {vs*100:.0f}%  |  Test {ts*100:.0f}%"
            )
            # color rules
            if ts < 0 or vs < 0 or (ts + vs) >= 1.0:
                split_pct_lbl.configure(fg="red")            # invalid
            elif tr < 0.10:
                split_pct_lbl.configure(fg="#D17D00")        # warning: tiny train set
            else:
                split_pct_lbl.configure(fg="black")          # ok

        # trigger on slider move + initialize
        test_split_var.trace_add("write", _update_split_labels)
        val_split_var.trace_add("write", _update_split_labels)
        _update_split_labels()
        # ---------------------------------------------------------------
        _settings_initializing = True #prevent debug print fires during button loading
        # Min class size (drop below this value)
        tk.Label(settings, text="Minimum samples per class").pack(pady=(12,0))
        tk.Spinbox(settings, from_=1, to=50, textvariable=min_samples_var, width=6).pack()
        # balancing block
        ttk.Checkbutton(settings, text="Enable class balancing (cap majority)", variable=use_balance_var).pack(pady=(12,0)) #to ttk
        row = tk.Frame(settings); row.pack(pady=2)
        tk.Label(row, text="Majority label:").pack(side="left")
        tk.Entry(row, textvariable=majority_label_var, width=16).pack(side="left", padx=6)
        row2 = tk.Frame(settings); row2.pack(pady=2)
        tk.Label(row2, text="Majority factor (×max minor):").pack(side="left")
        tk.Spinbox(row2, from_=1, to=20, textvariable=majority_factor_var, width=6).pack(side="left", padx=6)
        # stratify
        ttk.Checkbutton(settings, text="Stratify by label", variable=use_stratify_var).pack(pady=(12,0)) #tk to ttk
        # --- Random Forest options ---
        sep = ttk.Separator(settings, orient="horizontal"); sep.pack(fill="x", padx=10, pady=(12,6))
        tk.Label(settings, text="Random Forest Options").pack()

        row_rf = tk.Frame(settings); row_rf.pack(pady=2)
        tk.Label(row_rf, text="n_estimators (trees):").pack(side="left")
        tk.Spinbox(row_rf, from_=50, to=2000, increment=50, textvariable=n_estimators_var, width=7).pack(side="left", padx=6)

        ttk.Checkbutton(settings, text='Use class_weight = "balanced"', variable=class_weight_var).pack(pady=2)#tk to ttk
        ttk.Checkbutton(settings, text="Real-world test (cap training only)", variable=real_world_test_var).pack(pady=(4,0)) #tk to ttk

        sep3 = ttk.Separator(settings, orient="horizontal"); sep3.pack(fill="x", padx=10, pady=(12,6))
        tk.Label(settings, text="Prediction Thresholding").pack()

        ttk.Checkbutton(
            settings,
            text="Enable confidence threshold → fallback to Majority label",
            variable=enable_thresh_var).pack(pady=(2,2))
        
        row_thr = tk.Frame(settings); row_thr.pack(pady=2)
        tk.Label(row_thr, text="Threshold τ (0.00–0.99):").pack(side="left")
        tk.Spinbox(row_thr, from_=0.00, to=0.99, increment=0.01,
                textvariable=thresh_val_var, width=6).pack(side="left", padx=6)

        row_margin = tk.Frame(settings); row_margin.pack(pady=2)
        tk.Label(row_margin, text="Majority support margin δ (0.00–0.20):").pack(side="left")
        tk.Spinbox(row_margin, from_=0.00, to=0.20, increment=0.01,
                textvariable=margin_val_var, width=6).pack(side="left", padx=6)
        _settings_initializing = False
        # Extra json-like editor panel button
        tk.Button(settings, text="Edit ML Parameters…",
                command=open_ml_params_window).pack(pady=(8, 0))
        # Force update, should be useless after code review in 20260419 (validated everything works)
        tk.Button(settings, text="Force current params", command=_push_vars_into_state).pack(pady=(8,0))

        # push once so editors / preview reflect these values
        _on_change()
        # 20260419 code review: once user click the button to open train settings, show params regardless of user edit or not. Unsealed status notification.
        _update_ml_summary()

        tk.Button(settings, text="Close", command=_on_settings_close).pack(pady=8)
    # ---------- UPDATED: parameters window ----------

    def train_model(): #RF

        try:
            from sklearn.ensemble import RandomForestClassifier
            from sklearn.metrics import classification_report
            from sklearn.metrics import confusion_matrix
            from sklearn.preprocessing import LabelEncoder
            from sklearn.utils.multiclass import unique_labels
            import joblib
            import pretrain_normalizer as normalizer
            from sklearn.metrics import precision_recall_fscore_support, accuracy_score
        except ImportError:
            messagebox.showerror("Missing Dependencies",
                                 "Please install scikit-learn and joblib.")
            return
        #20251005
        use_mass_feature = bool(include_mass_train_var.get())

        # 20260419 code review: if user clicks Validate&Use, effective_ml_params gets not {}, cfg capture effective live params, otherwise call from tk values
        if effective_ml_params:
            cfg = effective_ml_params
            print(f"[ML-DEBUG] cfg from: effective_ml_params,n_estimators={cfg.get('model',{}).get('n_estimators')}")                             
        else:
            cfg = snapshot_train_vars() 
            print(f"[ML-DEBUG] cfg from: snapshot_train_vars(), n_estimators={cfg.get('model',{}).get('n_estimators')}")
        rf = cfg.get("model", {})
        #tr = cfg.get("train", {})  # Need to think about the params from advanced panel, current split doens't catch that tho
        # 20260418 Code review Debug prints
        print(f"[ML-DEBUG] cfg source: {'effective_ml_params' if effective_ml_params else 'BUILTIN_ML+snapshot'}")                                                             
        print(f"[ML-DEBUG] cfg keys: {list(cfg.keys())}")
        print(f"[ML-DEBUG] cfg.model: {cfg.get('model', {})}")                               
        print(f"[ML-DEBUG] cfg.train: {cfg.get('train', {})}")  

        # check if trainable csv to train is there
        if not train_csv_path:
            messagebox.showwarning("No File", "Please select a training CSV first.")
            return
        else:
            logger.log(f"[ML][train][info] Using {train_csv_path} to train model")

        label_col = label_dropdown.get().strip()
        if label_col not in ['Structure', 'IUPACname(optional)', 'Glycanannotation2', 'GlyToucan ID']:
            messagebox.showerror("Invalid Label", f"'{label_col}' is not a supported label column.")
            return

        # load
        try:
            df = pd.read_csv(train_csv_path)
        except Exception as e:
            logger.log(f"[ML][train][ERROR 1] Trainable csv can't be read:{e}")
            messagebox.showerror("Read Error", str(e)); return
        # if trainable csv don't have the selected label for training
        if label_col not in df.columns:
            logger.log(f"[ML][train][ERROR 1] Label column: {label_col} not in trainable csv")
            messagebox.showerror("[ML][train][ERROR 1]Missing Column", f"{label_col} not in trainable csv"); return

        #20250929 ver to accept SO3(s), HexA and PO3H(p)
        # convert the composition to manual-style ONLY when needed
        print("[ML][train][DEBUG] labels before:", sorted(set(df[label_col].astype(str)))[:10], "...")

        # Normalize ONLY when needed (do not strip A / s / p!), need to take a review in future at normalizer
        try:
            if label_col == "Structure":
                before = sorted(set(df[label_col].astype(str)))[:10]
                df[label_col] = df[label_col].map(
                    lambda x: normalize_structure_for_training(x, normalizer.parse_structure_to_manual)
                )
                print("[ML][train][DEBUG][\"Structure\"] labels preview before normalization:", before, "...")
            # For GlyToucan/IUPAC, leave as-is. 20260419 code review: add string conversion so at least in can be fed into RF trainer
            else:
                before = df[label_col][:10]
                df[label_col] = df[label_col].astype(str).str.strip()
                # 20260517 fix B-17: NaN-aware empty-label guard. pd.read_csv reads blanks as NaN; astype(str) turns them into the literal "nan", so a naive == "" check misses real empties. Cover the common pandas/CSV empty representations and warn-and-confirm before silently training on phantom-class rows (Codex Edit 3 flag).
                _b17_empty_mask = df[label_col].str.strip().str.lower().isin({"", "nan", "none", "null"})
                _b17_empty_count = int(_b17_empty_mask.sum())
                if _b17_empty_count > 0:
                    if not messagebox.askyesno(
                        "Empty Labels",
                        f"Label column '{label_col}' has {_b17_empty_count} of {len(df)} rows with empty values. "
                        f"These rows will be dropped before training. Continue?",
                    ):
                        logger.log(f"[ML][train] User cancelled training: {_b17_empty_count} empty {label_col} rows.")
                        return
                    df = df.loc[~_b17_empty_mask].copy()
                    if df.empty:
                        logger.log(f"[ML][train][ERROR 1] All rows dropped — no usable {label_col} labels remain.")
                        messagebox.showerror("No Labels", f"All rows had empty {label_col}. Nothing to train.")
                        return
                print("[ML][train][DEBUG] Other labels that will be turned into string for training:", before, "...")
        except Exception as e:
            logger.log(f"[ML][train][ERROR 1] structure normalization failed {e}")
            print("[ML][train][ERROR 1] structure normalization failed due to:", e)
            messagebox.showerror(
                "Structure normalization failed. See console/log for details.")
            return
        
        print("[ML][train][DEBUG] labels preview after normalization:",  sorted(set(df[label_col].astype(str)))[:10], "...")

        # encode label AFTER any string/tuple harmonization (if needed)
        y_labels = df[label_col].astype(str) # y_raw -> y_labels

        # run balancing + split, cfg getting bypassed if someone edits json-like with n_estimators, other values there may not been reflected (need investigation)
        # consider adding debug print in future tracing this potential design/logic weakness
        try:
            if use_balance_var.get():
                X_train, y_train, X_val, y_val, X_test, y_test, info = balance_and_split(
                    df.assign(**{label_col: y_labels}),
                    label_col=label_col,
                    majority_label=majority_label_var.get(),
                    min_count=int(min_samples_var.get()),
                    majority_factor=int(majority_factor_var.get()),
                    test_size=float(test_split_var.get()),
                    val_size=float(val_split_var.get()),
                    stratify=bool(use_stratify_var.get()),
                    cap_training_only=bool(real_world_test_var.get()),   # << NEW
                    random_state=42
                )
            else:
                # no balancing — just use helper with factor=0 & no cap
                X_train, y_train, X_val, y_val, X_test, y_test, info = balance_and_split(
                    df.assign(**{label_col: y_labels}),
                    label_col=label_col,
                    majority_label="__no_cap__",  # won't match in balance_and_split against labels existing → no cap
                    min_count=int(min_samples_var.get()),
                    majority_factor=1,
                    test_size=float(test_split_var.get()),
                    val_size=float(val_split_var.get()),
                    stratify=bool(use_stratify_var.get()),
                    cap_training_only=bool(real_world_test_var.get()),   # << NEW
                    random_state=42
                )
        except Exception as e:
            logger.log(f"[ML][train][ERROR 2] Split/Balancing Failed: {str(e)}")
            print(f"[ML][train][ERROR 2] Split/Balancing Failed: {str(e)}")
            messagebox.showerror("Split/Balancing Failed", str(e)); return

        #20251005
        # === MS1 feature handling (after split, before label encoding) ===
        # use_mass_feature already read above: use_mass_feature = bool(include_mass_train_var.get())

        def _prep_ms1_block(X, y, use_mass):

            X = X.copy()
            if use_mass:
                if "protonatedmass" in X.columns:
                    # make numeric and impute
                    X["protonatedmass"] = pd.to_numeric(X["protonatedmass"], errors="coerce")

                    # Non-glycan → 0.0 (use string labels here, before encoding)
                    if y is not None:
                        y_ser = pd.Series(list(y), index=X.index).astype(str)
                        X.loc[y_ser.eq("Non-glycan"), "protonatedmass"] = 0.0

                    # any remaining NaNs → 0.0 to satisfy RF
                    X["protonatedmass"] = X["protonatedmass"].fillna(0.0)
                else:
                    print("[Train] MS1 enabled but no 'protonatedmass' in X — continuing without it")
            else:
                # ensure we don't accidentally include it when toggle is off
                X = X.drop(columns=["protonatedmass"], errors="ignore")

            # absolute safety for other features too
            return X.fillna(0.0)

        X_train = _prep_ms1_block(X_train, y_train, use_mass_feature)
        X_val   = _prep_ms1_block(X_val,   y_val,   use_mass_feature)
        X_test  = _prep_ms1_block(X_test,  y_test,  use_mass_feature)
        # === end MS1 feature handling ===

        # LabelEncoder (consistent across splits)
        le = LabelEncoder()
        y_train_enc = le.fit_transform(y_train)
        y_val_enc   = le.transform(y_val) # 20260419 code review: in future curated libary having enough data, consider GridSearchCV or RandomizedSearchCV for looped optimization
        y_test_enc  = le.transform(y_test)


        # model (consume panel params if available; otherwise fallback to current UI/defaults)
        print(f"[ML][Train][finalcall] RF params: n_estimators={rf.get('n_estimators')},max_depth={rf.get('max_depth')}, class_weight={rf.get('class_weight')}")
        if not rf.get('n_estimators'):
            print("[ML][Train][WARN] n_estimators is empty/None — using hardcoded fallback (is cfg broken?)")
        model = RandomForestClassifier(
            n_estimators      = int(rf.get("n_estimators", 400)),  # for debug, change 400 to 666
            max_depth         = rf.get("max_depth", None),
            min_samples_split = int(rf.get("min_samples_split", 2)),
            min_samples_leaf  = int(rf.get("min_samples_leaf", 1)),
            random_state      = int(rf.get("random_state", 42)),
            class_weight      = rf.get("class_weight", ("balanced" if bool(class_weight_var.get()) else None)),
            n_jobs            = -1,
        )

        #20250929
        # --- NEW: sanitize features & lock column order ---
        X_train_num = _drop_meta_and_get_X(X_train)
        cols = list(X_train_num.columns)

        # align val/test to train's columns (fill missing with 0)
        X_val_num  = _drop_meta_and_get_X(X_val).reindex(columns=cols, fill_value=0)
        X_test_num = _drop_meta_and_get_X(X_test).reindex(columns=cols, fill_value=0)

        # replace originals, make sure only numeric columns enters the train process (or it will throw error)
        X_train, X_val, X_test = X_train_num, X_val_num, X_test_num

        # 20260419 code review: recovered the tags so now models show if it include protonatedmass as feature or not, not only showing the source csv status
        mass_suffix = "_withMass" if use_mass_feature else ""

        if X_train.isna().any().any():
            raise RuntimeError("[Train] NaNs remain in X_train after MS1 prep — unexpected")
        
        # Call RF to deal with numeric X feature matrices and Y (label mapped to integer)
        model.fit(X_train, y_train_enc)
        # Call RF to run prediction on test split, so we get performance evaluation later
        y_pred = model.predict(X_test)

        #extra binary glycan-vs-non-glycan report (print in terminal only)
        maj_idx = list(le.classes_).index(majority_label_var.get()) # find the major class index mapping integer
        y_true_bin = (y_test_enc != maj_idx).astype(int)   # Convert test labels -not majority(glycans) = 1, majority (Non-glycan) = 0
        y_pred_bin = (y_pred      != maj_idx).astype(int)  # Convert predetion from models - same logic

        p,r,f,_ = precision_recall_fscore_support(y_true_bin, y_pred_bin, average="binary", zero_division=0)
        acc = accuracy_score(y_true_bin, y_pred_bin)
        logger.log(f"[ML][train][report] Binary glycan vs non-glycan summary : Precison={p:.2f} Recall={r:.2f} F1 score={f:.2f} Accuracy={acc:.2f}")
        print(f"[ML][train][report] Binary glycan vs non-glycan summary : Precison={p:.2f} Recall={r:.2f} F1 score={f:.2f} Accuracy={acc:.2f}")

        # --- NEW: glycan-only confidence threshold -> fallback to majority label ---
        if enable_thresh_var.get():
            tau = float(thresh_val_var.get())
            #newly added
            #margin = 0.05  # δ: only demote if majority prob is within 0.05 of best predicted prob
            margin = float(margin_val_var.get())
            if hasattr(model, "predict_proba"):
                proba = model.predict_proba(X_test)      # [n_samples, n_classes]
                maxp = proba.max(axis=1)                 # best class prob per sample

                classes_ = list(le.classes_)
                try:
                    majority_idx = classes_.index(majority_label_var.get())
                except ValueError:
                    majority_idx = None

                if majority_idx is not None:
                    # predicted labels are encoded ints already
                    p_major = proba[:, majority_idx]
                    # Demote only when:
                    #   1) current pred is a glycan (not majority),
                    #   2) confidence below τ, and
                    #   3) majority prob is within δ of the best prob (model is "nearly indifferent")
                    demote = (y_pred != majority_idx) & (maxp < tau) & (p_major >= (maxp - margin))
                    if np.any(demote):
                        y_pred = y_pred.copy()
                        y_pred[demote] = majority_idx
                else:
                    print(f"[ML] Thresholding skipped: fallback label '{majority_label_var.get()}' not in training classes.")
            else:
                print("[ML] Thresholding skipped: classifier lacks predict_proba.")
            # --- end margin-aware demotion ---

        # --- END NEW ---
        
        # evaluate on test
        used = unique_labels(y_test_enc, y_pred)  # finds all class indices that appear in either the true labels OR the predictions
        used_names = [le.classes_[i] for i in used] #converts those integer indices back to human-readable labels using the encoder. So [0, 2, 5] becomes ["F1H5N4S2", "H5N4S2", "Non-glycan"].
        # quiet, deterministic handling of 0/0 cases 20250902
        report = classification_report(y_test_enc, y_pred, target_names=used_names, zero_division=0)
        #avoids cluttering the report with classes that have 0 samples in both true and predicted, which would show 0.00 across all metrics. 
        # (optional) flag classes with no predicted or no true samples
        labels_all = list(le.classes_)
        cm = confusion_matrix(y_test_enc, y_pred, labels=np.arange(len(labels_all)))
        no_pred = [labels_all[j] for j, s in enumerate(cm.sum(axis=0)) if s == 0]
        no_true = [labels_all[i] for i, s in enumerate(cm.sum(axis=1)) if s == 0]

        if no_pred: #never being assigned to these labels, since too similar to other labels (structures)
            logger.log(f"[ML][train][summary] No predicted samples for these classes: {no_pred}")
            print(f"[ML][train][summary] No predicted samples for these classes: {no_pred}")
        if no_true: #The class was in training but the random split put zero samples in the test set — happens with rare classes.
            logger.log(f"[ML][train][summary] (lost from splitting) No true samples for: {no_true}")
            print(f"[ML][train][summary] (lost from splitting) No true samples in test for: {no_true}")

        # Glycan only report by excluding the majority (Assume the majority label is Non-glycan)
        try:
            maj_idx = labels_all.index(majority_label_var.get())
        except ValueError:
            maj_idx = None

        if maj_idx is not None:
            mask_glycan_true = (y_test_enc != maj_idx)
            if np.any(mask_glycan_true):
                glycan_label_indices = [i for i, _ in enumerate(labels_all) if i != maj_idx]
                glycan_names = [labels_all[i] for i in glycan_label_indices]
                glycan_report = classification_report(
                    y_test_enc[mask_glycan_true],
                    y_pred[mask_glycan_true],
                    labels=glycan_label_indices,
                    target_names=glycan_names,
                    zero_division=0
                )
                print("\n[ML] Glycan-only report (excludes Non-glycan):\n", glycan_report)
        # --- end glycan-only report ---

        # notify
        # build a robust "cap" line that works for both modes
        cap_applied_to = info.get("majority_cap_applied_to") or info.get("train_majority_cap_applied_to")
        cap_value      = info.get("majority_cap") or info.get("train_majority_cap")
        cap_prefix     = "Capped training " if info.get("mode") == "cap_training_only" else "Capped "
        cap_line = (f"{cap_prefix}'{cap_applied_to}' at {cap_value}"
                    if cap_applied_to else "No majority cap applied.")

        msg = [
            "Random Forest trained successfully.",
            f"Mode: {'Cap training only' if real_world_test_var.get() else 'Cap before split'}",
            f"RF: {int(n_estimators_var.get())} trees, class_weight={'balanced' if class_weight_var.get() else 'none'}", 
            f"Thresholding: {'ON tau=' + format(thresh_val_var.get(), '.2f') + ', margin=' + format(margin_val_var.get(), '.2f') + ' -> ' + majority_label_var.get() if enable_thresh_var.get() else 'OFF'}",
            f"Classes kept: {len(info['kept_label_counts'])}",
            cap_line,
            f"Dropped tiny classes (< {min_samples_var.get()}):{sum(info['dropped_rare_counts'].values())}",
            "",                                                                                                                          
        ]
        messagebox.showinfo("Training Complete", "\n".join(msg))
        logger.log(f"[ML][Train][summary] {msg}")
        # save artifacts
        base = os.path.splitext(train_csv_path)[0] + mass_suffix # add mass_suffix, which is real value of if model get trained with mass or not
        # if user feels the double mass indicator annoying, or if we can wrap information into dataset explorer, we can make it tidied then.
        """
        for old_suffix in ("_withMass", "_noMass"):
            if base.endswith(old_suffix):
                base = base[:-len(old_suffix)]
                break
        base = base + mass_suffix
        """
        model_path  = base + "_rf_model.joblib"
        enc_path    = base + "_labelencoder.joblib"
        report_path = base + "_rf_performance.txt"
        features_path = base + "_features.json"    
        # 20260418 skops first, joblib as fallback, no more same model in 2 different extension
        try:                                                                                 
            from skops.io import dump as sk_dump
            sk_dump(model, base + "_rf_model.skops")                                         
            logger.log(f"[ML][Train] Model saved (skops): {base}_rf_model.skops")
        except ImportError:                                                                  
            joblib.dump(model, model_path)         
            logger.log(f"[ML][Train][fallback] Model saved (joblib): {base}_rf_model.joblib")  

        joblib.dump(le, enc_path)
        logger.log(f"[ML][Train] Model label encoder saved (joblib): {base}_labelencoder.joblib")  

        with open(report_path, "w", encoding="utf-8") as f:
            f.write(report)
            logger.log(f"[ML][Train] Model performance report saved (txt): {base}_rf_performance.txt")  
            if maj_idx is not None:
                f.write(glycan_report)
                logger.log(f"[ML][Train] Glycan-only (excluding the majority) report is appended.")  
            print(f"[ML][Train] Model performance report saved (txt): {base}_rf_performance.txt") 

        with open(features_path, "w", encoding="utf-8") as f:                                
            json.dump(list(X_train.columns), f, indent=2)                                    
        logger.log(f"[ML] Feature list saved: {features_path}")    

        # 20260418 code review to activate train environment json save
        try:
            snapshot_training_run(                                                           
                artifact_dir=os.path.dirname(base),           
                effective_params=cfg,                                                        
                column_order=list(X_train.columns),
                classes=list(le.classes_),                                                   
                inputs={"train_csv": train_csv_path, "linked_exp": linked_exp_json},
                hashes={"train_csv": file_sha256(train_csv_path) if train_csv_path else      
        None},                                                                               
                versions=collect_versions(),                                                 
            )         
            logger.log(f"[ML][Train] records save at {base} ")                                                                       
        except Exception as e:                                    
            logger.log(f"[ML] snapshot_training_run failed: {e}")
        # 20260419 code review: note that exp json may not be in the scope (nothing to do with current workflow)
        # persist params back to exp.json if present
        if linked_exp_json and os.path.exists(linked_exp_json):
            try:
                with open(linked_exp_json, "r", encoding="utf-8") as f:
                    exp_data = json.load(f)
                exp_data["train_parameters"] = {
                    "split_ratio_test": float(test_split_var.get()),
                    "split_ratio_val": float(val_split_var.get()),
                    "min_samples": int(min_samples_var.get()),
                    "balancing_enabled": bool(use_balance_var.get()),
                    "majority_label": majority_label_var.get(),
                    "majority_factor": int(majority_factor_var.get()),
                    "stratify": bool(use_stratify_var.get()),
                    "rf_n_estimators": int(n_estimators_var.get()),
                    "rf_class_weight": "balanced" if class_weight_var.get() else "none",
                    "real_world_test_cap_training_only": bool(real_world_test_var.get()),
                    "threshold_enabled": bool(enable_thresh_var.get()),
                    "threshold_tau": float(thresh_val_var.get()),
                    "threshold_margin": float(margin_val_var.get()),
                    "model_path": model_path,
                    "labelencoder_path": enc_path,
                    "report_path": report_path
                }
                with open(linked_exp_json, "w", encoding="utf-8") as f:
                    json.dump(exp_data, f, indent=4)
            except Exception as e:
                print("Failed to write training parameters to exp.json:", e)
    # end of RF model train section


    def select_model_file():
        nonlocal model_file_path
        path = filedialog.askopenfilename(title="Select Model File",filetypes=[("Model files", "*.joblib *.pkl *.skops"), ("All files", "*.*")])
        if path:
            model_file_path = os.path.abspath(path)
            messagebox.showinfo("Model Loaded", f"Model loaded from:\n{model_file_path}")

    def select_predict_input():
        nonlocal predict_input_path
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            predict_input_path = os.path.abspath(path)
            messagebox.showinfo("Prediction file selected", f"Data loaded from:\n{predict_input_path}")

    def run_prediction():

        try:
            import joblib
        except ImportError:
            messagebox.showerror("Missing Dependency", "joblib is not installed.")
            return

        if not model_file_path or not predict_input_path:
            print(f"[ML][Predict]: missing either model file {model_file_path} or prediction table {predict_input_path}")
            messagebox.showwarning("Missing Info", "Please select both a model file and an input CSV file.")
            return

        #20250909 add import for prediction report
        from prediction_report import (ReportParams, summarize_predictions,write_prediction_report, show_prediction_summary_popup)
        # 20260421 am 3:17 fixed the skops error (label encoder is not correctly parsed. Added loader condition to get correct le_path)
        try:
            model, loader = load_model_any(model_file_path)
            df = pd.read_csv(predict_input_path)
            logger.log(f"[ML][Predict] Model loaded via {loader}")
            # 20260521 fix B-33 (Branch 2 — tolerant cross-version notice, Codex-approved):
            # Surface a non-blocking info popup + logger line when the model loaded
            # successfully BUT trained sklearn version differs from current. Prediction
            # proceeds; user is informed of the version delta in case results look off.
            trained_sklearn = _read_trained_sklearn(model_file_path)
            current_sklearn = _current_sklearn_version()
            if trained_sklearn and trained_sklearn != current_sklearn:
                logger.log(f"[ML][Predict][version-notice] trained=sklearn-{trained_sklearn} current=sklearn-{current_sklearn} — load OK, proceeding")
                messagebox.showinfo(
                    "sklearn version notice",
                    f"Model trained with sklearn {trained_sklearn}; this machine has sklearn {current_sklearn}.\n\n"
                    f"Loaded successfully — prediction will proceed.\n"
                    f"If results look off, retrain on this machine for exact-version reproducibility."
                )
            # Try loading label encoder if available
            if loader == "joblib":
                logger.log(f"[ML][Predict] loading label encoder (joblib model)")
                try:
                    le_path = model_file_path.replace("_rf_model.joblib", "_labelencoder.joblib")
                    if os.path.exists(le_path):
                        le = joblib.load(le_path)
                    else:
                        le = None
                except:
                    logger.log(f"label encoder can't be loaded")
            elif loader == "skops":
                logger.log(f"[ML][Predict] loading label encoder (skops model)")
                try:
                    le_path = model_file_path.replace("_rf_model.skops", "_labelencoder.joblib")
                    if os.path.exists(le_path):
                        le = joblib.load(le_path)
                    else:
                        le = None
                except:
                    logger.log(f"label encoder can't be loaded")
        except Exception as e:
            logger.log(f"[ML][Predict][ERROR 0] exception in loading model {e}")
            # 20260521 fix B-33 (Branch 3 — version-aware guard rewrite, Codex-approved):
            # Replaces v14:7838-7853 prior guard. Three fixes vs prior code:
            # (A) operator precedence bug — prior `A and B or C` parsed as
            #     `(A and B) or C` so .skops dtype errors fell into the joblib branch.
            # (B) wrong message copy for .skops + dtype — prior message named
            #     "joblib model" and advised "re-export as .skops" which is wrong
            #     when the file already IS .skops and the dtype change is in sklearn
            #     internals (not the persistence format).
            # (C) unreachable .skops elif — couldn't fire for the dtype case because (A).
            # New: read training_run.json sidecar for trained sklearn version; compare
            # to current; surface BOTH versions in the message + logger line.
            err_text = str(e)
            is_version_error = "dtype" in err_text or "InconsistentVersionWarning" in err_text
            trained_sklearn = _read_trained_sklearn(model_file_path)
            current_sklearn = _current_sklearn_version()
            versions_differ = bool(trained_sklearn) and trained_sklearn != current_sklearn

            if is_version_error and versions_differ:
                logger.log(f"[ML][Predict][version-mismatch-fatal] trained=sklearn-{trained_sklearn} current=sklearn-{current_sklearn} err={err_text[:120]}")
                suffix = ".skops" if ".skops" in model_file_path else ".joblib"
                messagebox.showerror(
                    f"Model version mismatch ({suffix})",
                    f"This model cannot load on this machine.\n\n"
                    f"Trained with: sklearn {trained_sklearn}\n"
                    f"This machine: sklearn {current_sklearn}\n\n"
                    f"Fix: retrain on this machine.\n\n"
                    f"(Background: scikit-learn does not guarantee cross-version model loading. "
                    f"For RandomForest/DecisionTree models, sklearn 1.3 introduced tree-internal "
                    f"changes such as `missing_go_to_left`. Other estimator types and sklearn "
                    f"versions can have similar break points; retraining on this machine is the safe fix.)"
                )
            elif is_version_error:
                logger.log(f"[ML][Predict][load-fail-dtype-versions-match] trained=sklearn-{trained_sklearn or 'unknown'} current=sklearn-{current_sklearn} err={err_text[:120]}")
                messagebox.showerror(
                    "Model load failed",
                    f"Dtype/pickle error during load (sklearn versions appear to match: "
                    f"trained={trained_sklearn or 'unknown'}, current={current_sklearn}).\n\nReason:\n{e}"
                )
            elif ".skops" in model_file_path:
                logger.log(f"[ML][Predict][load-fail-skops] err={err_text[:120]}")
                messagebox.showerror("Model load failed", f"Failed to load .skops model.\n\nReason:\n{e}")
            else:
                logger.log(f"[ML][Predict][load-fail-other] err={err_text[:120]}")
                messagebox.showerror("Load Error", str(e))
            return
        # 20260419 code review: Stage 2 Feature Alignment
        #20251006 ver
        try:
            # 1) Get the model's training feature list first
            train_feats = get_training_features(model, model_file_path)
            if train_feats is None:
                logger.log(f"[ML][Predict][ERROR 0] Failed to determine training feature list. If you have *_features.json for this file, put in the same folder.")
                raise RuntimeError("Cannot determine training feature list. "
                                "Re-export with embedded features or provide *_features.json.")

            # 2) If the model expects MS1 but input lacks it, add placeholder NOW (before X)
            if "protonatedmass" in train_feats and "protonatedmass" not in df.columns:
                logger.log(f"[ML][Predict][WARN] This model supports mass feature, input is missing mass column so it will be set to 0 (lower prediction performance expected)")
                messagebox.showwarning(
                    "Model expects MS1",
                    "This model was trained with 'protonatedmass', but the input has no such column. "
                    "Proceeding with a zero-filled placeholder."
                )
                df["protonatedmass"] = 0.0  # float OK

            # (optional) if it expects delta_ppm too
            if "delta_ppm" in train_feats and "delta_ppm" not in df.columns:
                logger.log(f"[ML][Predict][WARN] This model supports delta_ppm, but it is missing in input file")
                df["delta_ppm"] = 0.0

            # 3) Now build X (after placeholders exist)
            X = df.drop(columns=["MS2scan_no"], errors="ignore")

            # 4) Align to training features
            missing = sorted(set(train_feats) - set(X.columns))
            if missing:
                logger.log(f"[ML][Predict][WARN] {len(missing)} features missing from input, filled with 0: {missing[:10]}{'...' if len(missing) > 10 else ''}")
            extra   = sorted(set(X.columns) - set(train_feats))
            if extra:
                logger.log(f"[ML][Predict][info] {len(extra)} extra columns in input (ignored): {extra[:10]}{'...' if len(extra) > 10 else ''}")

            for c in missing:
                X[c] = 0.0
            X = X[list(train_feats)]  # enforce order

            # 20260419 code review: Stage 3 predict (only this line)
            # Predict row by row, the order is not changed so we can make sure the predictions can be filled correctly to original MS2_scan or whatever
            y_pred = model.predict(X)

            # init proba for later demotion and reporter
            proba = None
            if hasattr(model, "predict_proba"):
                proba = model.predict_proba(X) #(keep this BEFORE gating and summarization)
            # Optional probabilities & margin
            #applying same filter to prediction model
            # --- optional: apply the same margin-aware demotion in prediction ---
            # 20260419 code review: tau and margin was applied always until 20260419
            # now it is bound with parameters in train/test parameter, need to add a button back here locally (probably changing the value itself and sync would be fine though)
            if enable_thresh_var.get():
                try:
                    tau = float(thresh_val_var.get())
                except Exception:
                    tau = 0.60  # sensible fallback
                #similar to tau, use same pattern on margin
                try:
                    margin = float(margin_val_var.get())
                except Exception:
                    margin = 0.05
                logger.log(f"[ML][Predict] Enable thresholding, tau = {tau}, margin = {margin}")

                # work in encoded-space (ints). If y_pred are strings because of a prior transform, re-encode temporarily.
                if le is not None and (len(y_pred) > 0 and isinstance(y_pred[0], str)):
                    y_pred_enc = le.transform(y_pred)
                    classes_ = list(le.classes_)
                else:
                    y_pred_enc = y_pred
                    classes_ = list(model.classes_)
                try:
                    majority_idx = classes_.index(majority_label_var.get())
                except ValueError:
                    majority_idx = None

                if majority_idx is not None:
                    maxp = proba.max(axis=1)
                    p_major = proba[:, majority_idx]
                    demote = (y_pred_enc != majority_idx) & (maxp < tau) & (p_major >= (maxp - margin))
                    if np.any(demote):
                        y_pred_enc = y_pred_enc.copy()
                        y_pred_enc[demote] = majority_idx
                    # push back to strings if we have a label encoder
                    if le is not None:
                        y_pred = le.inverse_transform(y_pred_enc)
                    else:
                        y_pred = y_pred_enc
                # --- end optional demotion at prediction ---
            else:
                logger.log(f"[ML][Predict] Thresholding is disabled in this prediction.")
            # If encoder available, decode label
            if le is not None: # LabelEncoder was loaded from the _labelencoder.joblib sidecar file.
                try:
                    y_pred = le.inverse_transform(y_pred) #convert integer predictions back to string labels (e.g., 2 → "F1H5N4S2")
                except:
                    pass
            elif hasattr(model, 'classes_'):
                y_pred = [model.classes_[i] if isinstance(i, int) else i for i in y_pred]
            else:
                logger.log("[ML][Predict][WARN] No label decoder available — predictions are raw integers") 
            # Append prediction labels back to original dataframe
            # FUTURE-API-HOOK: when Glycosmos API integration + dynamic prediction column naming lands post-freeze (Option A from B-17 item 4), rename column from `Predicted_Label` to `Predicted_<label_col>` (e.g., `Predicted_GlyToucan_ID`, `Predicted_WURCS`). Today static for freeze v1.10; reporter expects literal `pred_label` alias retained at the rename below.
            df['Predicted_Label'] = y_pred

            # 20260420 code review Stage 6 - make sure columns are there
            # --- ensure gate columns always exist --- (If we don't want this to show when gate is off, remove it.)
            if "ppm_precursor" not in df.columns:
                df["ppm_precursor"] = np.nan
            if "pred_ok" not in df.columns:
                df["pred_ok"] = True     # default: keep all rows

            # 20260420 code review Stage 7 - moved to stage 9
            # --- START: Prediction summary integration ---
            # 
            # 20260420 code review Stage 8 --- composition helpers (for precursor gate) ---
            # 20260517 fix B-01: KDN-first alternation (greedy); bare 'K' aliased to KDN for back-compat with legacy MAS-trained model labels.
            # Case-sensitive: uppercase F/H/N/S/G/A vs lowercase s/p are distinct semantic markers.
            # NOTE: future K-adduct work must use an adduct/metadata channel, not append 'K' to composition strings.
            _comp_pat = re.compile(r'(KDN|K|F|H|N|S|G|A|s|p)\s*([0-9]+)')
            def _canon_comp(s: str) -> str:
                """Turn any comp string into canonical 'F,H,N,S,G,KDN,A,s,p' order; omit zeros."""
                if not isinstance(s, str) or not s.strip():
                    return ""
                counts = {"F":0,"H":0,"N":0,"S":0,"G":0,"KDN":0,"A":0,"s":0,"p":0}
                for k,v in _comp_pat.findall(s):
                    if k == "K":  # 20260517 fix B-01: legacy single-char K aliased into KDN
                        k = "KDN"
                    counts[k] = counts.get(k, 0) + int(v)
                # canonical order
                parts = []
                for key in ("F","H","N","S","G","KDN","A","s","p"):
                    n = counts.get(key, 0)
                    if n > 0:
                        parts.append(f"{key}{n}")
                return "".join(parts)
            # ===============

            # --- START: reporter-friendly table (unchanged from v10) ---
            # 1) class names
            try:
                class_names = list(le.classes_) if le is not None else list(getattr(model, "classes_", []))
            except Exception:
                class_names = list(getattr(model, "classes_", []))

            
            pred_df = df.copy()
            pred_df = pred_df.rename(columns={"Predicted_Label": "pred_label"})  # reporter expects 'pred_label'

            # 2) proba  3way-monitoring and appends to pred_df: 1. return something 2. Python list/tuple (unlikely but defensive) 3. numpy array (what fires from sklearn)                               
            proba_cols = None
            if proba is not None and (isinstance(proba, (list, tuple)) or hasattr(proba, "shape")):
                try:
                    pred_df["proba_vector"] = [np.asarray(row, dtype=float) for row in proba]
                except Exception:
                    if class_names:
                        proba_cols = [f"proba_{c}" for c in class_names]
                    else:
                        proba_cols = [f"proba_{i}" for i in range(proba.shape[1])]
                    for j, col in enumerate(proba_cols):
                        pred_df[col] = proba[:, j]
            # --- END: reporter-friendly table ---

            # --- Packaging suffixes (so folders/files tell you MS1/gate) ---
            train_feats = train_feats or []
            ms1_used    = ("protonatedmass" in train_feats) or ("delta_ppm" in train_feats)
            ms1_suffix  = "_MS1feat" if ms1_used else "_noMS1"

            gate_on = bool(apply_pred_precursor_gate_var.get())
            try:
                gate_ppm = float(pred_precursor_ppm_var.get())  # 20260420 flood exceptions so we can log properly   #or 10)
                logger.log(f"[ML][Predict] Precursor gate ppm = {gate_ppm}")
                print(f"[ML][Predict] Precursor gate ppm = {gate_ppm}")
            except Exception:
                gate_ppm = 10.0
                logger.log(f"[ML][Predict][fallback] invalid Precursor gate ppm, set to default = {gate_ppm}")
                print(f"[ML][Predict][fallback] invalid Precursor gate ppm, set to default = {gate_ppm}")

            gate_suffix = f"_PG{int(gate_ppm)}ppm" if gate_on else "_noPG"
            method_suffix = ms1_suffix + gate_suffix

            method_base = _guess_method_basename_for_pack(predict_input_path, df)
            out_dir = Path(os.path.dirname(predict_input_path)) / f"{method_base}{method_suffix}"
            out_dir.mkdir(parents=True, exist_ok=True)

            # 20260420 code review: Stage 10 
            # --- Precursor gate (robust & optional) ---
            out_path_gate = None
            if gate_on:
                try:
                    # read in-silico and derive composition→[masses]
                    insilico_path = (insilico_csv_var.get() or "").strip()
                    if not insilico_path:
                        raise RuntimeError("No in-silico CSV selected.")
                    if logger: logger.log(f"[ML][Predict][info] Precursor gate uses in silico glycan library: {insilico_path}")
                    # 20260413 code review: change to existing main function call
                    # in silico file is not tsv, so prefer_tab is False (Claude Code suggests True). Need tests
                    lib = robust_read_csv(insilico_path, prefer_tab=False)
                    # 1) clean current headers
                    lib = clean_cols(lib) # line 7875: was lib.columns = clean_cols(lib.columns)
                    # 2) if it looks like a single-column read, try common delimiters
                    if len(lib.columns) == 1:
                        _one = lib.columns[0]
                        # log for debugging
                        if logger: logger.log(f"[ML][Predict][WARN] Detected in-silico looked single-col header: {_one!r}; reparsing…")
                        # try comma, semicolon (EU Excel), then tab
                        for _sep in (",", ";", "\t"):
                            try:
                                _tmp = pd.read_csv(insilico_path, sep=_sep, engine="python")
                                _tmp = clean_cols(_tmp) # line 7885: was _tmp.columns = clean_cols(_tmp.columns)
                                if len(_tmp.columns) > 1:
                                    lib = _tmp
                                    if logger: logger.log(f"[ML][Predict][WARN][auto] reparsed with sep='{_sep}': {lib.columns.tolist()}")
                                    break
                            except Exception as _e:
                                if logger: logger.log(f"[ML][Predict][ERROR 0] reparsed failed on '{_e}")
                                pass
                    if logger: logger.log(f"[ML][Predict][auto] in-silico columns (cleaned): {lib.columns.tolist()}")
                    # 3) build a case-insensitive lookup AFTER cleaning
                    # normalize headers
                    norm = {c.lower(): c for c in lib.columns}
                    logger.log(f"[ML][Predict][DEBUG][gate] norm keys: {list(norm.keys())}")                        
                    logger.log(f"[ML][Predict][DEBUG][gate] SO3 mapped to: {norm.get('so3')}, PO3H mapped to: {norm.get('po3h')}") 
                    # prefer string comp if present
                    comp_col = None
                    for cand in ("composition", "structure", "label", "comp_str"):
                        if cand in norm:
                            comp_col = norm[cand]; break
                    # pick mass column
                    mass_col = None
                    for cand in ("theoretical_mass", "mass", "Mass"):
                        if cand.lower() in norm:
                            mass_col = norm[cand.lower()]; break
                    if mass_col is None:
                        # tolerate Excel columns like 'Mass ' (trailing space) or different case
                        for c in lib.columns:
                            if c.strip().lower() in ("mass","theoretical_mass"):
                                mass_col = c; break
                    if mass_col is None:
                        if logger: logger.log(f"[ML][Predict][ERROR 2]In-silico CSV must have a Mass or theoretical_mass column.")
                        raise RuntimeError("In-silico CSV must have a Mass or theoretical_mass column.")
                    if logger: logger.log(f"[ML][Predict][info] pick the glycan library precursor mass column: {mass_col}")
                    # if no comp string, synthesize from counts
                    # If no composition string, synthesize one from monomer counts (FHN with modifiers A/s/p, plus G/K if present)
                                                            
                    # ADD: canonicalize library strings even if the CSV already has a composition column
                    def _canon_safe(x):
                        try:
                            return _canon_comp(str(x))
                        except Exception:
                            return str(x).strip()
                        
                    if comp_col is None:
                        def _get(row, key):
                            col = norm.get(key.lower()) #col = norm.get(key)
                            if col is None:
                                return 0
                            v = row.get(col, 0)
                            try:
                                return int(float(v))
                            except Exception:
                                return 0

                        def _counts_to_comp(row):
                            F  = _get(row, "Fuc")
                            H  = _get(row, "Hex")
                            N  = _get(row, "HexNAc")
                            NeuAc = _get(row, "NeuAc")
                            NeuGc = _get(row, "NeuGc")
                            KDN   = _get(row, "KDN")
                            HexA  = _get(row, "HexA")
                            SO3   = _get(row, "SO3") + _get(row, "s")
                            PO3H  = _get(row, "PO3H") + _get(row, "p")

                            parts = []
                            # Canonical order: F, H, N, (NeuAc/NeuGc/KDN), A (HexA), modifiers s/p
                            if F:     parts.append(f"F{F}")
                            if H:     parts.append(f"H{H}")
                            if N:     parts.append(f"N{N}")
                            if NeuAc: parts.append(f"S{NeuAc}")   # use 'S' for NeuAc if your label set uses it; otherwise remove this line
                            if NeuGc: parts.append(f"G{NeuGc}")
                            if KDN:   parts.append(f"KDN{KDN}")     # 20260517 fix B-01: 3-char KDN token (was "K{KDN}"); round-trips through _canon_comp
                            if HexA:  parts.append(f"A{HexA}")    # <-- HexA is 'A' in your schema
                            if SO3:   parts.append(f"s{SO3}")
                            if PO3H:  parts.append(f"p{PO3H}")
                            return "".join(parts)

                        lib = lib.copy()
                        lib["countedcomp"] = lib.apply(_counts_to_comp, axis=1)  # __comp -> countedcomp (for Python, it's not a class so no name mangling)
                        comp_col = "countedcomp"
                        logger.log(f"[ML][Predict][info] calculation finished. String composition is located at '{comp_col}', first 2:  {lib[comp_col].head(2).tolist()}") 
                    else:
                        # <<< NEW: canonicalize provided composition strings so they match prediction labels
                        lib = lib.copy()
                        lib[comp_col] = lib[comp_col].astype(str).map(_canon_safe)
                        logger.log(f"[Predict][gate] comp source: '{comp_col}', first 2:  {lib[comp_col].head(2).tolist()}") 

                    # build lookup composition → list of masses  (unchanged)
                    comp2masses = {}
                    for comp, mass in zip(lib[comp_col].astype(str), pd.to_numeric(lib[mass_col], errors="coerce")):
                        if pd.isna(mass):
                            continue
                        comp_key = _canon_safe(comp)
                        comp2masses.setdefault(comp_key, []).append(float(mass))

                    # apply gate (skip Non-glycan)
                    is_ng = pred_df["pred_label"].astype(str).eq("Non-glycan")  # ["pred_label"] becomes True for Non-glycan. Will be skipped later
                    pred_df["ppm_precursor"] = np.nan  #initialize the column to NaN always
                    pred_df["pred_ok"] = False # defaults to False, only when ppm passes, this column becomes true of that row.
                    
                    # ========================
                    # --- compute ppm for glycan rows only (robust: string LUT + counts fallback) ---
                    # 1) Canonical prediction strings, but keep the raw too for logging
                    pred_df["__canon_pred__"] = pred_df["pred_label"].astype(str).map(_canon_safe)
                    gly_mask = (~is_ng).to_numpy()
                    idx = np.nonzero(gly_mask)[0]
                    labels = pred_df.loc[gly_mask, "__canon_pred__"].astype(str).to_numpy()
                    obs    = pd.to_numeric(pred_df.loc[gly_mask, "protonatedmass"], errors="coerce").to_numpy()
                    best_ppm = np.full(labels.shape[0], np.nan, dtype=float)
                    keep     = np.zeros(labels.shape[0], dtype=bool)

                    # Create lib_counts for later comparison and extra report
                    lib_counts = lib.copy()
                    # Mass column (we already picked mass_col earlier)
                    lib_counts["__Mass__"] = pd.to_numeric(lib_counts[mass_col], errors="coerce")
                    # Parse a predicted label like "F1H4N3A1s1p1" → counts
                    def _parse_counts(label: str):
                        # 20260517 fix B-01: receives canonicalized labels from _canon_safe (A3) only; KDN-aware regex, no bare K (A2 already folded K -> KDN).
                        d = {"F":0,"H":0,"N":0,"S":0,"G":0,"KDN":0,"A":0,"s":0,"p":0}
                        for m in re.finditer(r'(KDN|[FHNSGAsp])(\d+)', str(label)):
                            ch, num = m.group(1), int(m.group(2))
                            if ch in d:
                                d[ch] += num
                        return d

                    #1) Normalize + numericize in-silico counts once
                    def _nz_int(colname):
                        if colname in lib_counts.columns:
                            lib_counts[colname] = pd.to_numeric(lib_counts[colname], errors="coerce").fillna(0).astype(int)
                        else:
                            lib_counts[colname] = 0

                    # standard headers to int, only on compositional field, DO NOT DO THIS ON MASS
                    _nz_int("Hex")
                    _nz_int("HexNAc")
                    _nz_int("Fuc")
                    _nz_int("NeuAc")
                    _nz_int("NeuGc")
                    _nz_int("KDN")
                    _nz_int("HexA")
                    _nz_int("SO3")
                    _nz_int("PO3H")

                    def _masses_by_counts(lbl: str):
                        cnt = _parse_counts(lbl)
                        m = lib_counts.loc[
                            (lib_counts["Fuc"]   == cnt["F"]) &
                            (lib_counts["Hex"]   == cnt["H"]) &
                            (lib_counts["HexNAc"]== cnt["N"]) &
                            (lib_counts["NeuAc"] == cnt["S"]) &
                            (lib_counts["NeuGc"] == cnt["G"]) &
                            (lib_counts["KDN"]   == cnt["KDN"]) &   # 20260517 fix B-01: dict key renamed K -> KDN (paired with _parse_counts change)
                            (lib_counts["HexA"]  == cnt["A"]) &
                            (lib_counts["SO3"]   == cnt["s"]) &
                            (lib_counts["PO3H"]  == cnt["p"]),
                            "__Mass__"
                        ].dropna()
                        return sorted(m.unique().tolist())
                    
                    lut_hits = 0
                    counts_hits = 0
                    misses = 0
                    gppm = float(gate_ppm)
                    for j, (lbl, pm) in enumerate(zip(labels, obs)):
                        if not lbl or np.isnan(pm):
                            continue
                        masses = comp2masses.get(lbl, [])
                        if masses:
                            lut_hits += 1
                        else:
                            masses = _masses_by_counts(lbl)
                            if masses:
                                counts_hits += 1
                            else:
                                misses += 1
                                continue

                        m = np.asarray(masses, dtype=float)
                        diffs = (pm - m) / m * 1e6  # PPM error against each candidate mass
                        k = np.nanargmin(np.abs(diffs)) # index of closest match
                        best_ppm[j] = float(diffs[k]) # store best PPM
                        passed = abs(best_ppm[j]) <= gppm # within tolerance?
                        keep[j] = passed # mark pass/fail

                    # 20260420 code review: update the report output as summary
                    logger.log(f"[ML][Predict][info] Precursor gate lookup summary: string_match={lut_hits}, counts_fallback={counts_hits}, no_match={misses}")
                    # write back only for glycan rows
                    pred_df.loc[idx, "ppm_precursor"] = best_ppm
                    pred_df.loc[idx, "pred_ok"]       = keep


                    # --- GATE DIAGNOSTICS ---
                    try:
                        # same canonicalizer for both sides
                        # build a quick "masses found?" map
                        # If you created comp2masses above, reuse it; otherwise reconstruct it here:
                        # comp2masses = {...}  # already built earlier

                        gly = pred_df[ pred_df["pred_label"].astype(str) != "Non-glycan" ].copy()
                        gly["canon_pred"]  = gly["pred_label"].astype(str).map(_canon_safe)
                        gly["masses_found"] = gly["canon_pred"].map(lambda c: len(comp2masses.get(c, [])))
                        gly["has_pm"]      = pd.to_numeric(gly["protonatedmass"], errors="coerce").notna()
                        gly["kept"]        = gly["pred_ok"].fillna(False)

                        kept  = int(gly["kept"].sum())
                        total = int(len(gly))
                        dropped = total - kept
                        no_lut = int((gly["masses_found"] == 0).sum())
                        no_pm  = int((~gly["has_pm"]).sum())
                        bad_ppm = int(((gly["masses_found"] > 0) & gly["has_pm"] & ~gly["kept"]).sum())

                        if logger:
                            logger.log(f"[Predict][gate] glycan rows total={total}, kept={kept}, dropped={dropped}, "
                                    f"no_lut={no_lut}, no_pm={no_pm}, ppm_exceeded={bad_ppm}")

                        # 20260420 code review: write both debug -> rename to dropped csv (1000 -> full) to reveal dropped columns (and we know reasons) and a pure passed csv
                        debug_cols = ["MS2scan_no","pred_label","canon_pred","protonatedmass",
                                    "ppm_precursor","masses_found","kept"]
                        gate_dropped = gly[~gly["kept"]][debug_cols]
                        gate_passed = gly[gly["kept"]][debug_cols]
                        dropped_path = (out_dir / (Path(predict_input_path).stem + "_gate_debug.csv")).as_posix()
                        gate_dropped.to_csv(dropped_path, index=False)
                        passed_path = (out_dir / (Path(predict_input_path).stem + "_gate_passed.csv")).as_posix()
                        gate_passed.to_csv(passed_path, index=False)
                        if logger: logger.log(f"[ML][Predict][info] After precursor gating, dropped rows summary is saved to → {dropped_path}")
                        if logger:logger.log(f"[ML][Predict][info] After precursor gating, passed rows summary is saved to → {passed_path} (count = {len(gate_passed)})")   
                    # --- EXTRA DIAGNOSTICS: for no_lut rows, fetch Mass candidates by counts-match ---
                        # 4) Take the first 30 no-lut rows and attach candidate Masses
                        no_lut_head = gly[gly["masses_found"] == 0].copy().head(30)
                        if not no_lut_head.empty:
                            no_lut_head["counts"] = no_lut_head["canon_pred"].map(_parse_counts)
                            no_lut_head["Mass_candidates"] = no_lut_head["canon_pred"].map(_masses_by_counts)   # 20260517 fix B-01 (Codex Finding 11): map over canon_pred strings, not the counts-dict column
                            # print to logger
                            if logger:
                                logger.log(f"[Predict][gate][no_lut] sample {len(no_lut_head)} rows with Mass candidates:")
                                for r in no_lut_head[["MS2scan_no","pred_label","canon_pred","protonatedmass","Mass_candidates"]].itertuples(index=False):
                                    logger.log("[no_lut] "
                                            f"scan={getattr(r,'MS2scan_no',None)}, "
                                            f"pred='{getattr(r,'pred_label',None)}' "
                                            f"canon='{getattr(r,'canon_pred',None)}' "
                                            f"pm={getattr(r,'protonatedmass',None)} "
                                            f"MassCandidates={getattr(r,'Mass_candidates',None)}")
                            # save a CSV so you can compare in Excel
                            nolut_bycounts_path = (out_dir / (Path(predict_input_path).stem + "_gate_nolut_bycounts_debug.csv")).as_posix()
                            no_lut_head[["MS2scan_no","pred_label","canon_pred","protonatedmass","Mass_candidates"]].to_csv(nolut_bycounts_path, index=False)
                            if logger: logger.log(f"[ML][Predict][gate] wrote by-counts no-lut debug → {nolut_bycounts_path}")

                    except Exception as _e:
                        if logger: logger.log(f"[ML][Predict][ERROR 0] Post-precursor gating diagnostic collection failed: {_e}")

                    # Add prediction ok or not (dropped or not) and ppm final values to full table
                    df["pred_ok"] = pred_df["pred_ok"].values
                    df["ppm_precursor"] = pred_df["ppm_precursor"].values 
                    
                    # Export Gate report
                    out_path_gate = (out_dir / (Path(predict_input_path).stem + f"_predicted{ms1_suffix}_PG{int(gate_ppm)}ppm_only.csv")).as_posix()
                    if out_path_gate: # 20260420 Export the Gate-passed only list, so no need to dig into full df output
                        gated_df = df[df["pred_ok"] == True].copy()
                        gated_df.to_csv(out_path_gate, index=False)
                        logger.log(f"[Predict][gate] Scans passed gated predictions is saved → {out_path_gate} ,( count = {len(gated_df)} )")
                    # For generating simple gated count table: composition, n (convenient passed class count)
                    pred_df_for_report = pred_df[pred_df["pred_ok"]].copy()
                    gate_counts = (pred_df_for_report
                                .groupby("pred_label", dropna=False)
                                .size()
                                .reset_index(name="count")
                                .sort_values("count", ascending=False))
                    gate_counts_path = (out_dir / (Path(predict_input_path).stem + "_PG_counts.csv")).as_posix()
                    gate_counts.to_csv(gate_counts_path, index=False)
                    if logger: logger.log(f"[ML][Predict][gate] Summary of passed class count is saved → {gate_counts_path}")

                except Exception as _e:
                    # Gate skipped → just report everything
                    if logger: logger.log(f"[Predict] Precursor gate skipped: {_e}")
                    pred_df_for_report = pred_df.copy()
            else:
                pred_df_for_report = pred_df.copy()

            # --- Harmonize probability vectors so all rows have same length ---
            import ast
            def _parse_vec(x):
                s = str(x).strip()
                if s in ("", "nan", "None"):
                    return []
                try:
                    try:
                        v = json.loads(s)
                    except Exception:
                        v = ast.literal_eval(s)
                    if isinstance(v, (list, tuple)):
                        return [float(z) for z in v]
                except Exception:
                    pass
                return []

            if "proba_vector" in pred_df_for_report.columns:
                vecs = pred_df_for_report["proba_vector"].map(_parse_vec)
                target = len(class_names) if (class_names and len(class_names) > 0) else max((len(v) for v in vecs), default=0)
                if target > 0:
                    vecs = vecs.map(lambda v: (v + [0.0]*target)[:target])
                    pred_df_for_report["proba_vector"] = vecs#.map(json.dumps)
                    proba_cols = None
            else:
                if proba_cols:
                    for c in proba_cols:
                        if c not in pred_df_for_report.columns:
                            pred_df_for_report[c] = 0.0

            # --- Reporter: summarize & write (keep run_sum as dict) ---
            if not enable_thresh_var.get():
                logger.log(f"[ML][Predict] Thresholding is disabled by default")
                print(f"[ML][Predict] Thresholding disabled by default")
            else:
                logger.log(f"[ML][Predict] Thresholding is enabled: tau {thresh_val_var.get()}, margin {margin_val_var.get()}")
                print(f"[ML][Predict] Thresholding is enabled: tau {thresh_val_var.get()}, margin {margin_val_var.get()}")
            params = ReportParams(
                tau=float(thresh_val_var.get()) if enable_thresh_var.get() else 0.0,
                margin=float(margin_val_var.get()) if enable_thresh_var.get() else 1.0,
                topk=5,
                sample_cols=("experiment_title", "sample_name")
            )   
            pred_rows, class_sum, by_sample_sum, run_sum = summarize_predictions(
                pred_df_for_report,
                class_names=class_names if class_names else None,
                proba_cols=proba_cols,   # None if using 'proba_vector'
                params=params,
            )
            if isinstance(run_sum, dict):
                run_sum["method_suffix"] = method_suffix

            arts = write_prediction_report(
                out_dir=out_dir,
                pred_rows=pred_rows,
                class_summary=class_sum,
                by_sample=by_sample_sum,
                run_summary=run_sum,
                params=params,
            )

            # Save tables (all rows; gated-only if present) + user popup
            # This one is intact prediction table
            out_path_all = (out_dir / (Path(predict_input_path).stem + f"_predicted{method_suffix}.csv")).as_posix()
            df.to_csv(out_path_all, index=False)
            if out_path_gate and "pred_ok" in df.columns and df["pred_ok"].any():
                # already written above when gate ran; keep this guard for safety
                pass
            ms1_line  = f"MS1 feature (expected by model): {'YES' if ms1_used else 'NO'}"
            gate_line = f"Precursor gate: {'ON' if gate_on else 'OFF'}" + (f" ({int(gate_ppm)} ppm)" if gate_on else "")
            try:
                show_prediction_summary_popup(root, run_summary=run_sum, artifacts=arts, class_summary_df=class_sum)
            except Exception as _e:
                print("[warn] Failed to show summary popup:", _e)
            extra = f"\nGated table:\n{out_path_gate}" if (out_path_gate and os.path.exists(out_path_gate)) else ""
            ## --- END: Prediction summary integration ---
            if logger: logger.log(f"[ML][Predict][info] Prediction Complete. {ms1_line} {gate_line} {out_dir}")
            messagebox.showinfo("Prediction Complete", f"{ms1_line}\n{gate_line}\n\nPacked into folder:\n{out_dir}\n\nMain table:\n{out_path_all}{extra}")

        except Exception as e:
            if logger: logger.log(f"[ML][Predict][Failed] Prediction Failed due to {str(e)}")
            messagebox.showerror("Prediction Failed", str(e))
    # Run prediction block ends
    # ==========================

    def _short_id(s: str) -> str:
        import hashlib
        return hashlib.sha1((s or "").encode("utf-8")).hexdigest()[:8]

    # 20260418 code review: possible for CLI implementation in future. Built 20250915 
    def create_unlabeled_from_method(method_path: str, default_ppm: int | None = None) -> str:
        import json, hashlib, pandas as pd, os
        with open(method_path, "r", encoding="utf-8") as f:
            m = json.load(f)

        # Resolve inputs
        parents   = m.get("parents") or {}
        ionblock  = m.get("ionlist") or {}
        raw_csv_s = parents.get("converted_csv") or ""
        ion_path_s= ionblock.get("path") or ""
        ion_sheet = ionblock.get("sheet") or "ionlist"
        ppm       = default_ppm or ionblock.get("ppm_tolerance") or 20

        # If converted_csv missing or file not found -> ask user and persist back
        raw_csv_p = pathcanon.to_native_path(raw_csv_s) if raw_csv_s else None
        if not raw_csv_p or not raw_csv_p.exists():
            picked = _pick_file_cli_or_gui("Select the *converted* CSV (ms2_*.csv)")
            if not picked:
                raise FileNotFoundError("No converted CSV provided.")
            raw_csv_p = pathcanon.to_native_path(picked)
            # write back to method.json in POSIX form
            parents["converted_csv"] = pathcanon.to_posix_str(raw_csv_p)
            m["parents"] = parents
            with open(method_path, "w", encoding="utf-8") as f:
                json.dump(m, f, indent=2, ensure_ascii=False)

        # Validate ion list path
        ion_path_p = pathcanon.to_native_path(ion_path_s)
        if not ion_path_p.exists():
            picked = _pick_file_cli_or_gui("Select ion list (CSV/XLSX)",
                                        (("CSV", "*.csv"), ("Excel", "*.xlsx;*.xls"), ("All files", "*.*")))
            if not picked:
                raise FileNotFoundError("No ion list provided.")
            ion_path_p = pathcanon.to_native_path(picked)
            ionblock["path"] = pathcanon.to_posix_str(ion_path_p)
            m["ionlist"] = ionblock
            with open(method_path, "w", encoding="utf-8") as f:
                json.dump(m, f, indent=2, ensure_ascii=False)

        # Load fragments (your existing smart loader)
        frags = extract_fragment_masses(ion_path_p.as_posix(), sheet_name=ion_sheet)

        # Use your existing extractor
        feature_df = extract_ion_intensities(raw_csv_p.as_posix(), frags, ppm=float(ppm))

        # Best-effort UID (no-op if columns missing)
        try:
            exp_id  = _short_id(m.get("experiment_title", ""))
            samp_id = _short_id(m.get("sample_name", ""))
            if "MS2scan_no" in feature_df.columns:
                s = feature_df["MS2scan_no"].astype(int).astype(str).str.zfill(6)
                feature_df = feature_df.copy()
                feature_df["UID"] = s.map(lambda x: f"{exp_id}:{samp_id}:{x}")
        except Exception:
            pass

        out_path = raw_csv_p.with_name(raw_csv_p.stem + f"_unlabeled_ppm{ppm}.csv")
        feature_df.to_csv(out_path, index=False)
        return out_path.as_posix()

    # 20260420 Code review
    # --- creating unlabeled datasets (PL or manual) ---

    def extract_fragment_masses(ion_path: str, sheet_name : str ="ionlist"):
        """
        Load a column of fragment masses from either:
        • Excel: prefers a sheet named 'ionlist' (or the provided sheet_name)
        • CSV:   uses a column named 'mass' (case-insensitive)
        Returns: list[float]
        """

        p = (ion_path or "").strip()
        if not p or not os.path.exists(p):
            raise FileNotFoundError(f"Ion list not found: {p}") 

        if p.lower().endswith((".xlsx", ".xls")):
            xls = pd.ExcelFile(p, engine="openpyxl")
            # try requested sheet first; else try common names; else first sheet that has a mass-like column
            sheet_to_use = None
            wanted = {sheet_name.lower(), "ionlist", "ions", "ion_list"}
            for s in xls.sheet_names:
                if s.strip().lower() in wanted:
                    sheet_to_use = s
                    break
            if sheet_to_use is None:
                # fallback: first sheet with something like 'mass'/'mz'
                for s in xls.sheet_names:
                    df = xls.parse(s)
                    cols = {str(c).strip().lower() for c in df.columns}
                    if any(k in cols for k in ("mass", "mz", "ion_mz", "m/z")):
                        sheet_to_use = s
                        break
            df = xls.parse(sheet_to_use or xls.sheet_names[0])
        else:
            # CSV ion list
            df = pd.read_csv(p, engine="python")

        # normalize to a 'mass' column
        colmap = {str(c).strip().lower(): c for c in df.columns}
        for k in ("mass", "mz", "ion_mz", "m/z"):
            if k in colmap:
                series = df[colmap[k]].dropna().astype(float)
                return series.tolist()
        return []


    def extract_ion_intensities(tsv_like_path, fragment_masses, ppm=20.0):
        """
        Read converted MS2 file (tab-separated but *.csv) and
        return a DataFrame with: protonatedmass, MS2scan_no, and one column per ion.
        """
        #20260413 code review: use main function for reading tsv. Need test
        df = robust_read_csv(tsv_like_path, prefer_tab=True) #read sep="\t" first 
        # tolerate either literal lists or python-literal strings in peaklist/peakintensity
        def _parse_list(x):
            if isinstance(x, (list, tuple)):
                return list(map(float, x))
            s = str(x).strip()
            # tolerate "[]" or "(...)" or "1,2,3"
            if s.startswith(("(", "[")) and s.endswith((")", "]")):
                try:
                    return [float(v) for v in s.strip("()[]").split(",") if v.strip() != ""]
                except Exception:
                    pass
            try:
                return [float(v) for v in s.split(",") if v.strip() != ""]
            except Exception:
                return []

        out_rows = []
        for _, row in df.iterrows():
            peaks = _parse_list(row.get("peaklist", ""))
            intensities = _parse_list(row.get("peakintensity", ""))
            if not peaks or not intensities:
                continue
            if len(peaks) != len(intensities):
                # skip malformed line
                continue

            feature_row = {
                "protonatedmass": row.get("protonatedmass", 0.0),
                "MS2scan_no": int(row.get("MS2scan_no", 0)),
            }
            peak_array = np.array(peaks, dtype=float)
            intensity_array = np.array(intensities, dtype=float)

            for target in fragment_masses:
                target = float(target)
                ppm_tol = target * float(ppm) / 1e6
                mask = np.abs(peak_array - target) <= ppm_tol
                if np.any(mask):
                    feature_row[str(target)] = float(np.log10(np.max(intensity_array[mask])) + 1.0) #fixed 20260122
                else:
                    feature_row[str(target)] = 1.0  # keep "1.0" default you were using

            out_rows.append(feature_row)
        return pd.DataFrame(out_rows)
    # ===
    # 20260420 code review, refactored by Claude Code with supervision from Henry.  
    def create_unlabeled_dataset():
        """
        Builds unlabeled feature CSVs for each eligible sample in an experiment JSON.
        Supports both v1 (samples as list) and legacy (samples as dict) formats.
        Eligibility per sample:
        - has a converted MS2 file (csv key or from method JSON inputs)
        - AND has an ion source (excel, ionlist_path, or from method JSON inputs)
        """
        import hashlib

        # --- Feature: Input source (experiment JSON via file dialog) ---
        exp_path = filedialog.askopenfilename(filetypes=[("Experiment JSON", "*.json")])
        if not exp_path:
            return

        try:
            with open(exp_path, "r", encoding="utf-8") as f:
                exp = json.load(f)
        except Exception as e:
            messagebox.showerror("Failed to Load JSON", str(e))
            return

        # --- Feature: PPM source (ask user explicitly) ---
        ppm_value = simpledialog.askinteger(
            "PPM Tolerance", "Enter PPM tolerance (e.g., 20):", minvalue=1, maxvalue=200
        )
        if ppm_value is None:
            return

        # --- Feature: JSON format (support both v1 list and legacy dict) ---
        raw_samples = exp.get("samples") or {}
        exp_base_dir = os.path.dirname(exp_path)
        # 20260521 confirmed list object
        if isinstance(raw_samples, list):
            # v1 format: list of dicts with "sample_name" and optionally "methods"
            for s in raw_samples:
                if not isinstance(s, dict):
                    continue
                sname = s.get("sample_name")
                if not sname:
                    continue
                # start with any flat keys present
                info = {}
                for k in ("csv", "excel", "ionlist_path", "ion_sheet"):
                    if s.get(k):
                        info[k] = s[k]
                # v1: try to resolve file paths from method JSON if flat keys are missing
                if (not info.get("csv") or not (info.get("excel") or info.get("ionlist_path"))):
                    methods = s.get("methods") or []
                    for m in methods:
                        mref = m.get("method_ref") or {}
                        mpath = mref.get("path") or m.get("path")
                        if not mpath or not isinstance(mpath, str):
                            continue
                        # resolve relative to experiment JSON folder
                        if not os.path.isabs(mpath):
                            mpath = os.path.normpath(os.path.join(exp_base_dir, mpath))
                        if not os.path.exists(mpath):
                            continue
                        try:
                            md = load_typed_json(mpath, expected_type=JSON_TYPE_METHOD, allow_legacy=True, context="[Unlabeled] ")
                            _, _, tree_entry, _ = normalize_method_json(md, mpath, base_dir=os.path.dirname(mpath))
                            if not info.get("csv") and tree_entry.get("csv"):
                                info["csv"] = tree_entry["csv"]
                            if not (info.get("excel") or info.get("ionlist_path")):
                                info["ionlist_path"] = tree_entry.get("ionlist_path")
                                info["excel"] = tree_entry.get("excel")
                            if tree_entry.get("ion_sheet"):
                                info["ion_sheet"] = tree_entry["ion_sheet"]
                        except Exception as e:
                            logger.log(f"[Unlabeled] Failed to load method for {sname}: {e}")
                        break  # use first valid method
                info["_v1_entry"] = s  # keep reference for writeback
                samples[sname] = info
        else:
            # legacy format: dict keyed by sample name
            samples = raw_samples

        exp_title = ""
        if isinstance(exp.get("experiment"), dict):
            exp_title = exp["experiment"].get("title", "")
        elif isinstance(exp.get("experiment"), str):
            exp_title = exp["experiment"]
        exp_id = _short_id(exp_title)

        # --- Batch processing ---
        made = 0
        missing = []
        seen_out = set()
        made_paths = []

        for sample_id, sample_info in samples.items():
            raw_csv   = sample_info.get("csv")
            ion_src   = sample_info.get("excel") or sample_info.get("ionlist_path")
            ion_sheet = sample_info.get("ion_sheet") or "ionlist"
            print(f"[DEBUG] sample_info: {sample_info}")
            # --- Feature: Missing file handling (skip and report) ---
            if not raw_csv or not ion_src:
                missing.append(f"{sample_id} (csv={'yes' if raw_csv else 'no'}, ion={'yes' if ion_src else 'no'})")
                logger.log(f"[Unlabeled] Skipped {sample_id}: csv={bool(raw_csv)}, ion={bool(ion_src)}")
                continue

            if not os.path.exists(raw_csv):
                missing.append(f"{sample_id} (csv not found: {os.path.basename(raw_csv)})")
                logger.log(f"[Unlabeled] Skipped {sample_id}: csv not found at {raw_csv}")
                continue

            if not os.path.exists(ion_src):
                missing.append(f"{sample_id} (ion list not found: {os.path.basename(ion_src)})")
                logger.log(f"[Unlabeled] Skipped {sample_id}: ion list not found at {ion_src}")
                continue

            try:
                # --- Feature: Ion loader (extract_fragment_masses — flexible) ---
                ion_list = extract_fragment_masses(ion_src, sheet_name=ion_sheet)
                if not ion_list:
                    raise ValueError(f"No 'mass' column found in {os.path.basename(ion_src)}")
                logger.log(f"[Unlabeled] {sample_id}: loaded {len(ion_list)} fragment masses")

                # --- Feature: Feature extractor (extract_ion_intensities) ---
                feats = extract_ion_intensities(raw_csv, ion_list, ppm=ppm_value)

                # --- Feature: UID generation (exp:sample:scan provenance) ---
                samp_id = _short_id(sample_id)
                if "MS2scan_no" in feats.columns:
                    feats = feats.copy()
                    scan_str = feats["MS2scan_no"].astype(str).str.zfill(6)
                    feats.insert(0, "UID", scan_str.map(lambda x: f"{exp_id}:{samp_id}:{x}"))

                # --- Feature: Output naming ({csv_stem}_unlabeled_ppm{ppm}.csv) ---
                out_dir = Path(raw_csv).resolve().parent
                base    = Path(raw_csv).stem
                out_path = out_dir / f"{base}_unlabeled_ppm{ppm_value}.csv"

                # --- Feature: Collision guard (append sample_id if path reused) ---
                if str(out_path) in seen_out:
                    safe_id = re.sub(r'[^A-Za-z0-9._-]+', '_', sample_id)
                    out_path = out_dir / f"{base}__{safe_id}_unlabeled_ppm{ppm_value}.csv"

                # write
                feats.to_csv(str(out_path), index=False, encoding="utf-8")
                logger.log(f"[Unlabeled] {sample_id}: saved → {out_path} ({len(feats)} rows)")

                # --- Feature: Write back to JSON (record unlabeled_csv path) ---
                sample_info["unlabeled_csv"] = str(out_path)
                # for v1 format, also update the original list entry
                if "_v1_entry" in sample_info:
                    sample_info["_v1_entry"]["unlabeled_csv"] = str(out_path)

                seen_out.add(str(out_path))
                made_paths.append(str(out_path))
                made += 1

            except Exception as e:
                missing.append(f"{sample_id} (error: {e})")
                logger.log(f"[Unlabeled] {sample_id} failed: {e}")

        # --- Feature: Completion report (messagebox + logger) ---
        if made == 0:
            messagebox.showwarning(
                "No Samples Processed",
                "No eligible samples were found in this experiment file.\n"
                "Tip: each sample needs a converted CSV + an ion list (Excel or CSV)."
            )
            logger.log(f"[Unlabeled] 0 samples processed from {exp_path}")
        else:
            # persist updated experiment file with recorded unlabeled paths
            try:
                # clean up internal keys before saving
                if isinstance(raw_samples, list):
                    for s in raw_samples:
                        if isinstance(s, dict):
                            s.pop("_v1_entry", None)
                with open(exp_path, "w", encoding="utf-8") as f:
                    json.dump(exp, f, indent=2, ensure_ascii=False)
                logger.log(f"[Unlabeled] Updated experiment JSON: {exp_path}")
            except Exception as e:
                logger.log(f"[Unlabeled] Failed to update experiment JSON: {e}")

            msg = [f"Created {made} unlabeled dataset(s)."]
            if missing:
                msg.append(f"\nSkipped ({len(missing)}):")
                for m in missing[:10]:
                    msg.append(f"  - {m}")
                if len(missing) > 10:
                    msg.append(f"  ... and {len(missing) - 10} more")
            messagebox.showinfo("Done", "\n".join(msg))
            logger.log(f"[Unlabeled] Done: {made} created, {len(missing)} skipped")

    # ML Analysis window init
    subwin = tk.Toplevel(root)
    subwin.title("ML Analysis")
    subwin.geometry("640x496")
    notebook = ttk.Notebook(subwin)
    notebook.pack(fill="both", expand=True)

    # --- Tab 1: Train Model ---
    train_tab = ttk.Frame(notebook)
    notebook.add(train_tab, text="Train Model")
    # Step 1: load trainable dataset
    tk.Label(train_tab, text="Step 1: Load Trainable Dataset (.csv)").grid(row=0, column=0, sticky="w", padx=10, pady=5)
    train_load_button = tk.Button(train_tab, text="Select CSV File", command=select_train_csv)
    train_load_button.grid(row=0, column=1, padx=5, pady=5)
    # options rows, init the value here
    is_pseudolabel_var = tk.BooleanVar(value=False) #default to False
    include_mass_train_var = tk.BooleanVar(value=False) #default to False
    # option rows: CGA dataset / Use mass feature for training
    options_row = ttk.Frame(train_tab)
    options_row.grid(row=1, column=0, columnspan=6, sticky="w", padx=10, pady=(0, 5))
    #CGA checker, does nothing for now
    tk.Checkbutton(options_row, text="CGA dataset?", variable=is_pseudolabel_var).pack(side="left", padx=(0, 16))
    #mass feature training
    tk.Checkbutton(options_row, text="Use protonated mass as a model feature (Non-glycan = 0)", variable=include_mass_train_var).pack(side="left", padx=(0, 16))
    # Step 2: choose label column
    tk.Label(train_tab, text="Step 2: Select Label Column").grid(row=2, column=0, sticky="w", padx=10, pady=5)
    label_dropdown = ttk.Combobox(train_tab, values=['Structure', 'IUPACname(optional)', 'Glycanannotation2', 'GlyToucan ID'])
    label_dropdown.set("Structure")
    label_dropdown.grid(row=2, column=1, padx=5, pady=5)
    # Step 3: choose classifier
    tk.Label(train_tab, text="Step 3: Choose Classifier").grid(row=3, column=0, sticky="w", padx=10, pady=5)
    classifier_dropdown = ttk.Combobox(train_tab, values=['Random Forest', ' Coming soon'])
    classifier_dropdown.set("Random Forest")
    classifier_dropdown.grid(row=3, column=1, padx=5, pady=5)
    # old tk button for model selection - buggy : I'd rather changing it to dropdown
    #classifier_var = tk.IntVar(value=0)
    #rf_button = tk.Radiobutton(train_tab, text="Random Forest (✔ functional)",variable=classifier_var, value=0)       
    #xgb_button = tk.Radiobutton(train_tab, text="XGBoost (placeholder)",variable=classifier_var, value=1) 
    #svm_button = ttk.Radiobutton(train_tab, text="SVM (placeholder)", variable=classifier_var, value="svm")
    #knn_button = ttk.Radiobutton(train_tab, text="KNN (placeholder)", variable=classifier_var, value="knn")
    # 20260418 code review: Set default button status. Tkinter can't automatically bind. Change from stringvar to intvar to see if it fixes
    #classifier_var.set(0)     
    #rf_button.grid(row=3, column=1, sticky="w")
    #xgb_button.grid(row=4, column=1, sticky="w")
    #svm_button.grid(row=5, column=1, sticky="w")
    #knn_button.grid(row=6, column=1, sticky="w")
    # Step 4: Train/Test Parameters
    tk.Label(train_tab, text="Step 4: Train/Test Parameters").grid(row=7, column=0, sticky="w", padx=10, pady=5)
    ttk.Label(train_tab, textvariable=ml_summary_var, wraplength=400).grid(row=8,column=0, columnspan=6, sticky="w", padx=10, pady=(0,6))  
    tk.Button(train_tab, text="Set Parameters / Train the Model", command=open_train_settings).grid(row=7, column=1, padx=5, pady=5)
    # Step 5: Train model
    train_button = tk.Button(train_tab, text="Train Model", command=train_model, bg="#CCE5FF")
    train_button.grid(row=9, column=0, columnspan=2, pady=10)
    # Info field, a intermediate plan before ML workflow is integrated to prepare dataset window (version 1.2+)
    tk.Label(train_tab, text="Trainable File Info (Origin Tracking)").grid(row=10, column=0, columnspan=2, sticky="w", padx=10, pady=(15, 5))
    origin_info = tk.Text(train_tab, height=4, width=70, state="disabled", wrap="word")
    origin_info.grid(row=11, column=0, columnspan=2, padx=10, pady=5)
    # Extra: pool datasets for training
    tk.Label(train_tab, text="Combine Datasets for Training").grid(row=12, column=0, columnspan=2, sticky="w", padx=10, pady=(15, 5))
    tk.Button(train_tab, text="Select Datasets", command=lambda: combine_trainable_datasets_ui(root)).grid(row=12, column=1, columnspan=2, padx=10, pady=5)

    # --- Tab 2: Predict ---
    predict_tab = ttk.Frame(notebook)
    notebook.add(predict_tab, text="Predict")

    tk.Label(predict_tab, text="Step 1: Load Trained Model (.skops/.joblib/.pkl)").grid(row=0, column=0, sticky="w", padx=10, pady=5)
    predict_model_button = tk.Button(predict_tab, text="Select Model File", command=select_model_file)
    predict_model_button.grid(row=0, column=1, padx=5, pady=5)

    tk.Label(predict_tab, text="Step 2: Prepare Unlabeled Input Dataset").grid(row=1, column=0, sticky="w", padx=10, pady=5)
    create_unlabeled_button = tk.Button(predict_tab, text="Create unlabeled dataset of certain experiment", command=create_unlabeled_dataset)
    create_unlabeled_button.grid(row=2, column=0, columnspan=2, padx=10, pady=5)

    predict_input_button = tk.Button(predict_tab, text="Load unlabeled dataset", command=select_predict_input)
    predict_input_button.grid(row=3, column=0, columnspan=2, padx=10, pady=5)

    # === PRED GATE: state vars (default ON) ===
    apply_pred_precursor_gate_var = tk.BooleanVar(value=True)
    pred_precursor_ppm_var        = tk.StringVar(value="10")
    insilico_csv_var              = tk.StringVar(value="")

    # === PRED GATE: UI frame ===
    _next = predict_tab.grid_size()[1]  # safe next free row
    pred_gate = ttk.LabelFrame(predict_tab, text="Precursor (MS1) Gate for Predictions", padding=(10,8))
    pred_gate.grid(row=_next, column=0, columnspan=3, sticky="we", padx=10, pady=(6,6))

    ttk.Checkbutton(pred_gate, text="Apply precursor gate", variable=apply_pred_precursor_gate_var)\
        .grid(row=0, column=0, sticky="w")

    ttk.Label(pred_gate, text="Tolerance (ppm):").grid(row=0, column=1, padx=(16,4), sticky="e")
    ttk.Entry(pred_gate, width=8, textvariable=pred_precursor_ppm_var).grid(row=0, column=2, sticky="w")

    ttk.Label(pred_gate, text="In-silico CSV (composition ↔ theoretical_mass):")\
        .grid(row=1, column=0, columnspan=2, sticky="w", pady=(6,0))
    ttk.Entry(pred_gate, width=60, textvariable=insilico_csv_var)\
        .grid(row=2, column=0, columnspan=2, sticky="we", pady=(2,6))

    def _pick_insilico_for_predict():
        p = filedialog.askopenfilename(
            title="Select in-silico CSV",
            filetypes=[("CSV","*.csv"), ("All files","*.*")]
        )
        if p:
            insilico_csv_var.set(p)

    ttk.Button(pred_gate, text="Select in-silico CSV file", command=_pick_insilico_for_predict)\
        .grid(row=2, column=2, padx=(8,0), sticky="w")

    pred_gate.columnconfigure(0, weight=1)
    # Run button goes AFTER the gate
    _next = predict_tab.grid_size()[1]
    predict_button = tk.Button(predict_tab, text="Run Prediction", bg="#D5F5E3", command=run_prediction)
    predict_button.grid(row=_next, column=0, columnspan=2, pady=10)
    # Combine section follows
    _next = predict_tab.grid_size()[1]
    # Future plan: for lazy guys doing 1 prediction on many datasets (need to share same feature matrix columns)
    #tk.Label(predict_tab, text="(🔜) Combine Datasets for Prediction").grid(row=_next, column=0, columnspan=2, sticky="w", padx=10, pady=(15, 5))
    #tk.Button(predict_tab, text="[Placeholder] Combine Datasets").grid(row=_next+1, column=0, columnspan=2, padx=10, pady=5)
    close_button = tk.Button(subwin, text="Close", command=subwin.destroy)
    close_button.pack(pady=5)

# main GUI
#icon
ico_path = os.path.join(os.path.dirname(__file__), 'GlycoMSPlogo.ico')
png_path = os.path.join(os.path.dirname(__file__), 'GlycoMSPlogo.png')

# Initialize GUI
root = tk.Tk()
if platform.system() == "Windows" and os.path.exists(ico_path):
    root.iconbitmap(default=ico_path)
elif platform.system() in ("Darwin", "Linux") and os.path.exists(png_path):
    icon_img = tk.PhotoImage(file=png_path)
    root.iconphoto(True, icon_img)

def on_closing():
    if messagebox.askokcancel("Quit", "Do you really want to quit?"):
        logger.log("Application closed by user.")
        root.destroy()  # Clean exit
    #if logger.entries:
        #logger.save("autosave.log", include_debug=True)
        #print("Log saved to autosave.log")

def save_log_to_file():
    filename = filedialog.asksaveasfilename(defaultextension=".log", filetypes=[("Log files", "*.log")])
    if filename:
        logger.save(filename, include_debug=True)
        messagebox.showinfo("Log Saved", f"Log saved to:\n{filename}")

# Fancy status bar but hasn't been wired to all places
def set_main_status(text, fg=None):
    status_var.set(text)
    if fg is not None:
        try:
            status_label.config(fg=fg)
        except Exception:
            pass

    
root.protocol("WM_DELETE_WINDOW", on_closing)
root.title("GlycoMSP File Manager GUI v1.10 Build 20260521 --finalcheck-CGA")
root.geometry("800x560")
root.minsize(800, 560)

# Buttons
button_frame = tk.Frame(root)
button_frame.pack(pady=10)

tk.Button(button_frame, text="Add Raw File", command=lambda: select_file("raw")).grid(row=0, column=0, padx=5)
tk.Button(button_frame, text="Add mzML File", command=lambda: select_file("mzml")).grid(row=0, column=1, padx=5)
tk.Button(button_frame, text="Clear", command=clear_files).grid(row=0, column=2, padx=5)
tk.Button(button_frame, text="About", command=open_about_window).grid(row=0, column=3, padx=5)
tk.Button(button_frame, text="Save Log", command=save_log_to_file).grid(row=0, column=4, padx=5)
#tk.Button(root, text="Save Log", command=save_log_to_file).pack(pady=5)
# Text widget to log selected files
text_widget = tk.Text(root, height=20, width=80)
text_widget.pack(pady=10)

# Status bar
status_var = tk.StringVar()
status_var.set("Idle")
status_label = tk.Label(root, textvariable=status_var, fg="blue")
status_label.pack(pady=5)

progress = ttk.Progressbar(root, orient="horizontal", mode="indeterminate", length=250)
progress.pack(pady=5)


# --- Analysis Tools Frame ---
analysis_frame = tk.LabelFrame(root, text="Analysis Tools", padx=10, pady=10)
analysis_frame.pack(padx=10, pady=10, fill="x")
convert_button = tk.Button(analysis_frame, text="Convert Raw to CSV",command=launch_metadata_batch)
convert_button.pack(side="left", padx=5)
tk.Button(analysis_frame, text="Prepare Dataset", command=open_prepare_dataset_window).pack(side="left", padx=5)
tk.Button(analysis_frame, text="Run ML Analysis", command=open_ml_analysis_window).pack(side="left", padx=5)


#import safe
if __name__ == "__main__":
    root.mainloop()
