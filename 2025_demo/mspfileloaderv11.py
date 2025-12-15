import os
version = "0.9999"
last_update = 20251208
import msprawextractor as mspext
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
import shutil
import traceback
import mspvalidator_merger as mspval
import pandas as pd
import platform
import msp_insilicomarker_withGPT as marker
# --- Add near top-level imports ---
import os, joblib
import compnewv4 as compv4  # assumes dev/test calls are guarded by if __name__ == "__main__"
# Ion mining feature (safe to import even if file is absent)
# 20250907 to solve ion suggest missing issue
import importlib, sys, os, traceback
from pathlib import Path
from pathcanon import to_posix_str, to_native_path, ensure_dir
from typing import Dict


# ---- GlycoMSP startup diagnostics: put this block at line 1 ----
import traceback, importlib.util

def _pause_if_no_tty():
    # If user double-clicked on Windows, stdin often isn't a real TTY.
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

def _quick_diag(note=None):
    print("=== GlycoMSP quick environment check ===")
    print(f"Python: {sys.version.split()[0]}  ({sys.executable})")
    print(f"Platform: {platform.platform()}")
    if note:
        print(note)

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
        print(f"{tag} {shown:<10} {info}")

    if os.name == "nt":
        ok, info = _check_mod("pymsfilereader")
        tag = "OK  " if ok else "MISS"
        print(f"{tag} pymsfilereader {info}  (needed only for RAW conversion)")

    if sys.version_info >= (3, 13):
        print("\n[Notice] Running on Python 3.13.",
              "If wheels for numpy/pandas/scikit-learn are missing, prefer Python 3.11/3.12 for now.")

def _excepthook(exctype, value, tb):
    print("\n[StartupError] Unhandled exception during launch/import:\n")
    traceback.print_exception(exctype, value, tb)
    print()
    _quick_diag()
    _pause_if_no_tty()
    # Swallow default hook to avoid duplicate tracebacks
sys.excepthook = _excepthook

# Allow a no-import diagnostics mode:
if "--diag" in sys.argv or "--diag" in sys.argv:
    _quick_diag(note="(ran with --diag; skipped importing the full GUI)")
    _pause_if_no_tty()
    sys.exit(0)

# ---- End diagnostics prelude ----


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

"""
try:
    from msp_ion_mining import export_ion_suggestions_csv, SuggestParams
    print("[dev] loaded msp_ion_mining")
except Exception:
    export_ion_suggestions_csv = None
    class SuggestParams:  # fallback stub
        def __init__(self, **kw): pass
"""
# v1.01? (future) fix the old macos crash issue due to malformed tkinter askopenfilename (see crash report analysis in GPT chat)
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

try: 
    import tkinter as tk
    from tkinter import filedialog, messagebox
except ImportError:
    raise ImportError('Please install tkinter to enable GUI-based file selection')
else:
    has_tkinter = True

selected_files = {}  # Dictionary to keep track of selected files


#probably same as pathcanon?
def _to_posix(p):
    return None if not p else Path(p).as_posix()

def _from_posix(p):
    return None if not p else os.path.normpath(p.replace('/', os.sep))

def _norm_for_display(p):
    return os.path.normpath(p) if p else p

#20250907 preventing key error in GUI #need real test to see behavior changes
FILETYPE_TO_KEY = {"csv": "csv", "excel": "excel", "json": "json",
                   "method": "json", "metadata": "json"}  # extend if you add new labels
def normalize_ftype(ft: str) -> str:
    return FILETYPE_TO_KEY.get(ft.lower().strip(), ft.lower().strip())
#
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


#logger
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

#20250822 replace composition window with pseudolabel window
# --- Pseudo-Labeling Setup window (replaces the old GlycanCompositionWindow) ---

class PseudoLabelingSetupWindow(tk.Toplevel):
    """
    Edit in-silico generation flags and preview essential metadata.
    - Defaults injected via `default_flags` (dict copied before editing)
    - Read-only metadata summary (Glycan Type / Charge / Derivatization)
    - Save/Load flags presets (JSON)
    - Calls `on_submit({"flags": ..., "metadata": {...}})` on Generate
    """
    def __init__(self, master,
             meta_json_path=None, meta_prefill=None, default_flags=None, on_submit=None,
             on_generate=None, on_link_existing=None, on_attach_ionlist=None, on_start=None,
             initial_insilico=None, initial_ionlist=None):

        super().__init__(master)
        self.title("Constraint-based Glycan Annotation Setup")
        self.geometry("620x640")
        self.resizable(True, True)

        self.meta_json_path = meta_json_path
        self.on_submit = on_submit
        self.flags = dict(default_flags or {})  # shallow copy; safe to mutate locally
        self.on_generate = on_generate
        self.on_link_existing = on_link_existing
        self.on_attach_ionlist = on_attach_ionlist
        self.on_start = on_start
        self.insilico_path_var = tk.StringVar(value=initial_insilico or "")
        self.ionlist_path_var  = tk.StringVar(value=initial_ionlist or "")


        """
        # --- Metadata summary (read-only) ---
        meta_frame = ttk.LabelFrame(self, text="Metadata Summary")
        meta_frame.pack(fill="x", padx=12, pady=(12, 6))
        self.meta_vars = {
            "Glycan Type": tk.StringVar(value=""),
            "Mass Analyzer charge mode": tk.StringVar(value=""),
            "Derivatization Type": tk.StringVar(value=""),
        }
        #debug lines
        print("[meta] GUI vars:", {k: v.get() for k, v in self.meta_vars.items()})
        # 1) prefill from dict if provided
        if meta_prefill:
            self.meta_vars["Glycan Type"].set(meta_prefill.get("Glycan Type", ""))
            self.meta_vars["Mass Analyzer charge mode"].set(meta_prefill.get("Mass Analyzer charge mode", ""))
            self.meta_vars["Derivatization Type"].set(meta_prefill.get("Derivatization Type", ""))

        # 2) then, if a path is given, try to read/overwrite from file
        if self.meta_json_path and os.path.exists(self.meta_json_path):
            try:
                with open(self.meta_json_path, "r", encoding="utf-8") as f:
                    meta = json.load(f)
                # only set if key exists, to avoid wiping prefill with blanks
                if "Glycan Type" in meta:
                    self.meta_vars["Glycan Type"].set(meta["Glycan Type"])
                if "Mass Analyzer charge mode" in meta:
                    self.meta_vars["Mass Analyzer charge mode"].set(meta["Mass Analyzer charge mode"])
                if "Derivatization Type" in meta:
                    self.meta_vars["Derivatization Type"].set(meta["Derivatization Type"])
            except Exception as e:
                messagebox.showwarning("Metadata", f"Could not read metadata:\n{self.meta_json_path}\n\n{e}")
        r = 0
        for label, var in self.meta_vars.items():
            ttk.Label(meta_frame, text=f"{label}:").grid(row=r, column=0, sticky="w", padx=10, pady=4)
            ttk.Label(meta_frame, textvariable=var).grid(row=r, column=1, sticky="w", padx=8, pady=4)
            r += 1
        """
        # --- Metadata (editable) ---
        meta_frame = ttk.LabelFrame(self, text="Metadata (can override here)")
        meta_frame.pack(fill="x", padx=12, pady=(12, 6))

        # choices (adjust if you have more)
        GLYCAN_CHOICES = ["N", "O"]
        CHARGE_CHOICES = ["+", "-"]
        DERIV_CHOICES  = ["PerMe(Reduced)", "PerMe(Freeend)", "None"]

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

        # render editable controls
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
        """
        # 2) then, if a path is given, try to read/overwrite from file (only if keys exist)
        file_meta = None
        if self.meta_json_path and os.path.exists(self.meta_json_path):
            try:
                with open(self.meta_json_path, "r", encoding="utf-8") as f:
                    file_meta = json.load(f)
                if "Glycan Type" in file_meta:
                    self.meta_vars["Glycan Type"].set(file_meta["Glycan Type"])
                if "Mass Analyzer charge mode" in file_meta:
                    self.meta_vars["Mass Analyzer charge mode"].set(file_meta["Mass Analyzer charge mode"])
                if "Derivatization Type" in file_meta:
                    self.meta_vars["Derivatization Type"].set(file_meta["Derivatization Type"])
            except Exception as e:
                messagebox.showwarning("Metadata", f"Could not read metadata:\n{self.meta_json_path}\n\n{e}"
        """
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
            #
            "arm_count","internal_minrep","internal_maxrep","topology",
            "compcheck","Hex_range","HexNAc_range","Neu5Ac_range","Neu5Gc_range","KDN_range","Fucose_range"
        }
        self.flag_vars = {}
        #newly added
        # Keys that belong to the "Advanced options" collapsible group
        ADVANCED_KEYS = {"debug", "dev", "force_exit", "topology", "perman"}
        # controls visibility of "Advanced options" panel
        self.show_advanced_var = tk.BooleanVar(value=False)
        

        left = ttk.Frame(flags_frame)
        right = ttk.Frame(flags_frame)
        #left.pack(side="left", fill="both", expand=True, padx=(10, 5), pady=8)
        #right.pack(side="left", fill="both", expand=True, padx=(5, 10), pady=8)
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
            ttk.Label(parent, text=key + ":").grid(row=row, column=0, sticky="w")
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
        
        """
        def _snapshot_flags():
            snap = {}
            # flags
            for k, w in self.flag_vars.items():
                if isinstance(w, tuple):
                    snap[k] = (int(w[0].get()), int(w[1].get()))
                elif isinstance(w, tk.BooleanVar):
                    snap[k] = bool(w.get())
                else:
                    snap[k] = int(w.get())
            # OG core types (added below)
            snap["_ogcore"] = {i: v.get() for i, v in self.og_core_vars.items()}
            return snap
        """
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
            for child in left.winfo_children(): child.destroy()
            for child in right.winfo_children(): child.destroy()
            for child in advanced_frame.winfo_children(): child.destroy()
            self.flag_vars.clear()
            # choose keys
            keys = [k for k in flag_spec.keys() if k in _current_keys()]
            main_keys = [k for k in keys if k not in ADVANCED_KEYS]
            adv_keys = [k for k in keys if k in ADVANCED_KEYS]
            half = (len(main_keys) + 1) // 2
            left_keys, right_keys = main_keys[:half], main_keys[half:]
            # (re)render
            def render_column(parent, keys_subset):
                r = 0
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
            # main columns
            render_column(left, left_keys)
            render_column(right, right_keys)
            # advanced panel (optional)
            if self.show_advanced_var.get() and adv_keys:
                if not advanced_frame.winfo_manager():
                    advanced_frame.pack(fill="x", padx=10, pady=(0, 8))
                render_column(advanced_frame, adv_keys)
            else:
                advanced_frame.pack_forget()
            # show/hide O-core panel according to glycan type
            _toggle_og_core_panel()
            # restore state
            _restore_flags(snap)   
        """
        def _rebuild_flag_panel(*_):
            # preserve state
            snap = _snapshot_flags()
            # clear frames
            for child in left.winfo_children(): child.destroy()
            for child in right.winfo_children(): child.destroy()
            self.flag_vars.clear()
            # choose keys
            keys = [k for k in flag_spec.keys() if k in _current_keys()]
            half = (len(keys) + 1) // 2
            left_keys, right_keys = keys[:half], keys[half:]
            # (re)render
            def render_column(parent, keys_subset):
                r = 0
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
            render_column(left, left_keys)
            render_column(right, right_keys)
            # show/hide O-core panel according to glycan type
            _toggle_og_core_panel()
            # restore state
            _restore_flags(snap)   
        """
        # --- O-glycan core types (multi-select) ---
        #core_frame = ttk.LabelFrame(self, text="O-glycan core types (select 1–4)")
        #core_frame.pack(fill="x", padx=12, pady=(0, 6))
        # Now a sub-frame inside "In-Silico Generation Flags"
        core_frame = ttk.LabelFrame(flags_frame, text="O-glycan core types (select 1–4)")

        self.og_core_vars = {i: tk.BooleanVar(value=False) for i in (0, 1, 2, 3, 4)}

        row = 0
        for i, label in [(0, "Tn Antigen"), (1, "Core 1"), (2, "Core 2"), (3, "Core 3"), (4, "Core 4")]:
            ttk.Checkbutton(core_frame, text=label, variable=self.og_core_vars[i]).grid(
                row=row // 2, column=row % 2, sticky="w", padx=8, pady=3
            )
            row += 1

        def _selected_coretypes():
            sel = [i for i, v in self.og_core_vars.items() if v.get()]
            return sel if sel else [1, 2, 3, 4]  # sensible default

        def _toggle_og_core_panel():
            # Show only for O-glycan
            gtype = self.meta_vars["Glycan Type"].get().strip().upper()
            core_frame.pack_forget()
            if gtype == "O":
                # if outside the block then use this old one
                #core_frame.pack(fill="x", padx=12, pady=(0, 6))
                # show inside the "In-Silico Generation Flags" block
                core_frame.pack(fill="x", padx=10, pady=(4, 6))

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

        """
        # --- Flag editor (selected subset; excludes termi_comp/internal_comp) ---
        flags_frame = ttk.LabelFrame(self, text="In-Silico Generation Flags")
        flags_frame.pack(fill="both", expand=True, padx=12, pady=6)

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
            # iteration logic
            "arm_count": ("int", (0, 8)),
            "internal_minrep": ("int", (0, 6)),
            "internal_maxrep": ("int", (0, 10)),
            "topology": ("bool", None),
            # core flags (NG)
            "corefuc": ("bool", None),
            "bicorefuc": ("bool", None),
            "highman": ("bool", None),
            "perman": ("bool", None),
            "hybrid": ("bool", None),  # note: hybrid components are derived later, not edited here
            # optional composition check
            "compcheck": ("bool", None),
            "Hex_range": ("range", (0, 20)),
            "HexNAc_range": ("range", (0, 20)),
            "Neu5Ac_range": ("range", (0, 10)),
            "Neu5Gc_range": ("range", (0, 10)),
            "KDN_range": ("range", (0, 10)),
            "Fucose_range": ("range", (0, 10)),
        }

        self.flag_vars = {}

        left = ttk.Frame(flags_frame)
        right = ttk.Frame(flags_frame)
        left.pack(side="left", fill="both", expand=True, padx=(10, 5), pady=8)
        right.pack(side="left", fill="both", expand=True, padx=(5, 10), pady=8)

        def _coerce_int_like(val, default=0):
            return int(val) if isinstance(val, (int, float, str)) and str(val).strip() != "" else int(default)

        def add_bool(parent, key, row):
            var = tk.BooleanVar(value=bool(self.flags.get(key, False)))
            ttk.Checkbutton(parent, text=key, variable=var).grid(row=row, column=0, sticky="w", pady=3)
            self.flag_vars[key] = var
            if key == "hybrid":
                ttk.Label(parent, text="(components auto-filled on submit)").grid(row=row, column=1, sticky="w")

        def add_int(parent, key, row, lo, hi):
            ttk.Label(parent, text=key + ":").grid(row=row, column=0, sticky="w")
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

        keys = list(flag_spec.keys())
        half = (len(keys) + 1) // 2
        left_keys, right_keys = keys[:half], keys[half:]

        def render_column(parent, keys_subset):
            r = 0
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

        render_column(left, left_keys)
        render_column(right, right_keys)
        """
        """
        btns = ttk.Frame(self)
        btns.pack(fill="x", padx=12, pady=(6, 12))
        ttk.Button(btns, text="Load Flags…", command=self.load_flags).pack(side="left", padx=4)
        ttk.Button(btns, text="Save Flags…", command=self.save_flags).pack(side="left", padx=4)
        ttk.Separator(btns, orient="vertical").pack(side="left", fill="y", padx=8)
        ttk.Button(btns, text="Generate In-Silico CSV", command=self.submit).pack(side="left", padx=4)
        ttk.Button(btns, text="Cancel", command=self.destroy).pack(side="right", padx=4)
        """
        # --- BELOW your metadata + flags UI ---

        # Small panel to show current links (insilico / ion list)
        links_frame = ttk.LabelFrame(self, text="Files selected for CGA analysis")
        links_frame.pack(fill="x", padx=12, pady=(6, 0))

        self.insilico_path_var = tk.StringVar(value=initial_insilico or "")
        self.ionlist_path_var  = tk.StringVar(value=initial_ionlist or "")

        ttk.Label(links_frame, text="In-silico glycan list (compositions):").grid(row=0, column=0, sticky="w", padx=10, pady=4)
        ttk.Label(links_frame, textvariable=self.insilico_path_var).grid(row=0, column=1, sticky="w", padx=8, pady=4)

        ttk.Label(links_frame, text="Fragmentation ion list (features):").grid(row=1, column=0, sticky="w", padx=10, pady=4)
        ttk.Label(links_frame, textvariable=self.ionlist_path_var).grid(row=1, column=1, sticky="w", padx=8, pady=4)

        # --- Footer buttons
        btns = ttk.Frame(self)
        btns.pack(fill="x", padx=12, pady=(6, 12))

        ttk.Button(btns, text="Load CGA settings", command=self.load_flags).pack(side="left", padx=4)
        ttk.Button(btns, text="Save CGA settings", command=self.save_flags).pack(side="left", padx=4)

        ttk.Separator(btns, orient="vertical").pack(side="left", fill="y", padx=8)
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

        """
        # Generate in-silico: call the provided on_generate/on_submit callback
        def _on_generate():
            payload = {"flags": self.collect_flags(),
                    "metadata": {
                        "Glycan Type": self.meta_vars["Glycan Type"].get(),
                        "Mass Analyzer charge mode": self.meta_vars["Mass Analyzer charge mode"].get(),
                        "Derivatization Type": self.meta_vars["Derivatization Type"].get(),
                        "_meta_json": self.meta_json_path or "",
                        "_overrides_applied": any(self.meta_vars[k].get() != self._meta_original.get(k, "")
                                                    for k in self._meta_original)
                    }}
            # Back-compat: accept either on_submit or on_generate
            cb = self.on_submit or self.on_generate
            if cb:
                # If the callback returns a path, reflect it in the UI
                maybe_path = cb(payload)
                if isinstance(maybe_path, str) and os.path.exists(maybe_path):
                    self.insilico_path_var.set(maybe_path)
        """
        ttk.Button(btns, text="Generate In-Silico CSV", command=_on_generate).pack(side="left", padx=4)

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

        ttk.Button(btns, text="Link in silico glycan list", command=_link_existing).pack(side="left", padx=4)

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

        ttk.Button(btns, text="Add Ion List", command=_attach_ionlist).pack(side="left", padx=4)

        ttk.Separator(btns, orient="vertical").pack(side="left", fill="y", padx=8)

        # Start pseudolabeling (enabled if converted CSV + insilico present)
        def _start():
            payload = {"flags": self.collect_flags(),
                    "metadata": {
                        "Glycan Type": self.meta_vars["Glycan Type"].get(),
                        "Mass Analyzer charge mode": self.meta_vars["Mass Analyzer charge mode"].get(),
                        "Derivatization Type": self.meta_vars["Derivatization Type"].get(),
                        "_meta_json": self.meta_json_path or "",
                        "_overrides_applied": any(self.meta_vars[k].get() != self._meta_original.get(k, "")
                                                    for k in self._meta_original)
                    },
                    "insilico_csv": self.insilico_path_var.get().strip(),
                    "ionlist_path": self.ionlist_path_var.get().strip(),
                    "coretype": _selected_coretypes()}  # <-- NEW
            cb = self.on_start
            if cb:
                cb(payload)

        ttk.Button(btns, text="Start CGA analysis", command=_start).pack(side="left", padx=4)

        ttk.Button(btns, text="Close", command=self.destroy).pack(side="right", padx=4)


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

    def submit(self):
        flags = self.collect_flags()

        # collect editable metadata (overrides)
        meta = {
            "Glycan Type": self.meta_vars["Glycan Type"].get(),
            "Mass Analyzer charge mode": self.meta_vars["Mass Analyzer charge mode"].get(),
            "Derivatization Type": self.meta_vars["Derivatization Type"].get(),
            "_meta_json": self.meta_json_path or "",
        }
        # mark whether anything was changed relative to the file/original values
        meta_changed = any(self.meta_vars[k].get() != self._meta_original.get(k, "") for k in self._meta_original)
        meta["_overrides_applied"] = bool(meta_changed)

        if self.on_submit:
            self.on_submit({"flags": flags, "metadata": meta})
        self.destroy()

#metadata class 

class MetadataEditorWindow:
    def __init__(self, parent, raw_file_list, on_each_metadata_ready_callback, on_finish=None, output_dir=None, skip_conversion=False):
        self.parent = parent
        self.raw_file_list = raw_file_list
        self.callback = on_each_metadata_ready_callback
        self.on_finish = on_finish
        self.current_index = 0
        self.last_metadata = None
        self.metadata = {}
        self.output_dir = output_dir
        self.skip_conversion = skip_conversion
        self.entries = {}
        self.files = {}   # ← add this line in 20250916
        fields = [
            "Experiment Title", "Experiment Description", "Author running this analysis",
            "Raw data acquired date", "Glycan Type", "Mass Analyzer charge mode", "Derivatization Type"
        ]
        #20250927
        self._tmp_paths = None

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

        if not self.skip_conversion:
            self.load_file(self.raw_file_list[self.current_index])

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

        # Buttons
        button_frame = tk.Frame(self.window)
        button_row = len(fields) + 1 #so the button adjusted itself in future if we add more contents into metadata
        button_frame.grid(row=button_row, column=0, columnspan=2, pady=10)
        tk.Button(button_frame, text="Import Metadata", command=self.import_metadata).grid(row=0, column=0, padx=5)
        tk.Button(button_frame, text="Use Last Metadata", command=self.use_last_metadata).grid(row=0, column=1, padx=5)
        tk.Button(button_frame, text="Clear", command=self.clear_fields).grid(row=0, column=2, padx=5)
        tk.Button(button_frame, text="Generate Dataset", command=self.generate).grid(row=0, column=3, padx=5)
        tk.Button(button_frame, text="Cancel", command=self.window.destroy).grid(row=0, column=4, padx=5)


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

        # Extract filename base (your existing line)
        self.rawfilename = os.path.splitext(os.path.basename(self.current_raw_file))[0]

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

        # Add status label if not already there
        if not hasattr(self, "status_label"):
            self.status_label = tk.Label(self.window, text="", fg="blue")
            self.status_label.grid(row=0, columnspan=2, pady=10)
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

                logger.log(f"[convert] Temp CSVs: {tmp_ms2}, {tmp_ms3}")

                self._convert_in_progress = False
                self._safe_after(0, self.on_conversion_complete)

            except Exception as e:
                self._convert_in_progress = False
                logger.log(f"[FATAL] Exception during raw file conversion: {e}")
                self._enter_error_mode("An error occurred. Please check the log.")

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
        """
        def background_conversion():
            if not self.skip_conversion:
                self._convert_in_progress = True  # NEW: mark running
                try:
                    tmp_ms2, tmp_ms3 = mspext.convert_raw_to_csv(
                        raw_file,
                        outdir=self.output_dir,
                        debug=True
                    )
                    # Remember where the temps are for finalize()
                    self.files["ms2tmp"] = tmp_ms2
                    self.files["ms3tmp"] = tmp_ms3

                    logger.log(f"[convert] Temp CSVs: {tmp_ms2}, {tmp_ms3}")

                    # mark done BEFORE we notify UI
                    self._convert_in_progress = False   # NEW
                    # notify UI safely (toplevel may be recreated/closed)
                    self._safe_after(0, self.on_conversion_complete)

                except Exception as e:
                    # mark done so UI can decide (batch may still have pending files)
                    self._convert_in_progress = False   # NEW
                    logger.log(f"[FATAL] Exception during raw file conversion: {e}")
                    #self._safe_set_status("Exception occurred. Please check the log.", fg="red")
                    self._enter_error_mode("An error occurred. Please check the log.")
                    # best-effort cleanup
                    try:
                        for p in [self.files.get("ms2tmp"), self.files.get("ms3tmp")]:
                            if p and os.path.exists(p):
                                os.remove(p)
                                logger.log(f"[Cleanup] Removed leftover temp file: {p}")
                    except Exception as ce:
                        logger.log(f"[WARNING] Cleanup failed: {ce}")

                    # show error (capture 'e' into lambda safely)
                    self._safe_after(
                        0,
                        lambda err=e: messagebox.showerror(
                            "Thermo Library Error",
                            f"Raw file could not be processed.\n\nDetails:\n{err}"
                        )
                    )

                    # update MAIN status to Idle (this file failed)
                    #self._safe_after(0, lambda: set_main_status("Idle"))
                    self._safe_after(0, self._update_main_status_for_batch_state)
                    # defer window close decision (safe close if not batching)
                    self._safe_after(500, self._attempt_close_after_error)  # NEW
                #20250917 what does finally do?
                finally:
                    # ensure flags are sane even if we bailed early
                    self._convert_in_progress = False
                    self._safe_after(0, self._update_main_status_for_batch_state)
            
            if not self.skip_conversion:
                try:
                    tmp_ms2, tmp_ms3 = mspext.convert_raw_to_csv(
                        raw_file,
                        outdir=self.output_dir,
                        debug=True
                    )
                    # Remember where the temps are for finalize()
                    self.files["ms2tmp"] = tmp_ms2
                    self.files["ms3tmp"] = tmp_ms3

                    # (Optional) sanity log
                    logger.log(f"[convert] Temp CSVs: {tmp_ms2}, {tmp_ms3}")

                    # notify UI
                    #self.window.after(0, self.on_conversion_complete)
                    self._safe_after(0, self.on_conversion_complete) #20250917 patch

                except Exception as e:
                    logger.log(f"[FATAL] Exception during raw file conversion: {e}")
                    self.status_label.config(text="Exception occurred. Please check the log.", fg="red")
                    # best-effort cleanup
                    try:
                        for p in [self.files.get("ms2tmp"), self.files.get("ms3tmp")]:
                            if p and os.path.exists(p):
                                os.remove(p)
                                logger.log(f"[Cleanup] Removed leftover temp file: {p}")
                    except Exception as ce:
                        logger.log(f"[WARNING] Cleanup failed: {ce}")

                    # show error (capture e into lambda)
                    self.window.after(
                        0,
                        lambda err=e: messagebox.showerror(
                            "Thermo Library Error",
                            f"Raw file could not be processed.\n\nDetails:\n{err}"
                        )
                    )
                    # 🔑 ensure status bar is reset
                    self.window.after(0, lambda: self.set_status("Idle"))    #this line is wrong. Set root status to Idle   
                    # destroy the window
                    self.window.destroy()    
            
            if not self.skip_conversion:
                try:
                    success = mspext.convert_raw_to_csv(raw_file, outdir=self.output_dir, debug=True)
                    if success:
                        self.window.after(0, self.on_conversion_complete)
                    else:
                        # Clean up temp files
                        rawfilename = os.path.splitext(os.path.basename(raw_file))[0]
                        temp_ms2 = f"ms2tmp_{rawfilename}.csv"
                        temp_ms3 = f"ms3tmp_{rawfilename}.csv"

                        for temp_file in [temp_ms2, temp_ms3]:
                            if os.path.exists(temp_file):
                                try:
                                    os.remove(temp_file)
                                    logger.log(f"[Cleanup] Removed leftover temp file: {temp_file}")
                                except Exception as e:
                                    logger.log(f"[WARNING] Failed to delete temp file: {temp_file} — {e}")
                        err_msg = str(success) if isinstance(success, Exception) else "Unknown error"
                        self.window.after(0, lambda: messagebox.showerror("Conversion Failed", f"The raw file could not be converted:\n{err_msg}"))
                        traceback.print_exc()
                except Exception as e:
                    logger.log(f"[FATAL] Exception during raw file conversion: {e}")
                    self.window.after(0, lambda: messagebox.showerror("Thermo Library or conversion Error",f"Raw file could not be processed.\n\nDetails:\n{e}"))
        """        
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
    
    #20250917
    def _widget_alive(self, w):
        try:
            return (w is not None) and w.winfo_exists()
        except Exception:
            return False

    def _safe_set_status(self, text, fg=None):
        try:
            if self._widget_alive(getattr(self, "status_label", None)):
                if fg is None:
                    self.status_label.config(text=text)
                else:
                    self.status_label.config(text=text, fg=fg)
        except Exception as e:
            logger.log(f"[UI] status update skipped: {e}")

    def _safe_after(self, ms, func, *args, **kwargs):
        host = getattr(self, "parent", None) or getattr(self, "window", None)
        if self._widget_alive(host):
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
            # main window status → batch still running
            set_main_status("Batch conversion running", fg="blue")
            # metadata window status + warning
            self._safe_set_status("Batch is running — cannot close this window now.", fg="orange")
            messagebox.showwarning(
                "Batch Running",
                "Batch conversion is still in progress. Please wait for this file to finish or cancel the batch."
            )
            try:
                self.window.lift()
                self.window.focus_force()
            except Exception:
                pass
            return  # block close
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
            if self._widget_alive(getattr(self, "window", None)):
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
            if self._widget_alive(getattr(self, "window", None)):
                try:
                    self.window.destroy()
                except Exception:
                    pass        

    def _update_main_status_for_batch_state(self):
        """Set main status based on whether batch is still active."""
        try:
            more_pending = self.current_index < (len(self.raw_file_list) - 1)
        except Exception:
            more_pending = False
        if getattr(self, "_convert_in_progress", False) or more_pending:
            set_main_status("Batch conversion running", fg="blue")
        else:
            set_main_status("Idle", fg="blue")

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

        extractdate = time.strftime("%Y%m%d")
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

                # fallback: construct likely locations
                raw_stem = os.path.splitext(os.path.basename(raw_path))[0] if raw_path != "not linked" else "no_raw"
                if not temp_ms2:
                    cand = os.path.join(self.output_dir, f"ms2tmp_{raw_stem}.csv")
                    temp_ms2 = cand if os.path.exists(cand) else os.path.abspath(f"ms2tmp_{raw_stem}.csv")
                if not temp_ms3:
                    cand = os.path.join(self.output_dir, f"ms3tmp_{raw_stem}.csv")
                    temp_ms3 = cand if os.path.exists(cand) else os.path.abspath(f"ms3tmp_{raw_stem}.csv")

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

    #20250917 ends
    """
    def generate(self):
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

        for field, entry in self.entries.items():
            self.metadata[field] = entry.get()

        extractdate = time.strftime("%Y%m%d")
        self.metadata["Parameters when GlycoMSP launched"] = ["Autofilled by GUI"]
        self.metadata["Date of file extracted from raw file"] = extractdate
        # --- Safely resolve raw file path ---
        raw_path = None
        # Priority: use field if present (linking mode)
        if hasattr(self, "raw_file_entry"):
            raw_path = self.raw_file_entry.get().strip()
        # Fallback: use internal attribute if available (conversion mode)
        if not raw_path:
            raw_path = getattr(self, "current_raw_file", None)
        # Default to not linked
        raw_path = raw_path or "not linked"
        # Save into metadata
        self.metadata["Original raw file path"] = raw_path
        self.metadata["Raw filename"] = os.path.splitext(os.path.basename(raw_path))[0] if raw_path != "not linked" else "(not linked)"

        mansavename = simpledialog.askstring("Project Name", "Enter a name for this batch/project:")
        if not mansavename:
            messagebox.showwarning("Missing Name", "You must enter a project/batch name.")
            return

        raw_base = os.path.splitext(os.path.basename(raw_path))[0] if raw_path != "not linked" else "no_raw"
        savename = f"{mansavename}_{extractdate}_{raw_base}"

        output_json_path = os.path.join(self.output_dir, savename + ".json")
        #if not self.skip_conversion:
        jsonfile = mspext.savemetadata(self.metadata, output_json_path)
        #else:
        #    print(f"[debug] not using mspext so no error happens if one tries to save metadata. Will metadata be saved?")
        logger.log(f"Metadata saved to {jsonfile}")
        
        # Store last metadata
        self.last_metadata = self.metadata.copy()


        #temp_ms2 = f"ms2tmp_{self.rawfilename}.csv"
        raw_stem = os.path.splitext(os.path.basename(raw_path))[0] if raw_path != "not linked" else "no_raw"
        final_ms2 = os.path.join(self.output_dir, f"ms2_{savename}.csv")
        final_ms3 = os.path.join(self.output_dir, f"ms3_{savename}.csv")
        success = True

        if self.skip_conversion:
            print("[debug] filling missing metadata only; no csv rename/check")
        else:
            try:
                Path(final_ms2).parent.mkdir(parents=True, exist_ok=True)

                # Prefer the exact temp paths from background_conversion
                # prefer paths captured during background_conversion
                temp_ms2 = self.files.get("ms2tmp")
                temp_ms3 = self.files.get("ms3tmp")
                #temp_ms2 = self.files.get("ms2tmp") or os.path.join(self.output_dir, f"ms2tmp_{raw_stem}.csv")
                #temp_ms3 = self.files.get("ms3tmp") or os.path.join(self.output_dir, f"ms3tmp_{raw_stem}.csv")

                raw_stem = os.path.splitext(os.path.basename(raw_path))[0] if raw_path != "not linked" else "no_raw"
                if not temp_ms2:
                    temp_ms2 = os.path.join(self.output_dir, f"ms2tmp_{raw_stem}.csv")
                if not temp_ms3:
                    temp_ms3 = os.path.join(self.output_dir, f"ms3tmp_{raw_stem}.csv")

                final_ms2 = os.path.join(self.output_dir, f"ms2_{savename}.csv")
                final_ms3 = os.path.join(self.output_dir, f"ms3_{savename}.csv")

                success = True
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
                
                if temp_ms2 and os.path.exists(temp_ms2):
                    os.replace(temp_ms2, final_ms2)  # atomic on same filesystem
                    logger.log(f"[finalize] Promoted {temp_ms2} → {final_ms2}")
                else:
                    logger.log(f"[WARNING] Missing temp MS2 file: {temp_ms2}")
                    success = False

                if temp_ms3 and os.path.exists(temp_ms3):
                    os.replace(temp_ms3, final_ms3)  # atomic on same filesystem
                    logger.log(f"[finalize] Promoted {temp_ms3} → {final_ms3}")
                else:
                    logger.log(f"[WARNING] Missing temp MS3 file: {temp_ms3}")
                    success = False

                # Update state with final paths
                if success:
                    self.files["ms2"] = final_ms2
                    self.files["ms3"] = final_ms3
                    self.files.pop("ms2tmp", None)
                    self.files.pop("ms3tmp", None)
                
            except Exception as e:
                logger.log(f"[FATAL] Finalize failed: {e}")
                success = False

            if success:
                mspext.finalize_extraction(self.current_raw_file, self.metadata, savename)
                logname = logger.save(os.path.join(self.output_dir, savename + ".log"), include_debug=True)
                logger.log(f"Saved log to: {logname}")

            self.current_index += 1
            if self.current_index < len(self.raw_file_list):
                self.load_file(self.raw_file_list[self.current_index])
            else:
                messagebox.showinfo("Metadata", "All metadata have been completed.")
                self.window.destroy()
                if self.on_finish:
                    self.on_finish()
            
            
        # Trigger callback after metadata is created
        if self.callback:
            self.callback(
                raw_path,             # resolved raw path
                self.metadata,        # metadata dict
                savename              # base name for saved file
            )
        """
        # Auto-close the window
        #self.window.destroy() #remove this to avoid the crash when dealing with batch processing

    """
        raw_stem = os.path.splitext(os.path.basename(raw_path))[0] if raw_path != "not linked" else "no_raw"
        temp_ms2 = f"ms2tmp_{raw_stem}.csv"
        temp_ms3 = f"ms3tmp_{raw_stem}.csv"
        final_ms2 = os.path.join(self.output_dir, f"ms2_{savename}.csv")
        final_ms3 = os.path.join(self.output_dir, f"ms3_{savename}.csv")
        success = True
        if self.skip_conversion:
            print(f"[debug]: probably filling missing metadata. No need to do csv rename and check")
            #pass
        else:
            try:
                Path(final_ms2).parent.mkdir(parents=True, exist_ok=True)

                if os.path.exists(temp_ms2):
                    os.replace(temp_ms2, final_ms2)  # atomic on same filesystem
                    logger.log(f"[finalize] Promoted {temp_ms2} → {final_ms2}")
                else:
                    logger.log(f"[WARNING] Missing temp MS2 file: {temp_ms2}")
                    success = False

                if os.path.exists(temp_ms3):
                    os.replace(temp_ms3, final_ms3)  # atomic on same filesystem
                    logger.log(f"[finalize] Promoted {temp_ms3} → {final_ms3}")
                else:
                    logger.log(f"[WARNING] Missing temp MS3 file: {temp_ms3}")
                    success = False
                
                if os.path.exists(temp_ms2):
                    shutil.move(temp_ms2, final_ms2)
                    logger.log(f"Renamed {temp_ms2} → {final_ms2} and move to {self.output_dir}")
                else:
                    logger.log(f"[WARNING] Missing temp MS2 file: {temp_ms2}")
                    success = False

                if os.path.exists(temp_ms3):
                    shutil.move(temp_ms3, final_ms3)
                    logger.log(f"Renamed {temp_ms3} → {final_ms3} and move to {self.output_dir}")
                else:
                    logger.log(f"[WARNING] Missing temp MS3 file: {temp_ms3}")
                    success = False
                
            except Exception as e:
                logger.log(f"[ERROR] Failed during file renaming: {e}")
                success = False

            except Exception as e:
                logger.log(f"[ERROR] Failed during file renaming: {e}")
                success = False
        """

def write_to_gui(message):
    text_widget.insert(tk.END, message + "\n")
    text_widget.see(tk.END)

logger.gui_writer = write_to_gui

#version info reader
#will add same one in extractor but may depreciate after moving metadata editing part to GUI
def get_version_info_for(module_name, manifest_path="project_version_manifest.ini"):
    config = configparser.ConfigParser()
    config.read(manifest_path)

    if module_name in config:
        version = config[module_name].get("version", "N/A")
        last_update = config[module_name].get("last_update", "N/A")
        return version, last_update
    return None, None

#newly added 20250824 for reading csv (is that essential?)
def _robust_read_csv(path, prefer_tab=False):
    """
    Try several parsing strategies (TSV first if prefer_tab=True).
    Returns a pandas DataFrame or raises the last error.
    """
    import pandas as pd, csv

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
            df = pd.read_csv(path, **kw)
            df.columns = [str(c).strip() for c in df.columns]
            return df
        except Exception as e:
            last_err = e
    raise last_err

#fix in 20251001
def _read_ion_df(path):
    import pandas as pd
    if not path:
        return None
    p = str(path).lower()
    ext = os.path.splitext(path)[1].lower() #need or not need
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
    elif p in (".csv", ".tsv"):
        sep = "\t" if ext == ".tsv" else ","
        return {"ionlist": pd.read_csv(path, sep=sep, encoding="utf-8")}
    else:
        try:
            df = pd.read_csv(path, engine="python")
        except:
            raise ValueError(f"Unsupported ion-list file: {path}")

    # normalize a mass column name
    col_map = {c.lower(): c for c in df.columns}
    for key in ("mass", "mz", "ion_mz", "m/z"):
        if key in col_map:
            if key != "mass":
                df = df.rename(columns={col_map[key]: "mass"})
            break
    return df[["mass"]].dropna() if "mass" in df.columns else None


# Function to create and hide Tkinter root window
# If future version supports GUI fully available, please reconstruct this part (not hiding main window, but think about what should be there)
def create_tkinter_root():
    root = tk.Tk()
    root.withdraw() # hide main window
    return root

#about tab showing the version info
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
    config.read(".\\2025_demo\\project_version_manifest.ini")
    #config.read("project_version_manifest.ini")

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

# Function to validate file path
def validate_file_path(file_path):
    if not os.path.exists(file_path):
        print(f"Error: The file {file_path} does not exist.")
        return False
    if not os.path.isfile(file_path):
        print(f"Error: The path {file_path} is not a file.")
        return False
    return True

def update_display():
    text_widget.delete(1.0, tk.END)
    for ftype, path in selected_files.items():
        text_widget.insert(tk.END, f"{ftype}: {path}\n")

def clear_files():
    # Placeholder for your implementation
    selected_files.clear()
    update_display()

def reset_main_status():
    status_var.set("Idle")
    progress.stop()
    convert_button.config(state="normal")


def on_metadata_ready(raw_file, metadata, savename):
    logger.log(f"Confirmed metadata for {raw_file}")
    # Pass to peak extractor
    #20250916 comment this section to see if duplicate conversion can be avoided
    #mspext.convert_raw_to_csv(raw_file, debug=False)

def launch_metadata_for_all(files):
    MetadataEditorWindow(root, files, on_metadata_ready, on_finish=reset_main_status)


def build_final_path(metadata: Dict, out_dir: str, kind: str) -> str:
    """
    Build final CSV path from metadata and type ('ms2' or 'ms3').
    Assumes you already have fields like Experiment Title / Sample / Date etc.
    """
    # Example: adapt to your naming scheme
    stem = metadata.get("Generated_Filename_Stem") or metadata.get("Experiment Title") or "untitled"
    stem = str(stem).strip().replace(" ", "_")
    return str(Path(out_dir) / f"{kind}_{stem}.csv")

from typing import Dict

def build_final_path(metadata: Dict, out_dir: str, kind: str) -> str:
    """
    Build final CSV path from metadata and type ('ms2' or 'ms3').
    Adjust this to your exact naming rules.
    """
    # Prefer your existing 'savename' logic; this is a fallback
    stem = metadata.get("Generated_Filename_Stem") \
        or metadata.get("Experiment Title") \
        or metadata.get("Raw filename") \
        or "untitled"
    stem = str(stem).strip().replace(" ", "_")
    return str(Path(out_dir) / f"{kind}_{stem}.csv")


def finalize_converted_files(temp_ms2: str, temp_ms3: str, metadata: Dict, out_dir: str) -> Dict[str, str]:
    final_ms2 = build_final_path(metadata, out_dir, kind="ms2")
    final_ms3 = build_final_path(metadata, out_dir, kind="ms3")
    Path(final_ms2).parent.mkdir(parents=True, exist_ok=True)
    os.replace(temp_ms2, final_ms2)
    os.replace(temp_ms3, final_ms3)
    return {"ms2": final_ms2, "ms3": final_ms3}

"""
def finalize_converted_files(temp_ms2: str, temp_ms3: str, metadata: Dict, out_dir: str) -> Dict[str, str]:
    #Promote temp CSVs to final names using atomic rename (no second write).
    final_ms2 = build_final_path(metadata, out_dir, kind="ms2")
    final_ms3 = build_final_path(metadata, out_dir, kind="ms3")
    Path(final_ms2).parent.mkdir(parents=True, exist_ok=True)
    os.replace(temp_ms2, final_ms2)
    os.replace(temp_ms3, final_ms3)
    return {"ms2": final_ms2, "ms3": final_ms3}
"""

# Function to select raw file for further pre-processing
def select_file(filetype):
    filetypes_dict = {
        "raw": [("Raw file", "*.raw")],
        "mzml": [("mzML file", "*.mzML")],
        "csv": [("CSV file", "*.csv")],
        "excel": [("Excel file", "*.xls *.xlsx")]
    }

    #filepath = filedialog.askopenfilename(filetypes=filetypes_dict.get(filetype, [("All files", "*.*")]))
    filepaths = filedialog.askopenfilenames(filetypes=filetypes_dict.get(filetype, [("All files", "*.*")]))
    if filepaths:
        valid_paths = [fp for fp in filepaths if validate_file_path(fp)]
        selected_files.setdefault(filetype, []).extend(valid_paths)
        update_display()
        #validation of values
        #file_contents_validation
    elif not filepaths:
        print(f"No valid {filetype} file selected.")


# generic select file (20250407)
def select_files_generic(filetype_key, allow_multiple=False, on_select_callback=None):
    filetypes_dict = {
        "csv": [("CSV files", "*.csv")],
        "excel": [("Excel files", "*.xls *.xlsx")],
        "json": [("JSON files", "*.json")],
        "all": [("All files", "*.*")]
    }

    filetypes = filetypes_dict.get(filetype_key, filetypes_dict["all"])

    if allow_multiple:
        filepaths = filedialog.askopenfilenames(filetypes=filetypes)
    else:
        filepath = filedialog.askopenfilename(filetypes=filetypes)
        filepaths = [filepath] if filepath else []

    if filepaths and on_select_callback:
        on_select_callback(filepaths)


def launch_metadata_batch():
    if not mspext.pymsreader:
        logger.log("[ERROR] Raw file conversion is not available — pymsfilereader missing.")
        messagebox.showerror(
            "Raw Conversion Not Available",
            "Thermo MSFileReader or pymsfilereader is not installed.\n"
            "Please run on a compatible Windows system with the required libraries."
        )
        return

    if "raw" not in selected_files or not selected_files["raw"]:
        messagebox.showwarning("No Raw File", "Please select at least one raw file.")
        return

    rawfilelist = selected_files["raw"]
    status_var.set("Converting raw file to csv...check the metadata assignment window.")
    progress.start()
    MetadataEditorWindow(root, rawfilelist, on_metadata_ready, on_finish=reset_main_status)





def on_batch_conversion_complete(converted_raws):
    progress.stop()
    convert_button.config(state="normal")
    status_var.set("Raw file conversions complete. Proceeding to metadata entry...")
    messagebox.showinfo("Conversion Complete", "Raw file conversion finished.\nNow entering metadata for each file.")
    
    launch_metadata_for_all(converted_raws)

#added for pseudolabeling
from pathlib import Path

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


def _same_drive(a: str, b: str) -> bool:
    da = os.path.splitdrive(os.path.abspath(a))[0].lower()
    db = os.path.splitdrive(os.path.abspath(b))[0].lower()
    return da == db

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

    if _same_drive(src, dst):
        os.replace(src, dst)
        logger.log(f"[finalize] Promoted (atomic) {src} → {dst}")
    else:
        shutil.copy2(src, dst)
        os.remove(src)
        logger.log(f"[finalize] Promoted (copy+delete) {src} → {dst} (cross-drive)")

def _score_metadata_candidate(meta: dict, sample_name: str, csv_path: str | None) -> int:
    """Content-based score: does this metadata look like it belongs to the selected sample?"""
    sample_tok = (sample_name or "").lower()
    csv_tok = ""
    if csv_path:
        csv_tok = Path(csv_path).stem.lower()
        if csv_tok.endswith(".raw"):  # normalize "...raw.csv" stems
            csv_tok = csv_tok[:-4]
    tokens = []
    tokens.append(os.path.basename((meta.get("Raw filename") or "")).lower())
    tokens.append(os.path.basename((meta.get("Original raw file path") or "")).lower())
    tokens.append((meta.get("Experiment Title") or "").lower())

    score = 0
    for t in tokens:
        if not t:
            continue
        if sample_tok and sample_tok in t:
            score += 1
        if csv_tok and csv_tok in t:
            score += 1
    return score


def _choose_metadata_dialog(parent, candidates):
    """
    candidates: list of tuples (path, meta_dict, score)
    Returns: (path, meta_dict) or (None, None) if cancelled.
    """
    import tkinter as tk
    from tkinter import ttk

    win = tk.Toplevel(parent)
    win.title("Select metadata JSON for this sample")
    win.transient(parent)
    win.grab_set()

    cols = ("file", "raw", "exp", "score")
    tree = ttk.Treeview(win, columns=cols, show="headings", height=min(10, len(candidates)))
    for c, w in zip(cols, (36, 20, 26, 6)):
        tree.heading(c, text=c.upper())
        tree.column(c, width=12 * w, anchor="w")
    for p, meta, s in candidates:
        tree.insert("", "end", values=(
            os.path.basename(p),
            os.path.basename(meta.get("Raw filename") or ""),
            (meta.get("Experiment Title") or ""),
            s,
        ))
    tree.pack(fill="both", expand=True, padx=10, pady=8)

    sel = {"idx": None}
    def _ok():
        cur = tree.selection()
        if cur:
            sel["idx"] = tree.index(cur[0])
        win.destroy()

    btns = ttk.Frame(win); btns.pack(fill="x", padx=10, pady=(0,10))
    ttk.Button(btns, text="OK", command=_ok).pack(side="right", padx=6)
    ttk.Button(btns, text="Cancel", command=win.destroy).pack(side="right")

    win.wait_window()
    if sel["idx"] is None:
        return None, None
    p, meta, _ = candidates[sel["idx"]]
    return p, meta

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

#newly added
def _fallback_simple_ion_scoring(matched_df, ion_df, ppm_value):
    """Minimal ion score using validator.findingions + marker.score_counter (anchors-aware)."""
    import numpy as np, pandas as pd
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
# --- Combine Trainable Datasets (Plan 1) ---
def _infer_uid_series(df: pd.DataFrame, src_base: str):
    """
    Returns a UID series. Prefer MS2scan_no if present; otherwise index-based.
    UID format: {src_base}#scan{MS2scan_no}  OR  {src_base}#row{n}
    """
    # normalize possible scan column variants
    scan_candidates = ["MS2scan_no", "MS2scan", "ScanNo", "scan", "scan_no", "ms2_scan_no"]
    scan_col = next((c for c in scan_candidates if c in df.columns), None)

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
            df = _robust_read_csv(p, prefer_tab=prefer_tab)
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

#end 20250929

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
    #20250917 start fixing json issues
    # --- BEGIN: exp.json save/load helpers (generic; PL-ready) ---
    from pathlib import Path

    def _norm(p):
        if not p: return None
        return os.path.normpath(os.path.expanduser(str(p)))

    def _rel(p, base_dir):
        p = _norm(p)
        if not p: return p
        try:
            return os.path.relpath(p, base_dir)
        except Exception:
            return p

    def _abs(p, base_dir):
        if not p: return p
        return os.path.normpath(p if os.path.isabs(p) else os.path.join(base_dir, p))

    def _current_exp_title():
        """Best-effort: get the Experiment currently selected in the tree, else the first one."""
        sel = tree.selection()
        if sel:
            node = sel[0]
            parent = tree.parent(node)
            exp_node = parent or node
            txt = tree.item(exp_node, "text")
        else:
            roots = tree.get_children()
            if not roots:
                return None
            txt = tree.item(roots[0], "text")
        return txt.replace("Experiment: ", "").split(" (")[0].strip()

    def export_experiment_json(exp_title, out_path):
        """Serialize everything we currently know for that experiment (all keys under a sample)."""
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        base = os.path.dirname(out_path)

        payload = {
            "experiment": exp_title,
            "generated_on": datetime.now().strftime("%Y-%m-%d %H:%M"),
            "samples": {}
        }
        samples = experiment_projects.get(exp_title, {}).get("samples", {})

        for sname, files in samples.items():
            # Shallow copy and relativize every string value
            record = {}
            for k, v in (files or {}).items():
                if isinstance(v, str) and v.strip():
                    record[k] = _to_posix(v)   # <— store POSIX in JSON
                    #record[k] = _rel(v, base)
            payload["samples"][sname] = record

        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)

        experiment_method_paths[exp_title] = out_path
        if exp_title in experiment_status_labels:
            experiment_status_labels[exp_title].config(text=f"EXP file: {out_path}")
        logger.log(f"[EXP] Saved → {out_path}")
        return out_path

    def import_experiment_json(in_path):
        base = os.path.dirname(in_path)
        with open(in_path, "r", encoding="utf-8") as f:
            data = json.load(f)

        exp_title = data.get("experiment") or "Unnamed Experiment"
        experiment_projects.setdefault(exp_title, {"samples": {}})
        dst = experiment_projects[exp_title]["samples"]

        # Back-compat: accept both a flat per-sample dict or a {'files':{...}} shape
        for sname, sample_blob in (data.get("samples") or {}).items():
            if isinstance(sample_blob, dict) and "files" in sample_blob:
                files = sample_blob.get("files") or {}
            else:
                files = sample_blob or {}
            # absolutize
            # resolved = {k: _abs(v, base) for k, v in files.items() if isinstance(v, str)}
            resolved = {k: _abs(_from_posix(v), base) for k, v in files.items() if isinstance(v, str) and v.strip()}
            dst[sname] = resolved

        experiment_method_paths[exp_title] = in_path
        logger.log(f"[EXP] Loaded ← {in_path}")
        refresh_tree()
        return exp_title
    # --- END: exp.json save/load helpers ---
    ###

    #20250917 future function: integrity check
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
            "This feature is in development.\n\n"
            "Planned: compute & store file hashes, check missing/changed files, "
            "and help you relink or bulk change roots."
        )
        # (Later you might call into mspik.writehashtojson(...) etc.)

    ttk.Button(toolbar, text="Check integrity", command=_check_integrity_clicked)\
       .pack(side="left", padx=4)

    #20251002
    include_mass_feat_var = tk.BooleanVar(value=False)
    #20250905 added for negative label
    # --- Negative sampling (Prepare Dataset) ---
    add_negatives_var   = tk.BooleanVar(value=False)
    neg_ratio_var       = tk.DoubleVar(value=3.0)   # max neg : pos
    #neg_scorecol_var    = tk.StringVar(value="score")
    #neg_scorethr_var    = tk.DoubleVar(value=0.05)  # keep scans with max(score) < thr
    neg_markercols_var  = tk.StringVar(value="")    # e.g. "core_marker_b,core_marker_y"
    neg_markermin_var   = tk.IntVar(value=1)        # require < this many marker hits
    min_hits_var      = tk.IntVar(value=3)      # keep scans with < min_hits hits
    ppm_tol_var       = tk.StringVar(value="10")# ppm tolerance; text with validation
    n_features_var    = tk.IntVar(value=0)      # total ion features (from ion_df)
    gate_hint_var     = tk.StringVar(value="Ion features not loaded yet")
    # --- Ion suggestion options ---
    ion_suggest_enable_var   = tk.BooleanVar(value=False)
    ion_suggest_ppm_var      = tk.DoubleVar(value=10.0)
    ion_suggest_dafloor_var  = tk.DoubleVar(value=0.03)
    ion_suggest_minsupp_var  = tk.IntVar(value=5)     # min glycan support per ion
    ion_suggest_topk_var     = tk.IntVar(value=60)    # how many to export
    last_suggest_csv_var = tk.StringVar(value="")

    # --- Treeview UI ---
    tree = ttk.Treeview(subwin)
    tree.heading("#0", text="Dataset Explorer", anchor="w")
    tree.pack(expand=True, fill="both", padx=10, pady=10)

    # --- Tree logic ---
    def on_tree_select(event):
        sel = tree.selection()
        if not sel:
            link_button.config(state="disabled")
            return

        item_text = tree.item(sel[0], "text")
        sample_name = clean_sample_name(item_text)

        parent_id = tree.parent(sel[0])
        exp_text = tree.item(parent_id, "text") if parent_id else ""
        exp_name = exp_text.replace("Experiment: ", "").split(" (")[0].strip()

        # Check if it's a valid sample in the experiment and if  it's okay to merge
        if exp_name in experiment_projects and sample_name in experiment_projects[exp_name]["samples"]:
            if (exp_name, sample_name) not in linked_validated_samples:
                link_button.config(state="normal")
                merge_button.config(state="disabled")
            else:
                link_button.config(state="disabled")
                merge_button.config(state="normal")

    tree.bind("<<TreeviewSelect>>", on_tree_select)

    def try_link_selected_sample():
        sel = tree.selection()
        if not sel:
            return
        sample_id = sel[0]
        exp_id = tree.parent(sample_id)
        sample_name = clean_sample_name(tree.item(sel[0], "text"))
        #sample_name = tree.item(sample_id, "text").replace("Sample: ", "").split(" (")[0]
        exp_name = tree.item(exp_id, "text").replace("Experiment: ", "")
        link_and_validate_sample(exp_name, sample_name)

    # --- Tree refresh logic ---
    def refresh_tree():
        tree.delete(*tree.get_children())

        for exp_title, exp_data in experiment_projects.items():
            exp_node = tree.insert("", "end", text=f"Experiment: {exp_title}", open=True)

            for sample_name, files in exp_data.get("samples", {}).items():
                has_csv = bool(files.get("csv"))
                has_excel = bool(files.get("excel"))
                has_json = bool(files.get("json"))

                # Decide how to display the sample label
                if (exp_title, sample_name) in linked_validated_samples:
                    sample_display = f"✅ Sample: {sample_name}"
                elif (exp_title, sample_name) in validation_failed_samples:
                    sample_display = f"⛔ Sample: {sample_name} (validation failed)"
                elif has_csv and has_excel and has_json:
                    sample_display = f"⚠️ Sample: {sample_name} (unvalidated)"
                elif not has_json:
                    sample_display = f"❌ Sample: {sample_name} (metadata missing)"
                else:
                    sample_display = f"Sample: {sample_name}"

                sample_node = tree.insert(exp_node, "end", text=sample_display, open=True)

                for ftype in ["csv", "excel", "json", "metadata", "insilico_csv", "ionlist_path", "pseudolabel_csv"]:
                    if files.get(ftype):
                        if ftype == "json":
                            label = "Method"
                        elif ftype == "metadata":
                            label = "Metadata"
                        elif ftype == "insilico_csv":
                            label = "In-silico CSV"
                        elif ftype == "ionlist_path":
                            label = "Ion List"
                        elif ftype == "pseudolabel_csv":
                            label = "Pseudolabel CSV"
                        else:
                            label = ftype.upper()
                        tree.insert(sample_node, "end", text=f"{label}: {os.path.basename(files[ftype])}")
    # --- statistics ---
    def show_experiment_summary(exp_name):
        if exp_name not in experiment_projects:
            return

        total = 0
        validated = 0

        for sname in experiment_projects[exp_name]["samples"]:
            total += 1
            if (exp_name, sname) in linked_validated_samples:
                validated += 1

        unvalidated = total - validated

        summary = (
            f"🧪 Experiment Summary: {exp_name}\n\n"
            f"• Total Samples: {total}\n"
            f"• ✅ Validated: {validated}\n"
            f"• ⚠️ Unvalidated: {unvalidated}"
        )
        messagebox.showinfo("Experiment Status", summary)
    # --- prevent symbols getting read --- v9 changed 
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
    # --- Add sample to selected experiment ---
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
        if sample_name in experiment_projects[exp_title]["samples"]:
            messagebox.showwarning("Duplicate Sample", f"Sample '{sample_name}' already exists.")
            return

        experiment_projects[exp_title]["samples"][sample_name] = {"csv": None, "excel": None, "json": None}
        refresh_tree()

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

    def store_metadata_back_to_sample(exp_name, sample_name, rawfile, metadata, savename):
        json_path = os.path.join(os.path.dirname(rawfile), savename + ".json")

        # Fallback if sample doesn't exist yet (e.g., metadata was created before file assignment)
        if sample_name not in experiment_projects[exp_name]["samples"]:
            experiment_projects[exp_name]["samples"][sample_name] = {"csv": None, "excel": None, "json": None}

        sample = experiment_projects[exp_name]["samples"][sample_name]
        sample["json"] = json_path

        # Try to rename sample to raw name (from metadata), only if different
        try:
            raw_base = os.path.splitext(os.path.basename(metadata.get("Raw filename", "")))[0]
        except Exception:
            raw_base = sample_name

        if raw_base != sample_name and raw_base not in experiment_projects[exp_name]["samples"]:
            experiment_projects[exp_name]["samples"][raw_base] = sample
            del experiment_projects[exp_name]["samples"][sample_name]
            sample_name = raw_base

        refresh_tree()
        write_method_file(exp_name)
    
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
            on_each_metadata_ready_callback=lambda rf, md, name: store_metadata_back_to_sample(exp_name, sample_name, rf, md, name),
            output_dir=outdir,
            skip_conversion=True
        )

        # Delay execution of load_file until window is fully initialized
        if not editor.skip_conversion:
            editor.window.after(10, lambda: editor.load_file(raw_file_path))
        else:
            print(f"[debug] skipping loading raw file")

    #live update of ion hits

    def current_selected_files():
        """Return the files dict for the currently selected sample,
        even if a child file node is selected."""
        sel = tree.selection()
        if not sel:
            return None

        node = sel[0]
        text = tree.item(node, "text")

        # If a file row like 'CSV: ...' is selected, go up to the sample row
        if ":" in text:
            node = tree.parent(node)
            text = tree.item(node, "text")

        # Now `node` should be a sample row; derive sample & experiment names
        sample_name = clean_sample_name(text)
        exp_node    = tree.parent(node)
        if not exp_node:
            return None

        exp_text = tree.item(exp_node, "text") or ""
        exp_name = exp_text.replace("Experiment: ", "").split(" (")[0].strip()

        try:
            return experiment_projects[exp_name]["samples"][sample_name]
        except KeyError:
            return None

    def _fetch_ion_count_for_dialog():
        # Try: current sample’s Excel; else let user pick one.
        from tkinter import filedialog, messagebox
        files = current_selected_files()
        excel_path = None
        try:
            # replace this with how you access the selected sample's file map
            # e.g., files = get_selected_sample_files()
            #files = experiment_projects["samples"]#current_selected_files()  # <-- implement this small accessor
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

    # added 20250905 place here bc validate+merge use this function
    def open_negative_options_dialog():
        import re
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

        # live hint
        tk.Label(dlg, textvariable=gate_hint_var, fg="#555").grid(row=4, column=0, columnspan=2,
                                                                sticky="w", padx=10, pady=(6,10))

        tk.Button(dlg, text="Load ion list now…",
          command=_fetch_ion_count_for_dialog).grid(row=4, column=0, columnspan=2,
                                                            sticky="w", padx=10, pady=(6,0))
        # move the hint to next row
        tk.Label(dlg, textvariable=gate_hint_var, fg="#555").grid(row=5, column=0, columnspan=2,
                                                                sticky="w", padx=10, pady=(6,10))

        # react when min-hits changes
        def _on_hits_change(*_): _update_gate_hint()
        min_hits_var.trace_add("write", _on_hits_change)

        # initialize text
        _update_gate_hint()

        tk.Button(dlg, text="Close", command=dlg.destroy).grid(row=5, column=0, columnspan=2, pady=(4,10))

    """ #score version
    def open_negative_options_dialog():
        dlg = tk.Toplevel(root)  # use the same parent you use elsewhere
        dlg.title("Negative Sampling Options")
        dlg.resizable(False, False)
        dlg.grab_set()

        tk.Checkbutton(dlg, text="Add Non-glycan entries (easy negatives)",
                    variable=add_negatives_var).grid(row=0, column=0, columnspan=2,
                                                        sticky="w", padx=10, pady=(10,4))

        tk.Label(dlg, text="Max ratio (neg:pos):").grid(row=1, column=0, sticky="e", padx=10)
        tk.Spinbox(dlg, from_=0.0, to=10.0, increment=0.5, width=6,
                textvariable=neg_ratio_var).grid(row=1, column=1, sticky="w", padx=6, pady=2)

        #tk.Label(dlg, text="Score column:").grid(row=2, column=0, sticky="e", padx=10)
        #tk.Entry(dlg, textvariable=neg_scorecol_var, width=18).grid(row=2, column=1, sticky="w", padx=6, pady=2)

        #tk.Label(dlg, text="Score threshold:").grid(row=3, column=0, sticky="e", padx=10)
        #tk.Spinbox(dlg, from_=0.00, to=1.00, increment=0.01, width=6,
        #        textvariable=neg_scorethr_var).grid(row=3, column=1, sticky="w", padx=6, pady=2)

        tk.Label(dlg, text="Marker cols (comma):").grid(row=4, column=0, sticky="e", padx=10)
        tk.Entry(dlg, textvariable=neg_markercols_var, width=28).grid(row=4, column=1, sticky="w", padx=6, pady=2)

        tk.Label(dlg, text="Max allowed marker hits:").grid(row=5, column=0, sticky="e", padx=10)
        tk.Spinbox(dlg, from_=0, to=5, increment=1, width=6,
                textvariable=neg_markermin_var).grid(row=5, column=1, sticky="w", padx=6, pady=(2,10))

        tk.Button(dlg, text="Close", command=dlg.destroy).grid(row=6, column=0, columnspan=2, pady=(4,10))
        """
    #20250906 add ion suggestion method function
    def open_ion_suggest_dialog():
        dlg = tk.Toplevel(root)
        dlg.title("Ion Suggestions")
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

    # --- link and validate the grouped sample ---
    def link_and_validate_sample(exp_name, sample_name):

        sample_name = clean_sample_name(sample_name)
        if sample_name not in experiment_projects[exp_name]["samples"]:
            logger.log(f"[ERROR] Cleaned sample name '{sample_name}' not found under '{exp_name}'")
            messagebox.showerror("Invalid Sample", f"Sample not found in experiment: {sample_name}")
            return
        sample = experiment_projects[exp_name]["samples"][sample_name]

        # Check file presence
        if not sample.get("csv") or not sample.get("excel"):
            messagebox.showerror("Missing File", "Sample must have both a CSV and Excel file before linking.")
            return
        

        # Validate CSV format
        #try:
        if not mspval.validate_csv_structure(sample["csv"]): #when it returns False
        #except Exception as e:
            
            messagebox.showerror("CSV Validation Failed")#, f"{sample['csv']}")
            logger.log(f"[Validation] CSV failed: {sample['csv']}\n")
            #raise ValueError("csv read error: file may be corrupted or unreadable.")
            validation_failed_samples.add((exp_name, sample_name))
            refresh_tree()
            return  # ← this prevents continuing to metadata
        else:
            logger.log(f"[Validation] CSV successful: {sample['csv']}\n")
        #Validate Excel format (assume sheet 'MSlist' exists and ion list is valid) 

        try:
            xl = pd.read_excel(sample["excel"], sheet_name=None)
            if not isinstance(xl, dict):
                raise ValueError("Excel read error: file may be corrupted or unreadable.")
            if "MSlist" not in xl:
                raise ValueError("Sheet 'MSlist' not found in Excel.")
        except Exception as e:
            messagebox.showerror("Excel Validation Failed", f"{sample['excel']}\n\n{e}")
            logger.log(f"[Validation] Excel failed: {sample['excel']}\n{e}")
            validation_failed_samples.add((exp_name, sample_name))
            refresh_tree()
            return  # ← this prevents continuing to metadata

        #try:
        anno = pd.ExcelFile(sample["excel"])
        MSlistdf = pd.read_excel(anno, sheet_name="MSlist")
        ionlistdf = pd.read_excel(anno, sheet_name="ionlist")
        if not mspval.validate_annotation_structure(MSlistdf):
        #except Exception as e:
            messagebox.showerror("Excel Validation Failed")#), f"{sample['excel']}\n")
            logger.log(f"[Validation] Excel failed: {sample['excel']}\n")
            #raise ValueError("Excel read error: file may be corrupted or unreadable.")
            validation_failed_samples.add((exp_name, sample_name))
            refresh_tree()
            return  # ← this prevents continuing to metadata
        else:
            logger.log(f"[Validation] Excel successful: {sample['csv']}\n")
        ion_df = ionlistdf[["mass"]] 
        if ion_df["mass"].dtype != "float64":
            logger.log(f"[Validation] ion list has invalid values: {sample['csv']}\n")
            return

        # Check metadata
        if not sample.get("json"):
            proceed = messagebox.askyesno("Metadata Missing", "No metadata found. Would you like to create it now?")
            if not proceed:
                return
            open_metadata_editor_for_sample(exp_name, sample_name)
            return

        # All checks passed → mark as validated
        linked_validated_samples.add((exp_name, sample_name))
        write_method_file(exp_name)
        messagebox.showinfo("Validated", f"Sample '{sample_name}' under '{exp_name}' is now validated.")
        refresh_tree()

    from pathlib import Path

    def _resolve_method_json(files: dict) -> str | None:
        """
        Return an absolute path to the sample's method JSON.
        Works whether files['json'] is absolute, relative, or just a basename.
        Searches (in order): as-given, CSV folder, Excel folder, CWD.
        """
        cand = (files.get("json") or "").strip()
        if not cand:
            return None

        p = Path(cand)
        if p.is_absolute() and p.exists():
            return str(p.resolve())

        # Use the basename (user may have stored only the filename)
        name = p.name if p.name else cand

        # search roots: where users most often keep the JSON
        roots = []
        csvp   = files.get("csv")
        excelp = files.get("excel")
        if csvp:   roots.append(Path(csvp).parent)
        if excelp: roots.append(Path(excelp).parent)
        roots.append(Path.cwd())

        for r in roots:
            try:
                q = (r / name)
                if q.exists():
                    return str(q.resolve())
            except Exception:
                pass
        return None

    #
        #20250919?
    # 20250919 move from prepare dataset
    #20250918
    # ---- ML params context glue for the panel ----
    def _current_exp_title():
        """Best-effort experiment currently selected (you already use this pattern above)."""
        sel = tree.selection()
        if sel:
            node = sel[0]
            parent = tree.parent(node)
            exp_node = parent or node
            txt = tree.item(exp_node, "text")
        else:
            roots = tree.get_children()
            if not roots:
                return None
            txt = tree.item(roots[0], "text")
        return txt.replace("Experiment: ", "").split(" (")[0].strip()

    def _get_ml_context():
        """
        Returns the three things the panel needs:
        - exp_json: path to the current experiment's .exp.json (if any)
        - method_json: path to the currently selected sample's method .json (if any)
        - current_params: optional in-memory set (we're not keeping one, so None)
        """
        exp_title = _current_exp_title()
        exp_json_path = experiment_method_paths.get(exp_title)

        files = current_selected_files() or {}
        # You already have this helper in the same function:
        method_json_path = _resolve_method_json(files)

        return {
            "exp_json": exp_json_path,
            "method_json": method_json_path,
            "current_params": None,
        }


    def try_merge_selected_sample():
        sel = tree.selection()
        if not sel:
            return

        sample_node = sel[0]
        sample_name = clean_sample_name(tree.item(sample_node, "text"))
        exp_node = tree.parent(sample_node)
        exp_name = tree.item(exp_node, "text").replace("Experiment: ", "").split(" (")[0].strip()

        files = experiment_projects[exp_name]["samples"][sample_name]
        if not all([files.get("csv"), files.get("excel"), files.get("json")]):
            messagebox.showerror("Error", "Sample is missing required files.")
            return

        # Ask user for output folder
        outdir = filedialog.askdirectory(title="Select output folder to save merged dataset")
        if not outdir:
            return

        # Ask user for ion sheet name (later we will allow external ion file)
        #ion_sheet_name = simpledialog.askstring("Ion Sheet", "Enter ion sheet name (in Excel):", initialvalue="core_OG")
        
        #20250905 quick fix on path issue
        #metadata_path=files["json"]
        metadata_path = _resolve_method_json(files)
        if not metadata_path:
            messagebox.showerror(
                "Merge Failed",
                "Could not locate the sample method JSON.\n\n"
                f"Original value: {files.get('json')}\n"
                "Tried: CSV folder, Excel folder, and the current working directory."
            )
            return

        print(f"[DEBUG] Method JSON resolved to: {metadata_path}")
        #print(f"debug: metadata file: {metadata_path}")
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
                
            
            today = datetime.now().strftime("%Y%m%d")
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
                    annotated_scans = pre_df["MS2Scan_no"],
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
            #another old version
            """
            if add_negatives_var.get():
                try:
                    neg_df = mspval.sample_real_negatives(
                        raw_tsv_path = files["csv"],
                        annotated_scans = pre_df["MS2scan_no"],
                        ion_df = ion_df,               # <-- use the same ion list you export with
                        ppm_tol = 10.0,                # match your exporter tolerance (you use 10 here) 
                        min_hits = int(neg_markermin_var.get()),  # e.g., 1 → keep only scans with 0 hits
                        max_neg_ratio = float(neg_ratio_var.get())
                    )
                    if len(neg_df):
                        pre_df = pd.concat([pre_df, neg_df], ignore_index=True)
                        print(f"[Prepare] Added {len(neg_df)} Non-glycan negatives (cap {neg_ratio_var.get():.1f}×).")
                    else:
                        print("[Prepare] No eligible negatives with current gates.")
                except Exception as e:
                    print(f"[Prepare] Negative sampling skipped due to error: {e}")
            """
            """ #old version using score
            if add_negatives_var.get():
                try:
                    marker_cols = [c.strip() for c in neg_markercols_var.get().split(",") if c.strip()] or None
                    score_col = (neg_scorecol_var.get().strip() or None)
                    neg_df = mspval.sample_real_negatives(
                        raw_tsv_path = files["csv"],
                        annotated_scans = pre_df["MS2scan_no"],
                        score_col = score_col,
                        score_thr = float(neg_scorethr_var.get()),
                        marker_cols = marker_cols,
                        marker_min = int(neg_markermin_var.get()),
                        max_neg_ratio = float(neg_ratio_var.get()),
                    )
                    added = len(neg_df)
                    if added:
                        pre_df = pd.concat([pre_df, neg_df], ignore_index=True)
                        print(f"[Prepare] Added {added} Non-glycan negatives (cap {neg_ratio_var.get():.1f}×).")
                    else:
                        print("[Prepare] No eligible negatives found with current gates.")
                except Exception as e:
                    print(f"[Prepare] Negative sampling skipped due to error: {e}")
            """

            #add ion suggestion
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
                        print(f"[Prepare] Ion suggestions saved: {suggest_csv}")
                        messagebox.showinfo("Ion suggestions",
                            f"Suggested ions written to:\n{os.path.basename(suggest_csv)}")
                        last_suggest_csv_var.set(suggest_csv)
                    except Exception as e:
                        messagebox.showwarning("Ion suggestions", f"Suggestion failed:\n{e}")


            """
            if ion_suggest_enable_var.get():
                if export_ion_suggestions_csv is None:
                    messagebox.showwarning("Ion suggestions",
                                        "Module msp_ion_mining.py not found; skipping.")
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
                        print(f"[Prepare] Ion suggestions saved: {suggest_csv}")
                        messagebox.showinfo("Ion suggestions",
                            f"Suggested ions written to:\n{os.path.basename(suggest_csv)}")
                        last_suggest_csv_var.set(suggest_csv)
                    except Exception as e:
                        messagebox.showwarning("Ion suggestions", f"Suggestion failed:\n{e}")
            """
            mspval.createnormailzedionlistcsv(iondfindex, pre_df,ion_df, outpath)
            messagebox.showinfo("Merge Complete", f"Dataset saved:\n{os.path.basename(os.path.basename(outpath))}")
        except Exception as e:
            messagebox.showerror("Merge Failed", f"Error:\n{str(e)}")

    #20250911 move convert pseudolabeled data to trainable csv functionality to here
    # ---------- PSEUDOLABEL → TRAINABLE (one-pass) ----------
    import ast, re, math
    from datetime import datetime

    _TUPLE_LIKE_RE = re.compile(r"^\(\s*\d+\s*,\s*\d+\s*,\s*\d+\s*,\s*\d+\s*,\s*\d+\s*,\s*\d+\s*\)$")
    _FHNSKDN_RE    = re.compile(r"^(F\d+)?(H\d+)?(N\d+)?(S\d+)?(G\d+)?(KDN\d+)?$")
    _NONMASS = {
        "entry_no","MS1scan_no","MS1_isolationmass","MS1_monoisolationmass","chargeState",
        "protonatedmass","MS2scan_no","label","Structure","composition","theoretical_mass",
        "observed_mass","ppm_error","ion score","ion hit count","ion hits m/z",
        "ion hits intensity","ion hits logI","ion hits relI","ID","Source","unique_ID",
        "IUPACname(optional)","Glycanannotation2","GlyToucan ID"
    }


    #20250912@mark fix MS2Scan_no missing in pseudo -> trainable csv
    SCAN_CANDIDATES = ("MS2scan_no", "MS2Scan_no", "ScanNum", "scan", "Scan", "unique_ID")
    def _canonicalize_scan_column(df):
        """
        Ensure there is a canonical 'MS2scan_no' string column in df.
        If a candidate exists, copy/rename it; otherwise leave df unchanged.
        Drops duplicate candidate columns after promoting.
        """
        import pandas as pd
        if df is None or df.empty:
            return df
        src = next((c for c in SCAN_CANDIDATES if c in df.columns), None)
        if not src:
            return df
        # create/overwrite canonical
        df = df.copy()
        df["MS2Scan_no"] = df[src].astype(str)
        # drop other candidates except the canonical
        for c in SCAN_CANDIDATES:
            if c in df.columns and c != "MS2scan_no":
                df.drop(columns=[c], inplace=True)
        return df

    def _coerce_comp_to_tuple(x):
        if isinstance(x, (list, tuple)) and len(x) == 6:
            return tuple(int(v) for v in x)
        s = str(x).strip()
        if _TUPLE_LIKE_RE.match(s):
            try:
                t = ast.literal_eval(s)
                if isinstance(t, (list, tuple)) and len(t) == 6:
                    return tuple(int(v) for v in t)
            except Exception:
                return None
        return None

    def _infer_ion_masses_from_wide_df(df):
        masses = []
        for c in df.columns:
            if c in _NONMASS: continue
            try:
                masses.append(float(c))
            except Exception:
                continue
        return sorted(set(masses))
    #20250912 fix MS2scan missing issue
    SCAN_CANDIDATES = ("MS2scan_no","MS2Scan_no","ScanNum","scan","Scan","unique_ID")

    def _ensure_scan(df, fallback=None):
        """
        Guarantee a canonical string column 'MS2scan_no' exists in df.
        If a candidate exists, copy→cast→drop dups.
        Else, if fallback (Series/array) matches length, insert it.
        Returns df (copy).
        """
        import pandas as pd
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
        print("[PL→Train][debug] _ensure_scan failed; df cols:", list(df.columns)[:20], "len=", len(df))
        raise KeyError("MS2scan_no")

    #20250912 fix bad calls (the function looks not identical to score_counter in marker)

    def _clean_cols(df):
        if df is None or df.empty: 
            return df
        df = df.copy()
        df.columns = [str(c).strip().replace("\ufeff","") for c in df.columns]
        return df

    def _comp_tuple_to_label(x):
        """
        Accepts a tuple/list/str and returns compact label like F1H4N2S3 (omit zeros).
        Order: (Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc) → H,N,S,G,K,F in label
        If you use F,H,N,S,G,K as your canonical, adjust the mapping below.
        """
        import ast, json
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
                x = [int(float(p)) for p in parts] if parts else []
        if isinstance(x, (list, tuple)):
            # Assume library order: Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc
            # Label order (O/N both): F, H, N, S(=NeuAc), G(=NeuGc), K(=KDN)
            if len(x) < 6:
                x = list(x) + [0]*(6-len(x))
            hex_, hexc, neuac, neugc, kdn, fuc = [int(v) for v in x[:6]]
            parts = []
            if fuc:   parts.append(f"F{fuc}")
            if hex_:  parts.append(f"H{hex_}")
            if hexc:  parts.append(f"N{hexc}")
            if neuac: parts.append(f"S{neuac}")
            if neugc: parts.append(f"G{neugc}")
            if kdn:   parts.append(f"K{kdn}")
            return "".join(parts) if parts else None
        return None

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

    def build_trainable_from_pseudolabels(
        sample_name: str,
        pseudo_path: str,
        ion_file_path: str | None,
        ion_sheet_name: str | None,
        salvage_path: str | None,
        thresholds: dict,
        neg_opts: dict,
        feature_mode: str,          # "rebuild" or "reuse"
        wide_feat_csv: str | None,  # used when feature_mode == "reuse"
        output_path: str | None,
        logger=None,
        precursor_gate_ppm: float | None = None,   # e.g., 10.0 to enable, None to keep legacy behavior
        converted_csv_path: str | None = None,     # converted TSV with MS2scan_no + protonatedmass
        include_mass_feature: bool = False,               # <— NEW #20251002
        ):
        """
        One-pass builder:
        - Load pseudolabeled long TSV
        - Normalize labels (compact FHNSGKDN)
        - Threshold selection (min ion score, max |ppm|, Top-N)
        - Optional salvage override
        - Build/Reuse features + optional negatives
        - Save trainable CSV
        Returns: (outpath, summary_dict)
        """
        import pandas as pd, os

        log = (logger.log if logger else print)
        log(f"[PL→Train] starting for sample={sample_name}")

        # 1) Load pseudo TSV (long form)
        pl = _robust_read_csv(pseudo_path, prefer_tab=True) #always tsv
        #new
        pl = _clean_cols(pl)
        if pl is None or pl.empty:
            raise RuntimeError("Pseudolabeled TSV is empty or unreadable.")
        # unify case/aliases early
        aliases = {c.lower(): c for c in pl.columns}
        def has(col): return col in pl.columns
        def has_lower(col): return col.lower() in aliases

        #20251002 probably include mass to trainable datasets
        # --- [MS1 attach] make sure precursor mass is present when gate or feature needs it
        # detect PL scan column now (we'll reuse it)
        _scan_candidates = ("MS2scan_no","MS2Scan_no","ScanNum","scan","Scan","unique_ID")
        pl_scan_col = next((c for c in _scan_candidates if c in pl.columns), None)
        if not pl_scan_col and any(k.lower() in aliases for k in _scan_candidates):
            # resolve via lowercase aliases if needed
            pl_scan_col = aliases[next(k for k in (s.lower() for s in _scan_candidates) if k in aliases)]
        if not pl_scan_col:
            raise RuntimeError("Pseudolabels lack an MS2 scan column (e.g., 'MS2scan_no').")
        pl[pl_scan_col] = pl[pl_scan_col].astype(str)

        #debug lines
        print("[PL→Train][dbg] scan_col in PL:", pl_scan_col)
        print("[PL→Train][dbg] include_mass_feature (arg):", include_mass_feature)
        print("[PL→Train][dbg] columns in PL (first 20):", pl.columns.tolist()[:20])
        print("[PL→Train][dbg] 'protonatedmass' in PL before MS1 attach:", "protonatedmass" in pl.columns)
        print("[PL→Train][dbg] converted_csv_path:", converted_csv_path)
        #####

        need_ms1 = bool(include_mass_feature) or (precursor_gate_ppm is not None)
        if need_ms1 and "protonatedmass" not in pl.columns:
            if not converted_csv_path or not os.path.exists(converted_csv_path):
                raise RuntimeError("MS1 mass required (gate/feature), but converted MS2 CSV/TSV is missing.")
            conv = _robust_read_csv(converted_csv_path, prefer_tab=True)
            # find scan + mass columns in converted file
            conv_scan = next((c for c in _scan_candidates if c in conv.columns), None)
            conv_mz   = next((c for c in ("protonatedmass","ProtonatedMass","precursor_mass","mz","MZ") if c in conv.columns), None)
            if not conv_scan or not conv_mz:
                raise RuntimeError("Converted CSV/TSV must include MS2scan_no and protonatedmass-like columns.")
            conv = conv[[conv_scan, conv_mz]].dropna()
            conv[conv_scan] = conv[conv_scan].astype(str)
            pl = pl.merge(conv.rename(columns={conv_mz: "protonatedmass"}), left_on=pl_scan_col, right_on=conv_scan, how="left").drop(columns=[conv_scan])

        #debug lines
        print("[PL→Train][dbg] after MS1 attach: 'protonatedmass' in PL =", "protonatedmass" in pl.columns)
        if "protonatedmass" in pl.columns:
            print("[PL→Train][dbg] protonatedmass head:", pl["protonatedmass"].head(5).tolist())

        #add here? to apply precursor filter
        import numpy as np, pandas as pd

        if precursor_gate_ppm is not None:
            # Need the library mass per composition; accept any available column name
            lib_mass_col = next((c for c in ("theoretical_mass","Mass","mass","TheoMass") if c in pl.columns), None)
            if lib_mass_col is None:
                raise RuntimeError("Pseudolabels lack a theoretical mass column to compare against.")

            with np.errstate(divide="ignore", invalid="ignore"):
                pl["ppm_precursor"] = 1e6 * (pl["protonatedmass"] - pl[lib_mass_col]) / pl[lib_mass_col]

            before = len(pl)
            pl = pl[pl["ppm_precursor"].abs() <= float(precursor_gate_ppm)]
            after = len(pl)
            if logger: logger.log(f"[PL→Train] precursor gate {precursor_gate_ppm} ppm: kept {after}/{before} rows")
        if "ppm_precursor" in pl.columns:
            print("[PL→Train][dbg] after gate: kept rows:", len(pl))
        """
        if precursor_gate_ppm is not None:
            if not converted_csv_path or not os.path.exists(converted_csv_path):
                raise RuntimeError("Precursor gate requested, but converted_csv_path is missing.")
            conv = _robust_read_csv(converted_csv_path, prefer_tab=True)  # your converted MS2 table
            # Find scan + mass columns robustly
            scan_cand = [c for c in ("MS2scan_no","MS2Scan_no","ScanNum","scan","Scan","unique_ID") if c in conv.columns]
            mass_cand = [c for c in ("protonatedmass","ProtonatedMass","precursor_mass","mz","MZ") if c in conv.columns]
            if not scan_cand or not mass_cand:
                raise RuntimeError("Converted CSV lacks MS2scan_no or protonatedmass-like columns.")
            sc_col, mz_col = scan_cand[0], mass_cand[0]

            # Normalize types/headers
            conv = conv[[sc_col, mz_col]].dropna()
            conv[sc_col] = conv[sc_col].astype(str)
            pl[sc_col]   = pl[sc_col].astype(str)

            # Attach protonatedmass to PL rows by scan
            pl = pl.merge(conv.rename(columns={mz_col: "protonatedmass"}), on=sc_col, how="left")

            # Need the library mass per composition; accept any available column name
            lib_mass_col = None
            for cand in ("theoretical_mass","Mass","mass","TheoMass"):
                if cand in pl.columns:
                    lib_mass_col = cand; break
            if lib_mass_col is None:
                raise RuntimeError("Pseudolabels lack a theoretical mass column to compare against.")

            # Compute ppm between precursor and library mass
            with np.errstate(divide="ignore", invalid="ignore"):
                pl["ppm_precursor"] = 1e6 * (pl["protonatedmass"] - pl[lib_mass_col]) / pl[lib_mass_col]
            # Gate
            before = len(pl)
            pl = pl[pl["ppm_precursor"].abs() <= float(precursor_gate_ppm)]
            after = len(pl)
            if logger: logger.log(f"[PL→Train] precursor gate {precursor_gate_ppm} ppm: kept {after}/{before} rows")

            # If you keep doing a second ppm filter on an older 'ppm_error', that’s fine;
            # the precursor gate is just an additional hard filter on top.
            """

        # preferred canonical names
        scan_candidates = ("MS2scan_no","MS2Scan_no","ScanNum","scan","Scan","unique_ID")
        scan_col = next((c for c in scan_candidates if has(c) or has_lower(c)), None)
        if scan_col and scan_col not in pl.columns and has_lower(scan_col):
            scan_col = aliases[scan_col.lower()]

        # make/standardize 'composition'
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
            print("[PL→Train][debug] pseudolabel headers:", pl.columns.tolist()[:30])
            raise RuntimeError("No composition/comp_tuple column in pseudolabels.")        



        # 2) Thresholding / selection
        min_score = float(thresholds.get("min_ion_score", 0.0))
        max_abs_ppm = float(thresholds.get("max_abs_ppm", 20.0))
        topn = int(thresholds.get("topn", 1))
        use_score = ("ion score" in pl.columns)
        if "ppm_error" not in pl.columns:
            pl["ppm_error"] = 9e9  # fallback if missing

        sel = pl.copy()
        if use_score:
            sel = sel[sel["ion score"] >= min_score]
        sel = sel[sel["ppm_error"].abs() <= max_abs_ppm]
        order_cols = [scan_col] + (["ion score"] if use_score else ["ppm_error"])
        ascending  = [True] + ([False] if use_score else [True])
        sel = sel.sort_values(order_cols, ascending=ascending).groupby(scan_col, as_index=False).head(topn)
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
        sel[scan_col] = sel[scan_col].astype(str)

        #debug track
        print("[PL→Train][dbg] sel cols:", [scan_col, "Structure"], 
            "null_Struct:", int(sel["Structure"].isna().sum()))
        # 3) Ion list
        ion_df = _read_ion_df(ion_file_path) if ion_file_path else None
        if ion_df is not None and ion_sheet_name:
            # (optional) if your _read_ion_df can select sheet, pass it there; otherwise ignore
            pass

        # 4) Features
        # If REBUILD: need long-form peaks and a numeric ion mass list
        # If REUSE : use the provided wide feature CSV and infer numeric masses (for negatives later)
        import numpy as np
        from ml_ng_utils import build_features_from_peaks_log10_plus1

        ion_masses = None
        if feature_mode == "rebuild":
            if ion_df is None or "mass" not in ion_df.columns:
                raise RuntimeError("Ion list with a 'mass' column is required to rebuild features.")
            ion_masses = ion_df["mass"].astype(float).tolist()

            # require long-form peaks to build features
            if not {"peaklist","peakintensity"}.issubset(pl.columns):
                raise RuntimeError("Pseudolabeled TSV must contain 'peaklist' and 'peakintensity' to rebuild features.")
            
            #20250912@mark 
            #pos_scans_df = sel[[scan_col,"peaklist","peakintensity"]].rename(columns={scan_col:"MS2scan_no"})
            #pos_scans_df["MS2scan_no"] = pos_scans_df["MS2scan_no"].astype(str)


            #pos_feat = build_features_from_peaks_log10_plus1(
            #    pos_scans_df, ion_masses, ppm=float(thresholds.get("ion_ppm", 10.0))
            #)
            #pos_feat = _ensure_scan(pos_feat, fallback=pos_scans_df["MS2scan_no"])
            # keep a private copy of scan ids for safety
            #pos_feat["_scan_fallback"] = pos_scans_df["MS2scan_no"].values  # NEW
            #20250912@mark fix MS2Scan missing
            #pos_feat = _canonicalize_scan_column(pos_feat)
            #20251030 fix MS2scan_no error
            pos_scans_df = sel[[scan_col,"peaklist","peakintensity"]].rename(columns={scan_col:"MS2scan_no"})
            pos_scans_df["MS2scan_no"] = pos_scans_df["MS2scan_no"].astype(str)

            # Build features
            pos_feat = build_features_from_peaks_log10_plus1(
                pos_scans_df, ion_masses, ppm=float(thresholds.get("ion_ppm", 10.0))
            )

            # --- NEW: robust scan handling even when features are empty ---
            scan_series = pos_scans_df["MS2scan_no"].astype("string")
            if pos_feat is None or pos_feat.empty:
                # Create a minimal shell so downstream code can proceed gracefully
                import pandas as pd
                pos_feat = pd.DataFrame({"MS2scan_no": scan_series})
                pos_feat["_scan_fallback"] = scan_series
            else:
                # Ensure canonical scan column then record fallback for later merges
                pos_feat = _ensure_scan(pos_feat, fallback=scan_series)
                pos_feat["_scan_fallback"] = pos_feat["MS2scan_no"].astype(str).values
            # --------------------------------------------------------------
            print("[PL→Train][dbg] ion_masses count:", 0 if ion_masses is None else len(ion_masses))
            """
            pos_scans_df = sel[[scan_col,"peaklist","peakintensity"]].rename(columns={scan_col:"MS2scan_no"})
            pos_scans_df["MS2scan_no"] = pos_scans_df["MS2scan_no"].astype(str)

            pos_feat = build_features_from_peaks_log10_plus1(pos_scans_df, ion_masses, ppm=float(thresholds.get("ion_ppm", 10.0)))
            # positives
            pos_feat = _ensure_scan(pos_feat, fallback=pos_scans_df["MS2scan_no"])
            pos_feat["_scan_fallback"] = pos_feat["MS2scan_no"].astype(str).values
            """
            #debug print
            print("[PL→Train][dbg] pos_feat pre-merge has:", 
            [c for c in ("MS2scan_no","Structure","_scan_fallback") if c in pos_feat.columns])

            # attach label next; do NOT subset columns yet
            pos_feat = pos_feat.merge(sel[[scan_col,"Structure"]],
                                    left_on="MS2scan_no", right_on=scan_col,
                                    how="left").drop(columns=[scan_col])    

            # === attach MS1 (protonatedmass) + delta_ppm into wide features when requested ===
            print("[PL→Train][dbg] rebuild: include_mass_feature", include_mass_feature)

            # === Include precursor mass as a feature (rebuild) ===========================
            if include_mass_feature:
                print("[PL→Train][dbg] rebuild: include_mass_feature =", include_mass_feature)

                # Build a (scan, protonatedmass) table from PL, renaming scan to MS2scan_no
                pm_src_scan = pl_scan_col  # detected earlier from PL
                if "protonatedmass" not in pl.columns:
                    print("[PL→Train][dbg][mass] WARNING: PL lacks 'protonatedmass' after MS1 attach; skipping mass feature")
                else:
                    pm = (pl[[pm_src_scan, "protonatedmass"]]
                            .dropna()
                            .drop_duplicates(pm_src_scan)
                            .rename(columns={pm_src_scan: "MS2scan_no"}))
                    pm["MS2scan_no"] = pm["MS2scan_no"].astype(str)
                    print("[PL→Train][dbg][mass] pm rows:", len(pm), 
                        "head:", pm.head(3).to_dict("records"))

                    # Make sure pos_feat actually has MS2scan_no; if not, rebuild it
                    if "MS2scan_no" not in pos_feat.columns:
                        print("[PL→Train][dbg][mass] pos_feat missing MS2scan_no; rebuilding from fallback…")
                        pos_feat = _ensure_scan(pos_feat, fallback=pos_feat.get("_scan_fallback"))
                        print("[PL→Train][dbg][mass] after ensure_scan, has MS2scan_no?",
                            "MS2scan_no" in pos_feat.columns)

                    # Final guard (fail soft if still missing)
                    if "MS2scan_no" in pos_feat.columns:
                        # Cast to string to avoid dtype merge issues
                        pos_feat["MS2scan_no"] = pos_feat["MS2scan_no"].astype(str)
                        # Merge
                        pos_feat = pos_feat.merge(pm, on="MS2scan_no", how="left")
                        print("[PL→Train][dbg][mass] merged mass → pos_feat cols now:",
                            [c for c in pos_feat.columns[:25]])
                        print("[PL→Train][dbg][mass] protonatedmass non-null count:",
                            int(pos_feat["protonatedmass"].notna().sum()))
                    else:
                        print("[PL→Train][dbg][mass] ABORT mass merge: no MS2scan_no in pos_feat")
            # ============================================================================
            """
            if include_mass_feature:
                if "protonatedmass" not in pl.columns:
                    print("[PL→Train][dbg][WARN] rebuild: PL still lacks protonatedmass; skip MS1 feature merge")
                else:
                    pm = (pl[[pl_scan_col, "protonatedmass"]]
                            .dropna()
                            .drop_duplicates(pl_scan_col)
                            .rename(columns={pl_scan_col: "MS2scan_no"}))
                    before_cols = pos_feat.columns.tolist()
                    pos_feat = pos_feat.merge(pm, on="MS2scan_no", how="left")
                    print("[PL→Train][dbg] rebuild: merged MS1 → new cols:",
                        [c for c in pos_feat.columns if c not in before_cols])
                    # optional delta_ppm if theo mass exists in PL
                    _lib_mass_pl = next((c for c in ("theoretical_mass","Mass","mass","TheoMass") if c in pl.columns), None)
                    if _lib_mass_pl:
                        tm = (pl[[pl_scan_col, _lib_mass_pl]]
                                .dropna()
                                .drop_duplicates(pl_scan_col)
                                .rename(columns={pl_scan_col: "MS2scan_no"}))
                        pos_feat = pos_feat.merge(tm, on="MS2scan_no", how="left")
                        with np.errstate(divide="ignore", invalid="ignore"):
                            pos_feat["delta_ppm"] = 1e6 * (pos_feat["protonatedmass"] - pos_feat[_lib_mass_pl]) / pos_feat[_lib_mass_pl]
                        print("[PL→Train][dbg] rebuild: added delta_ppm; head:",
                            pos_feat["delta_ppm"].head(5).tolist())
                    print("[PL→Train][dbg] rebuild: has protonatedmass in pos_feat?",
                        "protonatedmass" in pos_feat.columns)
                    if "protonatedmass" in pos_feat.columns:
                        print("[PL→Train][dbg] rebuild: protonatedmass head:",
                            pos_feat["protonatedmass"].head(5).tolist())
            
            # attach MS1 to wide table if user wants mass as a feature
            if include_mass_feature and "protonatedmass" in pl.columns and "MS2scan_no" in pos_feat.columns:
                pm = pl[[pl_scan_col, "protonatedmass"]].dropna().drop_duplicates(pl_scan_col)
                pos_feat = pos_feat.merge(pm.rename(columns={pl_scan_col: "MS2scan_no"}), on="MS2scan_no", how="left")

                # also add delta_ppm if theoretical mass exists in PL
                _lib_mass_pl = next((c for c in ("theoretical_mass","Mass","mass","TheoMass") if c in pl.columns), None)
                if _lib_mass_pl:
                    tm = pl[[pl_scan_col, _lib_mass_pl]].dropna().drop_duplicates(pl_scan_col)
                    pos_feat = pos_feat.merge(tm.rename(columns={pl_scan_col: "MS2scan_no"}), on="MS2scan_no", how="left")
                    with np.errstate(divide="ignore", invalid="ignore"):
                        pos_feat["delta_ppm"] = 1e6 * (pos_feat["protonatedmass"] - pos_feat[_lib_mass_pl]) / pos_feat[_lib_mass_pl]
            """
            #pos_feat = _ensure_scan(pos_feat, fallback=pos_scans_df["MS2scan_no"])
            #pos_feat["_scan_fallback"] = pos_scans_df["MS2scan_no"].values  # KEEP THIS

            #pos_scans_df["MS2scan_no"] = pos_scans_df["MS2scan_no"].astype(str)
            #pos_feat = pos_feat.merge(
            #    sel[[scan_col, "Structure"]],
            #    left_on="MS2scan_no",
            #    right_on=scan_col,
            #    how="left"
            #).drop(columns=[scan_col])

            # ensure scan column exists (from the earlier patch) then cast to str
            if "MS2scan_no" not in pos_feat.columns and len(pos_feat) == len(pos_scans_df):
                pos_feat.insert(0, "MS2scan_no", pos_scans_df["MS2scan_no"].values)
            pos_feat["MS2scan_no"] = pos_feat["MS2scan_no"].astype(str)

            if "Structure" not in pos_feat.columns:
                pos_feat = pos_feat.merge(
                    sel[[scan_col, "Structure"]],
                    left_on="MS2scan_no",
                    right_on=scan_col,
                    how="left",
                    suffixes=("", "_pl")  # avoid _x/_y confusion
                ).drop(columns=[scan_col])
            else:
                # Normalize any legacy duplicates from previous runs
                if "Structure_x" in pos_feat.columns and "Structure_y" in pos_feat.columns:
                    pos_feat["Structure"] = pos_feat["Structure_x"].fillna(pos_feat["Structure_y"])
                    pos_feat.drop(columns=["Structure_x","Structure_y"], inplace=True, errors="ignore")
                elif "Structure_pl" in pos_feat.columns:
                    pos_feat["Structure"] = pos_feat.get("Structure").fillna(pos_feat["Structure_pl"])
                    pos_feat.drop(columns=["Structure_pl"], inplace=True, errors="ignore")
            #debug print
            print("[PL→Train][dbg] pos_feat post-merge has:", 
                [c for c in ("MS2scan_no","Structure","Structure_pl","Structure_x","Structure_y","_scan_fallback") 
                if c in pos_feat.columns])

        else:  # feature_mode == "reuse"
            if not wide_feat_csv:
                raise RuntimeError("Provide a wide feature CSV when feature_mode='reuse'.")
            print("[PL→Train][dbg][reuse branch] pos_feat pre-merge has:", 
                [c for c in ("MS2scan_no","Structure","_scan_fallback") if c in pos_feat.columns])
            # 1) load the wide feature matrix
            feat = _robust_read_csv(wide_feat_csv)  # CSV (wide)
            if feat is None or feat.empty:
                raise RuntimeError("Wide feature CSV unreadable or empty.")

            # 2) detect its scan column *after* feat exists
            scan_feat = next((c for c in ("MS2scan_no","ScanNum","scan","Scan","unique_ID") if c in feat.columns), None)
            if not scan_feat:
                raise RuntimeError("No scan column found in wide feature CSV.")

            # 3) infer numeric ion masses from wide headers for later (negs)
            ion_masses = _infer_ion_masses_from_wide_df(feat)
            if not ion_masses:
                raise RuntimeError("No numeric ion masses inferred from wide feature CSV headers.")
            #20250912@mark
            feat[scan_feat] = feat[scan_feat].astype(str)
            # positives
            pos_feat = _ensure_scan(pos_feat, fallback=pos_scans_df["MS2scan_no"])
            pos_feat["_scan_fallback"] = pos_feat["MS2scan_no"].astype(str).values
            #debug print
            print("[PL→Train][dbg] pos_feat pre-merge has:", 
            [c for c in ("MS2scan_no","Structure","_scan_fallback") if c in pos_feat.columns])

            #20251001 fix by GPT?
            if "protonatedmass" in pl.columns and "MS2scan_no" in pos_feat.columns:
                pm = pl[["MS2scan_no","protonatedmass"]].dropna().drop_duplicates("MS2scan_no")
                pos_feat = pos_feat.merge(pm, on="MS2scan_no", how="left")
                print("[PL→Train][dbg] 20251001fix")
                if include_mass_feature:
                    _lib_mass_pl = next((c for c in ("theoretical_mass","Mass","mass","TheoMass") if c in pl.columns), None)
                    if _lib_mass_pl:
                        tm = pl[[pl_scan_col, _lib_mass_pl]].dropna().drop_duplicates(pl_scan_col)
                        pos_feat = pos_feat.merge(tm.rename(columns={pl_scan_col: "MS2scan_no"}), on="MS2scan_no", how="left")
                        with np.errstate(divide="ignore", invalid="ignore"):
                            pos_feat["delta_ppm"] = 1e6 * (pos_feat["protonatedmass"] - pos_feat[_lib_mass_pl]) / pos_feat[_lib_mass_pl]

            # attach label next; do NOT subset columns yet
            pos_feat = pos_feat.merge(sel[[scan_col,"Structure"]],
                                    left_on="MS2scan_no", right_on=scan_col,
                                    how="left").drop(columns=[scan_col])

            # attach label if missing
            if "Structure" not in pos_feat.columns:
                pos_feat = pos_feat.merge(
                    sel[[scan_col, "Structure"]],
                    left_on="MS2scan_no",
                    right_on=scan_col,
                    how="left",
                    suffixes=("", "_pl")  # avoid _x/_y confusion
                ).drop(columns=[scan_col])
            else:
                # Normalize any legacy duplicates from previous runs
                if "Structure_x" in pos_feat.columns and "Structure_y" in pos_feat.columns:
                    pos_feat["Structure"] = pos_feat["Structure_x"].fillna(pos_feat["Structure_y"])
                    pos_feat.drop(columns=["Structure_x","Structure_y"], inplace=True, errors="ignore")
                elif "Structure_pl" in pos_feat.columns:
                    pos_feat["Structure"] = pos_feat.get("Structure").fillna(pos_feat["Structure_pl"])
                    pos_feat.drop(columns=["Structure_pl"], inplace=True, errors="ignore")
            #debug print
            print("[PL→Train][dbg] pos_feat post-merge has:", 
            [c for c in ("MS2scan_no","Structure","Structure_pl","Structure_x","Structure_y","_scan_fallback") 
            if c in pos_feat.columns])
 
        # 5) Optionally add negatives (easy non-glycan) using the same ion set
        pos_ids_set = set(sel[scan_col].astype(str))
        add_negs = bool(neg_opts.get("enable", False))
        final_df = pos_feat.copy()
        if add_negs:
            # --- params ---
            min_hits  = int(neg_opts.get("min_hits", 3))
            ng_ppm    = float(neg_opts.get("ion_ppm", thresholds.get("ion_ppm", 10.0)))
            max_ratio = float(neg_opts.get("max_ratio", 3.0))
            mode      = str(neg_opts.get("sampling", "random")).lower()   # "random" | "first"
            seed      = int(neg_opts.get("seed", 42))
            NEG_LABEL = "Non-glycan"  # keep consistent with your ML prep

            # --- sanity: ion masses source ---
            if ion_df is None or ion_df.empty:
                ion_df = pd.DataFrame({"mass": ion_masses})

            # --- collect candidate negatives from long pseudolabel TSV ---
            #   (scans with < min_hits glycan-ion matches and not already positive)
            neg_rows = []
            for r in pl[[scan_col, "peaklist", "peakintensity"]].dropna().itertuples(index=False):
                sid = str(getattr(r, scan_col))
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
                print("[PL→Train][dbg][neg] kept rows:", len(neg_df),
                    "first scans:", neg_df["MS2scan_no"].head(5).tolist())

                # --- build features for negatives (only kept rows) ---
                neg_feat = build_features_from_peaks_log10_plus1(neg_df, ion_masses, ppm=ng_ppm)
                neg_feat = _ensure_scan(neg_feat, fallback=neg_df["MS2scan_no"])
                neg_feat["Structure"] = NEG_LABEL
                neg_feat["_scan_fallback"] = neg_feat["MS2scan_no"].astype(str).values

                print("[PL→Train][dbg][neg] after features:",
                    [c for c in ("MS2scan_no","Structure","_scan_fallback") if c in neg_feat.columns],
                    "first scans:", neg_feat["MS2scan_no"].head(5).tolist())
        # 6) Finalize + save
        if output_path is None or output_path.strip() == "":
            outdir = os.path.dirname(pseudo_path)
            # suffix reflects whether MS1 is part of the feature set
            mass_suffix = "_withMass" if include_mass_feature else "_noMass"
            outname = f"{sample_name}_trainable_fromPL_{datetime.now().strftime('%Y%m%d')}{mass_suffix}.csv"
            output_path = os.path.join(outdir, outname)
            if logger:
                logger.log(f"[PL→Train] include_mass_feature={include_mass_feature} → {os.path.basename(output_path)}")

        # assemble
        # --- assemble ---
        final_df = pos_feat if 'neg_feat' not in locals() else \
                pd.concat([pos_feat, neg_feat], ignore_index=True, sort=False)

        # Debug: check for protonatedmass survival before ensure_scan
        print("[PL→Train][dbg] after concat: has protonatedmass?","protonatedmass" in final_df.columns)
        if "protonatedmass" in final_df.columns:
            print("[PL→Train][dbg] final protonatedmass non-null:",
                int(final_df["protonatedmass"].notna().sum()))
        print("[PL→Train][dbg] after concat: cols=", final_df.columns.tolist()[:20])
        if "protonatedmass" in final_df.columns:
            print("[PL→Train][dbg] protonatedmass head:", final_df["protonatedmass"].head().tolist())

        # --- make sure we can reconstruct the scan column ---
        combined_fallback = final_df["_scan_fallback"].astype(str) if "_scan_fallback" in final_df.columns else None

        final_df = _ensure_scan(final_df, fallback=combined_fallback)

        # Debug after ensure_scan
        print("[PL→Train][dbg] after ensure_scan: cols=", final_df.columns.tolist()[:20])
        if "protonatedmass" in final_df.columns:
            print("[PL→Train][dbg] protonatedmass head:", final_df["protonatedmass"].head().tolist())

        # DEBUG (safe now)
        #print("[PL→Train][dbg] after ensure_scan: has_scan=", "MS2scan_no" in final_df.columns,
        #    "has_structure=", "Structure" in final_df.columns)
        #print("[PL→Train][dbg] first scans:", final_df["MS2scan_no"].head(5).tolist())

        # drop helper only after ensure_scan
        final_df.drop(columns=["_scan_fallback"], errors="ignore", inplace=True)

        print("[PL→Train][dbg] end-before-order: cols=", final_df.columns.tolist()[:30])

        if not include_mass_feature:
            final_df.drop(columns=["protonatedmass","delta_ppm","ppm_precursor","precursor_gate_comp"],
                        errors="ignore", inplace=True)
            print("[PL→Train][dbg] end: MS1 feature disabled → dropped MS1 columns")

        # order columns for saving
        feature_cols = [c for c in final_df.columns if c not in ("MS2scan_no","Structure")]
        if include_mass_feature and "protonatedmass" in final_df.columns:
            # ensure mass stays in feature columns
            if "protonatedmass" not in feature_cols:
                feature_cols.insert(0, "protonatedmass")
        else:
            # explicitly drop if user disabled
            feature_cols = [c for c in feature_cols if c != "protonatedmass"]
        # order columns for saving
        #feature_cols = [c for c in final_df.columns if c not in ("MS2scan_no","Structure")]
        
        final_df = final_df[["MS2scan_no","Structure"] + feature_cols]

        #do we need this really?
        #if not include_mass_feature:
        #    final_df.drop(columns=["protonatedmass","delta_ppm","ppm_precursor","precursor_gate_comp"],
        #                errors="ignore", inplace=True)

        final_df.to_csv(output_path, index=False)
        log(f"[PL→Train] saved: {output_path}")

        # summary
        classes = final_df["Structure"].value_counts().to_dict()
        summary = {
            "rows": int(len(final_df)),
            "cols": int(len(final_df.columns)),
            "classes": classes,
            "negatives_added": int((final_df["Structure"] == "None").sum()),
            "feature_mode": feature_mode,
            "ion_masses": len(ion_masses or []),
            "thresholds": thresholds,
            "neg_opts": neg_opts
        }
        return output_path, summary
    # ---------- end PSEUDOLABEL → TRAINABLE ----------



    #added 20250906 ion suggestion window?
    def open_ion_suggestions_viewer():
        import os, pandas as pd
        from tkinter import filedialog, messagebox, ttk

        path = last_suggest_csv_var.get().strip()
        if not path or not os.path.exists(path):
            # let user pick if we don't have a saved path yet
            path = filedialog.askopenfilename(
                title="Open ion suggestions CSV",
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

    #newly added
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

    def run_pseudolabeling(sample_name: str, files: dict, meta_overrides: dict, parent=None,
                        ppm_value: float = 20.0, keep_top_n_per_scan: int = 3,
                        ion_ppm: float = 10.0, anchors_required: int = 2):
        import pandas as pd, numpy as np, traceback, os
        from datetime import datetime
        # use your module helpers
        #label_style: str = "short"):
        """
        label_style: "short" -> A/s/p   ;  "long" -> A/Sul/Phos
        """

        csv_path = files.get("csv")
        ins_path = files.get("insilico_csv")
        ion_path = files.get("ionlist_path")

        if not csv_path or not os.path.exists(csv_path):
            messagebox.showwarning("Converted CSV missing", "Link a converted CSV for this sample.")
            return
        if not ins_path or not os.path.exists(ins_path):
            messagebox.showwarning("In-silico CSV missing", "Generate or link an in-silico CSV first.")
            return

        # 1) Load inputs
        df  = _robust_read_csv(csv_path, prefer_tab=True)   # converted TSV
        lib = _robust_read_csv(ins_path)                    # in-silico CSV

        # column heuristics (converted)
        scan_col = next((c for c in ["MS2scan_no","unique_ID","ScanNum","scan","Scan"] if c in df.columns), None)
        mass_col = next((c for c in ["protonatedmass","ProtonatedMass","precursor_mass","mz","MZ"] if c in df.columns), None)
        if not scan_col or not mass_col:
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

        libn = marker.normalize_insilico(
            lib,
            comp_cols=COMP_COLS_BASE,   # keep tuple = 6-core only
            mass_col="Mass",
            add_legacy_repr=True,
        ).copy()

        # IMPORTANT: keep modifier columns on the normalized table
        # normalize_insilico typically returns a row-per-entry frame preserving order,
        # so we can attach auxiliary columns directly by index alignment:
        for _c in MOD_COLS:
            if _c not in libn.columns and _c in lib.columns:
                libn[_c] = lib[_c].values

        libn = libn.sort_values("Mass").reset_index(drop=True)
        #old ones
        """
        # 2) Normalize in-silico library → comp_tuple/comp_str + sorted Mass
        libn = marker.normalize_insilico(
            lib,
            comp_cols=("Hex","HexNAc","NeuAc","NeuGc","KDN","Fuc"),
            mass_col="Mass",
            add_legacy_repr=True,
        )
        libn = libn.sort_values("Mass").reset_index(drop=True)
        """
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
        """
        for r in df[[scan_col, mass_col]].dropna().itertuples(index=False):
            scan, obs = getattr(r, scan_col), float(getattr(r, mass_col))
            hits = marker.find_compositions_for_mass(
                obs_mass=obs,
                insilico_sorted=libn,
                mass_col="Mass",
                ppm=ppm_value,
                mass_transform=None,  # Mass already protonated in your lib
                comp_cols=("Hex","HexNAc","NeuAc","NeuGc","KDN","Fuc"),
            )
            if not hits.empty:
                hits = hits.reindex(hits["ppm_error"].abs().sort_values().index)
                if keep_top_n_per_scan:
                    hits = hits.head(keep_top_n_per_scan)
                hits.insert(0, "MS2scan_no", scan)
                rows.append(hits[["MS2scan_no","observed_mass","theoretical_mass","ppm_error","comp_str","comp_tuple"]])

        matched = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
            columns=["MS2scan_no","observed_mass","theoretical_mass","ppm_error","comp_str","comp_tuple"]
        )
        """
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
        import math

        # 20250920 Try to add ion mining to PL
        # --- helper: build pre_df (positives from PL + optional sampled negatives) ---
        def _pre_df_for_ion_suggest_from_pl(matched, csv_path, ion_df):
            import pandas as pd
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
            print(f"[label] failed to build Pseudolabeling; reason={e}")
            matched["composition"] = ""

        # If downstream expects 'composition', mirror it (optional) #changed already above
        #matched["composition"] = matched["Predicted_Label"]

        """
        def _tuple_to_FHNSGKDN(t):
            # t: (Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc)
            h, n, s, g, kdn, f = [int(x) for x in t]
            parts = []
            if f:   parts.append(f"F{f}")
            if h:   parts.append(f"H{h}")
            if n:   parts.append(f"N{n}")
            if s:   parts.append(f"S{s}")
            if g:   parts.append(f"G{g}")
            if kdn: parts.append(f"KDN{kdn}")
            return "".join(parts) or "Non-glycan"

        if "comp_tuple" in matched.columns:
            # Build compact labels from tuples (most reliable)
            matched["composition"] = matched["comp_tuple"].apply(_tuple_to_FHNSGKDN)
        else:
            # Fallback: if composition already looks compact, keep it; else leave as-is
            import re
            _FHNSKDN_RE = re.compile(r"^(F\d+)?(H\d+)?(N\d+)?(S\d+)?(G\d+)?(KDN\d+)?$")
            looks_compact = matched["composition"].astype(str).str.match(_FHNSKDN_RE).all()
            if not looks_compact:
                print("[label] 'comp_tuple' missing; keeping composition as-is")
        """

        if matched.empty:
            messagebox.showinfo("Pseudolabeling", "No precursor matches within tolerance.")
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

        # --- 4.5) Ion suggestions from PL (requires peaklist/peakintensity now present) ---
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
                    print(f"[PL→Trainable] Ion suggestions saved: {suggest_csv}")
                    last_suggest_csv_var.set(suggest_csv)
                    messagebox.showinfo("Ion suggestions",
                                        f"Suggested ions written to:\n{os.path.basename(suggest_csv)}")
                except Exception as e:
                    messagebox.showwarning("Ion suggestions", f"Suggestion failed:\n{e}")
            else:
                messagebox.showwarning("Ion suggestions",
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
                    print("[ion] scorer returned no ion columns; falling back to simple score_counter")
                    raise RuntimeError("no_ion_columns")

            except Exception as e:
                import traceback; traceback.print_exc()
                # 5b) fallback — always produce basic ion columns with score_counter
                try:
                    matched = _fallback_simple_ion_scoring(matched, ion_df, ion_ppm)
                    ion_scoring_status = "ok(fallback)"
                    n_with_scores = int((matched["ion hit count"] > 0).sum())
                except Exception:
                    traceback.print_exc()
                    ion_scoring_status = "failed"

        # 6) Merge back onto the original converted file (one row per composition match)
        #enrich_cols = ["MS2scan_no","composition","theoretical_mass","ppm_error","observed_mass"]
        #for extra in ("ion score","ion hit count","ion hits m/z"):
        #    if extra in matched.columns: enrich_cols.append(extra)
        #out = df.merge(matched[enrich_cols], left_on=scan_col, right_on="MS2scan_no", how="left")
        scan_right = "MS2scan_no"
        if scan_right not in matched.columns:
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
       #20250930 fix win11 issue
       # 7) Save TSV next to converted CSV  (ABSOLUTE + explicit encoding)
        outdir  = os.path.dirname(os.path.abspath(csv_path))
        outname = f"{sample_name}_pseudolabels_{datetime.now().strftime('%Y%m%d')}.tsv"
        outpath = os.path.abspath(os.path.join(outdir, outname))
        out.to_csv(outpath, index=False, sep="\t", encoding="utf-8")

        # store absolute path so refresh_tree can always find it
        files["pseudolabel_csv"] = outpath

        """
        # 7) Save TSV next to converted CSV
        outdir  = os.path.dirname(csv_path)
        outname = f"{sample_name}_pseudolabels_{datetime.now().strftime('%Y%m%d')}.tsv"
        outpath = os.path.join(outdir, outname)
        out.to_csv(outpath, index=False, sep="\t")

        files["pseudolabel_csv"] = outpath
        """
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
                "Pseudolabeling complete", f"Saved and linked:\n{p}"
            ))
        subwin.after(0, refresh_tree)

        """
        if parent:
            messagebox.showinfo("Pseudolabeling complete",
                    "Saved and linked:\n" + _norm_for_display(outpath))
            logger.log(f"[PL] Saved and linked: {_norm_for_display(outpath)}")
            #messagebox.showinfo("Pseudolabeling complete", f"Saved and linked:\n{outpath}")
        try:
            refresh_tree()
        except Exception:
            pass
        """
    # --- Assign file to experiment/sample ---
    def assign_file(filetype, filepath, exp_title="Unassigned", sample_name="Unassigned"):
        if exp_title not in experiment_projects:
            experiment_projects[exp_title] = {"samples": {}}
        if sample_name not in experiment_projects[exp_title]["samples"]:
            experiment_projects[exp_title]["samples"][sample_name] = {"csv": None, "excel": None, "json": None}
        experiment_projects[exp_title]["samples"][sample_name][filetype] = filepath
        refresh_tree()

    def extract_rawname_from_metadata(json_path):
        try:
            with open(json_path, "r") as f:
                meta = json.load(f)
            rawbase = os.path.basename(meta.get("Raw filename", ""))
            return os.path.splitext(rawbase)[0]
        except:
            return None

    # --- File handlers ---
    def handle_csv_selection(filepaths):
        for path in filepaths:
            sample_id = os.path.splitext(os.path.basename(path))[0]
            assign_file("csv", path, "Unassigned", sample_id)

    def handle_excel_selection(filepaths):
        for path in filepaths:
            sample_id = os.path.splitext(os.path.basename(path))[0]
            assign_file("excel", path, "Unassigned", sample_id)

    def handle_json_selection(filepaths):
        for path in filepaths:
            title = extract_title_from_metadata(path) or "Unassigned"
            rawname = extract_rawname_from_metadata(path) or "Unassigned"

            assign_file("json", path, title, "Unassigned")

            sample_name = "Unassigned"

            # Auto-rename sample (move from 'Unassigned' to rawname)
            if title in experiment_projects:
                samples = experiment_projects[title]["samples"]
                if "Unassigned" in samples:
                    if rawname in samples:
                        messagebox.showwarning("Sample Exists", f"Sample '{rawname}' already exists. Skipping rename.")
                    else:
                        samples[rawname] = samples.pop("Unassigned")
                        sample_name = rawname  # <- use updated sample name
                        refresh_tree()

            # After rename, auto-write .method.json if csv + excel exist
            sample = experiment_projects[title]["samples"].get(sample_name)
            if sample and sample.get("csv") and sample.get("excel"):
                from datetime import datetime
                entry = {
                    "csv": os.path.basename(sample["csv"]),
                    "excel": os.path.basename(sample["excel"]),
                    "metadata": os.path.basename(path),
                    "raw_file": rawname,
                    "validated": False
                }
                sample_method = {
                    "experiment": title,
                    "samples": {sample_name: entry},
                    "generated_on": datetime.now().strftime("%Y-%m-%d %H:%M")
                }
                try:
                    sample_path = os.path.join(os.path.dirname(path), f"{sample_name}.method.json")
                    with open(sample_path, "w") as f:
                        json.dump(sample_method, f, indent=4)
                    logger.log(f"[Method] Auto-saved method file for sample: {sample_path}")
                    messagebox.showinfo("Method File Saved", f"A sample method file was saved:\n{os.path.basename(sample_path)}")

                except Exception as e:
                    logger.log(f"[WARNING] Failed to save method file for '{sample_name}': {e}")
    def extract_title_from_metadata(json_path):
        try:
            with open(json_path, "r") as f:
                meta = json.load(f)
            return meta.get("Experiment Title")
        except:
            return None

    def select_files_generic(filetype_key, allow_multiple=False, on_select_callback=None):
        filetypes_dict = {
            "csv": [("CSV files", "*.csv")],
            "excel": [("Excel files", "*.xls *.xlsx")],
            "json": [("JSON files", "*.json")],
            "all": [("All files", "*.*")]
        }
        filetypes = filetypes_dict.get(filetype_key, filetypes_dict["all"])
        if allow_multiple:
            filepaths = filedialog.askopenfilenames(filetypes=filetypes)
        else:
            filepath = filedialog.askopenfilename(filetypes=filetypes)
            filepaths = [filepath] if filepath else []
        if filepaths and on_select_callback:
            on_select_callback(filepaths)

    def update_status_display(exp_name):
        # Update experiment-level path display
        path_label, sample_label = experiment_status_labels.get(exp_name, (None, None))
        if path_label:
            path = experiment_method_paths.get(exp_name)
            path_label.set(f"Experiment method file: {path if path else 'None'}")

        if sample_label:
            sample_label.set(f"Sample method will be saved at: {sample_method_folder if sample_method_folder else 'None'}")

    #drag
    def on_drag_start(event):
        item_id = tree.identify_row(event.y)
        print("[DEBUG] Drag start:", drag_data)

        if not item_id:
              return
        item_text = tree.item(item_id, "text")
    
        if ":" in item_text:  # This is a file node
            sample_id = tree.parent(item_id)
            exp_id = tree.parent(sample_id)
            if not sample_id or not exp_id:
                return  # Avoid broken context
        
            ft_raw = item_text.split(":")[0].strip()
            drag_data["filetype"] = normalize_ftype(ft_raw)              # <— was .lower()
            drag_data["filename"] = item_text.split(":")[1].strip()
            drag_data["item"] = item_id
            #drag_data["filetype"] = item_text.split(":")[0].strip().lower()
            #drag_data["filename"] = item_text.split(":")[1].strip()
            drag_data["from_sample"] = clean_sample_name(tree.item(sample_id, "text"))
            drag_data["from_exp"] = tree.item(exp_id, "text").replace("Experiment: ", "")
        else:
            drag_data["item"] = None

            
    def on_drag_release(event):
        dest_id = tree.identify_row(event.y)
        if not dest_id or not drag_data["item"]:
            return

        dest_text = tree.item(dest_id, "text")
        if not dest_text.startswith("Sample:") and not any(dest_text.startswith(sym + " Sample:") for sym in ["✅", "⚠️", "❌", "⛔"]):
            return  # Only allow drop into sample

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


    # --- Right-click move logic ---
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

        filetype_raw = selected_text.split(":")[0].strip().lower()
        filetype = normalize_ftype(filetype_raw)                     # <— normalize before use
        filename = selected_text.split(":")[1].strip()
        #filetype = selected_text.split(":")[0].strip().lower()
        #filename = selected_text.split(":")[1].strip()

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
        if selected_text.startswith("Sample: "):
            parent_id = tree.parent(item_id)
            sample_name = selected_text.replace("Sample: ", "").split(" (")[0]
            exp_name = tree.item(parent_id, "text").replace("Experiment: ", "")
            menu.add_command(
                label="Link and Validate Sample",
                command=lambda: link_and_validate_sample(exp_name, sample_name)
            )
        filetype = selected_text.split(":")[0].strip().lower()
        filename = selected_text.split(":")[1].strip()
        # Add remove option if it's a valid file
        menu.add_command(label=f"Remove {filetype.upper()}",command=lambda: remove_file(exp_name, clean_sample_name(sample_name), filetype))  
        menu.post(event.x_root, event.y_root)
          



    def move_file(from_exp, from_sample, ftype, filename, to_exp, to_sample):
        
        entry = experiment_projects[from_exp]["samples"][from_sample][ftype]
        if entry and os.path.basename(entry) == filename:
            experiment_projects[from_exp]["samples"][from_sample][ftype] = None
            assign_file(ftype, entry, to_exp, to_sample)

    # --- write to method ---
    def write_method_file(exp_name, auto=False):
        method = {
            "experiment": exp_name,
            "generated_on": datetime.now().strftime("%Y-%m-%d %H:%M"),
            "samples": {}
        }

        for sname, files in experiment_projects[exp_name]["samples"].items():
            if not all([files.get("csv"), files.get("excel"), files.get("json")]):
                continue
            try:
                with open(files["json"], "r") as f:
                    meta = json.load(f)
            except Exception:
                meta = {}

            entry = {
                "csv": os.path.abspath(os.path.normpath(files["csv"])),#files["csv"],  #os.path.basename(files["csv"]),
                "excel": os.path.abspath(os.path.normpath(files["excel"])),#files["excel"], #os.path.basename(files["excel"]),
                "metadata": os.path.abspath(os.path.normpath(files["json"])),#files["json"], #os.path.basename(files["json"]),
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
                    sample_path = os.path.join(os.path.dirname(files["json"]), f"{sname}.method.json")

                with open(sample_path, "w") as sf:
                    json.dump(sample_method, sf, indent=4)
                logger.log(f"[Method] Per-sample method saved: {sample_path}")
            except Exception as e:
                logger.log(f"[WARNING] Failed to save per-sample method for {sname}: {e}")

        # Save experiment-level .exp.json
        if auto and exp_name in experiment_method_paths:
            output_path = experiment_method_paths[exp_name]
        else:
            output_path = filedialog.asksaveasfilename(
                defaultextension=".exp.json",
                initialfile=f"{exp_name}.exp.json",
                filetypes=[("Experiment Method JSON", "*.exp.json")]
            )
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

    # --- load method file ---

    # --- safe loader ---
    def safe_get_field(d, key, fallback="(not linked)"):
        val = d.get(key)
        return val if isinstance(val, str) and val.strip() else fallback
    #def check_relative_location():
    def load_method_file(paths=None):
        if not paths:
            paths = filedialog.askopenfilenames(
                title="Select One or More Sample Method Files",
                filetypes=[("Sample Method JSON", "*.method.json"), ("JSON files", "*.json")]
            )
        if not paths:
            return

        loaded = 0
        skipped = 0
        for path in paths:
            try:
                with open(path, "r") as f:
                    method = json.load(f)
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load method file:\n{path}\n{e}")
                continue

            exp_name = method.get("experiment", "Recovered")
            samples = method.get("samples", {})
            if not samples:
                skipped += 1
                continue

            if exp_name not in experiment_projects:
                experiment_projects[exp_name] = {"samples": {}}

            for sample_name, files in samples.items():
                sample_name = clean_sample_name(sample_name)  # Normalize sample name

                if sample_name in experiment_projects[exp_name]["samples"]:
                    skipped += 1
                    continue

                base = os.path.dirname(path)
                files_resolved = {
                    "csv": os.path.join(base, files.get("csv")),
                    "excel": os.path.join(base, files.get("excel")),
                    "json": os.path.join(base, files.get("metadata"))
                }
                experiment_projects[exp_name]["samples"][sample_name] = files_resolved

                if files.get("validated"):
                    linked_validated_samples.add((exp_name, sample_name))
                    
                loaded += 1

        refresh_tree()
        msg = f"Imported {loaded} sample(s) successfully.\nSkipped: {skipped}"
        if skipped > 0:
            msg += "\n(Skipped files were already loaded or invalid)"
        messagebox.showinfo("Method Import", msg)

    def change_experiment_method_path():
        path = filedialog.asksaveasfilename(
            title="Select path to save experiment method file",
            defaultextension=".exp.json",
            filetypes=[("Experiment Method JSON", "*.exp.json")]
        )
        if path:
            experiment_method_paths["Unassigned"] = path  # Replace key if you're in a real experiment context
            update_status_display("Unassigned")

    def change_sample_method_folder():
        global sample_method_folder
        path = filedialog.askdirectory(title="Select folder to save sample method files")
        if path:
            sample_method_folder = path
            update_status_display("Unassigned")

    #20250917 to avoid hard-fixing the items in exp json so the PL workflow can be saved as well
    def load_experiment_method_file():
        """Repurposed: load an experiment (.exp.json) and rebuild the tree."""
        path = filedialog.askopenfilename(
            title="Open Experiment (.exp.json)",
            filetypes=[("Experiment JSON", "*.exp.json"), ("JSON", "*.json")]
        )
        if not path:
            return
        exp_title = import_experiment_json(path)
        if exp_title in experiment_status_labels:
            experiment_status_labels[exp_title].config(text=f"EXP file: {path}")

    """
    def load_experiment_method_file():
        path = filedialog.askopenfilename(
            title="Load experiment method file",
            filetypes=[("Experiment Method JSON", "*.exp.json")]
        )
        if not path:
            return

        try:
            with open(path, "r") as f:
                method = json.load(f)
        except Exception as e:
            messagebox.showerror("Load Failed", f"Could not load file:\n{e}")
            return

        exp_name = method.get("experiment", "Recovered")
        if exp_name not in experiment_projects:
            experiment_projects[exp_name] = {"samples": {}}

        experiment_method_paths[exp_name] = path
        for sample_name, files in method.get("samples", {}).items():
            sample = {
                "csv": os.path.join(os.path.dirname(path), files["csv"]),
                "excel": os.path.join(os.path.dirname(path), files["excel"]),
                "json": os.path.join(os.path.dirname(path), files["metadata"])
            }
            if "ion_sheet" in files:
                sample["ion_sheet"] = files["ion_sheet"]
            if "ion_sheet_file" in files:
                sample["ion_sheet_file"] = files["ion_sheet_file"]
            experiment_projects[exp_name]["samples"][sample_name] = sample

            if files.get("validated") is True:
                linked_validated_samples.add((exp_name, sample_name))
            elif files.get("validated") is False:
                validation_failed_samples.add((exp_name, sample_name))

        experiment_status_labels[exp_name] = (exp_status_var, sample_status_var)
        refresh_tree()
        update_status_display(exp_name)
        messagebox.showinfo("Loaded", f"Method file loaded for experiment: {exp_name}")
    """
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
            path = filedialog.asksaveasfilename(
                title="Save Experiment (.exp.json)",
                defaultextension=".exp.json",
                initialfile=f"{safe}.exp.json",
                filetypes=[("Experiment JSON", "*.exp.json"), ("JSON", "*.json")]
            )
            if not path:
                return
        export_experiment_json(exp, path)
    #20250917 sarabada
    """ 
    def save_current_experiment_method():
        sel = tree.selection()
        if not sel:
            messagebox.showwarning("No Experiment Selected", "Please select an experiment to save.")
            return

        sample_id = sel[0]
        exp_id = tree.parent(sample_id)
        if not exp_id:  # If we're already on an experiment node
            exp_id = sample_id

        exp_name = tree.item(exp_id, "text").replace("Experiment: ", "").split(" (")[0].strip()
        write_method_file(exp_name, auto=True)
    """

    def bind_right_click(widget, callback):
        # Universal right-click binding for macOS, Windows, Linux
        widget.bind("<Button-3>", callback)  # Windows & Linux
        widget.bind("<Control-Button-1>", callback)  # macOS trackpad

    def get_selected_csv_path():
        selected = tree.focus()
        print(selected)
        node_info = experiment_projects.get(selected)
        print(f"node info: {node_info}")
        if node_info and node_info["type"] == "csv":
            exp = node_info["exp"]
            sample = node_info["sample"]
            path = experiment_projects[exp]["samples"][sample]["csv"]
            return path
        else:
            return None

        # Reuse the metadata stored in TreeView tags or descriptions
        #item_info = tree.item(selected_item)
        #values = item_info.get("values", [])
        
        # Sample: ["File", "Unassigned", "Sample_XYZ", "csv"]
        #if len(values) >= 4:
        #    node_type = values[0].lower()
        #    exp = values[1]
        #    sample = values[2]
        #    filetype = values[3].lower()
        #    
        #    if node_type in ("sample", "file") and filetype == "csv":
        #        return experiment_projects[exp]["samples"][sample].get("csv")

        return None

    # -- pseudo labeling --
    def launch_pseudo_labeling():
        #import compnewv4 as compv4  # assumes dev/test calls are guarded by if __name__ == "__main__"
        from datetime import datetime

        sel = tree.selection()
        if not sel:
            messagebox.showwarning("No Selection", "Select a sample in the tree first.")
            return

        node = sel[0]
        text = tree.item(node, "text")

        # If user clicked a file row like "CSV: ..." or "Metadata: ...", go up to the sample row
        if ":" in text:
            node = tree.parent(node)

        sample_name = clean_sample_name(tree.item(node, "text"))
        exp_node = tree.parent(node)
        if not exp_node:
            messagebox.showerror("Invalid Selection", "Please select a sample under an experiment.")
            return

        exp_text = tree.item(exp_node, "text")
        exp_name = exp_text.replace("Experiment: ", "").split(" (")[0].strip()

        files = experiment_projects.get(exp_name, {}).get("samples", {}).get(sample_name, {})
        csv_path = files.get("csv")
        meta_path = files.get("json")

        if not csv_path or not meta_path:
            messagebox.showerror("Missing Files", "This sample must have both CSV and Metadata (.json) linked.")
            return
        
        print("[launch] files:", files)
        print("[launch] using meta:", meta_path)

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
                    messagebox.showwarning("Metadata not valid",
                                        "Selected file lacks Glycan Type / Mass Analyzer charge mode.")
                    meta_path, meta_dict = None, None
            else:
                messagebox.showwarning(
                    "Metadata not found",
                    "Could not locate metadata for this sample. You can proceed, but defaults may be wrong."
                )

        print("[launch] using meta (final):", meta_path)


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
            outname = f"{sample_name}_insilico_{datetime.now().strftime('%Y%m%d')}.csv"
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
                refresh_tree()
                return outpath  # so the window can show it immediately

            except Exception as e:
                traceback.print_exc()
                messagebox.showerror("Generation failed", str(e))
                return None

        def on_link_existing(path):
            files["insilico_csv"] = path
            append_runlog(files, {
                "ts": datetime.now().isoformat(timespec="seconds"),
                "action": "insilico_link_existing",
                "sample": sample_name,
                "output": {"insilico_csv": path},
            })
            refresh_tree()

        def on_attach_ionlist(path):
            files["ionlist_path"] = path
            refresh_tree()

        def on_start(payload):
            # Use linked paths & run the wrapper
            run_pseudolabeling(
                sample_name=sample_name,
                files=files,
                meta_overrides=payload.get("metadata", {}),
                parent=root
            )

        # Open the setup window with the resolved metadata
        PseudoLabelingSetupWindow(
            root,
            meta_json_path=meta_path,
            meta_prefill=meta_dict or {},
            default_flags=compv4.NG_flags,
            on_generate=on_generate,                 # new
            on_link_existing=on_link_existing,       # new
            on_attach_ionlist=on_attach_ionlist,     # new
            on_start=on_start,                       # new
            initial_insilico=files.get("insilico_csv"),
            initial_ionlist=files.get("ionlist_path")
        )

    #20250911 
    def open_pl_to_trainable_modal(root, sample_name, files, logger):
        import tkinter as tk
        from tkinter import ttk, filedialog, messagebox

        win = tk.Toplevel(root)
        win.title("Pseudolabel → Trainable (one-pass)")
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
        converted_csv_var        = tk.StringVar(value=files.get("converted_csv",""))

        def browse(var, exts=(("All","*.*"),)):
            p = filedialog.askopenfilename(filetypes=exts)
            if p: var.set(p)

        row=0
        ttk.Label(frm, text="1) Pseudolabeled TSV/CSV (long):").grid(row=row, column=0, sticky="w"); 
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

        mode = tk.StringVar(value="rebuild")
        ttk.Radiobutton(feat_box, text="Rebuild ALL features from long-form peaks (log10(1+I))",
                        variable=mode, value="rebuild").grid(row=0, column=0, columnspan=3, sticky="w")
        ttk.Radiobutton(feat_box, text="Reuse existing wide features CSV",
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

        """
        # --- Features
        feat_box = ttk.LabelFrame(frm, text="6) Feature building")
        feat_box.grid(row=row, column=0, columnspan=3, sticky="we", pady=(4,8))
        mode = tk.StringVar(value="rebuild")
        ttk.Radiobutton(feat_box, text="Rebuild ALL features from long-form peaks (log10(1+I))", variable=mode, value="rebuild").grid(row=0, column=0, columnspan=3, sticky="w")
        ttk.Radiobutton(feat_box, text="Reuse existing wide features CSV", variable=mode, value="reuse").grid(row=1, column=0, columnspan=3, sticky="w")
        ttk.Entry(feat_box, textvariable=wide_var, width=70).grid(row=2, column=0, sticky="we")
        ttk.Button(feat_box, text="Choose…", command=lambda: browse(wide_var,(("CSV","*.csv"),))).grid(row=2, column=1, sticky="w")
        # --- Extra
        mass_box = ttk.LabelFrame(frm, text="7) Include mass as feature")
        mass_box.grid(row=row, column=0, columnspan=3, sticky="we", pady=(4,12))
        ttk.Checkbutton(mass_box, text="Include precursor mass as feature",
                variable=include_mass_feat_var).grid(row=0, column=0, sticky="we")
        """
        #consider the option of passing precursor mass to the trainable csv
        """
        # --- Precursor gate (NEW) ---
        gate_frm = ttk.LabelFrame(frm, text="Precursor (MS1) Gate", padding=(10,8))
        gate_frm.grid(row=row, column=0, columnspan=3, sticky="we", pady=(6,6))
        row += 1

        chk = ttk.Checkbutton(gate_frm, text="Apply precursor gate", variable=apply_precursor_gate_var)
        chk.grid(row=0, column=0, sticky="w")

        ttk.Label(gate_frm, text="Tolerance (ppm):").grid(row=0, column=1, padx=(16,4), sticky="e")
        ppm_ent = ttk.Entry(gate_frm, width=8, textvariable=precursor_ppm_var)
        ppm_ent.grid(row=0, column=2, sticky="w")

        ttk.Label(gate_frm, text="Converted CSV/TSV (has MS2scan_no + protonatedmass):").grid(
            row=1, column=0, columnspan=2, sticky="w", pady=(6,0)
        )
        conv_ent = ttk.Entry(gate_frm, width=60, textvariable=converted_csv_var)
        conv_ent.grid(row=2, column=0, columnspan=2, sticky="we", pady=(2,6))
        def _pick_converted():
            p = filedialog.askopenfilename(
                title="Select converted CSV/TSV",
                filetypes=[("Table files","*.csv *.tsv *.txt"), ("All files","*.*")]
            )
            if p:
                converted_csv_var.set(p)
        ttk.Button(gate_frm, text="Browse…", command=_pick_converted).grid(row=2, column=2, padx=(8,0), sticky="w")
        gate_frm.columnconfigure(0, weight=1)
        """
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
                print("[PL→Train][dbg] include_mass_feature UI:", include_mass_feat_var.get())
                outpath, summary = build_trainable_from_pseudolabels(
                    sample_name=sample_name,
                    pseudo_path=pseudo_var.get().strip(),
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
            except Exception as e:
                import traceback; traceback.print_exc()
                messagebox.showerror("Failed", str(e))

        ttk.Button(frm, text="Build Trainable CSV", command=run_once).grid(row=row+2, column=1, pady=8)

    def try_pl_to_trainable():
        sel = tree.selection()
        if not sel:
            messagebox.showerror("No Selection", "Please select a sample first.")
            return
        sample_node = sel[0]
        sample_name = clean_sample_name(tree.item(sample_node, "text"))
        exp_node = tree.parent(sample_node)
        exp_name = tree.item(exp_node, "text").replace("Experiment: ", "").split(" (")[0].strip()

        files = experiment_projects[exp_name]["samples"][sample_name]
        open_pl_to_trainable_modal(root, sample_name, files, logger)

    # --- Button panel ---
    button_frame = tk.Frame(subwin)
    button_frame.pack(pady=5)

    tk.Button(button_frame, text="Add MS2 CSV(s)", command=lambda: select_files_generic("csv", True, handle_csv_selection)).grid(row=0, column=0, padx=5)
    tk.Button(button_frame, text="Add Excels (ion list/man annotation)", command=lambda: select_files_generic("excel", True, handle_excel_selection)).grid(row=0, column=1, padx=5)
    tk.Button(button_frame, text="Add Sample using metadata", command=lambda: select_files_generic("json", True, handle_json_selection)).grid(row=0, column=2, padx=5)
    tk.Button(button_frame, text="Add Sample", command=add_sample).grid(row=1, column=0, padx=5)
    tk.Button(button_frame, text="Clean up empty unassigned sample tags", command=clean_unassigned_samples).grid(row=1, column=1, padx=5)
    link_button = tk.Button(button_frame, text="Link Sample", state="disabled", command=lambda: try_link_selected_sample())
    link_button.grid(row=1, column=2, padx=5) #why it was gone?
    merge_button = tk.Button(button_frame, text="Merge Sample", state="disabled", command=lambda: try_merge_selected_sample())
    merge_button.grid(row=1, column=3, padx=5)
    tk.Button(button_frame, text="Negative options…",
          command=open_negative_options_dialog).grid(row=2, column=0, padx=5)
    tk.Button(button_frame, text="Ion suggestions…",
          command=open_ion_suggest_dialog).grid(row=2, column=1, padx=5)
    tk.Button(button_frame, text="View suggestions…",
          command=open_ion_suggestions_viewer).grid(row=2, column=2, padx=5)
    tk.Button(button_frame, text="Load Method", command=load_method_file).grid(row=2, column=3, padx=5)
    ttk.Button(button_frame, text="CGA manager", command=lambda:launch_pseudo_labeling()).grid(row=3, column=0, padx=5, pady=5) 
    #GPT said without () it only passes the function, and work only if clicked
    tk.Button(button_frame, text="CGA → Trainable",
          command=try_pl_to_trainable).grid(row=3, column=1, padx=5)
    #ttk.Button(button_frame, text="Build trainable CSV from pseudolabeled TSV",
    #       command=lambda: open_pl_to_trainable_modal(root, sample_name, files, logger)).pack(pady=6)

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


    # --- Right-click bind ---
    #tree.bind("<Button-3>", on_right_click)
    tree.bind("<ButtonPress-1>", on_drag_start)
    tree.bind("<ButtonRelease-1>", on_drag_release)
    bind_right_click(tree, on_right_click)


def open_ml_analysis_window():
    import json
    import os
    import pandas as pd
    from tkinter import filedialog, messagebox
    from tkinter import ttk
    import tkinter as tk
    #20250901 add split
    from sklearn.model_selection import train_test_split
    import copy

    #v0.9923~0.9929 utilities
    from ml_ng_utils import (
        collect_ng_candidates,
        cap_non_glycan,
        predict_with_threshold,
        build_features_from_peaks_log10_plus1,  # optional if you need it directly
    )
    from ml_ng_utils_extras import resample_by_strategy, tau_sweep_summary
    #_ensure_ml_state()
    #20250919 try to let live param changes apply to save/loadable field
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
                "test_size": 0.25,
                "val_size": 0.10,
                "stratify": True,
                "real_world_test": False,
                "threshold": {"enabled": False, "tau": 0.65, "margin": 0.05},
            },
            "min_samples_per_class": 5,
            "balance": {"enabled": True, "majority_label": "Non-glycan", "majority_factor": 3},
        },
        # These are set by open_ml_params_window so Train/Test can push live updates:
        "editor_txt": None,
        "preview_txt": None,
    }

    def _render_effective_from_state():
        if _ml_state["preview_txt"] is None:
            return
        # If you already have a merge routine, use it; otherwise just echo current:
        eff = _ml_state["current"]  # or: _merge_ml(BUILTIN_ML, files_layer, _ml_state["current"])
        _ml_state["preview_txt"].delete("1.0", "end")
        _ml_state["preview_txt"].insert("1.0", json.dumps(eff, indent=2))

    def _write_editor_from_state():
        """Refresh the JSON text panes if they exist (no-op otherwise)."""
        try:
            import json
        except Exception:
            json = None

        ed, pv = _get_editor_preview_widgets()

        editor_blob   = _ml_state.get("editor") or {}
        effective_blob = _ml_state.get("effective_params") or {}

        s_editor    = json.dumps(editor_blob, indent=2) if json else str(editor_blob)
        s_effective = json.dumps(effective_blob, indent=2) if json else str(effective_blob)

        if ed:
            try:
                ed.configure(state="normal")
                ed.delete("1.0", "end")
                ed.insert("1.0", s_editor)
            except Exception:
                pass

        if pv:
            try:
                pv.configure(state="normal")
                pv.delete("1.0", "end")
                pv.insert("1.0", s_effective)
            except Exception:
                pass

    """
    _ml_state = {
        "editor": {   # what the user edits (left JSON)
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
                "test_size": 0.25,
                "val_size": 0.10,
                "stratify": True,
                "real_world_test": False,
                "threshold": {"enabled": False, "tau": 0.65, "margin": 0.05},
                "min_samples_per_class": 5,
            },
            "balance": {"enabled": True, "majority_label": "Non-glycan", "majority_factor": 3}
        },
        # optional: layers coming from files/experiment; if none, keep empty list
        "file_layers": [],
        # widgets hooks are filled later when windows are created
        "widgets": {"editor_txt": None, "preview_txt": None}
    }
    

    def _write_train_ui_into_state(
        test_split, val_split, min_samples, use_balance,
        majority_label, majority_factor, stratify,
        n_estimators, class_weight_balanced, real_world_test,
        th_enable, th_tau, th_margin
    ):
        ed = _ml_state["editor"]

        # model
        ed.setdefault("model", {})
        ed["model"]["type"] = "RandomForest"
        ed["model"]["n_estimators"] = int(n_estimators)
        ed["model"]["class_weight"] = "balanced" if class_weight_balanced else None

        # train
        ed.setdefault("train", {})
        ed["train"]["test_size"] = float(test_split)
        ed["train"]["val_size"] = float(val_split)
        ed["train"]["stratify"] = bool(stratify)
        ed["train"]["real_world_test"] = bool(real_world_test)
        ed["train"]["min_samples_per_class"] = int(min_samples)
        ed["train"]["threshold"] = {
            "enabled": bool(th_enable),
            "tau": float(th_tau),
            "margin": float(th_margin),
        }

        # balance
        ed["balance"] = {
            "enabled": bool(use_balance),
            "majority_label": majority_label,
            "majority_factor": int(majority_factor),
        }

        _refresh_json_textboxes()
    """
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
                t2.insert("1.0", json.dumps(_ml_state.get("effective", _ml_state["editor"]), indent=2))
                t2.configure(state="normal")
        except Exception:
            pass

    def _write_train_ui_into_state( 
        test_size: float,
        val_size: float,
        min_per_class: int,
        use_balance: bool,
        majority_label: str,
        majority_factor: int,
        stratify: bool,
        n_estimators: int,
        use_class_weight: bool,
        real_world_test: bool,
        thr_enable: bool,
        thr_tau: float,
        thr_margin: float,
    ):
        st = _ensure_ml_state()
        ed  = st["editor"]
        tr  = ed.setdefault("train", {})
        mdl = ed.setdefault("model", {})
        bal = ed.setdefault("balance", {})

        # train/test
        tr["test_size"] = float(test_size)
        tr["val_size"]  = float(val_size)
        tr["stratify"]  = bool(stratify)
        tr["real_world_test"] = bool(real_world_test)
        tr["threshold"] = {"enabled": bool(thr_enable), "tau": float(thr_tau), "margin": float(thr_margin)}

        # per-class minimum
        ed["min_samples_per_class"] = int(min_per_class)

        # balance
        bal["enabled"] = bool(use_balance)
        bal["majority_label"]  = str(majority_label).strip() if majority_label is not None else "Non-glycan"
        bal["majority_factor"] = int(majority_factor)

        # model (RF)
        mdl["type"]          = "RandomForest"
        mdl["n_estimators"]  = int(n_estimators)
        mdl["class_weight"]  = "balanced" if bool(use_class_weight) else None

        # recompute effective if you keep one
        st["effective"] = copy.deepcopy(ed)

    def _deep_merge(a, b):
        if isinstance(a, dict) and isinstance(b, dict):
            out = dict(a)
            for k, v in b.items():
                out[k] = _deep_merge(out.get(k), v) if k in out else v
            return out
        return b if b is not None else a

    def _effective_params():
        eff = BUILTIN_ML
        for layer in _ml_state.get("file_layers", []):
            eff = _deep_merge(eff, layer)
        eff = _deep_merge(eff, _ml_state["editor"])
        return eff

    def _refresh_json_textboxes():
        ed = json.dumps(_ml_state["editor"], indent=2)
        eff = json.dumps(_effective_params(), indent=2)
        t1 = _ml_state["widgets"]["editor_txt"]
        t2 = _ml_state["widgets"]["preview_txt"]
        if t1:
            t1.config(state="normal"); t1.delete("1.0", "end"); t1.insert("1.0", ed); t1.config(state="normal")
        if t2:
            t2.config(state="normal"); t2.delete("1.0", "end"); t2.insert("1.0", eff); t2.config(state="normal")

    #20250929
    def _drop_meta_and_get_X(df: pd.DataFrame) -> pd.DataFrame:
        """
        Returns a numeric-only feature matrix X.
        Drops obviously non-feature metadata columns like UID and Origin_*.
        Also attempts to coerce object columns to numeric if they are mostly numeric.
        """
        import numpy as np
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

        # try to coerce mostly numeric object columns
        for c in list(df2.columns):
            if df2[c].dtype == object:
                coerced = pd.to_numeric(df2[c], errors="coerce")
                # keep coerced if it's largely numeric (>=95% not NaN)
                if coerced.notna().mean() >= 0.95:
                    df2[c] = coerced

        # keep strictly numeric columns for X
        X = df2.select_dtypes(include=[np.number])
        return X    

    import re

    # Pure manual: FHNSGKDN order, no neg-mode extras
    _MANUAL_PURE = re.compile(r"^(?:F\d+)?H\d+N\d+(?:S\d+)?(?:G\d+)?(?:KDN\d+)?$", re.IGNORECASE)

    # Manual + PL extensions allowed:
    #   optional F..., then H...N..., optional S/G/KDN, optional A..., optional trailing s/p suffixes (lowercase)
    _MANUAL_EXT = re.compile(
        r"^(?:F\d+)?H\d+N\d+(?:S\d+)?(?:G\d+)?(?:KDN\d+)?(?:A\d+)?(?:[sp]\d+)?$",
        re.IGNORECASE
    )

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


    #20250909
    def load_model_any(model_path: str):
        ext = os.path.splitext(model_path)[1].lower()
        if ext in (".joblib", ".pkl"):
            model = joblib.load(model_path)
            loader = "joblib"
        elif ext == ".skops":
            try:
                from skops.io import load as sk_load
            except Exception as e:
                raise ImportError(
                    "This model is a .skops file but 'skops' is not installed. "
                    "Install it in this environment: pip install skops"
                ) from e
            model = sk_load(model_path, trusted=True)
            loader = "skops"
        else:
            raise ValueError(f"Unsupported model file extension: {ext}")
        if not hasattr(model, "predict"):
            raise TypeError(
                f"Loaded object is {type(model).__name__} and has no .predict(). "
                "Did you select the *_labelencoder.joblib by mistake?"
            )
        return model, loader

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


    def balance_and_split(
        df: pd.DataFrame,
        label_col="Structure",
        majority_label="Non-glycan",
        min_count=12,
        majority_factor=3,
        feature_exclude=("MS2scan_no","ID","Source","IUPACname(optional)","Glycanannotation2","GlyToucan ID","unique_ID"),
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
        from pathlib import Path
        folder = os.path.dirname(predict_input_path)

        # 1) explicit hint column
        for col in df.columns:
            if col.lower() in ("method_json", "method_path", "method", "json_path"):
                try:
                    cand = str(df[col].dropna().iloc[0]).strip()
                    if cand and cand.lower().endswith(".method.json") and os.path.exists(cand):
                        return Path(cand).stem
                except Exception:
                    pass

        # 2) neighbors
        candidates = sorted(glob.glob(os.path.join(folder, "*.method.json")))
        if len(candidates) == 1:
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
                    return Path(c).stem

            # 4) ask user explicitly
            try:
                from tkinter import filedialog
                sel = filedialog.askopenfilename(
                    title="Select method.json for packaging",
                    initialdir=folder,
                    filetypes=[("Method JSON", "*.method.json")]
                )
                if sel:
                    return Path(sel).stem
            except Exception:
                pass  # fall through to #5

        # 5) last resort
        return Path(predict_input_path).stem
    
    #new safe method json & raw csv
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





    #20250918 ML method export utilities
    # =========================
    # ML PARAM PRESET UTILITIES
    # =========================
    #from __future__ import annotations
    import json, os, hashlib, datetime, platform, sys
    from typing import Optional, Dict, Any

    # Sensible built-ins if neither exp nor method has overrides
    BUILTIN_ML: Dict[str, Any] = {
        "model": {"type": "RandomForest", "n_estimators": 500, "max_depth": None, "class_weight": "balanced"},
        "split": {"test_size": 0.2, "val_size": 0.0, "random_state": 42, "stratified": True},
        "filters": {"min_samples_per_class": 15, "drop_rare_in_test": True},
        "negatives": {"enabled": True, "max_ratio": 3.0},
        "features": {"exclude_cols": ["MS2scan_no", "protonatedmass"], "use_ion_suggestions": False},
        "thresholds": {"tau": 0.60, "margin": 0.05},
    }

    def _safe_read_json(path: str) -> Optional[Dict[str, Any]]:
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return None

    def _safe_write_json(path: str, data: Dict[str, Any]) -> None:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
            f.write("\n")

    def read_ml_defaults_from_exp(exp_json_path: str) -> Optional[Dict[str, Any]]:
        data = _safe_read_json(exp_json_path)
        if not data: return None
        return data.get("ml_defaults")

    def write_ml_defaults_to_exp(exp_json_path: str, params: Dict[str, Any]) -> bool:
        data = _safe_read_json(exp_json_path) or {}
        data["ml_defaults"] = params
        try:
            _safe_write_json(exp_json_path, data)
            return True
        except Exception:
            return False

    def read_ml_overrides_from_method(method_json_path: str) -> Optional[Dict[str, Any]]:
        data = _safe_read_json(method_json_path)
        if not data: return None
        return data.get("ml_overrides")

    def write_ml_overrides_to_method(method_json_path: str, params: Dict[str, Any]) -> bool:
        data = _safe_read_json(method_json_path) or {}
        data["ml_overrides"] = params
        try:
            _safe_write_json(method_json_path, data)
            return True
        except Exception:
            return False

    def load_ml_params_from_json(path: str) -> Optional[Dict[str, Any]]:
        return _safe_read_json(path)

    def save_ml_params_to_json(path: str, params: Dict[str, Any]) -> bool:
        try:
            _safe_write_json(path, params)
            return True
        except Exception:
            return False

    def _deep_merge(a: Dict[str, Any], b: Dict[str, Any]) -> Dict[str, Any]:
        """Return new dict = a merged with b (b overrides)."""
        out = dict(a)
        for k, v in (b or {}).items():
            if isinstance(v, dict) and isinstance(out.get(k), dict):
                out[k] = _deep_merge(out[k], v)
            else:
                out[k] = v
        return out

    def merge_ml_params(builtins: Dict[str, Any],
                        exp_defaults: Optional[Dict[str, Any]],
                        method_overrides: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        return _deep_merge(_deep_merge(builtins, exp_defaults or {}), method_overrides or {})

    # ---- ML state bootstrapper -----------------------------------------------
    def _ensure_ml_state():
        """Guarantee _ml_state exists with the keys we expect."""
        global _ml_state
        try:
            _ml_state  # noqa: F401
        except NameError:
            _ml_state = {}

        if not isinstance(_ml_state, dict):
            _ml_state = {}

        # BUILTIN_ML must already exist here
        if "editor" not in _ml_state:
            _ml_state["editor"] = copy.deepcopy(BUILTIN_ML)
        if "effective" not in _ml_state:
            _ml_state["effective"] = copy.deepcopy(_ml_state["editor"])
        if "source_files" not in _ml_state:
            _ml_state["source_files"] = []  # paths you load/merge from

        return _ml_state

    def file_sha256(path: str) -> Optional[str]:
        try:
            h = hashlib.sha256()
            with open(path, "rb") as f:
                for chunk in iter(lambda: f.read(1024*1024), b""):
                    h.update(chunk)
            return h.hexdigest()
        except Exception:
            return None

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
        _safe_write_json(out, payload)
        return out

    # ==================================
    # TK PANEL: ML Parameter Editor Pane
    # ==================================
    def build_ml_params_panel(parent,
                            get_current_context,
                            on_effective_params_ready):
        """
        parent: a tk/ttk container in your Train tab.
        get_current_context(): callable -> dict with keys:
            {
            "exp_json": str|None,
            "method_json": str|None,
            "current_params": dict|None   # if you already have a working set in memory
            }
        on_effective_params_ready(params_dict): called when user clicks "Validate & Use"
        """
        import tkinter as tk
        from tkinter import ttk, filedialog, messagebox
        from tkinter.scrolledtext import ScrolledText

        ctx = get_current_context() or {}
        exp_json_path = ctx.get("exp_json")
        method_json_path = ctx.get("method_json")
        working_params = ctx.get("current_params") or {}

        frm = ttk.LabelFrame(parent, text="ML Parameters", padding=8)
        frm.grid_columnconfigure(0, weight=1)
        frm.grid_columnconfigure(1, weight=1)
        frm.grid(row=0, column=0, sticky="nsew", padx=4, pady=4)

        # Editable JSON text
        ttk.Label(frm, text="Editable (this panel's params):").grid(row=0, column=0, sticky="w")
        edit_box = ScrolledText(frm, height=14, wrap="none")
        edit_box.grid(row=1, column=0, sticky="nsew", padx=(0,4))
        # Effective (preview after merge)
        ttk.Label(frm, text="Effective (Builtins ← Exp ← Method):").grid(row=0, column=1, sticky="w")
        eff_box = ScrolledText(frm, height=14, wrap="none", state="disabled")
        eff_box.grid(row=1, column=1, sticky="nsew")

        def _set_edit_box_from_dict(d):
            edit_box.delete("1.0", "end")
            edit_box.insert("1.0", json.dumps(d or {}, indent=2, ensure_ascii=False))

        def _set_eff_box_from_dict(d):
            eff_box.configure(state="normal")
            eff_box.delete("1.0", "end")
            eff_box.insert("1.0", json.dumps(d or {}, indent=2, ensure_ascii=False))
            eff_box.configure(state="disabled")

        def _current_edit_dict():
            try:
                return json.loads(edit_box.get("1.0", "end").strip() or "{}")
            except Exception as e:
                messagebox.showerror("Invalid JSON", f"Could not parse parameters:\n{e}")
                return None

        def _refresh_effective_preview():
            exp_def = read_ml_defaults_from_exp(exp_json_path) if exp_json_path else None
            meth_ov = read_ml_overrides_from_method(method_json_path) if method_json_path else None
            # The editor holds the *active layer* user is editing. Choose which layer?
            # By default, we treat editor content as what will be *applied* somewhere later.
            # For preview, just show how it would look if applied at the method level:
            try:
                edited = _current_edit_dict() or {}
            except Exception:
                edited = {}
            effective = merge_ml_params(BUILTIN_ML, exp_def, edited if meth_ov is None else _deep_merge(meth_ov, edited))
            _set_eff_box_from_dict(effective)

        # Initialize editor with either method overrides, else exp defaults, else builtins
        initial = read_ml_overrides_from_method(method_json_path) if method_json_path else None
        if initial is None:
            initial = read_ml_defaults_from_exp(exp_json_path) if exp_json_path else None
        if initial is None:
            initial = working_params or BUILTIN_ML
        _set_edit_box_from_dict(initial)
        _refresh_effective_preview()

        # ---- Buttons row
        btns = ttk.Frame(frm); btns.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(6,0))
        for i in range(6): btns.grid_columnconfigure(i, weight=1)

        def _load_from_file():
            path = filedialog.askopenfilename(title="Load ML params JSON",
                                            filetypes=[("JSON", "*.json"), ("All files","*.*")])
            if not path: return
            d = load_ml_params_from_json(path)
            if d is None:
                messagebox.showerror("Load failed", "Could not read JSON parameters.")
                return
            _set_edit_box_from_dict(d); _refresh_effective_preview()

        def _save_to_file():
            path = filedialog.asksaveasfilename(title="Save ML params JSON",
                                                defaultextension=".json",
                                                filetypes=[("JSON", "*.json"), ("All files","*.*")])
            if not path: return
            d = _current_edit_dict()
            if d is None: return
            ok = save_ml_params_to_json(path, d)
            if not ok:
                messagebox.showerror("Save failed", "Could not save JSON parameters.")
            else:
                messagebox.showinfo("Saved", f"Saved parameters to:\n{path}")

        def _pull_from_exp():
            if not exp_json_path:
                messagebox.showwarning("No experiment file", "No .exp.json is active.")
                return
            d = read_ml_defaults_from_exp(exp_json_path)
            if d is None:
                messagebox.showinfo("No defaults", "This experiment has no ml_defaults yet.")
                return
            _set_edit_box_from_dict(d); _refresh_effective_preview()

        def _pull_from_method():
            if not method_json_path:
                messagebox.showwarning("No method file", "No .method.json is active.")
                return
            d = read_ml_overrides_from_method(method_json_path)
            if d is None:
                messagebox.showinfo("No overrides", "This sample has no ml_overrides yet.")
                return
            _set_edit_box_from_dict(d); _refresh_effective_preview()

        def _apply_to_exp():
            if not exp_json_path:
                messagebox.showwarning("No experiment file", "No .exp.json is active.")
                return
            d = _current_edit_dict()
            if d is None: return
            if write_ml_defaults_to_exp(exp_json_path, d):
                messagebox.showinfo("Applied", f"Updated ml_defaults in:\n{exp_json_path}")
            else:
                messagebox.showerror("Failed", "Could not write to experiment JSON.")
            _refresh_effective_preview()

        def _apply_to_method():
            if not method_json_path:
                messagebox.showwarning("No method file", "No .method.json is active.")
                return
            d = _current_edit_dict()
            if d is None: return
            if write_ml_overrides_to_method(method_json_path, d):
                messagebox.showinfo("Applied", f"Updated ml_overrides in:\n{method_json_path}")
            else:
                messagebox.showerror("Failed", "Could not write to method JSON.")
            _refresh_effective_preview()

        def _reset_builtins():
            _set_edit_box_from_dict(BUILTIN_ML); _refresh_effective_preview()

        def _validate_and_use():
            d = _current_edit_dict()
            if d is None: return
            # Merge for runtime (method layer gets editor content if a method exists)
            exp_def = read_ml_defaults_from_exp(exp_json_path) if exp_json_path else None
            meth_ov = read_ml_overrides_from_method(method_json_path) if method_json_path else None
            effective = merge_ml_params(BUILTIN_ML, exp_def, _deep_merge(meth_ov or {}, d))
            on_effective_params_ready(effective)
            messagebox.showinfo("Parameters ready", "Effective ML parameters are validated and ready to use.")
            _set_eff_box_from_dict(effective)



        ttk.Button(btns, text="Load…", command=_load_from_file).grid(row=0, column=0, sticky="ew", padx=2)
        ttk.Button(btns, text="Save…", command=_save_to_file).grid(row=0, column=1, sticky="ew", padx=2)
        ttk.Button(btns, text="Pull from Experiment", command=_pull_from_exp).grid(row=0, column=2, sticky="ew", padx=2)
        ttk.Button(btns, text="Pull from Sample", command=_pull_from_method).grid(row=0, column=3, sticky="ew", padx=2)
        ttk.Button(btns, text="Apply to Experiment", command=_apply_to_exp).grid(row=0, column=4, sticky="ew", padx=2)
        ttk.Button(btns, text="Apply to Sample", command=_apply_to_method).grid(row=0, column=5, sticky="ew", padx=2)

        btns2 = ttk.Frame(frm); btns2.grid(row=3, column=0, columnspan=2, sticky="ew", pady=(6,0))
        btns2.grid_columnconfigure(0, weight=1); btns2.grid_columnconfigure(1, weight=1)
        ttk.Button(btns2, text="Reset to Built-ins", command=_reset_builtins).grid(row=0, column=0, sticky="ew", padx=2)
        ttk.Button(btns2, text="Validate & Use", command=_validate_and_use).grid(row=0, column=1, sticky="ew", padx=2)

        # Recompute preview on edits
        def _on_edit(*_):
            _refresh_effective_preview()
        edit_box.bind("<<Modified>>", lambda e: (edit_box.edit_modified(False), _on_edit()))

        return frm  # in case you want to pack/place/grid differently from caller

    predict_input_path = None
    train_csv_path = None
    linked_exp_json = None
    # ----- ML params helpers (scoped to ML window) -----
    import json, os
    from tkinter import filedialog, messagebox
    from tkinter.scrolledtext import ScrolledText

    BUILTIN_ML = {
        "model": {"type": "RandomForest", "n_estimators": 400, "max_depth": None,
                  "min_samples_split": 2, "min_samples_leaf": 1, "random_state": 42,
                  "class_weight": "balanced"},
    }

    def _read_ml_defaults_from_exp(path):
        try:
            with open(path, "r", encoding="utf-8") as f: data = json.load(f)
            return data.get("ml_defaults")
        except Exception:
            return None

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

    def _merge_ml(builtins: dict, exp_defaults: dict|None, editor: dict|None) -> dict:
        out = copy.deepcopy(builtins)
        for layer in (exp_defaults or {}, editor or {}):
            for k,v in layer.items():
                if isinstance(v, dict) and isinstance(out.get(k), dict):
                    out[k].update(v)
                else:
                    out[k] = v
        return out

    def build_ml_params_panel(parent, get_current_context, on_effective_params_ready):
        box = ttk.LabelFrame(parent, text="ML Parameters"); box.grid_columnconfigure(0, weight=1); box.grid_columnconfigure(1, weight=1)

        ttk.Label(box, text="Editor (JSON):").grid(row=0, column=0, sticky="w")
        ttk.Label(box, text="Effective (Builtins ← Exp ← Editor):").grid(row=0, column=1, sticky="w")

        editor = ScrolledText(box, height=12, wrap="none"); editor.grid(row=1, column=0, sticky="nsew", padx=(0,6))
        effbox = ScrolledText(box, height=12, wrap="none", state="disabled"); effbox.grid(row=1, column=1, sticky="nsew")




        def _set_eff(d):
            effbox.configure(state="normal"); effbox.delete("1.0","end")
            effbox.insert("1.0", json.dumps(d or {}, indent=2, ensure_ascii=False))
            effbox.configure(state="disabled")

        def _get_edit():
            try:
                txt = editor.get("1.0","end").strip() or "{}"
                return json.loads(txt)
            except Exception as e:
                messagebox.showerror("Invalid JSON", f"Editor JSON parse error:\n{e}")
                return None

        def _refresh():
            ctx = get_current_context() or {}
            exp_path = ctx.get("exp_json_path")
            exp_defs = _read_ml_defaults_from_exp(exp_path) if exp_path else None
            ed = _get_edit() or {}
            eff = _merge_ml(BUILTIN_ML, exp_defs, ed)
            _set_eff(eff)

        # initial load (use exp defaults if available, else builtins)
        ctx0 = get_current_context() or {}
        init_defs = _read_ml_defaults_from_exp(ctx0.get("exp_json_path")) if ctx0.get("exp_json_path") else None
        editor.insert("1.0", json.dumps(init_defs or BUILTIN_ML, indent=2, ensure_ascii=False))
        _refresh()

        # buttons
        row = 2
        btns = ttk.Frame(box); btns.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(6,0))
        for i in range(5): btns.grid_columnconfigure(i, weight=1)

        def _pull_from_exp():
            ctx = get_current_context() or {}
            exp_path = ctx.get("exp_json_path")
            if not exp_path:
                messagebox.showwarning("No experiment", "No .exp.json linked."); return
            defs = _read_ml_defaults_from_exp(exp_path)
            if defs is None:
                messagebox.showinfo("No defaults", "This experiment has no ml_defaults yet.")
                return
            editor.delete("1.0","end"); editor.insert("1.0", json.dumps(defs, indent=2, ensure_ascii=False))
            _refresh()

        def _apply_to_exp():
            ctx = get_current_context() or {}
            exp_path = ctx.get("exp_json_path")
            if not exp_path:
                messagebox.showwarning("No experiment", "No .exp.json linked."); return
            d = _get_edit()
            if d is None: return
            if _write_ml_defaults_to_exp(exp_path, d):
                messagebox.showinfo("Saved", f"Updated ml_defaults in:\n{exp_path}")
            _refresh()

        def _load_preset():
            path = filedialog.askopenfilename(title="Load params JSON", filetypes=[("JSON","*.json"), ("All files","*.*")])
            if not path: return
            try:
                with open(path, "r", encoding="utf-8") as f: d = json.load(f)
            except Exception as e:
                messagebox.showerror("Load failed", str(e)); return
            # accept whole-file payloads (exp files) or plain params
            d2 = d.get("ml_defaults", d)
            editor.delete("1.0","end"); editor.insert("1.0", json.dumps(d2, indent=2, ensure_ascii=False))
            _refresh()

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

        def _validate_use():
            ctx = get_current_context() or {}
            exp_path = ctx.get("exp_json_path")
            exp_defs = _read_ml_defaults_from_exp(exp_path) if exp_path else None
            d = _get_edit()
            if d is None: return
            eff = _merge_ml(BUILTIN_ML, exp_defs, d)
            on_effective_params_ready(eff)
            _set_eff(eff)
            messagebox.showinfo("Ready", "Effective ML parameters are validated and ready to use for training.")

        ttk.Button(btns, text="Load preset…", command=_load_preset).grid(row=0, column=0, sticky="ew", padx=2)
        ttk.Button(btns, text="Save preset as…", command=_save_preset_as).grid(row=0, column=1, sticky="ew", padx=2)
        ttk.Button(btns, text="Pull from experiment", command=_pull_from_exp).grid(row=0, column=2, sticky="ew", padx=2)
        ttk.Button(btns, text="Apply to experiment", command=_apply_to_exp).grid(row=0, column=3, sticky="ew", padx=2)
        ttk.Button(btns, text="Validate & Use", command=_validate_use).grid(row=0, column=4, sticky="ew", padx=2)

        # live preview on edit
        editor.bind("<<Modified>>", lambda e: (editor.edit_modified(False), _refresh()))
        return box
    #

    # ---------- END NEW helper ----------

    # state
    model_file_path = None
    predict_input_path = None
    train_csv_path = None
    linked_exp_json = None

    # ---------- NEW: GUI vars for balancing/stratify ----------
    n_estimators_var = tk.IntVar(value=400)     # default trees
    class_weight_var = tk.BooleanVar(value=True)  # use "balanced"  
    test_split_var   = tk.DoubleVar(value=0.20)   # test %
    val_split_var    = tk.DoubleVar(value=0.10)   # val %
    min_samples_var  = tk.IntVar(value=5)         # tiny class threshold
    use_balance_var  = tk.BooleanVar(value=True)  # enable balancing pipeline
    majority_label_var = tk.StringVar(value="Non-glycan")
    majority_factor_var = tk.IntVar(value=3)      # cap = factor * max(minor)
    use_stratify_var = tk.BooleanVar(value=True)
    # Cap only the training fold (leave Val/Test uncapped)
    real_world_test_var = tk.BooleanVar(value=False)
    # Confidence thresholding (post-prediction) 20250902, final addition
    enable_thresh_var = tk.BooleanVar(value=True)   # default ON
    thresh_val_var    = tk.DoubleVar(value=0.65)    # τ in [0,1]
    margin_val_var = tk.DoubleVar(value=0.05)  # δ for majority-support check
    #read current values
    def _on_change(*_):
        _write_train_ui_into_state(
            test_split_var.get(),
            val_split_var.get(),
            min_samples_var.get(),
            use_balance_var.get(),
            majority_label_var.get(),
            majority_factor_var.get(),
            use_stratify_var.get(),
            n_estimators_var.get(),
            class_weight_var.get(),
            real_world_test_var.get(),
            enable_thresh_var.get(),
            thresh_val_var.get(),
            margin_val_var.get(),
        )
    for v in (test_split_var, val_split_var, min_samples_var, use_balance_var,
          majority_label_var, majority_factor_var, use_stratify_var,
          n_estimators_var, class_weight_var, real_world_test_var,
          enable_thresh_var, thresh_val_var, margin_val_var):
        v.trace_add("write", _on_change)

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

    #move back to avoid ref before assignment
    # ===== Modal: ML Parameters Editor (Step 4) =====
    def _update_ml_summary():
        try:
            m = effective_ml_params.get("model", {})
            ml_summary_var.set(
                f'Params: RF n_estimators={m.get("n_estimators", 400)}, '
                f'max_depth={m.get("max_depth", None)}, '
                f'min_split={m.get("min_samples_split", 2)}, '
                f'min_leaf={m.get("min_samples_leaf", 1)}, '
                f'class_weight={m.get("class_weight", "balanced")}'
            )
        except Exception:
            ml_summary_var.set("Params: (using built-ins)")

    # ---- ML state (analysis window scope) ----
    BUILTIN_ML = {
        "model": {"type": "RandomForest", "n_estimators": 400, "max_depth": None,
                "min_samples_split": 2, "min_samples_leaf": 1, "random_state": 42,
                "class_weight": "balanced"},
        "train": {"test_size": 0.25, "val_size": 0.10, "stratify": True,
                "real_world_test": False,
                "threshold": {"enabled": False, "tau": 0.5, "margin": 0.02}}
    }

    effective_ml_params = {}  # what Train uses
    ml_summary_var = tk.StringVar(value="Params: (using built-ins)")

    def _widget_alive(w):
        try:
            return (w is not None) and int(w.winfo_exists()) == 1
        except Exception:
            return False

    def _get_editor_preview_widgets():
        ed = _ml_state.get("editor_txt")
        pv = _ml_state.get("preview_txt")
        if not _widget_alive(ed):
            ed = None
        if not _widget_alive(pv):
            pv = None
        return ed, pv

    def _update_ml_summary():
        m = (effective_ml_params or {}).get("model", {})
        ml_summary_var.set(
            f'Params: RF n_estimators={m.get("n_estimators", 400)}, '
            f'max_depth={m.get("max_depth", None)}, '
            f'min_split={m.get("min_samples_split", 2)}, '
            f'min_leaf={m.get("min_samples_leaf", 1)}, '
            f'class_weight={m.get("class_weight", "balanced")}'
        )
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

    def _merge_ml(*layers):
        out = {}
        for layer in layers:
            if not layer: continue
            for k, v in layer.items():
                if isinstance(v, dict) and isinstance(out.get(k), dict):
                    out[k] = {**out[k], **v}
                else:
                    out[k] = v
        return out  
    # inside open_ml_analysis_window(), before you define open_ml_params_window()
    def _get_ml_context():
        # If you later want to thread in better defaults from Prepare Dataset,
        # fill these two fields there and keep this signature.
        return {"method_json": None, "exp_json": None}

    # --- Shared state for ML (lives inside ML window scope) ---
    effective_ml_params = {}                     # final dict trainer will consume
    ml_summary_var = tk.StringVar(value="Params: (using built-ins)")

    # train/test + RF tk variables (shared between settings dialog and the rest)
    test_split_var       = tk.DoubleVar(value=0.25)
    val_split_var        = tk.DoubleVar(value=0.10)
    min_samples_var      = tk.IntVar(value=5)
    use_balance_var      = tk.BooleanVar(value=True)
    majority_label_var   = tk.StringVar(value="Non-glycan")
    majority_factor_var  = tk.IntVar(value=3)
    use_stratify_var     = tk.BooleanVar(value=True)
    real_world_test_var  = tk.BooleanVar(value=False)
    class_weight_var     = tk.BooleanVar(value=True)  # => "balanced" if True else None
    n_estimators_var     = tk.IntVar(value=400)

    # thresholding
    enable_thresh_var    = tk.BooleanVar(value=False)
    thresh_val_var       = tk.DoubleVar(value=0.65)
    margin_val_var       = tk.DoubleVar(value=0.05)

    # default model + train (used when nothing else provided)
    BUILTIN_ML = {
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
            "test_size": 0.25,
            "val_size": 0.10,
            "stratify": True,
            "real_world_test": False,
            "threshold": {"enabled": False, "tau": 0.50, "margin": 0.02},
        },
    }

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

    def _update_ml_summary():
        m = (effective_ml_params or BUILTIN_ML).get("model", {})
        ml_summary_var.set(
            f'Params: RF n_estimators={m.get("n_estimators")}, '
            f'class_weight={m.get("class_weight")}, '
            f'test={ (effective_ml_params or BUILTIN_ML).get("train",{}).get("test_size") }, '
            f'val={ (effective_ml_params or BUILTIN_ML).get("train",{}).get("val_size") }'
        )


    def open_ml_params_window():
        win = tk.Toplevel(root)
        win.title("ML Parameters")
        win.geometry("820x520")
        win.transient(root)
        win.grab_set()

        # --- left editor ---
        left = ttk.LabelFrame(win, text="Editor (JSON)")
        left.pack(side="left", fill="both", expand=True, padx=(10,5), pady=10)
        editor_txt = tk.Text(left, wrap="none")
        editor_txt.pack(fill="both", expand=True, padx=8, pady=8)

        seed = _ensure_ml_state().get("editor") or BUILTIN_ML
        editor_txt.insert("1.0", json.dumps(seed, indent=2))

        # seed the editor with current effective OR builtins
        #seed = effective_ml_params if effective_ml_params else BUILTIN_ML
        #editor_txt.insert("1.0", json.dumps(seed, indent=2))

        # --- right preview (effective after merge) ---
        right = ttk.LabelFrame(win, text="Effective (Builtins ← Files ← Editor)")
        right.pack(side="left", fill="both", expand=True, padx=(5,10), pady=10)
        preview_txt = tk.Text(right, wrap="none", state="disabled")
        preview_txt.pack(fill="both", expand=True, padx=8, pady=8)
        debug_var = tk.StringVar(value="")
        ttk.Label(right, textvariable=debug_var).pack(anchor="w", padx=8, pady=(0,6))

        #_ml_state["widgets"]["editor_txt"] = editor_txt   # left textbox
        #_ml_state["widgets"]["preview_txt"] = preview_txt # right textbox
        #_refresh_json_textboxes()

        # after building the two Text widgets:
        #_ml_state["editor_txt"] = editor_txt    # <--- your left text widget
        #_ml_state["preview_txt"] = preview_txt  # <--- your right text widget

        st = _ensure_ml_state()
        st["editor_txt"] = editor_txt
        st["preview_txt"] = preview_txt
        _emit_ml_state_changed_ui_refresh()

        def _on_params_close():
            # Drop dead widget refs so future updates become no-ops
            _ml_state["editor_txt"]  = None
            _ml_state["preview_txt"] = None
            win.destroy()

        win.protocol("WM_DELETE_WINDOW", _on_params_close)

        #_write_editor_from_state()
        # seed these panes from the current state when the window opens

        # --- local helpers that use the widgets above ---
        def _get_editor_json_or_empty():
            try:
                return json.loads(editor_txt.get("1.0", "end").strip() or "{}")
            except Exception:
                return {}

        def _read_ml_from_json(path):
            if not path or not os.path.exists(path): 
                return None
            try:
                d = json.load(open(path, "r", encoding="utf-8")) or {}
                ml = d.get("ml") or {}
                return ml.get("parameters") or ml
            except Exception:
                return None

        def _get_ml_context():
            # If you already wrote a context getter elsewhere, you can call it here.
            return {"exp_json": None, "method_json": None}

        def _refresh_effective():
            # ← the key: include snapshot_train_vars() in the merge
            ctx = _get_ml_context()
            exp_defs  = _read_ml_from_json(ctx.get("exp_json"))
            meth_defs = _read_ml_from_json(ctx.get("method_json"))
            editor    = _get_editor_json_or_empty()
            eff = _merge_ml(BUILTIN_ML,exp_defs,meth_defs,editor,snapshot_train_vars())               # editor first…snapshot_train_vars() # …then UI snapshot overrides
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

        def _validate_and_use():
            nonlocal effective_ml_params
            ctx = _get_ml_context()
            exp_defs  = _read_ml_from_json(ctx.get("exp_json"))
            meth_defs = _read_ml_from_json(ctx.get("method_json"))
            editor    = _get_editor_json_or_empty()
            eff = _merge_ml(BUILTIN_ML, exp_defs, meth_defs, snapshot_train_vars(), editor)
            effective_ml_params = eff
            _ensure_ml_state()["editor"] = editor        # persist what’s in the left pane
            _ensure_ml_state()["effective_params"] = eff # optional: keep a copy
            _update_ml_summary()
            win.destroy()
        """
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
            effective_ml_params = _merge_ml(BUILTIN_ML, exp_defs, meth_defs, ed)
            _update_ml_summary()
            _refresh_preview()
            messagebox.showinfo("Ready", "Effective ML parameters set for training.")
        """
        def _save_effective_to_method():
            # pick current method or let user create one
            path = _get_ml_context().get("method_json") or filedialog.asksaveasfilename(
                title="Save or choose method JSON", defaultextension=".json",
                filetypes=[("JSON","*.json")]
            )
            if not path: return

            # ensure effective is up-to-date with what’s in the editor
            try:
                ed = json.loads(editor_txt.get("1.0", "end"))
            except Exception:
                ed = {}
            ctx = _get_ml_context()
            exp_defs  = _read_ml_from_json(ctx.get("exp_json"))
            meth_defs = _read_ml_from_json(ctx.get("method_json"))
            editor = _get_editor_json_or_empty()
            eff = _merge_ml(BUILTIN_ML, exp_defs, meth_defs, editor, snapshot_train_vars())
            os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

            # write back to method JSON
            data = {}
            if os.path.exists(path):
                try:
                    with open(path, "r", encoding="utf-8") as f:
                        txt = f.read().strip()
                        data = json.loads(txt) if txt else {}
                except Exception:
                    data = {}
            data.setdefault("ml", {})
            data["ml"]["updated"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
            #data["ml"]["updated"] = datetime.now().strftime("%Y-%m-%d %H:%M")
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
        # bind Close button to _on_close
        # live preview on typing
        def _on_key(_evt=None): _refresh_effective()
        editor_txt.bind("<KeyRelease>", _on_key)
        # --- footer buttons ---
        btns = ttk.Frame(win)
        btns.pack(fill="x", padx=10, pady=(0,10))

        
        ttk.Button(btns, text="Load from method…",
                command=_load_from_method).pack(side="left", padx=4)

        ttk.Button(btns, text="Validate & Use",
                command=_validate_and_use).pack(side="left", padx=4)

        ttk.Button(btns, text="Save to method…",
                command=_save_effective_to_method).pack(side="left", padx=12)

        ttk.Button(btns, text="Close",
                command=_on_params_close).pack(side="right", padx=4)
        # live preview while typing in the left editor
        editor_txt.bind("<KeyRelease>", lambda _=None: _refresh_effective())

        # seed the right pane immediately on window open
        _refresh_effective()
        """
        # footer buttons (all local to the modal)
        btns = ttk.Frame(win); btns.pack(fill="x", padx=10, pady=(0,10))
        ttk.Button(btns, text="Load from method…", command=_load_from_method).pack(side="left", padx=4)
        ttk.Button(btns, text="Validate & Use",    command=_validate_and_use).pack(side="left", padx=4)
        # live preview
        editor_txt.bind("<KeyRelease>", lambda _=None: _refresh_effective())
        _refresh_effective()
        ttk.Button(btns, text="Save to method…",   command=_save_effective_to_method).pack(side="left", padx=12)
        ttk.Button(btns, text="Close",             command=win.destroy).pack(side="right", padx=4)

        _refresh_effective()
        """
        # 20250919 new version of ML exp json
        # ----- defaults you already use -----
        """
        BUILTIN_ML = {
            "model": {"type": "RandomForest", "n_estimators": 400, "max_depth": None,
                    "min_samples_split": 2, "min_samples_leaf": 1, "random_state": 42,
                    "class_weight": "balanced"},
            "train": {"test_size": 0.25, "val_size": 0.10, "stratify": True,
                    "real_world_test": False,
                    "threshold": {"enabled": False, "tau": 0.5, "margin": 0.02}}
        }

        effective_ml_params = {}
        ml_summary_var = tk.StringVar(value="Params: (using built-ins)")

        def _update_ml_summary():
            m = (effective_ml_params or {}).get("model", {})
            ml_summary_var.set(
                f'Params: RF n_estimators={m.get("n_estimators", 400)}, '
                f'max_depth={m.get("max_depth", None)}, '
                f'min_split={m.get("min_samples_split", 2)}, '
                f'min_leaf={m.get("min_samples_leaf", 1)}, '
                f'class_weight={m.get("class_weight", "balanced")}'
            )

        def _read_ml_from_json(path):
            if not path or not os.path.exists(path):
                return None
            try:
                with open(path, "r", encoding="utf-8") as f:
                    d = json.load(f) or {}
                ml = d.get("ml") or {}
                return ml.get("parameters") or ml  # accept legacy
            except Exception:
                return None

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
        """
        def _ensure_effective_from_editor_if_empty():
            """If user forgot Validate, synthesize effective from current editor + jsons."""
            nonlocal effective_ml_params
            if effective_ml_params:
                return
            ctx = _get_ml_context()  # you already have this getter
            exp_defs   = _read_ml_from_json(ctx.get("exp_json"))
            meth_defs  = _read_ml_from_json(ctx.get("method_json"))
            editor_now = _get_editor_json_or_empty()  # defined inside the ML modal; fallback to {}
            effective_ml_params = _merge_ml(BUILTIN_ML, exp_defs, meth_defs, editor_now)
            _update_ml_summary()    

        # --- editor widgets you already have ---
        # editor_txt = tk.Text(...); preview_txt = tk.Text(...)

        def _get_editor_json_or_empty():
            try:
                return json.loads(editor_txt.get("1.0", "end").strip() or "{}")
            except Exception:
                return {}

        def _refresh_effective_preview():
            _refresh_effective()  # keep one code path
            """
            ctx = _get_ml_context()
            exp_defs  = _read_ml_from_json(ctx.get("exp_json"))
            meth_defs = _read_ml_from_json(ctx.get("method_json"))
            ed        = _get_editor_json_or_empty()
            eff = _merge_ml(BUILTIN_ML, exp_defs, meth_defs, ed)
            preview_txt.config(state="normal")
            preview_txt.delete("1.0", "end")
            preview_txt.insert("1.0", json.dumps(eff, indent=2))
            preview_txt.config(state="disabled")
            """

        def _load_from_method():
            # prefer currently-linked method; else let user pick
            p = _get_ml_context().get("method_json") or filedialog.askopenfilename(
                title="Choose method JSON", filetypes=[("JSON","*.json")])
            if not p:
                return
            defs = _read_ml_from_json(p)
            if not defs:
                messagebox.showwarning("ML Parameters", "No ML block found in that JSON.")
                return
            editor_txt.delete("1.0", "end")
            editor_txt.insert("1.0", json.dumps(defs, indent=2))
            _refresh_effective_preview()

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
            effective_ml_params = _merge_ml(BUILTIN_ML, exp_defs, meth_defs, ed)
            _update_ml_summary()
            _refresh_effective_preview()
            messagebox.showinfo("Ready", "Effective ML parameters set for training.")

        def _save_effective_to_method():
            ctx = _get_ml_context()
            path = ctx.get("method_json")
            if not path:
                path = filedialog.asksaveasfilename(
                    title="Save or choose method JSON",
                    defaultextension=".json",
                    filetypes=[("JSON","*.json")])
            if not path:
                return
            # ensure we have something to save
            _ensure_effective_from_editor_if_empty()
            data = {}
            if os.path.exists(path):
                try:
                    with open(path, "r", encoding="utf-8") as f:
                        txt = f.read().strip()
                        data = json.loads(txt) if txt else {}
                except Exception:
                    data = {}
            data.setdefault("ml", {})
            data["ml"]["updated"] = datetime.now().strftime("%Y-%m-%d %H:%M")
            data["ml"]["parameters"] = effective_ml_params
            with open(path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            messagebox.showinfo("Saved", f"ML parameters written to:\n{os.path.basename(path)}")

        # ===== Step 4 (existing controls) =====
    # ... your Step 4 labels/entries ...

    # summary (read-only) shows what will be used if set
    #ml_summary_var = tk.StringVar(value="Params: (using built-ins)")
    #def _update_ml_summary():
    #    try:
    #        m = effective_ml_params.get("model", {})
    ###        ml_summary_var.set(
    #            f'Params: RF n_estimators={m.get("n_estimators", 400)}, '
    #            f'max_depth={m.get("max_depth", None)}, '
    #            f'min_split={m.get("min_samples_split", 2)}, '
    #            f'min_leaf={m.get("min_samples_leaf", 1)}, '
    #            f'class_weight={m.get("class_weight", "balanced")}'
    #        )
    #    except Exception:
    #        ml_summary_var.set("Params: (using built-ins)")    


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

        def _commit_and_close():
            # Generate a fresh effective preview by merging BUILTIN + live train vars
            # (Editor JSON / files will still be merged inside the ML-params window)
            nonlocal effective_ml_params
            eff_local = _merge_ml(BUILTIN_ML, snapshot_train_vars())
            # Don’t overwrite global effective here; we only show a summary hint.
            # The real commit is done by ML Parameters → “Validate & Use”.
            _update_ml_summary()
            settings.destroy()

        def _pull_state_into_vars():
            s = _ml_state["current"]
            t = s.get("train", {})
            b = s.get("balance", {})
            m = s.get("model", {})

            # sliders/spinboxes/checks/entries you already have:
            test_split_var.set(float(t.get("test_size", 0.20)))
            val_split_var.set(float(t.get("val_size", 0.0)))
            use_stratify_var.set(1 if t.get("stratify", True) else 0)
            real_world_test_var.set(1 if t.get("real_world_test", False) else 0)

            thr = t.get("threshold", {})
            enable_thresh_var.set(1 if thr.get("enabled", False) else 0)
            thresh_val_var.set(float(thr.get("tau", 0.65)))
            margin_val_var.set(float(thr.get("margin", 0.05)))

            min_samples_var.set(int(s.get("min_samples_per_class", 5)))

            use_balance_var.set(1 if b.get("enabled", True) else 0)
            majority_label_var.set(b.get("majority_label", "Non-glycan"))
            majority_factor_var.set(int(b.get("majority_factor", 3)))

            n_estimators_var.set(int(m.get("n_estimators", 400)))
            class_weight_var.set(1 if (m.get("class_weight") == "balanced") else 0)

        def _push_vars_into_state(*_):
            _write_train_ui_into_state(
                test_size       = float(test_split_var.get()),
                val_size        = float(val_split_var.get()),
                min_per_class   = int(min_samples_var.get()),
                use_balance     = bool(use_balance_var.get()),
                majority_label  = str(majority_label_var.get()),
                majority_factor = int(majority_factor_var.get()),
                stratify        = bool(use_stratify_var.get()),
                n_estimators    = int(n_estimators_var.get()),
                use_class_weight= bool(class_weight_var.get()),
                real_world_test = bool(real_world_test_var.get()),
                thr_enable      = bool(enable_thresh_var.get()),
                thr_tau         = float(thresh_val_var.get()),
                thr_margin      = float(margin_val_var.get()),
            )
            _emit_ml_state_changed_ui_refresh()  # safe no-op if editor window isn’t open

        # 1) seed the controls when the dialog opens
        _pull_state_into_vars()

        # 2) wire all variables for live updates
        for var in (
            test_split_var, val_split_var, use_stratify_var, real_world_test_var,
            enable_thresh_var, thresh_val_var, margin_val_var,
            min_samples_var, use_balance_var, majority_label_var, majority_factor_var,
            n_estimators_var, class_weight_var,
        ):
            var.trace_add("write", _push_vars_into_state)

        # 3) also bind Scales so dragging updates continuously
        def _scale_cb(_val):
            _push_vars_into_state()

        for sc in (test_split_var, val_split_var):  # use your actual Scale objects
            try:
                sc.configure(command=_scale_cb)
            except Exception:
                pass

        def _update_split_labels(*_):
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

        # tiny class
        tk.Label(settings, text="Minimum samples per class").pack(pady=(12,0))
        tk.Spinbox(settings, from_=1, to=50, textvariable=min_samples_var, width=6).pack()

        # balancing block
        tk.Checkbutton(settings, text="Enable class balancing (cap majority)", variable=use_balance_var).pack(pady=(12,0))
        row = tk.Frame(settings); row.pack(pady=2)
        tk.Label(row, text="Majority label:").pack(side="left")
        tk.Entry(row, textvariable=majority_label_var, width=16).pack(side="left", padx=6)
        row2 = tk.Frame(settings); row2.pack(pady=2)
        tk.Label(row2, text="Majority factor (×max minor):").pack(side="left")
        tk.Spinbox(row2, from_=1, to=20, textvariable=majority_factor_var, width=6).pack(side="left", padx=6)

        # stratify
        tk.Checkbutton(settings, text="Stratify by label", variable=use_stratify_var).pack(pady=(12,0))

        # --- Random Forest options ---
        sep = ttk.Separator(settings, orient="horizontal"); sep.pack(fill="x", padx=10, pady=(12,6))
        tk.Label(settings, text="Random Forest Options").pack()

        row_rf = tk.Frame(settings); row_rf.pack(pady=2)
        tk.Label(row_rf, text="n_estimators (trees):").pack(side="left")
        tk.Spinbox(row_rf, from_=50, to=2000, increment=50, textvariable=n_estimators_var, width=7).pack(side="left", padx=6)

        tk.Checkbutton(settings, text='Use class_weight = "balanced"', variable=class_weight_var).pack(pady=2)
        tk.Checkbutton(settings, text="Real-world test (cap training only)",
               variable=real_world_test_var).pack(pady=(4,0))
        

        sep3 = ttk.Separator(settings, orient="horizontal"); sep3.pack(fill="x", padx=10, pady=(12,6))
        tk.Label(settings, text="Prediction Thresholding").pack()

        tk.Checkbutton(
            settings,
            text="Enable confidence threshold → fallback to Majority label",
            variable=enable_thresh_var
        ).pack(pady=(2,2))

        row_thr = tk.Frame(settings); row_thr.pack(pady=2)
        tk.Label(row_thr, text="Threshold τ (0.00–0.99):").pack(side="left")
        tk.Spinbox(row_thr, from_=0.00, to=0.99, increment=0.01,
                textvariable=thresh_val_var, width=6).pack(side="left", padx=6)

        row_margin = tk.Frame(settings); row_margin.pack(pady=2)
        tk.Label(row_margin, text="Majority support margin δ (0.00–0.20):").pack(side="left")
        tk.Spinbox(row_margin, from_=0.00, to=0.20, increment=0.01,
                textvariable=margin_val_var, width=6).pack(side="left", padx=6)

        # ... all your Step 4 widgets ...
        tk.Button(settings, text="Edit ML Parameters…",
                command=open_ml_params_window).pack(pady=(8, 0))
        
        tk.Button(settings, text="Apply to ML Param[if not reflect]", command=_push_vars_into_state).pack(pady=(8,0))

        # (a) define tk variables: test_split_var, val_split_var, ... margin_val_var
        # (b) define _on_change(...) that writes the current widget values into _ml_state["editor"]
        # (c) build all the widgets that bind to those variables
        # ---- PASTE THE SEED BLOCK RIGHT HERE ----
        # (d) Close button, settings.mainloop/deiconify/return
        # ---- seed UI from current editor state ----
        ed = _ml_state.get("editor", {})
        tv = ed.get("train", {}) or {}
        bv = ed.get("balance", {}) or {}
        mv = ed.get("model", {}) or {}

        test_split_var.set(tv.get("test_size", 0.20))
        val_split_var.set(tv.get("val_size", 0.0))
        min_samples_var.set(tv.get("min_samples_per_class", 5))

        use_balance_var.set(bv.get("enabled", True))
        majority_label_var.set(bv.get("majority_label", "Non-glycan"))
        majority_factor_var.set(bv.get("majority_factor", 3))

        use_stratify_var.set(tv.get("stratify", True))

        n_estimators_var.set(mv.get("n_estimators", 400))
        class_weight_var.set((mv.get("class_weight") or "balanced") == "balanced")
        #real world change to True as default #20251006
        real_world_test_var.set(tv.get("real_world_test", True))#False))

        th = tv.get("threshold", {}) or {}
        enable_thresh_var.set(th.get("enabled", False))
        thresh_val_var.set(th.get("tau", 0.65))
        margin_val_var.set(th.get("margin", 0.05))

        # push once so editors / preview reflect these values
        _on_change()


        tk.Button(settings, text="Close", command=settings.destroy).pack(pady=8)

    #extra train settings save/load
    #row4 = train_tab.grid_size()[1]
    #step4_container = ttk.Frame(train_tab)
    ##step4_container.grid(row=row4, column=0, columnspan=3,
    #                    sticky="ew", padx=10, pady=(6, 0))

    #step4_bar = tk.Frame(step4_container); step4_bar.pack(fill="x")
    #tk.Button(step4_bar, text="Train/Test Parameters…",
    #        command=open_train_settings).pack(side="left")
    #tk.Button(step4_bar, text="Edit ML Parameters…",
    #        command=open_ml_params_window).pack(side="left", padx=6)
    #tk.Label(step4_bar, textvariable=ml_summary_var).pack(side="left", padx=12)
    # -----------------------------------------------

    def train_model():
        try:
            from sklearn.ensemble import RandomForestClassifier
            from sklearn.metrics import classification_report
            from sklearn.preprocessing import LabelEncoder
            from sklearn.utils.multiclass import unique_labels
            import joblib
            import pretrain_normalizer as normalizer
        except ImportError:
            messagebox.showerror("Missing Dependencies",
                                 "Please install scikit-learn and joblib.")
            return
        #20251005
        use_mass_feature = bool(include_mass_train_var.get())

        if not effective_ml_params:
            ctx = _get_ml_context()
            exp_defs  = _read_ml_from_json(ctx.get("exp_json"))
            meth_defs = _read_ml_from_json(ctx.get("method_json"))
            effective = _merge_ml(BUILTIN_ML, exp_defs, meth_defs, {})
        else:
            effective = effective_ml_params
        cfg = effective_ml_params or _merge_ml(BUILTIN_ML, snapshot_train_vars())
        rf = cfg.get("model", {})
        tr = cfg.get("train", {})

        #rf = effective.get("model", {})
        #train_cfg = effective.get("train", {})

        if not train_csv_path:
            messagebox.showwarning("No File", "Please select a training CSV first.")
            return

        label_col = label_dropdown.get().strip()
        if label_col not in ['Structure', 'IUPACname(optional)', 'Glycanannotation2', 'GlyToucan ID']:
            messagebox.showerror("Invalid Label", f"'{label_col}' is not a supported label column.")
            return

        # load
        try:
            df = pd.read_csv(train_csv_path)
        except Exception as e:
            messagebox.showerror("Read Error", str(e)); return
        if label_col not in df.columns:
            messagebox.showerror("Missing Column", f"{label_col} not found."); return

        #20250929 ver to accept SO3(s), HexA and PO3H(p)
        # convert the composition to manual-style ONLY when needed
        print("[MLdebug] labels before:", sorted(set(df[label_col].astype(str)))[:10], "...")

        # Normalize ONLY when needed (do not strip A / s / p!)
        try:
            if label_col == "Structure":
                before = sorted(set(df[label_col].astype(str)))[:10]
                df[label_col] = df[label_col].map(
                    lambda x: normalize_structure_for_training(x, normalizer.parse_structure_to_manual)
                )
                after  = sorted(set(df[label_col].astype(str)))[:10]
                print("[MLdebug] labels before:", before, "...")
                print("[MLdebug] labels after:",  after,  "...")
            # For GlyToucan/IUPAC, leave as-is.
        except Exception as e:
            print("[dev] structure normalization warning:", e)
            messagebox.showerror(
                "Structure normalization failed",
                "Can't normalize labels. See console for details."
            )
            return
        print("[MLdebug] labels after:",  sorted(set(df[label_col].astype(str)))[:10], "...")
        """
        #convert the composition to man-readable one
        try:
            df[label_col] = df[label_col].apply(normalizer.parse_structure_to_manual)
        except Exception as e:
            print(f"[dev]Encounter structure conversion error")
            messagebox.showerror("Can't convert the structure(composition) to man-like format", str(e)); return
        """
        # encode label AFTER any string/tuple harmonization (if needed)
        y_raw = df[label_col].astype(str)

        # run balancing + split
        try:
            if use_balance_var.get():
                X_train, y_train, X_val, y_val, X_test, y_test, info = balance_and_split(
                    df.assign(**{label_col: y_raw}),
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
                    df.assign(**{label_col: y_raw}),
                    label_col=label_col,
                    majority_label="__no_cap__",  # won't match → no cap
                    min_count=int(min_samples_var.get()),
                    majority_factor=1,
                    test_size=float(test_split_var.get()),
                    val_size=float(val_split_var.get()),
                    stratify=bool(use_stratify_var.get()),
                    cap_training_only=bool(real_world_test_var.get()),   # << NEW
                    random_state=42
                )
        except Exception as e:
            messagebox.showerror("Split/Balancing Failed", str(e)); return

        #20251005
        # === MS1 feature handling (after split, before label encoding) ===
        # use_mass_feature already read above: use_mass_feature = bool(include_mass_train_var.get())

        def _prep_ms1_block(X, y, use_mass):
            #import pandas as _pd
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
        y_val_enc   = le.transform(y_val)
        y_test_enc  = le.transform(y_test)


        # model (consume panel params if available; otherwise fallback to current UI/defaults)
        rf = effective_ml_params.get("model", {}) if effective_ml_params else {}
        model = RandomForestClassifier(
            n_estimators      = int(rf.get("n_estimators", 400)),
            max_depth         = rf.get("max_depth", None),
            min_samples_split = int(rf.get("min_samples_split", 2)),
            min_samples_leaf  = int(rf.get("min_samples_leaf", 1)),
            random_state      = int(rf.get("random_state", 42)),
            class_weight      = rf.get("class_weight", ("balanced" if bool(class_weight_var.get()) else None)),
            n_jobs            = -1,
        )
        # rf_cfg = {}
        #try:
        #    if effective_ml_params:
        #        if effective_ml_params.get("model", {}).get("type", "RandomForest") == "RandomForest":
        #            rf_cfg = effective_ml_params.get("model", {})
        #except Exception:
        #    rf_cfg = {}

        #model = RandomForestClassifier(
        #    n_estimators = int(rf_cfg.get("n_estimators", 400)),
        #    max_depth    = rf_cfg.get("max_depth", None),
        #    min_samples_split = int(rf_cfg.get("min_samples_split", 2)),
        #    min_samples_leaf  = int(rf_cfg.get("min_samples_leaf", 1)),
        #    random_state = int(rf_cfg.get("random_state", 42)),
        #    class_weight = rf_cfg.get("class_weight", ("balanced" if bool(class_weight_var.get()) else None)),
        #    n_jobs = -1
        #)

        # model
        #model = RandomForestClassifier(
        #    n_estimators=int(n_estimators_var.get()),
        #    class_weight=("balanced" if class_weight_var.get() else None),
        #    max_features="sqrt",
        #    random_state=42,
        #    n_jobs=-1
        #)

        #20250929
        # --- NEW: sanitize features & lock column order ---
        X_train_num = _drop_meta_and_get_X(X_train)
        cols = list(X_train_num.columns)

        # align val/test to train's columns (fill missing with 0)
        X_val_num  = _drop_meta_and_get_X(X_val).reindex(columns=cols, fill_value=0)
        X_test_num = _drop_meta_and_get_X(X_test).reindex(columns=cols, fill_value=0)

        # replace originals
        X_train, X_val, X_test = X_train_num, X_val_num, X_test_num

        # --- END NEW ---
        train_feats = list(X_train.columns)

        mass_suffix = "_withMass" if use_mass_feature else ""

        if X_train.isna().any().any():
            raise RuntimeError("[Train] NaNs remain in X_train after MS1 prep — unexpected")
        
        #model = RandomForestClassifier(n_estimators=100, random_state=42)
        model.fit(X_train, y_train_enc)
        from sklearn.metrics import classification_report, confusion_matrix
        from sklearn.utils.multiclass import unique_labels
        import numpy as np
        # evaluate with proba
        # raw predictions
        y_pred = model.predict(X_test)

        #extra binary glycan-vs-non-glycan report (print in terminal only)
        from sklearn.metrics import precision_recall_fscore_support, accuracy_score

        maj_idx = list(le.classes_).index(majority_label_var.get())
        y_true_bin = (y_test_enc != maj_idx).astype(int)   # 1=glycan, 0=non
        y_pred_bin = (y_pred      != maj_idx).astype(int)

        p,r,f,_ = precision_recall_fscore_support(y_true_bin, y_pred_bin, average="binary", zero_division=0)
        acc = accuracy_score(y_true_bin, y_pred_bin)
        print(f"[ML] Binary glycan-vs-non :: P={p:.2f} R={r:.2f} F1={f:.2f} Acc={acc:.2f}")

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

            """
            if hasattr(model, "predict_proba"):
                proba = model.predict_proba(X_test)  # [n_samples, n_classes]
                maxp = proba.max(axis=1)

                # indices in the fitted LabelEncoder
                classes_ = list(le.classes_)
                try:
                    majority_idx = classes_.index(majority_label_var.get())
                except ValueError:
                    majority_idx = None

                if majority_idx is not None:
                    import numpy as np
                    # apply threshold ONLY when the predicted class is a glycan (not the majority)
                    mask = (y_pred != majority_idx) & (maxp < tau)
                    if np.any(mask):
                        y_pred = y_pred.copy()
                        y_pred[mask] = majority_idx
                else:
                    print(f"[ML] Thresholding skipped: fallback label '{majority_label_var.get()}' not in training classes.")
            else:
                print("[ML] Thresholding skipped: classifier lacks predict_proba.")
                """
        # --- END NEW ---
        # --- NEW: confidence threshold -> fallback to majority label
        """
        if enable_thresh_var.get():
            tau = float(thresh_val_var.get())
            if hasattr(model, "predict_proba"):
                proba = model.predict_proba(X_test)  # shape: [n_samples, n_classes]
                maxp = proba.max(axis=1)

                # find fallback class index in the fitted LabelEncoder
                try:
                    fallback_idx = list(le.classes_).index(majority_label_var.get())
                except ValueError:
                    fallback_idx = None

                if fallback_idx is not None:
                    import numpy as np
                    mask = (maxp < tau)
                    if mask.any():
                        # y_pred is already encoded ints
                        y_pred = y_pred.copy()
                        y_pred[mask] = fallback_idx
                else:
                    print(f"[ML] Thresholding skipped: fallback label '{majority_label_var.get()}' not in training classes.")
            else:
                print("[ML] Thresholding skipped: classifier lacks predict_proba.")
        """
        #report = classification_report(y_test_enc, y_pred, target_names=used_names, zero_division=0)

        # --- END NEW ---
        
        # evaluate on test
        #y_pred = model.predict(X_test)
        used = unique_labels(y_test_enc, y_pred)
        used_names = [le.classes_[i] for i in used]
        # quiet, deterministic handling of 0/0 cases 20250902
        report = classification_report(y_test_enc, y_pred, target_names=used_names, zero_division=0)
        # (optional) flag classes with no predicted or no true samples
        labels_all = list(le.classes_)
        cm = confusion_matrix(y_test_enc, y_pred, labels=np.arange(len(labels_all)))
        no_pred = [labels_all[j] for j, s in enumerate(cm.sum(axis=0)) if s == 0]
        no_true = [labels_all[i] for i, s in enumerate(cm.sum(axis=1)) if s == 0]

        if no_pred or no_true:
            print(f"[ML] No predicted samples for: {no_pred}")
            print(f"[ML] No true samples in test for: {no_true}")
        #report = classification_report(y_test_enc, y_pred, target_names=used_names)

        #glycan only report (maybe wont export as report, need screenshot?)
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
        # build a robust “cap” line that works for both modes
        cap_applied_to = info.get("majority_cap_applied_to") or info.get("train_majority_cap_applied_to")
        cap_value      = info.get("majority_cap") or info.get("train_majority_cap")
        cap_prefix     = "Capped training " if info.get("mode") == "cap_training_only" else "Capped "
        cap_line = (f"{cap_prefix}'{cap_applied_to}' at {cap_value}"
                    if cap_applied_to else "No majority cap applied.")

        msg = [
            "Random Forest trained successfully.",
            f"Classes kept: {len(info['kept_label_counts'])}",
            cap_line,
            f"Dropped tiny classes (< {min_samples_var.get()}): {sum(info['dropped_rare_counts'].values())}",
            "",
            report,
        ]
        """
        msg = [
            "Random Forest trained successfully.",
            f"Classes kept: {len(info['kept_label_counts'])}",
            (f"Capped '{info['majority_cap_applied_to']}' at {info['majority_cap']} "
             if info['majority_cap_applied_to'] else "No majority cap applied."),
            f"Dropped tiny classes (< {min_samples_var.get()}): {sum(info['dropped_rare_counts'].values())}",
            "",
            report
        ]
        """
        #new lines, comment if I feel it annoying
        msg.insert(1, f"RF: {int(n_estimators_var.get())} trees, class_weight="
               f"{'balanced' if class_weight_var.get() else 'none'}")
        msg.insert(1, f"Mode: {'Cap training only' if real_world_test_var.get() else 'Cap before split'}")
        msg.insert(1, f"Thresholding: {'ON τ=' + format(thresh_val_var.get(), '.2f') + ' → ' + majority_label_var.get() if enable_thresh_var.get() else 'OFF'}")
        msg.insert(1, f"Thresholding: {'ON τ=' + format(thresh_val_var.get(), '.2f') + ', δ=' + format(margin_val_var.get(), '.2f') + ' → ' + majority_label_var.get() if enable_thresh_var.get() else 'OFF'}")
        messagebox.showinfo("Training Complete", "\n".join(msg))

        # save artifacts
        base = os.path.splitext(train_csv_path)[0]
        model_path  = base + "_rf_model.joblib"
        enc_path    = base + "_labelencoder.joblib"
        report_path = base + "_rf_performance.txt"
        joblib.dump(model, model_path)
        joblib.dump(le, enc_path)
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(report)

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

    def train_with_optional_ng(
    positives_df,            # trainable (≥0.07 or ≥0.06 or salvage)
    df_pseudo_full,          # long-form pseudolabel TSV
    ion_masses, ppm=20.0,
    include_ng=True,
    low_score_col="ion score",
    low_score_cut=0.03,
    ng_strategy="cap_ng",    # or 'undersample' or 'none'
    max_ng_ratio=1.0,
    class_weight_balanced=True,):
        train_df = positives_df.copy()

        if include_ng:
            ng_df = collect_ng_candidates(
                df_pseudo_full=df_pseudo_full,
                df_pseudo_filtered_pos=positives_df,
                ion_masses=ion_masses,
                ppm=ppm,
                low_score_col=low_score_col,
                low_score_cut=low_score_cut,
                use_low_score=True,
                use_no_hit=True,
                manual_unknown_df=None,  # or your MSlist slice if available
                scan_col="MS2scan_no",
            )
            if ng_strategy in ("cap_ng", "undersample"):
                comb = pd.concat([train_df, ng_df], ignore_index=True)
                train_df = resample_by_strategy(
                    comb, label_col="Structure", strategy=ng_strategy, max_ng_ratio=max_ng_ratio
                )
            else:
                train_df = pd.concat([train_df, ng_df], ignore_index=True)

    def select_model_file():
        nonlocal model_file_path
        path = filedialog.askopenfilename(title="Select Model File",filetypes=[("Model files", "*.joblib *.pkl *.skops"), ("All files", "*.*")])
        #path = filedialog.askopenfilename(filetypes=[("Model files", "*.joblib *.pkl")])
        if path:
            model_file_path = os.path.abspath(path)
            messagebox.showinfo("Model Loaded", f"Model loaded from:\n{model_file_path}")

    def select_predict_input():
        nonlocal predict_input_path
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            predict_input_path = os.path.abspath(path)
            messagebox.showinfo("Input File Selected", f"Data loaded from:\n{predict_input_path}")

    #def create_unlabeled_dataset():
    #    messagebox.showinfo("Not Yet Implemented", "This feature will allow you to select an experiment and automatically create a feature-matched dataset from its annotation and early raw-converted CSV.")
    
    def _read_any_table(p):
        # try TSV first, then CSV with common settings
        try:
            return pd.read_csv(p, sep="\t", engine="python")
        except Exception:
            return pd.read_csv(p, engine="python")

    import re
    _comp_pat = re.compile(r'(KDN|F|H|N|S|G|A|s|p)\s*([0-9]+)', re.I)
    def _canon_comp(s: str) -> str:
        counts = {"F":0,"H":0,"N":0,"S":0,"G":0,"KDN":0,"A":0,"s":0,"p":0}
        for k,v in _comp_pat.findall(s or ""):
            k = "KDN" if k.upper()=="KDN" else k
            counts[k] = counts.get(k,0) + int(v)
        parts=[]
        for key in ("F","H","N","S","G","KDN","A","s","p"):
            n=counts.get(key,0)
            if n>0: parts.append(f"{key}{n}")
        return "".join(parts)
    #20251001 safeguard
    out_gate = None

    def run_prediction():
        import numpy as np
        try:
            import joblib
        except ImportError:
            messagebox.showerror("Missing Dependency", "joblib is not installed.")
            return

        if not model_file_path or not predict_input_path:
            messagebox.showwarning("Missing Info", "Please select both a model file and an input CSV file.")
            return

        #20250909 add import for prediction report
        from prediction_report import (ReportParams, summarize_predictions,write_prediction_report, show_prediction_summary_popup)

        try:
            model, loader = load_model_any(model_file_path)
            #model = joblib.load(model_file_path)
            df = pd.read_csv(predict_input_path)

            # Try loading label encoder if available
            le_path = model_file_path.replace("_rf_model.joblib", "_labelencoder.joblib")
            if os.path.exists(le_path):
                le = joblib.load(le_path)
            else:
                le = None
        except Exception as e:
                # Suggestion to user if it's a version mismatch crash
            if ".joblib" in model_file_path and "InconsistentVersionWarning" in str(e) or "dtype" in str(e):
                messagebox.showerror(
                    "Model version mismatch",
                    "This joblib model was trained with a different scikit-learn version.\n\n"
                    "Fix: re-export as .skops on the training machine (or retrain), then load the .skops here."
                )
                return
            else:
                messagebox.showerror("Load Error", str(e))
            return
        #in case error occurs, add the debug lines
        """
        import sys, sklearn, os
        print("[env] python:", sys.executable)
        print("[env] sklearn:", sklearn.__version__)
        print("[model] file:", model_path)
        print("[model] type:", type(model).__name__)

        # Before aligning:
        print("[features] input cols:", len(df.columns))
        X = df.drop(columns=["MS2scan_no", "protonatedmass"], errors="ignore")
        print("[features] after drop:", len(X.columns))

        # After aligning:
        missing = sorted(set(train_feats) - set(X.columns))
        extra   = sorted(set(X.columns) - set(train_feats))
        for c in missing: X[c] = 0.0
        X = X[list(train_feats)]
        print(f"[features] +{len(missing)} filled, -{len(extra)} dropped, final={X.shape[1]}")
        """

        #202509ver
        """
        try:
            # Build X, drop traceability cols
            X = df.drop(columns=["MS2scan_no"], errors= "ignore")#, "protonatedmass"], errors="ignore")

            # Align features to training set
            train_feats = get_training_features(model, model_file_path)
            #20251006
            # detect if the loaded model expects MS1 features
            ms1_used = False
            if train_feats:
                tf = set(train_feats)
                ms1_used = ("protonatedmass" in tf) or ("delta_ppm" in tf)

            # optional: warn if the model expects mass but the input lacks it
            if ms1_used and "protonatedmass" not in df.columns:
                messagebox.showwarning(
                    "Model expects MS1",
                    "This model was trained with 'protonatedmass', but the input table has no such column. "
                    "Prediction will proceed with a zero-filled placeholder."
                )
                df["protonatedmass"] = 0.0

            if train_feats is None:
                raise RuntimeError("Cannot determine training feature list. "
                                "Re-export the model with embedded feature names or provide *_features.json.")

            missing = sorted(set(train_feats) - set(X.columns))
            extra   = sorted(set(X.columns) - set(train_feats))
            for c in missing: X[c] = 0.0
            X = X[list(train_feats)]  # enforce order
        """


        #20251006 ver
        try:
            # 1) Get the model's training feature list first
            train_feats = get_training_features(model, model_file_path)
            if train_feats is None:
                raise RuntimeError("Cannot determine training feature list. "
                                "Re-export with embedded features or provide *_features.json.")

            # 2) If the model expects MS1 but input lacks it, add placeholder NOW (before X)
            if "protonatedmass" in train_feats and "protonatedmass" not in df.columns:
                messagebox.showwarning(
                    "Model expects MS1",
                    "This model was trained with 'protonatedmass', but the input has no such column. "
                    "Proceeding with a zero-filled placeholder."
                )
                df["protonatedmass"] = 0.0  # float OK

            # (optional) if it expects delta_ppm too
            if "delta_ppm" in train_feats and "delta_ppm" not in df.columns:
                df["delta_ppm"] = 0.0

            # 3) Now build X (after placeholders exist)
            X = df.drop(columns=["MS2scan_no"], errors="ignore")

            # 4) Align to training features
            missing = sorted(set(train_feats) - set(X.columns))
            extra   = sorted(set(X.columns) - set(train_feats))
            for c in missing:
                X[c] = 0.0
            X = X[list(train_feats)]  # enforce order

            # 5) (optional) detect for UI
            #ms1_used = ("protonatedmass" in train_feats) or ("delta_ppm" in train_feats)

            # Predict
            y_pred = model.predict(X)

            # Optional probabilities & margin
            """
            if hasattr(model, "predict_proba"):
                import numpy as np
                proba = model.predict_proba(X)
                top1 = proba.max(axis=1)
                top2 = np.partition(proba, -2, axis=1)[:, -2] if proba.shape[1] > 1 else np.zeros(len(top1))
                df["proba_top1"], df["proba_top2"] = top1, top2
                df["margin"] = top1 - top2
            """
            #X = df.drop(columns=["MS2scan_no"], errors="ignore")
            #y_pred = model.predict(X)

            #applying same filter to prediction model
            # --- optional: apply the same margin-aware demotion in prediction ---
            try:
                tau = float(thresh_val_var.get())
            except Exception:
                tau = 0.60  # sensible fallback
            #similar to tau, use same pattern on margin
            try:
                margin = float(margin_val_var.get())
            except Exception:
                margin = 0.05

            proba = None
            if hasattr(model, "predict_proba"):
                proba = model.predict_proba(X)
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

            # If encoder available, decode
            if le is not None:
                try:
                    y_pred = le.inverse_transform(y_pred)
                except:
                    pass
            elif hasattr(model, 'classes_'):
                y_pred = [model.classes_[i] if isinstance(i, int) else i for i in y_pred]

            df['Predicted_Label'] = y_pred

            #safe guard
            # --- ensure gate columns always exist ---
            import numpy as _np
            if "ppm_precursor" not in df.columns:
                df["ppm_precursor"] = _np.nan
            if "pred_ok" not in df.columns:
                df["pred_ok"] = True     # default: keep all rows

            #add prediction reports
            # --- START: Prediction summary integration ---

            # 1) Collect class names (for readable labels in the report)
            try:
                class_names = list(le.classes_) if le is not None else list(getattr(model, "classes_", []))
            except Exception:
                class_names = list(getattr(model, "classes_", []))

            # 2) If available, get probabilities for margins/entropy
            #already declaired before, with more functionality?
            """
            y_proba = None
            if hasattr(model, "predict_proba"):
                try:
                    y_proba = model.predict_proba(X)  # shape: [n_samples, n_classes]
                except Exception:
                    y_proba = None
            """
            # --- composition helpers (for precursor gate) ---
            import re
            _comp_pat = re.compile(r'(KDN|F|H|N|S|G|A|s|p)\s*([0-9]+)', re.I)

            def _canon_comp(s: str) -> str:
                """Turn any comp string into canonical 'F,H,N,S,G,KDN,A,s,p' order; omit zeros."""
                if not isinstance(s, str) or not s.strip():
                    return ""
                # accumulate counts (case-insensitive; keep KDN token intact)
                counts = {"F":0,"H":0,"N":0,"S":0,"G":0,"KDN":0,"A":0,"s":0,"p":0}
                for k,v in _comp_pat.findall(s):
                    k = "KDN" if k.upper() == "KDN" else k  # keep KDN as a unit
                    counts[k] = counts.get(k, 0) + int(v)
                # canonical order
                parts = []
                for key in ("F","H","N","S","G","KDN","A","s","p"):
                    n = counts.get(key, 0)
                    if n > 0:
                        parts.append(f"{key}{n}")
                return "".join(parts)

            #new precursor gate maybe
            import re


            def _find_col(df: pd.DataFrame, candidates):
                """Return the first existing column from candidates (case-insensitive)."""
                cols = {c.lower(): c for c in df.columns}
                for name in candidates:
                    if name.lower() in cols:
                        return cols[name.lower()]
                return None

            _LABEL_RE = re.compile(r"""
                (?:F(?P<F>\d+))?      # Fuc
                (?:H(?P<H>\d+))?      # Hex
                (?:N(?P<N>\d+))?      # HexNAc
                (?:S(?P<S>\d+))?      # NeuAc (Sialic acid)
                (?:G(?P<G>\d+))?      # NeuGc
                (?:KDN(?P<KDN>\d+))?  # KDN (explicit token)
            """, re.VERBOSE)

            def _canon_label(s: str) -> str:
                """Normalize labels to FHNSGKDN order; keep only nonzero terms.
                Accepts variants like 'H5N2', 'F1H5N2S1', 'H5N2KDN1', etc."""
                if not isinstance(s, str) or not s:
                    return ""
                m = _LABEL_RE.fullmatch(s)
                if not m:
                    # Try to expand bare 'K1' -> 'KDN1', or tolerate lowercase
                    s2 = s.upper().replace("K", "KDN")
                    m = _LABEL_RE.fullmatch(s2)
                    if not m:
                        return ""
                parts = []
                F = int(m.group("F") or 0)
                H = int(m.group("H") or 0)
                N = int(m.group("N") or 0)
                S = int(m.group("S") or 0)
                G = int(m.group("G") or 0)
                KDN = int(m.group("KDN") or 0)
                if F:   parts.append(f"F{F}")
                if H:   parts.append(f"H{H}")
                if N:   parts.append(f"N{N}")
                if S:   parts.append(f"S{S}")
                if G:   parts.append(f"G{G}")
                if KDN: parts.append(f"KDN{KDN}")
                return "".join(parts)

            def _label_from_counts(row: pd.Series) -> str:
                """Construct FHNSGKDN label from count columns if no label string exists."""
                # Common column names (case-insensitive)
                def g(name): 
                    for cand in [name, name.capitalize()]:
                        if cand in row:
                            return int(row[cand]) if pd.notna(row[cand]) else 0
                    return 0
                F   = g("F") or g("Fuc")
                H   = g("H") or g("Hex")
                N   = g("N") or g("HexNAc")
                S   = g("S") or g("NeuAc")
                G   = g("G") or g("NeuGc")
                KDN = g("KDN")
                parts = []
                if F:   parts.append(f"F{F}")
                if H:   parts.append(f"H{H}")
                if N:   parts.append(f"N{N}")
                if S:   parts.append(f"S{S}")
                if G:   parts.append(f"G{G}")
                if KDN: parts.append(f"KDN{KDN}")
                return "".join(parts)

            def apply_mass_gate(df_pred: pd.DataFrame,
                                df_lib: pd.DataFrame,
                                ppm_threshold: float = 10.0,
                                logger=None) -> pd.DataFrame:
                """
                Adds columns:
                - 'ppm_precursor' (float)
                - 'pred_ok' (bool)
                Requires df_pred to have observed mass in one of:
                ['protonatedmass','precursor_mz','mz']
                and predicted label in one of:
                ['Predicted','Predicted_Label','label','pred_label']
                Library must provide theoretical mass via one of:
                ['Mass','mass','theoretical_mass'] and a label column or counts.
                """
                log = (logger.info if logger else print)

                # --- find columns (case-insensitive) ---
                pred_lab_col = _find_col(df_pred, ["Predicted","Predicted_Label","label","pred_label"])
                obs_mz_col   = _find_col(df_pred, ["protonatedmass","precursor_mz","mz"])
                if not pred_lab_col or not obs_mz_col:
                    log("[gate] ERROR: missing required columns in predictions "
                        f"(label? {bool(pred_lab_col)}, mass? {bool(obs_mz_col)}). Skipping gate.")
                    df_pred = df_pred.copy()
                    df_pred["ppm_precursor"] = pd.NA
                    df_pred["pred_ok"] = False
                    return df_pred

                theo_mass_col = _find_col(df_lib, ["Mass","mass","theoretical_mass"])
                lib_label_col = _find_col(df_lib, ["Label","label","Predicted","composition","Comp","Glycan"])
                if not theo_mass_col:
                    log("[gate] ERROR: library has no theoretical mass column (Mass/mass/theoretical_mass). Skipping gate.")
                    df_pred = df_pred.copy()
                    df_pred["ppm_precursor"] = pd.NA
                    df_pred["pred_ok"] = False
                    return df_pred

                # --- prepare library: ensure a canonical label column ---
                lib = df_lib.copy()
                if lib_label_col:
                    lib["__label__"] = lib[lib_label_col].astype(str).map(_canon_label)
                else:
                    lib["__label__"] = lib.apply(_label_from_counts, axis=1).map(_canon_label)

                # Drop library rows that still have empty labels or missing masses
                lib = lib[pd.notna(lib[theo_mass_col])]
                lib = lib[lib["__label__"] != ""]
                lib = lib.drop_duplicates(subset="__label__", keep="first")

                # --- prepare predictions: canonicalize labels ---
                out = df_pred.copy()
                out["__label__"] = out[pred_lab_col].astype(str).map(_canon_label)

                # --- join ---
                merged = out.merge(lib[["__label__", theo_mass_col]],
                                on="__label__", how="left", suffixes=("", "_lib"))

                # compute ppm only where theoretical mass is available
                theo = merged[theo_mass_col]
                obs  = merged[obs_mz_col]
                with pd.option_context("mode.use_inf_as_na", True):
                    ppm = (obs - theo) / theo * 1e6
                merged["ppm_precursor"] = ppm.where(pd.notna(theo))

                # determine pass/fail
                merged["pred_ok"] = merged["ppm_precursor"].abs() <= ppm_threshold

                # warn if nothing matched
                n_labels_matched = merged[theo_mass_col].notna().sum()
                if n_labels_matched == 0:
                    uniq_preds = sorted(set(out["__label__"]) - {""})
                    log(f"[gate] WARNING: 0 compositions matched between predictions and library. "
                        f"Pred labels example: {uniq_preds[:8]} ... "
                        "Check label format and library columns (Label/F/H/N/S/G/KDN + Mass).")

                # clean up temp column
                merged.drop(columns=["__label__"], inplace=True)
                return merged            


            # --- Precursor gate (optional) ---
            out_path_gate = None                    # make sure this exists for later UI messages
            gate_on = bool(apply_pred_precursor_gate_var.get())  # <- use the SAME var name you used in the UI
            # --- START: reporter-friendly table (unchanged from v10) ---
            # 1) class names
            try:
                class_names = list(le.classes_) if le is not None else list(getattr(model, "classes_", []))
            except Exception:
                class_names = list(getattr(model, "classes_", []))

            # 2) proba (keep this BEFORE gating and summarization)
            proba = None
            if hasattr(model, "predict_proba"):
                proba = model.predict_proba(X)

            pred_df = df.copy()
            pred_df = pred_df.rename(columns={"Predicted_Label": "pred_label"})  # reporter expects 'pred_label'

            proba_cols = None
            if proba is not None and (isinstance(proba, (list, tuple)) or hasattr(proba, "shape")):
                try:
                    import numpy as _np
                    pred_df["proba_vector"] = [_np.asarray(row, dtype=float) for row in proba]
                except Exception:
                    if class_names:
                        proba_cols = [f"proba_{c}" for c in class_names]
                    else:
                        proba_cols = [f"proba_{i}" for i in range(proba.shape[1])]
                    for j, col in enumerate(proba_cols):
                        pred_df[col] = proba[:, j]
            # --- END: reporter-friendly table ---

            # --- Packaging suffixes (so folders/files tell you MS1/gate) ---
            train_feats = get_training_features(model, model_file_path) or []
            ms1_used    = ("protonatedmass" in train_feats) or ("delta_ppm" in train_feats)
            ms1_suffix  = "_MS1feat" if ms1_used else "_noMS1"

            gate_on = bool(apply_pred_precursor_gate_var.get())
            try:
                gate_ppm = float(pred_precursor_ppm_var.get() or 10)
            except Exception:
                gate_ppm = 10.0
            gate_suffix = f"_PG{int(gate_ppm)}ppm" if gate_on else "_noPG"
            method_suffix = ms1_suffix + gate_suffix

            from pathlib import Path
            method_base = _guess_method_basename_for_pack(predict_input_path, df)
            out_dir = Path(os.path.dirname(predict_input_path)) / f"{method_base}{method_suffix}"
            out_dir.mkdir(parents=True, exist_ok=True)

            # --- Precursor gate (robust & optional) ---
            out_path_gate = None
            import pandas as _pd
            import numpy as np
            if gate_on:
                try:

                    # read in-silico and derive composition→[masses]
                    insilico_path = (insilico_csv_var.get() or "").strip()
                    if not insilico_path:
                        raise RuntimeError("No in-silico CSV selected.")
                    lib = _read_any_table(insilico_path)

                    # --- robust reparse / header clean for in-silico CSV ---
                    def _clean_cols(cols):
                        out = []
                        for c in cols:
                            s = str(c).replace("\uFEFF", "").strip()  # strip BOM + spaces
                            out.append(s)
                        return out

                    # 1) clean current headers
                    lib.columns = _clean_cols(lib.columns)

                    # 2) if it looks like a single-column read, try common delimiters
                    if len(lib.columns) == 1:
                        _one = lib.columns[0]
                        # log for debugging
                        if logger: logger.log(f"[Predict] in-silico looked single-col header: {_one!r}; reparsing…")
                        # try comma, semicolon (EU Excel), then tab
                        for _sep in (",", ";", "\t"):
                            try:
                                _tmp = _pd.read_csv(insilico_path, sep=_sep, engine="python")
                                _tmp.columns = _clean_cols(_tmp.columns)
                                if len(_tmp.columns) > 1:
                                    lib = _tmp
                                    if logger: logger.log(f"[Predict] reparsed with sep='{_sep}': {lib.columns.tolist()}")
                                    break
                            except Exception as _e:
                                pass

                    # 3) build a case-insensitive lookup AFTER cleaning
                    norm = {c.lower(): c for c in lib.columns}

                    # pick mass column robustly
                    mass_col = None
                    for key in ("theoretical_mass", "mass"):
                        if key in norm:
                            mass_col = norm[key]
                            break
                    # extra tolerance (e.g., stray spaces/case)
                    if mass_col is None:
                        for c in lib.columns:
                            if c.strip().lower() in ("mass", "theoretical_mass"):
                                mass_col = c
                                break

                    if logger:
                        logger.log(f"[Predict] in-silico columns (cleaned): {lib.columns.tolist()}")
                        logger.log(f"[Predict] chosen mass_col: {mass_col}")


                    # normalize headers
                    norm = {c.lower().strip(): c for c in lib.columns}
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
                        raise RuntimeError("In-silico CSV must have a Mass or theoretical_mass column.")

                    # if no comp string, synthesize from counts
                    # If no composition string, synthesize one from monomer counts (FHN with modifiers A/s/p, plus G/K if present)
                    if comp_col is None:
                        def _get(row, key):
                            col = norm.get(key)
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
                            if KDN:   parts.append(f"K{KDN}")
                            if HexA:  parts.append(f"A{HexA}")    # <-- HexA is 'A' in your schema
                            if SO3:   parts.append(f"s{SO3}")
                            if PO3H:  parts.append(f"p{PO3H}")
                            return "".join(parts)

                        lib = lib.copy()
                        lib["__comp"] = lib.apply(_counts_to_comp, axis=1)
                        comp_col = "__comp"

                    #new patch
                                        
                    # ADD: canonicalize library strings even if the CSV already has a composition column
                    def _canon_safe(x):
                        try:
                            return _canon_comp(str(x))
                        except Exception:
                            return str(x).strip()

                    if comp_col is None:
                        # (your existing counts→string code stays as-is)
                        ...
                        lib = lib.copy()
                        lib["__comp"] = lib.apply(_counts_to_comp, axis=1)
                        comp_col = "__comp"
                    else:
                        # <<< NEW: canonicalize provided composition strings so they match prediction labels
                        lib = lib.copy()
                        lib[comp_col] = lib[comp_col].astype(str).map(_canon_safe)

                    # build lookup composition → list of masses  (unchanged)
                    comp2masses = {}
                    for comp, mass in zip(lib[comp_col].astype(str), pd.to_numeric(lib[mass_col], errors="coerce")):
                        if pd.isna(mass):
                            continue
                        comp_key = _canon_safe(comp)
                        comp2masses.setdefault(comp_key, []).append(float(mass))

                    # build lookup composition → list of masses
                    #comp2masses = {}
                    #for comp, mass in zip(lib[comp_col].astype(str), pd.to_numeric(lib[mass_col], errors="coerce")):
                    #    if pd.isna(mass) or not comp:
                    #        continue
                    #    comp2masses.setdefault(comp, []).append(float(mass))

                    # apply gate (skip Non-glycan)
                    is_ng = pred_df["pred_label"].astype(str).eq("Non-glycan")
                    pred_df["ppm_precursor"] = np.nan
                    # CHANGE: default to False; only set True for glycans that pass
                    pred_df["pred_ok"] = False
                    #pred_df["pred_ok"] = True

                    #debug lines
                    # --- GATE DIAGNOSTICS ---
                    try:
                        # same canonicalizer for both sides
                        def _canon(x):
                            try:    return _canon_comp(x)
                            except: return str(x).strip()

                        # build a quick “masses found?” map
                        # If you created comp2masses above, reuse it; otherwise reconstruct it here:
                        # comp2masses = {...}  # already built earlier

                        gly = pred_df[ pred_df["pred_label"].astype(str) != "Non-glycan" ].copy()
                        gly["canon_pred"]  = gly["pred_label"].astype(str).map(_canon)
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

                        # write a small debug file of dropped glycans
                        debug_cols = ["MS2scan_no","pred_label","canon_pred","protonatedmass",
                                    "ppm_precursor","masses_found","kept"]
                        gate_dbg = gly[~gly["kept"]][debug_cols].head(1000)  # cap size
                        dbg_path = (out_dir / (Path(predict_input_path).stem + "_gate_debug.csv")).as_posix()
                        gate_dbg.to_csv(dbg_path, index=False)
                        if logger: logger.log(f"[Predict][gate] wrote debug rows → {dbg_path}")
                    except Exception as _e:
                        if logger: logger.log(f"[Predict][gate] diagnostics failed: {_e}")

                    # --- EXTRA DIAGNOSTICS: for no_lut rows, fetch Mass candidates by counts-match ---
                    try:
                        import re

                        # 1) Normalize + numericize in-silico counts once
                        lib_counts = lib.copy()

                        def _nz_int(colname):
                            if colname in lib_counts.columns:
                                lib_counts[colname] = pd.to_numeric(lib_counts[colname], errors="coerce").fillna(0).astype(int)
                            else:
                                lib_counts[colname] = 0

                        # standard headers we saw in your file
                        _nz_int("Hex")
                        _nz_int("HexNAc")
                        _nz_int("Fuc")
                        _nz_int("NeuAc")
                        _nz_int("NeuGc")
                        _nz_int("KDN")
                        _nz_int("HexA")
                        _nz_int("SO3")
                        _nz_int("PO3H")

                        # Mass column (we already picked mass_col earlier)
                        lib_counts["__Mass__"] = pd.to_numeric(lib_counts[mass_col], errors="coerce")

                        # 2) Parse a predicted label like "F1H4N3A1s1p1" → counts
                        def _parse_counts(label: str):
                            d = {"F":0,"H":0,"N":0,"S":0,"G":0,"K":0,"A":0,"s":0,"p":0}
                            for m in re.finditer(r'([FHNSGKAsp])(\d+)', str(label)):
                                ch, num = m.group(1), int(m.group(2))
                                if ch in d:
                                    d[ch] += num
                            return d
                        #fix GPT's loop, you forgot sialic acid, that's mandatory, not optional okay?
                        # 3) Build a finder: from counts → candidate Mass list
                        def _find_masses_by_counts(cnt):
                            m = lib_counts.loc[
                                (lib_counts["Fuc"]   == cnt["F"]) &
                                (lib_counts["Hex"]   == cnt["H"]) &
                                (lib_counts["HexNAc"]== cnt["N"]) &
                                (lib_counts["NeuAc"]   == cnt["S"]) &
                                (lib_counts["NeuGc"]   == cnt["G"]) &
                                (lib_counts["KDN"]   == cnt["K"]) &
                                (lib_counts["HexA"]  == cnt["A"]) &
                                (lib_counts["SO3"]   == cnt["s"]) &
                                (lib_counts["PO3H"]  == cnt["p"]),
                                "__Mass__"
                            ].dropna()
                            return sorted(m.unique().tolist())

                        # 4) Take the first 30 no-lut rows and attach candidate Masses
                        no_lut_head = gly[gly["masses_found"] == 0].copy().head(30)
                        if not no_lut_head.empty:
                            no_lut_head["counts"] = no_lut_head["canon_pred"].map(_parse_counts)
                            no_lut_head["Mass_candidates"] = no_lut_head["counts"].map(_find_masses_by_counts)

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
                            if logger: logger.log(f"[Predict][gate] wrote by-counts no-lut debug → {nolut_bycounts_path}")

                    except Exception as _e:
                        if logger: logger.log(f"[Predict][gate] counts-mass diagnostics failed: {_e}")

                    """
                    idx = (~is_ng).to_numpy().nonzero()[0]
                    if len(idx):
                        # canonicalize predicted comps (uses your _canon_comp if present)
                        def _canon(x):
                            ...
                        comps = pred_df.loc[idx, "pred_label"].astype(str).map(_canon)
                        obs   = pd.to_numeric(pred_df.loc[idx, "protonatedmass"], errors="coerce")
                        best  = np.full(len(idx), np.nan)
                        keep  = np.zeros(len(idx), dtype=bool)
                        for j, (c, pm) in enumerate(zip(comps, obs)):
                            ...
                        pred_df.loc[idx, "ppm_precursor"] = best
                        pred_df.loc[idx, "pred_ok"] = keep
                    
                    if len(idx):
                        # canonicalize predicted comps (uses your _canon_comp if present)
                        def _canon(x):
                            try:
                                return _canon_comp(x)
                            except Exception:
                                return str(x)
                        comps = pred_df.loc[idx, "pred_label"].astype(str).map(_canon)
                        obs   = pd.to_numeric(pred_df.loc[idx, "protonatedmass"], errors="coerce")
                        best  = np.full(len(idx), np.nan)
                        keep  = np.zeros(len(idx), dtype=bool)
                        for j, (c, pm) in enumerate(zip(comps, obs)):
                            if not c or np.isnan(pm):
                                continue
                            masses = comp2masses.get(c, [])
                            if not masses:
                                continue
                            best[j] = min(abs((pm - tm)/tm*1e6) for tm in masses)
                            keep[j] = (best[j] <= gate_ppm)
                        pred_df.loc[idx, "ppm_precursor"] = best
                        pred_df.loc[idx, "pred_ok"] = keep
                    """
                    # --- compute ppm for glycan rows only (robust: string LUT + counts fallback) ---
                    import numpy as _np
                    import pandas as _pd

                    # 1) Canonical prediction strings, but keep the raw too for logging
                    pred_df["__canon_pred__"] = pred_df["pred_label"].astype(str).map(_canon_safe)

                    gly_mask = (~is_ng).to_numpy()
                    idx = _np.nonzero(gly_mask)[0]

                    labels = pred_df.loc[gly_mask, "__canon_pred__"].astype(str).to_numpy()
                    obs    = _pd.to_numeric(pred_df.loc[gly_mask, "protonatedmass"], errors="coerce").to_numpy()

                    best_ppm = _np.full(labels.shape[0], _np.nan, dtype=float)
                    keep     = _np.zeros(labels.shape[0], dtype=bool)

                    # Reuse your counts parser and finder from diagnostics
                    def _parse_counts(label: str):
                        import re
                        d = {"F":0,"H":0,"N":0,"S":0,"G":0,"K":0,"A":0,"s":0,"p":0}
                        for m in re.finditer(r'([FHNSGKAsp])(\d+)', str(label)):
                            ch, num = m.group(1), int(m.group(2))
                            if ch in d:
                                d[ch] += num
                        return d

                    def _masses_by_counts(lbl: str):
                        cnt = _parse_counts(lbl)
                        m = lib_counts.loc[
                            (lib_counts["Fuc"]   == cnt["F"]) &
                            (lib_counts["Hex"]   == cnt["H"]) &
                            (lib_counts["HexNAc"]== cnt["N"]) &
                            (lib_counts["NeuAc"] == cnt["S"]) &
                            (lib_counts["NeuGc"] == cnt["G"]) &
                            (lib_counts["KDN"]   == cnt["K"]) &
                            (lib_counts["HexA"]  == cnt["A"]) &
                            (lib_counts["SO3"]   == cnt["s"]) &
                            (lib_counts["PO3H"]  == cnt["p"]),
                            "__Mass__"
                        ].dropna()
                        return sorted(m.unique().tolist())

                    def _masses_for_label(lbl: str):
                        # try string LUT first
                        m = comp2masses.get(lbl, [])
                        if m:
                            return m
                        # fallback to counts-based match
                        return _masses_by_counts(lbl)

                    gppm = float(gate_ppm)
                    for j, (lbl, pm) in enumerate(zip(labels, obs)):
                        if not lbl or _np.isnan(pm):
                            continue
                        masses = _masses_for_label(lbl)
                        if not masses:
                            # verbose log once for the first few misses
                            if logger and j < 5:
                                logger.log(f"[gate][miss] no library masses for label '{lbl}' (pm={pm:.6f})")
                            continue

                        m = _np.asarray(masses, dtype=float)
                        diffs = (pm - m) / m * 1e6
                        k = _np.nanargmin(_np.abs(diffs))
                        best_ppm[j] = float(diffs[k])

                        passed = abs(best_ppm[j]) <= gppm
                        keep[j] = passed

                        # DEBUG: log every hit (so we can see matches in the terminal)
                        try:
                            scan_no = pred_df.loc[idx[j], "MS2scan_no"]
                        except Exception:
                            scan_no = None
                        if logger:
                            logger.log(
                                f"[gate][hit] scan={scan_no} label='{lbl}' pm={pm:.6f} "
                                f"bestTheo={m[k]:.6f} ppm={best_ppm[j]:.2f} → {'PASS' if passed else 'FAIL'}"
                            )

                    # write back only for glycan rows
                    pred_df.loc[idx, "ppm_precursor"] = best_ppm
                    pred_df.loc[idx, "pred_ok"]       = keep

                    """
                    # --- compute ppm for glycan rows only (explicit protonatedmass, canonical labels) ---
                    pred_df["__canon_pred__"] = pred_df["pred_label"].astype(str).map(_canon_safe)

                    gly_mask = (~is_ng).to_numpy()
                    idx = np.nonzero(gly_mask)[0]

                    labels = pred_df.loc[gly_mask, "__canon_pred__"].astype(str).to_numpy()
                    obs    = pd.to_numeric(pred_df.loc[gly_mask, "protonatedmass"], errors="coerce").to_numpy()

                    best_ppm = np.full(labels.shape[0], np.nan, dtype=float)
                    keep     = np.zeros(labels.shape[0], dtype=bool)

                    for j, (lbl, pm) in enumerate(zip(labels, obs)):
                        if not lbl or np.isnan(pm):
                            continue
                        masses = comp2masses.get(lbl, [])
                        if not masses:
                            continue
                        m = np.asarray(masses, dtype=float)
                        diffs = (pm - m) / m * 1e6
                        k = np.nanargmin(np.abs(diffs))
                        best_ppm[j] = float(diffs[k])
                        keep[j] = (abs(best_ppm[j]) <= float(gate_ppm))  # gate_ppm is your existing threshold

                    # write back only for glycan rows
                    pred_df.loc[idx, "ppm_precursor"] = best_ppm
                    pred_df.loc[idx, "pred_ok"] = keep
                    """
                    # gated copy for disk & for report
                    out_path_gate = (out_dir / (Path(predict_input_path).stem + f"_predicted{ms1_suffix}_PG{int(gate_ppm)}ppm_only.csv")).as_posix()
                    df["pred_ok"] = pred_df["pred_ok"].values  # keep in the full table too
                    #now ppm falls to csv, not in logger or terminal only
                    df["ppm_precursor"] = pred_df["ppm_precursor"].values   # <<< add this
                    df.to_csv(out_dir / (Path(predict_input_path).stem + f"_predicted{method_suffix}.csv"), index=False)
                    pred_df_for_report = pred_df[pred_df["pred_ok"]].copy()
                    #convenient passed class count
                    # Save a simple gated count table: composition, n
                    gate_counts = (pred_df_for_report
                                .groupby("pred_label", dropna=False)
                                .size()
                                .reset_index(name="count")
                                .sort_values("count", ascending=False))
                    gate_counts_path = (out_dir / (Path(predict_input_path).stem + "_PG_counts.csv")).as_posix()
                    gate_counts.to_csv(gate_counts_path, index=False)
                    if logger: logger.log(f"[Predict][gate] wrote gated counts → {gate_counts_path}")


                except Exception as _e:
                    # Gate skipped → just report everything
                    if logger: logger.log(f"[Predict] Precursor gate skipped: {_e}")
                    pred_df_for_report = pred_df.copy()
            else:
                pred_df_for_report = pred_df.copy()

            #
            # --- Harmonize probability vectors so all rows have same length ---
            import json, ast
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
            params = ReportParams(tau=0.60, margin=0.05, topk=5, sample_cols=("experiment_title", "sample_name"))
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
            messagebox.showinfo("Prediction Complete", f"{ms1_line}\n{gate_line}\n\nPacked into folder:\n{out_dir}\n\nMain table:\n{out_path_all}{extra}")
            messagebox.showinfo("Prediction Complete",
                f"MS1 feature (expected by model): {'YES' if ms1_used else 'NO'}\n"
                f"Precursor gate: {'ON' if gate_on else 'OFF'}{(' ('+str(gate_ppm)+' ppm)') if gate_on else ''}\n\n"
                f"Packed into folder:\n{out_dir}\n\nMain table:\n{out_path_all}{extra}"
            )

            ## --- END: Prediction summary integration ---
            #out_path = os.path.splitext(predict_input_path)[0] + "_predicted.csv"
            #df.to_csv(out_path, index=False)
            messagebox.showinfo("Prediction Complete", f"Predictions saved to:\n{out_path_all}")
        except Exception as e:
            messagebox.showerror("Prediction Failed", str(e))

    #20250915 replace extract_fragment_masses for csv-only ion list
    # --- replace extract_fragment_masses(...) with this ---
    def read_fragment_masses_any(ion_path: str, sheet_name: str = "ionlist"):
        """
        Load fragment masses from either CSV/TSV (expects a 'mass' column)
        or from an Excel sheet (default 'ionlist').
        """
        if not ion_path or not os.path.exists(ion_path):
            raise FileNotFoundError(f"Ion list not found: {ion_path}")

        ext = os.path.splitext(ion_path)[1].lower()
        if ext in (".csv", ".tsv", ".txt"):
            # try TSV if it looks like one; else default to CSV
            sep = "\t" if ext == ".tsv" else ","
            df = pd.read_csv(ion_path, sep=sep)
            if "mass" not in df.columns:
                raise ValueError(f"Ion list '{ion_path}' must contain a 'mass' column.")
            return df["mass"].dropna().astype(float).tolist()
        else:
            xls = pd.ExcelFile(ion_path, engine="openpyxl")
            ion_df = xls.parse(sheet_name or "ionlist")
            if "mass" not in ion_df.columns:
                raise ValueError(f"Ion sheet '{sheet_name}' missing 'mass' column.")
            return ion_df["mass"].dropna().astype(float).tolist()

    def _short_id(s: str) -> str:
        import hashlib
        return hashlib.sha1((s or "").encode("utf-8")).hexdigest()[:8]

    #20250915 
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
        raw_csv_p = to_native_path(raw_csv_s) if raw_csv_s else None
        if not raw_csv_p or not raw_csv_p.exists():
            picked = _pick_file_cli_or_gui("Select the *converted* CSV (ms2_*.csv)")
            if not picked:
                raise FileNotFoundError("No converted CSV provided.")
            raw_csv_p = to_native_path(picked)
            # write back to method.json in POSIX form
            parents["converted_csv"] = to_posix_str(raw_csv_p)
            m["parents"] = parents
            with open(method_path, "w", encoding="utf-8") as f:
                json.dump(m, f, indent=2, ensure_ascii=False)

        # Validate ion list path
        ion_path_p = to_native_path(ion_path_s)
        if not ion_path_p.exists():
            picked = _pick_file_cli_or_gui("Select ion list (CSV/XLSX)",
                                        (("CSV", "*.csv"), ("Excel", "*.xlsx;*.xls"), ("All files", "*.*")))
            if not picked:
                raise FileNotFoundError("No ion list provided.")
            ion_path_p = to_native_path(picked)
            ionblock["path"] = to_posix_str(ion_path_p)
            m["ionlist"] = ionblock
            with open(method_path, "w", encoding="utf-8") as f:
                json.dump(m, f, indent=2, ensure_ascii=False)

        # Load fragments (your existing smart loader)
        frags = read_fragment_masses_any(ion_path_p.as_posix(), sheet_name=ion_sheet)

        # Use your existing extractor
        feature_df = extract_ion_intensities(raw_csv_p.as_posix(), frags, ppm=float(ppm))

        # Best-effort UID (no-op if columns missing)
        try:
            def _short_id(s: str) -> str:
                return hashlib.sha1((s or "").encode("utf-8")).hexdigest()[:8]
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

    """
    def create_unlabeled_from_method(method_path: str, default_ppm: int | None = None) -> str:
        with open(method_path, "r", encoding="utf-8") as f:
            m = json.load(f)

        # Resolve inputs from method.json (PL has no manual sheet)
        raw_csv = (m.get("parents") or {}).get("converted_csv")
        ion_path = (m.get("ionlist") or {}).get("path")
        ion_sheet = (m.get("ionlist") or {}).get("sheet") or "ionlist"
        ppm = default_ppm or (m.get("ionlist") or {}).get("ppm_tolerance") or 20

        if not raw_csv or not os.path.exists(raw_csv):
            raise FileNotFoundError(f"Converted CSV not found (method): {raw_csv}")
        if not ion_path or not os.path.exists(ion_path):
            raise FileNotFoundError(f"Ion list not found (method): {ion_path}")

        frags = read_fragment_masses_any(ion_path, sheet_name=ion_sheet)
        feat_df = extract_ion_intensities(raw_csv, frags, ppm=float(ppm))

        # Add UID if possible
        try:
            exp_id = _short_id(m.get("experiment_title", ""))
            samp_id = _short_id(m.get("sample_name", ""))
            feat_df = feat_df.copy()
            feat_df["UID"] = (
                feat_df["MS2scan_no"].astype(int).astype(str).str.zfill(6)
                .map(lambda s: f"{exp_id}:{samp_id}:{s}")
            )
        except Exception:
            pass

        out_path = os.path.splitext(raw_csv)[0] + f"_unlabeled_ppm{ppm}.csv"
        feat_df.to_csv(out_path, index=False)
        return out_path
    """
    #20250915 end here?    
    """
    #creating unlabeled datasets
    def extract_fragment_masses(excel_path, sheet_name="ionlist"):
        xls = pd.ExcelFile(excel_path, engine="openpyxl")
        #sheet_name = ion_sheet_name #or xls.sheet_names[1]
        ion_df = xls.parse(sheet_name)
        fragment_masses  = ion_df["mass"].dropna().astype(float).tolist()#df.iloc[:, 0].dropna().astype(float).tolist()
        return fragment_masses

    def extract_ion_intensities(tsv_path, fragment_masses, ppm=20.0):
        df = pd.read_csv(tsv_path, sep="\t")
        result = []
        #base_cols = df[["protonatedmass", "MS2scan_no"]].copy()
        import numpy as np
        for idx, row in df.iterrows():
            try:
                peaks = list(eval(row["peaklist"]))
                intensities = list(eval(row["peakintensity"]))
            except Exception:
                continue

            peak_array = np.array(peaks)
            intensity_array = np.array(intensities)
            feature_row = {
                "protonatedmass": row["protonatedmass"],
                "MS2scan_no": int(row["MS2scan_no"])
            }

            for target in fragment_masses:
                ppm_tol = target * ppm / 1e6
                mask = np.abs(peak_array - target) <= ppm_tol
                if np.any(mask):
                    intensity = float(np.max(intensity_array[mask]))
                    feature_row[str(target)] = np.log10(intensity + 1)
                else:
                    feature_row[str(target)] = 1.0 #0.0

            result.append(feature_row)

        return pd.DataFrame(result)
        """
    #20250929 ver
    # --- creating unlabeled datasets (PL or manual) ---

    def extract_fragment_masses(path, sheet_name="ionlist"):
        """
        Load a column of fragment masses from either:
        • Excel: prefers a sheet named 'ionlist' (or the provided sheet_name)
        • CSV:   uses a column named 'mass' (case-insensitive)
        Returns: list[float]
        """
        import pandas as pd
        p = (path or "").strip()
        if not p:
            return []

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
        import numpy as np
        import pandas as pd

        df = pd.read_csv(tsv_like_path, sep="\t")  # your converter writes TSV
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
                    feature_row[str(target)] = float(np.log10(np.max(intensity_array[mask]) + 1.0))
                else:
                    feature_row[str(target)] = 1.0  # keep “1.0” default you were using

            out_rows.append(feature_row)

        return pd.DataFrame(out_rows)

    #20251001 extra fix of insilico csv gate
    # --- [PREDICTION: precursor-gate helpers] ------------------------------------
    import re
    import pandas as _pd
    import numpy as _np

    def _ppm_delta(obs, theo):
        if theo == 0 or _np.isnan(theo) or _np.isnan(obs): 
            return _np.nan
        return (float(obs) - float(theo)) / float(theo) * 1e6

    # Parse flexible composition strings e.g. "H7N2F1S1A1s1p0", "F1H7N2S1A1", case-insensitive,
    # order-agnostic; missing parts -> 0. Supports A (HexA), s (SO3), p (PO3H).
    _comp_pat = re.compile(r'(F|H|N|S|G|KDN|A|s|p)\s*([0-9]+)', re.I)

    def _canonicalize_comp_string(s: str) -> str:
        if not isinstance(s, str) or not s.strip():
            return ""
        counts = {"F":0,"H":0,"N":0,"S":0,"G":0,"KDN":0,"A":0,"s":0,"p":0}
        for k,v in _comp_pat.findall(s):
            k = "KDN" if k.upper()=="KDN" else k  # keep KDN token
            counts[k] = counts.get(k, 0) + int(v)
        # canonical short form: F,H,N,S,G,KDN,A,s,p (omit zeros)
        parts = []
        for key in ("F","H","N","S","G","KDN","A","s","p"):
            n = counts.get(key,0)
            if n>0: parts.append(f"{key}{n}")
        return "".join(parts)

    def _load_insilico_for_gate(path: str) -> pd.DataFrame:
        """
        Accepts:
        - New wide CSV: Hex,HexNAc,NeuAc,NeuGc,KDN,Fuc[,HexA,SO3,PO3H], Mass
        - Legacy 2-col: comp_str, Mass
        - Already-canonical: composition, theoretical_mass
        Returns a DataFrame with at least: composition (canonical str), theoretical_mass (float)
        """
        df = pd.read_csv(path)
        cols = {c.lower(): c for c in df.columns}

        # Case 1: already has 'composition' + 'theoretical_mass' (or 'mass')
        if "composition" in cols and ("theoretical_mass" in cols or "mass" in cols):
            comp_col = cols["composition"]
            mass_col = cols.get("theoretical_mass", cols.get("mass"))
            out = pd.DataFrame({
                "composition": df[comp_col].astype(str).map(_canonicalize_comp_string),
                "theoretical_mass": pd.to_numeric(df[mass_col], errors="coerce")
            })
            return out.dropna(subset=["theoretical_mass"])

        # Case 2: wide numeric columns
        wide_keys = ["hex","hexnac","neuac","neugc","kdn","fuc"]
        if all(k in cols for k in wide_keys) and ("mass" in cols or "theoretical_mass" in cols):
            mass_col = cols.get("theoretical_mass", cols.get("mass"))
            # optional modifiers
            hexA = df[cols.get("hexa","HexA")] if "hexa" in cols else 0
            so3  = df[cols.get("so3","SO3")]   if "so3"  in cols else 0
            po3h = df[cols.get("po3h","PO3H")] if "po3h" in cols else 0
            # build canonical label
            def _build_row(row):
                parts = []
                if row["Fuc"]>0:    parts.append(f"F{int(row['Fuc'])}")
                if row["Hex"]>0:    parts.append(f"H{int(row['Hex'])}")
                if row["HexNAc"]>0: parts.append(f"N{int(row['HexNAc'])}")
                if row["NeuAc"]>0:  parts.append(f"S{int(row['NeuAc'])}")
                if row["NeuGc"]>0:  parts.append(f"G{int(row['NeuGc'])}")
                if row["KDN"]>0:    parts.append(f"KDN{int(row['KDN'])}")
                A = int(row.get("HexA",0)); s = int(row.get("SO3",0)); p = int(row.get("PO3H",0))
                if A>0: parts.append(f"A{A}")
                if s>0: parts.append(f"s{s}")
                if p>0: parts.append(f"p{p}")
                return "".join(parts)
            tmp = df.rename(columns={
                cols["hex"]:"Hex", cols["hexnac"]:"HexNAc", cols["neuac"]:"NeuAc",
                cols["neugc"]:"NeuGc", cols["kdn"]:"KDN", cols["fuc"]:"Fuc"
            }).copy()
            tmp["HexA"] = pd.to_numeric(hexA, errors="coerce").fillna(0).astype(int)
            tmp["SO3"]  = pd.to_numeric(so3,  errors="coerce").fillna(0).astype(int)
            tmp["PO3H"] = pd.to_numeric(po3h, errors="coerce").fillna(0).astype(int)
            out = pd.DataFrame({
                "composition": tmp.apply(_build_row, axis=1),
                "theoretical_mass": pd.to_numeric(tmp[mass_col], errors="coerce")
            })
            return out.dropna(subset=["theoretical_mass"])

        # Case 3: legacy 2-column like "(H,N,S,G,KDN,F),Mass"
        if len(df.columns) == 2:
            comp_col, mass_col = df.columns[0], df.columns[1]
            comp = df[comp_col].astype(str).map(_canonicalize_comp_string)
            mass = pd.to_numeric(df[mass_col], errors="coerce")
            return pd.DataFrame({"composition": comp, "theoretical_mass": mass}).dropna(subset=["theoretical_mass"])

        # Fallback → raise a clear message the UI can show
        raise ValueError("In-silico CSV must have composition/theoretical_mass OR standard wide columns.")


    #20251001 fix unbound error and etc
    def create_unlabeled_dataset():
        """
        Builds unlabeled feature CSVs for each eligible sample in an experiment JSON.
        Eligibility:
        - has a converted MS2 file at samples[*].csv
        - AND has either samples[*].excel OR samples[*].ionlist_path (+ optional ion_sheet)
        """
        import re
        from pathlib import Path

        exp_path = filedialog.askopenfilename(filetypes=[("Experiment JSON", "*.json")])
        if not exp_path:
            return

        def _sample_in_filename(sample_id, path):
            if not path:
                return False
            name = os.path.basename(path).lower()
            must_have = []
            for tok in sample_id.lower().split("_"):
                if tok in {"u937", "cells", "ng"}:
                    continue
                must_have.append(tok)
            return all(tok in name for tok in must_have)

        try:
            with open(exp_path, "r", encoding="utf-8") as f:
                exp = json.load(f)
        except Exception as e:
            messagebox.showerror("Failed to Load JSON", str(e))
            return

        ppm_value = simpledialog.askinteger(
            "PPM Tolerance", "Enter PPM tolerance (e.g., 20):", minvalue=1, maxvalue=200
        )
        if ppm_value is None:
            return

        made = 0
        missing = []
        seen_out = set()
        made_paths = []

        for sample_id, sample_info in (exp.get("samples") or {}).items():
            raw_csv   = sample_info.get("csv")
            ion_src   = sample_info.get("excel") or sample_info.get("ionlist_path")
            ion_sheet = sample_info.get("ion_sheet") or "ionlist"

            # 1) presence check first (avoids calling basename on None)
            if not raw_csv or not ion_src:
                missing.append(sample_id)
                continue

            # 2) filename sanity check (optional)
            if not _sample_in_filename(sample_id, raw_csv):
                missing.append(f"{sample_id} (csv looks mismatched: {os.path.basename(raw_csv)})")
                continue

            try:
                # fragment list
                ion_list = extract_fragment_masses(ion_src, sheet_name=ion_sheet)
                if not ion_list:
                    raise ValueError(f"No 'mass' column found in {os.path.basename(ion_src)}")

                # features
                feats = extract_ion_intensities(raw_csv, ion_list, ppm=ppm_value)

                # 3) output next to CSV (define out_dir/base *inside* the loop)
                out_dir = Path(raw_csv).resolve().parent
                base    = Path(raw_csv).stem
                out_path = out_dir / f"{base}_unlabeled_ppm{ppm_value}.csv"

                # collision guard (two samples may share same CSV)
                if str(out_path) in seen_out:
                    safe_id = re.sub(r'[^A-Za-z0-9._-]+', '_', sample_id)
                    out_path = out_dir / f"{base}__{safe_id}_unlabeled_ppm{ppm_value}.csv"

                # write
                feats.to_csv(str(out_path), index=False, encoding="utf-8")

                # remember for this sample and for duplicate guard
                sample_info["unlabeled_csv"] = str(out_path)
                seen_out.add(str(out_path))
                made_paths.append(str(out_path))
                made += 1

            except Exception as e:
                print(f"[ERROR] Unlabeled build failed for '{sample_id}': {e}")

        # tell user what happened
        if made == 0:
            messagebox.showwarning(
                "No Samples Processed",
                "No eligible samples were found in this experiment file.\n"
                "Tip: each sample needs 'csv' + ('excel' or 'ionlist_path')."
            )
        else:
            # persist updated experiment file with recorded unlabeled paths
            try:
                with open(exp_path, "w", encoding="utf-8") as f:
                    json.dump(exp, f, indent=2, ensure_ascii=False)
            except Exception:
                pass
            msg = [f"Created {made} unlabeled dataset(s)."]
            if missing:
                msg.append(f"Skipped (missing/mismatch): {', '.join(missing)}")
            messagebox.showinfo("Done", "\n".join(msg))

    """
    #20250915 ver
    def create_unlabeled_dataset():
        # Let the user pick either an experiment .exp.json OR a single sample .method.json
        exp_or_method = filedialog.askopenfilename(filetypes=[("JSON", "*.json")])
        if not exp_or_method:
            return

        # Try to detect which JSON it is
        try:
            with open(exp_or_method, "r", encoding="utf-8") as f:
                J = json.load(f)
        except Exception as e:
            messagebox.showerror("Failed to Load JSON", str(e))
            return

        # If it looks like a method.json (has ionlist + parents), do ONE sample now
        if isinstance(J, dict) and J.get("ionlist") and J.get("parents"):
            # PPM: prefer method value; if user inputs, it overrides
            ppm_default = (J.get("ionlist") or {}).get("ppm_tolerance")
            ppm_value = simpledialog.askinteger(
                "PPM Tolerance",
                f"Enter PPM tolerance (default {ppm_default or 20}):",
                minvalue=1, maxvalue=100
            )
            ppm_value = ppm_value or ppm_default or 20

            try:
                out_path = create_unlabeled_from_method(exp_or_method, default_ppm=ppm_value)
                messagebox.showinfo("Done", f"Unlabeled dataset saved:\n{out_path}")
            except Exception as e:
                messagebox.showerror("Unlabeled Build Failed", str(e))
            return

        # Else treat as experiment .exp.json (multi-sample)
        ppm_value = simpledialog.askinteger("PPM Tolerance", "Enter PPM tolerance (e.g. 20):", minvalue=1, maxvalue=100)
        if ppm_value is None:
            return

        # Expect exp["samples"] structure; each sample may have either:
        #  (A) csv + excel  → use Excel ion sheet (legacy manual route)
        #  (B) json (method.json) → use external ion list from the method (PL route)
        changed = False
        samples = J.get("samples", {}) if isinstance(J, dict) else {}

        #20250917 ver
        # inside create_unlabeled_dataset(), in the loop over samples:
        for sample_id, sample_info in samples.items():
            try:
                raw_csv = sample_info.get("csv")
                annotation_excel = sample_info.get("excel")
                method_json = sample_info.get("json") or sample_info.get("method_json")
                out_path = None

                if annotation_excel and raw_csv:
                    # (A) legacy manual: Excel ion sheet
                    ion_list = read_fragment_masses_any(annotation_excel, sheet_name="ionlist")
                    feature_df = extract_ion_intensities(raw_csv, ion_list, ppm=float(ppm_value))
                    out_path = Path(raw_csv).with_name(Path(raw_csv).stem + f"_unlabeled_ppm{ppm_value}.csv")

                elif method_json and os.path.exists(method_json):
                    # (B) PL route: method.json
                    out_path = create_unlabeled_from_method(method_json, default_ppm=ppm_value)

                elif sample_info.get("ionlist_path") and raw_csv:
                    # (C) NEW: converted CSV + external CSV ion list (no Excel/method)
                    ion_path = to_native_path(sample_info["ionlist_path"])
                    if not ion_path.exists():
                        # ask once if missing, then persist back
                        picked = _pick_file_cli_or_gui("Select ion list (CSV/XLSX)",
                                                    (("CSV","*.csv"),("Excel","*.xlsx;*.xls"),("All","*.*")))
                        if not picked:
                            raise FileNotFoundError("ionlist_path not found and no replacement provided.")
                        ion_path = to_native_path(picked)
                        sample_info["ionlist_path"] = to_posix_str(ion_path)

                    raw_csv_p = to_native_path(raw_csv)
                    if not raw_csv_p.exists():
                        picked = _pick_file_cli_or_gui("Select the *converted* CSV (ms2_*.csv)")
                        if not picked:
                            raise FileNotFoundError("csv not found and no replacement provided.")
                        raw_csv_p = to_native_path(picked)
                        sample_info["csv"] = to_posix_str(raw_csv_p)

                    # load fragments (CSV ion list has a 'mass' column)
                    frags = read_fragment_masses_any(ion_path.as_posix(), sheet_name="ionlist")
                    feature_df = extract_ion_intensities(raw_csv_p.as_posix(), frags, ppm=float(ppm_value))

                    # Optional: embed method pointer column if present
                    if method_json:
                        feature_df.insert(0, "method_json", Path(method_json).resolve().as_posix())

                    out_path = raw_csv_p.with_name(raw_csv_p.stem + f"_unlabeled_ppm{ppm_value}.csv")

                    # Optional (recommended): write a tiny method.json for reproducibility/packing
                    if not method_json:
                        mpath = raw_csv_p.with_suffix("").with_name(f"{sample_id}.fromexp.method.json")
                        mini = {
                            "version": "1.0",
                            "dataset_type": "pseudolabel",
                            "experiment_title": J.get("experiment", ""),
                            "sample_name": sample_id,
                            "parents": {
                                "converted_csv": to_posix_str(raw_csv_p),
                                "pseudolabels_tsv": "(unknown)",
                                "trainable_csv": "(unknown)"
                            },
                            "ionlist": {
                                "source": "external",
                                "path": to_posix_str(ion_path),
                                "sheet": "",
                                "sha1": "",
                                "ppm_tolerance": int(ppm_value)
                            }
                        }
                        with open(mpath, "w", encoding="utf-8") as f:
                            json.dump(mini, f, indent=2, ensure_ascii=False)
                        sample_info["method_json"] = mpath.as_posix()

                    feature_df.to_csv(out_path, index=False)

                if out_path:
                    sample_info["unlabeled_dataset"] = out_path.as_posix()
                    changed = True
                else:
                    print(f"[SKIP] sample '{sample_id}' had no excel/method/ionlist_path combo")

            except Exception as e:
                print(f"[ERROR] Failed on sample '{sample_id}': {e}")

        if changed:
            # Persist ppm to the experiment JSON for traceability
            J["prediction_parameters"] = {"ppm": int(ppm_value)}
            try:
                with open(exp_or_method, "w", encoding="utf-8") as f:
                    json.dump(J, f, indent=4)
            except Exception as e:
                messagebox.showerror("Failed to Save JSON", str(e))
                return

            messagebox.showinfo("Done", "Unlabeled datasets saved and experiment file updated.")
        else:
            messagebox.showwarning("No Samples Processed", "No eligible samples were found in this experiment file.")
    """
    """
    def create_unlabeled_dataset():
        exp_path = filedialog.askopenfilename(filetypes=[("Experiment JSON", "*.json")])
        if not exp_path:
            return
        try:
            with open(exp_path, "r", encoding="utf-8") as f:
                exp = json.load(f)
        except Exception as e:
            messagebox.showerror("Failed to Load JSON", str(e))
            return

        ppm_value = simpledialog.askinteger("PPM Tolerance", "Enter PPM tolerance (e.g., 20):", minvalue=1, maxvalue=100)
        if ppm_value is None:
            return

        for sample_id, sample_info in exp.get("samples", {}).items():
            raw_csv = sample_info.get("csv")
            annotation_excel = sample_info.get("excel")
            if not (raw_csv and annotation_excel):
                continue
            try:
                ion_list = extract_fragment_masses(annotation_excel)
                feature_df = extract_ion_intensities(raw_csv, ion_list, ppm=ppm_value)
                out_path = os.path.splitext(raw_csv)[0] + f"_unlabeled_ppm{ppm_value}.csv"
                feature_df.to_csv(out_path, index=False)
                sample_info["unlabeled_dataset"] = out_path
            except Exception as e:
                print(f"[ERROR] Failed on sample '{sample_id}':", e)

        exp["prediction_parameters"] = {"ppm": ppm_value}
        try:
            with open(exp_path, "w", encoding="utf-8") as f:
                json.dump(exp, f, indent=4)
        except Exception as e:
            messagebox.showerror("Failed to Save JSON", str(e))
            return
        messagebox.showinfo("Done", f"Unlabeled datasets saved and experiment file updated.")
    """
    # Predict tab state
    model_file_path = None
    predict_input_path = None



    train_csv_path = None
    linked_exp_json = None
    test_split_var = tk.DoubleVar(value=0.25)
    min_samples_var = tk.IntVar(value=5)

    subwin = tk.Toplevel(root)
    subwin.title("ML Analysis")
    subwin.geometry("680x560")

    notebook = ttk.Notebook(subwin)
    notebook.pack(fill="both", expand=True)

    # ... [no changes below this line: GUI layout remains as-is]


    
    # --- Tab 1: Train Model ---
    train_tab = ttk.Frame(notebook)
    notebook.add(train_tab, text="Train Model")

    tk.Label(train_tab, text="Step 1: Load Trainable Dataset (.csv)").grid(row=0, column=0, sticky="w", padx=10, pady=5)
    train_load_button = tk.Button(train_tab, text="Select CSV File", command=select_train_csv)
    train_load_button.grid(row=0, column=1, padx=5, pady=5)
    # flag: manual vs pseudo-labeled
    is_pseudolabel_var = tk.BooleanVar(value=False)
    pseudo_check = tk.Checkbutton(
        train_tab,
        text="CGA dataset?",
        variable=is_pseudolabel_var
    )
    pseudo_check.grid(row=0, column=2, padx=8, pady=5, sticky="w")
    # === Train tab: mass-as-feature option ===
    include_mass_train_var = tk.BooleanVar(value=False)

    mass_train_frame = ttk.LabelFrame(train_tab, text="Optional: include precursor mass")
    mass_train_frame.grid(row=0, column=3, columnspan=3, sticky="we", padx=10, pady=(4,8))

    ttk.Checkbutton(
        mass_train_frame,
        text="Use protonated mass as a model feature (Non-glycan = 0)",
        variable=include_mass_train_var
    ).grid(row=0, column=4, sticky="w")



    tk.Label(train_tab, text="Step 2: Select Label Column").grid(row=1, column=0, sticky="w", padx=10, pady=5)
    label_dropdown = ttk.Combobox(train_tab, values=['Structure', 'IUPACname(optional)', 'Glycanannotation2', 'GlyToucan ID'])
    label_dropdown.set("Structure")
    label_dropdown.grid(row=1, column=1, padx=5, pady=5)

    tk.Label(train_tab, text="Step 3: Choose Classifier").grid(row=2, column=0, sticky="w", padx=10, pady=5)
    classifier_var = tk.StringVar(value="rf")
    rf_button = tk.Radiobutton(train_tab, text="Random Forest (✔ functional)", variable=classifier_var, value="rf")
    xgb_button = tk.Radiobutton(train_tab, text="XGBoost (placeholder)", variable=classifier_var, value="xgb")
    svm_button = tk.Radiobutton(train_tab, text="SVM (placeholder)", variable=classifier_var, value="svm")
    knn_button = tk.Radiobutton(train_tab, text="KNN (placeholder)", variable=classifier_var, value="knn")

    rf_button.grid(row=2, column=1, sticky="w")
    xgb_button.grid(row=3, column=1, sticky="w")
    svm_button.grid(row=4, column=1, sticky="w")
    knn_button.grid(row=5, column=1, sticky="w")






   # ---- ML PARAM PANEL (paste under Step 3, before Step 4) ----
    effective_ml_params = {}  # nonlocal capture so train_model() can see it if desired

    def _get_current_context():
        # We prefer the .exp.json that was auto-detected when the user picked the CSV
        exp_json_path = linked_exp_json  # may be None until a CSV has been picked

        # Derive a friendly experiment title if we have the exp json
        exp_title = None
        try:
            if exp_json_path and os.path.exists(exp_json_path):
                with open(exp_json_path, "r", encoding="utf-8") as _f:
                    _exp = json.load(_f)
                exp_title = _exp.get("experiment")
        except Exception:
            pass

        # We don’t strictly need a method_json here; keep None (panel can still load from any file)
        method_json_path = None

        return {
            "experiment_title": exp_title,
            "exp_json_path": exp_json_path,
            "method_json_path": method_json_path,
        }

    def _on_effective_params_ready(merged):
        # Capture for use by train_model() or anywhere else in this window
        nonlocal effective_ml_params
        effective_ml_params = merged
        print("[ML PARAMS] effective parameters now in memory:", effective_ml_params)

    # Build the panel into the Train tab (it creates its own labeled frame)
    build_ml_params_panel(
        parent=train_tab,
        get_current_context=_get_current_context,
        on_effective_params_ready=_on_effective_params_ready,
    )

    # place this right AFTER the classifier radio buttons, BEFORE "Step 4" label
    #ml_frame = build_ml_params_panel(train_tab, _get_current_context, _on_effective_params_ready)
    #ml_frame.grid_configure(row=6, column=0, columnspan=3, sticky="nsew", padx=10, pady=(6, 10))

    # ---- Step 4: bar with 2 buttons + summary (Train tab) ----
    # (place this where your Step 4 label/button currently lives)
    # --- Step 4 bar (goes in open_ml_analysis_window, on the Train tab) ---
    ml_summary_var = tk.StringVar(value="Params: (using built-ins)")

    """
    # choose the next free grid row on train_tab
    row4 = train_tab.grid_size()[1]   # number of rows currently used

    
    # 3) ONLY AFTER train_tab exists and has content, add the Step-4 bar
    def _add_step4_bar():
        # next free grid row on train_tab
        next_row = train_tab.grid_size()[1]

        step4_container = ttk.Frame(train_tab)
        step4_container.grid(row=next_row, column=0, columnspan=3,
                            sticky="ew", padx=10, pady=(6, 0))

        # inside the container we can pack
        bar = tk.Frame(step4_container)
        bar.pack(fill="x")

        tk.Button(bar, text="Train/Test Parameters…",
                command=open_train_settings).pack(side="left")

        tk.Button(bar, text="Edit ML Parameters…",
                command=open_ml_params_window).pack(side="left", padx=6)

        tk.Label(bar, textvariable=ml_summary_var).pack(side="left", padx=12)

    # call it right here, after your Step 1/2/3 widgets:
    _add_step4_bar()
    """
    # a container that is gridded into train_tab
    #step4_container = ttk.Frame(train_tab)
    #step4_container.grid(row=row4, column=0, columnspan=3,
    #                    sticky="ew", padx=10, pady=(6, 0))

    # now INSIDE the container you can use pack freely
    #step4_bar = tk.Frame(step4_container)
    #step4_bar.pack(fill="x")

    ##tk.Button(step4_bar, text="Train/Test Parameters…", command=open_train_settings)\
    #.pack(side="left")

    #tk.Button(step4_bar, text="Edit ML Parameters…", command=open_ml_params_window)\
    #.pack(side="left", padx=6)



    #tk.Label(step4_bar, textvariable=ml_summary_var)\
    #.pack(side="left", padx=12)




    # button to launch the modal
    #ttk.Button(train_tab, text="Edit ML Parameters…", command=open_ml_params_window)\
    #    .grid(row=STEP4_ROW, column=0, sticky="w", padx=10, pady=(4,0))  # set STEP4_ROW to match your layout

    #ttk.Label(train_tab, textvariable=ml_summary_var)\
    #    .grid(row=STEP4_ROW, column=1, columnspan=2, sticky="w")

    # then bump your existing "Step 4: Train/Test Parameters" and below down to start at row=7 or 8
    
    tk.Label(train_tab, text="Step 4: Train/Test Parameters").grid(row=6, column=0, sticky="w", padx=10, pady=5)
    ttk.Label(train_tab, textvariable=ml_summary_var).grid(row=7, column=0, sticky="w", padx=10, pady=(0,6))
    tk.Button(train_tab, text="Set Parameters / Train the Model", command=open_train_settings).grid(row=6, column=1, padx=5, pady=5)

    train_button = tk.Button(train_tab, text="Train Model", command=train_model, bg="#CCE5FF")
    train_button.grid(row=8, column=0, columnspan=2, pady=10)

    tk.Label(train_tab, text="Trainable File Info (Origin Tracking)").grid(row=9, column=0, columnspan=2, sticky="w", padx=10, pady=(15, 5))
    origin_info = tk.Text(train_tab, height=4, width=70, state="disabled", wrap="word")
    origin_info.grid(row=10, column=0, columnspan=2, padx=10, pady=5)
    # -- Training tab and prediction tab UI (end reminder buttons) --
    tk.Label(train_tab, text="Combine Datasets for Training").grid(row=11, column=0, columnspan=2, sticky="w", padx=10, pady=(15, 5))
    tk.Button(train_tab, text="Select Datasets", command=lambda: combine_trainable_datasets_ui(root)).grid(row=11, column=1, columnspan=2, padx=10, pady=5)
    # --- Tab 2: Predict ---
    predict_tab = ttk.Frame(notebook)
    notebook.add(predict_tab, text="Predict")

    tk.Label(predict_tab, text="Step 1: Load Trained Model (.joblib/.pkl)").grid(row=0, column=0, sticky="w", padx=10, pady=5)
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
    tk.Label(predict_tab, text="(🔜) Combine Datasets for Prediction").grid(row=_next, column=0, columnspan=2, sticky="w", padx=10, pady=(15, 5))
    tk.Button(predict_tab, text="[Placeholder] Combine Datasets").grid(row=_next+1, column=0, columnspan=2, padx=10, pady=5)
    """
    # --- Tab 3 Build Trainable from Pseudolabels ---
    build_tab = ttk.Frame(notebook)
    # let the rightmost column expand so long paths aren’t cramped
    build_tab.columnconfigure(2, weight=1)

    # tiny status label above the Build button
    build_status_var = tk.StringVar(value="")
    tk.Label(build_tab, textvariable=build_status_var, fg="gray").grid(
        row=6, column=0, columnspan=3, sticky="w", padx=10, pady=(0, 4)
    )
    notebook.add(build_tab, text="From Pseudolabels")
    ng_frame = ttk.LabelFrame(build_tab, text="Include Non-glycan (NG) candidates")
    ng_frame.grid(row=3, column=0, columnspan=3, sticky="ew", padx=10, pady=(0,10))

    use_low_score_ng = tk.BooleanVar(value=True)
    use_nohit_ng     = tk.BooleanVar(value=True)
    low_score_cut    = tk.DoubleVar(value=0.07)

    ttk.Checkbutton(ng_frame, text="Include low-score bin as NG", variable=use_low_score_ng).grid(row=0, column=0, padx=8, pady=6, sticky="w")
    ttk.Label(ng_frame, text="low-score cut").grid(row=0, column=1, padx=(12,4), pady=6, sticky="e")
    tk.Spinbox(ng_frame, from_=0.0, to=0.5, increment=0.01, textvariable=low_score_cut, width=6).grid(row=0, column=2, padx=4, pady=6, sticky="w")

    ttk.Checkbutton(ng_frame, text="Include 'no-hit' spectra as NG", variable=use_nohit_ng).grid(row=0, column=3, padx=12, pady=6, sticky="w")

    # ---- Normalization / feature build options ----
    norm_frame = ttk.LabelFrame(build_tab, text="Feature normalization & ion list")
    norm_frame.grid(row=4, column=0, columnspan=3, sticky="ew", padx=10, pady=(0,10))

    rebuild_all_features = tk.BooleanVar(value=False)
    ttk.Checkbutton(norm_frame, text="Rebuild ALL features from long-form (ensure log10(I)+1 baseline=1.0)", 
                    variable=rebuild_all_features).grid(row=0, column=0, padx=8, pady=6, sticky="w")

    tk.Label(norm_frame, text="Ion list (CSV/XLSX with a 'mass' column)").grid(row=1, column=0, padx=8, pady=4, sticky="w")
    ionlist_path_var2 = tk.StringVar(value="")
    def _pick_ionlist():
        p = filedialog.askopenfilename(title="Select Ion List", 
                                    filetypes=[("CSV/XLSX", "*.csv;*.xlsx;*.xls"), ("All files", "*.*")])
        if p: ionlist_path_var2.set(p)
    ttk.Button(norm_frame, text="Choose…", command=_pick_ionlist).grid(row=1, column=1, padx=6, pady=4, sticky="w")
    ttk.Label(norm_frame, textvariable=ionlist_path_var2, wraplength=360).grid(row=1, column=2, columnspan=2, padx=6, pady=4, sticky="w")


    # State
    pseudo_path_var = tk.StringVar(value="")
    features_path_var = tk.StringVar(value="")
    out_path_var = tk.StringVar(value="(auto, next to features)")

    # Row 0: selectors
    tk.Label(build_tab, text="1) Pseudolabeled TSV/CSV (long form)").grid(row=0, column=0, sticky="w", padx=10, pady=6)
    def _sel_pseudo():
        p = filedialog.askopenfilename(title="Select pseudolabel TSV/CSV",
                                       filetypes=[("TSV/CSV", "*.tsv *.csv"), ("All files", "*.*")])
        if p: pseudo_path_var.set(os.path.abspath(p))
    tk.Button(build_tab, text="Choose…", command=_sel_pseudo).grid(row=0, column=1, padx=6, pady=6, sticky="w")
    tk.Label(build_tab, textvariable=pseudo_path_var, wraplength=360, anchor="w", justify="left").grid(row=0, column=2, sticky="w", padx=6, pady=6)

    tk.Label(build_tab, text="2) Feature/Merged CSV (wide, same sample)").grid(row=1, column=0, sticky="w", padx=10, pady=6)
    def _sel_features():
        p = filedialog.askopenfilename(title="Select feature/merged CSV (wide)",
                                       filetypes=[("CSV", "*.csv"), ("All files", "*.*")])
        if p: features_path_var.set(os.path.abspath(p))
    tk.Button(build_tab, text="Choose…", command=_sel_features).grid(row=1, column=1, padx=6, pady=6, sticky="w")
    tk.Label(build_tab, textvariable=features_path_var, wraplength=360, anchor="w", justify="left").grid(row=1, column=2, sticky="w", padx=6, pady=6)

    # Row 2: thresholds
    thr_frame = ttk.LabelFrame(build_tab, text="3) Thresholds / selection")
    thr_frame.grid(row=2, column=0, columnspan=3, sticky="ew", padx=10, pady=(4,10))
    thr_frame.columnconfigure(4, weight=1)

    tk.Label(thr_frame, text="min ion score").grid(row=0, column=0, padx=8, pady=6, sticky="w")
    ion_score_min = tk.DoubleVar(value=0.07)
    tk.Spinbox(thr_frame, from_=0.00, to=1.00, increment=0.01, textvariable=ion_score_min, width=6).grid(row=0, column=1, padx=6, pady=6)

    tk.Label(thr_frame, text="max |ppm_error|").grid(row=0, column=2, padx=12, pady=6, sticky="w")
    ppm_abs_max = tk.DoubleVar(value=20.0)
    tk.Spinbox(thr_frame, from_=0.0, to=200.0, increment=1.0, textvariable=ppm_abs_max, width=6).grid(row=0, column=3, padx=6, pady=6)

    tk.Label(thr_frame, text="Top N per scan (by ion score)").grid(row=0, column=4, padx=12, pady=6, sticky="w")
    topn_var = tk.IntVar(value=1)
    tk.Spinbox(thr_frame, from_=1, to=5, textvariable=topn_var, width=4).grid(row=0, column=5, padx=6, pady=6)

    # Row 3: label mapping
    map_frame = ttk.LabelFrame(build_tab, text="4) Map composition → label column")
    map_frame.grid(row=3, column=0, columnspan=3, sticky="ew", padx=10, pady=(0,10))
    tk.Label(map_frame, text="Target label column").grid(row=0, column=0, padx=8, pady=6, sticky="w")
    label_target = ttk.Combobox(map_frame,
        values=['Structure', 'IUPACname(optional)', 'Glycanannotation2', 'GlyToucan ID'])
    label_target.set('Structure')
    label_target.grid(row=0, column=1, padx=6, pady=6, sticky="w")

    fill_unlabeled_as_none = tk.BooleanVar(value=False)
    ttk.Checkbutton(map_frame, text="Fill non-matched scans as 'None' (negative class)",
                    variable=fill_unlabeled_as_none).grid(row=0, column=2, padx=12, pady=6, sticky="w")
    
    # Output line → row 5 (moved down so it doesn’t overlap)
    tk.Label(build_tab, text="5) Output (auto-named):").grid(row=5, column=0, sticky="w", padx=10, pady=6)
    tk.Label(build_tab, textvariable=out_path_var, anchor="w", justify="left").grid(row=5, column=1, columnspan=2, sticky="w", padx=6, pady=6)


    # ---- ML Parameters panel (optional but handy to keep with Prepare Dataset) ----
    ml_box = ttk.LabelFrame(subwin, text="ML Parameters")
    ml_box.pack(fill="both", expand=False, padx=10, pady=(0, 10))

    # This function was included in the helpers you pasted earlier:
    # build_ml_params_panel(parent, get_current_context, on_effective_params_ready)
    #build_ml_params_panel(ml_box, _get_ml_context, _on_params_ready)

    # Helpers


    def _choose_scan_col(df):
        for c in ["MS2scan_no","unique_ID","ScanNum","scan","Scan"]:
            if c in df.columns: return c
        return None

    def _save_and_autoload(df_out, base_features_path):
        base, ext = os.path.splitext(base_features_path)
        out_path = base + "_PLabeled.csv"
        df_out.to_csv(out_path, index=False)
        out_path_var.set(out_path)

        # auto-load into Train tab
        nonlocal train_csv_path, linked_exp_json
        train_csv_path = os.path.abspath(out_path)
        origin_info.configure(state="normal")
        origin_info.delete(1.0, "end")
        origin_info.insert("end", f"Loaded file: {os.path.basename(out_path)}\n")
        origin_info.insert("end", f"Path: {out_path}")
        origin_info.configure(state="disabled")

        messagebox.showinfo("Done", f"Trainable CSV saved:\n{out_path}\n\nYou can now click “Train Model” in the first tab.")

    # 20250910 convert tuple to canonical label FxHxNxSxGxKDNx
    import re, ast
    def _tuple_to_FHNSGKDN(t):
        # comp tuple order is (Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc)
        h, n, s, g, kdn, f = map(int, t)
        parts = []
        if f:   parts.append(f"F{f}")
        if h:   parts.append(f"H{h}")
        if n:   parts.append(f"N{n}")
        if s:   parts.append(f"S{s}")
        if g:   parts.append(f"G{g}")
        if kdn: parts.append(f"KDN{kdn}")
        return "".join(parts) or "Non-glycan"

    _TUPLE_LIKE_RE = re.compile(r"^\(\s*\d+\s*,\s*\d+\s*,\s*\d+\s*,\s*\d+\s*,\s*\d+\s*,\s*\d+\s*\)$")
    _FHNSKDN_RE    = re.compile(r"^(F\d+)?(H\d+)?(N\d+)?(S\d+)?(G\d+)?(KDN\d+)?$")

    def _coerce_comp_to_tuple(x):
        # return a (H,N,S,G,KDN,F) tuple or None
        if isinstance(x, (list, tuple)) and len(x) == 6:
            return tuple(int(v) for v in x)
        s = str(x).strip()
        if _TUPLE_LIKE_RE.match(s):
            try:
                t = ast.literal_eval(s)
                if isinstance(t, (list, tuple)) and len(t) == 6:
                    return tuple(int(v) for v in t)
            except Exception:
                return None
        return None
    # --- end helpers ---


    # 20250911 added for solving ion list required issue when rebuild feature is not needed 
    import math

    def _load_ion_masses_from_file(path):
        if not path:
            return None
        try:
            # supports .csv or .xlsx (pandas can read both if engine installed)
            df = _robust_read_csv(path)
            # common column names: 'mass', 'mz'
            for col in ("mass", "Mass", "mz", "m/z", "MZ"):
                if col in df.columns:
                    vals = []
                    for v in df[col].tolist():
                        try:
                            vals.append(float(v))
                        except Exception:
                            pass
                    return sorted(set(vals))
        except Exception:
            pass
        return None

    _EXCLUDE_NONMASS = {
        "entry_no","MS1scan_no","MS1_isolationmass","MS1_monoisolationmass",
        "chargeState","protonatedmass","MS2scan_no","label","Structure",
        "composition","theoretical_mass","observed_mass","ppm_error",
        # add any other metadata headers you know appear in your feature CSV
    }

    def _infer_ion_masses_from_feature_df(feat_df):
        masses = []
        for c in feat_df.columns:
            if c in _EXCLUDE_NONMASS:
                continue
            try:
                masses.append(float(c))
            except Exception:
                # ignore headers that aren't pure numeric (e.g., 'HCD_energy')
                continue
        return sorted(set(masses))

    """
    #remove funtion to avoid GPT getting distracted

    close_button = tk.Button(subwin, text="Close", command=subwin.destroy)
    close_button.pack(pady=5)












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


#icon
ico_path = os.path.join(os.path.dirname(__file__), 'GlycoMSPlogo.ico')
png_path = os.path.join(os.path.dirname(__file__), 'GlycoMSPlogo.png')




# Initialize GUI
root = tk.Tk()
#root.iconbitmap(default=ico_path)
if platform.system() == "Windows" and os.path.exists(ico_path):
    root.iconbitmap(default=ico_path)
elif platform.system() in ("Darwin", "Linux") and os.path.exists(png_path):
    icon_img = tk.PhotoImage(file=png_path)
    root.iconphoto(True, icon_img)
    
root.protocol("WM_DELETE_WINDOW", on_closing)
root.title("GlycoMSP File Manager GUI v0.6 Build 20251010 core 0.9998")
root.geometry("840x600")
root.minsize(840, 600)



# Buttons
button_frame = tk.Frame(root)
button_frame.pack(pady=10)

tk.Button(button_frame, text="Select Raw File", command=lambda: select_file("raw")).grid(row=0, column=0, padx=5)
#tk.Button(button_frame, text="Select mzML File", command=lambda: select_file("mzml")).grid(row=0, column=1, padx=5)
tk.Button(button_frame, text="Select CSV File", command=lambda: select_file("csv")).grid(row=0, column=2, padx=5)
tk.Button(button_frame, text="Select Excel File", command=lambda: select_file("excel")).grid(row=0, column=3, padx=5)
tk.Button(button_frame, text="Clear All", command=clear_files).grid(row=0, column=4, padx=5)
tk.Button(root, text="Save Log", command=save_log_to_file).pack(pady=5)
tk.Button(root, text="About", command=open_about_window).pack(pady=5)
# Text widget to log selected files
text_widget = tk.Text(root, height=15, width=80)
text_widget.pack(pady=10)


#status bar?
status_var = tk.StringVar()
status_var.set("Idle")
status_label = tk.Label(root, textvariable=status_var, fg="blue")
status_label.pack(pady=5)
#progress bar?
progress = ttk.Progressbar(root, orient="horizontal", mode="indeterminate", length=250)
progress.pack(pady=5)

# NEW: helper to update main window status safely
def set_main_status(text, fg=None):
    status_var.set(text)
    if fg is not None:
        try:
            status_label.config(fg=fg)
        except Exception:
            pass


# --- Add Analysis Tools Frame ---

analysis_frame = tk.LabelFrame(root, text="Analysis Tools", padx=10, pady=10)
analysis_frame.pack(padx=10, pady=10, fill="x")

convert_button = tk.Button(analysis_frame, text="Convert Raw to CSV",command=launch_metadata_batch)
convert_button.pack(side="left", padx=5)
tk.Button(analysis_frame, text="Prepare Dataset", command=open_prepare_dataset_window).pack(side="left", padx=5)
tk.Button(analysis_frame, text="Run ML Analysis", command=open_ml_analysis_window).pack(side="left", padx=5)


# Run the GUI
#root.mainloop()

#import safe
if __name__ == "__main__":
    root.mainloop()

#20250915 import safe and no side effects by loading too many modules, consider activate it in future
"""
# --- at top of mspfileloaderv10.py ---
root = None  # set up a module-global you can reference

def build_gui():
    global root
    import tkinter as tk
    from tkinter import ttk
    root = tk.Tk()
    root.title("GlycoMSParser v10")

    # ... all your widget creation & menu wiring here ...
    # e.g. TreeView, buttons, callbacks, etc.

    return root

# keep helpers importable for CLI/REPL
__all__ = [
    # your non-GUI helpers:
    "read_fragment_masses_any",
    "create_unlabeled_from_method",
    # (add others you want to script against)
]

if __name__ == "__main__":
    app = build_gui()
    app.mainloop()
"""