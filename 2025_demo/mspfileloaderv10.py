import os
version = "0.9963"
last_update = 20250910
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


#20250907 preventing key error in GUI #need real test to see behavior changes
FILETYPE_TO_KEY = {"csv": "csv", "excel": "excel", "json": "json",
                   "method": "json", "metadata": "json"}  # extend if you add new labels
def normalize_ftype(ft: str) -> str:
    return FILETYPE_TO_KEY.get(ft.lower().strip(), ft.lower().strip())
#

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
        self.title("Pseudo-Labeling Setup")
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
            "Hex_range": ("range", (0, 20)),
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
            "arm_count","internal_minrep","internal_maxrep","topology",
            "corefuc","bicorefuc","highman","perman","hybrid",
            "compcheck","Hex_range","HexNAc_range","Neu5Ac_range","Neu5Gc_range","KDN_range","Fucose_range"
        }
        OG_KEYS = {
            "debug","dev","force_exit",
            "allow5ac","allow5gc","allowkdn","allowfuc","allowpsa",
            "arm_count","internal_minrep","internal_maxrep","topology",
            "compcheck","Hex_range","HexNAc_range","Neu5Ac_range","Neu5Gc_range","KDN_range","Fucose_range"
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

        def _current_keys():
            gtype = self.meta_vars["Glycan Type"].get().strip().upper()
            return NG_KEYS if gtype == "N" else OG_KEYS

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

        # --- O-glycan core types (multi-select) ---
        core_frame = ttk.LabelFrame(self, text="O-glycan core types (select 1–4)")
        core_frame.pack(fill="x", padx=12, pady=(0, 6))

        self.og_core_vars = {i: tk.BooleanVar(value=False) for i in (1, 2, 3, 4)}

        row = 0
        for i, label in [(1, "Core 1"), (2, "Core 2"), (3, "Core 3"), (4, "Core 4")]:
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
                core_frame.pack(fill="x", padx=12, pady=(0, 6))

        # After creating the Glycan Type combobox (named via self.meta_vars["Glycan Type"])
        gly_cb = meta_frame.grid_slaves(row=0, column=1)[0]  # the Combobox you just created
        gly_cb.bind("<<ComboboxSelected>>", _rebuild_flag_panel)

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
        links_frame = ttk.LabelFrame(self, text="Linked files (optional)")
        links_frame.pack(fill="x", padx=12, pady=(6, 0))

        self.insilico_path_var = tk.StringVar(value=initial_insilico or "")
        self.ionlist_path_var  = tk.StringVar(value=initial_ionlist or "")

        ttk.Label(links_frame, text="In-silico CSV:").grid(row=0, column=0, sticky="w", padx=10, pady=4)
        ttk.Label(links_frame, textvariable=self.insilico_path_var).grid(row=0, column=1, sticky="w", padx=8, pady=4)

        ttk.Label(links_frame, text="Ion list (optional):").grid(row=1, column=0, sticky="w", padx=10, pady=4)
        ttk.Label(links_frame, textvariable=self.ionlist_path_var).grid(row=1, column=1, sticky="w", padx=8, pady=4)

        # --- Footer buttons
        btns = ttk.Frame(self)
        btns.pack(fill="x", padx=12, pady=(6, 12))

        ttk.Button(btns, text="Load Flags…", command=self.load_flags).pack(side="left", padx=4)
        ttk.Button(btns, text="Save Flags…", command=self.save_flags).pack(side="left", padx=4)

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

        ttk.Button(btns, text="Link Existing…", command=_link_existing).pack(side="left", padx=4)

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

        ttk.Button(btns, text="Attach Ion List…", command=_attach_ionlist).pack(side="left", padx=4)

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

        ttk.Button(btns, text="Start Pseudolabeling", command=_start).pack(side="left", padx=4)

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
        fields = [
            "Experiment Title", "Experiment Description", "Author running this analysis",
            "Raw data acquired date", "Glycan Type", "Mass Analyzer charge mode", "Derivatization Type"
        ]

        self.window = tk.Toplevel(parent)
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


    def on_conversion_complete(self):
        self.status_label.config(text="Conversion complete. Ready for metadata.", fg="green")

        # Extract filename base for use in renaming
        self.rawfilename = os.path.splitext(os.path.basename(self.current_raw_file))[0]

        # Re-enable entries and buttons
        for entry in self.entries.values():
            entry.config(state="normal")
        for child in self.window.winfo_children():
            if isinstance(child, tk.Button):
                child.config(state="normal")

        # Pre-fill if last used
        if self.last_metadata:
            self.fill_fields_from_metadata(self.last_metadata)

        self.window.title(f"Metadata for: {os.path.basename(self.current_raw_file)}")

    def load_file(self, raw_file):
        self.current_raw_file = raw_file
        self.window.title(f"Converting: {os.path.basename(raw_file)}")
        self.clear_fields()

        # Add status label if not already there
        if not hasattr(self, "status_label"):
            self.status_label = tk.Label(self.window, text="", fg="blue")
            self.status_label.grid(row=0, columnspan=2, pady=10)
        self.status_label.config(text="Converting raw file... please wait.", fg="blue")

        # Disable form and buttons temporarily
        for entry in self.entries.values():
            entry.config(state="disabled")
        for child in self.window.winfo_children():
            if isinstance(child, tk.Button):
                child.config(state="disabled")

        def background_conversion():
            if not self.skip_conversion:
                try:
                    success = mspext.convert_raw_to_csv(raw_file, debug=True)
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
                    self.window.after(0, lambda: messagebox.showerror("Thermo Library Error",f"Raw file could not be processed.\n\nDetails:\n{e}"))
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

        # Auto-close the window
        self.window.destroy()

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

def _read_ion_df(path):
    import pandas as pd
    if not path:
        return None
    p = str(path).lower()

    if p.endswith((".xlsx", ".xls")):
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
    else:
        df = pd.read_csv(path, engine="python")

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
    mspext.convert_raw_to_csv(raw_file, debug=False)

def launch_metadata_for_all(files):
    MetadataEditorWindow(root, files, on_metadata_ready, on_finish=reset_main_status)





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

        # 2) Normalize in-silico library → comp_tuple/comp_str + sorted Mass
        libn = marker.normalize_insilico(
            lib,
            comp_cols=("Hex","HexNAc","NeuAc","NeuGc","KDN","Fuc"),
            mass_col="Mass",
            add_legacy_repr=True,
        )
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
        matched = matched.rename(columns={"comp_str": "composition"})
        import math

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
        """
        peaks_cols = [c for c in ("peaklist","peakintensity") if c in df.columns]
        if len(peaks_cols) == 2:
            if scan_col == "MS2scan_no":
                # same key name on both sides → use 'on=' to avoid suffixes
                matched = matched.merge(
                    df[[scan_col] + peaks_cols],
                    on="MS2scan_no", how="left"
                )
            else:
                matched = matched.merge(
                    df[[scan_col] + peaks_cols],
                    left_on="MS2scan_no", right_on=scan_col, how="left"
                ).drop(columns=[scan_col], errors="ignore")
        """
        #if len(peaks_cols) == 2:
        #    matched = matched.merge(df[[scan_col] + peaks_cols], left_on="MS2scan_no", right_on=scan_col, how="left") \
        #                    .drop(columns=[scan_col])  


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

        """
        # 5) Ion scoring (optional; compact columns appended)
        ion_scoring_status, n_with_scores = "skipped", 0
        ion_df = _read_ion_df(ion_path) if ion_path else None

        #temp debug
        print("[ion] path:", ion_path,
        "ion_df_rows:", 0 if (ion_df is None) else len(ion_df),
        "has_peaks:", bool("peaklist" in matched.columns and "peakintensity" in matched.columns))
        if ion_df is not None and not ion_df.empty and all(c in matched.columns for c in ("peaklist","peakintensity")):
            try:
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
                    if src in matched.columns and dst not in matched.columns:
                        matched.rename(columns={src: dst}, inplace=True)
                #temp change to a later version
                
                matched = marker.attach_ion_score_on_matched(
                    matched_df=matched,
                    ion_df=ion_df,
                    ppm_value=ion_ppm,
                    scan_col="MS2scan_no",
                    ion_mass_col="mass",
                    score_col="ion score",
                    hitcount_col="ion hit count",
                    hitlist_col="ion hits m/z",
                )
                ion_scoring_status = "ok"
                n_with_scores = int((matched["ion hit count"] > 0).sum())
                
            except Exception:
                traceback.print_exc()
                ion_scoring_status = "failed"
            """

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
        # 7) Save TSV next to converted CSV
        outdir  = os.path.dirname(csv_path)
        outname = f"{sample_name}_pseudolabels_{datetime.now().strftime('%Y%m%d')}.tsv"
        outpath = os.path.join(outdir, outname)
        out.to_csv(outpath, index=False, sep="\t")

        files["pseudolabel_csv"] = outpath
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

        if parent:
            messagebox.showinfo("Pseudolabeling complete", f"Saved and linked:\n{outpath}")
        try:
            refresh_tree()
        except Exception:
            pass

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



        def on_submit(payload):
            # 1) Start from NG_flags defaults, overlay UI edits
            flags = dict(compv4.NG_flags)
            flags.update(payload.get("flags", {}))

            # 2) Determine glycan type from metadata (prefer window metadata; fallback to file)
            glycan_type = (payload["metadata"].get("Glycan Type") or "").strip().upper()
            charge_mode = payload["metadata"].get("Mass Analyzer charge mode")  # if needed downstream
            # Example normalization at the GUI edge (optional but nice):
            deriv_raw = (payload["metadata"].get("Derivatization Type") or "").strip().lower()
            reduced = False
            if "reduced" in deriv_raw:   # e.g., "Reduced Permethylation"
                reduced = True
            flags["derivatization_type"] = payload["metadata"].get("Derivatization Type", "PerMe")
            flags["reduced"] = reduced
            #deriv_type  = payload["metadata"].get("Derivatization Type")        # if needed downstream
            if payload["metadata"].get("_overrides_applied"):
                logger.log("[WARN] Using metadata overrides from UI; consider updating the metadata JSON.")
            # 3) Choose output path
            outdir = filedialog.askdirectory(title="Select output folder for in-silico CSV")
            if not outdir:
                return
            outname = f"{sample_name}_insilico_{datetime.now().strftime('%Y%m%d')}.csv"
            outpath = os.path.join(outdir, outname)

            try:
                # 4) Launch NG or OG generator (no hybrid component dicts yet)
                if glycan_type == "N":
                    compv4.NGlaunch(user_flags=flags, filename=outpath, debug=flags.get("debug", False))
                elif glycan_type == "O":
                    # OGbranch: we can extend with coretype/keep_topology later
                    compv4.OGlaunch(user_flags=flags, filename=outpath, debug=flags.get("debug", False))
                else:
                    messagebox.showwarning("Glycan Type Missing",
                                        "Metadata lacks a valid 'Glycan Type' (expected 'N' or 'O').")
                    return

                files["insilico_csv"] = outpath
                messagebox.showinfo("In-silico CSV generated", f"Saved and linked:\n{outpath}")
                refresh_tree()
                # (Optional) attach outpath to the sample record here so the next step can find it
                # files["insilico_csv"] = outpath
                refresh_tree()

            except Exception as e:
                traceback.print_exc()
                messagebox.showerror("Generation failed", str(e))

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

            outdir = filedialog.askdirectory(title="Select output folder for in-silico CSV")
            if not outdir:
                return None
            outname = f"{sample_name}_insilico_{datetime.now().strftime('%Y%m%d')}.csv"
            outpath = os.path.join(outdir, outname)

            try:
                if glycan_type == "N":
                    compv4.NGlaunch(user_flags=flags, filename=outpath, debug=flags.get("debug", False))
                elif glycan_type == "O":
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



    # --- Button panel ---
    button_frame = tk.Frame(subwin)
    button_frame.pack(pady=5)

    tk.Button(button_frame, text="Select MS2 CSV(s)", command=lambda: select_files_generic("csv", True, handle_csv_selection)).grid(row=0, column=0, padx=5)
    tk.Button(button_frame, text="Select Excel", command=lambda: select_files_generic("excel", True, handle_excel_selection)).grid(row=0, column=1, padx=5)
    tk.Button(button_frame, text="Add sample (metadata file required)", command=lambda: select_files_generic("json", True, handle_json_selection)).grid(row=1, column=0, padx=5)
    tk.Button(button_frame, text="Add Sample (missing metadata)", command=add_sample).grid(row=1, column=1, padx=5)
    tk.Button(button_frame, text="Clean up empty unassigned sample tags", command=clean_unassigned_samples).grid(row=1, column=2, padx=5)
    link_button = tk.Button(button_frame, text="Link Sample", state="disabled", command=lambda: try_link_selected_sample())
    link_button.grid(row=2, column=0, padx=5) #why it was gone?
    merge_button = tk.Button(button_frame, text="Merge Sample", state="disabled", command=lambda: try_merge_selected_sample())
    merge_button.grid(row=2, column=1, padx=5)
    tk.Button(button_frame, text="Negative options…",
          command=open_negative_options_dialog).grid(row=2, column=2, padx=5)
    tk.Button(button_frame, text="Ion suggestions…",
          command=open_ion_suggest_dialog).grid(row=2, column=3, padx=5)
    tk.Button(button_frame, text="View suggestions…",
          command=open_ion_suggestions_viewer).grid(row=2, column=4, padx=5)
    tk.Button(button_frame, text="Load Method", command=load_method_file).grid(row=3, column=1, padx=5)
    ttk.Button(button_frame, text="Assign Pseudo-Labels by Glycan Composition", command=lambda:launch_pseudo_labeling()).grid(row=3, column=0, padx=5, pady=5) 
    #GPT said without () it only passes the function, and work only if clicked
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


    #v0.9923~0.9929 utilities
    from ml_ng_utils import (
        collect_ng_candidates,
        cap_non_glycan,
        predict_with_threshold,
        build_features_from_peaks_log10_plus1,  # optional if you need it directly
    )
    from ml_ng_utils_extras import resample_by_strategy, tau_sweep_summary

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

    """
    # ---------- NEW: balancing helper ---------- 20250901
    def balance_and_split(
        df: pd.DataFrame,
        label_col="Structure",
        majority_label="Non-glycan",
        min_count=5,                  # tiny-class filter
        majority_factor=3,            # cap majority to factor * max(minor count)
        feature_exclude=("MS2scan_no","ID","Source","IUPACname(optional)","Glycanannotation2","GlyToucan ID","unique_ID"),
        test_size=0.2,
        val_size=0.1,
        stratify=True,
        random_state=42
    ):
        from sklearn.model_selection import train_test_split

        # A) drop tiny classes
        counts = df[label_col].value_counts()
        keep = counts[counts >= min_count].index
        dropped_rare = df[~df[label_col].isin(keep)]
        df1 = df[df[label_col].isin(keep)].copy()

        # B) cap the majority class to majority_factor * max(minor)
        cap_applied = None
        if majority_label in df1[label_col].unique():
            minor_counts = df1[df1[label_col] != majority_label][label_col].value_counts()
            max_minor = int(minor_counts.max()) if not minor_counts.empty else 0
            cap = max(majority_factor * max_minor, 1)
            majority_df = df1[df1[label_col] == majority_label]
            if len(majority_df) > cap:
                majority_df = majority_df.sample(n=cap, random_state=random_state)
                df1 = pd.concat([majority_df, df1[df1[label_col] != majority_label]], axis=0)
                cap_applied = cap

        # inside balance_and_split(...) right before the split section
        # added 20250902, to solve the validation set to 0 error occurred in 20250901
        feat_cols = [c for c in df1.columns if c not in set(feature_exclude) | {label_col}]
        X, y = df1[feat_cols], df1[label_col]
        if y.nunique() < 2:
            raise ValueError("After filtering, fewer than 2 classes remain. Loosen filters.")

        # --- NEW: robust split math ---
        ts = float(test_size)
        vs = float(val_size)
        if ts < 0 or vs < 0:
            raise ValueError("Test/validation splits must be >= 0.")
        if ts == 0 and vs == 0:
            raise ValueError("At least one of test or validation must be > 0.")
        holdout = ts + vs
        if holdout >= 0.999:
            raise ValueError(f"Test + Validation ({holdout:.2f}) must be < 1.0.")

        strat = y if stratify else None
        from sklearn.model_selection import train_test_split
        X_train, X_tmp, y_train, y_tmp = train_test_split(
            X, y, test_size=holdout, stratify=strat, random_state=random_state
        )

        if vs > 0:
            # second split between val and test
            rel_test = ts / holdout
            # clamp away from 0 and 1 for sklearn
            rel_test = min(max(rel_test, 1e-6), 1 - 1e-6)
            strat_tmp = y_tmp if stratify else None
            X_val, X_test, y_val, y_test = train_test_split(
                X_tmp, y_tmp, test_size=rel_test, stratify=strat_tmp, random_state=random_state
            )
        else:
            # no validation set requested
            X_val, y_val = X_tmp.iloc[0:0], y_tmp.iloc[0:0]  # empty
            X_test, y_test = X_tmp, y_tmp
        # --- END NEW ---
        
        # C) split (optionally stratified) →  train/val/test
        feat_cols = [c for c in df1.columns if c not in set(feature_exclude) | {label_col}]
        X, y = df1[feat_cols], df1[label_col]
        if y.nunique() < 2:
            raise ValueError("After filtering, fewer than 2 classes remain. Loosen filters.")

        strat = y if stratify else None
        X_train, X_tmp, y_train, y_tmp = train_test_split(
            X, y, test_size=(test_size + val_size), random_state=random_state, stratify=strat
        )
        rel_test = test_size / (test_size + val_size) if (test_size + val_size) > 0 else 0.0
        strat_tmp = y_tmp if stratify else None
        X_val, X_test, y_val, y_test = train_test_split(
            X_tmp, y_tmp, test_size=rel_test, random_state=random_state, stratify=strat_tmp
        )
        
        info = {
            "kept_label_counts": y.value_counts().to_dict(),
            "dropped_rare_counts": dropped_rare[label_col].value_counts().to_dict(),
            "majority_cap_applied_to": majority_label if cap_applied else None,
            "majority_cap": int(cap_applied) if cap_applied else None,
            "feature_cols": feat_cols,
        }
        return X_train, y_train, X_val, y_val, X_test, y_test, info
        """
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

    # ---------- UPDATED: parameters window ----------
    def open_train_settings():
        settings = tk.Toplevel()
        settings.title("Train/Test Parameters")
        settings.geometry("360x360")

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

        tk.Button(settings, text="Close", command=settings.destroy).pack(pady=12)
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

        #convert the composition to man-readable one
        try:
            df[label_col] = df[label_col].apply(normalizer.parse_structure_to_manual)
        except Exception as e:
            print(f"[dev]Encounter structure conversion error")
            messagebox.showerror("Can't convert the structure(composition) to man-like format", str(e)); return

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

        # LabelEncoder (consistent across splits)
        le = LabelEncoder()
        y_train_enc = le.fit_transform(y_train)
        y_val_enc   = le.transform(y_val)
        y_test_enc  = le.transform(y_test)

        # model
        model = RandomForestClassifier(
            n_estimators=int(n_estimators_var.get()),
            class_weight=("balanced" if class_weight_var.get() else None),
            max_features="sqrt",
            random_state=42,
            n_jobs=-1
        )
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


    """
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
                try:
                    exp_data["train_csv_is_pseudolabel"] = bool(is_pseudolabel_var.get())    
                except:
                    print("[dev]No pseudolabel info in metadata (no impact)")
                exp_data["train_csv"] = train_csv_path
                with open(exp_json_path, "w") as f:
                    json.dump(exp_data, f, indent=4)
            except Exception as e:
                print(f"[ML] Failed to update experiment JSON: {e}")

        origin_info.configure(state="normal")
        origin_info.delete(1.0, "end")
        origin_info.insert("end", f"\nDataset type: {'Pseudo-labeled' if is_pseudolabel_var.get() else 'Manual'}")
        origin_info.insert("end", f"Loaded file: {os.path.basename(path)}\n")
        if exp_json_path:
            origin_info.insert("end", f"Linked .exp.json: {os.path.basename(exp_json_path)}\n")
        origin_info.insert("end", f"Path: {path}")
        origin_info.configure(state="disabled")

    def open_train_settings():
        settings = tk.Toplevel()
        settings.title("Train/Test Parameters")
        settings.geometry("300x200")

        tk.Label(settings, text="Test Split Ratio (e.g., 0.25 = 75/25)").pack(pady=(10, 0))
        split_slider = tk.Scale(settings, from_=0.1, to=0.5, resolution=0.05, orient="horizontal", variable=test_split_var)
        split_slider.pack(pady=5)

        tk.Label(settings, text="Minimum Samples per Class").pack(pady=(15, 0))
        min_sample_spin = tk.Spinbox(settings, from_=1, to=50, textvariable=min_samples_var, width=5)
        min_sample_spin.pack(pady=5)

        tk.Button(settings, text="Close", command=settings.destroy).pack(pady=10)

    def train_model():
        try:
            from sklearn.ensemble import RandomForestClassifier
            from sklearn.model_selection import train_test_split
            from sklearn.metrics import classification_report
            from sklearn.preprocessing import LabelEncoder
            import joblib
            from sklearn.utils.multiclass import unique_labels

        except ImportError as e:
            messagebox.showerror("Missing Dependencies", "Required package not found.\n\nPlease make sure scikit-learn and joblib are installed.")
            return

        if not train_csv_path:
            messagebox.showwarning("No File", "Please select a training CSV file first.")
            return

        label_col = label_dropdown.get().strip()
        if label_col not in ['Structure', 'IUPACname(optional)', 'Glycanannotation2', 'GlyToucan ID']:
            messagebox.showerror("Invalid Label", f"'{label_col}' is not a supported label column.")
            return

        try:
            df = pd.read_csv(train_csv_path)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read CSV file:\n{e}")
            return

        if label_col not in df.columns:
            messagebox.showerror("Missing Column", f"'{label_col}' not found in the dataset.")
            return

        min_samples = min_samples_var.get()
        class_counts = df[label_col].value_counts()
        valid_classes = class_counts[class_counts >= min_samples].index
        df = df[df[label_col].isin(valid_classes)].copy()

        if is_pseudolabel_var.get():
            print("[ML] Training with pseudo-labeled dataset")
        else:
            print("[ML] Training with manual dataset")

        drop_cols = ['ID', 'Source', 'IUPACname(optional)', 'Glycanannotation2', 'GlyToucan ID', 'unique_ID']
        X = df.drop(columns=[col for col in drop_cols if col in df.columns] + [label_col], errors='ignore')
        y = df[label_col]

        if y.dtype == 'object':
            le = LabelEncoder()
            y = le.fit_transform(y)
            # Save encoder for later decoding
            le_path = os.path.splitext(train_csv_path)[0] + "_labelencoder.joblib"
            joblib.dump(le, le_path)
        try:
            test_ratio = test_split_var.get()
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_ratio, random_state=42)
            model = RandomForestClassifier(n_estimators=100, random_state=42)
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            #class_names = le.classes_ #so the encoded classes are restored
            #report = classification_report(y_test, y_pred, target_names=class_names)
            used_labels = unique_labels(y_test, y_pred)
            used_class_names = [le.classes_[i] for i in used_labels]
            report = classification_report(y_test, y_pred, target_names=used_class_names)
            #report = classification_report(y_test, y_pred)

            messagebox.showinfo("Training Complete", f"Random Forest trained successfully.\n\n{report}")

            base = os.path.splitext(train_csv_path)[0]
            model_path = base + "_rf_model.joblib"
            report_path = base + "_rf_performance.txt"
            joblib.dump(model, model_path)
            with open(report_path, "w") as f:
                f.write(report)

            if linked_exp_json and os.path.exists(linked_exp_json):
                try:
                    with open(linked_exp_json, "r") as f:
                        exp_data = json.load(f)
                    exp_data["train_parameters"] = {
                        "split_ratio": test_ratio,
                        "min_samples": min_samples,
                        "model_path": model_path,
                        "report_path": report_path,
                        "is_pseudolabeled": bool(is_pseudolabel_var.get()),
                    }
                    with open(linked_exp_json, "w") as f:
                        json.dump(exp_data, f, indent=4)
                except Exception as e:
                    print("Failed to write training parameters to exp.json:", e)

        except Exception as e:
            messagebox.showerror("Training Failed", str(e))
        """
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

        try:
            # Build X, drop traceability cols
            X = df.drop(columns=["MS2scan_no", "protonatedmass"], errors="ignore")

            # Align features to training set
            train_feats = get_training_features(model, model_file_path)
            if train_feats is None:
                raise RuntimeError("Cannot determine training feature list. "
                                "Re-export the model with embedded feature names or provide *_features.json.")

            missing = sorted(set(train_feats) - set(X.columns))
            extra   = sorted(set(X.columns) - set(train_feats))
            for c in missing: X[c] = 0.0
            X = X[list(train_feats)]  # enforce order

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
            # 3) Build a compact DataFrame for the reporter
            pred_df = df.copy()
            pred_df = pred_df.rename(columns={"Predicted_Label": "pred_label"})  # reporter expects 'pred_label'

            # Prefer a single vector column to keep CSV lean
            # adjust y_proba (duplicated) to proba that holding same meanings
            proba_cols = None
            if proba is not None and isinstance(proba, (list, tuple)) or hasattr(proba, "shape"):
                try:
                    import numpy as _np
                    pred_df["proba_vector"] = [ _np.asarray(row, dtype=float) for row in proba ]
                except Exception:
                    # fallback: expand into columns if needed
                    if class_names:
                        proba_cols = [f"proba_{c}" for c in class_names]
                    else:
                        proba_cols = [f"proba_{i}" for i in range(proba.shape[1])]
                    for j, col in enumerate(proba_cols):
                        pred_df[col] = proba[:, j]
            """
            if y_proba is not None and isinstance(y_proba, (list, tuple)) or hasattr(y_proba, "shape"):
                try:
                    import numpy as _np
                    pred_df["proba_vector"] = [ _np.asarray(row, dtype=float) for row in y_proba ]
                except Exception:
                    # fallback: expand into columns if needed
                    if class_names:
                        proba_cols = [f"proba_{c}" for c in class_names]
                    else:
                        proba_cols = [f"proba_{i}" for i in range(y_proba.shape[1])]
                    for j, col in enumerate(proba_cols):
                        pred_df[col] = y_proba[:, j]
            """
            # 4) Parameters (use the same thresholds you apply in your pipeline)
            params = ReportParams(tau=0.60, margin=0.05, topk=5, sample_cols=("experiment_title", "sample_name"))

            # 5) Choose output folder (same folder as input unlabeled CSV)
            from pathlib import Path
            out_dir = Path(os.path.dirname(predict_input_path))

            # 6) Run summarization + write artifacts
            pred_rows, class_sum, by_sample_sum, run_sum = summarize_predictions(
                pred_df,
                class_names=class_names if class_names else None,
                proba_cols=proba_cols,  # None if using 'proba_vector'
                params=params,
            )

            arts = write_prediction_report(
                out_dir=out_dir,
                pred_rows=pred_rows,
                class_summary=class_sum,
                by_sample=by_sample_sum,
                run_summary=run_sum,
                params=params,
            )

            # 7) Show popup (we're on the Tk main thread here; if you later thread this, wrap with root.after)
            try:
                show_prediction_summary_popup(root, run_summary=run_sum, artifacts=arts, class_summary_df=class_sum)
            except Exception as _e:
                print("[warn] Failed to show summary popup:", _e)

            # --- END: Prediction summary integration ---

            out_path = os.path.splitext(predict_input_path)[0] + "_predicted.csv"
            df.to_csv(out_path, index=False)
            messagebox.showinfo("Prediction Complete", f"Predictions saved to:\n{out_path}")
        except Exception as e:
            messagebox.showerror("Prediction Failed", str(e))


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
        text="This is a pseudo-labeled dataset",
        variable=is_pseudolabel_var
    )
    pseudo_check.grid(row=0, column=2, padx=10, pady=5, sticky="w")


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

    tk.Label(train_tab, text="Step 4: Train/Test Parameters").grid(row=6, column=0, sticky="w", padx=10, pady=5)
    tk.Button(train_tab, text="Set Parameters / Train the Model", command=open_train_settings).grid(row=6, column=1, padx=5, pady=5)

    train_button = tk.Button(train_tab, text="Train Model", command=train_model, bg="#CCE5FF")
    train_button.grid(row=7, column=0, columnspan=2, pady=10)

    tk.Label(train_tab, text="Trainable File Info (Origin Tracking)").grid(row=8, column=0, columnspan=2, sticky="w", padx=10, pady=(15, 5))
    origin_info = tk.Text(train_tab, height=4, width=70, state="disabled", wrap="word")
    origin_info.grid(row=9, column=0, columnspan=2, padx=10, pady=5)
    # -- Training tab and prediction tab UI (end reminder buttons) --
    tk.Label(train_tab, text="(🔜) Combine Datasets for Training").grid(row=10, column=0, columnspan=2, sticky="w", padx=10, pady=(15, 5))
    tk.Button(train_tab, text="[Placeholder] Combine Datasets").grid(row=11, column=0, columnspan=2, padx=10, pady=5)
    
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

    predict_button = tk.Button(predict_tab, text="Run Prediction", bg="#D5F5E3", command=run_prediction)
    predict_button.grid(row=4, column=0, columnspan=2, pady=10)

    # -- Training tab and prediction tab UI (end reminder buttons) --
    tk.Label(predict_tab, text="(🔜) Combine Datasets for Prediction").grid(row=5, column=0, columnspan=2, sticky="w", padx=10, pady=(15, 5))
    tk.Button(predict_tab, text="[Placeholder] Combine Datasets").grid(row=6, column=0, columnspan=2, padx=10, pady=5)

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


    """
    # Row 4: output
    tk.Label(build_tab, text="5) Output (auto-named):").grid(row=4, column=0, sticky="w", padx=10, pady=6)
    tk.Label(build_tab, textvariable=out_path_var, anchor="w", justify="left").grid(row=4, column=1, columnspan=2, sticky="w", padx=6, pady=6)
    """
    # Helpers
    def _read_any_table(p):
        # try TSV first, then CSV with common settings
        try:
            return pd.read_csv(p, sep="\t", engine="python")
        except Exception:
            return pd.read_csv(p, engine="python")

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


    # Action
    def _build_from_pseudolabels():
        p_path = pseudo_path_var.get().strip()      # long-form pseudolabels (full TSV/CSV)
        f_path = features_path_var.get().strip()    # existing wide feature CSV (optional if rebuilding all)
        if not p_path or not os.path.exists(p_path):
            messagebox.showwarning("Missing pseudolabels", "Please choose a pseudolabeled TSV/CSV (long form).")
            return
        if not rebuild_all_features.get() and (not f_path or not os.path.exists(f_path)):
            messagebox.showwarning("Missing features", "Choose a wide feature CSV (or enable 'Rebuild ALL features').")
            return

        try:
            df_full = _read_any_table(p_path)     # full long-form pseudolabels
            df_feat = pd.read_csv(f_path) if (f_path and os.path.exists(f_path)) else None
        except Exception as e:
            messagebox.showerror("Read error", str(e)); return

        # --- scan column detection ---
        scan_col = _choose_scan_col(df_full) or "MS2scan_no"
        if df_feat is not None:
            scan_feat = _choose_scan_col(df_feat) or "MS2scan_no"

        # --- positives (keep best composition per scan after your score/ppm filters) ---
        pos = df_full.copy()
        if "ion score" in pos.columns:
            pos = pos[pos["ion score"].astype(float) >= float(ion_score_min.get())]
        if "ppm_error" in pos.columns:
            pos = pos[pos["ppm_error"].abs().astype(float) <= float(ppm_abs_max.get())]
        # rank by score (or |ppm|)
        use_score = "ion score" in pos.columns
        sort_cols = [scan_col] + (["ion score"] if use_score else ["ppm_error"])
        ascending = [True] + ([False] if use_score else [True])
        pos = pos.sort_values(by=sort_cols, ascending=ascending).groupby(scan_col, as_index=False).head(int(topn_var.get()))
        # map composition -> chosen label col
        labcol = label_target.get().strip()

        if "composition" in pos.columns:
            # DEBUG: print dtype and a few sample types/values
            try:
                sample_vals = pos["composition"].head(5).tolist()
                sample_types = [type(v).__name__ for v in sample_vals]
                logger.log(f"[PL] composition dtype={pos['composition'].dtype}; sample types={sample_types}; samples={sample_vals}")
                print("[11111debug]")
                print({pos['composition'].dtype})
            except Exception:
                pass

            # Prefer comp_tuple if present (most reliable)
            if "comp_tuple" in pos.columns:
                logger.log("[PL] Using 'comp_tuple' to build FHNSGKDN labels")
                pos[labcol] = pos["comp_tuple"].apply(_tuple_to_FHNSGKDN)

            else:
                # Try to coerce 'composition' to tuple
                tuples = pos["composition"].apply(_coerce_comp_to_tuple)
                n_tuples = tuples.notna().sum()
                logger.log(f"[PL] tuple coercion from 'composition' → {n_tuples} rows")

                if n_tuples > 0:
                    pos[labcol] = tuples.apply(lambda t: _tuple_to_FHNSGKDN(t) if t is not None else "")
                else:
                    # If already in FHNSGKDN format, just use it; otherwise fallback to str
                    looks_compact = pos["composition"].astype(str).str.match(_FHNSKDN_RE).all()
                    if looks_compact:
                        logger.log("[PL] 'composition' already in FHNSGKDN format; using as-is")
                        pos[labcol] = pos["composition"].astype(str)
                    else:
                        logger.log("[PL] 'composition' not tuple-like; falling back to string cast")
                        pos[labcol] = pos["composition"].astype(str)

        elif labcol not in pos.columns:
            messagebox.showerror("No composition",
                                "Neither 'composition' nor the chosen label column exist in the pseudolabels.")
            return

        pos_labels = pos[[scan_col, labcol]].drop_duplicates()        
        
        """
        # map composition -> chosen label col
        labcol = label_target.get().strip()
        if "composition" in pos.columns:
            # DEBUG: print dtype and a few sample types/values
            try:
                sample_vals = pos["composition"].head(5).tolist()
                sample_types = [type(v).__name__ for v in sample_vals]
                logger.log(f"[PL] composition dtype={pos['composition'].dtype}; sample types={sample_types}; samples={sample_vals}")
            except Exception:
                pass

            # Prefer comp_tuple if present (most reliable)
            if "comp_tuple" in pos.columns:
                logger.log("[PL] Using 'comp_tuple' to build FHNSGKDN labels")
                pos[labcol] = pos["comp_tuple"].apply(_tuple_to_FHNSGKDN)

            else:
                # Try to coerce 'composition' to tuple
                tuples = pos["composition"].apply(_coerce_comp_to_tuple)
                n_tuples = tuples.notna().sum()
                logger.log(f"[PL] tuple coercion from 'composition' → {n_tuples} rows")

                if n_tuples > 0:
                    pos[labcol] = tuples.apply(lambda t: _tuple_to_FHNSGKDN(t) if t is not None else "")
                else:
                    # If already in FHNSGKDN format, just use it; otherwise fallback to str
                    looks_compact = pos["composition"].astype(str).str.match(_FHNSKDN_RE).all()
                    if looks_compact:
                        logger.log("[PL] 'composition' already in FHNSGKDN format; using as-is")
                        pos[labcol] = pos["composition"].astype(str)
                    else:
                        logger.log("[PL] 'composition' not tuple-like; falling back to string cast")
                        pos[labcol] = pos["composition"].astype(str)

        elif labcol not in pos.columns:
            messagebox.showerror("No composition",
                                "Neither 'composition' nor the chosen label column exist in the pseudolabels.")
            return
        """
        """
        if "composition" in pos.columns:
            print(f"[debug]:type of composition")
            print(pos["composition"].dtype)
            if pos["composition"].dtype == "object" and any(isinstance(x, (tuple, list)) for x in pos["composition"].head(5)):
                # convert tuple → compact string like F1H4N2S3
                pos[labcol] = pos["composition"].apply(_tuple_to_FHNSGKDN)
            else:
                pos[labcol] = pos["composition"].astype(str)
        elif labcol not in pos.columns:
            messagebox.showerror("No composition", "Neither 'composition' nor the chosen label column exist in the pseudolabels.")
            return
            """
        #if "composition" in pos.columns:
        #    pos[labcol] = pos["composition"].astype(str)
        #elif labcol not in pos.columns:
        #    messagebox.showerror("No composition", "Neither 'composition' nor the chosen label column exist in the pseudolabels."); return
        #pos_labels = pos[[scan_col, labcol]].drop_duplicates()

        #20250911 
        ion_path = ionlist_path_var2.get().strip()
        try:
            ion_masses = _load_ion_masses_from_file(ion_path)
        except:
            ion_masses = None
        if not ion_masses:
            ion_masses = _infer_ion_masses_from_feature_df(df_feat)

        print(f"[features] using {len(ion_masses)} ion masses "
            f"(first 8: {ion_masses[:8] if ion_masses else []})")
        if not ion_masses:
            messagebox.showerror("No ion masses",
                "Could not obtain ion masses from ion list file or feature CSV headers.")
            return
        if ion_masses is None:
            if rebuild_all_features.get():
                # we will rebuild ALL features from long-form peaks using an ion list
                ion_df = _read_ion_df(ionlist_path_var2.get().strip())
                if ion_df is None or "mass" not in ion_df.columns:
                    messagebox.showerror("Ion list required", "Provide an ion list with a 'mass' column to rebuild features."); return
                ion_masses = ion_df["mass"].astype(float).tolist()
            else:
                # use the existing wide matrix's feature columns as the ion set (drop admin/label)
                drop_cols = {
                    labcol, "ID", "Source", "IUPACname(optional)", "Glycanannotation2",
                    "GlyToucan ID", "unique_ID", scan_feat, "MS2scan_no", "protonatedmass",
                    "theoretical_mass", "observed_mass", "ppm_error", "ion score",
                    "ion hit count", "ion hits m/z", "ion hits intensity", "ion hits logI", "ion hits relI"
                }
                ion_masses = []
                for c in df_feat.columns:
                    if c in drop_cols:
                        continue
                    try:
                        ion_masses.append(float(c))
                    except Exception:
                        # ignore non-numeric headers
                        continue

        """
        # --- decide feature source & ion masses ---
        ion_masses = None
        if rebuild_all_features.get():
            # we will rebuild ALL features from long-form peaks using an ion list
            ion_df = _read_ion_df(ionlist_path_var2.get().strip())
            if ion_df is None or "mass" not in ion_df.columns:
                messagebox.showerror("Ion list required", "Provide an ion list with a 'mass' column to rebuild features."); return
            ion_masses = ion_df["mass"].astype(float).tolist()
        else:
            # use the existing wide matrix’ feature columns as the ion set (drop admin/label)
            drop_cols = {labcol, "ID", "Source", "IUPACname(optional)", "Glycanannotation2", "GlyToucan ID", "unique_ID"}
            ion_masses = [c for c in df_feat.columns if c not in drop_cols and c != scan_feat]
        """

        # --- build NG from FULL long-form (complement of positives) ---
        # 1) unique scans from the long TSV (carry only peaks we need to build features)
        all_scans_unique = (df_full.drop_duplicates(subset=[scan_col])
                            [[scan_col, "peaklist", "peakintensity"]]
                            .copy())

        # 2) complement = all scans not in the positive set
        pos_set = set(pd.to_numeric(pos_labels[scan_col], errors="coerce").astype("Int64").dropna().tolist())
        neg_scans_df = all_scans_unique[~pd.to_numeric(all_scans_unique[scan_col], errors="coerce")
                                        .astype("Int64").isin(pos_set)].reset_index(drop=True)

        # 3) build NG features with the same normalization (log10(I)+1; miss=1.0)
        #20250911 ver
        ng_feat = build_features_from_peaks_log10_plus1(neg_scans_df, ion_masses, ppm=float(ppm_abs_max.get()))

        #ng_feat = build_features_from_peaks_log10_plus1(neg_scans_df, ion_masses, ppm=float(ppm_abs_max.get()))

        # 4) assemble NG table (rename scan col and set label)
        ng_df = pd.concat([neg_scans_df[[scan_col]].reset_index(drop=True), ng_feat], axis=1)
        ng_df = ng_df.rename(columns={scan_col: "MS2scan_no"})
        ng_df["Structure"] = "Non-glycan"
        if labcol != "Structure":
            ng_df = ng_df.rename(columns={"Structure": labcol})

        """
        # --- build NG candidates (low-score + no-hit) with reference normalization ---
        ng_df = pd.DataFrame(columns=[scan_col])
        try:
            from ml_ng_utils import collect_ng_candidates, build_features_from_peaks_log10_plus1
        except Exception:
            messagebox.showerror("Import error", "Could not import ml_ng_utils helpers."); return

        if use_low_score_ng.get() or use_nohit_ng.get():
            ng_df = collect_ng_candidates(
                df_pseudo_full=df_full,
                df_pseudo_filtered_pos=pos,
                ion_masses=ion_masses,
                ppm=float(ppm_abs_max.get()),
                low_score_col="ion score",
                low_score_cut=float(low_score_cut.get()),
                use_low_score=bool(use_low_score_ng.get()),
                use_no_hit=bool(use_nohit_ng.get()),
                manual_unknown_df=None,
                scan_col=scan_col,
            )
            # ng_df contains features (rebuilt with log10(I)+1) and Structure='Non-glycan'
            # ensure label column name matches UI selection
            if "Structure" in ng_df.columns and labcol != "Structure":
                ng_df = ng_df.rename(columns={"Structure": labcol})

        # --- assemble final wide matrix ---
        
        if rebuild_all_features.get():
            # positives: rebuild features directly from long-form to match NG normalization
            feat_pos = build_features_from_peaks_log10_plus1(df_full.set_index(scan_col).loc[pos_labels[scan_col]].reset_index(),
                                                            ion_masses, ppm=float(ppm_abs_max.get()))
            feat_pos.insert(0, "MS2scan_no", pos_labels[scan_col].values)
            feat_pos[labcol] = pos_labels[labcol].values
            wide = feat_pos
            if not ng_df.empty:
                wide = pd.concat([wide, ng_df], ignore_index=True)
        else:
            # keep the existing wide features for positives + merge labels; then append NG (already wide)
            wide = df_feat.copy()
            if labcol not in wide.columns:
                wide[labcol] = None
            # join labels into wide
            wide = wide.merge(pos_labels, how="left", left_on=scan_feat, right_on=scan_col, suffixes=("", "_pl"))
            wide[labcol] = wide[labcol + "_pl"].where(wide[labcol + "_pl"].notna(), wide[labcol])
            wide = wide.drop(columns=[c for c in [labcol + "_pl", scan_col] if c in wide.columns])

            # append NG rows that aren’t already in wide
            if not ng_df.empty:
                already = set(wide[scan_feat].astype(int).tolist())
                ng_add = ng_df[~ng_df["MS2scan_no"].astype(int).isin(already)].copy()
                # align columns
                for c in wide.columns:
                    if c not in ng_add.columns:
                        ng_add[c] = None
                ng_add = ng_add[wide.columns]
                wide = pd.concat([wide, ng_add], ignore_index=True)
        
        # optional: fill unlabeled as 'None'
        if fill_unlabeled_as_none.get():
            wide[labcol] = wide[labcol].fillna("None")
        """

        # positives: rebuild features directly from long-form to match NG normalization 20250830
        # 1) make one row per positive scan (carry peaks only)
        scans_unique = (df_full.drop_duplicates(subset=[scan_col])
                        [[scan_col, "peaklist", "peakintensity"]]
                        .copy())

        # ensure consistent dtypes for the join
        scans_unique[scan_col] = pd.to_numeric(scans_unique[scan_col], errors="coerce").astype("Int64")
        pos_labels[scan_col]   = pd.to_numeric(pos_labels[scan_col],   errors="coerce").astype("Int64")

        # pull exactly the positive scans, preserving order of pos_labels
        scans_for_pos = scans_unique.set_index(scan_col).loc[pos_labels[scan_col]].reset_index()

        # 2) build features (log10(I)+1; miss=1.0) with the agreed ion set
        feat_pos = build_features_from_peaks_log10_plus1(
            scans_for_pos, ion_masses, ppm=float(ppm_abs_max.get())
        )

        # 3) insert aligned scan id and label (sizes match)
        feat_pos.insert(0, "MS2scan_no", scans_for_pos[scan_col].to_numpy())
        feat_pos[labcol] = pos_labels.set_index(scan_col).loc[scans_for_pos[scan_col], labcol].to_numpy()

        wide = feat_pos
        if not ng_df.empty:
            wide = pd.concat([wide, ng_df], ignore_index=True)


        # --- status readout: counts ---
        try:
            # positives = unique scans that survived score/ppm and were labeled
            pos_count = int(pos_labels[scan_col].nunique()) if 'pos_labels' in locals() and not pos_labels.empty else 0
        except Exception:
            pos_count = 0

        ng_count = int(len(ng_df)) if 'ng_df' in locals() and ng_df is not None and not ng_df.empty else 0
        total_count = int(len(wide))

        # IMPORTANT: NG is derived from the full long-form table, excluding positive scans.
        # This ensures NG is NOT a subset of the 2,481 filtered rows, but from the complement in df_full.
        build_status_var.set(f"Built dataset → positives={pos_count}, NG={ng_count}, total={total_count}")

        _save_and_autoload(wide, f_path or p_path)

    #tk.Button(build_tab, text="Build Trainable CSV", command=_build_from_pseudolabels, bg="#E6FFE6").grid(row=5, column=0, columnspan=3, pady=12)

    
    # Build button → row 7 (status label sits on row 6)
    tk.Button(build_tab, text="Build Trainable CSV", command=_build_from_pseudolabels, bg="#E6FFE6").grid(
        row=7, column=0, columnspan=3, pady=12
    )



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
root.title("GlycoMSP File Manager GUI v0.5")
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


# --- Add Analysis Tools Frame ---

analysis_frame = tk.LabelFrame(root, text="Analysis Tools", padx=10, pady=10)
analysis_frame.pack(padx=10, pady=10, fill="x")

convert_button = tk.Button(analysis_frame, text="Convert Raw to CSV",command=launch_metadata_batch)
convert_button.pack(side="left", padx=5)
tk.Button(analysis_frame, text="Prepare Dataset", command=open_prepare_dataset_window).pack(side="left", padx=5)
tk.Button(analysis_frame, text="Run ML Analysis", command=open_ml_analysis_window).pack(side="left", padx=5)


# Run the GUI
root.mainloop()



