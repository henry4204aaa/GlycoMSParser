import os, json, pandas as pd
from mspfileloaderv10 import extract_ion_intensities

def read_fragment_masses_any(ion_path: str, sheet_name: str = "ionlist"):
    if not ion_path or not os.path.exists(ion_path):
        raise FileNotFoundError(f"Ion list not found: {ion_path}")
    ext = os.path.splitext(ion_path)[1].lower()
    if ext in (".csv", ".tsv", ".txt"):
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

def create_unlabeled_from_method(method_path: str, default_ppm: int | None = None) -> str:
    with open(method_path, "r", encoding="utf-8") as f:
        m = json.load(f)

    raw_csv  = (m.get("parents") or {}).get("converted_csv")
    ion_path = (m.get("ionlist") or {}).get("path")
    ion_sheet= (m.get("ionlist") or {}).get("sheet") or "ionlist"
    ppm      = default_ppm or (m.get("ionlist") or {}).get("ppm_tolerance") or 20

    if not raw_csv or not os.path.exists(raw_csv):
        raise FileNotFoundError(f"Converted CSV not found (method): {raw_csv}")
    if not ion_path or not os.path.exists(ion_path):
        raise FileNotFoundError(f"Ion list not found (method): {ion_path}")

    frags = read_fragment_masses_any(ion_path, sheet_name=ion_sheet)

    # IMPORTANT: reuse your existing extractor used by the GUI
    # If it's a function in this file, call it here. Example:
    #   feature_df = extract_ion_intensities(raw_csv, frags, ppm=float(ppm))
    # If it's a class method, pull it out to module scope or make a thin wrapper.
    feature_df = extract_ion_intensities(raw_csv, frags, ppm=float(ppm))  # <- your existing function

    # Try to add UID if columns are present (safe no-op otherwise)
    try:
        import hashlib
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

    out_path = os.path.splitext(raw_csv)[0] + f"_unlabeled_ppm{ppm}.csv"
    feature_df.to_csv(out_path, index=False)
    return out_path

# make it importable
__all__ = ["read_fragment_masses_any", "create_unlabeled_from_method"]


"""
import sys, importlib
sys.path.insert(0, r"G:/其他電腦/My Computer/GlycoMSParser/2025_demo")
import mspfileloaderv10 as v10
importlib.reload(v10)

method = r"C:/Users/Sakazuki/Desktop/Khoolab_2025data/U937cells_NG_new/OG/train_fillUID/U937_COKOST1OE_OGneu.method.json"
print("Building from:", method)
out = v10.create_unlabeled_from_method(method, default_ppm=20)
print("Wrote:", out)
"""