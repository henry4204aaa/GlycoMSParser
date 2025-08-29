
import argparse
import ast
import math
import os
import re
from typing import List, Optional, Tuple, Union

import pandas as pd
import numpy as np

STRUCT_CANDIDATES = [
    "Structure","structure","ShortStructure","short_structure",
    "composition","Composition","Struct","label","Predicted_Label"
]
SCAN_CANDIDATES = [
    "MS2scan_no","ms2scan_no","ScanNo","scan","scan_no","Scan","scan_number"
]

MetaColsPreferred = [
    "MS2scan_no", "protonatedmass", "Structure", "Source",
    "score", "pseudolabel_score"
]

def _is_short_structure(s: str) -> bool:
    return bool(re.match(r"^[A-Za-z]+[0-9]+", s))

def parse_tuple_structure(s: Union[str, tuple, list]) -> Optional[Tuple[int,int,int,int,int,int]]:
    if isinstance(s, (tuple, list)) and len(s) == 6:
        try:
            return tuple(int(v) for v in s)  # type: ignore
        except Exception:
            return None
    if isinstance(s, str):
        s = s.strip()
        if _is_short_structure(s):
            return None
        try:
            parsed = ast.literal_eval(s)
            if isinstance(parsed, (tuple, list)) and len(parsed) == 6:
                return tuple(int(v) for v in parsed)
        except Exception:
            return None
    return None

def tuple_to_short_composition(t: Tuple[int,int,int,int,int,int]) -> str:
    H, N, S, G, KDN, F = t
    parts = []
    if F > 0:   parts.append(f"F{F}")
    if H > 0:   parts.append(f"H{H}")
    if N > 0:   parts.append(f"N{N}")
    if S > 0:   parts.append(f"S{S}")
    if G > 0:   parts.append(f"G{G}")
    if KDN > 0: parts.append(f"KDN{KDN}")
    return "".join(parts) if parts else "Non-glycan"

def normalize_structure_series(series: pd.Series) -> pd.Series:
    out = []
    for v in series.fillna("").tolist():
        tpl = parse_tuple_structure(v)
        if tpl is None:
            sv = str(v).strip()
            out.append(sv if sv else "Non-glycan")
        else:
            out.append(tuple_to_short_composition(tpl))
    return pd.Series(out, index=series.index, name="Structure")

def find_first_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None

def read_ion_masses(path: str, sheet_name: str="ionlist") -> List[float]:
    p = path.lower()
    if p.endswith(".xlsx") or p.endswith(".xls"):
        xls = pd.ExcelFile(path)
        if sheet_name not in xls.sheet_names:
            if "ionlist" in [s.lower() for s in xls.sheet_names]:
                sheet_name = xls.sheet_names[[s.lower() for s in xls.sheet_names].index("ionlist")]
            else:
                sheet_name = xls.sheet_names[0]
        ion_df = xls.parse(sheet_name)
    else:
        ion_df = pd.read_csv(path)
    col = None
    for c in ion_df.columns:
        if c.strip().lower() == "mass":
            col = c
            break
    if col is None:
        raise ValueError("Ion sheet must contain a 'mass' column.")
    masses = ion_df[col].dropna().astype(float).tolist()
    return masses

def build_features_from_peaks(df: pd.DataFrame, ion_masses: List[float], ppm: float=20.0, unit_norm: bool=False) -> pd.DataFrame:
    if "peaklist" not in df.columns or "peakintensity" not in df.columns:
        raise ValueError("Rebuild features requested, but 'peaklist' and/or 'peakintensity' not found in pseudolabeled table.")
    feature_cols = [str(m) for m in ion_masses]
    out = np.zeros((len(df), len(ion_masses)), dtype=float)
    for i, (_, row) in enumerate(df.iterrows()):
        try:
            peaks = np.array(ast.literal_eval(str(row["peaklist"])), dtype=float)
            intens = np.array(ast.literal_eval(str(row["peakintensity"])), dtype=float)
        except Exception:
            continue
        for j, target in enumerate(ion_masses):
            tol = target * ppm / 1e6
            mask = np.abs(peaks - target) <= tol
            if mask.any():
                I = float(np.max(intens[mask]))
                out[i, j] = math.log10(I) + 1
            else:
                out[i, j] = 1.0
        if unit_norm:
            norm = np.linalg.norm(out[i, :])
            if norm > 0:
                out[i, :] /= norm
    feat_df = pd.DataFrame(out, index=df.index, columns=feature_cols)
    return feat_df

def collect_existing_feature_columns(df: pd.DataFrame, meta_cols: Optional[List[str]]=None) -> List[str]:
    if meta_cols is None:
        meta_cols = MetaColsPreferred
    feats = []
    for c in df.columns:
        if c in meta_cols:
            continue
        lc = c.lower()
        if any(k in lc for k in ["structure", "source", "uid", "unique_id", "iupac", "glycanannotation", "glytoucan"]):
            continue
        if lc in ["protonatedmass", "ms2scan_no", "peaklist", "peakintensity"]:
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            feats.append(c)
            continue
        try:
            df[c].astype(float)
            feats.append(c)
        except Exception:
            pass
    return feats

def salvage_labels_with_manual(df: pd.DataFrame, manual_csv: str, scan_col_in_df: str) -> pd.DataFrame:
    try:
        mdf = pd.read_csv(manual_csv)
    except Exception as e:
        print(f"[salvage] Failed to read manual CSV: {e}")
        return df
    scan_col = find_first_col(mdf, SCAN_CANDIDATES)
    if scan_col is None:
        print("[salvage] Manual CSV lacks a scan number column; skip salvage.")
        return df
    struct_col = find_first_col(mdf, STRUCT_CANDIDATES)
    if struct_col is None:
        print("[salvage] Manual CSV lacks a Structure-like column; skip salvage.")
        return df
    mdf_norm = mdf[[scan_col, struct_col]].copy()
    mdf_norm.columns = ["MS2scan_no_manual", "Structure_manual"]
    mdf_norm["Structure_manual"] = normalize_structure_series(mdf_norm["Structure_manual"])
    out = df.merge(mdf_norm, left_on=scan_col_in_df, right_on="MS2scan_no_manual", how="left")
    out["Structure"] = out["Structure_manual"].combine_first(out.get("Structure"))
    out = out.drop(columns=["MS2scan_no_manual","Structure_manual"])
    return out

def make_trainable(
    pseudotab: str,
    out_csv: str,
    score_threshold: Optional[float]=None,
    salvage_manual_csv: Optional[str]=None,
    keep_source: bool=False,
    ion_sheet: Optional[str]=None,
    ion_sheet_name: str="ionlist",
    ppm: float=20.0,
    unit_norm: bool=False,
) -> str:
    df = pd.read_csv(pseudotab, sep='\t' if pseudotab.lower().endswith(".tsv") else ",")
    scan_col = find_first_col(df, SCAN_CANDIDATES)
    if scan_col is None:
        raise ValueError(f"Input missing scan column (tried: {SCAN_CANDIDATES}). Columns present: {list(df.columns)}")
    struct_col = find_first_col(df, STRUCT_CANDIDATES)
    if struct_col is None:
        df["Structure"] = "Non-glycan"
        struct_col = "Structure"
    df["Structure"] = normalize_structure_series(df[struct_col])
    sc = None
    for c in ("score","pseudolabel_score","pl_score","Score", "ion score"):
        if c in df.columns:
            sc = c; break
    if score_threshold is not None and sc is not None:
        df = df[df[sc].astype(float) >= float(score_threshold)].copy()
    elif score_threshold is not None and sc is None:
        print("[warn] No score column found; threshold ignored.")
    if salvage_manual_csv:
        df = salvage_labels_with_manual(df, salvage_manual_csv, scan_col_in_df=scan_col)
    keep_cols = []
    if "protonatedmass" in df.columns:
        keep_cols.append("protonatedmass")
    keep_cols.append(scan_col)
    keep_cols.append("Structure")
    if keep_source and "Source" in df.columns:
        keep_cols.append("Source")
    if ion_sheet is not None:
        ion_masses = read_ion_masses(ion_sheet, sheet_name=ion_sheet_name)
        feat_df = build_features_from_peaks(df, ion_masses, ppm=ppm, unit_norm=unit_norm)
        out_df = pd.concat([df[keep_cols].copy(), feat_df], axis=1)
    else:
        feature_cols = collect_existing_feature_columns(df, meta_cols=list(set(MetaColsPreferred + keep_cols)))
        out_df = df[keep_cols + feature_cols].copy()
        for c in feature_cols:
            out_df[c] = pd.to_numeric(out_df[c], errors="coerce").fillna(0.0)
    out_df = out_df.rename(columns={scan_col: "MS2scan_no"})
    out_df.to_csv(out_csv, index=False)
    return out_csv

def convert_pseudolabels(
    pseudotab: str,
    out_prefix: str,
    manual_csv: Optional[str]=None,
    keep_source: bool=False,
    ion_sheet: Optional[str]=None,
    ion_sheet_name: str="ionlist",
    ppm: float=20.0,
    unit_norm: bool=False,
) -> dict:
    res = {}
    a_path = f"{out_prefix}_score007.csv"
    b_path = f"{out_prefix}_score006.csv"
    c_path = f"{out_prefix}_score007_salvage.csv"
    res["score007"] = make_trainable(
        pseudotab, a_path, score_threshold=0.07, salvage_manual_csv=None,
        keep_source=keep_source, ion_sheet=ion_sheet, ion_sheet_name=ion_sheet_name,
        ppm=ppm, unit_norm=unit_norm
    )
    res["score006"] = make_trainable(
        pseudotab, b_path, score_threshold=0.06, salvage_manual_csv=None,
        keep_source=keep_source, ion_sheet=ion_sheet, ion_sheet_name=ion_sheet_name,
        ppm=ppm, unit_norm=unit_norm
    )
    res["score007_salvage"] = make_trainable(
        pseudotab, c_path, score_threshold=0.07, salvage_manual_csv=manual_csv,
        keep_source=keep_source, ion_sheet=ion_sheet, ion_sheet_name=ion_sheet_name,
        ppm=ppm, unit_norm=unit_norm
    )
    return res

def main():
    ap = argparse.ArgumentParser(description="Convert pseudolabeled TSV to trainable CSVs with ion sheet, thresholds, and salvage (v2.1 tolerant columns).")
    ap.add_argument("--pseudotab", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--manual-csv", default=None, help="Manual MSlist or full CSV including scan + Structure.")
    ap.add_argument("--keep-source", action="store_true")
    ap.add_argument("--ion-sheet", default=None, help="Excel/CSV containing an ion 'mass' column.")
    ap.add_argument("--ion-sheet-name", default="ionlist", help="Excel sheet name for ion masses (if Excel).")
    ap.add_argument("--ppm", type=float, default=20.0, help="PPM tolerance for feature extraction.")
    ap.add_argument("--unit-norm", action="store_true", help="Per-spectrum unit L2 normalization after log1p.")
    args = ap.parse_args()

    res = convert_pseudolabels(
        pseudotab=args.pseudotab,
        out_prefix=args.out_prefix,
        manual_csv=args.manual_csv,
        keep_source=args.keep_source,
        ion_sheet=args.ion_sheet,
        ion_sheet_name=args.ion_sheet_name,
        ppm=args.ppm,
        unit_norm=args.unit_norm,
    )
    print("Generated files:")
    for k, v in res.items():
        print(f"  {k}: {v}")

if __name__ == "__main__":
    main()
