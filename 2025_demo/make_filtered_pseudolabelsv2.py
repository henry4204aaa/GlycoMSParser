#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, ast, math, os
from typing import Optional, Tuple, List, Dict

import numpy as np
import pandas as pd

PROTON_MASS = 1.007276466812

# ---------------- I/O helpers ----------------
def read_table_auto(path: str) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext in [".xls", ".xlsx"]:
        return pd.read_excel(path)
    if ext in [".tsv", ".txt"]:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)

# ---------------- parsing helpers ----------------
def parse_tuple_like(x) -> Optional[Tuple[int, ...]]:
    if pd.isna(x): return None
    s = str(x).strip().replace("[","(").replace("]",")")
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1]
    parts = [p.strip() for p in s.split(",") if p.strip()]
    try:
        return tuple(int(round(float(p))) for p in parts)
    except Exception:
        try:
            v = ast.literal_eval(str(x))
            if isinstance(v, (list, tuple)):
                return tuple(int(round(float(p))) for p in v)
        except Exception:
            return None
    return None

import re

# Convert strings like "F1H3N5KDN2S1" to (H, N, S, G, KDN, F)
def parse_shorthand_glycan(s: str) -> Optional[tuple]:
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return None
    text = str(s).strip().upper().replace(" ", "")
    if not text:
        return None

    # Capture tokens: KDN first (3 letters), then single-letter symbols
    pattern = re.compile(r"(KDN|F|H|N|S|G)(\d*)", re.IGNORECASE)
    counts = {"H":0, "N":0, "S":0, "G":0, "KDN":0, "F":0}

    matched_any = False
    for m in pattern.finditer(text):
        key = m.group(1).upper()
        num = m.group(2)
        val = int(num) if num else 1
        if key in counts:
            counts[key] += val
            matched_any = True

    if not matched_any:
        return None

    # Final order: (H, N, S, G, KDN, F)
    return (counts["H"], counts["N"], counts["S"], counts["G"], counts["KDN"], counts["F"])

def coerce_float(x):
    try: return float(x)
    except Exception: return np.nan

def ppm_error(x: float, ref: float) -> float:
    return abs(x - ref) / ref * 1e6

def neutral_from_mz_charge(mz: float, z: int) -> Optional[float]:
    if pd.isna(mz) or pd.isna(z) or z == 0: return None
    try: return mz * z - z * PROTON_MASS
    except Exception: return None

def isotope_neutral_variants(neutral_mass: float, include_mminus1: bool) -> List[Tuple[str, float]]:
    out = [("M", neutral_mass), ("M+1", neutral_mass + 1.003355)]
    if include_mminus1:
        out.append(("M-1", neutral_mass - 1.003355))
    return out

# ---------------- column alias resolver ----------------
def norm(s: str) -> str:
    return s.strip().lower().replace(" ", "_").replace("-", "_")

def build_colmap(cols: List[str]) -> Dict[str, str]:
    # map normalized -> original
    return {norm(c): c for c in cols}

def resolve_or_fail(colmap: Dict[str,str], aliases: List[str], kind: str) -> str:
    for a in aliases:
        k = norm(a)
        if k in colmap:
            return colmap[k]
    raise SystemExit(f"[ERROR] required column not found; tried aliases: {aliases} (kind={kind})")

def resolve_optional(colmap: Dict[str,str], aliases: List[str]) -> Optional[str]:
    for a in aliases:
        k = norm(a)
        if k in colmap:
            return colmap[k]
    return None

# Canonical alias lists
ALIASES = {
    "scan": ["MS2scan_no","ms2scan_no","ms2_scan","ms2scan","scan","MS2"],
    "composition": ["composition","tuple","pseudolabel_tuple","pseudolabel_comp","predicted_comp","composition_tuple_str"],
    "ion_score": ["ion_score","ion score","ionscore","ion_similarity","score","ion_match_score"],
    "comp_mass": ["comp_neutral_mass","composition_neutral_mass","comp_mass","composition_mass","comp mass"],
    "prec_neutral": ["precursor_neutral_mass","protonatedmass","precursor_mass","neutral_mass","neutral mass"],
    "prec_mz": ["precursor_mz","precursor mz","mz_precursor"],
    "charge": ["charge","z","precursor_charge"],
    "manual_tuple": ["manual_tuple","manual_comp","manual composition","manual_tuple_str"]
}

def pick_precursor_neutral(row) -> Optional[float]:
    for k in ("precursor_neutral_mass","protonatedmass","precursor_mass"):
        if k in row and pd.notna(row[k]):
            return coerce_float(row[k])
    z = row.get("charge", row.get("z", None))
    mz = row.get("precursor_mz", None)
    try: z = int(z) if pd.notna(z) else None
    except Exception: z = None
    mz = coerce_float(mz)
    if z and mz and not math.isclose(z,0.0):
        return neutral_from_mz_charge(mz, z)
    return None

# ---------------- overrides loader ----------------
def load_override_csv(path: str, score_min: float) -> pd.DataFrame:
    df = read_table_auto(path)
    df.columns = [c.strip() for c in df.columns]
    cmap = build_colmap(df.columns)

    comp_col = resolve_or_fail(cmap, ALIASES["composition"], "override:composition")
    mass_col = resolve_or_fail(cmap, ALIASES["comp_mass"], "override:comp_neutral_mass")
    ion_col = resolve_optional(cmap, ALIASES["ion_score"])

    out = pd.DataFrame({
        "composition": df[comp_col],
        "comp_neutral_mass": pd.to_numeric(df[mass_col], errors="coerce"),
    })
    out["composition_tuple"] = out["composition"].apply(parse_tuple_like)
    if ion_col is None:
        out["ion_score"] = score_min
    else:
        out["ion_score"] = pd.to_numeric(df[ion_col], errors="coerce").fillna(score_min)
    return out[["composition","composition_tuple","comp_neutral_mass","ion_score"]].copy()

# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--manual", required=False)
    ap.add_argument("--override-comps", required=False)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--ppm", type=float, default=20.0)
    ap.add_argument("--score-min", type=float, default=0.07)
    ap.add_argument("--big-mass", type=float, default=2500.0)
    ap.add_argument("--enable-mminus1", action="store_true")
    ap.add_argument("--manual-sheet", required=False, help="Excel sheet name for --manual (e.g., 'MSlist')")
    ap.add_argument("--manual-scan-col", required=False, help="Manual scan column name override (e.g., 'MS2scan_no')")
    ap.add_argument("--manual-tuple-col", required=False, help="Manual tuple column name override (e.g., 'manual_tuple')")
    ap.add_argument("--manual-shorthand-col", required=False, help="Manual shorthand column (e.g., 'F1H3N5KDN1'); will be converted to (H,N,S,G,KDN,F)")

    args = ap.parse_args()

    # --- load candidates (CSV/TSV/XLSX ok) ---
    df = read_table_auto(args.candidates)
    df.columns = [c.strip() for c in df.columns]
    cmap = build_colmap(df.columns)

    scan_col  = resolve_or_fail(cmap, ALIASES["scan"], "candidates:scan")
    comp_col  = resolve_or_fail(cmap, ALIASES["composition"], "candidates:composition")
    score_col = resolve_or_fail(cmap, ALIASES["ion_score"], "candidates:ion_score")

    compmass_col   = resolve_optional(cmap, ALIASES["comp_mass"])
    precN_col      = resolve_optional(cmap, ALIASES["prec_neutral"])
    precmz_col     = resolve_optional(cmap, ALIASES["prec_mz"])
    charge_col     = resolve_optional(cmap, ALIASES["charge"])

    # Coerce
    df["_scan"] = df[scan_col]
    df["_composition"] = df[comp_col]
    df["_composition_tuple"] = df["_composition"].apply(parse_tuple_like)
    df["_ion_score"] = pd.to_numeric(df[score_col], errors="coerce")

    # Optional masses
    if compmass_col is not None:
        df["_comp_neutral"] = pd.to_numeric(df[compmass_col], errors="coerce")
    else:
        df["_comp_neutral"] = np.nan

    # Precursor neutral: directly if present, otherwise derive from mz+z
    df["_precursor_neutral"] = np.nan
    if precN_col is not None:
        df["_precursor_neutral"] = pd.to_numeric(df[precN_col], errors="coerce")
    else:
        # try derive if we have mz and charge
        if precmz_col is not None:
            df["precursor_mz"] = pd.to_numeric(df[precmz_col], errors="coerce")
        else:
            df["precursor_mz"] = np.nan
        if charge_col is not None:
            df["charge"] = pd.to_numeric(df[charge_col], errors="coerce")
        else:
            df["charge"] = np.nan
        df["_precursor_neutral"] = df.apply(pick_precursor_neutral, axis=1)

    # --- manual (optional; CSV/TSV/XLSX) ---
    df_manual = None
    if args.manual:
        ext = os.path.splitext(args.manual)[1].lower()
        if ext in [".xls", ".xlsx"]:
            m = pd.read_excel(args.manual, sheet_name=(args.manual_sheet or 0))
        else:
            m = read_table_auto(args.manual)

        m.columns = [c.strip() for c in m.columns]
        m_map = build_colmap(m.columns)

        # Scan column (either explicit or via aliases)
        if args.manual_scan_col:
            m_scan = args.manual_scan_col
            if m_scan not in m.columns:
                raise SystemExit(f"[ERROR] manual scan column '{m_scan}' not found in sheet '{args.manual_sheet or 'default'}'")
        else:
            m_scan = resolve_or_fail(m_map, ALIASES["scan"], "manual:scan")

        # Choose source for manual composition:
        # 1) If user provided a tuple-col, use & parse it as tuple-like.
        # 2) Else if user provided a shorthand-col, convert it with parse_shorthand_glycan.
        # 3) Else try aliases for a tuple-like column; if fails, try to detect a shorthand-like column by content.
        manual_tuple_series = None

        if args.manual_tuple_col and args.manual_tuple_col in m.columns:
            manual_tuple_series = m[args.manual_tuple_col].apply(parse_tuple_like)

        elif args.manual_shorthand_col and args.manual_shorthand_col in m.columns:
            manual_tuple_series = m[args.manual_shorthand_col].apply(parse_shorthand_glycan)

        else:
            # Try alias-based tuple-like column first
            try:
                m_tuple_col = resolve_or_fail(m_map, ALIASES["manual_tuple"], "manual:manual_tuple")
                manual_tuple_series = m[m_tuple_col].apply(parse_tuple_like)
            except SystemExit:
                # Fallback: try to auto-detect a shorthand column by looking for F/H/N/S/G/KDN tokens
                # Pick the first column that yields at least one non-None parse
                candidate = None
                for col in m.columns:
                    sample = m[col].head(50).apply(parse_shorthand_glycan)
                    if sample.notna().any():
                        candidate = col
                        break
                if candidate is None:
                    raise
                manual_tuple_series = m[candidate].apply(parse_shorthand_glycan)

        df_manual = pd.DataFrame({
            "MS2scan_no": m[m_scan],
            "manual_tuple": manual_tuple_series
        })

    # --- MASS GATING over candidates ---
    kept = []
    for _, r in df.iterrows():
        scan = r["_scan"]; comp = r["_composition"]; comp_tuple = r["_composition_tuple"]
        score = r["_ion_score"]; comp_mass = r["_comp_neutral"]; precN = r["_precursor_neutral"]

        # pass-through if mass info missing (score filter later)
        if pd.isna(comp_mass) or pd.isna(precN):
            kept.append({"MS2scan_no": scan, "composition": comp, "composition_tuple": comp_tuple,
                         "ion_score": score, "delta_mass_ppm": np.nan,
                         "is_isotope_window": "NA_no_mass_info", "source": "input"})
            continue

        windows = [("M", precN)]
        if precN >= args.big_mass:
            windows = isotope_neutral_variants(precN, args.enable_mminus1)

        best_ppm = None; best_lab = None
        for lab, target in windows:
            ppm = ppm_error(comp_mass, target)
            if ppm <= args.ppm:
                if best_ppm is None or ppm < best_ppm:
                    best_ppm, best_lab = ppm, lab
        if best_lab is not None:
            kept.append({"MS2scan_no": scan, "composition": comp, "composition_tuple": comp_tuple,
                         "ion_score": score, "delta_mass_ppm": best_ppm,
                         "is_isotope_window": best_lab, "source": "input"})

    df_kept = pd.DataFrame(kept)

    # --- OVERRIDES (salvage) ---
    if args.override_comps:
        over = load_override_csv(args.override_comps, score_min=args.score_min)
        injected = []
        scans_prec = df[["_scan","_precursor_neutral"]].drop_duplicates()
        for _, rr in scans_prec.iterrows():
            scan = rr["_scan"]; precN = rr["_precursor_neutral"]
            if pd.isna(precN): continue
            windows = [("M", precN)]
            if precN >= args.big_mass:
                windows = isotope_neutral_variants(precN, args.enable_mminus1)
            for _, ov in over.iterrows():
                comp = ov["composition"]; comp_tuple = ov["composition_tuple"]
                comp_mass = ov["comp_neutral_mass"]; ov_score = ov["ion_score"]
                best_ppm = None; best_lab = None
                for lab, target in windows:
                    ppm = ppm_error(comp_mass, target)
                    if ppm <= args.ppm:
                        if best_ppm is None or ppm < best_ppm:
                            best_ppm, best_lab = ppm, lab
                if best_lab is not None:
                    injected.append({"MS2scan_no": scan, "composition": comp, "composition_tuple": comp_tuple,
                                     "ion_score": ov_score, "delta_mass_ppm": best_ppm,
                                     "is_isotope_window": best_lab, "source": "override"})
        if injected:
            df_kept = pd.concat([df_kept, pd.DataFrame(injected)], ignore_index=True)
            print(f"[INFO] Injected {len(injected)} override rows (score ≥ {args.score_min}).")

    # de-dupe within (scan, composition): keep highest score
    if not df_kept.empty:
        df_kept.sort_values(["MS2scan_no","composition","ion_score"], ascending=[True,True,False], inplace=True)
        df_kept = df_kept.drop_duplicates(subset=["MS2scan_no","composition"], keep="first").reset_index(drop=True)

    # score filter
    before = len(df_kept)
    df_kept = df_kept[df_kept["ion_score"] >= args.score_min].copy()
    print(f"[INFO] Score filter ≥ {args.score_min}: kept {len(df_kept)}/{before} rows.")

    # join manual
    if df_manual is not None and not df_kept.empty:
        df_kept = df_kept.merge(df_manual, on="MS2scan_no", how="left")
        df_kept["tuple_match"] = (df_kept["composition_tuple"].astype(str) == df_kept["manual_tuple"].astype(str))
        df_kept["manual_check"] = df_kept["ion_score"].notna()
    else:
        df_kept["manual_tuple"] = np.nan
        df_kept["tuple_match"] = np.nan
        df_kept["manual_check"] = df_kept["ion_score"].notna()

    detailed_out = f"{args.out_prefix}_detailed.csv"
    df_kept.to_csv(detailed_out, index=False)
    print(f"[OK] Wrote detailed rows → {detailed_out} (rows={len(df_kept)})")

    # per-scan summary
    if df_kept.empty:
        summary = pd.DataFrame(columns=["MS2scan_no","n_candidates","any_match","manual_tuple","best_candidate","best_score","best_is_isotope_window"])
    else:
        df_kept["score_rank"] = df_kept.groupby("MS2scan_no")["ion_score"].rank("dense", ascending=False)
        idx_best = (df_kept.sort_values(["MS2scan_no","ion_score"], ascending=[True,False])
                    .groupby("MS2scan_no").head(1).reset_index(drop=True))
        g = df_kept.groupby("MS2scan_no")
        summary = pd.DataFrame({
            "MS2scan_no": g.size().index,
            "n_candidates": g.size().values,
            "any_match": g["tuple_match"].any().values if "tuple_match" in df_kept.columns else np.nan,
            "manual_tuple": g["manual_tuple"].first().values
        })
        summary = summary.merge(
            idx_best[["MS2scan_no","composition","ion_score","is_isotope_window"]],
            on="MS2scan_no", how="left"
        ).rename(columns={"composition":"best_candidate","ion_score":"best_score","is_isotope_window":"best_is_isotope_window"})
    summary_out = f"{args.out_prefix}_summary.csv"
    summary.to_csv(summary_out, index=False)
    print(f"[OK] Wrote per-scan summary → {summary_out} (rows={len(summary)})")

    if "tuple_match" in df_kept.columns and not df_kept.empty:
        print("[INFO] Coverage (after filters):", df_kept["tuple_match"].value_counts(dropna=False).to_dict())

if __name__ == "__main__":
    main()