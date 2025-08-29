
import argparse
import sys
import pandas as pd
import numpy as np
import re

PPM_DEFAULT = 20.0
ORDER_T = ["Hex","HexNAc","NeuAc","NeuGc","KDN","Fuc"]  # target order for tuple

ALIASES = {
    "H":"Hex", "HEX":"Hex",
    "N":"HexNAc", "HEXNAC":"HexNAc",
    "S":"NeuAc", "NEU5AC":"NeuAc", "NEUAC":"NeuAc", "SA":"NeuAc",
    "G":"NeuGc", "NEU5GC":"NeuGc", "NEUGC":"NeuGc",
    "KDN":"KDN", "K":"KDN",
    "F":"Fuc", "FUC":"Fuc"
}
token_re = re.compile(r"([A-Za-z]+)\s*([+-]?\d+)", re.IGNORECASE)

def manual_to_tuple(s: str):
    counts = {k:0 for k in ORDER_T}
    if not isinstance(s, str) or not s.strip():
        return tuple([0]*6)
    s_clean = s.replace(" ", "").upper()
    for m in token_re.finditer(s_clean):
        raw_key, num = m.group(1).upper(), int(m.group(2))
        key = ALIASES.get(raw_key) or ALIASES.get(raw_key.replace("5","").replace("AC",""))
        if key:
            counts[key] += num
    return tuple(counts[k] for k in ORDER_T)

def autodetect_mass_column(df, prefer=None):
    cand = prefer or []
    cand += ["protonatedmass","observed_mass","precursor_mass","precursor m/z","m/z","mz","neutral_mass","Mass","mass"]
    for c in cand:
        if c in df.columns:
            s = pd.to_numeric(df[c], errors="coerce")
            if s.notna().any():
                return c
    for c in df.columns:
        cl = str(c).lower()
        if "mass" in cl or cl in ["mz","m/z"] or "precursor" in cl:
            s = pd.to_numeric(df[c], errors="coerce")
            if s.notna().any():
                return c
    return None

def read_insilico_map(insilico_csv, mass_col=None):
    df = pd.read_csv(insilico_csv)
    # resolve comp columns ignoring case
    cols = {c.lower(): c for c in df.columns}
    need = {}
    for k in ORDER_T:
        lk = k.lower()
        if lk in cols: need[k] = cols[lk]
    if "fuc" not in need and "fucose" in cols:
        need["Fuc"] = cols["fucose"]
    if "neuac" not in need and "neu5ac" in cols:
        need["NeuAc"] = cols["neu5ac"]
    if "neugc" not in need and "neu5gc" in cols:
        need["NeuGc"] = cols["neu5gc"]
    missing = [k for k in ORDER_T if k not in need]
    if missing:
        raise ValueError(f"In-silico CSV missing composition columns: {missing}")

    dfc = df[[need[k] for k in ORDER_T]].copy()
    dfc.columns = ORDER_T
    for c in ORDER_T:
        dfc[c] = pd.to_numeric(dfc[c], errors="coerce").fillna(0).astype(int)

    # resolve mass column
    mcol = mass_col or autodetect_mass_column(df)
    if mcol is None:
        # try to compute a representative mass by grouping existing numeric columns that look like mass
        num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        if len(num_cols) == 0:
            raise ValueError("No mass-like column found in in-silico CSV. Specify --insilico-mass-col.")
        # choose the numeric column with highest variance as a proxy
        mcol = max(num_cols, key=lambda c: float(df[c].var(skipna=True) or 0.0))

    masses = pd.to_numeric(df[mcol], errors="coerce")
    dfc["mass"] = masses
    # drop rows without mass
    dfc = dfc.dropna(subset=["mass"]).copy()
    dfc["tuple"] = dfc.apply(lambda r: (int(r["Hex"]),int(r["HexNAc"]),int(r["NeuAc"]),int(r["NeuGc"]),int(r["KDN"]),int(r["Fuc"])), axis=1)
    tupl_mass = dfc.groupby("tuple", as_index=False)["mass"].mean()  # average if duplicates
    return tupl_mass  # columns: tuple, mass

def ppm_delta(a, b):
    return (a - b) / b * 1e6

def mass_pair(manual_xlsx, pseudo_csv, insilico_csv=None, insilico_mass_col=None, ppm=PPM_DEFAULT, save_prefix=None):
    # Load manual
    xls = pd.ExcelFile(manual_xlsx)
    if "MSlist" not in xls.sheet_names:
        raise ValueError(f"'MSlist' sheet not found in {manual_xlsx}")
    manu = xls.parse("MSlist")
    if "Structure" not in manu.columns or "MS2scan_no" not in manu.columns:
        raise ValueError("Manual sheet must contain 'Structure' and 'MS2scan_no'.")
    m = manu[["Structure","MS2scan_no"]].dropna(subset=["MS2scan_no"]).copy()
    m["MS2scan_no"] = pd.to_numeric(m["MS2scan_no"], errors="coerce").astype("Int64")
    m["tuple"] = m["Structure"].astype(str).map(manual_to_tuple)

    # Manual mass: try to find a mass column
    man_mass_col = autodetect_mass_column(manu)
    if man_mass_col is not None and man_mass_col != "mass shift":
        m["man_mass"] = pd.to_numeric(manu.loc[m.index, man_mass_col], errors="coerce")
    else:
        m["man_mass"] = np.nan

    # If no manual mass, try derive from in-silico
    if m["man_mass"].isna().all():
        if insilico_csv is None:
            raise ValueError("No manual mass column found. Provide --insilico CSV so we can derive masses per composition.")
        t2m = read_insilico_map(insilico_csv, mass_col=insilico_mass_col)
        t2m = t2m.set_index("tuple")["mass"]
        m["man_mass"] = m["tuple"].map(t2m)

    # Load pseudo survivors
    p = pd.read_csv(pseudo_csv)
    # candidate pseudo mass columns
    pseudo_mass_col = autodetect_mass_column(p, prefer=["observed_mass","protonatedmass","theoretical_mass"])
    if pseudo_mass_col is None:
        raise ValueError("No mass-like column found in pseudolabel CSV (tried observed_mass/protonatedmass/theoretical_mass).")
    p["pseudo_mass"] = pd.to_numeric(p[pseudo_mass_col], errors="coerce")
    p["pseudo_tuple"] = p["composition"].astype(str).apply(lambda s: tuple(int(float(x)) for x in eval(s, {"__builtins__":{}})) if s else tuple([0]*6))

    # Pair each manual row to the nearest pseudo within ppm
    p_sub = p[["MS2scan_no","pseudo_mass","pseudo_tuple"]].dropna(subset=["pseudo_mass"]).copy()

    results = []
    for i, row in m.iterrows():
        man_mass = row["man_mass"]
        t = row["tuple"]
        if pd.isna(man_mass):
            results.append({**row.to_dict(), "best_pseudo_scan": pd.NA, "best_pseudo_mass": np.nan, "ppm_error": np.nan, "tuple_match": False, "mass_match": False})
            continue
        # Compute ppm to all pseudo
        deltas = (p_sub["pseudo_mass"] - man_mass) / man_mass * 1e6
        within = p_sub.loc[deltas.abs() <= ppm].copy()
        if within.empty:
            results.append({**row.to_dict(), "best_pseudo_scan": pd.NA, "best_pseudo_mass": np.nan, "ppm_error": np.nan, "tuple_match": False, "mass_match": False})
        else:
            j = within.iloc[deltas[within.index].abs().argmin()]
            ppm_err = float(((j["pseudo_mass"] - man_mass) / man_mass) * 1e6)
            tuple_match = (tuple(j["pseudo_tuple"]) == tuple(t))
            results.append({**row.to_dict(), "best_pseudo_scan": j["MS2scan_no"], "best_pseudo_mass": j["pseudo_mass"], "ppm_error": ppm_err, "tuple_match": tuple_match, "mass_match": True})

    out = pd.DataFrame(results)

    # Metrics from manual side
    mass_matched = out["mass_match"].fillna(False)
    tuple_matched = out["tuple_match"].fillna(False)
    n_manual = len(out)
    n_mass_matched = int(mass_matched.sum())
    n_tuple_matched = int(tuple_matched.sum())
    n_no_candidate = int((~mass_matched).sum())
    n_tuple_mismatch = int((mass_matched & (~tuple_matched)).sum())

    summary = {
        "ppm": ppm,
        "manual_rows": n_manual,
        "mass_matched_scans": n_mass_matched,
        "tuple_matched_scans": n_tuple_matched,
        "no_candidate_scans": n_no_candidate,
        "tuple_mismatch_scans": n_tuple_mismatch
    }

    if save_prefix:
        out_path = f"{save_prefix}_mass_pairing.csv"
        out.to_csv(out_path, index=False)
        print(f"[saved] {out_path}")

    # Print concise summary
    print(f"Manual rows: {n_manual}")
    print(f"Within ±{ppm} ppm mass match: {n_mass_matched}")
    print(f"Tuple matches (mass+composition): {n_tuple_matched}")
    print(f"No candidate within ppm: {n_no_candidate}")
    print(f"Mass match but tuple mismatch: {n_tuple_mismatch}")

    return out, summary

def main():
    ap = argparse.ArgumentParser(description="Mass-based pairing between manual MSlist and pseudolabel survivors using ±ppm tolerance.")
    ap.add_argument("--manual", required=True, help="Manual Excel with sheet MSlist (columns Structure, MS2scan_no).")
    ap.add_argument("--pseudo", required=True, help="Filtered pseudolabel CSV (e.g., score≥0.1 survivors).")
    ap.add_argument("--ppm", type=float, default=PPM_DEFAULT, help="Mass tolerance in ppm (default 20).")
    ap.add_argument("--insilico", help="Optional in-silico CSV to derive manual masses if manual has none.")
    ap.add_argument("--insilico-mass-col", help="Column name in in-silico CSV to use for mass (auto-detect if omitted).")
    ap.add_argument("--save-prefix", help="If set, save detailed pairing CSV with this prefix.")
    args = ap.parse_args()

    try:
        mass_pair(
            manual_xlsx=args.manual,
            pseudo_csv=args.pseudo,
            insilico_csv=args.insilico,
            insilico_mass_col=args.insilico_mass_col,
            ppm=args.ppm,
            save_prefix=args.save_prefix
        )
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
