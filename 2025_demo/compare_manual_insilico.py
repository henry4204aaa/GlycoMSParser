import argparse
import sys
import pandas as pd
import re
from typing import Tuple, Dict, List, Set

#for examining if the in silico csv matches the manual annotation if they have


# Canonical order used in outputs and matching
ORDER = ["Hex", "HexNAc", "NeuAc", "NeuGc", "KDN", "Fuc"]

# Aliases accepted in manual "Structure" strings
ALIASES = {
    # Hexose
    "H": "Hex", "HEX": "Hex",
    # HexNAc
    "N": "HexNAc", "HEXNAC": "HexNAc",
    # NeuAc (sialic acid)
    "S": "NeuAc", "NEU5AC": "NeuAc", "NEUAC": "NeuAc", "SA": "NeuAc",
    # NeuGc
    "G": "NeuGc", "NEU5GC": "NeuGc", "NEUGC": "NeuGc",
    # KDN
    "KDN": "KDN", "K": "KDN",
    # Fucose
    "F": "Fuc", "FUC": "Fuc", "FUCOSE": "Fuc"
}

token_re = re.compile(r"([A-Za-z]+)\s*([+-]?\d+)", re.IGNORECASE)

def _safe_int(x):
    try:
        return int(x)
    except Exception:
        try:
            return int(float(x))
        except Exception:
            return 0

def _tuple_to_dict_key(t: Tuple[int, int, int, int, int, int]) -> Dict[str, int]:
    return {k: t[i] for i, k in enumerate(ORDER)}

def _dict_to_tuple(d: Dict[str, int]) -> Tuple[int, int, int, int, int, int]:
    return tuple(int(d.get(k, 0)) for k in ORDER)

def _tuple_to_compact_string(t: Tuple[int, int, int, int, int, int]) -> str:
    # Compact like F2H6N5S1G0KDN0 (order: F,H,N,S,G,KDN) to be human-friendly in logs
    labels = [("F", t[5]), ("H", t[0]), ("N", t[1]), ("S", t[2]), ("G", t[3]), ("KDN", t[4])]
    return "".join([f"{lab}{val}" for lab, val in labels if val != 0] + ([f"KDN0"] if t[4]==0 else []))

def parse_manual_structure_to_tuple(s: str) -> Tuple[int, int, int, int, int, int]:
    """
    Parse manual shorthand like 'F2H6N5S1' into tuple (Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc).
    Accepts aliases defined in ALIASES.
    Returns a 6-int tuple.
    """
    counts = {k: 0 for k in ORDER}
    if not isinstance(s, str) or not s.strip():
        return _dict_to_tuple(counts)

    s_clean = s.replace(" ", "").upper()
    for m in token_re.finditer(s_clean):
        raw_key, num = m.group(1).upper(), int(m.group(2))
        # normalize aliases
        key = ALIASES.get(raw_key)
        if key is None:
            # lenient collapse (remove numerals & AC to catch NEU5AC → NEUAC, etc.)
            raw_try = raw_key.replace("5", "").replace("AC", "")
            key = ALIASES.get(raw_try)
        if key is None:
            # Unrecognized token, ignore silently
            continue
        counts[key] = counts.get(key, 0) + num

    return _dict_to_tuple(counts)

def read_insilico_compositions(csv_path: str) -> pd.DataFrame:
    """
    Read in-silico CSV and return unique compositions as tuples in ORDER.
    Accepts column variants for Fucose (Fuc/Fucose).
    """
    df = pd.read_csv(csv_path)
    cols = {c.lower(): c for c in df.columns}
    # Resolve expected columns (case-insensitive)
    needed = {}
    for name in ORDER:
        key = name.lower()
        if key in cols:
            needed[name] = cols[key]
        else:
            # special-case Fuc/Fucose
            if name == "Fuc":
                if "fuc" in cols: needed[name] = cols["fuc"]
                elif "fucose" in cols: needed[name] = cols["fucose"]
            # Accept NeuAc/Neu5Ac as NeuAc, etc.
            elif name == "NeuAc":
                if "neuac" in cols: needed[name] = cols["neuac"]
                elif "neu5ac" in cols: needed[name] = cols["neu5ac"]
            elif name == "NeuGc":
                if "neugc" in cols: needed[name] = cols["neugc"]
                elif "neu5gc" in cols: needed[name] = cols["neu5gc"]
            elif name == "HexNAc":
                if "hexnac" in cols: needed[name] = cols["hexnac"]
    missing = [k for k in ORDER if k not in needed]
    if missing:
        raise ValueError(f"In-silico CSV is missing required columns: {missing}. Found columns={list(df.columns)}")

    # Cast to int and deduplicate
    comp_df = df[[needed[k] for k in ORDER]].copy()
    comp_df.columns = ORDER
    for c in ORDER:
        comp_df[c] = comp_df[c].map(_safe_int)
    comp_df["tuple"] = comp_df.apply(lambda r: tuple(int(r[k]) for k in ORDER), axis=1)
    comp_unique = comp_df.drop_duplicates(subset=["tuple"]).reset_index(drop=True)
    return comp_unique

def read_manual_compositions(excel_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read manual Excel (sheet 'MSlist'), pull ['Structure','MS2scan_no'],
    parse Structure to tuple (Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc).
    Returns:
      - manual_scans: each row is a scan with its tuple (may duplicate tuples).
      - manual_unique: unique tuples with a 'count_scans' column.
    """
    xls = pd.ExcelFile(excel_path)
    if "MSlist" not in xls.sheet_names:
        raise ValueError(f"Sheet 'MSlist' not found. Sheets={xls.sheet_names}")
    df = xls.parse("MSlist")

    # Ensure columns exist (case-insensitive resolve)
    cols = {c.lower(): c for c in df.columns}
    if "structure" not in cols or "ms2scan_no" not in cols:
        raise ValueError(f"'MSlist' must contain columns ['Structure','MS2scan_no']. Found={list(df.columns)}")

    df2 = df[[cols["structure"], cols["ms2scan_no"]]].rename(columns={cols["structure"]: "Structure", cols["ms2scan_no"]: "MS2scan_no"})
    # Drop scans without scan id
    df2 = df2.dropna(subset=["MS2scan_no"]).copy()
    df2["MS2scan_no"] = pd.to_numeric(df2["MS2scan_no"], errors="coerce").astype("Int64")
    df2 = df2.dropna(subset=["MS2scan_no"]).copy()

    # Parse to tuples
    df2["tuple"] = df2["Structure"].astype(str).map(parse_manual_structure_to_tuple)

    # Unique composition list with scan counts
    manu_unique = df2.groupby("tuple", dropna=False).size().reset_index(name="count_scans")
    return df2, manu_unique

def comparing_annotation_against_insilico(insilicocsv: str, manannotationexcel: str):
    """
    Compare composition coverage:
      - Overlap compositions
      - Manual-only compositions
      - In-silico-only compositions
    Prints a concise summary and returns a dict with dataframes.
    """
    # Load inputs
    insi = read_insilico_compositions(insilicocsv)
    manu_scans, manu_unique = read_manual_compositions(manannotationexcel)

    insi_set: Set[Tuple[int, int, int, int, int, int]] = set(insi["tuple"].tolist())
    manu_set: Set[Tuple[int, int, int, int, int, int]] = set(manu_unique["tuple"].tolist())

    overlap = sorted(list(insi_set & manu_set))
    manu_only = sorted(list(manu_set - insi_set))
    insi_only = sorted(list(insi_set - manu_set))

    # Metrics
    n_insi = len(insi_set)
    n_manu_comp = len(manu_set)
    n_manu_scans = int(manu_scans.shape[0])
    n_overlap = len(overlap)
    n_manu_only = len(manu_only)
    n_insi_only = len(insi_only)

    # Pretty print
    def pretty_list(tuples: List[Tuple[int,int,int,int,int,int]], limit=50) -> str:
        if not tuples:
            return "[]"
        items = [ _tuple_to_compact_string(t) for t in tuples[:limit] ]
        suffix = "" if len(tuples) <= limit else f" ... (+{len(tuples)-limit} more)"
        return "[" + ", ".join(items) + "]" + suffix

    print(f"in silico csv: {n_insi} compositions")
    print(f"manual annotation: {n_manu_comp} compositions, {n_manu_scans} MS2scan_no")
    print(f"manual annotation found in in silico csv:  {n_overlap}/{n_manu_comp}")
    print(f"missing manual annotation: {pretty_list(manu_only)}")
    print(f"in silico only (not in manual): {pretty_list(insi_only)}")

    # Build dataframes for programmatic use
    df_overlap = pd.DataFrame(overlap, columns=ORDER).assign(source="both")
    df_manu_only = pd.DataFrame(manu_only, columns=ORDER).assign(source="manual_only")
    df_insi_only = pd.DataFrame(insi_only, columns=ORDER).assign(source="insilico_only")

    return {
        "overlap": df_overlap,
        "manual_only": df_manu_only,
        "insilico_only": df_insi_only,
        "metrics": {
            "insilico_compositions": n_insi,
            "manual_compositions": n_manu_comp,
            "manual_ms2scan_count": n_manu_scans,
            "overlap_count": n_overlap,
            "manual_only_count": n_manu_only,
            "insilico_only_count": n_insi_only,
        }
    }

def _save_optional_tables(dfs: Dict[str, pd.DataFrame], prefix: str):
    for key, df in dfs.items():
        if isinstance(df, pd.DataFrame):
            df.to_csv(f"{prefix}_{key}.csv", index=False)

def main():
    ap = argparse.ArgumentParser(description="Compare composition coverage between in-silico CSV and manual annotation Excel (MSlist).")
    ap.add_argument("--insilico", required=True, help="Path to in-silico CSV with columns: Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc (case-insensitive).")
    ap.add_argument("--manual", required=True, help="Path to manual annotation Excel containing sheet 'MSlist' with columns: Structure, MS2scan_no.")
    ap.add_argument("--save-prefix", default=None, help="If provided, save overlap/manual_only/insilico_only tables to CSV with this prefix.")
    args = ap.parse_args()

    try:
        results = comparing_annotation_against_insilico(args.insilico, args.manual)
        if args.save_prefix:
            _save_optional_tables(results, args.save_prefix)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

"""
# Basic summary to stdout
python compare_insilico_vs_manual.py \
  --insilico "G:\其他電腦\My Computer\GlycoMSParser\zf_sPerMeNG_intestine_insilico_20250824.csv.csv" \
  --manual "G:\其他電腦\My Computer\GlycoMSParser\src\20240922_temp_zf_intestine_1.xlsx"

# Also save tables to CSVs (prefix creates *_overlap.csv, *_manual_only.csv, *_insilico_only.csv)
python compare_insilico_vs_manual.py \
  --insilico /path/to/insilico.csv \
  --manual /path/to/20240922_temp_zf_intestine_1.xlsx \
  --save-prefix zf_intestine_compcheck

  """