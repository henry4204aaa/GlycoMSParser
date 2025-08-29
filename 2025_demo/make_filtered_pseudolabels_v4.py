#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, ast, math, os, re
from typing import Optional, Tuple, List, Dict

import numpy as np
import pandas as pd

PROTON_MASS = 1.007276466812

# ---------------- I/O helpers ----------------
def read_table_auto(path: str) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext in [".xls", ".xlsx"]:
        # caller can pass sheet_name separately where needed
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

#fix salvage scan count issue
def scan_to_str_series(s: pd.Series) -> pd.Series:
    # Normalize scan IDs so 17821, "17821", 17821.0 → "17821"
    s = s.astype(str).str.strip()
    s = s.str.replace(r"\.0$", "", regex=True)
    return s


def parse_shorthand_glycan(s: str) -> Optional[tuple]:
    # Shorthand F/H/N/S/G/KDN -> tuple order (H, N, S, G, KDN, F)
    if s is None or (isinstance(s, float) and np.isnan(s)): return None
    text = str(s).strip().upper().replace(" ", "")
    if not text: return None
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

# Canonical alias lists (extended)
ALIASES = {
    "scan": ["MS2scan_no","ms2scan_no","ms2_scan","ms2scan","scan","MS2"],
    "composition": ["composition","tuple","pseudolabel_tuple","pseudolabel_comp","predicted_comp","composition_tuple_str"],
    "ion_score": ["ion_score","ion score","ionscore","ion_similarity","score","ion_match_score"],
    # include theoretical mass variants here:
    "comp_mass": ["comp_neutral_mass","composition_neutral_mass","comp_mass","composition_mass","comp mass",
                  "theoretical_mass","theoretical mass"],
    "prec_neutral": ["precursor_neutral_mass","protonatedmass","precursor_mass","neutral_mass","neutral mass","protonated_mass","precursor_protonatedmass"],
    "prec_mz": ["precursor_mz","precursor mz","mz_precursor"],
    "charge": ["charge","z","precursor_charge"],
    "manual_tuple": [
        "manual_tuple","manual_comp","manual composition","manual_tuple_str",
        "manual","manual_comp_tuple","composition (manual)","composition_manual",
        "assigned_tuple","assigned_comp","annotation_tuple","annotation",
        "composition","composition_tuple","tuple"
    ]
}

def pick_precursor_neutral(row) -> Optional[float]:
    for k in ("precursor_neutral_mass","protonatedmass","precursor_mass","protonated_mass","precursor_protonatedmass"):
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

#fix override score
def scan_to_key(x):
    s = str(x).strip()
    if s.endswith(".0"):
        s = s[:-2]
    return s

# ---------------- overrides loader ----------------
def load_override_csv(path: str, score_min: float) -> pd.DataFrame:
    """
    Load a salvage/override CSV.

    Supported schemas:
    A) Old style: has a single 'composition' column like "(H,N,S,G,KDN,F)"
       + a composition mass column (any alias).
    B) Your style: separate monosaccharide columns:
         Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc, Mass
       -> we'll construct:
         - composition: "(H,N,S,G,KDN,F)"
         - composition_tuple: (H,N,S,G,KDN,F) as ints
         - comp_neutral_mass: from 'Mass'
    If 'ion_score' is absent, we fill it with score_min so the rows pass the filter.
    """
    df = read_table_auto(path)
    df.columns = [c.strip() for c in df.columns]
    cmap = build_colmap(df.columns)

    # --- Path B: monosaccharide columns present (case-insensitive) ---
    mono_keys = {
        "H": resolve_optional(cmap, ["Hex", "hex"]),
        "N": resolve_optional(cmap, ["HexNAc", "hexnac", "HexNac"]),
        "S": resolve_optional(cmap, ["NeuAc", "neuac", "sialic_acid", "SialicAcid"]),
        "G": resolve_optional(cmap, ["NeuGc", "neugc"]),
        "KDN": resolve_optional(cmap, ["KDN", "kdn"]),
        "F": resolve_optional(cmap, ["Fuc", "fuc", "Fucose", "fucose"]),
    }
    # mass aliases (your file uses 'Mass')
    mass_col = resolve_optional(cmap, ["Mass", "mass", "comp_neutral_mass", "composition_neutral_mass",
                                       "theoretical_mass", "theoretical mass", "comp_mass", "composition_mass", "comp mass"])
    have_mono = all(v is not None for v in mono_keys.values())

    if have_mono and mass_col is not None:
        # Build tuple order (H, N, S, G, KDN, F)
        H = pd.to_numeric(df[mono_keys["H"]], errors="coerce").fillna(0).astype(int)
        N = pd.to_numeric(df[mono_keys["N"]], errors="coerce").fillna(0).astype(int)
        S = pd.to_numeric(df[mono_keys["S"]], errors="coerce").fillna(0).astype(int)
        G = pd.to_numeric(df[mono_keys["G"]], errors="coerce").fillna(0).astype(int)
        K = pd.to_numeric(df[mono_keys["KDN"]], errors="coerce").fillna(0).astype(int)
        F = pd.to_numeric(df[mono_keys["F"]], errors="coerce").fillna(0).astype(int)

        comp_tuple_series = list(zip(H, N, S, G, K, F))
        comp_str_series = [f"({h}, {n}, {s}, {g}, {k}, {f})"
                           for h, n, s, g, k, f in comp_tuple_series]
        comp_mass_series = pd.to_numeric(df[mass_col], errors="coerce")

        out = pd.DataFrame({
            "composition": comp_str_series,
            "composition_tuple": comp_tuple_series,
            "comp_neutral_mass": comp_mass_series
        })

        # ion_score default
        ion_col = resolve_optional(cmap, ALIASES["ion_score"])
        if ion_col is None:
            out["ion_score"] = float(score_min)
        else:
            out["ion_score"] = pd.to_numeric(df[ion_col], errors="coerce").fillna(score_min)

        return out[["composition", "composition_tuple", "comp_neutral_mass", "ion_score"]].copy()

    # --- Path A: legacy single-column composition + mass aliases ---
    comp_col = resolve_optional(cmap, ALIASES["composition"])
    mass_col = resolve_optional(cmap, ALIASES["comp_mass"])

    if comp_col is None or mass_col is None:
        raise SystemExit(
            "[ERROR] override CSV not recognized.\n"
            "Expected either:\n"
            "  A) columns: composition + comp_neutral_mass (or aliases), or\n"
            "  B) columns: Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc, Mass\n"
            f"Found columns: {list(df.columns)}"
        )

    out = pd.DataFrame({
        "composition": df[comp_col],
        "composition_tuple": df[comp_col].apply(parse_tuple_like),
        "comp_neutral_mass": pd.to_numeric(df[mass_col], errors="coerce"),
    })
    ion_col = resolve_optional(cmap, ALIASES["ion_score"])
    if ion_col is None:
        out["ion_score"] = float(score_min)
    else:
        out["ion_score"] = pd.to_numeric(df[ion_col], errors="coerce").fillna(score_min)

    return out[["composition", "composition_tuple", "comp_neutral_mass", "ion_score"]].copy()

# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--manual", required=False)
    ap.add_argument("--manual-sheet", required=False, help="Excel sheet name for --manual (e.g., 'MSlist')")
    ap.add_argument("--manual-scan-col", required=False, help="Manual scan column override")
    ap.add_argument("--manual-tuple-col", required=False, help="Manual tuple column override (already 6-tuple)")
    ap.add_argument("--manual-shorthand-col", required=False, help="Manual shorthand column (e.g., F1H3N5KDN1)")
    ap.add_argument("--override-comps", required=False)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--ppm", type=float, default=20.0)
    ap.add_argument("--score-min", type=float, default=0.07)
    ap.add_argument("--big-mass", type=float, default=2500.0)
    ap.add_argument("--enable-mminus1", action="store_true")
    ap.add_argument("--enable-isotope-salvage", action="store_true",
                help="Enable isotope-aware candidate rescue at generation (adds M+1, optional M−1) "
                     "so compositions not matching strict precursor neutral mass can still enter. "
                     "Disabled by default.")
    ap.add_argument("--max-candidates-per-scan", type=int, default=20,
                help="Max candidates to keep per scan (ranked by ppm error). Default=20")
    ap.add_argument("--override-score-mode", choices=["floor", "lookup"], default="floor",
                help="How to score injected override rows. 'floor' sets ion_score=--score-min. "
                     "'lookup' tries to use the real per-scan candidate score if available, "
                     "and floors it to --score-min so it passes the filter.")
    args = ap.parse_args()

    # --- load candidates (CSV/TSV/XLSX ok) ---
    df = read_table_auto(args.candidates)
    df.columns = [c.strip() for c in df.columns]
    cmap = build_colmap(df.columns)

    scan_col  = resolve_or_fail(cmap, ALIASES["scan"], "candidates:scan")
    comp_col  = resolve_or_fail(cmap, ALIASES["composition"], "candidates:composition")
    score_col = resolve_or_fail(cmap, ALIASES["ion_score"], "candidates:ion_score")

    compmass_col = resolve_optional(cmap, ALIASES["comp_mass"])
    precN_col    = resolve_optional(cmap, ALIASES["prec_neutral"])
    precmz_col   = resolve_optional(cmap, ALIASES["prec_mz"])
    charge_col   = resolve_optional(cmap, ALIASES["charge"])

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

    df["_precursor_neutral"] = np.nan
    if precN_col is not None:
        df["_precursor_neutral"] = pd.to_numeric(df[precN_col], errors="coerce")
    else:
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

        # scan col
        if args.manual_scan_col:
            m_scan = args.manual_scan_col
            if m_scan not in m.columns:
                raise SystemExit(f"[ERROR] manual scan column '{m_scan}' not found in sheet '{args.manual_sheet or 'default'}'")
        else:
            m_scan = resolve_or_fail(m_map, ALIASES["scan"], "manual:scan")

        # composition source: tuple col → shorthand col → alias guess → autodetect shorthand
        manual_tuple_series = None
        if args.manual_tuple_col and args.manual_tuple_col in m.columns:
            manual_tuple_series = m[args.manual_tuple_col].apply(parse_tuple_like)
        elif args.manual_shorthand_col and args.manual_shorthand_col in m.columns:
            manual_tuple_series = m[args.manual_shorthand_col].apply(parse_shorthand_glycan)
        else:
            try:
                m_tuple_col = resolve_or_fail(m_map, ALIASES["manual_tuple"], "manual:manual_tuple")
                manual_tuple_series = m[m_tuple_col].apply(parse_tuple_like)
            except SystemExit:
                # autodetect shorthand
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

    # Capture manual scan set BEFORE any filtering (for coverage report)
    manual_scan_set = set(scan_to_str_series(df_manual["MS2scan_no"])) if df_manual is not None else set()


    # Capture raw-candidates scan set BEFORE gating (for coverage report)
    raw_candidate_scan_set = set(scan_to_str_series(df["_scan"]))

    # Cap-trimming counters
    cap_trimmed_main = 0
    cap_scans_main = 0
    # Cap-trimming counters for salvage
    cap_trimmed_over = 0
    cap_scans_over = 0
    # Consolidated totals (initialized to 0, will be updated later)
    total_trimmed = 0
    total_scans = 0

    # --- MASS GATING over candidates ---
    kept = []
    for _, r in df.iterrows():
        scan = r["_scan"]; comp = r["_composition"]; comp_tuple = r["_composition_tuple"]
        score = r["_ion_score"]; comp_mass = r["_comp_neutral"]; precN = r["_precursor_neutral"]

        if pd.isna(comp_mass) or pd.isna(precN):
            kept.append({
                "MS2scan_no": scan, "composition": comp, "composition_tuple": comp_tuple,
                "ion_score": score, "delta_mass_ppm": np.nan,
                "is_isotope_window": "NA_no_mass_info", "source": "input"
            })
            continue

        # by default: strict M only
        windows = [("M", precN)]

        # isotope-aware rescue (optional flag)
        if args.enable_isotope_salvage and precN >= args.big_mass:
            windows = isotope_neutral_variants(precN, args.enable_mminus1)

        best_ppm = None; best_lab = None
        matches = []
        for lab, target in windows:
            ppm = ppm_error(comp_mass, target)
            if ppm <= args.ppm:
                matches.append((ppm, lab))

        if matches:
            # sort by ppm error, pick top-K (default 20)
            matches.sort(key=lambda x: x[0])
            # cap report
            if len(matches) > args.max_candidates_per_scan:
                cap_trimmed_main += (len(matches) - args.max_candidates_per_scan)
                cap_scans_main += 1
            #
            for ppm, lab in matches[:args.max_candidates_per_scan]:
                kept.append({
                    "MS2scan_no": scan,
                    "composition": comp,
                    "composition_tuple": comp_tuple,
                    "ion_score": score,
                    "delta_mass_ppm": ppm,
                    "is_isotope_window": lab,
                    "source": "input"
                })

    df_kept = pd.DataFrame(kept)


    # Build a lookup of real scores per (scan, composition_tuple)
    # Use strings for scan keys to be robust to int/float/string mix
    df["_scan_key"] = df["_scan"].apply(scan_to_key)
    cand_score_map = (
        df[~df["_composition_tuple"].isna() & df["_ion_score"].notna()]
        .groupby(["_scan_key", "_composition_tuple"])["_ion_score"]
        .max()
        .to_dict()
    )
    # --- OVERRIDES (salvage) ---
    if args.override_comps:
        over = load_override_csv(args.override_comps, score_min=args.score_min)
        injected = []
        scans_prec = df[["_scan","_precursor_neutral"]].drop_duplicates()
         # override-score-mode=lookup stats
        n_over_score_found = 0
        n_over_score_lifted = 0
        n_over_score_missing = 0
        for _, rr in scans_prec.iterrows():
            scan = rr["_scan"]; precN = rr["_precursor_neutral"]
            if pd.isna(precN): continue
            windows = [("M", precN)]
            if precN >= args.big_mass:
                windows = isotope_neutral_variants(precN, args.enable_mminus1)
            for _, ov in over.iterrows():
                comp = ov["composition"]
                comp_tuple = ov["composition_tuple"]
                comp_mass = ov["comp_neutral_mass"]
                # default: keep floor (from loader) unless we switch to lookup
                ov_score = ov["ion_score"]
                score_origin = "floor"
                if args.override_score_mode == "lookup":
                    scan_key = scan_to_key(scan)
                    found = cand_score_map.get((scan_key, comp_tuple))
                    scan_best = (
                                df[df["_scan_key"] == scan_key]["_ion_score"].max()
                                if (df is not None) else None
                            )
                    if found is not None:
                        score_origin = "lookup_hit"
                        if float(found) < float(args.score_min):
                            n_over_score_lifted += 1
                            score_origin = "lookup_lifted"
                        else:
                            n_over_score_found += 1
                        ov_score = max(float(found), float(args.score_min))
                    else:
                        n_over_score_missing += 1
                        score_origin = "lookup_missing"
                best_ppm = None; best_lab = None
                matches = []
                for lab, target in windows:
                    ppm = ppm_error(comp_mass, target)
                    if ppm <= args.ppm:
                        matches.append((ppm, lab))
                if matches:
                    matches.sort(key=lambda x: x[0])
                    # cap report
                    if len(matches) > args.max_candidates_per_scan:
                        cap_trimmed_over += (len(matches) - args.max_candidates_per_scan)
                        cap_scans_over += 1
                    #
                    for ppm, lab in matches[:args.max_candidates_per_scan]:
                        injected.append({
                            "MS2scan_no": int(scan),
                            "composition": comp,
                            "composition_tuple": comp_tuple,
                            "ion_score": ov_score,
                            "delta_mass_ppm": ppm,
                            "is_isotope_window": lab,
                            "source": "override",
                            "override_score_origin": score_origin,   # <-- add this
                            "override_real_score": found,      # None if missing
                            "scan_best_score": scan_best           # for display only
                        })
        if injected:
            df_kept = pd.concat([df_kept, pd.DataFrame(injected)], ignore_index=True)
            print(f"[INFO] Injected {len(injected)} override rows (score ≥ {args.score_min}).")
            #failsafe of injected scan no
            df_kept["MS2scan_no"] = pd.to_numeric(df_kept["MS2scan_no"], errors="coerce").astype("Int64")
                # report how override scores were assigned in lookup mode
            if args.override_score_mode == "lookup":
                print(f"[INFO] Override scoring (lookup): used real scores for {n_over_score_found} rows, "
                    f"lifted {n_over_score_lifted} below floor to {args.score_min}, "
                    f"no-match fallback for {n_over_score_missing}.")

            total_trimmed = cap_trimmed_main + cap_trimmed_over
            # Cap-trimming counters
            total_scans = cap_scans_main + cap_scans_over
            if total_trimmed > 0:
                print(f"[INFO] Cap trimmed {total_trimmed} candidates across {total_scans} scans (top-{args.max_candidates_per_scan} kept).")
            else:
                print(f"[INFO] Cap did not trim any candidates (top-{args.max_candidates_per_scan} kept).")
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

    # ---------- WRITE DETAILED ----------
    detailed_out = f"{args.out_prefix}_detailed.csv"
    df_kept.to_csv(detailed_out, index=False)
    print(f"[OK] Wrote detailed rows → {detailed_out} (rows={len(df_kept)})")

    # ---------- SUMMARY ----------

    if df_kept.empty:
        summary = pd.DataFrame(columns=["MS2scan_no","n_candidates","any_match","manual_tuple","best_candidate","best_score","best_is_isotope_window"])
        #failsafe of scanno
        summary["MS2scan_no"] = pd.to_numeric(summary["MS2scan_no"], errors="coerce").astype("Int64")
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
    #failsafe of scanno
    summary["MS2scan_no"] = pd.to_numeric(summary["MS2scan_no"], errors="coerce").astype("Int64")
    summary.to_csv(summary_out, index=False)
    print(f"[OK] Wrote per-scan summary → {summary_out} (rows={len(summary)})")

    # ---------- COVERAGE REPORT (for manuscript) ----------
    # After-filters scan set:
    after_filters_scan_set = set(scan_to_str_series(summary["MS2scan_no"]))

    # Manual coverage vs raw candidates
    if manual_scan_set:
        n_manual_total = len(manual_scan_set)
        n_in_raw = len(manual_scan_set & raw_candidate_scan_set)
        n_in_after = len(manual_scan_set & after_filters_scan_set)

        # Missing lists
        missing_in_raw = sorted(list(manual_scan_set - raw_candidate_scan_set), key=lambda x: (len(x), x))
        missing_after = sorted(list(manual_scan_set - after_filters_scan_set), key=lambda x: (len(x), x))

        # Print concise report
        print("[COVERAGE] Manual scans total:", n_manual_total)
        print("[COVERAGE] Manual scans present in RAW candidates:", n_in_raw)
        print("[COVERAGE] Manual scans present AFTER FILTERS:", n_in_after)

        # Write CSVs for your records
        cov_dir_txt = f"{args.out_prefix}_manual_coverage_report.txt"
        with open(cov_dir_txt, "w", encoding="utf-8") as f:
            # Core coverage
            f.write(f"Manual scans total: {n_manual_total}\n")
            f.write(f"Manual scans in RAW candidates: {n_in_raw}\n")
            f.write(f"Manual scans AFTER FILTERS: {n_in_after}\n")

            # Filter settings (useful for reproducibility)
            f.write("\n[Settings]\n")
            f.write(f"ppm_window: {args.ppm}\n")
            f.write(f"score_min: {args.score_min}\n")
            f.write(f"big_mass_threshold: {args.big_mass}\n")
            f.write(f"isotope_salvage_enabled: {bool(args.enable_isotope_salvage)}\n")
            f.write(f"max_candidates_per_scan: {args.max_candidates_per_scan}\n")

            # Cap-trim stats
            f.write("\n[Cap Trim]\n")
            f.write(f"trimmed_main_candidates: {cap_trimmed_main}  (affected_scans={cap_scans_main})\n")
            f.write(f"trimmed_override_candidates: {cap_trimmed_over}  (affected_scans={cap_scans_over})\n")
            f.write(f"trimmed_total: {total_trimmed}  (affected_scans_total={total_scans})\n")
        #Also write a tiny CSV for cap stats (optional but handy)
        cap_csv = f"{args.out_prefix}_cap_stats.csv"
        pd.DataFrame([{
            "ppm_window": args.ppm,
            "score_min": args.score_min,
            "big_mass_threshold": args.big_mass,
            "isotope_salvage_enabled": bool(args.enable_isotope_salvage),
            "max_candidates_per_scan": args.max_candidates_per_scan,
            "trimmed_main_candidates": cap_trimmed_main,
            "trimmed_override_candidates": cap_trimmed_over,
            "trimmed_total": total_trimmed,
            "affected_scans_main": cap_scans_main,
            "affected_scans_override": cap_scans_over,
            "affected_scans_total": total_scans,
        }]).to_csv(cap_csv, index=False)
        print(f"[OK] Wrote cap stats → {cap_csv}")

        pd.DataFrame({"MS2scan_no_missing_in_raw": missing_in_raw}).to_csv(f"{args.out_prefix}_missing_in_raw.csv", index=False)
        pd.DataFrame({"MS2scan_no_missing_after_filters": missing_after}).to_csv(f"{args.out_prefix}_missing_after_filters.csv", index=False)
        print(f"[OK] Wrote coverage report → {cov_dir_txt}")
        print(f"[OK] Wrote missing scan lists → {args.out_prefix}_missing_in_raw.csv, {args.out_prefix}_missing_after_filters.csv")

    #if "tuple_match" in df_kept.columns and not df_kept.empty:
    #    print("[INFO] Coverage (after filters):", df_kept["tuple_match"].value_counts(dropna=False).to_dict())
    # ---------- EXTRA BREAKDOWN ----------
    if "tuple_match" in df_kept.columns:
        # manual scans that survived filters
        scans_after = set(after_filters_scan_set & manual_scan_set)

        # how many of those retained scans actually matched the manual tuple?
        retained_match = df_kept[(df_kept["MS2scan_no"].astype(str).isin(scans_after)) &
                                    (df_kept["tuple_match"] == True)]["MS2scan_no"].nunique()

        retained_total = len(scans_after)
        retained_mismatch = retained_total - retained_match
        lost_by_score = len(manual_scan_set) - retained_total

        print("[BREAKDOWN] Manual scans total:", len(manual_scan_set))
        print("[BREAKDOWN] Retained & matched:", retained_match)
        print("[BREAKDOWN] Retained but mismatched:", retained_mismatch)
        print("[BREAKDOWN] Lost by score/mass filters:", lost_by_score)

        # Write 1-row breakdown CSV
        breakdown_out = f"{args.out_prefix}_manual_breakdown.csv"
        pd.DataFrame([{
            "manual_total": len(manual_scan_set),
            "retained_match": retained_match,
            "retained_mismatch": retained_mismatch,
            "lost_by_filters": lost_by_score
        }]).to_csv(breakdown_out, index=False)
        print(f"[OK] Wrote breakdown → {breakdown_out}")

        # ---------- LIST THE RETAINED-BUT-MISMATCHED SCANS ----------
        # Use 'summary' (one row per scan) to pick scans that:
        #   - are in manual,
        #   - survived filters,
        #   - but did not match the manual composition
        summ2 = summary.copy()
        summ2["MS2scan_no_str"] = summ2["MS2scan_no"].astype(str)

        mismatched_scans = summ2[
            (summ2["MS2scan_no_str"].isin(manual_scan_set)) &
            (summ2["any_match"] == False)        # explicitly False, not NaN
        ].copy()

        # Add "top-3 candidates" from df_kept for convenience
        topk_rows = []
        kept_sorted = df_kept.sort_values(["MS2scan_no", "ion_score"], ascending=[True, False])

        mismatch_scan_set = set(mismatched_scans["MS2scan_no"])  # keep as native dtype (int/str)
        for scan, sub in kept_sorted.groupby("MS2scan_no"):
            if scan not in mismatch_scan_set:
                continue
            top3 = sub[["composition", "ion_score"]].head(3)
            top3_str = "; ".join([f"{r['composition']} (s={r['ion_score']:.4f})" for _, r in top3.iterrows()])
            topk_rows.append({"MS2scan_no": scan, "top3_candidates": top3_str})

        # Ensure topk_df always has the expected columns
        if topk_rows:
            topk_df = pd.DataFrame(topk_rows, columns=["MS2scan_no", "top3_candidates"])
        else:
            topk_df = pd.DataFrame({"MS2scan_no": pd.Series(dtype=mismatched_scans["MS2scan_no"].dtype),
                                    "top3_candidates": pd.Series(dtype="object")})

        # Merge top-3 into the mismatched table; keep helpful columns
        mismatched_scans = mismatched_scans.merge(topk_df, on="MS2scan_no", how="left")
        mismatched_scans = mismatched_scans[[
            "MS2scan_no",
            "manual_tuple",
            "best_candidate",
            "best_score",
            "best_is_isotope_window",
            "n_candidates",
            "top3_candidates"
        ]].sort_values("MS2scan_no")

        mism_out = f"{args.out_prefix}_retained_mismatched_scans.csv"
        mismatched_scans.to_csv(mism_out, index=False)
        print(f"[OK] Wrote retained-but-mismatched list → {mism_out} (rows={len(mismatched_scans)})")

    """
    # ---------- EXTRA BREAKDOWN ----------
    if "tuple_match" in df_kept.columns:
        # manual scans after filters
        scans_after = set(after_filters_scan_set & manual_scan_set)

        retained_match = df_kept[(df_kept["MS2scan_no"].astype(str).isin(scans_after)) &
                                    (df_kept["tuple_match"] == True)]["MS2scan_no"].nunique()

        retained_total = len(scans_after)
        retained_mismatch = retained_total - retained_match
        lost_by_score = len(manual_scan_set) - retained_total

        print("[BREAKDOWN] Manual scans total:", len(manual_scan_set))
        print("[BREAKDOWN] Retained & matched:", retained_match)
        print("[BREAKDOWN] Retained but mismatched:", retained_mismatch)
        print("[BREAKDOWN] Lost by score/mass filters:", lost_by_score)

        # Write breakdown to CSV
        breakdown_out = f"{args.out_prefix}_manual_breakdown.csv"
        pd.DataFrame([{
            "manual_total": len(manual_scan_set),
            "retained_match": retained_match,
            "retained_mismatch": retained_mismatch,
            "lost_by_filters": lost_by_score
        }]).to_csv(breakdown_out, index=False)
        print(f"[OK] Wrote breakdown → {breakdown_out}")
    """

if __name__ == "__main__":
    main()