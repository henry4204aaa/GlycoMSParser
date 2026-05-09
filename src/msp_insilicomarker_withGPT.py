from __future__ import annotations
#Python3.9 compatibility fix (not working haha)
#part of it was copied from GPT5 answer based on the discussion
version = "1.09"
last_update = 20260411
import numpy as np
import csv
import mspvalidator_merger as validator
import pandas as pd
import re
import json



#changelog
#v1.1 (future) fix existing know bugs if exist. Remove unneeded functions and imports
#v1.09 code review 1, clean up
#v1.08 support negative mode, HexA and Sulphate, Phosphate labeling
#v1.0 works in positive mode labeling



#to get rid of warning
def ppm(x): return x



#extract_ionmasslist in mspvalidator
def extract_ionmasslist(ionmass_sheet):
    ion_df = ionmass_sheet[["mass"]]
    iondfindex = ion_df.values.flatten().tolist()
    return iondfindex

#sort mass for later searching
def parse_composition(compfile):
    #only for single func test#
    #compfile1 = pd.read_csv(compfile)
    ##
    #print(compfile1.head())
    print("sort composition list by mass ascend")
    sort = compfile.sort_values('Mass').reset_index(drop=True)
    print(sort.head())
    return sort

###GPT generated blocks###
import ast

#202507 version of in_silico file
def normalize_insilico(
    df: pd.DataFrame,
    comp_cols=("Hex", "HexNAc", "NeuAc", "NeuGc", "KDN", "Fuc"),
    mass_col="Mass",
    add_legacy_repr: bool = True,
):
    """
    Ensure numeric dtypes, add a tuple + string representation of composition.
    Returns a copy with columns: comp_cols..., Mass, comp_tuple, comp_str
    """
    out = df.copy()

    # force numeric
    for c in comp_cols:
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("Int64")
    out[mass_col] = pd.to_numeric(out[mass_col], errors="coerce").astype(float)

    if add_legacy_repr:
        # tuple in the legacy order (adjust comp_cols order if your legacy tuple order differs)
        out["comp_tuple"] = list(out[list(comp_cols)].itertuples(index=False, name=None))
        out["comp_str"]   = out["comp_tuple"].apply(lambda t: f"({', '.join(str(int(x)) for x in t)})")

    # a clean sorted view (keep the unsorted original around if you want)
    sorted_out = out.sort_values(mass_col).reset_index(drop=True)
    return sorted_out

#20250929 negative support
# ===== Canonical label builders (HexA / SO3 / PO3H aware) =====

def _fmt_part(prefix: str, n: int) -> str:
    """Return 'prefixN' only if n>0."""
    try:
        n = int(n)
    except Exception:
        n = 0
    return f"{prefix}{n}" if n > 0 else ""

def canonical_label_tuple(comp6, *, hexA=0, so3=0, po3h=0, style="short") -> str:
    """
    Build canonical label from a 6-tuple (H,N,Ac,Gc,KDN,F) plus modifiers.

    style="short": uses 'A' (HexA), 's' (SO3), 'p' (PO3H)
    style="long" : uses 'A' (HexA), 'Sul' (SO3), 'Phos' (PO3H)

    Base order agreed: F, H, N, S(NeuAc), G(NeuGc), KDN, then A/s/p (or A/Sul/Phos).
    """
    H, N, Ac, Gc, KDN, F = map(int, comp6)

    parts = []
    parts.append(_fmt_part("F", F))
    parts.append(_fmt_part("H", H))
    parts.append(_fmt_part("N", N))
    parts.append(_fmt_part("S", Ac))   # NeuAc
    parts.append(_fmt_part("G", Gc))   # NeuGc
    parts.append(_fmt_part("KDN", KDN))

    if style == "short":
        parts.append(_fmt_part("A", hexA))   # HexA
        parts.append(_fmt_part("s", so3))    # sulfate (SO3)
        parts.append(_fmt_part("p", po3h))   # phosphate monoester (PO3H)
    else:
        parts.append(_fmt_part("A", hexA))
        parts.append(_fmt_part("Sul", so3))
        parts.append(_fmt_part("Phos", po3h))

    return "".join(p for p in parts if p)

def canonical_label_from_row(row, style="short") -> str:
    """
    Build canonical label from a pandas row/dict that may include:
      Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc, HexA, SO3, PO3H
    Missing columns are treated as 0.
    """
    H   = int(row.get("Hex", 0))
    N   = int(row.get("HexNAc", 0))
    Ac  = int(row.get("NeuAc", 0))
    Gc  = int(row.get("NeuGc", 0))
    K   = int(row.get("KDN", 0))
    F   = int(row.get("Fuc", 0))
    A   = int(row.get("HexA", 0))
    s   = int(row.get("SO3", 0))
    p   = int(row.get("PO3H", 0))
    return canonical_label_tuple((H, N, Ac, Gc, K, F), hexA=A, so3=s, po3h=p, style=style)



#legacy, or 2023 version of in silico
def legacy_two_col_to_new(path, comp_cols=("Hex","HexNAc","NeuAc","NeuGc","KDN","Fuc")):
    raw = pd.read_csv(path, header=None, names=["comp_str","Mass"])
    # comp_str looks like: "(0, 3, 2, 0, 0)"
    tup = raw["comp_str"].str.strip().apply(ast.literal_eval)
    comps = pd.DataFrame(tup.tolist(), columns=comp_cols)
    df = pd.concat([comps, raw["Mass"]], axis=1)
    return normalize_insilico(df, comp_cols=comp_cols, mass_col="Mass")

def _apply_mass_transform(theo_array, mass_transform):
    """
    Apply an offset / transform so library masses live in the same domain
    as your observed protonated masses. Accepts:
      - None          -> no change
      - "M+H","M+Na","M+K","M+NH4","neutral"
      - numeric       -> add that offset
      - callable      -> f(theo_array) -> array
      - "offset:<x>"  -> add <x>
    """
    PROTON = 1.007276466812
    OFFSETS = {
        "neutral": 0.0,
        "M+H": PROTON, "[M+H]+": PROTON,
        "M+Na": 22.989218, "[M+Na]+": 22.989218,
        "M+K": 38.963158, "[M+K]+": 38.963158,
        "M+NH4": 18.033823, "[M+NH4]+": 18.033823,
    }
    arr = np.asarray(theo_array, dtype=float)
    if mass_transform is None:
        return arr
    if callable(mass_transform):
        return np.asarray(mass_transform(arr), dtype=float)
    if isinstance(mass_transform, (int, float)):
        return arr + float(mass_transform)
    key = str(mass_transform).strip()
    if key in OFFSETS:
        return arr + OFFSETS[key]
    if key.startswith("offset:"):
        return arr + float(key.split(":", 1)[1])
    raise ValueError(f"Unsupported mass_transform: {mass_transform!r}")

def find_compositions_for_mass(
    obs_mass: float,
    insilico_sorted: pd.DataFrame,
    mass_col: str = "Mass",
    ppm: float = 20.0,
    mass_transform=None,  # e.g., "M+H", "M+Na", a numeric offset, or a callable
    comp_cols=("Hex", "HexNAc", "NeuAc", "NeuGc", "KDN", "Fuc"),
):
    """
    Return ALL in-silico rows within ±ppm of obs_mass, with deltas.
    mass_transform lets you align the library mass domain to your observed domain:
      - "M+H", "M+Na", "M+K", "M+NH4"
      - a numeric offset (e.g., 1.007276)
      - a callable: lambda theo: theo + <something>
    """
    # 1) Raw library masses → transformed masses
    #theo_raw = insilico_df[mass_col].to_numpy(dtype=float)
    theo_raw = insilico_sorted[mass_col].to_numpy(dtype=float)
    theo_x   = _apply_mass_transform(theo_raw, mass_transform)  # <-- THIS is theo_x

    # 2) Sort BY the transformed masses (needed for searchsorted)
    order        = np.argsort(theo_x)
    theo_sorted  = theo_x[order]
    tol = abs(obs_mass) * float(ppm) / 1e6
    lo = np.searchsorted(theo_x, obs_mass - tol, side="left")
    hi = np.searchsorted(theo_x, obs_mass + tol, side="right")

    if lo >= hi:
        # ✅ ensure all columns exist even when empty
        empty = insilico_sorted.iloc[0:0].copy()
        for c in list(comp_cols):
            if c not in empty.columns:
                empty[c] = pd.Series(dtype="Int64")
        return empty.assign(
            observed_mass=obs_mass,
            theoretical_mass=np.nan,
            delta_Da=np.nan,
            ppm_error=np.nan,
            comp_tuple=pd.Series(dtype=object),
            comp_str=pd.Series(dtype=object),
        )

    hits = insilico_sorted.iloc[lo:hi].copy()
    t = theo_x[lo:hi]
    delta = obs_mass - t
    hits["observed_mass"]   = obs_mass
    hits["theoretical_mass"] = t
    hits["delta_Da"]        = delta
    hits["ppm_error"]       = (delta / t) * 1e6

    if "comp_tuple" not in hits.columns:
        hits["comp_tuple"] = list(hits[list(comp_cols)].itertuples(index=False, name=None))
        hits["comp_str"]   = hits["comp_tuple"].apply(lambda u: f"({', '.join(str(int(x)) for x in u)})")

    front = ["observed_mass", "theoretical_mass", "delta_Da", "ppm_error", "comp_str"]
    rest  = [c for c in hits.columns if c not in front]
    return hits[front + rest]


def hits_to_legacy_list(hits: pd.DataFrame):
    # ✅ guard for empty & flexible mass column fallback
    if hits is None or hits.empty:
        return []
    tm_col = "theoretical_mass" if "theoretical_mass" in hits.columns else (
             "Mass" if "Mass" in hits.columns else None)
    if tm_col is None:
        # last-ditch: don’t crash; emit comp only
        return [{"comp": tup} for tup in hits.get("comp_tuple", [])]

    comps = hits["comp_tuple"] if "comp_tuple" in hits.columns else None
    if comps is None:
        # build on the fly if you didn’t keep comp_tuple
        comp_cols = [c for c in ("Hex","HexNAc","NeuAc","NeuGc","KDN","Fuc") if c in hits.columns]
        comps = [tuple(int(x) for x in row) for row in hits[comp_cols].to_numpy()] if comp_cols else [None]*len(hits)

    return [
        {"comp": tup, "mass": float(tm), "ppm_error": float(ppm)}
        for tup, tm, ppm in zip(comps, hits[tm_col], hits.get("ppm_error", [np.nan]*len(hits)))
    ]



# 1) For each row’s observed mass, collect all comp matches within ppm
def collect_row_hits(row):
    x = float(row["protonatedmass"])
    hits = find_compositions_for_mass(x, lib, ppm=20.0, mass_transform="M+H")
    return hits_to_legacy_list(hits)  # or return the DataFrame itself if you prefer




# or, if you need speed and have many rows, iterate in vectorized batches or with itertuples

# 2) Build matched_df

#matched_df.to_csv("U937NGST1OENGneu.csv")

#dealing with ion list

import numpy as np
import pandas as pd
from mspvalidator_merger import findingions, expandpeaklist
# assume: from somewhere import score_counter  # returns 0..1

def _ppm_value(val=20.0):
    try:
        return float(ppm(val))  # project helper if present
    except NameError:
        return float(val)


##self-made function
def score_counter_old(hits,ionlist):
    hitcount = len(hits)
    totalioncount = len(ionlist)
    score = hitcount/totalioncount
    return score

#GPT suggestion, avoid 0 division and value <0 or >1
def score_counter(hits, ionlist):
    hitcount = len(hits)
    totalioncount = len(ionlist) if len(ionlist) else 1  # avoid /0
    return min(1.0, max(0.0, hitcount / totalioncount))


def _prepare_ion_df(ion_df, mass_col="mass"):
    """Return (ionlist_mz_array, iondf_for_findingions, mass_col)."""
    if isinstance(ion_df, (list, tuple, np.ndarray)):
        arr = np.asarray(ion_df, dtype=float)
        return arr, pd.DataFrame({mass_col: arr}), mass_col
    if mass_col not in ion_df.columns:
        for c in ("mass","mz","ion_mz"):
            if c in ion_df.columns:
                mass_col = c; break
    arr = ion_df[mass_col].astype(float).to_numpy()
    return arr, ion_df[[mass_col]].copy(), mass_col

def _coerce_peaks(df):
    """Parse tuple-strings if needed."""
    if df.shape[0] and (isinstance(df["peaklist"].iloc[0], str) or isinstance(df["peakintensity"].iloc[0], str)):
        return expandpeaklist(df.copy())
    return df

def attach_ion_score_on_matched(
    matched_df: pd.DataFrame,
    ion_df: pd.DataFrame | list | np.ndarray,
    ppm_value: float = 20.0,
    scan_col: str = "MS2scan_no",
    ion_mass_col: str = "mass",
    score_col: str = "ion score",
    hitcount_col: str = "ion hit count",
    hitlist_col: str = "ion hits m/z",
):
    out = _coerce_peaks(matched_df.copy())
    ionlist_mz, iondf_for_find, ion_mass_col = _prepare_ion_df(ion_df, mass_col=ion_mass_col)

    def _score_row(row):
        hits_out = findingions(row, iondf_for_find, ppm_value)
        # Detect shape:
        #  (A) matched-only: list of mz or (mz, intensity) for matches
        #  (B) full listing:  list of (mz, logI_plus1_or_1) for all reference ions
        if not hits_out:
            return 0.0, 0, ""
        elem = hits_out[0]
        if isinstance(elem, (list, tuple)) and len(elem) == 2 and len(hits_out) == len(ionlist_mz):
            # (B) full listing from findingions
            matched = [(mz, v) for mz, v in hits_out if v > 1.0]
            score = sum(1 for _mz, _v in hits_out if _v > 1.0) / len(hits_out)
            hit_list = ";".join(f"{mz:.4f}" for mz, _ in matched)
            return float(score), len(matched), hit_list
        else:
            # (A) matched-only from peaks_ppm
            # normalize to m/z list for score_counter
            if isinstance(elem, (list, tuple)):
                matched_mz = [mz for mz, *_ in hits_out]     # [(mz,intensity),...] → [mz,...]
            else:
                matched_mz = list(hits_out)                  # [mz,...]
            score = score_counter(matched_mz, ionlist_mz)    # 0..1
            hit_list = ";".join(f"{mz:.4f}" for mz in matched_mz)
            return float(score), len(matched_mz), hit_list

    triples = out.apply(_score_row, axis=1)
    triples = pd.DataFrame(triples.tolist(), index=out.index,
                           columns=[score_col, hitcount_col, hitlist_col])
    return pd.concat([out, triples], axis=1)




# how many anchors in each scan (using your current hit list)
#def parse_hits(s): return [float(x) for x in s.split(";") if x]
#matched_df["anchors"] = matched_df["ion hits m/z"].map(parse_hits).map(has_anchors)
#matched_df["passes_anchors"] = matched_df["anchors"] >= 2
#matched_df["passes_anchors"].mean()  # fraction passing the gate



#extra function to add intensity information back for later ML-based training
def _intensity_from_logi1(logi1):
    """
    findingions returns log10(I + 1) for matches, and 1.0 for 'no hit'.
    Recover I (approx) by I = 10**logi1 - 1 when logi1 > 1.0, else 0.
    """
    if logi1 is None or logi1 <= 1.0:
        return 0.0
    return float(10.0 ** float(logi1)) - 1.0 #fixed 20260122

def attach_ion_hits_with_intensity_on_matched(
    matched_df: pd.DataFrame,
    ion_df: pd.DataFrame | list | np.ndarray,
    ppm_value: float = 20.0,
    scan_col: str = "MS2scan_no",
    ion_mass_col: str = "mass",
    compute_peak_mz: bool = False,
):
    """
    Append per-row JSON arrays of matched ion m/z and matched intensities (raw & rel),
    and return a long (scan × ion) table with match flags and intensities.

    Columns added to matched_df:
      - ion hits m/z        (JSON array of matched ion masses from the ion list)
      - ion hits intensity  (JSON array of raw intensities aligned to the above)
      - ion hits logI       (JSON array of log10(I+1) aligned to the above)
      - ion hits relI       (JSON array of intensity / basepeak aligned to the above)
      (We leave your existing "ion score" / "ion hit count" columns untouched.)

    Returns: updated_matched_df, long_df
      long_df columns:
        [MS2scan_no, ion_mz, matched, logI_plus1, intensity, rel_intensity, best_peak_mz (opt), ppm_error (opt)]
    """
    out = _coerce_peaks(matched_df.copy())
    ionlist_mz, iondf_for_find, ion_mass_col = _prepare_ion_df(ion_df, mass_col=ion_mass_col)

    long_rows = []

    for idx, row in out.iterrows():
        # Basepeak for relative intensity
        peak_ints = np.asarray(row.get("peakintensity", []) if not isinstance(row.get("peakintensity"), str)
                               else json.loads(row["peakintensity"]), dtype=float)
        basepeak = float(peak_ints.max()) if peak_ints.size else 1.0

        hits_out = findingions(row, iondf_for_find, ppm_value)
        matched_mz = []
        matched_I  = []
        matched_logI = []
        matched_relI = []

        # Case A: findingions returns one entry per reference ion: (ion_mz, logI+1)
        # We can still build a full long table (matched or not)
        if hits_out and isinstance(hits_out[0], (list, tuple)) and len(hits_out) == len(ionlist_mz):
            for (ion_mz, logi1), ref_mz in zip(hits_out, ionlist_mz):
                matched = bool(logi1 > 1.0)
                I = _intensity_from_logi1(logi1)
                relI = I / basepeak if basepeak > 0 else 0.0

                best_mz = np.nan
                ppm_err = np.nan
                if compute_peak_mz and matched:
                    # Optional: compute best matched peak m/z & ppm via local search
                    peaks = np.asarray(row["peaklist"], dtype=float)
                    ints  = np.asarray(row["peakintensity"], dtype=float)
                    tol = abs(ref_mz) * float(ppm_value) / 1e6
                    mask = (np.abs(peaks - ref_mz) <= tol) if peaks.size == ints.size and peaks.size > 0 else np.zeros(0, dtype=bool)
                    if np.any(mask):
                        k = int(np.argmax(ints[mask]))
                        best_mz = float(np.asarray(peaks[mask])[k])
                        ppm_err = (best_mz - ref_mz) / ref_mz * 1e6

                long_rows.append({
                    scan_col: row[scan_col],
                    "ion_mz": float(ref_mz),
                    "matched": matched,
                    "logI_plus1": float(logi1),
                    "intensity": float(I),
                    "rel_intensity": float(relI),
                    "best_peak_mz": best_mz,
                    "ppm_error": ppm_err,
                })

                if matched:
                    matched_mz.append(ref_mz)
                    matched_I.append(I)
                    matched_logI.append(float(logi1))
                    matched_relI.append(relI)

        else:
            # Case B: you passed through a matched-only list earlier.
            # We'll still compute intensities by searching peaks within ppm for each matched ion.
            # If hits_out looks like [mz,...] or [(mz,intensity),...], normalize to list of m/z:
            if hits_out:
                first = hits_out[0]
                if isinstance(first, (list, tuple)):
                    matched_only_mz = [float(m) for m, *_ in hits_out]
                else:
                    matched_only_mz = [float(m) for m in hits_out]
            else:
                matched_only_mz = []

            peaks = np.asarray(row["peaklist"], dtype=float)
            ints  = np.asarray(row["peakintensity"], dtype=float)

            # Build a full long row set (including non-matches) for ML
            ion_set = set(np.round(matched_only_mz, 6))
            for ref_mz in ionlist_mz:
                matched = (np.round(ref_mz, 6) in ion_set)
                best_mz = np.nan
                I = 0.0
                logi1 = 1.0
                ppm_err = np.nan
                if matched and peaks.size == ints.size and peaks.size > 0:
                    tol = abs(ref_mz) * float(ppm_value) / 1e6
                    mask = np.abs(peaks - ref_mz) <= tol
                    if np.any(mask):
                        k = int(np.argmax(ints[mask]))
                        best_mz = float(np.asarray(peaks[mask])[k])
                        I = float(np.asarray(ints[mask])[k])
                        logi1 = float(np.log10(I) + 1.0) #fixed 20260122
                        ppm_err = (best_mz - ref_mz) / ref_mz * 1e6
                relI = I / basepeak if basepeak > 0 else 0.0

                long_rows.append({
                    scan_col: row[scan_col],
                    "ion_mz": float(ref_mz),
                    "matched": matched,
                    "logI_plus1": float(logi1),
                    "intensity": float(I),
                    "rel_intensity": float(relI),
                    "best_peak_mz": best_mz,
                    "ppm_error": ppm_err,
                })

                if matched:
                    matched_mz.append(ref_mz)
                    matched_I.append(I)
                    matched_logI.append(logi1)
                    matched_relI.append(relI)

        # Store JSON arrays (compact & lossless alignment)
        out.at[idx, "ion hits m/z"]        = json.dumps([round(float(x), 6) for x in matched_mz])
        out.at[idx, "ion hits intensity"]  = json.dumps([float(x) for x in matched_I])
        out.at[idx, "ion hits logI"]       = json.dumps([float(x) for x in matched_logI])
        out.at[idx, "ion hits relI"]       = json.dumps([float(x) for x in matched_relI])

    long_df = pd.DataFrame(long_rows, columns=[
        scan_col, "ion_mz", "matched", "logI_plus1", "intensity", "rel_intensity",
        "best_peak_mz", "ppm_error"
    ])
    return out, long_df

#avoid call error from mspfileloaderv10.py
if __name__ == "__main__":
#fast converting fixed path
    converted_csv = r"C:\Users\Sakazuki\Desktop\Khoolab_2025data\U937cells_NG_new\ms2_U937_20250729_U937_ST1OE_NGneu.csv"
    composition_list = r"G:\其他電腦\My Computer\GlycoMSParser\U937NG_fix.csv" #this is fault file, mass wrong  #fix wrong mass before adding comp list file
    ion_sheet = r"C:\Users\Sakazuki\Downloads\ionlist_from_zebrafish_NG.csv"
    insilico_df = parse_composition(pd.read_csv(composition_list))
# insilico_df is your new table
# 0) Prepare library once
    lib = normalize_insilico(insilico_df)  # sorted by Mass
    df = pd.read_csv(converted_csv, sep='\t') 
    df = df.copy()
    df["pseudo compositions"] = df["protonatedmass"].apply(lambda m: collect_row_hits({"protonatedmass": m}))
    matched_df = df[df["pseudo compositions"].map(bool)].copy().reset_index()
    print(matched_df.head())

    ion_df = extract_ionmasslist(pd.read_csv(ion_sheet))

    matched_df2 = attach_ion_score_on_matched(
        matched_df,
        ion_df,                 # DataFrame or list/array
        ppm_value=20.0,
        scan_col="MS2scan_no",
        ion_mass_col="mass",
    )
    print(matched_df2.head())
    print("rows:", len(matched_df2))
    print("scans with ≥1 ion hit:", int((matched_df2["ion hit count"] > 0).sum()))
    print("median ion hit count:", matched_df2["ion hit count"].median())
    print("median ion score:", matched_df2["ion score"].median())


    matched_df2.to_csv("U937NGST1OENGneu_withzfion.csv", index=False)


    # top scans by ion score
    top = matched_df2.sort_values("ion score", ascending=False).head(50)

    matched_df3, ion_long = attach_ion_hits_with_intensity_on_matched(
        matched_df2,
        ion_df,                   # DataFrame or list/array
        ppm_value=20.0,
        scan_col="MS2scan_no",
        ion_mass_col="mass",
        compute_peak_mz=False     # set True if you want best_peak_mz + ppm_error (slower)
    )

    # Save both for inspection / ML
    matched_df3.to_csv("pseudolabel_with_intensities.csv", index=False)
    ion_long.to_csv("ion_hits_long_with_intensity.csv", index=False)



#v20250815 we have old version to use
#need to restrict numbers, currently set to float, could be better to limit 4 dec
def getppm(refmass, comparemass, ppm):
    if refmass == 0 or comparemass == 0:
        return np.nan
    elif refmass > comparemass:
        return ((refmass-comparemass)/refmass)* 1e6
    elif comparemass > refmass:
        return ((comparemass-refmass)/comparemass)* 1e6
    else:
        raise ValueError("Invalid input of mass or ppm value provided")




#test line
#parse_composition(composition_list)

#copied from CompositionFinder2.ipynb
def binarycompositionsearch(array, x, low, high, boundary):
    # Repeat until the pointers low and high meet each other
    while low <= high:
        mid = low + (high - low)//2
        #print(f"mid is {mid} now")
        #print(f"value:{array[mid]-x}")
        if abs((array[mid])-x)< boundary:
            return mid
        elif ((array[mid])-x)< -(boundary):
            #print("mid<x")
            low = mid + 1
            #print(f"now lower bound is {low} and it's {array[low]}")
        else:
            #print("mid>x")
            high = mid - 1
            #print(f"now higher bound is {high} and it's {array[high]}")
    return -1



#v20250817 we have old version, not sure where it was
def readunlabeledcsv(csv, ionlist, compositionlist):
    #validation of data format
    loadcsv = pd.read_csv(csv, sep='\t')
    valid_csv = validator.validate_csv_structure(loadcsv)
    if valid_csv:
        print("Call mspvalidator function expandpeaklist to deal with peaklist and intensity")
        peak_flattern_csv = validator.expandpeaklist(loadcsv)
    elif not valid_csv:
        print("probably the csv is truncated or not supported format")
        raise ValueError("the unlabel csv file is invalid for GlycoMSP!")
    else:
        raise Exception("The validation process failed for unknown reason")
    loadionlist = pd.read_csv(ionlist)
    flatternions = extract_ionmasslist(loadionlist)
    
    #tested work
    complist = pd.read_csv(compositionlist)
    #

    peak_flattern_csv.extend(['PseudoComposition'])
    print("haha")





import numpy as np
import pandas as pd

def ppm_qc(matches_long: pd.DataFrame, scan_col="MS2scan_no"):
    if matches_long.empty:
        return {"n_matches": 0}
    s = matches_long["ppm_error"].abs().dropna()
    out = {
        "n_matches": int(s.size),
        "median_abs_ppm": float(np.median(s)),
        "p90_abs_ppm": float(np.percentile(s, 90)),
        "p95_abs_ppm": float(np.percentile(s, 95)),
        "p99_abs_ppm": float(np.percentile(s, 99)),
        "max_abs_ppm": float(s.max()),
        "scans_with_many_hits": matches_long.groupby(scan_col).size().sort_values(ascending=False).head(10).to_dict(),
    }
    return out




#some extra suggestions after testing with GPT

#20260411 code review: need to confirm if below code blocks is active and list the chain of triggering below functions, especially has_anchors and score_counter

#1) Anchor rule (cheap & effective)
# example anchors (placeholder values – swap for your dataset’s anchors)
ANCHORS = [204.087, 366.140, 512.197]  # adjust for derivatization/adduct
ANCHOR_TOL = 0.02  # Da window

def has_anchors(hit_mz_list):
    return sum(any(abs(m - a) <= ANCHOR_TOL for m in hit_mz_list) for a in ANCHORS)

def score_counter(hits_mz, ionlist_mz):
    # gate: at least 2 anchor ions
    if has_anchors(hits_mz) < 2:
        return 0.0
    return len(hits_mz) / max(1, len(ionlist_mz))

#2) Trim the denominator (don’t penalize with rarely observed ions)
#curated = ion_df.query("mass >= 150 & mass <= 2000")  # plus your own whitelist/blacklist
# Use `curated` for scoring, keep full list only for logs


#3) Add intensity awareness (optional, still 0..1)
# If you can switch to an internal matcher that returns the matched peak intensity:
def intensity_fraction(matched_intensities, all_intensities, min_i=0.0):
    num = sum(i for i in matched_intensities if i > min_i)
    den = sum(i for i in all_intensities if i > min_i) or 1.0
    return num / den  # 0..1

# Hybrid score (anchors gate + size + intensity)
#score = 0.5 * (len(hits_mz)/len(ionlist_mz)) + 0.5 * intensity_fraction(...)

#4) Tighten ppm at low m/z (optional)
def ppm_dynamic(mz):
    return 10.0 if mz < 400 else 20.0  # example

#5) Target–decoy FDR (data-driven threshold)
import numpy as np

def decoy_ions(ionlist_mz, shift=50.0):
    # simple decoy: shift ions outside real windows
    return [m + shift for m in ionlist_mz]

# run your scorer twice per scan: real ions vs decoy ions → compare distributions
# choose score threshold where decoy pass-rate / real pass-rate ≈ desired FDR

import os

def extract_ionmasslist(ionmass_sheet, mass_col="mass"):
    if isinstance(ionmass_sheet, (str, os.PathLike)):
        ionmass_sheet = pd.read_csv(ionmass_sheet)  # or read_excel for .xlsx
    df = ionmass_sheet.copy()
    if mass_col not in df.columns:
        for c in ("mass","mz","ion_mz","m/z"):
            if c in df.columns: mass_col = c; break
    return pd.to_numeric(df[mass_col], errors="coerce").dropna().tolist()