# msp_ion_mining.py
version = "0.1"
last_update = 20250906
from __future__ import annotations
import math, re
from ast import literal_eval
from dataclasses import dataclass
from typing import Iterable, List, Tuple, Optional

import numpy as np
import pandas as pd
from bisect import bisect_left

#This file aims to figure out possible important fragmentation ions based on
#Manually annotated spectra information
#We can consider if we should include/implement pseudo-labeled glycan-like spectra into this method
#Currently to avoid bias getting amplified we will stick to manual approach first. 

#Have not tested in 20250906

# -------------------------
# Utilities
# -------------------------

def _to_list(x):
    """Parse peaklist/peakintensity cell to list[float]."""
    if isinstance(x, (list, tuple, np.ndarray)):
        return [float(v) for v in x]
    s = str(x).strip()
    if not s:
        return []
    try:
        vals = literal_eval(s)
        if isinstance(vals, (list, tuple)):
            return [float(v) for v in vals]
    except Exception:
        pass
    s = s.strip("()[]")
    return [float(z) for z in s.split(",") if z.strip()]

def _hybrid_tol(mz: float, ppm: float, da_floor: float) -> float:
    return max(ppm * mz * 1e-6, da_floor)

# -------------------------
# 1D clustering of m/z values (greedy, ppm-aware with DA floor)
# -------------------------

def cluster_mz(values: Iterable[float], ppm: float = 10.0, da_floor: float = 0.03,
               min_count: int = 3) -> np.ndarray:
    """
    Greedy single-pass clustering on sorted m/z values.
    Two points join a cluster if they lie within tol of the current cluster centroid,
    where tol = max(ppm*centroid/1e6, da_floor).
    Returns cluster centers (mean m/z of each cluster) with size >= min_count.
    """
    arr = np.asarray([float(v) for v in values if np.isfinite(v)], dtype=float)
    if arr.size == 0:
        return np.array([], dtype=float)
    arr.sort()
    centers = []
    start = 0
    while start < arr.size:
        centroid = arr[start]
        end = start + 1
        while end < arr.size:
            tol = _hybrid_tol(centroid, ppm, da_floor)
            if arr[end] - centroid <= tol:
                # expand cluster and update centroid
                centroid = (centroid * (end - start) + arr[end]) / (end - start + 1)
                end += 1
            else:
                break
        if (end - start) >= min_count:
            centers.append(centroid)
        start = end
    return np.array(centers, dtype=float)

# -------------------------
# Presence map & scoring
# -------------------------

def _presence_for_centers(mz_list: List[float], centers: np.ndarray,
                          ppm: float, da_floor: float) -> np.ndarray:
    """Return boolean presence vector over centers for a single spectrum."""
    if len(mz_list) == 0 or centers.size == 0:
        return np.zeros(centers.size, dtype=bool)
    m = np.asarray(sorted(mz_list), dtype=float)
    c = centers  # already sorted by construction
    pres = np.zeros(c.size, dtype=bool)
    for mz in m:
        idx = bisect_left(c, mz)
        for j in (idx-1, idx, idx+1):
            if 0 <= j < c.size:
                tol = _hybrid_tol(c[j], ppm, da_floor)
                if abs(mz - c[j]) <= tol:
                    pres[j] = True
    return pres

@dataclass
class SuggestParams:
    ppm: float = 10.0
    da_floor: float = 0.03
    min_cluster_count: int = 3          # min #peak hits to form a cluster center
    min_support_glycan: int = 5         # min #glycan spectra containing the ion
    top_k: int = 60

def suggest_ions_from_pre_df(
    pre_df: pd.DataFrame,
    ion_df: Optional[pd.DataFrame] = None,
    label_col: str = "Structure",
    structure_majority_label: str = "Non-glycan",
    params: SuggestParams = SuggestParams(),
) -> pd.DataFrame:
    """
    Build ion suggestions from pre_df that contains columns:
      - 'peaklist' (list-like), 'peakintensity' (list-like), 'MS2scan_no'
      - a label column, default 'Structure'
    ion_df: used only to exclude ions already present (expects a 'mass' column).
    Returns DataFrame with columns:
      mz, support_glycan, support_non, lift, odds_ratio, already_in_list, recommended
    """
    # 1) Extract all m/z from glycan spectra and cluster to candidate centers
    is_glycan = pre_df[label_col].astype(str).ne(structure_majority_label)
    peaks_all = []
    for mzs in pre_df.loc[is_glycan, "peaklist"].map(_to_list):
        peaks_all.extend(mzs)
    centers = cluster_mz(peaks_all, ppm=params.ppm,
                         da_floor=params.da_floor,
                         min_count=params.min_cluster_count)

    # 2) Presence per spectrum for each center
    n = centers.size
    if n == 0:
        return pd.DataFrame(columns=[
            "mz","support_glycan","support_non","lift","odds_ratio",
            "already_in_list","recommended"
        ])
    pres_g = np.zeros(n, dtype=int)
    pres_n = np.zeros(n, dtype=int)
    # Iterate once through all spectra
    for _, row in pre_df.iterrows():
        lst = _to_list(row.get("peaklist", []))
        pres = _presence_for_centers(lst, centers, params.ppm, params.da_floor)
        if str(row[label_col]) == structure_majority_label:
            pres_n += pres.astype(int)
        else:
            pres_g += pres.astype(int)

    # 3) Score: lift & odds ratio with Laplace smoothing
    g_tot = int(is_glycan.sum())
    n_tot = int((~is_glycan).sum())
    eps = 1.0
    p_g = (pres_g + eps) / (g_tot + 2*eps)
    p_n = (pres_n + eps) / (n_tot + 2*eps)
    lift = p_g / np.maximum(p_n, 1e-12)
    odds_ratio = ((pres_g + eps) * (n_tot - pres_n + eps)) / ((pres_n + eps) * (g_tot - pres_g + eps))

    # 4) Pack table and filter
    out = pd.DataFrame({
        "mz": centers,
        "support_glycan": pres_g,
        "support_non": pres_n,
        "lift": lift,
        "odds_ratio": odds_ratio,
    }).sort_values(["lift","support_glycan"], ascending=[False, False]).reset_index(drop=True)

    if ion_df is not None and "mass" in ion_df.columns:
        known = np.asarray(pd.to_numeric(ion_df["mass"], errors="coerce").dropna().values, dtype=float)
        known.sort()
        def _already(mz):
            idx = bisect_left(known, mz)
            for j in (idx-1, idx, idx+1):
                if 0 <= j < known.size:
                    tol = _hybrid_tol(mz, params.ppm, params.da_floor)
                    if abs(known[j] - mz) <= tol:
                        return True
            return False
        out["already_in_list"] = out["mz"].map(_already)
    else:
        out["already_in_list"] = False

    out = out[out["support_glycan"] >= params.min_support_glycan].reset_index(drop=True)
    out["recommended"] = ~out["already_in_list"]

    return out.head(params.top_k)

def export_ion_suggestions_csv(
    pre_df: pd.DataFrame, ion_df: Optional[pd.DataFrame],
    out_csv: str, params: SuggestParams = SuggestParams(),
    label_col: str = "Structure", majority_label: str = "Non-glycan"
) -> str:
    df = suggest_ions_from_pre_df(
        pre_df, ion_df, label_col=label_col,
        structure_majority_label=majority_label, params=params
    )
    df.to_csv(out_csv, index=False)
    return out_csv