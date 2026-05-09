
import pandas as pd
import numpy as np
from typing import Tuple, Dict

DEFAULTS = dict(
    min_ion_score=0.05,
    min_hit_count=2,
    max_ppm=20.0,
    keep_top_n=1,
    include_mass_cols=False
)

def _safe_eval_tuple_list(x):
    if isinstance(x, (list, tuple, np.ndarray)):
        return np.array(x, dtype=float)
    try:
        val = eval(str(x), {"__builtins__": {}})
        return np.array(val, dtype=float)
    except Exception:
        return np.array([], dtype=float)

def filter_pseudolabels(
    df: pd.DataFrame,
    min_ion_score: float = DEFAULTS["min_ion_score"],
    min_hit_count: int = DEFAULTS["min_hit_count"],
    max_ppm: float = DEFAULTS["max_ppm"],
    keep_top_n: int = DEFAULTS["keep_top_n"],
) -> pd.DataFrame:
    work = df.copy()
    if "ion score" not in work.columns or "ion hit count" not in work.columns or "ppm_error" not in work.columns:
        raise ValueError("Required columns missing: 'ion score', 'ion hit count', 'ppm_error'.")
    if "MS2scan_no" not in work.columns:
        work["MS2scan_no"] = work.get("entry_no", np.arange(len(work)))

    work["ion score"] = pd.to_numeric(work["ion score"], errors="coerce")
    work["ion hit count"] = pd.to_numeric(work["ion hit count"], errors="coerce")
    work["ppm_error"] = pd.to_numeric(work["ppm_error"], errors="coerce")

    mask = (
        (work["ion score"] >= min_ion_score) &
        (work["ion hit count"] >= min_hit_count) &
        (work["ppm_error"].abs() <= max_ppm)
    )
    filt = work.loc[mask].copy()

    filt.sort_values(
        by=["MS2scan_no", "ion score", "ppm_error", "ion hit count"],
        ascending=[True, False, True, False],
        inplace=True
    )
    topn = filt.groupby("MS2scan_no", as_index=False, group_keys=False).head(keep_top_n).copy()
    topn.reset_index(drop=True, inplace=True)
    return topn

def _extract_feature_row_from_peaks(peaks: np.ndarray, intens: np.ndarray, targets: np.ndarray, ppm: float):
    out = {}
    if len(peaks) == 0 or len(intens) == 0 or len(targets) == 0:
        for t in targets:
            out[str(float(t))] = 1.0
        return out

    for t in targets:
        tol = float(t) * ppm / 1e6
        mask = np.abs(peaks - t) <= tol
        if mask.any():
            intensity = float(np.max(intens[mask]))
            out[str(float(t))] = float(np.log10(intensity)+ 1.0) #fixed 20260122
        else:
            out[str(float(t))] = 1.0
    return out

def to_training_dataset(
    df_filtered: pd.DataFrame,
    ionlist_excel_path: str,
    ion_sheet: str = "ionlist",
    ppm: float = DEFAULTS["max_ppm"],
    include_mass_cols: bool = DEFAULTS["include_mass_cols"]
) -> pd.DataFrame:
    xls = pd.ExcelFile(ionlist_excel_path)
    if ion_sheet not in xls.sheet_names:
        raise ValueError(f"Sheet '{ion_sheet}' not found in {ionlist_excel_path}. Available: {xls.sheet_names}")
    ion_df = xls.parse(ion_sheet)
    if "mass" not in ion_df.columns:
        raise ValueError("Ion sheet must contain a 'mass' column.")
    targets = ion_df["mass"].dropna().astype(float).to_numpy()

    meta_cols = ["MS2scan_no", "Structure", "Source"]
    if include_mass_cols:
        meta_cols += ["protonatedmass", "observed_mass", "ppm_error"]
    feat_cols = [str(float(t)) for t in np.sort(targets)]
    if df_filtered is None or len(df_filtered) == 0:
        return pd.DataFrame(columns=meta_cols + feat_cols)

    rows = []
    for _, row in df_filtered.iterrows():
        peaks = _safe_eval_tuple_list(row.get("peaklist", ()))
        intens = _safe_eval_tuple_list(row.get("peakintensity", ()))
        feats = _extract_feature_row_from_peaks(peaks, intens, targets, ppm)

        base = {
            "MS2scan_no": int(row["MS2scan_no"]),
            "Structure": str(row.get("composition", "")),
            "Source": "pseudolabel"
        }
        if include_mass_cols:
            base["protonatedmass"] = float(row.get("protonatedmass", np.nan))
            base["observed_mass"] = float(row.get("observed_mass", np.nan))
            base["ppm_error"] = float(row.get("ppm_error", np.nan))

        base.update(feats)
        rows.append(base)

    out = pd.DataFrame(rows)
    for c in feat_cols:
        if c not in out.columns:
            out[c] = 0.0

    return out[meta_cols + feat_cols]

def compare_against_manual(
    df_filtered: pd.DataFrame,
    manual_excel_path: str,
    manual_sheet: str = "MSlist",
    keep_top_n: int = DEFAULTS["keep_top_n"]
):
    xls = pd.ExcelFile(manual_excel_path)
    if manual_sheet not in xls.sheet_names:
        raise ValueError(f"Sheet '{manual_sheet}' not found in {manual_excel_path}.")
    manual = xls.parse(manual_sheet)

    if "MS2scan_no" not in manual.columns or "Structure" not in manual.columns:
        raise ValueError("Manual sheet must contain 'MS2scan_no' and 'Structure' columns.")

    manual2 = manual[["MS2scan_no", "Structure"]].dropna(subset=["MS2scan_no"]).copy()
    manual2["MS2scan_no"] = pd.to_numeric(manual2["MS2scan_no"], errors="coerce").astype("Int64")
    manual2 = manual2.dropna(subset=["MS2scan_no"]).drop_duplicates(subset=["MS2scan_no"], keep="first")

    if df_filtered is None or len(df_filtered) == 0:
        joined = manual2.copy()
        joined["composition"] = pd.NA
        joined["ion score"] = pd.NA
        joined["ion hit count"] = pd.NA
        joined["ppm_error"] = pd.NA
        joined["match"] = False
        tp = 0
        fp = 0
        fn = int(len(joined))
        precision = 0.0
        recall = 0.0
        f1 = 0.0
        metrics = dict(TP=tp, FP=fp, FN=fn, precision=precision, recall=recall, f1=f1)
        return joined, metrics

    pseudo2 = df_filtered[["MS2scan_no", "composition", "ion score", "ion hit count", "ppm_error"]].copy()
    pseudo2["MS2scan_no"] = pd.to_numeric(pseudo2["MS2scan_no"], errors="coerce").astype("Int64")

    joined = manual2.merge(pseudo2, on="MS2scan_no", how="left", suffixes=("_manual", "_pseudo"))
    joined["match"] = (joined["Structure"].astype(str).str.strip() == joined["composition"].astype(str).str.strip())

    tp = int((joined["match"] == True).sum())
    fp = int((joined["match"] == False).sum())
    fn = int(joined["composition"].isna().sum())

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    metrics = dict(TP=tp, FP=fp, FN=fn, precision=precision, recall=recall, f1=f1)
    return joined, metrics
