import ast
import math
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


def build_features_from_peaks_log10_plus1(
    df: pd.DataFrame,
    ion_masses: Sequence[float],
    ppm: float = 20.0,
) -> pd.DataFrame:
    """
    Build a feature matrix from 'peaklist' and 'peakintensity' columns.

    Normalization EXACTLY follows your reference behavior:
        - If a target ion has ≥1 hit within ±ppm: feature = log10(I_max) + 1
        - If no hit:                                    feature = 1.0

    Notes
    -----
    - Expects df['peaklist'] and df['peakintensity'] to be literal-evaluable
      sequences of floats (e.g., strings like "[100.1, 200.2, ...]")
    - ion_masses are the target m/z values (floats), 1 feature per mass.
    """
    if "peaklist" not in df.columns or "peakintensity" not in df.columns:
        raise ValueError("Missing required columns 'peaklist' and/or 'peakintensity'.")

    ion_masses = list(map(float, ion_masses))
    n = len(df)
    m = len(ion_masses)
    out = np.full((n, m), 1.0, dtype=float)  # baseline = 1.0 for 'miss'

    for i, (_, row) in enumerate(df.iterrows()):
        try:
            peaks = np.array(ast.literal_eval(str(row["peaklist"])), dtype=float)
            intens = np.array(ast.literal_eval(str(row["peakintensity"])), dtype=float)
        except Exception:
            # keep baseline row of 1.0s
            continue

        for j, target in enumerate(ion_masses):
            tol = target * ppm / 1e6
            mask = np.abs(peaks - target) <= tol
            if mask.any():
                I = float(np.max(intens[mask]))
                if I < 0:
                    I = 0.0  # guard, matches your reference note
                # IMPORTANT: log10(I) + 1 (NOT log10(I+1))
                out[i, j] = math.log10(I) + 1.0

    cols = [str(mz) for mz in ion_masses]
    return pd.DataFrame(out, index=df.index, columns=cols)


def _is_nohit_row(features_row: np.ndarray) -> bool:
    """
    True if all features equal the baseline (1.0), i.e., no target ions matched.
    """
    if features_row.size == 0:
        return True
    return np.allclose(features_row, 1.0)


def collect_ng_candidates(
    df_pseudo_full: pd.DataFrame,
    df_pseudo_filtered_pos: pd.DataFrame,
    ion_masses: Sequence[float],
    ppm: float = 20.0,
    low_score_col: Optional[str] = "ion score",
    low_score_cut: float = 0.03,
    use_low_score: bool = True,
    use_no_hit: bool = True,
    manual_unknown_df: Optional[pd.DataFrame] = None,
    scan_col: str = "MS2scan_no",
) -> pd.DataFrame:
    """
    Build a Non-glycan candidate table with matched features and label.

    Sources (configurable):
      1) Low-score bin in df_pseudo_full (below low_score_cut).
      2) No-hit spectra (features all baseline after transform).
      3) Manual unknowns (MSlist-only slice), if provided.

    Returns a DataFrame with at least:
      ['MS2scan_no', <feature columns...>, 'Structure']
    where Structure == 'Non-glycan'.
    """
    # normalize scan col name in positives
    if scan_col not in df_pseudo_full.columns:
        raise ValueError(f"scan_col '{scan_col}' not found in df_pseudo_full")
    pos_scans = set(df_pseudo_filtered_pos.get(scan_col, pd.Series([], dtype=int)).astype(int))

    pools = []

    # (1) low-score pool
    if use_low_score and low_score_col and (low_score_col in df_pseudo_full.columns):
        low = df_pseudo_full[df_pseudo_full[low_score_col] < low_score_cut].copy()
        low = low[~low[scan_col].astype(int).isin(pos_scans)]
        pools.append(low)

    # (2) no-hit pool (we will detect after feature building)
    if use_no_hit:
        nh = df_pseudo_full[~df_pseudo_full[scan_col].astype(int).isin(pos_scans)].copy()
        nh["_candidate_nohit"] = True
        pools.append(nh)

    # (3) manual unknowns
    if manual_unknown_df is not None and not manual_unknown_df.empty:
        mu = manual_unknown_df.copy()
        # be tolerant about scan column aliases
        scan_alt = None
        for cand in (scan_col, "MS2scan_no", "ms2scan_no", "ScanNo", "scan", "scan_no", "scan_number"):
            if cand in mu.columns:
                scan_alt = cand
                break
        if scan_alt is None:
            raise ValueError("manual_unknown_df lacks a scan number column.")
        mu = mu[~mu[scan_alt].astype(int).isin(pos_scans)].copy()
        # normalize to unified scan_col
        mu = mu.rename(columns={scan_alt: scan_col})
        pools.append(mu)

    if not pools:
        return pd.DataFrame(columns=[scan_col, "Structure"])

    ng_raw = pd.concat(pools, ignore_index=True).drop_duplicates(subset=[scan_col])

    # Build features (consistent with training)
    feat_df = build_features_from_peaks_log10_plus1(ng_raw, ion_masses, ppm=ppm)

    # If marking no-hit: keep only rows that are truly near-baseline across all features
    if "_candidate_nohit" in ng_raw.columns:
        mask_nohit = (np.isclose(feat_df.values, 1.0)).all(axis=1)
        # NEW: coerce the flag to boolean (NaN→False) before bitwise NOT
        flag = pd.Series(ng_raw["_candidate_nohit"], index=ng_raw.index).astype("boolean").fillna(False)
        keep = (~flag) | mask_nohit
        ng_raw = ng_raw[keep].drop(columns=["_candidate_nohit"])
        feat_df = feat_df.loc[ng_raw.index]

    # Final NG table
    out = pd.concat([ng_raw[[scan_col]].reset_index(drop=True),
                     feat_df.reset_index(drop=True)], axis=1)
    out["Structure"] = "Non-glycan"
    # Canonicalize scan column name
    out = out.rename(columns={scan_col: "MS2scan_no"})
    return out


def cap_non_glycan(
    ng_df: pd.DataFrame,
    train_df: pd.DataFrame,
    label_col: str = "Structure",
    max_ratio: float = 1.0,
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Cap the NG class size to avoid swamping rare classes.
    max_ratio defines the cap relative to the majority glycan class size.
    """
    if ng_df.empty:
        return ng_df
    pos = train_df[train_df[label_col] != "Non-glycan"]
    if pos.empty:
        return ng_df
    majority = pos[label_col].value_counts().max()
    cap = int(max_ratio * majority)
    if len(ng_df) <= cap:
        return ng_df
    return ng_df.sample(n=cap, random_state=random_state)


def predict_with_threshold(
    model,
    X: np.ndarray,
    label_encoder,
    tau: float = 0.55,
    ng_label: str = "Non-glycan",
):
    """
    Hybrid inference:
      - If model already has NG in its classes, keep that prediction.
      - Regardless, apply a max_proba threshold; if max_proba < tau => NG.
    Returns: (y_pred_labels, max_proba, proba, abstained_mask)
    """
    proba = model.predict_proba(X)            # (n, C)
    idx = np.argmax(proba, axis=1)            # argmax class index
    maxp = proba[np.arange(len(idx)), idx]
    y_hat = label_encoder.inverse_transform(idx)

    if ng_label in label_encoder.classes_:
        y_final = np.where(maxp < tau, ng_label, y_hat)
    else:
        y_final = np.where(maxp < tau, ng_label, y_hat)

    abstained = (y_final == ng_label) & (y_hat != ng_label)
    return y_final, maxp, proba, abstained