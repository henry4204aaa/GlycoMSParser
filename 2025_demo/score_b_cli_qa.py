from __future__ import annotations

import ast
import re
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from score_b_skeleton_v3 import (
    ObservedIonHit,
    ScoreBConfig,
    ScoreBResult,
    SpectrumEvidence,
    WorkbookValidationError,
    evaluate_candidate_score_b,
    evaluate_global_motif_evidence,
    load_scoreb_workbook,
    print_candidate_debug,
    summarize_rule_evaluations,
    summarize_target_evidence,
)


COMP_PAT = re.compile(r"(KDN|F|H|N|S|G|A|s|p)\s*([0-9]+)", re.I)


def load_cga_tsv(tsv_path: str | Path) -> pd.DataFrame:
    return pd.read_csv(tsv_path, sep="\t")


def get_labeled_rows(df: pd.DataFrame) -> pd.DataFrame:
    comp = df["composition"].astype("string")
    return df[comp.notna() & comp.str.strip().ne("")].copy()


def show_candidates_for_scan(df: pd.DataFrame, ms2scan_no: int) -> pd.DataFrame:
    rows = df[df["MS2scan_no"] == ms2scan_no].copy()
    cols = [
        "entry_no",
        "MS2scan_no",
        "composition",
        "protonatedmass",
        "theoretical_mass",
        "ppm_error",
        "ion score",
        "ion hit count",
    ]
    cols = [c for c in cols if c in rows.columns]
    return rows[cols]


def select_candidate_row(df: pd.DataFrame, ms2scan_no: int, composition: str) -> pd.Series:
    mask = (
        (df["MS2scan_no"] == ms2scan_no)
        & (df["composition"].astype("string").str.strip() == composition.strip())
    )
    rows = df[mask]
    if rows.empty:
        raise KeyError(f"No candidate row found for MS2scan_no={ms2scan_no}, composition={composition!r}")
    if len(rows) > 1:
        raise ValueError(
            f"Multiple rows found for MS2scan_no={ms2scan_no}, composition={composition!r}; "
            f"use a more specific selector if needed."
        )
    return rows.iloc[0]


def parse_compact_composition_label(label: str) -> tuple[int, int, int, int, int, int]:
    """
    Parse compact composition labels into internal tuple order (H, N, S, G, KDN, F).

    Examples:
        F3H5N4      -> (5, 4, 0, 0, 0, 3)
        H4N4G1KDN1  -> (4, 4, 0, 1, 1, 0)
        H5N4S1      -> (5, 4, 1, 0, 0, 0)
    """
    if not isinstance(label, str) or not label.strip():
        raise ValueError(f"Invalid composition label: {label!r}")

    counts = {"F": 0, "H": 0, "N": 0, "S": 0, "G": 0, "KDN": 0, "A": 0, "s": 0, "p": 0}
    for token, value in COMP_PAT.findall(label):
        token = "KDN" if token.upper() == "KDN" else token
        counts[token] = counts.get(token, 0) + int(value)

    return (
        counts.get("H", 0),
        counts.get("N", 0),
        counts.get("S", 0),
        counts.get("G", 0),
        counts.get("KDN", 0),
        counts.get("F", 0),
    )


def parse_numeric_tuple_string(text: Any) -> list[float]:
    if text is None or (isinstance(text, float) and pd.isna(text)):
        return []
    if isinstance(text, (list, tuple)):
        return [float(x) for x in text]
    if not isinstance(text, str):
        raise ValueError(f"Expected tuple-style string, got {type(text).__name__}: {text!r}")

    stripped = text.strip()
    if not stripped:
        return []

    try:
        value = ast.literal_eval(stripped)
    except Exception as exc:
        raise ValueError(f"Could not parse tuple-style numeric string: {text!r}") from exc

    if isinstance(value, (int, float)):
        return [float(value)]
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"Parsed value is not list/tuple-like: {value!r}")
    return [float(x) for x in value]


def build_peak_pairs_from_row(row: pd.Series) -> list[tuple[float, float]]:
    mzs = parse_numeric_tuple_string(row.get("peaklist"))
    intensities = parse_numeric_tuple_string(row.get("peakintensity"))

    if len(mzs) != len(intensities):
        raise ValueError(
            f"peaklist / peakintensity length mismatch for MS2scan_no={row.get('MS2scan_no')}: "
            f"{len(mzs)} != {len(intensities)}"
        )
    return list(zip(mzs, intensities))


def ppm_error(observed_mz: float, theoretical_mz: float) -> float:
    if theoretical_mz == 0:
        raise ZeroDivisionError("theoretical_mz cannot be zero")
    return (observed_mz - theoretical_mz) / theoretical_mz * 1_000_000.0


def extract_observed_ion_hits_from_peaks(
    peak_pairs: list[tuple[float, float]],
    cfg: ScoreBConfig,
    ppm_tolerance: float = 20.0,
) -> list[ObservedIonHit]:
    """
    Match raw peak pairs against workbook ion definitions.

    Current behavior:
      - at most one best peak is kept per ion_struct_id
      - best is defined by smallest absolute ppm error within tolerance
    """
    observed_hits: list[ObservedIonHit] = []

    for ion_struct_id, ion_def in cfg.ions_by_struct_id.items():
        best: Optional[tuple[float, float, float]] = None  # mz, intensity, ppm
        for obs_mz, obs_intensity in peak_pairs:
            err = ppm_error(obs_mz, ion_def.fragmentation_mass)
            if abs(err) <= ppm_tolerance:
                if best is None or abs(err) < abs(best[2]):
                    best = (obs_mz, obs_intensity, err)
        if best is not None:
            obs_mz, obs_intensity, err = best
            observed_hits.append(
                ObservedIonHit(
                    ion_struct_id=ion_struct_id,
                    mz=obs_mz,
                    intensity=obs_intensity,
                    ppm_error=err,
                )
            )

    return observed_hits


def score_cga_candidate_row(
    row: pd.Series,
    cfg: ScoreBConfig,
    candidate_flags: Optional[dict[str, bool]] = None,
    ppm_tolerance: float = 20.0,
    *,
    charge_mode: str = "any",
    derivatization: str = "any",
) -> tuple[list[ObservedIonHit], SpectrumEvidence, Any, ScoreBResult]:
#) -> tuple[SpectrumEvidence, Any, ScoreBResult]:
    composition_label = str(row.get("composition", "")).strip()
    if not composition_label:
        raise ValueError("Selected row does not have an assigned composition")

    candidate_composition = parse_compact_composition_label(composition_label)
    peak_pairs = build_peak_pairs_from_row(row)
    observed_hits = extract_observed_ion_hits_from_peaks(peak_pairs, cfg, ppm_tolerance=ppm_tolerance)

    spectrum = SpectrumEvidence(
        observed_hits=observed_hits,
        charge_mode=charge_mode,
        derivatization=derivatization,
    )
    global_evidence = evaluate_global_motif_evidence(cfg, spectrum)
    result = evaluate_candidate_score_b(
        cfg,
        global_evidence,
        candidate_composition=candidate_composition,
        candidate_flags=candidate_flags,
    )
    return observed_hits, spectrum, global_evidence, result


def _filter_df(df: pd.DataFrame, *, target_ids: Optional[list[str]] = None, row_indices: Optional[list[int]] = None) -> pd.DataFrame:
    out = df
    if target_ids and "target_id" in out.columns:
        out = out[out["target_id"].isin(target_ids)]
    if row_indices and "row_idx" in out.columns:
        out = out[out["row_idx"].isin(row_indices)]
    return out


def print_row_summary(row: pd.Series) -> None:
    keep_cols = [
        "entry_no",
        "MS1scan_no",
        "MS2scan_no",
        "chargeState",
        "protonatedmass",
        "composition",
        "theoretical_mass",
        "ppm_error",
        "ion score",
        "ion hit count",
    ]
    data = {k: row.get(k) for k in keep_cols if k in row.index}
    print("\n=== Candidate row summary ===")
    for k, v in data.items():
        print(f"{k}: {v}")

#newly added 0325
def summarize_observed_ion_hits(
    observed_hits: list[ObservedIonHit],
    cfg: ScoreBConfig,
) -> pd.DataFrame:
    """
    Summarize extracted ObservedIonHit entries with workbook mass reference.

    Output columns:
        - ion_struct_id
        - theoretical_mz
        - matched_mz
        - ppm_error
        - intensity
    """
    ion_mass_map: dict[str, float] = {}
    for ion_struct_id, ion_def in cfg.ions_by_struct_id.items():
        ion_mass_map[ion_struct_id] = float(ion_def.fragmentation_mass)

    rows: list[dict[str, Any]] = []
    for hit in observed_hits:
        rows.append(
            {
                "ion_struct_id": hit.ion_struct_id,
                "theoretical_mz": ion_mass_map.get(hit.ion_struct_id),
                "matched_mz": hit.mz,
                "ppm_error": hit.ppm_error,
                "intensity": hit.intensity,
            }
        )

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    return df.sort_values(
        by=["ion_struct_id", "ppm_error", "matched_mz"],
        ascending=[True, True, True],
    ).reset_index(drop=True)

def summarize_observed_ion_hits_by_intensity(
    observed_hits: list[ObservedIonHit],
    cfg: ScoreBConfig,
) -> pd.DataFrame:
    df = summarize_observed_ion_hits(observed_hits, cfg)
    if df.empty:
        return df
    return df.sort_values(
        by=["intensity", "ppm_error"],
        ascending=[False, True],
    ).reset_index(drop=True)

def print_observed_ion_hits(
    observed_hits: list[ObservedIonHit],
    cfg: ScoreBConfig,
) -> None:
    print("\n=== Extracted observed ion hits ===")
    hit_df = summarize_observed_ion_hits(observed_hits, cfg)
    if hit_df.empty:
        print("(no extracted observed ion hits)")
    else:
        print(hit_df.to_string(index=False))

#################################

def run_cga_tsv_candidate_case(
    tsv_path: str | Path,
    workbook_path: str | Path,
    *,
    ms2scan_no: int,
    composition: str,
    candidate_flags: dict[str, bool],
    ppm_tolerance: float = 20.0,
    focus_target_ids: Optional[list[str]] = None,
    focus_row_indices: Optional[list[int]] = None,
    charge_mode: str = "any",
    derivatization: str = "any",
) -> None:
    cfg = load_scoreb_workbook(workbook_path)
    df = get_labeled_rows(load_cga_tsv(tsv_path))
    row = select_candidate_row(df, ms2scan_no, composition)

    print_row_summary(row)
    parsed_comp = parse_compact_composition_label(str(row["composition"]))
    print("parsed_composition_tuple:", parsed_comp)

    observed_hits, spectrum, evidence, result = score_cga_candidate_row(
        row,
        cfg,
        candidate_flags=candidate_flags,
        ppm_tolerance=ppm_tolerance,
        charge_mode=charge_mode,
        derivatization=derivatization,
    )
    print_observed_ion_hits(observed_hits, cfg)
    print("extracted_observed_ion_hit_count:", len(observed_hits))

    pass1_df = summarize_rule_evaluations(evidence.pass1_rule_evaluations)
    pass1_df = _filter_df(pass1_df, target_ids=focus_target_ids, row_indices=focus_row_indices)
    #print("\n=== pass1 evals ===")
    #print(pass1_df.to_string(index=False) if not pass1_df.empty else "(no matching pass1 rows)")

    pass2_df = summarize_rule_evaluations(evidence.pass2_rule_evaluations)
    pass2_df = _filter_df(pass2_df, target_ids=focus_target_ids, row_indices=focus_row_indices)
    #print("\n=== pass2 evals ===")
    #print(pass2_df.to_string(index=False) if not pass2_df.empty else "(no matching pass2 rows)")

    target_df = summarize_target_evidence(evidence)
    target_df = _filter_df(target_df, target_ids=focus_target_ids)
    #print("\n=== target evidence ===")
    #print(target_df.to_string(index=False) if not target_df.empty else "(no matching target rows)")

    print_candidate_debug(result)


if __name__ == "__main__":
    WORKBOOK_PATH = Path(r"G:\其他電腦\My Computer\GlycoMSParser\2025_demo\GlycoMSP_scoring_update_example_v7.xlsx")
    #TSV_PATH = Path(r"C:\Users\Sakazuki\Desktop\zf_finalcompare\zf_sPerMeNG_formerge_ovary_pseudolabels_20251007_extratest.tsv")
    TSV_PATH = Path(r"C:\Users\Sakazuki\Desktop\zf_finalcompare\zf_sPerMeNG_formerge_intestine_pseudolabels_20251007_extratest.tsv")
    #TSV_PATH = Path(r"C:\Users\Sakazuki\Desktop\zf_finalcompare\zf_sPerMeNG_formerge_intestine_pseudolabels_20251007.tsv")
    #TSV_PATH = Path(r"C:\Users\Sakazuki\Desktop\zf_finalcompare\zf_sPerMeNG_formerge_ovary_pseudolabels_20251007.tsv")
    #TSV_PATH = Path(r"C:\Users\Sakazuki\Desktop\zf_finalcompare\zf_sPerMeNG_formerge_brain_pseudolabels_20251007.tsv")

    df = get_labeled_rows(load_cga_tsv(TSV_PATH))
    print("=== Candidate rows for MS2scan_no 26313  ===")
    print(show_candidates_for_scan(df, 26313).to_string(index=False))

    run_cga_tsv_candidate_case(
        TSV_PATH,
        WORKBOOK_PATH,
        ms2scan_no=20877,
        composition="F2H5N5KDN1",
        candidate_flags={
            "Allow5Ac": True,
            "AllowLewis": True,
            "Allow5Gc": True,
            "AllowLDNC": True, #True for ovary sample
            "AllowBG": False,
            "AllowFuc": True,
            "Forbid5Gc": False,
            "Forbid5Ac": False,
            "LeA_penalize_LeX": False,
            "LeX_penalize_LeA": False,
            "LeY_penalize_LeB": False,
            "LeB_penalize_LeY": False,
        },
        ppm_tolerance=20.0,
        #focus_target_ids=["Neu5Ac", "Neu5Gc", "SialylAcLacNAc", "SialylAcLDNC", "sLeX", "sLeA"],
    )
