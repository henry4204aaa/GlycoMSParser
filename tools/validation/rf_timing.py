#!/usr/bin/env python3
"""RF-TIMING: one reproducible measurement of GlycoMSP Random Forest training time.

Usage
-----
    python tools/validation/rf_timing.py --trainable <trainable CSV> --alias <name>
                                         [--repeats 3] [--trees 400] [--out <folder>]

What it does
------------
Builds the training matrix exactly as the GUI's "Train" button does for a `_withMass`
trainable (label column `Structure`, "Use protonated mass as a model feature" ON, every
other Train-panel setting at its default) and times, `--repeats` times with
`time.perf_counter()`:

  (a) the `RandomForestClassifier.fit()` call alone, and
  (b) the whole train step, from `pd.read_csv` to the fitted model.

It writes `rf_timing_<alias>.txt` into `--out` (default: the folder holding the CSV).
The report carries only `--alias`: the input path is never written into it.
Nothing in `src/` is imported except `pretrain_normalizer` (no Tk, no GUI start-up) and
nothing in the product is changed.

Where the logic comes from
--------------------------
The GUI keeps its training code as closures inside `train_model()` in
`src/mspfileloaderv14.py` (checked at base commit 26c14da), so it cannot be imported
without starting Tk. The functions below are line-for-line re-implementations; the
sources were located by grep in this session:

  normalize_structure_for_training + _MANUAL_PURE/_MANUAL_EXT   mspfileloaderv14.py:8054-8088
  _drop_meta_and_get_X                                          mspfileloaderv14.py:8020-8049
  balance_and_split (cap-training-only branch, val=0)           mspfileloaderv14.py:8169-8319
  _prep_ms1_block (mass ON)                                     mspfileloaderv14.py:9080-9102
  LabelEncoder fit/transform                                    mspfileloaderv14.py:9109-9113
  RandomForestClassifier(...) parameters                        mspfileloaderv14.py:9120-9128
  model.fit(X_train, y_train_enc)                               mspfileloaderv14.py:9149
  binary glycan-vs-non-glycan accuracy (log line)               mspfileloaderv14.py:9154-9160
  GUI default settings (Tk vars)                                mspfileloaderv14.py:7865-7875
  `pretrain_normalizer.parse_structure_to_manual` is imported from src/ (product code).

GUI defaults reproduced (all stock): test 0.20, val 0.00, min 5 spectra per class,
balancing ON with majority label "Non-glycan" x3, stratified split, random_state 42,
"real world test" ON (the Non-glycan cap is applied to the TRAINING fold only, i.e.
the `cap_training_only` branch), threshold OFF, class_weight "balanced".
Because that cap shrinks the training fold, the report lists the training-fold size
before and after it; the fit time depends on the size after it.

Random Forest parameters as the product sets them: n_estimators=400 (overridable with
--trees), max_depth=None, min_samples_split=2, min_samples_leaf=1, random_state=42,
class_weight="balanced", and n_jobs=-1 (the product passes n_jobs=-1 explicitly at :9127,
so this script does too).

Not timed: the GUI's file hash, prediction, report writing and model saving.
"""
from __future__ import annotations

import argparse
import os
import platform
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
import pretrain_normalizer as normalizer  # noqa: E402  (product module; no Tk)

LABEL_COL = "Structure"
MAJORITY_LABEL = "Non-glycan"
MAJORITY_FACTOR = 3
MIN_COUNT = 5
TEST_SIZE = 0.20
RANDOM_STATE = 42
# balance_and_split's own default feature_exclude (mspfileloaderv14.py:8175)
FEATURE_EXCLUDE = ("MS2scan_no", "ID", "Source", "IUPACname(optional)", "Glycanannotation2",
                   "GlyToucan ID", "WURCS", "unique_ID")

_ALIAS_OK = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")

# --- label normalisation (mspfileloaderv14.py:8054-8088) -------------------------------
_MANUAL_PURE = re.compile(r"^(?:F\d+)?H\d+N\d+(?:S\d+)?(?:G\d+)?(?:KDN\d+)?$", re.IGNORECASE)
_MANUAL_EXT = re.compile(
    r"^(?:F\d+)?H\d+N\d+(?:S\d+)?(?:G\d+)?(?:KDN\d+)?(?:A\d+)?(?:[sp]\d+)?$", re.IGNORECASE)


def _looks_tuplelike_for_structure(x) -> bool:
    if isinstance(x, (tuple, list)):
        return True
    if isinstance(x, str) and re.fullmatch(r"\(\s*\d+(?:\s*,\s*\d+){5}\s*\)", x):
        return True
    return False


def normalize_structure_for_training(val, parse_tuple_to_manual):
    s = str(val).strip()
    if _looks_tuplelike_for_structure(val):
        return parse_tuple_to_manual(val)
    if _MANUAL_EXT.match(s) or _MANUAL_PURE.match(s):
        return s
    if re.fullmatch(r"(?:[A-Za-z]+?\d+)+", s):
        try:
            return parse_tuple_to_manual(s)
        except Exception:
            return s
    return s


# --- numeric feature matrix (mspfileloaderv14.py:8020-8049) -----------------------------
def _drop_meta_and_get_X(df: pd.DataFrame) -> pd.DataFrame:
    meta_exact = {"uid", "predicted_label", "label_str", "class_str",
                  "origin_file", "origin_basename"}
    meta_prefix = ("origin_",)
    drop_cols = []
    for c in df.columns:
        lc = str(c).strip().lower()
        if lc in meta_exact or any(lc.startswith(p) for p in meta_prefix):
            drop_cols.append(c)
    df2 = df.drop(columns=drop_cols, errors="ignore").copy()
    for c in list(df2.columns):
        if df2[c].dtype == object:
            coerced = pd.to_numeric(df2[c], errors="coerce")
            if coerced.notna().mean() >= 0.95:
                df2[c] = coerced
    return df2.select_dtypes(include=[np.number])


# --- balance + split, cap-training-only branch, val=0 (mspfileloaderv14.py:8169-8319) ---
def balance_and_split(df, label_col=LABEL_COL, majority_label=MAJORITY_LABEL,
                      min_count=MIN_COUNT, majority_factor=MAJORITY_FACTOR,
                      feature_exclude=FEATURE_EXCLUDE, test_size=TEST_SIZE,
                      stratify=True, random_state=RANDOM_STATE):
    counts = df[label_col].value_counts()
    keep = counts[counts >= min_count].index
    dropped_rare = df[~df[label_col].isin(keep)]
    df1 = df[df[label_col].isin(keep)].copy()

    def _feat_cols(dff):
        return [c for c in dff.columns if c not in set(feature_exclude) | {label_col}]

    def _cap_majority_df(dff, cap):
        if majority_label not in dff[label_col].unique():
            return dff, None
        maj = dff[dff[label_col] == majority_label]
        if len(maj) <= cap:
            return dff, None
        maj_keep = maj.sample(n=int(cap), random_state=random_state)
        dff2 = pd.concat([maj_keep, dff[dff[label_col] != majority_label]], axis=0)
        return dff2, int(cap)

    feat_cols = _feat_cols(df1)
    X_all, y_all = df1[feat_cols], df1[label_col]
    if y_all.nunique() < 2:
        raise ValueError("After filtering, fewer than 2 classes remain.")

    strat = y_all if stratify else None
    try:
        X_train, X_tmp, y_train, y_tmp = train_test_split(
            X_all, y_all, test_size=test_size, stratify=strat, random_state=random_state)
    except ValueError:
        X_train, X_tmp, y_train, y_tmp = train_test_split(
            X_all, y_all, test_size=test_size, stratify=None, random_state=random_state)
    X_val, y_val = X_tmp.iloc[0:0], y_tmp.iloc[0:0]
    X_test, y_test = X_tmp, y_tmp

    train_before_cap = len(X_train)
    if majority_label in y_train.unique():
        minor_counts_train = y_train[y_train != majority_label].value_counts()
        max_minor_train = int(minor_counts_train.max()) if not minor_counts_train.empty else 0
        cap_train = max(majority_factor * max_minor_train, 1)
        train_df = df1.loc[X_train.index].copy()
        train_df_capped, cap_applied = _cap_majority_df(train_df, cap_train)
        X_train = train_df_capped[feat_cols]
        y_train = train_df_capped[label_col]
    else:
        cap_applied = None

    info = {
        "kept_label_counts": y_all.value_counts().to_dict(),
        "dropped_rare_counts": dropped_rare[label_col].value_counts().to_dict(),
        "rows_after_filter": int(len(df1)),
        "train_before_cap": int(train_before_cap),
        "cap_applied": cap_applied,
    }
    return X_train, y_train, X_val, y_val, X_test, y_test, info


# --- MS1 block with the mass feature ON (mspfileloaderv14.py:9080-9102) -----------------
def _prep_ms1_block(X, y, use_mass=True):
    X = X.copy()
    if use_mass:
        if "protonatedmass" in X.columns:
            X["protonatedmass"] = pd.to_numeric(X["protonatedmass"], errors="coerce")
            if y is not None:
                y_ser = pd.Series(list(y), index=X.index).astype(str)
                X.loc[y_ser.eq("Non-glycan"), "protonatedmass"] = 0.0
            X["protonatedmass"] = X["protonatedmass"].fillna(0.0)
    else:
        X = X.drop(columns=["protonatedmass"], errors="ignore")
    return X.fillna(0.0)


def train_step(csv_path: str, n_trees: int):
    """CSV load -> fitted model, as train_model() does it. Returns (state, fit_seconds)."""
    df = pd.read_csv(csv_path)
    if LABEL_COL not in df.columns:
        raise ValueError(f"label column '{LABEL_COL}' not in trainable csv")
    n_rows_total = int(len(df))
    df[LABEL_COL] = df[LABEL_COL].map(
        lambda x: normalize_structure_for_training(x, normalizer.parse_structure_to_manual))
    y_labels = df[LABEL_COL].astype(str)

    X_train, y_train, X_val, y_val, X_test, y_test, info = balance_and_split(
        df.assign(**{LABEL_COL: y_labels}))

    X_train = _prep_ms1_block(X_train, y_train)
    X_test = _prep_ms1_block(X_test, y_test)

    le = LabelEncoder()
    y_train_enc = le.fit_transform(y_train)
    y_test_enc = le.transform(y_test)

    model = RandomForestClassifier(
        n_estimators=int(n_trees), max_depth=None, min_samples_split=2, min_samples_leaf=1,
        random_state=42, class_weight="balanced", n_jobs=-1)

    X_train_num = _drop_meta_and_get_X(X_train)
    cols = list(X_train_num.columns)
    X_test_num = _drop_meta_and_get_X(X_test).reindex(columns=cols, fill_value=0)
    if X_train_num.isna().any().any():
        raise RuntimeError("NaNs remain in X_train after MS1 prep")

    t0 = time.perf_counter()
    model.fit(X_train_num, y_train_enc)
    fit_s = time.perf_counter() - t0

    state = {
        "model": model, "le": le, "X_test": X_test_num, "y_test_enc": y_test_enc,
        "rows_total": n_rows_total, "rows_after_filter": info["rows_after_filter"],
        "n_features": len(cols), "mass_included": "protonatedmass" in cols,
        "classes_kept": len(info["kept_label_counts"]),
        "classes_dropped": len(info["dropped_rare_counts"]),
        "train_before_cap": info["train_before_cap"], "train_rows": int(len(X_train_num)),
        "test_rows": int(len(X_test_num)), "cap_applied": info["cap_applied"],
    }
    return state, fit_s


def _cpu_model() -> str:
    """CPU model string; platform.processor() is the Windows source (as the ticket asks)."""
    sysname = platform.system()
    try:
        if sysname == "Darwin":
            out = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                                 capture_output=True, text=True, timeout=5).stdout.strip()
            if out:
                return out
        elif sysname == "Linux":
            with open("/proc/cpuinfo", encoding="utf-8") as f:
                for line in f:
                    if line.lower().startswith("model name"):
                        return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return platform.processor() or "unknown"


def _fmt(values):
    return ", ".join(f"{v:.3f}" for v in values)


def build_report(alias, trees, repeats, state, fit_times, full_times, acc, acc_bin) -> str:
    cap = state["cap_applied"]
    lines = [
        "RF training timing report",
        f"alias: {alias}",
        f"generated_utc: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "[data]",
        f"rows_total: {state['rows_total']}",
        f"rows_after_eligibility_filter: {state['rows_after_filter']}",
        f"features: {state['n_features']}",
        f"protonatedmass_feature_included: {'yes' if state['mass_included'] else 'no'}",
        f"classes_kept: {state['classes_kept']}",
        f"classes_dropped_below_5_spectra: {state['classes_dropped']}",
        f"train_rows_before_majority_cap: {state['train_before_cap']}",
        f"majority_cap_applied_to_training: {cap if cap else 'no cap applied'}",
        f"train_rows: {state['train_rows']}",
        f"test_rows: {state['test_rows']}",
        "",
        "[random forest, as the product sets it]",
        f"n_estimators: {trees}",
        "class_weight: balanced",
        "random_state: 42",
        "max_depth: None",
        "min_samples_split: 2",
        "min_samples_leaf: 1",
        "n_jobs: -1 (the product passes n_jobs=-1 explicitly; all logical cores)",
        "split: 80/20 stratified, random_state=42; Non-glycan capped x3 in the training fold only",
        "",
        "[environment]",
        f"python: {platform.python_version()}",
        f"scikit-learn: {sklearn.__version__}",
        f"numpy: {np.__version__}",
        f"pandas: {pd.__version__}",
        f"platform: {platform.platform()}",
        f"processor_string: {platform.processor() or 'unknown'}",
        f"cpu_model: {_cpu_model()}",
        f"logical_cores: {os.cpu_count()}",
        "",
        f"[timing, {repeats} repeat(s), seconds, time.perf_counter()]",
        f"fit_seconds_per_repeat: {_fmt(fit_times)}",
        f"full_step_seconds_per_repeat (CSV load -> fitted model): {_fmt(full_times)}",
        f"fit_seconds_mean: {sum(fit_times) / len(fit_times):.3f}",
        f"fit_seconds_max: {max(fit_times):.3f}",
        f"full_step_seconds_mean: {sum(full_times) / len(full_times):.3f}",
        f"full_step_seconds_max: {max(full_times):.3f}",
        "",
        "[held-out accuracy of the last fitted model, threshold off]",
        f"multiclass_accuracy: {acc:.4f}   (the 'accuracy' row of the GUI's rf_performance.txt)",
        "glycan_vs_nonglycan_binary_accuracy: "
        + (f"{acc_bin:.4f}   (the GUI's 'Binary glycan vs non-glycan' log line)"
           if acc_bin is not None else "n/a (no Non-glycan class kept)"),
        "",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Time the GUI's Random Forest training step.")
    ap.add_argument("--trainable", required=True, help="trainable CSV (_withMass)")
    ap.add_argument("--alias", required=True, help="short name used in the report and file name")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--trees", type=int, default=400)
    ap.add_argument("--out", default=None, help="output folder (default: next to the CSV)")
    args = ap.parse_args(argv)

    if not _ALIAS_OK.match(args.alias):
        print("error: --alias must be 1-64 characters of letters, digits, '_', '.', '-'",
              file=sys.stderr)
        return 2
    if args.repeats < 1 or args.trees < 1:
        print("error: --repeats and --trees must be >= 1", file=sys.stderr)
        return 2
    csv_path = Path(args.trainable)
    if not csv_path.is_file():
        print("error: --trainable is not a file", file=sys.stderr)
        return 2
    out_dir = Path(args.out) if args.out else csv_path.resolve().parent
    out_dir.mkdir(parents=True, exist_ok=True)

    fit_times, full_times, state = [], [], None
    for _ in range(args.repeats):
        t0 = time.perf_counter()
        state, fit_s = train_step(str(csv_path), args.trees)
        full_times.append(time.perf_counter() - t0)
        fit_times.append(fit_s)

    # held-out accuracy of the last fitted model (outside the timed region)
    y_pred = state["model"].predict(state["X_test"])
    y_true = state["y_test_enc"]
    acc = float(accuracy_score(y_true, y_pred))
    classes = list(state["le"].classes_)
    acc_bin = None
    if MAJORITY_LABEL in classes:
        maj_idx = classes.index(MAJORITY_LABEL)
        acc_bin = float(accuracy_score((y_true != maj_idx).astype(int),
                                       (y_pred != maj_idx).astype(int)))

    report = build_report(args.alias, args.trees, args.repeats, state,
                          fit_times, full_times, acc, acc_bin)
    report_path = out_dir / f"rf_timing_{args.alias}.txt"
    report_path.write_text(report, encoding="utf-8")
    print(report)
    print(f"report written: rf_timing_{args.alias}.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
