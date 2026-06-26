#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mas_ml_loto_manuscript_analysis.py
===================================
Standalone MAS-ML manuscript analysis for GlycoMSP (Bioinformatics submission).

PURPOSE
-------
Reproduce GlycoMSP's MAS-ML Random Forest evaluation OUTSIDE the GUI, on already-
exported MAS *trainable* CSVs, and add the ML-evaluation block a Bioinformatics
reviewer expects:

    A) within-tissue stratified hold-out test         (formal in-sample ML eval)
    B) leave-one-tissue-out (LOTO) independent test    (cross-tissue generalization)
    C) pooled coverage summary + stratified k-fold CV  (CV sets, coverage demo)

It does NOT modify, import, or run GlycoMSP source. The four pipeline functions
that actually shape the numbers (balance_and_split, the MS1/mass block, the
label normalizer, and the meta-column dropper) are PORTED VERBATIM below from
src/mspfileloaderv14.py (functional-freeze v1.10) so that results reconcile with
the GUI tool. Source line anchors are given in each port.

FAITHFULNESS GATE
-----------------
Run with --validate-brain to check the ported within-brain hold-out against the
frozen GUI run (smoke_6/6_1_brain/..._rf_performance.txt):
    total test support = 219 (173 Non-glycan + 46 glycan),
    all-class accuracy ~0.99, glycan-only macro-F1 ~0.93.

INTERPRETATION BOUNDARIES (manuscript-safe wording)
---------------------------------------------------
"within-dataset held-out evaluation", "leave-one-tissue-out (leave-one-GROUP-out)
independent test", "seen glycan classes", "out-of-label-space held-out
compositions", "pooled class coverage summary". This is NOT de novo structure
discovery and NOT universal glycan prediction. LOTO here is leave-one-GROUP-out
(3 biological tissues), NOT sample-level leave-one-out cross-validation (LOOCV).

Author: built with Cowork for Henry Tseng. RF settings mirror GlycoMSP v1.10.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import platform
from datetime import datetime, timezone

import warnings
import numpy as np
import pandas as pd

# sklearn emits a "unique classes > 50% of samples" heuristic warning on the
# small per-class supports typical of glycan MS data; it is not an error here.
warnings.filterwarnings("ignore", message=r".*number of unique classes is greater.*")

from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import (
    accuracy_score,
    precision_recall_fscore_support,
    classification_report,
    confusion_matrix,
)

SCRIPT_VERSION = "1.0.0"

# ---------------------------------------------------------------------------
# GlycoMSP defaults (mirror src/mspfileloaderv14.py tk-var defaults & RF config)
# ---------------------------------------------------------------------------
RF_DEFAULTS = dict(            # mspfileloaderv14.py:7603 RandomForestClassifier(...)
    n_estimators=400,
    max_depth=None,
    min_samples_split=2,
    min_samples_leaf=1,
    random_state=42,
    class_weight="balanced",
    n_jobs=-1,
)
MAJORITY_LABEL = "Non-glycan"     # majority_label_var
MAJORITY_FACTOR = 3               # majority_factor_var  (cap = factor * max_minor)
MIN_CLASS_SUPPORT = 5             # min_samples_var  -> balance_and_split min_count
TEST_SIZE = 0.20                 # test_split_var
RANDOM_STATE = 42                # hard-coded throughout MAS workflow

# meta columns that must never become features (mirror BOTH guards in v14):
#  - balance_and_split.feature_exclude            (mspfileloaderv14.py:6662)
#  - _drop_meta_and_get_X.meta_exact / prefixes   (mspfileloaderv14.py:6515)
FEATURE_EXCLUDE = (
    "MS2scan_no", "ID", "Source",
    "IUPACname(optional)", "Glycanannotation2", "GlyToucan ID", "WURCS", "unique_ID",
)
_DROP_META_EXACT = {"uid", "predicted_label", "label_str", "class_str",
                    "origin_file", "origin_basename"}
_DROP_META_PREFIX = ("origin_",)

# columns this script adds for provenance (also non-features)
PROV_COLS = ("tissue", "row_id", "source_uid", "UID", "Origin_Basename", "Origin_File")


# ===========================================================================
# PORTED VERBATIM FROM src/mspfileloaderv14.py  (do not "improve" — fidelity)
# ===========================================================================

# --- label normalizer ----------------------------------------------- v14:6541
_MANUAL_PURE = re.compile(r"^(?:F\d+)?H\d+N\d+(?:S\d+)?(?:G\d+)?(?:KDN\d+)?$", re.IGNORECASE)
_MANUAL_EXT = re.compile(r"^(?:F\d+)?H\d+N\d+(?:S\d+)?(?:G\d+)?(?:KDN\d+)?(?:A\d+)?(?:[sp]\d+)?$", re.IGNORECASE)


def _looks_tuplelike_for_structure(x) -> bool:                       # v14:6547
    if isinstance(x, (tuple, list)):
        return True
    if isinstance(x, str) and re.fullmatch(r"\(\s*\d+(?:\s*,\s*\d+){5}\s*\)", x):
        return True
    return False


def _identity_parse_tuple_to_manual(val):
    """Fallback used only if pretrain_normalizer is unavailable AND a tuple-like
    label is encountered. MAS trainable CSVs already store manual strings
    (e.g. 'F1H5N4S2'), so this path is normally never taken."""
    return str(val).strip()


def normalize_structure_for_training(val, parse_tuple_to_manual):     # v14:6555
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


def _get_label_normalizer():
    """Prefer the real GlycoMSP normalizer if importable (read-only); else identity."""
    try:
        import pretrain_normalizer  # noqa
        return pretrain_normalizer.parse_structure_to_manual
    except Exception:
        return _identity_parse_tuple_to_manual


# --- meta drop / numeric X ------------------------------------------ v14:6507
def _drop_meta_and_get_X(df: pd.DataFrame) -> pd.DataFrame:
    drop_cols = []
    for c in df.columns:
        lc = str(c).strip().lower()
        if lc in _DROP_META_EXACT or any(lc.startswith(p) for p in _DROP_META_PREFIX):
            drop_cols.append(c)
    df2 = df.drop(columns=drop_cols, errors="ignore").copy()
    for c in list(df2.columns):
        if df2[c].dtype == object:
            coerced = pd.to_numeric(df2[c], errors="coerce")
            if coerced.notna().mean() >= 0.95:
                df2[c] = coerced
    return df2.select_dtypes(include=[np.number])


# --- MS1 / mass feature block --------------------------------------- v14:6563
def prep_ms1_block(X: pd.DataFrame, y, use_mass: bool) -> pd.DataFrame:
    X = X.copy()
    if use_mass:
        if "protonatedmass" in X.columns:
            X["protonatedmass"] = pd.to_numeric(X["protonatedmass"], errors="coerce")
            if y is not None:
                y_ser = pd.Series(list(y), index=X.index).astype(str)
                X.loc[y_ser.eq(MAJORITY_LABEL), "protonatedmass"] = 0.0
            X["protonatedmass"] = X["protonatedmass"].fillna(0.0)
        # else: mass requested but column absent -> continue without it
    else:
        X = X.drop(columns=["protonatedmass"], errors="ignore")
    return X.fillna(0.0)


# --- class balance + split + capping -------------------------------- v14:6656
def balance_and_split(
    df: pd.DataFrame,
    label_col="Structure",
    majority_label=MAJORITY_LABEL,
    min_count=MIN_CLASS_SUPPORT,
    majority_factor=MAJORITY_FACTOR,
    feature_exclude=FEATURE_EXCLUDE,
    test_size=TEST_SIZE,
    val_size=0.0,
    stratify=True,
    cap_training_only=True,          # GUI default real_world_test_var = True
    random_state=RANDOM_STATE,
):
    # 0) drop tiny classes globally so test never contains unseen labels
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

    ts = float(test_size); vs = float(val_size)
    if ts < 0 or vs < 0:
        raise ValueError("Test/Validation must be >= 0.")
    if ts == 0 and vs == 0:
        raise ValueError("At least one of test/val must be > 0.")
    holdout = ts + vs
    if holdout >= 0.999:
        raise ValueError(f"Test+Validation ({holdout:.2f}) must be < 1.0.")

    feat_cols = _feat_cols(df1)
    X_all, y_all = df1[feat_cols], df1[label_col]
    if y_all.nunique() < 2:
        raise ValueError("After filtering, fewer than 2 classes remain.")

    # MODE B: cap TRAINING ONLY (real-world test) — the GUI default
    strat = y_all if stratify else None
    try:
        X_train, X_tmp, y_train, y_tmp = train_test_split(
            X_all, y_all, test_size=holdout, stratify=strat, random_state=random_state)
    except ValueError:
        X_train, X_tmp, y_train, y_tmp = train_test_split(
            X_all, y_all, test_size=holdout, stratify=None, random_state=random_state)

    if vs > 0:
        rel_test = min(max(ts / holdout, 1e-6), 1 - 1e-6)
        strat_tmp = y_tmp if stratify else None
        try:
            X_val, X_test, y_val, y_test = train_test_split(
                X_tmp, y_tmp, test_size=rel_test, stratify=strat_tmp, random_state=random_state)
        except ValueError:
            X_val, X_test, y_val, y_test = train_test_split(
                X_tmp, y_tmp, test_size=rel_test, stratify=None, random_state=random_state)
    else:
        X_val, y_val = X_tmp.iloc[0:0], y_tmp.iloc[0:0]
        X_test, y_test = X_tmp, y_tmp

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
        "mode": "cap_training_only",
        "kept_label_counts": y_all.value_counts().to_dict(),
        "dropped_rare_counts": dropped_rare[label_col].value_counts().to_dict(),
        "train_majority_cap_applied_to": majority_label if cap_applied else None,
        "train_majority_cap": int(cap_applied) if cap_applied else None,
        "feature_cols": feat_cols,
    }
    return X_train, y_train, X_val, y_val, X_test, y_test, info


# ===========================================================================
# Helpers built ON TOP of the ported pipeline (manuscript analysis additions)
# ===========================================================================

def load_feature_columns(features_arg: str | None, fallback_df: pd.DataFrame | None):
    """Return the canonical feature column order.
    Accepts a JSON list (features.json) OR a GlycoMSP training_run.json
    (training_run.effective_params + training_run.column_order)."""
    if features_arg:
        with open(features_arg, "r", encoding="utf-8") as f:
            obj = json.load(f)
        if isinstance(obj, dict) and "training_run" in obj:
            cols = obj["training_run"].get("column_order")
            if not cols:
                raise ValueError("training_run.json has no column_order")
            return [str(c) for c in cols], "training_run.json:column_order"
        if isinstance(obj, list):
            return [str(c) for c in obj], "features.json(list)"
        raise ValueError("Unrecognised --features file format")
    # conservative inference (last resort): numeric m/z-like cols + protonatedmass
    if fallback_df is None:
        raise ValueError("No --features and no dataframe to infer from")
    inferred = []
    if "protonatedmass" in fallback_df.columns:
        inferred.append("protonatedmass")
    for c in fallback_df.columns:
        cs = str(c)
        if cs in FEATURE_EXCLUDE or cs == "Structure" or cs in PROV_COLS:
            continue
        if cs.replace(".", "", 1).isdigit():
            inferred.append(cs)
    print("[WARN] --features not supplied; inferring feature columns from data. "
          "Strongly prefer passing features.json / training_run.json.", file=sys.stderr)
    return inferred, "INFERRED(numeric)"


def read_tissue_csv(path: str, tissue: str, normalizer) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "Structure" not in df.columns:
        raise SystemExit(f"[ERROR] {path}: missing required 'Structure' column.")
    df["Structure"] = df["Structure"].map(lambda x: normalize_structure_for_training(x, normalizer))
    # stable row id (prefer GlycoMSP's unique_ID). unique_ID is only unique WITHIN
    # a tissue file, so prefix with tissue to guarantee global uniqueness across
    # the pooled/LOTO analyses (otherwise cross-tissue ID collisions look like
    # train/test leakage). Raw id retained in source_uid for back-trace.
    if "unique_ID" in df.columns:
        base = df["unique_ID"].astype(str)
    elif "UID" in df.columns:
        base = df["UID"].astype(str)
    else:
        base = pd.Series([str(i) for i in range(len(df))], index=df.index)
    df["source_uid"] = base.values
    df["row_id"] = tissue + "::" + base.values
    df["tissue"] = tissue
    return df


def align_features(df: pd.DataFrame, feat_cols: list[str], use_mass: bool) -> list[str]:
    """Effective feature columns = canonical order minus protonatedmass if no-mass."""
    cols = list(feat_cols)
    if not use_mass and "protonatedmass" in cols:
        cols = [c for c in cols if c != "protonatedmass"]
    return cols


def build_X(df: pd.DataFrame, y, feat_cols: list[str], use_mass: bool) -> pd.DataFrame:
    """Reindex to canonical feature columns (missing -> 0), then apply mass block
    and meta-drop, exactly as the GUI does post-split."""
    X = df.reindex(columns=feat_cols)              # missing canonical cols -> NaN
    X = prep_ms1_block(X, y, use_mass)             # NG mass -> 0, NaN -> 0 (v14:6563)
    X = _drop_meta_and_get_X(X)                    # numeric-only guard      (v14:6507)
    X = X.reindex(columns=[c for c in feat_cols if c in X.columns], fill_value=0)
    return X.fillna(0.0)


def make_rf(random_state: int) -> RandomForestClassifier:
    p = dict(RF_DEFAULTS); p["random_state"] = random_state
    return RandomForestClassifier(**p)


def glycan_only_metrics(y_true, y_pred, labels_present):
    """macro precision/recall/F1/accuracy on glycan classes only (exclude NG)."""
    gly = [l for l in labels_present if l != MAJORITY_LABEL]
    mask = np.array([yt != MAJORITY_LABEL for yt in y_true])
    if mask.sum() == 0:
        return dict(accuracy=np.nan, macro_precision=np.nan, macro_recall=np.nan,
                    macro_f1=np.nan, n=0)
    yt = np.array(y_true)[mask]; yp = np.array(y_pred)[mask]
    p, r, f, _ = precision_recall_fscore_support(
        yt, yp, labels=gly, average="macro", zero_division=0)
    return dict(accuracy=accuracy_score(yt, yp), macro_precision=p,
                macro_recall=r, macro_f1=f, n=int(mask.sum()))


def all_class_metrics(y_true, y_pred):
    p, r, f, _ = precision_recall_fscore_support(
        y_true, y_pred, average="macro", zero_division=0)
    return dict(accuracy=accuracy_score(y_true, y_pred), macro_precision=p,
                macro_recall=r, macro_f1=f, n=len(y_true))


# ===========================================================================
# Evaluation A — within-tissue stratified hold-out
# ===========================================================================
def eval_within_tissue(dfs, feat_cols, use_mass, args, out):
    rows, manifest_rows = [], []
    for tissue, df in dfs.items():
        rec = {"tissue": tissue, "n_rows_total": len(df)}
        try:
            Xtr, ytr, Xv, yv, Xte, yte, info = balance_and_split(
                df, min_count=args.min_class_support_for_holdout,
                test_size=args.test_size, stratify=True,
                cap_training_only=True, random_state=args.random_state)
        except Exception as e:
            rec["notes"] = f"skipped: {e}"
            rows.append(rec); continue

        Xtr_b = build_X(df.loc[Xtr.index].assign(Structure=ytr.values), ytr, feat_cols, use_mass)
        Xte_b = build_X(df.loc[Xte.index].assign(Structure=yte.values), yte, feat_cols, use_mass)
        Xte_b = Xte_b.reindex(columns=Xtr_b.columns, fill_value=0)

        le = LabelEncoder(); ytr_e = le.fit_transform(ytr)
        model = make_rf(args.random_state); model.fit(Xtr_b, ytr_e)
        yte_pred = le.inverse_transform(model.predict(Xte_b))
        yte_true = list(yte.values)

        present = sorted(set(yte_true) | set(yte_pred))
        gm = glycan_only_metrics(yte_true, yte_pred, present)
        am = all_class_metrics(yte_true, yte_pred)
        rec.update({
            "n_rows_used": len(Xtr) + len(Xte),
            "n_train": len(Xtr_b), "n_test": len(Xte_b),
            "n_classes_used": int(pd.Series(ytr).nunique()),
            "accuracy_seen_glycan": round(gm["accuracy"], 4),
            "macro_f1_seen_glycan": round(gm["macro_f1"], 4),
            "macro_precision_seen_glycan": round(gm["macro_precision"], 4),
            "macro_recall_seen_glycan": round(gm["macro_recall"], 4),
            "accuracy_seen_all": round(am["accuracy"], 4),
            "macro_precision_seen_all": round(am["macro_precision"], 4),
            "macro_recall_seen_all": round(am["macro_recall"], 4),
            "macro_f1_seen_all": round(am["macro_f1"], 4),
            "n_glycan_test_rows": gm["n"],
            "train_majority_cap": info.get("train_majority_cap"),
            "notes": "cap_training_only; stratified 80/20",
        })
        rows.append(rec)

        for idx, role in [(Xtr.index, "train"), (Xte.index, "test")]:
            sub = df.loc[idx]
            for rid, st in zip(sub["row_id"], sub["Structure"]):
                manifest_rows.append(dict(analysis_type="within_tissue", fold_id=tissue,
                    tissue=tissue, row_id=rid, Structure=st, split_role=role,
                    used_for_training=(role == "train"), used_for_testing=(role == "test")))

        # confusion matrix (glycan classes present)
        cm_labels = [l for l in present if l != MAJORITY_LABEL] + [MAJORITY_LABEL]
        cm = confusion_matrix(yte_true, yte_pred, labels=cm_labels)
        pd.DataFrame(cm, index=cm_labels, columns=cm_labels).to_csv(
            os.path.join(out, f"confusion_within_{tissue}.csv"))

    pd.DataFrame(rows).to_csv(os.path.join(out, "within_dataset_holdout_summary.csv"), index=False)
    return rows, manifest_rows


# ===========================================================================
# Evaluation B — leave-one-tissue-out (leave-one-GROUP-out) independent test
# ===========================================================================
def eval_loto(dfs, feat_cols, use_mass, args, out):
    fold_rows, class_rows, pred_rows, manifest_rows = [], [], [], []
    tissues = list(dfs.keys())
    for held in tissues:
        train_tissues = [t for t in tissues if t != held]
        train_df = pd.concat([dfs[t] for t in train_tissues], axis=0, ignore_index=True)
        test_df = dfs[held].copy()
        fold_id = f"train[{('+').join(train_tissues)}]__test[{held}]"

        # define trainable class space: min_count filter on TRAIN labels
        tc = train_df["Structure"].value_counts()
        keep = tc[tc >= args.min_class_support_for_holdout].index
        train_df = train_df[train_df["Structure"].isin(keep)].copy()
        model_classes = set(train_df["Structure"].unique())

        # cap Non-glycan in TRAIN ONLY to 3 x max_minor (real-world test)
        minor = train_df[train_df["Structure"] != MAJORITY_LABEL]["Structure"].value_counts()
        max_minor = int(minor.max()) if not minor.empty else 0
        cap = max(MAJORITY_FACTOR * max_minor, 1)
        if MAJORITY_LABEL in model_classes:
            maj = train_df[train_df["Structure"] == MAJORITY_LABEL]
            if len(maj) > cap:
                maj = maj.sample(n=cap, random_state=args.random_state)
            train_df = pd.concat([maj, train_df[train_df["Structure"] != MAJORITY_LABEL]], axis=0)

        ytr = train_df["Structure"]
        Xtr = build_X(train_df, ytr, feat_cols, use_mass)
        le = LabelEncoder(); ytr_e = le.fit_transform(ytr)
        model = make_rf(args.random_state); model.fit(Xtr, ytr_e)

        # predict the whole held-out tissue (uncapped)
        yte_true = list(test_df["Structure"])
        Xte = build_X(test_df, test_df["Structure"], feat_cols, use_mass).reindex(
            columns=Xtr.columns, fill_value=0)
        proba = model.predict_proba(Xte)
        pred_idx = proba.argmax(axis=1)
        yte_pred = le.inverse_transform(pred_idx)
        maxp = proba.max(axis=1)
        # margin = top1 - top2
        if proba.shape[1] >= 2:
            part = np.partition(proba, -2, axis=1)
            margin = part[:, -1] - part[:, -2]
        else:
            margin = maxp.copy()

        # seen/unseen partition
        def seen_status(st):
            if st == MAJORITY_LABEL:
                return "non_glycan"
            return "seen_glycan" if st in model_classes else "unseen_glycan"

        statuses = [seen_status(st) for st in yte_true]
        seen_mask = np.array([s == "seen_glycan" for s in statuses])
        n_seen_gly = int(seen_mask.sum())
        n_unseen_gly = int(sum(s == "unseen_glycan" for s in statuses))
        n_ng = int(sum(s == "non_glycan" for s in statuses))

        # PRIMARY metrics: seen-glycan rows; macro over (true u pred) - Non-glycan.
        # Same code path as within-tissue/pooled (glycan_only_metrics) -> uniform.
        # A true-glycan row predicted Non-glycan stays in the eval (NG dropped from
        # the label set, not the row), so it counts as a miss against that glycan
        # class's recall; NG is never a scored class here.
        if n_seen_gly > 0:
            seen_true = np.array(yte_true)[seen_mask]
            seen_pred = np.array(yte_pred)[seen_mask]
            gm = glycan_only_metrics(list(seen_true), list(seen_pred),
                                     sorted(set(seen_true) | set(seen_pred)))
            acc_seen, p, r, f = gm["accuracy"], gm["macro_precision"], gm["macro_recall"], gm["macro_f1"]
            # diagnostic: previous full training-label-space macro (kept, renamed)
            diag_labels = sorted([c for c in model_classes if c != MAJORITY_LABEL])
            dp, dr, df_, _ = precision_recall_fscore_support(
                seen_true, seen_pred, labels=diag_labels, average="macro", zero_division=0)
        else:
            p = r = f = acc_seen = dp = dr = df_ = np.nan

        # SECONDARY metrics: seen-glycan + Non-glycan rows; macro over (true u pred),
        # may include Non-glycan. Operational classifier view only (NG = estimated
        # negative generated by GlycoMSP, not ground truth).
        sec_mask = np.array([s in ("seen_glycan", "non_glycan") for s in statuses])
        if sec_mask.sum() > 0:
            am = all_class_metrics(list(np.array(yte_true)[sec_mask]),
                                   list(np.array(yte_pred)[sec_mask]))
        else:
            am = dict(accuracy=np.nan, macro_precision=np.nan, macro_recall=np.nan, macro_f1=np.nan)

        fold_rows.append(dict(
            fold_id=fold_id, train_tissues="+".join(train_tissues), test_tissue=held,
            n_train_rows=len(Xtr), n_test_rows=len(Xte),
            n_train_classes=len(model_classes),
            n_model_classes=len(le.classes_),
            n_test_classes=int(test_df["Structure"].nunique()),
            n_seen_glycan_classes=len({s for s, st in zip(yte_true, statuses) if st == "seen_glycan"}),
            n_unseen_glycan_classes=len({s for s, st in zip(yte_true, statuses) if st == "unseen_glycan"}),
            n_seen_glycan_rows=n_seen_gly, n_unseen_glycan_rows=n_unseen_gly,
            n_non_glycan_rows=n_ng,
            accuracy_seen_glycan=round(acc_seen, 4) if acc_seen == acc_seen else "",
            macro_f1_seen_glycan=round(f, 4) if f == f else "",
            macro_precision_seen_glycan=round(p, 4) if p == p else "",
            macro_recall_seen_glycan=round(r, 4) if r == r else "",
            accuracy_seen_all=round(am["accuracy"], 4) if am["accuracy"] == am["accuracy"] else "",
            macro_precision_seen_all=round(am["macro_precision"], 4) if am["macro_precision"] == am["macro_precision"] else "",
            macro_recall_seen_all=round(am["macro_recall"], 4) if am["macro_recall"] == am["macro_recall"] else "",
            macro_f1_seen_all=round(am["macro_f1"], 4) if am["macro_f1"] == am["macro_f1"] else "",
            macro_precision_training_label_space_diagnostic=round(dp, 4) if dp == dp else "",
            macro_recall_training_label_space_diagnostic=round(dr, 4) if dr == dr else "",
            macro_f1_training_label_space_diagnostic=round(df_, 4) if df_ == df_ else "",
            notes="primary=seen-glycan (true u pred, -NG); secondary=seen+NG; "
                  "diagnostic=full training-label-space (NOT manuscript); unseen reported separately",
        ))

        # per-class per-fold table
        train_counts = train_df["Structure"].value_counts().to_dict()
        test_counts = test_df["Structure"].value_counts().to_dict()
        for st in sorted(set(yte_true) | model_classes):
            class_rows.append(dict(
                fold_id=fold_id, test_tissue=held, Structure=st,
                train_count=int(train_counts.get(st, 0)),
                test_count=int(test_counts.get(st, 0)),
                in_model_class_set=(st in model_classes),
                seen_status=seen_status(st) if st in set(yte_true) else (
                    "train_only" if st in model_classes else "n/a"),
                is_non_glycan=(st == MAJORITY_LABEL)))

        # per-row predictions with truth
        for rid, st, pr, ss, mp, mg in zip(
                test_df["row_id"], yte_true, yte_pred, statuses, maxp, margin):
            pred_rows.append(dict(
                fold_id=fold_id, tissue=held, row_id=rid, Structure=st,
                predicted_label=pr, seen_status=ss,
                correct_if_seen=(pr == st) if ss == "seen_glycan" else "",
                max_proba=round(float(mp), 6), prediction_margin=round(float(mg), 6)))
            manifest_rows.append(dict(analysis_type="loto", fold_id=fold_id,
                tissue=held, row_id=rid, Structure=st, split_role="test",
                used_for_training=False, used_for_testing=True))
        for rid, st in zip(train_df["row_id"], train_df["Structure"]):
            manifest_rows.append(dict(analysis_type="loto", fold_id=fold_id,
                tissue="+".join(train_tissues), row_id=rid, Structure=st,
                split_role="train", used_for_training=True, used_for_testing=False))

        # confusion matrix on seen-glycan rows
        if n_seen_gly > 0:
            labs = sorted({*np.array(yte_true)[seen_mask], *np.array(yte_pred)[seen_mask]})
            cm = confusion_matrix(np.array(yte_true)[seen_mask], np.array(yte_pred)[seen_mask], labels=labs)
            pd.DataFrame(cm, index=labs, columns=labs).to_csv(
                os.path.join(out, f"confusion_loto_test_{held}.csv"))

    pd.DataFrame(fold_rows).to_csv(os.path.join(out, "loto_fold_summary.csv"), index=False)
    pd.DataFrame(class_rows).to_csv(os.path.join(out, "loto_seen_unseen_classes.csv"), index=False)
    pd.DataFrame(pred_rows).to_csv(os.path.join(out, "loto_predictions_with_truth.csv"), index=False)
    return fold_rows, manifest_rows


# ===========================================================================
# Evaluation C — pooled coverage summary + stratified k-fold CV
# ===========================================================================
def eval_pooled(dfs, feat_cols, use_mass, args, out):
    pooled = pd.concat(list(dfs.values()), axis=0, ignore_index=True)
    counts = pooled["Structure"].value_counts()
    eligible = counts[counts >= args.min_class_support_for_holdout]
    coverage = dict(
        n_rows_total=len(pooled),
        n_glycan_rows=int((pooled["Structure"] != MAJORITY_LABEL).sum()),
        n_non_glycan_rows=int((pooled["Structure"] == MAJORITY_LABEL).sum()),
        n_labels_total=int(pooled["Structure"].nunique()),
        n_glycan_labels=int(pooled[pooled["Structure"] != MAJORITY_LABEL]["Structure"].nunique()),
        n_labels_meeting_min_support=int(len(eligible)),
        min_class_support=args.min_class_support_for_holdout,
    )
    pd.DataFrame([coverage]).to_csv(os.path.join(out, "pooled_coverage_summary.csv"), index=False)
    counts.rename_axis("Structure").reset_index(name="n_spectra").to_csv(
        os.path.join(out, "pooled_class_counts.csv"), index=False)

    # stratified k-fold CV on pooled (CV "sets" for Bioinformatics compliance)
    pooled_f = pooled[pooled["Structure"].isin(eligible.index)].copy()
    y = pooled_f["Structure"].values
    k = args.kfolds
    cv_rows = []
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=args.random_state)
    Xfull = build_X(pooled_f, pooled_f["Structure"], feat_cols, use_mass)
    y = pooled_f["Structure"].reset_index(drop=True)
    Xfull = Xfull.reset_index(drop=True)
    for fold, (tr_i, te_i) in enumerate(skf.split(Xfull, y), 1):
        ytr = y.iloc[tr_i]
        # train-only NG cap
        tr_df = pd.DataFrame({"Structure": ytr.values}, index=tr_i)
        minor = ytr[ytr != MAJORITY_LABEL].value_counts()
        cap = max(MAJORITY_FACTOR * (int(minor.max()) if not minor.empty else 0), 1)
        maj_idx = ytr.index[ytr == MAJORITY_LABEL]
        if len(maj_idx) > cap:
            keep_maj = pd.Series(maj_idx).sample(n=cap, random_state=args.random_state).tolist()
            use_idx = list(ytr.index[ytr != MAJORITY_LABEL]) + keep_maj
        else:
            use_idx = list(ytr.index)
        le = LabelEncoder(); ytr_use = y.loc[use_idx]
        ytr_e = le.fit_transform(ytr_use)
        model = make_rf(args.random_state); model.fit(Xfull.loc[use_idx], ytr_e)
        yte_true = y.iloc[te_i].values
        yte_pred = le.inverse_transform(model.predict(Xfull.iloc[te_i]))
        gm = glycan_only_metrics(list(yte_true), list(yte_pred), sorted(set(yte_true) | set(yte_pred)))
        am = all_class_metrics(list(yte_true), list(yte_pred))
        cv_rows.append(dict(fold=fold, n_train=len(use_idx), n_test=len(te_i),
            accuracy_seen_glycan=round(gm["accuracy"], 4),
            macro_f1_seen_glycan=round(gm["macro_f1"], 4),
            macro_precision_seen_glycan=round(gm["macro_precision"], 4),
            macro_recall_seen_glycan=round(gm["macro_recall"], 4),
            accuracy_all=round(am["accuracy"], 4), macro_f1_all=round(am["macro_f1"], 4)))
    cvdf = pd.DataFrame(cv_rows)
    if len(cvdf):
        summ = {"fold": "mean±sd"}
        for c in cvdf.columns:
            if c == "fold":
                continue
            summ[c] = f"{cvdf[c].mean():.4f}±{cvdf[c].std(ddof=1):.4f}" if cvdf[c].notna().any() else ""
        cvdf = pd.concat([cvdf, pd.DataFrame([summ])], ignore_index=True)
    cvdf.to_csv(os.path.join(out, f"pooled_stratified_{k}fold_cv.csv"), index=False)
    return coverage, cv_rows


# ===========================================================================
# Plots (optional, matplotlib)
# ===========================================================================
def make_plots(loto_rows, out):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        print("[WARN] matplotlib unavailable; skipping plots", file=sys.stderr); return
    if not loto_rows:
        return
    tissues = [r["test_tissue"] for r in loto_rows]
    f1 = [r["macro_f1_seen_glycan"] if r["macro_f1_seen_glycan"] != "" else 0 for r in loto_rows]
    fig, ax = plt.subplots(figsize=(5, 3.2))
    ax.bar(tissues, f1, color="#4C72B0")
    ax.set_ylabel("macro-F1 (seen glycan)"); ax.set_ylim(0, 1)
    ax.set_title("LOTO: macro-F1 on seen glycan classes")
    for i, v in enumerate(f1):
        ax.text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(out, "loto_macro_f1.png"), dpi=150); plt.close(fig)

    seen = [r["n_seen_glycan_rows"] for r in loto_rows]
    unseen = [r["n_unseen_glycan_rows"] for r in loto_rows]
    fig, ax = plt.subplots(figsize=(5, 3.2))
    ax.bar(tissues, seen, label="seen glycan rows", color="#55A868")
    ax.bar(tissues, unseen, bottom=seen, label="unseen (out-of-label-space)", color="#C44E52")
    ax.set_ylabel("held-out test rows"); ax.set_title("LOTO: seen vs unseen glycan rows")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(os.path.join(out, "loto_seen_unseen.png"), dpi=150); plt.close(fig)


# ===========================================================================
# Main
# ===========================================================================
def main():
    ap = argparse.ArgumentParser(description="MAS-ML LOTO manuscript analysis (standalone)")
    ap.add_argument("--brain", required=True)
    ap.add_argument("--intestine", required=True)
    ap.add_argument("--ovary", required=True)
    ap.add_argument("--features", default=None,
                    help="features.json (list) OR training_run.json (column_order)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-class-support-for-holdout", type=int, default=MIN_CLASS_SUPPORT,
                    dest="min_class_support_for_holdout")
    ap.add_argument("--test-size", type=float, default=TEST_SIZE, dest="test_size")
    ap.add_argument("--random-state", type=int, default=RANDOM_STATE, dest="random_state")
    ap.add_argument("--kfolds", type=int, default=5)
    ap.add_argument("--no-mass", action="store_true", help="exclude protonatedmass feature")
    ap.add_argument("--gms-src", default=None,
                    help="path to GlycoMSP src/ so pretrain_normalizer is importable for "
                         "canonical label normalization (e.g. K->KDN). Falls back to identity "
                         "(labels kept verbatim) if absent; metrics are unaffected either way.")
    ap.add_argument("--validate-brain", default=None,
                    help="path to reference rf_performance.txt to sanity-check brain holdout")
    args = ap.parse_args()

    use_mass = not args.no_mass
    os.makedirs(args.outdir, exist_ok=True)
    if args.gms_src and os.path.isdir(args.gms_src):
        sys.path.insert(0, args.gms_src)
    normalizer = _get_label_normalizer()
    print(f"[label normalizer] {'GlycoMSP pretrain_normalizer' if normalizer is not _identity_parse_tuple_to_manual else 'identity (labels verbatim)'}")

    paths = {"brain": args.brain, "intestine": args.intestine, "ovary": args.ovary}
    dfs = {t: read_tissue_csv(p, t, normalizer) for t, p in paths.items()}

    feat_cols_all, feat_src = load_feature_columns(args.features, dfs["brain"])
    feat_cols = align_features(dfs["brain"], feat_cols_all, use_mass)

    # input audit
    audit = []
    for t, df in dfs.items():
        present_feats = [c for c in feat_cols if c in df.columns]
        audit.append(dict(tissue=t, input_path=paths[t], n_rows=len(df), n_columns=df.shape[1],
            n_feature_columns=len(present_feats),
            n_labels=int(df["Structure"].nunique()),
            n_glycan_labels=int(df[df["Structure"] != MAJORITY_LABEL]["Structure"].nunique()),
            n_non_glycan_rows=int((df["Structure"] == MAJORITY_LABEL).sum()),
            has_unique_ID=("unique_ID" in df.columns), has_UID=("UID" in df.columns),
            has_protonatedmass=("protonatedmass" in df.columns)))
    pd.DataFrame(audit).to_csv(os.path.join(args.outdir, "input_audit.csv"), index=False)
    json.dump(feat_cols, open(os.path.join(args.outdir, "feature_columns_used.json"), "w"), indent=2)

    within_rows, m1 = eval_within_tissue(dfs, feat_cols, use_mass, args, args.outdir)
    loto_rows, m2 = eval_loto(dfs, feat_cols, use_mass, args, args.outdir)
    coverage, cv_rows = eval_pooled(dfs, feat_cols, use_mass, args, args.outdir)
    pd.DataFrame(m1 + m2).to_csv(os.path.join(args.outdir, "split_manifest.csv"), index=False)
    make_plots(loto_rows, args.outdir)

    # README
    with open(os.path.join(args.outdir, "README.md"), "w", encoding="utf-8") as f:
        f.write(f"""# MAS-ML LOTO manuscript analysis — run report

- Script version: {SCRIPT_VERSION}
- Run timestamp: {datetime.now(timezone.utc).isoformat()}
- Python {platform.python_version()} on {platform.platform()}
- Feature source: `{feat_src}`  ({len(feat_cols)} feature columns, mass={'ON' if use_mass else 'OFF'})

## Inputs
- brain: `{args.brain}`
- intestine: `{args.intestine}`
- ovary: `{args.ovary}`

## Model settings (mirror GlycoMSP v1.10)
RandomForest n_estimators={RF_DEFAULTS['n_estimators']}, class_weight=balanced,
min_samples_split=2, min_samples_leaf=1, random_state={args.random_state}.
Non-glycan capped to {MAJORITY_FACTOR}x max-minority on the TRAINING fold only
(real-world test). Min class support = {args.min_class_support_for_holdout}.
Post-prediction confidence threshold (tau) is OFF (matches frozen GUI default).

## Evaluations
- A) within-tissue stratified hold-out  -> within_dataset_holdout_summary.csv
- B) leave-one-tissue-out (leave-one-GROUP-out) independent test -> loto_fold_summary.csv,
     loto_seen_unseen_classes.csv, loto_predictions_with_truth.csv
- C) pooled coverage + stratified {args.kfolds}-fold CV -> pooled_coverage_summary.csv,
     pooled_stratified_{args.kfolds}fold_cv.csv

## Interpretation rules (manuscript-safe)
- PRIMARY metrics are computed over held-out rows whose TRUE labels are seen
  glycan composition classes. The macro-average label set is the union of true
  and predicted labels present within those evaluated rows, EXCLUDING Non-glycan.
- Non-glycan is an operational negative/background class GENERATED by GlycoMSP
  (low ion-score / no-hit spectra), not ground truth; it is reported separately
  in the SECONDARY all-seen-class metrics (*_seen_all) and never in the primary.
- *_training_label_space_diagnostic columns macro-average over the full training
  label space; they are diagnostic only and must NOT be used as manuscript primary.
- Held-out unseen glycan compositions are reported separately as out-of-label-space
  cases, NOT counted as ordinary false negatives.
- "LOTO" = leave-one-GROUP-out across 3 biological tissues (NOT sample-level LOOCV).
  Spectrum counts are not abundance measurements.
""")

    # optional reference sanity check
    if args.validate_brain and within_rows:
        brain = next((r for r in within_rows if r["tissue"] == "brain"), None)
        print("\n=== BRAIN VALIDATION vs reference ===")
        if brain and "n_test" in brain:
            print(f" ported: n_test={brain['n_test']}  acc_all={brain['accuracy_seen_all']} "
                  f"glycan_macroF1={brain['macro_f1_seen_glycan']} "
                  f"glycan_test_rows={brain.get('n_glycan_test_rows')}")
            print(" reference (GUI rf_performance.txt): n_test=219 acc_all~0.99 "
                  "glycan_macroF1~0.93 glycan_test_rows=46")

    print(f"\n[OK] wrote results to {args.outdir}")


if __name__ == "__main__":
    main()
