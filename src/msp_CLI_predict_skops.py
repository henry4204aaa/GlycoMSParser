"""
Predict with a .skops model on a CSV feature matrix.

Usage:
  python predict_with_skops.py --model path/to/model.skops \
                               --input path/to/unlabeled.csv \
                               --out path/to/predictions.csv \
                               [--drop-cols MS2scan_no protonatedmass] \
                               [--class-col Predicted_Label]

Behavior:
- Loads the model with skops (version-agnostic).
- Aligns feature columns to model.feature_names_in_ if available; otherwise
  tries to use `feature_names` from metadata (if exported).
- Adds proba_top1, proba_top2, margin if the estimator supports predict_proba.
"""

import argparse
import json
import os
import numpy as np
import pandas as pd
from skops.io import load

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, help="Path to .skops model file")
    ap.add_argument("--input", required=True, help="CSV with features (unlabeled)")
    ap.add_argument("--out", required=True, help="Output CSV with predictions")
    ap.add_argument("--drop-cols", nargs="*", default=["MS2scan_no", "protonatedmass"], help="Columns to drop before prediction")
    ap.add_argument("--class-col", default="Predicted_Label", help="Output column name for predicted class")
    args = ap.parse_args()

    # Load model (trusted=True because this is your own artifact)
    model = load(args.model, trusted=True)

    df = pd.read_csv(args.input)
    X = df.drop(columns=args.drop_cols, errors="ignore")

    # Determine training feature order
    train_feats = getattr(model, "feature_names_in_", None)
    if train_feats is None:
        # Try metadata embedded in the .skops file
        meta = getattr(model, "__skops_metadata__", None)
        if meta and isinstance(meta, dict):
            train_feats = meta.get("feature_names") or meta.get("feature_names_in_")
    if train_feats is None:
        raise RuntimeError("No feature_names_in_ on model and no feature metadata embedded. "
                           "Re-export with --features or ensure the estimator stores feature_names_in_.")

    train_feats = list(map(str, train_feats))

    # Reconcile columns: add missing as 0.0, drop extras, enforce order
    missing = sorted(set(train_feats) - set(X.columns))
    extra   = sorted(set(X.columns) - set(train_feats))
    for c in missing:
        X[c] = 0.0
    X = X[train_feats]

    # Predict
    y_pred = model.predict(X)
    df[args.class_col] = y_pred

    # Optional probabilities
    if hasattr(model, "predict_proba"):
        proba = model.predict_proba(X)
        top1 = proba.max(axis=1)
        top2 = np.partition(proba, -2, axis=1)[:, -2] if proba.shape[1] > 1 else np.zeros(len(top1))
        df["proba_top1"] = top1
        df["proba_top2"] = top2
        df["margin"] = top1 - top2

    df.to_csv(args.out, index=False)
    print(f"Wrote predictions: {args.out}")
    if "proba_top1" in df.columns:
        cov = (df["proba_top1"] >= 0.65).mean()
        print(f"(Info) Coverage @ tau=0.65: {cov:.1%}")

if __name__ == "__main__":
    main()