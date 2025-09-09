"""
Export a scikit-learn model saved via joblib to a portable .skops file.

Usage:
  python export_model_to_skops.py --joblib path/to/model_rf_model.joblib \
                                  --out path/to/model.skops \
                                  [--features path/to/_features.json]

Notes:
- Run this on a machine that can load the original joblib (i.e., with the same sklearn version used for training).
- The .skops file can then be loaded across different sklearn versions:
    from skops.io import load
    model = load("model.skops", trusted=True)
"""

import argparse
import os
import json
import joblib
import sklearn
from skops.io import dump

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--joblib", required=True, help="Path to joblib-saved sklearn model")
    ap.add_argument("--out", help="Output .skops path (default: same base name)")
    ap.add_argument("--features", help="Optional JSON with the exact training feature names (list[str])")
    args = ap.parse_args()

    in_path = args.joblib
    out_path = args.out or os.path.splitext(in_path)[0] + ".skops"

    # Load the joblib model (must be on a compatible sklearn version)
    model = joblib.load(in_path)

    # Build metadata
    metadata = {
        "sklearn_version": sklearn.__version__,
        "estimator_class": type(model).__name__,
    }
    if args.features and os.path.exists(args.features):
        with open(args.features, "r", encoding="utf-8") as f:
            try:
                feats = json.load(f)
                if isinstance(feats, list):
                    metadata["feature_names"] = feats
            except Exception:
                pass
    # If model exposes feature_names_in_, store it too
    feats_attr = getattr(model, "feature_names_in_", None)
    if feats_attr is not None:
        metadata["feature_names_in_"] = list(map(str, feats_attr))

    # Write .skops
    dump(model, out_path, metadata=metadata)
    print(f"Exported to: {out_path}")
    print(f"Metadata: {json.dumps(metadata, ensure_ascii=False)}")

if __name__ == "__main__":
    main()