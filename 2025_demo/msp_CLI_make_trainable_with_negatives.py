import argparse, sys
import pandas as pd

import mspvalidator_merger as mspval  # must be importable (PYTHONPATH)

def main():
    p = argparse.ArgumentParser(description="Build trainable CSV with optional Non-glycan negatives and decoys.")
    p.add_argument("--excel", required=True, help="Manual annotation Excel (contains MSlist & ionlist)")
    p.add_argument("--raw", required=True, help="Raw TSV converted from extractor")
    p.add_argument("--deriv", default="PerMe(Reduced)", help="Derivatization type string")
    p.add_argument("--out", required=True, help="Output trainable CSV")
    # real negatives
    p.add_argument("--add-neg", action="store_true", help="Add real Non-glycan negatives from raw TSV")
    p.add_argument("--neg-ratio", type=float, default=3.0, help="Max negatives per positive (default 3.0)")
    p.add_argument("--score-col", default="score", help="Score column in raw TSV (or 'none')")
    p.add_argument("--score-thr", type=float, default=0.05, help="Keep scans with max(score) < thr")
    p.add_argument("--marker-cols", default="", help="Comma-separated marker columns (optional)")
    p.add_argument("--marker-min", type=int, default=1, help="Require < this many marker hits")
    # feature-space decoys
    p.add_argument("--add-decoys", action="store_true", help="Append feature-space decoy negatives to the wide CSV")
    p.add_argument("--decoy-method", choices=["permute","sparse_mask","mixup"], default="permute")
    p.add_argument("--decoy-ratio", type=float, default=0.5, help="Decoy rows per positive row")
    args = p.parse_args()

    # 1) Merge manual + peaks → pre_df, ion_index, ion_df
    pre_df, ion_index, ion_df = mspval.directassign_files(args.excel, args.raw, args.deriv, debug=False)

    # 2) Optional: add real negatives
    if args.add_neg:
        from mspvalidator_merger import sample_real_negatives
        marker_cols = [c.strip() for c in args.marker_cols.split(",") if c.strip()] or None
        score_col = None if args.score_col.lower() == "none" else args.score_col
        neg_df = sample_real_negatives(
            args.raw, annotated_scans=pre_df["MS2scan_no"],
            score_col=score_col, score_thr=args.score_thr,
            marker_cols=marker_cols, marker_min=args.marker_min,
            max_neg_ratio=args.neg_ratio
        )
        print(f"[cli] Added real negatives: {len(neg_df)}")
        pre_df = pd.concat([pre_df, neg_df], ignore_index=True)

    # 3) Export wide trainable CSV (manual-ion features)
    mspval.createnormailzedionlistcsv(ion_index, pre_df, ion_df, args.out)
    print(f"[cli] Wrote trainable CSV: {args.out}")

    # 4) Optional: append feature-space decoys directly to the wide CSV file
    if args.add_decoys:
        from mspfileloaderv10 import append_feature_space_decoys  # or place it in a shared utils file
        wide = pd.read_csv(args.out)
        wide2 = append_feature_space_decoys(
            wide, method=args.decoy_method, ratio=args.decoy_ratio
        )
        wide2.to_csv(args.out, index=False)
        print(f"[cli] Appended decoys via '{args.decoy_method}', new rows: {len(wide2)-len(wide)}")

if __name__ == "__main__":
    sys.exit(main())