#from __future__ import annotations
import json, math, os, sys, webbrowser
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ---------- small utils ----------
def _entropy_base2(p: np.ndarray) -> float:
    # stable entropy (sum p log p), allow zeros
    p = np.clip(p, 1e-12, 1.0)
    return float(-np.sum(p * np.log2(p)))

def _open_path_cross_platform(path: Path) -> None:
    try:
        if sys.platform.startswith("darwin"):
            os.system(f'open "{path}"')
        elif os.name == "nt":
            os.startfile(str(path))  # type: ignore[attr-defined]
        else:
            os.system(f'xdg-open "{path}"')
    except Exception:
        try:
            webbrowser.open_new_tab(path.as_uri())
        except Exception:
            pass


# ---------- core API ----------
@dataclass
class ReportParams:
    tau: float = 0.60
    margin: float = 0.05
    topk: int = 5
    sample_cols: Tuple[str, str] = ("experiment_title", "sample_name")

@dataclass
class ReportArtifacts:
    pred_with_conf: Path
    class_summary: Path
    by_sample_summary: Optional[Path]
    review_low_conf: Path
    review_low_margin: Path
    review_ambiguous: Path
    run_summary_json: Path
    run_summary_md: Path

def summarize_predictions(
    pred_df: pd.DataFrame,
    class_names: Optional[List[str]] = None,
    proba_cols: Optional[List[str]] = None,
    params: ReportParams = ReportParams(),
) -> Tuple[pd.DataFrame, pd.DataFrame, Optional[pd.DataFrame], Dict]:
    """
    Inputs
    ------
    pred_df:
        Must contain:
          - 'pred_label'  (top1 label str)
          - either:
              a) probability columns per class in `proba_cols` (floats, 0..1), OR
              b) a single 'proba_vector' column containing list/np.array per row
        Recommended to also contain: 'MS2scan_no', 'protonatedmass'
        Optional: params.sample_cols ('experiment_title','sample_name')

    Returns
    -------
    (pred_rows_df, class_summary_df, by_sample_df_or_None, run_summary_dict)
    """
    df = pred_df.copy()

    # --- get probability matrix ---
    if proba_cols is not None:
        proba_mat = df[proba_cols].to_numpy(dtype=float)
        classes = class_names if class_names is not None else proba_cols
    else:
        # expect a column 'proba_vector' containing array-like
        if "proba_vector" not in df.columns:
            raise ValueError("Either provide proba_cols or include 'proba_vector' column.")
        seqs = df["proba_vector"].apply(lambda x: np.asarray(x, dtype=float))
        # normalize lengths
        lens = seqs.map(len)
        if lens.nunique() != 1:
            raise ValueError("All probability vectors must have the same length.")
        proba_mat = np.stack(seqs.values, axis=0)
        classes = class_names if class_names is not None else [f"class_{i}" for i in range(proba_mat.shape[1])]

    n, C = proba_mat.shape
    if class_names is not None and len(classes) != C:
        raise ValueError("len(class_names) must match number of probability columns/vector length.")

    # --- top2, margin, entropy ---
    top1_idx = np.argmax(proba_mat, axis=1)
    top1 = proba_mat[np.arange(n), top1_idx]
    # Get second best
    # argsort ascending; last two are top2
    sort_idx = np.argsort(proba_mat, axis=1)
    sec_idx = sort_idx[:, -2] if C >= 2 else top1_idx
    sec = proba_mat[np.arange(n), sec_idx]
    margin = top1 - sec

    # entropy per row
    ent = np.apply_along_axis(_entropy_base2, 1, proba_mat)

    df["max_proba"] = top1
    df["second_label"] = [classes[i] for i in sec_idx]
    df["second_proba"] = sec
    df["margin"] = margin
    df["entropy"] = ent

    # acceptance
    tau, m = params.tau, params.margin
    low_conf = df["max_proba"] < tau
    low_margin = df["margin"] < m
    df["accepted"] = ~(low_conf | low_margin)
    def _rej_reason(lc, lm):
        if lc and lm: return "low_conf_and_margin"
        if lc: return "low_conf"
        if lm: return "low_margin"
        return ""
    df["reject_reason"] = [ _rej_reason(lc, lm) for lc, lm in zip(low_conf, low_margin) ]

    # optional topk as CSV strings (compact)
    if params.topk and params.topk > 2:
        K = min(params.topk, C)
        # note: reversed order for top-k descending
        topk_idx = np.flip(np.argsort(proba_mat, axis=1)[:, -K:], axis=1)
        df["topk_labels"] = [ ",".join(classes[j] for j in row) for row in topk_idx ]
        df["topk_probas"] = [ ",".join(f"{proba_mat[i, j]:.4f}" for j in topk_idx[i]) for i in range(n) ]

    # class summary
    grp = df.groupby("pred_label", dropna=False)
    class_summary = grp.agg(
        pred_count=("pred_label", "count"),
        accepted_count=("accepted", "sum"),
        mean_max_proba=("max_proba", "mean"),
        median_max_proba=("max_proba", "median"),
        mean_margin=("margin", "mean"),
        median_margin=("margin", "median"),
        mean_entropy=("entropy", "mean"),
        median_entropy=("entropy", "median"),
    ).reset_index()
    class_summary["accepted_rate_within_class"] = (
        class_summary["accepted_count"] / class_summary["pred_count"].replace({0: np.nan})
    ).fillna(0.0)

    # by-sample (if columns exist)
    by_sample = None
    s_exp, s_sample = params.sample_cols
    if s_exp in df.columns or s_sample in df.columns:
        # fill missing columns if needed
        if s_exp not in df.columns: df[s_exp] = ""
        if s_sample not in df.columns: df[s_sample] = ""
        by_sample = (
            df.groupby([s_exp, s_sample], dropna=False)
              .agg(
                  n_total=("pred_label", "count"),
                  n_accepted=("accepted", "sum"),
                  mean_max_proba=("max_proba", "mean"),
                  mean_margin=("margin", "mean"),
                  mean_entropy=("entropy", "mean"),
              ).reset_index()
        )
        by_sample["n_rejected"] = by_sample["n_total"] - by_sample["n_accepted"]
        by_sample["accepted_rate"] = by_sample["n_accepted"] / by_sample["n_total"].replace({0: np.nan})
        # top3 classes per sample
        # compute distribution per sample via value_counts normalized
        def _top3(sub: pd.DataFrame) -> str:
            vc = sub["pred_label"].value_counts(normalize=True)
            items = [f"{k}:{v:.2f}" for k, v in vc.head(3).items()]
            return "; ".join(items)
        by_sample["top3_classes"] = (
            df.groupby([s_exp, s_sample], dropna=False)
              .apply(_top3)
              .reindex(pd.MultiIndex.from_frame(by_sample[[s_exp, s_sample]]))
              .values
        )

    # run summary dictionary
    def _five(x: np.ndarray) -> Dict[str, float]:
        if len(x) == 0: return {"mean": 0, "median": 0, "p25": 0, "p75": 0}
        return {
            "mean": float(np.mean(x)),
            "median": float(np.median(x)),
            "p25": float(np.percentile(x, 25)),
            "p75": float(np.percentile(x, 75)),
        }

    run_summary = {
        "n_total": int(len(df)),
        "n_accepted": int(df["accepted"].sum()),
        "n_rejected": int((~df["accepted"]).sum()),
        "tau": float(tau),
        "margin": float(m),
        "max_proba_stats": _five(df["max_proba"].to_numpy(np.float64)),
        "margin_stats": _five(df["margin"].to_numpy(np.float64)),
        "entropy_stats": _five(df["entropy"].to_numpy(np.float64)),
        "top_classes_by_count": [
            {"label": k, "count": int(v), "frac": float(v/len(df) if len(df) else 0.0)}
            for k, v in df["pred_label"].value_counts().items()
        ],
        "top_classes_by_accepted": [
            {"label": k, "accepted_count": int(v), "accepted_frac_of_total": float(v/len(df) if len(df) else 0.0)}
            for k, v in df[df["accepted"]]["pred_label"].value_counts().items()
        ],
        "ambiguous_count": int((df["margin"] < m).sum()),
        "low_conf_count": int((df["max_proba"] < tau).sum()),
        "class_names": list(classes),
    }

    # add class_summary into dict (without example ids)
    run_summary["class_summary"] = class_summary.to_dict(orient="records")

    return df, class_summary, by_sample, run_summary


def write_prediction_report(
    out_dir: Path,
    pred_rows: pd.DataFrame,
    class_summary: pd.DataFrame,
    by_sample: Optional[pd.DataFrame],
    run_summary: Dict,
    params: ReportParams = ReportParams(),
) -> ReportArtifacts:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    p1 = out_dir / "predictions_with_conf.csv"
    p2 = out_dir / "class_summary.csv"
    p3 = out_dir / "by_sample_summary.csv"
    r1 = out_dir / "to_review_low_conf.csv"
    r2 = out_dir / "to_review_low_margin.csv"
    r3 = out_dir / "to_review_ambiguous_topk.csv"
    j1 = out_dir / "run_summary.json"
    m1 = out_dir / "run_summary.md"

    # save CSVs
    pred_rows.to_csv(p1, index=False)
    class_summary.to_csv(p2, index=False)
    if by_sample is not None and len(by_sample):
        by_sample.to_csv(p3, index=False)
    else:
        p3 = None  # type: ignore[assignment]

    # reviewers
    tau, m = params.tau, params.margin
    pred_rows[pred_rows["max_proba"] < tau].to_csv(r1, index=False)
    pred_rows[pred_rows["margin"] < m].to_csv(r2, index=False)

    # ambiguous among accepted (runner-up close)
    amb_mask = (pred_rows["accepted"]) & ((pred_rows["second_proba"] > (pred_rows["max_proba"] - m)))
    pred_rows[amb_mask].to_csv(r3, index=False)

    # JSON summary
    j1.write_text(json.dumps(run_summary, indent=2), encoding="utf-8")

    # Markdown summary (compact)
    lines = []
    lines.append("# Prediction run summary\n")
    lines.append(f"- Total: {run_summary['n_total']}")
    lines.append(f"- Accepted: {run_summary['n_accepted']}  (Rejected: {run_summary['n_rejected']})")
    lines.append(f"- τ (confidence): {run_summary['tau']}  ·  margin: {run_summary['margin']}\n")
    mp = run_summary["max_proba_stats"]; mg = run_summary["margin_stats"]; en = run_summary["entropy_stats"]
    lines.append(f"- max_proba — mean {mp['mean']:.3f}, median {mp['median']:.3f}, p25 {mp['p25']:.3f}, p75 {mp['p75']:.3f}")
    lines.append(f"- margin — mean {mg['mean']:.3f}, median {mg['median']:.3f}, p25 {mg['p25']:.3f}, p75 {mg['p75']:.3f}")
    lines.append(f"- entropy — mean {en['mean']:.3f}, median {en['median']:.3f}, p25 {en['p25']:.3f}, p75 {en['p75']:.3f}\n")
    topc = run_summary["top_classes_by_count"][:8]
    if topc:
        lines.append("## Top classes by count")
        for it in topc:
            lines.append(f"- {it['label']}: {it['count']} ({it['frac']:.2%})")
        lines.append("")
    m1.write_text("\n".join(lines), encoding="utf-8")

    return ReportArtifacts(
        pred_with_conf=p1,
        class_summary=p2,
        by_sample_summary=p3,
        review_low_conf=r1,
        review_low_margin=r2,
        review_ambiguous=r3,
        run_summary_json=j1,
        run_summary_md=m1,
    )


# ---------- Tkinter popup (non-blocking) ----------
def show_prediction_summary_popup(
    root,
    run_summary: Dict,
    artifacts: ReportArtifacts,
    class_summary_df: pd.DataFrame,
    max_rows: int = 5,
):
    """
    Small non-blocking dialog with key stats and quick-open buttons.
    """
    import tkinter as tk
    from tkinter import ttk

    top = tk.Toplevel(root)
    top.title("Prediction summary")
    top.attributes("-topmost", True)
    top.resizable(False, False)

    # header
    n_total = run_summary.get("n_total", 0)
    n_acc = run_summary.get("n_accepted", 0)
    n_rej = run_summary.get("n_rejected", 0)
    tau = run_summary.get("tau", 0.0)
    margin = run_summary.get("margin", 0.0)
    acc_rate = (n_acc / n_total) if n_total else 0.0

    hdr = ttk.Frame(top, padding=10)
    hdr.pack(fill="x")
    ttk.Label(hdr, text=f"Total: {n_total}").grid(row=0, column=0, sticky="w", padx=(0,20))
    ttk.Label(hdr, text=f"Accepted: {n_acc}   Rejected: {n_rej}   (rate {acc_rate:.1%})").grid(row=0, column=1, sticky="w")
    ttk.Label(hdr, text=f"τ={tau} · margin={margin}").grid(row=1, column=0, sticky="w", pady=(4,0))

    # small table: top classes
    tbl_frame = ttk.Labelframe(top, text="Top classes (by predictions)", padding=8)
    tbl_frame.pack(fill="x", padx=10, pady=8)

    cols = ("label", "pred_count", "accepted_count", "accepted_rate")
    tree = ttk.Treeview(tbl_frame, columns=cols, show="headings", height=min(max_rows, 8))
    for c in cols:
        tree.heading(c, text=c)
        tree.column(c, width=130 if c!="label" else 180, anchor="center")
    tree.pack(fill="x")

    # derive accepted_rate within class and sort by pred_count
    tmp = class_summary_df.copy()
    if "accepted_rate_within_class" in tmp.columns:
        tmp["accepted_rate"] = tmp["accepted_rate_within_class"]
    else:
        tmp["accepted_rate"] = tmp["accepted_count"] / tmp["pred_count"].replace({0: np.nan})
    tmp = tmp.sort_values("pred_count", ascending=False).head(max_rows)

    for _, r in tmp.iterrows():
        tree.insert("", "end", values=(
            r["pred_label"],
            int(r["pred_count"]),
            int(r["accepted_count"]),
            f"{(r['accepted_rate'] if pd.notnull(r['accepted_rate']) else 0.0):.1%}",
        ))

    # buttons
    btns = ttk.Frame(top, padding=10)
    btns.pack(fill="x")
    def _open(p: Path): _open_path_cross_platform(p)

    ttk.Button(btns, text="Open class_summary.csv",
               command=lambda: _open(artifacts.class_summary)).pack(side="left", padx=(0,8))
    ttk.Button(btns, text="Open low_conf.csv",
               command=lambda: _open(artifacts.review_low_conf)).pack(side="left", padx=(0,8))
    ttk.Button(btns, text="Open run_summary.md",
               command=lambda: _open(artifacts.run_summary_md)).pack(side="left", padx=(0,8))
    ttk.Button(btns, text="Open predictions_with_conf.csv",
               command=lambda: _open(artifacts.pred_with_conf)).pack(side="left")

    # close
    ttk.Button(top, text="Close", command=top.destroy).pack(pady=(0,10))