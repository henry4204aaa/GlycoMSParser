import numpy as np
import pandas as pd
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from sklearn.metrics import f1_score, accuracy_score
from sklearn.utils.class_weight import compute_class_weight

# Reuse predict_with_threshold from the main utils module
from ml_ng_utils import predict_with_threshold


def compute_class_weight_dict(labels: Sequence[str]) -> Dict[str, float]:
    """
    Compute class weights (balanced) as a mapping {label: weight}.
    """
    classes = np.unique(labels)
    weights = compute_class_weight(class_weight="balanced", classes=classes, y=np.array(labels))
    return {lbl: float(w) for lbl, w in zip(classes, weights)}


def resample_by_strategy(
    df: pd.DataFrame,
    label_col: str = "Structure",
    strategy: str = "none",
    max_ng_ratio: float = 1.0,
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Simple resampling strategies to complement/replace class_weight='balanced'.

    strategy:
      - 'none'         : no resampling
      - 'undersample'  : undersample all classes to the minority count
      - 'cap_ng'       : cap 'Non-glycan' to <= majority glycan size * max_ng_ratio
    """
    if strategy == "none":
        return df

    vc = df[label_col].value_counts()
    if strategy == "undersample":
        target = int(vc.min())
        parts = []
        for cls, n in vc.items():
            take = min(n, target)
            parts.append(df[df[label_col] == cls].sample(n=take, random_state=random_state))
        return pd.concat(parts, ignore_index=True)

    if strategy == "cap_ng":
        if "Non-glycan" not in vc.index:
            return df
        glycan_vc = vc.drop("Non-glycan", errors="ignore")
        if glycan_vc.empty:
            return df
        cap = int(glycan_vc.max() * max_ng_ratio)
        ng = df[df[label_col] == "Non-glycan"]
        keep_ng = ng.sample(n=min(len(ng), cap), random_state=random_state)
        rest = df[df[label_col] != "Non-glycan"]
        return pd.concat([rest, keep_ng], ignore_index=True)

    # Fallback
    return df


def tau_sweep_summary(
    model,
    X: np.ndarray,
    y_true_labels: Sequence[str],
    label_encoder,
    taus: Iterable[float] = (0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8),
    ng_label: str = "Non-glycan",
) -> pd.DataFrame:
    """
    Sweep τ and summarize:
      - abstention_rate: % of predictions mapped to NG by threshold
      - macro_f1_glycan: macro-F1 over glycan classes only (exclude NG)
      - f1_ng (if NG in y_true): F1 for the NG class
      - accuracy_on_glycan: accuracy restricted to rows whose true label is NOT NG

    Returns a DataFrame with one row per τ.
    """
    y_true_labels = np.array(y_true_labels)
    glycan_mask = (y_true_labels != ng_label)
    glycan_labels_unique = np.unique(y_true_labels[glycan_mask]) if glycan_mask.any() else np.array([])

    rows = []
    for tau in taus:
        y_pred, maxp, proba, abstained = predict_with_threshold(model, X, label_encoder, tau=tau, ng_label=ng_label)

        abstention_rate = float(np.mean(y_pred == ng_label))

        # Macro-F1 over glycan classes only
        if glycan_labels_unique.size > 0:
            macro_f1_glycan = f1_score(
                y_true_labels[glycan_mask],
                np.array(y_pred)[glycan_mask],
                average="macro",
                labels=glycan_labels_unique
            )
            accuracy_on_glycan = accuracy_score(
                y_true_labels[glycan_mask],
                np.array(y_pred)[glycan_mask],
            )
        else:
            macro_f1_glycan = np.nan
            accuracy_on_glycan = np.nan

        # F1 for NG if present in ground truth
        if np.any(y_true_labels == ng_label):
            f1_ng = f1_score(y_true_labels == ng_label, np.array(y_pred) == ng_label)
        else:
            f1_ng = np.nan

        rows.append({
            "tau": tau,
            "abstention_rate": abstention_rate,
            "macro_f1_glycan": macro_f1_glycan,
            "accuracy_on_glycan": accuracy_on_glycan,
            "f1_ng": f1_ng,
        })

    return pd.DataFrame(rows)