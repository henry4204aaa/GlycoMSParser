# pretrain_normalizer.py
from __future__ import annotations
import re
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split


import re
import numpy as np

COMP_ORDER = ["H","N","S","G","K", "F"]  # Hex, HexNAc, Fuc, NeuAc, NeuGc, KDN

def _tuple_like(x):
    if isinstance(x, (tuple, list, np.ndarray)):
        return True
    if isinstance(x, str) and re.fullmatch(r"\s*\(\s*\d+(?:\s*,\s*\d+){5}\s*\)\s*", x):
        return True
    return False

def _parse_tuple_like(x):
    if isinstance(x, (tuple, list, np.ndarray)):
        vals = list(map(int, list(x)))
    else:
        vals = list(map(int, re.findall(r"\d+", x)))
    if len(vals) != 6:
        raise ValueError("Composition must have 6 integers (H,N,F,S,G,K).")
    return vals  # [H,N,F,S,G,K] by convention

def _compose_label_from_counts(counts, order=COMP_ORDER):
    return "".join(f"{sym}{cnt}" for sym, cnt in zip(order, counts) if cnt > 0)

def _canon_structure_label(x, order=COMP_ORDER):
    # 1) tuple-like → canonical string
    if _tuple_like(x):
        counts = _parse_tuple_like(x)  # [H,N,F,S,G,K]
        return _compose_label_from_counts(counts, order=order)

    # 2) string-like paths handled as before (your previous normalization rules) …
    s = str(x).strip()
    # already H#N#F#...
    if re.fullmatch(r'(?:[HNFSGK]\d+)+', s):
        return s
    z = re.sub(r'[\s,_\-:;]+', '', s).lower()
    token_map = {
        "hexnac":"N","n":"N","gna":"N",
        "hex":"H","h":"H",
        "neuac":"S","sia":"S","s":"S",
        "neugc":"G","g":"G",
        "kdn":"K","k":"K",
        "fuc":"F","f":"F",
    }
    for k,v in sorted(token_map.items(), key=lambda kv: -len(kv[0])):
        z = z.replace(k, v.lower())
    parts = re.findall(r'([hnfsgk])(\d+)', z)
    counts = {k:0 for k in COMP_ORDER}
    for letter,num in parts:
        counts[letter.upper()] += int(num)
    if any(counts.values()):
        return "".join(f"{k}{counts[k]}" for k in order if counts[k] > 0)
    return s



# --- 1) Label canonicalizer ---------------------------------------------------
# Example rule set: turn composition-ish strings into a single form H#N#F#S#G#K (Hex, HexNAc, Fuc, NeuAc, NeuGc, KDN)
TOKEN_MAP = {
    "hexnac": "N", "n": "N", "gna": "N",
    "hex": "H", "h": "H",
    "fuc": "F", "f": "F",
    "neuac": "S", "sia": "S", "s": "S",
    "neugc": "G", "g": "G",
    "kdn": "K", "k": "K",
}
COMP_ORDER = ["H", "N", "F", "S", "G", "K"]

def _canon_structure_label(x: str) -> str:
    if pd.isna(x):
        return x
    s = str(x).strip()

    # fast path: already like H3N2F1...
    if re.fullmatch(r'(?:[HNFSGK]\d+)+', s):
        return s

    # normalize text
    z = re.sub(r'[\s,_\-:;]+', '', s).lower()

    # replace tokens with single letters
    for k, v in sorted(TOKEN_MAP.items(), key=lambda kv: -len(kv[0])):  # longest first
        z = z.replace(k, v.lower())

    # extract letter+integer pairs
    parts = re.findall(r'([hnfsgk])(\d+)', z)
    counts = {k.upper(): 0 for k in COMP_ORDER}
    for letter, num in parts:
        counts[letter.upper()] += int(num)

    if any(v > 0 for v in counts.values()):
        return ''.join(f"{k}{counts[k]}" for k in COMP_ORDER if counts[k] > 0)

    # fallback: return original (e.g., a named class like "Non-glycan")
    return s

# --- 2) UID builder -----------------------------------------------------------
def make_uid(df: pd.DataFrame,
             uid_cols=("experiment","sample","raw_basename","MS2scan_no")) -> pd.Series:
    cols_present = [c for c in uid_cols if c in df.columns]
    if cols_present:
        parts = []
        for c in cols_present:
            if c == "raw_basename":
                parts.append(df[c].astype(str).apply(lambda p: Path(p).stem))
            else:
                parts.append(df[c].astype(str))
        uid = pd.Series(['|'.join(row) for row in zip(*parts)], index=df.index)
    elif "MS2scan_no" in df.columns:
        uid = df["MS2scan_no"].astype(str)
    else:
        uid = df.index.astype(str)
    return uid

# --- 3) Downsample majority class (e.g., "Non-glycan") ------------------------
def cap_majority(df: pd.DataFrame, label_col: str,
                 majority_label: str, cap: int, rng: int = 42) -> pd.DataFrame:
    if majority_label not in df[label_col].unique():
        return df
    maj = df[df[label_col] == majority_label]
    if len(maj) <= cap:
        return df
    keep_idx = maj.sample(n=cap, random_state=rng).index
    return pd.concat([df.loc[keep_idx], df[df[label_col] != majority_label]], axis=0)

# --- 4) Main entry ------------------------------------------------------------
def prepare_for_training(
    df: pd.DataFrame,
    feature_cols: list[str],
    label_col: str = "Structure",
    score_col: str | None = None,       # e.g., pseudolabel score; optional
    majority_label: str | None = "Non-glycan",
    majority_cap_strategy: dict | None = None,  # e.g., {"mode":"x2_of_max_minor", "hard_cap": 10000}
    test_size: float = 0.2,
    val_size: float = 0.1,
    random_state: int = 42,
):
    out = df.copy()

    # canonicalize label
    out[label_col] = out[label_col].apply(_canon_structure_label)

    # add uid (no impact on features; just for auditing)
    out["uid"] = make_uid(out)

    # optional: cap majority
    if majority_label is not None and majority_label in out[label_col].unique():
        if majority_cap_strategy is None:
            # default: cap majority to 2x of the largest minority class
            counts = out[out[label_col] != majority_label][label_col].value_counts()
            max_minor = counts.max() if len(counts) else 0
            cap = max(2 * max_minor, 1)
        else:
            counts = out[out[label_col] != majority_label][label_col].value_counts()
            max_minor = counts.max() if len(counts) else 0
            mode = majority_cap_strategy.get("mode", "x2_of_max_minor")
            hard_cap = majority_cap_strategy.get("hard_cap", 10000)
            if mode == "x2_of_max_minor":
                cap = min(2 * max_minor, hard_cap)
            elif mode == "fixed":
                cap = majority_cap_strategy.get("cap", 5000)
            else:
                cap = min(2 * max_minor, hard_cap)
        out = cap_majority(out, label_col, majority_label, cap, rng=random_state)

    # build sample_weight (optional)
    sample_weight = None
    if score_col and score_col in out.columns:
        s = out[score_col].astype(float).clip(lower=0, upper=1)
        # map to sqrt weighting to soften extremes
        sample_weight = np.sqrt(s.values)

    # stratified train/val/test
    X = out[feature_cols].copy()
    y = out[label_col].copy()

    X_train, X_tmp, y_train, y_tmp, idx_train, idx_tmp = train_test_split(
        X, y, out.index, test_size=(test_size + val_size), stratify=y, random_state=random_state
    )
    relative_test = test_size / (test_size + val_size) if (test_size + val_size) > 0 else 0.0
    X_val, X_test, y_val, y_test, idx_val, idx_test = train_test_split(
        X_tmp, y_tmp, idx_tmp, test_size=relative_test, stratify=y_tmp, random_state=random_state
    )

    sw_train = None
    if sample_weight is not None:
        sw_train = pd.Series(sample_weight, index=out.index).loc[idx_train].values

    # packs
    packs = dict(
        train=dict(X=X_train, y=y_train, idx=idx_train, sample_weight=sw_train),
        val=dict(X=X_val, y=y_val, idx=idx_val),
        test=dict(X=X_test, y=y_test, idx=idx_test),
        cleaned=out,  # full cleaned dataframe with uid
    )
    return packs



"""
#usage

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

df = pd.read_csv("your_pseudolabeled.csv")

feature_cols = [c for c in df.columns if c not in ("Structure","MS2scan_no","score","uid")]
packs = prepare_for_training(
    df,
    feature_cols=feature_cols,
    label_col="Structure",
    score_col="score",                    # omit if you don't have it
    majority_label="Non-glycan",          # set None to skip capping
    majority_cap_strategy={"mode":"x2_of_max_minor","hard_cap":10000},
    test_size=0.2, val_size=0.1, random_state=7
)

rf = RandomForestClassifier(
    n_estimators=400,
    max_features="sqrt",
    class_weight="balanced",              # important!
    random_state=7,
    n_jobs=-1
)
rf.fit(packs["train"]["X"], packs["train"]["y"], sample_weight=packs["train"]["sample_weight"])

y_pred = rf.predict(packs["test"]["X"])
print(classification_report(packs["test"]["y"], y_pred, digits=2))

"""