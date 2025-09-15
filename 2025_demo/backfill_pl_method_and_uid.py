# backfill_pl_method_and_uid.py
import argparse, os, json, hashlib, glob, sys
import pandas as pd
from datetime import datetime

REQUIRED_META_KEYS = ("Glycan Type", "Mass Analyzer charge mode")

def short_id(s: str) -> str:
    return hashlib.sha1((s or "").encode("utf-8")).hexdigest()[:8]

def sha1_file(path):
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

def looks_like_metadata(d: dict) -> bool:
    return isinstance(d, dict) and all(k in d for k in REQUIRED_META_KEYS)

def try_load_json(p):
    try:
        with open(p, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None

def infer_experiment_title(root, exp_json_path=None, metadata_json=None, fallback="Unassigned"):
    # priority: explicit .exp.json -> metadata.json -> folder name
    if exp_json_path and os.path.exists(exp_json_path):
        exp = try_load_json(exp_json_path)
        if isinstance(exp, dict) and exp.get("experiment"):
            return str(exp["experiment"])
    if metadata_json and looks_like_metadata(metadata_json):
        t = str(metadata_json.get("Experiment Title") or "").strip()
        if t:
            return t
    return os.path.basename(os.path.abspath(root)) or fallback

def infer_sample_name_from_filename(stem: str) -> str:
    # strip common tags
    for tag in ("_trainable_", "_trainable", "_pseudolabels_", "_pseudolabels", "_unlabeled_", "_unlabeled"):
        pos = stem.lower().find(tag)
        if pos > 0:
            return stem[:pos]
    return stem

def ensure_uid(df: pd.DataFrame, exp_title: str, sample_name: str) -> pd.DataFrame:
    # Accept both MS2scan_no and MS2Scan_no (plus common alternates)
    scan_cols = [c for c in df.columns if c.lower() in {"ms2scan_no", "ms2scan", "scan", "scannum"} or c == "MS2Scan_no"]
    scan_col = scan_cols[0] if scan_cols else None
    if scan_col is None:
        return df  # cannot make UID without a scan number
    exp_id    = short_id(exp_title)
    sample_id = short_id(sample_name)
    s = df[scan_col].astype(str).str.replace(r"\D", "", regex=True)
    s = s.str.zfill(6).where(s.str.len() > 0, other="000000")
    df = df.copy()
    if "UID" not in df.columns:
        df["UID"] = s.map(lambda x: f"{exp_id}:{sample_id}:{x}")
    return df

def build_method(experiment_title, sample_name, parents, ion_path, ion_sheet, ppm):
    payload = {
        "version": "1.0",
        "dataset_type": "pseudolabel",
        "experiment_title": experiment_title,
        "sample_name": sample_name,
        "parents": parents,
        "ionlist": {
            "source": "external",
            "path": ion_path,
            "sheet": ion_sheet or "",
            "sha1": sha1_file(ion_path) if ion_path and os.path.exists(ion_path) else "",
            "ppm_tolerance": ppm,
        },
        "uid": {
            "namespace": "GLYCOMSP",
            "rule": "UID = exp_id:sample_id:MS2scan_no (zfill 6)"
        },
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds")
    }
    return payload

def main():
    ap = argparse.ArgumentParser(description="Backfill PL method.json and UID for trainable/unlabeled CSVs.")
    ap.add_argument("--root", required=True, help="Folder to scan (recurses).")
    ap.add_argument("--ionlist", required=True, help="Ion list CSV/XLSX used for PL.")
    ap.add_argument("--ion-sheet", default="", help="Sheet name if ion list is Excel.")
    ap.add_argument("--ppm", type=int, default=20, help="PPM tolerance for ion matching.")
    ap.add_argument("--exp-json", default="", help="Optional .exp.json to infer experiment title.")
    ap.add_argument("--write-suffix", default="_uidfixed", help="Suffix for rewritten CSV (if UID added).")
    args = ap.parse_args()

    root = args.root
    ion_path = args.ionlist
    if not os.path.exists(root):
        print(f"[ERR] root not found: {root}")
        sys.exit(2)
    if not os.path.exists(ion_path):
        print(f"[ERR] ionlist not found: {ion_path}")
        sys.exit(2)

    # Try to load any metadata JSON in the root to help infer exp title
    meta_json = None
    for cand in glob.glob(os.path.join(root, "*.json")):
        d = try_load_json(cand)
        if looks_like_metadata(d):
            meta_json = d
            break

    experiment_title = infer_experiment_title(root, args.exp_json, meta_json)

    # Find candidate CSVs
    csvs = []
    for pat in ("*trainable*.csv", "*_trainable.csv", "*_unlabeled*.csv"):
        csvs.extend(glob.glob(os.path.join(root, "**", pat), recursive=True))
    if not csvs:
        print("[INFO] no trainable/unlabeled csvs found; nothing to do.")
        return

    for csv_path in sorted(set(csvs)):
        stem = os.path.splitext(os.path.basename(csv_path))[0]
        sample_name = infer_sample_name_from_filename(stem)

        # Nearby pseudolabels (optional)
        tsv = ""
        for pat in (f"{sample_name}*pseudolabel*.tsv", f"{sample_name}*pseudolabel*.csv", "*pseudolabel*.tsv", "*pseudolabel*.csv"):
            hits = glob.glob(os.path.join(os.path.dirname(csv_path), pat))
            if hits:
                tsv = hits[0]
                break

        # Nearby converted CSV (optional)
        conv = ""
        for pat in (f"ms2_{sample_name}*.csv", f"{sample_name}*.raw.csv", "*ms2_*.csv"):
            hits = glob.glob(os.path.join(os.path.dirname(csv_path), pat))
            if hits:
                conv = hits[0]
                break

        # Backfill UID
        df = pd.read_csv(csv_path)
        had_uid = "UID" in df.columns
        df2 = ensure_uid(df, experiment_title, sample_name)
        if ("UID" not in df.columns) and ("UID" in df2.columns):
            out_csv = os.path.join(os.path.dirname(csv_path), stem + args.write_sufix + ".csv") if hasattr(args, "write_sufix") else os.path.join(os.path.dirname(csv_path), stem + args.write_suffix + ".csv")
            df2.to_csv(out_csv, index=False)
            print(f"[UID] wrote {out_csv}")
            csv_final = out_csv
        else:
            csv_final = csv_path
            print(f"[UID] ok: {os.path.basename(csv_path)} (had_uid={had_uid})")

        # Write method.json
        parents = {
            "converted_csv": conv or "(unknown)",
            "pseudolabels_tsv": tsv or "(unknown)",
            "trainable_csv": csv_final
        }
        method = build_method(experiment_title, sample_name, parents, ion_path, args.ion_sheet, args.ppm)
        method_path = os.path.join(os.path.dirname(csv_final), f"{sample_name}.method.json")
        with open(method_path, "w", encoding="utf-8") as f:
            json.dump(method, f, indent=2, ensure_ascii=False)
        print(f"[METHOD] wrote {method_path}")

if __name__ == "__main__":
    main()



