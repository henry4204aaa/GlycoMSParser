# normalize_paths_in_json.py
import json, os, sys, glob
from pathcanon import to_posix_str, strip_trailing_sep

def fix_dict(d: dict) -> bool:
    changed = False
    if "parents" in d and isinstance(d["parents"], dict):
        for k, v in list(d["parents"].items()):
            if isinstance(v, str) and v not in ("(unknown)", ""):
                newv = to_posix_str(strip_trailing_sep(v))
                if newv != v:
                    d["parents"][k] = newv; changed = True
    if "ionlist" in d and isinstance(d["ionlist"], dict):
        v = d["ionlist"].get("path")
        if isinstance(v, str) and v:
            newv = to_posix_str(strip_trailing_sep(v))
            if newv != v:
                d["ionlist"]["path"] = newv; changed = True
    # exp.json shape: samples → fields like csv/excel/json/unlabeled_dataset
    if "samples" in d and isinstance(d["samples"], dict):
        for s, info in d["samples"].items():
            if isinstance(info, dict):
                for key in ("csv","excel","json","unlabeled_dataset","trainable_csv","predicted_csv"):
                    v = info.get(key)
                    if isinstance(v, str) and v:
                        newv = to_posix_str(strip_trailing_sep(v))
                        if newv != v:
                            info[key] = newv; changed = True
    return changed

def main(root):
    paths = glob.glob(os.path.join(root, "**", "*.json"), recursive=True)
    fixed = 0
    for p in paths:
        try:
            with open(p, "r", encoding="utf-8") as f:
                d = json.load(f)
            if isinstance(d, dict) and fix_dict(d):
                with open(p, "w", encoding="utf-8") as f:
                    json.dump(d, f, indent=2, ensure_ascii=False)
                print("[fixed]", p); fixed += 1
        except Exception as e:
            print("[skip]", p, "->", e)
    print("Done. Files fixed:", fixed)

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")