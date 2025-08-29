import re, pandas as pd, numpy as np

ORDER = ["H","N","S","G","KDN","F"]
ALIASES = {
    "H":"H","HEX":"H","N":"N","HEXNAC":"N","S":"S","NEU5AC":"S","NEUAC":"S","SA":"S",
    "G":"G","NEU5GC":"G","NEUGC":"G","KDN":"KDN","K":"KDN","F":"F","FUC":"F"
}
token_re = re.compile(r"([A-Za-z]+)\s*([+-]?\d+)", re.IGNORECASE)

def manual_to_tuple(s):
    c = {k:0 for k in ORDER}
    if not isinstance(s,str) or not s.strip(): return None
    for m in token_re.finditer(s.replace(" ","").upper()):
        k = ALIASES.get(m.group(1).upper()) or ALIASES.get(m.group(1).upper().replace("5","").replace("AC",""))
        if k: c[k] += int(m.group(2))
    return tuple(c[k] for k in ORDER)

def pseudo_to_tuple(s):
    if s is None or (isinstance(s,float) and np.isnan(s)): return None
    try:
        vals = list(eval(str(s), {"__builtins__":{}}))
    except Exception:
        return None
    vals = [int(float(x)) if pd.notna(x) else 0 for x in vals] + [0]*6
    return tuple(vals[:6])

# === paths ===
SURVIVORS = r"G:\其他電腦\My Computer\GlycoMSParser\zf_intestine_pseudolabels_filtered_score0p1.csv"
MANUAL    = r"G:\其他電腦\My Computer\GlycoMSParser\src\20240922_temp_zf_intestine_1.xlsx"
OUT       = SURVIVORS.replace(".csv", "_vs_manual_tuple_by_scan.csv")

# Load
pseudo = pd.read_csv(SURVIVORS)[["MS2scan_no","composition"]].copy()
pseudo["MS2scan_no"] = pd.to_numeric(pseudo["MS2scan_no"], errors="coerce").astype("Int64")
pseudo["pseudo_tuple"] = pseudo["composition"].apply(pseudo_to_tuple)

ms = pd.ExcelFile(MANUAL).parse("MSlist")[["Structure","MS2scan_no"]].dropna(subset=["MS2scan_no"])
ms["MS2scan_no"] = pd.to_numeric(ms["MS2scan_no"], errors="coerce").astype("Int64")
ms["manual_tuple"] = ms["Structure"].apply(manual_to_tuple)

# Merge and score
j = ms.merge(pseudo, on="MS2scan_no", how="left")
has_pseudo = j["pseudo_tuple"].apply(lambda x: isinstance(x, tuple))
is_match   = j.apply(lambda r: isinstance(r["manual_tuple"], tuple) and isinstance(r["pseudo_tuple"], tuple)
                               and r["manual_tuple"] == r["pseudo_tuple"], axis=1)

TP = int(is_match.sum())
FN = int((~has_pseudo).sum())
FP = int((has_pseudo & (~is_match)).sum())

print(f"Same-scan tuple match — TP={TP}, FP={FP}, FN={FN}")
j["tuple_match"] = is_match
j.to_csv(OUT, index=False)
print(f"Saved -> {OUT}")