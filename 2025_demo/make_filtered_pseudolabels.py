import pandas as pd
import importlib.util

# --- load pseudolabel_utils.py even if it's not installed ---
UTIL_PATH = r".\2025_demo\pseudolabel_utils.py"  # <- change if located elsewhere
spec = importlib.util.spec_from_file_location("pseudolabel_utils", UTIL_PATH)
plu = importlib.util.module_from_spec(spec); spec.loader.exec_module(plu)
filter_pseudolabels = plu.filter_pseudolabels

# --- paths ---
IN_TSV  = r"G:\其他電腦\My Computer\GlycoMSParser\zf_sPerMeNG_intestine_pseudolabels_20250825_freeend.tsv"
OUT_CSV = r"G:\其他電腦\My Computer\GlycoMSParser\zf_intestine_pseudolabels_filtered_score0p1.csv"

# --- params (your requested defaults with score=0.1) ---
MIN_SCORE   = 0.1
MIN_HITS    = 2
MAX_PPM     = 20.0
TOP_N       = 1

# --- stream + filter ---
chunks = []
for chunk in pd.read_csv(IN_TSV, sep="\t", chunksize=50000):
    chunks.append(
        filter_pseudolabels(
            chunk,
            min_ion_score=MIN_SCORE,
            min_hit_count=MIN_HITS,
            max_ppm=MAX_PPM,
            keep_top_n=TOP_N,
        )
    )
filtered = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
filtered.to_csv(OUT_CSV, index=False)
print(f"{len(filtered)} rows -> {OUT_CSV}")


"""
# START with a long-form table that has one row per (scan, composition) candidate
# BEFORE any collapsing:
# cols needed (example): MS2scan_no, composition, ion_score, neutral_mass, z, manual_tuple

def isotope_neutral_variants(neutral_mass):
    return [("M", neutral_mass),
            ("M+1", neutral_mass + 1.003355),
            ("M-1", neutral_mass - 1.003355)]

PPM = 30  # example
SCORE_MIN = 0.07
BIG_MASS = 2500.0

def within_ppm(mass, target):
    return abs(mass - target) / target * 1e6

# assume df_raw has one row per (scan, composition) w/ candidate neutral mass from comp
rows = []
for _, r in df_raw.iterrows():
    scan = r.MS2scan_no
    comp = r.composition  # tuple or string
    cand_mass = r.comp_neutral_mass
    neutral = r.neutral_mass  # precursor neutral (from spectrum, i.e., (mz - proton_mass)*z)
    score = r.ion_score

    windows = [("M", neutral)]
    if neutral >= BIG_MASS:
        windows = isotope_neutral_variants(neutral)

    hit = False
    iso_flag = None
    best_ppm = None
    for label, target in windows:
        ppm = within_ppm(cand_mass, target)
        if ppm <= PPM:
            hit = True
            iso_flag = label
            best_ppm = ppm if (best_ppm is None or ppm < best_ppm) else best_ppm
    if not hit:
        continue  # drop by mass gate

    rows.append({
        "MS2scan_no": scan,
        "composition": comp,
        "ion_score": score,
        "delta_mass_ppm": best_ppm,
        "is_isotope_window": iso_flag,
        # optional: n_ions_matched, sum_intensity_matched
    })

df_candidates = pd.DataFrame(rows)

# dedupe only within (scan, composition), keep highest score
df_candidates.sort_values(["MS2scan_no","composition","ion_score"], ascending=[True,True,False], inplace=True)
df_candidates = df_candidates.drop_duplicates(subset=["MS2scan_no","composition"], keep="first")

# SCORE FILTER (tunable)
df_candidates = df_candidates[df_candidates["ion_score"] >= SCORE_MIN].copy()

# JOIN manual for comparison WITHOUT collapsing candidates
df_out = df_candidates.merge(df_manual[["MS2scan_no","manual_tuple"]],
                             on="MS2scan_no", how="left")
df_out["tuple_match"] = (df_out["composition"].astype(str) == df_out["manual_tuple"].astype(str))
df_out["manual_check"] = df_out["ion_score"].notna()  # same semantics you used

# OPTIONAL: per-scan summary (separate artifact)
agg = (df_out
       .assign(score_rank=df_out.groupby("MS2scan_no")["ion_score"].rank("dense", ascending=False))
       .sort_values(["MS2scan_no","score_rank"])
       .groupby("MS2scan_no")
       .agg(
           manual_tuple=("manual_tuple","first"),
           any_match=("tuple_match","any"),
           best_candidate=("composition","first"),
           best_score=("ion_score","max"),
           best_is_isotope_window=("is_isotope_window","first"),
           n_candidates=("composition","size"),
       )
       .reset_index())
"""