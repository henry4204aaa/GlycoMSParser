import joblib, pandas as pd, numpy as np
#fix the model prediction using (0, 1, 2, 3, 4, 5) as structure label and decode to proper string
#use this as a salvage method if you mess up the model
enc = joblib.load("W:\GlycoMSP_temp\MLgeneratetest1_score007_salvage_PLabeled_labelencoder_20test0val_balanced_realcap_tau062_margin008.joblib")
df  = pd.read_csv("G:\其他電腦\My Computer\GlycoMSParser\ms2_zebrafish_Yann_glycome_20240922_zf_sPerMeNG_intestine.raw_unlabeled_ppm20_predicted_pseudo.csv")

# ints -> original objects (likely tuples)
y = enc.inverse_transform(df["Predicted_Label"].astype(float).astype(int).to_numpy())

# canonicalize to strings
def canonicalize_label(lbl):
    if isinstance(lbl, str):
        s = lbl.strip()
        return "Non-glycan" if s.lower().startswith("non") else s
    if isinstance(lbl, (list, tuple)) and lbl and isinstance(lbl[0], (list, tuple)) and len(lbl[0]) == 2:
        lbl = dict(lbl)
    if isinstance(lbl, dict):
        order = ["F","H","N","G","KDN","S"]
        parts = [f"{k}{int(lbl.get(k,0))}" for k in order if int(lbl.get(k,0)) > 0]
        return "Non-glycan" if not parts else "".join(parts)
    if isinstance(lbl, (list, tuple)) and all(isinstance(x, (int, float)) for x in lbl):
        keys = ["F","H","N","G","KDN","S"]
        parts = [f"{k}{int(v)}" for k, v in zip(keys, lbl) if int(v) > 0]
        return "Non-glycan" if not parts else "".join(parts)
    return str(lbl)

df["Predicted_Label"] = [canonicalize_label(v) for v in y]
df.to_csv("W:\GlycoMSP_temp\MLgeneratetest1_score007_salvage_PLabeled_labelencoder_20test0val_balanced_realcap_tau062_margin008_ms2predicted_pseudo_decoded.csv", index=False)
print("Wrote *_decoded.csv with canonical string labels.")