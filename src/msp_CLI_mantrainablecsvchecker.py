import pandas as pd
df = pd.read_csv("<<<your_merged_trainable.csv>>>")

# 1) Label counts
print(df["Structure"].value_counts().head(20))
print("Non-glycan total:", (df["Structure"]=="Non-glycan").sum())

# 2) No duplicate scan ids
uid_col = "unique_ID" if "unique_ID" in df.columns else "MS2scan_no"
print("Duplicate unique_IDs:", df[uid_col].duplicated().sum())

# 3) How sparse are negatives w.r.t. ion features?
meta = {"Structure","IUPACname(optional)","Glycanannotation2","protonatedmass",uid_col}
feat_cols = [c for c in df.columns if c not in meta]
neg_hits = (df.loc[df.Structure=="Non-glycan", feat_cols] > 0).sum(axis=1)
print("Neg hit-count summary:\n", neg_hits.describe())  # should be mostly < 3 by design

# 4) Compare positive vs negative average feature activity
pos_hits = (df.loc[df.Structure!="Non-glycan", feat_cols] > 0).sum(axis=1)
print("Avg pos hits:", pos_hits.mean(), " | Avg neg hits:", neg_hits.mean())

# 5) Spot-check a few Non-glycan rows (features aren’t all zeros)
print(df.loc[df.Structure=="Non-glycan", feat_cols].head(3))