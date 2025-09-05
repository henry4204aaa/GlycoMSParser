import mspvalidator_merger as mspval
import pandas as pd

excel = "path/to/your_manual_annotation.xlsx"   # has MSlist & ionlist
raw_tsv = "path/to/your_raw_converted.csv"      # tab-separated from extractor
deriv = "PerMe(Reduced)"                        # or value from metadata dialog

# 1) Parse + align manual annotations with peaks
pre_df, ion_index, ion_df = mspval.directassign_files(excel, raw_tsv, deriv, debug=False)

# 2) (Optional) keep N-glycan only using the helper column carried through
pre_df_N = pre_df[pre_df["Glycanannotation2"].eq("N")].copy()

# 3) Emit a trainable wide CSV (features + labels)
out_csv = "my_experiment_N_only_merged.csv"
mspval.createnormailzedionlistcsv(ion_index, pre_df_N, ion_df, out_csv)
print("Saved:", out_csv)