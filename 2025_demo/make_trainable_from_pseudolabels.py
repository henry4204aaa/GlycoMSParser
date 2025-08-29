# make_trainable_from_pseudolabels.py
import pandas as pd
import importlib.util

# --- load utils (adjust the path if needed) ---
UTIL_PATH = r".\2025_demo\pseudolabel_utils.py"
spec = importlib.util.spec_from_file_location("pseudolabel_utils", UTIL_PATH)
plu = importlib.util.module_from_spec(spec); spec.loader.exec_module(plu)
to_training_dataset = plu.to_training_dataset

# --- inputs (edit paths as needed) ---
SURVIVORS_CSV = r"G:\其他電腦\My Computer\GlycoMSParser\zf_intestine_pseudolabels_filtered_score0p1.csv"
ANNOT_XLSX    = r"G:\其他電腦\My Computer\GlycoMSParser\src\20240922_temp_zf_intestine_1.xlsx"

# --- params ---
ION_SHEET     = "ionlist"
PPM           = 20.0
INCLUDE_MASS  = False   # set True if you want protonatedmass/ppm_error columns included

# --- run ---
df = pd.read_csv(SURVIVORS_CSV)
train_df = to_training_dataset(df, ionlist_excel_path=ANNOT_XLSX, ion_sheet=ION_SHEET, ppm=PPM, include_mass_cols=INCLUDE_MASS)

OUT = SURVIVORS_CSV.replace(".csv", "_trainable.csv")
train_df.to_csv(OUT, index=False)
print(f"Trainable saved -> {OUT} ; rows={len(train_df)}, features={train_df.shape[1]-3 if not INCLUDE_MASS else train_df.shape[1]-6}")