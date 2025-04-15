# prepare_dataset_module.py

import pandas as pd
import ast
import os
import numpy as np
from mspvalidator_merger import calc_protonated_mass, calc_observed_mass  # Assuming you're using a separate module

# ---------------------
# Validation Functions
# ---------------------

def validate_csv_structure(csv_df):
    required_columns = {"MS2scan_no", "peaklist", "peakintensity"}
    if not required_columns.issubset(csv_df.columns):
        raise ValueError(f"CSV file missing required columns: {required_columns - set(csv_df.columns)}")



def validate_annotation_structure(anno_df):
    required_columns = {"MS2scan_no", "Stripped_sequence", "Glycan_structure", "Adduct", "Charge"}
    if not required_columns.issubset(anno_df.columns):
        raise ValueError(f"Annotation sheet missing required columns: {required_columns - set(anno_df.columns)}")

# ------------------------
# Glycan Type Setup
# ------------------------

def typeselection(glycan_type):
    if glycan_type == "NG":
        return 203.07937, 1.007276, 0.0
    elif glycan_type == "OG":
        return 203.07937, 1.007276, 18.01056
    else:
        raise ValueError("Glycan type must be 'NG' or 'OG'")

# ------------------------
# Data Preparation
# ------------------------

def add_mass_columns(df, glycan_type):
    derivatization_mass, proton_mass, water_mass = typeselection(glycan_type)
    df["protonated_mass"] = df.apply(lambda row: calc_protonated_mass(row, derivatization_mass), axis=1)
    df["observed_mass"] = df.apply(lambda row: calc_observed_mass(row, derivatization_mass, water_mass, proton_mass), axis=1)
    return df

# ------------------------
# Merge Logic
# ------------------------

def extractannotation(anno_df, csv_df):
    return pd.merge(anno_df, csv_df, on="MS2scan_no", how="inner")

def expand_peaklists(df):
    df["peaklist"] = df["peaklist"].apply(ast.literal_eval)
    df["peakintensity"] = df["peakintensity"].apply(ast.literal_eval)
    return df

# ------------------------
# Ion Matching
# ------------------------

def extract_ion_list(ion_df):
    return ion_df["mz"].tolist(), ion_df.columns[1:].tolist()  # ion masses, labels

def match_ions(peaklist, peakintensity, ion_masses, ppm=20):
    matched = []
    for target in ion_masses:
        found = False
        for mz, inten in zip(peaklist, peakintensity):
            if abs(mz - target) / target * 1e6 <= ppm:
                matched.append(inten)
                found = True
                break
        if not found:
            matched.append(0)
    return matched

# ------------------------
# Final Dataset Generator
# ------------------------

def prepare_dataset(anno_df, csv_df, ion_df, glycan_type="NG"):
    validate_csv_structure(csv_df)
    validate_annotation_structure(anno_df)
    combined = extractannotation(anno_df, csv_df)
    combined = add_mass_columns(combined, glycan_type)
    combined = expand_peaklists(combined)

    ion_masses, ion_labels = extract_ion_list(ion_df)

    feature_vectors = combined.apply(
        lambda row: match_ions(row["peaklist"], row["peakintensity"], ion_masses), axis=1
    )
    feature_matrix = pd.DataFrame(feature_vectors.tolist(), columns=ion_labels)

    output_df = pd.concat([combined[["MS2scan_no", "Stripped_sequence", "Glycan_structure"]], feature_matrix], axis=1)
    return output_df

# ------------------------
# Optional Postprocessing
# ------------------------

def add_unique_id_column(df, prefix="sample"):
    df = df.copy()
    df.insert(0, "unique_ID", [f"{prefix}_{i+1}" for i in range(len(df))])
    return df
