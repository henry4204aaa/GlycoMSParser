# prepare_and_export_dataset.py

import os
import json
import pandas as pd
from datetime import datetime

import os
import json
from datetime import datetime
from tkinter import filedialog, messagebox


from .prepare_dataset_module import (
    prepare_dataset,
    validate_csv_structure,
    validate_annotation_structure,
    extract_ion_list,
    add_unique_id_column
)

method_save_paths = {}  # New: Maps experiment → last-used method file path

# --------------------------
# Utility: Load Metadata
# --------------------------

def load_metadata(metadata_path):
    with open(metadata_path, "r") as f:
        metadata = json.load(f)
    return metadata

# --------------------------
# Main Wrapper Function
# --------------------------

def prepare_and_export_dataset(csv_path, excel_path, ion_sheet_name, metadata_path, output_dir, sample_name=None):
    # Load files
    csv_df = pd.read_csv(csv_path)
    excel_data = pd.read_excel(excel_path, sheet_name=None)
    metadata = load_metadata(metadata_path)

    # Extract required sheets
    anno_df = excel_data.get("MSlist")
    ion_df = excel_data.get(ion_sheet_name)

    if anno_df is None or ion_df is None:
        raise ValueError("Excel file must contain 'MSlist' and the specified ion list sheet")

    # Validate structures
    validate_csv_structure(csv_df)
    validate_annotation_structure(anno_df)

    # Determine glycan type from metadata
    glycan_type = metadata.get("Glycan Type", "OG")

    # Run main merge logic
    output_df = prepare_dataset(anno_df, csv_df, ion_df, glycan_type)

    # Apply optional ID prefix
    prefix = sample_name or os.path.splitext(os.path.basename(csv_path))[0]
    output_df = add_unique_id_column(output_df, prefix=prefix)

    # Save merged output
    timestamp = datetime.now().strftime("%Y%m%d")
    out_filename = f"{prefix}_merged_{timestamp}.csv"
    out_path = os.path.join(output_dir, out_filename)
    output_df.to_csv(out_path, index=False)

    # Save updated method file
    method = {
        "experiment": metadata.get("Experiment Title", "Unassigned"),
        "generated_on": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "samples": {
            prefix: {
                "csv": os.path.basename(csv_path),
                "excel": os.path.basename(excel_path),
                "metadata": os.path.basename(metadata_path),
                "raw_file": metadata.get("Raw filename", "not linked"),
                "validated": True,
                "ion_sheet": ion_sheet_name
            }
        }
    }

    exp_name = method["experiment"]
    if exp_name in method_save_paths:
        method_path = method_save_paths[exp_name]
    else:
        method_path = os.path.join(output_dir, f"{prefix}.method.json")
        method_save_paths[exp_name] = method_path

    with open(method_path, "w") as f:
        json.dump(method, f, indent=4)

    return out_path, method_path, output_df
