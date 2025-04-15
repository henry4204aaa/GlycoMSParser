import os
version = "0.91"
last_update = 20250415
import msprawextractor as mspext
import threading
from tkinter import ttk
import time #for testing
import random #for random failure
#import datetime
from datetime import datetime
import configparser #load config
import json
import tkinter.simpledialog as simpledialog
from tkinter import filedialog
import shutil
import traceback
import mspvalidator_merger as mspval
import pandas as pd
# v0.99: demo version before cleaning code
# v0.91: ML added. Lacking combining data and include "Non-glycan labels for training"
# v0.9: add ML window
# v0.81: adding derivatization flag
# v0.8: able to export merged trainable datasets
# v0.65: introduce real validation function
# v0.64c: pseudo-link (no real validation) but link the files
# v0.64b: load method file and skip unneeded file conversion for filling missing metadata part
# v0.64a: record current assignment and link information (in one method file or update json and keep sample json in method?)
# v0.63: dataset management + link and validate data (placeholder)
# v0.61: build functional dataset management window (if new version failed, rollback to this)
# v0.6: support saving converted files by selecting user defined folders
# for trace errors
# v0.58: UI refine. (v0.57 works perfectly)
# v0.56: trace errors and also delete temp csv while metadata assingment and rename csv process failed
# v0.56  change the logic to avoid main window freeze during conversion
# v0.55: refine the logic. Go to simulation before real run.
# version: 0.54
# add persistent metadata edit window
# version: 0.53
# add metadata class
# version: 0.52
# add quit destroy function
# version: 0.51
# add logger
# test of adding tkinter windows with the help of GPT4o
# add GUI, add progress bar and making buttons work

# version: 0.4
# added mzml support function (waiting mzml extractor implementation - will be noted as v0.5 then, while MSP will be v0.6)
# version: 0.3
# added excel reader for annotated file. Will add validation here or in annotationreader.py
# version: 0.2
# date: 20240730
# about this file: original trytolistoutheaders.py and GlycoMSP_demo.py modulated to support future standardized processing workflow
# import tkinter part
#mind that we plan to add GUI support in project managing, the tkinter detection may be moved to funtion: GUI in MSPinit.py in future
# re-organized the code structure and readability using copilot and ChatGPT4-o
# no sensitive contents were sent to the server for this part

try: 
    import tkinter as tk
    from tkinter import filedialog, messagebox
except ImportError:
    raise ImportError('Please install tkinter to enable GUI-based file selection')
else:
    has_tkinter = True

selected_files = {}  # Dictionary to keep track of selected files


#logger
class AppLogger:
    def __init__(self):
        self.entries = []
        self.debug_logs = []
        self.gui_writer = None  # Optional live display hook

    def log(self, message):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        entry = f"[{timestamp}] {message}"
        print(entry)
        self.entries.append(entry)
        if self.gui_writer:
            self.gui_writer(entry)

    def debug(self, message):
        self.debug_logs.append(f"[DEBUG] {message}")

    def save(self, filename, include_debug=False):
        with open(filename, "w", encoding="utf-8") as f: #add utf-8 encode to get arrow supported
            for line in self.entries:
                f.write(line + "\n")
            if include_debug and self.debug_logs:
                f.write("-----\ndebug enabled-----\n")
                for line in self.debug_logs:
                    f.write(line + "\n")
                f.write("-----\ndebug log ends-----\n")
        return filename

logger = AppLogger()


#metadata class 

class MetadataEditorWindow:
    def __init__(self, parent, raw_file_list, on_each_metadata_ready_callback, on_finish=None, output_dir=None, skip_conversion=False):
        self.parent = parent
        self.raw_file_list = raw_file_list
        self.callback = on_each_metadata_ready_callback
        self.on_finish = on_finish
        self.current_index = 0
        self.last_metadata = None
        self.metadata = {}
        self.output_dir = output_dir
        self.skip_conversion = skip_conversion
        self.entries = {}
        fields = [
            "Experiment Title", "Experiment Description", "Author running this analysis",
            "Raw data acquired date", "Glycan Type", "Mass Analyzer charge mode", "Derivatization Type"
        ]

        self.window = tk.Toplevel(parent)
        self.window.title("Enter Experiment Metadata")
        self.window.geometry("600x500")
        # Row 0: Status label placeholder (spans both columns)
        self.status_label = tk.Label(self.window, text="", fg="blue")
        self.status_label.grid(row=0, column=0, columnspan=2, pady=(10, 5))

        if not self.skip_conversion:
            self.load_file(self.raw_file_list[self.current_index])

        for i, field in enumerate(fields, start=1):
            label = tk.Label(self.window, text=field)
            label.grid(row=i, column=0, padx=10, pady=5, sticky="w")

            if field == "Mass Analyzer charge mode":
                var = tk.StringVar()
                entry = ttk.Combobox(self.window, textvariable=var, values=["+", "-"], state="readonly", width=10)
                entry.current(0)
            elif field == "Derivatization Type":
                var = tk.StringVar()
                entry = ttk.Combobox(
                    self.window, textvariable=var,
                    values=["PerMe(Freeend)", "PerMe(Reduced)", "Others"],  # Add more if needed
                    state="readonly", width=20
                )
                entry.current(0)

            elif field == "Glycan Type":
                var = tk.StringVar()
                entry = ttk.Combobox(self.window, textvariable=var, values=["N", "O", "N+O"], state="readonly", width=10) #add the values when going to support new glycan types
                entry.current(0)
            else:
                entry = tk.Entry(self.window, width=50)

            entry.grid(row=i, column=1, padx=10, pady=5)
            self.entries[field] = entry

        # Buttons
        button_frame = tk.Frame(self.window)
        button_row = len(fields) + 1 #so the button adjusted itself in future if we add more contents into metadata
        button_frame.grid(row=button_row, column=0, columnspan=2, pady=10)
        tk.Button(button_frame, text="Import Metadata", command=self.import_metadata).grid(row=0, column=0, padx=5)
        tk.Button(button_frame, text="Use Last Metadata", command=self.use_last_metadata).grid(row=0, column=1, padx=5)
        tk.Button(button_frame, text="Clear", command=self.clear_fields).grid(row=0, column=2, padx=5)
        tk.Button(button_frame, text="Generate Dataset", command=self.generate).grid(row=0, column=3, padx=5)
        tk.Button(button_frame, text="Cancel", command=self.window.destroy).grid(row=0, column=4, padx=5)


        label = tk.Label(self.window, text="Raw File (optional):")
        label.grid(row=button_row + 1, column=0, padx=10, pady=5, sticky="w")

        raw_file_frame = tk.Frame(self.window)
        raw_file_frame.grid(row=button_row + 1, column=1, padx=10, pady=5, sticky="w")

        self.raw_file_entry = tk.Entry(raw_file_frame, width=40)
        self.raw_file_entry.pack(side=tk.LEFT)

        def select_raw_file():
            path = filedialog.askopenfilename(title="Select RAW File", filetypes=[("RAW files", "*.raw")])
            if path:
                self.raw_file_entry.delete(0, tk.END)
                self.raw_file_entry.insert(0, path)

        tk.Button(raw_file_frame, text="Browse", command=select_raw_file).pack(side=tk.LEFT, padx=5)


    def on_conversion_complete(self):
        self.status_label.config(text="Conversion complete. Ready for metadata.", fg="green")

        # Extract filename base for use in renaming
        self.rawfilename = os.path.splitext(os.path.basename(self.current_raw_file))[0]

        # Re-enable entries and buttons
        for entry in self.entries.values():
            entry.config(state="normal")
        for child in self.window.winfo_children():
            if isinstance(child, tk.Button):
                child.config(state="normal")

        # Pre-fill if last used
        if self.last_metadata:
            self.fill_fields_from_metadata(self.last_metadata)

        self.window.title(f"Metadata for: {os.path.basename(self.current_raw_file)}")

    def load_file(self, raw_file):
        self.current_raw_file = raw_file
        self.window.title(f"Converting: {os.path.basename(raw_file)}")
        self.clear_fields()

        # Add status label if not already there
        if not hasattr(self, "status_label"):
            self.status_label = tk.Label(self.window, text="", fg="blue")
            self.status_label.grid(row=0, columnspan=2, pady=10)
        self.status_label.config(text="Converting raw file... please wait.", fg="blue")

        # Disable form and buttons temporarily
        for entry in self.entries.values():
            entry.config(state="disabled")
        for child in self.window.winfo_children():
            if isinstance(child, tk.Button):
                child.config(state="disabled")

        def background_conversion():
            if not self.skip_conversion:
                try:
                    success = mspext.convert_raw_to_csv(raw_file, debug=True)
                    if success:
                        self.window.after(0, self.on_conversion_complete)
                    else:
                        # Clean up temp files
                        rawfilename = os.path.splitext(os.path.basename(raw_file))[0]
                        temp_ms2 = f"ms2tmp_{rawfilename}.csv"
                        temp_ms3 = f"ms3tmp_{rawfilename}.csv"

                        for temp_file in [temp_ms2, temp_ms3]:
                            if os.path.exists(temp_file):
                                try:
                                    os.remove(temp_file)
                                    logger.log(f"[Cleanup] Removed leftover temp file: {temp_file}")
                                except Exception as e:
                                    logger.log(f"[WARNING] Failed to delete temp file: {temp_file} — {e}")
                        err_msg = str(success) if isinstance(success, Exception) else "Unknown error"
                        self.window.after(0, lambda: messagebox.showerror("Conversion Failed", f"The raw file could not be converted:\n{err_msg}"))
                        traceback.print_exc()
                except Exception as e:
                    logger.log(f"[FATAL] Exception during raw file conversion: {e}")
                    self.window.after(0, lambda: messagebox.showerror("Thermo Library Error",f"Raw file could not be processed.\n\nDetails:\n{e}"))
        threading.Thread(target=background_conversion, daemon=True).start()

    def import_metadata(self):
        path = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
        if path:
            with open(path, "r") as f:
                self.metadata = json.load(f)
            self.fill_fields_from_metadata(self.metadata)

    def use_last_metadata(self):
        if self.last_metadata:
            self.fill_fields_from_metadata(self.last_metadata)
        else:
            messagebox.showinfo("No previous metadata", "You haven't submitted any metadata yet.")

    def fill_fields_from_metadata(self, meta):
        for field, entry in self.entries.items():
            entry.delete(0, tk.END)
            entry.insert(0, meta.get(field, ""))

    def clear_fields(self):
        for entry in self.entries.values():
            entry.delete(0, tk.END)

    def generate(self):
        if self.output_dir is None:
            selected_dir = filedialog.askdirectory(
                title="Select folder to save metadata and output files",
                initialdir=os.getcwd()
            )
            if not selected_dir:
                messagebox.showwarning("Cancelled", "You must select an output folder.")
                return
            self.output_dir = selected_dir
            logger.log(f"[Metadata] Output will be saved to: {self.output_dir}")

        for field, entry in self.entries.items():
            self.metadata[field] = entry.get()

        extractdate = time.strftime("%Y%m%d")
        self.metadata["Parameters when GlycoMSP launched"] = ["Autofilled by GUI"]
        self.metadata["Date of file extracted from raw file"] = extractdate
        # --- Safely resolve raw file path ---
        raw_path = None
        # Priority: use field if present (linking mode)
        if hasattr(self, "raw_file_entry"):
            raw_path = self.raw_file_entry.get().strip()
        # Fallback: use internal attribute if available (conversion mode)
        if not raw_path:
            raw_path = getattr(self, "current_raw_file", None)
        # Default to not linked
        raw_path = raw_path or "not linked"
        # Save into metadata
        self.metadata["Original raw file path"] = raw_path
        self.metadata["Raw filename"] = os.path.splitext(os.path.basename(raw_path))[0] if raw_path != "not linked" else "(not linked)"

        mansavename = simpledialog.askstring("Project Name", "Enter a name for this batch/project:")
        if not mansavename:
            messagebox.showwarning("Missing Name", "You must enter a project/batch name.")
            return

        raw_base = os.path.splitext(os.path.basename(raw_path))[0] if raw_path != "not linked" else "no_raw"
        savename = f"{mansavename}_{extractdate}_{raw_base}"

        output_json_path = os.path.join(self.output_dir, savename + ".json")
        #if not self.skip_conversion:
        jsonfile = mspext.savemetadata(self.metadata, output_json_path)
        #else:
        #    print(f"[debug] not using mspext so no error happens if one tries to save metadata. Will metadata be saved?")
        logger.log(f"Metadata saved to {jsonfile}")
        
        # Store last metadata
        self.last_metadata = self.metadata.copy()


        #temp_ms2 = f"ms2tmp_{self.rawfilename}.csv"

        raw_stem = os.path.splitext(os.path.basename(raw_path))[0] if raw_path != "not linked" else "no_raw"
        temp_ms2 = f"ms2tmp_{raw_stem}.csv"
        temp_ms3 = f"ms3tmp_{raw_stem}.csv"
        final_ms2 = os.path.join(self.output_dir, f"ms2_{savename}.csv")
        final_ms3 = os.path.join(self.output_dir, f"ms3_{savename}.csv")
        success = True
        if self.skip_conversion:
            print(f"[debug]: probably filling missing metadata. No need to do csv rename and check")
            #pass
        else:
            try:
                if os.path.exists(temp_ms2):
                    shutil.move(temp_ms2, final_ms2)
                    logger.log(f"Renamed {temp_ms2} → {final_ms2} and move to {self.output_dir}")
                else:
                    logger.log(f"[WARNING] Missing temp MS2 file: {temp_ms2}")
                    success = False

                if os.path.exists(temp_ms3):
                    shutil.move(temp_ms3, final_ms3)
                    logger.log(f"Renamed {temp_ms3} → {final_ms3} and move to {self.output_dir}")
                else:
                    logger.log(f"[WARNING] Missing temp MS3 file: {temp_ms3}")
                    success = False
            except Exception as e:
                logger.log(f"[ERROR] Failed during file renaming: {e}")
                success = False

            except Exception as e:
                logger.log(f"[ERROR] Failed during file renaming: {e}")
                success = False

            if success:
                mspext.finalize_extraction(self.current_raw_file, self.metadata, savename)
                logname = logger.save(os.path.join(self.output_dir, savename + ".log"), include_debug=True)
                logger.log(f"Saved log to: {logname}")

            self.current_index += 1
            if self.current_index < len(self.raw_file_list):
                self.load_file(self.raw_file_list[self.current_index])
            else:
                messagebox.showinfo("Metadata", "All metadata have been completed.")
                self.window.destroy()
                if self.on_finish:
                    self.on_finish()
        # Trigger callback after metadata is created
        if self.callback:
            self.callback(
                raw_path,             # resolved raw path
                self.metadata,        # metadata dict
                savename              # base name for saved file
            )

        # Auto-close the window
        self.window.destroy()

def write_to_gui(message):
    text_widget.insert(tk.END, message + "\n")
    text_widget.see(tk.END)

logger.gui_writer = write_to_gui

#version info reader
#will add same one in extractor but may depreciate after moving metadata editing part to GUI
def get_version_info_for(module_name, manifest_path="project_version_manifest.ini"):
    config = configparser.ConfigParser()
    config.read(manifest_path)

    if module_name in config:
        version = config[module_name].get("version", "N/A")
        last_update = config[module_name].get("last_update", "N/A")
        return version, last_update
    return None, None


# Function to create and hide Tkinter root window
# If future version supports GUI fully available, please reconstruct this part (not hiding main window, but think about what should be there)
def create_tkinter_root():
    root = tk.Tk()
    root.withdraw() # hide main window
    return root

#about tab showing the version info
def open_about_window():
    about_win = tk.Toplevel(root)
    about_win.title("About GlycoMSP")
    about_win.geometry("500x400")

    tk.Label(about_win, text="About GlycoMSP GUI", font=("Arial", 14, "bold")).pack(pady=10)
    tk.Label(about_win, text="Made by Huan-Chuan TSENG", font=("Arial", 10)).pack(pady=5)
    
    # Create a scrollable text box
    about_text = tk.Text(about_win, wrap=tk.WORD, height=20, width=60)
    about_text.pack(padx=10, pady=5, fill="both", expand=True)

    scrollbar = tk.Scrollbar(about_text)
    scrollbar.pack(side="right", fill="y")
    about_text.config(yscrollcommand=scrollbar.set)
    scrollbar.config(command=about_text.yview)

    # Load version info from manifest, mind the path issues (on demo folder now)
    config = configparser.ConfigParser()
    config.read(".\\2025_demo\\project_version_manifest.ini")
    #config.read("project_version_manifest.ini")

    if "version_check" in config:
        ts = config["version_check"].get("checked_on", "Unknown")
        about_text.insert(tk.END, f"Version manifest generated on: {ts}\n\n")

    for section in config.sections():
        if section != "version_check":
            version = config[section].get("version", "N/A")
            updated = config[section].get("last_update", "N/A")
            about_text.insert(tk.END, f"[{section}]\n  Version: {version}\n  Last Update: {updated}\n\n")

    about_text.config(state=tk.DISABLED)

    tk.Label(about_win, text="Component version listed above", font=("Arial", 10)).pack(pady=5)

    tk.Button(about_win, text="Close", command=about_win.destroy).pack(pady=10)

# Function to validate file path
def validate_file_path(file_path):
    if not os.path.exists(file_path):
        print(f"Error: The file {file_path} does not exist.")
        return False
    if not os.path.isfile(file_path):
        print(f"Error: The path {file_path} is not a file.")
        return False
    return True

def update_display():
    text_widget.delete(1.0, tk.END)
    for ftype, path in selected_files.items():
        text_widget.insert(tk.END, f"{ftype}: {path}\n")

def clear_files():
    # Placeholder for your implementation
    selected_files.clear()
    update_display()

def reset_main_status():
    status_var.set("Idle")
    progress.stop()
    convert_button.config(state="normal")


def on_metadata_ready(raw_file, metadata, savename):
    logger.log(f"Confirmed metadata for {raw_file}")
    # Pass to peak extractor
    mspext.convert_raw_to_csv(raw_file, debug=False)

def launch_metadata_for_all(files):
    MetadataEditorWindow(root, files, on_metadata_ready, on_finish=reset_main_status)





# Function to select raw file for further pre-processing
def select_file(filetype):
    filetypes_dict = {
        "raw": [("Raw file", "*.raw")],
        "mzml": [("mzML file", "*.mzML")],
        "csv": [("CSV file", "*.csv")],
        "excel": [("Excel file", "*.xls *.xlsx")]
    }

    #filepath = filedialog.askopenfilename(filetypes=filetypes_dict.get(filetype, [("All files", "*.*")]))
    filepaths = filedialog.askopenfilenames(filetypes=filetypes_dict.get(filetype, [("All files", "*.*")]))
    if filepaths:
        valid_paths = [fp for fp in filepaths if validate_file_path(fp)]
        selected_files.setdefault(filetype, []).extend(valid_paths)
        update_display()
        #validation of values
        #file_contents_validation
    elif not filepaths:
        print(f"No valid {filetype} file selected.")


# generic select file (20250407)
def select_files_generic(filetype_key, allow_multiple=False, on_select_callback=None):
    filetypes_dict = {
        "csv": [("CSV files", "*.csv")],
        "excel": [("Excel files", "*.xls *.xlsx")],
        "json": [("JSON files", "*.json")],
        "all": [("All files", "*.*")]
    }

    filetypes = filetypes_dict.get(filetype_key, filetypes_dict["all"])

    if allow_multiple:
        filepaths = filedialog.askopenfilenames(filetypes=filetypes)
    else:
        filepath = filedialog.askopenfilename(filetypes=filetypes)
        filepaths = [filepath] if filepath else []

    if filepaths and on_select_callback:
        on_select_callback(filepaths)


def launch_metadata_batch():
    if not mspext.pymsreader:
        logger.log("[ERROR] Raw file conversion is not available — pymsfilereader missing.")
        messagebox.showerror(
            "Raw Conversion Not Available",
            "Thermo MSFileReader or pymsfilereader is not installed.\n"
            "Please run on a compatible Windows system with the required libraries."
        )
        return

    if "raw" not in selected_files or not selected_files["raw"]:
        messagebox.showwarning("No Raw File", "Please select at least one raw file.")
        return

    rawfilelist = selected_files["raw"]
    status_var.set("Converting raw file to csv...check the metadata assignment window.")
    progress.start()
    MetadataEditorWindow(root, rawfilelist, on_metadata_ready, on_finish=reset_main_status)





def on_batch_conversion_complete(converted_raws):
    progress.stop()
    convert_button.config(state="normal")
    status_var.set("Raw file conversions complete. Proceeding to metadata entry...")
    messagebox.showinfo("Conversion Complete", "Raw file conversion finished.\nNow entering metadata for each file.")
    
    launch_metadata_for_all(converted_raws)


def open_prepare_dataset_window():
    subwin = tk.Toplevel(root)
    subwin.title("Prepare Dataset")
    subwin.geometry("800x550")

    # --- Nested project data ---
    experiment_projects = {}  # {experiment: {"samples": {sample_name: {csv, excel, json}}}}
    linked_validated_samples = set() #validated samples
    validation_failed_samples = set() #failed sample
    #added 20250412
    experiment_method_paths = {}  # Store .exp.json path per experiment
    sample_method_folder = None  # Global path for saving per-sample method.json files
    experiment_status_labels = {}  # GUI labels for status display, indexed by experiment
    # --- Treeview UI ---
    tree = ttk.Treeview(subwin)
    tree.heading("#0", text="Dataset Explorer", anchor="w")
    tree.pack(expand=True, fill="both", padx=10, pady=10)

    # --- Tree logic ---
    def on_tree_select(event):
        sel = tree.selection()
        if not sel:
            link_button.config(state="disabled")
            return

        item_text = tree.item(sel[0], "text")
        sample_name = clean_sample_name(item_text)

        parent_id = tree.parent(sel[0])
        exp_text = tree.item(parent_id, "text") if parent_id else ""
        exp_name = exp_text.replace("Experiment: ", "").split(" (")[0].strip()

        # Check if it's a valid sample in the experiment and if  it's okay to merge
        if exp_name in experiment_projects and sample_name in experiment_projects[exp_name]["samples"]:
            if (exp_name, sample_name) not in linked_validated_samples:
                link_button.config(state="normal")
                merge_button.config(state="disabled")
            else:
                link_button.config(state="disabled")
                merge_button.config(state="normal")

    tree.bind("<<TreeviewSelect>>", on_tree_select)

    def try_link_selected_sample():
        sel = tree.selection()
        if not sel:
            return
        sample_id = sel[0]
        exp_id = tree.parent(sample_id)
        sample_name = clean_sample_name(tree.item(sel[0], "text"))
        #sample_name = tree.item(sample_id, "text").replace("Sample: ", "").split(" (")[0]
        exp_name = tree.item(exp_id, "text").replace("Experiment: ", "")
        link_and_validate_sample(exp_name, sample_name)

    # --- Tree refresh logic ---
    def refresh_tree():
        tree.delete(*tree.get_children())

        for exp_title, exp_data in experiment_projects.items():
            exp_node = tree.insert("", "end", text=f"Experiment: {exp_title}", open=True)

            for sample_name, files in exp_data.get("samples", {}).items():
                has_csv = bool(files.get("csv"))
                has_excel = bool(files.get("excel"))
                has_json = bool(files.get("json"))

                # Decide how to display the sample label
                if (exp_title, sample_name) in validation_failed_samples:
                    sample_display = f"⛔ Sample: {sample_display} (validation failed)"
                if (exp_title, sample_name) in linked_validated_samples:
                    sample_display = f"✅ Sample: {sample_name}"
                elif has_csv and has_excel and has_json:
                    sample_display = f"⚠️ Sample: {sample_name} (unvalidated)"
                elif not has_json:
                    sample_display = f"❌ Sample: {sample_name} (metadata missing)"
                else:
                    sample_display = f"Sample: {sample_name}"

                sample_node = tree.insert(exp_node, "end", text=sample_display, open=True)

                for ftype in ["csv", "excel", "json"]:
                    if files.get(ftype):
                        label = ftype.upper() if ftype != "json" else "Metadata"
                        tree.insert(sample_node, "end", text=f"{label}: {os.path.basename(files[ftype])}")
    # --- statistics ---
    def show_experiment_summary(exp_name):
        if exp_name not in experiment_projects:
            return

        total = 0
        validated = 0

        for sname in experiment_projects[exp_name]["samples"]:
            total += 1
            if (exp_name, sname) in linked_validated_samples:
                validated += 1

        unvalidated = total - validated

        summary = (
            f"🧪 Experiment Summary: {exp_name}\n\n"
            f"• Total Samples: {total}\n"
            f"• ✅ Validated: {validated}\n"
            f"• ⚠️ Unvalidated: {unvalidated}"
        )
        messagebox.showinfo("Experiment Status", summary)
    # --- prevent symbols getting read --- v9 changed 
    def clean_sample_name(text):
        for symbol in ["✅", "⚠️", "❌", "⛔"]:
            if text.startswith(symbol):
                text = text[len(symbol):].strip()

        if text.startswith("Sample:"):
            text = text[len("Sample:"):].strip()

        # Remove validation tag suffix
        if " (metadata missing)" in text:
            text = text.replace(" (metadata missing)", "")
        if " (validation failed)" in text:
            text = text.replace(" (validation failed)", "")
        if " (unvalidated)" in text:
            text = text.replace(" (unvalidated)", "")

        return text.strip()
    # --- Add sample to selected experiment ---
    def add_sample():
        sel = tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select an experiment to add a sample.")
            return

        item_id = sel[0]
        item_text = tree.item(item_id, "text")
        if not item_text.startswith("Experiment:"):
            messagebox.showinfo("Invalid Selection", "You must select an experiment, not a sample or file.")
            return

        exp_title = item_text.replace("Experiment: ", "")
        sample_name = simpledialog.askstring("New Sample", "Enter sample name:")
        if not sample_name:
            return
        if sample_name in experiment_projects[exp_title]["samples"]:
            messagebox.showwarning("Duplicate Sample", f"Sample '{sample_name}' already exists.")
            return

        experiment_projects[exp_title]["samples"][sample_name] = {"csv": None, "excel": None, "json": None}
        refresh_tree()

    def clean_unassigned_samples():
        exp = "Unassigned"
        if exp not in experiment_projects:
            return
        samples = experiment_projects[exp]["samples"]
        to_delete = []

        for s_name, files in samples.items():
            if not files.get("csv") and not files.get("excel") and not files.get("json"):
                to_delete.append(s_name)

        for s in to_delete:
            del samples[s]

        if not samples:
            del experiment_projects[exp]  # if no samples left, remove experiment too
        refresh_tree()

    def store_metadata_back_to_sample(exp_name, sample_name, rawfile, metadata, savename):
        json_path = os.path.join(os.path.dirname(rawfile), savename + ".json")

        # Fallback if sample doesn't exist yet (e.g., metadata was created before file assignment)
        if sample_name not in experiment_projects[exp_name]["samples"]:
            experiment_projects[exp_name]["samples"][sample_name] = {"csv": None, "excel": None, "json": None}

        sample = experiment_projects[exp_name]["samples"][sample_name]
        sample["json"] = json_path

        # Try to rename sample to raw name (from metadata), only if different
        try:
            raw_base = os.path.splitext(os.path.basename(metadata.get("Raw filename", "")))[0]
        except Exception:
            raw_base = sample_name

        if raw_base != sample_name and raw_base not in experiment_projects[exp_name]["samples"]:
            experiment_projects[exp_name]["samples"][raw_base] = sample
            del experiment_projects[exp_name]["samples"][sample_name]
            sample_name = raw_base

        refresh_tree()
        write_method_file(exp_name)
    
    def open_metadata_editor_for_sample(exp_name, sample_name):
        from tkinter import simpledialog

        sample = experiment_projects[exp_name]["samples"][sample_name]

        # Let user select raw file if needed (optional)
        raw_file_path = filedialog.askopenfilename(
            title="Select RAW file for this sample (optional)",
            filetypes=[("Raw files", "*.raw")]
        )
        raw_file_path = raw_file_path or "(not linked)"
        print(f"[debug] set to not linked and need to check if we can add one back in metadata editor window")

        # Ask user for output folder
        outdir = filedialog.askdirectory(title="Select output folder to save metadata")
        if not outdir:
            messagebox.showwarning("Cancelled", "Metadata creation cancelled.")
            return

        editor = MetadataEditorWindow(
            parent=root,
            raw_file_list=[raw_file_path],
            on_each_metadata_ready_callback=lambda rf, md, name: store_metadata_back_to_sample(exp_name, sample_name, rf, md, name),
            output_dir=outdir,
            skip_conversion=True
        )

        # Delay execution of load_file until window is fully initialized
        if not editor.skip_conversion:
            editor.window.after(10, lambda: editor.load_file(raw_file_path))
        else:
            print(f"[debug] skipping loading raw file")


    # --- link and validate the grouped sample ---
    def link_and_validate_sample(exp_name, sample_name):

        sample_name = clean_sample_name(sample_name)
        if sample_name not in experiment_projects[exp_name]["samples"]:
            logger.log(f"[ERROR] Cleaned sample name '{sample_name}' not found under '{exp_name}'")
            messagebox.showerror("Invalid Sample", f"Sample not found in experiment: {sample_name}")
            return
        sample = experiment_projects[exp_name]["samples"][sample_name]

        # Check file presence
        if not sample.get("csv") or not sample.get("excel"):
            messagebox.showerror("Missing File", "Sample must have both a CSV and Excel file before linking.")
            return
        

        # Validate CSV format
        #try:
        if not mspval.validate_csv_structure(sample["csv"]): #when it returns False
        #except Exception as e:
            
            messagebox.showerror("CSV Validation Failed")#, f"{sample['csv']}")
            logger.log(f"[Validation] CSV failed: {sample['csv']}\n")
            #raise ValueError("csv read error: file may be corrupted or unreadable.")
            validation_failed_samples.add((exp_name, sample_name))
            refresh_tree()
            return  # ← this prevents continuing to metadata
        else:
            logger.log(f"[Validation] CSV successful: {sample['csv']}\n")
        #Validate Excel format (assume sheet 'MSlist' exists and ion list is valid) 

        try:
            xl = pd.read_excel(sample["excel"], sheet_name=None)
            if not isinstance(xl, dict):
                raise ValueError("Excel read error: file may be corrupted or unreadable.")
            if "MSlist" not in xl:
                raise ValueError("Sheet 'MSlist' not found in Excel.")
        except Exception as e:
            messagebox.showerror("Excel Validation Failed", f"{sample['excel']}\n\n{e}")
            logger.log(f"[Validation] Excel failed: {sample['excel']}\n{e}")
            validation_failed_samples.add((exp_name, sample_name))
            refresh_tree()
            return  # ← this prevents continuing to metadata

        #try:
        anno = pd.ExcelFile(sample["excel"])
        MSlistdf = pd.read_excel(anno, sheet_name="MSlist")
        ionlistdf = pd.read_excel(anno, sheet_name="ionlist")
        if not mspval.validate_annotation_structure(MSlistdf):
        #except Exception as e:
            messagebox.showerror("Excel Validation Failed")#), f"{sample['excel']}\n")
            logger.log(f"[Validation] Excel failed: {sample['excel']}\n")
            #raise ValueError("Excel read error: file may be corrupted or unreadable.")
            validation_failed_samples.add((exp_name, sample_name))
            refresh_tree()
            return  # ← this prevents continuing to metadata
        else:
            logger.log(f"[Validation] Excel successful: {sample['csv']}\n")
        ion_df = ionlistdf[["mass"]] 
        if ion_df["mass"].dtype != "float64":
            logger.log(f"[Validation] ion list has invalid values: {sample['csv']}\n")
            return

        # Check metadata
        if not sample.get("json"):
            proceed = messagebox.askyesno("Metadata Missing", "No metadata found. Would you like to create it now?")
            if not proceed:
                return
            open_metadata_editor_for_sample(exp_name, sample_name)
            return

        # All checks passed → mark as validated
        linked_validated_samples.add((exp_name, sample_name))
        write_method_file(exp_name)
        messagebox.showinfo("Validated", f"Sample '{sample_name}' under '{exp_name}' is now validated.")
        refresh_tree()

    def try_merge_selected_sample():
        sel = tree.selection()
        if not sel:
            return

        sample_node = sel[0]
        sample_name = clean_sample_name(tree.item(sample_node, "text"))
        exp_node = tree.parent(sample_node)
        exp_name = tree.item(exp_node, "text").replace("Experiment: ", "").split(" (")[0].strip()

        files = experiment_projects[exp_name]["samples"][sample_name]
        if not all([files.get("csv"), files.get("excel"), files.get("json")]):
            messagebox.showerror("Error", "Sample is missing required files.")
            return

        # Ask user for output folder
        outdir = filedialog.askdirectory(title="Select output folder to save merged dataset")
        if not outdir:
            return

        # Ask user for ion sheet name (later we will allow external ion file)
        #ion_sheet_name = simpledialog.askstring("Ion Sheet", "Enter ion sheet name (in Excel):", initialvalue="core_OG")
        metadata_path=files["json"]
        print(f"debug: metadata file: {metadata_path}")
        try:
            with open(metadata_path, "r") as f:
                metadata = json.load(f)
            derivatization_type = metadata.get("Derivatization Type", "Others")
            if derivatization_type == "Others":
                value = simpledialog.askfloat("Custom Derivatization Mass","Enter custom derivatization mass (e.g., 50.1234):")
                if value is None:
                    return  # user cancelled
                else:
                    derivatization_type = value
                    # pass it into directassign_files or your protonated mass logic
                    # need to save the info back later
                
            
            today = datetime.now().strftime("%Y%m%d")
            outname = f"{sample_name}_merged_{today}.csv"
            outpath = os.path.join(outdir, outname)
            
            pre_df, iondfindex, ion_df = mspval.directassign_files(files["excel"], files["csv"],derivatization_type, debug = False)
            mspval.createnormailzedionlistcsv(iondfindex, pre_df,ion_df, outpath)
            messagebox.showinfo("Merge Complete", f"Dataset saved:\n{os.path.basename(os.path.basename(outpath))}")
        except Exception as e:
            messagebox.showerror("Merge Failed", f"Error:\n{str(e)}")

    # --- Assign file to experiment/sample ---
    def assign_file(filetype, filepath, exp_title="Unassigned", sample_name="Unassigned"):
        if exp_title not in experiment_projects:
            experiment_projects[exp_title] = {"samples": {}}
        if sample_name not in experiment_projects[exp_title]["samples"]:
            experiment_projects[exp_title]["samples"][sample_name] = {"csv": None, "excel": None, "json": None}
        experiment_projects[exp_title]["samples"][sample_name][filetype] = filepath
        refresh_tree()

    def extract_rawname_from_metadata(json_path):
        try:
            with open(json_path, "r") as f:
                meta = json.load(f)
            rawbase = os.path.basename(meta.get("Raw filename", ""))
            return os.path.splitext(rawbase)[0]
        except:
            return None

    # --- File handlers ---
    def handle_csv_selection(filepaths):
        for path in filepaths:
            sample_id = os.path.splitext(os.path.basename(path))[0]
            assign_file("csv", path, "Unassigned", sample_id)

    def handle_excel_selection(filepaths):
        for path in filepaths:
            sample_id = os.path.splitext(os.path.basename(path))[0]
            assign_file("excel", path, "Unassigned", sample_id)

    def handle_json_selection(filepaths):
        for path in filepaths:
            title = extract_title_from_metadata(path) or "Unassigned"
            rawname = extract_rawname_from_metadata(path) or "Unassigned"

            assign_file("json", path, title, "Unassigned")

            sample_name = "Unassigned"

            # Auto-rename sample (move from 'Unassigned' to rawname)
            if title in experiment_projects:
                samples = experiment_projects[title]["samples"]
                if "Unassigned" in samples:
                    if rawname in samples:
                        messagebox.showwarning("Sample Exists", f"Sample '{rawname}' already exists. Skipping rename.")
                    else:
                        samples[rawname] = samples.pop("Unassigned")
                        sample_name = rawname  # <- use updated sample name
                        refresh_tree()

            # After rename, auto-write .method.json if csv + excel exist
            sample = experiment_projects[title]["samples"].get(sample_name)
            if sample and sample.get("csv") and sample.get("excel"):
                from datetime import datetime
                entry = {
                    "csv": os.path.basename(sample["csv"]),
                    "excel": os.path.basename(sample["excel"]),
                    "metadata": os.path.basename(path),
                    "raw_file": rawname,
                    "validated": False
                }
                sample_method = {
                    "experiment": title,
                    "samples": {sample_name: entry},
                    "generated_on": datetime.now().strftime("%Y-%m-%d %H:%M")
                }
                try:
                    sample_path = os.path.join(os.path.dirname(path), f"{sample_name}.method.json")
                    with open(sample_path, "w") as f:
                        json.dump(sample_method, f, indent=4)
                    logger.log(f"[Method] Auto-saved method file for sample: {sample_path}")
                    messagebox.showinfo("Method File Saved", f"A sample method file was saved:\n{os.path.basename(sample_path)}")

                except Exception as e:
                    logger.log(f"[WARNING] Failed to save method file for '{sample_name}': {e}")
    def extract_title_from_metadata(json_path):
        try:
            with open(json_path, "r") as f:
                meta = json.load(f)
            return meta.get("Experiment Title")
        except:
            return None

    def select_files_generic(filetype_key, allow_multiple=False, on_select_callback=None):
        filetypes_dict = {
            "csv": [("CSV files", "*.csv")],
            "excel": [("Excel files", "*.xls *.xlsx")],
            "json": [("JSON files", "*.json")],
            "all": [("All files", "*.*")]
        }
        filetypes = filetypes_dict.get(filetype_key, filetypes_dict["all"])
        if allow_multiple:
            filepaths = filedialog.askopenfilenames(filetypes=filetypes)
        else:
            filepath = filedialog.askopenfilename(filetypes=filetypes)
            filepaths = [filepath] if filepath else []
        if filepaths and on_select_callback:
            on_select_callback(filepaths)

    def update_status_display(exp_name):
        # Update experiment-level path display
        path_label, sample_label = experiment_status_labels.get(exp_name, (None, None))
        if path_label:
            path = experiment_method_paths.get(exp_name)
            path_label.set(f"Experiment method file: {path if path else 'None'}")

        if sample_label:
            sample_label.set(f"Sample method will be saved at: {sample_method_folder if sample_method_folder else 'None'}")


    # --- Right-click move logic ---
    def on_right_click(event):
        item_id = tree.identify_row(event.y)
        if not item_id:
            return

        selected_text = tree.item(item_id, "text")
        sample_id = tree.parent(item_id)
        exp_id = tree.parent(sample_id)
        sample_name = clean_sample_name(tree.item(sample_id, "text"))#.replace("Sample: ", "").split(" (")[0]
        exp_name = tree.item(exp_id, "text").replace("Experiment: ", "")

        if ":" not in selected_text:
            return  # not a file entry
        
        # Right-click menu option: Remove file
        def remove_file(exp_name, sample_name, filetype):
            entry = experiment_projects[exp_name]["samples"][sample_name]
            entry[filetype] = None

            # Clear validation state if needed
            linked_validated_samples.discard((exp_name, sample_name))
            validation_failed_samples.discard((exp_name, sample_name))

            refresh_tree()

        filetype = selected_text.split(":")[0].strip().lower()
        filename = selected_text.split(":")[1].strip()

        menu = tk.Menu(subwin, tearoff=0)
        move_menu = tk.Menu(menu, tearoff=0)

        for e_name, e_data in experiment_projects.items():
            for s_name in e_data["samples"]:
                move_menu.add_command(
                    label=f"{e_name} → {s_name}",
                    command=lambda en=e_name, sn=s_name: move_file(from_exp=exp_name, from_sample=clean_sample_name(sample_name), ftype=filetype, filename=filename, to_exp=en, to_sample=clean_sample_name(sn))
                )
        menu.add_cascade(label="Move to...", menu=move_menu)
        
        # Right-click on sample: enable linking
        if selected_text.startswith("Sample: "):
            parent_id = tree.parent(item_id)
            sample_name = selected_text.replace("Sample: ", "").split(" (")[0]
            exp_name = tree.item(parent_id, "text").replace("Experiment: ", "")
            menu.add_command(
                label="Link and Validate Sample",
                command=lambda: link_and_validate_sample(exp_name, sample_name)
            )
        filetype = selected_text.split(":")[0].strip().lower()
        filename = selected_text.split(":")[1].strip()
        # Add remove option if it's a valid file
        menu.add_command(label=f"Remove {filetype.upper()}",command=lambda: remove_file(exp_name, clean_sample_name(sample_name), filetype))  
        menu.post(event.x_root, event.y_root)
          



    def move_file(from_exp, from_sample, ftype, filename, to_exp, to_sample):
        
        entry = experiment_projects[from_exp]["samples"][from_sample][ftype]
        if entry and os.path.basename(entry) == filename:
            experiment_projects[from_exp]["samples"][from_sample][ftype] = None
            assign_file(ftype, entry, to_exp, to_sample)

    # --- write to method ---
    def write_method_file(exp_name, auto=False):
        method = {
            "experiment": exp_name,
            "generated_on": datetime.now().strftime("%Y-%m-%d %H:%M"),
            "samples": {}
        }

        for sname, files in experiment_projects[exp_name]["samples"].items():
            if not all([files.get("csv"), files.get("excel"), files.get("json")]):
                continue
            try:
                with open(files["json"], "r") as f:
                    meta = json.load(f)
            except Exception:
                meta = {}

            entry = {
                "csv": os.path.abspath(os.path.normpath(files["csv"])),#files["csv"],  #os.path.basename(files["csv"]),
                "excel": os.path.abspath(os.path.normpath(files["excel"])),#files["excel"], #os.path.basename(files["excel"]),
                "metadata": os.path.abspath(os.path.normpath(files["json"])),#files["json"], #os.path.basename(files["json"]),
                "raw_file": os.path.abspath(os.path.normpath(meta.get("Raw filename", "not linked"))),#meta.get("Raw filename", "not linked"),
                "validated": (exp_name, sname) in linked_validated_samples
            }
            if "ion_sheet" in files:
                entry["ion_sheet"] = files["ion_sheet"]
            if "ion_sheet_file" in files:
                entry["ion_sheet_file"] = files["ion_sheet_file"]

            method["samples"][sname] = entry

            # Save per-sample method.json
            try:
                sample_method = {
                    "experiment": exp_name,
                    "samples": {sname: entry},
                    "generated_on": method["generated_on"]
                }
                if sample_method_folder:
                    sample_path = os.path.join(sample_method_folder, f"{sname}.method.json")
                else:
                    sample_path = os.path.join(os.path.dirname(files["json"]), f"{sname}.method.json")

                with open(sample_path, "w") as sf:
                    json.dump(sample_method, sf, indent=4)
                logger.log(f"[Method] Per-sample method saved: {sample_path}")
            except Exception as e:
                logger.log(f"[WARNING] Failed to save per-sample method for {sname}: {e}")

        # Save experiment-level .exp.json
        if auto and exp_name in experiment_method_paths:
            output_path = experiment_method_paths[exp_name]
        else:
            output_path = filedialog.asksaveasfilename(
                defaultextension=".exp.json",
                initialfile=f"{exp_name}.exp.json",
                filetypes=[("Experiment Method JSON", "*.exp.json")]
            )
            if not output_path:
                return
            experiment_method_paths[exp_name] = output_path

        try:
            with open(output_path, "w") as f:
                json.dump(method, f, indent=4)
            logger.log(f"[Method] Experiment method saved to: {output_path}")
        except Exception as e:
            messagebox.showerror("Save Failed", f"Could not save method file:\n{e}")

        update_status_display(exp_name)

    # --- load method file ---

    # --- safe loader ---
    def safe_get_field(d, key, fallback="(not linked)"):
        val = d.get(key)
        return val if isinstance(val, str) and val.strip() else fallback
    #def check_relative_location():
    def load_method_file(paths=None):
        if not paths:
            paths = filedialog.askopenfilenames(
                title="Select One or More Sample Method Files",
                filetypes=[("Sample Method JSON", "*.method.json"), ("JSON files", "*.json")]
            )
        if not paths:
            return

        loaded = 0
        skipped = 0
        for path in paths:
            try:
                with open(path, "r") as f:
                    method = json.load(f)
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load method file:\n{path}\n{e}")
                continue

            exp_name = method.get("experiment", "Recovered")
            samples = method.get("samples", {})
            if not samples:
                skipped += 1
                continue

            if exp_name not in experiment_projects:
                experiment_projects[exp_name] = {"samples": {}}

            for sample_name, files in samples.items():
                sample_name = clean_sample_name(sample_name)  # Normalize sample name

                if sample_name in experiment_projects[exp_name]["samples"]:
                    skipped += 1
                    continue

                base = os.path.dirname(path)
                files_resolved = {
                    "csv": os.path.join(base, files.get("csv")),
                    "excel": os.path.join(base, files.get("excel")),
                    "json": os.path.join(base, files.get("metadata"))
                }
                experiment_projects[exp_name]["samples"][sample_name] = files_resolved

                if files.get("validated"):
                    linked_validated_samples.add((exp_name, sample_name))
                    
                loaded += 1

        refresh_tree()
        msg = f"Imported {loaded} sample(s) successfully.\nSkipped: {skipped}"
        if skipped > 0:
            msg += "\n(Skipped files were already loaded or invalid)"
        messagebox.showinfo("Method Import", msg)

    def change_experiment_method_path():
        path = filedialog.asksaveasfilename(
            title="Select path to save experiment method file",
            defaultextension=".exp.json",
            filetypes=[("Experiment Method JSON", "*.exp.json")]
        )
        if path:
            experiment_method_paths["Unassigned"] = path  # Replace key if you're in a real experiment context
            update_status_display("Unassigned")

    def change_sample_method_folder():
        global sample_method_folder
        path = filedialog.askdirectory(title="Select folder to save sample method files")
        if path:
            sample_method_folder = path
            update_status_display("Unassigned")

    def load_experiment_method_file():
        path = filedialog.askopenfilename(
            title="Load experiment method file",
            filetypes=[("Experiment Method JSON", "*.exp.json")]
        )
        if not path:
            return

        try:
            with open(path, "r") as f:
                method = json.load(f)
        except Exception as e:
            messagebox.showerror("Load Failed", f"Could not load file:\n{e}")
            return

        exp_name = method.get("experiment", "Recovered")
        if exp_name not in experiment_projects:
            experiment_projects[exp_name] = {"samples": {}}

        experiment_method_paths[exp_name] = path
        for sample_name, files in method.get("samples", {}).items():
            sample = {
                "csv": os.path.join(os.path.dirname(path), files["csv"]),
                "excel": os.path.join(os.path.dirname(path), files["excel"]),
                "json": os.path.join(os.path.dirname(path), files["metadata"])
            }
            if "ion_sheet" in files:
                sample["ion_sheet"] = files["ion_sheet"]
            if "ion_sheet_file" in files:
                sample["ion_sheet_file"] = files["ion_sheet_file"]
            experiment_projects[exp_name]["samples"][sample_name] = sample

            if files.get("validated") is True:
                linked_validated_samples.add((exp_name, sample_name))
            elif files.get("validated") is False:
                validation_failed_samples.add((exp_name, sample_name))

        experiment_status_labels[exp_name] = (exp_status_var, sample_status_var)
        refresh_tree()
        update_status_display(exp_name)
        messagebox.showinfo("Loaded", f"Method file loaded for experiment: {exp_name}")


    def save_current_experiment_method():
        sel = tree.selection()
        if not sel:
            messagebox.showwarning("No Experiment Selected", "Please select an experiment to save.")
            return

        sample_id = sel[0]
        exp_id = tree.parent(sample_id)
        if not exp_id:  # If we're already on an experiment node
            exp_id = sample_id

        exp_name = tree.item(exp_id, "text").replace("Experiment: ", "").split(" (")[0].strip()
        write_method_file(exp_name, auto=True)

    # --- Button panel ---
    button_frame = tk.Frame(subwin)
    button_frame.pack(pady=5)

    tk.Button(button_frame, text="Select MS2 CSV(s)", command=lambda: select_files_generic("csv", True, handle_csv_selection)).grid(row=0, column=0, padx=5)
    tk.Button(button_frame, text="Select Excel", command=lambda: select_files_generic("excel", True, handle_excel_selection)).grid(row=0, column=1, padx=5)
    tk.Button(button_frame, text="Add sample (metadata file required)", command=lambda: select_files_generic("json", True, handle_json_selection)).grid(row=1, column=0, padx=5)
    tk.Button(button_frame, text="Add Sample (missing metadata)", command=add_sample).grid(row=1, column=1, padx=5)
    tk.Button(button_frame, text="Clean up empty unassigned sample tags", command=clean_unassigned_samples).grid(row=1, column=2, padx=5)
    link_button = tk.Button(button_frame, text="Link Sample", state="disabled", command=lambda: try_link_selected_sample())
    link_button.grid(row=2, column=0, padx=5)
    merge_button = tk.Button(button_frame, text="Merge Sample", state="disabled", command=lambda: try_merge_selected_sample())
    merge_button.grid(row=2, column=1, padx=5)
    tk.Button(button_frame, text="Load Method", command=load_method_file).grid(row=2, column=2, padx=5)
    tk.Button(subwin, text="Close", command=subwin.destroy).pack(pady=10)


    #status label
    status_frame = tk.Frame(subwin, relief=tk.SUNKEN, borderwidth=1)
    status_frame.pack(fill=tk.X, padx=10, pady=(0, 5))
    exp_status_var = tk.StringVar(value="Experiment method file: None")
    sample_status_var = tk.StringVar(value="Sample method will be saved at: None")
    tk.Label(status_frame, textvariable=exp_status_var, anchor="w").pack(fill=tk.X, pady=2)
    tk.Label(status_frame, textvariable=sample_status_var, anchor="w").pack(fill=tk.X, pady=2)
    # Store label references for this experiment
    experiment_status_labels["Unassigned"] = (exp_status_var, sample_status_var)

    btn_frame = tk.Frame(subwin)
    btn_frame.pack(fill=tk.X)

    tk.Button(btn_frame, text="Set Experiment Method File", command=change_experiment_method_path).pack(side=tk.LEFT, padx=10, pady=5)
    tk.Button(btn_frame, text="Set Sample Method Folder", command=change_sample_method_folder).pack(side=tk.LEFT, padx=10, pady=5)
    #bottom place for "global" exp method file
    tk.Button(btn_frame, text="Load .exp.json", command=load_experiment_method_file).pack(side=tk.LEFT, padx=10)
    tk.Button(btn_frame, text="Save .exp.json", command=save_current_experiment_method).pack(side=tk.LEFT)

    # --- Right-click bind ---
    tree.bind("<Button-3>", on_right_click)


def open_ml_analysis_window():
    import json
    import os
    import pandas as pd
    from tkinter import filedialog, messagebox
    from tkinter import ttk
    import tkinter as tk

    def select_train_csv():
        nonlocal train_csv_path, linked_exp_json
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if not path:
            return
        train_csv_path = os.path.abspath(path)

        exp_json_path = None
        base_dir = os.path.dirname(path)
        for fname in os.listdir(base_dir):
            if fname.endswith(".exp.json"):
                exp_json_path = os.path.join(base_dir, fname)
                break

        linked_exp_json = exp_json_path

        if exp_json_path and os.path.exists(exp_json_path):
            try:
                with open(exp_json_path, "r") as f:
                    exp_data = json.load(f)
                exp_data["train_csv"] = train_csv_path
                with open(exp_json_path, "w") as f:
                    json.dump(exp_data, f, indent=4)
            except Exception as e:
                print(f"[ML] Failed to update experiment JSON: {e}")

        origin_info.configure(state="normal")
        origin_info.delete(1.0, "end")
        origin_info.insert("end", f"Loaded file: {os.path.basename(path)}\n")
        if exp_json_path:
            origin_info.insert("end", f"Linked .exp.json: {os.path.basename(exp_json_path)}\n")
        origin_info.insert("end", f"Path: {path}")
        origin_info.configure(state="disabled")

    def open_train_settings():
        settings = tk.Toplevel()
        settings.title("Train/Test Parameters")
        settings.geometry("300x200")

        tk.Label(settings, text="Test Split Ratio (e.g., 0.25 = 75/25)").pack(pady=(10, 0))
        split_slider = tk.Scale(settings, from_=0.1, to=0.5, resolution=0.05, orient="horizontal", variable=test_split_var)
        split_slider.pack(pady=5)

        tk.Label(settings, text="Minimum Samples per Class").pack(pady=(15, 0))
        min_sample_spin = tk.Spinbox(settings, from_=1, to=50, textvariable=min_samples_var, width=5)
        min_sample_spin.pack(pady=5)

        tk.Button(settings, text="Close", command=settings.destroy).pack(pady=10)

    def train_model():
        try:
            from sklearn.ensemble import RandomForestClassifier
            from sklearn.model_selection import train_test_split
            from sklearn.metrics import classification_report
            from sklearn.preprocessing import LabelEncoder
            import joblib
            from sklearn.utils.multiclass import unique_labels

        except ImportError as e:
            messagebox.showerror("Missing Dependencies", "Required package not found.\n\nPlease make sure scikit-learn and joblib are installed.")
            return

        if not train_csv_path:
            messagebox.showwarning("No File", "Please select a training CSV file first.")
            return

        label_col = label_dropdown.get().strip()
        if label_col not in ['Structure', 'IUPACname(optional)', 'Glycanannotation2', 'GlyToucan ID']:
            messagebox.showerror("Invalid Label", f"'{label_col}' is not a supported label column.")
            return

        try:
            df = pd.read_csv(train_csv_path)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read CSV file:\n{e}")
            return

        if label_col not in df.columns:
            messagebox.showerror("Missing Column", f"'{label_col}' not found in the dataset.")
            return

        min_samples = min_samples_var.get()
        class_counts = df[label_col].value_counts()
        valid_classes = class_counts[class_counts >= min_samples].index
        df = df[df[label_col].isin(valid_classes)].copy()

        drop_cols = ['ID', 'Source', 'IUPACname(optional)', 'Glycanannotation2', 'GlyToucan ID', 'unique_ID']
        X = df.drop(columns=[col for col in drop_cols if col in df.columns] + [label_col], errors='ignore')
        y = df[label_col]

        if y.dtype == 'object':
            le = LabelEncoder()
            y = le.fit_transform(y)
            # Save encoder for later decoding
            le_path = os.path.splitext(train_csv_path)[0] + "_labelencoder.joblib"
            joblib.dump(le, le_path)
        try:
            test_ratio = test_split_var.get()
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_ratio, random_state=42)
            model = RandomForestClassifier(n_estimators=100, random_state=42)
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            #class_names = le.classes_ #so the encoded classes are restored
            #report = classification_report(y_test, y_pred, target_names=class_names)
            used_labels = unique_labels(y_test, y_pred)
            used_class_names = [le.classes_[i] for i in used_labels]
            report = classification_report(y_test, y_pred, target_names=used_class_names)
            #report = classification_report(y_test, y_pred)

            messagebox.showinfo("Training Complete", f"Random Forest trained successfully.\n\n{report}")

            base = os.path.splitext(train_csv_path)[0]
            model_path = base + "_rf_model.joblib"
            report_path = base + "_rf_performance.txt"
            joblib.dump(model, model_path)
            with open(report_path, "w") as f:
                f.write(report)

            if linked_exp_json and os.path.exists(linked_exp_json):
                try:
                    with open(linked_exp_json, "r") as f:
                        exp_data = json.load(f)
                    exp_data["train_parameters"] = {
                        "split_ratio": test_ratio,
                        "min_samples": min_samples,
                        "model_path": model_path,
                        "report_path": report_path
                    }
                    with open(linked_exp_json, "w") as f:
                        json.dump(exp_data, f, indent=4)
                except Exception as e:
                    print("Failed to write training parameters to exp.json:", e)

        except Exception as e:
            messagebox.showerror("Training Failed", str(e))


    def select_model_file():
        nonlocal model_file_path
        path = filedialog.askopenfilename(filetypes=[("Model files", "*.joblib *.pkl")])
        if path:
            model_file_path = os.path.abspath(path)
            messagebox.showinfo("Model Loaded", f"Model loaded from:\n{model_file_path}")

    def select_predict_input():
        nonlocal predict_input_path
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            predict_input_path = os.path.abspath(path)
            messagebox.showinfo("Input File Selected", f"Data loaded from:\n{predict_input_path}")

    def create_unlabeled_dataset():
        messagebox.showinfo("Not Yet Implemented", "This feature will allow you to select an experiment and automatically create a feature-matched dataset from its annotation and early raw-converted CSV.")

    def run_prediction():
        try:
            import joblib
        except ImportError:
            messagebox.showerror("Missing Dependency", "joblib is not installed.")
            return

        if not model_file_path or not predict_input_path:
            messagebox.showwarning("Missing Info", "Please select both a model file and an input CSV file.")
            return

        try:
            model = joblib.load(model_file_path)
            df = pd.read_csv(predict_input_path)

            # Try loading label encoder if available
            le_path = model_file_path.replace("_rf_model.joblib", "_labelencoder.joblib")
            if os.path.exists(le_path):
                le = joblib.load(le_path)
            else:
                le = None
        except Exception as e:
            messagebox.showerror("Load Error", str(e))
            return

        try:
            X = df.drop(columns=["MS2scan_no"], errors="ignore")
            y_pred = model.predict(X)

            # If encoder available, decode
            if le is not None:
                try:
                    y_pred = le.inverse_transform(y_pred)
                except:
                    pass
            elif hasattr(model, 'classes_'):
                y_pred = [model.classes_[i] if isinstance(i, int) else i for i in y_pred]

            df['Predicted_Label'] = y_pred
            out_path = os.path.splitext(predict_input_path)[0] + "_predicted.csv"
            df.to_csv(out_path, index=False)
            messagebox.showinfo("Prediction Complete", f"Predictions saved to:\n{out_path}")
        except Exception as e:
            messagebox.showerror("Prediction Failed", str(e))


    #creating unlabeled datasets
    def extract_fragment_masses(excel_path, sheet_name="ionlist"):
        xls = pd.ExcelFile(excel_path, engine="openpyxl")
        #sheet_name = ion_sheet_name #or xls.sheet_names[1]
        ion_df = xls.parse(sheet_name)
        fragment_masses  = ion_df["mass"].dropna().astype(float).tolist()#df.iloc[:, 0].dropna().astype(float).tolist()
        return fragment_masses

    def extract_ion_intensities(tsv_path, fragment_masses, ppm=20.0):
        df = pd.read_csv(tsv_path, sep="\t")
        result = []
        #base_cols = df[["protonatedmass", "MS2scan_no"]].copy()
        import numpy as np
        for idx, row in df.iterrows():
            try:
                peaks = list(eval(row["peaklist"]))
                intensities = list(eval(row["peakintensity"]))
            except Exception:
                continue

            peak_array = np.array(peaks)
            intensity_array = np.array(intensities)
            feature_row = {
                "protonatedmass": row["protonatedmass"],
                "MS2scan_no": int(row["MS2scan_no"])
            }

            for target in fragment_masses:
                ppm_tol = target * ppm / 1e6
                mask = np.abs(peak_array - target) <= ppm_tol
                if np.any(mask):
                    intensity = float(np.max(intensity_array[mask]))
                    feature_row[str(target)] = np.log10(intensity + 1)
                else:
                    feature_row[str(target)] = 0.0

            result.append(feature_row)

        return pd.DataFrame(result)

    def create_unlabeled_dataset():
        exp_path = filedialog.askopenfilename(filetypes=[("Experiment JSON", "*.json")])
        if not exp_path:
            return
        try:
            with open(exp_path, "r", encoding="utf-8") as f:
                exp = json.load(f)
        except Exception as e:
            messagebox.showerror("Failed to Load JSON", str(e))
            return

        ppm_value = simpledialog.askinteger("PPM Tolerance", "Enter PPM tolerance (e.g., 20):", minvalue=1, maxvalue=100)
        if ppm_value is None:
            return

        for sample_id, sample_info in exp.get("samples", {}).items():
            raw_csv = sample_info.get("csv")
            annotation_excel = sample_info.get("excel")
            if not (raw_csv and annotation_excel):
                continue
            try:
                ion_list = extract_fragment_masses(annotation_excel)
                feature_df = extract_ion_intensities(raw_csv, ion_list, ppm=ppm_value)
                out_path = os.path.splitext(raw_csv)[0] + f"_unlabeled_ppm{ppm_value}.csv"
                feature_df.to_csv(out_path, index=False)
                sample_info["unlabeled_dataset"] = out_path
            except Exception as e:
                print(f"[ERROR] Failed on sample '{sample_id}':", e)

        exp["prediction_parameters"] = {"ppm": ppm_value}
        try:
            with open(exp_path, "w", encoding="utf-8") as f:
                json.dump(exp, f, indent=4)
        except Exception as e:
            messagebox.showerror("Failed to Save JSON", str(e))
            return
        messagebox.showinfo("Done", f"Unlabeled datasets saved and experiment file updated.")
    # Predict tab state
    model_file_path = None
    predict_input_path = None



    train_csv_path = None
    linked_exp_json = None
    test_split_var = tk.DoubleVar(value=0.25)
    min_samples_var = tk.IntVar(value=5)

    subwin = tk.Toplevel(root)
    subwin.title("ML Analysis")
    subwin.geometry("650x520")

    notebook = ttk.Notebook(subwin)
    notebook.pack(fill="both", expand=True)

    # ... [no changes below this line: GUI layout remains as-is]


    # --- Tab 1: Train Model ---
    train_tab = ttk.Frame(notebook)
    notebook.add(train_tab, text="Train Model")

    tk.Label(train_tab, text="Step 1: Load Trainable Dataset (.csv)").grid(row=0, column=0, sticky="w", padx=10, pady=5)
    train_load_button = tk.Button(train_tab, text="Select CSV File", command=select_train_csv)
    train_load_button.grid(row=0, column=1, padx=5, pady=5)

    tk.Label(train_tab, text="Step 2: Select Label Column").grid(row=1, column=0, sticky="w", padx=10, pady=5)
    label_dropdown = ttk.Combobox(train_tab, values=['Structure', 'IUPACname(optional)', 'Glycanannotation2', 'GlyToucan ID'])
    label_dropdown.set("Structure")
    label_dropdown.grid(row=1, column=1, padx=5, pady=5)

    tk.Label(train_tab, text="Step 3: Choose Classifier").grid(row=2, column=0, sticky="w", padx=10, pady=5)
    classifier_var = tk.StringVar(value="rf")
    rf_button = tk.Radiobutton(train_tab, text="Random Forest (✔ functional)", variable=classifier_var, value="rf")
    xgb_button = tk.Radiobutton(train_tab, text="XGBoost (placeholder)", variable=classifier_var, value="xgb")
    svm_button = tk.Radiobutton(train_tab, text="SVM (placeholder)", variable=classifier_var, value="svm")
    knn_button = tk.Radiobutton(train_tab, text="KNN (placeholder)", variable=classifier_var, value="knn")

    rf_button.grid(row=2, column=1, sticky="w")
    xgb_button.grid(row=3, column=1, sticky="w")
    svm_button.grid(row=4, column=1, sticky="w")
    knn_button.grid(row=5, column=1, sticky="w")

    tk.Label(train_tab, text="Step 4: Train/Test Parameters").grid(row=6, column=0, sticky="w", padx=10, pady=5)
    tk.Button(train_tab, text="Set Parameters / Train the Model", command=open_train_settings).grid(row=6, column=1, padx=5, pady=5)

    train_button = tk.Button(train_tab, text="Train Model", command=train_model, bg="#CCE5FF")
    train_button.grid(row=7, column=0, columnspan=2, pady=10)

    tk.Label(train_tab, text="Trainable File Info (Origin Tracking)").grid(row=8, column=0, columnspan=2, sticky="w", padx=10, pady=(15, 5))
    origin_info = tk.Text(train_tab, height=4, width=70, state="disabled", wrap="word")
    origin_info.grid(row=9, column=0, columnspan=2, padx=10, pady=5)
    # -- Training tab and prediction tab UI (end reminder buttons) --
    tk.Label(train_tab, text="(🔜) Combine Datasets for Training").grid(row=10, column=0, columnspan=2, sticky="w", padx=10, pady=(15, 5))
    tk.Button(train_tab, text="[Placeholder] Combine Datasets").grid(row=11, column=0, columnspan=2, padx=10, pady=5)

    # --- Tab 2: Predict ---
    predict_tab = ttk.Frame(notebook)
    notebook.add(predict_tab, text="Predict")

    tk.Label(predict_tab, text="Step 1: Load Trained Model (.joblib/.pkl)").grid(row=0, column=0, sticky="w", padx=10, pady=5)
    predict_model_button = tk.Button(predict_tab, text="Select Model File", command=select_model_file)
    predict_model_button.grid(row=0, column=1, padx=5, pady=5)

    tk.Label(predict_tab, text="Step 2: Prepare Unlabeled Input Dataset").grid(row=1, column=0, sticky="w", padx=10, pady=5)
    create_unlabeled_button = tk.Button(predict_tab, text="Create unlabeled dataset of certain experiment", command=create_unlabeled_dataset)
    create_unlabeled_button.grid(row=2, column=0, columnspan=2, padx=10, pady=5)

    predict_input_button = tk.Button(predict_tab, text="Load unlabeled dataset", command=select_predict_input)
    predict_input_button.grid(row=3, column=0, columnspan=2, padx=10, pady=5)

    predict_button = tk.Button(predict_tab, text="Run Prediction", bg="#D5F5E3", command=run_prediction)
    predict_button.grid(row=4, column=0, columnspan=2, pady=10)

    # -- Training tab and prediction tab UI (end reminder buttons) --
    tk.Label(predict_tab, text="(🔜) Combine Datasets for Prediction").grid(row=5, column=0, columnspan=2, sticky="w", padx=10, pady=(15, 5))
    tk.Button(predict_tab, text="[Placeholder] Combine Datasets").grid(row=6, column=0, columnspan=2, padx=10, pady=5)

    close_button = tk.Button(subwin, text="Close", command=subwin.destroy)
    close_button.pack(pady=5)












def on_closing():
    if messagebox.askokcancel("Quit", "Do you really want to quit?"):
        logger.log("Application closed by user.")
        root.destroy()  # Clean exit
    #if logger.entries:
        #logger.save("autosave.log", include_debug=True)
        #print("Log saved to autosave.log")

def save_log_to_file():
    filename = filedialog.asksaveasfilename(defaultextension=".log", filetypes=[("Log files", "*.log")])
    if filename:
        logger.save(filename, include_debug=True)
        messagebox.showinfo("Log Saved", f"Log saved to:\n{filename}")


#icon
ico_path = os.path.join(os.path.dirname(__file__), 'GlycoMSPlogo.ico')

# Initialize GUI
root = tk.Tk()
root.iconbitmap(default=ico_path)
root.protocol("WM_DELETE_WINDOW", on_closing)
root.title("GlycoMSP File Manager GUI v0.4a")
root.geometry("840x600")
root.minsize(840, 600)



# Buttons
button_frame = tk.Frame(root)
button_frame.pack(pady=10)

tk.Button(button_frame, text="Select Raw File", command=lambda: select_file("raw")).grid(row=0, column=0, padx=5)
#tk.Button(button_frame, text="Select mzML File", command=lambda: select_file("mzml")).grid(row=0, column=1, padx=5)
tk.Button(button_frame, text="Select CSV File", command=lambda: select_file("csv")).grid(row=0, column=2, padx=5)
tk.Button(button_frame, text="Select Excel File", command=lambda: select_file("excel")).grid(row=0, column=3, padx=5)
tk.Button(button_frame, text="Clear All", command=clear_files).grid(row=0, column=4, padx=5)
tk.Button(root, text="Save Log", command=save_log_to_file).pack(pady=5)
tk.Button(root, text="About", command=open_about_window).pack(pady=5)
# Text widget to log selected files
text_widget = tk.Text(root, height=15, width=80)
text_widget.pack(pady=10)


#status bar?
status_var = tk.StringVar()
status_var.set("Idle")
status_label = tk.Label(root, textvariable=status_var, fg="blue")
status_label.pack(pady=5)
#progress bar?
progress = ttk.Progressbar(root, orient="horizontal", mode="indeterminate", length=250)
progress.pack(pady=5)


# --- Add Analysis Tools Frame ---

analysis_frame = tk.LabelFrame(root, text="Analysis Tools", padx=10, pady=10)
analysis_frame.pack(padx=10, pady=10, fill="x")

convert_button = tk.Button(analysis_frame, text="Convert Raw to CSV",command=launch_metadata_batch)
convert_button.pack(side="left", padx=5)
tk.Button(analysis_frame, text="Prepare Dataset", command=open_prepare_dataset_window).pack(side="left", padx=5)
tk.Button(analysis_frame, text="Run ML Analysis", command=open_ml_analysis_window).pack(side="left", padx=5)


# Run the GUI
root.mainloop()



