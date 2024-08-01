import os
version = 0.2
last_update = 20240730

# version: 0.2
# date: 20240730
# about this file: original trytolistoutheaders.py and GlycoMSP_demo.py modulated to support future standardized processing workflow
# import tkinter part
#mind that we plan to add GUI support in project managing, the tkinter detection may be moved to funtion: GUI in MSPinit.py in future
# re-organized the code structure and readability using copilot and ChatGPT4-o
# no sensitive contents were sent to the server for this part

try: 
    import tkinter as tk
    from tkinter import filedialog
except ImportError:
    raise ImportError('Please install tkinter to enable GUI-based file selection')
else:
    has_tkinter = True

# Function to create and hide Tkinter root window
# If future version supports GUI fully available, please reconstruct this part (not hiding main window, but think about what should be there)
def create_tkinter_root():
    root = tk.Tk()
    root.withdraw() # hide main window
    return root

# Function to validate file path
def validate_file_path(file_path):
    if not os.path.exists(file_path):
        print(f"Error: The file {file_path} does not exist.")
        return False
    if not os.path.isfile(file_path):
        print(f"Error: The path {file_path} is not a file.")
        return False
    return True

# Function to select raw file for further pre-processing
def select_raw_file():
    root = create_tkinter_root()
    rawfile_path = filedialog.askopenfilename(filetypes=[("RAW files", "*.raw")])

    if rawfile_path:
        if validate_file_path(rawfile_path):
            print(f"Selected raw file: {rawfile_path}")
        else:
            rawfile_path = None
    else:
        print("No raw file has been selected!")

    return rawfile_path

# Function to select csv files (annotation file OR processed/analyse-in-progress project file stored as csv)
def select_project_file():
    root = create_tkinter_root()
    csvfile_path = filedialog.askopenfilename(filetypes=[("Comma-separated values", "*.csv")])

    if csvfile_path:
        if validate_file_path(csvfile_path):
            print(f"Selected csv file: {csvfile_path}")
        else:
            csvfile_path = None
    else:
        print("No csv file has been selected!")

    return csvfile_path

# Dealing with loading functions
def fileloader(filetype):
    if not has_tkinter:
        print("Please install tkinter for GUI-supported file selection")
        filepath = None
    else:
        if filetype == "raw":
            filepath = select_raw_file()
        elif filetype == "csv":
            filepath = select_project_file()
    return filepath


'''
usage from other files:
import mspfileloader
mspfileloader.fileloader(type)
type could be string "raw" for Thermo raw file or "csv" for pre-processed file 

fileloader(raw):
                will return raw file path for further data pre-processing
fileloader(csv):
                will return csv file path for further analyse progress
when fileloader(any) returns None:
                tkinter not installed. Catch the None attribute as condition.
                
need to be validated:
DataType(filepath) when it is a path. And should I import os.path or pathlib            
'''


