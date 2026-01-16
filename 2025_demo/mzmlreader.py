import os
import numpy as np
version = 0.5
last_update = 20260116
#data-loading and saving
import glob
import os
import pathlib
from pathlib import Path
import os
import pathlib
import tempfile
import time
import json

import pymzml

#changelog:
#v0.5 functional mzmlreader (proto) - may not be flexible on different parameters during extraction

# Load an mzML file
#run = pymzml.run.Reader(r"C:\Users\Sakazuki\Downloads\zf_sPerMeNG_brain.mzML")


"""
#for debug
import pydoc
doc = pydoc.render_doc(type(spectrum), "Help on %s")
with open("spectrum_help.txt", "w", encoding="utf-8") as f:
    f.write(doc)
print("Wrote spectrum_help.txt")
#help(spectrum)
"""

"""
# Iterate through the spectra
t = 0
print("test function here")
for spectrum in run:
    if t > 0:
        break
    else:
        if spectrum['ms level'] == 2:
            continue
            print("ID:", spectrum.ID)
            #print("selected_precursors:", spectrum.selected_precursors)
            id = spectrum.ID
            precursor_list = spectrum.selected_precursors 

            ms1mz = precursor_list['mz']
            ms1intensity = precursor_list['i']
            charge = precursor_list['charge']
            ms1scan = precursor_list['precursor id']

            empty_containerpeaks = []
            empty_containerintensity = []
            print(f"Getting peak list")
            p_peaklist = spectrum.mz
            p_intensity = spectrum.i
            #print(f"type of peaklist is {type(p_peaklist)} and the output is {p_peaklist}")
            for peaks in p_peaklist:
                #print(f"peaks are in type of {type(peaks)} and the data inside is {peaks}")
                #peaks are in type of <class 'numpy.float32'> and the data inside is 97.00896453857422
                empty_containerpeaks.append(np.round(peaks, 4))
            print(f"peaks after pre-processing rounding is {empty_containerpeaks}")
            for intensities in p_intensity:
                empty_containerintensity.append(np.round(intensities, 4))
            print(f"intensity after pre-processing rounding is {empty_containerintensity}")


            #peaklist = p_peaklist.split()
            #intensiy = p_intensity.split()
            #print(f"peaklist is {peaklist} and intensity is {intensiy}")
            t += 1
        elif spectrum['ms level'] == 3:
            continue
                #if spectrum.ms_level == 3:
            print("=== MS3 spectrum found ===")
            print("ID:", spectrum.ID)
            print("index:", spectrum.index)
            print(f"Getting peak list")
            p_peaklist = spectrum.mz
            p_intensity = spectrum.i
            print("peak",p_peaklist)
            print("intensity",p_intensity)
            print("precursors",spectrum.selected_precursors)
            break
        #print("\n--- Public attributes / methods ---")
        #print([n for n in dir(spectrum) if not n.startswith("_")])
"""

#####functional blocks imported from msprawrxtactor.py#####
###consider isolate the metadata function to another pyhton file###
###copied 20260115, rawextractor version = 0.53, last_update = 20250927###
def fillexpinfo(mzmlfilepath, loadfile=None, debug=False):
    if not loadfile:
        print("You imported experiment description file")
    else:
        if debug:
            print("[debug]Start filling new experiment information for metadata")
    expTitle = input("Enter title of this experiment. For example 'human T cells'.\n")
    expDescription = input("Enter details of this experiment. For example 'SiaT KO test'\n")
    expAuthor = input("Enter your name so other knows who did the analysis.\n")
    expRawdate = input("Enter the date when the raw file was obtained in YYYY/MM/DD format. For example '20240229'.\n")
    expglycantype = input("Enter glycan type of analysis. (N/O/GL/N+O/others).\n")
    expmsmode = input("Enter mass spec detection mode (+/-)")
    #expmsinstrument = input("Enter mass spec instrument name")

    #autofill
    parameters = ["[debug]Autofill the version info"]
    extractdate = time.strftime("%Y%m%d") #20001105 for example
    rawfile = mzmlfilepath #filename.raw (str)
    rawfilename =  os.path.splitext(mzmlfilepath)[0] # filename (str)

    #define metadata
    metadata = {
        "Experiment Title": expTitle,
        "Experiment Description": expDescription,
        "Author running this analysis": expAuthor,
        "Raw data acquired date": expRawdate,
        "Glycan Type": expglycantype,
        "Mass Analyzer charge mode": expmsmode,
        "Parameters when GlycoMSP launched":parameters,
        "Date of file extracted from raw file": extractdate,
        "Original raw file path":rawfile,
        "Raw filename":rawfilename
        }
    return metadata
    
    #define project name for saving metadata and log file #
def generatefilename(rawfilepath):
    mansavename = input("Enter the name you want to save for this project or batch of analysis.\n")
    rawfilename =  os.path.basename(rawfilepath) # filename (str) #os.path.splitext(rawfilepath)[0]
    extractdate = time.strftime("%Y%m%d") #20001105 for example
    savename = f"{mansavename}_{extractdate}_{rawfilename}"
    return savename

def savemetadata(metadata, filepath):
    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=4)
    return filepath

def save_log(filepath, log_entries, debug_logs=None):
    with open(filepath, 'w', encoding='utf-8') as f:
        for entry in log_entries:
            f.write(f"{entry}\n")
        if debug_logs:
            f.write("-----\ndebug enabled-----\n")
            for logs in debug_logs:
                f.write(f"{logs}\n")
            f.write("-----\ndebug log ends-----\n")
    return filepath

def calcproton(isolatedmass, charge):
    if charge == 1 :   #change abs(charge) to charge to negative can be calculated too, 
        #A tiny nuance to keep in mind as we go forward: Thermo’s API typically returns the magnitude of charge (often positive) 
        # even for negative-mode data, depending on method/instrument
        return isolatedmass
    elif charge >= 2:
        conv = (isolatedmass * charge - (charge - 1) * 1.00784)
        return conv
    elif charge == 0:
        return 0
    elif charge < -1:
        negconv = (isolatedmass * charge + (charge - 1) * 1.00784)
        return negconv

try:
    import xml.etree.ElementTree as ET
    def mzml_spectrum_count(mzml_path: str) -> int | None:
        # iterparse is streaming; we exit early once spectrumList is found
        for event, elem in ET.iterparse(mzml_path, events=("start",)):
            # elem.tag may include a namespace like {http://psi.hupo.org/ms/mzml}spectrumList
            if elem.tag.endswith("spectrumList"):
                cnt = elem.attrib.get("count")
                return int(cnt) if cnt is not None else None
        return None
except:
    def mzml_spectrum_count(mzml_path):
        return None

def mzml_file_examination(mzmlfile, auto = True, time = 0.0,debug=False):
    if not auto:
        print(f'This raw file has {mzml_spectrum_count(mzmlfile)} spectra')
        scan_number = int(input('enter spectrum no you want to summarize'))
        maxn = int(mzml_spectrum_count(mzmlfile))
    else:
        scan_number = int(mzml_spectrum_count(mzmlfile))
        maxn = scan_number
    return maxn 


def peaklist_validation(peaklist, peakintensity, debug = False):

    if len(peaklist) == len(peakintensity):
        return True, None
    else:
        if debug:
            peaklist_length = [len(peaklist), len(peakintensity)]
            print(f"[debug] peaklist and intensity are not matched, {peaklist_length}")
            return False, peaklist_length
        return False, None    

def save_to_csv(csvname, data, debug=False):
    """
    Atomically write a tab-separated CSV:
    - Write to a temp file in the same directory, then os.replace to final.
    - Keeps the legacy behavior of using a *stem* (no .csv) for csvname.
    """
    exportcsv = f"{csvname}.csv"
    parent = os.path.dirname(exportcsv) or "."
    os.makedirs(parent, exist_ok=True)

    perf_esti2 = time.time()
    # Create a temp file in the same directory so os.replace is atomic on the same volume
    with tempfile.NamedTemporaryFile("w", delete=False, dir=parent, suffix=".tmp") as tmp:
        # Write rows as tab-separated (matches legacy)
        # NOTE: data is a list of row-like items, first row may be header tuple
        for row in data:
            tmp.write('\t'.join(map(str, row)) + '\n')
        tmp.flush()
        os.fsync(tmp.fileno())
        tmp_path = tmp.name

    # Atomic promotion to final name (no second write)
    os.replace(tmp_path, exportcsv)

    if debug:
        csvwritetime = time.time() - perf_esti2
        print(f"[debug] Finished exported {exportcsv} from raw data")
        print(f"[debug] Time spent for extraction: {csvwritetime:.3f} seconds.")
    return exportcsv

#############################
#########end of copy#########

def extract_mzML(mzmlfilepath=None, outdir=None, round=False, debug = False):
    if debug:
        print("[debug] need to create debug logs when debug is set to True.")
        print("[debug] Also need a readable log that records needed information and the operation can be reproduced when loading that record.")
        print("[debug add debug log and record log in debug block]")
        print("[debug]: start calculate time spent for converting data")
        perf_esti1 = time.time() #first time stamp for evaluating performance
    try:
        mzmlfile = pymzml.run.Reader(mzmlfilepath)
    except Exception as e:
        print(f"[Error] MSFileReader failed: {e}")
        return False
    # Default to the raw file’s folder, else CWD
    if not outdir:
        outdir = os.path.dirname(mzmlfilepath) or os.getcwd()

    #init copied from msprawextractor
    scan_number = mzml_file_examination(mzmlfilepath)
    #init variables for storing parsed data
    ms1list, ms2list, ms3list, errlist , unmatchlist, logs = [], [], [], [], [], []
    #from v0.4 named Tuple of (peaklist,intensity) were changed to 2 separated lists. Include length validation 
    ms2data, ms3data = [], []  #b for MS2, c for MS3 < replacing
    fixms2header = ('entry_no','MS1scan_no', 'MS1_isolationmass', 'MS1_monoisolationmass','chargeState','protonatedmass',
           'MS2scan_no', 'peaklist', 'peakintensity')
    ms2data.append(fixms2header)
    fixms3header =  ('entry no','MS3scan_no', 'MS2_isolationmass', 'MS2_monoisomass', 'MS2scan_no', 'peaklist', 'peakintensity')
    ms3data.append(fixms3header)
    ms2count, ms3count = 1,1  
    print(f"[dev] scan number is {scan_number}")
    for spectrum in mzmlfile:
        if spectrum['ms level'] == 1:
            if debug:
                ms1list.append(spectrum.ID)
        elif spectrum['ms level'] == 2:
            if debug:
                ms2list.append(spectrum.ID)
            #print("ID:", spectrum.ID)
            #print("selected_precursors:", spectrum.selected_precursors)
            id = spectrum.ID
            precursor_list = spectrum.selected_precursors[0]
            ms1mz = precursor_list['mz']
            ms1intensity = precursor_list.get("i",0.0)#precursor_list['i']
            charge = precursor_list['charge']
            ms1scan = precursor_list['precursor id']
            empty_containerpeaks = []
            empty_containerintensity = []
            #print(f"Getting MS2 peak list")
            p_peaklist = spectrum.mz
            p_intensity = spectrum.i
            #print(f"type of peaklist is {type(p_peaklist)} and the output is {p_peaklist}")
            if round:
                print("[dev]Perform round to 4 digits in data extraction")
                print("[dev]Need to save log of this operation")
                for peaks in p_peaklist:
                    #print(f"peaks are in type of {type(peaks)} and the data inside is {peaks}")
                    #peaks are in type of <class 'numpy.float32'> and the data inside is 97.00896453857422
                    empty_containerpeaks.append(np.round(peaks, 4))
                #relocate back, update the peaklist
                #warning, check the DataType
                p_peaklist = empty_containerpeaks
                print(f"peaks after pre-processing rounding is {empty_containerpeaks}")
                for intensities in p_intensity:
                    empty_containerintensity.append(np.round(intensities, 4))
                print(f"intensity after pre-processing rounding is {empty_containerintensity}")
                p_intensity = empty_containerintensity
            #mass post calculation
            protonatedmass = calcproton(ms1mz,charge)
            ms2spectrumdata = (
                ms2count,
                ms1scan,#parent scan
                ms1mz,#isolation mass
                0, #monoisotopic isolation mass if present
                charge,
                protonatedmass,
                spectrum.ID,
                p_peaklist,
                p_intensity
            )
            ms2data.append(ms2spectrumdata)
            ms2count += 1
            peakcheck, failedlist = peaklist_validation(p_peaklist, p_intensity, debug = True)
            if not peakcheck: #if the peaklist and intensity aren't matched, not False = True
                reportunmatch = [spectrum.ID, failedlist]
                unmatchlist.append(reportunmatch)
            #print(f"[debug] ms2data = {ms2data}, will break later")   
        elif spectrum['ms level'] == 3:
            id = spectrum.ID
            if debug:
                ms3list.append(id)
            precursor_list = spectrum.selected_precursors[0]
            ms2mz = precursor_list['mz']
            ms2intensity = precursor_list.get('i', 0.0)  # likely None for MS3
            charge = precursor_list.get('charge', 0)
            ms2scan = precursor_list['precursor id']
            empty_containerpeaks = []
            empty_containerintensity = []
            #print(f"Getting MS3 peak list")
            p_peaklist = spectrum.mz
            p_intensity = spectrum.i
            #print(f"type of peaklist is {type(p_peaklist)} and the output is {p_peaklist}")
            if round:
                print("[dev]Perform round to 4 digits in data extraction")
                print("[dev]Need to save log of this operation")
                for peaks in p_peaklist:
                    #print(f"peaks are in type of {type(peaks)} and the data inside is {peaks}")
                    #peaks are in type of <class 'numpy.float32'> and the data inside is 97.00896453857422
                    empty_containerpeaks.append(np.round(peaks, 4))
                #relocate back, update the peaklist
                #warning, check the DataType
                p_peaklist = empty_containerpeaks
                print(f"peaks after pre-processing rounding is {empty_containerpeaks}")
                for intensities in p_intensity:
                    empty_containerintensity.append(np.round(intensities, 4))
                print(f"intensity after pre-processing rounding is {empty_containerintensity}")
                p_intensity = empty_containerintensity
            ms3spectrumdata = (
                ms3count,
                ms2scan,
                ms2mz,
                0,
                spectrum.ID,
                p_peaklist,
                p_intensity
                )
            ms3data.append(ms3spectrumdata)
            ms3count +=1
        else:
            errlist.append(spectrum.ID)   
    print("Extraction finished.")
    #if possible, delete the mzml in memory here
    if unmatchlist is not []:
        logs.append("Unmatched peaklists are recorded as [[spectrum no], [peaklist entries, peak intensity entries]]")
        logs.append(unmatchlist)

    if debug:
        debug_headers = "Debug enabled. Debug log information: [MS1 list, MS2 list, MS3 list, exception list]"
        debug_logs = [debug_headers, ms1list, ms2list, ms3list, errlist]
    else:
        debug_logs = None
    if debug:
        extractiontime = time.time() - perf_esti1  #timestamp1: time spent for raw file extraction
        print("[debug]Raw file has been closed.")
        print("[debug]Preparing information for writing...")
        print(f"[debug]Time spent for extraction: {int(extractiontime)} seconds.")

    filename = os.path.splitext(os.path.basename(mzmlfilepath))[0]
    ms2output = "ms2tmp_" + filename
    ms3output = "ms3tmp_" + filename

    rawstem = pathlib.Path(mzmlfilepath).stem
    #v0.53
    ms2tmp_path = os.path.join(outdir, f"ms2tmp_{rawstem}.csv")
    ms3tmp_path = os.path.join(outdir, f"ms3tmp_{rawstem}.csv")
    # save_to_csv still takes a *stem*; pass the stem without .csv but with full dir
    save_to_csv(ms2tmp_path[:-4], ms2data, debug=debug)
    save_to_csv(ms3tmp_path[:-4], ms3data, debug=debug)
    versioninfo = ["extractor info", version, last_update]
    logs.append(versioninfo)
    print("Conversion finished from mspextractor.")

    return {"ms2tmp": ms2tmp_path, "ms3tmp": ms3tmp_path}

def finalize_extraction(rawfile, metadata, filename, debug=False):
    print(f"[Finalizing extraction] For: {filename}")
    print(f"[Metadata provided] Sample name: {metadata.get('Experiment Title', 'Unknown')}")

    # Optional: log the fact that finalization completed
    if debug:
        print("[Debug] finalize_extraction called successfully")

#testing line
#extract_mzML(mzmlfilepath=r"C:\Users\Sakazuki\Downloads\zf_sPerMeNG_brain.mzML", round=False, debug = False)

'''
for spectrum in run:
    if t > 3:
        break
    elif spectrum['ms level'] == 2:
        print(f"type of this spectrum is {type(spectrum)}")
        print(f"testing pymzml functionality")
        print(f"Spectrum ID: {spectrum.ID}, {type(spectrum.ID)}")
        print(f"MS Level: {spectrum.ms_level}, {type(spectrum.ms_level)}")
        print(f"Retention Time: {spectrum.scan_time_in_minutes()} minutes, {type(spectrum.scan_time_in_minutes)}")
        print(f"Precursor m/z: {spectrum.selected_precursors}, {type(spectrum.selected_precursors)}")
        print(f"Polarity: {spectrum.get('scan polarity')}, {type(spectrum.get('scan polarity'))}")
        print(f"Instrument Type: {spectrum.get('instrument model')}, {type(spectrum.get('instrument model'))}")
        # Access other metadata as needed
        print(spectrum.mz, spectrum.i, type(spectrum.mz), type(spectrum.i))  
        t+=1
#print("testing certain ms number")
#spectrum_with_id_311 = run[ 311 ]
#print(spectrum_with_id_311())
print("testing run itself")
spectrum = next(run)
print(dir(spectrum))
'''
print(f"Testing mzml reader finished")

#docs are summarized by GPT4o (20241104) I can't find clear instructions on internet 
