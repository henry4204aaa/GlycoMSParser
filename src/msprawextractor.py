######version info#####

version= "0.61"
last_update = 20261023

#v0.61 READER-DIAG R2 (Codex round 2, 20261023): a CSV is named in the diagnostics only after it is written;
#      Ctrl-C/SystemExit also finish the diagnostics (re-raised unchanged); stage 'iterate' starts once the file is
#      open; the open-failure message is redacted; every exported MS2 status row records mass_mz_source (was only
#      rows without a reported charge). CSVs unchanged.
#v0.60 READER-DIAG CP-D: every conversion also writes <stem>.scan_status.tsv (one row per MS2/MS3 scan) and
#      <stem>.reader_diag.json/.txt (outcome, exit code, warnings, counts, coverage, failure) next to the CSVs, on
#      success and on failure (msp_spectrum_record v0.4). Optional stem= (output name; the CLI passes the alias) and
#      entry=. The return dict gains outcome, exit_code, warnings, diagnostics. CSVs unchanged. The MS order is
#      read once per scan (was up to three times), and MS3 scans read their filter for the polarity.

#v0.59 READER-DIAG CP-C.5: optional charge sources for MS2 scans without a valid reported charge (both off by default):
#      charge_comb (MS1 isotope comb, msp_charge_comb) then dataset_charge (one user value). A reported charge is never
#      overwritten. The mass uses the isolation m/z when the file has no monoisotopic m/z (was 0). Summary gains "charge".
#v0.58 READER-DIAG: calcproton takes the scan polarity; negative mode gives the [M-H]- equivalent (was the positive formula, (z-1)*2.015 Da low).
#v0.57 READER-DIAG: MS2/MS3 scans without peaks are not written to the CSV (counted in the summary).
#v0.56 READER-DIAG: MS2/MS3 peaks come from read_scan_peaks (FTMS scans: label data as before;
#      other analyzers: mass list, vendor-centroided when the scan is profile). MS3 rows read their own peaks.
#      MS3 header names fixed: column 2 is the parent MS2 scan, column 5 the MS3 scan (names were swapped).
#v0.55 code review 1. Remove unused blocks
#v0.53 fix tmp file not found issue when directly click the v10 python file (cwd not specified -> fall to default)
#v0.52 check and prepare for negative mode support
#v0.51 failsafe in case the Thermo Library or foundation is missing
#v0.5 trying to solve duplicated metadata and log issue
#v0.4 move version control out to prevent circular import for making metadata
#dynamic version control will be added soon (independently)
#load other files in this package
#import mspfileloader

#combining trytolistoutheaders.py and comparepeaklist

import configparser

import pandas as pd
import json

#handling mass peak list (~GlycoMSP v0.3/ Will be removed before 1.0)
#about to change the way storing it
from collections import namedtuple

#generate timestamp and perfermance evaluation
import time

#data-loading and saving
import glob
import os
import pathlib
from pathlib import Path
import os
import pathlib
import tempfile
import re

#suppress warning during data hanlding process
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#v0.59 READER-DIAG CP-C.5: shared charge rules and the MS1 isotope comb (new modules, standard library only)
import msp_charge_comb
from msp_spectrum_record import (normalize_charge, resolve_charge, validate_dataset_charge, CHARGE_SOURCES,
                                 CHARGE_SOURCE_REPORTED, CHARGE_SOURCE_NONE, ReaderDiagnostics)

#specific raw file extraction key module
try:
    from pymsfilereader import MSFileReader
    pymsreader = True
except:
    print("Please run this code on a Windows-based OS with pymsfilereader and Thermo library installed")
    pymsreader = False


config = configparser.ConfigParser()
config.read("project_version_manifest.ini")

#########import finished########

    #file metadata generation
    #some information can be retrieved when raw file is opened
def fillexpinfo(rawfilepath, loadfile=None, debug=False):
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
    rawfile = rawfilepath #filename.raw (str)
    rawfilename =  os.path.splitext(rawfilepath)[0] # filename (str)

    #define metadata
    metadata = {
        "json_type": "glycomsp.metadata",
        "schema_version": "1.0.0",
        #above are fixes for next version json update
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


    #returns True when a user decides to export as arff
def exportasarff(format):
    if format.lower() == 'a' or format.lower() == 'arff':
        arff = True
    else:
        arff = False
    return arff

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


#work work 
def save_to_arff():#(arffname, data, headers, attributes, debug=False):
    print("Waiting to be built...")
    #supposed to be added in 0.4 after csv files are validated okay.
#work work 


#calculates and returns protonated mass for precursor in each spetrum
#v0.58 READER-DIAG: negative-mode precursor mass. The in-silico library (compnewv4, mode '-') stores [M-H]-,
#and Thermo reports the charge magnitude in negative mode too, so the scan polarity decides the formula.
PROTON_MASS = 1.007276

def calcproton(isolatedmass, charge, polarity=None):
    if polarity == '-' or charge < 0:
        #[M-zH]z- -> [M-H]- equivalent: m/z * z + (z - 1) * proton. A negative charge value is read as negative mode.
        z = abs(charge)
        if z == 0:
            return 0
        return isolatedmass * z + (z - 1) * PROTON_MASS
    if charge == 1 :   #change abs(charge) to charge to negative can be calculated too, 
        #A tiny nuance to keep in mind as we go forward: Thermo’s API typically returns the magnitude of charge (often positive) 
        # even for negative-mode data, depending on method/instrument
        return isolatedmass
    elif charge >= 2:
        conv = (isolatedmass * charge - (charge - 1) * 1.00784)
        return conv
    elif charge == 0:
        return 0
    

#peak extractor parts, need modification

#READER-DIAG 20261013: where a scan's peaks come from.
#GetLabelData returns the FT label stream; ion-trap (ITMS) scans have none and come back empty,
#although their points are available from the mass list. The pilot LTQ probe confirmed this, and the
#vendor-centroided mass list is identical to the MSConvert "Thermo/Xcalibur peak picking" output.
PEAK_SOURCE_LABEL_DATA = "label_data"
PEAK_SOURCE_MASS_LIST = "mass_list"
PEAK_SOURCE_MASS_LIST_CENTROID = "mass_list_vendor_centroid"
_FT_ANALYZER_RE = re.compile(r"\bFTMS\b")
_POLARITY_RE = re.compile(r"(?:^|\s)([+-])\s+[pc](?:\s|$)")   # "ITMS - p ESI ...", "FTMS {1,1} + c NSI ..."

def scan_polarity(scan_filter):
    """'+', '-' or None, from the polarity token of a Thermo scan filter."""
    match = _POLARITY_RE.search(scan_filter or "")
    return match.group(1) if match else None

def read_scan_peaks(rawfile, scan, scan_filter=None):
    """Return (peaklist, peakintensity, peak_source) for one native scan.

    FTMS scans keep the label-data path unchanged. Any other analyzer reads the mass list with no
    filter, no intensity cutoff and no peak limit; a profile scan is centroided by the vendor library
    (default peak width), a centroid scan is returned as stored.
    """
    if scan_filter is None:
        scan_filter = rawfile.GetFilterForScanNum(scan) or ""
    if _FT_ANALYZER_RE.search(scan_filter):
        labels = rawfile.GetLabelData(scan)[0]
        return labels[0], labels[1], PEAK_SOURCE_LABEL_DATA
    centroid = bool(rawfile.IsProfileScanForScanNum(scan))
    masslist = rawfile.GetMassListFromScanNum(scan, scanFilter='', intensityCutoffType=0, intensityCutoffValue=0,
                                              maxNumberOfPeaks=0, centroidResult=centroid, centroidPeakWidth=0.0)[0]
    return masslist[0], masslist[1], (PEAK_SOURCE_MASS_LIST_CENTROID if centroid else PEAK_SOURCE_MASS_LIST)

def _new_peak_summary():
    return {"ms2": {"by_source": {}, "exported": 0, "empty": []}, "ms3": {"by_source": {}, "exported": 0, "empty": []}}

def _tally_peaks(block, scan, peaklist, peakintensity, source):
    """Count one scan; return False when it has no peaks and is not to be exported (as MSConvert drops it)."""
    block["by_source"][source] = block["by_source"].get(source, 0) + 1
    if len(peaklist) == 0 and len(peakintensity) == 0:
        block["empty"].append(scan)
        return False
    block["exported"] += 1
    return True

#READER-DIAG CP-C.5 (20261015): charge for MS2 scans without a valid reported charge.
#Order (Henry): reported charge, never overwritten -> comb on the parent MS1 profile (charge_comb=True)
#-> the user's dataset charge -> unknown (0, as before). The comb runs only on non-FT profile MS1, the one
#condition it was tried on; each parent MS1 is read once.
COMB_NOT_RUN = "not_run"                 # charge reported, or the comb is off
COMB_NO_ISOLATION_MZ = "no_isolation_mz"
COMB_NO_PARENT_MS1 = "no_parent_ms1"
COMB_NO_PROFILE_MS1 = "no_profile_ms1"   # parent MS1 stored as centroids
COMB_FT_MS1_NOT_TRIED = "ft_ms1_not_tried"
MASS_MZ_MONOISOTOPIC = "monoisotopic"
MASS_MZ_ISOLATION = "isolation"

def _parent_ms1_profile(rawfile, parent, cache):
    """(mz, intensity) of the parent MS1 profile, or a comb status string. Read once per parent scan."""
    if parent in cache:
        return cache[parent]
    try:
        scan = int(parent)
        if scan <= 0 or rawfile.GetMSOrderForScanNum(scan) != 1:
            result = COMB_NO_PARENT_MS1
        elif _FT_ANALYZER_RE.search(rawfile.GetFilterForScanNum(scan) or ""):
            result = COMB_FT_MS1_NOT_TRIED
        elif not rawfile.IsProfileScanForScanNum(scan):
            result = COMB_NO_PROFILE_MS1
        else:
            masslist = rawfile.GetMassListFromScanNum(scan, scanFilter='', intensityCutoffType=0, intensityCutoffValue=0,
                                                      maxNumberOfPeaks=0, centroidResult=False, centroidPeakWidth=0.0)[0]
            result = (masslist[0], masslist[1])
    except Exception:
        result = COMB_NO_PARENT_MS1
    cache[parent] = result
    return result

def _new_charge_summary(charge_comb, dataset_charge):
    return {"settings": {"charge_comb": charge_comb, "dataset_charge": dataset_charge,
                         "comb": msp_charge_comb.parameters() if charge_comb else None},
            "sources": {k: 0 for k in CHARGE_SOURCES}, "comb_status": {}, "comb_vs_dataset_disagree": 0,
            "mass_mz_source": {MASS_MZ_MONOISOTOPIC: 0, MASS_MZ_ISOLATION: 0},
            "scans": []}   # one entry per exported MS2 scan without a valid reported charge

def _resolve_scan_charge(rawfile, scan, info, charge_comb, dataset_charge, ms1_cache, charge_summary):
    """Return (charge cell, mass m/z, mass m/z source) for one MS2 scan and record the decision in charge_summary."""
    isolation_mz, mono_mz, reported_cell, parent = info[0], info[1], info[2], info[3]
    mass_mz, mz_source = (mono_mz, MASS_MZ_MONOISOTOPIC) if (mono_mz or 0) > 0 else (isolation_mz, MASS_MZ_ISOLATION)
    charge_summary["mass_mz_source"][mz_source] += 1
    reported = normalize_charge(reported_cell)
    comb, status = None, COMB_NOT_RUN
    if reported is None and charge_comb:
        if not (isolation_mz or 0) > 0:
            status = COMB_NO_ISOLATION_MZ
        else:
            profile = _parent_ms1_profile(rawfile, parent, ms1_cache)
            if isinstance(profile, str):
                status = profile
            else:
                comb = msp_charge_comb.infer_charge(profile[0], profile[1], isolation_mz)
                status = comb.status
    comb_z = comb.z if comb is not None else None
    charge, source = resolve_charge(reported, comb_z, dataset_charge)
    charge_summary["sources"][source] += 1
    charge_summary["comb_status"][status] = charge_summary["comb_status"].get(status, 0) + 1
    if comb_z is not None and dataset_charge is not None and comb_z != dataset_charge:
        charge_summary["comb_vs_dataset_disagree"] += 1
    if source != CHARGE_SOURCE_REPORTED:
        entry = {"scan": scan, "charge_source": source, "charge": charge or 0, "mass_mz_source": mz_source,
                 "comb_status": status}
        if comb is not None and comb.score is not None:
            entry.update(comb_z=comb_z, comb_best_z=max(comb.scores, key=lambda z: comb.scores[z][0]),
                         comb_score=round(comb.score, 4), comb_margin=round(comb.margin, 4),
                         comb_anchor=round(comb.anchor, 4))
        charge_summary["scans"].append(entry)
    # A reported charge, or none at all, is written exactly as the file gave it (byte-identical rows).
    cell = reported_cell if source in (CHARGE_SOURCE_REPORTED, CHARGE_SOURCE_NONE) else charge
    return cell, mass_mz, mz_source

#return scan number limit
#define number of spectra to parse (TODO: make time exclusion convertible and need to be added to record file.)

def raw_file_examination(rawfile, auto = True, time = 0.0,debug=False):
    #time and debug are placeholders
    if not auto:
        print(f'This raw file has {rawfile.GetNumSpectra()} spectra')
        scan_number = int(input('enter spectrum no you want to summarize'))
        maxn = int(rawfile.GetNumSpectra())
    else:
        scan_number = int(rawfile.GetNumSpectra())
        maxn = scan_number
    return maxn 
#maxn defines as upper limit for spectrum extraction

####construction####

def peaklist_validation(peaklist, peakintensity, debug = False):

    if len(peaklist) == len(peakintensity):
        return True, None
    else:
        if debug:
            peaklist_length = [len(peaklist), len(peakintensity)]
            print(f"[debug] peaklist and intensity are not matched, {peaklist_length}")
            return False, peaklist_length
        return False, None

#READER-DIAG CP-D (20261020): instrument model and software version for the diagnostics (fail soft). The instrument
#name, serial number and file creator are not read.
def _raw_instrument(rawfile):
    info = {}
    for key, call in (("model", "GetInstModel"), ("software_version", "GetInstSoftwareVersion")):
        try:
            value = getattr(rawfile, call)()
            info[key] = value if isinstance(value, (str, int, float)) or value is None else str(value)
        except Exception:
            info[key] = None
    return info

def _backend_info():
    info = {"name": "pymsfilereader", "available": bool(pymsreader)}
    for pkg in ("pymsfilereader", "comtypes"):
        try:
            from importlib import metadata
            info[pkg] = metadata.version(pkg)
        except Exception:
            info[pkg] = None
    return info

def _status_row(ms_level, scan, scan_filter, peaklist, peakintensity, peaksource, exported, **fields):
    """Status-file row for one RAW scan. Exported peaks are always centroids (label data or the mass list)."""
    if not exported:
        arrays = "empty"
    elif len(peaklist) != len(peakintensity):
        arrays = "length_mismatch"
    else:
        arrays = "ok"
    return {"ms_level": ms_level, "scan": scan, "exported": exported, "arrays_status": arrays,
            "n_peaks": len(peaklist), "polarity": scan_polarity(scan_filter), "representation": "centroid",
            "peak_source": peaksource, **fields}

#usage: peak_extractor(pathofrawfile, raw_file_examination(rawfile))
#def peak_extractor(rawfileinput, metadata=None, filename=None, auto=True, debug=False):
def convert_raw_to_csv(rawfileinput, outdir=None, debug=False, charge_comb=False, dataset_charge=None, stem=None,
                       entry=None): #split peak_extractor to 2 functions
    """Convert one RAW file to the MS2/MS3 temp CSVs; also writes the status file and diagnostics files (CP-D).

    `stem` names every output (default: the input file stem, as before); `entry` is recorded in the diagnostics
    (the CLI passes "cli-convert"). Returns the legacy dict plus outcome, exit_code, warnings and diagnostics
    paths, or False when the file cannot be opened. An exception during conversion is recorded and re-raised.
    """
    #v0.59 CP-C.5: charge_comb / dataset_charge fill MS2 charges the file does not report (both off by default)
    dataset_charge = validate_dataset_charge(dataset_charge)
    if not outdir:
        outdir = os.path.dirname(rawfileinput) or os.getcwd()
    stem = stem or pathlib.Path(rawfileinput).stem
    diag = ReaderDiagnostics("msprawextractor", version, last_update, rawfileinput, stem, "raw", entry=entry,
                             settings={"charge_comb": bool(charge_comb), "dataset_charge": dataset_charge,
                                       "comb": msp_charge_comb.parameters() if charge_comb else None},
                             backend=_backend_info())
    try:
        result = _convert_raw_to_csv(rawfileinput, outdir, debug, charge_comb, dataset_charge, stem, diag)
    except BaseException as exc:   #v0.61 R2-F6: Ctrl-C / SystemExit finish the diagnostics too; re-raised unchanged
        diag.finish(outdir, exc)
        raise
    done = diag.finish(outdir)
    if isinstance(result, dict):
        result.update(done)
    return result

def _convert_raw_to_csv(rawfileinput, outdir, debug, charge_comb, dataset_charge, stem, diag):
    try:
        rawfile = MSFileReader(rawfileinput)
    except Exception as e:
        print(f"[Error] MSFileReader failed: {diag.redact(e)}")   #v0.61 R2-F3: redacted
        diag.fail(e)
        return False
    diag.stage = "iterate"   #v0.61 R2-F7: only MSFileReader(...) above counts as "could not be opened"
    diag.instrument = _raw_instrument(rawfile)
    # Default to the raw file’s folder, else CWD
    if not outdir:
        outdir = os.path.dirname(rawfileinput) or os.getcwd()
    
    scan_number = raw_file_examination(rawfile, auto = True, time = 0.0,debug=False)
    diag.declared = scan_number
    #init variables for storing parsed data
    ms1list, ms2list, ms3list, errlist , unmatchlist, logs = [], [], [], [], [], []
    #from v0.4 named Tuple of (peaklist,intensity) were changed to 2 separated lists. Include length validation 
    ms2data, ms3data = [], []  #b for MS2, c for MS3 < replacing
    fixms2header = ('entry_no','MS1scan_no', 'MS1_isolationmass', 'MS1_monoisolationmass','chargeState','protonatedmass',
           'MS2scan_no', 'peaklist', 'peakintensity')
    ms2data.append(fixms2header)
    fixms3header =  ('entry no','MS2scan_no', 'MS2_isolationmass', 'MS2_monoisomass', 'MS3scan_no', 'peaklist', 'peakintensity')
    ms3data.append(fixms3header)
    #MS2peaklist = namedtuple('MS2peaklist', ('dMass', 'dIntensity'))
    #MS3peaklist = namedtuple('MS3peaklist', ('dMass', 'dIntensity'))
    ms2count, ms3count = 1, 1 
    peak_summary = _new_peak_summary()
    charge_summary = _new_charge_summary(bool(charge_comb), dataset_charge)
    diag.charge = charge_summary
    ms1_cache = {}

    #development area and exception prevention
    if debug:
        print("[debug] need to create debug logs when debug is set to True.")
        print("[debug] Also need a readable log that records needed information and the operation can be reproduced when loading that record.")
        print("[debug add debug log and record log in debug block]")
        print("[debug]: start calculate time spent for converting data")
        perf_esti1 = time.time() #first time stamp for evaluating performance
    
    
    #perform extraction in defined spectra range
    for i in range(scan_number):
        #i starts from 0, using j for spectra count
        j = i+1
        ms_order = rawfile.GetMSOrderForScanNum(j)   #v0.60: read once per scan
        diag.step(j, ms_order)
        #do nothing for ms1 spectrum
        if ms_order == 1:
            if debug:
                ms1list.append(j)
                #trace back ms1 spectrum list
            pass
        #extract ms2 label data (processed by Thermo library)
        elif ms_order == 2:
            if debug:
                ms2list.append(j)
                #trace back ms2 spectrum list
            #Extract peaks (label data for FTMS scans, mass list otherwise)
            scan_filter = rawfile.GetFilterForScanNum(j) or ""
            peaklist, peakintensity, peaksource = read_scan_peaks(rawfile, j, scan_filter)
            if not _tally_peaks(peak_summary["ms2"], j, peaklist, peakintensity, peaksource):
                diag.add(_status_row(2, j, scan_filter, peaklist, peakintensity, peaksource, False))
                continue
            precursorinfo = rawfile.GetPrecursorInfoFromScanNum(j)
            #v0.59: charge in the agreed order; mass from the mono m/z, or the isolation m/z when the file has none
            chargestate, massmz, mzsource = _resolve_scan_charge(rawfile, j, precursorinfo, charge_comb,
                                                                 dataset_charge, ms1_cache, charge_summary)
            protonatedmass = calcproton(massmz, chargestate, polarity=scan_polarity(scan_filter)) #Get protonated mass using calcproton(isolatedmass, charge)
            scans = charge_summary["scans"]
            entry = scans[-1] if scans and scans[-1]["scan"] == j else {"charge_source": CHARGE_SOURCE_REPORTED}
            diag.add(_status_row(2, j, scan_filter, peaklist, peakintensity, peaksource, True,
                                 parent_scan=precursorinfo[3], precursor_mz=precursorinfo[0], charge=chargestate,
                                 mass=protonatedmass, mass_mz_source=mzsource,
                                 **{k: v for k, v in entry.items() if k not in ("scan", "charge", "mass_mz_source")}))

            #Arrange extracted data
            ms2spectrumdata = (
                ms2count,
                precursorinfo[3],#parentScanNo, MS1
                precursorinfo[0],#Isolation mass
                precursorinfo[1],#monoIsotopic Isolationmass
                chargestate,#chargeState
                protonatedmass, 
                j, #MS2scan no
                peaklist,
                peakintensity
                )
            ms2data.append(ms2spectrumdata)
            ms2count +=1
            #check if peaklist length unmatch
            peakcheck, failedlist = peaklist_validation(peaklist, peakintensity, debug = True)
            if not peakcheck: #if the peaklist and intensity aren't matched, not False = True
                reportunmatch = [j, failedlist]
                unmatchlist.append(reportunmatch)
            #ms2 extraction ends
        #extract ms3 label data as well (Note that position of data is not exactly the same) 
        #no validation of peaklists here since MS3 has really few peaks
        elif ms_order == 3:
            if debug:
                ms3list.append(j)
                #trace back ms3 spectrum list
            #v0.56: read this scan's own peaks (earlier versions wrote the previous MS2 scan's arrays here)
            scan_filter = rawfile.GetFilterForScanNum(j) or ""   #v0.60: read here for the polarity (was inside read_scan_peaks)
            peaklist, peakintensity, peaksource = read_scan_peaks(rawfile, j, scan_filter)
            if not _tally_peaks(peak_summary["ms3"], j, peaklist, peakintensity, peaksource):
                diag.add(_status_row(3, j, scan_filter, peaklist, peakintensity, peaksource, False))
                continue
            ms3spectrumdata = (
                ms3count,
                rawfile.GetPrecursorInfoFromScanNum(j)[3],#parent MS2 scan no
                rawfile.GetPrecursorInfoFromScanNum(j)[0],#Isolation mass
                rawfile.GetPrecursorInfoFromScanNum(j)[1],#monoisotopic isolationmass
                j, #MS3scan no
                peaklist,
                peakintensity 
                )
            ms3data.append(ms3spectrumdata)
            ms3count +=1
            diag.add(_status_row(3, j, scan_filter, peaklist, peakintensity, peaksource, True,
                                 parent_scan=ms3spectrumdata[1], precursor_mz=ms3spectrumdata[2]))
        else:
            errlist.append(j)        
        i+=1
    print("Extraction finished.")
    for level in ("ms2", "ms3"):
        block = peak_summary[level]
        if block["by_source"]:
            print(f"[peaks] {level.upper()} peak sources: {block['by_source']}; exported {block['exported']}, skipped without peaks: {len(block['empty'])}")
    print(f"[charge] MS2 charge sources: {charge_summary['sources']}; comb: {charge_summary['comb_status']}; "
          f"mass m/z: {charge_summary['mass_mz_source']}")
    rawfile.Close()
    #close reading raw file to release memory after extracting labeled data
    
    #logs are version information and exceptions returned from the program
    #debug_logs are only appended when we have debug=True
    if unmatchlist is not []:
        logs.append("Unmatched peaklists are recorded as [[spectrum no], [peaklist entries, peak intensity entries]]")
        logs.append(unmatchlist)

    if debug:
        debug_headers = "Debug enabled. Debug log information: [MS1 list, MS2 list, MS3 list, exception list]"
        debug_logs = [debug_headers, ms1list, ms2list, ms3list, errlist]
    else:
        debug_logs = None

    #estimate time spent for running the extraction
    if debug:
        extractiontime = time.time() - perf_esti1  #timestamp1: time spent for raw file extraction
        print("[debug]Raw file has been closed.")
        print("[debug]Preparing information for writing...")
        print(f"[debug]Time spent for extraction: {int(extractiontime)} seconds.")
    #call function fillexpinfo to complete metadata
    #if metadata is None or filename is None:
    #    raise ValueError("Metadata and filename must be provided by the GUI.")
    #ask for output format (csv as default or arff)
    #temp refrain the askformat section
    askformat = "c"
    #askformat = input("Please decide format for saving extracted file. (A/a/arff) for arff format, (C/c/csv) or leave it blank will generate csv file.")
    #make filename
    #os.path.splitext(rawfile)[0]
    #filename = str(generatefilename(rawfileinput)) #mind that rawfile is MSLoader.raw, not path like object

    #20250908: observed ms2, ms3 tmp persist issue, I/O performed twice, consider if we can simplify it to one

    filename = os.path.splitext(os.path.basename(rawfileinput))[0]
    ms2output = "ms2tmp_" + filename
    ms3output = "ms3tmp_" + filename

    #check if export in arff (future)
    if exportasarff(askformat):
        save_to_arff()
        print("Will be supported soon. Please use csv for now instead.")
    else:
        #20250916 try to fix duplicate tmp files issue
    # New: extractor only writes TEMP files; GUI will finalize/rename later.
        rawstem = stem   #v0.60 CP-D: the input file stem unless the caller names the outputs (CLI alias)
        #v0.53
        ms2tmp_path = os.path.join(outdir, f"ms2tmp_{rawstem}.csv")
        ms3tmp_path = os.path.join(outdir, f"ms3tmp_{rawstem}.csv")
        # save_to_csv still takes a *stem*; pass the stem without .csv but with full dir
        diag.stage = "write"
        save_to_csv(ms2tmp_path[:-4], ms2data, debug=debug)
        diag.wrote("ms2", ms2tmp_path)   #v0.61 R2-F2: named in the diagnostics once written
        save_to_csv(ms3tmp_path[:-4], ms3data, debug=debug)
        diag.wrote("ms3", ms3tmp_path)
        return {"ms2tmp": ms2tmp_path, "ms3tmp": ms3tmp_path, "summary": {"peaks": peak_summary, "charge": charge_summary}}
    #export log files...
    versioninfo = ["extractor info", version, last_update]
    logs.append(versioninfo)
    print("Conversion finished from mspextractor.")
    return True


def finalize_extraction(rawfile, metadata, filename, debug=False):
    print(f"[Finalizing extraction] For: {filename}")
    print(f"[Metadata provided] Sample name: {metadata.get('Experiment Title', 'Unknown')}")

    # Optional: log the fact that finalization completed
    if debug:
        print("[Debug] finalize_extraction called successfully")


'''
TODO list in future:
0. support arff conversion after validating everything works and documentations are done (URGENT)
1. read metadata from previous studies (another file) 
    and let it autofill when editing metadata <- this file
2. support other formats (another file) and convert to same format (maybe another file too) - mzml prototype finished
3. convert above things into GUI-based stuff <- this file?
'''
