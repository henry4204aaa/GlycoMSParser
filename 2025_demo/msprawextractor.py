######version info#####

version= "0.53"
last_update = 20250927
#import msprawextractor   -> msprawextractor.version 
import versioninfo
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
#importing modules
#should be loaded somewhere for general data processing
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

#suppress warning during data hanlding process
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

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


    #save the metadata#
#def savemetadata(metadata,savename):
#    savename = f"{savename}.json"
#    with open(savename, 'w') as f:
#        json.dump(metadata, f, indent=4)
#    return savename
def savemetadata(metadata, filepath):
    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=4)
    return filepath


    #save the log file#
#def save_log(savename, log_entries,debug_logs = None):
#    logname = f"{savename}.log"
#    with open(logname, 'w') as f:
#        for entry in log_entries:
#            f.write(f"{entry}\n")
#        if debug_logs is not None:
#            f.write(f"-----\ndebug enabled-----")
#            for logs in debug_logs:
#                f.write(f"{logs}\n")
#            f.write(f"-----\ndebug log ends-----")
#    return logname
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

"""
def save_to_csv(csvname,data,debug =False):
    #export spectra file in csv format
    exportcsv =  f"{csvname}.csv"
    perf_esti2 = time.time()
    with open(exportcsv, 'wt') as c:
        print("writing to csv...")
        for d in range(len(data)):
            print('\t'.join(map(str, (data[d]))), file = c)
    if debug:
        csvwritetime = time.time() - perf_esti2
        print(f"[debug]Finished exported {exportcsv} from raw data")
        print(f"[debug]Time spent for extraction: {float(csvwritetime)} seconds.")
    return exportcsv
"""
#work work 
def save_to_arff():#(arffname, data, headers, attributes, debug=False):
    print("Waiting to be built...")
    #supposed to be added in 0.4 after csv files are validated okay.
#work work 


#calculates and returns protonated mass for precursor in each spetrum
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
    

#peak extractor parts, need modification

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

#usage: peak_extractor(pathofrawfile, raw_file_examination(rawfile))
#def peak_extractor(rawfileinput, auto = True, debug = False): #old version
#def peak_extractor(rawfileinput, metadata=None, filename=None, auto=True, debug=False):
def convert_raw_to_csv(rawfileinput, outdir=None, debug=False): #split peak_extractor to 2 functions
    try:
        rawfile = MSFileReader(rawfileinput)
    except Exception as e:
        print(f"[Error] MSFileReader failed: {e}")
        return False
    # Default to the raw file’s folder, else CWD
    if not outdir:
        outdir = os.path.dirname(rawfileinput) or os.getcwd()
    
    scan_number = raw_file_examination(rawfile, auto = True, time = 0.0,debug=False)
    #init variables for storing parsed data
    ms1list, ms2list, ms3list, errlist , unmatchlist, logs = [], [], [], [], [], []
    #from v0.4 named Tuple of (peaklist,intensity) were changed to 2 separated lists. Include length validation 
    ms2data, ms3data = [], []  #b for MS2, c for MS3 < replacing
    fixms2header = ('entry_no','MS1scan_no', 'MS1_isolationmass', 'MS1_monoisolationmass','chargeState','protonatedmass',
           'MS2scan_no', 'peaklist', 'peakintensity')
    ms2data.append(fixms2header)
    fixms3header =  ('entry no','MS3scan_no', 'MS2_isolationmass', 'MS2_monoisomass', 'MS2scan_no', 'peaklist', 'peakintensity')
    ms3data.append(fixms3header)
    #MS2peaklist = namedtuple('MS2peaklist', ('dMass', 'dIntensity'))
    #MS3peaklist = namedtuple('MS3peaklist', ('dMass', 'dIntensity'))
    ms2count, ms3count = 1, 1 

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
        #do nothing for ms1 spectrum
        if rawfile.GetMSOrderForScanNum(j) == 1:
            if debug:
                ms1list.append(j)
                #trace back ms1 spectrum list
            pass
        #extract ms2 label data (processed by Thermo library)
        elif rawfile.GetMSOrderForScanNum(j) == 2:
            if debug:
                ms2list.append(j)
                #trace back ms2 spectrum list
            #GetLabelData(spectrum number) will return plenty of data
            #GetLabelData(n)[0] will return lists (need validation) [peaklist], [peak intensity]
            #peaklist = MS2peaklist((rawfile.GetLabelData(j)[0][0]), (rawfile.GetLabelData(j)[0][1]))  @oldmethod
            #Extract labeled data
            peaklist = rawfile.GetLabelData(j)[0][0]
            peakintensity = rawfile.GetLabelData(j)[0][1]
            isolationmass = rawfile.GetPrecursorInfoFromScanNum(j)[1]
            chargestate = rawfile.GetPrecursorInfoFromScanNum(j)[2]
            protonatedmass = calcproton(isolationmass,chargestate) #Get protonated mass using calcproton(isolatedmass, charge)

            #Arrange extracted data
            ms2spectrumdata = (
                ms2count,
                rawfile.GetPrecursorInfoFromScanNum(j)[3],#parentScanNo, MS1
                rawfile.GetPrecursorInfoFromScanNum(j)[0],#Isolation mass
                isolationmass,#monoIsotopic Isolationmass
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
        elif rawfile.GetMSOrderForScanNum(j) == 3:
            if debug:
                ms3list.append(j)
                #trace back ms3 spectrum list
            #peaklist = MS3peaklist((rawfile.GetMassListFromScanNum(j)[0][0]), (rawfile.GetMassListFromScanNum(j)[0][1]))
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
        else:
            errlist.append(j)        
        i+=1
    print("Extraction finished.")
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
        rawstem = pathlib.Path(rawfileinput).stem
        #v0.53
        ms2tmp_path = os.path.join(outdir, f"ms2tmp_{rawstem}.csv")
        ms3tmp_path = os.path.join(outdir, f"ms3tmp_{rawstem}.csv")
        # save_to_csv still takes a *stem*; pass the stem without .csv but with full dir
        save_to_csv(ms2tmp_path[:-4], ms2data, debug=debug)
        save_to_csv(ms3tmp_path[:-4], ms3data, debug=debug)
        return {"ms2tmp": ms2tmp_path, "ms3tmp": ms3tmp_path}
        #v0.52
        #ms2tmp = f"ms2tmp_{rawstem}.csv"
        #ms3tmp = f"ms3tmp_{rawstem}.csv"
        #save_to_csv(ms2tmp.replace(".csv",""), ms2data, debug=debug)
        #save_to_csv(ms3tmp.replace(".csv",""), ms3data, debug=debug)

        # tmp_ms2 = Path(outdir) / f"ms2tmp_{rawstem}.csv"
        #tmp_ms3 = Path(outdir) / f"ms3tmp_{rawstem}.csv"

        # Ensure parent exists
        #tmp_ms2.parent.mkdir(parents=True, exist_ok=True)
        #tmp_ms3.parent.mkdir(parents=True, exist_ok=True)

        # Write temp CSVs (no finals here)
        #save_to_csv(str(tmp_ms2).replace(".csv",""), ms2data, debug=debug)
        #save_to_csv(str(tmp_ms3).replace(".csv",""), ms3data, debug=debug)

        #if debug:
        #    print(f"[extractor] temp CSVs written: {tmp_ms2.name}, {tmp_ms3.name}")

        # Return the TEMP paths for the GUI to finalize.
        # Return filename of the temp files
        #return os.path.abspath(ms2tmp), os.path.abspath(ms3tmp)
        #export to csv by default
        #ms2done = save_to_csv(ms2output, ms2data, debug=debug)
        #print(f"MS2 export {ms2done} has finished")
        #ms3done = save_to_csv(ms3output, ms3data, debug=debug)
        #print(f"MS3 export {ms3done} has finished")
    #export log files...
    versioninfo = ["extractor info", version, last_update]
    logs.append(versioninfo)
    print("Conversion finished from mspextractor.")
    return True
    #save metada to json /added in 20240731
    #savemetadata(metadata,filename)

    #save logs and add debug logs only if debug=True
    #if debug:
    #    save_log(filename, logs, debug_logs = debug_logs)
    #    print(f"[debug] Export debug logs...")
    #    #info in logs: datetime, filename, file path, scan number, version of GlycoMSP and EACH components
    #else:
    #    save_log(filename, logs)
    #print("Extraction finished")
    #print("from version 0.4 the metadata part is separated for GUI UI/UX flow")

    ############

def finalize_extraction(rawfile, metadata, filename, debug=False):
    print(f"[Finalizing extraction] For: {filename}")
    print(f"[Metadata provided] Sample name: {metadata.get('Experiment Title', 'Unknown')}")

    # Optional: log the fact that finalization completed
    if debug:
        print("[Debug] finalize_extraction called successfully")
    #try:
    #        # Just finalize logs and metadata
    #        jsonfile = savemetadata(metadata, filename)
    #        print(f"Metadata saved to: {jsonfile}")###

    #        logs = [f"File processed: {rawfile}", f"Saved metadata: {jsonfile}"]
    #        debug_logs = ["Debug mode ON"] if debug else []
    #        logname = save_log(filename, logs, debug_logs)
    #        print(f"Log saved to: {logname}")

    #except Exception as e:
    #    print(f"[ERROR] Finalization failed: {e}")
#read file
#rawfileinput = mspfileloader.fileloader('raw')  #filename.raw (str)
#set scan number outside
#run extractor function
#peak_extractor(rawfileinput, auto = True, debug = True) #debug = True for testing
##running code, save metadata and logs
#filename = generatefilename(rawfilenameraw)
#metadata = fillexpinfo(rawfilenameraw, loadfile=None, debug=False)
#metajson = savemetadata(metadata,filename)

#from loader v0.51 call from loader
'''
TODO list in future:
0. support arff conversion after validating everything works and documentations are done (URGENT)
1. read metadata from previous studies (another file) 
    and let it autofill when editing metadata <- this file
2. support other formats (another file) and convert to same format (maybe another file too)
3. convert above things into GUI-based stuff <- this file?
'''

'''
debug records:
why I can't define scan number?
'''
