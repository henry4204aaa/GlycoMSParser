from collections import namedtuple
import time
import glob
import pandas as pd
import os
import pathlib

#copied from comparepeaklist.py
def calcproton(isolatedmass, charge):
    if abs(charge) == 1 :
        return isolatedmass
    elif charge >= 2:
        conv = (isolatedmass * charge - (charge - 1) * 1.00784)
        return conv
    elif charge == 0:
        return 0
    elif charge < -1:
        negconv = (isolatedmass * charge + (charge - 1) * 1.00784)
        return negconv
    

def peak_extractor(rawfile, auto = True, debug = False):

    print(f'This raw file has {rawfile.GetNumSpectra()} spectra')
    #define number of spectra to parse (TODO: make time exclusion convertible and need to be added to record file.)
    if not auto:
        scan_number = int(input('enter spectrum no you want to summarize'))
        maxn = int(rawfile.GetNumSpectra())
    else:
        scan_number = int(rawfile.GetNumSpectra())
        maxn = scan_number

    #init variables for storing parsed data
    ms1list, ms2list, ms3list, errlist = {}, {}, {}, {}
    ms1fromms2 = []
    b, c = [], []  #b for MS2, c for MS3
    header = ('entry no','MS1scan no', 'MS1Isolation mass', 'MS1monoIsomass','chargeState','in [H+]',
          'intensity','Structure', 'MS2 Scan no', 'peaklist')
    b.append(header)
    ms3header =  ('entry no','MS3scan no', 'MS2Isolation mass', 'MS2monoIsomass', 'MS2 Scan no', 'peaklist')
    c.append(ms3header)
    MS2peaklist = namedtuple('MS2peaklist', ('dMass', 'dIntensity'))
    MS3peaklist = namedtuple('MS3peaklist', ('dMass', 'dIntensity'))
    ms2count, ms3count = 1, 1, 

    #development area and exception prevention
    if debug:
        print("[debug] need to create debug logs when debug is set to True.")
        print("[debug]Also need a readable log that records needed information and the operation can be reproduced when loading that record.")
        print("[debug add debug log and record log in debug block]")
        print("[debug]: start calculate time spent for converting data")
        perf_esti1 = time.time()
    
    if scan_number > maxn:
        print(r"you're attempting to run a number more than the spectra this file has.")
        scan_number = maxn
        print('set scan_number to', maxn)
    #
    
    #perform file conversion in defined spectra range
    for i in range(scan_number):
        j = i+1
        if i > maxn:
            break
        elif rawfile.GetMSOrderForScanNum(j) == 1:
            #ms1list[i] = j
            pass
            #shooud I keep this? the list is totally referred from MS2
        elif rawfile.GetMSOrderForScanNum(j) == 2:
            peaklist = MS2peaklist((rawfile.GetLabelData(j)[0][0]), (rawfile.GetLabelData(j)[0][1]))
            isolationmass = rawfile.GetPrecursorInfoFromScanNum(j)[1]
            chargestate = rawfile.GetPrecursorInfoFromScanNum(j)[2]
            inHmass = calcproton(isolationmass,chargestate)
            a = (ms2count,
                rawfile.GetPrecursorInfoFromScanNum(j)[3],#parentScanNo, MS1
                rawfile.GetPrecursorInfoFromScanNum(j)[0],#Isolation mass
                isolationmass,#monoIsomass
                chargestate,#chargeState
                inHmass, #calculate calcproton(isolatedmass, charge)
                'ext from peak list',
                'structure na',
                j, #MS2scan no
                peaklist
                )
            b.append(a)
            ms2count +=1
        elif rawfile.GetMSOrderForScanNum(j) == 3:
            ms3list[i] = j
            peaklist = MS3peaklist((rawfile.GetMassListFromScanNum(j)[0][0]), (rawfile.GetMassListFromScanNum(j)[0][1]))
            d = (ms3count,
                rawfile.GetPrecursorInfoFromScanNum(j)[3],#parentScanNo
                rawfile.GetPrecursorInfoFromScanNum(j)[0],#Isolation mass
                rawfile.GetPrecursorInfoFromScanNum(j)[1],#monoIsomass
                j, #MS2scan no
                peaklist #surely will have error
                )
            c.append(d)
            ms3count +=1
        else:
            errlist[i] = j        
        i+=1
    print("Extraction finished.")
    rawfile.Close()

    if debug:
        extractiontime = time.time() - perf_esti1
        print("[debug]Raw file has been closed.")
        print("[debug]Preparing information for writing...")
        print(f"[debug]Time spent for extraction: {int(extractiontime)} seconds.")

    timestamp = time.strftime("%Y%m%d-%H%M%S") #datetime.datetime.now().strftime("%Y%m%d-%H%M%S")

    #export MS2 spectra file in csv format
    save_ms2filename= input("Please enter the file name or leave it blank to generate a filename with datetime")
    if len(save_ms2filename) == 0:
        fileext = rawname + ' MS2 summary from ' + str(scan_number) + " at " + timestamp + '.csv'
    else:
        fileext = rawname + save_ms2filename + " at " + timestamp + ".csv"
    print(f"Filename output:{fileext}")
    perf_esti2 = time.time()
    with open(fileext, 'wt') as g:
        print("writing to csv...")
        for q in range(len(b)):
            print('\t'.join(map(str, (b[q]))), file = g)
    if debug:
        ms2writetime = time.time() - perf_esti2
        print("[debug]Finished writing MS2 information")
        print(f"[debug]Time spent for extraction: {int(ms2writetime)} seconds.")

    #export MS3 spectra file in csv format
    timestamp = time.strftime("%Y%m%d-%H%M%S") 
    save_ms3filename= input("Please enter the MS3 file name or leave it blank to generate a filename with datetime")
    if len(save_ms3filename) ==0 :
        fileext = rawname + ' MS3 summary from ' + str(scan_number) + " at " + timestamp + '.csv'
    else:
        fileext = rawname + save_ms3filename + " at " + timestamp + ".csv"
    print(f"Filename output:{fileext}")
    perf_esti3 = time.time()
    with open(fileext, 'wt') as h:
        print("writing to csv...")
        for q in range(len(c)):
            print('\t'.join(map(str, (c[q]))), file = h)
    if debug:
        ms3writetime = time.time() - perf_esti3
        print("[debug]Finished writing MS3 information")
        print(f"[debug]Time spent for extraction: {int(ms3writetime)} seconds.")


    #export log files...
    if debug:
        print("Write logs")
        print(f"file name: {fileext}")
        print(f"datetime (latest): {timestamp}")
        print(f"file name: {rawname}")
        print("file path is not possible in current context.")
        print(f"scan number set in this experiment: {scan_number}")
        print(f"GlycoMSP main version: {GlycoMSP_version}")
        #info in logs: datetime, filename, file path, scan number, version of GlycoMSP and EACH components
    
    #report errors
    print('err list')
    for key,value in errlist.items():
        print(value)

    print("dump finished")




def peak_exthandler(rawfile, filetype=None, debug=False):
    if filetype == "raw":
        peak_extractor(rawfile, auto= True, debug=debug)#or debug=debug can work
    elif filetype == "mzML":
        print("mzML file will be supported in v1.1")
    else:
        print("unsupported format")


def runextractor(debug=False):
    if GlycoMSPinit:
        try:
            import glob
            import pandas as pd
            import os
            import pathlib
            import numpy as np
            import time  #for calculating efficiency
            #from collections import namedtuple
            import warnings
            warnings.simplefilter(action='ignore', category=FutureWarning) ###to supress warning
        except:
            print("Missing essential package for following analysis") #try to list out the package missing
        if rawfile:
            print("Start raw file extration process...")
            filetype="raw"
            peak_exthandler(rawfile, filetype, debug)

runextractor(debug=True)


'''
Start raw file extration process...
This raw file has 44218 spectra
[debug] need to create debug logs when debug is set to True.
[debug]Also need a readable log that records needed information and the operation can be reproduced when loading that record.
[debug add debug log and record log in debug block]
[debug]: start calculate time spent for converting data
Extraction finished.
[debug]Raw file has been closed.
[debug]Preparing information for writing...
[debug]Time spent for extraction: 288 seconds.
Please enter the file name or leave it blank to generate a filename with datetime
Filename output:G:/zf_sPerMeNG_brain.raw MS2 summary from 44218 at 20231009-203347.csv
writing to csv...
[debug]Finished writing MS2 information
[debug]Time spent for extraction: 4 seconds.
Please enter the MS3 file name or leave it blank to generate a filename with datetime
Filename output:G:/zf_sPerMeNG_brain.raw MS3 summary from 44218 at 20231009-222454.csv
writing to csv...
[debug]Finished writing MS3 information
[debug]Time spent for extraction: 0 seconds.
Write logs
file name: G:/zf_sPerMeNG_brain.raw MS3 summary from 44218 at 20231009-222454.csv
datetime (latest): 20231009-222454
file name: G:/zf_sPerMeNG_brain.raw
file path is not possible in current context.
scan number set in this experiment: 44218
GlycoMSP main version: 0.3
err list
dump finished

'''