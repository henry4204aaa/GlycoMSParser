import glob
import pandas as pd
import os
import pathlib
import numpy as np
from pymsfilereader import MSFileReader
from collections import namedtuple
#for calculating efficiency
import time  
###to supress warning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
#####
try:
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
except ImportError:
    raise ImportError('Please install tkinter')

### initialize

def openfile():
    print('Going to select a file')
    Tk().withdraw()
    filename = askopenfilename()
    print('Selected', filename)
    if filename == "":
        filename = os.path.abspath("G:\\zf_sPerMeNG_brain.raw")
        print("default", filename, "is used")
    rawfile = MSFileReader(filename)
    return rawfile

def listoutheaderstest(rawfile):
    print('test function listoutheaderstest')


def listoutheaders(rawfile):    
    print('this raw file has', rawfile.GetNumSpectra(), 'spectra')
    scan_number = int(input('enter spectrum no you want to summarize'))
    maxn = int(rawfile.GetNumSpectra())
    ms1list = {}
    ms1fromms2 = []
    ms2list = {}
    ms3list = {}
    errlist = {}
    b = []
    #headtitle = ('MS2scan no', 'Isolation mass', 'monoIsomass','chargeState','parentScanNo')
    #b.append(headtitle)
    header = ('entry no','MS1scan no', 'MS1Isolation mass', 'MS1monoIsomass','chargeState','in [H+]',
              'intensityN\/A','StructureN\/A', 'MS2 Scan no', 'peaklist')
    b.append(header)
    ms1no = 1
    MS2peaklist = namedtuple('MS2peaklist', ('dMass', 'dIntensity'))

    #running search throughout the spectra you called
    for i in range(scan_number):
        j = i+1
        #prevent error
        if scan_number > maxn:
            print('you\'re trying to get more than what the file has')
            scan_number = maxn
            print('set scan_number to,', maxn)
        elif i > maxn:
            break
        elif rawfile.GetMSOrderForScanNum(j) == 1:
            ms1list[i] = j
            #shooud I keep this? the list is totally referred from MS2
        elif rawfile.GetMSOrderForScanNum(j) == 2:
            peaklist = MS2peaklist((rawfile.GetLabelData(j)[0][0]), (rawfile.GetLabelData(j)[0][1]))
            a = (ms1no,
                 rawfile.GetPrecursorInfoFromScanNum(j)[3],#parentScanNo
                 rawfile.GetPrecursorInfoFromScanNum(j)[0],#Isolation mass
                 rawfile.GetPrecursorInfoFromScanNum(j)[1],#monoIsomass
                 rawfile.GetPrecursorInfoFromScanNum(j)[2],#chargeState
                 'in H', #calculate
                 'ext from peak list',
                 'structure na',
                 j, #MS2scan no
                 peaklist #surely will have error
                 )
            b.append(a)
            ms1no +=1
        elif rawfile.GetMSOrderForScanNum(j) == 3:
            ms3list[i] = j
        else:
            errlist[i] = j        
        i+=1
    rawfile.Close()

    print('MS1 list is skipped,see next')
    print('MS2')
    #extract MS2 info into csv
    fileext = 'the MS2 summary from' + str(scan_number) +  '.csv'
    print(b, 'is also saved in ', fileext)
    with open(fileext, 'wt') as g:
        for q in range(len(b)):
            print('\t'.join(map(str, (b[q]))), file = g)
    print("MS1 from MS2")

    print('MS1 list from MS2')
    '''
    print('MS2 list')
    for key,value in ms2list.items():
        print(value)
    '''
    print('MS3 list')
    for key,value in ms3list.items():
        print(value)
    print('err list')
    for key,value in errlist.items():
        print(value)

#native python pseudocode

#test OG in 100 (only to ms2)
#for i in range(spectrumNumber):
#    if MSLevel = 1:
#        ms1list = [rawfile function]
#        dict.append(mslist) (does it work?)
#    elif MSLevel = 2:
#        list =[rawfile function]
#        dictforMS2.append(list) (does it work?)
#    elif MSlevel = 3:
#          # here the functions for MS3 will be different, be aware of it
#        listms3 = [rawdile function]
#    else:
#        errlist.append[i]
# monitor memory use and time spend on the function
#write list/dicts to csv
#MS1
#for j in range(len(ms1list)):
#    write MS1 dicts to tmpMS1.csv
#for k in range(len(ms2list)):
#    write MS2 dicts to tmpMS2.csv
#for m in range(len(errlist)):
#    

