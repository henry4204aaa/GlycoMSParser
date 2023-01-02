import glob
import pandas as pd
import os
import pathlib
import numpy as np
from pymsfilereader import MSFileReader
import time  #for calculating efficiency
from collections import namedtuple
rawfile = MSFileReader("zf_sPerMeOG_intestine")

###to supress warning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
#####

### initialize

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



def extractionlist():
    #see the implemention in another file
    pass

#running search throughout the spectra you called
for i in range(scan_number):
    j = i+1
    #prevent error
    if scan_number > maxn:
        print('you\'re trying to get more than what the file has')
        scan_number = maxn
        print('set scan_number to,' maxn)
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
    ############## useless
    #elif rawfile.GetMSOrderForScanNum(j) == 1:
    #    ms1list[i] = j

    #headtitle = ('MS2scan no', 'Isolation mass', 'monoIsomass','chargeState','parentScanNo')
    ''' block A
    elif rawfile.GetMSOrderForScanNum(j) == 2:
        a = ( j,rawfile.GetPrecursorInfoFromScanNum(j)[0],
              rawfile.GetPrecursorInfoFromScanNum(j)[1],
              rawfile.GetPrecursorInfoFromScanNum(j)[2],
              rawfile.GetPrecursorInfoFromScanNum(j)[3])
        b.append(a)
    '''
    ###new
rawfile.Close()

print('MS1 list is skipped,see next')
'''
for key,value in ms1list.items():
    print(value)
'''
print('MS2')
#extract MS2 info into csv
fileext = 'the MS2 summary from' + str(scan_number) +  '.csv'
print(b, 'is also saved in ', fileext)
with open(fileext, 'wt') as g:
    for q in range(len(b)):
        print('\t'.join(map(str, (b[q]))), file = g)
print("MS1 from MS2")

#exclude the first items "parentScanno"
''' this works for blockA 
for r in range(1, len(b)):
    ms1fromms2.append(b[r][-1])
print(ms1fromms2)
'''
###this work
'''
for q in range(len(b)):
    print(q , 'is', b[q])
    for s in range(len(b[q])):
        print('sub', (b[q])[s])
''' 

'''


'''
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




#below using pandas runs really slow, consider using python itself

'''
headertmp = pd.DataFrame([])
headertmp1 = pd.DataFrame([])
perf1 = 0
perf2 = 0
#headertmp1 = pd.Series([])
tmp = []
dfheader = pd.Series([],dtype=float)
for i in range(scan_number):
    j = i+1
    if scan_number > maxn:
        scan_number = maxn
    elif i > maxn:
        print("too large")
        break
    else:
        #print(i, j)
        #for MS1
        if rawfile.GetMSOrderForScanNum(j) == 1:
            continue
            #print("oh it's MS1")
        #for MS2
        elif rawfile.GetMSOrderForScanNum(j) == 2:
            dfheader = pd.Series(data=(j,
                          rawfile.GetPrecursorInfoFromScanNum(j)[0],
                          rawfile.GetPrecursorInfoFromScanNum(j)[1],
                          rawfile.GetPrecursorInfoFromScanNum(j)[2],
                          rawfile.GetPrecursorInfoFromScanNum(j)[3]),
                          index=["MS2scan no","Isolation mass","monoIsomass","chargeState","parentScanNo"]) 
            #warn but work
            start = time.time()
            headertmp = headertmp.append(dfheader, ignore_index=True)
            end = time.time()
            a = end -start
            perf1 += a
            #work but wrong format
            #tmp.append(dfheader)
            start1 = time.time()
            headertmp1 = pd.concat([headertmp1,dfheader], axis=1, join='outer', ignore_index=True)
            end1 = time.time()
            b = end1 -start1
            perf2 += b
        else:
            continue
            #print('lmao')
    i+=1
#headertmp1 = pd.concat(tmp, axis=0, join='outer', )
c = time.time()
headertmp2 = headertmp1.T
d = time.time()
e = d-c
f = perf2 + e
print('time for running on append', perf1)
print('time for running on concat', f, 'while concat spend', perf2, 'sec, transpose spend', e)
#FutureWarning: The frame.append method is deprecated and will
#be removed from pandas in a future version. Use pandas.concat instead.
'''
