import glob
import pandas as pd
import os
import pathlib
import numpy as np
from pymsfilereader import MSFileReader
#import datetime
import time  #for calculating efficiency
from collections import namedtuple
rawfile = MSFileReader("zf_sPerMeNG_brain")
rawname = "zf_sPerMeNG_brain"
###to supress warning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
#####

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
#ms3
c = []
#headtitle = ('MS2scan no', 'Isolation mass', 'monoIsomass','chargeState','parentScanNo')
#b.append(headtitle)
header = ('entry no','MS1scan no', 'MS1Isolation mass', 'MS1monoIsomass','chargeState','in [H+]',
          'intensityN\/A','StructureN\/A', 'MS2 Scan no', 'peaklist')
b.append(header)
ms3header =  ('entry no','MS3scan no', 'MS2Isolation mass', 'MS2monoIsomass', 'MS2 Scan no', 'peaklist')
c.append(ms3header)
MS2peaklist = namedtuple('MS2peaklist', ('dMass', 'dIntensity'))
MS3peaklist = namedtuple('MS3peaklist', ('dMass', 'dIntensity'))
ms1no = 1
ms2no = 1

#prevent error
if scan_number > maxn:
    print('you\'re trying to get more than what the file has')
    scan_number = maxn
    print('set scan_number to', maxn)

#running search throughout the spectra you called
for i in range(scan_number):
    j = i+1
    if i > maxn:
        break
    elif rawfile.GetMSOrderForScanNum(j) == 1:
        ms1list[i] = j
        #shooud I keep this? the list is totally referred from MS2
    elif rawfile.GetMSOrderForScanNum(j) == 2:
        peaklist = MS2peaklist((rawfile.GetLabelData(j)[0][0]), (rawfile.GetLabelData(j)[0][1]))
        isolationmass = rawfile.GetPrecursorInfoFromScanNum(j)[1]
        chargestate = rawfile.GetPrecursorInfoFromScanNum(j)[2]
        inHmass = calcproton(isolationmass,chargestate)
        a = (ms1no,
             rawfile.GetPrecursorInfoFromScanNum(j)[3],#parentScanNo
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
        ms1no +=1
    elif rawfile.GetMSOrderForScanNum(j) == 3:
        ms3list[i] = j
        peaklist = MS3peaklist((rawfile.GetMassListFromScanNum(j)[0][0]), (rawfile.GetMassListFromScanNum(j)[0][1]))
        d = (ms2no,
             rawfile.GetPrecursorInfoFromScanNum(j)[3],#parentScanNo
             rawfile.GetPrecursorInfoFromScanNum(j)[0],#Isolation mass
             rawfile.GetPrecursorInfoFromScanNum(j)[1],#monoIsomass
             j, #MS2scan no
             peaklist #surely will have error
             )
        c.append(d)
        ms2no +=1
    else:
        errlist[i] = j        
    i+=1

rawfile.Close()

print('MS1 list is skipped,see next')
'''
for key,value in ms1list.items():
    print(value)
'''
print('MS2 detection')
#extract MS2 info into csv

print('MS3 detection')
#extract MS2 info into csv
timestamp = time.strftime("%Y%m%d-%H%M%S") #datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
save_ms2filename= input("Please enter the file name or leave it blank to generate a filename with datetime")
if len(save_ms2filename) ==0:
    fileext = rawname + ' MS2 summary from ' + str(scan_number) + " at " + timestamp + '.csv'
else:
    fileext = rawname + save_ms2filename + " at " + timestamp + ".csv"
print(f"Filename output:{fileext}")
with open(fileext, 'wt') as g:
    print("writing to csv...")
    for q in range(len(b)):
        print('\t'.join(map(str, (b[q]))), file = g)
print(b, 'is saved in ', fileext)
print("MS1 will be quiried from MS2")
timestamp = time.strftime("%Y%m%d-%H%M%S") 
save_ms3filename= input("Please enter the MS3 file name or leave it blank to generate a filename with datetime")
if len(save_ms3filename) ==0 :
    fileext = rawname + ' MS3 summary from ' + str(scan_number) + " at " + timestamp + '.csv'
else:
    fileext = rawname + save_ms3filename + " at " + timestamp + ".csv"
print(f"Filename output:{fileext}")
with open(fileext, 'wt') as g:
    print("writing to csv...")
    for q in range(len(c)):
        print('\t'.join(map(str, (c[q]))), file = g)
print(c, 'is saved in ', fileext)
print("MS3 from MS2")


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
'''
print('MS2 list')
for key,value in ms2list.items():
    print(value)
'''
#print('MS3 list')
#for key,value in ms3list.items():
#    print(value)
print('err list')
for key,value in errlist.items():
    print(value)

print("dump finished")
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
