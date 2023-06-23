import csv
import math
from decimal import Decimal,getcontext, ROUND_HALF_UP
from collections import namedtuple
import re
import pandas as pd
import comparepeaklist
import numpy as np

list1 = (344.1698)
roundedlist = []
roundedlist2 = []
tmplist = []
MS2peaklist = namedtuple('MS2peaklist', ('dMass', 'dIntensity'))



def setcomparelist():
    comparelist = []
    while True:
        a = input("enter fragments for searching")
        try:
            comparelist.append(float(a))
        except:
            print("You entered one isn't a number")
            break
    return comparelist
def setcomparelistfromitem():
    comparelist = []
    #pseudocode
    #list out a dict (fragment name, mass) and allow
    #adding multiple entries at one time
    #delete added fragment after adding (allow adding until user exit)
    
#return the tuple of both peak mass (rounded 0.01) and float intensity
# as a dict or (m, i) or in separated list (better for future design but unreadable by manwork
# get the whole list info and add simple annotation (referring previous research results)
#question: cant recognize correctly. See glypick results and extract certain spectrum for testing
def readms2csvpdreg():
    df = pd.read_csv('zf_sPerMeNG_intestine MS2 summary from 48127 at 20230611-182546 20230611-194107 predictedcomp.csv', sep='\t') #'debug1.csv the MS2 summary from44218.csv'
    #remove unmaned index
    df.drop(columns=df.columns[0], axis=1, inplace=True) 
    print(df.dtypes)
    hitspectrumlist = []
    err = 0
    notfound = 0
    hits = 0
    ms2scanno= 0
    subdf = []
    tmphitreturnlist = []
    ppm = comparepeaklist.ppmcalculation()
    cc = []
    ncc = []
    tmpdf = pd.DataFrame()
    for i in range(len(df)):
        #initialize the placeholder for input peak information (to be parsed)
        inputlist = []
        inputintensitylist = []
        #parse peaklist
        #load no. i row's peaklist of the dataframe
        tmp = df.loc[i, "peaklist"]
        #the peak intensity and peak mass are included as well.
        inputlist1 = tmp.split('(')[2].split(')')[0].split(',') #peakmass
        d = tmp.split('(')[3].split(')')[0].split(',') #peakintensity
        #print('inputlist =',type(inputlist), inputlist1, 'intensity=',type(d), d)
        #print('going to modify the number')
        for j in range(len(inputlist1)):
            #should be 0.0001 or adjustable in the future
            k = Decimal((inputlist1[j])).quantize(Decimal("0.001"),rounding=ROUND_HALF_UP) 
            k = float(k)
            l = Decimal((d[j])).quantize(Decimal("0.001"),rounding=ROUND_HALF_UP)
            l = float(l)
            inputlist.append(k)
            inputintensitylist.append(l)
        #print('modified list=', inputlist, inputintensitylist)
        #loss one OHCH3= 32.0419
        #call user input and if no vaild input, set it to default list
        comparelist = [246.1336, 260.1492, 303.1438, 335.1700 ,344.1704, 374.1810 ,376.1966, 406.2072, 432.2071, 450.2334, 464.2490, 539.2698, 580.2964, 621.3229, 638.3382, 793.3808, 825.4227, 842.4380]
        #comparelist = []#setcomparelist(), should be included in somewhere of argument of this function
        if comparelist is None:
            #here is pre-defined list [(Neu5Ac-OMe B ion, 344.17), (Neu5Gc-OMe B ion, 374.18), (Neu5AcB ion, 376.20), (Neu5Gc B ion, 406.21), (LacNAc-OMe B ion, 432.22),
            #(LacNAc-BY ion, 450.23), (LacNAc-B ion, 464.25), (Neu5AcHex B ion, 580.30), (SiaLacNAc-OMe B ion, 793.38), (SiaLacNAc B ion, 825.42)]
            #comparelist = [344.17, 374.18, 376.20, 406.21, 432.22, 450.23, 464.25, 580.30, 793.38, 825.42]
            #comparelist = [344.1547, 374.1653, 376.1966, 406.2072, 432.2071, 450.2334, 464.2490, 580.2964, 793.3808, 825.4227]
            comparelist = [246.1336, 260.1492, 303.1438, 335.1700 ,344.1704, 374.1810 ,376.1966, 406.2072, 432.2071, 450.2334, 464.2490, 539.2698, 580.2964, 621.3229, 638.3382, 793.3808, 825.4227, 842.4380]
        #print('compare list=', comparelist)
        result = comparepeaklist.comparepeaklistppm(inputlist, comparelist, ppm, inputintensitylist)
        #print(result)
        #######pseudocode
        # comparelist=input("fragments for searching in spectra")
        # check if all input is float type
        # if the user input is invalid, return default list or trap in a editing loop done
        if result[0] > 0:
            df.loc[i, "in [H+]"] = comparepeaklist.calcproton(df.loc[i, "MS1Isolation mass"],df.loc[i, "chargeState"])
            ms2scan = df.loc[i, "MS2 Scan no"]
            #hitspectrumlist.append([ms2scan, result[0], result[1], result[3], result[4]])
            hits+=1
            #before append calculate the charge conversion
            protonateass = comparepeaklist.calcproton(df.loc[i, "MS1Isolation mass"], df.loc[i, "chargeState"])
            #I think I should add them back at thise stage
            #And theoretical mass
            #And prediction of composition if prediction flag is set to True)
            addinfo = pd.Series([result[1], result[3], result[4]], index=['hitpeaklist', 'selectedpeakintensity','normalizedselectedin']) #, result[5], 'maxintensity'
            tmpdf = pd.concat([df.loc[i], addinfo], axis = 0)
            subdf.append(tmpdf)
            #debug1 = int(df.loc[i, "MS2 Scan no"])
            #if debug1 == int(14857):    #test why 376 and the relative intensity value is wrong
            #    print(type(addinfo), type(subdf))
            #    print(addinfo, subdf)
            #    break
        elif result[0] == 0:
            #print('no hits')
            notfound+=1
        else:
            print('unexpected error')
            err+=1
            #do nothing
        #need to catch the returned value and decide if this spectra is ok
    print('hits', hits, 'notfound', notfound, 'err', err) #'hits', hitspectrumlist)
    w = pd.DataFrame()
    #subdf.to_csv("ms2filteredtotallistwithHtest.csv", index=False)
    with open('ms2filteredtotallistwithH_NGintestine20230612_50ppm.csv', 'w', newline='') as zz:
        roww = csv.writer(zz)
        roww.writerow(["entry no", "MS1scan no", "MS1Isolation mass", "MS1monoIsomass", "chargeState", "in [H+]", "intensity", "Structure", "MS2 Scan no", "peaklist", "hitpeaklist", "selectedpeakintensity", "normalizedselectedintensity"])
        roww.writerows(subdf)
    print("Finished exporting")
    #return hitspectrumlist
    #with open('ms2searchsummary.csv', 'w', newline='') as z:
    #    roww = csv.writer(z)
    #    roww.writerows(hitspectrumlist)
    '''
    for q in range(len(hitspectrumlist)):
        print('spectrum no', int((hitspectrumlist[q][0])), 'list',
                (hitspectrumlist[q][1]))
    '''
    


def readms2csvpd():
    df = pd.read_csv('the MS2 summary from350.csv', sep='\t')
    print(df.dtypes)
    for i in range(len(df)):
        print(df.loc[i, "entry no"], df.loc[i, "peaklist"])
        #parse peaklist
        tmp = df.loc[i, "peaklist"]
        tmp = tmp.replace('(', '!').replace(')', '!').split('!')
        #drop named tuple value
        print(tmp[0],'\<zero \> first', tmp[1])
        #tmp = tmp[2]
        print(type(tmp), tmp)
        peaklist = tmp[2:-2]
        print(peaklist)
    #df1 = df["peaklist"]
    print(df)
