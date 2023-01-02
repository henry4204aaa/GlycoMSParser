import csv
import math
from collections import namedtuple
import re
import pandas as pd
import numpy as np
import ast
from decimal import *

def parsecompositiontomatrix0():
    #this function is created for converting the structure dictionary
    #into independent column and fill missing value with 0
    #complete = "Annotated_20221130.csv"
    complete = "revised_zfNGbrain_20221129_1219_raw.csv"
    liteversion = "Annotatedlite20221130.csv"
    dfa = pd.read_csv(complete)#, sep='\t')
    dfb = pd.read_csv(liteversion)
    dfa.sort_index(inplace=True)
    dfb.sort_index(inplace=True)
    tempdfsel = input(r"enter 'lite' for dfb while others will be dfa")
    if tempdfsel == "lite":
        tempdf = dfb
    else:
        tempdf = dfa
    #init the empty list
    tmpopt = []
    tmpseries = ()
    tmp = ()
    for i in range(len(tempdf)):
        #use .copy to avoid chained settingwithcopywarning
        tmp = tempdf.iloc[i].copy()
        dictparser = tmp['Structure dict']
        parsed = ast.literal_eval(dictparser)
        tmpkeys = list(parsed.keys())
        tmpvalues = list(parsed.values())
        if i == 1:
            print("Hello this is line 1, and we'll get dict keys")
            print(type(parsed), parsed)
            #Store to list to avoid not subscriptable error
            print(dictparser)
            print(tmpkeys)
        #print(tmpvalues)
        for j in range(len(tmpkeys)):
            tmp[tmpkeys[j]] = tmpvalues[j]
        tmpopt.append(tmp)
    tmpopt = pd.DataFrame(tmpopt)
    tmpopt = tmpopt.reset_index(drop=True)
    tmpopt = tmpopt.fillna("0")
    file=input("Enter the file name")
    file = file + ".csv"
    tmpopt.to_csv(file)
    print("Finish exporting", file)
    
def parseintensitytomatrix0():
    #this function is created for converting the input and relative intensity
    #into independent column and fill missing value with 0
    #filename = "liteforPCA20221213v2.csv"
    filename = "revised_zfNGbrain_20221129_1219_proc1.csv"
    df = pd.read_csv(filename)
    df.sort_index(inplace=True)
    #init the empty list
    tmpopt = []
    tmpseries = ()
    tmp = ()
    for i in range(len(df)):
        #use .copy to avoid chained settingwithcopywarning
        tmp = df.iloc[i].copy()
        peaklistparser = tmp['hitpeaklist']
        intensityparser = tmp['normalizedselectedintensity']
        parsedpeaklistheader = list(ast.literal_eval(peaklistparser))
        parsedintensity = list(ast.literal_eval(intensityparser))
        #print("peaklist", type(parsedpeaklistheader), parsedpeaklistheader)
        #print("intensity", type(parsedintensity), parsedintensity)
        for j in range(len(parsedpeaklistheader)):
            #print(parsedpeaklistheader[j])
            #print(parsedintensity[j])
            tmp[parsedpeaklistheader[j]] = parsedintensity[j]
        tmpopt.append(tmp)
    tmpopt = pd.DataFrame(tmpopt)
    tmpopt = tmpopt.reset_index(drop=True)
    tmpopt = tmpopt.fillna(0)
    #need sum tho those cells are filled with 0 so it should be fine
    tmpopt = tmpopt.groupby(level=0, axis=1).sum()
    file=input("Enter the file name")
    file = file + ".csv"
    tmpopt.to_csv(file)
    print("Finish exporting", file)

def fixhitlistbugmethodA():
    #oroginal hitlist records original instrument readout
    #resulting in small mass difference making columns not able to merge
    #any possible functions to merge "similar column names?"
    #consider adding this parallel with the bug fixes
        #this function is created for converting the input and relative intensity
    #into independent column and fill missing value with 0
    filename = "revised_zfNGbrain_20221129_1219_proc2.csv"
    #filename = "liteforPCA20221213v2.csv"
    df = pd.read_csv(filename)
    df.sort_index(inplace=True)
    #init the empty list
    tmpopt = []
    tmpseries = ()
    tmp = ()
    for i in range(len(df)):
        #use .copy to avoid chained settingwithcopywarning
        tmp = df.iloc[i].copy()
        peaklistparser = tmp['hitpeaklist']
        intensityparser = tmp['normalizedselectedintensity']
        parsedpeaklistheader = list(ast.literal_eval(peaklistparser))
        parsedintensity = list(ast.literal_eval(intensityparser))
        #print("peaklist", type(parsedpeaklistheader), parsedpeaklistheader)
        #print("intensity", type(parsedintensity), parsedintensity)
        for j in range(len(parsedpeaklistheader)):
            #print(parsedpeaklistheader[j])
            #print(parsedintensity[j])
            tmp[parsedpeaklistheader[j]] = parsedintensity[j]
        tmpopt.append(tmp)
    tmpopt = pd.DataFrame(tmpopt)
    tmpopt = tmpopt.reset_index(drop=True)
    tmpopt = tmpopt.fillna(0)
    a = list(tmpopt.columns)
    b = []
    for atoint in a:
        try:
            atoint = int(atoint)
        except:
            atoint = str(atoint)
        b.append(atoint)
    #comparelist = [344.1547, 374.1653, 376.1966, 406.2072,
    #432.2071, 450.2334, 464.2490, 580.2964, 793.3808, 825.4227]
    c = []
    for bconv in b:
        if bconv == int(344):
            bconv = 344.1547
        elif bconv == int(374):
            bconv = 374.1653
        elif bconv == int(376):
            bconv = 376.1966
        elif bconv == int(406):
            bconv = 406.2072
        elif bconv == int(432):
            bconv = 432.2071
        elif bconv == int(450):
            bconv = 450.2334
        elif bconv == int(464):
            bconv = 464.2490
        elif bconv == int(580):
            bconv = 580.2964
        elif bconv == int(793):
            bconv = 793.3808
        elif bconv == int(825):
            bconv = 825.4227
        c.append(bconv)
    print('modified?', c)
    #replace the columns
    tmpopt.set_axis(c, axis=1, inplace=True)#Warning
    #FutureWarning: DataFrame.set_axis 'inplace' keyword is deprecated and will be removed in a future version. Use `obj = obj.set_axis(..., copy=False)` instead
    #tmpopt.set_axis(c, axis=1, copy=False)
    print('modified index', tmpopt.columns)
    tmpopt = tmpopt.groupby(level=0, axis=1).sum()
    tmpopt = tmpopt.reset_index(drop=True)
    cols = tmpopt.columns.tolist()
    cols = ['MS1scan no','in [H+]', 'Theoretical mass', 'Structure', 'Structure dict'
            , 'MS2 Scan no', 'hitpeaklist', 'normalizedselectedintensity'
            , 'F', 'H' , 'N', 'S', 'G', 'Na', 'Questionable'
            , 344.1547, 374.1653, 376.1966, 406.2072, 432.2071
            , 450.2334, 464.2490, 580.2964, 793.3808, 825.4227]
    tmpopt = tmpopt[cols]
    tmpopt = tmpopt.reset_index(drop=True)
    print(tmpopt.columns)
    file=input("Enter the file name")
    file = file + ".csv"
    tmpopt.to_csv(file, index=False)
    print("Finish exporting", file)
    
    '''
convert column name from string to float and use round to generate
same name and continue merging
OR
clone the ppm comparison code and modify the mechanics inside it
(more calculation)
    '''
    
def fixhitlistbugmethodB():
    #methodA: convert to integer and assign correct value
    #methodB: ppm comparison
    #Don't forget to check if the output in comparepeaklist outputs 'compare value'
    #not the real one tho it may be needed in future.
    
    filename = "liteforPCA20221213v2.csv"
    df = pd.read_csv(filename)
    df.sort_index(inplace=True)
    #init the empty list
    tmpopt = []
    tmpseries = ()
    tmp = ()
    for i in range(len(df)):
        #use .copy to avoid chained settingwithcopywarning
        tmp = df.iloc[i].copy()
        peaklistparser = tmp['hitpeaklist']
        intensityparser = tmp['normalizedselectedintensity']
        parsedpeaklistheader = list(ast.literal_eval(peaklistparser))
        parsedintensity = list(ast.literal_eval(intensityparser))
        #print("peaklist", type(parsedpeaklistheader), parsedpeaklistheader)
        #print("intensity", type(parsedintensity), parsedintensity)
        for j in range(len(parsedpeaklistheader)):
            #print(parsedpeaklistheader[j])
            #print(parsedintensity[j])
            tmp[parsedpeaklistheader[j]] = parsedintensity[j]
        tmpopt.append(tmp)
    tmpopt = pd.DataFrame(tmpopt)
    tmpopt = tmpopt.reset_index(drop=True)
    tmpopt = tmpopt.fillna(0)
    a = list(tmpopt.columns)
    b = []
    for atoint in a:
        try:
            atoint = int(atoint)
        except:
            atoint = str(atoint)
        b.append(atoint)
    #comparelist = [344.1547, 374.1653, 376.1966, 406.2072,
    #432.2071, 450.2334, 464.2490, 580.2964, 793.3808, 825.4227]
    for bconv in b:
        if bconv == 1:
            pass
    print(tmpopt.columns)
    print('modified?', b)
