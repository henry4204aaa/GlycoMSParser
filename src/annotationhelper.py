import csv
import math
from collections import namedtuple
import re
import pandas as pd
import numpy as np
import ast
from decimal import *

def fillmissinginfo():
    #askopenfile 2 times if user doesn't pass the filename
    #funcname(datawithnewinfo, originaldata)
    newfile = "revised_zfNGbrain_20221129_1.csv"
    annotated = "ms2filteredselectedannotatedlist.csv"
    dfa = pd.read_csv(annotated, sep='\t')
    dfb = pd.read_csv(newfile, sep='\t')
    dfa.sort_index(inplace=True)
    dfb.sort_index(inplace=True)
    modified = []
    jred = 0
    #pretty slow since it run through the whole table
    #change the search column to user input variable
    for i in range(len(dfa)):
        tmp = dfb.iloc[i]
        #print(tmp)
        for j in range(len(dfb)):
            #do something to prevent j starts from 0
            #write this test function outside
            tmp2 = dfa.iloc[j]
            if tmp.loc['MS1scan no'] == tmp2.loc['MS1scan no']:
                tmp3 = tmp2['Structure']
                print('tmp3=',tmp3)
                print(tmp.name, 'find')
                dfb.loc[i, 'Structure'] = tmp3
                print('now the tmp is',dfb.loc[i])
                #modified.append(tmp)
                #here testing block, assigning jred to jth elements which has finished comparing
                jred = j
                break
    dfb.to_csv('revised_zfNGbrain_20221129_1_updatedlist.csv', sep='\t')
    '''
    for i in range(len(dfa)):
        tmp = dfa.iloc[i]
        for j in range(len(dfb)):
            #do something to prevent j starts from 0
            #write this test function outside
            tmp2 = dfb.iloc[j]
            if tmp.loc['MS1scan no'] == tmp2.loc['MS1scan no']:
                tmp3 = tmp2[-3:]
                #print(tmp3)
                #print(tmp.name, 'find')
                tmp = pd.concat([tmp, tmp3])
                modified.append(tmp)
                #here testing block, assigning jred to jth elements which has finished comparing
                jred = j
                break
    '''
    #output = pd.DataFrame(modified)

    #output.to_csv('updatedlist1.csv', sep='\t')
    #print(finished)
    '''
    dfa.loc[(dfb['entry no'] == dfa['entry no']), 'hitpeaklist'] = dfb['hitpeaklist']
    dfa.loc[(dfb['entry no'] == dfa['entry no']), 'selectedpeakintensity'] = dfb['selectedpeakintensity']
    dfa.loc[(dfb['entry no'] == dfa['entry no']), 'normalizedselectedintensity'] = dfb['normalizedselectedintensity']
    '''

def structuretodictparser():
    #should be selectable in future
    newdf = pd.read_csv("revised_zfNGbrain_20221129_1_updatedlist.csv", sep='\t', index_col=[0])
    exportdf = []
    #handling the data in a df
    newdf.insert(8, 'Structure dict',"")
    for i in range(len(newdf)):
        parser = newdf.iloc[i]['Structure']
        strucdict = {}
        strucdictcom = pd.DataFrame()
        tmpseries = ()
        #start from F
        fuc = re.search('F\d+', parser)
        hexose = re.search('H\d+', parser)
        hexnac = re.search('N\d+', parser)
        neu5ac = re.search('S\d+', parser)
        neu5gc = re.search('G\d+', parser)
        #ambigious strings, get nothing if not present
        Naadduct = re.search('Na', parser)   # if not Naadduct: do nothing ; else add an adduct
        sus = re.search('\?', parser)
        #not included in the first test data
        kdn = re.search('Kdn\d+', parser)
        NH3adduct = re.search('NH3', parser) #may be removed if we calculate this in other function
        Kadduct = re.search('K', parser)
        #conbine the  info into single list in each loop
        if not fuc:
            strucdict.update({'F': 0})
        else:
            strucdict.update({'F': (fuc.group()[1:])})
        if not hexose:
            strucdict.update({'H': 0})
        else:
            strucdict.update({'H': (hexose.group()[1:])})
        if not hexnac:
            strucdict.update({'N': 0})
        else:
            strucdict.update({'N': (hexnac.group()[1:])})
        if not neu5ac:
            strucdict.update({'S': 0})
        else:
            strucdict.update({'S': (neu5ac.group()[1:])})
        if not neu5gc:
            strucdict.update({'G': 0})
        else:
            strucdict.update({'G': (neu5gc.group()[1:])})
        # I won't add them as default since they are exceptions.
        if not Naadduct:
            pass
            #strucdict.update({'Na': 0) if the loop on the list throws an error, try this
        else:
            strucdict.update({'Na': 1})
        if not sus:
            pass
            #strucdict.update({'Questionable': 0) if the loop on the list throws an error, try this
        else:
            strucdict.update({'Questionable': 1})
        tmpseries = newdf.iloc[i][1:]
        tmpseries['Structure dict'] = strucdict
        exportdf.append(tmpseries)
    exportdf = pd.DataFrame(exportdf)
    exportdf = exportdf.reset_index(drop=True)
    print(exportdf)
    exportdf.to_csv("revised_zfNGbrain_20221129_1_updatedlist_parsed.csv", sep='\t')
    #replace the iloc[i]['Structure']  <- SettingWithCopyWarning and the df is not modified.
    
#it's not correct
  
'''
    if parser.find('F') == -1 :
        strucdict['F'] = 0
    else:
        strucdict['F'] = parser[(parser.index('F')+1)]
    print strucdict
'''



def deftheoreticalms1():
    #should be an element class holding these values
    df = pd.read_csv("revised_zfNGbrain_20221129_1_updatedlist_parsed.csv", sep='\t', index_col=[0])
    #askopenfile?
    theomasslist = []
    newdf = pd.DataFrame()
    #print(df.dtypes)
    getcontext().prec = 8
    H = Decimal(1.007825) 
    C = Decimal(12.000)
    O = Decimal(15.9994)
    #N = 14
    Nadif = Decimal(22.989769) - Decimal(1.007825)   #replace H with Na
    mperMedHex = Decimal(174.08921)
    mperMeHex = Decimal(204.09977)
    mperMeHexNAc = Decimal(245.12632)
    mperMeNeu5Ac = Decimal(361.17367)
    mperMeNeu5Gc = Decimal(391.18423)
    mperMefreeend = Decimal(46.04186)  #CH3COH3
    mperMeterminal = Decimal(15.023475) #CH3, in positive mode (H+)
    #read the df
    #for testing: obsmass = {'F': 1, 'H':3, 'N': 2, 'S': 0, 'G': 0 }
    #get the structure dict
    for i in range(len(df)):
        tmp = df.iloc[i]['Structure dict']
        obsmass = ast.literal_eval(tmp)
        theomass = 0
        #calculate it
        #print(obsmass, obsmass['F'], type(obsmass['F']))
        Fuc = Decimal(obsmass['F'])
        Hex = Decimal(obsmass['H'])
        HexNAc = Decimal(obsmass['N'])
        Ac = Decimal(obsmass['S'])
        Gc = Decimal(obsmass['G'])
        keys = list(obsmass)
        #print('keys', keys[-1])
        if (keys[-1] == 'Na') or (keys[-2] == 'Na'):
            #print('find Na')
            sodium = 1
        else:
            sodium = 0
        if(keys[-1] == 'Questionable'):
            print('it has a mass difference from theoretical value')
        theomass = str((Fuc * mperMedHex) + (Hex * mperMeHex) + (HexNAc * mperMeHexNAc) + (Ac * mperMeNeu5Ac) + (Gc * mperMeNeu5Gc) + mperMefreeend + H + (Nadif * sodium))
        theomasslist.append(theomass)
    print(theomasslist)
    df.insert(5, 'Theoretical mass' , theomasslist, allow_duplicates = True)
    print(df)
    df.to_csv('Annotated_20221130.csv', sep='\t', index= False)
    print('finish writing csv')
    #add a new column and insert


def fix221115():
    #quick load  df = pd.read_csv("Annotated_20221114.csv", sep='\t')
    annotated = "Annotated_20221114.csv"
    dfa = pd.read_csv(annotated, sep='\t')
    dfa.sort_index(inplace=True)
    print(dfa)


'''
df = pd.read_csv("Annotated_20221114.csv", sep='\t')
                     
df
                     
     entry no  ...                        normalizedselectedintensity
0        8002  ...  [32.815755904120806, 6.984713114752358, 2.8406...
1        8213  ...            [6.604417678607369, 32.815755904120806]
2        8214  ...            [6.604417678607369, 32.815755904120806]
3        8362  ...  [5.163043743480633, 0.6416186619881706, 1.3002...
4        8363  ...  [5.163043743480633, 0.6416186619881706, 1.3002...
..        ...  ...                                                ...
259     13018  ...  [1.1388082011049414, 2.203262544260023, 62.952...
260     13139  ...  [1.1643946131449447, 1.262355235614242, 9.4110...
261     13235  ...  [1.262355235614242, 4.893756787743044, 6.76365...
262     13418  ...  [1.3950007899305432, 5.163043743480633, 1.2576...
263     13551  ...  [4.893756787743044, 9.41103570475239, 1.281006...

[264 rows x 14 columns]
df.columns
                     
Index(['entry no', 'MS1scan no', 'MS1Isolation mass', 'MS1monoIsomass',
       'chargeState', 'Theoretical mass', 'in [H+]', 'intensity', 'Structure',
       'MS2 Scan no', 'peaklist', 'hitpeaklist', 'selectedpeakintensity',
       'normalizedselectedintensity'],
      dtype='object')
del df['MS1Isolation mass']
                     
df.columns
                     
Index(['entry no', 'MS1scan no', 'MS1monoIsomass', 'chargeState',
       'Theoretical mass', 'in [H+]', 'intensity', 'Structure', 'MS2 Scan no',
       'peaklist', 'hitpeaklist', 'selectedpeakintensity',
       'normalizedselectedintensity'],
      dtype='object')
del df['MS1monoIsomass']
                     
del df['chargeState']
                     
del df['intensity']
                     
del df['in [H+]']
                     
del df['peaklist']
                     
del df['selectedpeakintensity']
                     
df
                     
     entry no  ...                        normalizedselectedintensity
0        8002  ...  [32.815755904120806, 6.984713114752358, 2.8406...
1        8213  ...            [6.604417678607369, 32.815755904120806]
2        8214  ...            [6.604417678607369, 32.815755904120806]
3        8362  ...  [5.163043743480633, 0.6416186619881706, 1.3002...
4        8363  ...  [5.163043743480633, 0.6416186619881706, 1.3002...
..        ...  ...                                                ...
259     13018  ...  [1.1388082011049414, 2.203262544260023, 62.952...
260     13139  ...  [1.1643946131449447, 1.262355235614242, 9.4110...
261     13235  ...  [1.262355235614242, 4.893756787743044, 6.76365...
262     13418  ...  [1.3950007899305432, 5.163043743480633, 1.2576...
263     13551  ...  [4.893756787743044, 9.41103570475239, 1.281006...

[264 rows x 7 columns]
df.columns
                     
Index(['entry no', 'MS1scan no', 'Theoretical mass', 'Structure',
       'MS2 Scan no', 'hitpeaklist', 'normalizedselectedintensity'],
      dtype='object')

df.to_csv("Annotatedlite20221115.csv", index=False)
''' 


