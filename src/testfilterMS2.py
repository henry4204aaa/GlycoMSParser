import csv
import math
from decimal import Decimal,getcontext, ROUND_HALF_UP
from collections import namedtuple
import re
import pandas as pd
import comparepeaklist

list1 = (344.1698)
roundedlist = []
roundedlist2 = []
tmplist = []
MS2peaklist = namedtuple('MS2peaklist', ('dMass', 'dIntensity'))

#this works on defined list
'''
with open('no23083isMS2.tsv', newline='') as file1:
    reader = csv.reader(file1, delimiter='\t', quotechar='|')
    next(file1) #skip header
    for row in reader:
        #a = row[1]
        #b = Decimal(a
        #print(type(a))
        #print(row[1])
        
        b = Decimal((row[1])).quantize(Decimal(".0001"),rounding=ROUND_HALF_UP)
        c = float(b)
        roundedlist2.append(c)
print(roundedlist2)
'''
###



#try to extract info 
def readms2csv():
    with open('the MS2 summary from350.csv', newline='') as file1:
        reader = csv.reader(file1, delimiter='\t', quotechar='|')
        next(file1) #skip header
        rows = list(map(tuple, reader))
        print(rows)
        print(type(rows))
        c = rows
        print('c is', c[0])
        '''for row in reader:
            MS2peaklist = tuple(row[9])
            print(MS2peaklist[1])
            print('next line')
        '''
        #b = Decimal((row[1])).quantize(Decimal(".0001"),rounding=ROUND_HALF_UP)
        #c = float(b)
        #roundedlist.append(c)
#print(roundedlist2)

#return the tuple of both peak mass (rounded 0.01) and float intensity
# as a dict or (m, i) or in separated list (better for future design but unreadable by manwork
# get the whole list info and add simple annotation (referring previous research results)
#question: cant recognize correctly. See glypick results and extract certain spectrum for testing
def readms2csvpdreg():
    df = pd.read_csv('the MS2 summary from44218.csv', sep='\t')
    print(df.dtypes)
    hitspectrumlist = []
    err = 0
    notfound = 0
    hits = 0
    ms2scanno= 0
    for i in range(len(df)):
        inputlist = []
        #print(df.loc[i, "entry no"], df.loc[i, "peaklist"])
        #parse peaklist
        tmp = df.loc[i, "peaklist"]
        #tmp = tmp[2]
        #print(type(tmp), tmp)
        #the peak intensity and peak mass are included as well.
        inputlist1 = tmp.split('(')[2].split(')')[0].split(',') #peakmass
        d = tmp.split('(')[3].split(')')[0].split(',') #peakintensity
        #print('inputlist =',type(inputlist), inputlist, 'd=',type(d), d)
        #print('going to modify the number')
        for j in range(len(inputlist1)):
            k = Decimal((inputlist1[j])).quantize(Decimal("0.01"),rounding=ROUND_HALF_UP)
            #k = Decimal((inputlist1[j])).quantize(Decimal(".0001"),rounding=ROUND_HALF_UP)
            k = float(k)
            inputlist.append(k)
        #print('modified list=', inputlist)
        #loss one OHCH3= 32.0419
        comparelist = [344.17, 374.18, 376.20, 406.21, 432.22, 450.23, 464.25, 580.30, 793.38, 825.42]
        #print('compare list=', comparelist)
        result = comparepeaklist.comparepeaklist(inputlist,comparelist)
        #print(result)
        #######pseudocode
        # comparelist=input("fragments for searching in spectra")
        # check if all input is float type
        # if the user input is invalid, return default list or trap in a editing loop
        #
        if result[0] > 0:
            #print('this spectrum', i+1, 'doesnt have any hits')
            ms2scan = df.loc[i, "MS2 Scan no"]
            hitspectrumlist.append([ms2scan, result[0], result[1]])
            hits+=1
            print('hits')
            
        elif result[0] == 0:
            #print('no hits')
            notfound+=1
        else:
            print('unexpected error')
            err+=1
            #do nothing
        #need to catch the returned value and decide if this spectra is ok
    print('hits', hits, 'notfound', notfound, 'err', err, 'hits', hitspectrumlist)
        
    #return hitspectrumlist
    with open('ms2searchsummary.csv', 'w', newline='') as z:
        roww = csv.writer(z)
        roww.writerows(hitspectrumlist)
    for q in range(len(hitspectrumlist)):
        print('spectrum no', int((hitspectrumlist[q][0])), 'list',
                (hitspectrumlist[q][1]))



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
