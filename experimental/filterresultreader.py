#this file is made for parsing v1 output search results
import pandas as pd
import numpy
import csv
import math
#from decimal import Decimal,getcontext, ROUND_HALF_UP
from collections import namedtuple
import re
from pathlib import Path
#pseudocode



def readsearchresults():
    #select file
    try:
        from tkinter import Tk
        from tkinter.filedialog import askopenfilename
        Tk().withdraw()
        print('Please select a csv search results')
        filename = askopenfilename()
        print('Selected', filename)
        print('You are now on', filename)
        if not filename :
            print('You haven\'t select a file. Use default one')
            d = Path(__file__).resolve().parents[1]
            filename = d / 'ms2filteredtotallistwithHmanuallyannotation.csv'
            print(filename)
            df = pd.read_csv(filename, sep='\t')
        else:
            df = pd.read_csv(filename, sep='\t')
    except ImportError:
        raise ImportError('Please install tkinter')
    #default
    #search function
    print('the csv file you selected has following data structure')
    print(df.dtypes)
    print('You can try to search in certain scans or find out certain matching texts')
    #trap you in the search function. It can be another function in the future
    editflag = False
    searchflag = True
    searched = False
    dumpdestination = ""
    intflag = False
    while searchflag == True :
        n = input("enter s for searching, \n d for dumping the files, \n e for edit, \n q for exit \n")
        print(n, 'you entered')
        if n == 's' or n == 'S':
            #searching mode
            print('on searching mode')
            n = ''
            print('column labels are', df.columns, '\n with', df.index, 'rows in this data')
            c = input("\n Enter the number of row number which starts from 1, or search from column name")
            try:
                #this hits certain row
                c = int(c)
                intflag = True
                c = c-1
                if c > len(df.index):
                    print("You entered a number over the rows it has")
                    #will return to the loop top
                    continue
                else:
                    #successful search
                    find = df.iloc[c]
                    print('c is', c)
                    print(find)
                    searched = True
                    continue
            except ValueError:
                pass
            if intflag == False:
                #not int, try to find the column name
                if c in {'entry no', 'MS1scan no', 'MS1Isolation mass', 'MS1monoIsomass', 'chargeState',
                         'in [H+]', 'intensity', 'Structure', 'MS2 Scan no', 'peaklist'}:
                    #store all rows for that column going to search
                    tmp = df.loc[:,c]
                    print('now you can search in', c)
                    d = input('enter what you want to search')
                    if df.loc[:,c].dtype == "object":
                        print("it\'s object but we'll treat it as string")
                        objectformat = "string"
                    else:
                        print("It's not an object")
                    #check datatype
                    while type(d) != type(objectformat):
                        print(tmp.name,'you entered', type(d), 'while',df.loc[:,c].dtype, 'is the datatype')
                        d = input('type not matching')
                    find = df.loc[df[c] == d]
                    print('find',find)
                    searched = True
                else:
                    print('no related name found')
                    continue
            else:
                print('You should already have searched on rows')
        elif n == 'd' or n == 'D':
            existed = False
            print('on dumping mode')
            n = ''
            z = ''
            if searched == True:
                print('you have', find, 'for dumping')
                print(type(find))
                print('select a file or create a new one')
                if dumpdestination != "":
                    print('you are on', dumpdestination)
                    con = input("enter for appending on it, n for new file, o for open a file")
                    if con == "":
                        #these are redundent code, I think I should manage them as classes?
                        if isinstance(find, pd.DataFrame):
                            print('Writing dataframe to csv')
                            find.to_csv(dumpdestination, mode='a', sep='\t', header = False, index = False)
                        elif isinstance(find, pd.Series):
                            print('Writing series to csv')
                            #why series can be serious on the row/column structure? /solved
                            finddf = pd.DataFrame(find)
                            finddf = pd.DataFrame.transpose(finddf)
                            finddf.to_csv(dumpdestination, mode='a', sep='\t', header = False, index = False)
                        else:
                            with open(dumpdestination, 'w', newline='') as searchresults:
                                iterrow = csv.writer(searchresults)
                                iterrow.writerow(["entry no", "MS1scan no", "MS1Isolation mass", "MS1monoIsomass", "chargeState", "in [H+]", "intensity", "Structure", "MS2 Scan no", "peaklist"])
                                iterrow.writerows(find)
                        print("Finished dumping")
                        #close search tags to avoid you storing value twice
                        searched = False
                        z = 'done'
                if z == '':
                    z = input("n for new file, o for open a file")
                if z == 'n' or z == 'N':
                    #new file
                    dumpdestination = input('enter the file name you want')
                    dumpdestination = dumpdestination + ".csv"
                    if isinstance(find, pd.DataFrame):
                        print('Writing dataframe to csv')
                        find.to_csv(dumpdestination, sep='\t', index = False)
                    elif isinstance(find, pd.Series):
                        print('Writing series to csv')
                        finddf = pd.DataFrame(find)
                        finddf = pd.DataFrame.transpose(finddf)
                        finddf.to_csv(dumpdestination, sep='\t', index = False)
                    else:
                        with open(dumpdestination, 'w', newline='') as searchresults:
                            iterrow = csv.writer(searchresults)
                            iterrow.writerow(["entry no", "MS1scan no", "MS1Isolation mass", "MS1monoIsomass", "chargeState", "in [H+]", "intensity", "Structure", "MS2 Scan no", "peaklist"])
                            iterrow.writerows(find)
                    print("Finished dumping")
                    #close search tags to avoid you storing value twice
                    searched = False
                    z = ''
                elif z == 'o' or z == 'O':
                    #open file
                    dumpdestination = askopenfilename("csv formats", "*.csv")
                    print('Selected', dumpdestination)
                    if not dumpdestination:
                        print('You haven\'t select a file. Use default one')
                        dumpdestination = "searchresults.csv"
                        print(dumpdestination)
                    if isinstance(find, pd.DataFrame):
                        print('Writing dataframe to csv')
                        find.to_csv(dumpdestination, mode='a', sep='\t', header = False, index = False)
                    elif isinstance(find, pd.Series):
                        print('Writing series to csv')
                        finddf = pd.DataFrame(find)
                        finddf = pd.DataFrame.transpose(finddf)
                        finddf.to_csv(dumpdestination, mode='a', sep='\t', header = False, index = False)
                    else:
                        with open(dumpdestination, 'w', newline='') as searchresults:
                            iterrow = csv.writer(searchresults)
                            iterrow.writerow(["entry no", "MS1scan no", "MS1Isolation mass", "MS1monoIsomass", "chargeState", "in [H+]", "intensity", "Structure", "MS2 Scan no", "peaklist"])
                            iterrow.writerows(find)
                    print("Finished dumping")
                    #close search tags to avoid you storing value twice
                    searched = False
                    z = ""
                elif z == 'done':
                    pass
                else:
                    print('Error')
            else:
                print("You haven't search anything")
                continue
        elif n == 'e' or n == 'E':
            print('Haven\'t done')
            searchflag = False
            editflag = True
            break
        elif n == 'q' or n == 'Q':
            print('exit')
            searchflag = False
            break 
        else:
            n = input('please enter valid string')
    else:
        print('exit the search')
    while editflag == True :
        print('in editing mode')
        
    #output search result to a user-defined file
    
    #block to add certain annotation to the file(before the save function)
    
#def reader():
    
def checkifcsvisvalid():
    print('Hello')
