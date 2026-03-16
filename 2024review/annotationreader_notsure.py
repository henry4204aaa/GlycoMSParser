#for reading annotation file and convert it to annotated data
import os
import mspfileloader as mspload
#import corner
try:
    print("import  openpyxl to read excel files")
    import openpyxl
    #if we have pandas, use pandas for data processing
    supportxls = True
except ImportError:
    supportxls = False
    print("Unable to read excel files. Use csv annotation sheet instead.")

try:
    print("import pandas for data processing")
    import pandas as pd 
    dataprocess = True
except ImportError:
        print("You have not install packages for data processing")
        print("Using native openpyxl command to extract annotation data")
        dataprocess = False


#v0.4 GlycoMSP headers for annotated excel file:
MS2listheaders = ["Structure", "MS2scan_no", "peak", "charge", "mass shift", "adduct", "diff profile (integer for groups if different from same annotation)", "note"]


#validate if the MS2list has correct column names. We only check headers, the values are NOT validated
def sheetheadervalidation(inputsheet, correctheaders):
    try:
        actual_headers = list(inputsheet.columns)            
        if actual_headers == correctheaders:
            print("Headers are valid.")
            return True
        else:
            print(f"Invalid headers. Expected {correctheaders}, but got {actual_headers}.")
            return False
    except Exception as e:
        print(f"[Dev-report] the input sheet is not a pandas dataframe {e}.")
        return False



def annotation_excelcheck(annotatedfile):
    #examine if file loaded is a valid excel file
    #return False when we cannot extract the MS2list as wished 
    try:
        excelfile = pd.ExcelFile(annotatedfile)  #saw that read excel will read through whole file, while this only pass an Excel object
    except Exception as e:
        print(f"Not an excel file error: {e}")
        return False
    #check if the sheet MS2list exists 
    if "MS2list" in excelfile.sheet_names:
        df = pd.read_excel(excelfile, sheet_name="MSlist")
    #    return df
        if sheetheadervalidation(df, MS2listheaders):
            return df
        else:
            return False
    else:
        print("No MS2list sheet found in the file. Please choose correct file or use single csv sheet instead.")
        return False
    #for checking annotation file data integrity. Write another function for it


    #print("Check annotated data format")
    #a = os.path.basename(annotatedfile)
    #please suggest me how to validate if a file is a excel file

def temp_annotationtest(annotatedfile, debug=True):
    test_df = annotation_excelcheck(annotatedfile)
    print(f"excel readed as {type(test_df)}")
    print(test_df.head)    



def extractdata_formatcheck():
    #read converted data and check if format and contents are correct
    print("check extracted data")
def annotation_validation():
    #confirm composition matches ms2 no (before extracting spectra for annotation)
    print("validate annotated data before extract & combining to training dataset")

def openpyxl_annotation_excel(excel_file):
    wb = openpyxl.load_workbook(excel_file)
    ws = wb["MSlist"]
    print(ws)
    #https://qiita.com/YukiYamam/items/d0e38921f5b51a950c40
    #https://openpyxl.readthedocs.io/en/stable/tutorial.html

def read_annotation_excel(excel_file, supportxls):
    print("read excel")
    if supportxls is False:
        print("You dont have packages for reading excel files. Please install openpyxl for reading excel")
        print("You can read csv annotation sheet only to solve this issue if you don't want to add new packages...")
        print("[FeaturePending]: ask for a csv annotation file here")
    elif dataprocess is False:
        print("Please install pandas to process data. You can install it by running 'pip install pandas'")
        openpyxl_annotation_excel(excel_file)
        print("Using Openpyxl to get sheet data. No more processing in this version.")
        #Exactly we can still deal with it. But not for now... making useless efforts here isnt good
    else:
        print("[DEBUG?] validate if the file is an excel file")
        valdf = annotation_excelcheck(excel_file)
        if valdf:
            print(f"MS2list found and headers were validated, now perform data comparison and calculate the glycan mass")
            print(valdf.head())
            
        else:
            print("MS2list not found or the headers were not matched ")
        

def read_csv(csv_file):
    print("Haha read csv")
    return csv_file


def selectprojectfile():
    print("Selecet 1 annotation file and one extracted data file")
    print("Select annotation file first")
    annotatefile = mspload.fileloader("excel")
    annotatesheet = read_annotation_excel(annotatefile, supportxls)
    print("Select extracted data file")
    extractfile = mspload.fileloader("csv")
    extractsheet = read_csv(extractfile)
    return annotatesheet, extractsheet
    

def extractannotationdata(annotationdata, extractedms2data):
    print("main part")


file = selectprojectfile()
print(type(file[0]))
print("annotate file ---- csv file")
print(type(file[1]))
#temp_annotationtest(file[0], debug=True)

########
'''mass information from annotationreader.py
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
    
Hex1Neu5Ac1HexNAc1 = 204.09977+361.17367+245.12632 = 810.39976

original OCH2-OH end -> OCH2-OCH3 +H 
reduced open ring OH (1st C-> 5' O) + R-CH=O  ->  OCH3 + R-CH2OCH3   + H   -> 2O 3C 9H
= 2*15.9994 + 3*12.000 + 9*1.007825 = 77.069225bb
+NAdif = 84 approved (should be correct. Better checkup later)

810+64+21 = 895 correct


468 = reducing end HexNAc dHex (Y ion)
dHex 174 + CH2 14 + H+ = 189
HexNAc redend = 245 + (O-> OCH3 delta = 15) + (CHOH-> CHOCH3 delta = 14 )

nothing added = 482

189+245+29+18 (from H2O) + H+ =482

general formula
CH3 (non-reducing end)  + glycan units + CH3 + CH2 + H2O + H

15 + m + 15 + 14 + 18+ 1 = the mass of reduced glycan (protonated)
for non-reduced...

CH3+ glycan units + CH3 + CH3 + H?
15 + m + 15 + 15 + 1? = 46 = the mass of non-reduced glycan?

810+46=856+21= 877

'''

########