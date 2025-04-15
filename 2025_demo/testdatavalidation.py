version = 0.2
last_update = 20250411
#v0.2: validated working, going to inherit back to annotation_org.py
#functions inherited from annotation_org.py

import os
import re
import pandas as pd
from ast import literal_eval
import numpy as np


'''
PS G:\其他電腦\My Computer\GlycoMSParser> python .\2025_demo\mspdatavalidation.py
the columns here are ['Structure', 'MS2scan_no', 'peak', 'charge', 'mass shift', 'adduct', 'diff profile (integer for groups if different from same annotation)', 'IUPACname(optional)', 'Glycanannotation2', 'note'] and the types are <class 'list'>
Header version is v3
annotation list validated
Column 'mass' is float64.
raw converted csv is valid
PS G:\其他電腦\My Computer\GlycoMSParser> 
'''



def validate_csv_structure(csv_path):
    expected_headers = ['entry_no', 'MS1scan_no', 'MS1_isolationmass', 'MS1_monoisolationmass', 'chargeState', 'protonatedmass', 'MS2scan_no', 'peaklist', 'peakintensity']
    expected_dtypes = {
    'entry_no': 'int64',
    'MS1scan_no': 'int64',
    'MS1_monoisolationmass': 'float64',
    'MS1_isolationmass': 'float64',
    'chargeState': 'int64',
    'protonatedmass': 'float64',
    'MS2scan_no': 'int64',
    'peaklist': 'object',
    'peakintensity': 'object'
    }
    try:
        df = pd.read_csv(csv_path,index_col=None, sep='\t')

        # Check headers match
        if list(df.columns) != expected_headers:
            print("Header mismatch:")
            print("Expected:", expected_headers)
            print("Found:", list(df.columns))
            return False

        # Check dtypes
        for col, expected_type in expected_dtypes.items():
            if col not in df.columns:
                print(f"Missing column: {col}")
                return False
            actual_type = df[col].dtype.name
            if actual_type != expected_type:
                print(f"Type mismatch for column '{col}': expected {expected_type}, got {actual_type}")
                return False

        #print("CSV structure is valid.")
        return True

    except Exception as e:
        print(f"Error during validation: {e}")
        return False



def is_string_column(series):
    return series.map(lambda x: isinstance(x, str) or pd.isna(x)).all()

def validate_annotation_structure(excel_df):#extra arguments may be asked in future

    v4_expected_headers = ['Structure', 'MS2scan_no', 'peak', 'charge', 'mass shift', 'adduct', 'diff profile (integer for groups if different from same annotation)', 'IUPACname(optional)', 'Glycanannotation2', 'note', 'GlyToucan ID']
    v3_expected_headers = ['Structure', 'MS2scan_no', 'peak', 'charge', 'mass shift', 'adduct', 'diff profile (integer for groups if different from same annotation)', 'IUPACname(optional)', 'Glycanannotation2', 'note']
    basic_expected_headers = ['Structure', 'MS2scan_no', 'charge', 'mass shift', 'adduct', 'note']
    #all str should call is_string_column to convert object to str
    #add new rules both on headers and datatype condition statement in this function
    #don't check the note column
    #ok to leave blank: any field not in needcheck_expected_dtypes, they are all labels and can be used in ML training as well
    #future-feature: validate those labels and allow user to train by that column. THE METHOD should be logged as records
    needcheck_expected_dtypes = {
    'Structure': 'str' ,#F1H5N3S2 #LABEL#ESSENTIAL
    'MS2scan_no': 'int64', #12345 #IDENTIFIER (trackable only in same sample, need another id when merging multiple datasets)
    'charge': 'int64', #2
    'mass shift': 'int64', #1
    'adduct': 'str', #2H, NH4...
    }
    v3_added_dtypes = {
    'peak': 'int64', #1 #LABEL, but purely man-based, error-prone
    'diff profile (integer for groups if different from same annotation)': 'int64', #1 #LABEL#
    'IUPACname(optional)': 'str', #the one complete annotation  #LABEL#
    }
    v4_added_dtypes ={
    'Glycanannotation2': 'str', #placeholder, put some formats to test  #LABEL#
    'GlyToucan ID': 'str' #glycan in database, mind if there are stem/parent ID use that instead to avoid undefined linkage #LABEL#
    }
    dv3 = {**needcheck_expected_dtypes, **v3_added_dtypes}
    dv4 = {**needcheck_expected_dtypes, **v4_added_dtypes}
    version_flags = 0 #default is 0
    print(f"the columns here are {list(excel_df.columns)} and the types are {type(list(excel_df.columns))}")
    try:
        #df = pd.read_excel(excel_path, sheet_name="MSlist")

        # Check headers match, always start from latest version
        if list(excel_df.columns) == v4_expected_headers: 
            print("Header version is v4")
            version_flags = 4
        elif   list(excel_df.columns) == v3_expected_headers:
            print("Header version is v3")
            version_flags = 3
        elif list(excel_df.columns) == basic_expected_headers:
            print("Header version is v1, prototype")
            version_flags = 1
        else:
            print("Header mismatched, the file is invalid")
            return False

        # Check certain dtypes
        #if version_flags = 1:  #when other version flag or label validation functions are implemented, move this condition to the last part. Default minimal validation 
        for col, expected_type in needcheck_expected_dtypes.items():

            if col not in excel_df.columns:
                print(f"Missing column: {col}")
                return False

            col_data = excel_df[col]
            
            if col_data.isna().all():
                print(f"Column '{col}' is entirely empty – allowed.")
                continue  # pass for entirely empty column

            if col_data.isna().any():
                print(f"Column '{col}' has missing values – not allowed.")
                return False  # fail for partially missing column

            actual = col_data.dtype.name
            if expected_type == 'str':
                if actual != 'object' or not col_data.map(lambda x: isinstance(x, str)).all():
                    print(f"Column '{col}' is not all strings.")
                    return False
            else:
                if actual != expected_type:
                    print(f"Type mismatch for column '{col}': expected {expected_type}, got {actual}")
                    return False
        #if version_flags = 3 or IUPACname_label is True:  #enable IUPAC name based training
        ####for col, expected_type in dv3.items(): ...
        #if version_flags = 4 or IUPACname_label is True:  #enable IUPAC name based training
        ####for col, expected_type in dv4.items(): ...

        #print("CSV structure is valid.")
        return True

    except Exception as e:
        print(f"Error during validation: {e}")
        return False



#no merging df version of directassign_files
def raw_validation(annotation_file, raw_csv, debug = False):
    anno = pd.ExcelFile(annotation_file)
    notfounderror = []
    MSlistcheck, ionlistcheck, csvfilecheck = None, None, None
    #Still missing MSlist and ionlist validation 
    if "MSlist" in anno.sheet_names:
        df = pd.read_excel(anno, sheet_name="MSlist")
        if validate_annotation_structure(df):
            print("annotation list validated")
            MSlistcheck = True
        else:
            print("Annotation list invalid")
            MSlistcheck = False
        #if debug:
        #    print(f"annotation MSlist head{df.head()}, types = {df.dtypes}")
    else:
        notfounderror.append[1]
        print("No MSlist in sheet names")
        MSlistcheck = False
    #read ion list information
    if "ionlist" in anno.sheet_names:
        ion_df = pd.read_excel(anno, sheet_name="ionlist")
        ion_df = ion_df[["mass"]] #will also return this
        if ion_df["mass"].dtype == "float64":
            print("Column 'mass' is float64.")
            ionlistcheck = True
        else:
            print(f"Column 'mass' has dtype {ion_df['mass'].dtype}, expected float64.")
        ionlistcheck = False
    else:
        notfounderror.append[2]
        print("No ion list in sheet names")
        ionlistcheck = False
    #TODO raw file validation
    if validate_csv_structure(raw_csv):
        print("raw converted csv is valid")
        df2 = pd.read_csv(raw_csv, sep='\t')
    else:
        csvfilecheck = False
        print("No raw csv file valid in sheet names")
    #converting annotation file
    #add 2 mass information to dataframe
    #df = adding_protonated_and_observed_mass(df)

    #merged_df = extractannotation(df, df2, debug = True)
    #pre_df = slice_combined_df(merged_df)
    #pre_df = expandpeaklist(pre_df)
    #iondfindex = extract_ionmasslist(ion_df)
    #return pre_df, iondfindex, ion_df


##1	1	536.1658936	536.1658936	1	536.1658936	3	(96.3175277709961, 98.88307189941406, 101.24893951416016, 102.06212615966797, 102.79436492919922, 113.02568054199219, 113.75145721435547, 116.54693603515625, 120.58672332763672, 124.80850982666016, 124.99066925048828, 134.31565856933594, 135.7213592529297, 137.3096160888672, 142.7156524658203, 146.1102294921875, 147.30699157714844, 147.73255920410156, 148.69737243652344, 148.7539520263672, 150.69403076171875, 167.0554656982422, 172.38291931152344, 175.46035766601562, 206.3819580078125, 223.06370544433594, 239.09483337402344, 267.0886535644531, 281.052490234375, 283.03009033203125, 293.6178894042969, 299.0614013671875, 300.0619201660156, 307.49407958984375, 359.0290832519531, 415.037109375, 482.83709716796875, 503.1072692871094, 519.138427734375, 536.1648559570312, 537.1666259765625)	(1632.2470703125, 1846.3802490234375, 1482.30126953125, 1618.1602783203125, 1541.8468017578125, 1868.7392578125, 1689.5087890625, 1597.5814208984375, 2158.280517578125, 2136.55029296875, 1864.4212646484375, 2147.241455078125, 1912.2576904296875, 2128.126708984375, 2334.9931640625, 2246.912353515625, 2409.70361328125, 2016.370849609375, 1925.5081787109375, 2243.49853515625, 2522.318603515625, 6009.876953125, 2500.128173828125, 2772.1328125, 3015.68701171875, 3632.06005859375, 5836.669921875, 2570.6044921875, 3980.724365234375, 4580.166015625, 2757.16064453125, 98228.59375, 4958.79931640625, 2900.84716796875, 9016.7119140625, 17808.96484375, 2650.1845703125, 84050.7265625, 72264.5703125, 201270.125, 8773.1953125)

#input a read excel sheet that is ionlist, export the column names for normalized ion list
def extract_ionmasslist(ionmass_sheet):
    ion_df = ionmass_sheet[["mass"]]
    iondfindex = ion_df.values.flatten().tolist()
    iondfindex.extend(['Structure', 'IUPACname(optional)', 'Glycanannotation2'])
    iondfindex.insert(0, 'protonatedmass')
    return iondfindex

#validate_csv_structure("G:\\其他電腦\\My Computer\\GlycoMSParser\\ms2_zebrafish_Yann_glycome_20241008_zf_sPerMeOG_intestine.raw.csv")
raw_validation("G:\\其他電腦\\My Computer\\GlycoMSParser\\src\\OG_int_online_annotation_revised_2024_glycan_zf_OG.xlsx", "G:\\其他電腦\\My Computer\\GlycoMSParser\\ms2_zebrafish_Yann_glycome_20241008_zf_sPerMeOG_intestine.raw.csv", debug = True)