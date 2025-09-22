version = 0.93
last_update = 20250922
#v0.93 add negative label (non-glycans) back to manual annotation workflow
#v0.9 workable file and awaiting to be merged to main workflow. Validation prototype built and tested.
#v0.8 workable file (w/o validation)


#Refactor version of annotationhelper.py

#functions awaiting
#calculate proper mass for sugar based on derivatization, adducts
#export intermediate file with proper filename
#auto-exporting combined file with proper filename
#GUI (tkinter) that can - view ion list, sample information, metadata and
#Batch processing for files using same protocols (need to bind csv to excel before processing)
#save logs, information and ion list export/import 

import os
import re
import pandas as pd
from ast import literal_eval
import numpy as np
#add R before path string in Windows environment


# 20250411 shall we add this information in metadata?
# save mass and saccharide information in another file seems to be a good idea?
# import saccharideinfo as sugarinfo
def import_sugarinfo(type, mode, extramass = 0):
    if type == "PerMe":
        print("Set config to Perme Profile")
        #monosaccharide_masses = metadata.permethylated_masses
        #adducts_masses = sugarinfo.adducts_masses  #always the same
        #elemental_masses = sugarinfo.elemental_masses #always the same
        #etc_masses = sugarinfo.etc_masses
    elif type == "Native":
        #monosaccharide_masses = sugarinfo.native_masses
        print("Set config to Native Profile")
    if mode == "reduced":
        print("Set reducing end profile to reduced")
        #terminalmass = sugarinfo.reduded_mass
    #return monosaccharide_masses, adducts_masses, elemental_masses, etc_masses, terminalmass, extramass


#####PERMETHYLATION TEMPLATE#####
#find out CHARACTERS (caps only) have following digits  
pattern = re.compile(r'([A-Z]+)(\d+)')
patternex = re.compile(r'([A-Z]+)(\d*)')  #allow no digits
monosaccharide_masses = {
'F': 174.08921,
'H': 204.09977,
'N': 245.12632,
'S': 361.17367,
'G': 391.18423,
'KDN': 320.1465,
'K': 320.1465
}
adducts_masses = {
'H': 1.007825,
'Na': 22.9898,
'NH3': 17.0272,
#K,SO3-, H-, ...
}

elemental_masses = {
'H': 1.007825,
'Na': 22.9898,
'N': 14.0037,
'O': 15.9995
}

etc_masses = {
'CH3': 15.0235,
'CH2': 14.0157,
'H2O': 18.0151,
'PerMe(Reduced)': 62.0778,
'PerMe(Freeend)': 46.0465, #not sure if it's 46.0465 or 46.0783
'Others': 0 #should raise dialogbox to ask user enter a value
}
#reduced 15(non-reducing end) + 15(CH3) + 14(CH2) + 18(H2O)
#nonreduced 15(non-reducing end) + 15 (CH3) + 31 (CH2OH?)


# 20250921 GPT supported addition [attention, may be not what I want]
# --- add/extend mass tables near existing dicts ---
monosaccharide_masses.update({
    # GlcA as HexA under PerMe; tune if you refine later
    'HexA': 218.09502,   # placeholder PerMe-HexA; adjust when you finalize
})
# Functional groups that add to the neutral framework (covalent)
functional_group_masses = {
    'SO3': 79.956815,    # sulfate (SO3)
    'PO3H': 79.966331,   # HPO3 ; adjust if you prefer PO3 or H2PO3
}

# negative-mode friendly adducts (in addition to H/Na/NH3)
adducts_masses.update({
    'Cl': 34.968853,     # [M+Cl]− common in ESI−
    'Ac': 59.013851,     # acetate [M+CH3COO]− if present
})

H_MASS = 1.007825

def _parse_functional_groups(row):
    """
    Optional: parse functional groups count if you store them (e.g., 'SO3', 'PO3H').
    For now, read from columns 'SO3' and 'PO3H' if present; else 0.
    """
    so3 = int(row['SO3']) if 'SO3' in row and pd.notna(row['SO3']) else 0
    pgr = int(row['PO3H']) if 'PO3H' in row and pd.notna(row['PO3H']) else 0
    return so3, pgr

def calc_neutral_mass(row, deri_mode_or_name='PerMe'):
    """
    Sum monosaccharides + reducing end (etc) + functional groups to get NEUTRAL mass.
    Replaces the old 'calc_protonated_mass' baseline.
    """
    annotation = row['Structure']
    total_mass = 0.0
    composition = {}

    # parse composition like F1H5N2S1...
    matches = pattern.findall(annotation)
    for comp, count in matches:
        composition[comp] = int(count)

    # core monosaccharides (PerMe masses by default)
    for comp, mass in monosaccharide_masses.items():
        count = composition.get(comp, 0)
        total_mass += count * mass

    # functional groups, if any
    so3, pgr = _parse_functional_groups(row)
    total_mass += so3 * functional_group_masses['SO3']
    total_mass += pgr * functional_group_masses['PO3H']

    # reducing end / derivatization term (existing logic kept)
    etc = row.get('derivatization', None)  # may be 'PerMe(Reduced)' or 'PerMe(Freeend)' etc
    if isinstance(etc, float):
        total_mass += etc
    else:
        total_mass += etc_masses.get(etc, 0.0)

    return round(total_mass, 4)

def calc_precursor_mz_from_neutral(row, mode='+'):
    """
    Convert neutral mass -> theoretical precursor m/z using charge/adduct/mass_shift.
    mode: '+' or '-'
    """
    charge = int(row['charge'])
    mass_shift = float(row.get('mass shift', 0))
    adduction = str(row.get('adduct', '')).strip()

    # parse adduct string like "Na", "2H", "Cl", "H2", etc.
    adduct_mass = 0.0
    matches = patternex.findall(adduction)
    for elem, count in matches:
        n = int(count) if count else 1
        adduct_mass += n * elemental_masses.get(elem, adducts_masses.get(elem, 0.0))

    neutral = float(row['theo neutral mass'])

    # number of protons implicitly added/removed (default 1 if adduct string doesn’t include H explicitly)
    # keep behavior close to your current: +mode assumes +H, −mode assumes −H
    nH = 1 if ('H' not in adduction and 'h' not in adduction) else 0

    if mode == '+':
        precursor = neutral + adduct_mass + nH * H_MASS + mass_shift
        denom = max(1, abs(charge))
    else:  # negative mode
        precursor = neutral + adduct_mass - nH * H_MASS + mass_shift
        denom = max(1, abs(charge))

    return round(precursor / denom, 4)

def adding_masses_mode_aware(df, derivatization_tag, mode='+'):
    """
    Replacement for adding_protonated_and_observed_mass():
    - writes 'theo neutral mass'
    - writes 'theo protonated mass' (if +) or 'theo deprotonated mass' (if −) for debug
    - writes 'obs' as precursor m/z
    """
    # put derivatization string into a column if not present (keeps your etc path)
    if 'derivatization' not in df.columns:
        df = df.copy()
        df['derivatization'] = derivatization_tag

    df['theo neutral mass'] = df.apply(lambda r: calc_neutral_mass(r, derivatization_tag), axis=1)

    if mode == '+':
        df['theo protonated mass'] = df['theo neutral mass'] + H_MASS  # debug/compat
    else:
        df['theo deprotonated mass'] = df['theo neutral mass'] - H_MASS  # debug

    df['obs'] = df.apply(lambda r: calc_precursor_mz_from_neutral(r, mode), axis=1)
    return df




#20241011: add type, replaced etc_masses.get('reduced')
#need to add mode (positive or negative) and adducts in future

def calc_protonated_mass(row,etc,debug = False):
    annotation = row['Structure']
    total_mass = 0.00
    composition = {}
    #search sugar units from annotation and record in a dict
    matches = pattern.findall(annotation)
    for comp, count in matches:
        composition[comp] = int(count)
    for comp, mass in monosaccharide_masses.items():
        count = composition.get(comp, 0)  # return 0 if there is no such an element
        total_mass += count * mass
    #add adduct and reducing end option, while latter one should be able to manage somewhere
    #if debug:
    #    print(f"current total mass after glycan unit calculation is {total_mass}")
        #print(f"adducts is: {adduction} and the mass is {adducts_masses.get(adduction)}")
    #    print(f"reduced glycan extra mass is {etc_masses.get('reduced')}")
    # 20250413 added custom
    if type(etc) is float:
        total_mass +=  1.007325 + etc
    else:
        total_mass += 1.007325 + etc_masses.get(etc)#('reduced')#get(type)
    total_mass = round(total_mass, 4) #note that we have many digits as float... should be ok
    return total_mass

#always apply this function after calculating protonated mass
def calc_observed_mass(row):
    charge = row['charge'] 
    mass_shift = row['mass shift'] 
    adduction = str(row['adduct'])
    adduct_mass = 0
    obs_mass = 0
    adducts = {}
    #search adduct units from annotation and record in a dict
    matches = patternex.findall(adduction)
    for elements, count in matches:
        if count:
            adducts[elements] = int(count)
        else:
            adducts[elements] = 1
    for elements, mass in elemental_masses.items():
        count = adducts.get(elements, 0)  # return -1 if there is no such an element
        adduct_mass += count * mass    
    #print(f"[dev]: adduct mass is {adduct_mass} while charge is {charge} and mass shift is {mass_shift}")
    
    #now we have adduct mass, calculate observed mass.
    try:
        inHmass =  row['theo protonated mass'] 
        #print(f"[dev] protonated mass is {inHmass}.")
        obs_mass = (inHmass - 1.007825 + adduct_mass + float(mass_shift))/float(charge)
        obs_mass = round(obs_mass, 4) #note that we have many digits as float... should be ok
        return obs_mass
    except Exception as e:
        print(f"[DEBUG-report]please calculate protonated mass first {e}")
        return "NaN"


def adding_protonated_and_observed_mass(df, derivatization): #type should be solved here
    #df['theo protonated mass'] = df.apply(calc_protonated_mass, axis = 1)
    df['theo protonated mass'] = df.apply(lambda row: calc_protonated_mass(row, derivatization), axis=1) #GPT suggestion to use lambda
    df['obs'] = df.apply(calc_observed_mass, axis = 1)
    return df

def extractannotation(annotationdf, extractiondf, debug = False):
    if debug:
        print(extractiondf.head())
        print(annotationdf.head())
        print(f"merged dafaframe by MS2scan_no{extractiondf.columns}")
        print(f"annotationdf columns are {annotationdf.columns}")
    extracted_df = pd.merge(annotationdf, extractiondf, on="MS2scan_no", how="inner")
    #no error handling part now
    if debug:
        #temporarily save the merged file for debugging. Complete this part with writing information to log file as well.
        #extracted_df.to_csv("annotated_intermediate1.csv", index=False)
        print("temporarily disabled intermediate file export after debug finished v0.8")        
    return extracted_df

#slice the combined dataframe to the columns we need
#added MS2scan_no for further id tracking (developing)
def slice_combined_df(df):
    return df[["protonatedmass","peaklist","peakintensity","Structure", "IUPACname(optional)","Glycanannotation2", "MS2scan_no"]] 


def intensity_normalization(ion, intensity):
    #added to avoid errors
    if intensity < 0:
        intensity = 0 #need to output a warning
    intensity = np.log10(intensity) + 1  #plus 1 shift after finished filling zero values
    return [ion, intensity]


def peaks_ppm(inputpeaklist, peakintensity, reflist, setppm):
    #print(f"The items in reference list is {len(reflist)}")
    hitlist = []
    peaklist = inputpeaklist.copy()  # 建立一個副本來進行操作，避免改動原來的列表
    for ref in reflist:
        for peaks in peaklist:  # 遍歷整個 peaklist
            hitpeakinfo = []  # 清空 hitpeakinfo
            peakindex = 0  # 初始化 index
            # 檢查 peaks 是否符合條件
            if peaks <= ref:
                ppm = 1000000 * (ref - peaks) / ref
                if ppm <= setppm:  # 符合 ppm 條件
                    peakindex = peaklist.index(peaks)
                    hitpeakinfo = [ref, peakintensity[peakindex]]
                    hitlist.append(hitpeakinfo)
            elif peaks > ref:
                ppm = 1000000 * (peaks - ref) / peaks
                if ppm < setppm:  # 符合 ppm 條件
                    peakindex = peaklist.index(peaks)
                    hitpeakinfo = [ref, peakintensity[peakindex]]
                    hitlist.append(hitpeakinfo)
                else:
                    # 不再切割 peaklist，直接跳過這次的 ref
                    break
    #print(f"[dev] hit list is {hitlist}")
    nhitlist = []
    #peaks_ppm.counthelp +=1
    #print(f"[dev] count is {peaks_ppm.counthelp}")
    for i in range(len(hitlist)):
        #test without normalization
        #appendlist = (hitlist[i][0], hitlist[i][1])
        appendlist = intensity_normalization(hitlist[i][0], hitlist[i][1])
        nhitlist.append(appendlist)
    #print(f"[dev] normalized hit list is {nhitlist}")
    zeroionlist = []
    for logi in reflist:
        # 檢查 nhitlist 裡所有元素的第一個值是否匹配 logi
        if logi not in [item[0] for item in nhitlist]:
            zeroions = [logi, 1]
            nhitlist.append(zeroions)
            zeroionlist.append(zeroions)
        else:
            pass
    #print(f"[dev] after adding 0 value, hit list is {nhitlist}")
    #print(f"[dev] check zero ion list: {zeroionlist}")
    #filling extra value +1 and fill blank with 0+1 to avoid division error
    #print(f"[dev] normalized hit list is {nhitlist}")
    return nhitlist
#peaks_ppm.counthelp = 0

def findingions(dfrow, iondf, ppm):
    peaks = list(dfrow['peaklist'])
    intensity = list(dfrow['peakintensity'])
    #print(type(peaks))
    #print(f"peaks are {peaks[0:5]}")
    ions = []
    ions = iondf.values.flatten().tolist() #values in object
    return peaks_ppm(peaks, intensity, ions,ppm)

#wrap below codes into a function
#need variables: ion_df, pre_df, setppm
#which resolves the issue of string to list conversion in peaklist and peakintensity. Modify the dataframe itself.
def expandpeaklist(df):
    df = df.copy()
    df.loc[:, 'peaklist'] = df['peaklist'].apply(lambda x: literal_eval(x.strip('()')))
    df.loc[:, 'peakintensity'] = df['peakintensity'].apply(lambda x: literal_eval(x.strip('()')))
    return df

#input a read excel sheet that is ionlist, export the column names for normalized ion list
def extract_ionmasslist(ionmass_sheet):
    ion_df = ionmass_sheet[["mass"]]
    iondfindex = ion_df.values.flatten().tolist()
    iondfindex.extend(['Structure', 'IUPACname(optional)', 'Glycanannotation2'])
    iondfindex.insert(0, 'protonatedmass')
    return iondfindex

#read excel file 
#already existed in annotationreader.py: annotation_excelcheck(excel_file)

#####copied from mspdatavalidation.py #####

#convert object to string
def is_string_column(series):
    return series.map(lambda x: isinstance(x, str) or pd.isna(x)).all()
#validate csv file from raw#
#update columns if we change the file conversion logic#  
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

#validate annotation file and ion list in same Excel#
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
    #print(f"the columns here are {list(excel_df.columns)} and the types are {type(list(excel_df.columns))}")
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
            print(f"[debug] mismatch {excel_df.columns}")
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


#####
def directassign_files(annotation_file, raw_csv, derivatizationtags, debug = False):
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
    #df = adding_protonated_and_observed_mass(df,derivatizationtags)

    #20250921 GPT supported
    mode = '+'
    if isinstance(derivatizationtags, str) and '-' in derivatizationtags:
        mode = '-'
    elif isinstance(derivatizationtags, dict):
        mode = '-' if str(derivatizationtags.get('Mass Analyzer charge mode','+')).strip() == '-' else '+'

    df = adding_masses_mode_aware(df, derivatizationtags, mode=mode)

    merged_df = extractannotation(df, df2, debug = True)
    pre_df = slice_combined_df(merged_df)
    pre_df = expandpeaklist(pre_df)
    iondfindex = extract_ionmasslist(ion_df)
    return pre_df, iondfindex, ion_df

#xls = pd.ExcelFile(excel_file)
#print(xls.sheet_names)
#for checking annotation file data integrity. Write another function for it
#df = pd.read_excel(xls, sheet_name="MSlist")

#adding_protonated_and_observed_mass(df)  if we did not return df, we can use this
#df = adding_protonated_and_observed_mass(df)


#df2 = pd.read_csv(csv_file, sep='\t')
#print(df2.head())
#merged_df = extractannotation(df, df2, debug = True)
#pre_df = slice_combined_df(merged_df)

#print(pre_df.head())
#print(pre_df.info())
#print("Load ion filter and convert the df to proper shape and then export it. Do it.")
#read same excel file. Need to consider when the ion list is in another file
#ion_df = pd.read_excel(xls, sheet_name="ionlist")
#ion_df = ion_df[["mass"]]


#pre_df.loc[:, 'peaklist'] = pre_df['peaklist'].apply(lambda x: literal_eval(x.strip('()')))
#pre_df.loc[:, 'peakintensity'] = pre_df['peakintensity'].apply(lambda x: literal_eval(x.strip('()')))

#iondfindex = extract_ionmasslist(ion_df)

#pre_df = expandpeaklist(pre_df)

#20250905 added by GPT
#to read ion list feature counts before proceeding to merge
def peek_ion_feature_count(excel_path: str):
    """Read 'ionlist' sheet and count masses; return (count, masses_series)."""
    import pandas as pd
    xls = pd.ExcelFile(excel_path)
    # be lenient with sheet name
    cand = [s for s in xls.sheet_names if s.lower().replace("_","") in ("ionlist","ionlists")]
    if not cand:
        raise ValueError("No 'ionlist' sheet found")
    df = pd.read_excel(excel_path, sheet_name=cand[0])
    # be lenient with column name
    mass_col = next((c for c in df.columns if str(c).strip().lower() == "mass"), None)
    if mass_col is None:
        raise ValueError("No 'mass' column in ionlist")
    masses = pd.to_numeric(df[mass_col], errors="coerce").dropna()
    return int(masses.size), masses


#20250905 added by GPT
#notice that our feature extraction approach will make decoy method no effect at all. Consider adding this back in future if we have other approach
def append_feature_space_decoys(
    wide_df,
    *,
    label_col: str = "Structure",
    exclude_cols: tuple[str, ...] = (
        "Structure", "IUPACname(optional)", "Glycanannotation2",
        "MS2scan_no", "unique_ID", "protonatedmass"
    ),
    method: str = "permute",     # "permute" | "sparse_mask" | "mixup"
    ratio: float = 0.5,          # decoys per *positive* row (0.5 = half as many)
    sparse_keep: float = 0.2,    # for sparse_mask: keep 20% of nonzeros
    noise_sd: float = 0.05,      # small feature noise
    random_state: int = 42
):
    """
    Return a new DataFrame with extra Non-glycan decoys created in *feature space*.
    """
    import numpy as np, pandas as pd
    rs = np.random.RandomState(random_state)
    df = wide_df.copy()

    feat_cols = [c for c in df.columns if c not in exclude_cols]
    if not feat_cols:
        return df

    # base pool: use existing Non-glycan if present; otherwise positives
    base_neg = df[df[label_col] == "Non-glycan"]
    base = base_neg if len(base_neg) else df[df[label_col] != "Non-glycan"]
    if base.empty:
        return df

    n_decoys = int(max(1, ratio * len(df[df[label_col] != "Non-glycan"])))
    src = base.sample(n=min(len(base), n_decoys), random_state=random_state, replace=(len(base) < n_decoys))

    rows = []
    for _, r in src.iterrows():
        v = r[feat_cols].astype(float).values.copy()

        if method == "permute":
            # permute a random slice of columns
            k = max(1, int(0.5 * len(feat_cols)))
            idx = rs.choice(len(feat_cols), size=k, replace=False)
            rs.shuffle(v[idx])

        elif method == "sparse_mask":
            # keep a small subset of non-zero features, zero the rest
            nz = np.where(v != 0)[0]
            if len(nz) > 0:
                keep = rs.choice(nz, size=max(1, int(sparse_keep * len(nz))), replace=False)
                mask = np.ones_like(v, dtype=bool)
                mask[keep] = False
                v[mask] = 0.0

        elif method == "mixup":
            # average with another negative row and add small noise
            r2 = base.sample(n=1, random_state=rs.randint(0, 1_000_000)).iloc[0]
            v = 0.5 * (v + r2[feat_cols].astype(float).values)

        # light noise
        if noise_sd > 0:
            v = v + rs.normal(0.0, noise_sd, size=v.shape)

        new_r = r.copy()
        new_r[feat_cols] = v
        new_r[label_col] = "Non-glycan"
        rows.append(new_r)

    decoys = pd.DataFrame(rows, columns=df.columns)
    return pd.concat([df, decoys], ignore_index=True)

def sample_real_negatives(
    raw_tsv_path: str,
    annotated_scans,                  # iterable of MS2scan_no already in pre_df
    *,
    ion_df=None,                      # pass the ionlist df you already loaded
    ppm_tol: float = 10.0,            # mass tolerance for ion hits
    min_hits: int = 1,                # keep scans with < min_hits glycan-ion matches
    max_neg_ratio: float = 3.0,       # cap negatives to ratio × positives
    random_state: int = 42
):
    """
    Build Non-glycan rows from raw TSV scans NOT annotated.
    Gate by glycan-ion 'hits' using ion_df['mass'] and ppm tolerance.
    Returns rows aligned with your exporter expectations:
      MS2scan_no, protonatedmass, peaklist, peakintensity, Structure, IUPACname(optional), Glycanannotation2
    """
    import numpy as np, pandas as pd
    from ast import literal_eval

    raw = pd.read_csv(raw_tsv_path, sep="\t")
    # Schema is validated elsewhere: peaklist/peakintensity objects exist. 

    # Parse list-like columns if they are strings
    def to_list_safe(x):
        if isinstance(x, (list, tuple, np.ndarray)): return list(x)
        s = str(x).strip()
        try:
            return list(literal_eval(s))
        except Exception:
            # handle "(1,2,3)"-style or comma-only strings
            s = s.strip("()[]")
            return [float(z) for z in s.split(",") if z.strip()]

    raw["peaklist"] = raw["peaklist"].apply(to_list_safe)
    raw["peakintensity"] = raw["peakintensity"].apply(to_list_safe)
    raw["MS2scan_no"] = raw["MS2scan_no"].astype(int)

    annotated_set = set(map(int, annotated_scans))
    all_scans = set(raw["MS2scan_no"].unique())
    cand_scans = list(all_scans - annotated_set)
    if not cand_scans:
        return pd.DataFrame(columns=[
            "MS2scan_no","protonatedmass","peaklist","peakintensity",
            "Structure","IUPACname(optional)","Glycanannotation2"
        ])

    # Prepare ion masses
    ion_masses = None
    if ion_df is not None and "mass" in ion_df.columns:
        ion_masses = ion_df["mass"].astype(float).values

    # Quick helper: count matches to ionlist within ppm
    def count_hits(mz_list):
        if ion_masses is None or len(mz_list) == 0:
            return 0
        mz = np.asarray(mz_list, dtype=float)
        # vectorized: for each ion mass, check any mz within ppm window
        hits = 0
        for mi in ion_masses:
            tol = mi * ppm_tol * 1e-6
            if np.any((mz >= mi - tol) & (mz <= mi + tol)):
                hits += 1
        return hits

    # Summarize per-scan and gate
    rows = []
    sub = raw[raw["MS2scan_no"].isin(cand_scans)]
    g = sub.groupby("MS2scan_no", sort=False)
    for scan_id, gdf in g:
        # Each row already holds one scan; if multiple rows, take the first (they should be identical by schema)
        r = gdf.iloc[0]
        hits = count_hits(r["peaklist"])
        if hits < min_hits:
            rows.append({
                "MS2scan_no": int(scan_id),
                "protonatedmass": float(r.get("protonatedmass", np.nan)),
                "peaklist": list(map(float, r["peaklist"])),
                "peakintensity": list(map(float, r["peakintensity"])),
                "Structure": "Non-glycan",
                "IUPACname(optional)": "",
                "Glycanannotation2": ""
            })

    neg = pd.DataFrame(rows)
    if neg.empty:
        return neg

    # Cap by ratio
    n_pos = max(1, len(annotated_set))
    keep = min(len(neg), int(max_neg_ratio * n_pos))
    if len(neg) > keep:
        neg = neg.sample(n=keep, random_state=random_state).reset_index(drop=True)
    return neg


#added 20250905 old version of adding negative label (but rely on score - which is not presented in manual anno workflow)
"""
def sample_real_negatives(
    raw_tsv_path: str,
    annotated_scans,                  # iterable of MS2scan_no already in pre_df
    *,
    scan_col: str = "MS2scan_no",
    precursor_col: str = "precursor_mz",
    frag_mz_col: str = "fragment_mz",
    frag_int_col: str = "fragment_intensity",
    score_col: str = "score",         # set to None if you don't have one
    score_thr: float = 0.05,          # keep scans with max(score) < score_thr
    marker_cols: list[str] | None = None,  # e.g., ["core_marker_b", "core_marker_y"]
    marker_min: int = 1,              # require < marker_min markers to keep (0 or 1)
    max_neg_ratio: float = 3.0,       # cap negatives to ratio × positives
    random_state: int = 42
):
    """
"""
    Build Non-glycan rows from raw TSV scans that are NOT annotated.
    Gating:
      - if score_col present: keep scans with max(score) < score_thr
      - if marker_cols present: keep scans with sum(marker_flags) < marker_min
    Returns a DataFrame aligned to pre_df's expected columns:
      MS2scan_no, protonatedmass, peaklist, peakintensity, Structure, IUPACname(optional), Glycanannotation2
"""
"""
    import numpy as np, pandas as pd

    raw = pd.read_csv(raw_tsv_path, sep="\t")
    if scan_col not in raw.columns:
        raise ValueError(f"'{scan_col}' not found in {raw_tsv_path}")
    raw = raw.copy()
    raw[scan_col] = raw[scan_col].astype(int)
    annotated_set = set(map(int, annotated_scans))

    # 1) candidate scans = not annotated
    all_scans = set(raw[scan_col].unique())
    cand_scans = all_scans - annotated_set
    if not cand_scans:
        return pd.DataFrame(columns=[
            "MS2scan_no", "protonatedmass", "peaklist", "peakintensity",
            "Structure", "IUPACname(optional)", "Glycanannotation2"
        ])

    # 2) score gate
    if score_col and score_col in raw.columns:
        per_scan_score = raw.groupby(scan_col)[score_col].max()
        cand_scans &= set(per_scan_score.index[per_scan_score < score_thr])

    # 3) marker-ion absence gate (optional)
    if marker_cols:
        missing = [c for c in marker_cols if c not in raw.columns]
        if missing:
            print(f"[mspval] marker columns missing in raw TSV, skipping: {missing}")
        else:
            flags = raw.groupby(scan_col)[marker_cols].max()  # assumes 0/1
            allowed = flags.index[(flags.sum(axis=1) < marker_min)]
            cand_scans &= set(allowed)

    if not cand_scans:
        return pd.DataFrame(columns=[
            "MS2scan_no", "protonatedmass", "peaklist", "peakintensity",
            "Structure", "IUPACname(optional)", "Glycanannotation2"
        ])

    # 4) build one row per scan
    df = raw[raw[scan_col].isin(cand_scans)].copy()
    if frag_mz_col not in df.columns or frag_int_col not in df.columns:
        raise ValueError(f"'{frag_mz_col}'/'{frag_int_col}' not found in raw TSV")

    g = df.groupby(scan_col)

    def _one_scan(gdf):
        pmass = float(gdf[precursor_col].iloc[0]) if precursor_col in gdf.columns else float("nan")
        plist = gdf[frag_mz_col].astype(float).to_list()
        pint  = gdf[frag_int_col].astype(float).to_list()
        return pd.Series({
            "MS2scan_no": int(gdf.name),
            "protonatedmass": pmass,
            "peaklist": plist,
            "peakintensity": pint,
            "Structure": "Non-glycan",
            "IUPACname(optional)": "",
            "Glycanannotation2": ""
        })

    neg = g.apply(_one_scan).reset_index(drop=True)

    # 5) cap by ratio
    n_pos = max(1, len(annotated_set))
    n_keep = min(len(neg), int(max_neg_ratio * n_pos))
    neg = neg.sample(n=n_keep, random_state=random_state) if len(neg) > n_keep else neg

    return neg
"""

#autofilled by copilot. Need manual validation
def createnormailzedionlistcsv(ionindex, converted_df, ion_df, filename):
    niondf = pd.DataFrame(columns=ionindex)
    for i in range(len(converted_df)):
        annotationlist = [['protonatedmass',converted_df.iloc[i]['protonatedmass']],['Structure',converted_df.iloc[i]['Structure']], ['IUPACname(optional)' ,converted_df.iloc[i]['IUPACname(optional)']],['Glycanannotation2',converted_df.iloc[i]['Glycanannotation2']],['unique_ID',converted_df.iloc[i]['MS2scan_no']]]
        tempions = findingions(converted_df.iloc[i], ion_df, 10)
        #print(f"[DEBUG] tempions for row {i}: {tempions}")
        #tempions1 = sorted(tempions, key=lambda x: x[0])
        #print(f"[DEBUG] tempions for row after sorting {i}: {tempions}")
        data_dict = {}
        for item in tempions:
            key, value = item[0], item[1]
            # 只添加當前的 ion 當它不存在於 data_dict 中，或者它不是 1
            if key not in data_dict or (key in data_dict and value != 1):
                #print(f"[DEBUG] Adding {key}: {value} to data_dict")
                data_dict[key] = value
        #print(f"[DEBUG] data_dict for row {i} BEFORE ADDING ANNOTATION: {data_dict}")
        for item in annotationlist:
            if item[0] not in data_dict:
                data_dict[item[0]] = item[1] #fix the bug that existed ions are overwritten by annotationlist
        #print(f"[DEBUG] data_dict for row {i}: {data_dict}")
        new_row = pd.Series(data_dict)
        niondf = pd.concat([niondf, new_row.to_frame().T],ignore_index=True)
    print(f"[DEBUG] niondf after row {i}: {niondf.head()}")

    niondf.to_csv(filename, index=False)

 
def typeselection(type):
    if type == "NG" or type == "non-reduced":
        return 'nonreduced'
    elif type == "OG" or type == "reduced":
        return 'reduced'


#very primitive way:
#type = typeselection("OG")#("NG")#("OG")
#excel_file = R"G:\其他電腦\My Computer\GlycoMSParser\src\OG_int_online_annotation_revised_2024_glycan_zf_OG.xlsx"  #test validation
#excel_file = R"G:\其他電腦\My Computer\GlycoMSParser\src\20240922_temp_zf_intestine_1.xlsx"
#R"G:\其他電腦\My Computer\GlycoMSParser\src\202401012_fixed_temp_zf_ovary.xlsx"
#R"G:\其他電腦\My Computer\GlycoMSParser\src\zf_NG_brain_annotation_20240922.xlsx"
#R"C:\Users\Sakazuki\Documents\online_annotation_revised_2024_glycan_zf_OG_brain_20240817temp.xlsx"
#"G:\其他電腦\My Computer\GlycoMSParser\src\OG_int_online_annotation_revised_2024_glycan_zf_OG.xlsx"
#R"G:\其他電腦\My Computer\GlycoMSParser\src\zf_OG_ovary_annotation_202410.xlsx"
#"/Users/hnstseng/Downloads/online_annotation_revised_2024_glycan_zf_OG.xlsx"#
#R"C:\Users\Sakazuki\Documents\online_annotation_revised_2024_glycan_zf_OG_brain_20240817temp.xlsx" #0812temp is older one #R
#csv_file = R"G:\其他電腦\My Computer\GlycoMSParser\ms2_zebrafish_Yann_glycome_20241008_zf_sPerMeOG_intestine.raw.csv"
#csv_file = R"G:\其他電腦\My Computer\GlycoMSParser\ms2_zebrafish_Yann_glycome_20240922_zf_sPerMeNG_intestine.raw.csv"
#R"G:\其他電腦\My Computer\GlycoMSParser\ms2_zebrafish_Yann_glycome_20240922_20240922_zf_sPerMeNG_ovary.raw.csv"
#R"G:\其他電腦\My Computer\GlycoMSParser\ms2_zebrafish_Yann_glycome_20240921_zf_sPerMeNG_brain_20240921_zf_sPerMeNG_brain.raw.csv" 
#R"C:\Users\Sakazuki\Downloads\ms2_zebrafish_Yann_glycome_20240731_zf_PerMeOG_brain.raw.csv"
#R"G:\其他電腦\My Computer\GlycoMSParser\ms2_zebrafish_Yann_glycome_20241008_zf_sPerMeOG_intestine.raw.csv"
#R"/Users/hnstseng/Downloads/ms2_zebrafish_Yann_glycome_20240731_zf_PerMeOG_brain.raw.csv"#
#R"C:\Users\Sakazuki\Downloads\ms2_zebrafish_Yann_glycome_20240731_zf_PerMeOG_brain.raw.csv" #R 
#preprocessed_df, ionlist, ion_df = directassign_files(excel_file, csv_file, debug = False)

#createnormailzedionlistcsv(ionlist, preprocessed_df,ion_df, R"G:\其他電腦\My Computer\GlycoMSParser\zfOG_int_passedvalidation.csv")


def add_unique_id(csv, extname=None):
    basename, ext = os.path.splitext(csv)
    if extname is None:
        extname = str(input("please enter the extension name to apply to unique id"))
    tdf = pd.read_csv(csv, header= 0)
    print(tdf.head)
    tdf['unique_ID'] = tdf['unique_ID'].astype(str) + (extname)
    filename = basename + "_adduid" + ext
    tdf.to_csv(filename, index=False)
    return "adding temp id finished"

#copy paste the csv in  "createnormailzedionlistcsv" and remember 21st of the September... don't forget to change extname or you will be locked in a shoebox. Crazy? I was crazy once."
#add_unique_id(R"G:\其他電腦\My Computer\GlycoMSParser\zfOG_int_passedvalidation.csv", extname="intestine")