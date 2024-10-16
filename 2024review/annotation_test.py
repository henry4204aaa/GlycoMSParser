version = 0.7
#this file is a test file to deal with excel file and csv file handling
#those functions finished tests will be added to proper python file afterwards.
#changelog: 20240822 v0.7: general test finished. Save a new file for organizing code.


print("Finished Add current code into annotationreader.py")

import re
import pandas as pd
#add R before path in Windows environment

from ast import literal_eval
import numpy as np
#excel_file = R"C:\Users\Sakazuki\Documents\online_annotation_revised_2024_glycan_zf_OG_brain_20240817temp.xlsx" #0812temp is older one #R"/Users/hnstseng/Downloads/online_annotation_revised_2024_glycan_zf_OG.xlsx"#
#20240922 for brain NG
#excel_file = R"G:\其他電腦\My Computer\GlycoMSParser\src\zf_NG_brain_annotation_20240922.xlsx"
#20240922 for ovary NG (new rule, deleted repeat )
#excel_file = R"G:\其他電腦\My Computer\GlycoMSParser\src\20240922_temp_zf_ovary.xlsx"
#20240922 for intestine NG (2024reviewed sheet, only select ok, ok(10e5). No peak information. Charge information is referred from extracted data)
excel_file = R"G:\其他電腦\My Computer\GlycoMSParser\src\OG_int_online_annotation_revised_2024_glycan_zf_OG.xlsx"
xls = pd.ExcelFile(excel_file)
print(xls.sheet_names)
#for checking annotation file data integrity. Write another function for it
df = pd.read_excel(xls, sheet_name="MSlist")
print(df.head())
print(df.dtypes)

#find out CHARACTERS (caps only) have following digits  
pattern = re.compile(r'([A-Z]+)(\d+)')
patternex = re.compile(r'([A-Z]+)(\d*)')  #allow no digits
monosaccharide_masses = {
'F': 174.08921,
'H': 204.09977,
'N': 245.12632,
'S': 361.17367,
'G': 391.18423,
'KDN': 320.1465
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
'reduced': 62.0778,
}
#reduced 15(non-reducing end) + 15(CH3) + 14(CH2) + 18(H2O)

try:
    df1 = pd.read_excel(xls, sheet_name="fakesheet")
except Exception as e:
    print(f"no such a sheet exists {e}")

if "MSlist" in xls.sheet_names:
    df2 = pd.read_excel(xls, sheet_name="MSlist")
    print(df2.head())
else:
    print("No MSlist in sheet names")

if "fakesheet" in xls.sheet_names:
    df3 = pd.read_excel(xls, sheet_name="fakesheet")
    print(df2.head())
else:
    print("No fakesheet in sheet names")



def calculate_protonated_mass(row):
    annotation = str(row['Structure'])
    
    # check if 
    F_value = int(annotation.split('F')[1][0]) if 'F' in annotation else 0
    H_value = int(annotation.split('H')[1][0]) if 'H' in annotation else 0
    N_value = int(annotation.split('N')[1][0]) if 'N' in annotation else 0
    S_value = int(annotation.split('S')[1][0]) if 'S' in annotation else 0
    G_value = int(annotation.split('G')[1][0]) if 'G' in annotation else 0
    KDN_value = int(annotation.split('KDN')[1][0]) if 'KDN' in annotation else 0
    mass_shift = row['mass shift']
    #print(f"f value is {F_value} and it's type is {type(F_value)}")
    protonated_mass = int(F_value)*174 + int(H_value)*204 + int(N_value)*245 + int(S_value)*361 + int(G_value)*391 + int(KDN_value)*320 + 15 + 15 + 14 + 18+ 1 + mass_shift
    return protonated_mass


#c = df.iloc[2]
#print(c)
#print("test")
#d = calculate_protonated_mass(c)
#print(d)
# 使用apply方法計算新欄位
#df['protonatedmass'] = df.apply(calculate_protonated_mass, axis = 1)
#print("re-calculated")
#print(df.head())


def v2_calc_protonated_mass(row):
    annotation = row['Structure']
    #mass_shift = row['mass shift']   #no need for "calculating theoretic mass"
    #adduction = str(row['adduct'])
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
    #print(f"current total mass after glycan unit calculation is {total_mass}")
    #print(f"adducts is: {adduction} and the mass is {adducts_masses.get(adduction)}")
    #print(f"reduced glycan extra mass is {etc_masses.get('reduced')}")
    total_mass += 1.007325 + etc_masses.get('reduced')
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


#TODO: the adduct, charge and mass shift will be used to add : observed mass @done 20240813
 

#c = df.iloc[2]
#print(c)
#print("test")
#d = calc_observed_mass(c)
#print(d)
# 使用apply方法計算新欄位


def ppm_error1(observed_mass, protonated_mass, ppm):
    if ppm > (observed_mass - protonated_mass) / protonated_mass * 1e6:
        return True
    else:
        return False


df['theo protonated mass'] = df.apply(v2_calc_protonated_mass, axis = 1)
df['obs'] = df.apply(calc_observed_mass, axis = 1)
#print(f"calculated dataframe is {df.head()}")

#testmode
csv_file = R"C:\Users\Sakazuki\Downloads\ms2_zebrafish_Yann_glycome_20240731_zf_PerMeOG_brain.raw.csv" #R"/Users/hnstseng/Downloads/ms2_zebrafish_Yann_glycome_20240731_zf_PerMeOG_brain.raw.csv"#
#20240922 brain NG
#csv_file = R"G:\其他電腦\My Computer\GlycoMSParser\ms2_zebrafish_Yann_glycome_20240921_zf_sPerMeNG_brain_20240921_zf_sPerMeNG_brain.raw.csv" 
#20240922 ovary NG
#csv_file = R"G:\其他電腦\My Computer\GlycoMSParser\ms2_zebrafish_Yann_glycome_20240922_20240922_zf_sPerMeNG_ovary.raw.csv"
#20240922 intestine NG
#csv_file = "G:\其他電腦\My Computer\GlycoMSParser\ms2_zebrafish_Yann_glycome_20240922_zf_sPerMeNG_intestine.raw.csv"
#20241008 ovary OG
#csv_file = "G:\其他電腦\My Computer\GlycoMSParser\ms2_zebrafish_Yann_glycome_20241008_zf_sPerMeOG_ovary.raw.csv"
#csv_file = "G:\其他電腦\My Computer\GlycoMSParser\ms2_zebrafish_Yann_glycome_20241008_zf_sPerMeOG_intestine.raw.csv"
df2 = pd.read_csv(csv_file, sep='\t')
#print(df2.head())

matchlist = []

merged_df = pd.merge(df, df2, on="MS2scan_no", how="inner")

print("DEBUG in progress")
print(merged_df)
print(f"merged dafaframe by MS2scan_no{merged_df.columns}")
#print("store the merged file when debug is enabled")
#print("set a new dataframe for only needed information")
merged_df.to_csv("C:\Users\Sakazuki\Downloads\debugannotationstep1.csv", index=False)

def convert_to_float_list(s):
    return [float(x) for x in s.split(", ")]

pre_df = merged_df[["protonatedmass","peaklist","peakintensity","Structure", "IUPACname(optional)","Glycanannotation2"]]
#print(pre_df.head())
print(pre_df.info())
print("Load ion filter and convert the df to proper shape and then export it. Do it.")

ion_df = pd.read_excel(xls, sheet_name="ionlist")
ion_df = ion_df[["mass"]]
#print(ion_df.info())
#print(f"type of ion list mass is {type(ion_df)}")
# print(ion_df)


#####pre_df["peaklist"] = pre_df['peaklist'].apply(convert_to_float_list)
#####pre_df["peakintensity"] = pre_df['peakintensity'].apply(convert_to_float_list)
pre_df.loc[:, 'peaklist'] = pre_df['peaklist'].apply(lambda x: literal_eval(x.strip('()')))
pre_df.loc[:, 'peakintensity'] = pre_df['peakintensity'].apply(lambda x: literal_eval(x.strip('()')))

#warning#pre_df['peaklist'] = pre_df['peaklist'].apply(lambda x: literal_eval(x.strip('()')))
#warning#pre_df['peakintensity'] = pre_df['peakintensity'].apply(lambda x: literal_eval(x.strip('()')))

def intensity_normalization(ion, intensity):
    intensity = np.log10(intensity) + 1  #plus 1 shift after finished filling zero values
    return [ion, intensity]

def peaks_ppm(peaklist, peakintensity, reflist,setppm):
    print(f"The items in reference list is {len(reflist)}")
    hitlist = []
    for ref in reflist:
        #ref loop
        for peaks in peaklist:
            hitpeakinfo = [] #clear the value every time
            peakindex = 0 #reset index to zero, mind the index may change since we delete the range we've searched
            #peak loop
            if peaks <= ref:
                ppm = 1000000*(ref-peaks)/ref
                if ppm > setppm:
                    pass
                else:
                    #print(f"[dev]record peak. peak < ref")
                    peakindex = peaklist.index(peaks) 
                    hitpeakinfo = [ref, peakintensity[peakindex]]
                    hitlist.append(hitpeakinfo)
            elif peaks > ref:
                ppm = 1000000*(peaks-ref)/peaks
                if ppm < setppm:
                    #print(f"[dev]record peak. peak < ref")
                    peakindex = peaklist.index(peaks) 
                    hitpeakinfo = [ref, peakintensity[peakindex]]
                    hitlist.append(hitpeakinfo)
                else:
                    #print(f"[dev] peak >> ref. Delete peaklist until the index -1")
                    peakindex = peaklist.index(peaks)
                    #not sure if this works
                    peaklist = peaklist[peakindex:]
                    break #so the search start from next reference
    #print(f"[dev] hit list is {hitlist}")
    nhitlist = []
    for i in range(len(hitlist)):
        appendlist = intensity_normalization(hitlist[i][0], hitlist[i][1])
        nhitlist.append(appendlist)
    for logi in reflist:
        if logi not in nhitlist[:][0]:
            zeroions = [logi, 1]
            nhitlist.append(zeroions)
        else:
            pass

    #filling extra value +1 and fill blank with 0+1 to avoid division error
    #print(f"[dev] normalized hit list is {nhitlist}")
    return nhitlist






def findingions(dfrow, iondf, ppm):
    peaks = list(dfrow['peaklist'])
    intensity = list(dfrow['peakintensity'])
    #print(type(peaks))
    #print(f"peaks are {peaks[0:5]}")
    ions = []
    ions = iondf.values.flatten().tolist() #values in object
    #print(ions)
    #for i in range(5):
    #    print(f"ion type {type(ions[i])} and value is {ions[i]}")
    #    ion = float(ions[i])
    #    print(f"ion type {type(ion)} and value is {ion}")
    return peaks_ppm(peaks, intensity, ions,ppm)

iondfindex = ion_df.values.flatten().tolist()
iondfindex.extend(['Structure', 'IUPACname(optional)', 'Glycanannotation2'])
iondfindex.insert(0, 'protonatedmass')
niondf = pd.DataFrame(columns=iondfindex)
print(niondf.head())
print(f"length of pre_df = {len(pre_df)}")
#niondf['Structure'],niondf['IUPACname(optional)'],niondf['Glycanannotation2']  = '','', ''
#print(f"added annotation columns {niondf.head()}")
#print(niondf.iloc[2])
######
print(f"the size (entries) of original assigned list is {len(pre_df)} ")
for i in range(len(pre_df)):#range(5):#
    print(f"index is {i}")
    #init a series of current datarow. it's always throwing me KeyError: 0 or IndexError: single positional indexer is out-of-bounds
    annotationlist = [['protonatedmass',pre_df.iloc[i]['protonatedmass']],['Structure',pre_df.iloc[i]['Structure']], ['IUPACname(optional)' ,pre_df.iloc[i]['IUPACname(optional)']],['Glycanannotation2',pre_df.iloc[i]['Glycanannotation2']]]
    tempions = findingions(pre_df.iloc[i], ion_df, 10)
    data_dict = {item[0]: item[1] for item in tempions + annotationlist}
    new_row = pd.Series(data_dict)
    #print(f"tempions is {tempions}")
    #addseries = pd.Series(data = tempions, index = iondfindex)
    #print(f"annotation info from {i} entry of preprocessed df {annotationlist}")
    #print(f"Series to append:{new_row} ")
    niondf = pd.concat([niondf, new_row.to_frame().T],ignore_index=True)
    #print(f"normalized df {niondf.iloc[i]['Structure']}  and source df {pre_df.iloc[i]['Structure']}")
    #niondf.iloc[i]['Structure'] =  pre_df.iloc[i]['Structure']
    #niondf.loc[i]['IUPACname(optional)'] =  pre_df.iloc[i]['IUPACname(optional)']
    #niondf.loc[i]['Glycanannotation2'] =  pre_df.iloc[i]['Glycanannotation2']
#    tempions = findingions(pre_df.iloc[i], ion_df, 10)
    #for addions, ionvalue in enumerate(tempions):
    #    for col, value in ionvalue:
    #        niondf.iloc[i, col] = value
    # 将未被填充的列填为0
    #niondf.iloc[i, niondf.loc[i].isna()] = 0
#will need to fill not found cells with 0 and apply extra plus 1 to all values
#which one is faster? give 1 when searching? or add 1 to all values?
#for i in range(iondfindex):
#    df[i] +=1
#print(f"modified df: {niondf.head()}")
#peaks_ppm(pre_df.iloc[2], peakintensity, reflist,setppm)
#print(f"now niondf is {niondf.head()}")


#niondf.to_csv("zfOGbrain_converted.csv", index=False)
#20240922
#niondf.to_csv("zfNGbrain_converted.csv", index=False)
#20240922_2
#niondf.to_csv("zfNGovary_converted.csv", index=False)
#20240922_3
#niondf.to_csv("zfNGintestine_converted.csv", index=False)
#20241008
#niondf.to_csv("zfOGovary_converted.csv", index=False)
niondf.to_csv("zfOGintestine_converted.csv", index=False)
#df3 = (df2[df2["MS2scan_no"].isin(df["MS2scan_no"])] | (ppm_error(df2["protonatedmass"],df["protonated mass"],10)))
#print(f"df3 is {df3.head()}")

def ppm_error2(row, df, ppm):
    matching_rows = df.loc[df['MS2scan_no'] == row['MS2scan_no'], 'theo protonated mass']
    if matching_rows.empty:
        return False  # 如果找不到匹配项，则跳过这一行
    protonated_mass = matching_rows.values[0]
    return ppm > abs(row['protonatedmass'] - protonated_mass) / protonated_mass * 1e6
# 使用 apply 方法在 df2 的每一行上调用 ppm_error 函数
#df3 = df2[df2.apply(lambda row: ppm_error2(row, df, 10), axis=1)]
#df4 = df3.drop(columns=['entry_no', 'MS1_isolationmass', 'MS1_monoisolationmass', 'chargeState', ], inplace=True)
#df4 = df[df.apply(lambda row: ppm_error(row, df2, 10), axis=1)]
#print(f'changed df3 is \n {df3.head()}')