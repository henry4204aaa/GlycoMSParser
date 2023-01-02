import glob
import pandas as pd
import os
import pathlib
#test w/o saving spectrum to csv
#from DetermineMSLevel import ExtractMSlevelauto 


#consider loading csv from tmp
ms1 = "no23043isMS1.tsv"

'''read one csv'''
dfms1 = pd.dataframe(pd.read_csv(ms1), sep='\t')

#Unicode error in windows
#dfs = {f.stem: pd.read_csv(f) for f in pathlib.Path('.NewFolder').glob('*.csv')}





#monitor memory use to prevent crash


#read first set of MS1+2 tables and clean

#MS1
#keep No [0][0],  mass[0][1], intensity? [0][2], charge [0][6]

#pandas add empty column for labeling

#Ms2
#read header and keep only precursor and...? monoisotopic mass

#think about the execution time of calling scanhearinfo and extrac only those you need
#fast way to drop unneeded data:
#for i in peaklist (find any peaks in wishlist <- Glypick may work in this logic)


#relate intensity% isn't export directlly, need to calculate before drop? Needed?
#partial maximum (largest peak in list = 1.00
#global maximum (while the maximum peak may be dropped)



#copied from stackovrrflow, not working

#
'''
print('test gathering data from pandas')




###read multiple files###
path = r'C:\Users\hnstseng\Downloads' 
all_files = glob.glob(os.path.join(path , "*.tsv"))

li = []

print(all_files)

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    li.append(df)

frame = pd.concat(li, axis=0, ignore_index=True)
##############################
'''
