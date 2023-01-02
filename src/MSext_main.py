from pymsfilereader import MSFileReader
import glob
import pandas as pd
import os
#import needed file


#load file  -- it has been moved to DetermineMSLevel.py #####
'''
try:
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
    Tk().withdraw()
    filename = askopenfilename
    print('Selected', filename)
    rawfile = MSFileReader(filename)
except ImportError:
    raise ImportError('Please install tkinter')

'''
#############################################################
try:
    import testfunc
except ImportError:
    raise ImportError('testfunc is not found')
#no need to test it now
#print(testfunc.test1())

print('testfunc loaded')

try:
    import DetermineMSLevel
except ImportError:
    raise ImportError('testfunc is not found')
#test for the findingMSlevel function
#print(DetermineMSLevel.findingMSlevel())
print('DetermineMSLevel finished')


#This part will load zf_sPerMeOG_intestine specifically in Downloads
'''
print('main')
try:
    import DetermineMSLeveltest
except ImportError:
    raise ImportError('DetermineMSLevel is not found')
print('test 1 finished')
'''
#####################################################################



#menu
#considering using tkinter to build GUI for cross-platform
#current lib is windows-dependent, try to find out how thermo nupkg works

#for testing function#
print('self test')
#print(testfunc.test1())
print(testfunc.platformtest())

print('###################\n#   Menu     #\n###################')
print('select which you want to do.\n 1. Check MSlevel 2. Extract MS1 data')
c = input('')

if c == '1':
    print('check MSLevel')
    print(DetermineMSLevel.findingMSlevel())

elif c == '2':
    print('test calling ExtractMS')
    print(DetermineMSLevel.ExtractMSlevel())
elif c == '3':
    print('go to next test section, grabbing tsv')
    pass
else:
    print('don\'t enter other number.')
print 


###function for pandas

###read one file and check the way to handle df properly.
#keep only charge/ m/z in MS1 with headers (put at where)
#and MS2 peaklists

### more aggressive way: delete all of unneeded MS1 and MS2 don't own useful ions
### how to? check the code difference between Boston univ's and glypick

'''
print('test gathering data from pandas')
path = r'C:\Users\hnstseng\Downloads'

all_files = glob.glob(os.path.join(path , "/*.tsv"))

li = []

for filename in all_files:
    df = pd.read_csv(filename,sep='\t')
    li.append(df)

print(li[0])
'''
