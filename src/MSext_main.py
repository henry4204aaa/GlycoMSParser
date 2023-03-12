'''
from pymsfilereader import MSFileReader
import glob
import pandas as pd
import os
'''
#import needed file


#2023/2/15

#load file and run tkinter, if no then enter command line interface
def missingpackage(name):
    def init(self, *args, **kwargs):
        raise ImportError(
            f'Missing {name} package so the program cannot be executed\n'
            f'Please install it and run this tool again.'
        )
    return type(name, (), {'__init__':init})

try:
    import glob
except ImportError:
    globerr = missingpackage('glob')
try:
    import os
except ImportError:
    globerr = missingpackage('os')
try:
    import pandas as pd
except ImportError:
    globerr = missingpackage('pandas')


try:
    from pymsfilereader import MSFileReader
except ImportError:
    print('You need pymsfilereader installed to extract raw files')
    msglyco = False

try:
    import tkinter as tk
    root = tk.Tk()
    root.title('GlycoMSParser GUI version 0.1')
    root.geometry('640x360')
    root.gameloop()
except ImportError:
    print('Install tkinter to activate GUI function or ')
    gui_available = False
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
