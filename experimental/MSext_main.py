#global variables
class Flagcheck():
    global pyms_activated, gui_start, gui_support
    pyms_activated = False
    gui_start = False
    gui_support = False

#####init block####

print("#" * 20, "Initializing", "#" * 20)
import sys

try:
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
    print('tkinter initialized')
    gui_support = True
except ImportError:
    raise ImportError('Please install tkinter')

def loadrawextractor_pyms():
    try:
        from pymsfilereader import MSFileReader
        pyms_activated = True
    except ImportError:
        raise ImportError('Please install pymsfilereader')

def loadrawextractor_myextractor():
    #Build a new one that's universal for all platform...
    print('Will have it in future')


def load_generalpackage():
    try:
        import pandas as pd
        import glob
        import os
    except ImportError as errors:
        '''
        # From https://blog.airbrake.io/blog/python/importerror-and-modulenotfounderror
        # Output expected ImportErrors.
        Logging.log_exception(error)
        # Include the name and path attributes in output.
        Logging.log(f'error.name: {error.name}')
        Logging.log(f'error.path: {error.path}')
        except Exception as exception:
        # Output unexpected Exceptions.
        Logging.log_exception(exception, False)
        '''
        raise ImportError("some of the general packages aren't loaded")

def printimportedmodules():  
    print("you need to have pymsfilereader, pandas and tkinter installed")
    for module_name in sys.modules:
        module = sys.modules[module_name]
        print(module_name, getattr(module, '__version__', 'n/a'))


#############################################################
#Import GlycoMSParser dependencies
def loadFilterMS2():
    try:
        import FilterMS2
    except ImportError:
        raise ImportError('testfunc is not found')
    #test for the findingMSlevel function
    #print(DetermineMSLevel.findingMSlevel())
    print('FilterMS2 loaded')

def loadcomparepeaklist():
    try:
        import comparepeaklist
    except ImportError:
        raise ImportError('comparepeaklist is not found')
    #test for the findingMSlevel function
    #print(DetermineMSLevel.findingMSlevel())
    print('comparepeaklist loaded')

#####################################################################



#menu
#considering using tkinter to build GUI for cross-platform
#current lib is windows-dependent, try to find out how thermo nupkg works
Menu = True
while Menu:
    print('###################\n#   Menu     #\n###################')
    print('select which you want to do.')
    menus = ["Load RAWfileextractor", "Load FilterMS2", "Load annotationhelper", "Load testfiles", "Loadlibraries", "Load all files"]
    for i in range(len(menus)):
        a = str(i+1)+ "."
        print(a, menus[i])
    print("0. for printimportedmodules")
    print('###################')
    c = input('Select the number')
    if c == '1':
        print('load functions')
        import trytolistoutheaders
    elif c == '2':
        loadFilterMS2()
    elif c == '3':
        import annotationhelper
        pass    
    elif c == '4':
        print('import comparepeaklist')
        import comparepeaklist
    elif c == '5':
        print('import library')
        load_generalpackage()
    elif c == '6':
        print("load all functions but how")
        import trytolistoutheaders
        loadFilterMS2()
        import annotationhelper
        import comparepeaklist
    elif c == '0':
        printimportedmodules()
    else:
        print('don\'t enter other number.')
        break


