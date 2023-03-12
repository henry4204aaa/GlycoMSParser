import sys

#global use#
gui_available = None
gui_autostart = False
msglyco = None



#start argument settings#
def startparam(gui):
    global gui_autostart
    if (str(gui) == "-GUI") or (str(gui) == "-g"):
        gui_autostart = True
    else:
        gui_autostart = False
if __name__ == "__main__":
    if len(sys.argv) == 2:
        startparam(str(sys.argv[1]))
    else:
        pass

#when missing package, throw a warning
#THIS DOESNT STOP THE CODE, please set other flags to prevent errors
def missingpackage(name):
    def init(self, *args, **kwargs):
        raise ImportError(
            f'Missing {name} package so the program cannot be executed.\nPlease install it and run this tool again.'
        )
    return type(name, (), {'__init__':init})


#import needed core components

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
    msglyco = True
except ImportError:
    print('You need pymsfilereader installed to extract raw files')
    msglyco = False


#call rawextractor.py
def rawextractor():
    if msglyco is True:
        try:
            import trytolistoutheadersv2
        except ImportError:
            f"Missing core files 'trytolistoutheaders' to extract the spectra. Package corrupted?"
        try:
            import numpy
        except ImportError:
            globerr = missingpackage('numpy')
        #print("Raw extractor function part")
    else:
        print(f"Glycan related extrator function has been disabled since you didn't install pymsfilereader.")

#call 
def filterms2():
    #print("FilterMS2 part")
    pass

#call
def annotationhelper():
    #print("annotation helper part")
    pass

#call
def miscellaneous():
    #print("Placeholder for other developing functions")
    pass



def GUIinit(gui_autostart=False):
    print("GUI init")
    try:
        import tkinter as tk
        global gui_available
        gui_available = True
        if gui_autostart is True:
            root = tk.Tk()
            root.title('GlycoMSParser GUI version 0.1')
            root.geometry('640x360')
            devinfo = tk.Label(root, text="Under development, please close this window", font=('Arial', 24))
            devinfo.pack
            root.mainloop()
            return False
            main_loop()
        else:
            main_loop()
    except ImportError:
        print('Install tkinter to activate GUI function or ')
        gui_autostart = False
        gui_available = False
        main_loop()



def callfunc(n):
    if n in function_dlist:
        if n.lower() == "q":
            print(function_dlist['q'])
            return False
        elif n.lower() == "g":
            print('g')
            GUIinit(gui_autostart = True)
            #when exit it will run again (that's question)
        elif n.isdigit():
            if int(n) == 1:
                print('call RAWheaderextractor - still testing')
                print('this function also requires numpy in beta version')
                function_dlist[n]     
            elif int(n) == 2:
                print('call function2')
                function_dlist[n]
            elif int(n) == 3:
                print('call function3')
                function_dlist[n]
            elif int(n) == 4:
                print('call function4')
                function_dlist[n]
        else:
            print("Something error occurred")
    else:
        print("Invaild command")
        return





#build up essentials for the mainloop
def main_loop():
    loopinit = True
    Datahandle = False
    #need a loop here
    while loopinit:
        global gui_autostart
        while gui_autostart is True:
            GUIinit(gui_autostart=True)
            break
        print("-"*20)
        print("GlycoMSParser CLI v0.1")
        print(f"GUI start= {gui_autostart}")
        print(f"GUI status None is not checked= {gui_available}")
        print(f"Glycoextractor availibility = {msglyco}")
        print("-"*20)
        print("Mode selection loop")
        print("select the functions you want to execute\n")
        print("1 --- Rawextractor \n2 --- MS2 Filtering \n3 --- Annotation helper\n4 --- others\nq --- quit \ng start GUI")
        n = input("")
        if callfunc(n) == False:
            break




function_dlist = {'1': rawextractor(), '2': filterms2(), '3': annotationhelper(), '4': "miscellaneous in DEV", 'q': "quit", 'g': 'GUI'}



main_loop()

#code structure
#set global
#set arguments when starting the app
#define actions when missing packages
#define CLI calling functions
#define CLI itself
#main block
#try to import needed modules
#call tkinter if possible/GUI is set to true

