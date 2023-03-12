#thisispseudocode for extracting roughly extracted csv (merge it to listoutheadersv2 later)
import os
import glob




#func pickheader
def checkheaderlist(checklist):
    #get csv header
    #compare csv header with currently existed
    #keep those matched column
    #make a list or dict to allow user to add preset or load presets, and add any cooumn they want
    pass
    #consider the integration in the future, we may need to build a id array to make the analysis pipeline clearer
    #return pickedlist

#load csv file
def loadcsvfile(defaultcsv=None):
    #load or list out available csv and show the selected csv and ask user before start picking
    #this may be a menu or outside the function, chexk and see which runs faster
    path = "./"
    ext = "*.csv"
    csvlist = list(filter(lambda x: '.csv' in x, os.listdir(path)))
    #print(csvlist)   #this print out all csvs in this directory
    contain_def = "MS2 summary"
    defcsvformat = []
    num_csvlist=0
    for i in csvlist:
        if contain_def in i:
            num_csvlist+=1
            defcsvformat.append(i)
            print(f"no. {num_csvlist}:{i} infoblahblahblah")
    #print(defcsvformat)
    while True:#addind a file selection loop
        selcsv = input("The csv you're going to type from pre-processed csv")
        print("This is file selection loop")
        if selcsv.isdigit():
            try:
                parsecsv = defcsvformat[int(selcsv)-1]
                print(f"You selected {parsecsv}")
                return parsecsv
                break
            except IndexError:
                print("You're asking a file out of range")
        else:
            if selcsv in defcsvformat:
                parsecsv = selcsv
                print(f"You selected {parsecsv}")
                return parsecsv
                break
            elif selcsv in csvlist:
                print(f"{selcsv} is not in standard filename. It may not work in next processing. Select it anyway?")
                usrans = input("y/n")
                if usrans.lower() == "y":
                    parsecsv = selcsv
                    print(f"You selected {parsecsv}")
                    return parsecsv
                    break
                else:
                    print("Restart the selection")
            else:
                usrans = input("You're not giving any proper csv file, enter q to quit or re-select other one")
                if usrans.lower() == "q":
                    print("exit selection, nothing has been processed")
                    break



def extractpickeddata(ext):
    #extract with defined paramteters and export the file
    #should we have headers or unique id label in filename (cab be decoded when parsing filename)
    #shared list or dict may also be suitable
    pass
