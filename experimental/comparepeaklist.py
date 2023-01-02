#for test: defined list input
#in real: the list determined prior to each run
#leave logs for tracking what users have defined

#init
#move this to init.py
try:
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
    Tk().withdraw()
    filename = askopenfilename()
    print('Selected', filename)
    
except ImportError:
    raise ImportError('Please install tkinter')

inputlist = []
comparelist = []
inputintensitylist = []
ninputintensitylist = []
#force int type after validate it's not NoneType
def ppmcalculation(userdefinedppm=None):
    defineppmflag = True
    #if the ppm is not set or not an integer, this loop will ask you continuously
    while (userdefinedppm is None) or (defineppmflag is True):
        userdefinedppm = input('please enter the ppm between 0 and 100000')
        try:
            userdefinedppm = int(userdefinedppm)
            if ((0 < userdefinedppm) and (userdefinedppm < 100000 ) is True):
                defineppmflag = False
            else:
                print('please enter the correct ppm between 0 and 100000')
        except:
            print('please enter an integer')
    return userdefinedppm


def comparepeaklist(inputlist, comparelist):
    #pseudocode
    #minimize calculation
    #in defined list
    #evaluate the efficiency of which loop will run faster
    #though I think it depends on the data (no of i should >>j)
    #for i in inputlist:
    #    for j in comparelist:
    #            if j < i:
    #                check next j (pass or continue)
    #            elif j == i:  (< defined delta ppm)
    #                hitcount +=1
    #                break
    #            elif j > i:
    #                print('no this ion found, move next')
    #                 break
    #            else:
    hitcount = 0
    hitlist =[]
    for i in inputlist:
        #debug
        #print('input is',i)
        for j in comparelist:
            #print('comparing with', j)
            if i > j:
                #print('input > comparepeak, go up')
                del comparelist[0]
                #continue
            elif i == j:
                #print('should be the ppm calculation')
                ##### pseudocode
                # compare input[i] as a and comparelist[j] as b
                # if (abs((a-b)/b)*1000000) < userdefinedppm:
                #     add the input[i] into hit list
                ## using abs will make comparing input and compare unavailable
                ## set?
                hitcount +=1
                hitlist.append(i)
                del comparelist[0]
            elif i < j:
                #print('input < comparepeak')
                break #byebye, check next i
            else:
                print('exception')
            
    #print('hitcount=', hitcount)
    #print('ions, ', hitlist)
    return hitcount, hitlist
#Think about glee senpai said if built-in set datatype
#can run faster (not sure if the numbers will be compared
#as string or float number

def comparepeaklistppm(inputlist, comparelist, ppm, inputintensitylist):
    #pre-assign input list and compare list
    if (inputlist == []):
        inputlist = [344, 344.16, 344.19, 344.25, 345, 376, 432, 464]
        print('use default inputlist', inputlist)
    else:
        #print('inputlist is', inputlist)
        pass
    if (comparelist == []):
        comparelist = [344.17, 374.18, 376.20, 406.21, 432.22, 450.23, 464.25, 580.30, 793.38, 825.42]
        #print('use default comparelist', comparelist)
    else:
        #print('comparelist is', comparelist)
        pass
    if (inputintensitylist == []):
        inputintensitylist = [500, 1000, 34567, 43210, 666666, 1, 810, 114514, 100000, 100000000000]
        #print('use default comparelist', comparelist)
    else:
        #print('comparelist is', comparelist)
        pass
    #ppm = ppmcalculation()
    #normalize from inputintensitylist
    #pseudocode:
    #find largest value in inputintensitylist = maxin
    #ninputintensitylist = (inputintensitylist/maxin)*100 for i in range(inputintensitylist)
    maxin = max(inputintensitylist[:])
    #print('max', maxin)
    ninputintensitylist = []
    for i in range(len(inputintensitylist)):
        ninputintensitylist.append((inputintensitylist[i]/maxin)*100)
    #for test the normalized bug (fixed)
    '''
    if str(maxin)== str(3798419.75):
        print('hey', ninputintensitylist)
    else:
        pass
    '''
    hitcount = 0
    hitlist =[]
    location = []
    hitintensitylist = []
    nhitintensitylist = []
    for i in inputlist:
        #print('input', i)
        for j in comparelist:
            #print('comparing with', j)
            #find out which is larger
            if i > j:
                k = (i-j)/i
                ioverj = True
            elif i < j:
                k = (j-i)/j
                ioverj = False
            elif i == j:
                #test why 376 isn't appearing
                '''
                if i == 376.20:
                    print(i, 'found')
                ioverj = False
                '''  
                k = 0
            else:
                print('error in comparing i j')

            if (k * 1000000 < ppm):
                #location for finding intensity info
                loc = inputlist.index(i)
                location.append(loc)
                htmp = inputintensitylist[loc]
                ntmp = ninputintensitylist[loc]
                hitintensitylist.append(htmp)
                nhitintensitylist.append(ntmp)
                hitcount +=1
                hitlist.append(i)
                #print('it\'s hit')
                #print("now going to delete", comparelist[comparelist.index(j)])
                del comparelist[comparelist.index(j)]
                break
            elif (ioverj is True) and (k * 1000000 > ppm):
                #print('input > comparepeak, go up')
                #delete those never been compared in the coming iteration
                #print("now going to delete", comparelist[comparelist.index(j)])
                del comparelist[comparelist.index(j)]
                #print('delete those comparing ions smaller than input', comparelist)
            elif (ioverj is False) and (k *1000000 > ppm) :
                #print('input < comparepeak')
                break #byebye, check next i
            else:
                print('exception')
    return hitcount, hitlist, location, hitintensitylist, nhitintensitylist #, maxin
#Think about glee senpai said if built-in set datatype
#can run faster (not sure if the numbers will be compared
#as string or float number
def calcproton(isolatedmass, charge):
    if abs(charge) == 1 :
        return isolatedmass
    elif charge >= 2:
        conv = (isolatedmass * charge - (charge - 1) * 1.00784)
        return conv
    elif charge == 0:
        return 0
    elif charge < -1:
        negconv = (isolatedmass * charge + (charge - 1) * 1.00784)
        return negconv
    #### pseudocode
    #convert datatype into float and int
    # if charge == 1:
    # return isolatedmass
    # elif charge == 2:
    # [m+2H]/2 -> m+H
    # return (2*isolatedmass - H)
    # elif charge == 3:
    # return (3*isolateass - 2 * H)
    #... (seems the code can be condensed)
    # elif charge == 0 or None:
    # return 0    (confusing C people?)
    
    

