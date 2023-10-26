def comparepeaklistppm(inputlist, comparelist, ppm, inputintensitylist):
    #pre-assign input list and compare list
    if (inputlist == []):
        inputlist = [344, 344.16, 344.19, 344.25, 345, 376, 432, 464]
        print('use default inputlist', inputlist)
    else:
        print('inputlist is', inputlist)
        pass
    if (comparelist == []):
        comparelist = [344.17, 374.18, 376.20, 406.21, 432.22, 450.23, 464.25, 580.30, 793.38, 825.42]
        print('use default comparelist', comparelist)
    else:
        print('comparelist is', comparelist)
        pass
    if (inputintensitylist == []):
        inputintensitylist = [500, 1000, 34567, 43210, 666666, 1, 810, 114514, 100000, 100000000000]
        #print('use default comparelist', comparelist)
    else:
        print('inputintensitylist is', comparelist)
        pass
    a = "testing function"
    return a