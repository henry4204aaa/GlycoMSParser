def structuretodictparser():
    #should be selectable in future
    newdf = pd.read_csv("revised_zfNGbrain_20221129_1_updatedlist.csv", sep='\t', index_col=[0])
    exportdf = []
    #handling the data in a df
    newdf.insert(8, 'Structure dict',"")
    for i in range(len(newdf)):
        parser = newdf.iloc[i]['Structure']
        strucdict = {}
        strucdictcom = pd.DataFrame()
        tmpseries = ()
        #start from F
        fuc = re.search('F\d+', parser)
        hexose = re.search('H\d+', parser)
        hexnac = re.search('N\d+', parser)
        neu5ac = re.search('S\d+', parser)
        neu5gc = re.search('G\d+', parser)
        #ambigious strings, get nothing if not present
        Naadduct = re.search('Na', parser)   # if not Naadduct: do nothing ; else add an adduct
        sus = re.search('\?', parser)
        #not included in the first test data
        kdn = re.search('Kdn\d+', parser)
        NH3adduct = re.search('NH3', parser) #may be removed if we calculate this in other function
        Kadduct = re.search('K', parser)
        #conbine the  info into single list in each loop
        if not fuc:
            strucdict.update({'F': 0})
        else:
            strucdict.update({'F': (fuc.group()[1:])})
        if not hexose:
            strucdict.update({'H': 0})
        else:
            strucdict.update({'H': (hexose.group()[1:])})
        if not hexnac:
            strucdict.update({'N': 0})
        else:
            strucdict.update({'N': (hexnac.group()[1:])})
        if not neu5ac:
            strucdict.update({'S': 0})
        else:
            strucdict.update({'S': (neu5ac.group()[1:])})
        if not neu5gc:
            strucdict.update({'G': 0})
        else:
            strucdict.update({'G': (neu5gc.group()[1:])})
        if not kdn:
            strucdict.update({'KDN': 0})
        else:
            strucdict.update({'KDN': (kdn.group()[1:])})
        # I won't add them as default since they are exceptions.
        if not Naadduct:
            pass
            #strucdict.update({'Na': 0) if the loop on the list throws an error, try this
        else:
            strucdict.update({'Na': 1})
        if not sus:
            pass
            #strucdict.update({'Questionable': 0) if the loop on the list throws an error, try this
        else:
            strucdict.update({'Questionable': 1})
        tmpseries = newdf.iloc[i][1:]
        tmpseries['Structure dict'] = strucdict
        exportdf.append(tmpseries)
    exportdf = pd.DataFrame(exportdf)
    exportdf = exportdf.reset_index(drop=True)
    print(exportdf)
    exportdf.to_csv("revised_zfNGbrain_20221129_1_updatedlist_parsed.csv", sep='\t')
    #replace the iloc[i]['Structure']  <- SettingWithCopyWarning and the df is not modified.
    
#it's not correct
  
'''
    if parser.find('F') == -1 :
        strucdict['F'] = 0
    else:
        strucdict['F'] = parser[(parser.index('F')+1)]
    print strucdict
'''



def deftheoreticalms1():
    #should be an element class holding these values
    df = pd.read_csv("revised_zfNGbrain_20221129_1_updatedlist_parsed.csv", sep='\t', index_col=[0])
    #askopenfile?
    theomasslist = []
    newdf = pd.DataFrame()
    #print(df.dtypes)
    getcontext().prec = 8
    H = Decimal(1.007825) 
    C = Decimal(12.000)
    O = Decimal(15.9994)
    #N = 14
    Nadif = Decimal(22.989769) - Decimal(1.007825)   #replace H with Na
    mperMedHex = Decimal(174.08921)
    mperMeHex = Decimal(204.09977)
    mperMeHexNAc = Decimal(245.12632)
    mperMeNeu5Ac = Decimal(361.17367)
    mperMeNeu5Gc = Decimal(391.18423)
    mperMefreeend = Decimal(46.04186)  #CH3COH3
    mperMeterminal = Decimal(15.023475) #CH3, in positive mode (H+)
    #read the df
    #for testing: obsmass = {'F': 1, 'H':3, 'N': 2, 'S': 0, 'G': 0 }
    #get the structure dict
    for i in range(len(df)):
        tmp = df.iloc[i]['Structure dict']
        obsmass = ast.literal_eval(tmp)
        theomass = 0
        #calculate it
        #print(obsmass, obsmass['F'], type(obsmass['F']))
        Fuc = Decimal(obsmass['F'])
        Hex = Decimal(obsmass['H'])
        HexNAc = Decimal(obsmass['N'])
        Ac = Decimal(obsmass['S'])
        Gc = Decimal(obsmass['G'])
        keys = list(obsmass)
        #print('keys', keys[-1])
        if (keys[-1] == 'Na') or (keys[-2] == 'Na'):
            #print('find Na')
            sodium = 1
        else:
            sodium = 0
        if(keys[-1] == 'Questionable'):
            print('it has a mass difference from theoretical value')
        theomass = str((Fuc * mperMedHex) + (Hex * mperMeHex) + (HexNAc * mperMeHexNAc) + (Ac * mperMeNeu5Ac) + (Gc * mperMeNeu5Gc) + mperMefreeend + H + (Nadif * sodium))
        theomasslist.append(theomass)
    print(theomasslist)
    df.insert(5, 'Theoretical mass' , theomasslist, allow_duplicates = True)
    print(df)
    df.to_csv('Annotated_20221130.csv', sep='\t', index= False)
    print('finish writing csv')
    #add a new column and insert