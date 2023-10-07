
def predcomplsttodict(complist):
    predcomp = {}
    if len(complist) >= 5:
        predcomp["F"] = complist[0]
        predcomp["H"] = complist[1]
        predcomp["N"] = complist[2]
        predcomp["S"] = complist[3]
        predcomp["G"] = complist[4]
        if len(complist) == 6:
            predcomp["KDN"] = complist[5]
    print(f"predicted composition dict is {predcomp}")
    return(predcomp)

def checkifpossiblededuce(invalflag): #for debugging
    if (invalflag is True):
        print("This structural enumeration is not possible")
    elif(invalflag is False):
        print("Applying calculations...")
    else:
        raise ValueError("you haven't check if the enumeration is possible")
        
def compactcompcalc(calccomp, enumerateepitope, looplevel=0): #separate the elumarating METHOD outside as a independent function
        invalflag = None #init
        remaining_unit = [] #flag
        remaining_log = [] #Logs
        print(f"iteration on {enumerateepitope[0]}, nested loop depth= {looplevel}")
        valdict = {units: calccomp[units] - enumerateepitope[1].get(units, 0) for units in calccomp}
        #print(f"current valdict in compactcompcalc {valdict}")
        for unittest in valdict.values():
            if unittest < 0:
                invalflag = True
                #print("caught negative value")
                remaining_unit.append(False)
            else:
                remaining_unit.append(True)
                pass
        if (False in remaining_unit) is False:  #why this statement works?
            invalflag = False
        #checkifpossiblededuce(invalflag)  #for debug
        if invalflag is False:
            #print("add to logs")
            remaining_log.append(enumerateepitope[0]) #need to combine
            return True, remaining_log, valdict #added valdict
        elif invalflag is True:
            #print("do nothing. Test next glycotope")
            return False, "", valdict
        else:
            raise ValueError("The core assignment is somhow not executed.")

def unittest(valdict):
    remaining_unit = []
    invalflag = None
    for unittest in valdict.values():
        if unittest < 0:
            invalflag = True
            #print("caught negative value")
            remaining_unit.append(False)
        else:
            remaining_unit.append(True)
            pass
    if invalflag:  
        return False 
    else:
        #print("no negative values found")
        return True
            
def testzero(dictinput):
    for unittest in dictinput.values():
        if (unittest != 0): #set it invalid and leave no logs recorded
            return False
    return True

def NGcorearms(spectrapredcomp):
    corededuce = {'H': 3, 'N': 2}
    coreelements = {units: spectrapredcomp[units] - corededuce.get(units, 0) for units in spectrapredcomp} #for further calculation
    #try 4 arms
    biante = {'H': 2, 'N': 2}
    triante = {'H': 3, 'N': 3}
    tetraante = {'H': 4, 'N': 4}
    bisect = {'H': 2, 'N': 3}
    hybrid = {'H': 7, 'N': 3}
    tmplist = [biante,triante,tetraante,bisect,hybrid]
    possiblearms = []
    i = 0
    #print(f"after subtract corededuce = {coreelements}")
    for arms in tmplist:
        arm = {units: coreelements[units] - arms.get(units, 0) for units in coreelements}
        #print(f"i = {i} and current calc is {arm}")
        #print(f"unit test on arm :{unittest(arm)}")
        if unittest(arm):
            if i == 0: #"biante":
                possiblearms.append(2)
                print(f"now possible arms is biante")
            elif i == 1: #arms == "triante":
                possiblearms.append(3)
                print(f"now possible arms is triante")
            elif i == 2: #arms == "tetraante":
                possiblearms.append(4)
                print(f"now possible arms is tetraante")
            elif i == 3: #arms == "bisect":
                possiblearms.append(2)
                print(f"now possible arms is bisect")
            elif i == 4: #arms == "hybrid":
                possiblearms.append(1)
                print(f"now possible arms is hybrid")
        else:
            possiblearms.append(-1) #0 -> -1 to avoid conflicts
        i+=1
    #print(possiblearms)
    if possiblearms is []:
        possiblearms.append(-2) #-2 means error
    else:
        pass
    #print(f"possible arms {possiblearms}")
    #arms = 0
    arms = max(possiblearms)
    #print(f"max possible arms is {arms} in function")
    #consider returning all arms to reduce complexity in deducer
    return arms, coreelements 

#testing block
#print("testing NGcore function")    
#testcoomp = {'H': 6, 'N': 5}
#NGcorearms(testcoomp)
#print("Finished testing. Need to pass arms and composition")
#print("testing NGcore function2")    
#testcoomp1 = {'H': 4, 'N': 4}
#NGcorearms(testcoomp1)
#print("Finished testing. Need to pass arms and composition")
#print("testing NGcore function3")    
#testcoomp2 = {'H': 2, 'N': 2}
#NGcorearms(testcoomp2)
#testing block ends

def deducer(specific_epitopes, spectrapredcomp):
    #print("Try to assign NG core decreasing")
    calc_comp = spectrapredcomp #store to another dict to avoid updating on original data needed later

    #print("importing grouped epitope definition...")
    #group_NGcore = sorted_coreelements  #first level
    group_terminal = sorted_terelements #second level
    group_extensive = sorted_extelements #final level
    #how to calculate sugar unit changes: #cdict = {key: calc[key] - bdict.get(key, 0) for key in adict}
    
    #print("adding log registering space")
    complog = [] #when to collect and clean complog?
    finalresults = []
    armsno, forcalcomp= NGcorearms(calc_comp)
    #print(f"Testing NGcore decuce and the current maxiumn arm is {armsno} and comp for further deduction is {forcalcomp}.")
    ###core###
    for armdeduce in corearms:
        #get arms number first to reducde unneeded calculation...
        #arms= -1 means the composition is less than H5N4
        #currently corrupted NG (such as H4N3+sth) isn't calculated correctly
        #maybe we still need corrupted units enumeration + terminal if possible... but how? 
        if armsno == -1: #means only NG core may be possible
            #print("confirm if the current dict left zero")
            if unittest(forcalcomp): #check negative value first, and then start terminal addition
                loopoverflag = True
                compforloop = forcalcomp
                flag1111 = 0
                flag2222 = 0
                while loopoverflag: #loops until all units get zeroed or any of negative value appears
                    #debug
                    #print(f"Current comp {compforloop} and first input comp {forcalcomp}")
                    #print(f"debug flag 1 = {flag1111}")
                    flag1111+=1
                    #debug ends
                    for leftovers in cor_epitopes:
                        #print(f"debug flag 2 = {flag2222}")
                        flag2222+=1
                        #print("start checking terminal addition")
                        compforloop1 = {units: compforloop[units] - leftovers[1].get(units, 0) for units in compforloop}
                        compforloop.update(compforloop1)
                        if unittest(compforloop):
                            print(f"add epitopes log {leftovers}")
                            complog.append(leftovers)
                            #continues calculation but break if all remaining units = 0 or negative numbers
                            if testzero(compforloop):
                                loopoverflag = False
                                complog.append("exits with 0")
                                #print("reach all-zero in arms = -1")
                                #return complog
                                break
                        else:
                            loopoverflag = False
                            #print("reach negative value in arms = -1")
                            #return complog
                            break
                print(f"how many times here {complog}")
                break
            else:
                #print("the composition couldn't afford even simplest NG core.")
                break
                #return complog #get an empty []
        elif armsno == 0:
            #print("Reached the arms no = 0 after descending enumerates arms no")
            break
        else:
            ########################################################################################
            print("need to assign forcalcomp everytime in the beginning: get original composition and start dealing with ..")
            #print("Start enumeration from largest arms")
            #print(f"arms no: {armsno}")
            #print(f"current arms enumeration : {armdeduce[-1]}")
            #print("Hybrid and maybe bisect type need unique armsno")
            if armsno < armdeduce[-1]:
                #print("skipping impossible arms enumeration")
                continue
            elif armsno == armdeduce[-1]:
                #get the composition left fitting current arms no, for example biantennary takes H3N2+ H2N2 = H5N4
                valdict = {units: forcalcomp[units] - armdeduce[1].get(units, 0) for units in forcalcomp}
                print(f"arms testing and now the valdict is{valdict} with enumerating on {armdeduce[1]}")
                print("Reset and run the following enumeration for maximumarms -1 time")
                #check if the composition gets zero. If so, record and exit. Otherwise start terminal enumeration
                validflag = unittest(valdict)
                if validflag: #invalflag is False:
                    if testzero(valdict): #假如結構在N-core就剛好完全用掉，就應該停止計算
                        complog.append(armdeduce[0]) #armdeduce ->corededuce #log ->complog
                        complog.append("end w/ only core+arms.")
                        #print(f"core arms fits zero. Continue to skip this iteration")
                        armsno -=1
                        continue
                    else:
                        ########################################################################################
                        #armsno = global one; armsforepitopes = local one for enumerating all epitopes
                        print("Continue on extensivble epitope enumeration based on selected core type")
                        #try to add glycotopes from as much arms as possible. when not possible, -1 -> -2 arm
                        #assign temp composition storing space and isolate the calculating composition
                        while armsno > 0: #global one
                            #when armsno -1 -> new round of enumeration begins
                            tempcomp = [] #store epitopes when arms = 4#1 7 6 4 
                            armsforepitopes = armsno 
                            print(f"now temp looping arms = {armsforepitopes}")
                            print("will reset and start new enumeration from armsno -1...until it gets 0")
                            ff = forcalcomp #reset the composition to original - H3N2
                            valdict1 = valdict #reset the composition to original - (core+arms)
                            #do sth on terminal calculation
                                #SET arms no = max from armsno (for example 4)
                                #if adding terminal is possible
                                    #if meets zero:
                                        #record the epitope [arms+terminal+(debug)end indicator], CONTINUE to next terminal test
                                    #if not meeting zero:
                                        #do sth on extensive calculation
                                            #if if adding extensive is possible:
                                                #if meets zero:
                                                    #record the epitope [arms+ter+ext+(debug)indicator]
                                                #else terminal not meets zero:
                                                    #update the dict and CONTINUE on other extensive epitope?
                                #else
                                    #CONTINUE to next terminal
                            while armsforepitopes > 0:
                                print(f"current arms for enumeration LIMITED {armsforepitopes} arms")
                                tmplog = [] #store terminal and extensive informatiom, reset to None in each "loop"
                                avgunit = {}
                                for reminame, reminno in ff.items(): # valdict1 -> forcalcomp, so arms are getting calculated as well
                                    avgunit[reminame] = (reminno/armsforepitopes)
                                print(f"arms {armsforepitopes} and average units every time arms -1 before terminal{avgunit}")
                                #try to add terminal epitopes on "arms" (LIMITED)
                                reminflag = []
                                for name, reminno in avgunit.items():  #use average unit to get correct calc answer
                                    if reminno < 1:
                                        reminflag.append(False)
                                        #print(f"there is no integer unit left for adding on {armsforepitopes} arm(s)")
                                    elif reminno >= 1.0:
                                        #print(f"integer unit left, able for adding terminal structure first")
                                        reminflag.append(True)
                                #if any of units has more than one, try to add it on all arms first
                                if (True in reminflag): 
                                    print("possible for adding terminal units...doing")
                                    #get current max arms for enumeration
                                    #enuarms = armsforepitopes 
                                    for i in group_terminal:
                                        print("debug mode, consider adding units from 4->1 arms while keep origin arm no info")
                                    #for example max = 4 arms and if 4 arms has no units available, start from 3 and so does here
                                        
                                        a, b, secvaldict = compactcompcalc(avgunit, i, 2) #calculating on "avg unit"
                                        if a:
                                            #if the terminal epitope can be added do the next 
                                            c = "".join(b) #convert epitope log list to str
                                            print("recovering remaining dict")
                                            for updatecomp, updatecompno in secvaldict.items():
                                                valdict1[updatecomp] = (updatecompno * armsforepitopes)
                                                print("Iterate on extensive units")
                                                print(f"debug: now the dict remaining is {secvaldict}")
                                                print(f"debug: now the updated original dict is {valdict1}")
                                                #suppose "each type" of terminal structure only needs to be iterate once
                                                if testzero(secvaldict):
                                                    #means the core+arm+terminal equals to composition, done
                                                    #add arm information and appended terminal epitopes
                                                    tmplog.append(armdeduce[0])
                                                    tmplog.append(c)
                                                    tmplog.append("terminal units fits zero.")
                                                    tempcomp.append(tmplog)
                                                    print(f"debug: logs are {tmplog} right now, continuing on next terminal epitope")
                                                    #print(f"if there's only one arm left and we get many units possible, how?")
                                                    continue #go to next terminal iterations
                                                #else are skipped since we have continue^, if not >> one level right
                                                for j in group_extensive:
                                                    aa, bb, finvaldict = compactcompcalc(secvaldict, j, 3)
                                                    print("trying to add extensive epitopes")
                                                    if aa:
                                                        cc = "".join(bb)
                                                        finalresult = testzero(finvaldict)
                                                        if finalresult:
                                                            #means the core+arm+terminal+extension equals to composition, done
                                                            print(f"return logs: {complog}")
                                                            tmplog.append(armdeduce[0])
                                                            tmplog.append(c)
                                                            tmplog.append(cc)
                                                            tmplog.append("terminal units fits zero.")
                                                            tempcomp.append(tmplog)

                                                            break #not sure if this words
                                                        else:
                                                            #print("No logs possible")
                                                            pass
                                                        #added to hey FAQ
                                                        tmplog.append(armdeduce[0]) #newly added
                                                        tmplog.append(c)
                                                        print("Add extensive glycotopes...")
                                                        #tmplog.append(cc)
                                                        print("Finished testing on extensive units")
                                                        complog.append(cc)

                                                    else:
                                                        #print("The extensive units aren't able to be added")
                                                        pass
                                            #print(f"testing {j}")
                                            #compactcompcalc(valdict, j) #need to add logs if successfully predict
                                            ff.update(valdict1) #do this to update back to the enumeration loop
                                    #F2H5N4S2 here becomes F2S2, if 2arms: get one F,S off
                                    #forcalcomp provides F2H2N2S2, avg gets 1111          
                                    #update current dict?
                                    #for reminame, reminno in valdict1.items():
                                    #    avgunit[reminame] = (reminno/enuarms)
                                    #print(f"arms {enuarms} and average units every time arms -1 before terminal{avgunit}")
                                    #a= True/False b= epitope log c= remaining composition
                                        
                                    #when we're assigning "items" all on 1 arm, need to consider original arms no
                                    #and give all permutations of "how the get added on all locations"
                                    #coreFuc should be removed from terminal?
                                print("after one turn of enumeration, LIMITED arms-1")
                                armsforepitopes -= 1
                            #######################################################################
                            print("after all calculation in this arm no, arms-1")          
                            armsno -= 1
                            print("append tempcomp")
                            complog.append(tempcomp)
                
                            #######################################################################
                            #calculate average unit numbers. Do this every time when arms no updates (iterated -1)
                            # 1 7 6 4 - 0 7 6 0 = 1 0 0 4 at valdict1
                            

                            #print("If none of units showing more than one, add all of them to first arm")

                                    else:
                                        print("not a = the terminal isn't possible, go to the next one")
                                        pass #or continue?
                                    #print(f"current logs after core-terminal-extensive nested loop \n{complog}")
                                    print("think about if extensive units need to be added multiple times")
                                #level for i in group_terminal:
                                    tempcomp.append(tmplog)
                                    print(f"appending logs {tmplog}")
                            else:
                                #print("reminflag is False, should be no further calculation")
                                pass
                            #print("The terminal and extensive unit calculation is done")
                            #print("armsforepitopes - 1 here after enumeration")
                            armsforepitopes -=1
                            print(f"current total logs are {tempcomp}")
                        armsno -= 1
                        #print("armsno - 1 here after enumeration")
                else:
                    armsno -= 1 #reduce one arm to see if less arms results in non-negative value?
                    print("Should not entering this exception part...")
                    #pass
                    #print("this arms isn't possible...could that be happening here?")
            else:
                print("error when doing tricks on armsno...")
        #print(f"logs are {complog}") 
        finalresults.append(complog)
    print(f"ending...")
    return finalresults


    
def deducecompositiona(specific_epitopes, speccompositions):
    print("developing, now testing without enzyme and fragments")
    #create index first
    index = 0
    #print(f"len of input compositions: {len(speccompositions)}")  #this calculated all strs so it will give like 6~10 len of each comp
    print("Input pd.Series")
    strlist1 = speccompositions["predictedcomp"]
    list2 = ast.literal_eval(strlist1)
    #add log container
    jjlog = {}
    for j in range(len(list2)):
        #print(type(list2[j][1]))
        #print(list2[j][1])
        complist = ast.literal_eval(list2[j][1]) #convert string to tuple
        spectrapredcomp = predcomplsttodict(complist) #convert tuple to dict (add KDN if structure has "6" components)
        #this extracted tuple value is what we're going to calculate...or manipulate? in real.
        #deducer(specific_epitopes, spectrapredcomp)
        jjlog[j] = deducer(demo_epitopes, spectrapredcomp)  #using demo_epitopes, not the specific epitopes
    return jjlog

#calculate time spend
start = time.time()
#read 
sample1=pd.read_csv('Annotated_revised_zfNGintenstine_20230612_cloud.csv', sep='\t')

deduceinput = findpossibleepitope(c, zfdemo)
speccompositions = sample1.head(1) #24 for one spectrum has multiple assignments
tmpseries = sample1.iloc[23]
#print(type(sample1.iloc[2])) #series
strlist123 = tmpseries["predictedcomp"]
#print(f"the strlist is {strlist123}")
#print(tmpseries)
#deducecompositiona(deduceinput, speccompositions)  #for df
final_output = deducecompositiona(deduceinput, tmpseries)
print(f"the enumeration output is {final_output}")
#print(type(speccompositions))  #it's pandas df
timer1 = time.time()-start
timer2 = time.time()
print('running through the load csv to strucutre prediction', timer1, 'seconds.')
#https://stackoverflow.com/questions/36459969/how-to-convert-a-list-to-a-dictionary-with-indexes-as-values isn't working

#second prompt
tmpseries1 = sample1.iloc[2]
final_output2 = deducecompositiona(deduceinput, tmpseries1)
print(f"the enumeration output is {final_output2}")
timer3 = time.time()-start #from beginning
timer4 = time.time()-timer2 #diff from previous
print('running through the load csv to strucutre prediction', timer3, 'seconds.')
print('delta time from previous analysis is',  timer4, 'seconds more.')
