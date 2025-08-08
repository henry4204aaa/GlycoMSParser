import mspvalidator_merger as mspval
import json
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN
import pandas as pd
import csv
import math
import copy
version = "0.25"
last_update = 20250808

#v0.5: code review and need functions for better integration. Need pre-filling in subwindow (not this file) when it detect N/O/GSL glycans etc
#v0.4: add pseudo-label to converted csv (stop here and test the accuracy if possible)
#v0.3: generate in silico csv in decodecompinput and pass it to fit_composition 
#v0.25 supress precursormassv2 message
#v0.22: test version that collects all information needed properly.
#v0.2: debugging test (use the code from compositioncalc.py)
#v0.1: placeholder



#editing for future support#
def precursormassv2(composition, deri="PerMe", debug=False): #general version
    a, b, c, d, e = composition[0], composition[1], composition[2], composition[3], composition[4]
    if composition[5]: #KDN
        f = composition[5]
    else:
        f = 0
    if deri == "PerMe":
        M =  a * 174.08921 + b * 204.09977 + c * 245.12632 + d * 361.17367 + e * 391.18423 + f * 335.1700 + 46.04186 + 1.0073
        # M = a x Fucose (triangle) + b x Hexose (circle) + c x HexNAc (Square) + d * Neu5Ac (purple diamond) + e * Neu5Gc (aqua diamond) + reducing end + [H+]
        # a, b, c, d, e >=0, is integer
        if debug:
            print("precursormassv2 calculating permetylation precursor mass")
    else:
        M = a * 162 + b * 180 + c * 204 + d * 245 + e * 275 + 46.04186 + 1.0073
        if debug:
            print("precursormassv2 calculating native precursor mass (not percise mass)")
    return(M)

def compositionbuilder(f=0, h=0, n=0, s=0, g=0, k=0):
    #compositionlist=[]
    
    #for F in range (f):
    #    for H in range(h):
    #        for N in range(n):
    #            for S in range(s):
    #                for G in range(G):
    #                    for K in range(k):
    #                        compositionlist.append((F, H, N, S, G, K))
    #return (f,h,n,s,g,k)#compositionlist
    return [f,h,n,s,g,k] #for copying current glycan + fucose, tuple is immutable so we can't copy&change the value

def glycancompositionrestraints(glycantype, arms=2, options=None, profiler_version = 1):
    #no arguments needed. write the limitations down here
    #OG core = N1
    if glycantype == "NG" and profiler_version == 1: #v1: try to search, list 
        print("Run through: Man5-9 by deafult, if poly?mannose is selected [add a button in GUI], allow search for yeast-like NG, max N2H??? 20?")
        print("Run based on arm number maximum [if not defined or no button, set maximum to 4] and follow the maximum cage calculation for each")
        print("Consider testing time spent via passing composition list versus directlt generate the csv within here <- need to rename if we stick to this")
        hmcompositionlist = []#high-mannose
        hbcompositionlist = []#hybrid-type (placeholder, skip it)
        cpcompositionlist = []#complex type

        #mass = []
        #high mannose enumeration
        if "high-mannose" in options:
            if "hyperman" in options:
                max_hex = 20
            else:
                max_hex = 9
            if "mancore2fuc" in options:
                max_fuc = 2
            elif "mancorefuc" in options:
                print("[debug]: 1 core-fucose set")
                max_fuc = 1
            else:
                print("[debug]: no core-fucose set")
                max_fuc = 0
            if "truncatedman" in options:
                min_hex = 3
            else:
                min_hex = 5
            print("Start highman_calculation")   
            #produce high-mannose glycans based on rules defined
            for hex in range(min_hex,max_hex+1):
                hmcompositionlist.append(compositionbuilder(0, hex, 2))
            #if core fucose is allowed
            if max_fuc != 0:
                copylist = copy.deepcopy(hmcompositionlist)
                for fuc in range(1,max_fuc+1):
                    for elements in copylist:
                        elements[0] = fuc
                        hmcompositionlist.append(elements)
            print(f"High mannose output is {hmcompositionlist}")
        else:
            print("No high-mannose selected from options")
        #complext type enumeration
        #NG core = H3N2
        #NG arm limits = 4
        #NG extension limits = 2
        #NG allowed maximum each "cage": 5 sugars, HexNac + Hex = 3, Fuc + Neu = 2 , all Neu should be 1 in one arm
        NGcore = [0,3,2,0,0,0]
        results = []
        armslog = []
        count = 0
        #default arms no = 2
        if not arms:
            arms = 2 #default
        for i in range(1,arms+1):
            for hex in range(3):  # max a+b ≤ 2
                for hexnac in range(3 - hex):  # ensures a + b ≤ 2
                    for neuac in range(2):  # max c+d+e ≤ 1
                        for neugc in range(2 - neuac):  # ensures c + d ≤ 1
                            for kdn in range(2 - neuac - neugc):  # ensures c + d + e ≤ 1
                                for fuc in range(3 - neuac - neugc - kdn):  # ensures c+d+e+f ≤ 2
                                    total = hex + hexnac + neuac + neugc + kdn + fuc
                                    if total <= 5:
                                        results.append([fuc, hex+3, hexnac+2, neuac, neugc, kdn])
                                        count+=1
                                        armslog.append([arms, count])
        #default extension = 0, if add extension, copy the complex above and add extension based on arms no

        print(f"Total combinations: {armslog}, and detail is {results}")



        
    else:
        print("The profiler version is not defined, use default in silico production")
        if glycantype == "NG":
            allcomp = (5, 9, 10, 4, 4) #Max F=5,N=9, H=10, S=4, G=4, KDN= None
            return allcomp
        elif glycantype == "OG":
            allcomp = (4, 5, 6, 4, 4) 
        else:
            print("Not NG or OG selected, returning None to raise Errors")
            return None
        
        return allcomp

    NG_limitations = {"0":0}

def decodecompinput(config, infofrommetadata):
    #metadata needs: derivatization type (such as PerMethylated), charge mode (positive/negative), reduced/free-end
    hex = config.get("Hex", 0)
    hexnac= config.get("HexNAc", 0)
    fuc = config.get("Fuc", 0)
    Neu5Ac = config.get("Neu5Ac", 0)
    Neu5Gc = config.get("Neu5Gc", 0)
    KDN = config.get("KDN", 0)
    
    if infofrommetadata:
        with open(infofrommetadata, "r") as f:
            metadata = json.load(f)
    else:
        print("No metadata selected")
    glycantype = metadata.get("Glycan Type", None) 
    mode = metadata.get("Mass Analyzer charge mode", None)
    derivatization = metadata.get("Derivatization Type", None) #values=["PerMe(Freeend)", "PerMe(Reduced)", "Others"],  # Add more if needed
    
    #detect N or O glycan and automatically import default settings if possible
    #N-glycans: HexNAc >= 2, Hex >= 3, with other minor rules included
    #O-glycans: HexNAc >= 1, with other minor rules included

    #mode check
    if mode == "+":
        mode_protonated = 1
    elif mode == "-":
        mode_protonated = -1
    else:
        print("mode is out of range")
        mode_protonated = 0
    #derivatization check
    deri_mass = mspval.etc_masses.get(derivatization)

    if "PerMe(Freeend)" == derivatization or "PerMe(Reduced)" == derivatization:
        permesugar = mspval.monosaccharide_masses.copy()
    elif "Native" == derivatization:
        print("Native glycan setting is not implemented")
        permesugar = None
    elif derivatization is None:
        print("Warning, derivatization information missing, use PerMe as default")
        permesugar = mspval.monosaccharide_masses.copy()
    else:
        raise ValueError("Derivatization data cannot be retrieved by exceptions in derivatization values")

    #from mspvalidator_merger.py etc_masses
    decodeconfig = {"Hex": hex, "HexNAc": hexnac, "Fuc": fuc, "Neu5Ac": Neu5Ac, "Neu5Gc": Neu5Gc, "KDN": KDN, "Glycan type (for pre-filling)": glycantype,
                    "mode": mode, "Derivatization (for proper sugar mass)": derivatization, "Monosaccharide unit mass": permesugar}
    return decodeconfig
#  File "G:\其他電腦\My Computer\GlycoMSParser\2025_demo\mspcomposition.py", line 11, in decodecompinput
#    glycantype = infofrommetadata.get("Glycan Type", None)
#AttributeError: 'str' object has no attribute 'get'


#- code adapted from compositioncalc.py

#for insilicover2
def glycancompositionver2(a, b, c, d, e, f):
    M =  a * 174.08921 + b * 204.09977 + c * 245.12632 + d * 361.17367 + e * 391.18423 + f * 335.1700 +  46.04186 + 1.0073
    # M = a x Fucose (triangle) + b x Hexose (circle) + c x HexNAc (Square) + d * Neu5Ac (purple diamond) + 
    # e * Neu5Gc (aqua diamond) + f * KDN (green diamond) reducing end + [H+]
    # a, b, c, d, e, f>=0, is integer
    return(M)

#F0~6, H3~10, N2~11, S0~5, G0~5, K0~4 in silico glycan list#
def insilicover2():
    compositionlist = []
    mass = []
    for F in range (6):
        for H in range(3,10):
            for N in range(2,11):
                for S in range(5):
                    for G in range(5):
                        for K in range(4):
                            compositionlist.append((F, H, N, S, G, K))
                            calcmass = glycancompositionver2(F, H, N, S, G, K)
                            mass.append(Decimal(str(calcmass)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
    with open('glycaninsilicowithkdn_fixed.csv', 'w', encoding = 'utf-8', newline='') as z:
        rowlists = zip(tuple(compositionlist), mass)
        writer = csv.writer(z)
        for row in rowlists:
            writer.writerow(row)
#

def fit_composition (config, selected_csv, metadata=None, debug=True):
    print(f"config is {config}")
    print(f"selected csv is {selected_csv}")
    print(f"metadata is {metadata}")
    decodedconfig = decodecompinput(config, infofrommetadata=metadata)
    print(f"collected unzipped config and metadata info:{decodedconfig}")


#test run code

glycancompositionrestraints("NG", arms=None, options=["high-mannose","mancorefuc"], profiler_version = 1)