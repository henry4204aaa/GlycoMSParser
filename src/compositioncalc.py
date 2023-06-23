from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN
import pandas as pd
import csv
import math

def glycancomposition(a, b, c, d, e):
    M =  a * 174.08921 + b * 204.09977 + c * 245.12632 + d * 361.17367 + e * 391.18423 + 46.04186 + 1.0073
    # M = a x Fucose (triangle) + b x Hexose (circle) + c x HexNAc (Square) + d * Neu5Ac (purple diamond) + e * Neu5Gc (aqua diamond) + reducing end + [H+]
    # a, b, c, d, e >=0, is integer
    return(M)
#F1H5N4S1G



allcomp = input("enter the composition")
if not allcomp:
    allcomp = (5, 9, 10, 4, 4)
    a, b, c, d, e = allcomp[0],allcomp[1],allcomp[2],allcomp[3],allcomp[4]
    print(f"a = {a}, b = {b}, c= {c}, d= {d}, e= {e}")
    #core 2N3H + 0-1F
    #max arm no. = 4
    #each arm + 1N 1-2H 0-1F 0-1S
    compositionlist = []
    mass = []
    glycanbuild = ["highman", "hybrid", "bian", "trian", "tetraan"]
    for i in glycanbuild:
        #high man = 2H + 5~9N
        if i == "highman":
            F = 0
            S = 0
            G = 0
            N = 2
            for H in range(5,10):
                compositionlist.append((F, H, N, S, G))
                calcmass = glycancomposition(F, H, N, S, G)
                mass.append(Decimal(str(calcmass)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
        #hybrid = 2N + 5H + 0~1F + 1~2 arm (1N, 0~2H, 0~1F, 0~1S, 0~1G)    
        else:
            print("finished test")
    print(compositionlist)
    print(mass)
'''
        elif i == "hybrid":
            for coreFuc in range(0,2):
                for arm in range (1,3):
                    N = 2 + arm
                    for Hex in range (0,3):
                        if Hex == 0:
                            F = coreFuc
                            S = 0
                            G = 0
                            H = 5
                            compositionlist.append((F, H, N, S, G))
                            calcmass = glycancomposition(F, H, N, S, G)
                            mass.append(Decimal(str(calcmass)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
'''       


         
def insilico():
    compositionlist = []
    mass = []
    for F in range (6):
        for H in range(3,10):
            for N in range(2,11):
                for S in range(5):
                    for G in range(5):
                        compositionlist.append((F, H, N, S, G))
                        calcmass = glycancomposition(F, H, N, S, G)
                        mass.append(Decimal(str(calcmass)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))
    with open('glycaninsilico.csv', 'w', encoding = 'utf-8', newline='') as z:
        rowlists = zip(tuple(compositionlist), mass)
        writer = csv.writer(z)
        for row in rowlists:
            writer.writerow(row)



def glycancompositionver2(a, b, c, d, e, f):
    M =  a * 174.08921 + b * 204.09977 + c * 245.12632 + d * 361.17367 + e * 391.18423 + f * 335.1700 +  46.04186 + 1.0073
    # M = a x Fucose (triangle) + b x Hexose (circle) + c x HexNAc (Square) + d * Neu5Ac (purple diamond) + 
    # e * Neu5Gc (aqua diamond) + f * KDN (green diamond) reducing end + [H+]
    # a, b, c, d, e, f>=0, is integer
    return(M)

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
    with open('glycaninsilicowithkdn.csv', 'w', encoding = 'utf-8', newline='') as z:
        rowlists = zip(tuple(compositionlist), mass)
        writer = csv.writer(z)
        for row in rowlists:
            writer.writerow(row)
#F1H5N4S1G
