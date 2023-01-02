from pymsfilereader import MSFileReader

rawfile = MSFileReader("zf_sPerMeOG_intestine.raw")

print('Version', rawfile.Version())
print('GetFileName', rawfile.GetFileName())
print('test')


print(type(rawfile.GetMassListFromScanNum(1,scanFilter="",
                               intensityCutoffType=1,
                               intensityCutoffValue=100,
                               maxNumberOfPeaks=4,
                               centroidResult=False,
                               centroidPeakWidth=0.0)))

print(len(rawfile.GetMassListFromScanNum(1,scanFilter="",
                               intensityCutoffType=1,
                               intensityCutoffValue=100,
                               maxNumberOfPeaks=4,
                               centroidResult=False,
                               centroidPeakWidth=0.0)))

a = rawfile.GetMassListFromScanNum(1,scanFilter="",
                               intensityCutoffType=1,
                               intensityCutoffValue=5000,
                               maxNumberOfPeaks=0,
                               centroidResult=False,
                               centroidPeakWidth=0.0)[0]

b = a[0]
c = a[1]

#peaklist
'''
with open('testtest1.tsv', 'wt') as f:
    print('\t'.join(map(str, ('peak num',
                              'peak intensity'
                                ))), file=f)

    for i in range(len(b)):
        print('\t'.join(map(str, (i,
                                b[i],
                                c[i],  
                                    ))), file=f)
'''

de = rawfile.GetLabelData(1)
print(type(de))
print(len(rawfile.GetLabelData(1)))

fe = de[0]
ff = de[1]

print(fe[0][0])
print(fe[1][0])
print(fe[2][0])
print(fe[3][0])
print(fe[4][0])
print(fe[5][0])
#temp

with open('testtestlabeldatafe0.tsv', 'wt') as g:
    print('\t'.join(map(str, ('dMass',
                              'dIntensity',
                              'fResolution',
                              'fBase',
                              'fNoise',
                              'charge'
                              
                                ))), file=g)

    for j in range(len(fe[0])):
        print('\t'.join(map(str, (j,   
                                  fe[0][j],
                                  fe[1][j],
                                  fe[2][j],
                                  fe[3][j],
                                  fe[4][j],
                                  fe[5][j],
                                 ))), file=g)

#temp end

with open('ff.tsv', 'wt') as g:
    print('\t'.join(map(str, (ff[0]))), file=g)

'''
with open('testtestlabeldata.tsv', 'wt') as g:
    print('\t'.join(map(str, ('peak mass',
                              'labels'
                                ))), file=g)

    for i in range(len(de)):
        print('\t'.join(map(str, (i,
                                fe[i], 
                                    ))), file=g)

'''

rawfile.Close()
