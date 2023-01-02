from pymsfilereader import MSFileReader

#assign file in code
'''
rawfile = MSFileReader("zf_sPerMeOG_intestine.raw")
print('test')
print('File choosed', rawfile.Version())
print('GetFileName', rawfile.GetFileName())
'''

#choose file
#please add file extension filter later
try:
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
    Tk().withdraw()
    filename = askopenfilename()
    print('Selected', filename)
    rawfile = MSFileReader(filename)
    print('test')
    print('File choosed', rawfile.Version())
    print('GetFileName', rawfile.GetFileName())
except ImportError:
    raise ImportError('Please install tkinter')



def findingMSlevel():
    k = input('enter the number of spectrum you want to test')
    for i in range(rawfile.FirstSpectrumNumber, rawfile.LastSpectrumNumber + 1):
        if rawfile.GetMSOrderForScanNum(i)==1:
            print(i, 'spectrum is MS1.')
            print('Do something on MS1')
        elif rawfile.GetMSOrderForScanNum(i)==2:
            print(i, 'spectrum is MS2.')
            print('Do something on MS2.')
        else:
            print(i, 'spectrum is MS3.')
        if i == int(k):
            print('stop processing after',k,'spectrum')
            break
    print('No information exported')
    rawfile.Close()
    print('rawfileclosed')

#add one more function which process the extraction automatically
#donot ask user input

#this one run from the first to defined spctrum number
def ExtractMSlevel():
    k = input('enter the number of spectrum you want to test')
    for i in range(rawfile.FirstSpectrumNumber, rawfile.LastSpectrumNumber + 1):
        if rawfile.GetMSOrderForScanNum(i)==1:
            print(i, 'spectrum is MS1.')
            de = rawfile.GetLabelData(i)
            fe = de[0]
            ff = de[1]
            #print(fe[0][0])
            #print(fe[1][0])
            #print(fe[2][0])
            #print(fe[3][0])
            #print(fe[4][0])
            #print(fe[5][0])
            fileext = 'no' + str(i) +'isMS1' +  '.tsv'
            with open(fileext, 'wt') as g:
                print('\t'.join(map(str, ('No of entry',
                                          'dMass',
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
            
        elif rawfile.GetMSOrderForScanNum(i)==2:
            print(i, 'spectrum is MS2.')
            print('Do something on MS2.')
            de = rawfile.GetLabelData(i)
            print(type(de))
            #print(len(rawfile.GetLabelData(i)))
            fe = de[0]
            ff = de[1]
            #print(fe[0][0])
            #print(fe[1][0])
            #print(fe[2][0])
            #print(fe[3][0])
            #print(fe[4][0])
            #print(fe[5][0])
            fileext = 'no' + str(i) +'isMS2' +  '.tsv'
            with open(fileext, 'wt') as g:
                print('\t'.join(map(str, ('No of entry',
                                          'dMass',
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
            

        else:
            print(i, 'spectrum is MS3. Nothing has been done')
        if i == int(k):
            print('stop processing after',k,'spectrum')
            break
    print('see extracted tsv files')
    rawfile.Close()
    print('rawfileclosed')


#for specific
def SExtractMSlevel():
    i = int(input('enter the specific spectrum you want to test'))
    if i < rawfile.LastSpectrumNumber:
        print('valid spectrum')
        if rawfile.GetMSOrderForScanNum(i)==1:
            print(i, 'spectrum is MS1.')
            de = rawfile.GetLabelData(i)
            fe = de[0]
            ff = de[1]
            #print(fe[0][0])
            #print(fe[1][0])
            #print(fe[2][0])
            #print(fe[3][0])
            #print(fe[4][0])
            #print(fe[5][0])
            fileext = 'no' + str(i) +'isMS1' +  '.tsv'
            with open(fileext, 'wt') as g:
                print('\t'.join(map(str, ('No of entry',
                                          'dMass',
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
            print('see extracted tsv file')
        elif rawfile.GetMSOrderForScanNum(i)==2:
            print(i, 'spectrum is MS2.')
            print('Do something on MS2.')
            de = rawfile.GetLabelData(i)
            print(type(de))
            #print(len(rawfile.GetLabelData(i)))
            fe = de[0]
            ff = de[1]
            #print(fe[0][0])
            #print(fe[1][0])
            #print(fe[2][0])
            #print(fe[3][0])
            #print(fe[4][0])
            #print(fe[5][0])
            fileext = 'no' + str(i) +'isMS2' +  '.tsv'
            with open(fileext, 'wt') as g:
                print('\t'.join(map(str, ('No of entry',
                                          'dMass',
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
            

        else:
            print(i, 'spectrum is MS3. Nothing has been done')
    else:
        print('invalid spectra %i')
   
    print('Enter rawfile.Close() after ending processing')
