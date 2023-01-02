from pymsfilereader import MSFileReader

#rawfile = MSFileReader("zf_sPerMeOG_intestine.raw")
rawfile = MSFileReader("zf_sPerMeNG_brain")

print('Version', rawfile.Version())
print('GetFileName', rawfile.GetFileName())

#scan_number = 3


#block for assigning spectrum no
intcheck = False

while intcheck == False:
    scan_number = input('enter spectrum no you want to show')
    print('The input is', scan_number)
    try:
        scan_number = int(scan_number)
        break
    except ValueError:
        print ('please enter a int')

#end
print('enter test(n) to run it')
print('close() before leave')
        

'''
print('############################################## XCALIBUR INTERFACE BEGIN')
print('GetScanHeaderInfoForScanNum',
        rawfile.GetScanHeaderInfoForScanNum(scan_number))  # "View/Scan header", upper part
print('GetTrailerExtraForScanNum', rawfile.GetTrailerExtraForScanNum(scan_number))  # "View/Scan header", lower part

print("MSOrder", rawfile.GetMSOrderForScanNum(scan_number),'\n',
      "Scan number", rawfile.ScanNumFromRT(rawfile.RTFromScanNum(scan_number)),'\n',
      "Scan filter", rawfile.GetFilterForScanNum(scan_number),'\n',
      "Precursor", rawfile.GetPrecursorInfoFromScanNum(scan_number),'\n')
'''
def test1(scan_number):
    print("MSOrder", rawfile.GetMSOrderForScanNum(scan_number),'\n',
          "Scan number", rawfile.ScanNumFromRT(rawfile.RTFromScanNum(scan_number)),'\n',
          "Scan filter", rawfile.GetFilterForScanNum(scan_number),'\n',
          "Precursor", rawfile.GetPrecursorInfoFromScanNum(scan_number),'\n')
    print('#############')
    print('GetScanHeaderInfoForScanNum',
            rawfile.GetScanHeaderInfoForScanNum(scan_number))  # "View/Scan header", upper part
    print('#############')
    print('GetTrailerExtraForScanNum', rawfile.GetTrailerExtraForScanNum(scan_number))  # "View/Scan header", lower part


def test2(scan_number):
    print('scan number=', scan_number)
    print('MSOrder', rawfile.GetMSOrderForScanNum(scan_number),'\n',)
    print("Scan filter", rawfile.GetFilterForScanNum(scan_number),'\n')
    print('Isolation mass', rawfile.GetPrecursorInfoFromScanNum(scan_number)[0], '\n')
    print('monoIsomass', rawfile.GetPrecursorInfoFromScanNum(scan_number)[1], '\n')
    print('chargeState', rawfile.GetPrecursorInfoFromScanNum(scan_number)[2], '\n')
    print('parentScanNo', rawfile.GetPrecursorInfoFromScanNum(scan_number)[3], '\n')
    
def close():
    rawfile.Close()

