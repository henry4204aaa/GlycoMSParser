from pymsfilereader import MSFileReader

rawfile = MSFileReader("zf_sPerMeOG_intestine.raw")

print('Version', rawfile.Version())
print('GetFileName', rawfile.GetFileName())
print('test')

print('i = spectrum number')
print('j = loop for MS1 spec')

k = input('enter the number of spectrum you want to test')

for i in range(rawfile.FirstSpectrumNumber, rawfile.LastSpectrumNumber + 1):
    print 

    if rawfile.GetMSOrderForScanNum(i)==1:
        print(i, 'spectrum is MS1.')
        print('Do something on MS1')
    elif rawfile.GetMSOrderForScanNum(i)==2:
        print(i, 'spectrum is MS2.')
        print('Do something on MS2 \nOK?')
    else:
        print(i, 'spectrum is MS3.')
    if i == int(k):
        print('stop processing after',k,'spectrum')
        break
