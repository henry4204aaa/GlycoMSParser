import glob
import pandas as pd
import os
import pathlib

from pymsfilereader import MSFileReader

#consider loading csv from tmp

#consider make a df containing all essential spectra info
#pseudo code
# for i in range(1,spectrumno+1)
# if MSlevel(i) = 1:
# append (extracted header) to dfms1index
# elif MSlevel(i) = 2:
# append (extracted header) to dfms2index
# elif MSlevel(i) = 3:
# append (extracted header) to dfms3index
# else:
# append (extracted header) to dfothersindex
#



#will be serialized in batch processing
dfms1 = pd.DataFrame(pd.read_table('no23043isMS1.tsv'))
dfms1b = pd.DataFrame(pd.read_table('no22998isMS1.tsv'))

#same
dfms2 = pd.DataFrame(pd.read_table('no23083isMS2.tsv'))
dfms2b = pd.DataFrame(pd.read_table('no23045isMS2.tsv'))

#read headerinfo and contact it for reference

rawfile = MSFileReader("zf_sPerMeNG_brain")
#mind that this loads RAW file and you should close them whatever
scan_number = int(input('enter spectrum no you want to show'))

#Q: can I store different type in a series?
dfheader = pd.Series([scan_number,
                      rawfile.GetPrecursorInfoFromScanNum(scan_number)[0],
                      rawfile.GetPrecursorInfoFromScanNum(scan_number)[1],
                      rawfile.GetPrecursorInfoFromScanNum(scan_number)[2],
                      rawfile.GetPrecursorInfoFromScanNum(scan_number)[3]],
                      index=["MS2scan no","Isolation mass","monoIsomass","chargeState","parentScanNo"])
rawfile.Close()
