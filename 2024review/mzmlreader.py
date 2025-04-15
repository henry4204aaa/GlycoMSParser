import os
import numpy as np
version = 0.1
last_update = 20241112


import pymzml

# Load an mzML file
run = pymzml.run.Reader(f"G:\zf_sPerMeNG_brain.mzML")

# Iterate through the spectra
t = 0

for spectrum in run:
    if t > 0:
        break
    else:
        if spectrum['ms level'] == 2:
            empty_containerpeaks = []
            empty_containerintensity = []
            print(f"Getting peak list")
            p_peaklist = spectrum.mz
            p_intensity = spectrum.i
            #print(f"type of peaklist is {type(p_peaklist)} and the output is {p_peaklist}")
            for peaks in p_peaklist:
                print(f"peaks are in type of {type(peaks)} and the data inside is {peaks}")
                #peaks are in type of <class 'numpy.float32'> and the data inside is 97.00896453857422
                empty_containerpeaks.append(np.round(peaks, 4))
            print(f"peaks after pre-processing rounding is {empty_containerpeaks}")
            for intensities in p_intensity:
                empty_containerintensity.append(np.round(intensities, 4))
            print(f"peaks after pre-processing rounding is {empty_containerintensity}")


            #peaklist = p_peaklist.split()
            #intensiy = p_intensity.split()
            #print(f"peaklist is {peaklist} and intensity is {intensiy}")
            t += 1
'''
for spectrum in run:
    if t > 3:
        break
    elif spectrum['ms level'] == 2:
        print(f"type of this spectrum is {type(spectrum)}")
        print(f"testing pymzml functionality")
        print(f"Spectrum ID: {spectrum.ID}, {type(spectrum.ID)}")
        print(f"MS Level: {spectrum.ms_level}, {type(spectrum.ms_level)}")
        print(f"Retention Time: {spectrum.scan_time_in_minutes()} minutes, {type(spectrum.scan_time_in_minutes)}")
        print(f"Precursor m/z: {spectrum.selected_precursors}, {type(spectrum.selected_precursors)}")
        print(f"Polarity: {spectrum.get('scan polarity')}, {type(spectrum.get('scan polarity'))}")
        print(f"Instrument Type: {spectrum.get('instrument model')}, {type(spectrum.get('instrument model'))}")
        # Access other metadata as needed
        print(spectrum.mz, spectrum.i, type(spectrum.mz), type(spectrum.i))  
        t+=1
#print("testing certain ms number")
#spectrum_with_id_311 = run[ 311 ]
#print(spectrum_with_id_311())
print("testing run itself")
spectrum = next(run)
print(dir(spectrum))
'''
print(f"Testing mzml reader finished")

#docs are summarized by GPT4o (20241104) I can't find clear instructions on internet 
