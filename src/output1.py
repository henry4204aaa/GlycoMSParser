Python 3.10.7 (tags/v3.10.7:6cc6b13, Sep  5 2022, 14:08:36) [MSC v.1933 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license()" for more information.

============= RESTART: C:\Users\hnstseng\Downloads\MSFileReader.py =============
Version 3.0.31.0
GetFileName zf_sPerMeOG_intestine.raw
GetCreatorID Fusion
GetVersionNumber 66
GetCreationDate time.struct_time(tm_year=1970, tm_mon=1, tm_mday=1, tm_hour=12, tm_min=6, tm_sec=33, tm_wday=3, tm_yday=1, tm_isdst=0)
IsError False
IsNewFile False
IsThereMSData True
HasExpMethod True
InAcquisition False
GetErrorCode 0
GetErrorMessage 
GetWarningMessage 
RefreshViewOfFile None
GetNumberOfControllers 2
GetNumberOfControllersOfType('No device') 0
GetNumberOfControllersOfType('MS') 1
GetNumberOfControllersOfType('Analog') 0
GetNumberOfControllersOfType('A/D card') 0
GetNumberOfControllersOfType('PDA') 0
GetNumberOfControllersOfType('UV') 1
GetControllerType('MS') MS
GetCurrentController() (0, 1)
GetCurrentController() (0, 1)
GetCurrentController() (0, 1)
GetExpectedRunTime() 90.0
GetMaxIntegratedIntensity() 4356860928.0
GetMaxIntensity() 0
GetInletID() 0
GetErrorFlag() 0
GetFlags() 
GetAcquisitionFileName() 
GetOperator() 
GetComment1() 
GetComment2() 
GetFilters() 
GetMassTolerance() (False, 500.0, 0)
rawfile.SetMassTolerance(userDefined=True, massTolerance=555.0, units=2) None
GetMassTolerance() (True, 555.0, 2)
rawfile.SetMassTolerance(userDefined=False, massTolerance=500.0, units=0) None
GetMassResolution 0.5
GetNumTrailerExtra 46125
GetLowMass 90.0
GetHighMass 2000.0
GetStartTime 0.0011264189333333335
GetEndTime 90.0013702568
GetNumSpectra 46125
GetFirstSpectrumNumber 1
GetLastSpectrumNumber 46125
GetAcquisitionDate 
GetUniqueCompoundNames ()
############################################## INSTRUMENT BEGIN
GetInstrumentDescription C#000
GetInstrumentID 0
GetInstSerialNumber FSN10176
GetInstName Orbitrap Fusion
GetInstModel Orbitrap Fusion
GetInstSoftwareVersion 3.0.2041
GetInstHardwareVersion 
GetInstFlags 
GetInstNumChannelLabels 0
IsQExactive False
############################################## INSTRUMENT END
############################################## XCALIBUR INTERFACE BEGIN
GetScanHeaderInfoForScanNum OrderedDict([('numPackets', 61304), ('StartTime', 0.0011264189333333335), ('LowMass', 500.0), ('HighMass', 2000.0), ('TIC', 5333820.5), ('BasePeakMass', 536.1658935546875), ('BasePeakIntensity', 242910.546875), ('numChannels', 0), ('uniformTime', 0), ('Frequency', 0.0)])
GetTrailerExtraForScanNum OrderedDict([('', '\t\t\t\t\t\t\t\t\t         \t\t\t\t\t\t\t\t\t\t\t\t\t\t'), ('Scan Description', '                '), ('AGC', 'On          '), ('Micro Scan Count', 1.0), ('Ion Injection Time (ms)', 50.0), ('Elapsed Scan Time (sec)', 0.324), ('Average Scan by Inst', 'No'), ('Orbitrap Resolution', 120000.0), ('Access ID', 0.0), ('API Process Delay', -1.0), ('Dependency Type', 0.0), ('Multi Inject Info', '                                '), ('Master Scan Number', -1.0), ('Monoisotopic M/Z', -1.0), ('Charge State', 0.0), ('HCD Energy', -1.0), ('MS2 Isolation Width', -1.0), ('SPS Masses', '                                                                                                                                '), ('SPS Masses Continued', '                                                                                                                                '), ('Conversion Parameter I', 0.0), ('Conversion Parameter A', 0.0), ('Conversion Parameter B', 211820001.740694), ('Conversion Parameter C', -33697650.20725), ('Conversion Parameter D', 0.0), ('Conversion Parameter E', 0.0), ('Temperature Comp. (ppm)', -5.12), ('RF Comp. (ppm)', 0.0), ('Space Charge Comp. (ppm)', -1.27), ('Resolution Comp. (ppm)', -0.05), ('Number of LM Found', 0.0), ('LM Correction (ppm)', 0.0), ('RawOvFtT', 296606.2), ('Injection t0', -0.063), ('Reagent Ion Injection Time (ms)', 0.0)])
GetNumTuneData 1
GetTuneData(0) 
GetNumInstMethods 2
GetInstMethodNames ('TNG-Calcium', 'Dionex Chromatography MS Link')
-------------------------------------------------------------------------------
Orbitrap Fusion Method Summary
  
Creator: Fusion-PC\Fusion          Last Modified: 5/7/2019 11:44:39 AM by Fusion-PC\Fusion
  
Global Settings
	Use Ion Source Settings from Tune = True
	Method Duration (min)= 90
	Pressure Mode = Standard
	Default Charge State = 2
	Advanced Precursor Determination = False
Experiment 1
	Start Time (min) = 0
	End Time (min) = 90
	Cycle Time (sec) = 3
		


		Scan MasterScan
			MSn Level = 1
			Use Wide Quad Isolation = True
			Detector Type = Orbitrap
			Orbitrap Resolution = 120K
			Mass Range = Normal
			Scan Range (m/z) = 500-2000
			Maximum Injection Time (ms) = 50
			AGC Target = 400000
			Microscans = 1
			RF Lens (%) = 60
			Use ETD Internal Calibration = False
			DataType = Profile
			Polarity = Positive
			Source Fragmentation = False
			Scan Description = 

		Filter MIPS
			MIPS Mode = Peptide

		Filter ChargeState
			Include undetermined charge states = False
			Include charge state(s) = 1-4
			Include charge states 25 and higher = False

		Filter DynamicExclusion
			Exclude after n times = 1
			Exclusion duration (s) = 10
			Mass Tolerance = ppm
			Mass tolerance low = 10
			Mass tolerance high = 10
			Exclude isotopes = True
			Perform dependent scan on single charge state per precursor only = False

		Filter IntensityThreshold
			Intensity Filter Type = IntensityThreshold
			Maximum Intensity = 1E+20
			Minimum Intensity = 50000
			Relative Intensity Threshold = 0

		Data Dependent Properties
			Data Dependent Mode= Cycle Time
Scan Event 1
		


		Scan ddMSnScan
			MSn Level = 2
			Isolation Mode = Quadrupole
			Isolation Offset = Off
			Isolation Window = 1.6
			Reported Mass = Offset Mass
			Multi-notch Isolation = False
			Scan Range Mode = Auto Normal
			FirstMass = 90
			Scan Priority= 1
			ActivationType = HCD
			Is Stepped Collision Energy On = True
			Stepped Collision Energy (%) = 5
			Collision Energy (%) = 15
			Detector Type = Orbitrap
			Orbitrap Resolution = 30K
			Maximum Injection Time (ms) = 54
			AGC Target = 50000
			Inject ions for all available parallelizable time = False
			Microscans = 1
			Use ETD Internal Calibration = False
			DataType = Profile
			Polarity = Positive
			Source Fragmentation = False
			Scan Description = 

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------

Program for Dionex Chromatography MS Link

       ColumnOven.TempCtrl =         On
       ColumnOven.Temperature.Nominal =  50.0 [°C]
       ColumnOven.Temperature.LowerLimit =  20.0 [°C]
       ColumnOven.Temperature.UpperLimit =  75.0 [°C]
       EquilibrationTime =           0.1 [min]
       ColumnOven.ReadyTempDelta =   5.0 [°C]
       Sampler.TempCtrl =            On
       Sampler.Temperature.Nominal = 6.0 [°C]
       Sampler.Temperature.LowerLimit =  4.0 [°C]
       Sampler.Temperature.UpperLimit =  45.0 [°C]
       Sampler.ReadyTempDelta =      2.0 [°C]
       LoadingPump.Pressure.LowerLimit =  0 [psi]
       LoadingPump.Pressure.UpperLimit =  9000 [psi]
       LoadingPump.MaximumFlowRampDown =  10 [µl/min²]
       LoadingPump.MaximumFlowRampUp =  10 [µl/min²]
       LoadingPump.%A.Equate =       "%A"
       LoadingPump.%B.Equate =       "%B"
       %C.Equate =                   "%C"
       NC_Pump.Pressure.LowerLimit = 0 [psi]
       NC_Pump.Pressure.UpperLimit = 11600 [psi]
       NC_Pump.MaximumFlowRampDown = 0.500 [µl/min²]
       NC_Pump.MaximumFlowRampUp =   0.500 [µl/min²]
       NC_Pump.%A.Equate =           "%A"
       NC_Pump.%B.Equate =           "%B"
       DrawSpeed =                   50 [nl/s]
       DrawDelay =                   5000 [ms]
       DispSpeed =                   2000 [nl/s]
       DispenseDelay =               5000 [ms]
       WasteSpeed =                  4000 [nl/s]
       WashSpeed =                   4000 [nl/s]
       LoopWashFactor =              2.000
       SampleHeight =                2.500 [mm]
       PunctureDepth =               10.000 [mm]
       WashVolume =                  50.000 [µl]
       RinseBetweenReinjections =    Yes
       LowDispersionMode =           Off
       InjectMode =                  ulPickUp
       FlushVolume =                 5.000 [µl]
       NC_Pump_Pressure.Step =       0.01 [s]
       NC_Pump_Pressure.Average =    Off
       ValveRight =                  1_2
       FirstTransportVial =          R5
       LastTransportVial =           R5
       TransportVialCapacity =       99999
       TransLiquidHeight =           5.000 [mm]
       TransVialPunctureDepth =      10.000 [mm]
       LoadingPump.Flow =            0.000 [µl/min]
       LoadingPump.%B =              0.0 [%]
       %C =                          0.0 [%]

 0.000 Wait                          LoadingPump.Ready and NC_Pump.Ready and ColumnOven.Ready and Sampler.Ready and PumpModule.Ready
       ;Chromeleon sets this property to signal to Xcalibur that it is ready to start a run.
       ReadyToRun =                  1
       ;Xcalibur sets this property to start the run or injection.
       Wait                          StartRun
       NC_Pump.Flow =                0.500 [µl/min]
       NC_Pump.%B =                  1.0 [%]
       NC_Pump.Curve =               5
       Wait                          LoadingPump.Ready and NC_Pump.Ready and ColumnOven.Ready and Sampler.Ready and PumpModule.Ready
       Inject
       NC_Pump_Pressure.AcqOn
       ;Chromeleon sets this property to signal the injection to Xcalibur.
       InjectResponse =              1
       ;Depending on your system configuration it might be necessary to manually insert
       ;a "Relay" command below in order to send the start signal to the MS.
       ;Typical syntaxes:
       ;Pump_Relay_1.Closed  Duration = 2.00
       ;UM3PUMP_Relay1.On    Duration = 2.00
       Relay_4.State Off
       NC_Pump.Flow =                0.500 [µl/min]
       NC_Pump.%B =                  1.0 [%]
       NC_Pump.Curve =               5

 0.010 Relay_4.State On

 1.000 Relay_4.State Off
       NC_Pump.Flow =                0.500 [µl/min]
       NC_Pump.%B =                  1.0 [%]
       NC_Pump.Curve =               5

 1.010 NC_Pump.Flow =                0.500 [µl/min]
       NC_Pump.%B =                  30.0 [%]
       NC_Pump.Curve =               5

45.000 NC_Pump_Pressure.AcqOff

70.000 NC_Pump.Flow =                0.500 [µl/min]
       NC_Pump.%B =                  80.0 [%]
       NC_Pump.Curve =               5

70.100 NC_Pump.Flow =                0.500 [µl/min]
       NC_Pump.%B =                  90.0 [%]
       NC_Pump.Curve =               5

80.000 NC_Pump.Flow =                0.500 [µl/min]
       NC_Pump.%B =                  90.0 [%]
       NC_Pump.Curve =               5

80.100 NC_Pump.Flow =                0.500 [µl/min]
       NC_Pump.%B =                  1.0 [%]
       NC_Pump.Curve =               5

90.000 NC_Pump.Flow =                0.500 [µl/min]
       NC_Pump.%B =                  1.0 [%]
       NC_Pump.Curve =               5
       InjectResponse =              0
       End

-------------------------------------------------------------------------------
rawfile.ExtractInstMethodFromRaw None
GetVialNumber 0
GetInjectionVolume 0.0
GetInjectionAmountUnits 
GetSampleVolume 0.0
GetSampleVolumeUnits 
GetSampleWeight 0.0
GetSampleAmountUnits 
GetSeqRowNumber 8
GetSeqRowSampleType Unknown
GetSeqRowDataPath 
GetSeqRowRawFileName F:\IBC_LAB603\Hsiao\20190508\zf_sPerMeOG_intestine.raw
GetSeqRowSampleName 
GetSeqRowSampleID RA1
GetSeqRowComment 
GetSeqRowLevelName 
GetSeqRowUserText 
GetSeqRowUserText 
GetSeqRowUserText 
GetSeqRowUserText 
GetSeqRowUserText 
GetSeqRowInstrumentMethod F:\IBC_LAB603\Hsiao\20190508\ON_Glycan-90min_HCDpdCID_zfOG.meth
GetSeqRowProcessingMethod 
GetSeqRowCalibrationFile 
GetSeqRowVial BD5
GetSeqRowInjectionVolume 2.0
GetSeqRowSampleWeight 0.0
GetSeqRowSampleVolume 0.0
GetSeqRowISTDAmount 0.0
GetSeqRowDilutionFactor 1.0
GetSeqRowUserLabel Study
GetSeqRowUserLabel Client
GetSeqRowUserLabel Laboratory
GetSeqRowUserLabel Company
GetSeqRowUserLabel Phone
GetSeqRowUserTextEx 
GetSeqRowUserTextEx 
GetSeqRowUserTextEx 
GetSeqRowUserTextEx 
GetSeqRowUserTextEx 
GetSeqRowBarcode 
GetSeqRowBarcodeStatus 0
GetNumStatusLog 2838
GetStatusLogForScanNum
Traceback (most recent call last):
  File "C:\Users\hnstseng\Downloads\MSFileReader.py", line 143, in <module>
    pprint(rawfile.GetStatusLogForScanNum(scan_number))
NameError: name 'pprint' is not defined. Did you mean: 'print'?

============= RESTART: C:\Users\hnstseng\Downloads\MSFileReader.py =============
Version 3.0.31.0
GetFileName zf_sPerMeOG_intestine.raw
GetCreatorID Fusion
GetVersionNumber 66
GetCreationDate time.struct_time(tm_year=1970, tm_mon=1, tm_mday=1, tm_hour=12, tm_min=6, tm_sec=33, tm_wday=3, tm_yday=1, tm_isdst=0)
IsError False
IsNewFile False
IsThereMSData True
HasExpMethod True
InAcquisition False
GetErrorCode 0
GetErrorMessage 
GetWarningMessage 
RefreshViewOfFile None
GetNumberOfControllers 2
GetNumberOfControllersOfType('No device') 0
GetNumberOfControllersOfType('MS') 1
GetNumberOfControllersOfType('Analog') 0
GetNumberOfControllersOfType('A/D card') 0
GetNumberOfControllersOfType('PDA') 0
GetNumberOfControllersOfType('UV') 1
GetControllerType('MS') MS
GetCurrentController() (0, 1)
GetCurrentController() (0, 1)
GetCurrentController() (0, 1)
GetExpectedRunTime() 90.0
GetMaxIntegratedIntensity() 4356860928.0
GetMaxIntensity() 0
GetInletID() 0
GetErrorFlag() 0
GetFlags() 
GetAcquisitionFileName() 
GetOperator() 
GetComment1() 
GetComment2() 
GetFilters() 
GetMassTolerance() (False, 500.0, 0)
rawfile.SetMassTolerance(userDefined=True, massTolerance=555.0, units=2) None
GetMassTolerance() (True, 555.0, 2)
rawfile.SetMassTolerance(userDefined=False, massTolerance=500.0, units=0) None
GetMassResolution 0.5
GetNumTrailerExtra 46125
GetLowMass 90.0
GetHighMass 2000.0
GetStartTime 0.0011264189333333335
GetEndTime 90.0013702568
GetNumSpectra 46125
GetFirstSpectrumNumber 1
GetLastSpectrumNumber 46125
GetAcquisitionDate 
GetUniqueCompoundNames ()
############################################## INSTRUMENT BEGIN
GetInstrumentDescription C#000
GetInstrumentID 0
GetInstSerialNumber FSN10176
GetInstName Orbitrap Fusion
GetInstModel Orbitrap Fusion
GetInstSoftwareVersion 3.0.2041
GetInstHardwareVersion 
GetInstFlags 
GetInstNumChannelLabels 0
IsQExactive False
############################################## INSTRUMENT END
############################################## XCALIBUR INTERFACE BEGIN
GetScanHeaderInfoForScanNum OrderedDict([('numPackets', 61304), ('StartTime', 0.0011264189333333335), ('LowMass', 500.0), ('HighMass', 2000.0), ('TIC', 5333820.5), ('BasePeakMass', 536.1658935546875), ('BasePeakIntensity', 242910.546875), ('numChannels', 0), ('uniformTime', 0), ('Frequency', 0.0)])
GetTrailerExtraForScanNum OrderedDict([('', '\t\t\t\t\t\t\t\t\t         \t\t\t\t\t\t\t\t\t\t\t\t\t\t'), ('Scan Description', '                '), ('AGC', 'On          '), ('Micro Scan Count', 1.0), ('Ion Injection Time (ms)', 50.0), ('Elapsed Scan Time (sec)', 0.324), ('Average Scan by Inst', 'No'), ('Orbitrap Resolution', 120000.0), ('Access ID', 0.0), ('API Process Delay', -1.0), ('Dependency Type', 0.0), ('Multi Inject Info', '                                '), ('Master Scan Number', -1.0), ('Monoisotopic M/Z', -1.0), ('Charge State', 0.0), ('HCD Energy', -1.0), ('MS2 Isolation Width', -1.0), ('SPS Masses', '                                                                                                                                '), ('SPS Masses Continued', '                                                                                                                                '), ('Conversion Parameter I', 0.0), ('Conversion Parameter A', 0.0), ('Conversion Parameter B', 211820001.740694), ('Conversion Parameter C', -33697650.20725), ('Conversion Parameter D', 0.0), ('Conversion Parameter E', 0.0), ('Temperature Comp. (ppm)', -5.12), ('RF Comp. (ppm)', 0.0), ('Space Charge Comp. (ppm)', -1.27), ('Resolution Comp. (ppm)', -0.05), ('Number of LM Found', 0.0), ('LM Correction (ppm)', 0.0), ('RawOvFtT', 296606.2), ('Injection t0', -0.063), ('Reagent Ion Injection Time (ms)', 0.0)])
GetNumTuneData 1
GetTuneData(0) 
GetNumInstMethods 2
GetInstMethodNames ('TNG-Calcium', 'Dionex Chromatography MS Link')
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
rawfile.ExtractInstMethodFromRaw None
GetVialNumber 0
GetInjectionVolume 0.0
GetInjectionAmountUnits 
GetSampleVolume 0.0
GetSampleVolumeUnits 
GetSampleWeight 0.0
GetSampleAmountUnits 
GetSeqRowNumber 8
GetSeqRowSampleType Unknown
GetSeqRowDataPath 
GetSeqRowRawFileName F:\IBC_LAB603\Hsiao\20190508\zf_sPerMeOG_intestine.raw
GetSeqRowSampleName 
GetSeqRowSampleID RA1
GetSeqRowComment 
GetSeqRowLevelName 
GetSeqRowUserText 
GetSeqRowUserText 
GetSeqRowUserText 
GetSeqRowUserText 
GetSeqRowUserText 
GetSeqRowInstrumentMethod F:\IBC_LAB603\Hsiao\20190508\ON_Glycan-90min_HCDpdCID_zfOG.meth
GetSeqRowProcessingMethod 
GetSeqRowCalibrationFile 
GetSeqRowVial BD5
GetSeqRowInjectionVolume 2.0
GetSeqRowSampleWeight 0.0
GetSeqRowSampleVolume 0.0
GetSeqRowISTDAmount 0.0
GetSeqRowDilutionFactor 1.0
GetSeqRowUserLabel Study
GetSeqRowUserLabel Client
GetSeqRowUserLabel Laboratory
GetSeqRowUserLabel Company
GetSeqRowUserLabel Phone
GetSeqRowUserTextEx 
GetSeqRowUserTextEx 
GetSeqRowUserTextEx 
GetSeqRowUserTextEx 
GetSeqRowUserTextEx 
GetSeqRowBarcode 
GetSeqRowBarcodeStatus 0
GetNumStatusLog 2838
GetStatusLogForScanNum
(2.321666592308702e-08, [('RFC1 Detected RF', '1.19677'), ('RFC1 Modulation', '1.06126'), ('RFC1 Mass DAC', '1.19371188'), ('RFC1 DAC Error', '-0.00229010'), ('RFC1 5V Reference', '4.99'), ('RFC1 Diode A Temp', '59.4'), ('RFC1 Diode B Temp', '59.9'), ('RFC1 Ambient Temp', '45.9'), ('RFC2 Detected RF', '0.00778'), ('RFC2 Modulation', '-0.51225'), ('RFC2 Mass DAC', '0.00076171'), ('RFC2 DAC Error', '0.08316042'), ('RFC2 5V Reference', '5.00'), ('RFC2 Diode A Temperature', '59.9'), ('RFC2 Diode B Temperature', '59.5'), ('RFC2 Temperature', '38.5'), ('RFA1 Modulation', '1.1'), ('RFA1 Temperature', '43.4'), ('RFA1 +48 V PS', '47.9'), ('RFA1 + 15 V PS', '15.1'), ('RFA1 - 15 V PS', '-14.6'), ('RFA1 +48 V Current', '3.11'), ('RFA1 RF Voltage', '19.1'), ('RFA1 RF Current', '0.24'), ('RFA1 VSWR', '0.1'), ('RFA1 SWR Fault', '1.0'), ('RFA1 MOSFET 1 Bias', '3.8'), ('RFA1 MOSFET 2 Bias', '3.6'), ('RFA1 MOSFET 3 Bias', '3.7'), ('RFA1 MOSFET 4 Bias', '3.6'), ('RFA2 Modulation', '-0.5'), ('RFA2 Temperature', '40.2'), ('RFA2 +48 V PS', '48.2'), ('RFA2 + 15 V PS', '15.1'), ('RFA2 - 15 V PS', '-14.8'), ('RFA2 +48 V Current', '1.92'), ('RFA2 RF Voltage', '1.5'), ('RFA2 RF Current', '0.1'), ('RFA2 VSWR', '0.0'), ('RFA2 SWR Fault', '0.0'), ('RFA2 MOSFET 1 Bias', '3.8'), ('RFA2 MOSFET 2 Bias', '3.8'), ('RFA2 MOSFET 3 Bias', '3.8'), ('RFA2 MOSFET 4 Bias', '3.8'), ('+24V', '23.7'), ('+15V', '15.0'), ('+5V', '5.2'), ('-15V', '-15.0'), ('X Rods', '-27.5'), ('Y Rods', '25.5'), ('RD +1000 V PS', '1002.8'), ('RD –1000 V PS', '-1001.7'), ('RD + 24 V Current', '0.14'), ('RD Temperature', '38.6'), ('RD Detected RF', '-1.2'), ('SRIG RF FB Loop', '0.0'), ('SRIG RF', '1.1'), ('Sweep Gas', '247.0'), ('Aux Gas P', '179.7'), ('Sheath Gas', '179.8'), ('Immediate CID P', '-10.9'), ('Sweep Gas', '0.5'), ('Aux Gas', '0.1'), ('Sheath Gas', '0.1'), ('Ion Gauge Pressure', '5.03e-05'), ('Collision Pressure', '0.008'), ('Source Pressure', '1.649'), ('Ion Transfer Tube Temp', '275.1'), ('Vaporizer Temp', '0.0'), ('Dynode', '-12144'), ('Multiplier', '-1892'), ('Spray Current', '0.0'), ('Spray Voltage', '1671.4'), ('Gate', '-11.8'), ('G2B', '-11.8'), ('G1A', '-65.7'), ('S-lens DC', '26.2'), ('MP0 Drag', '-43.8'), ('TK1', '5.0'), ('L0', '8.3'), ('MP1', '4.4'), ('MP0', '5.8'), ('MP00', '11.1'), ('HP Front Section', '0.0'), ('LP Front Section', '-0.0'), ('LP Back Section', '-0.0'), ('TK2', '-63.0'), ('L2', '-188.0'), ('LP Center Section', '-23.0'), ('HP Back Section', '0.0'), ('HP Center Section', '-15.0'), ('Front Lens', '-3.7'), ('Center Lens', '-1.1'), ('Back Lens', '35.0'), ('L32', '-50.6'), ('L31', '42.4'), ('MP3', '-4.0'), ('HCD Drag', '63.9'), ('HCD Offset', '5.8'), ('IGRF3 Detected RF', '810.6'), ('IGRF3 +24 V Current', '0.23'), ('IGRF3 Command Amp', '823.5'), ('IGRF3 Modulation', '-2.8'), ('IGRF3 +24 V PS', '23.6'), ('IGRF2 Detected RF', '109.2'), ('IGRF2 +24 V Current', '0.1'), ('IGRF2 Command Amp', '-63.3'), ('IGRF2 Modulation', '-0.0'), ('IGRF2 +24 V PS', '23.7'), ('IGRF1 Detected RF', '8.4'), ('IGRF1 +24 V Current', '0.1'), ('IGRF1 Command Amp', '3.8'), ('IGRF1 Modulation', '0.0'), ('IGRF1 +24 V PS', '23.7'), ('WFG + 10 V Ref. Volatage', '10.0'), ('WFG - 2.5 V Ref. Volatage', '-2.6'), ('WFG + 2.5 V Ref. Volatage', '2.5'), ('WFG + 2.5 V PS', '2.5'), ('WFG + 3.3 V PS', '3.3'), ('WFG + 5 V DAC', '5.0'), ('WFG + 5 V ADC', '5.0'), ('WFG + 5 V PS', '5.1'), ('WFG Amplitude', '-0.0'), ('WFG Temperature', '27.3'), ('WFG + 15 V PS', '14.7'), ('WFG + 48 V PS', '47.6'), ('WFG - 15 V PS', '-15.1'), ('FT ADAP FPGA Version', '11'), ('FT Handshake Trigger', '0'), ('FT Mirror', '59402'), ('FT OSPI Board Detect', '33685504'), ('FT ADAP SW Version', '66103'), ('FT Timer', '-397.66742'), ('FT Trigger Timeout', '0'), ('FT Fan Current', '3.2'), ('FT DAQ FPGA Version', '11'), ('FT - 5 V PS', '-5.0'), ('FT Oven Alarm Voltage', '5.0'), ('FT + 15 V PS', '15.0'), ('FT + 1 V PS', '1.2'), ('FT + 3 V PS', '3.2'), ('FT DAQ + 5V PS', '5.0'), ('FT Missing Trigger', '0'), ('FT Fan Fail', '0'), ('FT PLL Lock', '1.0'), ('Gate Lens', '-59.3'), ('Trap Lens', '-59.4'), ('C-Trap Offset (Inject)', '2302.9'), ('C-Trap Offset (Pass)', '0.1'), ('C-Trap Push', '388.3'), ('C-Trap Pull', '-420.9'), ('C-Trap Detected RF', '2386.1'), ('FT RF Current', '0.4'), ('FT RF Frequency', '3.2570'), ('FT RF Overload', '0'), ('CE Inject (Cation)', '-3713.5'), ('CE Inject (Anion)', '3726.4'), ('CE Measure (Cation)', '-5001.1'), ('CE Measure (Anion)', '5001.0'), ('DE Inject', '2.7'), ('DE Measure', '728.8'), ('FT Lens 6', '1185.0'), ('FT Lens 3', '-358.5'), ('Analyzer Temp', '29.6'), ('FT Analyzer Overload', '0.0'), ('FT CE Temp (Neg)', '38.8'), ('FFT CE Temp (Pos)', '39.0'), ('FT HV PS', '5538.4'), ('FT Opt Overload', '0'), ('FT Peltier Temperature', '37.0'), ('FT CPU Load', '1088.0'), ('FT PC Free RAM', '1153.4'), ('UHV Pressure', '1.2e-10'), ('FT Vacuum SW Version', '65812'), ('FT Vacuum Up Time', '49.17'), ('FT CPU Throttle Events', '0'), ('CPU Board Temperature', '30'), ('CPU Die Temperature', '42'), ('CPU Free Ram', '690'), ('ETD - 15V PS', '-14.9'), ('ETD + 24 V PS', '23.7'), ('ETD + 15V PS', '14.9'), ('ETD + 5V PS', '5.2'), ('ETD - 2kV PS', '-2027.3'), ('Discharge Current', '0.00'), ('Discharge Voltage', '-1.9'), ('ETD HS Temperature', '1818.2'), ('Oven 1 Temp', '75.0'), ('ETD Oven 1 Current', '0.12'), ('ETD Oven 1 Voltage', '13.8'), ('ETD Oven 2 Temperature', '1818.2'), ('ETD Oven 2 Current', '-0.00'), ('ETD Oven 2 Voltage', '0.2'), ('Split Temp', '130.1'), ('ETD Split Current', '2.26'), ('ETD Split Voltage', '13.4'), ('ETD RIS Voltage', '18.9'), ('RIS Heater Current', '0.11'), ('Valve 1', '15.1'), ('ETD Valve 2', '80.0'), ('TURBO PUMP 1', ''), ('Status', 'Running'), ('Speed (Hz)', '800.0'), ('Temperature (°C)', '41.00'), ('Life Time (hours)', '45825'), ('Power (Watts)', '71'), ('TURBO PUMP 2', ''), ('Status', 'Running'), ('Speed (Hz)', '1200.0'), ('Temperature (°C)', '31.00'), ('Life Time (hours)', '45670'), ('Power (Watts)', '17'), ('TURBO PUMP 3', ''), ('Status', 'Connected'), ('Speed (Hz)', '999.0'), ('Temperature (°C)', '35.00'), ('Life Time (hours)', '45694'), ('Power (Watts)', '12')])
GetStatusLogForPos(position=0) 
GetStatusLogForPos(position=1) 
GetStatusLogPlottableIndex() 
GetNumErrorLog 25
GetErrorLogItem 0 ('Possible spray instability was detected in the full scan near 0.253 min, scan number 38 \n', 0.841193675994873)
GetErrorLogItem 1 ('Possible spray instability was detected in the full scan near 1.688 min, scan number 250 \n', 1.692299485206604)
GetErrorLogItem 2 ('Possible spray instability was detected in the full scan near 2.609 min, scan number 384 \n', 3.44126033782959)
GetErrorLogItem 3 ('Possible spray instability was detected in the full scan near 4.057 min, scan number 581 \n', 4.325737953186035)
GetErrorLogItem 4 ('Possible spray instability was detected in the full scan near 4.877 min, scan number 694 \n', 5.191875457763672)
GetErrorLogItem 5 ('Possible spray instability was detected in the full scan near 5.795 min, scan number 825 \n', 6.058018207550049)
GetErrorLogItem 6 ('Possible spray instability was detected in the full scan near 6.318 min, scan number 900 \n', 6.9258246421813965)
GetErrorLogItem 7 ('Possible spray instability was detected in the full scan near 11.299 min, scan number 1591 \n', 11.374929428100586)
GetErrorLogItem 8 ('Possible spray instability was detected in the full scan near 12.026 min, scan number 1693 \n', 12.224346160888672)
GetErrorLogItem 9 ('Possible spray instability was detected in the full scan near 12.342 min, scan number 1738 \n', 13.108858108520508)
GetErrorLogItem 10 ('Possible spray instability was detected in the full scan near 14.736 min, scan number 2085 \n', 14.89116382598877)
GetErrorLogItem 11 ('Possible spray instability was detected in the full scan near 15.253 min, scan number 2176 \n', 15.775678634643555)
GetErrorLogItem 12 ('Possible spray instability was detected in the full scan near 75.767 min, scan number 41596 \n', 76.19174194335938)
GetErrorLogItem 13 ('Possible spray instability was detected in the full scan near 78.194 min, scan number 42599 \n', 78.1911392211914)
GetErrorLogItem 14 ('Possible spray instability was detected in the full scan near 79.714 min, scan number 43000 \n', 80.04195404052734)
GetErrorLogItem 15 ('Possible spray instability was detected in the full scan near 80.962 min, scan number 43372 \n', 80.95814514160156)
GetErrorLogItem 16 ('Possible spray instability was detected in the full scan near 81.311 min, scan number 43487 \n', 81.87438201904297)
GetErrorLogItem 17 ('Possible spray instability was detected in the full scan near 82.563 min, scan number 43783 \n', 82.7922592163086)
GetErrorLogItem 18 ('Possible spray instability was detected in the full scan near 83.700 min, scan number 44084 \n', 83.70851135253906)
GetErrorLogItem 19 ('Possible spray instability was detected in the full scan near 84.064 min, scan number 44171 \n', 84.64144897460938)
GetErrorLogItem 20 ('Possible spray instability was detected in the full scan near 85.458 min, scan number 44578 \n', 85.5743637084961)
GetErrorLogItem 21 ('Possible spray instability was detected in the full scan near 86.501 min, scan number 44930 \n', 86.50894927978516)
GetErrorLogItem 22 ('Possible spray instability was detected in the full scan near 87.227 min, scan number 45182 \n', 87.45853424072266)
GetErrorLogItem 23 ('Possible spray instability was detected in the full scan near 87.606 min, scan number 45316 \n', 88.39146423339844)
GetErrorLogItem 24 ('Possible spray instability was detected in the full scan near 88.901 min, scan number 45773 \n', 89.37445068359375)
############################################## XCALIBUR INTERFACE END
GetMassListFromScanNum 
GetMassListRangeFromScanNum 
GetSegmentedMassListFromScanNum 
GetAverageMassList 
GetAveragedMassSpectrum 
GetSummedMassSpectrum 
GetLabelData 
GetAveragedLabelData 
GetAllMSOrderData (AllMSOrderData_Labels(mass=(), intensity=(), resolution=(), baseline=(), noise=(), charge=()), AllMSOrderData_Flags(activation_type=(), is_precursor_range_valid=()), 0)
GetChroData 
GetFullMSOrderPrecursorDataFromScanNum(scan_number,0) FullMSOrderPrecursorData(precursorMass=50.0, isolationWidth=1.0, collisionEnergy=25.0, collisionEnergyValid=5e-324, rangeIsValid=0.0, firstPrecursorMass=0.0, lastPrecursorMass=0.0, isolationWidthOffset=-9.4060009631604e-311)
GetFullMSOrderPrecursorDataFromScanNum(scan_number,1) FullMSOrderPrecursorData(precursorMass=50.0, isolationWidth=1.0, collisionEnergy=25.0, collisionEnergyValid=5e-324, rangeIsValid=0.0, firstPrecursorMass=0.0, lastPrecursorMass=0.0, isolationWidthOffset=-1.65401491286514e-310)
GetPrecursorInfoFromScanNum(scan_number,1) None
