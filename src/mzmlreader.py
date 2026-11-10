import os
import numpy as np
version = 0.94
last_update = 20261026
#data-loading and saving
import glob
import os
import pathlib
from pathlib import Path
import os
import pathlib
import tempfile
import time
import json
import re
import gzip
try:
    import pymzml
    mzmlavailable = True
except:
    mzmlavailable = False
from msp_spectrum_record import (SpectrumRecord, normalize_charge, derive_mass, finish_metadata_flags,
                                 check_arrays, legacy_ms2_row, legacy_ms3_row, new_summary, tally,
                                 count_level, resolve_charge, validate_dataset_charge,
                                 PRECURSOR_METADATA_ERROR, DECODE_ERROR, ARRAYS_LENGTH_MISMATCH,
                                 CHARGE_SOURCES, CHARGE_SOURCE_REPORTED, ReaderDiagnostics, status_row,
                                 polarity_from_filter)
import msp_charge_comb

#changelog:
#v0.93 [READER-DIAG R2 (Codex round 2), 20261023]: polarity read from the filter string when the spectrum has no polarity term (Henry D5; the negative-mode mass formula then applies, CSV masses of such files change); the XML precursor fallback reads one selectedIon only (no charge borrowed from another ion); the comb matches its parent by the full spectrumRef; a CSV is named in the diagnostics only after it is written; Ctrl-C/SystemExit also finish the diagnostics; stage 'iterate' starts once the file is open; the open-failure message is redacted.
#v0.94 [READER-DIAG R3 (Codex round 3), 20261026]: R3-F3: the precursor m/z, intensity and charge always come from the first selectedIon of the first precursor in the spectrum XML, and the parent from that precursor, also when pymzml's selected_precursors is filled (pymzml pairs these by position across all ions of all precursors, so one ion's m/z could get another ion's or another precursor's charge); pymzml's dicts are used only when the spectrum has no XML element. CSV changes (R3C-F1): a first ion that lacks a charge no longer takes one from a later ion or precursor (charge and mass stay 0); a spectrumRef whose scan number is empty gives parent 0 instead of a blank cell. Complete first ions are unchanged.
#v0.92 [READER-DIAG CP-D, 20261020]: every conversion also writes <stem>.scan_status.tsv (one row per MS2/MS3 spectrum) and <stem>.reader_diag.json/.txt (outcome, exit code, warnings, counts, coverage, failure) next to the CSVs, on success and on failure (msp_spectrum_record v0.4). Optional stem= (output name; the CLI passes the alias) and entry=. The return dict gains outcome, exit_code, warnings, diagnostics. CSVs unchanged.
#v0.91 [READER-DIAG CP-C.6, 20261017]: optional charge_comb (default off): the MS1 isotope comb (msp_charge_comb) fills MS2 charges the file does not report, before the dataset charge. Only the latest MS1 is held; the comb runs only on a profile MS1 named by the MS2's spectrumRef. summary["charge"] now has the RAW route's keys.
#v0.9 [READER-DIAG CP-C.5, 20261015]: optional dataset_charge (one user charge for the whole file, default off). It fills MS2 spectra without a valid reported charge; a reported charge is never overwritten. The source is counted in summary['ms2']['charge_sources']. The MS1 isotope comb is RAW-only for now (this reader keeps no MS1 peaks).
#v0.8 [READER-DIAG CP-C.4, 20261013]: calcproton takes the spectrum polarity; negative mode gives the [M-H]- equivalent (was the positive formula, (z-1)*2.015 Da low). Same change in msprawextractor v0.58.
#v0.8 [READER-DIAG CP-C.3, 20261013]: spectra without peaks are no longer written to the CSV (counted as 'empty' in the summary; msp_spectrum_record v0.2).
#v0.8 [READER-DIAG CP-C.2, 20261013]: MS3 header names fixed (column 2 is the parent MS2 scan, column 5 the MS3 scan; they were swapped). Same change in msprawextractor.
#v0.8 [READER-DIAG CP-A, 20261013]: MS2/MS3 spectra go through read_mzml_spectrum (shared rules in msp_spectrum_record): missing charge no longer raises KeyError and is never substituted; peaks are kept when precursor metadata is missing; undecodable spectra are counted and skipped; return dict gains 'summary'.
#v0.7 [s7-mzml-fold / B-23, 20260518]: CSV serialization fix — peaklist/peakintensity tuple-cast (matches RAW path repr; eliminates str(np.array) line-wrap that broke tab-row structure); CV-term-aware monoisotopic m/z extraction with ms1mz fallback (was hardcoded 0).
#v0.6 code review 1.
#v0.5 functional mzmlreader (proto) - may not be flexible on different parameters during extraction

# Load an mzML file
#run = pymzml.run.Reader(r"C:\Users\Sakazuki\Downloads\zf_sPerMeNG_brain.mzML")


#####functional blocks imported from msprawrxtactor.py#####
###consider isolate the metadata function to another pyhton file###
###copied 20260115, rawextractor version = 0.53, last_update = 20250927###
def fillexpinfo(mzmlfilepath, loadfile=None, debug=False):
    if not loadfile:
        print("You imported experiment description file")
    else:
        if debug:
            print("[debug]Start filling new experiment information for metadata")
    expTitle = input("Enter title of this experiment. For example 'human T cells'.\n")
    expDescription = input("Enter details of this experiment. For example 'SiaT KO test'\n")
    expAuthor = input("Enter your name so other knows who did the analysis.\n")
    expRawdate = input("Enter the date when the raw file was obtained in YYYY/MM/DD format. For example '20240229'.\n")
    expglycantype = input("Enter glycan type of analysis. (N/O/GL/N+O/others).\n")
    expmsmode = input("Enter mass spec detection mode (+/-)")
    #expmsinstrument = input("Enter mass spec instrument name")

    #autofill
    parameters = ["[debug]Autofill the version info"]
    extractdate = time.strftime("%Y%m%d") #20001105 for example
    rawfile = mzmlfilepath #filename.raw (str)
    rawfilename =  os.path.splitext(mzmlfilepath)[0] # filename (str)

    #define metadata
    metadata = {
        "json_type": "glycomsp.metadata",
        "schema_version": "1.0.0",
        #above are fixes for next version json update
        "Experiment Title": expTitle,
        "Experiment Description": expDescription,
        "Author running this analysis": expAuthor,
        "Raw data acquired date": expRawdate,
        "Glycan Type": expglycantype,
        "Mass Analyzer charge mode": expmsmode,
        "Parameters when GlycoMSP launched":parameters,
        "Date of file extracted from raw file": extractdate,
        "Original raw file path":rawfile,
        "Raw filename":rawfilename
        }
    return metadata
    
    #define project name for saving metadata and log file #
def generatefilename(rawfilepath):
    mansavename = input("Enter the name you want to save for this project or batch of analysis.\n")
    rawfilename =  os.path.basename(rawfilepath) # filename (str) #os.path.splitext(rawfilepath)[0]
    extractdate = time.strftime("%Y%m%d") #20001105 for example
    savename = f"{mansavename}_{extractdate}_{rawfilename}"
    return savename

def savemetadata(metadata, filepath):
    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=4)
    return filepath

def save_log(filepath, log_entries, debug_logs=None):
    with open(filepath, 'w', encoding='utf-8') as f:
        for entry in log_entries:
            f.write(f"{entry}\n")
        if debug_logs:
            f.write("-----\ndebug enabled-----\n")
            for logs in debug_logs:
                f.write(f"{logs}\n")
            f.write("-----\ndebug log ends-----\n")
    return filepath

#v0.8 CP-C.4 READER-DIAG: negative-mode precursor mass. The in-silico library (compnewv4, mode '-') stores [M-H]-,
#and Thermo reports the charge magnitude in negative mode too, so the scan polarity decides the formula.
PROTON_MASS = 1.007276

def calcproton(isolatedmass, charge, polarity=None):
    if polarity == '-' or charge < 0:
        #[M-zH]z- -> [M-H]- equivalent: m/z * z + (z - 1) * proton. A negative charge value is read as negative mode.
        z = abs(charge)
        if z == 0:
            return 0
        return isolatedmass * z + (z - 1) * PROTON_MASS
    if charge == 1 :   #change abs(charge) to charge to negative can be calculated too, 
        #A tiny nuance to keep in mind as we go forward: Thermo’s API typically returns the magnitude of charge (often positive) 
        # even for negative-mode data, depending on method/instrument
        return isolatedmass
    elif charge >= 2:
        conv = (isolatedmass * charge - (charge - 1) * 1.00784)
        return conv
    elif charge == 0:
        return 0

try:
    import xml.etree.ElementTree as ET
    def mzml_spectrum_count(mzml_path: str) -> int | None:
        # iterparse is streaming; we exit early once spectrumList is found
        for event, elem in ET.iterparse(mzml_path, events=("start",)):
            # elem.tag may include a namespace like {http://psi.hupo.org/ms/mzml}spectrumList
            if elem.tag.endswith("spectrumList"):
                cnt = elem.attrib.get("count")
                return int(cnt) if cnt is not None else None
        return None
except:
    def mzml_spectrum_count(mzml_path):
        return None

def mzml_file_examination(mzmlfile, auto = True, time = 0.0,debug=False):
    if not auto:
        print(f'This raw file has {mzml_spectrum_count(mzmlfile)} spectra')
        scan_number = int(input('enter spectrum no you want to summarize'))
        maxn = int(mzml_spectrum_count(mzmlfile))
    else:
        scan_number = int(mzml_spectrum_count(mzmlfile))
        maxn = scan_number
    return maxn 


def peaklist_validation(peaklist, peakintensity, debug = False):

    if len(peaklist) == len(peakintensity):
        return True, None
    else:
        if debug:
            peaklist_length = [len(peaklist), len(peakintensity)]
            print(f"[debug] peaklist and intensity are not matched, {peaklist_length}")
            return False, peaklist_length
        return False, None

# 20260518 [B-23] Walk precursor XML for PSI/Thermo monoisotopic-m/z, else fall back.
# Why: many converters (incl. the one that produced our smoke #2a fixture) drop the
# monoisotopic CV term; hardcoding 0 silently zeros the column. Fallback to ms1mz
# matches RAW-path semantics for the common case monoisotopic == isolation.
def _extract_monoiso_mz(precursor_dict, fallback_mz):
    element = precursor_dict.get('element') if isinstance(precursor_dict, dict) else None
    if element is not None:
        for sub in element.iter():
            tag = sub.tag.rsplit('}', 1)[-1]
            if tag == 'cvParam' and sub.attrib.get('accession') == 'MS:1002817':
                try:
                    return float(sub.attrib.get('value'))
                except (TypeError, ValueError):
                    pass
            elif tag == 'userParam':
                name = (sub.attrib.get('name') or '').lower()
                if 'monoisotopic' in name and 'm/z' in name:
                    try:
                        return float(sub.attrib.get('value'))
                    except (TypeError, ValueError):
                        pass
    try:
        return float(fallback_mz)
    except (TypeError, ValueError):
        return 0

def save_to_csv(csvname, data, debug=False):
    """
    Atomically write a tab-separated CSV:
    - Write to a temp file in the same directory, then os.replace to final.
    - Keeps the legacy behavior of using a *stem* (no .csv) for csvname.
    """
    exportcsv = f"{csvname}.csv"
    parent = os.path.dirname(exportcsv) or "."
    os.makedirs(parent, exist_ok=True)

    perf_esti2 = time.time()
    # Create a temp file in the same directory so os.replace is atomic on the same volume
    with tempfile.NamedTemporaryFile("w", delete=False, dir=parent, suffix=".tmp") as tmp:
        # Write rows as tab-separated (matches legacy)
        # NOTE: data is a list of row-like items, first row may be header tuple
        for row in data:
            tmp.write('\t'.join(map(str, row)) + '\n')
        tmp.flush()
        os.fsync(tmp.fileno())
        tmp_path = tmp.name

    # Atomic promotion to final name (no second write)
    os.replace(tmp_path, exportcsv)

    if debug:
        csvwritetime = time.time() - perf_esti2
        print(f"[debug] Finished exported {exportcsv} from mzml data")
        print(f"[debug] Time spent for extraction: {csvwritetime:.3f} seconds.")
    return exportcsv

#############################
#########end of copy#########

def _mzml_cv_choice(element, choices):
    """Return the label of the first accession in `choices` present on the spectrum, else None."""
    if element is None:
        return None
    present = {sub.attrib.get('accession') for sub in element.iter() if sub.tag.rsplit('}', 1)[-1] == 'cvParam'}
    for accession, label in choices:
        if accession in present:
            return label
    return None

def _mzml_cv_value(element, accession):
    """The value of the first cvParam with this accession on the spectrum, else None."""
    if element is None:
        return None
    for sub in element.iter():
        if sub.tag.rsplit('}', 1)[-1] == 'cvParam' and sub.attrib.get('accession') == accession:
            return sub.attrib.get('value')
    return None

#READER-DIAG CP-C.6 (20261017): the comb on the mzML route, same order and status strings as msprawextractor v0.59.
#Only the latest MS1 is held (data-dependent runs: the parent precedes its MS2s). The comb runs only when the MS2's
#spectrumRef names that MS1, the MS1 is a profile spectrum, and its Thermo filter string (if any) is not FTMS.
COMB_NOT_RUN = "not_run"                   # charge reported, or the comb is off
COMB_NO_ISOLATION_MZ = "no_isolation_mz"
COMB_NO_PARENT_MS1 = "no_parent_ms1"       # no spectrumRef, or it names an MS1 other than the latest one
COMB_NO_PROFILE_MS1 = "no_profile_ms1"     # parent MS1 stored as centroids, or not stated
COMB_FT_MS1_NOT_TRIED = "ft_ms1_not_tried"
COMB_MS1_NOT_DECODABLE = "ms1_not_decodable"
_FT_ANALYZER_RE = re.compile(r"\bFTMS\b")

class _MzmlChargeComb:
    """Holds the latest MS1 for the comb and keeps summary["charge"] (keys as the RAW route, without mass_mz_source)."""

    def __init__(self, enabled, dataset_charge):
        self.enabled = bool(enabled)
        self.dataset_charge = dataset_charge
        self.held_id, self.held, self.held_native = None, None, None
        self.pending = None
        self.summary = {"settings": {"charge_comb": self.enabled, "dataset_charge": dataset_charge,
                                     "comb": msp_charge_comb.parameters() if self.enabled else None},
                        "sources": {k: 0 for k in CHARGE_SOURCES}, "comb_status": {}, "comb_vs_dataset_disagree": 0,
                        "scans": []}   # one entry per exported MS2 spectrum without a valid reported charge

    def hold_ms1(self, spectrum):
        """Replace the held MS1: its (m/z, intensity) arrays, or the status that stops the comb on it."""
        if not self.enabled:
            return
        element = getattr(spectrum, 'element', None)
        self.held_id = str(spectrum.ID)
        self.held_native = element.get('id') if element is not None else None   #R2-F8: full native ID
        if _mzml_cv_choice(element, (('MS:1000127', 'centroid'), ('MS:1000128', 'profile'))) != 'profile':
            self.held = COMB_NO_PROFILE_MS1
        elif _FT_ANALYZER_RE.search(_mzml_cv_value(element, 'MS:1000512') or ''):
            self.held = COMB_FT_MS1_NOT_TRIED
        else:
            try:
                mz, intensity = spectrum.mz, spectrum.i
                self.held = (mz, intensity) if len(mz) == len(intensity) else COMB_MS1_NOT_DECODABLE
            except Exception:
                self.held = COMB_MS1_NOT_DECODABLE

    def infer(self, rec):
        """Comb charge for an MS2 record whose charge is unknown, or None. Called before the charge is resolved."""
        comb, status = None, COMB_NOT_RUN
        if self.enabled:
            if not (rec.precursor_mz or 0) > 0:
                status = COMB_NO_ISOLATION_MZ
            elif not self._is_held_parent(rec):
                status = COMB_NO_PARENT_MS1
            elif isinstance(self.held, str):
                status = self.held
            else:
                comb = msp_charge_comb.infer_charge(self.held[0], self.held[1], rec.precursor_mz)
                status = comb.status
        self.pending = (rec.scan, comb, status)
        return comb.z if comb is not None else None

    def _is_held_parent(self, rec):
        """True when the MS2's parent is the held MS1: by the full spectrumRef when both IDs exist (R2-F8: two
        controllers can share a scan number), else by pymzml's reduced scan number as before."""
        if rec.parent_scan is None or self.held_id is None:
            return False
        if rec.parent_ref and self.held_native:
            return rec.parent_ref == self.held_native
        return str(rec.parent_scan) == self.held_id

    def record(self, rec):
        """Count one exported MS2 record after its charge was resolved."""
        comb, status = None, COMB_NOT_RUN
        if self.pending is not None and self.pending[0] == rec.scan:
            comb, status = self.pending[1], self.pending[2]
        self.pending = None
        cs = self.summary
        cs["sources"][rec.charge_source] += 1
        cs["comb_status"][status] = cs["comb_status"].get(status, 0) + 1
        comb_z = comb.z if comb is not None else None
        if comb_z is not None and self.dataset_charge is not None and comb_z != self.dataset_charge:
            cs["comb_vs_dataset_disagree"] += 1
        if rec.charge_source != CHARGE_SOURCE_REPORTED:
            entry = {"scan": rec.scan, "charge_source": rec.charge_source, "charge": rec.charge or 0,
                     "comb_status": status}
            if comb is not None and comb.score is not None:
                entry.update(comb_z=comb_z, comb_best_z=max(comb.scores, key=lambda z: comb.scores[z][0]),
                             comb_score=round(comb.score, 4), comb_margin=round(comb.margin, 4),
                             comb_anchor=round(comb.anchor, 4))
            cs["scans"].append(entry)

def _first_precursor(element):
    """The first <precursor> element of a spectrum, or None."""
    if element is None:
        return None
    for sub in element.iter():
        if sub.tag.rsplit('}', 1)[-1] == 'precursor':
            return sub
    return None

def _precursor_from_element(element):
    """Read the first precursor straight from the spectrum XML.

    The precursor source for every spectrum with an XML element (R3-F3, 20261026). pymzml's selected_precursors
    collects the m/z, intensity and charge terms of all selected ions into separate lists and pairs them by
    position, so an ion without a charge can get the next ion's charge; it also loses the parent reference of a
    precursor without m/z.
    Same keys as pymzml ('mz', 'i', 'charge', 'precursor id', 'element'); absent values are omitted.
    R2-F9 (20261023): m/z, intensity and charge come from the first selectedIon only, all together; a field that
    ion lacks stays absent (never borrowed from a second ion).
    """
    precursor = _first_precursor(element)
    if precursor is None:
        return {}
    found = {'element': precursor}
    ref = precursor.get('spectrumRef')
    if ref is not None:
        match = pymzml.regex_patterns.SPECTRUM_ID_PATTERN.search(ref)
        if match and match.group(1):
            found['precursor id'] = match.group(1)
    ion = None
    for sub in precursor.iter():
        if sub.tag.rsplit('}', 1)[-1] == 'selectedIon':
            ion = sub
            break
    for sub in (ion.iter() if ion is not None else ()):
        if sub.tag.rsplit('}', 1)[-1] != 'cvParam':
            continue
        key = {'MS:1000744': 'mz', 'MS:1000042': 'i', 'MS:1000041': 'charge'}.get(sub.attrib.get('accession'))
        if key and key not in found:
            try:
                found[key] = float(sub.attrib.get('value'))
            except (TypeError, ValueError):
                pass
    return found

def _round4(values):
    rounded = []
    for value in values:
        rounded.append(np.round(value, 4))
    return rounded

# 20261013 [READER-DIAG CP-A] Precursor metadata and peak arrays are read independently.
# Why: pymzml omits the 'charge' key when the selected ion has no charge term, and the old
# direct index raised KeyError for the whole file before any peaks were read.
def read_mzml_spectrum(spectrum, ms_level, round=False, dataset_charge=None, comb=None):
    """Build a SpectrumRecord for one MS2/MS3 pymzml spectrum. Never raises for missing metadata.

    For an MS2 spectrum without a reported charge, `comb` (an _MzmlChargeComb) is tried first, then
    `dataset_charge` (validated by the caller). A reported charge is never replaced.
    """
    rec = SpectrumRecord(ms_level=ms_level, scan=spectrum.ID)
    element = getattr(spectrum, 'element', None)
    if element is not None:
        rec.native_id = element.get('id')
    rec.polarity = _mzml_cv_choice(element, (('MS:1000130', '+'), ('MS:1000129', '-')))
    if rec.polarity is not None:
        rec.polarity_source = 'cv'
    else:   #v0.93 R2 D5 (Henry 20261023): no polarity term -> the Thermo filter string's polarity token, if any
        rec.polarity = polarity_from_filter(_mzml_cv_value(element, 'MS:1000512'))
        rec.polarity_source = 'filter' if rec.polarity is not None else None
    rec.representation = _mzml_cv_choice(element, (('MS:1000127', 'centroid'), ('MS:1000128', 'profile')))

    try:
        precursors = spectrum.selected_precursors or []
    except Exception as e:
        precursors = []
        rec.flag(PRECURSOR_METADATA_ERROR)
        rec.note_error("precursor", e)
    # v0.94 R3-F3: one selectedIon's fields together; pymzml's dicts only when there is no XML element
    precursor = _precursor_from_element(element) if element is not None else (precursors[0] if precursors else {})
    rec.precursor_mz = precursor.get('mz')
    rec.precursor_intensity = precursor.get('i')
    rec.parent_scan = precursor.get('precursor id')
    first = _first_precursor(element)
    rec.parent_ref = first.get('spectrumRef') if first is not None else None   #R2-F8
    rec.charge = normalize_charge(precursor.get('charge'))
    if ms_level == 2:
        comb_z = comb.infer(rec) if (comb is not None and rec.charge is None) else None
        rec.charge, rec.charge_source = resolve_charge(rec.charge, comb_z, dataset_charge)
    if rec.precursor_mz is not None:
        rec.monoiso_mz = _extract_monoiso_mz(precursor, rec.precursor_mz)  # 20260518 [B-23] CV/userParam-aware
    rec.derived_mass = derive_mass(rec.precursor_mz, rec.charge,
                                   lambda mz, z: calcproton(mz, z, polarity=rec.polarity))
    finish_metadata_flags(rec)

    try:
        p_peaklist = spectrum.mz
        p_intensity = spectrum.i
        if round:
            print("[dev]Perform round to 4 digits in data extraction")
            print("[dev]Need to save log of this operation")
            p_peaklist = _round4(p_peaklist)
            print(f"peaks after pre-processing rounding is {p_peaklist}")
            p_intensity = _round4(p_intensity)
            print(f"intensity after pre-processing rounding is {p_intensity}")
        # 20260518 [B-23] cast to tuple of py floats; str(np.array) wraps at 75 chars and breaks tab rows
        rec.mz = tuple(float(x) for x in p_peaklist)
        rec.intensity = tuple(float(x) for x in p_intensity)
    except Exception as e:
        rec.mz = rec.intensity = None
        rec.flag(DECODE_ERROR)
        rec.note_error("decode", e)
    check_arrays(rec)
    return rec

#READER-DIAG CP-D (20261020): instrument terms for the diagnostics, read from the mzML header only. cvParam names
#and software versions; never cvParam values, the run id or the sourceFile names and locations.
_HEADER_SECTIONS = {"referenceableParamGroupList": "instrument", "instrumentConfigurationList": "instrument",
                    "softwareList": "software", "dataProcessingList": "processing"}

def _mzml_header_info(mzmlfilepath, limit=4 * 1024 * 1024):
    try:
        opener = gzip.open if str(mzmlfilepath).lower().endswith(".gz") else open
        with opener(mzmlfilepath, "rb") as fh:
            head = fh.read(limit).decode("utf-8", "replace")
    except Exception:
        return {}
    head = head.split("<run", 1)[0]
    info = {"instrument": [], "software": [], "processing": []}
    for section, key in _HEADER_SECTIONS.items():
        for block in re.findall(r"<%s\b.*?</%s>" % (section, section), head, re.S):
            if key == "software":
                for attrs, body in re.findall(r"<software\b([^>]*)>(.*?)</software>", block, re.S):
                    ver = re.search(r'\bversion="([^"]*)"', attrs)
                    for name in re.findall(r'<cvParam\b[^>]*\bname="([^"]*)"', body):
                        entry = {"name": name, "version": ver.group(1) if ver else None}
                        if entry not in info["software"]:
                            info["software"].append(entry)
            else:
                for name in re.findall(r'<cvParam\b[^>]*\bname="([^"]*)"', block):
                    if name not in info[key]:
                        info[key].append(name)
    return info

def _backend_info():
    return {"name": "pymzml", "version": getattr(pymzml, "__version__", None) if mzmlavailable else None}

def extract_mzML(mzmlfilepath=None, outdir=None, round=False, mzmlavailable=False, debug = False, dataset_charge=None,
                 charge_comb=False, stem=None, entry=None):
    """Convert one mzML file to the MS2/MS3 temp CSVs; also writes the status file and diagnostics files (CP-D).

    `stem` names every output (default: the input file stem, as before); `entry` is recorded in the diagnostics
    (the CLI passes "cli-convert"). Returns the legacy dict plus outcome, exit_code, warnings and diagnostics
    paths, or False when the file cannot be opened. An exception during conversion is recorded and re-raised.
    """
    #v0.9 CP-C.5: one user charge for the whole file; fills only spectra without a valid reported charge
    #v0.91 CP-C.6: charge_comb tries the MS1 isotope comb first (profile MS1 only)
    dataset_charge = validate_dataset_charge(dataset_charge)
    if not outdir:
        outdir = os.path.dirname(mzmlfilepath or "") or os.getcwd()
    stem = stem or (pathlib.Path(mzmlfilepath).stem if mzmlfilepath else "unnamed")
    comb = _MzmlChargeComb(charge_comb, dataset_charge)
    diag = ReaderDiagnostics("mzmlreader", version, last_update, mzmlfilepath, stem, "mzml", entry=entry,
                             settings={"round": bool(round), **comb.summary["settings"]}, backend=_backend_info())
    try:
        result = _extract_mzML(mzmlfilepath, outdir, round, debug, dataset_charge, comb, stem, diag)
    except BaseException as exc:   #v0.93 R2-F6: Ctrl-C / SystemExit finish the diagnostics too; re-raised unchanged
        diag.charge = comb.summary
        diag.finish(outdir, exc)
        raise
    diag.charge = comb.summary
    done = diag.finish(outdir)
    if isinstance(result, dict):
        result.update(done)
    return result

def _extract_mzML(mzmlfilepath, outdir, round, debug, dataset_charge, comb, stem, diag):
    if debug:
        print("[debug] need to create debug logs when debug is set to True.")
        print("[debug] Also need a readable log that records needed information and the operation can be reproduced when loading that record.")
        print("[debug add debug log and record log in debug block]")
        print("[debug]: start calculate time spent for converting data")
        perf_esti1 = time.time() #first time stamp for evaluating performance
    try:
        mzmlfile = pymzml.run.Reader(mzmlfilepath)
    except Exception as e:
        print(f"[Error] MSFileReader failed: {diag.redact(e)}")   #v0.93 R2-F3: redacted
        diag.fail(e)
        return False
    diag.stage = "iterate"   #v0.93 R2-F7: only the Reader(...) call above counts as "could not be opened"
    diag.instrument = _mzml_header_info(mzmlfilepath)
    # Default to the raw file’s folder, else CWD
    if not outdir:
        outdir = os.path.dirname(mzmlfilepath) or os.getcwd()

    #init copied from msprawextractor
    scan_number = mzml_file_examination(mzmlfilepath)
    diag.declared = scan_number
    #init variables for storing parsed data
    ms1list, ms2list, ms3list, errlist , unmatchlist, logs = [], [], [], [], [], []
    #from v0.4 named Tuple of (peaklist,intensity) were changed to 2 separated lists. Include length validation 
    ms2data, ms3data = [], []  #b for MS2, c for MS3 < replacing
    fixms2header = ('entry_no','MS1scan_no', 'MS1_isolationmass', 'MS1_monoisolationmass','chargeState','protonatedmass',
           'MS2scan_no', 'peaklist', 'peakintensity')
    ms2data.append(fixms2header)
    fixms3header =  ('entry no','MS2scan_no', 'MS2_isolationmass', 'MS2_monoisomass', 'MS3scan_no', 'peaklist', 'peakintensity')
    ms3data.append(fixms3header)
    ms2count, ms3count = 1,1  
    print(f"[dev] scan number is {scan_number}")
    summary = new_summary()
    for spectrum in mzmlfile:
        ms_level = spectrum['ms level']
        diag.step(spectrum.ID, ms_level)
        if ms_level == 1:
            if debug:
                ms1list.append(spectrum.ID)
            count_level(summary, 1)
            comb.hold_ms1(spectrum)
        elif ms_level == 2:
            if debug:
                ms2list.append(spectrum.ID)
            rec = read_mzml_spectrum(spectrum, 2, round=round, dataset_charge=dataset_charge, comb=comb)
            tally(summary, rec)
            if not rec.exportable:
                comb.pending = None
                diag.add(status_row(rec))
                continue
            comb.record(rec)
            scans = comb.summary["scans"]
            diag.add(status_row(rec, scans[-1] if scans and scans[-1]["scan"] == rec.scan else None))
            ms2data.append(legacy_ms2_row(ms2count, rec))
            ms2count += 1
            peakcheck, failedlist = peaklist_validation(rec.mz, rec.intensity, debug = True)
            if not peakcheck: #if the peaklist and intensity aren't matched, not False = True
                reportunmatch = [spectrum.ID, failedlist]
                unmatchlist.append(reportunmatch)
        elif ms_level == 3:
            if debug:
                ms3list.append(spectrum.ID)
            rec = read_mzml_spectrum(spectrum, 3, round=round)
            tally(summary, rec)
            diag.add(status_row(rec))
            if not rec.exportable:
                continue
            ms3data.append(legacy_ms3_row(ms3count, rec))
            ms3count +=1
        else:
            count_level(summary, ms_level)
            errlist.append(spectrum.ID)   
    print("Extraction finished.")
    #if possible, delete the mzml in memory here
    if unmatchlist is not []:
        logs.append("Unmatched peaklists are recorded as [[spectrum no], [peaklist entries, peak intensity entries]]")
        logs.append(unmatchlist)

    if debug:
        debug_headers = "Debug enabled. Debug log information: [MS1 list, MS2 list, MS3 list, exception list]"
        debug_logs = [debug_headers, ms1list, ms2list, ms3list, errlist]
    else:
        debug_logs = None
    if debug:
        extractiontime = time.time() - perf_esti1  #timestamp1: time spent for raw file extraction
        print("[debug]mzml file conversion finished.")
        print("[debug]Preparing information for writing...")
        print(f"[debug]Time spent for extraction: {int(extractiontime)} seconds.")

    filename = os.path.splitext(os.path.basename(mzmlfilepath))[0]
    ms2output = "ms2tmp_" + filename
    ms3output = "ms3tmp_" + filename

    rawstem = stem   #v0.92 CP-D: the input file stem unless the caller names the outputs (CLI alias)
    #v0.53
    ms2tmp_path = os.path.join(outdir, f"ms2tmp_{rawstem}.csv")
    ms3tmp_path = os.path.join(outdir, f"ms3tmp_{rawstem}.csv")
    # save_to_csv still takes a *stem*; pass the stem without .csv but with full dir
    diag.stage = "write"
    save_to_csv(ms2tmp_path[:-4], ms2data, debug=debug)
    diag.wrote("ms2", ms2tmp_path)   #v0.93 R2-F2: named in the diagnostics once written
    save_to_csv(ms3tmp_path[:-4], ms3data, debug=debug)
    diag.wrote("ms3", ms3tmp_path)
    versioninfo = ["extractor info", version, last_update]
    logs.append(versioninfo)
    ms2s = summary["ms2"]
    summary["charge"] = comb.summary
    print(f"[mzml] MS2 exported {ms2s['exported']}/{ms2s['seen']} (charge unknown {ms2s['charge_unknown']}, "
          f"skipped without peaks {ms2s['empty']}, not decodable {ms2s['failed']})")
    print(f"[mzml] MS2 charge sources: {ms2s['charge_sources']} (dataset charge: {dataset_charge}); "
          f"comb: {comb.summary['comb_status']}")
    print("Conversion finished from mspextractor.")

    return {"ms2tmp": ms2tmp_path, "ms3tmp": ms3tmp_path, "summary": summary}

def finalize_extraction(rawfile, metadata, filename, debug=False):
    print(f"[Finalizing extraction] For: {filename}")
    print(f"[Metadata provided] Sample name: {metadata.get('Experiment Title', 'Unknown')}")

    # Optional: log the fact that finalization completed
    if debug:
        print("[Debug] finalize_extraction called successfully")

if mzmlavailable:
    print(f"mzml support confirmed")
elif not mzmlavailable:
    print(f"mzml not installed")

#docs are summarized by GPT4o (20241104) I can't find clear instructions on internet 
