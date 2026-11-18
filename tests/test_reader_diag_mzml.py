"""READER-DIAG CP-A: mzML reader missing-metadata states (2026-10-13).

Synthetic fixtures only (tests/_synth_mzml.py). The optional pilot test rebuilds a small mzML
from the hand-cut excerpt under handoff/reader_diag_inputs/ when that local file exists; the
excerpt and the rebuilt file never leave tmp_path.

CSV values asserted for unknown metadata (chargeState/protonatedmass 0) are the INTERIM
serialisation, pending the CP-B export contract.
"""
import csv
import os
import re
import sys

import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, '..', 'src'))
sys.path.insert(0, HERE)

pymzml = pytest.importorskip("pymzml")

import mzmlreader  # noqa: E402
import msp_spectrum_record as R  # noqa: E402
import _synth_mzml as S  # noqa: E402

PILOT_EXCERPT = os.path.join(HERE, '..', 'handoff', 'reader_diag_inputs', 'header_and_6785.mzML')

MS1 = lambda idx, scan: S.spectrum_xml(idx, scan, 1, [400.0, 500.0, 600.0], [10.0, 20.0, 30.0])


def _records(path, ms_level=2):
    out = {}
    for spectrum in pymzml.run.Reader(path):
        if spectrum['ms level'] == ms_level:
            out[spectrum.ID] = mzmlreader.read_mzml_spectrum(spectrum, ms_level)
    return out


MS3_HEADER = ['entry no', 'MS2scan_no', 'MS2_isolationmass', 'MS2_monoisomass', 'MS3scan_no',
              'peaklist', 'peakintensity']


def _rows(path):
    with open(path, newline='') as fh:
        return list(csv.reader(fh, delimiter='\t'))


@pytest.fixture
def mixed_file(tmp_path):
    spectra = [
        MS1(0, 10),
        S.spectrum_xml(1, 11, 2, [150.0, 200.5], [1.0, 2.0],
                       precursor={"parent": 10, "mz": 500.25, "charge": 2, "intensity": 99.0}),
        S.spectrum_xml(2, 12, 2, [150.0, 210.0, 300.25], [3.0, 4.0, 5.0],
                       precursor={"parent": 10, "mz": 512.5}),                       # no charge
        S.spectrum_xml(3, 13, 2, [160.0], [1.0],
                       precursor={"parent": 10, "mz": 300.0, "charge": 0}),          # charge 0
        S.spectrum_xml(4, 14, 2, [170.0], [1.0], precursor={"parent": 10}),          # no m/z
        S.spectrum_xml(5, 15, 2, [180.0], [1.0], precursor={"mz": 450.0, "charge": 2}),  # no parent
        S.spectrum_xml(6, 16, 2, [], [], precursor={"parent": 10, "mz": 520.0}),     # empty arrays
        S.spectrum_xml(7, 17, 2, [190.0], [1.0], precursor={"parent": 10, "mz": 530.0},
                       mz_binary="eJwBPADD/notzlib/AAAA"),                            # decode error
        S.spectrum_xml(8, 18, 2, [195.0, 196.0], [1.0, 2.0],
                       precursor={"parent": 10, "mz": 540.0, "charge": 3}, polarity="+", centroid=False),
    ]
    return S.write_mzml(tmp_path / "synthetic.mzML", spectra)


def test_root_cause_pymzml_omits_charge_key(mixed_file):
    """pymzml only adds 'charge' when the selected ion carries MS:1000041 (the v0.7 KeyError)."""
    for spectrum in pymzml.run.Reader(mixed_file):
        if spectrum.ID == 12:
            precursor = spectrum.selected_precursors[0]
            assert 'charge' not in precursor
            assert precursor['precursor id'] == '10'
            return
    pytest.fail("scan 12 not found")


def test_known_charge_record(mixed_file):
    rec = _records(mixed_file)[11]
    assert rec.charge == 2 and rec.polarity == '-'
    assert rec.derived_mass == pytest.approx(500.25 * 2 + 1.007276)   # negative mode: [M-H]- equivalent
    assert rec.parent_scan == '10' and rec.precursor_mz == 500.25
    assert rec.mz == (150.0, 200.5) and rec.intensity == (1.0, 2.0)
    assert rec.issues == [] and rec.mass_eligible


@pytest.mark.parametrize("scan", [12, 13])
def test_missing_or_zero_charge_stays_unknown(mixed_file, scan):
    rec = _records(mixed_file)[scan]
    assert rec.charge is None                      # never 1, -1 or a method default
    assert rec.derived_mass is None
    assert R.CHARGE_UNKNOWN in rec.issues and R.DERIVED_MASS_UNAVAILABLE in rec.issues
    assert rec.parent_scan == '10' and rec.precursor_mz is not None
    assert len(rec.mz) > 0 and rec.exportable and not rec.mass_eligible


def test_missing_precursor_mz_keeps_parent_and_peaks(mixed_file):
    rec = _records(mixed_file)[14]
    assert R.PRECURSOR_MZ_MISSING in rec.issues
    assert rec.parent_scan == '10'                 # pymzml alone would lose it
    assert rec.mz == (170.0,) and not rec.mass_eligible


def test_missing_parent_reference(mixed_file):
    rec = _records(mixed_file)[15]
    assert R.PARENT_REF_MISSING in rec.issues
    assert rec.parent_scan is None
    assert rec.charge == 2 and rec.mass_eligible   # mass does not need the parent


def test_empty_arrays_are_readable_not_errors(mixed_file):
    rec = _records(mixed_file)[16]
    assert rec.mz == () and rec.intensity == ()
    assert R.ARRAYS_EMPTY in rec.issues and R.DECODE_ERROR not in rec.issues
    assert rec.decoded and not rec.exportable       # readable, but not written (D2, 2026-10-13)


def test_decode_error_is_distinct(mixed_file):
    rec = _records(mixed_file)[17]
    assert R.DECODE_ERROR in rec.issues and R.ARRAYS_EMPTY not in rec.issues
    assert rec.mz is None and not rec.exportable
    assert rec.errors and rec.errors[0][0] == "decode"
    assert rec.precursor_mz == 530.0               # metadata still read


def test_polarity_representation_native_id(mixed_file):
    recs = _records(mixed_file)
    assert recs[11].polarity == '-' and recs[11].representation == 'centroid'
    assert recs[18].polarity == '+' and recs[18].representation == 'profile'
    assert recs[11].native_id == "controllerType=0 controllerNumber=1 scan=11"


def test_extract_mzml_mixed_file(mixed_file, tmp_path):
    res = mzmlreader.extract_mzML(mixed_file, outdir=str(tmp_path))
    rows = _rows(res["ms2tmp"])
    assert rows[0] == ['entry_no', 'MS1scan_no', 'MS1_isolationmass', 'MS1_monoisolationmass',
                       'chargeState', 'protonatedmass', 'MS2scan_no', 'peaklist', 'peakintensity']
    by_scan = {r[6]: r for r in rows[1:]}
    assert sorted(by_scan) == ['11', '12', '13', '14', '15', '18']   # 16 empty, 17 not decodable
    assert [r[0] for r in rows[1:]] == [str(i) for i in range(1, 7)]  # entry_no stays contiguous
    # Known charge: exact legacy text.
    assert by_scan['11'] == ['1', '10', '500.25', '500.25', '2', str(500.25 * 2 + 1 * 1.007276), '11',
                             '(150.0, 200.5)', '(1.0, 2.0)']                  # negative mode
    assert by_scan['18'][4:6] == ['3', str(540.0 * 3 - (3 - 1) * 1.00784)]   # positive mode, legacy formula
    # Unknown charge: interim sentinel, peaks kept.
    assert by_scan['12'][4:6] == ['0', '0'] and by_scan['12'][7] == '(150.0, 210.0, 300.25)'

    s = res["summary"]
    assert s["spectra_by_ms_level"] == {"1": 1, "2": 8}
    ms2 = s["ms2"]
    assert {k: v for k, v in ms2.items() if k not in ("issues", "errors")} == {
        "seen": 8, "exported": 6, "empty": 1, "failed": 1,
        "charge_known": 3, "charge_unknown": 3, "mass_eligible": 3,
        "charge_sources": {"reported": 3, "comb": 0, "dataset_value": 0, "none": 3}}
    assert s["charge"]["settings"] == {"charge_comb": False, "dataset_charge": None, "comb": None}
    assert s["charge"]["comb_status"] == {"not_run": 6}             # exported MS2 only
    assert ms2["issues"][R.ARRAYS_EMPTY] == [16]
    assert ms2["issues"][R.CHARGE_UNKNOWN] == [12, 13, 14, 16, 17]   # 17 also has no charge
    assert ms2["issues"][R.DECODE_ERROR] == [17]
    assert ms2["errors"][0]["scan"] == 17 and ms2["errors"][0]["stage"] == "decode"


def test_ms3_missing_metadata_does_not_stop_conversion(tmp_path):
    spectra = [MS1(0, 10),
               S.spectrum_xml(1, 11, 2, [150.0], [1.0], precursor={"parent": 10, "mz": 500.0, "charge": 2}),
               S.spectrum_xml(2, 12, 3, [90.0, 95.0], [1.0, 2.0], precursor={"mz": 150.0})]
    res = mzmlreader.extract_mzML(S.write_mzml(tmp_path / "ms3.mzML", spectra), outdir=str(tmp_path))
    rows = _rows(res["ms3tmp"])
    assert rows[0] == MS3_HEADER
    assert rows[1] == ['1', '0', '150.0', '150.0', '12', '(90.0, 95.0)', '(1.0, 2.0)']


def test_ms3_header_names_match_values(tmp_path):
    spectra = [MS1(0, 10),
               S.spectrum_xml(1, 11, 2, [150.0], [1.0], precursor={"parent": 10, "mz": 500.0, "charge": 2}),
               S.spectrum_xml(2, 12, 3, [90.0], [1.0], precursor={"parent": 11, "mz": 150.0})]
    res = mzmlreader.extract_mzML(S.write_mzml(tmp_path / "ms3h.mzML", spectra), outdir=str(tmp_path))
    row = dict(zip(*_rows(res["ms3tmp"])[:2]))
    assert row["MS2scan_no"] == "11" and row["MS3scan_no"] == "12"


@pytest.mark.parametrize("value, expected", [
    (None, None), (0, None), ("0", None), ("", None), ("nan", None), (float("nan"), None),
    (2.5, None), (True, None), (2, 2), ("2", 2), (2.0, 2), (-2, -2), (3, 3),
])
def test_normalize_charge(value, expected):
    assert R.normalize_charge(value) == expected


def _rebuild_excerpt(text):
    """Turn the hand-cut excerpt into a small valid non-indexed mzML (complete spectra only)."""
    head = re.search(r'<spectrumList [^>]*>', text)
    header = re.sub(r'<indexedmzML[^>]*>\s*', '', text[:head.start()])
    body_start = text.index('<spectrum index="6753"')
    body = text[body_start:text.rindex('</spectrum>') + len('</spectrum>')]
    counter = iter(range(1000))
    body = re.sub(r'<spectrum index="\d+"', lambda m: f'<spectrum index="{next(counter)}"', body)
    n = body.count('<spectrum index=')
    spectrum_list = re.sub(r'count="\d+"', f'count="{n}"', head.group(0))
    return header + spectrum_list + "\n" + body + "\n</spectrumList>\n</run>\n</mzML>\n"


@pytest.mark.skipif(not os.path.exists(PILOT_EXCERPT), reason="local pilot excerpt not present")
def test_pilot_excerpt_rebuilt(tmp_path):
    path = tmp_path / "pilot_LTQ_native_excerpt.mzML"
    with open(PILOT_EXCERPT, encoding="utf-8") as fh:
        path.write_text(_rebuild_excerpt(fh.read()), encoding="utf-8")
    recs = _records(str(path))
    assert {scan: len(r.mz) for scan, r in recs.items()} == {6781: 7, 6783: 10, 6785: 15}
    rec = recs[6785]
    assert rec.parent_scan == '6784' and rec.charge is None and rec.derived_mass is None
    assert rec.polarity == '-' and rec.representation == 'centroid'
    assert not rec.mass_eligible and rec.exportable
    res = mzmlreader.extract_mzML(str(path), outdir=str(tmp_path))
    assert res["summary"]["ms2"]["charge_unknown"] == 3 and res["summary"]["ms2"]["mass_eligible"] == 0


@pytest.mark.parametrize("mz, charge, polarity, expected", [
    (1038.8733, 2, '-', 1038.8733 * 2 + 1.007276),        # [M-2H]2- -> [M-H]-
    (1367.4973, 2, '-', 1367.4973 * 2 + 1.007276),
    (900.0, 3, '-', 900.0 * 3 + 2 * 1.007276),
    (1235.4418, 1, '-', 1235.4418),                       # z=1 negative: m/z is already [M-H]-
    (1038.8733, -2, None, 1038.8733 * 2 + 1.007276),      # negative charge value = negative mode
    (700.0, -1, None, 700.0),
    (700.0, 0, '-', 0),
    (500.25, 2, '+', 500.25 * 2 - 1.00784),               # positive mode unchanged
    (500.25, 2, None, 500.25 * 2 - 1.00784),              # unknown polarity: legacy formula
    (500.25, 1, None, 500.25),
])
def test_calcproton_polarity(mz, charge, polarity, expected):
    assert mzmlreader.calcproton(mz, charge, polarity=polarity) == pytest.approx(expected)


# READER-DIAG CP-C.5 (2026-10-15): one user charge for the whole file fills only unreported MS2 charges.
def test_dataset_charge_fills_only_unreported(mixed_file, tmp_path):
    res = mzmlreader.extract_mzML(mixed_file, outdir=str(tmp_path), dataset_charge=2)
    by_scan = {r[6]: r for r in _rows(res["ms2tmp"])[1:]}
    assert by_scan['12'][4:6] == ['2', str(512.5 * 2 + 1 * 1.007276)]     # no charge -> 2, negative mode
    assert by_scan['13'][4] == '2'                                         # reported 0 = not determined -> 2
    assert by_scan['14'][4:6] == ['2', '0']                                # no precursor m/z: still no mass
    assert by_scan['11'][4] == '2' and by_scan['18'][4] == '3'             # reported charges kept
    ms2 = res["summary"]["ms2"]
    assert ms2["charge_sources"] == {"reported": 3, "comb": 0, "dataset_value": 3, "none": 0}
    assert ms2["issues"][R.CHARGE_UNKNOWN] == []
    ch = res["summary"]["charge"]
    assert ch["settings"] == {"charge_comb": False, "dataset_charge": 2, "comb": None}
    assert ch["sources"] == {"reported": 3, "comb": 0, "dataset_value": 3, "none": 0}


def test_dataset_charge_off_is_byte_identical(mixed_file, tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    a = mzmlreader.extract_mzML(mixed_file, outdir=str(tmp_path / "a"))
    b = mzmlreader.extract_mzML(mixed_file, outdir=str(tmp_path / "b"), dataset_charge=None)
    assert open(a["ms2tmp"], "rb").read() == open(b["ms2tmp"], "rb").read()


@pytest.mark.parametrize("bad", [0, 2.5, False, "1"])
def test_invalid_dataset_charge_raises(mixed_file, tmp_path, bad):
    with pytest.raises(ValueError, match="dataset charge"):
        mzmlreader.extract_mzML(mixed_file, outdir=str(tmp_path), dataset_charge=bad)


@pytest.mark.parametrize("reported, comb, dataset, expected", [
    (2, 3, 1, (2, R.CHARGE_SOURCE_REPORTED)),
    (None, 3, 1, (3, R.CHARGE_SOURCE_COMB)),
    (None, None, 1, (1, R.CHARGE_SOURCE_DATASET)),
    (None, None, None, (None, R.CHARGE_SOURCE_NONE)),
    (-2, None, 1, (-2, R.CHARGE_SOURCE_REPORTED)),
])
def test_resolve_charge_order(reported, comb, dataset, expected):
    assert R.resolve_charge(reported, comb, dataset) == expected
