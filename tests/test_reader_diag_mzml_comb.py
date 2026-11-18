"""READER-DIAG CP-C.6: the MS1 isotope comb on the mzML route (2026-10-17).

Synthetic mzML only (tests/_synth_mzml.py); isotope patterns as in test_reader_diag_charge_comb.py.
Order as on the RAW route: reported charge (never replaced) -> comb -> dataset charge -> 0. Only the
latest MS1 is held, and the comb runs only on a profile MS1 named by the MS2's spectrumRef.
"""
import csv
import json
import math
import os
import sys

import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, '..', 'src'))
sys.path.insert(0, HERE)

pytest.importorskip("pymzml")

import mzmlreader  # noqa: E402
import _synth_mzml as S  # noqa: E402

PROTON = 1.007276


def _profile(mono, z, rel=(1.0, 0.85, 0.45, 0.18), sigma=0.08):
    mz = [mono - 4.2 + 0.031 + 0.0833 * i for i in range(101)]
    return mz, [sum(r * 1000 * math.exp(-(x - (mono + k * 1.00336 / z)) ** 2 / (2 * sigma ** 2))
                    for k, r in enumerate(rel)) for x in mz]


def _ms1(idx, scan, profile, **kw):
    return S.spectrum_xml(idx, scan, 1, profile[0], profile[1], **kw)


def _ms2(idx, scan, parent=None, mz=1039.1, charge=None):
    return S.spectrum_xml(idx, scan, 2, [300.1, 450.2], [5.0, 9.0],
                          precursor={"parent": parent, "mz": mz, "charge": charge})


Z2 = _profile(1038.87, 2)
LONE = _profile(799.9, 1, rel=(1.0,))
Z2_AND_LONE = (LONE[0] + Z2[0], LONE[1] + Z2[1])


@pytest.fixture
def comb_file(tmp_path):
    spectra = [
        _ms1(0, 10, Z2_AND_LONE, centroid=False),
        _ms2(1, 11, parent=10),                                   # comb z = 2
        _ms2(2, 12, parent=10, charge=3),                         # reported: kept, comb not run
        _ms2(3, 13, parent=10, mz=800.0),                         # lone peak: below margin
        _ms2(4, 14),                                              # no spectrumRef
        _ms1(5, 20, Z2),                                          # centroid MS1
        _ms2(6, 21, parent=20),
        _ms2(7, 22, parent=10),                                   # names an MS1 that is no longer held
        _ms1(8, 30, Z2, centroid=False, filter_string="FTMS - p ESI Full ms [380.00-2000.00]"),
        _ms2(9, 31, parent=30),
        _ms1(10, 40, Z2, centroid=False, filter_string="ITMS - p ESI Full ms [380.00-2000.00]"),
        _ms2(11, 41, parent=40),                                  # comb z = 2
        _ms2(12, 42, parent=40, mz=None),                         # no precursor m/z
        _ms1(13, 50, Z2, representation=False),                   # centroid/profile not stated
        _ms2(14, 51, parent=50),
    ]
    return S.write_mzml(tmp_path / "comb.mzML", spectra)


def _convert(path, out, **kw):
    res = mzmlreader.extract_mzML(path, outdir=str(out), **kw)
    with open(res["ms2tmp"], newline='') as fh:
        rows = list(csv.reader(fh, delimiter='\t'))[1:]
    return res, {r[6]: r for r in rows}


def test_comb_off_is_byte_identical(comb_file, tmp_path):
    a, _ = _convert(comb_file, tmp_path / "a")
    b, _ = _convert(comb_file, tmp_path / "b", charge_comb=False)
    assert open(a["ms2tmp"], "rb").read() == open(b["ms2tmp"], "rb").read()
    assert a["summary"]["charge"]["comb_status"] == {"not_run": 10}


def test_comb_fills_unreported_charge_and_never_overwrites(comb_file, tmp_path):
    res, rows = _convert(comb_file, tmp_path, charge_comb=True)
    assert rows['11'][4:6] == ['2', str(1039.1 * 2 + PROTON)]         # negative mode, selected-ion m/z
    assert rows['41'][4] == '2'
    assert rows['12'][4] == '3'
    for scan in ('13', '14', '21', '22', '31', '42', '51'):
        assert rows[scan][4] == '0'
    ch = res["summary"]["charge"]
    assert ch["settings"]["charge_comb"] is True and ch["settings"]["comb"]["margin"] == 0.1
    assert ch["sources"] == {"reported": 1, "comb": 2, "dataset_value": 0, "none": 7}
    assert res["summary"]["ms2"]["charge_sources"] == ch["sources"]
    assert ch["comb_status"] == {"called": 2, "not_run": 1, "below_margin": 1, "no_parent_ms1": 2,
                                 "no_profile_ms1": 2, "ft_ms1_not_tried": 1, "no_isolation_mz": 1}
    by_scan = {e["scan"]: e for e in ch["scans"]}
    assert sorted(by_scan) == [11, 13, 14, 21, 22, 31, 41, 42, 51]        # unreported spectra only
    assert by_scan[11]["comb_z"] == 2 and by_scan[11]["comb_margin"] >= 0.1
    assert by_scan[13]["comb_status"] == "below_margin" and by_scan[13]["comb_z"] is None
    assert by_scan[22]["comb_status"] == "no_parent_ms1"                  # only the latest MS1 is held


def test_dataset_charge_fills_what_the_comb_leaves(comb_file, tmp_path):
    res, rows = _convert(comb_file, tmp_path, charge_comb=True, dataset_charge=1)
    assert rows['11'][4] == '2' and rows['41'][4] == '2' and rows['12'][4] == '3'
    assert [rows[s][4] for s in ('13', '14', '21', '22', '31', '42', '51')] == ['1'] * 7
    ch = res["summary"]["charge"]
    assert ch["sources"] == {"reported": 1, "comb": 2, "dataset_value": 7, "none": 0}
    assert ch["comb_vs_dataset_disagree"] == 2


def test_cli_convert_mzml_with_comb(comb_file, tmp_path):
    import msp_extract_cli as CLI
    out = tmp_path / "cli"
    args = ["convert", "--mzml", comb_file, "--out", str(out), "--alias", "mzml_alias", "--charge-comb"]
    assert CLI.main(args) == CLI.EXIT_WARNINGS          # 7 MS2 charges stay unknown (CP-D)
    record = json.loads((out / "mzml_alias.convert.json").read_text())
    assert record["options"] == {"charge_comb": True, "dataset_charge": None}
    assert record["module"]["name"] == "mzmlreader" and record["module"]["version"] == str(mzmlreader.version)
    assert record["summary"]["charge"]["sources"]["comb"] == 2
