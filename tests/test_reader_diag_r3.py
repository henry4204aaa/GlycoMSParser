"""READER-DIAG round 3 (Codex R3-F3, 2026-10-26): the precursor fields of an MS2 come from one selectedIon.

pymzml's selected_precursors collects the m/z, intensity and charge terms of all selected ions into separate lists
and pairs them by position, so an ion without a charge got the next ion's charge and a mass was derived from it.
The reader now reads the first selectedIon from the spectrum XML on the normal path too (R2-F9 fixed only the
fallback). Synthetic mzML (tests/_synth_mzml.py) only; no third-party data.
Spec: handoff/codex_READER_DIAG_r3_review_20261026_handoff.md (R3-F3); fix record in handoff/session_report.md.
"""
import os
import sys

import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, '..', 'src'))
sys.path.insert(0, HERE)

pymzml = pytest.importorskip("pymzml")

import mzmlreader  # noqa: E402
import _synth_mzml as S  # noqa: E402

MS1 = ([800.0, 800.5, 801.0], [10.0, 6.0, 3.0])
PEAKS = ([204.0867, 366.1395, 528.1923], [1500.5, 320.25, 77.125])


def _file(tmp_path, name, precursor):
    path = tmp_path / name
    path.write_text(S.build_mzml([S.spectrum_xml(0, 1, 1, *MS1, polarity="+"),
                                  S.spectrum_xml(1, 2, 2, *PEAKS, polarity="+", precursor=precursor)]),
                    encoding="utf-8")
    return str(path)


# Ion A: m/z 500, no charge, no intensity. Ion B: m/z 600, charge 3, intensity 42.
TWO_IONS_A_UNCHARGED = {"parent": 1, "mz": 500.0, "more_ions": [{"mz": 600.0, "charge": 3, "intensity": 42.0}]}


def _ms2_spectrum(path):
    for spectrum in pymzml.run.Reader(path):
        if spectrum.ms_level == 2:
            return spectrum
    pytest.fail("no MS2 spectrum")


def _csv_rows(path):
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    cols = lines[0].split("\t")
    return [dict(zip(cols, line.split("\t"))) for line in lines[1:]]


def test_root_cause_pymzml_pairs_selected_ion_fields_by_position(tmp_path):
    """The backend behaviour R3-F3 works around (pymzml 2.6): ion A's m/z is paired with ion B's charge."""
    precursor = _ms2_spectrum(_file(tmp_path, "a.mzML", TWO_IONS_A_UNCHARGED)).selected_precursors[0]
    assert precursor["mz"] == 500.0 and precursor["charge"] == 3 and precursor["i"] == 42.0


def test_normal_path_never_borrows_charge_or_intensity_from_another_ion(tmp_path):
    rec = mzmlreader.read_mzml_spectrum(_ms2_spectrum(_file(tmp_path, "a.mzML", TWO_IONS_A_UNCHARGED)), 2)
    assert rec.precursor_mz == 500.0
    assert rec.charge is None and rec.precursor_intensity is None
    assert rec.derived_mass is None and not rec.mass_eligible
    assert rec.parent_scan == "1"


def test_extract_writes_unknown_charge_and_mass_as_zero(tmp_path):
    res = mzmlreader.extract_mzML(_file(tmp_path, "a.mzML", TWO_IONS_A_UNCHARGED), outdir=str(tmp_path / "o"))
    (row,) = _csv_rows(res["ms2tmp"])
    assert (row["MS1_isolationmass"], row["chargeState"], row["protonatedmass"]) == ("500.0", "0", "0")


def test_later_ions_do_not_change_a_complete_first_ion(tmp_path):
    """CSV bytes for a complete first ion are the same with or without further ions (one-ion files unchanged)."""
    one = mzmlreader.extract_mzML(_file(tmp_path, "one.mzML", {"parent": 1, "mz": 500.0, "charge": 2}),
                                  outdir=str(tmp_path / "o1"), stem="x")
    two = mzmlreader.extract_mzML(_file(tmp_path, "two.mzML", {"parent": 1, "mz": 500.0, "charge": 2,
                                                               "more_ions": [{"mz": 600.0, "charge": 3}]}),
                                  outdir=str(tmp_path / "o2"), stem="x")
    with open(one["ms2tmp"], "rb") as a, open(two["ms2tmp"], "rb") as b:
        assert a.read() == b.read()
    assert _csv_rows(one["ms2tmp"])[0]["chargeState"] == "2"


def _with_second_precursor(spectrum_xml, mz, charge):
    """Append a second <precursor> (one selectedIon) to a synthetic spectrum's precursorList."""
    second = (f'<precursor spectrumRef="controllerType=0 controllerNumber=1 scan=1"><selectedIonList count="1">'
              f'<selectedIon>{S._cv("MS:1000744", "selected ion m/z", mz)}{S._cv("MS:1000041", "charge state", charge)}'
              f'</selectedIon></selectedIonList></precursor></precursorList>')
    return spectrum_xml.replace('<precursorList count="1">', '<precursorList count="2">').replace(
        '</precursorList>', second)


def test_no_charge_borrowed_from_a_second_precursor(tmp_path):
    """R3C-F1: two precursors with one ion each; pymzml paired the first ion's m/z with the second's charge."""
    path = tmp_path / "p.mzML"
    ms2 = _with_second_precursor(S.spectrum_xml(1, 2, 2, *PEAKS, polarity="+", precursor={"parent": 1, "mz": 500.0}),
                                 600.0, 3)
    path.write_text(S.build_mzml([S.spectrum_xml(0, 1, 1, *MS1, polarity="+"), ms2]), encoding="utf-8")
    assert _ms2_spectrum(str(path)).selected_precursors[0]["charge"] == 3   # the backend behaviour
    res = mzmlreader.extract_mzML(str(path), outdir=str(tmp_path / "o"))
    (row,) = _csv_rows(res["ms2tmp"])
    assert (row["MS1scan_no"], row["chargeState"], row["protonatedmass"]) == ("1", "0", "0")


def test_empty_scan_number_in_spectrum_ref_gives_parent_0(tmp_path):
    """R3C-F1 (documented CSV change): pymzml read `scan=` as an empty parent (blank cell); now unknown, written 0."""
    res = mzmlreader.extract_mzML(_file(tmp_path, "e.mzML", {"parent": "", "mz": 500.0, "charge": 2}),
                                  outdir=str(tmp_path / "o"))
    (row,) = _csv_rows(res["ms2tmp"])
    assert (row["MS1scan_no"], row["chargeState"]) == ("0", "2")


class _NoElementSpectrum:
    """A spectrum without an XML element: pymzml's first precursor dict is the only source left."""
    ID = 2
    element = None
    selected_precursors = [{"mz": 650.0, "charge": 2, "precursor id": "1"}]
    mz = (204.0867,)
    i = (10.0,)


def test_without_xml_element_pymzml_dict_is_used():
    rec = mzmlreader.read_mzml_spectrum(_NoElementSpectrum(), 2)
    assert rec.precursor_mz == 650.0 and rec.charge == 2 and rec.parent_scan == "1"
