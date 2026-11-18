"""READER-DIAG round 2 (Codex R2-F1 ... R2-F14, Henry's decisions D1-D5, 2026-10-23).

Synthetic mzML (tests/_synth_mzml.py) and fake MSFileReader backends only; no third-party data.
Spec: handoff/READER_DIAG_R2_fix_spec_20261023.md.
"""
import json
import math
import os
import re
import sys
import xml.etree.ElementTree as ET

import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, '..', 'src'))
sys.path.insert(0, HERE)

pytest.importorskip("pymzml")

import mzmlreader  # noqa: E402
import msprawextractor as X  # noqa: E402
import msp_extract_cli as CLI  # noqa: E402
import msp_raw_probe  # noqa: E402
import msp_charge_comb  # noqa: E402
import msp_spectrum_record as R  # noqa: E402
import _synth_mzml as S  # noqa: E402

MS1 = ([800.0, 800.5, 801.0], [10.0, 6.0, 3.0])
PEAKS_A = ([204.0867, 366.1395, 528.1923], [1500.5, 320.25, 77.125])
PEAKS_B = ([138.0545, 274.0921], [44.0, 9.5])
MS3_HEADER = "entry no\tMS2scan_no\tMS2_isolationmass\tMS2_monoisomass\tMS3scan_no\tpeaklist\tpeakintensity\n"


def _write(path, spectra):
    path.write_text(S.build_mzml(spectra), encoding="utf-8")
    return str(path)


def _ms1(idx=0, scan=1, **kw):
    return S.spectrum_xml(idx, scan, 1, *MS1, **kw)


def _ms2(idx, scan, charge=2, mz=800.0, peaks=PEAKS_B, **kw):
    return S.spectrum_xml(idx, scan, 2, *peaks, precursor={"parent": 1, "mz": mz, "charge": charge}, **kw)


def _ms3(idx, scan):
    return S.spectrum_xml(idx, scan, 3, *PEAKS_B, precursor={"parent": 2, "mz": 366.1395})


def _diag(outdir, stem):
    with open(os.path.join(outdir, f"{stem}.reader_diag.json"), encoding="utf-8") as fh:
        return json.load(fh)


def _status(outdir, stem):
    with open(os.path.join(outdir, f"{stem}.scan_status.tsv"), encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    cols = lines[0].split("\t")
    return [dict(zip(cols, line.split("\t"))) for line in lines[1:]]


def _all_text(folder):
    text = ""
    for root, _, files in os.walk(folder):
        for f in files:
            text += f + "\n" + open(os.path.join(root, f), encoding="utf-8", errors="replace").read()
    return text


def _work_folders(out):
    return [d for d in os.listdir(out) if d.startswith("_work_")]


# ---- R2-F14: golden CSV bytes from the base 72935af readers ------------------------------------------------------
# Generated once (2026-10-23) by running `git show 72935af:src/mzmlreader.py` / `msprawextractor.py` on the inputs
# below. Domain: positive mode, every charge and m/z reported, no empty spectra (spec section 4). The MS3 header
# differs by the approved CP-C.2 name swap, so only its rows are compared.

GOLDEN_MZML_MS2 = (
    'entry_no\tMS1scan_no\tMS1_isolationmass\tMS1_monoisolationmass\tchargeState\tprotonatedmass\tMS2scan_no\tpeaklist\t'
    'peakintensity\n'
    '1\t1\t800.3412\t800.3412\t1\t800.3412\t2\t(204.08670043945312, 366.1394958496094, 528.1923217773438)\t'
    '(1500.5, 320.25, 77.125)\n'
    '2\t1\t650.7781\t650.7781\t2\t1300.54836\t3\t(138.05450439453125, 274.09210205078125)\t(44.0, 9.5)\n'
    '3\t1\t512.2501\t512.2501\t3\t1534.73462\t5\t(204.08670043945312, 366.1394958496094, 528.1923217773438)\t'
    '(1500.5, 320.25, 77.125)\n')
GOLDEN_MZML_MS3_ROWS = '1\t3\t366.1395\t366.1395\t4\t(138.05450439453125, 274.09210205078125)\t(44.0, 9.5)\n'
GOLDEN_RAW_MS2 = (
    'entry_no\tMS1scan_no\tMS1_isolationmass\tMS1_monoisolationmass\tchargeState\tprotonatedmass\tMS2scan_no\tpeaklist\t'
    'peakintensity\n'
    '1\t1\t800.3412\t800.3398\t1\t800.3398\t2\t[204.0867, 366.1395, 528.1923]\t[1500.5, 320.25, 77.125]\n'
    '2\t1\t650.7781\t650.7764\t2\t1300.54496\t3\t[138.0545, 274.0921]\t[44.0, 9.5]\n'
    '3\t1\t512.2501\t512.2489\t3\t1534.7310200000002\t4\t[204.0867, 366.1395, 528.1923]\t[1500.5, 320.25, 77.125]\n')


def _golden_mzml_spectra():
    return [S.spectrum_xml(0, 1, 1, *MS1, polarity="+"),
            S.spectrum_xml(1, 2, 2, *PEAKS_A, polarity="+", precursor={"parent": 1, "mz": 800.3412, "charge": 1}),
            S.spectrum_xml(2, 3, 2, *PEAKS_B, polarity="+", precursor={"parent": 1, "mz": 650.7781, "charge": 2}),
            S.spectrum_xml(3, 4, 3, *PEAKS_B, polarity="+", precursor={"parent": 3, "mz": 366.1395}),
            S.spectrum_xml(4, 5, 2, *PEAKS_A, polarity="+", precursor={"parent": 1, "mz": 512.2501, "charge": 3})]


class GoldenFTRaw:
    """Positive-mode FTMS scans with complete metadata (the base reader's label-data path)."""
    scans = {
        1: (1, "FTMS + p NSI Full ms [400.00-2000.00]", None, MS1),
        2: (2, "FTMS + c NSI d Full ms2 800.34@hcd30.00 [110.00-1610.00]", [800.3412, 800.3398, 1, 1], PEAKS_A),
        3: (2, "FTMS + c NSI d Full ms2 650.78@hcd30.00 [110.00-1310.00]", [650.7781, 650.7764, 2, 1], PEAKS_B),
        4: (2, "FTMS + c NSI d Full ms2 512.25@hcd30.00 [110.00-1550.00]", [512.2501, 512.2489, 3, 1], PEAKS_A),
    }
    fail_at = None

    def __init__(self, path):
        pass

    def GetNumSpectra(self):
        return len(self.scans)

    def GetMSOrderForScanNum(self, scan):
        if scan == self.fail_at:
            raise RuntimeError("backend stopped")
        return self.scans[scan][0]

    def GetFilterForScanNum(self, scan):
        return self.scans[scan][1]

    def GetPrecursorInfoFromScanNum(self, scan):
        return list(self.scans[scan][2])

    def GetLabelData(self, scan):
        return ((self.scans[scan][3][0], self.scans[scan][3][1], (), ()), ())

    def IsProfileScanForScanNum(self, scan):
        return 0

    def GetInstModel(self):
        return "Orbitrap"

    def Close(self):
        pass


@pytest.fixture
def ftraw(tmp_path, monkeypatch):
    monkeypatch.setattr(X, "MSFileReader", GoldenFTRaw, raising=False)
    monkeypatch.setattr(GoldenFTRaw, "fail_at", None)
    raw = tmp_path / "g.raw"
    raw.write_bytes(b"")
    return str(raw)


def _csv_bytes(path):
    with open(path, "rb") as fh:
        return fh.read()


def _expected_bytes(text):
    """R3-F14 (2026-10-26): compare bytes, not newline-translated text. The readers write their CSVs in text mode,
    so every line ends in os.linesep (LF on macOS/Linux, CRLF on Windows); the goldens above are written with LF."""
    return text.replace("\n", os.linesep).encode("utf-8")


def test_golden_mzml_csv_equals_base_reader(tmp_path):
    path = _write(tmp_path / "g.mzML", _golden_mzml_spectra())
    res = mzmlreader.extract_mzML(path, outdir=str(tmp_path / "o"))
    assert _csv_bytes(res["ms2tmp"]) == _expected_bytes(GOLDEN_MZML_MS2)
    assert _csv_bytes(res["ms3tmp"]) == _expected_bytes(MS3_HEADER + GOLDEN_MZML_MS3_ROWS)


def test_golden_raw_csv_equals_base_reader(tmp_path, ftraw):
    res = X.convert_raw_to_csv(ftraw, outdir=str(tmp_path / "o"))
    assert _csv_bytes(res["ms2tmp"]) == _expected_bytes(GOLDEN_RAW_MS2)
    assert _csv_bytes(res["ms3tmp"]) == _expected_bytes(MS3_HEADER)
    assert res["exit_code"] == R.EXIT_OK


# ---- R2-F1: the CLI never deletes what it did not create (Henry: CLI only, v14-style promotion) -----------------

def test_cli_keeps_existing_files_and_reused_work_folder_names(tmp_path):
    out = tmp_path / "out"
    old_work = out / "_work_x"
    old_work.mkdir(parents=True)
    (old_work / "notes.txt").write_text("keep")
    inp = _write(old_work / "input.mzML", [_ms1(), _ms2(1, 2)])     # the input inside the old folder name
    (out / "other.txt").write_text("keep")
    assert CLI.main(["convert", "--mzml", inp, "--out", str(out), "--alias", "x"]) == R.EXIT_OK
    assert (old_work / "notes.txt").read_text() == "keep" and os.path.isfile(inp)
    assert (out / "other.txt").read_text() == "keep"
    assert _work_folders(str(out)) == ["_work_x"]                     # this run's own folder is gone
    record = json.loads((out / "x.convert.json").read_text())
    assert "work_folder_left" not in record and record["ms2"]["file"] == "ms2_x.csv"


# ---- R2-F2: a CSV written before a failure is kept (Henry D1); header-only MS3 stays a normal result -------------

def _fail_second_write(monkeypatch, exc=None):
    real, calls = mzmlreader.save_to_csv, []

    def flaky(stem, data, debug=False):
        calls.append(stem)
        if len(calls) == 2:
            raise exc or OSError(28, "No space left on device")
        return real(stem, data, debug=debug)
    monkeypatch.setattr(mzmlreader, "save_to_csv", flaky)


def test_cli_keeps_ms2_when_ms3_write_fails(tmp_path, monkeypatch):
    inp = _write(tmp_path / "in.mzML", [_ms1(), _ms2(1, 2), _ms3(2, 3)])
    out = tmp_path / "out"
    _fail_second_write(monkeypatch)
    assert CLI.main(["convert", "--mzml", inp, "--out", str(out), "--alias", "x"]) == R.EXIT_FAILED
    record = json.loads((out / "x.convert.json").read_text())
    assert record["outcome"] == "failed" and record["exception"] == "OSError"
    assert record["ms2"]["file"] == "ms2_x.csv" and record["ms2"]["rows"] == 1 and record["ms3"] is None
    assert (out / "ms2_x.csv").is_file() and not (out / "ms3_x.csv").exists()
    diag = _diag(str(out), "x")
    assert diag["outcome"] == "failed" and diag["failure"]["stage"] == "write"
    assert diag["outputs"]["ms2_csv"] == "ms2_x.csv" and diag["outputs"].get("ms3_csv") is None
    rows = {r["scan"]: r for r in _status(str(out), "x")}
    assert rows["2"]["exported"] == "yes" and rows["3"]["exported"] == "no"           # Henry D3
    assert "not_written" in rows["3"]["issues"].split(",")
    assert diag["counts"]["ms3"]["not_written"] == 1 and diag["counts"]["ms2"]["exported"] == 1
    assert _work_folders(str(out)) == []


def test_file_without_ms3_still_writes_header_only_ms3_and_exits_0(tmp_path):
    inp = _write(tmp_path / "in.mzML", [_ms1(), _ms2(1, 2)])
    out = tmp_path / "out"
    assert CLI.main(["convert", "--mzml", inp, "--out", str(out), "--alias", "x"]) == R.EXIT_OK
    assert (out / "ms3_x.csv").read_text() == MS3_HEADER


def test_cli_finalize_failure_is_failed_in_both_records(tmp_path, monkeypatch):
    inp = _write(tmp_path / "in.mzML", [_ms1(), _ms2(1, 2)])
    out = tmp_path / "out"
    real = CLI._describe

    def describe(path):
        if os.path.basename(path).startswith("ms3_"):
            raise PermissionError("locked")
        return real(path)
    monkeypatch.setattr(CLI, "_describe", describe)
    assert CLI.main(["convert", "--mzml", inp, "--out", str(out), "--alias", "x"]) == R.EXIT_FAILED
    record = json.loads((out / "x.convert.json").read_text())
    assert record["finalize_error"]["stage"] == "finalize" and record["outcome"] == "failed"
    diag = _diag(str(out), "x")
    assert diag["outcome"] == "failed" and diag["exit_code"] == R.EXIT_FAILED
    assert diag["failure"]["stage"] == "finalize"
    assert "Outcome: failed" in (out / "x.reader_diag.txt").read_text()


def test_cli_csv_that_cannot_be_moved_stays_in_the_work_folder(tmp_path, monkeypatch):
    inp = _write(tmp_path / "in.mzML", [_ms1(), _ms2(1, 2)])
    out = tmp_path / "out"
    real = CLI._promote

    def promote(src, dst):
        if os.path.basename(src).startswith("ms3tmp_"):
            raise PermissionError("locked")
        return real(src, dst)
    monkeypatch.setattr(CLI, "_promote", promote)
    assert CLI.main(["convert", "--mzml", inp, "--out", str(out), "--alias", "x"]) == R.EXIT_FAILED
    record = json.loads((out / "x.convert.json").read_text())
    left = record["work_folder_left"]
    assert left["entries"] == ["ms3tmp_x.csv"] and (out / left["folder"] / "ms3tmp_x.csv").is_file()
    assert (out / "ms2_x.csv").is_file() and (out / "x.reader_diag.json").is_file()


# ---- R2-F6: the diagnostics never replace the conversion's exception; interruptions are recorded ----------------

def test_diagnostics_build_error_keeps_original_exception(tmp_path, monkeypatch):
    path = _write(tmp_path / "in.mzML", [_ms1(), _ms2(1, 2)])
    original = RuntimeError("conversion broke")

    def boom(*a, **kw):
        raise original
    monkeypatch.setattr(mzmlreader, "_extract_mzML", boom)
    monkeypatch.setattr(R, "_count_rows", lambda rows: 1 / 0)
    with pytest.raises(RuntimeError) as info:
        mzmlreader.extract_mzML(path, outdir=str(tmp_path / "o"))
    assert info.value is original


def test_diagnostics_error_on_success_gives_cli_exit_5(tmp_path, monkeypatch):
    inp = _write(tmp_path / "in.mzML", [_ms1(), _ms2(1, 2)])
    monkeypatch.setattr(R, "write_diagnostics", lambda *a, **k: 1 / 0)
    res = mzmlreader.extract_mzML(inp, outdir=str(tmp_path / "o"))
    assert res["exit_code"] == R.EXIT_OK and res["diagnostics"] == {}
    assert res["diagnostics_error"].startswith("ZeroDivisionError")
    out = tmp_path / "out"
    assert CLI.main(["convert", "--mzml", inp, "--out", str(out), "--alias", "x"]) == R.EXIT_FAILED   # Henry D2
    record = json.loads((out / "x.convert.json").read_text())
    assert record["diagnostics_error"] and record["ms2"]["file"] == "ms2_x.csv"


def test_keyboard_interrupt_finishes_diagnostics_and_is_reraised(tmp_path, monkeypatch):
    inp = _write(tmp_path / "in.mzML", [_ms1(), _ms2(1, 2), _ms2(2, 3)])
    stop = KeyboardInterrupt()
    real, n = mzmlreader.read_mzml_spectrum, []

    def reader(*a, **kw):
        n.append(1)
        if len(n) == 2:
            raise stop
        return real(*a, **kw)
    monkeypatch.setattr(mzmlreader, "read_mzml_spectrum", reader)
    out = tmp_path / "out"
    with pytest.raises(KeyboardInterrupt) as info:
        CLI.main(["convert", "--mzml", inp, "--out", str(out), "--alias", "x"])
    assert info.value is stop
    diag = _diag(str(out), "x")
    assert diag["failure"]["exception"] == "KeyboardInterrupt" and diag["failure"]["last_scan_started"] == 3
    assert json.loads((out / "x.convert.json").read_text())["outcome"] == "failed"
    assert _work_folders(str(out)) == []


# ---- R2-F7: stage and counts after a failure ---------------------------------------------------------------------

def test_failure_at_first_scan_is_not_an_open_failure(tmp_path, ftraw, monkeypatch):
    monkeypatch.setattr(GoldenFTRaw, "fail_at", 1)
    with pytest.raises(RuntimeError):
        X.convert_raw_to_csv(ftraw, outdir=str(tmp_path / "o"))
    diag = _diag(str(tmp_path / "o"), "g")
    assert diag["failure"]["stage"] == "iterate" and "could not be opened" not in diag["explanation"]


def test_mid_run_failure_marks_read_rows_not_written(tmp_path, monkeypatch):
    path = _write(tmp_path / "in.mzML", [_ms1(), _ms2(1, 2), _ms2(2, 3)])
    real, n = mzmlreader.read_mzml_spectrum, []

    def reader(*a, **kw):
        n.append(1)
        if len(n) == 2:
            raise RuntimeError("boom")
        return real(*a, **kw)
    monkeypatch.setattr(mzmlreader, "read_mzml_spectrum", reader)
    with pytest.raises(RuntimeError):
        mzmlreader.extract_mzML(path, outdir=str(tmp_path / "o"))
    diag = _diag(str(tmp_path / "o"), "in")
    assert diag["counts"]["ms2"]["exported"] == 0 and diag["counts"]["ms2"]["not_written"] == 1
    assert diag["failure"]["last_scan_started"] == 3 and diag["failure"]["last_ms_level"] == 2
    assert any("not_written" in w for w in diag["warnings"])


# ---- R2-F3/F4/F5: privacy -----------------------------------------------------------------------------------------

def test_no_name_or_folder_in_any_output_or_stdout(tmp_path, monkeypatch, capsys):
    folder = tmp_path / "privfolder"
    folder.mkdir()
    inp = folder / "secretname.mzML"
    text = S.build_mzml([_ms1(), _ms2(1, 2)]).replace(
        'id="controllerType=0 controllerNumber=1 scan=2"', 'id="file=secretname.mzML scan=2"')
    inp.write_text(text, encoding="utf-8")
    real = mzmlreader.read_mzml_spectrum

    def reader(*a, **kw):
        rec = real(*a, **kw)
        rec.errors.append(("precursor", "ValueError", f"bad value in {inp}"))
        return rec
    monkeypatch.setattr(mzmlreader, "read_mzml_spectrum", reader)
    out = folder / "out"                                  # --out inside the input's folder
    CLI.main(["convert", "--mzml", str(inp), "--out", str(out), "--alias", "x"])
    produced = _all_text(str(out)) + capsys.readouterr().out
    assert "secretname" not in produced and "privfolder" not in produced
    assert "file=x scan=2" in (out / "x.scan_status.tsv").read_text()


def test_redactor_path_forms_case_short_folder_and_root():
    red = R.Redactor(r"C:\Data\Secret.raw", "x")
    for message in (repr(r"C:\Data\Secret.raw"), "c:\\data\\secret.RAW", "C:/Data/Secret.raw",
                    json.dumps({"p": r"C:\Data\Secret.raw"}), "folder C:\\Data missing"):
        assert "Data" not in red(message) and "data" not in red(message).lower().replace("x", ""), message
        assert "Secret" not in red(message), message
    assert R.Redactor("/ab/secret.raw", "x")("cannot list /ab") == "cannot list x"
    assert R.Redactor("/secret.raw", "x")("/ is the root") == "/ is the root"


def test_redactor_stem_whole_word_and_idempotent():
    red = R.Redactor("/tmp/q/a.raw", "public")
    assert red("Charge State: 2 for a.raw and a") == "Charge State: 2 for public and public"
    red = R.Redactor("/tmp/q/bad.raw", "bad_raw")
    once = red("RAWfile /tmp/q/bad.raw failed; bad.raw; bad")
    assert once == "RAWfile bad_raw failed; bad_raw; bad_raw" and red(once) == once


def test_probe_keeps_trailer_labels_with_a_short_stem():
    red = msp_raw_probe.Redactor("/tmp/q/a.raw", "public")
    trailer = msp_raw_probe._jsonable({"Charge State:": 2, "Master Scan Number:": 7, "Note": "/tmp/q/a.raw"}, red)
    assert msp_raw_probe._trailer_value(trailer, "charge state") == 2
    assert msp_raw_probe._trailer_value(trailer, "master scan number") == 7
    assert trailer["Note"] == "public"


# ---- R2-F8/F9/F10/F11 and D5 --------------------------------------------------------------------------------------

def _comb_status(tmp_path, ms1_id):
    profile = S.spectrum_xml(0, 1, 1, [799.0, 800.0, 800.5, 801.0, 801.5], [1.0, 9.0, 4.0, 2.0, 1.0],
                             centroid=False)
    text = S.build_mzml([profile, _ms2(1, 2, charge=None)])
    text = text.replace('id="controllerType=0 controllerNumber=1 scan=1"', f'id="{ms1_id}"')
    path = tmp_path / f"c{len(ms1_id)}.mzML"
    path.write_text(text, encoding="utf-8")
    mzmlreader.extract_mzML(str(path), outdir=str(tmp_path), charge_comb=True, stem="c")
    return _status(str(tmp_path), "c")[0]["comb_status"]


def test_comb_parent_matched_by_full_native_id(tmp_path):
    assert _comb_status(tmp_path, "controllerType=0 controllerNumber=2 scan=1") == "no_parent_ms1"
    assert _comb_status(tmp_path, "controllerType=0 controllerNumber=1 scan=1") != "no_parent_ms1"


def test_xml_fallback_reads_one_selected_ion():
    el = ET.fromstring(
        '<spectrum xmlns="http://psi.hupo.org/ms/mzml"><precursorList><precursor spectrumRef="scan=1">'
        '<selectedIonList><selectedIon><cvParam accession="MS:1000744" value="500"/></selectedIon>'
        '<selectedIon><cvParam accession="MS:1000744" value="600"/><cvParam accession="MS:1000041" value="3"/>'
        '</selectedIon></selectedIonList></precursor></precursorList></spectrum>')
    found = mzmlreader._precursor_from_element(el)
    assert found["mz"] == 500.0 and "charge" not in found


def test_non_finite_inputs():
    assert R.normalize_charge(float("inf")) is None and R.normalize_charge(float("-inf")) is None
    for mz, inten in (([800.0, 800.5], [float("nan"), 1.0]), ([800.0, float("inf")], [1.0, 1.0]),
                      ([800.0, 800.5, 801.0], [1.0, 2.0])):
        res = msp_charge_comb.infer_charge(mz, inten, 800.0)
        assert res.status == msp_charge_comb.STATUS_INVALID_PROFILE and res.z is None and res.margin is None


def test_mass_eligibility_is_one_rule_and_the_explanation_names_the_mz():
    rec = R.SpectrumRecord(ms_level=2, scan=2, precursor_mz=0.0, charge=1, derived_mass=0.0, mz=(1.0,),
                           intensity=(1.0,))
    assert not rec.mass_eligible
    row = R._complete_row(R.status_row(rec))
    assert row["mass_status"] == "unavailable"
    counts, rep, pol, _ = R._count_rows([row])
    _, code, warnings, _, explanation = R.classify_outcome(counts, rep, None, pol)
    assert code == R.EXIT_WARNINGS and explanation == R.NO_MASS_ELIGIBLE_MZ_EXPLANATION


def test_raw_rows_record_mass_mz_source_for_reported_charges(tmp_path, ftraw):
    X.convert_raw_to_csv(ftraw, outdir=str(tmp_path / "o"))
    rows = [r for r in _status(str(tmp_path / "o"), "g") if r["ms_level"] == "2"]
    assert [r["mass_mz_source"] for r in rows] == ["monoisotopic"] * 3
    assert all(r["charge_source"] == "reported" for r in rows)


def _without_polarity(spectrum_xml):
    return re.sub(r'<cvParam[^>]*accession="MS:10001(?:29|30)"[^>]*/>', "", spectrum_xml)


def test_polarity_from_the_filter_string_when_the_term_is_missing(tmp_path):
    """Henry D5 (2026-10-23): no polarity term -> the filter string's token; the negative formula then applies."""
    ms2 = _without_polarity(_ms2(1, 2, charge=2, filter_string="ITMS - c ESI d w Full ms2 800.00@cid35.00"))
    ms2b = _without_polarity(_ms2(2, 3, charge=2))
    path = _write(tmp_path / "p.mzML", [_ms1(), ms2, ms2b])
    res = mzmlreader.extract_mzML(path, outdir=str(tmp_path / "o"))
    lines = open(res["ms2tmp"], encoding="utf-8").read().splitlines()
    assert math.isclose(float(lines[1].split("\t")[5]), 800.0 * 2 + 1.007276)          # [M-H]- equivalent
    assert math.isclose(float(lines[2].split("\t")[5]), 800.0 * 2 - 1.00784)            # unknown: legacy positive
    rows = {r["scan"]: r for r in _status(str(tmp_path / "o"), "p")}
    assert (rows["2"]["polarity"], rows["2"]["polarity_source"]) == ("-", "filter")
    assert (rows["3"]["polarity"], rows["3"]["polarity_source"]) == ("", "")
    assert any("state no polarity" in n for n in _diag(str(tmp_path / "o"), "p")["notes"])
