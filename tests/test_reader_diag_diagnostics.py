"""READER-DIAG CP-D: per-scan status file, reader diagnostics files, outcome and exit codes (2026-10-20).

Synthetic mzML (tests/_synth_mzml.py) and a fake MSFileReader only. Schema: READER-DIAG handoff B3, B6, B7.
Henry 2026-10-20: both readers write the files (GUI included); exit codes as B6; spectra skipped without peaks
are a note, not a warning.
"""
import csv
import json
import os
import sys

import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, '..', 'src'))
sys.path.insert(0, HERE)

pytest.importorskip("pymzml")

import mzmlreader  # noqa: E402
import msprawextractor as X  # noqa: E402
import msp_spectrum_record as R  # noqa: E402
import _synth_mzml as S  # noqa: E402

PEAKS = ([300.1, 450.2], [5.0, 9.0])


def _ms1(idx, scan):
    return S.spectrum_xml(idx, scan, 1, [800.0, 801.0], [10.0, 5.0])


def _ms2(idx, scan, charge=None, peaks=PEAKS, **kw):
    return S.spectrum_xml(idx, scan, 2, peaks[0], peaks[1],
                          precursor={"parent": 1, "mz": 800.0, "charge": charge}, **kw)


def _write(path, spectra, header=""):
    text = S.build_mzml(spectra)
    if header:
        text = text.replace("<run id=", header + "<run id=")
    path.write_text(text, encoding="utf-8")
    return str(path)


def _files(outdir, stem):
    with open(os.path.join(outdir, f"{stem}.reader_diag.json"), encoding="utf-8") as fh:
        diag = json.load(fh)
    with open(os.path.join(outdir, f"{stem}.scan_status.tsv"), newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    text = open(os.path.join(outdir, f"{stem}.reader_diag.txt"), encoding="utf-8").read()
    return diag, rows, text


# ---- mzML route -------------------------------------------------------------------------------------------------

def test_gui_style_call_complete_exit_0(tmp_path):
    """No stem given (as v14 calls it): files named after the input stem, next to the CSVs."""
    path = _write(tmp_path / "run_a.mzML", [_ms1(0, 1), _ms2(1, 2, charge=2), _ms2(2, 3, charge=1)])
    res = mzmlreader.extract_mzML(path, outdir=str(tmp_path / "out"))
    assert (res["outcome"], res["exit_code"], res["warnings"]) == ("complete", 0, [])
    assert os.path.basename(res["ms2tmp"]) == "ms2tmp_run_a.csv"
    diag, rows, text = _files(str(tmp_path / "out"), "run_a")
    assert diag["json_type"] == "glycomsp.reader_diagnostics" and diag["schema_version"] == "1.1.0"   # R2: not_written, failure scan
    assert diag["tool"]["reader"] == "mzmlreader" and diag["tool"]["entry"] == "direct"
    assert diag["counts"]["spectra_by_ms_level"] == {"1": 1, "2": 2}
    assert diag["counts"]["ms2"]["mass_eligible"] == 2 and diag["counts"]["ms2"]["charge_sources"]["reported"] == 2
    assert diag["coverage"] == {"spectra_declared": 3, "spectra_processed": 3, "last_scan_processed": 3,
                                "aborted": False}
    assert diag["outputs"]["ms2_csv"] == "ms2tmp_run_a.csv" and diag["failure"] is None
    assert [r["scan"] for r in rows] == ["2", "3"]
    assert rows[0]["charge"] == "2" and rows[0]["charge_status"] == "known" and rows[0]["mass_status"] == "eligible"
    assert rows[0]["charge_source"] == "reported" and rows[0]["comb_status"] == "not_run"
    assert rows[0]["representation"] == "centroid" and rows[0]["polarity"] == "-" and rows[0]["n_peaks"] == "2"
    assert "exit code 0" in text and "3 complete with warnings" in text


def test_csv_unchanged_by_stem(tmp_path):
    path = _write(tmp_path / "run_b.mzML", [_ms1(0, 1), _ms2(1, 2, charge=2), _ms2(2, 3)])
    a = mzmlreader.extract_mzML(path, outdir=str(tmp_path / "a"))
    b = mzmlreader.extract_mzML(path, outdir=str(tmp_path / "b"), stem="alias_b")
    assert os.path.basename(b["ms2tmp"]) == "ms2tmp_alias_b.csv"
    for key in ("ms2tmp", "ms3tmp"):
        assert open(a[key], "rb").read() == open(b[key], "rb").read()


def test_unknown_charge_is_a_warning_and_empty_spectra_a_note(tmp_path):
    path = _write(tmp_path / "run_c.mzML", [_ms1(0, 1), _ms2(1, 2, charge=2), _ms2(2, 3),
                                            _ms2(3, 4, charge=2, peaks=([], []))])
    res = mzmlreader.extract_mzML(path, outdir=str(tmp_path), stem="c")
    assert (res["outcome"], res["exit_code"]) == ("complete", R.EXIT_WARNINGS)
    diag, rows, text = _files(str(tmp_path), "c")
    assert len(diag["warnings"]) == 1 and diag["warnings"][0].startswith("1 of 2 exported MS2 spectra have no charge")
    assert diag["notes"] == ["1 MS2 and 0 MS3 spectra without peaks were not exported (as MSConvert does)."]
    assert diag["explanation"] is None
    by_scan = {r["scan"]: r for r in rows}
    assert by_scan["3"]["charge"] == "" and by_scan["3"]["charge_status"] == "unknown"
    assert by_scan["3"]["mass_status"] == "unavailable" and by_scan["3"]["charge_source"] == "none"
    assert "charge_unknown" in by_scan["3"]["issues"].split(",")
    assert by_scan["4"]["exported"] == "no" and by_scan["4"]["arrays_status"] == "empty"
    assert by_scan["4"]["charge_status"] == "" and by_scan["4"]["n_peaks"] == "0"
    assert diag["affected_scans"]["ms2"]["charge_unknown"] == [3]
    assert "complete with warnings - exit code 3" in text


def test_no_mass_eligible_gives_the_explanation(tmp_path):
    path = _write(tmp_path / "run_d.mzML", [_ms1(0, 1), _ms2(1, 2), _ms2(2, 3)])
    res = mzmlreader.extract_mzML(path, outdir=str(tmp_path), stem="d")
    assert res["exit_code"] == R.EXIT_WARNINGS
    diag, _, _ = _files(str(tmp_path), "d")
    assert diag["explanation"] == R.NO_MASS_ELIGIBLE_EXPLANATION
    assert R.NO_MASS_ELIGIBLE_EXPLANATION in diag["warnings"]


def test_comb_and_dataset_sources_are_notes(tmp_path):
    path = _write(tmp_path / "run_e.mzML", [_ms1(0, 1), _ms2(1, 2), _ms2(2, 3, charge=1)])
    res = mzmlreader.extract_mzML(path, outdir=str(tmp_path), stem="e", dataset_charge=2)
    assert res["exit_code"] == R.EXIT_OK
    diag, rows, _ = _files(str(tmp_path), "e")
    assert diag["notes"] == ["1 MS2 charges came from the dataset charge (status file: charge_source dataset_value)."]
    assert diag["settings"]["dataset_charge"] == 2 and diag["charge"]["sources"]["dataset_value"] == 1
    assert rows[0]["charge_source"] == "dataset_value" and rows[0]["charge"] == "2"


def test_profile_ms2_is_a_warning(tmp_path):
    path = _write(tmp_path / "run_f.mzML", [_ms1(0, 1), _ms2(1, 2, charge=1, centroid=False)])
    res = mzmlreader.extract_mzML(path, outdir=str(tmp_path), stem="f")
    assert res["exit_code"] == R.EXIT_WARNINGS
    diag, rows, _ = _files(str(tmp_path), "f")
    assert diag["representation"]["ms2"] == {"centroid": 0, "profile": 1, "unknown": 0}
    assert "profile data" in diag["warnings"][0] and "peakPicking vendor msLevel=2-" in diag["warnings"][0]
    assert rows[0]["representation"] == "profile"


def test_undecodable_spectrum_is_partial(tmp_path):
    path = _write(tmp_path / "run_g.mzML", [_ms1(0, 1), _ms2(1, 2, charge=1),
                                            _ms2(2, 3, charge=1, mz_binary="!!not-base64!!")])
    res = mzmlreader.extract_mzML(path, outdir=str(tmp_path), stem="g")
    assert (res["outcome"], res["exit_code"]) == ("partial", R.EXIT_PARTIAL)
    diag, rows, _ = _files(str(tmp_path), "g")
    assert diag["counts"]["ms2"]["failed"] == 1 and diag["outputs"]["partial"] is True
    assert {r["scan"]: r["arrays_status"] for r in rows} == {"2": "ok", "3": "decode_error"}


def test_no_ms2_with_peaks_is_failed(tmp_path):
    path = _write(tmp_path / "run_h.mzML", [_ms1(0, 1), _ms2(1, 2, charge=1, peaks=([], []))])
    res = mzmlreader.extract_mzML(path, outdir=str(tmp_path), stem="h")
    assert (res["outcome"], res["exit_code"]) == ("failed", R.EXIT_FAILED)
    diag, _, _ = _files(str(tmp_path), "h")
    assert diag["explanation"] == "No MS2 spectrum with peaks was exported."


def test_open_failure_returns_false_and_records_it(tmp_path):
    bad = tmp_path / "SECRET_NAME_bad.mzML"
    bad.write_text("not xml", encoding="utf-8")
    assert mzmlreader.extract_mzML(str(bad), outdir=str(tmp_path / "o"), stem="alias_x") is False
    diag, rows, _ = _files(str(tmp_path / "o"), "alias_x")
    assert diag["outcome"] == "failed" and diag["exit_code"] == R.EXIT_FAILED and rows == []
    assert diag["failure"]["stage"] == "open" and diag["coverage"]["aborted"] is False
    assert all(":" in frame and "/" not in frame and "\\" not in frame for frame in diag["failure"]["traceback"])


def test_exception_while_reading_is_recorded_and_reraised(tmp_path, monkeypatch):
    path = _write(tmp_path / "SECRET_NAME_run.mzML", [_ms1(0, 1), _ms2(1, 2, charge=1), _ms2(2, 3, charge=1)])
    real = mzmlreader.read_mzml_spectrum

    def boom(spectrum, ms_level, **kw):
        if str(spectrum.ID) == "3":
            raise RuntimeError(f"cannot read {path}")
        return real(spectrum, ms_level, **kw)
    monkeypatch.setattr(mzmlreader, "read_mzml_spectrum", boom)
    with pytest.raises(RuntimeError):
        mzmlreader.extract_mzML(path, outdir=str(tmp_path / "o"), stem="alias_y")
    diag, rows, text = _files(str(tmp_path / "o"), "alias_y")
    assert diag["outcome"] == "failed" and diag["failure"]["stage"] == "iterate"
    assert diag["coverage"]["aborted"] is True and diag["coverage"]["last_scan_processed"] == 3
    assert len(rows) == 1
    for name in os.listdir(tmp_path / "o"):
        assert "SECRET_NAME" not in open(tmp_path / "o" / name, encoding="utf-8").read()
    assert "SECRET_NAME" not in diag["failure"]["message"] and "alias_y" in diag["failure"]["message"]


HEADER = (
    '<fileDescription><sourceFileList count="1"><sourceFile id="SECRET_SRC" name="SECRET_NAME.raw" '
    'location="file:///C:/SECRET_DIR"><cvParam cvRef="MS" accession="MS:1000569" name="SHA-1" value="abc"/>'
    '</sourceFile></sourceFileList></fileDescription>'
    '<referenceableParamGroupList count="1"><referenceableParamGroup id="CommonInstrumentParams">'
    '<cvParam cvRef="MS" accession="MS:1000447" name="LTQ" value=""/>'
    '<cvParam cvRef="MS" accession="MS:1000529" name="instrument serial number" value="SECRET_SERIAL"/>'
    '</referenceableParamGroup></referenceableParamGroupList>'
    '<softwareList count="1"><software id="pwiz" version="3.0.1">'
    '<cvParam cvRef="MS" accession="MS:1000615" name="ProteoWizard software" value=""/></software></softwareList>'
    '<instrumentConfigurationList count="1"><instrumentConfiguration id="IC1">'
    '<cvParam cvRef="MS" accession="MS:1000133" name="electrospray ionization" value=""/>'
    '</instrumentConfiguration></instrumentConfigurationList>'
    '<dataProcessingList count="1"><dataProcessing id="pwiz_Reader_conversion"><processingMethod order="0" '
    'softwareRef="pwiz"><cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML" value=""/>'
    '</processingMethod></dataProcessing></dataProcessingList>\n')


def test_instrument_terms_only_names_from_the_header(tmp_path):
    path = _write(tmp_path / "run_i.mzML", [_ms1(0, 1), _ms2(1, 2, charge=1)], header=HEADER)
    mzmlreader.extract_mzML(path, outdir=str(tmp_path), stem="i")
    diag, _, _ = _files(str(tmp_path), "i")
    assert diag["instrument"] == {"instrument": ["LTQ", "instrument serial number", "electrospray ionization"],
                                  "software": [{"name": "ProteoWizard software", "version": "3.0.1"}],
                                  "processing": ["Conversion to mzML"]}
    text = open(tmp_path / "i.reader_diag.json", encoding="utf-8").read()
    for secret in ("SECRET", "synthetic_run", "SHA-1"):
        assert secret not in text


# ---- RAW route (fake MSFileReader) ------------------------------------------------------------------------------

class FakeRaw:
    """Scan -> (ms order, filter, precursor info [iso, mono, charge, parent], peaks)."""
    scans = {
        1: (1, "ITMS - c ESI Full ms [380.00-2000.00]", None, PEAKS),
        2: (2, "ITMS - c ESI d w Full ms2 800.00@cid35.00 [215.00-1610.00]", [800.0, 0.0, 2, 1], PEAKS),
        3: (2, "ITMS - c ESI d w Full ms2 900.00@cid35.00 [245.00-1810.00]", [900.0, 0.0, 0, 1], PEAKS),
        4: (2, "ITMS - c ESI d w Full ms2 950.00@cid35.00 [260.00-1910.00]", [950.0, 0.0, 1, 1], ([], [])),
        5: (3, "ITMS - c ESI d w Full ms3 800.00@cid35.00 400.00@cid35.00 [110.00-810.00]", [400.0, 0.0, 0, 2],
            PEAKS),
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

    def IsProfileScanForScanNum(self, scan):
        return 0

    def GetPrecursorInfoFromScanNum(self, scan):
        return list(self.scans[scan][2])

    def GetMassListFromScanNum(self, scanNumber, scanFilter='', intensityCutoffType=0, intensityCutoffValue=0,
                               maxNumberOfPeaks=0, centroidResult=False, centroidPeakWidth=0.0):
        return (self.scans[scanNumber][3], ())

    def GetInstModel(self):
        return "LTQ"

    def Close(self):
        pass


def _raw(tmp_path, monkeypatch, name="run.raw"):
    monkeypatch.setattr(X, "MSFileReader", FakeRaw, raising=False)
    monkeypatch.setattr(FakeRaw, "fail_at", None)
    raw = tmp_path / name
    raw.write_bytes(b"")
    return str(raw)


def test_raw_status_rows_and_warnings(tmp_path, monkeypatch):
    raw = _raw(tmp_path, monkeypatch)
    res = X.convert_raw_to_csv(raw, outdir=str(tmp_path))
    assert (res["outcome"], res["exit_code"]) == ("complete", R.EXIT_WARNINGS)
    diag, rows, _ = _files(str(tmp_path), "run")
    assert diag["tool"]["reader"] == "msprawextractor" and diag["input"]["format"] == "raw"
    assert diag["instrument"] == {"model": "LTQ", "software_version": None}
    assert diag["counts"]["spectra_by_ms_level"] == {"1": 1, "2": 3, "3": 1}
    assert diag["counts"]["ms2"] == {"seen": 3, "exported": 2, "empty": 1, "failed": 0, "not_written": 0, "charge_known": 1,
                                     "charge_unknown": 1, "mass_eligible": 1,
                                     "charge_sources": {"reported": 1, "comb": 0, "dataset_value": 0, "none": 1}}
    by_scan = {r["scan"]: r for r in rows}
    assert by_scan["2"]["charge_status"] == "known" and by_scan["2"]["peak_source"] == "mass_list"
    assert by_scan["2"]["parent_scan"] == "1" and by_scan["2"]["precursor_mz"] == "800.0"
    assert by_scan["3"]["charge"] == "" and by_scan["3"]["mass_status"] == "unavailable"
    assert by_scan["3"]["issues"] == "charge_unknown,derived_mass_unavailable"
    assert by_scan["4"]["exported"] == "no" and by_scan["4"]["arrays_status"] == "empty"
    assert by_scan["5"]["ms_level"] == "3" and by_scan["5"]["exported"] == "yes" and by_scan["5"]["charge_status"] == ""
    assert diag["notes"] == ["1 MS2 and 0 MS3 spectra without peaks were not exported (as MSConvert does)."]


def test_raw_open_failure(tmp_path, monkeypatch):
    raw = _raw(tmp_path, monkeypatch)

    def broken(path):
        raise OSError(f"cannot open {path}")
    monkeypatch.setattr(X, "MSFileReader", broken)
    assert X.convert_raw_to_csv(raw, outdir=str(tmp_path / "o"), stem="alias_r") is False
    diag, _, _ = _files(str(tmp_path / "o"), "alias_r")
    assert diag["failure"]["stage"] == "open" and diag["exit_code"] == R.EXIT_FAILED
    assert str(tmp_path) not in json.dumps(diag)


def test_raw_exception_mid_run_is_recorded_and_reraised(tmp_path, monkeypatch):
    raw = _raw(tmp_path, monkeypatch)
    monkeypatch.setattr(FakeRaw, "fail_at", 3)
    with pytest.raises(RuntimeError, match="backend stopped"):
        X.convert_raw_to_csv(raw, outdir=str(tmp_path))
    diag, rows, _ = _files(str(tmp_path), "run")
    assert diag["failure"]["stage"] == "iterate" and diag["coverage"]["last_scan_processed"] == 2
    assert diag["coverage"]["aborted"] is True and diag["outcome"] == "failed" and len(rows) == 1


# ---- CLI --------------------------------------------------------------------------------------------------------

def test_cli_records_outcome_and_moves_the_files(tmp_path, monkeypatch):
    import msp_extract_cli as CLI
    path = _write(tmp_path / "SECRET_NAME_run.mzML", [_ms1(0, 1), _ms2(1, 2, charge=1), _ms2(2, 3)])
    out = tmp_path / "cli"
    assert CLI.main(["convert", "--mzml", path, "--out", str(out), "--alias", "al"]) == CLI.EXIT_WARNINGS
    assert sorted(p.name for p in out.iterdir()) == ["al.convert.json", "al.reader_diag.json", "al.reader_diag.txt",
                                                     "al.scan_status.tsv", "ms2_al.csv", "ms3_al.csv"]
    record = json.loads((out / "al.convert.json").read_text())
    assert record["ok"] is True and record["outcome"] == "complete" and record["exit_code"] == 3
    assert record["diagnostics"] == {"json": "al.reader_diag.json", "txt": "al.reader_diag.txt",
                                     "status_tsv": "al.scan_status.tsv"}
    diag, _, text = _files(str(out), "al")
    assert diag["tool"]["entry"] == "cli-convert"
    assert diag["outputs"]["ms2_csv"] == "ms2_al.csv" and "ms2_al.csv" in text
    assert not any("SECRET_NAME" in p.read_text() for p in out.iterdir())


def test_cli_failed_and_partial_exit_codes(tmp_path):
    import msp_extract_cli as CLI
    bad = tmp_path / "bad.mzML"
    bad.write_text("not xml", encoding="utf-8")
    out = tmp_path / "cli"
    assert CLI.main(["convert", "--mzml", str(bad), "--out", str(out), "--alias", "b"]) == CLI.EXIT_FAILED
    record = json.loads((out / "b.convert.json").read_text())
    assert record["ok"] is False and record["outcome"] == "failed" and record["exit_code"] == 5
    assert (out / "b.reader_diag.txt").is_file()
    part = _write(tmp_path / "p.mzML", [_ms1(0, 1), _ms2(1, 2, charge=1),
                                        _ms2(2, 3, charge=1, mz_binary="!!not-base64!!")])
    assert CLI.main(["convert", "--mzml", part, "--out", str(out), "--alias", "p"]) == CLI.EXIT_PARTIAL


# ---- CP-D.1: one-pass redaction (Lab PC 2026-10-20: "bad.raw" with alias "bad_raw" gave "bad_raw_raw") ----------

def test_redaction_never_rescans_the_inserted_alias(tmp_path, monkeypatch):
    import msp_extract_cli as CLI
    import msp_raw_probe
    path = str(tmp_path / "bad.raw")
    message = f"RAWfile {path} could not be opened; bad.raw; bad"
    assert R.Redactor(path, "bad_raw")(message) == "RAWfile bad_raw could not be opened; bad_raw; bad_raw"
    assert msp_raw_probe.Redactor(path, "bad_raw")(message) == "RAWfile bad_raw could not be opened; bad_raw; bad_raw"
    open(path, "wb").close()

    def broken(p):
        raise OSError(f"RAWfile {p} could not be opened")
    monkeypatch.setattr(X, "MSFileReader", broken, raising=False)
    out = tmp_path / "cli"
    assert CLI.main(["convert", "--raw", path, "--out", str(out), "--alias", "bad_raw"]) == CLI.EXIT_FAILED
    diag, _, text = _files(str(out), "bad_raw")
    assert diag["failure"]["message"] == "RAWfile bad_raw could not be opened"
    assert "bad_raw_raw" not in text and str(tmp_path) not in text


def test_cli_exception_message_redacted_once(tmp_path, monkeypatch):
    import msp_extract_cli as CLI
    path = _write(tmp_path / "run.mzML", [_ms1(0, 1), _ms2(1, 2, charge=1)])

    def boom(*a, **kw):
        raise RuntimeError(f"cannot read {path} (run)")
    monkeypatch.setattr(mzmlreader, "_extract_mzML", boom)
    out = tmp_path / "cli"
    assert CLI.main(["convert", "--mzml", path, "--out", str(out), "--alias", "run_x"]) == CLI.EXIT_FAILED
    record = json.loads((out / "run_x.convert.json").read_text())
    assert record["message"] == "cannot read run_x (run_x)"
