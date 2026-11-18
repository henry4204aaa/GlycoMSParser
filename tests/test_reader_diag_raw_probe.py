"""READER-DIAG CP-C: read-only RAW probe and `inspect` CLI, against a fake backend (2026-10-13).

The real Thermo backend is Windows-only and is exercised in Henry's Lab Windows step. These tests
pin the probe's own behaviour: independent requests, charge sources, parent window, alias redaction.
"""
import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import msp_raw_probe as P  # noqa: E402
import msp_extract_cli as CLI  # noqa: E402

SECRET = "SECRETNAME_20170101"


class BackendError(Exception):
    pass


class FakeReader:
    """Profile ion-trap file: MS1 100 -> MS2 101 (charge 0), MS2 102 (centroiding fails)."""
    ms1 = [(499.6, 5.0), (499.75, 10.0), (500.25, 100.0), (500.75, 40.0), (501.25, 12.0), (510.0, 7.0)]
    ms2 = [(150.0, 3.0), (150.1, 4.0), (355.1, 100.0)]

    def __init__(self, path):
        self.path = path
        self.closed = False

    def Version(self):
        return "fake 1.0"

    def GetNumSpectra(self):
        return 3

    def GetInstModel(self):
        return "FakeTrap"

    def GetMSOrderForScanNum(self, scan):
        return {100: 1, 101: 2, 102: 2}[scan]

    def GetFilterForScanNum(self, scan):
        return "ITMS - p ESI Full ms [380.00-2000.00]" if scan == 100 else \
            "ITMS - p ESI d w Full ms2 500.25@cid35.00 [135.00-1015.00]"

    def IsProfileScanForScanNum(self, scan):
        return True

    def IsCentroidScanForScanNum(self, scan):
        return False

    def GetPrecursorInfoFromScanNum(self, scan):
        if scan == 102:
            raise BackendError(f"COM failure reading {self.path}")
        return [500.25, 500.25, 0, 100]

    def GetTrailerExtraForScanNum(self, scan):
        return {"Charge State:": 0, "Monoisotopic M/Z:": 0.0, "Master Scan Number:": 100}

    def GetScanHeaderInfoForScanNum(self, scan):
        return {"StartTime": 1.5, "LowMass": 135.0}

    def GetLabelData(self, scan):
        return (((), (), (), (), (), ()), ())

    def GetMassListFromScanNum(self, scanNumber, scanFilter='', intensityCutoffType=0, intensityCutoffValue=0,
                               nMaxNumberOfPeaks=0, centroidResult=False, centroidPeakWidth=0):
        if centroidResult and scanNumber == 102:
            raise BackendError("Failed to centroid")
        points = self.ms1 if scanNumber == 100 else self.ms2
        if centroidResult:
            points = points[::2]
        return ((tuple(m for m, _ in points), tuple(i for _, i in points)), (0,) * len(points))

    def Close(self):
        self.closed = True


class PositionalReader(FakeReader):
    def GetMassListFromScanNum(self, *args):
        assert args[1:] in (("", 0, 0, 0, False, 0.0), ("", 0, 0, 0, True, 0.0))
        return FakeReader.GetMassListFromScanNum(self, args[0], centroidResult=args[5])


@pytest.fixture
def raw(tmp_path):
    path = tmp_path / f"{SECRET}.raw"
    path.write_bytes(b"not a real raw file")
    return str(path)


def _run(raw, tmp_path, scans=(101, 102), factory=FakeReader, **kw):
    out = tmp_path / "diag"
    report = P.inspect_raw_scans(raw, list(scans), "pilot_alias", str(out), reader_factory=factory, **kw)
    return report, out


def test_requests_are_independent(raw, tmp_path):
    report, _ = _run(raw, tmp_path)
    by_scan = {e["scan"]: e for e in report["scans"]}
    r101 = by_scan[101]["requests"]
    assert r101["label_data"]["ok"] and r101["label_data"]["n_mz"] == 0
    assert r101["mass_list_centroid_off"]["n_mz"] == 3
    assert r101["mass_list_centroid_on"]["n_mz"] == 2
    r102 = by_scan[102]["requests"]
    assert not r102["mass_list_centroid_on"]["ok"]
    assert r102["mass_list_centroid_on"]["exception"] == "BackendError"
    assert r102["mass_list_centroid_off"]["ok"] and r102["label_data"]["ok"]
    assert report["outcome"] == "complete"


def test_charge_sources_and_parent_window(raw, tmp_path):
    report, out = _run(raw, tmp_path)
    e = {x["scan"]: x for x in report["scans"]}[101]
    assert e["charge_sources"]["from_precursor_info"]["charge"] == 0
    assert e["charge_sources"]["from_trailer"]["charge_state"] == 0
    assert e["charge_sources"]["from_filter"]["precursor_mz"] == 500.25
    assert e["parent"]["scan"] == 100 and e["parent"]["found_by"] == "precursor_info"
    win = e["parent"]["requests"]["mass_list_centroid_off"]
    assert win["n_mz"] == 6 and win["n_in_window"] == 5          # 510.0 is outside +/-3
    rows = (out / win["window_points_file"]).read_text().splitlines()
    assert rows[0] == "mz\tintensity" and len(rows) == 6


def test_parent_falls_back_when_precursor_info_fails(raw, tmp_path):
    report, _ = _run(raw, tmp_path)
    e = {x["scan"]: x for x in report["scans"]}[102]
    assert e["charge_sources"]["from_precursor_info"] is None
    assert e["parent"]["scan"] == 100 and e["parent"]["found_by"] == "trailer_master_scan"


def test_alias_redaction_everywhere(raw, tmp_path):
    report, out = _run(raw, tmp_path)
    names = os.listdir(out)
    assert "pilot_alias.inspect.json" in names and "pilot_alias.inspect.txt" in names
    for name in names:
        assert SECRET not in name
        assert SECRET not in (out / name).read_text(), name
    msg = {x["scan"]: x for x in report["scans"]}[102]["info_calls"]
    failed = [c for c in msg if not c["ok"]][0]
    assert "pilot_alias" in failed["message"] and failed["traceback"]


def test_signatures_recorded_and_positional_fallback(raw, tmp_path):
    report, _ = _run(raw, tmp_path, scans=(101,), factory=PositionalReader)
    assert report["environment"]["signatures"]["GetMassListFromScanNum"] == "(*args)"
    r = report["scans"][0]["requests"]["mass_list_centroid_off"]
    assert r["ok"] and r["kwargs"] == {} and r["args"] == [101, "", 0, 0, 0, False, 0.0]
    report, _ = _run(raw, tmp_path, scans=(101,))
    assert report["scans"][0]["requests"]["mass_list_centroid_on"]["kwargs"]["centroidResult"] is True


def test_open_failure_still_writes_report(raw, tmp_path):
    def broken(path):
        raise BackendError(f"cannot open {path}")
    report, out = _run(raw, tmp_path, factory=broken)
    assert report["outcome"] == "failed" and report["open"]["exception"] == "BackendError"
    saved = json.loads((out / "pilot_alias.inspect.json").read_text())
    assert saved["outcome"] == "failed" and SECRET not in json.dumps(saved)


def test_cli_exit_codes(raw, tmp_path, monkeypatch):
    out = str(tmp_path / "cli")
    assert CLI.main(["inspect", raw, "--alias", "bad/alias", "--scans", "1", "--out", out]) == CLI.EXIT_USAGE
    assert CLI.main(["inspect", raw + ".missing", "--alias", "a", "--scans", "1", "--out", out]) == CLI.EXIT_USAGE
    real = P.inspect_raw_scans
    monkeypatch.setattr(P, "inspect_raw_scans",
                        lambda *a, **k: real(*a, reader_factory=FakeReader, **k))
    assert CLI.main(["inspect", raw, "--alias", "a", "--scans", "101", "--out", out]) == CLI.EXIT_OK
    monkeypatch.setattr(P, "inspect_raw_scans",
                        lambda *a, **k: real(*a, reader_factory=lambda p: (_ for _ in ()).throw(OSError("x")), **k))
    assert CLI.main(["inspect", raw, "--alias", "a", "--scans", "101", "--out", out]) == CLI.EXIT_FAILED
