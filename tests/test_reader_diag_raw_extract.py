"""READER-DIAG: RAW conversion peak source per scan, against a fake MSFileReader (2026-10-13).

The real Thermo backend is Windows-only. These tests pin which request `convert_raw_to_csv` makes for
each scan type and that MS3 rows carry their own peaks.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import msprawextractor as X  # noqa: E402

FT_LABELS = ((300.1, 450.2, 600.3), (10.0, 20.0, 30.0))
MASS_LIST_CENTROID = ((149.01, 355.07), (4.5, 233.1))
MASS_LIST_STORED = ((201.5,), (7.0,))
MS3_CENTROID = ((180.0, 220.0), (3.0, 9.0))


class FakeRaw:
    """Scan table: number -> (ms order, filter, profile?)."""
    scans = {
        1: (1, "FTMS - p ESI Full ms [380.00-2000.00]", True),
        2: (2, "FTMS - p ESI d Full ms2 800.50@cid35.00 [215.00-1610.00]", True),
        3: (1, "ITMS - p ESI Full ms [380.00-2000.00]", True),
        4: (2, "ITMS - p ESI d w Full ms2 385.34@cid35.00 [95.00-785.00]", True),
        5: (2, "ITMS - c ESI d Full ms2 400.00@cid35.00 [100.00-815.00]", False),
        6: (3, "ITMS - p ESI d Full ms3 385.34@cid35.00 355.07@cid35.00 [95.00-370.00]", True),
        7: (2, "ITMS - p ESI d w Full ms2 1417.55@cid35.00 [380.00-2000.00]", True),
    }

    def __init__(self, path):
        self.calls = []
        FakeRaw.last = self

    def GetNumSpectra(self):
        return len(self.scans)

    def GetMSOrderForScanNum(self, scan):
        return self.scans[scan][0]

    def GetFilterForScanNum(self, scan):
        return self.scans[scan][1]

    def IsProfileScanForScanNum(self, scan):
        return 1 if self.scans[scan][2] else 0

    def GetPrecursorInfoFromScanNum(self, scan):
        return [385.34, 385.34, 0.0, scan - 1]

    def GetLabelData(self, scan):
        self.calls.append(("label", scan))
        if self.scans[scan][1].startswith("FTMS"):
            return (FT_LABELS + ((), (), (), ()), ())
        return (((), (), (), (), (), ()), ())

    def GetMassListFromScanNum(self, scanNumber, scanFilter='', intensityCutoffType=0, intensityCutoffValue=0,
                               maxNumberOfPeaks=0, centroidResult=False, centroidPeakWidth=0.0):
        self.calls.append(("masslist", scanNumber, scanFilter, intensityCutoffType, intensityCutoffValue,
                           maxNumberOfPeaks, centroidResult, centroidPeakWidth))
        if scanNumber == 7:
            return (((), ()), ())
        if scanNumber == 6:
            return (MS3_CENTROID, (0, 0))
        return ((MASS_LIST_CENTROID if centroidResult else MASS_LIST_STORED), (0,) * 2)

    def Close(self):
        pass


@pytest.fixture
def converted(tmp_path, monkeypatch):
    monkeypatch.setattr(X, "MSFileReader", FakeRaw, raising=False)
    raw = tmp_path / "synthetic.raw"
    raw.write_bytes(b"")
    res = X.convert_raw_to_csv(str(raw), outdir=str(tmp_path))
    ms2 = [line.split("\t") for line in open(res["ms2tmp"]).read().splitlines()[1:]]
    ms3_lines = [line.split("\t") for line in open(res["ms3tmp"]).read().splitlines()]
    ms3 = [dict(zip(ms3_lines[0], r)) for r in ms3_lines[1:]]
    return res, {int(r[6]): r for r in ms2}, {int(r["MS3scan_no"]): r for r in ms3}, FakeRaw.last


def test_ftms_scan_keeps_label_data(converted):
    res, ms2, _, raw = converted
    assert ms2[2][7] == str(FT_LABELS[0]) and ms2[2][8] == str(FT_LABELS[1])
    assert ("label", 2) in raw.calls
    assert not [c for c in raw.calls if c[0] == "masslist" and c[1] == 2]


def test_ion_trap_profile_uses_vendor_centroided_mass_list(converted):
    res, ms2, _, raw = converted
    assert ms2[4][7] == str(MASS_LIST_CENTROID[0]) and ms2[4][8] == str(MASS_LIST_CENTROID[1])
    assert ("masslist", 4, "", 0, 0, 0, True, 0.0) in raw.calls
    assert ("label", 4) not in raw.calls


def test_ion_trap_centroid_scan_returns_stored_list(converted):
    res, ms2, _, raw = converted
    assert ms2[5][7] == str(MASS_LIST_STORED[0])
    assert ("masslist", 5, "", 0, 0, 0, False, 0.0) in raw.calls


def test_ms3_row_has_its_own_peaks(converted):
    res, ms2, ms3, raw = converted
    assert ms3[6]["peaklist"] == str(MS3_CENTROID[0]) and ms3[6]["peakintensity"] == str(MS3_CENTROID[1])
    assert ms3[6]["peaklist"] != ms2[5][7]


def test_ms3_header_names_match_values(converted):
    res, _, ms3, _ = converted
    assert list(ms3[6]) == ['entry no', 'MS2scan_no', 'MS2_isolationmass', 'MS2_monoisomass', 'MS3scan_no',
                            'peaklist', 'peakintensity']
    assert ms3[6]["MS2scan_no"] == "5" and ms3[6]["MS3scan_no"] == "6"   # fake parent = scan - 1


def test_scan_without_peaks_is_skipped_and_reported(converted):
    res, ms2, _, _ = converted
    assert 7 not in ms2 and sorted(ms2) == [2, 4, 5]
    assert sorted(int(r[0]) for r in ms2.values()) == [1, 2, 3]       # entry_no stays contiguous
    peaks = res["summary"]["peaks"]
    assert peaks["ms2"]["by_source"] == {X.PEAK_SOURCE_LABEL_DATA: 1, X.PEAK_SOURCE_MASS_LIST_CENTROID: 2,
                                         X.PEAK_SOURCE_MASS_LIST: 1}
    assert peaks["ms2"]["exported"] == 3 and peaks["ms2"]["empty"] == [7]
    assert peaks["ms3"] == {"by_source": {X.PEAK_SOURCE_MASS_LIST_CENTROID: 1}, "exported": 1, "empty": []}


def test_ftms_detection_reads_the_analyzer_token():
    class R:
        def __init__(self, flt):
            self.flt = flt

        def GetFilterForScanNum(self, scan):
            return self.flt

        def GetLabelData(self, scan):
            return (FT_LABELS, ())

        def IsProfileScanForScanNum(self, scan):
            return 0

        def GetMassListFromScanNum(self, scan, **kw):
            return (MASS_LIST_STORED, ())

    assert X.read_scan_peaks(R("FTMS {1,2} + p NSI d Full ms2 800.50@hcd30.00 [100-1610]"), 1)[2] == \
        X.PEAK_SOURCE_LABEL_DATA
    assert X.read_scan_peaks(R("ITMS + c NSI d Full ms2 800.50@cid35.00 [215-1610]"), 1)[2] == \
        X.PEAK_SOURCE_MASS_LIST
    assert X.read_scan_peaks(R(""), 1)[2] == X.PEAK_SOURCE_MASS_LIST


@pytest.mark.parametrize("flt, expected", [
    ("ITMS - p ESI d w Full ms2 1235.29@cid35.00 [330.00-2000.00]", "-"),
    ("FTMS + p NSI Full ms [380.00-2000.00]", "+"),
    ("FTMS {1,1} - c NSI d Full ms2 800.50@hcd30.00 [100.00-1610.00]", "-"),
    ("", None),
    ("ms2 something without polarity", None),
])
def test_scan_polarity(flt, expected):
    assert X.scan_polarity(flt) == expected


def test_negative_mode_mass_uses_polarity(tmp_path, monkeypatch):
    class NegCharge2(FakeRaw):
        scans = {1: (1, "ITMS - p ESI Full ms [380.00-2000.00]", True),
                 2: (2, "ITMS - p ESI d w Full ms2 1039.10@cid35.00 [275.00-2000.00]", True),
                 3: (1, "FTMS + p NSI Full ms [380.00-2000.00]", True),
                 4: (2, "FTMS + p NSI d Full ms2 1039.10@hcd30.00 [100.00-2000.00]", True)}

        def GetPrecursorInfoFromScanNum(self, scan):
            return [1039.10, 1038.8733, 2, scan - 1]

    monkeypatch.setattr(X, "MSFileReader", NegCharge2, raising=False)
    raw = tmp_path / "neg.raw"
    raw.write_bytes(b"")
    res = X.convert_raw_to_csv(str(raw), outdir=str(tmp_path))
    rows = {int(r.split("\t")[6]): r.split("\t") for r in open(res["ms2tmp"]).read().splitlines()[1:]}
    assert float(rows[2][5]) == pytest.approx(1038.8733 * 2 + 1.007276)   # [M-H]- equivalent
    assert float(rows[4][5]) == pytest.approx(1038.8733 * 2 - 1.00784)    # positive scan: legacy formula
