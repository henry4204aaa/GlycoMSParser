"""READER-DIAG CP-C.5: charge order in the RAW converter, against a fake MSFileReader (2026-10-15).

Order: valid reported charge (never overwritten) -> MS1 isotope comb (charge_comb=True) -> dataset charge -> 0.
Synthetic profiles only (see test_reader_diag_charge_comb.py).
"""
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import msprawextractor as X  # noqa: E402

PROTON = 1.007276
MS2_PEAKS = ((300.1, 450.2), (5.0, 9.0))


def _profile(mono, z, rel=(1.0, 0.85, 0.45, 0.18), sigma=0.08):
    mz = [mono - 4.2 + 0.031 + 0.0833 * i for i in range(101)]
    return mz, [sum(r * 1000 * math.exp(-(x - (mono + k * 1.00336 / z)) ** 2 / (2 * sigma ** 2))
                    for k, r in enumerate(rel)) for x in mz]


Z2_PROFILE = _profile(1038.87, 2)
LONE_PROFILE = _profile(799.9, 1, rel=(1.0,))


class ChargeRaw:
    """Scan -> (ms order, filter, profile?, precursor info [iso, mono, charge, parent], MS1 profile)."""
    scans = {
        1: (1, "ITMS - p ESI Full ms [380.00-2000.00]", True, None, Z2_PROFILE),
        2: (2, "ITMS - p ESI d w Full ms2 1039.10@cid35.00 [275.00-2000.00]", True, [1039.10, 0.0, 0, 1], None),
        3: (2, "ITMS - p ESI d w Full ms2 1039.10@cid35.00 [275.00-2000.00]", True, [1039.10, 1038.87, 3, 1], None),
        4: (2, "ITMS - p ESI d w Full ms2 1039.20@cid35.00 [275.00-2000.00]", True, [1039.20, 0.0, 0, 1], None),
        5: (1, "ITMS - c ESI Full ms [380.00-2000.00]", False, None, None),
        6: (2, "ITMS - p ESI d w Full ms2 800.00@cid35.00 [215.00-1610.00]", True, [800.0, 0.0, 0.0, 5], None),
        7: (1, "FTMS - p ESI Full ms [380.00-2000.00]", True, None, Z2_PROFILE),
        8: (2, "ITMS - p ESI d w Full ms2 1039.10@cid35.00 [275.00-2000.00]", True, [1039.10, 0.0, 0, 7], None),
        9: (2, "ITMS - p ESI d w Full ms2 1039.10@cid35.00 [275.00-2000.00]", True, [1039.10, 0.0, 0, 0], None),
        10: (1, "ITMS - p ESI Full ms [380.00-2000.00]", True, None, LONE_PROFILE),
        11: (2, "ITMS - p ESI d w Full ms2 800.00@cid35.00 [215.00-1610.00]", True, [800.0, 0.0, 0, 10], None),
        12: (2, "ITMS - p ESI d w Full ms2 1235.29@cid35.00 [330.00-2000.00]", True, [1235.29, 0.0, 1, 10], None),
    }

    def __init__(self, path):
        self.calls = []
        ChargeRaw.last = self

    def GetNumSpectra(self):
        return len(self.scans)

    def GetMSOrderForScanNum(self, scan):
        return self.scans[scan][0]

    def GetFilterForScanNum(self, scan):
        return self.scans[scan][1]

    def IsProfileScanForScanNum(self, scan):
        return 1 if self.scans[scan][2] else 0

    def GetPrecursorInfoFromScanNum(self, scan):
        return list(self.scans[scan][3])

    def GetMassListFromScanNum(self, scanNumber, scanFilter='', intensityCutoffType=0, intensityCutoffValue=0,
                               maxNumberOfPeaks=0, centroidResult=False, centroidPeakWidth=0.0):
        self.calls.append(("masslist", scanNumber, centroidResult))
        if self.scans[scanNumber][0] == 1:
            return (self.scans[scanNumber][4], ())
        return (MS2_PEAKS, ())

    def Close(self):
        pass


def _convert(tmp_path, monkeypatch, **kw):
    monkeypatch.setattr(X, "MSFileReader", ChargeRaw, raising=False)
    raw = tmp_path / "synthetic.raw"
    raw.write_bytes(b"")
    res = X.convert_raw_to_csv(str(raw), outdir=str(tmp_path), **kw)
    lines = open(res["ms2tmp"]).read().splitlines()
    return res, {int(r.split("\t")[6]): r.split("\t") for r in lines[1:]}, ChargeRaw.last


def _ms1_reads(raw):
    return [c for c in raw.calls if ChargeRaw.scans[c[1]][0] == 1]


def test_defaults_off_keep_legacy_cells(tmp_path, monkeypatch):
    res, rows, raw = _convert(tmp_path, monkeypatch)
    assert rows[2][4:6] == ["0", "0"] and rows[6][4:6] == ["0.0", "0"]      # charge written as the file gave it
    assert rows[3][3:6] == ["1038.87", "3", str(1038.87 * 3 + 2 * PROTON)]  # reported charge, mono m/z
    assert _ms1_reads(raw) == []                                            # no MS1 read when the comb is off
    ch = res["summary"]["charge"]
    assert ch["settings"] == {"charge_comb": False, "dataset_charge": None, "comb": None}
    assert ch["sources"] == {"reported": 2, "comb": 0, "dataset_value": 0, "none": 6}
    assert ch["comb_status"] == {"not_run": 8}


def test_mass_uses_isolation_mz_when_the_file_has_no_mono(tmp_path, monkeypatch):
    res, rows, _ = _convert(tmp_path, monkeypatch)
    assert rows[12][2:6] == ["1235.29", "0.0", "1", "1235.29"]   # was 0.0: the mono m/z column stays as the file gave it
    assert res["summary"]["charge"]["mass_mz_source"] == {"monoisotopic": 1, "isolation": 7}


def test_comb_fills_unreported_charge_and_never_overwrites(tmp_path, monkeypatch):
    res, rows, raw = _convert(tmp_path, monkeypatch, charge_comb=True)
    assert rows[2][4:6] == ["2", str(1039.10 * 2 + PROTON)]                # comb z = 2, negative mode, isolation m/z
    assert rows[4][4] == "2"
    assert rows[3][4] == "3" and rows[12][4] == "1"                       # reported charges kept
    for scan in (6, 8, 9, 11):
        assert rows[scan][4:6] in (["0", "0"], ["0.0", "0"])
    ch = res["summary"]["charge"]
    assert ch["sources"] == {"reported": 2, "comb": 2, "dataset_value": 0, "none": 4}
    assert ch["comb_status"] == {"called": 2, "not_run": 2, "no_profile_ms1": 1, "ft_ms1_not_tried": 1,
                                 "no_parent_ms1": 1, "below_margin": 1}
    assert ch["settings"]["comb"]["margin"] == 0.1
    by_scan = {e["scan"]: e for e in ch["scans"]}
    assert sorted(by_scan) == [2, 4, 6, 8, 9, 11]                          # unreported scans only
    assert by_scan[2]["charge_source"] == "comb" and by_scan[2]["comb_z"] == 2 and by_scan[2]["comb_margin"] >= 0.1
    assert by_scan[11]["comb_status"] == "below_margin" and by_scan[11]["comb_z"] is None


def test_each_parent_ms1_is_read_once_as_profile(tmp_path, monkeypatch):
    _, _, raw = _convert(tmp_path, monkeypatch, charge_comb=True)
    reads = _ms1_reads(raw)
    assert sorted(reads) == [("masslist", 1, False), ("masslist", 10, False)]   # 2 and 4 share scan 1


def test_dataset_charge_fills_what_the_comb_leaves(tmp_path, monkeypatch):
    res, rows, _ = _convert(tmp_path, monkeypatch, charge_comb=True, dataset_charge=1)
    assert rows[2][4] == "2" and rows[4][4] == "2"                        # comb first
    assert [rows[s][4] for s in (6, 8, 9, 11)] == ["1", "1", "1", "1"]
    assert rows[6][5] == "800.0"                                          # z = 1: mass = m/z
    assert rows[3][4] == "3" and rows[12][4] == "1"
    ch = res["summary"]["charge"]
    assert ch["sources"] == {"reported": 2, "comb": 2, "dataset_value": 4, "none": 0}
    assert ch["comb_vs_dataset_disagree"] == 2


def test_dataset_charge_alone(tmp_path, monkeypatch):
    res, rows, raw = _convert(tmp_path, monkeypatch, dataset_charge=-2)    # sign ignored: polarity from the scan
    assert [rows[s][4] for s in (2, 4, 6, 8, 9, 11)] == ["2"] * 6
    assert rows[3][4] == "3"
    assert _ms1_reads(raw) == []
    assert res["summary"]["charge"]["settings"]["dataset_charge"] == 2


@pytest.mark.parametrize("bad", [0, 1.5, 2.0, True, "2"])
def test_invalid_dataset_charge_stops_before_reading(tmp_path, monkeypatch, bad):
    monkeypatch.setattr(X, "MSFileReader", ChargeRaw, raising=False)
    ChargeRaw.last = None
    with pytest.raises(ValueError, match="dataset charge"):
        X.convert_raw_to_csv(str(tmp_path / "synthetic.raw"), outdir=str(tmp_path), dataset_charge=bad)
    assert ChargeRaw.last is None


def test_cli_convert(tmp_path, monkeypatch):
    import json
    import msp_extract_cli as CLI
    monkeypatch.setattr(X, "MSFileReader", ChargeRaw, raising=False)
    raw = tmp_path / "SECRET_NAME_run.raw"
    raw.write_bytes(b"")
    out = tmp_path / "cli"
    args = ["convert", "--raw", str(raw), "--out", str(out), "--alias", "ltq_alias", "--charge-comb",
            "--dataset-charge", "1"]
    assert CLI.main(args) == CLI.EXIT_OK
    assert sorted(p.name for p in out.iterdir()) == ["ltq_alias.convert.json", "ltq_alias.reader_diag.json",
                                                     "ltq_alias.reader_diag.txt", "ltq_alias.scan_status.tsv",
                                                     "ms2_ltq_alias.csv", "ms3_ltq_alias.csv"]
    assert not any("SECRET_NAME" in p.read_text() for p in out.iterdir())
    text = (out / "ltq_alias.convert.json").read_text()
    record = json.loads(text)
    assert record["ok"] and record["options"] == {"charge_comb": True, "dataset_charge": 1}
    assert record["module"]["name"] == "msprawextractor" and record["module"]["version"] == X.version
    assert record["summary"]["charge"]["sources"] == {"reported": 2, "comb": 2, "dataset_value": 4, "none": 0}
    assert record["ms2"]["rows"] == 8
    assert record["outcome"] == "complete" and record["exit_code"] == CLI.EXIT_OK     # every charge filled


def test_cli_convert_usage_and_failure(tmp_path, monkeypatch):
    import msp_extract_cli as CLI
    raw = tmp_path / "a.raw"
    raw.write_bytes(b"")
    out = str(tmp_path / "cli")
    assert CLI.main(["convert", "--raw", str(raw), "--out", out, "--dataset-charge", "0"]) == CLI.EXIT_USAGE
    assert CLI.main(["convert", "--raw", str(raw), "--out", out, "--dataset-charge", "1.5"]) == CLI.EXIT_USAGE
    assert CLI.main(["convert", "--raw", str(raw) + ".missing", "--out", out]) == CLI.EXIT_USAGE

    def broken(path):
        raise OSError("cannot open")
    monkeypatch.setattr(X, "MSFileReader", broken, raising=False)
    assert CLI.main(["convert", "--raw", str(raw), "--out", out, "--alias", "x"]) == CLI.EXIT_FAILED
