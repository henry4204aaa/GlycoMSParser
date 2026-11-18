"""msp_CLI_precursor_mass (V14-FOLLOWUP CP-A.2): the CLI prints what compnewv4.precursormassv4 returns.

The CLI is a runnable statement of how GlycoMSP forms precursor masses; it must never hold a second mass
implementation. These tests compare its output with a direct precursormassv4 call on the same inputs, run the
CP-A.1 self-test, and check that it imports no GUI module.
"""
import csv
import json
import os
import subprocess
import sys

import pytest

SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, SRC)

import compnewv4  # noqa: E402
import msp_CLI_precursor_mass as cli  # noqa: E402

CLI = os.path.join(SRC, "msp_CLI_precursor_mass.py")


def _run(*args):
    return subprocess.run([sys.executable, CLI, *args], capture_output=True, text=True)


@pytest.mark.parametrize("label, comp, deri, reduced, kw", [
    ("H5N4F1S2", (5, 4, 2, 0, 0, 1), "native", False, {}),
    ("F1H5N4S2", (5, 4, 2, 0, 0, 1), "native", True, {}),
    ("H5N4", (5, 4, 0, 0, 0, 0), "perme", False, {}),
    ("H5N4", (5, 4, 0, 0, 0, 0), "perme", True, {}),
    ("H4N4KDN1G1", (4, 4, 0, 1, 1, 0), "native", False, {}),
    ("H5N4A1s1", (5, 4, 0, 0, 0, 0), "native", True, {"hexA": 1, "so3": 1}),
    ("H5N4p1", (5, 4, 0, 0, 0, 0), "perme", False, {"po3h": 1}),
])
def test_cli_json_equals_precursormassv4(label, comp, deri, reduced, kw):
    args = ["--composition", label, "--deri", deri, "--json"] + (["--reduced"] if reduced else [])
    proc = _run(*args)
    assert proc.returncode == 0, proc.stderr
    out = json.loads(proc.stdout)
    name = {"native": "Native", "perme": "PerMe"}[deri]
    counts = cli.parse_composition(label)
    for mode in "+-":
        direct = compnewv4.precursormassv4(comp, deri=name, reduced=reduced, mode=mode, **kw)
        assert out["ions"][cli.ion_label(counts, name, mode)] == direct   # exact: same function, full float
    assert out["function"] == "compnewv4.precursormassv4"


def test_group_flags_match_label_letters():
    a = json.loads(_run("--composition", "H5N4", "--so3", "1", "--hexa", "1", "--deri", "native", "--json").stdout)
    b = json.loads(_run("--composition", "H5N4A1s1", "--deri", "native", "--json").stdout)
    assert a["ions"] == b["ions"]


def test_mode_selects_one_ion():
    out = json.loads(_run("--composition", "H5N4", "--deri", "native", "--reduced", "--mode", "-", "--json").stdout)
    assert list(out["ions"]) == ["[M-H]-"]
    assert out["ions"]["[M-H]-"] == compnewv4.precursormassv4((5, 4, 0, 0, 0, 0), deri="Native", reduced=True, mode="-")


def test_selftest_passes():
    proc = _run("--selftest")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "12 of 12 rows ok" in proc.stdout


def test_selftest_detects_a_mismatch(monkeypatch):
    rows = list(cli.ACCEPTANCE_ROWS)
    rows[2] = rows[2][:5] + (rows[2][5] + 0.01,) + rows[2][6:]   # Native H5N4 reference off by 10 mDa
    monkeypatch.setattr(cli, "ACCEPTANCE_ROWS", rows)
    assert cli.main(["--selftest"]) == 1


def test_selftest_perme_rows_are_exact():
    for row in cli.ACCEPTANCE_ROWS:
        if row[1] == "PerMe":
            r = cli.evaluate(row)
            assert r["result"] == "identical"


def test_table_default_rows_csv_and_md(tmp_path):
    out_csv, out_md = tmp_path / "t.csv", tmp_path / "t.md"
    assert _run("--table", str(out_csv)).returncode == 0
    assert _run("--table", str(out_md)).returncode == 0
    rows = list(csv.DictReader(out_csv.open()))
    assert len(rows) == len(cli.ACCEPTANCE_ROWS)
    assert out_md.read_text().count("\n") == len(cli.ACCEPTANCE_ROWS) + 2


def test_table_custom_compositions(tmp_path):
    out = tmp_path / "c.csv"
    proc = _run("--table", str(out), "--composition", "H5N4", "--composition", "H3N2", "--deri", "native")
    assert proc.returncode == 0, proc.stderr
    rows = list(csv.DictReader(out.open()))
    assert [(r["composition"], r["mode"]) for r in rows] == [("H5N4", "+"), ("H5N4", "-"), ("H3N2", "+"), ("H3N2", "-")]
    direct = compnewv4.precursormassv4((3, 2, 0, 0, 0, 0), deri="Native", reduced=False, mode="-")
    assert rows[3]["ion_mz"] == f"{direct:.5f}"


def test_perme_with_groups_has_no_neutral():
    out = json.loads(_run("--composition", "H5N4s2", "--json").stdout)
    assert out["neutral_mass"] is None


@pytest.mark.parametrize("args", [
    ["--composition", "H5X4"],
    ["--composition", "H5H4"],
    ["--composition", "H5N4s1", "--so3", "2"],
    ["--composition", "H5N4", "--mode", "x"],
    ["--composition", "H5N4", "--deri", "2ab"],
    [],
])
def test_usage_errors_exit_2(args):
    assert _run(*args).returncode == 2


def test_cli_imports_no_gui():
    code = (f"import sys; sys.path.insert(0, {SRC!r}); import msp_CLI_precursor_mass as c; "
            "c.main(['--composition', 'H5N4']); "
            "bad = [m for m in ('mspfileloaderv14', 'tkinter', 'pandas', 'numpy') if m in sys.modules]; "
            "print('BAD', bad); sys.exit(1 if bad else 0)")
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr


# ---------------------------------------------------------------- Codex round 1 (V14F-F3, F4, F8)

def _table(tmp_path, *args):
    out = tmp_path / "t.csv"
    proc = _run("--table", str(out), *args)
    assert proc.returncode == 0, proc.stderr
    return list(csv.DictReader(out.open()))


@pytest.mark.parametrize("flags, label", [
    (["--so3", "1"], "H5N4s1"), (["--hexa", "1"], "H5N4A1"), (["--po3h", "1"], "H5N4p1"),
    (["--so3", "1", "--hexa", "2"], "H5N4A2s1"),
])
@pytest.mark.parametrize("deri", ["native", "perme"])
def test_table_group_flags_equal_label_form(tmp_path, flags, label, deri):
    """V14F-F3: --table with --composition applied no group flag before 0.2."""
    a = _table(tmp_path, "--composition", "H5N4", "--composition", "H3N2", "--deri", deri, *flags)
    b = _table(tmp_path, "--composition", label, "--composition", label.replace("H5N4", "H3N2"), "--deri", deri)
    assert a == b and a[0]["composition"] == label


def test_table_group_flags_with_mode(tmp_path):
    rows = _table(tmp_path, "--composition", "H5N4", "--deri", "native", "--so3", "1", "--mode", "-")
    direct = compnewv4.precursormassv4((5, 4, 0, 0, 0, 0), deri="Native", reduced=False, mode="-", so3=1)
    assert len(rows) == 1 and rows[0]["ion_mz"] == f"{direct:.5f}" and rows[0]["ion_mz"] == "1719.54168"


@pytest.mark.parametrize("args", [
    ["--composition", "H5N4s1", "--so3", "2"],              # conflicts with the label
    ["--composition", "H5N4", "--so3", "-1"],               # negative
    ["--composition", "H5N4s3"],                            # PerMe: no rule above 2
    ["--composition", "H5N4", "--po3h", "3"],
])
@pytest.mark.parametrize("table", [False, True])
def test_table_and_single_share_validation(tmp_path, args, table):
    extra = ["--table", str(tmp_path / "x.csv")] if table else []
    assert _run(*args, *extra).returncode == 2
    assert not (tmp_path / "x.csv").exists()


@pytest.mark.parametrize("args", [["--selftest", "--so3", "1"], ["--table", "t.md", "--hexa", "1"]])
def test_group_flags_refused_for_fixed_rows(args):
    assert _run(*args).returncode == 2


def test_native_has_no_group_limit():
    out = json.loads(_run("--composition", "H5N4s3", "--deri", "native", "--mode", "-", "--json").stdout)
    assert out["ions"] == {"[M-H]-": compnewv4.precursormassv4((5, 4, 0, 0, 0, 0), deri="Native", mode="-", so3=3)}


@pytest.mark.parametrize("label, labels", [
    ("H5N4", {"+": "[M+H]+", "-": "[M-H]-"}),
    ("H5N4s1", {"+": "m/z (z=+1)", "-": "m/z (z=-1)"}),
    ("H5N4p2", {"+": "m/z (z=+1)", "-": "m/z (z=-1)"}),
    ("H5N4s2", {"+": "[M+2H]2+", "-": "[M-2H]2-"}),
])
def test_perme_group_labels(label, labels):
    """V14F-F4 (Henry 2026-10-27): PerMe SO3/PO3H values are an m/z with fixed corrections; two SO3 is z = 2."""
    out = json.loads(_run("--composition", label, "--json").stdout)
    assert list(out["ions"]) == [labels["+"], labels["-"]]
    for mode in "+-":
        counts = cli.parse_composition(label)
        assert cli.ion_label(counts, "PerMe", mode) == labels[mode]
        assert cli.ion_label(counts, "Native", mode) == {"+": "[M+H]+", "-": "[M-H]-"}[mode]


def test_perme_two_so3_value_is_unchanged():
    """The label changed in 0.2, the number did not (base output at 53f3b9b)."""
    out = json.loads(_run("--composition", "H5N4s2", "--mode", "-", "--json").stdout)
    assert out["ions"] == {"[M-2H]2-": compnewv4.precursormassv4((5, 4, 0, 0, 0, 0), mode="-", so3=2)}
    assert repr(out["ions"]["[M-2H]2-"]) == "1088.45686"


def test_selftest_row_without_reference_fails(monkeypatch):
    """V14F-F8: a Native row with no reference used to pass on ion > 0."""
    rows = list(cli.ACCEPTANCE_ROWS)
    rows[6] = rows[6][:4] + (None, None) + rows[6][6:]
    monkeypatch.setattr(cli, "ACCEPTANCE_ROWS", rows)
    assert cli.main(["--selftest"]) == 1


@pytest.mark.parametrize("i", [6, 7, 8])
def test_selftest_rows_7_to_9_detect_a_wrong_residue(monkeypatch, i):
    """V14F-F8: a wrong KDN / NeuGc / HexA contribution (here +0.01 Da) now fails its row."""
    real = compnewv4.precursormassv4
    monkeypatch.setattr(compnewv4, "precursormassv4", lambda *a, **k: real(*a, **k) + 0.01)
    assert cli.evaluate(cli.ACCEPTANCE_ROWS[i])["ok"] is False


# Elemental monoisotopic masses (C 12 by definition); independent of the residue constants in compnewv4.
_EL = {"C": 12.0, "H": 1.00782503207, "N": 14.0030740048, "O": 15.99491461956, "S": 31.97207100, "P": 30.97376163}
_ELECTRON = 0.00054857990946
_RESIDUE = {"Hex": "C6H10O5", "HexNAc": "C8H13NO5", "Fuc": "C6H10O4", "NeuAc": "C11H17NO8", "NeuGc": "C11H17NO9",
            "KDN": "C9H14O8", "HexA": "C6H8O6", "SO3": "SO3", "PO3H": "HPO3"}


def _formula_mass(formula):
    import re
    return sum(_EL[el] * int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula))


@pytest.mark.parametrize("label, reduced", [
    ("H5N4KDN1", False), ("H5N4G1", False), ("H5N4A1", False), ("F1H5N4S2", True), ("H5N4s1p1", False),
])
def test_native_neutral_from_elemental_formulas(label, reduced):
    counts = cli.parse_composition(label)
    expected = _formula_mass("H2O") + (_formula_mass("H2") if reduced else 0)
    expected += sum(n * _formula_mass(_RESIDUE[k]) for k, n in counts.items() if n)
    neutral, ions = cli.masses(counts, "Native", reduced)
    proton = _EL["H"] - _ELECTRON
    assert abs(neutral - expected) < 0.0005
    assert abs(ions["+"] - (expected + proton)) < 0.0005 and abs(ions["-"] - (expected - proton)) < 0.0005


def test_rows_7_to_9_references_agree_with_elemental_formulas():
    for row in cli.ACCEPTANCE_ROWS[6:9]:
        counts = cli.parse_composition(row[0])
        expected = _formula_mass("H2O") + sum(n * _formula_mass(_RESIDUE[k]) for k, n in counts.items() if n)
        assert abs(row[4] - expected) < 0.002, row


# ---------------------------------------------------------------- Codex round 2 / Henry 2026-10-27: charge rule

@pytest.mark.parametrize("args", [
    ["--composition", "H5N4s1p1"], ["--composition", "H5N4s2p1"], ["--composition", "H5N4s1p2"],
    ["--composition", "H5N4s2p2"], ["--composition", "H5N4s1", "--po3h", "1"], ["--composition", "H5N4p2", "--so3", "2"],
])
@pytest.mark.parametrize("mode", [[], ["--mode", "-"], ["--mode", "+"]])
def test_perme_mixed_so3_po3h_refused(tmp_path, args, mode):
    """Each SO3/PO3H carries one charge; compnewv4 models one group type at a time, so mixed is refused."""
    for extra in ([], ["--json"], ["--table", str(tmp_path / "m.csv")]):
        proc = _run(*args, *mode, *extra)
        assert proc.returncode == 2 and "both SO3 and PO3H" in proc.stderr
    assert not (tmp_path / "m.csv").exists()


def test_native_mixed_so3_po3h_still_flat():
    out = json.loads(_run("--composition", "H5N4s1p1", "--deri", "native", "--json").stdout)
    for mode, ion in (("+", "[M+H]+"), ("-", "[M-H]-")):
        assert out["ions"][ion] == compnewv4.precursormassv4((5, 4, 0, 0, 0, 0), deri="Native", mode=mode,
                                                             so3=1, po3h=1)


@pytest.mark.parametrize("label, kw, z", [
    ("H5N4s1", {"so3": 1}, 1), ("H5N4p1", {"po3h": 1}, 1), ("H5N4p2", {"po3h": 2}, 1), ("H5N4s2", {"so3": 2}, 2),
])
def test_perme_supported_groups_value_and_charge(label, kw, z):
    """The four supported PerMe group cases: frozen values, labelled with the charge Henry's rule gives."""
    out = json.loads(_run("--composition", label, "--json").stdout)
    want = {1: {"+": "m/z (z=+1)", "-": "m/z (z=-1)"}, 2: {"+": "[M+2H]2+", "-": "[M-2H]2-"}}[z]
    for mode in "+-":
        assert out["ions"][want[mode]] == compnewv4.precursormassv4((5, 4, 0, 0, 0, 0), mode=mode, **kw)
