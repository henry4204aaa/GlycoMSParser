"""Native (underivatized) precursor mass in compnewv4.precursormassv4 (V14-FOLLOWUP CP-A, 2026-10-27).

Before CP-A the Native branch returned 0, so a CGA library built for Native metadata had no usable
masses, and the GUI value "Native(Reduced)" fell through to the PerMe-reduced default.

Reference values: V14-FOLLOWUP kickoff (Cowork), residue constants checked against GlycoMod by Henry
2026-10-27. Pass rule from the kickoff: Native within 0.002 Da; PerMe identical to the base output.
Composition tuple order: (Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc).
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import compnewv4
from compnewv4 import precursormassv4, _normalize_derivatization

TOL = 0.002  # Da, kickoff pass rule for Native rows

H5N4 = (5, 4, 0, 0, 0, 0)
F1H5N4S2 = (5, 4, 2, 0, 0, 1)
H3N2 = (3, 2, 0, 0, 0, 0)


def _native(comp, reduced=False, mode="+", **kw):
    return precursormassv4(comp, deri="Native", reduced=reduced, mode=mode, **kw)


def _neutral(comp, deri, reduced, **kw):
    # the mode term is symmetric (+p / -p), so the neutral mass is the midpoint of the two ions
    return (precursormassv4(comp, deri=deri, reduced=reduced, mode="+", **kw)
            + precursormassv4(comp, deri=deri, reduced=reduced, mode="-", **kw)) / 2


@pytest.mark.parametrize("comp, reduced, mode, neutral_ref, ion_ref", [
    (H5N4, False, "+", 1640.592, 1641.599),
    (H5N4, True, "+", 1642.608, 1643.615),
    (F1H5N4S2, False, "+", 2368.841, 2369.848),
    (H3N2, False, "+", 910.328, 911.335),
    (H5N4, False, "-", 1640.592, 1639.585),
    (H5N4, True, "-", 1642.608, 1641.6005),
    (F1H5N4S2, True, "-", 2370.857, 2369.8493),
])
def test_native_reference_rows(comp, reduced, mode, neutral_ref, ion_ref):
    assert abs(_native(comp, reduced, mode) - ion_ref) <= TOL
    assert abs(_neutral(comp, "Native", reduced) - neutral_ref) <= TOL


def test_native_proton_term():
    # Henry 2026-10-27: Native uses 1.00728 (PerMe keeps 1.0073)
    assert _native(H5N4, mode="+") - _native(H5N4, mode="-") == pytest.approx(2 * 1.00728, abs=1e-9)


@pytest.mark.parametrize("comp, kw, residue", [
    ((5, 4, 0, 0, 1, 0), {}, 250.06887),     # KDN
    ((5, 4, 0, 1, 0, 0), {}, 307.09033),     # NeuGc
    ((5, 4, 1, 0, 0, 0), {}, 291.09542),     # NeuAc
    ((5, 4, 0, 0, 0, 1), {}, 146.05791),     # Fuc
    (H5N4, {"hexA": 1}, 176.03209),          # HexA
    (H5N4, {"so3": 1}, 79.95682),            # SO3
    (H5N4, {"po3h": 1}, 79.96633),           # HPO3
])
def test_native_increments(comp, kw, residue):
    assert _native(comp, **kw) - _native(H5N4) == pytest.approx(residue, abs=1e-9)


def test_native_groups_stay_singly_charged():
    # Henry 2026-10-27: flat additions, no PerMe so3 == 2 halving
    assert _native(H5N4, so3=2) - _native(H5N4) == pytest.approx(2 * 79.95682, abs=1e-9)
    assert _native(H5N4, po3h=2) - _native(H5N4) == pytest.approx(2 * 79.96633, abs=1e-9)


def test_reduced_offset():
    assert _native(H5N4, reduced=True) - _native(H5N4) == pytest.approx(2.01565, abs=1e-9)


def test_near_isobaric_pairs():
    # an ion trap cannot separate these: derivatization, reduction and polarity must come from metadata
    assert abs(_native(H5N4, True, "-") - _native(H5N4, False, "+")) == pytest.approx(0.0011, abs=0.0002)
    assert abs(_native(F1H5N4S2, True, "-") - _native(F1H5N4S2, False, "+")) == pytest.approx(0.0011, abs=0.0002)


def test_perme_unchanged():
    # base output at 53f3b9b (repr); see also tests/test_b05_kdn_mass_smoke.py
    assert repr(precursormassv4(H5N4, deri="PerMe", reduced=False, mode="+")) == "2048.0532900000003"
    assert repr(precursormassv4(H5N4, deri="PerMe", reduced=True, mode="+")) == "2064.0845900000004"


def test_unknown_deri_still_zero():
    assert precursormassv4(H5N4, deri="2-AB") == 0


@pytest.mark.parametrize("value, reduced_flag, expected", [
    ("Native", False, ("Native", False)),
    ("Native(Reduced)", True, ("Native", True)),
    ("Native(Reduced)", None, ("Native", True)),
    ("PerMe(Freeend)", False, ("PerMe", False)),
    ("PerMe(Reduced)", True, ("PerMe", True)),
    ("None", False, ("PerMe", False)),
    ("", False, ("PerMe", False)),
])
def test_gui_values_normalize(value, reduced_flag, expected):
    # flags as mspfileloaderv14 builds them for NGlaunch/OGlaunch (derivatization_type + reduced)
    flags = {"derivatization_type": value}
    if reduced_flag is not None:
        flags["reduced"] = reduced_flag
    assert _normalize_derivatization(flags) == expected


def test_native_library_csv_has_masses(tmp_path):
    out = tmp_path / "native_lib"
    compnewv4.save_glycan_pseudocomp_to_csv([H5N4, F1H5N4S2], filename=str(out), derivatization="Native",
                                            reduced=True, mode="-")
    rows = (tmp_path / "native_lib.csv").read_text().splitlines()
    masses = [float(r.split(",")[-1]) for r in rows[1:]]
    assert len(masses) == 2 and all(m > 0 for m in masses)
    assert f"{_native(H5N4, True, '-'):.4f}" in rows[1] + rows[2]
