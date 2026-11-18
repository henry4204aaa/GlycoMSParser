"""Stage 3a: glycan composition enumeration in compnewv4 (TESTS-PUBLIC CP-B, 2026-11-18).

`NGlaunch` (N-glycans) and `OGlaunch` (O-glycans) build a composition library under user constraints and write it
as a CSV (Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc, HexA, Mode, SO3, PO3H, Mass).  The existing tests of this module cover
the mass formulas (`test_native_precursor_mass.py`, `test_b05_kdn_mass_smoke.py`); these cover WHICH compositions are
produced.  Small, fast constraint sets, run into a temporary folder; nothing is read from disk but the module itself.

The pinned sets (26 complex/hybrid-range compositions, the 12-composition box) are the behaviour of the v1.10 code
recorded as a regression guard; they are not a statement of what the biology allows.  Hybrid mode is not tested:
`hybrid=True` changed nothing in the ranges used here (the code comments call the hybrid arguments a placeholder).
O-glycans are covered only at the level of core selection and range filtering.
Composition tuple order everywhere: (Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc).
"""
import contextlib
import io
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import compnewv4 as C  # noqa: E402

KEYS = ["Hex", "HexNAc", "NeuAc", "NeuGc", "KDN", "Fuc"]
PROTON = 1.007276

# No sialic acid, KDN or fucose: only the Hex / HexNAc skeleton is enumerated
PLAIN = dict(Hex_range=[2, 12], HexNAc_range=[2, 6], Neu5Ac_range=[0, 0], Neu5Gc_range=[0, 0], KDN_range=[0, 0],
             Fucose_range=[0, 0], HexA_range=[0, 0], allow5ac=False, allow5gc=False, allowfuc=False, corefuc=False,
             derivatization="Native", reduced=False)
SIALYLATED = dict(PLAIN, Hex_range=[3, 8], HexNAc_range=[2, 5], Neu5Ac_range=[0, 2], allow5ac=True)

# complex-type compositions of PLAIN with two antennae and the default core (high-mannose on)
PLAIN_PINNED = {
    (3, 2), (3, 3), (3, 4), (4, 2), (4, 3), (4, 4), (4, 5), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6),
    (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (7, 2), (7, 4), (7, 5), (7, 6), (8, 2), (8, 5), (8, 6),
    (9, 2), (9, 6),
}


def _comp_set(df):
    return {tuple(int(x) for x in row) for row in df[KEYS].to_numpy()}


@pytest.fixture
def ng(tmp_path):
    counter = {"n": 0}

    def run(**flags):
        counter["n"] += 1
        name = tmp_path / f"lib{counter['n']}"
        with contextlib.redirect_stdout(io.StringIO()):           # the module prints progress lines
            C.NGlaunch(dict(flags, debug=False), filename=str(name))
        return pd.read_csv(f"{name}.csv")
    return run


@pytest.fixture
def og(tmp_path):
    counter = {"n": 0}

    def run(coretype, **flags):
        counter["n"] += 1
        name = tmp_path / f"og{counter['n']}"
        with contextlib.redirect_stdout(io.StringIO()):
            C.OGlaunch(dict(flags, debug=False), coretype=coretype, filename=str(name))
        return pd.read_csv(f"{name}.csv")
    return run


# ---------------------------------------------------------------- N-glycans: file layout and masses

def test_box_of_small_ranges_is_filled_completely_and_the_file_layout_is_fixed(ng):
    df = ng(**dict(PLAIN, Hex_range=[3, 6], HexNAc_range=[2, 4]))
    assert list(df.columns) == ["Hex", "HexNAc", "NeuAc", "NeuGc", "KDN", "Fuc", "HexA", "Mode", "SO3", "PO3H", "Mass"]
    assert len(df) == 12 and _comp_set(df) == {(h, n, 0, 0, 0, 0) for h in range(3, 7) for n in range(2, 5)}
    assert set(df["Mode"]) == {"+"} and set(df["SO3"]) == {0} and set(df["PO3H"]) == {0} and set(df["HexA"]) == {0}


def test_masses_in_the_file_are_precursormassv4_rounded_to_four_decimals(ng):
    df = ng(**dict(PLAIN, Hex_range=[3, 6], HexNAc_range=[2, 4]))
    for _, r in df.iterrows():
        comp = tuple(int(r[k]) for k in KEYS)
        assert r["Mass"] == pytest.approx(C.precursormassv4(comp, deri="Native", reduced=False, mode="+"), abs=6e-5)
    # independent anchors (Native, [M+H]+): H3N2 and H5N4 (H5N4: GlycoMod-checked 1641.599 in test_native_precursor_mass)
    by = {tuple(int(r[k]) for k in KEYS[:2]): r["Mass"] for _, r in df.iterrows()}
    assert by[(3, 2)] == pytest.approx(911.335, abs=1e-3)
    assert by[(5, 4)] == pytest.approx(1641.599, abs=2e-3)


def test_negative_mode_shifts_every_mass_by_two_protons_and_keeps_the_compositions(ng):
    pos = ng(**dict(PLAIN, Hex_range=[3, 6], HexNAc_range=[2, 4]))
    neg = ng(**dict(PLAIN, Hex_range=[3, 6], HexNAc_range=[2, 4], mode="-"))
    assert set(neg["Mode"]) == {"-"} and _comp_set(neg) == _comp_set(pos)
    merged = pos.merge(neg, on=KEYS, suffixes=("_pos", "_neg"))
    assert len(merged) == len(pos)
    assert ((merged["Mass_pos"] - merged["Mass_neg"]) - 2 * PROTON).abs().max() < 2e-4


def test_derivatization_changes_masses_not_compositions(ng):
    native = ng(**dict(PLAIN, Hex_range=[3, 6], HexNAc_range=[2, 4]))
    perme = ng(**dict(PLAIN, Hex_range=[3, 6], HexNAc_range=[2, 4], derivatization="PerMe", reduced=True))
    assert _comp_set(perme) == _comp_set(native)
    m = native.merge(perme, on=KEYS, suffixes=("_n", "_p"))
    assert (m["Mass_p"] > m["Mass_n"]).all()
    for _, r in perme.iterrows():
        comp = tuple(int(r[k]) for k in KEYS)
        assert r["Mass"] == pytest.approx(C.precursormassv4(comp, deri="PerMe", reduced=True, mode="+"), abs=6e-5)


def test_enumeration_is_deterministic(ng):
    a, b = ng(**PLAIN), ng(**PLAIN)
    assert a.equals(b)


# ---------------------------------------------------------------- N-glycans: constraints

def test_plain_skeleton_set_is_pinned(ng):
    df = ng(**PLAIN)
    assert {(h, n) for h, n, *_ in _comp_set(df)} == PLAIN_PINNED
    assert len(df) == 26 and not df.duplicated(KEYS).any()


def test_high_mannose_flag_adds_exactly_man6_to_man9(ng):
    with_hm = _comp_set(ng(**dict(PLAIN, highman=True)))
    without = _comp_set(ng(**dict(PLAIN, highman=False)))
    assert without < with_hm
    assert with_hm - without == {(h, 2, 0, 0, 0, 0) for h in (6, 7, 8, 9)}


def test_more_antennae_only_add_compositions(ng):
    sets = {n: _comp_set(ng(**dict(PLAIN, arm_count=n))) for n in (1, 2, 3)}
    assert (len(sets[1]), len(sets[2]), len(sets[3])) == (15, 26, 32)
    assert sets[1] < sets[2] < sets[3]


def test_ranges_filter_the_enumeration_and_tightening_only_removes(ng):
    wide = _comp_set(ng(**PLAIN))
    tight_flags = dict(PLAIN, Hex_range=[4, 7], HexNAc_range=[3, 5])
    tight = _comp_set(ng(**tight_flags))
    assert tight and all(4 <= h <= 7 and 3 <= n <= 5 for h, n, *_ in tight)
    assert tight == {c for c in wide if 4 <= c[0] <= 7 and 3 <= c[1] <= 5}


def test_sialic_acid_range_and_switch(ng):
    full = _comp_set(ng(**SIALYLATED))
    assert {c[2] for c in full} == {0, 1, 2} and {c[3] for c in full} == {0} and len(full) == 50
    up_to_one = _comp_set(ng(**dict(SIALYLATED, Neu5Ac_range=[0, 1])))
    assert up_to_one < full and {c[2] for c in up_to_one} == {0, 1} and len(up_to_one) == 36
    at_least_one = _comp_set(ng(**dict(SIALYLATED, Neu5Ac_range=[1, 2])))
    assert at_least_one == {c for c in full if c[2] >= 1} and len(at_least_one) == 30
    switched_off = _comp_set(ng(**dict(SIALYLATED, allow5ac=False)))
    assert {c[2] for c in switched_off} == {0}


def test_neu5gc_follows_its_own_switch_and_range(ng):
    gc = _comp_set(ng(**dict(SIALYLATED, Neu5Ac_range=[0, 0], Neu5Gc_range=[0, 1], allow5gc=True)))
    assert {c[3] for c in gc} == {0, 1} and {c[2] for c in gc} == {0} and len(gc) == 36


def test_fucose_range_and_switch(ng):
    base = dict(SIALYLATED, Neu5Ac_range=[0, 0])
    fuc = _comp_set(ng(**dict(base, Fucose_range=[0, 2], allowfuc=True, corefuc=True)))
    assert {c[5] for c in fuc} == {0, 1, 2} and len(fuc) == 56
    assert {c[5] for c in _comp_set(ng(**dict(base, Fucose_range=[0, 2], allowfuc=False, corefuc=False)))} == {0}


# ---------------------------------------------------------------- O-glycans

OG_FLAGS = dict(Hex_range=[0, 4], HexNAc_range=[1, 4], Neu5Ac_range=[0, 0], Neu5Gc_range=[0, 0], KDN_range=[0, 0],
                Fucose_range=[0, 1], HexA_range=[0, 0], derivatization="Native", reduced=False)


def test_o_glycan_core_selection_and_ranges(og):
    tn = _comp_set(og([0], **OG_FLAGS))
    assert (0, 1, 0, 0, 0, 0) in tn                                            # the Tn antigen itself: one HexNAc
    sub = _comp_set(og([1, 2], **OG_FLAGS))
    everything = _comp_set(og(None, **OG_FLAGS))                               # None = every core type
    assert tn and sub and tn | sub <= everything
    for c in everything:
        assert 0 <= c[0] <= 4 and 1 <= c[1] <= 4 and 0 <= c[5] <= 1 and c[2] == c[3] == c[4] == 0


def test_o_glycan_core_type_must_be_a_list(og):
    with pytest.raises(TypeError):
        og(1, **OG_FLAGS)
