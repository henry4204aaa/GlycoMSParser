"""READER-DIAG CP-C.5: the MS1 isotope comb (msp_charge_comb) on synthetic profiles (2026-10-15).

Synthetic isotope patterns only: Gaussian peaks at mono + k * 1.00336 / z on a 0.0833 m/z grid, as the
LTQ profile MS1 in the charge probe. No third-party spectra.
"""
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import msp_charge_comb as C  # noqa: E402

GLYCAN_ENVELOPE = (1.0, 0.85, 0.45, 0.18, 0.06)


def profile(mono, z, rel=GLYCAN_ENVELOPE, sigma=0.08, grid=0.0833, height=1000.0, noise=None):
    mz, inten = [], []
    x = mono - 4.2 + 0.031
    while x < mono + 4.2:
        y = sum(r * height * math.exp(-(x - (mono + k * 1.00336 / z)) ** 2 / (2 * sigma ** 2))
                for k, r in enumerate(rel))
        mz.append(x)
        inten.append(y + (noise(x) if noise else 0.0))
        x += grid
    return mz, inten


@pytest.mark.parametrize("z", [1, 2, 3])
@pytest.mark.parametrize("iso_offset", [0.0, 0.15, -0.2, 0.3])
def test_resolved_pattern_gives_its_charge(z, iso_offset):
    mz, inten = profile(1000.0, z)
    r = C.infer_charge(mz, inten, 1000.0 + iso_offset)
    assert (r.z, r.status) == (z, C.STATUS_CALLED)
    assert r.margin >= C.MARGIN


def test_unresolved_pattern_stays_unknown():
    # LTQ-like peak width: z = 3 isotopes 0.33 apart merge, so no charge is called (and none is wrong).
    mz, inten = profile(1000.0, 3, sigma=0.17)
    r = C.infer_charge(mz, inten, 1000.1)
    assert r.z is None and r.status == C.STATUS_BELOW_MARGIN


def test_lone_peak_has_no_charge():
    mz, inten = profile(1000.0, 1, rel=(1.0,))
    r = C.infer_charge(mz, inten, 1000.1)
    assert r.z is None and r.status == C.STATUS_BELOW_MARGIN


@pytest.mark.parametrize("mz, inten, iso", [
    ((), (), 1000.0),
    ((999.9, 1000.0, 1000.1), (0.0, 0.0, 0.0), 1000.0),
    ((990.0, 990.1), (5.0, 9.0), 1000.0),                     # nothing near the isolation m/z
])
def test_no_signal(mz, inten, iso):
    r = C.infer_charge(mz, inten, iso)
    assert r.z is None and r.status == C.STATUS_NO_SIGNAL


def test_point_order_does_not_matter():
    mz, inten = profile(1500.0, 2)
    a = C.infer_charge(mz, inten, 1500.2)
    b = C.infer_charge(list(reversed(mz)), list(reversed(inten)), 1500.2)
    assert (a.z, a.score, a.anchor) == (b.z, b.score, b.anchor)


def test_only_the_window_is_used():
    mz, inten = profile(1000.0, 2)
    far_mz, far_i = profile(1010.0, 1, height=50000.0)
    r = C.infer_charge(mz + far_mz, inten + far_i, 1000.1)
    assert r.z == 2


def test_score_formula_with_the_tooth_below_the_anchor():
    """Scored by hand on a 0.5 m/z grid (interpolation exact): on = a + k/z, k = 0..3; off = a + (k + 1/2)/z, k = -1..3."""
    mz = [998.0 + 0.5 * i for i in range(11)]                  # 998.0 .. 1003.0; anchors 999.5, 1000.0, 1000.5
    inten = [0.0] * 11
    inten[4] = 100.0                                           # 1000.0 = I_max
    inten[3] = 40.0                                            # 999.5 = 1000.0 - 1/2, the k = -1 off tooth for z = 1
    inten[6] = 60.0                                            # 1001.0 = 1000.0 + 1
    r = C.infer_charge(mz, inten, 1000.0)
    expected_z1 = ((100.0 + 60.0) / 4 - 40.0 / 5) / 100.0     # 0.32; without the k = -1 tooth it would be 0.40
    assert r.scores[1] == (pytest.approx(expected_z1), 1000.0)


def test_parameters_match_the_agreed_design():
    p = C.parameters()
    assert (p["window_half_width"], p["imax_half_width"], p["anchor_half_width"]) == (4.0, 1.0, 0.6)
    assert p["charges"] == [1, 2, 3] and p["margin"] == 0.1
    assert p["on_teeth"] == [0, 1, 2, 3] and p["off_teeth"] == [-1, 0, 1, 2, 3]
