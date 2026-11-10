"""Precursor charge from the parent MS1 isotope pattern: the anchor-scanned comb (READER-DIAG CP-C.5).

version 0.11, 20261023. Standard library only. (v0.11 R2-F10: a profile with unequal array lengths or a
non-finite value returns STATUS_INVALID_PROFILE; nothing is repaired.)

For each charge z the comb compares the MS1 profile at the positions where that charge's isotopes would
sit ("on": anchor + k/z, k = 0..3) with the positions halfway between them ("off": anchor + (k + 1/2)/z,
k = -1..3; the k = -1 tooth sits just below the anchor, where a monoisotopic peak has nothing). The
anchor is tried at every profile point within the isolation m/z +/- ANCHOR_HALF_WIDTH and the best
score is kept per z. Scores are divided by the highest profile point within the isolation m/z
+/- IMAX_HALF_WIDTH. A charge is called only when the best z beats the runner-up by MARGIN.

Scope and evidence: designed and tried on ion-trap profile MS1 (0.083 m/z grid) only, on 5 spectra whose
charges were known (4 called correctly, 1 below the margin, none wrong). The margin was set after seeing
them. The comb returns a charge, never a monoisotopic m/z: its best anchor may be another isotope.
Design: handoff/READER_DIAG_comb_design_20261015.md.
"""
import math
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field

version = "0.11"
last_update = 20261023

WINDOW_HALF_WIDTH = 4.0      # profile points kept around the isolation m/z
IMAX_HALF_WIDTH = 1.0        # I_max is the highest profile point within this
ANCHOR_HALF_WIDTH = 0.6      # anchors tried within this
CHARGES = (1, 2, 3)
ON_TEETH = (0, 1, 2, 3)      # anchor + k/z
OFF_TEETH = (-1, 0, 1, 2, 3)  # anchor + (k + 1/2)/z
MARGIN = 0.1

STATUS_CALLED = "called"
STATUS_BELOW_MARGIN = "below_margin"
STATUS_NO_SIGNAL = "no_signal"   # no profile point near the isolation m/z, or I_max is 0
STATUS_INVALID_PROFILE = "invalid_profile"   # arrays of unequal length, or a NaN/inf value (R2-F10)


def parameters():
    """The settings as a plain dict, for run records."""
    return {"version": version, "window_half_width": WINDOW_HALF_WIDTH, "imax_half_width": IMAX_HALF_WIDTH,
            "anchor_half_width": ANCHOR_HALF_WIDTH, "charges": list(CHARGES), "on_teeth": list(ON_TEETH),
            "off_teeth": list(OFF_TEETH), "margin": MARGIN}


@dataclass
class CombResult:
    z: object                 # the called charge (int), or None
    status: str
    score: object = None      # best z's score
    margin: object = None     # best minus runner-up
    anchor: object = None     # best z's anchor m/z
    i_max: object = None
    scores: dict = field(default_factory=dict)   # z -> (score, anchor)


class _Profile:
    """Sorted profile points with linear interpolation (0 outside the sampled range)."""

    def __init__(self, mz, intensity):
        pairs = sorted(zip((float(x) for x in mz), (float(y) for y in intensity)))
        self.mz = [p[0] for p in pairs]
        self.intensity = [p[1] for p in pairs]

    def between(self, lo, hi):
        i, j = bisect_left(self.mz, lo), bisect_right(self.mz, hi)
        return list(zip(self.mz[i:j], self.intensity[i:j]))

    def at(self, x):
        mz = self.mz
        if not mz or x < mz[0] or x > mz[-1]:
            return 0.0
        j = bisect_left(mz, x)
        if mz[j] == x:
            return self.intensity[j]
        x0, x1, y0, y1 = mz[j - 1], mz[j], self.intensity[j - 1], self.intensity[j]
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def infer_charge(mz, intensity, isolation_mz):
    """Run the comb on one parent MS1 profile (any length; only the window around isolation_mz is used)."""
    try:
        # one C-speed sum per array: NaN or inf anywhere makes the sum non-finite (an absurd finite overflow too)
        valid = (len(mz) == len(intensity) and math.isfinite(float(isolation_mz))
                 and math.isfinite(sum(map(float, mz))) and math.isfinite(sum(map(float, intensity))))
    except (TypeError, ValueError, OverflowError):
        valid = False
    if not valid:
        return CombResult(None, STATUS_INVALID_PROFILE)
    window = _Profile(mz, intensity).between(isolation_mz - WINDOW_HALF_WIDTH, isolation_mz + WINDOW_HALF_WIDTH)
    profile = _Profile([x for x, _ in window], [y for _, y in window])
    near = profile.between(isolation_mz - IMAX_HALF_WIDTH, isolation_mz + IMAX_HALF_WIDTH)
    i_max = max((y for _, y in near), default=0.0)
    anchors = [x for x, _ in profile.between(isolation_mz - ANCHOR_HALF_WIDTH, isolation_mz + ANCHOR_HALF_WIDTH)]
    if i_max <= 0 or not anchors:
        return CombResult(None, STATUS_NO_SIGNAL, i_max=i_max)

    scores = {}
    for z in CHARGES:
        best = None
        for a in anchors:
            on = sum(profile.at(a + k / z) for k in ON_TEETH) / len(ON_TEETH)
            off = sum(profile.at(a + (k + 0.5) / z) for k in OFF_TEETH) / len(OFF_TEETH)
            s = (on - off) / i_max
            if best is None or s > best[0]:
                best = (s, a)
        scores[z] = best
    ranked = sorted(CHARGES, key=lambda z: -scores[z][0])
    top, second = ranked[0], ranked[1]
    margin = scores[top][0] - scores[second][0]
    called = margin >= MARGIN
    return CombResult(top if called else None, STATUS_CALLED if called else STATUS_BELOW_MARGIN,
                      score=scores[top][0], margin=margin, anchor=scores[top][1], i_max=i_max, scores=scores)
