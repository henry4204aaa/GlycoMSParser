"""B-24 regression: branch detection + simple-fraction score_counter.

Locks the post-fix contract for `_score_row` and `score_counter`:
- hit count = number of (mz, logI) entries where float(logI) > 1.0
- score = hit_count / len(ionlist_mz), using the user-provided ion list length
  as denominator (NOT the deduped unique-mass count)

See handoff/codex_B24_scoreA_shadow_anchors_review.md §7-8.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
import numpy as np
import pandas as pd
import msp_insilicomarker_withGPT as marker


def test_score_counter_is_simple_fraction():
    """B-24-L2: score_counter is the simple fraction; no anchor gate."""
    # 5 matches out of 81 reference -> 5/81 regardless of which m/z they are
    hits = [100.0, 200.0, 300.0, 400.0, 500.0]
    ionlist = list(range(81))  # any 81 reference ions
    assert marker.score_counter(hits, ionlist) == pytest.approx(5 / 81)

    # zero matches -> 0
    assert marker.score_counter([], ionlist) == 0.0

    # score is bounded [0, 1]
    assert marker.score_counter(list(range(200)), ionlist) <= 1.0


def test_score_counter_ignores_anchor_membership():
    """B-24-L2: score does NOT depend on whether 204.087 / 366.140 / 512.197 are present.

    Pre-fix: hits without >=2 of these anchors returned 0.0. Post-fix: no such gate.
    """
    ionlist = list(range(81))
    # Five non-anchor m/z values
    no_anchor_hits = [100.0, 150.0, 250.0, 300.0, 400.0]
    # Five anchor m/z values (the deleted ANCHORS triplet + 2 extras)
    anchor_hits = [204.087, 366.140, 512.197, 600.0, 700.0]
    # Both should produce the SAME score post-fix (anchor membership is irrelevant)
    assert marker.score_counter(no_anchor_hits, ionlist) == pytest.approx(
        marker.score_counter(anchor_hits, ionlist)
    )


def test_anchor_helpers_removed():
    """B-24-L2: ANCHORS / ANCHOR_TOL / has_anchors are gone from the module."""
    assert not hasattr(marker, "ANCHORS")
    assert not hasattr(marker, "ANCHOR_TOL")
    assert not hasattr(marker, "has_anchors")


def test_score_row_filters_by_logI_regardless_of_length_match():
    """B-24-L1: branch detection uses element shape, not length equality.

    Pre-fix: when len(hits_out) != len(ionlist_mz) (e.g., findingions dedupes
    duplicate ion masses), execution fell into the 'matched-only' branch and
    counted miss placeholders [mz, 1] as hits. Post-fix: the 2-tuple shape
    triggers the logI > 1.0 filter regardless of length.
    """
    # Simulate findingions output: 72 entries (deduped), but only 3 real hits.
    # The other 69 are miss placeholders with logI = 1.0.
    fake_hits_out = [[100.0 + i, 1.0] for i in range(69)]  # 69 misses
    fake_hits_out += [[500.0, 2.5], [600.0, 3.1], [700.0, 2.0]]  # 3 hits

    # Simulate ionlist_mz of length 81 (raw, with 9 duplicate masses)
    ionlist_mz = np.array([100.0 + i for i in range(81)], dtype=float)

    matched_df = pd.DataFrame({
        "MS2scan_no": ["test_scan"],
        "peaklist": [[]],          # not used by mocked findingions
        "peakintensity": [[]],
    })
    ion_df = pd.DataFrame({"mass": ionlist_mz})

    # Patch findingions to return our fixture; call _score_row indirectly
    # via attach_ion_score_on_matched on a single-row DataFrame.
    original_findingions = marker.findingions
    marker.findingions = lambda row, df, ppm: fake_hits_out
    try:
        result = marker.attach_ion_score_on_matched(
            matched_df=matched_df, ion_df=ion_df, ppm_value=20.0,
            scan_col="MS2scan_no", ion_mass_col="mass",
        )
    finally:
        marker.findingions = original_findingions

    # B-24-L1 lock: hit count is 3 (real logI>1), NOT 72 (inflated by placeholders)
    assert int(result["ion hit count"].iloc[0]) == 3, \
        f"hit count should be 3 (logI>1 only), got {result['ion hit count'].iloc[0]}"

    # B-24 denominator contract: score = 3 / len(ionlist_mz) = 3 / 81
    # (NOT 3 / 72 -- the raw ion list is the scoring contract, not deduped)
    assert result["ion score"].iloc[0] == pytest.approx(3 / 81), \
        f"score should be 3/81={3/81:.6f}, got {result['ion score'].iloc[0]:.6f}"


def test_b24_denominator_uses_raw_ion_list_length_not_unique():
    """B-24 denominator contract (Henry-locked 2026-05-19): when the user provides
    81 ion rows with 72 unique masses, the denominator stays 81. The ion sheet
    is the scoring contract; deduplication is upstream's responsibility, not
    the scorer's.
    """
    # 81 entries, 9 duplicates
    ionlist_raw = [float(i) for i in range(72)] + [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    assert len(ionlist_raw) == 81
    assert len(set(ionlist_raw)) == 72  # confirm duplicates exist

    hits = [10.0, 20.0, 30.0]
    score = marker.score_counter(hits, ionlist_raw)
    # Score is 3/81, NOT 3/72
    assert score == pytest.approx(3 / 81)
    assert score != pytest.approx(3 / 72)
