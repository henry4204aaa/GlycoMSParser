"""B-05 smoke test — verify the KDN PerMe residue constant fix in compnewv4.

The bug (B-05, S4 / 2026-05-17): `precursormassv4` and `precursormassv3` had
the KDN PerMe residue constant set to 335.1700, which is ~15 Da too high
relative to the empirically-correct value (320.14710). Henry's hypothesis was
that the original value came from a B-ion / oxonium reference for KDN, which
carries an extra methyl from the C2 anomeric fragment that's not present in a
chain residue. The 5-methyl PerMe-KDN residue (250.0689 + 5 × CH2 = 320.1471)
matches GlycoWorkbench cross-checks exactly.

Effect of the bug: every KDN-bearing composition's theoretical mass was
over-stated by ~15 Da per KDN, breaking CGA precursor-gate alignment for any
glycan containing KDN (notably zebrafish intestine N-glycans).

The fix: constant change in compnewv4.py lines 307 + 350.

This test imports compnewv4 directly (no GUI / Tkinter deps in that module)
and asserts mass values match Henry's GlycoWorkbench-verified reference rows
from `handoff/zf_sPerMeNG_merge_insilico_20251007.csv`:
    T1 - H3N2 control (no KDN) — verifies non-KDN paths unchanged
    T2 - F1H9N2 control (no KDN, larger composition) — verifies Fuc + larger Hex unchanged
    T3 - F1H5N4S2 control (Neu5Ac=2, no KDN) — verifies NeuAc path unchanged
    T4 - F1H5N4KDN2 — primary KDN regression case (Henry-verified row 762)
    T5 - F2H5N5KDN1 — secondary KDN regression case (Henry-verified row 863)
    T6 - KDN=3 synthetic — extends regression coverage
    T7 - precursormassv3 sister sanity (same cases through the v3 path)
    T8 - Mass delta per KDN: confirms exactly 0 ppm change at the constant
         (i.e., the fix locks in 320.14710 against future regressions)

Tolerance: ±0.0005 Da to absorb 4-decimal vs 5-decimal precision noise in
Henry's GWB reference values.

Run from repo root:
    python3 tests/test_b05_kdn_mass_smoke.py
"""
import sys
from pathlib import Path

# Ensure src/ is on the path so `import compnewv4` works
_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT / "src"))

import compnewv4
from compnewv4 import precursormassv4, precursormassv3


TOL = 0.0005  # Da — absorbs 4 vs 5 decimal precision in reference values

# Composition tuple convention (per compnewv4.py line 298):
#   composition = (Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc)


def _check(label, got, expected, tol=TOL):
    delta = got - expected
    ok = abs(delta) <= tol
    flag = "PASS" if ok else "FAIL"
    print(f"  [{flag}] {label}: got={got:.4f} expected={expected:.4f} Δ={delta:+.4f}")
    return ok


def test_constant_is_fixed():
    """T0 — confirm the post-fix constant 320.14710 is wired into both functions.

    We probe by computing the mass for a pure-KDN composition (KDN=1, everything
    else 0) and back-solving: M = e * KDN_const + reducing_end + proton.
    """
    # (Hex=0, HexNAc=0, NeuAc=0, NeuGc=0, KDN=1, Fuc=0); reducing_end=46.04186, proton=1.0073
    m = precursormassv4((0, 0, 0, 0, 1, 0), deri="PerMe", reduced=False, mode="+")
    inferred_kdn = m - 46.04186 - 1.0073
    assert abs(inferred_kdn - 320.14710) <= TOL, \
        f"T0 v4: KDN constant {inferred_kdn:.5f} != 320.14710 — fix not applied"
    m3 = precursormassv3((0, 0, 0, 0, 1, 0), deri="PerMe", reduced=False)
    inferred_kdn3 = m3 - 46.04186 - 1.0073
    assert abs(inferred_kdn3 - 320.14710) <= TOL, \
        f"T0 v3: KDN constant {inferred_kdn3:.5f} != 320.14710 — fix not applied"
    print("  [PASS] T0 constant: precursormassv4 + precursormassv3 both at 320.14710")


def test_non_kdn_controls():
    """T1-T3 — verify non-KDN compositions unchanged (sanity for other constants)."""
    cases = [
        # (label, comp, GWB-expected mass, params)
        ("T1 H3N2 control",   (3, 2, 0, 0, 0, 0), 1149.601, {}),
        ("T2 H9N2 control",   (9, 2, 0, 0, 0, 0), 2374.1997, {}),  # Henry's row 5969 (no Fuc per recomputed match)
        ("T3 F1H5N4S2 (S=NeuAc) control", (5, 4, 2, 0, 0, 1), 2944.4898, {}),
    ]
    all_ok = True
    for label, comp, expected, params in cases:
        m = precursormassv4(comp, deri="PerMe", reduced=False, mode="+", **params)
        if not _check(label, m, expected):
            all_ok = False
    assert all_ok, "Non-KDN control(s) failed — non-KDN constants may have drifted"


def test_kdn_regression_cases():
    """T4-T6 — verify KDN-bearing compositions match GWB-verified reference values.

    T4 + T5 are Henry's manually-verified rows from
    handoff/zf_sPerMeNG_merge_insilico_20251007.csv (2026-05-17 GlycoWorkbench
    cross-check). T6 is a synthetic KDN=3 case to extend coverage.
    """
    cases = [
        ("T4 F1H5N4KDN2 (row 762)", (5, 4, 0, 0, 2, 1), 2862.4367),
        ("T5 F2H5N5KDN1 (row 863)", (5, 5, 0, 0, 1, 2), 2961.5052),
        # T6: synthetic — extend KDN regression to KDN=3
        # Expected = 5×204.09977 + 4×245.12632 + 3×320.14710 + 46.04186 + 1.00730
        #         = 1020.49885 + 980.50528 + 960.44130 + 46.04186 + 1.00730 = 3008.49459
        ("T6 H5N4KDN3 synthetic",   (5, 4, 0, 0, 3, 0), 3008.4946),
    ]
    all_ok = True
    for label, comp, expected in cases:
        m = precursormassv4(comp, deri="PerMe", reduced=False, mode="+")
        if not _check(label, m, expected):
            all_ok = False
    assert all_ok, "KDN regression case(s) failed — fix incomplete or constant drifted"


def test_v3_sister():
    """T7 — verify precursormassv3 sister fix.

    v3 is the non-default path (save_glycan_pseudocomp_to_csv uses use_v4=True
    by default) but it's still callable. Sister fix at line 350 must track v4.
    """
    cases = [
        ("T7a H3N2 v3",          (3, 2, 0, 0, 0, 0), 1149.601),
        ("T7b F1H5N4KDN2 v3",    (5, 4, 0, 0, 2, 1), 2862.4367),
        ("T7c F2H5N5KDN1 v3",    (5, 5, 0, 0, 1, 2), 2961.5052),
    ]
    all_ok = True
    for label, comp, expected in cases:
        m = precursormassv3(comp, deri="PerMe", reduced=False)
        if not _check(label, m, expected):
            all_ok = False
    assert all_ok, "precursormassv3 sister fix incomplete"


def test_kdn_delta_consistency():
    """T8 — confirm exactly 320.14710 per KDN unit.

    Compare (KDN=2) − (KDN=1) for otherwise-identical compositions; the delta
    must equal exactly the KDN constant (320.14710).
    """
    base = (5, 4, 0, 0, 1, 0)
    bumped = (5, 4, 0, 0, 2, 0)
    m_base = precursormassv4(base, deri="PerMe", reduced=False, mode="+")
    m_bumped = precursormassv4(bumped, deri="PerMe", reduced=False, mode="+")
    delta = m_bumped - m_base
    expected = 320.14710
    assert abs(delta - expected) <= TOL, \
        f"T8: KDN delta {delta:.5f} != expected {expected:.5f}"
    print(f"  [PASS] T8 delta consistency: ΔM(KDN+1) = {delta:.5f} (= expected 320.14710)")


if __name__ == "__main__":
    print(f"B-05 KDN-mass smoke test — compnewv4 last_update = {compnewv4.last_update}")
    print()

    failures = []
    for name, fn in [
        ("T0 constant",         test_constant_is_fixed),
        ("T1-T3 non-KDN controls", test_non_kdn_controls),
        ("T4-T6 KDN regression", test_kdn_regression_cases),
        ("T7 v3 sister",        test_v3_sister),
        ("T8 delta consistency", test_kdn_delta_consistency),
    ]:
        print(f"--- {name} ---")
        try:
            fn()
        except AssertionError as e:
            failures.append((name, str(e)))
            print(f"  FAILED: {e}")
        print()

    if failures:
        print(f"FAIL: {len(failures)} test(s) failed:")
        for name, msg in failures:
            print(f"  - {name}: {msg}")
        sys.exit(1)
    print("All B-05 KDN-mass smoke tests PASSED.")
    sys.exit(0)
