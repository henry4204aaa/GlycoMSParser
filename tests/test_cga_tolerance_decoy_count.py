"""Synthetic tests for tools/validation/cga_tolerance_decoy_count.py (CGA-TOLERANCE, 2026-10-31).

Run (from the repository root):  python3 -m pytest tests/test_cga_tolerance_decoy_count.py -q
"""
import importlib.util
import json
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)                                      # tests -> repository root
TOOLS = os.path.join(REPO, "tools", "validation")
_spec = importlib.util.spec_from_file_location("decoy", os.path.join(TOOLS, "cga_tolerance_decoy_count.py"))
D = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(D)

LIB = [911.3361, 1073.3889, 2223.7980]
IONS = [204.0867 + 30.0 * i for i in range(14)]


def _write(tmp, scans, lib=LIB, ions=IONS):
    conv = tmp / "secret_input_name_conv.csv"
    with open(conv, "w") as f:
        f.write("MS2scan_no\tprotonatedmass\tpeaklist\tpeakintensity\n")
        for i, (m, peaks, ints) in enumerate(scans, 1):
            f.write(f"{i}\t{m}\t{tuple(peaks)}\t{tuple(ints)}\n")
    libp = tmp / "secret_lib.csv"
    libp.write_text("Hex,HexNAc,Mass\n" + "".join(f"1,1,{m}\n" for m in lib))
    ionp = tmp / "secret_ions.csv"
    ionp.write_text("ion name,mass\n" + "".join(f",{m}\n" for m in ions))
    return str(conv), str(libp), str(ionp)


def test_parse_tolerances():
    assert D.parse_tolerances("20ppm, 0.3Da,300 PPM") == [(20.0, "ppm"), (0.3, "Da"), (300.0, "ppm")]
    for bad in ("20", "0ppm", "x Da", "-1ppm"):
        with pytest.raises(ValueError):
            D.parse_tolerances(bad)


def test_precursor_counts_ppm_da_and_decoy():
    lib = np.array(sorted(LIB))
    # -0.15 Da at 911 (~165 ppm); +3.37 Da exactly hits the decoy of 1073.3889; 0 = no mass
    scans = [(911.3361 - 0.15, np.array([])), (1073.3889 + 3.37, np.array([])), (0.0, np.array([]))]
    r = D.run(scans, lib, None, [(20.0, "ppm"), (300.0, "ppm"), (0.2, "Da")], [])
    p = {x["precursor"]: x for x in r["precursor"]}
    assert r["n_scans_with_mass"] == 2
    assert (p["20 ppm"]["real_scans"], p["20 ppm"]["decoy_scans"]) == (0, 1)
    assert (p["300 ppm"]["real_scans"], p["300 ppm"]["decoy_scans"]) == (1, 1)
    assert (p["0.2 Da"]["real_scans"], p["0.2 Da"]["decoy_scans"]) == (1, 1)


def test_fragment_ppm_window_keeps_the_v1_10_asymmetry():
    ion = np.array([500.0])
    f = 10 / 1e6
    def hits(pk):
        return D.fragment_hits(D.make_scan(1.0, [pk], [1000.0])[1], ion, (10.0, "ppm"))[0]
    assert hits(500.0 * (1 - f) + 1e-9) == 1           # below the ion, inside 10 ppm (<=)
    assert hits(500.0 / (1 - f) - 1e-9) == 1           # above, just under 10 ppm -> hit
    assert hits(500.0 / (1 - f) + 1e-9) == 0           # above, just over -> no hit


def test_fragment_da_window_symmetric_and_repeats_counted_two_ways():
    ions = np.array([300.0, 400.0])
    spec = D.make_scan(1.0, [299.75, 300.0, 300.25, 300.375, 400.5], [1000.0] * 5)[1]
    per_ion = D.fragment_hits(spec, ions, (0.25, "Da"))
    assert per_ion.tolist() == [3, 0]
    assert D.score_a(spec, ions, (0.25, "Da")) == (3, 1)        # glycomsp vs unique (B-39)


def test_end_to_end_cli_counts_only_and_no_paths(tmp_path, capsys):
    shifted = sorted(m + 0.1 for m in IONS[:10])
    scans = [(911.3361 - 0.15, shifted, [1000.0] * 10),
             (1073.3889 * (1 + 5e-6), sorted(IONS[10:12]), [1000.0, 1000.0]),
             (1073.3889 + 3.37, sorted(m + 3.37 for m in IONS[:3]), [1000.0] * 3),   # decoy precursor + decoy ions
             (1500.0, [300.0], [0.5])]                                             # low intensity, no match
    conv, lib, ions = _write(tmp_path, scans)
    out_csv = tmp_path / "counts.csv"
    rc = D.main(["--alias", "pilot", "--converted", conv, "--library", lib, "--ions", ions,
                 "--precursor", "20ppm,300ppm", "--fragment", "10ppm,0.3Da", "--out", str(out_csv), "--json"])
    assert rc == 0
    res = json.loads(capsys.readouterr().out)
    p = {x["precursor"]: x for x in res["precursor"]}
    assert (p["20 ppm"]["real_scans"], p["20 ppm"]["decoy_scans"]) == (1, 1)
    assert (p["300 ppm"]["real_scans"], p["300 ppm"]["decoy_scans"]) == (2, 1)
    fr = {(x["precursor"], x["fragment"]): x for x in res["fragments"]}
    w = fr[("300 ppm", "0.3 Da")]["glycomsp"]
    assert w["real_prec+real_ions"] == 2                 # scan 1: 10/14 at 0.3 Da; scan 2: 2 exact ions = 0.143
    # at 10 ppm scan 1's +0.1 Da fragments miss, so only scan 2 passes
    assert fr[("20 ppm", "10 ppm")]["glycomsp"]["real_prec+real_ions"] == 1
    assert fr[("300 ppm", "10 ppm")]["glycomsp"]["real_prec+real_ions"] == 1
    assert fr[("300 ppm", "0.3 Da")]["glycomsp"]["decoy_prec+decoy_ions"] == 1
    text = out_csv.read_text()
    assert "secret" not in text and str(tmp_path) not in text and "pilot" in text
    rc = D.main(["--alias", "pilot", "--converted", conv, "--library", lib, "--precursor", "20ppm"])
    printed = capsys.readouterr().out
    assert rc == 0 and "alias: pilot" in printed and "secret" not in printed and str(tmp_path) not in printed


def test_min_hits_rule_and_path_free_errors(tmp_path, capsys):
    scans = [(911.3361, sorted(IONS[:2]), [1000.0, 1000.0])]
    conv, lib, ions = _write(tmp_path, scans)
    D.main(["--alias", "a", "--converted", conv, "--library", lib, "--ions", ions, "--min-hits", "3", "--json"])
    res = json.loads(capsys.readouterr().out)
    assert res["pass_rule"] == "hits >= 3"
    assert res["fragments"][0]["glycomsp"]["real_prec+real_ions"] == 0
    with pytest.raises(SystemExit) as e:
        D.main(["--alias", "a", "--converted", str(tmp_path / "missing_secret.csv"), "--library", lib])
    assert "secret" not in str(e.value) and str(tmp_path) not in str(e.value)
    with pytest.raises(SystemExit):
        D.main(["--alias", "bad/alias", "--converted", conv, "--library", lib])


# ---------- CP-C.1 (Codex round 1: CGAT-F7 path-free errors, CGAT-F8 exact equivalence) ----------

import sys  # noqa: E402
sys.path.insert(0, os.path.join(REPO, "src"))


def _production_hits(peaks, ints, ions, tol):
    import mspvalidator_merger as M
    v, u = tol
    out = M.peaks_ppm(list(peaks), list(ints), list(ions), v, tol_da=(v if u == "Da" else None))
    return sum(1 for _, val in out if float(val) > 1.0)


@pytest.mark.parametrize("tol", [(10.0, "ppm"), (300.0, "ppm"), (0.3, "Da"), (0.05, "Da")])
def test_glycomsp_count_equals_production_on_random_spectra(tol):
    pytest.importorskip("mspvalidator_merger")
    rng = np.random.default_rng(11)
    ions = np.round(np.sort(rng.uniform(200, 900, 20)), 4)
    for trial in range(300):
        n = int(rng.integers(0, 40))
        peaks = list(np.round(rng.uniform(190, 910, n), 2))
        if trial % 3 == 0 and n:                      # repeated m/z with a different intensity
            peaks += [peaks[0], peaks[-1]]
        ints = list(rng.choice([0.0, 0.5, 1.0, 1.5, 1000.0], len(peaks)))
        if trial % 2 == 0:
            peaks_in = sorted(peaks)
            ints_in = [i for _, i in sorted(zip(peaks, ints), key=lambda t: t[0])]
        else:                                         # unsorted
            peaks_in, ints_in = peaks, ints
        _, scan = D.make_scan(500.0, peaks_in, ints_in)
        g, _ = D.score_a(scan, ions, tol)
        with np.errstate(divide="ignore"):
            assert g == _production_hits(peaks_in, ints_in, ions, tol), (trial, tol)


def test_codex_duplicate_peak_case():
    pytest.importorskip("mspvalidator_merger")
    ions = np.array([300.0] + [500.0 + i for i in range(13)])
    _, scan = D.make_scan(1000.0, [300.0, 300.0], [0.5, 1000.0])
    assert isinstance(scan, D.Irregular)
    g, u = D.score_a(scan, ions, (10.0, "ppm"))
    with np.errstate(divide="ignore"):
        assert g == u == 0 == _production_hits([300.0, 300.0], [0.5, 1000.0], ions, (10.0, "ppm"))


def test_no_path_or_content_in_any_error(tmp_path, capsys):
    scans = [(911.3361, sorted(IONS[:2]), [1000.0, 1000.0])]
    conv, lib, ions = _write(tmp_path, scans)
    secret_dir = tmp_path / "secret_dir_does_not_exist"
    # --out into a missing folder
    with pytest.raises(SystemExit) as e:
        D.main(["--alias", "a", "--converted", conv, "--library", lib, "--out", str(secret_dir / "x.csv")])
    assert "secret" not in str(e.value) and str(tmp_path) not in str(e.value)
    # a malformed peak token carrying identifying text
    bad = tmp_path / "secret_bad.csv"
    bad.write_text("MS2scan_no\tprotonatedmass\tpeaklist\tpeakintensity\n1\t911.3\t(1.0, secret_token)\t(1.0, 2.0)\n")
    with pytest.raises(SystemExit) as e:
        D.main(["--alias", "a", "--converted", str(bad), "--library", lib])
    assert "secret" not in str(e.value)
    # argparse must not echo a stray argument (an unquoted file name with a space)
    with pytest.raises(SystemExit):
        D.main(["--alias", "a", "--converted", conv, "--library", lib, "secret name.csv"])
    err = capsys.readouterr().err
    assert "secret" not in err and "not echoed" in err
    # a bad tolerance is reported without echoing it
    with pytest.raises(SystemExit) as e:
        D.main(["--alias", "a", "--converted", conv, "--library", lib, "--precursor", "secret"])
    assert "secret" not in str(e.value)



# ---------- CP-C.2 (Codex round 2, CGAT-R2-F4): float boundaries and near-one intensities ----------

def test_codex_r2_boundary_cases():
    pytest.importorskip("mspvalidator_merger")
    ions = np.array([300.0] + [600.0 + i for i in range(13)])
    # 300 + 0.3 rounds to the float 300.3, but production computes 300.3 - 300 = 0.30000000000001137 > 0.3
    spec = D.make_scan(1.0, [300.3], [1000.0])[1]
    with np.errstate(divide="ignore"):
        assert D.score_a(spec, ions, (0.3, "Da"))[0] == _production_hits([300.3], [1000.0], ions, (0.3, "Da")) == 0
    # first intensity 1.0000000000000002: log10(I) + 1 rounds to 1.0, so production counts neither duplicate
    irr = D.make_scan(1.0, [300.0, 300.0], [1.0000000000000002, 1000.0])[1]
    with np.errstate(divide="ignore"):
        assert D.score_a(irr, ions, (10.0, "ppm"))[0] == \
            _production_hits([300.0, 300.0], [1.0000000000000002, 1000.0], ions, (10.0, "ppm")) == 0


@pytest.mark.parametrize("tol", [(10.0, "ppm"), (300.0, "ppm"), (20.0, "ppm"), (0.3, "Da"), (0.25, "Da"), (0.1, "Da")])
def test_edges_equal_production(tol):
    """Peaks placed exactly at the window edges as float operations produce them, and one float either side;
    intensities around 1. Both the ascending path and the irregular (literal loop) path must equal production."""
    pytest.importorskip("mspvalidator_merger")
    v, u = tol
    refs = [300.0, 300.3, 512.3456, 1234.5678, 2223.798]
    near_one = [1.0, 1.0000000000000002, 0.9999999999999999, 1.0000000001, 1000.0, 0.0, -5.0]
    k = 0
    for ref in refs:
        if u == "Da":
            edges = [ref - v, ref + v]
        else:
            f = v / 1e6
            edges = [ref * (1 - f), ref / (1 - f), ref - ref * f, ref + ref * f]
        cand = []
        for e in edges:
            cand += [np.nextafter(e, -np.inf), e, np.nextafter(e, np.inf)]
        cand = sorted(set(float(c) for c in cand))
        for start in range(len(cand)):
            peaks = cand[start:] + cand[:start]                     # rotations: ascending and unsorted
            ints = [near_one[(k + j) % len(near_one)] for j in range(len(peaks))]
            k += 1
            for pk_in, it_in in ((peaks, ints), (sorted(peaks), [i for _, i in sorted(zip(peaks, ints))])):
                _, scan = D.make_scan(1.0, pk_in, it_in)
                with np.errstate(divide="ignore", invalid="ignore"):
                    assert D.score_a(scan, np.array(refs), tol)[0] == _production_hits(pk_in, it_in, refs, tol), \
                        (tol, ref, pk_in, it_in)


# ---------- CP-C.3 (Codex round 3: CGAT-R3-F3 extreme windows, CGAT-R3-F4 B-39 shared-peak counter) ----------

_spec_b = importlib.util.spec_from_file_location("b39shared", os.path.join(TOOLS, "B39_shared_peaks.py"))
B = importlib.util.module_from_spec(_spec_b)
_spec_b.loader.exec_module(B)


def test_codex_r3_extreme_ppm_case():
    pytest.importorskip("mspvalidator_merger")
    ions = np.array([300.0])
    pk, it, tol = [1499999999.9568815], [1000.0], (999999.8, "ppm")
    _, spec = D.make_scan(1.0, pk, it)
    with np.errstate(divide="ignore"):
        assert D.score_a(spec, ions, tol)[0] == _production_hits(pk, it, ions, tol) == 1


@pytest.mark.parametrize("seed", range(4))
def test_out_of_domain_windows_equal_production(seed):
    pytest.importorskip("mspvalidator_merger")
    rng = np.random.default_rng(100 + seed)
    tols = [(float(rng.uniform(4e5, 999999.99)), "ppm"), (float(rng.uniform(5e5, 5e5 + 1)), "ppm"),
            (float(rng.uniform(900, 5000)), "Da"), (1e-12, "Da")]
    for tol in tols:
        ions = np.sort(rng.uniform(100, 5e6, 6))
        for _ in range(40):
            peaks = sorted(set(float(x) for x in rng.uniform(1, 2e7, int(rng.integers(0, 12)))))
            ints = list(rng.choice([0.5, 1.0000000000000002, 1000.0], len(peaks)))
            _, spec = D.make_scan(1.0, peaks, ints)
            with np.errstate(divide="ignore", invalid="ignore"):
                assert D.score_a(spec, ions, tol)[0] == _production_hits(peaks, ints, ions, tol), (tol, peaks)


def test_b39_shared_counter_matches_glycomsp_rules():
    sh = lambda pk, it, ions, tol: B.shared_in_scan(D.make_scan(1.0, pk, it)[1], np.array(ions), tol)
    # Codex: 300.3 is outside both 0.3 Da windows with GlycoMSP's subtraction (0.30000000000001137 > 0.3)
    assert sh([300.3], [1000.0], [300.0, 300.6], (0.3, "Da")) == 0
    # repeated m/z: the first occurrence's intensity decides for both positions
    assert sh([300.0, 300.0], [0.5, 1000.0], [299.9, 300.1], (0.2, "Da")) == 0
    assert sh([300.0, 300.0], [1000.0, 0.5], [299.9, 300.1], (0.2, "Da")) == 2
    # unsorted: GlycoMSP stops at the first out-of-window peak above the ion (300.3), so nothing is counted
    assert sh([300.3, 299.9, 300.05], [1000.0] * 3, [300.0, 300.1], (0.2, "Da")) == 0
    assert sh(sorted([300.3, 299.9, 300.05]), [1000.0] * 3, [300.0, 300.1], (0.2, "Da")) == 1
    # per-ion hit counts behind the shared counter equal production
    pytest.importorskip("mspvalidator_merger")
    for pk, it in (([300.3, 299.9, 300.05], [1000.0] * 3), ([300.0, 300.0], [1000.0, 0.5])):
        for ion in (299.9, 300.0, 300.1, 300.6):
            got = len(D.fragment_hit_indices(D.make_scan(1.0, pk, it)[1], [ion], (0.2, "Da"))[0])
            with np.errstate(divide="ignore"):
                assert got == _production_hits(pk, it, [ion], (0.2, "Da"))


# ---------------------------------------------------------------- Codex round 4 (CGAT-R4-F1): ion mass <= 0, ppm

def test_codex_r4_negative_ion_ppm_case():
    """Codex: ion -300, peak -500, 10 ppm. GlycoMSP's lower-side ppm is 1e6 * 200 / -300 < 0 <= 10 -> 1 hit."""
    pytest.importorskip("mspvalidator_merger")
    _, spec = D.make_scan(1.0, [-500.0], [1000.0])
    assert isinstance(spec, D.Spectrum)                        # a regular spectrum ...
    assert not D._fast_ok(spec, np.array([-300.0]), (10.0, "ppm"))   # ... but routed to the literal loop
    assert D.score_a(spec, np.array([-300.0]), (10.0, "ppm"))[0] == _production_hits([-500.0], [1000.0], [-300.0],
                                                                                       (10.0, "ppm")) == 1
    # positive ions keep the fast path even with zero / negative peaks; Da windows keep it for any sign
    _, spec2 = D.make_scan(1.0, [-0.5, 0.0, 300.0], [1000.0] * 3)
    assert D._fast_ok(spec2, np.array([300.0]), (10.0, "ppm"))
    assert D._fast_ok(spec, np.array([-300.0]), (0.3, "Da"))


def test_codex_r4_zero_ion_raises_like_glycomsp():
    """Ion 0 with a peak <= 0: GlycoMSP divides by zero; the helper raises the same error instead of returning 0."""
    pytest.importorskip("mspvalidator_merger")
    _, spec = D.make_scan(1.0, [-0.5], [1000.0])
    with pytest.raises(ZeroDivisionError):
        _production_hits([-0.5], [1000.0], [0.0], (10.0, "ppm"))
    with pytest.raises(ZeroDivisionError):
        D.fragment_hit_indices(spec, [0.0], (10.0, "ppm"))
    with pytest.raises(ZeroDivisionError):
        B.shared_in_scan(spec, np.array([0.0, 300.0]), (10.0, "ppm"))
    # ion 0 with positive peaks only: GlycoMSP's first peak above 0 is 1e6 ppm away -> stop, 0 hits; helper agrees
    _, spec_pos = D.make_scan(1.0, [100.0, 200.0], [1000.0, 1000.0])
    assert D.score_a(spec_pos, np.array([0.0]), (10.0, "ppm"))[0] == _production_hits([100.0, 200.0], [1000.0] * 2,
                                                                                      [0.0], (10.0, "ppm")) == 0


@pytest.mark.parametrize("seed", range(4))
def test_mixed_sign_masses_equal_production(seed):
    """Random ions and peaks of both signs (no exact zero ion), ppm and Da: helper == GlycoMSP, and the B-39 shared
    counter equals the count from GlycoMSP's literal per-ion positions."""
    pytest.importorskip("mspvalidator_merger")
    rng = np.random.default_rng(400 + seed)
    for tol in ((10.0, "ppm"), (300.0, "ppm"), (4e5, "ppm"), (0.3, "Da"), (50.0, "Da")):
        for _ in range(40):
            ions = rng.choice([-1, 1], 5) * rng.uniform(0.5, 2000, 5)
            peaks = sorted(set(float(x) for x in rng.uniform(-2000, 2000, int(rng.integers(0, 12)))))
            ints = list(rng.choice([0.5, 1000.0], len(peaks)))
            _, spec = D.make_scan(1.0, peaks, ints)
            with np.errstate(divide="ignore", invalid="ignore"):
                assert D.score_a(spec, ions, tol)[0] == _production_hits(peaks, ints, ions, tol), (tol, ions, peaks)
                lit = D._literal_hit_indices(peaks, ints, ions, tol)
            assert D.fragment_hit_indices(spec, ions, tol) == lit
            from collections import Counter
            c = Counter(k for h in lit for k in set(h))
            assert B.shared_in_scan(spec, ions, tol) == sum(1 for n in c.values() if n >= 2)
