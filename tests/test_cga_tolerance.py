"""CGA-TOLERANCE (2026-10-31): CGA precursor (T1) and Score A fragment (T2) tolerances in ppm or Da.

Defaults must reproduce v1.10 exactly; a Da window is symmetric (±value) and only applies when set.
Synthetic data only (no fixture), so the module always runs.

Run from repo root:  python3 -m pytest tests/test_cga_tolerance.py -q
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _v14_extract as X  # noqa: E402

sys.path.insert(0, X.SRC)
import msp_insilicomarker_withGPT as marker  # noqa: E402
import mspvalidator_merger as merger  # noqa: E402

COMP = ("Hex", "HexNAc", "NeuAc", "NeuGc", "KDN", "Fuc")


def _lib(masses):
    rows = [dict(zip(COMP, (i + 2, 2, 0, 0, 0, 0)), Mass=m) for i, m in enumerate(masses)]
    return pd.DataFrame(rows).sort_values("Mass").reset_index(drop=True)


# Native reduced [M-H]- (negative mode; H3N2 911.3361 from the ticket) and PerMe [M+Na]+-like (positive mode) masses
NEG_LIB = _lib([911.3361, 1073.3889, 1235.4417, 1438.5211, 2223.7980])
POS_LIB = _lib([1579.7835, 1783.8833, 2244.1080, 2605.2817])


def _brute(obs, lib, ppm=None, da=None):
    t = lib["Mass"].to_numpy(float)
    tol = abs(obs) * ppm / 1e6 if da is None else da
    return sorted(t[(t >= obs - tol) & (t <= obs + tol)].tolist())


def _found(obs, lib, **kw):
    return sorted(marker.find_compositions_for_mass(obs, lib, mass_col="Mass", **kw)["theoretical_mass"].tolist())


# ---------- find_compositions_for_mass ----------

@pytest.mark.parametrize("lib", [NEG_LIB, POS_LIB], ids=["neg", "pos"])
def test_default_is_the_v1_10_ppm_window(lib):
    rng = np.random.default_rng(7)
    for theo in lib["Mass"]:
        for off_ppm in rng.uniform(-40, 40, 25):
            obs = theo * (1 + off_ppm / 1e6)
            a = marker.find_compositions_for_mass(obs, lib, mass_col="Mass")
            b = marker.find_compositions_for_mass(obs, lib, mass_col="Mass", ppm=20.0, tol_da=None)
            pd.testing.assert_frame_equal(a, b)
            assert sorted(a["theoretical_mass"].tolist()) == _brute(obs, lib, ppm=20.0)


def test_ltq_like_offsets_need_the_wide_window():
    # z=1: -0.15 Da on 911.34 (~ -165 ppm); z=2: +0.40 Da on 2223.80 (~ +180 ppm)
    for theo, off in ((911.3361, -0.15), (2223.7980, +0.40)):
        obs = theo + off
        assert _found(obs, NEG_LIB) == []                                  # 20 ppm default: no match
        assert theo in _found(obs, NEG_LIB, ppm=300.0)                     # 300 ppm
        assert theo in _found(obs, NEG_LIB, tol_da=0.5)                    # 0.5 Da
    assert _found(911.3361 - 0.15, NEG_LIB, tol_da=0.3) == [911.3361]
    assert _found(2223.7980 + 0.40, NEG_LIB, tol_da=0.3) == []             # 0.40 Da > 0.3 Da window


@pytest.mark.parametrize("lib", [NEG_LIB, POS_LIB], ids=["neg", "pos"])
@pytest.mark.parametrize("da", [0.05, 0.3, 1.0])
def test_da_window_is_symmetric_and_matches_brute_force(lib, da):
    for theo in lib["Mass"]:
        for off in (-1.2 * da, -0.999 * da, -0.5 * da, 0.0, 0.5 * da, 0.999 * da, 1.2 * da):
            obs = theo + off
            assert _found(obs, lib, tol_da=da) == _brute(obs, lib, da=da)
            assert _found(obs, lib, ppm=1e9, tol_da=da) == _brute(obs, lib, da=da)   # ppm ignored with tol_da


def test_da_window_does_not_scale_with_mass():
    # 0.3 Da is ~330 ppm at 911 and ~135 ppm at 2224; 300 ppm is 0.27 / 0.67 Da
    assert _found(911.3361 + 0.29, NEG_LIB, tol_da=0.3) == [911.3361]
    assert _found(911.3361 + 0.29, NEG_LIB, ppm=300.0) == []
    assert _found(2223.7980 + 0.5, NEG_LIB, tol_da=0.3) == []
    assert _found(2223.7980 + 0.5, NEG_LIB, ppm=300.0) == [2223.7980]


def test_columns_and_errors_unchanged_with_da():
    obs = 1235.4417 + 0.2
    a = marker.find_compositions_for_mass(obs, NEG_LIB, mass_col="Mass", tol_da=0.3)
    assert list(a.columns[:5]) == ["observed_mass", "theoretical_mass", "delta_Da", "ppm_error", "comp_str"]
    assert a["delta_Da"].iloc[0] == pytest.approx(0.2)
    assert a["ppm_error"].iloc[0] == pytest.approx(0.2 / 1235.4417 * 1e6)
    empty = marker.find_compositions_for_mass(5000.0, NEG_LIB, mass_col="Mass", tol_da=0.3)
    assert empty.empty and "ppm_error" in empty.columns


# ---------- v14 module-level helpers (extracted, no Tk root) ----------

@pytest.fixture(scope="module")
def v14():
    try:
        return X.load_v14()
    except X.AnchorMissing as e:
        pytest.skip(f"v14 anchor moved (coverage LOST, update this test): {e}")


def test_parse_tolerance(v14):
    p = v14["_parse_tolerance"]
    assert p("20", "ppm") == {"value": 20.0, "unit": "ppm"}
    assert p(" 0.3 ", "Da") == {"value": 0.3, "unit": "Da"}
    for bad in ("", "abc", "0", "-1", "nan", "inf"):
        with pytest.raises(ValueError):
            p(bad, "ppm")
    with pytest.raises(ValueError):
        p("1", "mDa")


def test_tolerance_kwargs_default_is_empty(v14):
    k = v14["_cga_tolerance_kwargs"]
    assert k(None) == {} and k({}) == {}
    assert k({"precursor": {"value": 20.0, "unit": "ppm"}}) == {"ppm_value": 20.0}
    assert k({"precursor": {"value": 300.0, "unit": "ppm"}}) == {"ppm_value": 300.0}
    assert k({"precursor": {"value": 0.3, "unit": "Da"}}) == {"tol_da": 0.3}
    assert k({"fragment": {"value": 10.0, "unit": "ppm"}}) == {"ion_ppm": 10.0}
    assert k({"precursor": {"value": 300.0, "unit": "ppm"}, "fragment": {"value": 0.3, "unit": "Da"}}) == \
        {"ppm_value": 300.0, "ion_tol_da": 0.3}
    assert v14["CGA_FRAGMENT_TOL_DEFAULT"] == {"value": 10.0, "unit": "ppm"}
    assert v14["CGA_PRECURSOR_TOL_DEFAULT"] == {"value": 20.0, "unit": "ppm"}
    assert v14["_tolerance_label"]({"value": 300.0, "unit": "ppm"}) == "300 ppm"
    assert v14["_tolerance_label"]({"value": 0.3, "unit": "Da"}) == "0.3 Da"


# ---------- run_pseudolabeling end to end (extracted closure, synthetic LTQ-like input) ----------

# (scan, charge, observed precursor mass): z=1 -0.15 Da, z=2 +0.40 Da (LTQ isolation m/z), one +5 ppm control,
# one unrelated mass. Library: native reduced [M-H]- masses, negative mode.
LTQ_SCANS = ((11, 1, 911.3361 - 0.15), (12, 2, 2223.7980 + 0.40), (13, 1, 1073.3889 * (1 + 5e-6)), (14, 1, 1500.0))
LTQ_LIB = ((3, 2, 0, 911.3361), (4, 2, 0, 1073.3889), (5, 4, 1, 2223.7980))


def _write_ltq_inputs(d):
    conv = os.path.join(d, "ms2_ltq_like.csv")
    with open(conv, "w", encoding="utf-8") as f:
        f.write("entry_no\tMS1scan_no\tMS1_isolationmass\tMS1_monoisolationmass\tchargeState\tprotonatedmass"
                "\tMS2scan_no\tpeaklist\tpeakintensity\n")
        for i, (scan, z, m) in enumerate(LTQ_SCANS, 1):
            f.write(f"{i}\t{scan - 1}\t{(m + (z - 1) * 1.007276) / z}\t0\t{z}\t{m}\t{scan}"
                    f"\t(300.1, 400.2)\t(1000.0, 2000.0)\n")
    lib = os.path.join(d, "lib.csv")
    with open(lib, "w", encoding="utf-8") as f:
        f.write("Hex,HexNAc,NeuAc,NeuGc,KDN,Fuc,HexA,Mode,SO3,PO3H,Mass\n")
        for h, n, fu, m in LTQ_LIB:
            f.write(f"{h},{n},0,0,0,{fu},0,-,0,0,{m}\n")
    return conv, lib


@pytest.fixture(scope="module")
def cga_ns():
    try:
        return X.load_v14(extract=("run_pseudolabeling",))
    except X.AnchorMissing as e:
        pytest.skip(f"v14 anchor moved (coverage LOST, update this test): {e}")


def _run_cga(cga_ns, tmp_path, **kw):
    runlog = []
    cga_ns.update(X.cga_stubs(runlog))
    conv, lib = _write_ltq_inputs(str(tmp_path))
    out = cga_ns["run_pseudolabeling"]("S", {"csv": conv, "insilico_csv": lib}, {}, parent=None, **kw)
    assert out and os.path.exists(out)
    tsv = pd.read_csv(out, sep="\t")
    matched = sorted(tsv.loc[tsv["composition"].notna(), "MS2scan_no"].unique().tolist())
    return tsv, matched, runlog[-1], cga_ns["logger"].lines


def test_cga_defaults_keep_the_v1_10_result(cga_ns, tmp_path):
    tsv, matched, entry, log = _run_cga(cga_ns, tmp_path)
    assert matched == [13]                                   # only the +5 ppm control survives 20 ppm
    assert len(tsv) == len(LTQ_SCANS)                        # one row per scan (left merge)
    p = entry["params"]
    assert p["ppm_value"] == 20.0 and p["ion_ppm"] == 10.0
    assert p["precursor_tolerance"] == {"value": 20.0, "unit": "ppm"} and p["tolerance_non_default"] is False
    assert any("precursor match tolerance: 20 ppm (default)" in l for l in log)


@pytest.mark.parametrize("kw, expect", [
    ({"ppm_value": 300.0}, [11, 12, 13]),
    ({"tol_da": 0.3}, [11, 13]),
    ({"tol_da": 0.5}, [11, 12, 13]),
], ids=["300ppm", "0.3Da", "0.5Da"])
def test_cga_wide_windows_match_ltq_offsets(cga_ns, tmp_path, kw, expect):
    tsv, matched, entry, log = _run_cga(cga_ns, tmp_path, **kw)
    assert matched == expect
    p = entry["params"]
    assert p["tolerance_non_default"] is True
    if "tol_da" in kw:
        assert p["ppm_value"] is None and p["precursor_tolerance"] == {"value": kw["tol_da"], "unit": "Da"}
    else:
        assert p["ppm_value"] == 300.0
    assert any("NON-DEFAULT" in l for l in log)
    # ppm_error stays in ppm whatever the window unit
    row = tsv[tsv["MS2scan_no"] == 11].iloc[0]
    assert row["ppm_error"] == pytest.approx(-0.15 / 911.3361 * 1e6, rel=1e-6)


# ---------- T2: fragment window (merger peaks_ppm / findingions, Score A scorer, v14 fallback) ----------

def _hit_refs(nhitlist):
    return sorted(float(r) for r, v in nhitlist if float(v) > 1.0)


def test_peaks_ppm_default_path_unchanged():
    peaks = [100.0, 200.001, 300.0, 300.002, 450.0]
    ints = [10.0, 20.0, 30.0, 40.0, 50.0]
    refs = [200.0, 300.0, 450.004, 600.0]
    a = merger.peaks_ppm(peaks, ints, refs, 10)
    assert merger.peaks_ppm(peaks, ints, refs, 10, tol_da=None) == a
    assert _hit_refs(a) == [200.0, 300.0, 300.0, 450.004]   # v1.10 behaviour: one hit per peak in the window


def test_peaks_da_window_symmetric_inclusive():
    refs = [300.0]
    # binary-exact offsets so the boundary itself is tested (inclusive on both sides)
    for off, hit in ((-0.375, False), (-0.25, True), (-0.125, True), (0.0, True), (0.125, True), (0.25, True),
                     (0.375, False)):
        out = merger.peaks_ppm([300.0 + off], [5.0], refs, 10, tol_da=0.25)
        assert (_hit_refs(out) == [300.0]) is hit, off
    # setppm ignored with tol_da
    assert _hit_refs(merger.peaks_ppm([300.1], [5.0], refs, 1e9, tol_da=0.05)) == []


def test_peaks_da_same_output_shape_as_ppm():
    peaks = [199.95, 200.1, 350.0]
    ints = [100.0, 200.0, 300.0]
    refs = [200.0, 400.0, 400.0]                     # a duplicated ref is zero-filled once, as in the ppm path
    da = merger.peaks_ppm(peaks, ints, refs, 10, tol_da=0.2)
    pp = merger.peaks_ppm(peaks, ints, refs, 2000)   # 2000 ppm at 200 = 0.4 Da: same hits
    assert da == pp
    assert [r for r, v in da if v == 1] == [400.0]


def _frag_row(peaks):
    return pd.DataFrame({"MS2scan_no": [1], "peaklist": [str(tuple(peaks))],
                         "peakintensity": [str(tuple(1000.0 for _ in peaks))]})


ION_LIST14 = [204.0867 + 30.0 * i for i in range(14)]


def test_score_a_da_window():
    peaks = sorted([m + 0.1 for m in ION_LIST14[:10]] + ION_LIST14[10:12])     # 10 shifted +0.1 Da, 2 exact
    row = _frag_row(peaks)
    ions = pd.DataFrame({"mass": ION_LIST14})
    d = marker.attach_ion_score_on_matched(row, ions, ppm_value=10.0)
    assert d["ion hit count"].iloc[0] == 2 and d["ion score"].iloc[0] == pytest.approx(2 / 14)
    d0 = marker.attach_ion_score_on_matched(row, ions, ppm_value=10.0, tol_da=None)
    pd.testing.assert_frame_equal(d, d0)
    w = marker.attach_ion_score_on_matched(row, ions, ppm_value=10.0, tol_da=0.3)
    assert w["ion hit count"].iloc[0] == 12 and w["ion score"].iloc[0] == pytest.approx(12 / 14)
    n = marker.attach_ion_score_on_matched(row, ions, ppm_value=10.0, tol_da=0.05)
    assert n["ion hit count"].iloc[0] == 2


def test_fallback_scorer_follows_the_same_window(v14):
    peaks = sorted([m + 0.1 for m in ION_LIST14[:10]] + ION_LIST14[10:12])
    row = _frag_row(peaks)
    ions = pd.DataFrame({"mass": ION_LIST14})
    fb = v14["_fallback_simple_ion_scoring"]
    assert fb(row, ions, 10.0)["ion hit count"].iloc[0] == 2
    assert fb(row, ions, 10.0, tol_da=None)["ion hit count"].iloc[0] == 2
    assert fb(row, ions, 10.0, tol_da=0.3)["ion hit count"].iloc[0] == 12


def _write_ion_list(d):
    p = os.path.join(d, "ions14.csv")
    with open(p, "w", encoding="utf-8") as f:
        f.write("ion name,ion type,mass,notes\n")
        for m in ION_LIST14:
            f.write(f",,{m},\n")
    return p


def _run_cga_ions(cga_ns, tmp_path, **kw):
    """LTQ-like scans whose spectra hold the 14 ions: 10 at +0.1 Da, 2 exact."""
    runlog = []
    cga_ns.update(X.cga_stubs(runlog))
    conv, lib = _write_ltq_inputs(str(tmp_path))
    peaks = sorted([m + 0.1 for m in ION_LIST14[:10]] + ION_LIST14[10:12])
    df = pd.read_csv(conv, sep="\t")
    df["peaklist"] = str(tuple(peaks))
    df["peakintensity"] = str(tuple(1000.0 for _ in peaks))
    df.to_csv(conv, sep="\t", index=False)
    files = {"csv": conv, "insilico_csv": lib, "ionlist_path": _write_ion_list(str(tmp_path))}
    out = cga_ns["run_pseudolabeling"]("S", files, {}, parent=None, **kw)
    tsv = pd.read_csv(out, sep="\t")
    return tsv[tsv["composition"].notna()], runlog[-1], cga_ns["logger"].lines


def test_cga_default_fragment_window(cga_ns, tmp_path):
    m, entry, log = _run_cga_ions(cga_ns, tmp_path)
    assert m["MS2scan_no"].tolist() == [13]
    assert m["ion hit count"].tolist() == [2]
    p = entry["params"]
    assert p["ion_ppm"] == 10.0 and p["fragment_tolerance"] == {"value": 10.0, "unit": "ppm"}
    assert p["tolerance_non_default"] is False
    assert any("Score A fragment match tolerance: 10 ppm (default)" in l for l in log)


def test_cga_ltq_settings_300ppm_and_0_3_da(cga_ns, tmp_path):
    m, entry, log = _run_cga_ions(cga_ns, tmp_path, ppm_value=300.0, ion_tol_da=0.3)
    assert sorted(m["MS2scan_no"].tolist()) == [11, 12, 13]
    assert m["ion hit count"].tolist() == [12, 12, 12]
    assert m["ion score"].tolist() == pytest.approx([12 / 14] * 3)
    p = entry["params"]
    assert p["ion_ppm"] is None and p["fragment_tolerance"] == {"value": 0.3, "unit": "Da"}
    assert p["tolerance_non_default"] is True
    assert any("Score A fragment match tolerance: 0.3 Da (NON-DEFAULT" in l for l in log)


def test_cga_fragment_only_non_default_is_marked(cga_ns, tmp_path):
    m, entry, _ = _run_cga_ions(cga_ns, tmp_path, ion_ppm=500.0)
    assert m["ion hit count"].tolist() == [12]          # 0.1 Da at 204-594 m/z is 170-490 ppm
    assert entry["params"]["tolerance_non_default"] is True


# ---------- pinned baseline (CP-C.1, Codex round 1 Q10) ----------
# sha256 of the CGA TSV (CRLF normalised to LF) that the BASE build 7910068 (v1.10 CGA chain, v14 1.13) writes for the
# synthetic inputs above at default settings, computed 2026-10-31 by running the base run_pseudolabeling and base
# insilicomarker through tests/_v14_extract.py. A coordinated change of branch code and test can no longer pass
# unnoticed. If a pandas/numpy upgrade changes float text, re-pin from 7910068 the same way, never from the branch.
PINNED_BASE_7910068 = {
    "precursor_only": "c0fa89a00d005651b4a0cc057f9d9e8c3ea86814660ca3827145a810d711709b",
    "with_14_ions": "623a1267e245c5e685c3a226ef7e12a497747bf3a40de5b4dfde3e3e42e31981",
}


def _tsv_digest(path):
    import hashlib
    with open(path, "rb") as f:
        return hashlib.sha256(f.read().replace(b"\r\n", b"\n")).hexdigest()


def test_default_output_equals_pinned_base(cga_ns, tmp_path):
    cga_ns.update(X.cga_stubs([]))
    d1 = tmp_path / "p"; d1.mkdir()
    conv, lib = _write_ltq_inputs(str(d1))
    out = cga_ns["run_pseudolabeling"]("S", {"csv": conv, "insilico_csv": lib}, {}, parent=None)
    assert _tsv_digest(out) == PINNED_BASE_7910068["precursor_only"]
    d2 = tmp_path / "i"; d2.mkdir()
    conv, lib = _write_ltq_inputs(str(d2))
    peaks = sorted([m + 0.1 for m in ION_LIST14[:10]] + ION_LIST14[10:12])
    df = pd.read_csv(conv, sep="\t")
    df["peaklist"] = str(tuple(peaks))
    df["peakintensity"] = str(tuple(1000.0 for _ in peaks))
    df.to_csv(conv, sep="\t", index=False)
    cga_ns.update(X.cga_stubs([]))
    out = cga_ns["run_pseudolabeling"]("S", {"csv": conv, "insilico_csv": lib, "ionlist_path": _write_ion_list(str(d2))},
                                       {}, parent=None)
    assert _tsv_digest(out) == PINNED_BASE_7910068["with_14_ions"]


# ---------------------------------------------------------------- CP-D.2 CGA Setup layout (real Tk; skips without a display)
def test_cga_setup_buttons_stay_visible_when_window_is_short(v14):
    """2026-11-01 (Henry, Lab PC): the Mass tolerances frame pushed the window's requested height past its fixed
    800 px, and pack clips the last-packed widgets first, so Start/Close vanished unless the window was maximised.
    The bottom stack is now anchored to the bottom, so only the flags panel is clipped. The visual order is unchanged."""
    tk = pytest.importorskip("tkinter")
    try:
        root = tk.Tk()
    except tk.TclError as e:
        pytest.skip(f"no display for Tk: {e}")
    root.withdraw()
    try:
        import compnewv4
        w = v14["CGASetupWindow"](root, meta_prefill={"Glycan Type": "N"}, default_flags=compnewv4.NG_flags)
        w.minsize(1, 1)
        w.geometry("900x600+0+0")
        w.update()
        kids = w.pack_slaves()
        title = lambda c: c.cget("text") if "text" in c.keys() else c.winfo_children()[0].cget("text")
        by = {title(c): c for c in kids}
        flags = by["In-Silico Generation Flags"]
        for name in ("Files selected for CGA analysis", "Mass tolerances (CGA matching)",
                     "Load CGA settings", "Start CGA analysis"):
            c = by[name]
            assert c.winfo_height() == c.winfo_reqheight(), name     # fully shown
        assert flags.winfo_height() < flags.winfo_reqheight()         # the flags panel takes the clipping
        top_to_bottom = [title(c) for c in sorted(kids, key=lambda c: c.winfo_y())]
        assert top_to_bottom == ["Metadata (can override here)", "In-Silico Generation Flags",
                                 "Files selected for CGA analysis", "Mass tolerances (CGA matching)",
                                 "Load CGA settings", "Start CGA analysis"]
        w.destroy()
    finally:
        root.destroy()
