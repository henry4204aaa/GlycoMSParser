"""CGA-TOLERANCE CP-C (2026-10-31): the tolerances actually used are recorded in the CGA method JSON
`parameters.cga` (tree key `_cga_parameters`), carried to training_run.json, and prefill the prediction fields
only when they differ from the v1.10 defaults. Headless (no Tk root), synthetic files only.

Run from repo root:  python3 -m pytest tests/test_cga_tolerance_record.py -q
"""
import datetime as _dt
import json
import os
import sys
import tempfile
import typing

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _v14_extract as X  # noqa: E402

PPM20 = {"value": 20.0, "unit": "ppm"}
PPM10 = {"value": 10.0, "unit": "ppm"}
WIDE_P = {"value": 300.0, "unit": "ppm"}
WIDE_F = {"value": 0.3, "unit": "Da"}
TH_DEFAULT = {"min_ion_score": 0.07, "max_abs_ppm": 20.0, "topn": 1, "ion_ppm": 10.0}
NEG_DEFAULT = {"enable": True, "max_ratio": 3.0, "min_hits": 3, "ion_ppm": 10.0}


@pytest.fixture(scope="module")
def ns():
    try:
        n = X.load_v14(extract=("build_method_v1_from_tree",))
    except X.AnchorMissing as e:
        pytest.skip(f"v14 anchor moved (coverage LOST, update this test): {e}")
    n.setdefault("experiment_projects", {})
    n["_tmp"] = tempfile.mkdtemp()
    return n


@pytest.fixture(autouse=True)
def _clean_registry(ns):
    ns["_CGA_TRAINABLE_REGISTRY"].clear()
    yield
    ns["_CGA_TRAINABLE_REGISTRY"].clear()


def _sample(ns):
    d = tempfile.mkdtemp(dir=ns["_tmp"])
    files = {}
    for key, name in (("csv", "s.csv"), ("metadata", "s.json"), ("ionlist_path", "ion.xlsx"),
                      ("insilico_csv", "ins.csv"), ("excel", "mas.xlsx"), ("pseudolabel_csv", "s_CGA.tsv")):
        p = os.path.join(d, name)
        open(p, "w").write("x")
        files[key] = p
    ns["experiment_projects"].clear()
    ns["experiment_projects"]["E"] = {"samples": {"S": files}}
    return d, files


def test_no_cga_run_keeps_empty_parameters(ns):
    _, files = _sample(ns)
    v1 = ns["build_method_v1_from_tree"]("E", "S", force_family="CGA")
    assert v1["parameters"] == {"mas": {}, "cga": {}, "scoring": {}, "ml": {}}


def test_run_and_trainable_recorded_in_cga_method_only(ns):
    d, files = _sample(ns)
    ns["_cga_record_run"](files, files["pseudolabel_csv"], PPM20, PPM10)
    v1 = ns["build_method_v1_from_tree"]("E", "S", force_family="CGA")
    run = v1["parameters"]["cga"]["run"]
    assert run["precursor_tolerance"] == PPM20 and run["fragment_tolerance"] == PPM10
    assert run["tolerance_non_default"] is False and run["tsv"] == os.path.normpath(files["pseudolabel_csv"])
    assert ns["build_method_v1_from_tree"]("E", "S", force_family="MAS")["parameters"]["cga"] == {}
    # trainable from the same TSV carries the CGA run tolerances
    files["_method_v1_cache_CGA"] = v1
    tcsv = os.path.join(d, "t.csv"); open(tcsv, "w").write("trainable")
    tr = ns["_cga_record_trainable"](files, tcsv, files["pseudolabel_csv"], TH_DEFAULT, NEG_DEFAULT, 10.0)
    assert tr["csv_sha256"] == ns["file_sha256"](tcsv)
    assert tr["cga_run"]["precursor_tolerance"] == PPM20 and tr["tolerance_non_default"] is False
    assert tr["precursor_gate_ppm"] == 10.0 and tr["max_abs_ppm"] == 20.0 and tr["feature_ion_ppm"] == 10.0
    assert tr["min_ion_score"] == 0.07 and tr["negatives"] == {"enabled": True, "ion_ppm": 10.0}
    # the cached method (what the autosave writes verbatim) follows
    assert files["_method_v1_cache_CGA"]["parameters"]["cga"]["trainable"] == tr


def test_trainable_from_another_tsv_has_unknown_run(ns):
    d, files = _sample(ns)
    ns["_cga_record_run"](files, files["pseudolabel_csv"], WIDE_P, WIDE_F)
    other = os.path.join(d, "older_CGA.tsv")
    tr = ns["_cga_record_trainable"](files, os.path.join(d, "t.csv"), other, TH_DEFAULT, NEG_DEFAULT, 10.0)
    assert tr["cga_run"] is None and tr["tolerance_non_default"] is False
    tr2 = ns["_cga_record_trainable"](files, os.path.join(d, "t.csv"), files["pseudolabel_csv"],
                                      TH_DEFAULT, NEG_DEFAULT, 10.0)
    assert tr2["cga_run"]["tolerance_non_default"] is True and tr2["tolerance_non_default"] is True


def test_score_b_output_path_followed(ns):
    d, files = _sample(ns)
    ns["_cga_record_run"](files, files["pseudolabel_csv"], PPM20, PPM10)
    newp = os.path.join(d, "s_CGA_scoreB.tsv")
    ns["_cga_record_tsv"](files, newp)
    assert files["_cga_parameters"]["run"]["tsv"] == os.path.normpath(newp)


@pytest.mark.parametrize("gate, th, nd", [
    (None, TH_DEFAULT, True),                                   # gate off
    (300.0, TH_DEFAULT, True),
    (10.0, dict(TH_DEFAULT, max_abs_ppm=300.0), True),
    (10.0, dict(TH_DEFAULT, ion_ppm=500.0), True),
    (10.0, dict(TH_DEFAULT, min_ion_score=0.2), False),        # not a tolerance
])
def test_trainable_non_default_flag(ns, gate, th, nd):
    d, files = _sample(ns)
    tr = ns["_cga_record_trainable"](files, os.path.join(d, "t.csv"), None, th, NEG_DEFAULT, gate)
    assert tr["tolerance_non_default"] is nd


def test_load_restores_and_reexports_and_registers(ns):
    d, files = _sample(ns)
    ns["_cga_record_run"](files, files["pseudolabel_csv"], WIDE_P, WIDE_F)
    tcsv = os.path.join(d, "t.csv"); open(tcsv, "w").write("x")
    files["trainable_csv"] = tcsv
    ns["_cga_record_trainable"](files, tcsv, files["pseudolabel_csv"], dict(TH_DEFAULT, ion_ppm=500.0),
                                NEG_DEFAULT, 300.0)
    v1 = ns["build_method_v1_from_tree"]("E", "S", force_family="CGA")
    mp = os.path.join(d, "S.CGA.method.json")
    json.dump(v1, open(mp, "w"), indent=2)
    # fresh session: no tree key, no cache, empty registry
    files.pop("_cga_parameters"); files.pop("_method_v1_cache_CGA", None)
    ns["_CGA_TRAINABLE_REGISTRY"].clear()
    _, _, tree, _ = ns["normalize_method_json"](json.load(open(mp)), mp, base_dir=d)
    assert tree["_cga_parameters"] == v1["parameters"]["cga"]
    assert ns["_CGA_TRAINABLE_REGISTRY"] == {}                 # CP-C.1 F1: normalising alone registers nothing
    assert ns["_find_cga_record_for_trainable"](tcsv) == (None, "not found")
    for k, v in tree.items():
        files[k] = v
    ns["_cga_register_loaded"](tree, mp)                       # what both loaders do once the import is accepted
    assert ns["build_method_v1_from_tree"]("E", "S", force_family="CGA")["parameters"]["cga"] == v1["parameters"]["cga"]
    tr, src = ns["_find_cga_record_for_trainable"](tcsv)
    assert tr == v1["parameters"]["cga"]["trainable"] and src == os.path.normpath(mp)
    # an older method JSON (parameters.cga == {}) adds no key
    obj = json.load(open(mp)); obj["parameters"]["cga"] = {}
    _, _, tree2, _ = ns["normalize_method_json"](obj, mp, base_dir=d)
    assert "_cga_parameters" not in tree2


def test_lookup_scans_method_json_in_folder(ns, tmp_path):
    tcsv = tmp_path / "t.csv"; tcsv.write_text("x")
    rec = {"csv": os.path.normpath(str(tcsv)), "csv_sha256": ns["file_sha256"](str(tcsv)),
           "precursor_gate_ppm": 300.0, "max_abs_ppm": 20.0, "feature_ion_ppm": 500.0, "cga_run": None}
    (tmp_path / "a.CGA.method.json").write_text(json.dumps({"parameters": {"cga": {"trainable": dict(rec, csv="zz.csv")}}}))
    (tmp_path / "b.CGA.method.json").write_text(json.dumps({"parameters": {"cga": {"trainable": rec}}}))
    (tmp_path / "broken.method.json").write_text("{")
    tr, src = ns["_find_cga_record_for_trainable"](str(tcsv), registry={}, search_dirs=[str(tmp_path)])
    assert tr == rec and src.endswith("b.CGA.method.json")
    (tmp_path / "pooled.csv").write_text("p")
    assert ns["_find_cga_record_for_trainable"](str(tmp_path / "pooled.csv"), registry={},
                                                search_dirs=[str(tmp_path)]) == (None, "not found")


def test_prefill_only_non_default(ns):
    p = ns["_model_tolerance_prefill"]
    assert p(None) == {} and p({}) == {}
    assert p({"precursor_gate_ppm": 10.0, "feature_ion_ppm": 10.0}) == {}
    assert p({"precursor_gate_ppm": None, "feature_ion_ppm": 10.0}) == {"gate_on": False}   # CP-C.1 F3
    assert p({"feature_ion_ppm": 10.0}) == {}                                                 # key absent: unknown
    assert p({"precursor_gate_ppm": "bad", "feature_ion_ppm": float("nan")}) == {}            # CP-C.1 F5
    assert p({"precursor_gate_ppm": 300.0, "feature_ion_ppm": 500.0}) == {"gate_ppm": 300.0, "feature_ppm": 500.0}
    assert ns["_ppm_tag"](20) == "20" and ns["_ppm_tag"](20.0) == "20" and ns["_ppm_tag"](12.5) == "12.5"
    for v in (12.345678, 0.1 + 0.2, 499.9999999):                                          # CP-C.1 F9
        assert float(ns["_ppm_tag"](v)) == v


def test_training_run_json_carries_the_record(tmp_path):
    try:
        n = X.load_v14(extract=("snapshot_training_run",),
                       stubs={"datetime": _dt, "Dict": typing.Dict, "Any": typing.Any, "Optional": typing.Optional})
    except X.AnchorMissing as e:
        pytest.skip(f"v14 anchor moved (coverage LOST, update this test): {e}")
    snap = n["snapshot_training_run"]
    kw = dict(effective_params={}, column_order=["a"], classes=["Glycan"], inputs={}, hashes={}, versions={})
    out = snap(str(tmp_path / "a"), **kw)
    assert json.load(open(out))["training_run"]["cga_tolerances"] == {"status": "not found"}
    blk = {"status": "found", "source": "x", "trainable": {"precursor_gate_ppm": 300.0, "feature_ion_ppm": 500.0}}
    out = snap(str(tmp_path / "b"), cga_tolerances=blk, **kw)
    assert json.load(open(out))["training_run"]["cga_tolerances"] == blk
    # and the prediction side reads it back
    ns2 = X.load_v14()
    rec = ns2["_read_model_cga_record"](str(tmp_path / "b" / "m_rf_model.skops"))
    assert ns2["_model_tolerance_prefill"](rec.get("trainable")) == {"gate_ppm": 300.0, "feature_ppm": 500.0}
    assert ns2["_read_model_cga_record"](str(tmp_path / "none" / "m.skops")) == {}


def test_override_notices(ns):
    g, f = ns["_gate_override_note"], ns["_feature_override_note"]
    pre = {"gate_ppm": 300.0, "feature_ppm": 500.0}
    assert g({}, True, 10.0) is None and g({}, False, 10.0) is None          # no record: today's behaviour
    assert g(pre, True, 300.0) is None
    assert "300 ppm" in g(pre, True, 10.0) and "uses 10 ppm" in g(pre, True, 10.0)
    assert "uses off" in g(pre, False, 300.0)
    assert f({}, 20) is None and f(pre, 500.0) is None
    assert "500 ppm feature window" in f(pre, 20) and "use 20 ppm" in f(pre, 20)


# ---------- CP-C.1 (Codex round 1) ----------

def _built(ns, d, files, content="trainable", gate=300.0, ion_ppm=500.0):
    tcsv = os.path.join(d, "t.csv"); open(tcsv, "w").write(content)
    ns["_cga_record_run"](files, files["pseudolabel_csv"], WIDE_P, WIDE_F)
    tr = ns["_cga_record_trainable"](files, tcsv, files["pseudolabel_csv"], dict(TH_DEFAULT, ion_ppm=ion_ppm),
                                     NEG_DEFAULT, gate)
    return tcsv, tr


def test_f1_skipped_import_cannot_replace_the_record(ns):
    d, files = _sample(ns)
    tcsv, tr = _built(ns, d, files)
    # an older default-valued method for the same CSV path is normalised (and then skipped by the loader)
    old = ns["build_method_v1_from_tree"]("E", "S", force_family="CGA")
    old["parameters"]["cga"]["trainable"] = dict(tr, precursor_gate_ppm=10.0, feature_ion_ppm=10.0,
                                                 tolerance_non_default=False)
    mp = os.path.join(d, "old.CGA.method.json"); json.dump(old, open(mp, "w"))
    ns["normalize_method_json"](json.load(open(mp)), mp, base_dir=d)
    got, src = ns["_find_cga_record_for_trainable"](tcsv)
    assert got == tr and src == "built in this session"


def test_f1_rebuilt_csv_at_same_path_is_stale(ns, tmp_path):
    d, files = _sample(ns)
    tcsv, tr = _built(ns, d, files)
    open(tcsv, "w").write("rebuilt elsewhere with other settings")
    assert ns["_find_cga_record_for_trainable"](tcsv) == (None, "stale")
    # a method JSON in the folder still describing the old file is stale too, never used
    m = {"parameters": {"cga": {"trainable": tr}}}
    (tmp_path / "a.CGA.method.json").write_text(json.dumps(m))
    assert ns["_find_cga_record_for_trainable"](tcsv, registry={}, search_dirs=[str(tmp_path)]) == (None, "stale")


def test_f1_two_method_files_disagree_is_ambiguous(ns, tmp_path):
    d, files = _sample(ns)
    tcsv, tr = _built(ns, d, files)
    (tmp_path / "a.CGA.method.json").write_text(json.dumps({"parameters": {"cga": {"trainable": tr}}}))
    (tmp_path / "b.CGA.method.json").write_text(json.dumps({"parameters": {"cga": {"trainable": dict(tr, feature_ion_ppm=20.0)}}}))
    assert ns["_find_cga_record_for_trainable"](tcsv, registry={}, search_dirs=[str(tmp_path)]) == (None, "ambiguous")
    (tmp_path / "b.CGA.method.json").write_text(json.dumps({"parameters": {"cga": {"trainable": tr}}}))
    assert ns["_find_cga_record_for_trainable"](tcsv, registry={}, search_dirs=[str(tmp_path)])[0] == tr


def test_f4_null_parameters_in_cached_method_is_repaired(ns):
    d, files = _sample(ns)
    cache = ns["build_method_v1_from_tree"]("E", "S", force_family="CGA")
    cache["parameters"] = None
    files["_method_v1_cache_CGA"] = cache
    tcsv, tr = _built(ns, d, files)
    assert files["_method_v1_cache_CGA"]["parameters"]["cga"]["trainable"] == tr
    assert set(files["_method_v1_cache_CGA"]["parameters"]) == {"mas", "cga", "scoring", "ml"}


def test_f5_malformed_records_are_dropped_not_fatal(ns):
    d, files = _sample(ns)
    v1 = ns["build_method_v1_from_tree"]("E", "S", force_family="CGA")
    for bad in ({"run": {"tsv": files["pseudolabel_csv"]}},                         # no tolerances
                {"run": {"precursor_tolerance": {"value": "x", "unit": "ppm"}, "fragment_tolerance": PPM10}},
                {"trainable": {"csv": "t.csv", "precursor_gate_ppm": "bad"}},
                "not a dict"):
        v1["parameters"]["cga"] = bad
        _, _, tree, _ = ns["normalize_method_json"](json.loads(json.dumps(v1)), "m.json", base_dir=d)
        assert "_cga_parameters" not in tree
    v1["parameters"] = None
    _, _, tree, _ = ns["normalize_method_json"](json.loads(json.dumps(v1)), "m.json", base_dir=d)
    assert "_cga_parameters" not in tree
    # a malformed run already on the tree: the trainable is still recorded, with cga_run unknown
    files["_cga_parameters"] = {"run": {"tsv": files["pseudolabel_csv"]}}
    tcsv = os.path.join(d, "t.csv"); open(tcsv, "w").write("x")
    tr = ns["_cga_record_trainable"](files, tcsv, files["pseudolabel_csv"], TH_DEFAULT, NEG_DEFAULT, 10.0)
    assert tr["cga_run"] is None


def test_f2_f3_model_switch_keeps_user_values(ns):
    up = ns["_gate_prefill_update"]
    # fresh window, model A records 300 ppm: default "10" is replaced
    f, on, kept, auto = up({}, "10", True, {"gate_ppm": 300.0})
    assert (f, on, kept, auto) == ("300", None, [], {"field": "300"})
    # user typed 150, model B records 200: the user's value is kept, reported, and owned from now on
    f, on, kept, auto = up({"field": "300"}, "150", True, {"gate_ppm": 200.0})
    assert f is None and kept == ["precursor gate 150 ppm"] and auto == {"field_user": True}
    # untouched automatic value, model without record: back to the default
    f, on, kept, auto = up({"field": "300"}, "300", True, {})
    assert (f, on, kept, auto) == ("10", None, [], {})
    # model trained with the gate off: untick it; next default model ticks it again
    f, on, kept, auto = up({}, "10", True, {"gate_on": False})
    assert (f, on, auto) == (None, False, {"gate_on": False})
    f, on, kept, auto = up({"gate_on": False}, "10", False, {})
    assert (f, on) == (None, True)
    # user re-ticked the gate after an automatic untick: kept
    f, on, kept, auto = up({"gate_on": False}, "10", True, {"gate_on": False})
    assert on is None and kept == ["precursor gate on"] and auto == {"on_user": True}


def test_f3_gate_off_notice(ns):
    g = ns["_gate_override_note"]
    assert g({"gate_on": False}, False, 10.0) is None
    assert "gate off" in g({"gate_on": False}, True, 10.0) and "uses 10 ppm" in g({"gate_on": False}, True, 10.0)



# ---------- CP-C.2 (Codex round 2) ----------

def test_r2f1_conflicting_records_for_same_bytes_are_ambiguous(ns, tmp_path):
    d, files = _sample(ns)
    tcsv, tr = _built(ns, d, files)                                     # in-session build: 300 / 500
    other = dict(tr, feature_ion_ppm=20.0)
    ns["_cga_register_loaded"]({"_cga_parameters": {"trainable": other}}, "accepted.CGA.method.json")
    assert ns["_find_cga_record_for_trainable"](tcsv) == (None, "ambiguous")   # no precedence for the build
    # registry agrees with itself, but a folder method disagrees
    ns["_CGA_TRAINABLE_REGISTRY"].clear()
    ns["_cga_registry_add"](tr, "built in this session")
    (tmp_path / "x.CGA.method.json").write_text(json.dumps({"parameters": {"cga": {"trainable": other}}}))
    assert ns["_find_cga_record_for_trainable"](tcsv, search_dirs=[str(tmp_path)]) == (None, "ambiguous")
    # the autosaved method holding the same record is not a conflict
    (tmp_path / "x.CGA.method.json").write_text(json.dumps({"parameters": {"cga": {"trainable": tr}}}))
    assert ns["_find_cga_record_for_trainable"](tcsv, search_dirs=[str(tmp_path)]) == (tr, "built in this session")
    # the same record registered twice stays one entry
    ns["_cga_registry_add"](tr, "again")
    assert len(next(iter(ns["_CGA_TRAINABLE_REGISTRY"].values()))) == 1


def test_r2f2_digest_is_the_one_given(ns):
    d, files = _sample(ns)
    tcsv = os.path.join(d, "t.csv"); open(tcsv, "w").write("A")
    sha_a = ns["file_sha256"](tcsv)
    ns["_cga_record_run"](files, files["pseudolabel_csv"], WIDE_P, WIDE_F)
    open(tcsv, "w").write("B")                                          # replaced after the write, before recording
    tr = ns["_cga_record_trainable"](files, tcsv, files["pseudolabel_csv"], TH_DEFAULT, NEG_DEFAULT, 300.0,
                                     csv_sha256=sha_a)
    assert tr["csv_sha256"] == sha_a
    assert ns["_find_cga_record_for_trainable"](tcsv) == (None, "stale")           # file is B now
    assert ns["_find_cga_record_for_trainable"](tcsv, csv_sha256=sha_a)[0] == tr    # training read A


def test_r2f3_user_ownership_survives_several_switches(ns):
    up = ns["_gate_prefill_update"]
    auto = {}
    f, on, kept, auto = up(auto, "10", True, {"gate_ppm": 300.0}); assert f == "300"          # model A
    field = "10"                                                                              # user types 10
    f, on, kept, auto = up(auto, field, True, {"gate_ppm": 200.0}); assert f is None          # model B: kept
    f, on, kept, auto = up(auto, field, True, {"gate_ppm": 400.0})                            # model C
    assert f is None and kept == ["precursor gate 10 ppm"] and auto.get("field_user") is True
    f, on, kept, auto = up(auto, field, True, {}); assert f is None and kept == []            # no record: still user's
    # checkbox: A unticks automatically, user ticks, B keeps it, C (gate off again) keeps it too
    auto = {}
    f, on, kept, auto = up(auto, "10", True, {"gate_on": False}); assert on is False
    f, on, kept, auto = up(auto, "10", True, {"gate_on": False}); assert on is None
    f, on, kept, auto = up(auto, "10", True, {"gate_on": False}); assert on is None and auto.get("on_user") is True


def test_r2f5_gate_outcome_vs_model(ns):
    m = ns["_gate_outcome_matches_model"]
    assert m({}, True, 10.0) is None and m({}, False, 10.0) is None
    assert m({"gate_ppm": 300.0}, True, 300.0) is True
    assert m({"gate_ppm": 300.0}, False, 300.0) is False          # requested 300 but the gate was skipped
    assert m({"gate_ppm": 300.0}, True, 10.0) is False
    assert m({"gate_on": False}, False, 10.0) is True
    assert m({"gate_on": False}, True, 10.0) is False


def test_r2f6_oversized_numbers_are_unusable(ns):
    big = json.loads('{"precursor_gate_ppm": 1' + "0" * 400 + ', "feature_ion_ppm": 1' + "0" * 400 + '}')
    assert ns["_model_tolerance_prefill"](big) == {}
    assert ns["_pos_num"](10 ** 400) is None
    assert ns["_clean_cga_parameters"]({"run": {"precursor_tolerance": {"value": 10 ** 400, "unit": "ppm"},
                                                "fragment_tolerance": PPM10}}) == {}


# ---------- CP-C.3 (Codex round 3) ----------

def test_r3f1_failed_capture_is_never_rehashed(ns):
    d, files = _sample(ns)
    tcsv = os.path.join(d, "t.csv"); open(tcsv, "w").write("B (replacement written after the failed capture)")
    ns["_cga_record_run"](files, files["pseudolabel_csv"], WIDE_P, WIDE_F)
    tr = ns["_cga_record_trainable"](files, tcsv, files["pseudolabel_csv"], TH_DEFAULT, NEG_DEFAULT, 300.0,
                                     csv_sha256=None)                  # the hash right after the write failed
    assert tr["csv_sha256"] is None                                     # not replaced by a later hash of B
    assert ns["_CGA_TRAINABLE_REGISTRY"] == {}                          # no identity => not registered
    assert ns["_find_cga_record_for_trainable"](tcsv) == (None, "not found")
    assert "trainable" not in ns["_clean_cga_parameters"](files["_cga_parameters"])   # and not restorable
    # argument omitted (other callers): hashed now, as before
    tr2 = ns["_cga_record_trainable"](files, tcsv, files["pseudolabel_csv"], TH_DEFAULT, NEG_DEFAULT, 300.0)
    assert tr2["csv_sha256"] == ns["file_sha256"](tcsv)


class _TraceVar:
    """Tk-variable stand-in whose set() runs write traces, like tkinter."""
    def __init__(self, v):
        self.v, self.cbs = v, []

    def trace_add(self, mode, cb):
        assert mode == "write"
        self.cbs.append(cb)

    def get(self):
        return self.v

    def set(self, v):
        self.v = v
        for cb in self.cbs:
            cb("name", "", "write")


def test_r3f2_edit_and_restore_still_counts_as_the_users(ns):
    field, gate = _TraceVar("10"), _TraceVar(True)
    auto, guard = {}, {"programmatic": False}
    ns["_attach_gate_ownership_traces"](field, gate, auto, guard)
    up = ns["_gate_prefill_update"]

    def load_model(prefill):                       # what select_model_file does
        f, on, kept, new = up(auto, field.get(), gate.get(), prefill)
        guard["programmatic"] = True
        try:
            if f is not None:
                field.set(f)
            if on is not None:
                gate.set(on)
        finally:
            guard["programmatic"] = False
        auto.clear(); auto.update(new)
        return kept

    load_model({"gate_ppm": 300.0, "gate_on": False})
    assert field.get() == "300" and gate.get() is False and "field_user" not in auto and "on_user" not in auto
    field.set("150"); field.set("300")             # user edits, then restores the automatic value by hand
    gate.set(True); gate.set(False)                # same for the checkbox
    load_model({"gate_ppm": 200.0})
    assert field.get() == "300" and gate.get() is False       # both kept: the user owns them
    load_model({})
    assert field.get() == "300" and gate.get() is False and auto.get("field_user") and auto.get("on_user")
    # an untouched window still follows the models
    field2, gate2 = _TraceVar("10"), _TraceVar(True)
    auto, guard = {}, {"programmatic": False}
    ns["_attach_gate_ownership_traces"](field2, gate2, auto, guard)
    field, gate = field2, gate2
    load_model({"gate_ppm": 300.0}); assert field.get() == "300"
    load_model({}); assert field.get() == "10" and not auto


# ---------------------------------------------------------------- CP-C.4 (2026-11-01): CGA -> Trainable precursor gate field

def test_cp_c4_gate_field_values(ns):
    """The new field gives the same precursor_gate_ppm as v1.10's fixed float("10") at its default, None when off."""
    g = ns["_trainable_gate_ppm"]
    assert g(True, "10") == float("10") == 10.0 and type(g(True, "10")) is float
    assert g(True, " 300 ") == 300.0 and g(True, "12.5") == 12.5
    assert g(False, "10") is None and g(False, "not a number") is None      # off: the text is never read
    for bad in ("", "0", "-5", "abc", "inf", "nan"):
        with pytest.raises(ValueError, match="precursor gate"):
            g(True, bad)


def _modal_source():
    import ast
    with open(X.V14_PATH, "rb") as f:
        src = f.read().decode("utf-8")
    tree = ast.parse(src)
    fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "open_pl_to_trainable_modal")
    return ast, fn


def test_cp_c4_gate_widgets_wired_with_v1_10_defaults():
    ast, fn = _modal_source()
    calls = [n for n in ast.walk(fn) if isinstance(n, ast.Call)]

    def attr(c):
        return c.func.attr if isinstance(c.func, ast.Attribute) else getattr(c.func, "id", None)

    def kw(c, name):
        v = next((k.value for k in c.keywords if k.arg == name), None)
        return v.id if isinstance(v, ast.Name) else None

    def first_arg(c):
        return c.args[0].id if c.args and isinstance(c.args[0], ast.Name) else None

    chk = [c for c in calls if attr(c) == "Checkbutton" and kw(c, "variable") == "apply_precursor_gate_var"]
    ent = [c for c in calls if attr(c) == "Entry" and kw(c, "textvariable") == "precursor_ppm_var"]
    assert len(chk) == 1 and first_arg(chk[0]) == "thresh_box"
    assert len(ent) == 1 and first_arg(ent[0]) == "thresh_box"
    # defaults unchanged: BooleanVar(value=True), StringVar(value="10")
    assigns = {t.id: n.value for n in ast.walk(fn) if isinstance(n, ast.Assign)
               for t in n.targets if isinstance(t, ast.Name)}
    for name, ctor, val in (("apply_precursor_gate_var", "BooleanVar", True), ("precursor_ppm_var", "StringVar", "10")):
        c = assigns[name]
        assert attr(c) == ctor and c.keywords[0].arg == "value" and c.keywords[0].value.value == val
    # run_once reads the pair only through the checked helper, and hands the result to the builder and the recorder
    run_once = next(n for n in ast.walk(fn) if isinstance(n, ast.FunctionDef) and n.name == "run_once")
    helper = [c for c in ast.walk(run_once) if isinstance(c, ast.Call) and getattr(c.func, "id", None) == "_trainable_gate_ppm"]
    assert len(helper) == 1
    reads = [n for n in ast.walk(run_once) if isinstance(n, ast.Name) and n.id in ("apply_precursor_gate_var", "precursor_ppm_var")]
    assert len(reads) == 2                                               # the two helper arguments only
    builder = next(c for c in ast.walk(run_once) if isinstance(c, ast.Call) and getattr(c.func, "id", None) == "build_trainable_from_CGA")
    assert kw(builder, "precursor_gate_ppm") == "pg_ppm"
    rec = next(c for c in ast.walk(run_once) if isinstance(c, ast.Call) and getattr(c.func, "id", None) == "_cga_record_trainable")
    assert isinstance(rec.args[5], ast.Name) and rec.args[5].id == "pg_ppm"
    # a bad entry stops before the build (the error dialog, then return)
    tries = [t for t in ast.walk(run_once) if isinstance(t, ast.Try)
             and any(isinstance(h.type, ast.Name) and h.type.id == "ValueError" for h in t.handlers)
             and any(c in helper for s in t.body for c in ast.walk(s))]
    assert len(tries) == 1 and any(isinstance(s, ast.Return) for h in tries[0].handlers for s in h.body)
    assert tries[0].lineno < builder.lineno


def test_cp_c4_widened_gate_is_recorded_and_prefills(ns):
    """A widened gate from the field is non-default, so the model prefills the prediction gate from it."""
    d, files = _sample(ns)
    tr = ns["_cga_record_trainable"](files, os.path.join(d, "t.csv"), None, TH_DEFAULT, NEG_DEFAULT,
                                     ns["_trainable_gate_ppm"](True, "300"))
    assert tr["precursor_gate_ppm"] == 300.0 and tr["tolerance_non_default"] is True
    assert ns["_model_tolerance_prefill"](tr).get("gate_ppm") == 300.0
    tr10 = ns["_cga_record_trainable"](files, os.path.join(d, "t.csv"), None, TH_DEFAULT, NEG_DEFAULT,
                                       ns["_trainable_gate_ppm"](True, "10"))
    assert tr10["tolerance_non_default"] is False and "gate_ppm" not in ns["_model_tolerance_prefill"](tr10)
