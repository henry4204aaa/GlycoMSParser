"""S1 — retained regression tests for the v14 wiring that has no importable entry point (Codex round 2 ask).

mspfileloaderv14.py creates the Tk root at module level, so it cannot be imported. This test executes the
module source up to (not including) `root = tk.Tk()`, which defines every module-level function including
`normalize_method_json`, `_glytoucan_fill_df` and `_glytoucan_record`, then extracts the closure-local
`build_method_v1_from_tree` (a nested def inside `open_prepare_dataset_window` that only uses module-level
names) and executes it at module level. No GUI, no network.

Covers:
  F1  a loaded method JSON carrying artifacts.glytoucan_resolution restores `_glytoucan_resolution_<FAMILY>`
      on the tree entry, and the rebuilt method re-emits the identical block (Export / autosave paths)
  F2  a new run without a fill clears tree key + cached block; the rebuilt JSON has no block
  family keying (CGA block never in MAS method and vice-versa); block-less v1 adds no key
  option off -> helper returns None and the frame is untouched; missing DB -> 'db_missing'

Run from repo root:  python3 -m pytest tests/test_s1_v14_headless.py -q
"""
import json
import os
import sys
import tempfile

import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
sys.path.insert(0, SRC)

import msp_glytoucan_resolver as R  # noqa: E402


def _stub(url, body, timeout):
    items = json.loads(body.decode("utf-8"))
    out = []
    for it in items:
        echo = {k: v for k, v in it.items() if v != "0"}
        echo["wurcs"] = "WURCS=2.0/stub"
        echo["id"] = "G0000001"
        out.append(echo)
    return json.dumps({"data": out}).encode("utf-8")


def _find_anchors(lines):
    """Only a MOVED/RENAMED anchor is a legitimate skip (Codex R3-F4); everything else must fail loudly."""
    def first(pred, start=0):
        for i in range(start, len(lines)):
            if pred(lines[i]):
                return i
        return None
    root_i = first(lambda l: l.startswith("root = tk.Tk()"))
    bi = first(lambda l: l.startswith("    def build_method_v1_from_tree("))
    ei = first(lambda l: l.strip() == "return v1", bi or 0)
    missing = [n for n, v in (("root = tk.Tk()", root_i), ("def build_method_v1_from_tree", bi), ("return v1", ei)) if v is None]
    if missing:
        pytest.skip(f"v14 anchors not found (coverage LOST, update this test): {missing}")
    return root_i, bi, ei + 1


@pytest.fixture(scope="module")
def v14ns(tmp_path_factory):
    """Namespace with v14's module-level names + extracted build_method_v1_from_tree.

    Codex R3-F3: the GLYCOMSP_REFERENCE_DB override is scoped to this fixture with a MonkeyPatch context that
    is undone on success, failure and skip, restoring the exact prior state (including 'unset').
    Codex R3-F4: only missing anchors skip; compile/import/exec errors FAIL with their traceback."""
    with pytest.MonkeyPatch.context() as mp:
        tmp = tmp_path_factory.mktemp("v14headless")
        dbp = tmp / "ref.db"
        db = R.ReferenceDB(dbp)
        R.prebuild(["F1H5N4S2"], db, tmp / "l.jsonl", http_post=_stub, sleep=lambda s: None)
        db.close()
        mp.setenv("GLYCOMSP_REFERENCE_DB", str(dbp))
        path = os.path.join(SRC, "mspfileloaderv14.py")
        src = open(path, "rb").read().decode("utf-8")
        lines = src.splitlines(keepends=True)
        root_i, bi, ei = _find_anchors(lines)
        ns = {"__name__": "v14_headless", "__file__": path}
        old_hook = sys.excepthook
        mp.chdir(SRC)
        try:
            exec(compile("".join(lines[:root_i]), "v14_headless", "exec"), ns)   # any error here FAILS the module
            body = "".join((l[4:] if l.startswith("    ") else l) for l in lines[bi:ei])
            exec(compile(body, "build_method_v1_from_tree_extract", "exec"), ns)
        finally:
            sys.excepthook = old_hook
        ns.setdefault("experiment_projects", {})
        ns["_tmp"] = str(tmp)
        ns["_dbp"] = str(dbp)
        yield ns


def _sample(ns):
    d = tempfile.mkdtemp(dir=ns["_tmp"])
    files = {}
    for key, name in (("csv", "s.csv"), ("metadata", "s.json"), ("ionlist_path", "ion.xlsx"),
                      ("insilico_csv", "ins.csv"), ("excel", "mas.xlsx")):
        p = os.path.join(d, name)
        open(p, "w").write("x")
        files[key] = p
    files["pseudolabel_csv"] = files["csv"]
    ns["experiment_projects"].clear()
    ns["experiment_projects"]["E"] = {"samples": {"S": files}}
    return d, files


def test_option_off_is_noop_and_missing_db_is_soft(v14ns):
    fill, GF = v14ns["_glytoucan_fill_df"], v14ns["GLYTOUCAN_FILL"]
    GF["enabled"] = False
    df = pd.DataFrame({"composition": ["F1H5N4S2"], "GlyToucan ID": [""], "WURCS": [""]})
    assert fill(df, "composition") is None and df["WURCS"].tolist() == [""]
    GF["enabled"] = True
    try:
        with pytest.MonkeyPatch.context() as mp:
            missing = os.path.join(v14ns["_tmp"], "nope", "x.db")
            mp.setenv("GLYCOMSP_REFERENCE_DB", missing)
            blk = fill(pd.DataFrame({"composition": ["F1H5N4S2"]}), "composition")
            assert blk["status"] == "db_missing" and not os.path.exists(os.path.dirname(missing))
        assert os.environ["GLYCOMSP_REFERENCE_DB"] == v14ns["_dbp"]     # inner override undone
    finally:
        GF["enabled"] = False


def test_f2_fill_off_rerun_clears_block_and_family_keying(v14ns):
    fill, rec, build, GF = (v14ns[k] for k in ("_glytoucan_fill_df", "_glytoucan_record", "build_method_v1_from_tree", "GLYTOUCAN_FILL"))
    _, files = _sample(v14ns)
    GF["enabled"] = True
    try:
        blk = fill(pd.DataFrame({"composition": ["F1H5N4S2"], "GlyToucan ID": [""], "WURCS": [""]}), "composition")
        assert blk["status"] == "filled" and blk["db"]["n_entries"] == 1
        rec(files, "CGA", blk)
        v1 = build("E", "S", force_family="CGA")
        assert v1["artifacts"]["glytoucan_resolution"] == blk
        files["_method_v1_cache_CGA"] = v1
        assert "glytoucan_resolution" not in build("E", "S", force_family="MAS")["artifacts"]
        GF["enabled"] = False
        rec(files, "CGA", fill(pd.DataFrame({"composition": ["F1H5N4S2"]}), "composition"))   # H1 with option off
        assert "_glytoucan_resolution_CGA" not in files
        assert "glytoucan_resolution" not in files["_method_v1_cache_CGA"]["artifacts"]
        assert "glytoucan_resolution" not in build("E", "S", force_family="CGA")["artifacts"]
    finally:
        GF["enabled"] = False


def test_f1_load_restores_provenance_and_export_reemits(v14ns):
    fill, rec, build, norm, GF = (v14ns[k] for k in ("_glytoucan_fill_df", "_glytoucan_record", "build_method_v1_from_tree", "normalize_method_json", "GLYTOUCAN_FILL"))
    d, files = _sample(v14ns)
    GF["enabled"] = True
    try:
        blk = fill(pd.DataFrame({"composition": ["F1H5N4S2"], "GlyToucan ID": [""], "WURCS": [""]}), "composition")
    finally:
        GF["enabled"] = False
    rec(files, "CGA", blk)
    v1 = build("E", "S", force_family="CGA")
    mp = os.path.join(d, "S.CGA.method.json")
    json.dump(v1, open(mp, "w"), indent=2)
    # simulate a fresh session: no tree key, no cache
    files.pop("_glytoucan_resolution_CGA"); files.pop("_method_v1_cache_CGA", None)
    _, _, tree, v1obj = norm(json.load(open(mp)), mp, base_dir=d)
    assert tree["_glytoucan_resolution_CGA"] == v1["artifacts"]["glytoucan_resolution"]
    for k, v in tree.items():          # what both loaders do
        files[k] = v
    assert build("E", "S", force_family="CGA")["artifacts"]["glytoucan_resolution"] == v1["artifacts"]["glytoucan_resolution"]
    assert "glytoucan_resolution" not in build("E", "S", force_family="MAS")["artifacts"]
    # MAS-family file restores under the MAS key only
    v1m = build("E", "S", force_family="MAS"); v1m["artifacts"]["glytoucan_resolution"] = dict(blk)
    mp2 = os.path.join(d, "S.MAS.method.json"); json.dump(v1m, open(mp2, "w"))
    _, _, tree2, _ = norm(json.load(open(mp2)), mp2, base_dir=d)
    assert "_glytoucan_resolution_MAS" in tree2 and "_glytoucan_resolution_CGA" not in tree2
    # block-less v1 adds no key (v1.10 tree entry unchanged)
    obj = json.load(open(mp2)); obj["artifacts"] = {}
    _, _, tree3, _ = norm(obj, mp2, base_dir=d)
    assert not any(k.startswith("_glytoucan") for k in tree3)


def test_r5f1_cached_cga_method_links_new_trainable_with_its_block(v14ns, tmp_path):
    """Codex R5-F1: after a trainable write, the cached CGA method (reused verbatim by autosave) must carry the NEW
    trainable path together with the H5 block; without a fill the path still updates and the old block is kept."""
    fill, rec, build, link, GF = (v14ns[k] for k in ("_glytoucan_fill_df", "_glytoucan_record", "build_method_v1_from_tree", "_cga_cache_link_trainable", "GLYTOUCAN_FILL"))
    d, files = _sample(v14ns)
    # (a) cache built at CGA completion: no trainable artifact yet
    files.pop("trainable_csv", None)
    cache = build("E", "S", force_family="CGA")
    files["_method_v1_cache_CGA"] = cache
    assert "trainable_csv" not in cache["artifacts"]
    new1 = os.path.join(d, "t1.csv"); open(new1, "w").write("x")
    files["trainable_csv"] = new1
    link(files, new1)
    GF["enabled"] = True
    try:
        blk = fill(pd.DataFrame({"Structure": ["F1H5N4S2"], "GlyToucan ID": [""], "WURCS": [""]}), "Structure")
    finally:
        GF["enabled"] = False
    rec(files, "CGA", blk)
    saved = files["_method_v1_cache_CGA"]            # what save_method_v1_for_sample serialises
    assert saved["artifacts"]["trainable_csv"]["path"] == os.path.normpath(new1)
    assert saved["artifacts"]["glytoucan_resolution"] == blk
    # (b) older trainable path in the cache, then a second conversion with the option OFF:
    #     path must move to the new file, the previous block is kept (approved no-clear rule)
    new2 = os.path.join(d, "t2.csv"); open(new2, "w").write("x")
    files["trainable_csv"] = new2
    link(files, new2)
    assert saved["artifacts"]["trainable_csv"]["path"] == os.path.normpath(new2)
    assert saved["artifacts"]["glytoucan_resolution"] == blk
    # (c) no cache: helper is a no-op and the builder reads the tree on the next save
    files.pop("_method_v1_cache_CGA")
    link(files, new2)
    assert build("E", "S", force_family="CGA")["artifacts"]["trainable_csv"]["path"] == os.path.normpath(new2)
