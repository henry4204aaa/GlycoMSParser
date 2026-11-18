"""S1 (v1.20, T1) — GlyTouCan / WURCS resolver: store, client, prebuild, run-time fill.

Spec: handoff/WURCS_designspec20260921.md. No network: the HTTP layer is stubbed via the
`http_post` injection point, and `sleep` is captured so the 10 s pause is asserted, not waited.

Covers the module half of the done-when list:
  - canonicalisation + symbol mapping (incl. the S / s collision, K -> KDN alias)
  - store: create / lookup / missing / upsert / identity
  - client: batching <= 30, pause between batches, echo verification, fail-soft on
    transport error and on malformed response
  - prebuild: diff -> fetch -> upsert -> log; second run makes zero calls; --no-network diff
  - fill: empty cells only, MAS value wins, comp column untouched, db_missing path,
    artifact block shape
  - CLI: --no-network diff on a real CGA fixture (if present), --summary

Run from repo root:
    python3 -m pytest tests/test_s1_glytoucan_resolver.py -q
"""
import json
import os
import sys

import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import msp_glytoucan_resolver as R  # noqa: E402


# ---------------------------------------------------------------- helpers
def make_stub(calls, *, id_for=lambda item: "G" + "".join(f"{k}{v}" for k, v in sorted(item.items()))[:7]):
    """HTTP stub that echoes non-zero counts + a stub wurcs; accession only when neu5ac > 0."""
    def stub(url, body, timeout):
        items = json.loads(body.decode("utf-8"))
        calls.append(len(items))
        out = []
        for it in items:
            echo = {k: v for k, v in it.items() if v != "0"}
            echo["wurcs"] = "WURCS=2.0/stub/" + "".join(f"{k}{v}" for k, v in sorted(echo.items()))
            echo["id"] = id_for(echo) if it.get("neu5ac", "0") != "0" else ""
            out.append(echo)
        return json.dumps({"data": out}).encode("utf-8")
    return stub


@pytest.fixture
def tmpdb(tmp_path):
    db = R.ReferenceDB(tmp_path / "ref.db")
    yield db, tmp_path / "log.jsonl"
    db.close()


# ---------------------------------------------------------------- canonicalisation
@pytest.mark.parametrize("raw,expected", [
    ("F1H5N4S2", "F1H5N4S2"),
    ("H5N4F1S2", "F1H5N4S2"),          # order-insensitive
    (" H3 N2 ", "H3N2"),
    ("H5N4K1", "H5N4KDN1"),            # legacy K -> KDN
    ("H5N4KDN1", "H5N4KDN1"),
    ("H3N2s1p1", "H3N2s1p1"),
    ("Non-glycan", ""),
    ("", ""),
    (None, ""),
    (float("nan"), ""),
    ("Hex5HexNAc4", ""),               # not the GlycoMSP grammar
    ("G12345AB", ""),
])
def test_canonical_composition(raw, expected):
    assert R.canonical_composition(raw) == expected


def test_symbol_mapping_and_S_collision():
    assert R.composition_to_api_query("F1H5N4S2") == {"dhex": 1, "hex": 5, "hexnac": 4, "neu5ac": 2}
    q = R.composition_to_api_query("H3N2s1p1A1KDN1G1")
    assert q == {"hex": 3, "hexnac": 2, "S": 1, "P": 1, "hexa": 1, "kdn": 1, "neu5gc": 1}
    with pytest.raises(ValueError):
        R.composition_to_api_query("Non-glycan")
    body = R._api_body_item(q)
    assert set(body) == set(R.API_QUERY_KEYS) and body["Ac"] == "0" and body["neu5ac"] == "0"


# ---------------------------------------------------------------- store
def test_store_create_lookup_missing_upsert(tmpdb):
    db, _ = tmpdb
    assert db.count() == 0
    m = db.meta()
    assert m["schema_version"] == R.SCHEMA_VERSION and m["db_name"] == R.DEFAULT_DB_NAME
    assert "CC BY 4.0" in m["attribution"]
    assert db.missing(["F1H5N4S2", "H5N4F1S2", "Non-glycan", ""]) == ["F1H5N4S2"]
    e = R.Entry.from_api("F1H5N4S2", {"dhex": 1, "hex": 5, "hexnac": 4, "neu5ac": 2}, "WURCS=2.0/x", "G22310AV")
    assert e.status == R.STATUS_RESOLVED
    e2 = R.Entry.from_api("H5N4KDN1", {"hex": 5, "hexnac": 4, "kdn": 1}, "WURCS=2.0/y", "")
    assert e2.status == R.STATUS_WURCS_ONLY
    e3 = R.Entry.from_api("H9N2", {"hex": 9, "hexnac": 2}, None, None)
    assert e3.status == R.STATUS_UNRESOLVABLE
    assert db.upsert([e, e2, e3]) == 3
    assert db.missing(["F1H5N4S2", "H5N4K1", "H9N2", "H1N1"]) == ["H1N1"]  # any status counts as present
    got = db.lookup(["H5N4F1S2"])
    assert got["F1H5N4S2"].accession == "G22310AV"
    ident = db.identity()
    assert ident["n_entries"] == 3 and len(ident["sha256"]) == 64 and ident["version"] == R.DEFAULT_DB_VERSION
    s = db.summary()
    assert s["by_status"] == {"resolved": 1, "wurcs_only": 1, "unresolvable": 1}
    # upsert is idempotent and preserves added_utc
    added = db.get("F1H5N4S2").added_utc
    e.accession = "G22310AV"
    e.added_utc = "1999-01-01T00:00:00Z"
    db.upsert([e])
    assert db.count() == 3 and db.get("F1H5N4S2").added_utc == added


# ---------------------------------------------------------------- client
def test_fetch_batches_pause_and_echo_check():
    calls, slept = [], []
    keys = [f"H{i}N2" for i in range(1, 71)] + ["H5N4F1S2", "H5N4F1S2"]   # 71 distinct
    fr = R.fetch_compositions(keys, http_post=make_stub(calls), sleep=slept.append, batch_size=30)
    assert fr.status == "ok" and fr.n_requested == 71
    assert calls == [30, 30, 11] and slept == [10.0, 10.0]
    assert len(fr.entries) == 71 and fr.api_calls == 3
    by = {e.composition_key: e for e in fr.entries}
    assert by["F1H5N4S2"].status == R.STATUS_RESOLVED and by["H1N2"].status == R.STATUS_WURCS_ONLY


def test_fetch_batch_size_capped_at_30():
    calls = []
    R.fetch_compositions([f"H{i}N2" for i in range(1, 41)], http_post=make_stub(calls), sleep=lambda s: None, batch_size=500)
    assert calls == [30, 10]


def test_fetch_fail_soft_transport_error():
    def down(url, body, timeout):
        raise OSError("no network")
    fr = R.fetch_compositions(["H5N4", "H6N2"], http_post=down, sleep=lambda s: None)
    assert fr.status == "offline" and fr.entries == [] and sorted(fr.unsent) == ["H5N4", "H6N2"]
    assert fr.api_calls == 0 and fr.errors


def test_fetch_partial_on_second_batch_failure():
    calls = []
    ok = make_stub(calls)

    def flaky(url, body, timeout):
        if len(calls) >= 1:
            raise OSError("dropped")
        return ok(url, body, timeout)
    keys = [f"H{i}N2" for i in range(1, 41)]
    fr = R.fetch_compositions(keys, http_post=flaky, sleep=lambda s: None)
    assert fr.status == "partial" and len(fr.entries) == 30 and len(fr.unsent) == 10 and fr.api_calls == 1


def test_fetch_rejects_malformed_and_mismatched_echo():
    def bad_shape(url, body, timeout):
        return b'{"message":"500 Internal Server Error"}'
    fr = R.fetch_compositions(["H5N4"], http_post=bad_shape, sleep=lambda s: None)
    assert fr.status == "offline" and fr.unsent == ["H5N4"]

    def wrong_echo(url, body, timeout):
        return json.dumps({"data": [{"hex": "9", "hexnac": "9", "wurcs": "W", "id": "G1"}]}).encode()
    fr = R.fetch_compositions(["H5N4"], http_post=wrong_echo, sleep=lambda s: None)
    assert fr.entries == [] and fr.unsent == ["H5N4"] and fr.status == "partial"


# ---------------------------------------------------------------- prebuild
def test_prebuild_diff_fetch_log_and_second_run_zero_calls(tmpdb):
    db, log = tmpdb
    calls, slept = [], []
    src = ["F1H5N4S2", "H5N4F1S2", "H5N4", "H3N2s1", "H5N4K1", "Non-glycan", ""]
    r1 = R.prebuild(src, db, log, http_post=make_stub(calls), sleep=slept.append)
    assert r1.status == "ok" and r1.n_source == 4 and r1.n_already == 0 and r1.n_fetched == 4 and r1.api_calls == 1
    assert db.count() == 4 and r1.db["n_entries"] == 4
    lines = [json.loads(l) for l in open(log, encoding="utf-8")]
    assert len(lines) == 4 and {l["op"] for l in lines} == {"api_resolve"}
    assert {l["composition_key"] for l in lines} == {"F1H5N4S2", "H5N4", "H3N2s1", "H5N4KDN1"}
    assert all(l["endpoint"] == R.API_ENDPOINT and l["db_version"] == R.DEFAULT_DB_VERSION for l in lines)
    # second run: zero calls
    r2 = R.prebuild(src, db, log, http_post=make_stub(calls), sleep=slept.append)
    assert r2.status == "nothing_to_do" and r2.api_calls == 0 and calls == [4]
    # incremental: one new composition -> one call with one item
    r3 = R.prebuild(src + ["H9N2"], db, log, http_post=make_stub(calls), sleep=slept.append)
    assert r3.n_already == 4 and r3.n_fetched == 1 and calls == [4, 1]
    # no-network diff stores nothing
    r4 = R.prebuild(["H1N1"], db, log, allow_network=False)
    assert r4.status == "disabled" and r4.unsent == ["H1N1"] and db.get("H1N1") is None
    # offline stores nothing, so the miss is retried next time
    def down(url, body, timeout):
        raise OSError("no network")
    r5 = R.prebuild(["H1N1"], db, log, http_post=down, sleep=slept.append)
    assert r5.status == "offline" and db.get("H1N1") is None
    assert sum(1 for _ in open(log, encoding="utf-8")) == 5


# ---------------------------------------------------------------- run-time fill
def test_annotate_fills_empty_only_and_never_touches_comp(tmpdb):
    db, log = tmpdb
    R.prebuild(["F1H5N4S2", "H3N2s1"], db, log, http_post=make_stub([]), sleep=lambda s: None)
    df = pd.DataFrame({
        "composition": ["F1H5N4S2", "H5N4F1S2", "H3N2s1", "H9N2", "Non-glycan", ""],
        "GlyToucan ID": ["", "USER-ID", "", "", "", ""],
        "WURCS": ["", "", "", "", "", ""],
        "other": [1, 2, 3, 4, 5, 6],
    })
    before_comp = df["composition"].tolist()
    fr = R.annotate_dataframe(df, "composition", db)
    assert df["composition"].tolist() == before_comp
    assert df["GlyToucan ID"].tolist()[0].startswith("G") and df["WURCS"].tolist()[0].startswith("WURCS=2.0/")
    assert df["GlyToucan ID"].tolist()[1] == "USER-ID"                 # MAS value wins
    assert df["WURCS"].tolist()[1] == ""                               # obs. d: user accession != record -> WURCS left as-is
    assert df["GlyToucan ID"].tolist()[2] == "" and df["WURCS"].tolist()[2].startswith("WURCS")  # wurcs_only row
    assert df["GlyToucan ID"].tolist()[3] == "" and df["WURCS"].tolist()[3] == ""               # missing row
    assert df["GlyToucan ID"].tolist()[4] == "" and df["WURCS"].tolist()[5] == ""
    assert fr.status == "partial" and (fr.n_compositions, fr.n_filled, fr.n_missing) == (3, 2, 1)
    assert fr.n_cells_written == 3 and fr.n_user_accession == 1
    blk = fr.to_artifact_block()
    assert set(blk) == {"status", "route", "db", "n_compositions", "n_filled", "n_missing", "n_user_accession", "n_user_wurcs", "last_resolved_utc"}
    assert set(blk["db"]) == {"name", "version", "path", "sha256", "n_entries", "updated_utc"} and blk["route"] == "database"


def test_annotate_creates_columns_and_db_missing(tmpdb):
    db, _ = tmpdb
    df = pd.DataFrame({"composition": ["F1H5N4S2"]})
    fr = R.annotate_dataframe(df, "composition", None)
    assert fr.status == "db_missing" and list(df.columns) == ["composition", "GlyToucan ID", "WURCS"]
    assert df["GlyToucan ID"].tolist() == [""] and df["WURCS"].tolist() == [""]
    # empty database -> nothing filled, status 'partial' with n_filled 0
    fr2 = R.annotate_dataframe(df, "composition", db)
    assert fr2.status == "partial" and fr2.n_filled == 0 and fr2.n_missing == 1 and fr2.n_cells_written == 0
    # NaN cells count as empty
    df3 = pd.DataFrame({"Structure": ["Non-glycan"], "GlyToucan ID": [float("nan")], "WURCS": [None]})
    fr3 = R.annotate_dataframe(df3, "Structure", db)
    assert fr3.status == "none" and fr3.n_compositions == 0


# ---------------------------------------------------------------- CLI
def test_cli_no_network_diff_and_summary(tmp_path, capsys):
    import msp_CLI_glytoucan_prebuild as cli
    src = tmp_path / "src.tsv"
    pd.DataFrame({"composition": ["F1H5N4S2", "H5N4", "Non-glycan"], "x": [1, 2, 3]}).to_csv(src, sep="\t", index=False)
    dbp = tmp_path / "ref.db"
    rc = cli.main(["--source", str(src), "--db", str(dbp), "--log", str(tmp_path / "l.jsonl"), "--no-network", "--json"])
    out = json.loads(capsys.readouterr().out)
    assert rc == 0 and out["status"] == "disabled" and out["n_source"] == 2 and out["api_calls"] == 0
    rc = cli.main(["--summary", "--db", str(dbp), "--json"])
    out = json.loads(capsys.readouterr().out)
    assert rc == 0 and out["n_entries"] == 0 and out["db_version"] == R.DEFAULT_DB_VERSION


def test_compositions_from_file_autodetect(tmp_path):
    fixture = os.path.join(ROOT, "tests", "fixtures", "smoke_4", "zf_sPerMeNG_brain_CGA_20260522_180609.tsv")
    if os.path.exists(fixture):
        keys = R.compositions_from_file(fixture)
        assert keys and all(R.canonical_composition(k) == k for k in keys) and "F1H5N4S2" in keys
    xl = tmp_path / "mas.csv"
    pd.DataFrame({"Structure": ["H5N4F1S2", "Non-glycan", "H5N4F1S2"]}).to_csv(xl, index=False)
    assert R.compositions_from_file(xl) == ["F1H5N4S2"]
    pd.DataFrame({"a": [1]}).to_csv(tmp_path / "nocol.csv", index=False)
    with pytest.raises(ValueError):
        R.compositions_from_file(tmp_path / "nocol.csv")


# ---------------------------------------------------------------- H3 wiring (real mspval writer)
def test_h3_writer_enrich_fn_default_none_is_v110_and_fill_applies_once(tmp_path):
    """createnormailzedionlistcsv: enrich_fn=None -> identical output to the v1.10 call;
    enrich_fn given -> applied once, MAS sheet values win, Structure untouched, failure is soft."""
    import mspvalidator_merger as mspval
    ion_df = pd.DataFrame({"mass": [204.087, 366.140]})
    ionindex = mspval.extract_ionmasslist(ion_df)

    def frame():
        return pd.DataFrame([
            {"protonatedmass": 500.5, "peaklist": [204.087], "peakintensity": [100.0], "Structure": "F1H5N4S2",
             "IUPACname(optional)": "", "Glycanannotation2": "N", "GlyToucan ID": "", "WURCS": "", "MS2scan_no": 101},
            {"protonatedmass": 600.5, "peaklist": [366.140], "peakintensity": [50.0], "Structure": "H5N4F1S2",
             "IUPACname(optional)": "", "Glycanannotation2": "N", "GlyToucan ID": "SHEET", "WURCS": "", "MS2scan_no": 102},
            {"protonatedmass": 700.5, "peaklist": [], "peakintensity": [], "Structure": "Non-glycan",
             "IUPACname(optional)": "", "Glycanannotation2": "", "GlyToucan ID": "", "WURCS": "", "MS2scan_no": 103},
        ])

    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    mspval.createnormailzedionlistcsv(ionindex, frame(), ion_df, str(a))                    # v1.10 call
    mspval.createnormailzedionlistcsv(ionindex, frame(), ion_df, str(b), enrich_fn=None)    # explicit default
    assert a.read_bytes() == b.read_bytes()

    db = R.ReferenceDB(tmp_path / "ref.db")
    R.prebuild(["F1H5N4S2"], db, tmp_path / "l.jsonl", http_post=make_stub([]), sleep=lambda s: None)
    calls = []

    def enrich(df):
        calls.append(1)
        R.annotate_dataframe(df, "Structure", db)
        return df
    c = tmp_path / "c.csv"
    mspval.createnormailzedionlistcsv(ionindex, frame(), ion_df, str(c), enrich_fn=enrich)
    db.close()
    assert calls == [1]
    out = pd.read_csv(c, dtype=str).fillna("")
    assert out["Structure"].tolist() == ["F1H5N4S2", "H5N4F1S2", "Non-glycan"]        # never canonicalised/overwritten
    assert out["GlyToucan ID"].tolist()[0].startswith("G") and out["WURCS"].tolist()[0].startswith("WURCS")
    assert out["GlyToucan ID"].tolist()[1] == "SHEET" and out["WURCS"].tolist()[1] == ""   # obs. d rule
    assert out["GlyToucan ID"].tolist()[2] == "" and out["WURCS"].tolist()[2] == ""
    base = pd.read_csv(a, dtype=str).fillna("")
    assert list(out.columns) == list(base.columns)                                       # no new columns in MAS CSV

    def broken(df):
        raise RuntimeError("boom")
    d = tmp_path / "d.csv"
    mspval.createnormailzedionlistcsv(ionindex, frame(), ion_df, str(d), enrich_fn=broken)  # fail soft
    assert d.read_bytes() == a.read_bytes()


# ---------------------------------------------------------------- Codex S1 findings (2026-09-22)
def test_f3_native_insilico_library_count_schema(tmp_path):
    """F3: compnewv4.save_glycan_pseudocomp_to_csv writes residue counts, not a composition string."""
    hdr = "Hex,HexNAc,NeuAc,NeuGc,KDN,Fuc,HexA,Mode,SO3,PO3H,Mass\n"
    rows = ["3,2,0,0,0,0,0,+,0,0,1149.6011\n",
            "5,4,2,0,0,1,0,+,0,0,2792.4\n",     # F1H5N4S2
            "5,4,1,1,0,1,0,+,0,0,2800.0\n",     # F1H5N4S1G1
            "5,4,0,0,1,1,0,+,0,0,2700.0\n",     # F1H5N4KDN1
            "3,2,0,0,0,0,1,+,0,0,1300.0\n",     # H3N2A1
            "3,2,0,0,0,0,0,-,1,0,1230.0\n",     # H3N2s1
            "5,2,0,0,0,0,0,+,0,1,1600.0\n",     # H5N2p1
            "5,4,2,0,0,1,0,+,0,0,2792.4\n"]     # duplicate row -> one key
    lib = tmp_path / "lib_insilico.csv"
    lib.write_text(hdr + "".join(rows))
    keys = R.compositions_from_file(lib)
    assert keys == sorted(["H3N2", "F1H5N4S2", "F1H5N4S1G1", "F1H5N4KDN1", "H3N2A1", "H3N2s1", "H5N2p1"])
    # every key round-trips through the API mapping (incl. KDN/HexA/sulfate/phosphate)
    assert R.composition_to_api_query("F1H5N4KDN1")["kdn"] == 1 and R.composition_to_api_query("H3N2s1")["S"] == 1
    # no-network diff on the native library works end to end via the CLI
    import msp_CLI_glytoucan_prebuild as cli
    rc = cli.main(["--source", str(lib), "--db", str(tmp_path / "r.db"), "--log", str(tmp_path / "l.jsonl"), "--no-network", "--json"])
    assert rc == 0
    # a string column still wins when present; a frame with neither shape still raises
    with pytest.raises(ValueError):
        R.compositions_from_count_frame(pd.DataFrame({"Hex": [1]}))


def test_f4_malformed_item_is_soft_and_earlier_batches_persist(tmpdb):
    """F4: a non-string wurcs/id in one item must not abort the run nor lose earlier batches."""
    db, log = tmpdb
    calls = []
    good = make_stub(calls)

    def stub(url, body, timeout):
        raw = good(url, body, timeout)
        if len(calls) == 2:  # corrupt one item of the second batch only
            d = json.loads(raw)
            d["data"][0]["wurcs"] = 123
            d["data"][1]["id"] = {"unexpected": "object"}
            raw = json.dumps(d).encode()
        if len(calls) == 3:
            raise OSError("dropped after two batches")
        return raw
    keys = [f"H{i}N2" for i in range(1, 71)]   # 3 batches: 30 / 30 / 10
    r = R.prebuild(keys, db, log, http_post=stub, sleep=lambda s: None)
    assert r.status == "partial" and r.api_calls == 2
    assert r.n_fetched == 58 and db.count() == 58            # 30 + 28 persisted before the third batch failed
    assert sum(1 for _ in open(log, encoding="utf-8")) == 58
    b2 = sorted(keys)[30:32]                                  # prebuild sorts keys lexicographically
    assert len(r.unsent) == 12 and all(k in r.unsent for k in b2)
    assert any("malformed response item" in e for e in r.errors) and any("dropped" in e for e in r.errors)
    assert all(db.get(k) is None for k in b2)                 # malformed items stay retryable
    # from_api rejects bad types / non-WURCS strings outright
    with pytest.raises(TypeError):
        R.Entry.from_api("H1N1", {"hex": 1, "hexnac": 1}, 123, "")
    with pytest.raises(TypeError):
        R.Entry.from_api("H1N1", {"hex": 1, "hexnac": 1}, "not a wurcs", "G1")


def test_f5_runtime_readonly_never_writes(tmp_path):
    """F5: run-time fill opens mode=ro; bytes unchanged; empty/foreign/missing files -> db_missing."""
    import hashlib
    dbp = tmp_path / "ref.db"
    w = R.ReferenceDB(dbp)
    R.prebuild(["F1H5N4S2"], w, tmp_path / "l.jsonl", http_post=make_stub([]), sleep=lambda s: None)
    w.close()
    before = hashlib.sha256(dbp.read_bytes()).hexdigest()
    os.chmod(dbp, 0o444)
    try:
        ro = R.ReferenceDB(dbp, readonly=True)
        df = pd.DataFrame({"composition": ["F1H5N4S2", "H9N9"]})
        fr = R.annotate_dataframe(df, "composition", ro)
        assert fr.status == "partial" and df["WURCS"].tolist()[0].startswith("WURCS")
        with pytest.raises(Exception):
            ro.upsert([R.Entry.from_api("H9N9", {"hex": 9, "hexnac": 9}, "WURCS=2.0/x", "")])
        ro.close()
    finally:
        os.chmod(dbp, 0o644)
    assert hashlib.sha256(dbp.read_bytes()).hexdigest() == before
    # empty sqlite file: not initialised, reported as db_missing
    empty = tmp_path / "empty.db"
    empty.write_bytes(b"")
    fr = R.annotate_dataframe(pd.DataFrame({"composition": ["F1H5N4S2"]}), "composition", R.ReferenceDB(empty, readonly=True))
    assert fr.status == "db_missing" and empty.read_bytes() == b""
    # missing file: not created
    missing = tmp_path / "nope" / "ref.db"
    fr = R.annotate_dataframe(pd.DataFrame({"composition": ["F1H5N4S2"]}), "composition", R.ReferenceDB(missing, readonly=True))
    assert fr.status == "db_missing" and not missing.parent.exists()


def test_writer_default_path_equals_v110_baseline_from_git(tmp_path):
    """Codex note: compare the new default path against the BASELINE merger (git f27a872), not itself."""
    import subprocess, importlib.util
    try:
        src = subprocess.check_output(["git", "-C", ROOT, "show", "f27a872:src/mspvalidator_merger.py"])
    except Exception as ex:  # no git / no baseline in this checkout
        pytest.skip("v1.10 baseline commit f27a872 is not in this checkout's git history (expected in a public "
                    "clone); this check runs only in the development repository")
    base_path = tmp_path / "mspvalidator_merger_baseline.py"
    base_path.write_bytes(src)
    spec = importlib.util.spec_from_file_location("mspval_baseline", base_path)
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    import mspvalidator_merger as new
    ion_df = pd.DataFrame({"mass": [204.087, 366.140]})
    ionindex = new.extract_ionmasslist(ion_df)
    rows = [{"protonatedmass": 500.5, "peaklist": [204.087], "peakintensity": [100.0], "Structure": "F1H5N4S2",
             "IUPACname(optional)": "", "Glycanannotation2": "N", "GlyToucan ID": "G1", "WURCS": "WURCS=2.0/x", "MS2scan_no": 101},
            {"protonatedmass": 700.5, "peaklist": [], "peakintensity": [], "Structure": "Non-glycan",
             "IUPACname(optional)": "", "Glycanannotation2": "", "GlyToucan ID": "", "WURCS": "", "MS2scan_no": 103}]
    a, b = tmp_path / "base.csv", tmp_path / "new.csv"
    base.createnormailzedionlistcsv(ionindex, pd.DataFrame(rows), ion_df, str(a))
    new.createnormailzedionlistcsv(ionindex, pd.DataFrame(rows), ion_df, str(b))
    assert a.read_bytes() == b.read_bytes()


def test_mas_workbook_reads_mslist_sheet(tmp_path):
    """MAS workbooks carry annotations on 'MSlist', not the first sheet."""
    xl = tmp_path / "mas.xlsx"
    with pd.ExcelWriter(xl) as w:
        pd.DataFrame({"Glycan name/abbr": ["x"], "Composition": ["H9N9"]}).to_excel(w, sheet_name="Glycanlist", index=False)
        pd.DataFrame({"Structure": ["H5N4F1S2", "Non-glycan", "H3N2"], "MS2scan_no": [1, 2, 3]}).to_excel(w, sheet_name="MSlist", index=False)
    assert R.compositions_from_file(xl) == ["F1H5N4S2", "H3N2"]
    assert R.compositions_from_file(xl, sheet="Glycanlist") == ["H9N9"]


# ---------------------------------------------------------------- Henry manual test 2026-09-22, obs. c + d
def test_obs_d_user_accession_keeps_wurcs_untouched(tmpdb):
    """A row with a user accession that differs from the composition record must not receive the
    generic composition-level WURCS (run 9 shape: two H3N6S2 rows with different accessions + blanks)."""
    db, log = tmpdb
    R.prebuild(["H3N6S2"], db, log, http_post=make_stub([]), sleep=lambda s: None)
    generic = db.get("H3N6S2")
    df = pd.DataFrame({
        "Structure": ["H3N6S2", "H3N6S2", "H3N6S2", "H3N6S2"],
        "GlyToucan ID": ["G73319AZ", "G92654QI", "", generic.accession],
        "WURCS": ["", "WURCS=2.0/user-specific", "", ""],
    })
    fr = R.annotate_dataframe(df, "Structure", db)
    assert df["GlyToucan ID"].tolist() == ["G73319AZ", "G92654QI", generic.accession, generic.accession]
    assert df["WURCS"].tolist()[0] == ""                                  # user accession, blank WURCS stays blank
    assert df["WURCS"].tolist()[1] == "WURCS=2.0/user-specific"           # user accession + user WURCS untouched
    assert df["WURCS"].tolist()[2] == generic.wurcs                        # blank row: generic accession + WURCS
    assert df["WURCS"].tolist()[3] == generic.wurcs                        # same accession as record: WURCS filled
    assert fr.n_user_accession == 2 and fr.n_cells_written == 3
    assert fr.to_artifact_block()["n_user_accession"] == 2


def test_obs_c_updated_utc_and_n_entries_in_block(tmpdb):
    db, log = tmpdb
    ident0 = db.identity()
    assert ident0["updated_utc"] == db.meta()["created_utc"] and ident0["n_entries"] == 0
    R.prebuild(["F1H5N4S2"], db, log, http_post=make_stub([]), sleep=lambda s: None)
    m = db.meta()
    assert "updated_utc" in m and m["updated_utc"] >= m["created_utc"]
    fr = R.annotate_dataframe(pd.DataFrame({"composition": ["F1H5N4S2"]}), "composition", db)
    blk = fr.to_artifact_block()
    assert set(blk["db"]) == {"name", "version", "path", "sha256", "n_entries", "updated_utc"}
    assert blk["db"]["n_entries"] == 1 and blk["db"]["updated_utc"] == m["updated_utc"]


# ---------------------------------------------------------------- Codex round 2 (2026-09-23): R2-F1..F4 + D1
class _CommitFailsOnce:
    """Connection proxy whose first commit() raises, to simulate a lock/disk failure mid-transaction."""
    def __init__(self, conn):
        self._c = conn
        self.fail = True
    def commit(self):
        if self.fail:
            self.fail = False
            raise sqlite3.OperationalError("database is locked (simulated)")
        return self._c.commit()
    def __getattr__(self, name):
        return getattr(self._c, name)


import sqlite3  # noqa: E402


def test_r2f1_upsert_is_transactional_and_rolls_back(tmpdb):
    db, _ = tmpdb
    db._connect()
    real = db._conn
    db._conn = _CommitFailsOnce(real)
    entries = [R.Entry.from_api(f"H{i}N2", {"hex": i, "hexnac": 2}, "WURCS=2.0/x", "") for i in range(1, 6)]
    with pytest.raises(sqlite3.OperationalError):
        db.upsert(entries)
    assert not real.in_transaction                      # rolled back, nothing left open
    assert db.count() == 0 and "updated_utc" not in db.meta()
    # a later inspection must not make anything durable
    db.identity()
    assert db.count() == 0
    # next upsert works and is complete
    assert db.upsert(entries) == 5 and db.count() == 5 and "updated_utc" in db.meta()


def test_r2f1_persist_failure_propagates_and_stops(tmpdb):
    db, log = tmpdb
    calls = []
    good = make_stub(calls)
    orig_upsert = db.upsert
    def flaky_upsert(entries):
        if len(calls) == 2:
            raise sqlite3.OperationalError("disk I/O error (simulated)")
        return orig_upsert(entries)
    db.upsert = flaky_upsert
    keys = [f"H{i}N2" for i in range(1, 71)]            # 30 / 30 / 10
    r = R.prebuild(keys, db, log, http_post=good, sleep=lambda s: None)
    assert r.status == "persist_failed" and not r.storage_ok
    assert r.api_calls == 2                              # stopped after the failing batch
    assert r.n_fetched == 30 and db.count() == 30        # only committed rows are counted
    assert len(r.unpersisted) == 30 and len(r.unsent) == 10
    assert sum(1 for _ in open(log, encoding="utf-8")) == 30
    assert any("persistence failed" in e for e in r.errors)
    # retry: the unpersisted + unsent keys are still missing, so the next run fetches exactly them
    db.upsert = orig_upsert
    r2 = R.prebuild(keys, db, log, http_post=good, sleep=lambda s: None)
    assert r2.status == "ok" and r2.n_already == 30 and r2.n_fetched == 40 and db.count() == 70


def test_r2f2_log_failure_is_reported_and_rebuildable(tmp_path):
    dbp = tmp_path / "ref.db"
    db = R.ReferenceDB(dbp)
    bad_log = tmp_path / "logdir"
    bad_log.mkdir()                                       # a directory: open(..., "a") fails
    r = R.prebuild([f"H{i}N2" for i in range(1, 4)], db, bad_log, http_post=make_stub([]), sleep=lambda s: None)
    assert r.status == "log_failed" and not r.storage_ok
    assert r.n_fetched == 3 and db.count() == 3           # database committed regardless
    assert sorted(r.unlogged) == ["H1N2", "H2N2", "H3N2"] and any("could not append" in e for e in r.errors)
    # a second run finds nothing to fetch (DB is source of truth) — the log gap is NOT silently repaired ...
    r2 = R.prebuild([f"H{i}N2" for i in range(1, 4)], db, bad_log, http_post=make_stub([]), sleep=lambda s: None)
    assert r2.status == "nothing_to_do"
    # ... rebuild_log recovers every committed transaction as valid JSONL
    good_log = tmp_path / "log.jsonl"
    good_log.write_text('{"old": true}\n')
    n = R.rebuild_log(db, good_log)
    lines = [json.loads(l) for l in good_log.read_text(encoding="utf-8").splitlines()]
    assert n == 3 and len(lines) == 3 and all(l["rebuilt"] and l["op"] == "api_resolve" for l in lines)
    assert {l["composition_key"] for l in lines} == {"H1N2", "H2N2", "H3N2"}
    assert (tmp_path / "log.jsonl.bak").read_text() == '{"old": true}\n'
    # append_log itself reports the partial count
    with pytest.raises(R.LogWriteError) as ei:
        R.append_log(bad_log, [{"a": 1}])
    assert ei.value.written == 0
    db.close()


def test_r2f2_cli_exit_codes_and_rebuild_log(tmp_path, capsys):
    import msp_CLI_glytoucan_prebuild as cli
    src = tmp_path / "src.csv"
    pd.DataFrame({"composition": ["F1H5N4S2"]}).to_csv(src, index=False)
    dbp = tmp_path / "ref.db"
    bad_log = tmp_path / "logdir"; bad_log.mkdir()
    # network path stubbed at module level
    R_http = R._default_http_post
    R._default_http_post = make_stub([])
    try:
        rc = cli.main(["--source", str(src), "--db", str(dbp), "--log", str(bad_log), "--pause", "0", "--json"])
        out = json.loads(capsys.readouterr().out)
        assert rc == 2 and out["status"] == "log_failed" and out["n_unlogged"] == 1
        rc = cli.main(["--rebuild-log", "--db", str(dbp), "--log", str(tmp_path / "l.jsonl"), "--json"])
        out = json.loads(capsys.readouterr().out)
        assert rc == 0 and out["rebuilt_lines"] == 1
        rc = cli.main(["--source", str(src), "--db", str(dbp), "--log", str(tmp_path / "l.jsonl"), "--json"])
        assert rc == 0 and json.loads(capsys.readouterr().out)["status"] == "nothing_to_do"
    finally:
        R._default_http_post = R_http


@pytest.mark.parametrize("name", ["ref%20db.db", "ref db.db", "ref#1.db", "ref?q.db", "réf_ü.db", "ref%2Fx.db"])
def test_r2f4_readonly_uri_addresses_the_same_file(tmp_path, name):
    if os.name == "nt" and any(ch in name for ch in "?"):
        pytest.skip("'?' is not a legal Windows filename")
    dbp = tmp_path / name
    w = R.ReferenceDB(dbp)
    R.prebuild(["F1H5N4S2"], w, tmp_path / "l.jsonl", http_post=make_stub([]), sleep=lambda s: None)
    w.close()
    # decoy with the percent-decoded name must NOT be what read-only opens
    decoy = tmp_path / name.replace("%20", " ").replace("%2F", "_")
    if decoy != dbp:
        d = R.ReferenceDB(decoy)
        R.prebuild(["H9N9"], d, tmp_path / "l2.jsonl", http_post=make_stub([]), sleep=lambda s: None)
        d.close()
    ro = R.ReferenceDB(dbp, readonly=True)
    assert ro.get("F1H5N4S2") is not None and ro.get("H9N9") is None
    ident = ro.identity()
    assert ident["path"] == str(dbp) and ident["n_entries"] == 1
    import hashlib
    assert ident["sha256"] == hashlib.sha256(dbp.read_bytes()).hexdigest()
    ro.close()


def test_d1_user_wurcs_with_blank_accession_keeps_accession_blank(tmpdb):
    db, log = tmpdb
    R.prebuild(["H3N6S1G1", "F1H5N4S2"], db, log, http_post=make_stub([]), sleep=lambda s: None)
    rec = db.get("F1H5N4S2")
    df = pd.DataFrame({
        "Structure":    ["H3N6S1G1",                 "F1H5N4S2", "F1H5N4S2",    "F1H5N4S2"],
        "GlyToucan ID": ["",                          "",         rec.accession, ""],
        "WURCS":        ["WURCS=2.0/from-GlycanBuilder", "",     "WURCS=2.0/user-different", rec.wurcs],
    })
    fr = R.annotate_dataframe(df, "Structure", db)
    assert df["GlyToucan ID"].tolist()[0] == "" and df["WURCS"].tolist()[0] == "WURCS=2.0/from-GlycanBuilder"  # D1
    assert df["GlyToucan ID"].tolist()[1] == rec.accession and df["WURCS"].tolist()[1] == rec.wurcs               # blank row
    assert df["GlyToucan ID"].tolist()[2] == rec.accession and df["WURCS"].tolist()[2] == "WURCS=2.0/user-different"  # preserved pre-existing pair
    assert df["GlyToucan ID"].tolist()[3] == rec.accession and df["WURCS"].tolist()[3] == rec.wurcs               # WURCS equals record -> accession filled
    assert fr.n_user_wurcs == 1 and fr.n_user_accession == 0 and fr.n_cells_written == 3
    assert fr.to_artifact_block()["n_user_wurcs"] == 1


def test_r2f3_static_ordering_publish_after_write():
    """R2-F3 (static): in v14, provenance is recorded only after the output write in both H1 and H3."""
    src = open(os.path.join(ROOT, "src", "mspfileloaderv14.py"), "rb").read().decode("utf-8")
    h1 = src.index("def run_pseudolabeling(")
    h1_write = src.index('out.to_csv(outpath, index=False, sep="\\t"', h1)
    h1_record = src.index('_glytoucan_record(files, "CGA", _gt_block)', h1)
    assert h1_record > h1_write
    h3 = src.index("def try_merge_selected_sample(")
    h3_write = src.index("mspval.createnormailzedionlistcsv(iondfindex, pre_df,ion_df, outpath, enrich_fn=_gt_enrich)", h3)
    h3_record = src.index('_glytoucan_record(files, "MAS", _gt_holder.get("block"))', h3)
    assert h3_record > h3_write
    # and no early record/clear remains inside the enrich closure or before the writer
    assert '_glytoucan_record(_files, "MAS"' not in src and '_glytoucan_record(files, "MAS", None)' not in src


# ---------------------------------------------------------------- Codex round 3 (2026-09-23): R3-F1, F2, F5, F6
def _seed_db(tmp_path, name="ref.db", keys=("F1H5N4S2", "H3N2")):
    dbp = tmp_path / name
    db = R.ReferenceDB(dbp)
    R.prebuild(list(keys), db, tmp_path / "seed.jsonl", http_post=make_stub([]), sleep=lambda s: None)
    db.close()
    return dbp


def test_r3f1_rebuild_refuses_missing_and_foreign_sources(tmp_path):
    live = tmp_path / "log.jsonl"; live.write_text('{"keep": 1}\n')
    bak = tmp_path / "log.jsonl.bak"; bak.write_text('{"old": 1}\n')
    # missing database: refused, nothing created, log + bak untouched
    missing = tmp_path / "nope" / "ref.db"
    with pytest.raises(R.RebuildLogError):
        R.rebuild_log(R.ReferenceDB(missing, readonly=True), live)
    with pytest.raises(R.RebuildLogError):
        R.rebuild_log(R.ReferenceDB(missing), live)          # writable handle: still refused, still not created
    assert not missing.parent.exists() and live.read_text() == '{"keep": 1}\n' and bak.read_text() == '{"old": 1}\n'
    # empty / foreign sqlite file: refused, bytes untouched
    for name, payload in (("empty.db", b""), ("foreign.db", None)):
        f = tmp_path / name
        if payload is None:
            c = sqlite3.connect(f); c.execute("CREATE TABLE other(x)"); c.commit(); c.close()
        else:
            f.write_bytes(payload)
        before = f.read_bytes()
        with pytest.raises(R.RebuildLogError):
            R.rebuild_log(R.ReferenceDB(f, readonly=True), live)
        assert f.read_bytes() == before and live.read_text() == '{"keep": 1}\n'
    # valid database still works, and a deliberately EMPTY valid reference database is allowed (0 lines)
    good = _seed_db(tmp_path)
    assert R.rebuild_log(R.ReferenceDB(good, readonly=True), live) == 2
    emptydb = tmp_path / "valid_empty.db"; _d = R.ReferenceDB(emptydb); _d.count(); _d.close()   # connect lazily -> file exists
    assert R.rebuild_log(R.ReferenceDB(emptydb, readonly=True), live) == 0


def test_r3f1_cli_rebuild_never_creates_database(tmp_path, capsys):
    import msp_CLI_glytoucan_prebuild as cli
    live = tmp_path / "log.jsonl"; live.write_text('{"keep": 1}\n')
    missing = tmp_path / "ghost.db"
    rc = cli.main(["--rebuild-log", "--db", str(missing), "--log", str(live), "--json"])
    out = json.loads(capsys.readouterr().out)
    assert rc == 2 and out["status"] == "rebuild_failed" and not missing.exists()
    assert live.read_text() == '{"keep": 1}\n' and not (tmp_path / "log.jsonl.bak").exists()


def test_r3f2_backup_by_copy_abort_on_obstruction_and_failed_replace(tmp_path, monkeypatch):
    good = _seed_db(tmp_path)
    ro = R.ReferenceDB(good, readonly=True)
    live = tmp_path / "log.jsonl"; live.write_text('{"keep": 1}\n')
    # obstructed backup destination (a directory): abort, live untouched, no tmp left
    bakdir = tmp_path / "log.jsonl.bak"; bakdir.mkdir()
    with pytest.raises(R.RebuildLogError):
        R.rebuild_log(ro, live)
    assert live.read_text() == '{"keep": 1}\n' and bakdir.is_dir() and not (tmp_path / "log.jsonl.tmp").exists()
    bakdir.rmdir()
    # failure at the final replace: live untouched, backup present, tmp removed
    real_replace = os.replace
    def failing_replace(src, dst):
        if str(dst) == str(live):
            raise PermissionError("simulated replace failure")
        return real_replace(src, dst)
    monkeypatch.setattr(R.os, "replace", failing_replace)
    with pytest.raises(R.RebuildLogError):
        R.rebuild_log(ro, live)
    monkeypatch.undo()
    assert live.read_text() == '{"keep": 1}\n'
    assert (tmp_path / "log.jsonl.bak").read_text() == '{"keep": 1}\n' and not (tmp_path / "log.jsonl.tmp").exists()
    # success: live is the new projection, .bak holds the previous content, live path never absent
    n = R.rebuild_log(ro, live)
    lines = [json.loads(l) for l in live.read_text().splitlines()]
    assert n == 2 and len(lines) == 2 and (tmp_path / "log.jsonl.bak").read_text() == '{"keep": 1}\n'
    ro.close()


def test_r3f5_identity_failure_is_soft_on_every_return_path(tmpdb, monkeypatch):
    db, log = tmpdb
    R.prebuild(["F1H5N4S2"], db, log, http_post=make_stub([]), sleep=lambda s: None)
    def boom():
        raise PermissionError("hash file unreadable (simulated)")
    monkeypatch.setattr(db, "identity", boom)
    r1 = R.prebuild(["F1H5N4S2"], db, log, http_post=make_stub([]), sleep=lambda s: None)
    assert r1.status == "nothing_to_do" and r1.n_already == 1 and r1.db == {"path": str(db.path)} and any("identity() failed" in e for e in r1.errors)
    r2 = R.prebuild(["F1H5N4S2", "H9N9"], db, log, allow_network=False)
    assert r2.status == "disabled" and r2.unsent == ["H9N9"] and r2.db == {"path": str(db.path)} and any("identity() failed" in e for e in r2.errors)
    r3 = R.prebuild(["H8N8"], db, log, http_post=make_stub([]), sleep=lambda s: None)
    assert r3.status == "ok" and r3.n_fetched == 1 and r3.db == {"path": str(db.path)}


def test_r3f6_help_shows_exit_codes_and_log_failed_with_unsent(tmp_path, capsys):
    import msp_CLI_glytoucan_prebuild as cli
    with pytest.raises(SystemExit) as ei:
        cli.main(["--help"])
    assert ei.value.code == 0
    help_text = capsys.readouterr().out
    assert "exit codes" in help_text and "--rebuild-log" in help_text and "n_unsent" in help_text
    # committed first batch + log failure + later network failure -> log_failed AND unsent > 0
    dbp = tmp_path / "ref.db"; db = R.ReferenceDB(dbp)
    baddir = tmp_path / "logdir"; baddir.mkdir()
    calls = []; good = make_stub(calls)
    def flaky(url, body, timeout):
        if len(calls) >= 1:
            raise OSError("network dropped")
        return good(url, body, timeout)
    r = R.prebuild([f"H{i}N2" for i in range(1, 41)], db, baddir, http_post=flaky, sleep=lambda s: None)
    assert r.status == "log_failed" and not r.storage_ok
    assert r.n_fetched == 30 and len(r.unlogged) == 30 and len(r.unsent) == 10
    db.close()


# ---------------------------------------------------------------- Codex round 4 (2026-09-23): R4-F1, R4-F2
def _sqlite_sig(path):
    """(bytes hash, table names) so a test can prove a source was neither rewritten nor initialised."""
    import hashlib
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    try:
        c = sqlite3.connect(f"{path.resolve().as_uri()}?mode=ro", uri=True)
        tables = sorted(r[0] for r in c.execute("SELECT name FROM sqlite_master WHERE type='table'"))
        c.close()
    except Exception:
        tables = None
    return h, tables


@pytest.mark.parametrize("handle", ["default", "create_false", "readonly"])
@pytest.mark.parametrize("kind", ["empty", "foreign"])
def test_r4f1_library_recovery_never_initialises_invalid_source(tmp_path, handle, kind):
    f = tmp_path / f"{kind}.db"
    if kind == "empty":
        f.write_bytes(b"")
    else:
        c = sqlite3.connect(f); c.execute("CREATE TABLE other(x)"); c.commit(); c.close()
    live = tmp_path / "log.jsonl"; live.write_text('{"keep": 1}\n')
    before = _sqlite_sig(f)
    db = {"default": lambda: R.ReferenceDB(f),
          "create_false": lambda: R.ReferenceDB(f, create=False),
          "readonly": lambda: R.ReferenceDB(f, readonly=True)}[handle]()
    with pytest.raises(R.RebuildLogError):
        R.rebuild_log(db, live)
    assert _sqlite_sig(f) == before                       # not rewritten, no reference tables added
    assert live.read_text() == '{"keep": 1}\n' and not (tmp_path / "log.jsonl.bak").exists()
    assert not list(tmp_path.glob("log.jsonl.*.tmp"))
    # the caller's writable handle was never connected by recovery
    assert db._conn is None
    db.close()


def test_r4f1_valid_sources_still_work_through_any_handle(tmp_path):
    good = _seed_db(tmp_path)
    live = tmp_path / "log.jsonl"
    for db in (R.ReferenceDB(good), R.ReferenceDB(good, create=False), R.ReferenceDB(good, readonly=True)):
        assert R.rebuild_log(db, live) == 2
        db.close()
    emptydb = tmp_path / "valid_empty.db"; d = R.ReferenceDB(emptydb); d.count(); d.close()
    assert R.rebuild_log(R.ReferenceDB(emptydb), live) == 0


def test_r4f2_destination_aliasing_is_refused_before_mutation(tmp_path):
    good = _seed_db(tmp_path)
    sig = _sqlite_sig(good)
    # source == live log
    with pytest.raises(R.RebuildLogError):
        R.rebuild_log(R.ReferenceDB(good, readonly=True), good)
    assert _sqlite_sig(good) == sig and not (tmp_path / "ref.db.bak").exists()
    # source == backup path of a live log
    live2 = tmp_path / "other.jsonl"; live2.write_text('{"keep": 2}\n')
    src_at_bak = tmp_path / "other.jsonl.bak"
    import shutil; shutil.copy2(good, src_at_bak)
    sig2 = _sqlite_sig(src_at_bak)
    with pytest.raises(R.RebuildLogError):
        R.rebuild_log(R.ReferenceDB(src_at_bak, readonly=True), live2)
    assert _sqlite_sig(src_at_bak) == sig2 and live2.read_text() == '{"keep": 2}\n'
    # symlinked backup pointing at the source
    live3 = tmp_path / "third.jsonl"; live3.write_text('{"keep": 3}\n')
    link = tmp_path / "third.jsonl.bak"
    try:
        link.symlink_to(good)
    except (OSError, NotImplementedError):
        pytest.skip("symlinks not available")
    with pytest.raises(R.RebuildLogError):
        R.rebuild_log(R.ReferenceDB(good, readonly=True), live3)
    assert _sqlite_sig(good) == sig and live3.read_text() == '{"keep": 3}\n'
    link.unlink()
    # symlinked live log is refused too
    live4 = tmp_path / "fourth.jsonl"; live4.symlink_to(good)
    with pytest.raises(R.RebuildLogError):
        R.rebuild_log(R.ReferenceDB(good, readonly=True), live4)
    assert _sqlite_sig(good) == sig
    # hard link between live log and backup
    live5 = tmp_path / "fifth.jsonl"; live5.write_text('{"keep": 5}\n')
    try:
        os.link(live5, tmp_path / "fifth.jsonl.bak")
    except (OSError, NotImplementedError):
        pytest.skip("hard links not available")
    with pytest.raises(R.RebuildLogError):
        R.rebuild_log(R.ReferenceDB(good, readonly=True), live5)
    assert live5.read_text() == '{"keep": 5}\n'
    # temp files are unique and cleaned up; ordinary rebuild still fine
    live6 = tmp_path / "sixth.jsonl"; live6.write_text('{"keep": 6}\n')
    assert R.rebuild_log(R.ReferenceDB(good, readonly=True), live6) == 2
    assert not list(tmp_path.glob("sixth.jsonl.*.tmp")) and (tmp_path / "sixth.jsonl.bak").read_text() == '{"keep": 6}\n'


# ---------------------------------------------------------------- H5 (2026-09-23): CGA -> Trainable fill
def test_h5_static_fill_before_write_and_record_after_in_caller():
    src = open(os.path.join(ROOT, "src", "mspfileloaderv14.py"), "rb").read().decode("utf-8")
    b = src.index("def build_trainable_from_CGA(")
    fill = src.index('_gt_block = _glytoucan_fill_df(final_df, "Structure")', b)
    write = src.index("final_df.to_csv(output_path, index=False)", b)
    assert fill < write                                   # fill happens on the frame that is written
    assert '"glytoucan_resolution": _gt_block' in src[b:src.index("end CGA (PSEUDOLABEL)", b)]
    call = src.index("outpath, summary = build_trainable_from_CGA(")
    rec = src.index('_glytoucan_record(files, "CGA", summary["glytoucan_resolution"])', call)
    setp = src.index('files["trainable_csv"] = outpath', call)
    assert rec > setp                                     # recorded only after the CSV exists and is linked
    assert 'is not None' in src[rec-200:rec]              # option off -> no record, no clear
    linkc = src.index("_cga_cache_link_trainable(files, outpath)", call)
    assert setp < linkc < rec                             # R5-F1: cache path synced before the block is recorded
