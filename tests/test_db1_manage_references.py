"""DB1 (v1.12, T1-b "Manage References") — headless acceptance checks.

Spec: handoff/T1b_manage_references_spec_20260929.md §10 (headless items). No network, no Tk:
the HTTP layer is stubbed via `http_post`, the subprocess runner is driven through its queue.

Sections (added checkpoint by checkpoint):
  - CP-A  resolver: last_prebuild_* meta on every prebuild return path (writable handle only),
          summary() exposure, CLI source_label, G1 snapshot-consistent identity()          (§10 item 5)
  - CP-B  v14 prologue: settings.json load/save rules, session db_path honoured by the fill (§10 items 1-2)
  - CP-C  _GtCliRunner: exact argv, queue-driven runner against tests/fixtures/db1_stub_cli.py, exit/status mapping,
          terminate -> kill, real-CLI status line matched by the regex                        (§10 items 3-4)
  - CP-D  GUI: popup text builder, active-DB label, manager summary (read-only), shared fill flag, static audit of
          every new command=/bind/after/protocol target, rows 0-3 + menu order, optional real-Tk smoke (§10 items 6-7)

Run from repo root:
    python3 -m pytest tests/test_db1_manage_references.py -q
"""
import hashlib
import json
import os
import sqlite3
import sys

import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import msp_glytoucan_resolver as R  # noqa: E402


# ---------------------------------------------------------------- helpers
def make_stub(calls=None):
    """HTTP stub: echoes non-zero counts + a stub wurcs; accession only when neu5ac > 0 (as in test_s1)."""
    calls = [] if calls is None else calls

    def stub(url, body, timeout):
        items = json.loads(body.decode("utf-8"))
        calls.append(len(items))
        out = []
        for it in items:
            echo = {k: v for k, v in it.items() if v != "0"}
            echo["wurcs"] = "WURCS=2.0/stub/" + "".join(f"{k}{v}" for k, v in sorted(echo.items()))
            echo["id"] = "G" + str(len(out)).zfill(6) if it.get("neu5ac", "0") != "0" else ""
            out.append(echo)
        return json.dumps({"data": out}).encode("utf-8")
    return stub


def _seed_db(tmp_path, name="ref.db", keys=("F1H5N4S2", "H3N2")):
    dbp = tmp_path / name
    db = R.ReferenceDB(dbp)
    R.prebuild(list(keys), db, tmp_path / "seed.jsonl", http_post=make_stub(), sleep=lambda s: None)
    db.close()
    return dbp


def _file_sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def _raw_meta(dbp):
    c = sqlite3.connect(str(dbp))
    try:
        return dict(c.execute("SELECT key, value FROM meta").fetchall())
    finally:
        c.close()


# ================================================================ CP-A — resolver (§7)
def test_cpa_last_prebuild_meta_on_every_return_path(tmp_path):
    dbp = tmp_path / "ref.db"
    log = tmp_path / "log.jsonl"
    db = R.ReferenceDB(dbp)
    try:
        # normal path
        r = R.prebuild(["F1H5N4S2", "H5N4"], db, log, http_post=make_stub(), sleep=lambda s: None,
                       source_label="/some/dir/brain_CGA.tsv")
        m = db.meta()
        assert r.status == "ok"
        assert m["last_prebuild_source"] == "brain_CGA.tsv"          # basename only
        assert m["last_prebuild_status"] == "ok" and m["last_prebuild_utc"].endswith("Z")
        # nothing_to_do path
        r = R.prebuild(["H5N4"], db, log, http_post=make_stub(), sleep=lambda s: None, source_label="second.csv")
        m = db.meta()
        assert r.status == "nothing_to_do"
        assert (m["last_prebuild_source"], m["last_prebuild_status"]) == ("second.csv", "nothing_to_do")
        # disabled (--no-network) path
        r = R.prebuild(["H9N9"], db, log, allow_network=False, source_label="third.xlsx")
        m = db.meta()
        assert r.status == "disabled"
        assert (m["last_prebuild_source"], m["last_prebuild_status"]) == ("third.xlsx", "disabled")
        # network fail-soft path is a normal return too
        def down(url, body, timeout):
            raise OSError("no network")
        r = R.prebuild(["H8N8"], db, log, http_post=down, sleep=lambda s: None, source_label="fourth.tsv")
        assert r.status == "offline" and db.meta()["last_prebuild_status"] == "offline"
        # no label given -> empty source, still recorded
        r = R.prebuild(["H5N4"], db, log)
        assert db.meta()["last_prebuild_source"] == "" and db.meta()["last_prebuild_status"] == "nothing_to_do"
    finally:
        db.close()


def test_cpa_reported_hash_matches_file_after_meta_write(tmp_path):
    """The meta write happens before identity, so res.db describes the file as prebuild leaves it."""
    dbp = _seed_db(tmp_path)
    db = R.ReferenceDB(dbp)
    try:
        r = R.prebuild(["H3N2"], db, tmp_path / "l.jsonl", source_label="x.tsv")
        assert r.status == "nothing_to_do"
        assert r.db["sha256"] == _file_sha(dbp)
    finally:
        db.close()


def test_cpa_persist_failed_status_is_recorded(tmp_path, monkeypatch):
    dbp = tmp_path / "ref.db"
    db = R.ReferenceDB(dbp)
    try:
        def broken_upsert(entries):
            raise sqlite3.OperationalError("disk I/O error (simulated)")
        monkeypatch.setattr(db, "upsert", broken_upsert)
        r = R.prebuild(["H5N4"], db, tmp_path / "l.jsonl", http_post=make_stub(), sleep=lambda s: None,
                       source_label="s.tsv")
        assert r.status == "persist_failed"
        assert db.meta()["last_prebuild_status"] == "persist_failed"
    finally:
        db.close()


def test_cpa_readonly_handle_never_writes(tmp_path):
    dbp = _seed_db(tmp_path)
    before_sha, before_meta = _file_sha(dbp), _raw_meta(dbp)
    ro = R.ReferenceDB(dbp, readonly=True)
    try:
        r1 = R.prebuild(["H3N2"], ro, tmp_path / "l.jsonl", source_label="ro.tsv")
        r2 = R.prebuild(["H9N9"], ro, tmp_path / "l.jsonl", allow_network=False, source_label="ro.tsv")
        assert (r1.status, r2.status) == ("nothing_to_do", "disabled")
        assert not any("last_prebuild" in e for e in r1.errors + r2.errors)   # skipped, not attempted-and-failed
    finally:
        ro.close()
    assert _file_sha(dbp) == before_sha and _raw_meta(dbp) == before_meta
    with pytest.raises(sqlite3.OperationalError):
        R.ReferenceDB(dbp, readonly=True).set_meta_many({"x": "y"})


def test_cpa_meta_write_failure_is_soft(tmp_path, monkeypatch):
    dbp = _seed_db(tmp_path)
    db = R.ReferenceDB(dbp)
    try:
        def boom(items):
            raise sqlite3.OperationalError("database is locked (simulated)")
        monkeypatch.setattr(db, "set_meta_many", boom)
        r = R.prebuild(["H3N2"], db, tmp_path / "l.jsonl", source_label="x.tsv")
        assert r.status == "nothing_to_do"
        assert any(e.startswith("last_prebuild meta write failed: OperationalError") for e in r.errors)
        assert r.db.get("n_entries") == 2                       # identity still taken
    finally:
        db.close()


def test_cpa_set_meta_many_is_all_or_nothing(tmp_path):
    dbp = _seed_db(tmp_path)
    db = R.ReferenceDB(dbp)
    try:
        class Bad:
            def __str__(self):
                raise RuntimeError("unrenderable")
        with pytest.raises(Exception):
            db.set_meta_many({"k1": "v1", "k2": Bad()})
        assert "k1" not in db.meta()
    finally:
        db.close()


def test_cpa_summary_exposes_last_prebuild_and_updated(tmp_path):
    dbp = tmp_path / "ref.db"
    db = R.ReferenceDB(dbp)
    try:
        s0 = db.summary()
        assert not any(k.startswith("last_prebuild") for k in s0)   # absent until a prebuild ran
        assert {"n_entries", "by_status", "by_source", "db_name", "db_version", "schema_version"} <= set(s0)
        R.prebuild(["H5N4"], db, tmp_path / "l.jsonl", http_post=make_stub(), sleep=lambda s: None,
                   source_label="a.tsv")
        s1 = db.summary()
        assert s1["last_prebuild_source"] == "a.tsv" and s1["last_prebuild_status"] == "ok"
        assert s1["last_prebuild_utc"] and s1["updated_utc"]
        assert "attribution" not in s1 and "endpoint" not in s1      # filter still a filter
    finally:
        db.close()


def test_cpa_cli_passes_source_basename(tmp_path, capsys):
    import msp_CLI_glytoucan_prebuild as cli
    src = tmp_path / "my source 数据.tsv"
    pd.DataFrame({"composition": ["F1H5N4S2", "H5N4"]}).to_csv(src, sep="\t", index=False)
    dbp = tmp_path / "ref.db"
    rc = cli.main(["--source", str(src), "--db", str(dbp), "--log", str(tmp_path / "l.jsonl"), "--no-network"])
    capsys.readouterr()
    assert rc == 0
    m = _raw_meta(dbp)
    assert m["last_prebuild_source"] == "my source 数据.tsv" and m["last_prebuild_status"] == "disabled"
    rc = cli.main(["--summary", "--db", str(dbp), "--json"])
    out = json.loads(capsys.readouterr().out)
    assert rc == 0 and out["last_prebuild_source"] == "my source 数据.tsv"


# ---------------------------------------------------------------- G1 — snapshot-consistent identity()
def _consistent(ident, states):
    """ident must match exactly one recorded (updated_utc, n_entries) state."""
    return (ident["updated_utc"], ident["n_entries"]) in states


def test_g1_identity_shape_unchanged(tmp_path):
    dbp = _seed_db(tmp_path)
    ro = R.ReferenceDB(dbp, readonly=True)
    try:
        ident = ro.identity()
    finally:
        ro.close()
    assert set(ident) == {"name", "version", "schema_version", "path", "sha256", "n_entries", "updated_utc"}
    assert ident["n_entries"] == 2 and ident["sha256"] == _file_sha(dbp)


def test_g1_commit_between_meta_and_count_cannot_split_the_snapshot(tmp_path):
    """A second writable connection tries to commit (new row + new updated_utc) right when identity() reads the
    count. Either the commit is held off by the read transaction, or it lands before the snapshot — the
    reported updated_utc / n_entries / sha256 must describe one state either way."""
    dbp = _seed_db(tmp_path)
    state_a = (_raw_meta(dbp)["updated_utc"], 2)
    state_b = ("2099-01-01T00:00:00Z", 3)
    ro = R.ReferenceDB(dbp, readonly=True)
    outcome = {}

    def trace(sql):
        if "COUNT(*) FROM entries" in sql and not outcome:
            w = sqlite3.connect(str(dbp), timeout=0.05)
            try:
                w.execute("INSERT INTO entries(composition_key, query_json, status, source, added_utc, last_checked_utc) "
                          "VALUES ('H9N9', '{}', 'resolved', 'api', 'x', 'x')")
                w.execute("INSERT OR REPLACE INTO meta(key, value) VALUES ('updated_utc', ?)", (state_b[0],))
                w.commit()
                outcome["committed"] = True
            except sqlite3.OperationalError:
                w.rollback()
                outcome["committed"] = False
            finally:
                w.close()
    try:
        ro._connect().set_trace_callback(trace)
        ident = ro.identity()
    finally:
        ro.close()
    assert outcome, "injection did not fire"
    assert _consistent(ident, {state_a, state_b}), (ident, outcome)
    assert ident["sha256"] == _file_sha(dbp)


def test_g1_commits_while_hashing_are_held_off(tmp_path, monkeypatch):
    """Codex DB1 F1: a writer trying to commit while identity() hashes — once or twice — cannot change the file
    before the hash is taken, so updated_utc / n_entries / sha256 describe one committed state."""
    dbp = _seed_db(tmp_path)
    state_a = (_raw_meta(dbp)["updated_utc"], 2)
    ro = R.ReferenceDB(dbp, readonly=True)
    real_sha = ro.sha256
    outcomes = []

    def try_commit(i):
        w = sqlite3.connect(str(dbp), timeout=0.05)
        try:
            w.execute("INSERT INTO entries(composition_key, query_json, status, source, added_utc, last_checked_utc) "
                      f"VALUES ('H9N{i}', '{{}}', 'resolved', 'api', 'x', 'x')")
            w.execute("INSERT OR REPLACE INTO meta(key, value) VALUES ('updated_utc', ?)", (f"2099-0{i}-01T00:00:00Z",))
            w.commit()
            outcomes.append(True)
        except sqlite3.OperationalError:
            w.rollback()
            outcomes.append(False)
        finally:
            w.close()

    def sha_with_commits():
        try_commit(1)
        h = real_sha()
        try_commit(2)
        return h
    monkeypatch.setattr(ro, "sha256", sha_with_commits)
    try:
        ident = ro.identity()
    finally:
        ro.close()
    assert outcomes == [False, False]                       # both held off by the read transaction
    assert (ident["updated_utc"], ident["n_entries"]) == state_a
    assert ident["sha256"] == _file_sha(dbp)                 # file unchanged since


def test_f6_summary_aggregates_come_from_one_state(tmp_path):
    """Codex DB1 F6: a commit attempted between summary()'s aggregate queries cannot split the totals."""
    dbp = _seed_db(tmp_path)
    ro = R.ReferenceDB(dbp, readonly=True)
    outcome = {}

    def trace(sql):
        if "GROUP BY source" in sql and not outcome:
            w = sqlite3.connect(str(dbp), timeout=0.05)
            try:
                w.execute("INSERT INTO entries(composition_key, query_json, status, source, added_utc, last_checked_utc) "
                          "VALUES ('H9N9', '{}', 'resolved', 'api', 'x', 'x')")
                w.commit()
                outcome["committed"] = True
            except sqlite3.OperationalError:
                w.rollback()
                outcome["committed"] = False
            finally:
                w.close()
    try:
        ro._connect().set_trace_callback(trace)
        s = ro.summary()
    finally:
        ro.close()
    assert outcome == {"committed": False}
    assert s["n_entries"] == sum(s["by_status"].values()) == sum(s["by_source"].values()) == 2


def test_g1_identity_leaves_a_writers_open_transaction_alone(tmp_path):
    """On a writable handle mid-transaction, identity() must neither commit nor roll back that transaction."""
    dbp = _seed_db(tmp_path)
    db = R.ReferenceDB(dbp)
    try:
        conn = db._connect()
        conn.execute("BEGIN")
        conn.execute("INSERT INTO entries(composition_key, query_json, status, source, added_utc, last_checked_utc) "
                     "VALUES ('H7N7', '{}', 'resolved', 'api', 'x', 'x')")
        db.identity()
        assert conn.in_transaction                     # still open
        conn.rollback()
        assert db.count() == 2                         # the unfinished write was never made durable
    finally:
        db.close()


# ================================================================ CP-B — v14 session state + settings (§5)
SRC = os.path.join(ROOT, "src")
V14 = os.path.join(SRC, "mspfileloaderv14.py")


def _prologue_source():
    """v14 source up to (not including) the first `import pandas as pd`: GLYTOUCAN_FILL, the settings helpers, the
    import-time settings load and `_glytoucan_fill_df` / `_glytoucan_record`. Only a moved anchor skips."""
    lines = open(V14, "rb").read().decode("utf-8").splitlines(keepends=True)
    end = next((i for i, l in enumerate(lines) if l.startswith("import pandas as pd")), None)
    if end is None:
        pytest.skip("v14 anchor `import pandas as pd` not found (coverage LOST, update this test)")
    return "".join(lines[:end])


@pytest.fixture
def home(tmp_path, monkeypatch):
    """Isolated per-user data dir (HOME on POSIX, LOCALAPPDATA on Windows); env overrides unset."""
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setenv("HOME", str(h))
    monkeypatch.setenv("LOCALAPPDATA", str(h))
    monkeypatch.delenv("GLYCOMSP_REFERENCE_DB", raising=False)
    monkeypatch.delenv("GLYCOMSP_GLYTOUCAN_LOG", raising=False)
    monkeypatch.chdir(SRC)
    return h


def _data_dir():
    return R.default_data_dir()


def _exec_prologue(capsys=None):
    ns = {"__name__": "v14_prologue", "__file__": V14}
    old_hook = sys.excepthook
    try:
        exec(compile(_prologue_source(), "v14_prologue", "exec"), ns)
    finally:
        sys.excepthook = old_hook
    out = capsys.readouterr().out if capsys is not None else ""
    return ns, out


def _write_settings(obj_or_text):
    d = _data_dir()
    d.mkdir(parents=True, exist_ok=True)
    p = d / "settings.json"
    p.write_text(obj_or_text if isinstance(obj_or_text, str) else json.dumps(obj_or_text), encoding="utf-8")
    return p


def test_cpb_first_launch_defaults_and_no_file_created(home, capsys):
    ns, out = _exec_prologue(capsys)
    assert ns["GLYTOUCAN_FILL"] == {"enabled": False, "db_path": None, "log_path": None}
    assert ns["GLYTOUCAN_REMEMBER"] == {"on": False}
    assert "[settings]" not in out
    assert not (_data_dir() / "settings.json").exists()           # §9 invariant 2
    assert not list(home.rglob("settings.json*"))


def test_cpb_round_trip_applies_only_when_remember_true(home, tmp_path, capsys):
    ns, _ = _exec_prologue(capsys)
    dbp = tmp_path / "dir with space" / "my_ref.db"
    assert ns["_glytoucan_save_settings"](remember=True, fill_on_run=True, db_path=dbp) is None
    sp = _data_dir() / "settings.json"
    data = json.loads(sp.read_text(encoding="utf-8"))
    assert data["schema_version"] == "1.0.0" and data["updated_utc"].endswith("Z")
    assert data["glytoucan"] == {"remember": True, "fill_on_run": True, "db_path": str(dbp).replace("\\", "/")}
    assert not list(_data_dir().glob("*.tmp"))                    # atomic replace leaves no temp file
    ns2, out2 = _exec_prologue(capsys)
    GF = ns2["GLYTOUCAN_FILL"]
    assert GF["enabled"] is True and str(GF["db_path"]) == str(dbp)          # applied even though it does not exist
    assert str(GF["log_path"]) == str(dbp.with_name("my_ref.jsonl"))       # §3.3: log follows a non-default DB
    assert ns2["GLYTOUCAN_REMEMBER"]["on"] is True and "[settings]" not in out2
    # remember false => nothing applied, file untouched by the load
    _write_settings({"schema_version": "1.0.0", "glytoucan": {"remember": False, "fill_on_run": True, "db_path": str(dbp)}})
    before = sp.read_bytes()
    ns3, _ = _exec_prologue(capsys)
    assert ns3["GLYTOUCAN_FILL"] == {"enabled": False, "db_path": None, "log_path": None}
    assert ns3["GLYTOUCAN_REMEMBER"]["on"] is False and sp.read_bytes() == before
    # remember must be exactly true (not a truthy string)
    _write_settings({"glytoucan": {"remember": "yes", "fill_on_run": True}})
    ns4, _ = _exec_prologue(capsys)
    assert ns4["GLYTOUCAN_FILL"]["enabled"] is False


@pytest.mark.parametrize("content", ["{not json", "[1, 2]", '{"glytoucan": "on"}', "﻿\x00"])
def test_cpb_corrupt_or_unreadable_file_keeps_defaults_with_warning(home, capsys, content):
    p = _write_settings(content)
    before = p.read_bytes()
    ns, out = _exec_prologue(capsys)
    assert ns["GLYTOUCAN_FILL"] == {"enabled": False, "db_path": None, "log_path": None}
    assert "[settings][WARN]" in out and "using defaults" in out
    assert p.read_bytes() == before                               # load never writes


def test_cpb_unknown_keys_preserved_on_save(home, capsys):
    _write_settings({"schema_version": "1.0.0", "future_consumer": {"x": [1, 2]},
                     "glytoucan": {"remember": False, "fill_on_run": False, "db_path": None, "future_key": "keep"}})
    ns, _ = _exec_prologue(capsys)
    assert ns["_glytoucan_save_settings"](fill_on_run=True, remember=True) is None
    data = json.loads((_data_dir() / "settings.json").read_text(encoding="utf-8"))
    assert data["future_consumer"] == {"x": [1, 2]} and data["glytoucan"]["future_key"] == "keep"
    assert data["glytoucan"]["fill_on_run"] is True and data["glytoucan"]["db_path"] is None
    # "Set as default" alone touches only db_path
    assert ns["_glytoucan_save_settings"](db_path="/tmp/x.db") is None
    g = json.loads((_data_dir() / "settings.json").read_text(encoding="utf-8"))["glytoucan"]
    assert g == {"remember": True, "fill_on_run": True, "db_path": "/tmp/x.db", "future_key": "keep"}


def test_cpb_env_db_takes_precedence_over_settings(home, tmp_path, monkeypatch, capsys):
    _write_settings({"glytoucan": {"remember": True, "fill_on_run": True, "db_path": "/remembered/ref.db"}})
    monkeypatch.setenv("GLYCOMSP_REFERENCE_DB", str(tmp_path / "env.db"))
    ns, out = _exec_prologue(capsys)
    GF = ns["GLYTOUCAN_FILL"]
    assert GF["enabled"] is True and GF["db_path"] is None and GF["log_path"] is None
    assert out.count("GLYCOMSP_REFERENCE_DB is set; remembered database path ignored") == 1
    # Remember-save while env governs keeps the remembered path instead of overwriting it with 'default'
    assert ns["_glytoucan_save_remembered"]() is None
    g = json.loads((_data_dir() / "settings.json").read_text(encoding="utf-8"))["glytoucan"]
    assert g["db_path"] == "/remembered/ref.db" and g["remember"] is True


def test_cpb_save_failure_is_reported_not_raised(home, tmp_path, capsys):
    ns, _ = _exec_prologue(capsys)
    blocker = tmp_path / "afile"
    blocker.write_text("x")
    err = ns["_glytoucan_save_settings"](remember=True, path=blocker / "settings.json")
    assert isinstance(err, str) and err


def test_f2_failed_replace_removes_its_own_temp_file(home, tmp_path, capsys):
    """Codex DB1 F2: unique temp per save; a failed os.replace leaves no temp file and the old file intact."""
    ns, _ = _exec_prologue(capsys)
    d = tmp_path / "s"
    (d / "settings.json").mkdir(parents=True)                      # replace onto a directory fails
    err = ns["_glytoucan_save_settings"](remember=True, path=d / "settings.json")
    assert err and not list(d.glob("*.tmp")) and (d / "settings.json").is_dir()


def test_f2_concurrent_saves_use_distinct_temp_files(home, capsys, monkeypatch):
    ns, _ = _exec_prologue(capsys)
    import tempfile as _tf
    seen = []
    real = _tf.mkstemp

    def spy(*a, **k):
        fd, p = real(*a, **k)
        seen.append(p)
        return fd, p
    monkeypatch.setattr(_tf, "mkstemp", spy)
    assert ns["_glytoucan_save_settings"](remember=True) is None
    assert ns["_glytoucan_save_settings"](fill_on_run=True) is None
    assert len(seen) == 2 and seen[0] != seen[1] and all(p.endswith(".json.tmp") for p in seen)
    assert json.loads((_data_dir() / "settings.json").read_text(encoding="utf-8"))["glytoucan"]["fill_on_run"] is True


def test_f4_failed_application_leaves_defaults(home, capsys):
    """Codex DB1 F4: a remembered db_path whose log cannot be derived ("/" has no file name) applies nothing."""
    _write_settings({"glytoucan": {"remember": True, "fill_on_run": True, "db_path": "/"}})
    ns, out = _exec_prologue(capsys)
    assert ns["GLYTOUCAN_FILL"] == {"enabled": False, "db_path": None, "log_path": None}
    assert ns["GLYTOUCAN_REMEMBER"]["on"] is False
    assert "[settings][WARN] could not apply" in out and "using defaults" in out


def test_cpb_log_for_db_default_and_other(home, tmp_path, capsys):
    ns, _ = _exec_prologue(capsys)
    f = ns["_glytoucan_log_for_db"]
    assert f(None) is None and f(R.default_db_path()) is None
    assert str(f(tmp_path / "a.b" / "ref_v2.db")) == str(tmp_path / "a.b" / "ref_v2.jsonl")


def test_cpb_fill_df_honours_session_db_path(home, tmp_path, capsys):
    ns, _ = _exec_prologue(capsys)
    fill, GF = ns["_glytoucan_fill_df"], ns["GLYTOUCAN_FILL"]
    dbp = _seed_db(tmp_path, name="session.db", keys=("F1H5N4S2",))
    df = pd.DataFrame({"composition": ["F1H5N4S2"], "GlyToucan ID": [""], "WURCS": [""]})
    # flag off => None, frame untouched, even with a session DB set
    GF["db_path"] = dbp
    assert fill(df, "composition") is None and df["WURCS"].tolist() == [""]
    GF["enabled"] = True
    blk = fill(df, "composition")
    assert blk["status"] == "filled" and blk["db"]["path"] == str(dbp) and df["WURCS"].iloc[0].startswith("WURCS=")
    # missing session DB => db_missing, nothing created (default path is NOT silently used instead)
    missing = tmp_path / "nope" / "x.db"
    GF["db_path"] = missing
    blk = fill(pd.DataFrame({"composition": ["F1H5N4S2"]}), "composition")
    assert blk["status"] == "db_missing" and not missing.parent.exists()
    assert not R.default_db_path().exists()
    # None => resolver default (env -> data dir), still read-only
    GF["db_path"] = None
    blk = fill(pd.DataFrame({"composition": ["F1H5N4S2"]}), "composition")
    assert blk["status"] == "db_missing" and not R.default_db_path().exists()


# ================================================================ CP-C — subprocess runner (§6)
import subprocess  # noqa: E402
import time  # noqa: E402

STUB = os.path.join(ROOT, "tests", "fixtures", "db1_stub_cli.py")
CLI = os.path.join(SRC, "msp_CLI_glytoucan_prebuild.py")


@pytest.fixture
def v14p(home, capsys):
    ns, _ = _exec_prologue(capsys)
    return ns


def _drive(runner, timeout=15.0):
    """Poll like root.after(100) would, without Tk. Returns when poll() reports the process finished."""
    t0 = time.time()
    while runner.poll():
        assert time.time() - t0 < timeout, "runner did not finish"
        time.sleep(0.02)


def _run_stub(ns, *stub_args, op="prebuild"):
    r = ns["_GtCliRunner"]()
    got, done = [], []
    r.start([sys.executable, STUB, *stub_args], op, on_line=got.append,
            on_exit=lambda text, level, rc: done.append((text, level, rc)))
    _drive(r)
    assert len(done) == 1 and not r.running
    return got, done[0]


def test_cpc_argv_prebuild_exact(v14p, tmp_path):
    r = v14p["_GtCliRunner"](python="/py/python3", script="/app/src/msp_CLI_glytoucan_prebuild.py")
    src, db, log = tmp_path / "a b%c 数据.tsv", tmp_path / "r.db", tmp_path / "r.jsonl"
    assert r.argv("prebuild", db=db, log=log, source=src) == [
        "/py/python3", "/app/src/msp_CLI_glytoucan_prebuild.py",
        f"--source={src}", f"--db={db}", f"--log={log}", "--batch-size=30", "--pause=10.0"]
    assert r.argv("prebuild", db=db, log=log, source=src, column="Structure", sheet="MSlist",
                  no_network=True, batch_size=12, pause=2) == [
        "/py/python3", "/app/src/msp_CLI_glytoucan_prebuild.py",
        f"--source={src}", f"--db={db}", f"--log={log}", "--column=Structure", "--sheet=MSlist",
        "--no-network", "--batch-size=12", "--pause=2.0"]
    assert "--json" not in r.argv("prebuild", db=db, log=log, source=src)
    with pytest.raises(ValueError):
        r.argv("prebuild", db=db, log=log)


def test_cpc_argv_rebuild_exact_and_defaults_always_explicit(v14p):
    r = v14p["_GtCliRunner"](python="py", script="cli.py")
    assert r.argv("rebuild_log", db="/d/r.db", log="/d/r.jsonl") == [
        "py", "cli.py", "--rebuild-log", "--db=/d/r.db", "--log=/d/r.jsonl", "--json"]
    a = r.argv("rebuild_log", db=None, log=None)          # session None => resolved resolver defaults, still explicit
    assert f"--db={R.default_db_path()}" in a and f"--log={R.default_log_path()}" in a
    with pytest.raises(ValueError):
        r.argv("summary", db=None, log=None)


def test_cpc_argv_is_accepted_by_the_real_cli_parser(v14p, tmp_path):
    """The built argv must parse with the real CLI (no drift between builder and argparse)."""
    import msp_CLI_glytoucan_prebuild as cli
    src = tmp_path / "s.tsv"
    pd.DataFrame({"composition": ["H5N4"]}).to_csv(src, sep="\t", index=False)
    r = v14p["_GtCliRunner"]()
    a = r.argv("prebuild", db=tmp_path / "r.db", log=tmp_path / "r.jsonl", source=src, no_network=True,
               batch_size=5, pause=0)
    assert cli.main(a[2:]) == 0 and _raw_meta(tmp_path / "r.db")["last_prebuild_source"] == "s.tsv"


def test_f3_leading_dash_values_reach_the_real_cli(v14p, tmp_path, monkeypatch):
    """Codex DB1 F3: relative file names and a column starting with '-' are passed as values, not options."""
    import msp_CLI_glytoucan_prebuild as cli
    monkeypatch.chdir(tmp_path)
    pd.DataFrame({"-comp": ["H5N4", "F1H5N4S2"]}).to_csv("-weird.tsv", sep="\t", index=False)
    r = v14p["_GtCliRunner"]()
    a = r.argv("prebuild", db="-ref.db", log="-ref.jsonl", source="-weird.tsv", column="-comp", no_network=True)
    assert cli.main(a[2:]) == 0
    m = _raw_meta(tmp_path / "-ref.db")
    assert m["last_prebuild_source"] == "-weird.tsv" and m["last_prebuild_status"] == "disabled"
    a = r.argv("rebuild_log", db="-ref.db", log="-ref.jsonl")
    assert cli.main(a[2:]) == 0 and (tmp_path / "-ref.jsonl").exists()


def test_cpc_availability_reasons(v14p, tmp_path):
    R_ = v14p["_GtCliRunner"]
    assert R_().available, R_().reason                                   # real interpreter + CLI beside v14
    assert "interpreter" in R_(python=str(tmp_path / "nope")).reason
    assert "CLI not found" in R_(script=str(tmp_path / "nope.py")).reason
    assert R_().script == CLI


def test_cpc_lines_arrive_in_order_and_exit_0_mapped(v14p):
    got, (text, level, rc) = _run_stub(v14p, "--unicode")
    assert got[:3] == ["[prebuild] source stub.tsv: 3 distinct compositions",
                       "[prebuild] 3 compositions in source; 0 already in database; 3 to fetch",
                       "[prebuild] batch 1/1 (3 compositions) — H5N4 数据"]
    assert got[3].startswith("[prebuild] status=ok ") and len(got) == 5
    assert (text, level, rc) == ("done: status=ok committed=3 api_calls=1", "ok", 0)


def test_cpc_exit_2_is_storage_problem(v14p):
    _, (text, level, rc) = _run_stub(v14p, "--status", "log_failed", "--rc", "2")
    assert (text, level, rc) == ("storage problem: log_failed; see log pane", "error", 2)


def test_cpc_fail_soft_and_missing_status_line(v14p):
    _, (text, level, _) = _run_stub(v14p, "--status", "offline")
    assert text == "network fail-soft — rerun to resume (n_unsent=2)" and level == "warn"
    _, (text, level, rc) = _run_stub(v14p, "--no-status", "--rc", "1")
    assert (text, level, rc) == ("process ended without a status line, exit 1", "error", 1)


def test_cpc_rebuild_json_mapping(v14p):
    ok = json.dumps({"status": "ok", "rebuilt_lines": 7, "log": "/x/l.jsonl"}, indent=2)
    _, (text, level, _) = _run_stub(v14p, "--json-result", ok, op="rebuild_log")
    assert text == "rebuilt 7 lines → /x/l.jsonl (previous kept as .bak)" and level == "ok"
    bad = json.dumps({"status": "rebuild_failed", "error": "rebuild-log refused: not a GlycoMSP reference database"})
    _, (text, level, rc) = _run_stub(v14p, "--json-result", bad, "--rc", "2", op="rebuild_log")
    assert (text, level, rc) == ("rebuild-log refused: not a GlycoMSP reference database", "error", 2)
    _, (text, level, _) = _run_stub(v14p, "--no-status", op="rebuild_log")
    assert text.startswith("process ended without a JSON result") and level == "error"


def _start_hanging(ns, *extra, kill_after=3.0):
    r = ns["_GtCliRunner"](kill_after=kill_after)
    got, done = [], []
    r.start([sys.executable, STUB, "--hang", *extra], "prebuild", on_line=got.append,
            on_exit=lambda text, level, rc: done.append((text, level, rc)))
    t0 = time.time()
    while "[prebuild] ready" not in got:                     # stub is past signal setup
        r.poll()
        assert time.time() - t0 < 15
        time.sleep(0.02)
    return r, got, done


def test_cpc_cancel_mid_stream_yields_cancelled_result(v14p):
    r, got, done = _start_hanging(v14p)
    assert r.running and r.cancel() is True
    _drive(r)
    assert done and done[0][0] == "cancelled — batches already committed are kept; rerun to resume"
    assert done[0][1] == "warn" and done[0][2] != 0 and r.proc is None
    assert r.cancel() is False                                # nothing left to cancel


@pytest.mark.skipif(sys.platform == "win32", reason="SIGTERM cannot be ignored on Windows (terminate == kill)")
def test_cpc_terminate_escalates_to_kill(v14p):
    r, _, done = _start_hanging(v14p, "--ignore-term", kill_after=0.3)
    t0 = time.time()
    r.cancel()
    _drive(r)
    assert done and done[0][1] == "warn" and done[0][2] == -9 and time.time() - t0 < 10


def test_cpc_shutdown_stops_a_live_child(v14p):
    r, _, _ = _start_hanging(v14p, kill_after=2.0)
    proc = r.proc
    r.shutdown()
    assert proc.poll() is not None


def test_cpc_second_start_refused_while_running(v14p):
    r, _, _ = _start_hanging(v14p)
    try:
        with pytest.raises(RuntimeError):
            r.start([sys.executable, STUB], "prebuild")
    finally:
        r.cancel()
        _drive(r)


def test_cpc_status_regex_matches_the_real_cli_line(v14p, tmp_path):
    """Run the real CLI (no network) as the runner does and map its genuine output."""
    src = tmp_path / "s.tsv"
    pd.DataFrame({"composition": ["H5N4", "F1H5N4S2"]}).to_csv(src, sep="\t", index=False)
    r = v14p["_GtCliRunner"]()
    got, done = [], []
    r.start(r.argv("prebuild", db=tmp_path / "r.db", log=tmp_path / "r.jsonl", source=src, no_network=True),
            "prebuild", on_line=got.append, on_exit=lambda *a: done.append(a))
    _drive(r, timeout=60)
    assert any(l.startswith("[prebuild] status=disabled ") for l in got), got
    assert done[0] == ("done: status=disabled committed=0 api_calls=0", "ok", 0)


def test_cpc_result_mapper_pure_cases(v14p):
    f = v14p["_gt_cli_result"]
    line = "[prebuild] status=partial committed=5 already=1 unsent=25 unpersisted=0 unlogged=0 api_calls=2"
    assert f("prebuild", 0, ["x", line, "[prebuild] database: p (6 entries, v)"]) == \
        ("network fail-soft — rerun to resume (n_unsent=25)", "warn")
    assert f("prebuild", 0, [line], cancelled=True)[1] == "warn"            # finished before the cancel landed
    assert f("prebuild", -15, [], cancelled=True)[0].startswith("cancelled")
    assert f("prebuild", 1, [line.replace("partial", "ok")])[1] == "error"   # unexpected non-zero exit
    assert f("prebuild", 0, ["[prebuild] status=ok but free text"])[0].startswith("process ended without")


# ================================================================ CP-D — GUI helpers, popup text, static audit (§3, §4)
import ast  # noqa: E402


def _block(status="filled", n_user_wurcs=0, n_missing=0, fam_db="glycomsp_reference"):
    return {"status": status, "route": "database",
            "db": {"name": fam_db, "version": "v0.01.20260921", "path": "/x/ref.db", "sha256": "ab" * 32,
                   "n_entries": 2230, "updated_utc": "2026-09-23T09:21:44Z"},
            "n_compositions": 10, "n_filled": 10 - n_missing, "n_missing": n_missing, "n_user_accession": 1,
            "n_user_wurcs": n_user_wurcs, "last_resolved_utc": "2026-09-23T09:21:44Z"}


def test_cpd_reference_info_text_cga_mas_both_none(v14p):
    f = v14p["_reference_info_text"]
    assert f({}) == "No reference fill recorded for this sample."
    assert f(None) == "No reference fill recorded for this sample."
    assert f({"csv": "/a.csv", "_glytoucan_resolution_XYZ": _block()}) == "No reference fill recorded for this sample."
    cga = f({"_glytoucan_resolution_CGA": _block(n_missing=2)})
    assert cga.startswith("CGA: reference fill filled") and "MAS:" not in cga
    assert "compositions: 10   filled: 8   not in database: 2" in cga
    assert "path: /x/ref.db" in cga and "entries: 2230   updated: 2026-09-23T09:21:44Z" in cga
    assert "last resolved: 2026-09-23T09:21:44Z" in cga
    mas = f({"_glytoucan_resolution_MAS": _block(status="partial")})
    assert mas.startswith("MAS: reference fill partial") and "CGA:" not in mas
    both = f({"_glytoucan_resolution_CGA": _block(), "_glytoucan_resolution_MAS": _block()})
    assert both.index("CGA:") < both.index("MAS:") and "\n\n" in both


def test_cpd_reference_info_text_p1_local_annotation_wording(v14p):
    """P-1: a user WURCS with a blank accession is a local annotation — never 'missing', 'unresolved' or wurcs_only."""
    txt = v14p["_reference_info_text"]({"_glytoucan_resolution_CGA": _block(n_user_wurcs=3)})
    assert "user WURCS kept (local annotation): 3" in txt
    low = txt.lower()
    assert "wurcs_only" not in low and "unresolved" not in low and "missing" not in low
    assert "regist" not in low                                  # no registration nudge


def test_cpd_reference_info_text_tolerates_malformed_block(v14p):
    txt = v14p["_reference_info_text"]({"_glytoucan_resolution_CGA": {"status": "db_missing", "db": "oops"}})
    assert txt.startswith("CGA: reference fill db_missing") and "last resolved: —" in txt


def test_cpd_db_label_text_states(v14p, tmp_path):
    ns = v14p
    lab, GF = ns["_glytoucan_db_label_text"], ns["GLYTOUCAN_FILL"]
    GF["db_path"] = tmp_path / "none.db"
    assert lab() == f"Active DB: not found ({tmp_path / 'none.db'})" and not (tmp_path / "none.db").exists()
    GF["db_path"] = _seed_db(tmp_path)
    assert lab() == f"Active DB: glycomsp_reference {R.DEFAULT_DB_VERSION} · 2 entries"
    foreign = tmp_path / "foreign.db"
    foreign.write_bytes(b"not sqlite at all")
    GF["db_path"] = foreign
    assert lab() == f"Active DB: unreadable ({foreign})" and foreign.read_bytes() == b"not sqlite at all"

    class Busy:
        running = True
    ns["_GT_UI"]["runner"] = Busy()
    try:
        assert lab() == "Active DB: (updating…)"
    finally:
        ns["_GT_UI"]["runner"] = None


def test_cpd_manager_summary_is_read_only_and_complete(v14p, tmp_path):
    ns = v14p
    summ, GF = ns["_glytoucan_manager_summary"], ns["GLYTOUCAN_FILL"]
    # missing DB: nothing created
    GF["db_path"], GF["log_path"] = tmp_path / "m" / "none.db", tmp_path / "m" / "none.jsonl"
    s = summ()
    assert s["db_state"] == "missing" and s["name"].startswith("not found") and s["log"].endswith("not found")
    assert not (tmp_path / "m").exists()
    # good DB with log + .bak
    dbp = _seed_db(tmp_path, keys=("F1H5N4S2", "H3N2", "H5N4"))
    log = tmp_path / "seed.jsonl"
    R.rebuild_log(R.ReferenceDB(dbp, readonly=True), log)      # makes the log (3 lines); a 2nd run makes the .bak
    R.rebuild_log(R.ReferenceDB(dbp, readonly=True), log)
    GF["db_path"], GF["log_path"] = dbp, log
    before = _file_sha(dbp), log.read_bytes()
    s = summ()
    assert s["db_state"] == "ok" and s["path"] == str(dbp)
    assert s["name"] == f"glycomsp_reference    Version {R.DEFAULT_DB_VERSION}" and s["schema"] == ""
    assert s["entries"].startswith("3   resolved ") and "by source api 3" in s["entries"]
    assert s["sha_full"] == _file_sha(dbp) and s["hash"] == s["sha_full"][:12] + "…"
    assert s["last"].endswith("· ok") and "3 lines" in s["log"] and ".bak present" in s["log"]
    assert (_file_sha(dbp), log.read_bytes()) == before       # §9 invariant 3
    # schema mismatch shown; foreign file shows the resolver's own error text
    w = R.ReferenceDB(dbp)
    w.set_meta("schema_version", "9.9.9")
    w.close()
    assert "schema 9.9.9" in summ()["schema"]
    foreign = tmp_path / "f.db"
    c = sqlite3.connect(str(foreign))
    c.execute("CREATE TABLE other(x)")
    c.commit()
    c.close()
    GF["db_path"] = foreign
    s = summ()
    assert s["db_state"] == "error" and "not a GlycoMSP reference database" in s["name"]


def test_cpd_set_enabled_saves_only_with_remember_and_syncs_vars(v14p, tmp_path):
    ns = v14p

    class Var:
        def __init__(self):
            self.v = None

        def set(self, v):
            self.v = v
    alive, dead = object(), object()
    ns["widget_alive"] = lambda w: w is alive                  # defined later in v14; the light prologue lacks it
    lv, fv, dfv = Var(), Var(), Var()
    ns["_GT_UI"]["prepare"][:] = [{"widget": alive, "label_var": lv, "fill_var": fv},
                                  {"widget": dead, "label_var": None, "fill_var": dfv}]
    sp = _data_dir() / "settings.json"
    assert ns["_glytoucan_set_enabled"](True) is None
    assert ns["GLYTOUCAN_FILL"]["enabled"] is True and fv.v is True and lv.v.startswith("Active DB:")
    assert dfv.v is None and len(ns["_GT_UI"]["prepare"]) == 1          # dead window pruned
    assert not sp.exists()                                              # Remember off => no file (§9 inv. 2)
    ns["GLYTOUCAN_REMEMBER"]["on"] = True
    assert ns["_glytoucan_set_enabled"](False) is None
    assert json.loads(sp.read_text(encoding="utf-8"))["glytoucan"]["fill_on_run"] is False and fv.v is False


# ---------------------------------------------------------------- static audit (§10 item 7; UI static audit rule)
def _v14_tree():
    return ast.parse(open(V14, "rb").read().decode("utf-8"))


def _func(tree, name):
    for n in ast.walk(tree):
        if isinstance(n, ast.FunctionDef) and n.name == name:
            return n
    raise AssertionError(f"{name} not found in v14")


def _command_targets(node):
    """Names passed as command= / bound via .bind(..., <callable>) / .after(ms, <callable>) / protocol(..., <callable>)."""
    out = []
    for n in ast.walk(node):
        if isinstance(n, ast.Call):
            for kw in n.keywords:
                if kw.arg == "command" and isinstance(kw.value, ast.Name):
                    out.append(kw.value.id)
            f = n.func
            if isinstance(f, ast.Attribute) and f.attr in ("bind", "after", "protocol") and len(n.args) >= 2 \
                    and isinstance(n.args[1], ast.Name):
                out.append(n.args[1].id)
    return out


def test_cpd_static_every_new_ui_target_is_defined():
    tree = _v14_tree()
    module_defs = {n.name for n in tree.body if isinstance(n, ast.FunctionDef)}
    mgr = _func(tree, "open_manage_references_window")
    nested = {n.name for n in ast.walk(mgr) if isinstance(n, ast.FunctionDef)}
    targets = _command_targets(mgr)
    assert {"_select", "_set_default", "_copy_hash", "_on_fill", "_on_remember", "_prebuild_dialog", "_rebuild",
            "_refresh", "_cancel", "_close", "_pump", "_browse", "_start"} <= set(targets)
    missing = [t for t in targets if t not in nested and t not in module_defs]
    assert not missing, f"undefined UI callbacks in the manager: {missing}"
    # Prepare Dataset: the launcher + the popup target
    prep = _func(tree, "open_prepare_dataset_window")
    prep_nested = {n.name for n in ast.walk(prep) if isinstance(n, ast.FunctionDef)}
    assert "open_manage_references_window" in _command_targets(prep) and "open_manage_references_window" in module_defs
    assert "_show_reference_info" in prep_nested
    # every module-level DB1 helper the GUI calls exists
    for name in ("_glytoucan_set_enabled", "_glytoucan_refresh_labels", "_glytoucan_register_prepare",
                 "_glytoucan_db_label_text", "_glytoucan_manager_summary", "_reference_info_text",
                 "_glytoucan_shutdown_runner", "_glytoucan_save_settings", "_glytoucan_save_remembered",
                 "_glytoucan_log_for_db", "_glytoucan_settings_path", "_glytoucan_active_paths", "widget_alive"):
        assert name in module_defs, name
    # no self._* references in DB1 code (module functions, no class) — the static-audit rule is vacuous but checked
    src = ast.get_source_segment(open(V14, "rb").read().decode("utf-8"), mgr) or ""
    assert "self._" not in src
    # on_closing stops a live runner
    assert "_glytoucan_shutdown_runner" in {n.func.id for n in ast.walk(_func(tree, "on_closing"))
                                            if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}


def test_cpd_prepare_rows_0_to_3_unchanged_and_menu_order():
    """§9 invariants 5-6: rows 0-3 keep their widgets; the two sample items follow "Link and Validate Sample"."""
    src = open(V14, "rb").read().decode("utf-8").replace("\r\n", "\n")
    i = src.index('    button_frame = tk.Frame(subwin)\n')
    grid = src[i:src.index('tk.Button(subwin, text="Close", command=subwin.destroy)', i)]
    for text, row, col in (("Add MS2 CSV(s)", 0, 0), ("Add Sample", 1, 0), ("Negative Sampling Setting", 2, 0),
                           ("Load Method", 2, 3), ("CGA manager", 3, 0), ("CGA → Trainable", 3, 1)):
        j = grid.index(f'text="{text}"')
        assert f".grid(row={row}, column={col}" in grid[j:j + 400], text
    assert 'text="Fill GlyToucan ID / WURCS from reference DB (see Manage References)"' in grid
    assert 'text="Manage References\\u2026"' in grid and ".grid(row=4, column=0" in grid
    k = src.index('label="Link and Validate Sample"')
    a = src.index('label="Reference database info\\u2026"', k)
    b = src.index('label="Resolve missing references (T1-h)", state="disabled"', a)
    r = src.index('menu.add_command(label=f"Remove {filetype.upper()}"', k)
    assert k < a < b < r


# ---------------------------------------------------------------- optional real-Tk smoke (skips without a display)
@pytest.fixture
def v14full(home, capsys, monkeypatch):
    """Whole v14 module up to `root = tk.Tk()` (as test_s1_v14_headless), then a withdrawn real root."""
    tk = pytest.importorskip("tkinter")
    lines = open(V14, "rb").read().decode("utf-8").splitlines(keepends=True)
    root_i = next((i for i, l in enumerate(lines) if l.startswith("root = tk.Tk()")), None)
    if root_i is None:
        pytest.skip("v14 anchor `root = tk.Tk()` not found (coverage LOST, update this test)")
    try:
        root = tk.Tk()
    except tk.TclError as e:
        pytest.skip(f"no display for Tk: {e}")
    root.withdraw()
    ns = {"__name__": "v14_full", "__file__": V14}
    old_hook = sys.excepthook
    try:
        exec(compile("".join(lines[:root_i]), "v14_full", "exec"), ns)
    finally:
        sys.excepthook = old_hook
    ns["root"] = root
    # the manager and its colour table are defined after `root = tk.Tk()` (beside on_closing): exec them from the AST
    src = "".join(lines)
    tree = ast.parse(src)
    for node in tree.body:
        if (isinstance(node, ast.FunctionDef) and node.name == "open_manage_references_window") or \
                (isinstance(node, ast.Assign) and any(getattr(t, "id", "") == "_GT_LEVEL_FG" for t in node.targets)):
            exec(compile(ast.Module(body=[node], type_ignores=[]), "v14_full_tail", "exec"), ns)
    assert "open_manage_references_window" in ns and "_GT_LEVEL_FG" in ns
    capsys.readouterr()
    yield ns
    try:
        r = ns["_GT_UI"].get("runner")
        if r is not None:
            r.shutdown()
        root.destroy()
    except Exception:
        pass


def _pump_until(root, cond, timeout=20.0):
    t0 = time.time()
    while not cond():
        root.update()
        assert time.time() - t0 < timeout, "timed out"
        time.sleep(0.02)


def _find(widget, cls_name, text):
    for w in widget.winfo_children():
        try:
            if w.winfo_class() == cls_name and w.cget("text") == text:
                return w
        except Exception:
            pass
        hit = _find(w, cls_name, text)
        if hit is not None:
            return hit
    return None


def test_cpd_tk_manager_smoke_select_refuse_and_stub_run(v14full, tmp_path, monkeypatch):
    ns, root = v14full, v14full["root"]
    ns["GLYTOUCAN_FILL"]["db_path"] = _seed_db(tmp_path)
    win = ns["open_manage_references_window"]()
    root.update()
    assert ns["open_manage_references_window"]() is win                       # singleton
    assert _find(win, "Button", "Check again (T1-e)").cget("state") == "disabled"
    assert _find(win, "Button", "Cancel").cget("state") == "disabled"
    # Select a non-DB file: refused, previous kept, nothing created
    junk = tmp_path / "junk.db"
    junk.write_text("x")
    errors = []
    monkeypatch.setattr(ns["filedialog"], "askopenfilename", lambda **k: str(junk))
    monkeypatch.setattr(ns["messagebox"], "showerror", lambda *a, **k: errors.append(a))
    prev = ns["GLYTOUCAN_FILL"]["db_path"]
    _find(win, "Button", "Select database…").invoke()
    assert errors and ns["GLYTOUCAN_FILL"]["db_path"] == prev and junk.read_text() == "x"
    # Run the stub through the manager's own pump; buttons follow the busy rule
    runner = ns["_GT_UI"]["runner"]
    real_start = runner.start
    monkeypatch.setattr(runner, "start", lambda argv, op, **kw: real_start([sys.executable, STUB, "--hang"], op, **kw))
    monkeypatch.setattr(runner, "argv", lambda *a, **k: ["py", "cli", "--rebuild-log"])
    monkeypatch.setattr(ns["messagebox"], "askyesno", lambda *a, **k: True)
    _find(win, "Button", "Rebuild log").invoke()
    _pump_until(root, lambda: runner.running)
    assert _find(win, "Button", "Close").cget("state") == "disabled"
    assert _find(win, "Button", "Cancel").cget("state") == "normal"
    assert ns["_glytoucan_db_label_text"]() == "Active DB: (updating…)"
    _find(win, "Button", "Cancel").invoke()
    _pump_until(root, lambda: not runner.running)
    root.update()
    assert _find(win, "Button", "Close").cget("state") == "normal"
    _find(win, "Button", "Close").invoke()
    root.update()
    assert ns["_GT_UI"]["manager"] is None


def test_cpd_tk_prepare_row4_label_and_fill_sync_create_nothing(v14full, home):
    """Real Prepare Dataset window: row-4 launcher + active-DB label; the manager and Prepare Dataset checkboxes
    follow one flag; opening both windows creates no file or directory (§9 invariants 2-3)."""
    import tkinter as tk
    ns, root = v14full, v14full["root"]
    ns["open_prepare_dataset_window"]()
    root.update()
    prep = [w for w in root.winfo_children() if isinstance(w, tk.Toplevel) and w.title() == "Prepare Dataset"]
    assert prep, "Prepare Dataset window not created"
    btn = _find(prep[0], "Button", "Manage References…")
    assert btn is not None and (btn.grid_info()["row"], btn.grid_info()["column"]) == (4, 0)
    cb = _find(prep[0], "Checkbutton", "Fill GlyToucan ID / WURCS from reference DB (see Manage References)")
    assert cb is not None and (cb.grid_info()["row"], cb.grid_info()["column"]) == (3, 2)
    entries = ns["_GT_UI"]["prepare"]
    assert len(entries) == 1 and entries[0]["label_var"].get().startswith("Active DB: not found (")
    btn.invoke()                                                   # opens the manager
    root.update()
    assert ns["_GT_UI"]["manager"] is not None and len(ns["_GT_UI"]["prepare"]) == 2
    ns["_glytoucan_set_enabled"](True)
    assert [e["fill_var"].get() for e in ns["_GT_UI"]["prepare"]] == [True, True]
    ns["_glytoucan_set_enabled"](False)
    assert [e["fill_var"].get() for e in ns["_GT_UI"]["prepare"]] == [False, False]
    assert not any(home.rglob("*"))                                # nothing created under the data dir
