"""V14-FOLLOWUP CP-B: v14 uses the reader's result (ticket WP-A; D-2, D-4).

mspfileloaderv14.py builds its Tk root at module level, so (as in tests/test_s1_v14_headless.py) the module source
is executed up to `root = tk.Tk()`. That defines the module-level helpers and the MetadataEditorWindow class. The
class is then driven headless: Tk widgets, dialogs, the thread and the reader are replaced by fakes, and the real
load_file -> background_conversion -> generate path runs against files in tmp_path. The diagnostics files are
written by the real msp_spectrum_record.ReaderDiagnostics, so their schema is the reader's own.

Covers: dict / False / exception returns; complete, complete with warnings, partial, failed; stale, wrong-type and
missing diagnostics; diagnostics moved under the final save name with their output references rewritten; failed
never promoted; the batch's later files converted into the chosen folder. Codex round 1: malformed diagnostics
never change the reader's outcome (V14F-F1); a file not written by this call is not used (V14F-F2); EXDEV copy
fallback (V14F-F5); the generate() guard, mzML identity and a later batch failure tested directly (V14F-F7).
Codex round 2: a file that could not be read before the call is "unverifiable" (V14F-R2-F1); round 3: so is
one whose lstat failed with anything but FileNotFoundError (V14F-R3-F1).
"""
import glob
import json
import os
import sys
import time
import types

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
sys.path.insert(0, SRC)

import msprawextractor as real_mspext  # noqa: E402
from msp_spectrum_record import ReaderDiagnostics, write_diagnostics  # noqa: E402


@pytest.fixture(scope="module")
def ns():
    path = os.path.join(SRC, "mspfileloaderv14.py")
    lines = open(path, "rb").read().decode("utf-8").splitlines(keepends=True)
    root_i = next((i for i, l in enumerate(lines) if l.startswith("root = tk.Tk()")), None)
    if root_i is None:
        pytest.skip("v14 anchor 'root = tk.Tk()' not found (coverage LOST, update this test)")
    names = {"__name__": "v14_headless_cpb", "__file__": path}
    old_hook, old_cwd = sys.excepthook, os.getcwd()
    try:
        os.chdir(SRC)
        exec(compile("".join(lines[:root_i]), "v14_headless_cpb", "exec"), names)
    finally:
        sys.excepthook = old_hook
        os.chdir(old_cwd)
    return names


# ---------------------------------------------------------------- fixtures: reader files

def _write_diag(raw, outdir, outcome="complete", exit_code=0, warnings=(), explanation="", failure=None,
                csvs=True, reader="msprawextractor", finished_utc=None):
    """Write temp CSVs and the three diagnostics files as the reader would; return the reader-style dict."""
    stem = os.path.splitext(os.path.basename(raw))[0]
    d = ReaderDiagnostics(reader, "0.61", 20261026, raw, stem, "raw")
    res = {}
    if csvs:
        for level in ("ms2", "ms3"):
            p = os.path.join(outdir, f"{level}tmp_{stem}.csv")
            with open(p, "w") as fh:
                fh.write(f"{level}col\n1\n")
            d.wrote(level, p)
            res[f"{level}tmp"] = p
    if failure is not None:
        d.fail(failure)
    done = d.finish(outdir)
    jp = done["diagnostics"]["json"]
    diag = json.load(open(jp, encoding="utf-8"))
    if failure is None:
        diag.update(outcome=outcome, exit_code=exit_code, warnings=list(warnings), explanation=explanation)
    if finished_utc:
        diag["timing"]["finished_utc"] = finished_utc
    write_diagnostics(jp, diag)
    res.update(outcome=diag["outcome"], exit_code=diag["exit_code"], warnings=diag["warnings"],
               diagnostics=done["diagnostics"])
    return res


def _raw(tmp_path, name="s1.raw", folder="in"):
    d = tmp_path / folder
    d.mkdir(exist_ok=True)
    p = d / name
    p.write_bytes(b"x" * 10)
    return str(p)


# ---------------------------------------------------------------- pure helpers

def test_output_folder_rule(ns):
    f = ns["reader_output_folder"]
    assert f("/a/b/s.raw", None) == "/a/b"
    assert f("/a/b/s.raw", "/out") == "/out"
    assert f("s.raw", None) == os.getcwd()


def test_dict_complete_fresh(ns, tmp_path):
    raw = _raw(tmp_path)
    started = time.time()
    res = _write_diag(raw, os.path.dirname(raw))
    s = ns["summarize_reader_result"](res, None, raw, os.path.dirname(raw), started, "msprawextractor")
    assert s["outcome"] == "complete" and s["promote"] and s["dialog"] is None and s["diag_status"] == "fresh"
    assert set(s["paths"]) == {"json", "txt", "status_tsv"}
    assert ns["reader_result_status"](s) == ("Conversion complete. Ready for metadata.", "green")


def test_dict_warnings_and_partial(ns, tmp_path):
    raw = _raw(tmp_path)
    started = time.time()
    res = _write_diag(raw, os.path.dirname(raw), "complete", 3, ["12 of 40 exported MS2 spectra have no charge"])
    s = ns["summarize_reader_result"](res, None, raw, os.path.dirname(raw), started, "msprawextractor")
    lines = ns["reader_result_log_lines"](s)
    assert "[reader] outcome: complete with warnings (exit code 3)" in lines
    assert any("have no charge" in l for l in lines) and s["dialog"] is None
    assert ns["reader_result_status"](s)[1] == "orange"
    res = _write_diag(raw, os.path.dirname(raw), "partial", 4, ["2 spectra could not be decoded"])
    s = ns["summarize_reader_result"](res, None, raw, os.path.dirname(raw), started, "msprawextractor")
    assert s["promote"] and s["dialog"] == "partial"
    assert "partial conversion" in ns["reader_result_message"](s, raw)


def test_dict_failed_is_not_promoted(ns, tmp_path):
    raw = _raw(tmp_path)
    started = time.time()
    res = _write_diag(raw, os.path.dirname(raw), "failed", 5, explanation="No MS2 spectrum with peaks was exported.")
    s = ns["summarize_reader_result"](res, None, raw, os.path.dirname(raw), started, "msprawextractor")
    assert not s["promote"] and s["dialog"] == "failed"
    assert "No MS2 spectrum with peaks" in s["reason"]


def test_false_return_reason_from_fresh_json(ns, tmp_path):
    raw = _raw(tmp_path)
    started = time.time()
    _write_diag(raw, os.path.dirname(raw), failure=OSError("backend could not open the file"), csvs=False)
    s = ns["summarize_reader_result"](False, None, raw, os.path.dirname(raw), started, "msprawextractor")
    assert s["outcome"] == "failed" and not s["promote"]
    assert "could not be opened" in s["reason"] and "backend could not open" in s["reason"]
    assert "no current diagnostics" not in s["reason"]


def test_stale_json_never_overrides_a_failure(ns, tmp_path):
    raw = _raw(tmp_path)
    _write_diag(raw, os.path.dirname(raw), "complete", 0, finished_utc="2020-01-01T00:00:00Z")
    s = ns["summarize_reader_result"](False, None, raw, os.path.dirname(raw), time.time(), "msprawextractor")
    assert s["outcome"] == "failed" and not s["promote"] and s["diag_status"] == "stale"
    assert s["paths"] == {} and "stale" in s["reason"]
    s = ns["summarize_reader_result"](None, RuntimeError("boom"), raw, os.path.dirname(raw), time.time(),
                                      "msprawextractor")
    assert s["outcome"] == "failed" and "RuntimeError: boom" in s["reason"]


def test_stale_json_with_dict_return_is_not_promoted(ns, tmp_path):
    raw = _raw(tmp_path)
    res = _write_diag(raw, os.path.dirname(raw), "complete", 0, finished_utc="2020-01-01T00:00:00Z")
    s = ns["summarize_reader_result"](res, None, raw, os.path.dirname(raw), time.time(), "msprawextractor")
    assert s["outcome"] == "complete" and s["promote"] and s["paths"] == {}


@pytest.mark.parametrize("mutate, status", [
    (lambda d: d.update(json_type="glycomsp.metadata"), "wrong_type"),
    (lambda d: d.update(schema_version="2.0.0"), "wrong_type"),
    (lambda d: d["tool"].update(reader="mzmlreader"), "wrong_type"),
    (lambda d: d["input"].update(size_bytes=11), "stale"),
    (lambda d: d["input"].update(alias="other"), "stale"),
])
def test_diag_checks(ns, tmp_path, mutate, status):
    raw = _raw(tmp_path)
    started = time.time()
    res = _write_diag(raw, os.path.dirname(raw))
    jp = res["diagnostics"]["json"]
    d = json.load(open(jp))
    mutate(d)
    json.dump(d, open(jp, "w"))
    assert ns["load_reader_diagnostics"](raw, os.path.dirname(raw), started, "msprawextractor")[1] == status


def test_missing_and_unreadable(ns, tmp_path):
    raw = _raw(tmp_path)
    load = ns["load_reader_diagnostics"]
    assert load(raw, os.path.dirname(raw), time.time())[1] == "missing"
    open(os.path.join(os.path.dirname(raw), "s1.reader_diag.json"), "w").write("{not json")
    assert load(raw, os.path.dirname(raw), time.time())[1] == "unreadable"


def test_legacy_tuple_and_diagnostics_error(ns, tmp_path):
    raw = _raw(tmp_path)
    s = ns["summarize_reader_result"](("a.csv", "b.csv"), None, raw, os.path.dirname(raw), time.time())
    assert s["outcome"] is None and s["promote"]
    assert ns["reader_result_status"](s) == ("Conversion complete. Ready for metadata.", "green")
    res = {"ms2tmp": "a", "ms3tmp": "b", "outcome": "complete", "exit_code": 0, "warnings": [], "diagnostics": {},
           "diagnostics_error": "PermissionError: denied"}
    s = ns["summarize_reader_result"](res, None, raw, os.path.dirname(raw), time.time())
    assert s["promote"] and any("could not be written" in l for l in ns["reader_result_log_lines"](s))


def test_summarize_never_raises(ns):
    s = ns["summarize_reader_result"](False, None, None, None, time.time())
    assert s["outcome"] == "failed" and not s["promote"]


def test_promote_moves_and_rewrites(ns, tmp_path):
    raw = _raw(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    started = time.time()
    res = _write_diag(raw, os.path.dirname(raw))
    s = ns["summarize_reader_result"](res, None, raw, os.path.dirname(raw), started, "msprawextractor")
    log = types.SimpleNamespace(lines=[])
    log.log = log.lines.append
    moved = ns["promote_reader_diagnostics"](s["paths"], str(out), "proj_X_s1", str(out / "ms2_proj_X_s1.csv"),
                                             str(out / "ms3_proj_X_s1.csv"), log)
    assert sorted(os.listdir(out)) == ["proj_X_s1.reader_diag.json", "proj_X_s1.reader_diag.txt",
                                       "proj_X_s1.scan_status.tsv"]
    assert set(moved) == {"json", "txt", "status_tsv"}
    assert not glob.glob(os.path.join(os.path.dirname(raw), "*.reader_diag.*"))
    d = json.load(open(out / "proj_X_s1.reader_diag.json"))
    assert d["outputs"]["ms2_csv"] == "ms2_proj_X_s1.csv" and d["outputs"]["ms3_csv"] == "ms3_proj_X_s1.csv"
    assert d["outputs"]["status_tsv"] == "proj_X_s1.scan_status.tsv"
    txt = (out / "proj_X_s1.reader_diag.txt").read_text()
    assert "ms2_proj_X_s1.csv" in txt and "ms2tmp_" not in txt and "'s1.scan_status.tsv'" not in txt


def test_promote_tolerates_missing_and_bad_json(ns, tmp_path):
    log = types.SimpleNamespace(lines=[])
    log.log = log.lines.append
    assert ns["promote_reader_diagnostics"]({}, str(tmp_path), "x", "a", "b", log) == {}
    bad = tmp_path / "in.reader_diag.json"
    bad.write_text("{}")
    out = tmp_path / "o"
    moved = ns["promote_reader_diagnostics"]({"json": str(bad)}, str(out), "x", "a", "b", log)
    assert os.path.isfile(moved["json"]) and any("references not updated" in l for l in log.lines)


# ---------------------------------------------------------------- headless MetadataEditorWindow

class _W:
    """Stands in for a Tk widget; after() runs the callback at once."""
    def __init__(self, value=""):
        self.value, self.cfg, self.destroyed = value, {}, False
    def winfo_exists(self): return 0 if self.destroyed else 1
    def after(self, ms, func, *a, **k): func(*a, **k)
    def config(self, **k): self.cfg.update(k)
    def title(self, *a): pass
    def destroy(self): self.destroyed = True
    def winfo_children(self): return []
    def delete(self, *a): self.value = ""
    def insert(self, index, text): self.value = text
    def get(self): return self.value


@pytest.fixture
def gui(ns, tmp_path, monkeypatch):
    """Fakes for dialogs, thread, status and reader; returns a builder for an editor over a list of RAW files."""
    calls = {"dialogs": [], "log": [], "reader_outdir": [], "reader_kwargs": []}
    mb = types.SimpleNamespace(
        showerror=lambda t, m: calls["dialogs"].append(("error", t, m)),
        showwarning=lambda t, m: calls["dialogs"].append(("warning", t, m)),
        showinfo=lambda t, m: calls["dialogs"].append(("info", t, m)),
        askyesno=lambda *a: True)
    out = tmp_path / "chosen"
    out.mkdir()
    monkeypatch.setitem(ns, "messagebox", mb)
    monkeypatch.setitem(ns, "simpledialog", types.SimpleNamespace(askstring=lambda *a, **k: "proj"))
    monkeypatch.setitem(ns, "filedialog", types.SimpleNamespace(askdirectory=lambda **k: str(out)))
    monkeypatch.setitem(ns, "threading", types.SimpleNamespace(
        Thread=lambda target, daemon=None: types.SimpleNamespace(start=target)))
    monkeypatch.setitem(ns, "set_main_status", lambda *a, **k: None)
    monkeypatch.setitem(ns, "progress", types.SimpleNamespace(stop=lambda: None))
    monkeypatch.setattr(ns["logger"], "gui_writer", calls["log"].append)
    monkeypatch.setattr(ns["logger"], "entries", [])

    def build(raws, behaviour, input_kind="raw", conversion_options=None):
        rname = "mzmlreader" if input_kind == "mzml" else "msprawextractor"   # the identity each reader writes

        def convert_raw_to_csv(raw, outdir=None, debug=False, **kw):
            calls["reader_outdir"].append(outdir)
            calls["reader_kwargs"].append(kw)
            folder = outdir or os.path.dirname(raw)
            kind = behaviour(raw) if callable(behaviour) else behaviour
            if kind == "false":
                _write_diag(raw, folder, failure=OSError("cannot open"), csvs=False, reader=rname)
                return False
            if kind == "false_no_diag":
                return False
            if kind == "raise":
                _write_diag(raw, folder, failure=RuntimeError("decoder crashed"), reader=rname)
                raise RuntimeError("decoder crashed")
            if kind == "failed_malformed":   # Codex V14F-F1: an accepted JSON whose "failure" is not an object
                res = _write_diag(raw, folder, failure=RuntimeError("x"), reader=rname)
                jp = res["diagnostics"]["json"]
                d = json.load(open(jp))
                d["failure"] = "bad"
                json.dump(d, open(jp, "w"))
                return res
            spec = {"complete": ("complete", 0, []), "warn": ("complete", 3, ["5 of 9 MS2 have no charge"]),
                    "partial": ("partial", 4, ["1 spectrum could not be decoded"]),
                    "failed": ("failed", 5, [])}[kind]
            return _write_diag(raw, folder, *spec, explanation="no MS2 with peaks" if kind == "failed" else "",
                               reader=rname)
        monkeypatch.setitem(ns, "mspext", types.SimpleNamespace(
            convert_raw_to_csv=convert_raw_to_csv, savemetadata=real_mspext.savemetadata,
            finalize_extraction=lambda *a, **k: None))
        monkeypatch.setitem(ns, "mspmzmlext", types.SimpleNamespace(
            extract_mzML=lambda path, outdir=None, debug=False, **kw: convert_raw_to_csv(path, outdir, debug, **kw)))
        cls = ns["MetadataEditorWindow"]
        ed = cls.__new__(cls)
        fields = ["Experiment Title", "Experiment Description", "Author running this analysis",
                  "Raw data acquired date", "Glycan Type", "Mass Analyzer charge mode", "Derivatization Type"]
        ed.__dict__.update(parent=_W(), window=_W(), status_label=_W(), raw_file_entry=_W(raws[0]),
                           entries={f: _W("v") for f in fields}, raw_file_list=list(raws), current_index=0,
                           callback=None, on_finish=None, metadata={}, last_metadata=None, output_dir=None,
                           skip_conversion=False, input_kind=input_kind, files={}, _convert_in_progress=False,
                           _error_mode=False, _batch_halted=False,
                           conversion_options=dict(conversion_options or ns["CONVERSION_OPTIONS_OFF"]))
        return ed

    calls["out"] = out
    return build, calls


def _finals(out):
    return sorted(os.path.basename(p) for p in glob.glob(os.path.join(out, "*")))


def test_gui_complete_promotes_csvs_and_diagnostics(gui, tmp_path):
    build, calls = gui
    raw = _raw(tmp_path)
    ed = build([raw], "complete")
    ed.load_file(raw)
    assert ed.status_label.cfg == {"text": "Conversion complete. Ready for metadata.", "fg": "green"}
    tmp_ms2 = open(ed.files["ms2tmp"]).read()
    ed.generate()
    names = _finals(calls["out"])
    stems = {n.split(".")[0] for n in names if n.endswith(".reader_diag.json")}
    assert len(stems) == 1
    save = stems.pop()
    assert {f"ms2_{save}.csv", f"ms3_{save}.csv", f"{save}.reader_diag.json", f"{save}.reader_diag.txt",
            f"{save}.scan_status.tsv", f"{save}.json"} <= set(names)
    assert open(calls["out"] / f"ms2_{save}.csv").read() == tmp_ms2     # CSV promoted unchanged
    assert not glob.glob(os.path.join(os.path.dirname(raw), "*"))[1:]   # only the RAW is left in its folder
    d = json.load(open(calls["out"] / f"{save}.reader_diag.json"))
    assert d["outputs"]["ms2_csv"] == f"ms2_{save}.csv" and d["outputs"]["status_tsv"] == f"{save}.scan_status.tsv"
    assert not [x for x in calls["dialogs"] if x[0] in ("error", "warning")]
    assert any("[reader] outcome: complete (exit code 0)" in l for l in calls["log"])


def test_gui_warnings_logged_no_dialog(gui, tmp_path):
    build, calls = gui
    raw = _raw(tmp_path)
    ed = build([raw], "warn")
    ed.load_file(raw)
    assert ed.status_label.cfg["fg"] == "orange"
    assert any("[reader] warning: 5 of 9 MS2 have no charge" in l for l in calls["log"])
    assert not [x for x in calls["dialogs"] if x[0] in ("error", "warning")]


def test_gui_partial_promotes_with_warning(gui, tmp_path):
    build, calls = gui
    raw = _raw(tmp_path)
    ed = build([raw], "partial")
    ed.load_file(raw)
    assert [x[1] for x in calls["dialogs"]] == ["Partial conversion"]
    ed.generate()
    assert any(n.startswith("ms2_proj_") for n in _finals(calls["out"]))
    assert any("[finalize] partial conversion promoted" in l for l in calls["log"])


@pytest.mark.parametrize("kind, title, needle", [
    ("failed", "Conversion failed", "no MS2 with peaks"),
    ("false", "Conversion failed", "cannot open"),
    ("false_no_diag", "Conversion failed", "no current diagnostics file: missing"),
    ("raise", "Thermo Library Error", "Reader diagnostics:"),
    ("failed_malformed", "Conversion failed", "bad"),
])
@pytest.mark.parametrize("input_kind", ["raw", "mzml"])
def test_gui_failure_never_promoted(gui, tmp_path, kind, title, needle, input_kind):
    build, calls = gui
    raw = _raw(tmp_path, "s1.mzML" if input_kind == "mzml" else "s1.raw")
    ed = build([raw], kind, input_kind=input_kind)
    ed.load_file(raw)
    errors = [x for x in calls["dialogs"] if x[0] == "error"]
    assert len(errors) == 1 and errors[0][1] == title and needle in errors[0][2]
    assert ed._error_mode and ed.window.destroyed
    ed.generate()                      # even if a user reached it, nothing is promoted
    assert _finals(calls["out"]) == []
    assert any("[reader] outcome: failed" in l for l in calls["log"])


def test_gui_stale_json_does_not_speak_for_a_failure(gui, tmp_path):
    build, calls = gui
    raw = _raw(tmp_path)
    _write_diag(raw, os.path.dirname(raw), "complete", 0, finished_utc="2020-01-01T00:00:00Z")
    ed = build([raw], "false_no_diag")
    ed.load_file(raw)
    err = [x for x in calls["dialogs"] if x[0] == "error"][0]
    assert "unchanged" in err[2] and ed._reader_result["paths"] == {}   # not rewritten by this call (V14F-F2)


def test_gui_batch_later_files_use_chosen_folder(gui, tmp_path):
    build, calls = gui
    raws = [_raw(tmp_path, "a.raw"), _raw(tmp_path, "b.raw")]
    ed = build(raws, "complete")
    ed.load_file(raws[0])
    ed.generate()                      # sets output_dir, then converts b.raw into it
    ed.generate()
    assert calls["reader_outdir"] == [None, str(calls["out"])]
    names = _finals(calls["out"])
    assert sum(n.endswith(".reader_diag.json") for n in names) == 2
    assert not any(n.startswith(("a.", "b.")) for n in names)   # renamed to the save names


# ---------------------------------------------------------------- Codex round 1 (V14F-F1, F2, F5, F7)

@pytest.mark.parametrize("bad", ["bad", ["a", "b"], 7])
def test_malformed_failure_field_keeps_a_failed_dict_failed(ns, tmp_path, bad):
    """V14F-F1: a string/list "failure" made the judgement raise before the outcome was copied -> promote=True."""
    raw = _raw(tmp_path)
    started = time.time()
    res = _write_diag(raw, os.path.dirname(raw), failure=RuntimeError("boom"))
    jp = res["diagnostics"]["json"]
    d = json.load(open(jp))
    d["failure"] = bad
    json.dump(d, open(jp, "w"))
    s = ns["summarize_reader_result"](res, None, raw, os.path.dirname(raw), started, "msprawextractor")
    assert s["outcome"] == "failed" and s["promote"] is False and s["dialog"] == "failed" and s["reason"]


def test_malformed_diagnostics_never_fail_a_complete_dict(ns, tmp_path):
    raw = _raw(tmp_path)
    started = time.time()
    res = _write_diag(raw, os.path.dirname(raw))
    jp = res["diagnostics"]["json"]
    d = json.load(open(jp))
    d.update(failure="bad", explanation=["x"])
    json.dump(d, open(jp, "w"))
    s = ns["summarize_reader_result"](res, None, raw, os.path.dirname(raw), started, "msprawextractor")
    assert s["outcome"] == "complete" and s["promote"] is True and s["dialog"] is None


def test_enrichment_error_is_dropped_and_the_outcome_stands(ns, tmp_path):
    """Anything that raises while reading the diagnostics only drops them (here: warnings is not a list)."""
    raw = _raw(tmp_path)
    started = time.time()
    _write_diag(raw, os.path.dirname(raw), failure=OSError("cannot open"), csvs=False)
    jp = os.path.join(os.path.dirname(raw), "s1.reader_diag.json")
    d = json.load(open(jp))
    d["warnings"] = 5
    json.dump(d, open(jp, "w"))
    s = ns["summarize_reader_result"](False, None, raw, os.path.dirname(raw), started, "msprawextractor")
    assert s["outcome"] == "failed" and not s["promote"] and s["diag_status"].startswith("error")
    assert s["paths"] == {} and "could not be opened" in s["reason"]
    res = {"ms2tmp": "a", "ms3tmp": "b", "outcome": "complete", "exit_code": 0, "warnings": 5, "diagnostics": {}}
    s = ns["summarize_reader_result"](res, None, raw, os.path.dirname(raw), started)
    assert s["promote"] and s["outcome"] == "complete"


def test_diagnostics_error_never_falls_back_to_an_old_file(ns, tmp_path):
    """V14F-F2: the reader could not write its files; the earlier run's JSON at the usual path is not read."""
    raw = _raw(tmp_path)
    _write_diag(raw, os.path.dirname(raw), explanation="OLD RUN")
    started = time.time()
    res = {"outcome": "failed", "exit_code": 5, "warnings": [], "diagnostics": {},
           "diagnostics_error": "OSError: disk full"}
    s = ns["summarize_reader_result"](res, None, raw, os.path.dirname(raw), started, "msprawextractor")
    assert s["diag_status"] == "not_written" and s["paths"] == {} and "OLD RUN" not in s["reason"]
    assert not s["promote"] and any("could not be written" in l for l in ns["reader_result_log_lines"](s))


def test_same_second_file_not_rewritten_is_unchanged(ns, tmp_path):
    """V14F-F2: an earlier run's JSON from the same second passes the time check but was not written by this call."""
    raw = _raw(tmp_path)
    out = os.path.dirname(raw)
    _write_diag(raw, out, explanation="OLD RUN")
    started = time.time()
    before = ns["reader_diag_fingerprint"](raw, out)
    assert before
    s = ns["summarize_reader_result"](False, None, raw, out, started, "msprawextractor", before=before)
    assert s["diag_status"] == "unchanged" and "OLD RUN" not in s["reason"] and s["paths"] == {}
    s = ns["summarize_reader_result"](None, RuntimeError("x"), raw, out, started, "msprawextractor", before=before)
    assert s["diag_status"] == "unchanged" and "OLD RUN" not in s["reason"]


def test_rewritten_file_is_fresh_with_a_fingerprint(ns, tmp_path):
    raw = _raw(tmp_path)
    out = os.path.dirname(raw)
    _write_diag(raw, out, explanation="OLD RUN")
    started = time.time()
    before = ns["reader_diag_fingerprint"](raw, out)
    _write_diag(raw, out, failure=OSError("cannot open"), csvs=False)        # this call's file
    s = ns["summarize_reader_result"](False, None, raw, out, started, "msprawextractor", before=before)
    assert s["diag_status"] == "fresh" and "OLD RUN" not in s["reason"] and "cannot open" in s["reason"]


def test_previous_second_is_stale(ns, tmp_path):
    """V14F-F2: the reader stamps whole seconds with the same clock, so nothing before int(started) is current."""
    raw = _raw(tmp_path)
    started = time.time()
    prev = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(int(started) - 1))
    same = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(int(started)))
    load = ns["load_reader_diagnostics"]
    _write_diag(raw, os.path.dirname(raw), finished_utc=prev)
    assert load(raw, os.path.dirname(raw), started, "msprawextractor")[1] == "stale"
    _write_diag(raw, os.path.dirname(raw), finished_utc=same)
    assert load(raw, os.path.dirname(raw), started, "msprawextractor")[1] == "fresh"


def test_fingerprint_and_load_share_the_stem_for_mzml_gz(ns, tmp_path):
    raw = _raw(tmp_path, "s1.mzML.gz")
    out = os.path.dirname(raw)
    assert ns["reader_diag_fingerprint"](raw, out) is None
    started = time.time()
    _write_diag(raw, out, reader="mzmlreader")
    assert os.path.isfile(os.path.join(out, "s1.mzML.reader_diag.json"))
    assert ns["reader_diag_fingerprint"](raw, out)
    assert ns["load_reader_diagnostics"](raw, out, started, "mzmlreader")[1] == "fresh"


def test_gui_old_same_second_file_does_not_speak_for_a_failure(gui, tmp_path):
    build, calls = gui
    raw = _raw(tmp_path)
    _write_diag(raw, os.path.dirname(raw), explanation="OLD RUN")       # same second, not rewritten by the call
    ed = build([raw], "false_no_diag")
    ed.load_file(raw)
    err = [x for x in calls["dialogs"] if x[0] == "error"][0]
    assert "OLD RUN" not in err[2] and "unchanged" in err[2]


def test_generate_guard_blocks_a_failed_result_with_temps_present(gui, tmp_path):
    """V14F-F7: the explicit failed-result guard in generate(), reached with real temp CSVs on disk."""
    build, calls = gui
    raw = _raw(tmp_path)
    ed = build([raw], "complete")
    ed.load_file(raw)
    assert os.path.isfile(ed.files["ms2tmp"]) and os.path.isfile(ed.files["ms3tmp"])
    ed._reader_result = dict(ed._reader_result, outcome="failed", exit_code=5, promote=False,
                             reason="The reader reported a failed conversion.")
    ed.generate()
    assert _finals(calls["out"]) == []
    assert [x[1] for x in calls["dialogs"] if x[0] == "error"] == ["Conversion failed"]
    assert any("[generate] blocked: the reader reported a failed conversion" in l for l in calls["log"])
    assert os.path.isfile(ed.files["ms2tmp"])          # left where the reader wrote it


def test_gui_mzml_diagnostics_are_accepted_and_promoted(gui, tmp_path):
    build, calls = gui
    raw = _raw(tmp_path, "s1.mzML")
    ed = build([raw], "complete", input_kind="mzml")
    ed.load_file(raw)
    assert ed._reader_result["diag_status"] == "fresh"
    ed.generate()
    names = _finals(calls["out"])
    assert sum(n.endswith(".reader_diag.json") for n in names) == 1
    assert sum(n.endswith(".scan_status.tsv") for n in names) == 1


def test_gui_mzml_rejects_raw_reader_diagnostics(ns, tmp_path):
    raw = _raw(tmp_path, "s1.mzML")
    started = time.time()
    res = _write_diag(raw, os.path.dirname(raw), reader="msprawextractor")
    s = ns["summarize_reader_result"](res, None, raw, os.path.dirname(raw), started, "mzmlreader")
    assert s["diag_status"] == "wrong_type" and s["promote"] and s["paths"] == {}


def test_gui_batch_later_failure_after_a_success(gui, tmp_path):
    build, calls = gui
    raws = [_raw(tmp_path, "a.raw"), _raw(tmp_path, "b.raw")]
    ed = build(raws, lambda raw: "failed" if raw.endswith("b.raw") else "complete")
    ed.load_file(raws[0])
    ed.generate()                      # a.raw promoted; b.raw then fails in the chosen folder
    names = _finals(calls["out"])
    assert sum(n.startswith("ms2_proj_") for n in names) == 1
    assert [x[1] for x in calls["dialogs"] if x[0] == "error"] == ["Conversion failed"]
    ed.generate()
    assert sum(n.startswith("ms2_proj_") for n in _finals(calls["out"])) == 1


class _Log:
    def __init__(self):
        self.lines = []
    def log(self, m):
        self.lines.append(m)


def test_promote_exdev_falls_back_to_copy(ns, tmp_path, monkeypatch):
    """V14F-F5: on macOS/Linux two volumes share the empty drive, and os.replace raises EXDEV."""
    import errno
    src, dst = tmp_path / "ms2tmp_s1.csv", tmp_path / "out" / "ms2_s1.csv"
    src.write_text("data\n")
    def exdev(a, b):
        raise OSError(errno.EXDEV, "Invalid cross-device link")
    monkeypatch.setattr(os, "replace", exdev)
    log = _Log()
    ns["promote_temp_to_final"](str(src), str(dst), log)
    assert dst.read_text() == "data\n" and not src.exists()
    assert "(other filesystem)" in log.lines[-1]


def test_promote_exdev_copy_failure_keeps_the_source(ns, tmp_path, monkeypatch):
    import errno
    import shutil
    src, dst = tmp_path / "ms2tmp_s1.csv", tmp_path / "out" / "ms2_s1.csv"
    src.write_text("data\n")
    copied = []
    def full(a, b):
        copied.append((a, b))
        raise OSError(errno.ENOSPC, "disk full")
    monkeypatch.setattr(os, "replace", lambda a, b: (_ for _ in ()).throw(OSError(errno.EXDEV, "x")))
    monkeypatch.setattr(shutil, "copy2", full)
    with pytest.raises(OSError) as e:
        ns["promote_temp_to_final"](str(src), str(dst), _Log())
    assert e.value.errno == errno.ENOSPC and len(copied) == 1    # the copy was tried (V14F-R2-F2)
    assert src.read_text() == "data\n"


def test_promote_other_errors_still_raise(ns, tmp_path, monkeypatch):
    import errno
    src, dst = tmp_path / "ms2tmp_s1.csv", tmp_path / "out" / "ms2_s1.csv"
    src.write_text("data\n")
    monkeypatch.setattr(os, "replace", lambda a, b: (_ for _ in ()).throw(PermissionError(errno.EACCES, "denied")))
    with pytest.raises(PermissionError):
        ns["promote_temp_to_final"](str(src), str(dst), _Log())
    assert src.exists() and not dst.exists()


def test_promote_same_filesystem_is_unchanged(ns, tmp_path):
    src, dst = tmp_path / "ms2tmp_s1.csv", tmp_path / "out" / "ms2_s1.csv"
    src.write_text("data\n")
    log = _Log()
    ns["promote_temp_to_final"](str(src), str(dst), log)
    assert dst.read_text() == "data\n" and not src.exists() and "Promoted (atomic)" in log.lines[-1]


# ---------------------------------------------------------------- Codex round 2 (V14F-R2-F1)

def _unreadable_fingerprint(ns, monkeypatch, raw, out):
    """reader_diag_fingerprint while reading the existing file fails (e.g. a transient network-drive error)."""
    def boom(*a, **k):
        raise OSError("transient read error")
    monkeypatch.setitem(ns, "open", boom)
    try:
        return ns["reader_diag_fingerprint"](raw, out)
    finally:
        monkeypatch.delitem(ns, "open")


def test_fingerprint_three_states(ns, tmp_path, monkeypatch):
    raw = _raw(tmp_path)
    out = os.path.dirname(raw)
    assert ns["reader_diag_fingerprint"](raw, out) is None                      # no file
    _write_diag(raw, out)
    assert len(ns["reader_diag_fingerprint"](raw, out)) == 40                    # sha1
    assert _unreadable_fingerprint(ns, monkeypatch, raw, out) == ns["DIAG_FINGERPRINT_UNREADABLE"]
    assert len(ns["reader_diag_fingerprint"](raw, out)) == 40                    # the later read works again


@pytest.mark.parametrize("res, exc", [(False, None), (None, RuntimeError("decoder crashed")), (("a", "b"), None)])
def test_unreadable_before_is_unverifiable(ns, tmp_path, monkeypatch, res, exc):
    """V14F-R2-F1: before, a failed pre-call read gave None and skipped the comparison, so a same-second old
    file was read as fresh. Now freshness cannot be established and the file is not used."""
    raw = _raw(tmp_path)
    out = os.path.dirname(raw)
    _write_diag(raw, out, explanation="OLD RUN")
    started = time.time()
    before = _unreadable_fingerprint(ns, monkeypatch, raw, out)
    s = ns["summarize_reader_result"](res, exc, raw, out, started, "msprawextractor", before=before)
    assert s["diag_status"] == "unverifiable" and s["paths"] == {} and s["diag_txt"] is None
    assert "OLD RUN" not in s["reason"] and s["diag_reason"] == ""
    if res is False or exc is not None:
        assert s["outcome"] == "failed" and not s["promote"]
    else:
        assert s["outcome"] is None and s["promote"]                             # legacy tuple: as before


def test_named_paths_ignore_the_fingerprint(ns, tmp_path):
    """A dict names its own files; the pre-call state of the usual path does not apply to it."""
    raw = _raw(tmp_path)
    out = os.path.dirname(raw)
    started = time.time()
    res = _write_diag(raw, out)
    s = ns["summarize_reader_result"](res, None, raw, out, started, "msprawextractor",
                                      before=ns["DIAG_FINGERPRINT_UNREADABLE"])
    assert s["diag_status"] == "fresh" and set(s["paths"]) == {"json", "txt", "status_tsv"}


def test_identical_rewrite_is_unchanged_and_the_outcome_stands(ns, tmp_path):
    """Accepted policy: a byte-identical rewrite cannot be shown to be this call's; only the diagnostics are skipped."""
    raw = _raw(tmp_path)
    out = os.path.dirname(raw)
    _write_diag(raw, out, failure=OSError("cannot open"), csvs=False)
    jp = os.path.join(out, "s1.reader_diag.json")
    content = open(jp, "rb").read()
    started = time.time()
    before = ns["reader_diag_fingerprint"](raw, out)
    open(jp, "wb").write(content)                                               # the call writes the same bytes
    s = ns["summarize_reader_result"](False, None, raw, out, started, "msprawextractor", before=before)
    assert s["diag_status"] == "unchanged" and s["outcome"] == "failed" and not s["promote"]
    assert "could not be opened" in s["reason"]


def test_gui_unreadable_fingerprint_does_not_speak_for_a_failure(gui, ns, tmp_path, monkeypatch):
    build, calls = gui
    raw = _raw(tmp_path)
    _write_diag(raw, os.path.dirname(raw), explanation="OLD RUN")
    monkeypatch.setitem(ns, "reader_diag_fingerprint", lambda *a: ns["DIAG_FINGERPRINT_UNREADABLE"])
    ed = build([raw], "false_no_diag")
    ed.load_file(raw)
    err = [x for x in calls["dialogs"] if x[0] == "error"][0]
    assert "OLD RUN" not in err[2] and "unverifiable" in err[2]
    assert _finals(calls["out"]) == []


# ---------------------------------------------------------------- Codex round 3 (V14F-R3-F1)

def _lstat_failing_fingerprint(ns, monkeypatch, raw, out, err):
    """reader_diag_fingerprint while lstat of the existing file fails; access is restored afterwards."""
    import errno as _errno
    real = os.lstat
    target = os.path.join(out, "s1.reader_diag.json")
    def lstat(path, *a, **k):
        if os.fspath(path) == target:
            raise OSError(getattr(_errno, err), err)
        return real(path, *a, **k)
    with monkeypatch.context() as m:
        m.setattr(os, "lstat", lstat)
        return ns["reader_diag_fingerprint"](raw, out)


@pytest.mark.parametrize("err", ["EACCES", "EIO"])
@pytest.mark.parametrize("res, exc", [(False, None), (None, RuntimeError("decoder crashed"))])
def test_lstat_error_is_not_absence(ns, tmp_path, monkeypatch, err, res, exc):
    """V14F-R3-F1: os.path.lexists turned EACCES/EIO into "absent", so before=None and a same-second old file
    was read as fresh once access recovered."""
    raw = _raw(tmp_path)
    out = os.path.dirname(raw)
    _write_diag(raw, out, explanation="OLD RUN")
    started = time.time()
    before = _lstat_failing_fingerprint(ns, monkeypatch, raw, out, err)
    assert before == ns["DIAG_FINGERPRINT_UNREADABLE"]
    s = ns["summarize_reader_result"](res, exc, raw, out, started, "msprawextractor", before=before)
    assert s["diag_status"] == "unverifiable" and s["paths"] == {} and "OLD RUN" not in s["reason"]
    assert s["outcome"] == "failed" and not s["promote"]


def test_lstat_file_not_found_is_absence(ns, tmp_path, monkeypatch):
    raw = _raw(tmp_path)
    out = os.path.dirname(raw)
    with monkeypatch.context() as m:
        m.setattr(os, "lstat", lambda p, *a, **k: (_ for _ in ()).throw(FileNotFoundError(2, "gone")))
        assert ns["reader_diag_fingerprint"](raw, out) is None


def test_broken_link_is_unreadable(ns, tmp_path):
    raw = _raw(tmp_path)
    out = os.path.dirname(raw)
    try:
        os.symlink(os.path.join(out, "nowhere.json"), os.path.join(out, "s1.reader_diag.json"))
    except (OSError, NotImplementedError) as e:
        pytest.skip(f"no symlinks here: {e}")
    assert ns["reader_diag_fingerprint"](raw, out) == ns["DIAG_FINGERPRINT_UNREADABLE"]
