"""V14-FOLLOWUP CP-C: precursor-charge options for GUI conversion (ticket WP-B).

The "Conversion options" dialog (text approved by Henry 2026-10-27) runs before the metadata window. Both options
start off; with both off the readers are called with exactly the v1.10 arguments. The choice is passed to the
readers as msp_extract_cli does (charge_comb=, dataset_charge=) and recorded in the project metadata JSON under two
additive keys. Uses the headless v14 namespace and the fake editor of tests/test_v14_reader_result.py.
"""
import glob
import inspect
import json
import os
import types

import pytest

from test_v14_reader_result import _raw, gui, ns  # noqa: F401  (fixtures; also puts src/ on sys.path)

import msprawextractor  # noqa: E402
import mzmlreader  # noqa: E402


@pytest.mark.parametrize("text, value", [("", None), ("   ", None), ("1", 1), (" 2 ", 2), ("-2", 2), ("+3", 3)])
def test_parse_dataset_charge_ok(ns, text, value):
    assert ns["parse_dataset_charge_text"](text) == value


@pytest.mark.parametrize("text", ["0", "1.5", "abc", "1e2", "−1"])
def test_parse_dataset_charge_rejects(ns, text):
    with pytest.raises(ValueError) as e:
        ns["parse_dataset_charge_text"](text)
    assert str(e.value) == ns["DATASET_CHARGE_ERROR"]


def test_reader_kwargs(ns):
    kw = ns["conversion_reader_kwargs"]
    assert kw(None) == {} and kw(ns["CONVERSION_OPTIONS_OFF"]) == {}
    assert kw({"charge_comb": True, "dataset_charge": None}) == {"charge_comb": True}
    assert kw({"charge_comb": False, "dataset_charge": 1}) == {"dataset_charge": 1}
    assert kw({"charge_comb": True, "dataset_charge": 2}) == {"charge_comb": True, "dataset_charge": 2}


def test_reader_kwargs_match_both_readers(ns):
    names = set(ns["conversion_reader_kwargs"]({"charge_comb": True, "dataset_charge": 1}))
    for fn in (msprawextractor.convert_raw_to_csv, mzmlreader.extract_mzML):
        assert names <= set(inspect.signature(fn).parameters)


def test_metadata_entries(ns):
    me = ns["conversion_metadata_entries"]
    assert me(None) == {"Conversion charge comb": False, "Conversion dataset charge": None}
    assert me({"charge_comb": True, "dataset_charge": 1}) == {"Conversion charge comb": True,
                                                               "Conversion dataset charge": 1}


def _metadata(out):
    paths = [p for p in glob.glob(os.path.join(out, "*.json")) if not p.endswith(".reader_diag.json")]
    assert len(paths) == 1
    return json.load(open(paths[0]))


@pytest.mark.parametrize("kind", ["raw", "mzml"])
def test_options_off_call_readers_as_before(gui, tmp_path, kind):
    build, calls = gui
    raw = _raw(tmp_path, "s1.mzML" if kind == "mzml" else "s1.raw")
    ed = build([raw], "complete", input_kind=kind)
    ed.load_file(raw)
    ed.generate()
    assert calls["reader_kwargs"] == [{}]
    md = _metadata(calls["out"])
    assert md["Conversion charge comb"] is False and md["Conversion dataset charge"] is None
    assert len(glob.glob(os.path.join(calls["out"], "*.reader_diag.json"))) == 1   # this reader's file, accepted


@pytest.mark.parametrize("kind", ["raw", "mzml"])
def test_options_on_reach_reader_and_metadata(gui, tmp_path, kind):
    build, calls = gui
    raw = _raw(tmp_path, "s1.mzML" if kind == "mzml" else "s1.raw")
    ed = build([raw], "complete", input_kind=kind, conversion_options={"charge_comb": True, "dataset_charge": 1})
    ed.load_file(raw)
    ed.generate()
    assert calls["reader_kwargs"] == [{"charge_comb": True, "dataset_charge": 1}]
    md = _metadata(calls["out"])
    assert md["Conversion charge comb"] is True and md["Conversion dataset charge"] == 1
    assert len(glob.glob(os.path.join(calls["out"], "*.reader_diag.json"))) == 1


def test_options_apply_to_every_file_in_batch(gui, tmp_path):
    build, calls = gui
    raws = [_raw(tmp_path, "a.raw"), _raw(tmp_path, "b.raw")]
    ed = build(raws, "complete", conversion_options={"charge_comb": True, "dataset_charge": None})
    ed.load_file(raws[0])
    ed.generate()
    ed.generate()
    assert calls["reader_kwargs"] == [{"charge_comb": True}, {"charge_comb": True}]


@pytest.fixture
def launcher(ns, monkeypatch):
    """launch_metadata_batch with the main-window globals and the editor replaced by recorders."""
    rec = {"editor": [], "reset": 0, "status": []}
    monkeypatch.setitem(ns, "root", object())
    monkeypatch.setitem(ns, "status_var", types.SimpleNamespace(set=rec["status"].append))
    monkeypatch.setitem(ns, "progress", types.SimpleNamespace(start=lambda: None, stop=lambda: None))
    monkeypatch.setitem(ns, "reset_main_status", lambda: rec.__setitem__("reset", rec["reset"] + 1))
    monkeypatch.setitem(ns, "MetadataEditorWindow", lambda *a, **k: rec["editor"].append(k))
    monkeypatch.setitem(ns, "mspext", types.SimpleNamespace(pymsreader=True))
    return rec


@pytest.mark.parametrize("files", [{"raw": ["x.raw"]}, {"mzml": ["x.mzML"]}])
def test_launch_cancel_starts_nothing(ns, launcher, monkeypatch, files):
    monkeypatch.setitem(ns, "selected_files", files)
    monkeypatch.setitem(ns, "ask_conversion_options", lambda parent: None)
    ns["launch_metadata_batch"]()
    assert launcher["editor"] == [] and launcher["reset"] == 1 and launcher["status"] == []


@pytest.mark.parametrize("files", [{"raw": ["x.raw"]}, {"mzml": ["x.mzML"]}])
def test_launch_passes_options(ns, launcher, monkeypatch, files):
    opts = {"charge_comb": True, "dataset_charge": 2}
    monkeypatch.setitem(ns, "selected_files", files)
    monkeypatch.setitem(ns, "ask_conversion_options", lambda parent: dict(opts))
    ns["launch_metadata_batch"]()
    assert len(launcher["editor"]) == 1 and launcher["editor"][0]["conversion_options"] == opts


def test_skip_conversion_metadata_has_no_option_keys(gui, tmp_path):
    build, calls = gui
    raw = _raw(tmp_path)
    ed = build([raw], "complete")
    ed.skip_conversion = True
    ed.output_dir = str(calls["out"])
    ed.generate()
    md = _metadata(calls["out"])
    assert "Conversion charge comb" not in md and "Conversion dataset charge" not in md


# ---------------------------------------------------------------- the dialog itself (real Tk, skipped headless)

@pytest.fixture
def tkroot(ns):
    tk = ns["tk"]
    try:
        root = tk.Tk()
    except tk.TclError as e:
        pytest.skip(f"no display for Tk: {e}")
    root.withdraw()
    yield root
    root.destroy()


def _walk(w):
    for c in w.winfo_children():
        yield c
        yield from _walk(c)


def _drive(ns, root, monkeypatch, actions):
    """Run ask_conversion_options; `actions(dlg, widgets)` stands in for the user while the dialog waits."""
    errors = []
    monkeypatch.setitem(ns, "messagebox", types.SimpleNamespace(showerror=lambda t, m, **k: errors.append(m)))
    tk = ns["tk"]

    def fake_wait(dlg):
        ws = list(_walk(dlg))
        actions(dlg, {"entry": next(w for w in ws if isinstance(w, tk.Entry)),
                      "check": next(w for w in ws if isinstance(w, tk.Checkbutton)),
                      "buttons": {w.cget("text"): w for w in ws if isinstance(w, tk.Button)},
                      "labels": [w.cget("text") for w in ws if isinstance(w, tk.Label)]})
    monkeypatch.setattr(root, "wait_window", fake_wait)
    return ns["ask_conversion_options"](root), errors


def test_dialog_defaults_off(ns, tkroot, monkeypatch):
    seen = {}
    def act(dlg, w):
        seen.update(check=w["check"].cget("text"), labels=w["labels"])
        w["buttons"]["Convert"].invoke()
    res, errors = _drive(ns, tkroot, monkeypatch, act)
    assert res == {"charge_comb": False, "dataset_charge": None} and errors == []
    assert seen["check"] == "Calculate the charge from the MS1 isotope pattern (charge 1–3)"
    assert any(t.startswith("Dataset charge for the rest (blank = off):") for t in seen["labels"])


def test_dialog_values_invalid_then_cancel(ns, tkroot, monkeypatch):
    def act(dlg, w):
        w["check"].select()
        w["entry"].insert(0, "1")
        w["buttons"]["Convert"].invoke()
    res, _ = _drive(ns, tkroot, monkeypatch, act)
    assert res == {"charge_comb": True, "dataset_charge": 1}

    def act_bad(dlg, w):
        w["entry"].insert(0, "0")
        w["buttons"]["Convert"].invoke()        # refused: dialog stays open
        assert dlg.winfo_exists()
        w["buttons"]["Cancel"].invoke()
    res, errors = _drive(ns, tkroot, monkeypatch, act_bad)
    assert res is None and errors == [ns["DATASET_CHARGE_ERROR"]]
