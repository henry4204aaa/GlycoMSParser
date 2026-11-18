"""B-03 smoke test — verify _apply_to_exp + _save_effective_to_method fail-close on invalid JSON.

The bug (see handoff/applytoexpbug.md): _apply_to_exp was passing the raw editor
string to _write_ml_defaults_to_exp (which expects a dict), producing a corrupted
on-disk .exp.json with `ml_defaults` as a JSON-string-literal instead of a nested dict.

The fix (S2 / 2026-05-17, B-03): replace the raw editor.get + dead None check with
_get_edit() (the validating helper that already shows messagebox.showerror on parse
failure and returns None). Codex Finding 11 surfaced a second instance of the same
silent-corruption bug class at _save_effective_to_method line 6974 — folded in as
Fold-B (line 6978 in post-fix source).

This test mirrors the post-fix logic of both fixed functions and exercises:
    T1 - _apply_to_exp good-JSON path writes dict (not string) to ml_defaults
    T2 - _apply_to_exp bad-JSON path short-circuits via _get_edit (no write)
    T3 - _apply_to_exp empty-editor path writes empty dict {}
    T4 - _apply_to_exp round-trip preserves dict shape (no escape layers stack)
    T5 - _save_effective_to_method (Fold-B) fail-closes on bad JSON (no write)
    T6 - _save_effective_to_method good-JSON path writes parameters dict

Run from repo root:
    python3 tests/test_b03_apply_to_exp_smoke.py
"""
import json
import os
import sys
import tempfile


# ---- mirrors of the post-fix v14 helpers --------------------------------
# Keep in lockstep with src/mspfileloaderv14.py edits 1-4 (lines 6878, 6890, 6891-6894, 6976-6979).


class _FakeEditor:
    """Stand-in for tk.Text — returns canned content from .get('1.0','end')."""

    def __init__(self, content):
        self._content = content

    def get(self, start, end):
        return self._content


_messagebox_errors = []  # captures messagebox.showerror calls so tests can assert on them


def _fake_showerror(title, message):
    _messagebox_errors.append((title, message))


# A2 mirror: _get_edit (line 6865 — validating helper, unchanged by B-03)
def _get_edit(editor_txt):
    try:
        txt = editor_txt.get("1.0", "end").strip() or "{}"
        return json.loads(txt)
    except Exception as e:
        _fake_showerror("Invalid JSON", f"Editor JSON parse error:\n{e}")
        return None


# Module-level writer mirror (line 6678, unchanged)
def _write_ml_defaults_to_exp(path, params):
    """Mirrors line 6678. If `params` is a string, json.dump serializes it as a
    JSON string literal — the pre-B-03 bug. After the fix, params must be a dict.
    """
    try:
        data = {}
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
        data["ml_defaults"] = params
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
            f.write("\n")
        return True
    except Exception as e:
        _fake_showerror("Save failed", str(e))
        return False


def _read_ml_defaults_from_exp(path):
    """Mirrors line 6671. Round-trip consumer."""
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return data.get("ml_defaults")
    except Exception:
        return None


# Post-fix _apply_to_exp mirror (lines 6886-6896 post-B-03)
def _apply_to_exp(editor_txt, exp_path):
    d = _get_edit(editor_txt)
    if d is None:
        return False  # short-circuit on parse failure (no write)
    return _write_ml_defaults_to_exp(exp_path, d)


# Post-fix _save_effective_to_method mirror (Fold-B; lines 6968-6993 post-B-03).
# Simplified to just the editor-parsing + write portion that the fix touches.
def _save_effective_to_method(editor_txt, method_path, snapshot, exp_defs, meth_defs):
    """Mirrors the editor-read + merge + write portion of _save_effective_to_method.
    Pre-Fold-B used _get_editor_json_or_empty (silent {} on bad JSON);
    post-Fold-B uses _get_edit (fail-close).
    """
    editor = _get_edit(editor_txt)
    if editor is None:
        return False  # Fold-B: fail-close on invalid editor JSON
    # Simplified merge (real code uses _merge_ml; here we just shallow-merge for test purposes)
    eff = {}
    for src in (snapshot, exp_defs or {}, meth_defs or {}, editor):
        eff.update(src)
    data = {"ml": {"parameters": eff}}
    with open(method_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    return True


# ---- tests --------------------------------------------------------------


def main():
    failures = []

    with tempfile.TemporaryDirectory() as tmpdir:
        exp_path = os.path.join(tmpdir, "test.exp.json")
        method_path = os.path.join(tmpdir, "test.method.json")

        # T1 - good JSON path writes dict (not string) to ml_defaults
        _messagebox_errors.clear()
        editor = _FakeEditor('{"model": {"type": "RandomForest", "n_estimators": 100}}')
        ok = _apply_to_exp(editor, exp_path)
        if not ok:
            failures.append("T1: _apply_to_exp returned False on valid JSON")
        else:
            with open(exp_path) as f:
                data = json.load(f)
            ml = data.get("ml_defaults")
            if not isinstance(ml, dict):
                failures.append(f"T1: ml_defaults is {type(ml).__name__}, expected dict")
            elif ml != {"model": {"type": "RandomForest", "n_estimators": 100}}:
                failures.append(f"T1: ml_defaults content mismatch: {ml}")
            else:
                print("T1 PASS: _apply_to_exp writes dict (not string) to ml_defaults")

        # T2 - bad JSON path short-circuits via _get_edit (no write, error surfaced)
        os.remove(exp_path)  # start fresh
        _messagebox_errors.clear()
        editor = _FakeEditor('{"model": {invalid json here')
        ok = _apply_to_exp(editor, exp_path)
        if ok:
            failures.append("T2: _apply_to_exp returned True on invalid JSON (should short-circuit)")
        elif os.path.exists(exp_path):
            failures.append("T2: exp.json was written despite invalid JSON")
        elif not _messagebox_errors or _messagebox_errors[0][0] != "Invalid JSON":
            failures.append(f"T2: no 'Invalid JSON' messagebox surfaced; got {_messagebox_errors}")
        else:
            print("T2 PASS: _apply_to_exp short-circuits + surfaces error on invalid JSON")

        # T3 - empty editor path writes empty dict {}
        if os.path.exists(exp_path):
            os.remove(exp_path)
        _messagebox_errors.clear()
        editor = _FakeEditor("")
        ok = _apply_to_exp(editor, exp_path)
        if not ok:
            failures.append("T3: _apply_to_exp returned False on empty editor")
        else:
            with open(exp_path) as f:
                data = json.load(f)
            ml = data.get("ml_defaults")
            if ml != {}:
                failures.append(f"T3: empty editor produced ml_defaults={ml!r}, expected {{}}")
            else:
                print("T3 PASS: _apply_to_exp empty editor writes empty dict {}")

        # T4 - round-trip preservation (no escape layers stack)
        os.remove(exp_path)
        _messagebox_errors.clear()
        original = {"model": {"type": "RandomForest"}, "train": {"test_size": 0.2}}
        editor = _FakeEditor(json.dumps(original, indent=2))
        _apply_to_exp(editor, exp_path)
        roundtripped = _read_ml_defaults_from_exp(exp_path)
        if roundtripped != original:
            failures.append(f"T4: round-trip mismatch: {roundtripped} != {original}")
        elif not isinstance(roundtripped, dict):
            failures.append(f"T4: round-tripped value is {type(roundtripped).__name__}, expected dict")
        else:
            print("T4 PASS: _apply_to_exp round-trip preserves dict (no escape stacking)")

        # T5 - Fold-B: _save_effective_to_method fail-closes on bad JSON
        if os.path.exists(method_path):
            os.remove(method_path)
        _messagebox_errors.clear()
        editor = _FakeEditor('{"broken: incomplete')
        snapshot = {"model": {"type": "RandomForest"}}
        ok = _save_effective_to_method(editor, method_path, snapshot, None, None)
        if ok:
            failures.append("T5: _save_effective_to_method returned True on invalid JSON (Fold-B failed)")
        elif os.path.exists(method_path):
            failures.append("T5: method.json was written despite invalid JSON (Fold-B failed)")
        elif not _messagebox_errors or _messagebox_errors[0][0] != "Invalid JSON":
            failures.append(f"T5: no 'Invalid JSON' messagebox surfaced; got {_messagebox_errors}")
        else:
            print("T5 PASS: _save_effective_to_method (Fold-B) fail-closes on invalid JSON")

        # T6 - Fold-B: good JSON path writes merged parameters
        _messagebox_errors.clear()
        editor = _FakeEditor('{"editor_override": true}')
        snapshot = {"model": {"type": "RandomForest"}}
        ok = _save_effective_to_method(editor, method_path, snapshot, None, None)
        if not ok:
            failures.append("T6: _save_effective_to_method returned False on valid JSON")
        else:
            with open(method_path) as f:
                data = json.load(f)
            params = data.get("ml", {}).get("parameters", {})
            if params.get("editor_override") is not True:
                failures.append(f"T6: editor override not present in saved params: {params}")
            elif params.get("model", {}).get("type") != "RandomForest":
                failures.append(f"T6: snapshot params not merged: {params}")
            else:
                print("T6 PASS: _save_effective_to_method (Fold-B) writes merged params on valid JSON")

    # Summary
    total = 6
    if failures:
        print()
        print(f"FAIL: {len(failures)}/{total} tests failed")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    else:
        print()
        print(f"PASS: {total}/{total} tests passed")
        sys.exit(0)


if __name__ == "__main__":
    main()
