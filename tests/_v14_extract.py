"""Headless access to v14 closure-local functions (CGA-TOLERANCE, 2026-10-31).

mspfileloaderv14.py creates the Tk root at module level and keeps the CGA functions (`run_pseudolabeling`,
`build_trainable_from_CGA`, ...) nested inside `open_prepare_dataset_window`, so none of them can be imported.
Same technique as tests/test_s1_v14_headless.py: execute the module source up to (not including)
`root = tk.Tk()`, which defines every module-level name, then extract the requested nested defs with `ast`,
dedent them and execute them in that namespace. Closure names the extracted code needs (logger, Tk variables,
append_runlog, refresh_tree, subwin, ...) are supplied by the caller as `stubs`.

`source` may be any version of the file (e.g. `git show <base>:src/mspfileloaderv14.py`), so a test or a
harness can run the base build and the branch build on the same inputs and compare the outputs byte for byte.
"""
import ast
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
V14_PATH = os.path.join(SRC, "mspfileloaderv14.py")


class AnchorMissing(LookupError):
    """A requested anchor is gone (moved/renamed): the caller decides whether that is a skip."""


def _nested_defs(tree, names):
    found = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name in names and node.col_offset > 0 \
                and node.name not in found:
            found[node.name] = node
    return found


def load_v14(source=None, *, extract=(), stubs=None, v14_path=V14_PATH):
    """Return a namespace with v14's module-level names plus the `extract`ed nested defs.

    Raises AnchorMissing if `root = tk.Tk()` or a requested def is not found. Any other error propagates."""
    if source is None:
        with open(v14_path, "rb") as f:
            source = f.read().decode("utf-8")
    lines = source.splitlines(keepends=True)
    root_i = next((i for i, l in enumerate(lines) if l.startswith("root = tk.Tk()")), None)
    if root_i is None:
        raise AnchorMissing("root = tk.Tk()")
    tree = ast.parse(source)
    defs = _nested_defs(tree, set(extract))
    missing = [n for n in extract if n not in defs]
    if missing:
        raise AnchorMissing(", ".join(missing))

    if SRC not in sys.path:
        sys.path.insert(0, SRC)
    ns = {"__name__": "v14_headless", "__file__": v14_path}
    old_hook, old_cwd = sys.excepthook, os.getcwd()
    try:
        os.chdir(SRC)
        exec(compile("".join(lines[:root_i]), "v14_headless", "exec"), ns)
        ns.update(stubs or {})
        for name in extract:
            node = defs[name]
            start = node.lineno - 1 - len(node.decorator_list)
            seg = lines[start:node.end_lineno]
            pad = node.col_offset
            body = "".join(l[pad:] if l[:pad].strip() == "" else l for l in seg)
            exec(compile(body, f"v14_extract:{name}", "exec"), ns)
    finally:
        sys.excepthook = old_hook
        os.chdir(old_cwd)
    return ns


class Var:
    """Minimal stand-in for a Tk variable."""
    def __init__(self, value=None):
        self._v = value

    def get(self):
        return self._v

    def set(self, v):
        self._v = v


class Logger:
    def __init__(self):
        self.lines = []

    def log(self, msg):
        self.lines.append(str(msg))


class _After:
    def after(self, _ms, fn=None, *a):
        return None


class MessageboxStub:
    def __init__(self):
        self.calls = []

    def __getattr__(self, name):
        def _f(*a, **k):
            self.calls.append((name, a, k))
            return None
        return _f


def cga_stubs(runlog=None):
    """Closure stubs for run_pseudolabeling (ion mining and negatives off, no Tk)."""
    runlog = [] if runlog is None else runlog
    return {
        "logger": Logger(),
        "messagebox": MessageboxStub(),
        "append_runlog": lambda files, entry: runlog.append(entry),
        "refresh_tree": lambda *a, **k: None,
        "subwin": _After(),
        "add_negatives_var": Var(False),
        "ion_suggest_enable_var": Var(False),
        "ion_suggest_dafloor_var": Var(0.0),
        "ion_suggest_minsupp_var": Var(0.0),
        "ion_suggest_ppm_var": Var("10"),
        "ion_suggest_topk_var": Var(0),
        "last_suggest_csv_var": Var(""),
        "min_hits_var": Var(3),
        "neg_ratio_var": Var(3.0),
        "ppm_tol_var": Var("10"),
    }
