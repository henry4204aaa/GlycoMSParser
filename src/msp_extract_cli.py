"""GlycoMSP extraction CLI (READER-DIAG CP-C) - thin wrapper over the reader-side functions.

version 0.42, 20261023 (v0.2 READER-DIAG CP-C.5: `convert`; v0.3 CP-C.6: --charge-comb on mzML too; v0.4 CP-D: diagnostics
files and outcome-based exit codes; v0.41 CP-D.1: messages redacted in one pass; v0.42 R2: own work folder per run,
CSVs already written are kept after a failure, no recursive delete, whole record redacted). Commands in this build:

  inspect   read-only probe of selected native scans in a Thermo RAW file (Windows + MSFileReader)
  convert   normal extraction to the MS2/MS3 CSVs, with the optional charge sources

  python src/msp_extract_cli.py inspect INPUT.raw --alias NAME --scans 6781 6783 6785 --out DIR
         [--window 3.0] [--full-points]
  python src/msp_extract_cli.py convert (--raw INPUT.raw | --mzml INPUT.mzML) --out DIR [--alias NAME]
         [--charge-comb] [--dataset-charge N]

convert writes ms2_<name>.csv, ms3_<name>.csv, <name>.convert.json, <name>.scan_status.tsv and
<name>.reader_diag.json/.txt in DIR, where <name> is the alias, or the input file stem when no alias is given. With an alias the input name appears in no output. The CSVs
are the converter's usual 9-/7-column files. --charge-comb and --dataset-charge fill MS2 charges the file does
not report; a reported charge is never overwritten. The comb needs a profile parent MS1 (on mzML: an MS1
converted without peak picking).

Exit codes, convert (the meaning is also written in <name>.reader_diag.txt):
  0 complete: converted, nothing to check
  3 complete with warnings: converted; MS2 charges unknown, no usable precursor mass, or profile MS2
  4 partial: converted, but some spectra could not be decoded and were not exported
  5 failed: the input could not be opened, the conversion stopped, no MS2 spectrum with peaks was exported, or the
    diagnostics files could not be written. A CSV written before the failure is kept (the record says failed).
  2 usage error (including an invalid --dataset-charge)
Exit codes, inspect: 0 done (individual backend requests may still have failed; see the report), 2 usage error,
5 the input could not be opened or the backend is unavailable.
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import msp_raw_probe  # noqa: E402
from msp_spectrum_record import (EXIT_OK, EXIT_USAGE, EXIT_WARNINGS, EXIT_PARTIAL, EXIT_FAILED,  # noqa: E402,F401
                                 OUTCOME_FAILED, Redactor, redact_tree, update_diagnostics_outputs)
version = "0.42"
last_update = 20261023
_ALIAS_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")


def _alias(value):
    if not _ALIAS_RE.match(value) or value in (".", ".."):
        raise argparse.ArgumentTypeError("alias must be 1-64 characters: letters, digits, '_', '-', '.'")
    return value


def build_parser():
    parser = argparse.ArgumentParser(prog="msp_extract_cli", description=__doc__.split("\n\n")[0],
                                     epilog="Exit codes: see the module docstring or <name>.reader_diag.txt.")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("inspect", help="read-only probe of selected RAW scans")
    p.add_argument("input", help="Thermo .raw file")
    p.add_argument("--alias", required=True, type=_alias,
                   help="name used for every output file and inside the report instead of the input file name")
    p.add_argument("--scans", required=True, type=int, nargs="+", help="native scan numbers to probe")
    p.add_argument("--out", required=True, help="output folder (created if missing)")
    p.add_argument("--window", type=float, default=3.0,
                   help="m/z half-width kept around the precursor in the parent MS1 (default 3.0)")
    p.add_argument("--full-points", action="store_true",
                   help="also write every returned point of MS1 scans requested directly")
    c = sub.add_parser("convert", help="convert a RAW or mzML file to the MS2/MS3 CSVs")
    src = c.add_mutually_exclusive_group(required=True)
    src.add_argument("--raw", help="Thermo .raw file (Windows + MSFileReader)")
    src.add_argument("--mzml", help=".mzML file (pymzml)")
    c.add_argument("--out", required=True, help="output folder (created if missing)")
    c.add_argument("--alias", type=_alias, help="name used for every output file instead of the input file name")
    c.add_argument("--charge-comb", action="store_true",
                   help="infer missing MS2 charges from the parent MS1 isotope pattern (profile MS1 only, not FTMS)")
    c.add_argument("--dataset-charge", type=int,
                   help="one charge for every MS2 scan whose charge is still unknown (non-zero integer)")
    return parser


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _describe(path):
    with open(path, "rb") as fh:
        header = fh.readline().decode("utf-8", "replace").rstrip("\r\n")
        rows = sum(1 for _ in fh)
    return {"file": os.path.basename(path), "rows": rows, "header": header.split("\t"),
            "bytes": os.path.getsize(path), "sha256": _sha256(path)}


def _promote(src, dst):
    """Move one finished file into place as v14's promote_temp_to_final does: rename on the same drive, otherwise copy
    and then remove the source file only (R2-F1)."""
    try:
        os.replace(src, dst)
    except OSError:
        shutil.copy2(src, dst)
        os.remove(src)
    return dst


def _move_diagnostics(work, out, name):
    """Move the reader's status and diagnostics files from the work folder to `out`; {key: file name}."""
    moved = {}
    for key, fname in (("json", f"{name}.reader_diag.json"), ("txt", f"{name}.reader_diag.txt"),
                       ("status_tsv", f"{name}.scan_status.tsv")):
        src = os.path.join(work, fname)
        if os.path.isfile(src):
            _promote(src, os.path.join(out, fname))
            moved[key] = fname
    return moved


def _clean_work(work):
    """Remove this run's own work folder when it is empty, after deleting the reader's half-written `*.tmp` files only.
    Never recursive, never anything outside it (R2-F1); a CSV that could not be moved stays (Henry D1). Returns the
    names left behind (empty when the folder is gone)."""
    left = []
    for entry in os.scandir(work):
        try:
            if entry.name.endswith(".tmp") and entry.is_file(follow_symlinks=False):
                os.remove(entry.path)
                continue
        except OSError:
            pass
        left.append(entry.name)
    try:
        os.rmdir(work)
    except OSError:
        left = left or [os.path.basename(work)]
    return left


def _finalize(work, out, name, res, record):
    """Promote every CSV the reader wrote (also after a failure, Henry D1), move the diagnostics, set the outcome.

    Files are moved first and described afterwards, so a failure while describing still leaves every file in place."""
    promoted = {}
    for level in ("ms2", "ms3"):
        src = res.get(f"{level}tmp") if isinstance(res, dict) else None
        src = src or os.path.join(work, f"{level}tmp_{name}.csv")
        if os.path.isfile(src):
            promoted[level] = _promote(src, os.path.join(out, f"{level}_{name}.csv"))
    record["diagnostics"] = _move_diagnostics(work, out, name)
    for level in ("ms2", "ms3"):
        record[level] = _describe(promoted[level]) if level in promoted else None
    if isinstance(res, dict):
        record["ok"] = res["outcome"] != OUTCOME_FAILED
        record["outcome"], record["exit_code"], record["warnings"] = res["outcome"], res["exit_code"], res["warnings"]
        if res.get("diagnostics_error") or "json" not in record["diagnostics"]:   # Henry D2: the CLI promises them
            record["diagnostics_error"] = res.get("diagnostics_error") or "diagnostics files missing"
            record["ok"], record["outcome"], record["exit_code"] = False, OUTCOME_FAILED, EXIT_FAILED
    else:                       # the reader returned False or raised: the diagnostics say failed
        record["ok"] = False
        record["outcome"], record["exit_code"] = OUTCOME_FAILED, EXIT_FAILED
    if "json" in record["diagnostics"]:
        update_diagnostics_outputs(os.path.join(out, record["diagnostics"]["json"]),
                                   ms2_csv=(record["ms2"] or {}).get("file"), ms3_csv=(record["ms3"] or {}).get("file"))


def _exception_fields(exc, redact, stage=None):
    fields = {"exception": type(exc).__name__, "message": redact(str(exc)),
              "traceback": [f"{os.path.basename(f.filename)}:{f.lineno}"
                            for f in traceback.extract_tb(exc.__traceback__)]}
    if stage:
        fields["stage"] = stage
    return fields


def run_convert(args):
    """Run one conversion; return (exit code, record). Outputs and the record carry the alias, never the input name."""
    route, inp = ("raw", args.raw) if args.raw else ("mzml", args.mzml)
    if not os.path.isfile(inp):
        print("input file not found", file=sys.stderr)
        return EXIT_USAGE, None
    if args.dataset_charge == 0:
        print("--dataset-charge must be a non-zero integer", file=sys.stderr)
        return EXIT_USAGE, None
    inp = os.path.abspath(inp)
    stem = os.path.splitext(os.path.basename(inp))[0]
    name = args.alias or stem
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    work = tempfile.mkdtemp(prefix="_work_", dir=out)   # v0.42 R2-F1: a new folder owned by this run, never reused

    redact = Redactor(inp, name)   # v0.41: one pass (was sequential replaces: "bad.raw" -> "bad_raw_raw")

    record = {"json_type": "glycomsp_convert_record", "schema_version": 1, "name": name, "route": route,
              "cli_version": version,
              "options": {"charge_comb": bool(args.charge_comb), "dataset_charge": args.dataset_charge},
              "started": time.strftime("%Y-%m-%d %H:%M:%S")}
    t0 = time.time()
    res, interrupted = None, None
    try:
        if route == "raw":
            import msprawextractor as mod
            res = mod.convert_raw_to_csv(inp, outdir=work, charge_comb=args.charge_comb,
                                         dataset_charge=args.dataset_charge, stem=name, entry="cli-convert")
        else:
            import mzmlreader as mod
            res = mod.extract_mzML(inp, outdir=work, dataset_charge=args.dataset_charge,
                                   charge_comb=args.charge_comb, stem=name, entry="cli-convert")
        record["module"] = {"name": mod.__name__, "version": str(getattr(mod, "version", None)),
                            "sha256": _sha256(mod.__file__)}
        if isinstance(res, dict):
            record["summary"] = res.get("summary")
        else:
            record["returned"] = repr(res)
    except BaseException as exc:   # v0.42: Ctrl-C too - keep what was written, then re-raise it below
        record.update(_exception_fields(exc, redact))
        if not isinstance(exc, Exception):
            interrupted = exc
    record["seconds"] = round(time.time() - t0, 1)
    try:
        _finalize(work, out, name, res, record)
    except Exception as exc:       # v0.42 R2-F2: a failure while moving files is a failed conversion, in both records
        failure = _exception_fields(exc, redact, stage="finalize")
        record["finalize_error"] = failure
        record["ok"], record["outcome"], record["exit_code"] = False, OUTCOME_FAILED, EXIT_FAILED
        try:
            if not record.get("diagnostics"):   # the failure came before the diagnostics were moved
                record["diagnostics"] = _move_diagnostics(work, out, name)
        except Exception:
            pass
        try:
            if (record.get("diagnostics") or {}).get("json"):
                update_diagnostics_outputs(os.path.join(out, record["diagnostics"]["json"]), failure=failure)
        except Exception:
            pass
    record.setdefault("diagnostics", {})
    left = _clean_work(work)
    if left:
        record["work_folder_left"] = {"folder": os.path.basename(work), "entries": left}
    with open(os.path.join(out, f"{name}.convert.json"), "w", encoding="utf-8") as fh:
        json.dump(redact_tree(record, redact), fh, indent=2, default=str)
    if interrupted is not None:
        raise interrupted
    return record["exit_code"], record


def main(argv=None):
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        return EXIT_USAGE if exc.code else EXIT_OK
    if args.command == "inspect":
        if not os.path.isfile(args.input):
            print("input file not found", file=sys.stderr)
            return EXIT_USAGE
        report = msp_raw_probe.inspect_raw_scans(args.input, args.scans, args.alias, args.out,
                                                 window=args.window, full_points=args.full_points)
        print(msp_raw_probe.summary_text(report))
        print(f"report: {args.alias}.inspect.json (in the --out folder)")
        return EXIT_OK if report["outcome"] == "complete" else EXIT_FAILED
    if args.command == "convert":
        code, record = run_convert(args)
        if record is not None:
            charge = (record.get("summary") or {}).get("charge")
            print(json.dumps({k: record.get(k) for k in ("name", "route", "ok", "outcome", "exit_code", "seconds",
                                                         "exception", "message")}))
            if charge and "sources" in charge:
                print(f"charge sources: {charge['sources']}  comb: {charge['comb_status']}")
            for warning in record.get("warnings") or []:
                print(f"warning: {warning}")
            # v0.42 R2-F3: file names only; the --out folder can sit inside the input's folder
            if record["diagnostics"].get("txt"):
                print(f"diagnostics: {record['diagnostics']['txt']} (in the --out folder)")
            print(f"record: {record['name']}.convert.json (in the --out folder)")
        return code
    return EXIT_USAGE


if __name__ == "__main__":
    sys.exit(main())
