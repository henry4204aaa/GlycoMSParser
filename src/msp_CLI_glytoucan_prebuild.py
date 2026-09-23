#!/usr/bin/env python3
"""
msp_CLI_glytoucan_prebuild.py — prebuild / update the GlyTouCan-WURCS reference database
from a composition source, without the GUI. (v1.20 S1 / T1; spec: handoff/WURCS_designspec20260921.md)

This is the ONLY network path in v1.20 today. Run-time fill inside GlycoMSP reads the
database and never calls the API.

Usage:
    python src/msp_CLI_glytoucan_prebuild.py --source path/to/sample_CGA_*.tsv
    python src/msp_CLI_glytoucan_prebuild.py --source path/to/MAS_sheet.xlsx --column Structure
    python src/msp_CLI_glytoucan_prebuild.py --source path/to/insilico_library.csv --no-network   # diff only
    python src/msp_CLI_glytoucan_prebuild.py --summary                                          # database status

Sources: CGA TSV (`composition`), trainable CSV / MAS sheet (`Structure`), in-silico library
CSV (auto-detected; use --column to override). Batches of <=30 with a pause between batches
(GlyCosmos courtesy rule). Fail soft: a network failure keeps what was fetched and exits 0
with status 'partial'/'offline' printed; nothing is stored for compositions that were not answered.

Locations default to the per-user data directory; override with --db / --log or the
GLYCOMSP_REFERENCE_DB / GLYCOMSP_GLYTOUCAN_LOG environment variables.

Exit codes: 0 = done, nothing to do, or network fail-soft (offline/partial: rerun to resume);
            2 = storage problem — database write failed (rows retried next run) or the JSONL log
                could not be written (database is complete; run --rebuild-log).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import msp_glytoucan_resolver as gtr  # noqa: E402


EXIT_CODES_HELP = """\
exit codes:
  0  done, nothing to do, or network fail-soft (status offline/partial: rerun to resume — check n_unsent)
  2  storage problem — status persist_failed (DB write failed; rows in n_unpersisted are retried next run)
     or log_failed (fetched rows are committed but the JSONL log is incomplete: fix the log location, run
     --rebuild-log, and STILL check n_unsent/n_unpersisted — log_failed does not mean every requested
     composition is in the database); or --rebuild-log refused (missing/foreign database, backup/replace failure)
recovery:
  --rebuild-log opens the database READ-ONLY (never creates or initialises one), writes the full projection to
  <log>.tmp, copies the current log to <log>.bak, then replaces the live log atomically.
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1].strip(),
                                 epilog=EXIT_CODES_HELP,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", help="CGA TSV / trainable CSV / MAS sheet / in-silico library")
    ap.add_argument("--column", default=None, help="composition column (auto-detect if omitted)")
    ap.add_argument("--sheet", default=None, help="Excel sheet name (default: MSlist if present, else first sheet)")
    ap.add_argument("--db", default=None, help=f"reference database path (default: {gtr.default_db_path()})")
    ap.add_argument("--log", default=None, help=f"transaction log path (default: {gtr.default_log_path()})")
    ap.add_argument("--no-network", action="store_true", help="diff against the database only; no API calls")
    ap.add_argument("--batch-size", type=int, default=30, help="compositions per request (max 30)")
    ap.add_argument("--pause", type=float, default=10.0, help="seconds between batches")
    ap.add_argument("--timeout", type=float, default=20.0, help="per-request timeout (s)")
    ap.add_argument("--summary", action="store_true", help="print database summary and exit")
    ap.add_argument("--rebuild-log", action="store_true",
                    help="regenerate the JSONL transaction log from the database (recovery after a log write failure); previous log kept as .bak")
    ap.add_argument("--json", action="store_true", help="print the result as JSON (machine-readable)")
    args = ap.parse_args(argv)

    if not args.summary and not args.rebuild_log and not args.source:
        ap.error("--source is required unless --summary or --rebuild-log is given")

    if args.rebuild_log:
        # Codex R3-F1: recovery must never create/initialise a database — open read-only and validate first.
        dbp = args.db or str(gtr.default_db_path())
        logp = args.log or gtr.default_log_path()
        try:
            db_ro = gtr.ReferenceDB(dbp, readonly=True)
            n = gtr.rebuild_log(db_ro, logp)
            db_ro.close()
        except (gtr.RebuildLogError, FileNotFoundError, OSError) as ex:
            msg = f"rebuild-log refused: {ex}"
            print(json.dumps({"status": "rebuild_failed", "error": msg, "log": str(logp)}, indent=2) if args.json
                  else f"[prebuild][ERROR] {msg}", file=sys.stdout)
            return 2
        print(json.dumps({"status": "ok", "rebuilt_lines": n, "log": str(logp)}, indent=2) if args.json
              else f"[prebuild] rebuilt transaction log: {n} lines -> {logp} (previous kept as .bak)")
        return 0

    db = gtr.ReferenceDB(args.db)
    try:
        if args.summary:
            s = db.summary()
            s["path"] = str(db.path)
            print(json.dumps(s, indent=2) if args.json else
                  "\n".join(f"{k}: {v}" for k, v in s.items()))
            return 0

        keys = gtr.compositions_from_file(args.source, args.column, sheet=args.sheet)
        say = (lambda s: None) if args.json else (lambda s: print(f"[prebuild] {s}", flush=True))
        say(f"source {os.path.basename(args.source)}: {len(keys)} distinct compositions")
        res = gtr.prebuild(keys, db, args.log, allow_network=not args.no_network, progress=say,
                           batch_size=args.batch_size, pause_s=args.pause, timeout_s=args.timeout)
        out = {
            "status": res.status, "n_source": res.n_source, "n_already": res.n_already,
            "n_fetched": res.n_fetched, "api_calls": res.api_calls, "n_unsent": len(res.unsent),
            "n_unpersisted": len(res.unpersisted), "n_unlogged": len(res.unlogged),
            "errors": res.errors, "db": res.db,
        }
        if args.json:
            print(json.dumps(out, indent=2))
        else:
            for e in res.errors:
                print(f"[prebuild][WARN] {e}")
            print(f"[prebuild] status={res.status} committed={res.n_fetched} already={res.n_already} "
                  f"unsent={len(res.unsent)} unpersisted={len(res.unpersisted)} unlogged={len(res.unlogged)} "
                  f"api_calls={res.api_calls}")
            print(f"[prebuild] database: {res.db.get('path')} ({res.db.get('n_entries')} entries, "
                  f"{res.db.get('version')})")
            if res.status == "log_failed":
                print(f"[prebuild][WARN] fetched rows were committed but the transaction log is incomplete "
                      f"({len(res.unlogged)} unlogged) — fix the log location and run --rebuild-log. "
                      f"Also check unsent={len(res.unsent)} unpersisted={len(res.unpersisted)}: "
                      f"rerun prebuild if any requested composition is still missing.")
        # Exit codes: 0 = done or network fail-soft (offline/partial are retried by rerunning);
        #             2 = STORAGE problem (persist_failed / log_failed) — needs attention (Codex R2-F1/F2).
        return 0 if res.storage_ok else 2
    finally:
        db.close()


if __name__ == "__main__":
    sys.exit(main())
