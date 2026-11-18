#!/usr/bin/env python3
"""DB1 test stub standing in for msp_CLI_glytoucan_prebuild.py (tests/test_db1_manage_references.py, CP-C).

Prints three progress lines, then (unless suppressed) the CLI's own status line and database line, then exits.
No network, no database. Options:
  --rc N            exit code (default 0)
  --status S        value after "status=" (default ok)
  --no-status       omit the status line (simulates a crash before the summary)
  --hang            after the first line, sleep until terminated (for Cancel tests)
  --ignore-term     ignore SIGTERM (POSIX) so only kill() ends the process
  --json-result J   print J (a JSON object) instead of progress lines (rebuild-log mode)
  --unicode         include a non-ASCII line (encoding round trip)
"""
import argparse
import signal
import sys
import time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rc", type=int, default=0)
    ap.add_argument("--status", default="ok")
    ap.add_argument("--no-status", action="store_true")
    ap.add_argument("--hang", action="store_true")
    ap.add_argument("--ignore-term", action="store_true")
    ap.add_argument("--json-result", default=None)
    ap.add_argument("--unicode", action="store_true")
    a = ap.parse_args()
    if a.ignore_term and hasattr(signal, "SIGTERM") and sys.platform != "win32":
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
    if a.json_result is not None:
        print(a.json_result, flush=True)
        return a.rc
    print("[prebuild] source stub.tsv: 3 distinct compositions", flush=True)
    if a.hang:
        print("[prebuild] ready", flush=True)
        while True:
            time.sleep(0.05)
    print("[prebuild] 3 compositions in source; 0 already in database; 3 to fetch", flush=True)
    print("[prebuild] batch 1/1 (3 compositions)" + (" — H5N4 数据" if a.unicode else ""), flush=True)
    if not a.no_status:
        print(f"[prebuild] status={a.status} committed=3 already=0 unsent=2 unpersisted=0 unlogged=0 api_calls=1")
        print("[prebuild] database: /tmp/stub.db (3 entries, v0.01.20260921)")
    return a.rc


if __name__ == "__main__":
    sys.exit(main())
