"""B-04 smoke test — exercise _excepthook in isolation.

Spawns a tiny subprocess that defines the same diagnostics-prelude code as
src/mspfileloaderv14.py lines 1-127, then deliberately raises an uncaught
exception so _excepthook fires. Asserts:

1. Exit code is 1 (Fix 2).
2. Stdout contains the [StartupError] banner, the traceback text, and the
   environment-check header — all in the same single stream (Fix 1).
3. A log file was written under ~/.glycomsp/ (macOS / Linux path) with
   matching content (Fix 3).
4. Microsecond filename suffix is present (Codex Finding 2 fix).

Run from repo root:
    python3 tests/test_b04_excepthook_smoke.py
"""
import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path


def main():
    # Synthetic prelude — mirrors the production diagnostics block lines 1-127.
    # Stripped of unrelated module imports so the test runs in any clean Python env.
    snippet = textwrap.dedent("""
        import os, sys, importlib, traceback, io, platform
        from pathlib import Path
        from datetime import datetime

        def _pause_if_no_tty():
            try:
                if os.name == "nt" and not sys.stdin.isatty():
                    input("\\nPress Enter to close...")
            except Exception:
                pass

        def _check_mod(name, alt=None):
            modname = alt or name
            spec = importlib.util.find_spec(modname)
            if not spec:
                return False, "not found"
            try:
                m = __import__(modname)
                v = getattr(m, "__version__", "")
                return True, (v or "present")
            except Exception as e:
                return False, f"import failed: {e.__class__.__name__}: {e}"

        def _quick_diag(note=None, file=None):
            f = file if file is not None else sys.stdout
            print("=== GlycoMSP quick environment check ===", file=f)
            print(f"Python: {sys.version.split()[0]}  ({sys.executable})", file=f)
            print(f"Platform: {platform.platform()}", file=f)
            if note:
                print(note, file=f)

        def _excepthook(exctype, value, tb):
            log_path = None
            try:
                if os.name == "nt":
                    log_dir = Path(os.environ.get("LOCALAPPDATA", os.path.expanduser("~"))) / "GlycoMSP"
                else:
                    log_dir = Path(os.path.expanduser("~")) / ".glycomsp"
                log_dir.mkdir(parents=True, exist_ok=True)
                log_path = log_dir / f"startup_error_{datetime.now():%Y%m%d_%H%M%S_%f}.log"
            except Exception:
                pass

            buf = io.StringIO()
            banner = "\\n[StartupError] Unhandled exception during launch/import:\\n"
            print(banner, file=sys.stdout)
            buf.write(banner + "\\n")
            traceback.print_exception(exctype, value, tb, file=sys.stdout)
            traceback.print_exception(exctype, value, tb, file=buf)
            print(file=sys.stdout)
            buf.write("\\n")
            _quick_diag(file=sys.stdout)
            _quick_diag(file=buf)

            if log_path is not None:
                try:
                    log_path.write_text(buf.getvalue(), encoding="utf-8")
                    print(f"\\n[StartupError] Crash details saved to: {log_path}", file=sys.stdout)
                except Exception as e:
                    print(f"\\n[StartupError] Failed to save crash log: {e}", file=sys.stdout)

            _pause_if_no_tty()
            sys.exit(1)

        sys.excepthook = _excepthook
        raise RuntimeError("B-04 synthetic smoke test trigger")
    """).strip()

    # Snapshot existing log files so we can detect the new one.
    log_dir = Path.home() / ".glycomsp"
    before = set(log_dir.glob("startup_error_*.log")) if log_dir.exists() else set()

    proc = subprocess.run(
        [sys.executable, "-c", snippet],
        capture_output=True, text=True, timeout=15,
    )

    after = set(log_dir.glob("startup_error_*.log"))
    new_logs = after - before

    failures = []

    # Check 1 — exit code 1
    if proc.returncode != 1:
        failures.append(f"FAIL exit code: expected 1, got {proc.returncode}")
    else:
        print("PASS exit code: 1")

    # Check 2a — banner on stdout
    if "[StartupError] Unhandled exception" not in proc.stdout:
        failures.append("FAIL banner missing from stdout")
    else:
        print("PASS banner on stdout")

    # Check 2b — traceback on stdout (single-stream)
    if "RuntimeError: B-04 synthetic smoke test trigger" not in proc.stdout:
        failures.append("FAIL traceback text missing from stdout (Fix 1 — single stream)")
    else:
        print("PASS traceback on stdout (single stream)")

    # Check 2c — env check on stdout
    if "GlycoMSP quick environment check" not in proc.stdout:
        failures.append("FAIL env check missing from stdout")
    else:
        print("PASS env check on stdout")

    # Check 2d — order: banner before traceback before env-check
    banner_pos = proc.stdout.find("[StartupError] Unhandled exception")
    tb_pos = proc.stdout.find("RuntimeError: B-04 synthetic")
    env_pos = proc.stdout.find("GlycoMSP quick environment check")
    if not (banner_pos < tb_pos < env_pos):
        failures.append(f"FAIL stream order: banner@{banner_pos} tb@{tb_pos} env@{env_pos}")
    else:
        print("PASS stream order (banner -> traceback -> env)")

    # Check 3 — log file written
    if not new_logs:
        failures.append(f"FAIL no new log file written under {log_dir}")
    else:
        log_file = next(iter(new_logs))
        print(f"PASS log file written: {log_file}")
        # Check 4 — microsecond suffix
        if not re.match(r"^startup_error_\d{8}_\d{6}_\d{6}\.log$", log_file.name):
            failures.append(f"FAIL log filename pattern: {log_file.name}")
        else:
            print("PASS microsecond filename suffix")
        # Check 3b — log content matches stdout
        content = log_file.read_text(encoding="utf-8")
        if "[StartupError] Unhandled exception" not in content:
            failures.append("FAIL banner missing from log file")
        elif "RuntimeError: B-04 synthetic smoke test trigger" not in content:
            failures.append("FAIL traceback missing from log file")
        elif "GlycoMSP quick environment check" not in content:
            failures.append("FAIL env check missing from log file")
        else:
            print("PASS log file content matches stdout dump")
        # Check 3c — crash-path line on stdout
        if f"Crash details saved to:" not in proc.stdout:
            failures.append("FAIL 'Crash details saved to:' line missing from stdout")
        else:
            print("PASS 'Crash details saved to:' message on stdout")

    print()
    if failures:
        print("=" * 60)
        print("FAILURES:")
        for f in failures:
            print("  " + f)
        print("=" * 60)
        print("STDOUT was:")
        print(proc.stdout)
        print("STDERR was:")
        print(proc.stderr)
        sys.exit(1)
    else:
        print("=" * 60)
        print("B-04 smoke test: ALL CHECKS PASS")
        print("=" * 60)


if __name__ == "__main__":
    main()
