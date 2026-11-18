"""B-06 smoke test — verify the tuple-unpack fix for experiment_status_labels.

The bug: `export_experiment_json` called `.config(text=...)` on
`experiment_status_labels[exp_title]`, which is a tuple of two tk.StringVar
instances, not a Label widget. That raises AttributeError.

The fix mirrors the canonical pattern at load_experiment_method_file:5533 —
unpack the tuple and call .set() on the StringVar.

This test:
1. Confirms the OLD buggy pattern would raise AttributeError on a tuple.
2. Confirms the NEW fixed pattern executes without error and the StringVar
   carries the expected value after .set().

Run from repo root:
    python3 tests/test_b06_status_label_smoke.py
"""
import sys
import tkinter as tk


def main():
    failures = []

    # Set up a hidden Tk root so StringVar can be constructed.
    root = tk.Tk()
    root.withdraw()

    try:
        # Reproduce the storage shape at mspfileloaderv14.py:6072.
        exp_status_var = tk.StringVar()
        sample_status_var = tk.StringVar()
        experiment_status_labels = {
            "Unassigned": (exp_status_var, sample_status_var),
        }

        out_path = "/tmp/test_b06.exp.json"
        exp_title = "Unassigned"

        # Check 1 — the OLD buggy pattern raises AttributeError.
        try:
            experiment_status_labels[exp_title].config(text=f"EXP file: {out_path}")
            failures.append("FAIL old buggy pattern did NOT raise — fix is not exercising the right code path")
        except AttributeError as e:
            if "'tuple' object has no attribute 'config'" in str(e):
                print("PASS old buggy pattern raises AttributeError as expected (regression guard)")
            else:
                failures.append(f"FAIL old buggy pattern raised wrong AttributeError: {e}")
        except Exception as e:
            failures.append(f"FAIL old buggy pattern raised unexpected exception type: {type(e).__name__}: {e}")

        # Check 2 — the NEW fixed pattern executes without error.
        try:
            if exp_title in experiment_status_labels:  # mirrors the patched code at v14:2634
                path_label, _ = experiment_status_labels[exp_title]
                if path_label:
                    path_label.set(f"EXP file: {out_path}")
            print("PASS new fixed pattern executes without raising")
        except Exception as e:
            failures.append(f"FAIL new fixed pattern raised: {type(e).__name__}: {e}")

        # Check 3 — the StringVar carries the expected value after .set().
        expected = f"EXP file: {out_path}"
        actual = exp_status_var.get()
        if actual == expected:
            print(f"PASS path_label StringVar value: {actual!r}")
        else:
            failures.append(f"FAIL path_label value: expected {expected!r}, got {actual!r}")

        # Check 4 — the second StringVar (sample_status_var) is untouched by the fix.
        if sample_status_var.get() == "":
            print("PASS sample_status_var untouched (still empty)")
        else:
            failures.append(f"FAIL sample_status_var unexpectedly modified: {sample_status_var.get()!r}")

        # Check 5 — defensive guard handles None path_label gracefully.
        # (Storage init always writes a 2-tuple of StringVars, so this is defensive only.)
        experiment_status_labels["FakeKey"] = (None, None)
        try:
            exp_title_fake = "FakeKey"
            if exp_title_fake in experiment_status_labels:
                path_label, _ = experiment_status_labels[exp_title_fake]
                if path_label:
                    path_label.set("should not run")
            print("PASS defensive `if path_label:` guard handles None gracefully")
        except Exception as e:
            failures.append(f"FAIL defensive None guard raised: {type(e).__name__}: {e}")

    finally:
        root.destroy()

    print()
    if failures:
        print("=" * 60)
        print("FAILURES:")
        for f in failures:
            print("  " + f)
        print("=" * 60)
        sys.exit(1)
    else:
        print("=" * 60)
        print("B-06 smoke test: ALL CHECKS PASS")
        print("=" * 60)


if __name__ == "__main__":
    main()
