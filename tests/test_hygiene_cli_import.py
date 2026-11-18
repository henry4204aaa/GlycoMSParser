"""Hygiene H2: `msp_CLI_make_trainable_with_negatives --add-decoys` imports
`append_feature_space_decoys` from `mspvalidator_merger`, not from the GUI.

Pre-fix the branch did `from mspfileloaderv14 import append_feature_space_decoys`,
which does not exist there (ImportError) and pulled in the whole Tkinter GUI.

Runs the CLI `main()` in a fresh interpreter so other test modules that import
v14 cannot pollute `sys.modules`. The merge/export steps are stubbed; the
`--add-decoys` branch runs for real. See handoff/hygiene_LF1_kickoff_20260926.md (H2).
"""
import os
import subprocess
import sys
import textwrap

SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))


def test_add_decoys_branch_does_not_import_gui(tmp_path):
    out_csv = tmp_path / "trainable.csv"
    script = textwrap.dedent(f"""
        import sys
        sys.path.insert(0, {SRC!r})
        import pandas as pd
        import mspvalidator_merger as mv
        import msp_CLI_make_trainable_with_negatives as cli

        mv.directassign_files = lambda *a, **k: (pd.DataFrame({{"MS2scan_no": [1]}}), [], None)
        def _fake_create(ion_index, pre_df, ion_df, out):
            pd.DataFrame({{"Structure": ["H5N4", "H3N4"],
                          "f1": [1.0, 2.0], "f2": [0.5, 0.0]}}).to_csv(out, index=False)
        mv.createnormailzedionlistcsv = _fake_create

        sys.argv = ["cli", "--excel", "x.xlsx", "--raw", "y.tsv",
                    "--out", {str(out_csv)!r}, "--add-decoys"]
        cli.main()
        assert "mspfileloaderv14" not in sys.modules, "CLI imported the GUI module"
        print("ROWS", len(pd.read_csv({str(out_csv)!r})))
    """)
    proc = subprocess.run([sys.executable, "-c", script],
                          capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr
    assert "ROWS 3" in proc.stdout  # 2 positives + int(0.5 * 2) = 1 decoy
