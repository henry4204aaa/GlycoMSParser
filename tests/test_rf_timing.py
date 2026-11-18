"""RF-TIMING: tools/validation/rf_timing.py on a tiny synthetic trainable (no fixtures)."""
import importlib.util
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "tools" / "validation" / "rf_timing.py"

EXPECTED_FIELDS = (
    "alias:", "rows_total:", "rows_after_eligibility_filter:", "features:", "classes_kept:",
    "train_rows:", "test_rows:", "python:", "scikit-learn:", "numpy:", "cpu_model:",
    "processor_string:", "logical_cores:", "fit_seconds_per_repeat:",
    "full_step_seconds_per_repeat", "fit_seconds_mean:", "fit_seconds_max:",
    "full_step_seconds_mean:", "full_step_seconds_max:", "multiclass_accuracy:",
)


def _load_script():
    spec = importlib.util.spec_from_file_location("rf_timing_under_test", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _synthetic_trainable(path):
    rng = np.random.default_rng(0)
    labels = ["Non-glycan", "H3N2", "H5N2", "H5N4", "F1H5N4", "H5N4S2"]  # 6 classes x 5 rows
    rows = []
    for k, lab in enumerate(labels):
        for _ in range(5):
            r = {"Structure": lab, "MS2scan_no": len(rows) + 1}
            for j in range(10):
                r[f"ion_{204 + j * 40}"] = float(rng.random() * 3 + k)
            r["protonatedmass"] = 0.0 if lab == "Non-glycan" else 800.0 + 100 * k
            rows.append(r)
    pd.DataFrame(rows).to_csv(path, index=False)


def test_rf_timing_report(tmp_path, capsys):
    mod = _load_script()
    csv = tmp_path / "synthetic_trainable.csv"
    _synthetic_trainable(csv)

    t0 = time.perf_counter()
    rc = mod.main(["--trainable", str(csv), "--alias", "synth", "--repeats", "1",
                   "--trees", "10", "--out", str(tmp_path)])
    elapsed = time.perf_counter() - t0
    capsys.readouterr()

    assert rc == 0
    assert elapsed < 5.0
    report = tmp_path / "rf_timing_synth.txt"
    assert report.is_file()
    text = report.read_text(encoding="utf-8")
    for field in EXPECTED_FIELDS:
        assert field in text, field
    kept = int(next(l for l in text.splitlines() if l.startswith("classes_kept:")).split(":")[1])
    assert kept >= 5
    assert "protonatedmass_feature_included: yes" in text
    # the input path is never written into the report
    assert str(tmp_path) not in text and csv.name not in text
