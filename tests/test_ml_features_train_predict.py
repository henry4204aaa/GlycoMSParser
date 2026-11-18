"""Stages 5-6 of the workflow, end to end on synthetic spectra (TESTS-PUBLIC CP-B, 2026-11-18).

features (ml_ng_utils) -> Random Forest -> .skops file -> command-line prediction (msp_CLI_predict_skops)
-> tau confidence rule -> prediction report (prediction_report).  No fixtures, no network; every input is built
here from a fixed seed, and the whole file runs in a few seconds.

What is and is not product code here
* Product code under test: `build_features_from_peaks_log10_plus1`, `predict_with_threshold`, `tau_sweep_summary`,
  `msp_CLI_predict_skops.main`, `summarize_predictions`, `write_prediction_report`, and the margin-aware demotion
  expression that the GUI uses (read from the v14 source text, so a change there is noticed).
* The GUI's training step lives in closures inside `train_model()` and cannot be imported without starting Tk.
  The tests therefore build the estimator with the product's settings (RandomForestClassifier, labels encoded with
  LabelEncoder, fitted on a DataFrame) and check the artefact contract the later stages depend on.  Timing of the
  GUI's train step is covered by tools/validation/rf_timing.py.
* Two tests are marked CHARACTERISATION: they record what the command line does today so that a change is
  deliberate.  They are not a statement that the behaviour is desired (see their docstrings).
"""
import ast
import json
import math
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("skops")
import joblib  # noqa: E402
from sklearn.ensemble import RandomForestClassifier  # noqa: E402
from sklearn.preprocessing import LabelEncoder  # noqa: E402
import skops  # noqa: E402
from skops.io import dump as skops_dump, get_untrusted_types, load as skops_load  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import msp_CLI_predict_skops as cli  # noqa: E402
import prediction_report as PR  # noqa: E402
from ml_ng_utils import build_features_from_peaks_log10_plus1 as build_features  # noqa: E402
from ml_ng_utils import predict_with_threshold  # noqa: E402
from ml_ng_utils_extras import tau_sweep_summary  # noqa: E402

V14 = ROOT / "src" / "mspfileloaderv14.py"
SKOPS_REJECTS_TRUSTED_TRUE = tuple(int(x) for x in skops.__version__.split(".")[:2]) >= (0, 10)

ION_MASSES = [204.0867, 366.1395, 528.1924, 657.2349]
NG = "Non-glycan"
CLASS_A, CLASS_B = "H5N2", "F1H5N4S2"
MASS = {CLASS_A: 1235.44, CLASS_B: 2369.85, NG: 0.0}   # the GUI sets the mass feature of Non-glycan rows to 0.0


def _load(path):
    """Load a .skops file the way the GUI does (it passes the file's own untrusted types, v14 `_load_model`)."""
    return skops_load(path, trusted=get_untrusted_types(file=path))


# ---------------------------------------------------------------- synthetic spectra

def _spectra(n_per_class, seed):
    """Class A shows ions 0 and 1, class B ions 2 and 3, Non-glycan only noise away from the ions."""
    rng = np.random.default_rng(seed)
    rows = []
    for label, hit in ((CLASS_A, (0, 1)), (CLASS_B, (2, 3)), (NG, ())):
        for _ in range(n_per_class):
            peaks = [float(x) for x in rng.uniform(100.0, 150.0, 4)]        # noise below every ion
            ints = [float(x) for x in rng.uniform(10.0, 50.0, 4)]           # plain floats: numpy 2 reprs are not literals
            for j in hit:
                peaks.append(float(ION_MASSES[j] * (1 + rng.uniform(-5e-6, 5e-6))))  # within 5 ppm of the ion
                ints.append(float(rng.uniform(1e3, 1e5)))
            rows.append({"peaklist": str([round(p, 5) for p in peaks]),
                         "peakintensity": str([round(i, 3) for i in ints]),
                         "Structure": label})
    return pd.DataFrame(rows)


def _trainable(n_per_class, seed, with_mass):
    spectra = _spectra(n_per_class, seed)
    feats = build_features(spectra, ION_MASSES, ppm=10.0)
    if with_mass:
        feats["protonatedmass"] = spectra["Structure"].map(MASS).to_numpy()
    return feats, spectra["Structure"]


def _fit(X, labels):
    le = LabelEncoder().fit(labels)
    model = RandomForestClassifier(n_estimators=60, random_state=42, class_weight="balanced", n_jobs=1)
    model.fit(X, le.transform(labels))            # the GUI fits on encoded integers
    return model, le


@pytest.fixture(scope="module")
def chain(tmp_path_factory):
    """Train on 10 spectra per class, save as .skops (+ label-encoder sidecar), hold out 5 per class from another seed."""
    out = tmp_path_factory.mktemp("chain")
    X, y = _trainable(10, 1, with_mass=False)
    model, le = _fit(X, y)
    skops_path = out / "toy_rf_model.skops"
    skops_dump(model, skops_path)
    le_path = out / "toy_labelencoder.joblib"
    joblib.dump(le, le_path)
    Xh, yh = _trainable(5, 2, with_mass=False)
    return {"dir": out, "model": model, "le": le, "skops": skops_path, "le_path": le_path, "X": X, "Xh": Xh, "yh": yh}


# ---------------------------------------------------------------- features: log10(I) + 1, baseline 1.0

def _one(peaks, ints):
    return pd.DataFrame({"peaklist": [str(peaks)], "peakintensity": [str(ints)]})


def test_features_hit_is_log10_intensity_plus_one_and_miss_is_one():
    f = build_features(_one([204.0867, 400.0], [1000.0, 50.0]), ION_MASSES, ppm=10.0)
    assert list(f.columns) == [str(m) for m in ION_MASSES]
    assert f.iloc[0, 0] == pytest.approx(math.log10(1000.0) + 1.0)      # 4.0, NOT log10(1000 + 1)
    assert f.iloc[0, 0] == pytest.approx(4.0)
    assert list(f.iloc[0, 1:]) == [1.0, 1.0, 1.0]


def test_features_use_the_largest_intensity_inside_the_window():
    f = build_features(_one([204.0867, 204.0868, 204.0869], [10.0, 1000.0, 100.0]), [204.0867], ppm=10.0)
    assert f.iloc[0, 0] == pytest.approx(4.0)


def test_features_ppm_window_is_relative_to_the_target():
    target = 500.0
    inside = target * (1 + 9e-6)       # 9 ppm
    outside = target * (1 + 11e-6)     # 11 ppm
    assert build_features(_one([inside], [100.0]), [target], ppm=10.0).iloc[0, 0] == pytest.approx(3.0)
    assert build_features(_one([outside], [100.0]), [target], ppm=10.0).iloc[0, 0] == 1.0
    # the window is symmetric
    assert build_features(_one([target * (1 - 9e-6)], [100.0]), [target], ppm=10.0).iloc[0, 0] == pytest.approx(3.0)


def test_features_default_window_is_20_ppm():
    target = 500.0
    peak = target * (1 + 15e-6)
    assert build_features(_one([peak], [100.0]), [target]).iloc[0, 0] == pytest.approx(3.0)        # default ppm=20
    assert build_features(_one([peak], [100.0]), [target], ppm=10.0).iloc[0, 0] == 1.0


def test_features_unparseable_row_stays_at_baseline_and_index_is_kept():
    df = pd.DataFrame({"peaklist": ["[204.0867]", "not a list"], "peakintensity": ["[1000.0]", "[1.0]"]},
                      index=[7, 9])
    f = build_features(df, [204.0867], ppm=10.0)
    assert list(f.index) == [7, 9]
    assert f.loc[7, "204.0867"] == pytest.approx(4.0)
    assert f.loc[9, "204.0867"] == 1.0


def test_features_need_peaklist_and_peakintensity():
    with pytest.raises(ValueError, match="peaklist"):
        build_features(pd.DataFrame({"peaklist": ["[1.0]"]}), [1.0])


# ---------------------------------------------------------------- training and the .skops artefact

def test_forest_separates_the_synthetic_classes_on_held_out_spectra(chain):
    model, le = chain["model"], chain["le"]
    predicted = le.inverse_transform(model.predict(chain["Xh"]))
    assert (predicted == chain["yh"].to_numpy()).mean() >= 0.9
    assert list(le.classes_) == sorted([CLASS_A, CLASS_B, NG])


def test_skops_round_trip_keeps_predictions_probabilities_and_feature_names(chain):
    loaded = _load(chain["skops"])
    X = chain["Xh"]
    assert np.array_equal(loaded.predict(X), chain["model"].predict(X))
    assert np.array_equal(loaded.predict_proba(X), chain["model"].predict_proba(X))
    assert list(loaded.feature_names_in_) == list(X.columns)
    assert list(loaded.classes_) == list(chain["model"].classes_)


def test_label_encoder_sidecar_round_trip(chain):
    le = joblib.load(chain["le_path"])
    assert list(le.classes_) == list(chain["le"].classes_)
    assert list(le.inverse_transform([0, 1, 2])) == list(chain["le"].classes_)


# ---------------------------------------------------------------- command-line prediction

def _run_cli(monkeypatch, capsys, model, in_csv, out_csv, *extra):
    # The command line passes `trusted=True`, which skops >= 0.10 refuses (see the xfail test below); every other
    # test of the command line loads the model the way the GUI does so that the rest of its logic is checked.
    monkeypatch.setattr(cli, "load", lambda path, trusted=None: _load(path))
    monkeypatch.setattr(sys, "argv", ["msp_CLI_predict_skops", "--model", str(model), "--input", str(in_csv),
                                      "--out", str(out_csv), *extra])
    cli.main()
    return pd.read_csv(out_csv), capsys.readouterr().out


def _unlabeled(chain, tmp_path, name="unlabeled.csv"):
    df = chain["Xh"].copy()
    df.insert(0, "MS2scan_no", np.arange(1, len(df) + 1))
    path = tmp_path / name
    df.to_csv(path, index=False)
    return df, path


@pytest.mark.xfail(SKOPS_REJECTS_TRUSTED_TRUE, strict=True, raises=TypeError,
                   reason="KNOWN DEFECT (reported, src/ unchanged by this ticket): msp_CLI_predict_skops.main() calls "
                          "skops.io.load(model, trusted=True); skops >= 0.10 raises TypeError for a boolean, and "
                          "src/requirements.txt pins skops 0.13.0. The GUI loads with the file's untrusted types. "
                          "When the command line is fixed this test XPASSes (strict) - remove the marker then.")
def test_cli_loads_a_model_saved_with_the_pinned_skops(chain, tmp_path, monkeypatch, capsys):
    _, in_csv = _unlabeled(chain, tmp_path)
    monkeypatch.setattr(sys, "argv", ["msp_CLI_predict_skops", "--model", str(chain["skops"]), "--input", str(in_csv),
                                      "--out", str(tmp_path / "pred.csv")])
    cli.main()
    assert (tmp_path / "pred.csv").is_file()


def test_cli_adds_prediction_and_confidence_columns(chain, tmp_path, monkeypatch, capsys):
    df, in_csv = _unlabeled(chain, tmp_path)
    out, stdout = _run_cli(monkeypatch, capsys, chain["skops"], in_csv, tmp_path / "pred.csv")
    assert list(out.columns) == list(df.columns) + ["Predicted_Label", "proba_top1", "proba_top2", "margin"]
    assert len(out) == len(df) and list(out["MS2scan_no"]) == list(df["MS2scan_no"])
    proba = chain["model"].predict_proba(chain["Xh"])
    ordered = np.sort(proba, axis=1)
    assert np.allclose(out["proba_top1"], ordered[:, -1])
    assert np.allclose(out["proba_top2"], ordered[:, -2])
    assert np.allclose(out["margin"], ordered[:, -1] - ordered[:, -2])
    assert (out["proba_top1"] >= out["proba_top2"]).all()
    coverage = (out["proba_top1"] >= 0.65).mean()
    assert f"Coverage @ tau=0.65: {coverage:.1%}" in stdout
    assert "Wrote predictions" in stdout


def test_cli_characterisation_writes_the_models_class_values_not_decoded_labels(chain, tmp_path, monkeypatch, capsys):
    """CHARACTERISATION.  The GUI fits on LabelEncoder integers, so for a GUI-style model the command line writes
    those integers in `Predicted_Label`; the names are in the `*_labelencoder.joblib` sidecar, which the command
    line does not read.  Recorded so that a change (e.g. decoding with the sidecar) is a deliberate one."""
    _, in_csv = _unlabeled(chain, tmp_path)
    out, _ = _run_cli(monkeypatch, capsys, chain["skops"], in_csv, tmp_path / "pred.csv")
    assert set(out["Predicted_Label"]) <= {0, 1, 2}
    decoded = chain["le"].inverse_transform(out["Predicted_Label"].to_numpy())
    assert (decoded == chain["yh"].to_numpy()).mean() >= 0.9


def test_cli_class_col_option_renames_the_prediction_column(chain, tmp_path, monkeypatch, capsys):
    _, in_csv = _unlabeled(chain, tmp_path)
    out, _ = _run_cli(monkeypatch, capsys, chain["skops"], in_csv, tmp_path / "pred.csv", "--class-col", "pred_enc")
    assert "pred_enc" in out.columns and "Predicted_Label" not in out.columns


def test_cli_aligns_columns_to_the_training_order_and_ignores_extras(chain, tmp_path, monkeypatch, capsys):
    df, in_csv = _unlabeled(chain, tmp_path)
    base, _ = _run_cli(monkeypatch, capsys, chain["skops"], in_csv, tmp_path / "base.csv")
    shuffled = df[["MS2scan_no"] + list(reversed(list(chain["Xh"].columns)))].copy()
    shuffled["unrelated_column"] = "x"
    sh_csv = tmp_path / "shuffled.csv"
    shuffled.to_csv(sh_csv, index=False)
    out, _ = _run_cli(monkeypatch, capsys, chain["skops"], sh_csv, tmp_path / "shuffled_pred.csv")
    assert list(out["Predicted_Label"]) == list(base["Predicted_Label"])
    assert np.allclose(out["proba_top1"], base["proba_top1"])


def test_cli_fills_a_missing_feature_column_with_zero(chain, tmp_path, monkeypatch, capsys):
    df, _ = _unlabeled(chain, tmp_path)
    missing_col = str(ION_MASSES[1])
    without = df.drop(columns=[missing_col])
    without.to_csv(tmp_path / "without.csv", index=False)
    explicit = df.copy()
    explicit[missing_col] = 0.0
    explicit.to_csv(tmp_path / "explicit.csv", index=False)
    a, _ = _run_cli(monkeypatch, capsys, chain["skops"], tmp_path / "without.csv", tmp_path / "a.csv")
    b, _ = _run_cli(monkeypatch, capsys, chain["skops"], tmp_path / "explicit.csv", tmp_path / "b.csv")
    assert list(a["Predicted_Label"]) == list(b["Predicted_Label"])
    assert np.allclose(a["proba_top1"], b["proba_top1"])


def test_cli_fill_value_for_a_missing_feature_is_zero(tmp_path, monkeypatch, capsys):
    """Feature `z` is 0.0 for Non-glycan and 1.0 for the glycan, so 0.0 and 1.0 land on different sides of the split."""
    X = pd.DataFrame({"z": np.r_[np.zeros(10), np.ones(10)], "w": np.linspace(1.0, 1.1, 20)})
    model, le = _fit(X, np.r_[np.full(10, NG), np.full(10, CLASS_A)])
    skops_dump(model, tmp_path / "z_model.skops")
    pd.DataFrame({"MS2scan_no": [1, 2], "w": [1.05, 1.05]}).to_csv(tmp_path / "unl.csv", index=False)   # no `z` column
    out, _ = _run_cli(monkeypatch, capsys, tmp_path / "z_model.skops", tmp_path / "unl.csv", tmp_path / "o.csv")
    assert set(out["Predicted_Label"]) == {le.transform([NG])[0]}


def test_cli_characterisation_default_drop_cols_removes_the_mass_feature(tmp_path, monkeypatch, capsys):
    """CHARACTERISATION.  The default `--drop-cols` is `MS2scan_no protonatedmass`.  For a model trained WITH the
    mass feature (the production setting), the mass column is dropped from the input and then re-created as 0.0,
    which is also the value the GUI gives Non-glycan rows.  Here the class depends on the mass only, so the default
    call predicts Non-glycan for every row while `--drop-cols MS2scan_no` keeps the mass and predicts the glycan.
    Recorded so that a change of the default is a deliberate one."""
    rng = np.random.default_rng(3)
    n = 20
    mass = np.r_[np.full(n, 1235.44), np.zeros(n)]
    X = pd.DataFrame({"204.0867": rng.uniform(1.0, 1.2, 2 * n), "protonatedmass": mass})
    labels = np.r_[np.full(n, CLASS_A), np.full(n, NG)]
    model, le = _fit(X, labels)
    skops_path = tmp_path / "mass_model.skops"
    skops_dump(model, skops_path)
    unl = pd.DataFrame({"MS2scan_no": range(1, 6), "204.0867": np.full(5, 1.1), "protonatedmass": np.full(5, 1235.44)})
    unl.to_csv(tmp_path / "unl.csv", index=False)
    glycan_code, ng_code = le.transform([CLASS_A])[0], le.transform([NG])[0]

    direct = model.predict(unl[["204.0867", "protonatedmass"]])
    assert set(direct) == {glycan_code}                                         # the model itself uses the mass
    default, _ = _run_cli(monkeypatch, capsys, skops_path, tmp_path / "unl.csv", tmp_path / "d.csv")
    assert set(default["Predicted_Label"]) == {ng_code}                         # mass dropped, zero-filled
    kept, _ = _run_cli(monkeypatch, capsys, skops_path, tmp_path / "unl.csv", tmp_path / "k.csv",
                       "--drop-cols", "MS2scan_no")
    assert set(kept["Predicted_Label"]) == {glycan_code}


# ---------------------------------------------------------------- the tau rule

class _StubModel:
    """Returns fixed probabilities so the rule is checked on exact numbers."""
    def __init__(self, proba):
        self._p = np.asarray(proba, dtype=float)

    def predict_proba(self, X):
        return self._p


def _encoder():
    return LabelEncoder().fit([CLASS_B, NG])      # classes_ = ['F1H5N4S2', 'Non-glycan']


def test_predict_with_threshold_demotes_low_confidence_glycans_only():
    proba = [[0.90, 0.10],    # glycan, confident          -> kept
             [0.60, 0.40],    # glycan, below tau          -> Non-glycan, abstained
             [0.65, 0.35],    # glycan, exactly tau        -> kept (rule is maxp < tau)
             [0.30, 0.70],    # model already says NG      -> NG, not an abstention
             [0.40, 0.60]]    # NG with maxp below tau     -> NG, not an abstention
    y, maxp, p, abst = predict_with_threshold(_StubModel(proba), np.zeros((5, 1)), _encoder(), tau=0.65)
    assert list(y) == [CLASS_B, NG, CLASS_B, NG, NG]
    assert list(abst) == [False, True, False, False, False]
    assert np.allclose(maxp, [0.90, 0.60, 0.65, 0.70, 0.60])
    assert np.array_equal(p, np.asarray(proba))


def test_predict_with_threshold_default_tau_is_055():
    # 0.60 survives the default tau (0.55) but not tau 0.65
    proba = [[0.60, 0.40]]
    assert list(predict_with_threshold(_StubModel(proba), np.zeros((1, 1)), _encoder())[0]) == [CLASS_B]
    assert list(predict_with_threshold(_StubModel(proba), np.zeros((1, 1)), _encoder(), tau=0.65)[0]) == [NG]


def test_tau_sweep_abstention_rate_is_monotone_and_hand_checked():
    proba = [[0.95, 0.05], [0.85, 0.15], [0.70, 0.30], [0.55, 0.45], [0.20, 0.80]]
    y_true = [CLASS_B, CLASS_B, CLASS_B, CLASS_B, NG]
    s = tau_sweep_summary(_StubModel(proba), np.zeros((5, 1)), y_true, _encoder(), taus=(0.5, 0.65, 0.8, 0.9))
    assert list(s["tau"]) == [0.5, 0.65, 0.8, 0.9]
    # rows demoted to NG: tau 0.5 -> only the row the model already calls NG; 0.65 -> +0.55; 0.8 -> +0.70; 0.9 -> +0.85
    assert list(s["abstention_rate"]) == pytest.approx([0.2, 0.4, 0.6, 0.8])
    assert list(s["abstention_rate"]) == sorted(s["abstention_rate"])
    # accuracy on the four true glycans: kept ones are correct
    assert list(s["accuracy_on_glycan"]) == pytest.approx([1.0, 0.75, 0.5, 0.25])


def _v14_demote_expressions():
    text = V14.read_text(encoding="utf-8")
    found = re.findall(r"^\s*demote = (.+)$", text, flags=re.M)
    if len(found) != 2:
        pytest.skip(f"v14 anchor moved: expected 2 `demote = ...` lines, found {len(found)} (coverage LOST, update this test)")
    return found


# tau 0.5 and margin 0.25 are exact in binary floating point, so the edge rows are exact
TAU, MARGIN, MAJORITY = 0.5, 0.25, 2
DEMOTE_CASES = [
    # predicted class, max proba, proba of the majority class, expected demotion
    (0, 0.375, 0.25, True),     # glycan, below tau, majority within the margin
    (1, 0.375, 0.125, True),    # majority exactly maxp - margin: the comparison is >=
    (0, 0.375, 0.0625, False),  # majority too far below the best
    (0, 0.5, 0.4375, False),    # confidence exactly tau: not demoted
    (0, 0.875, 0.0625, False),  # confident glycan
    (2, 0.375, 0.375, False),   # already the majority class
]


@pytest.mark.parametrize("which", [0, 1], ids=["train_evaluation", "predict_tab"])
def test_gui_margin_aware_demotion_expression(which):
    """Runs the exact expression the GUI uses (train-time evaluation and Predict tab), read from the v14 source."""
    expr = _v14_demote_expressions()[which]
    ast.parse(expr, mode="eval")
    pred = np.array([c[0] for c in DEMOTE_CASES])
    ns = {"y_pred": pred, "y_pred_enc": pred, "majority_idx": MAJORITY, "tau": TAU, "margin": MARGIN,
          "maxp": np.array([c[1] for c in DEMOTE_CASES]), "p_major": np.array([c[2] for c in DEMOTE_CASES])}
    demote = eval(expr, {"__builtins__": {}}, ns)     # noqa: S307 - the expression comes from this repository's own source
    assert list(demote) == [c[3] for c in DEMOTE_CASES]


def test_gui_threshold_defaults_are_off_with_tau_065_and_margin_005():
    text = V14.read_text(encoding="utf-8")
    if '"threshold": {"enabled": False, "tau": 0.65, "margin": 0.05}' not in text:
        pytest.skip("v14 anchor moved: default threshold dictionary not found (coverage LOST, update this test)")
    assert '"thresholds": {"tau": 0.65, "margin": 0.05}' in text


# ---------------------------------------------------------------- prediction report

REPORT_CLASSES = ["A", "B", "NG"]
REPORT_PROBA = [[0.875, 0.125, 0.0],      # accepted
                [0.5, 0.375, 0.125],      # max == tau but margin 0.125 < 0.25 -> low_margin
                [0.4, 0.35, 0.25],        # below tau and margin 0.05          -> low_conf_and_margin
                [0.25, 0.625, 0.125],     # accepted
                [0.4375, 0.1875, 0.375],  # below tau and margin 0.0625        -> low_conf_and_margin
                [0.625, 0.375, 0.0]]      # margin exactly 0.25: not below the margin -> accepted


def _report_frame():
    df = pd.DataFrame(REPORT_PROBA, columns=REPORT_CLASSES)
    df["pred_label"] = ["A", "A", "A", "B", "A", "A"]
    df["MS2scan_no"] = [1, 2, 3, 4, 5, 6]
    df["experiment_title"] = "exp"
    df["sample_name"] = ["s1", "s1", "s1", "s2", "s2", "s2"]
    return df


def test_summarize_predictions_acceptance_rules_and_counts():
    params = PR.ReportParams(tau=0.5, margin=0.25, topk=3)
    rows, class_sum, by_sample, run = PR.summarize_predictions(_report_frame(), class_names=REPORT_CLASSES,
                                                               proba_cols=REPORT_CLASSES, params=params)
    assert list(rows["accepted"]) == [True, False, False, True, False, True]
    assert list(rows["reject_reason"]) == ["", "low_margin", "low_conf_and_margin", "", "low_conf_and_margin", ""]
    assert list(rows["second_label"])[0] == "B" and list(rows["second_label"])[3] == "A"
    assert rows["margin"].iloc[0] == pytest.approx(0.75)
    assert rows["entropy"].iloc[0] == pytest.approx(0.5435644, abs=1e-6)
    assert rows["topk_labels"].iloc[0] == "A,B,NG"

    cs = class_sum.set_index("pred_label")
    assert (cs.loc["A", "pred_count"], cs.loc["A", "accepted_count"]) == (5, 2)
    assert (cs.loc["B", "pred_count"], cs.loc["B", "accepted_count"]) == (1, 1)
    assert cs.loc["A", "accepted_rate_within_class"] == pytest.approx(0.4)

    bs = by_sample.set_index("sample_name")
    assert list(bs.loc["s1", ["n_total", "n_accepted", "n_rejected"]]) == [3, 1, 2]
    assert list(bs.loc["s2", ["n_total", "n_accepted", "n_rejected"]]) == [3, 2, 1]

    assert (run["n_total"], run["n_accepted"], run["n_rejected"]) == (6, 3, 3)
    assert (run["tau"], run["margin"]) == (0.5, 0.25)
    assert run["low_conf_count"] == 2 and run["ambiguous_count"] == 3
    assert run["class_names"] == REPORT_CLASSES


def test_summarize_predictions_defaults_are_tau_060_margin_005():
    p = PR.ReportParams()
    assert (p.tau, p.margin, p.topk) == (0.60, 0.05, 5)


def test_summarize_predictions_needs_probabilities():
    with pytest.raises(ValueError, match="proba_cols"):
        PR.summarize_predictions(_report_frame(), class_names=REPORT_CLASSES)
    with pytest.raises(ValueError, match="len\\(class_names\\)"):
        PR.summarize_predictions(_report_frame(), class_names=["A", "B"], proba_cols=REPORT_CLASSES)


def test_write_prediction_report_files_and_review_lists(tmp_path):
    params = PR.ReportParams(tau=0.5, margin=0.25, topk=3)
    rows, class_sum, by_sample, run = PR.summarize_predictions(_report_frame(), class_names=REPORT_CLASSES,
                                                               proba_cols=REPORT_CLASSES, params=params)
    arts = PR.write_prediction_report(tmp_path / "report", rows, class_sum, by_sample, run, params)
    for path in (arts.pred_with_conf, arts.class_summary, arts.by_sample_summary, arts.review_low_conf,
                 arts.review_low_margin, arts.review_ambiguous, arts.run_summary_json, arts.run_summary_md):
        assert Path(path).is_file()
    assert len(pd.read_csv(arts.pred_with_conf)) == 6
    assert list(pd.read_csv(arts.review_low_conf)["MS2scan_no"]) == [3, 5]            # max_proba < tau
    assert list(pd.read_csv(arts.review_low_margin)["MS2scan_no"]) == [2, 3, 5]       # margin < margin
    assert json.loads(Path(arts.run_summary_json).read_text(encoding="utf-8"))["n_total"] == 6
    assert "Total: 6" in Path(arts.run_summary_md).read_text(encoding="utf-8")


def test_report_has_no_file_names_from_the_environment(tmp_path):
    """The report writes the folder it was given only; no absolute path of the input appears in its files."""
    params = PR.ReportParams(tau=0.5, margin=0.25, topk=3)
    rows, class_sum, by_sample, run = PR.summarize_predictions(_report_frame(), class_names=REPORT_CLASSES,
                                                               proba_cols=REPORT_CLASSES, params=params)
    out = tmp_path / "report"
    PR.write_prediction_report(out, rows, class_sum, by_sample, run, params)
    for f in out.iterdir():
        assert str(tmp_path) not in f.read_text(encoding="utf-8")
