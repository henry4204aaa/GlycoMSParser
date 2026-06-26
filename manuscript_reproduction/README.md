# MAS-ML leave-one-tissue-out (LOTO) — manuscript reproduction

Standalone analysis used to produce the **MAS-ML evaluation in Result 2 / Methods 2.8** of the GlycoMSP manuscript (zebrafish N-glycan within-tissue hold-out, leave-one-tissue-out, and pooled cross-validation). It reads GlycoMSP-exported MAS trainable CSVs and reproduces the reported numbers **without modifying, importing, or running GlycoMSP core code**. The pipeline functions that shape the numbers (`balance_and_split`, the MS1/mass block, label normalization, meta-column dropping) are ported verbatim from frozen `mspfileloaderv14.py` so results reconcile with the GUI tool. Corresponds to Supplementary Table S2c, package ID = 7.

## Contents

- `mas_ml_loto_manuscript_analysis.py` — the analysis script.
- `MAS_ML_LOTO_technote.md` — design, inputs, evaluation logic, metric definitions, interpretation rules.
- `MAS_ML_LOTO_validation_report.md` — evidence that the script reproduces the GlycoMSP GUI exactly (and is environment-independent).
- `MAS_ML_LOTO_results/` — output bundle (summary tables, per-row predictions, confusion matrices, figures, run README).

## Requirements

Python 3.10+ with `pandas`, `numpy`, `scikit-learn` (and `matplotlib` for the optional figures). The reported numbers were generated with scikit-learn 1.2.0; outputs were verified byte-identical on scikit-learn 1.7.2 (see validation report).

## Usage

```bash
python mas_ml_loto_manuscript_analysis.py \
  --brain     zf_sPerMeNG_brain_merged.csv \
  --intestine zf_sPerMeNG_intestine_merged.csv \
  --ovary     zf_sPerMeNG_ovary_merged.csv \
  --features  training_run.json \
  --gms-src   /path/to/GlycoMSP/src \
  --outdir    MAS_ML_LOTO_results
```

- `--features` accepts a GlycoMSP `training_run.json` (uses `column_order`) or a `features.json` list; it fixes the canonical feature columns (`protonatedmass` + 282 ion m/z = 283).
- `--gms-src` is optional; pointing it at GlycoMSP `src/` lets the script use GlycoMSP's label normalizer so class names match the tool (e.g. `KDN` rather than `K`). Metric values are unaffected by this flag.
- Input CSVs and `training_run.json` are the MAS trainable exports provided in this package (and in the smoke_6 fixtures of the repository).

## Expected results (faithfulness anchors)

- **Within-tissue hold-out, glycan macro-F1:** brain 0.931, intestine 0.935, ovary 0.952. The brain run reproduces the GUI `rf_performance.txt` exactly (test support 219; all-class accuracy 0.99; glycan macro-F1 0.93).
- **Leave-one-tissue-out, seen-glycan macro-F1:** brain 0.635, intestine 0.599, ovary 0.832 (accuracy 0.607 / 0.743 / 0.908).
- **Pooled stratified 5-fold CV:** glycan accuracy 0.928 ± 0.029, glycan macro-F1 0.904 ± 0.038; pooling raises trainable glycan labels from 14–25 per tissue to 43.

## Notes

- Primary metrics are glycan-only (Non-glycan, an operational/estimated negative, is excluded and reported separately). Held-out tissue-specific compositions absent from training are reported as out-of-label-space, not as classification errors. See the tech note for the exact metric conventions.
- "LOTO" here is a leave-one-group-out evaluation across three biologically distinct tissues, not sample-level leave-one-out cross-validation.
- RF settings mirror GlycoMSP v1.10: 400 trees, `class_weight="balanced"`, `random_state=42`, `min_samples_split=2`, `min_samples_leaf=1`; class-eligibility filter ≥5 spectra; Non-glycan capped at 3× the majority glycan class on the training fold only.
- Distributed with GlycoMSP under the MIT License.
