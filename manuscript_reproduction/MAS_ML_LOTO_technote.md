# MAS-ML Leave-One-Tissue-Out Re-analysis — Technical Note

**Deliverable:** `mas_ml_loto_manuscript_analysis.py` (v1.0.0) + `MAS_ML_LOTO_results/`
**Context:** Bioinformatics (Original Paper) manuscript, Result 2 (MAS-ML).
**Core rule honored:** GlycoMSP v1.10 source is *not* modified, imported, or executed. The pipeline functions that shape the numbers are ported verbatim from `src/mspfileloaderv14.py` (line anchors in the script header) so results reconcile with the GUI tool.

---

## 1. Purpose

Provide the ML-evaluation block a Bioinformatics reviewer expects for MAS-ML, which the GUI did not previously emit: a within-dataset held-out test, an **independent** cross-tissue test, and cross-validation sets — computed on the already-exported MAS *trainable* CSVs, with numbers faithful to GlycoMSP. It replaces the prior "train and predict on the same file" demonstration, which is not an independent evaluation.

## 2. Input files

Three GlycoMSP MAS trainable (`*_merged.csv`) exports — zebrafish brain / intestine / ovary N-glycans — plus a feature definition. The script accepts either a `features.json` (list) or a GlycoMSP `training_run.json` (`column_order`). For this run the canonical feature set is `protonatedmass` + 282 fragment-ion m/z columns (283 features, mass ON), taken from `6_1_brain/training_run.json`. All three tissues share an identical 282-ion set, so cross-tissue alignment is exact.

Reference for validation: `6_1_brain/..._rf_performance.txt` (the frozen GUI run).

## 3. Preprocessing assumptions

- Feature columns in the merged CSV are **already** `log10(I)+1`-transformed by GlycoMSP at export (`ml_ng_utils.build_features_from_peaks_log10_plus1`); the script does **not** re-transform.
- Each tissue is reindexed to the canonical `column_order` (missing columns → 0), mirroring the GUI's column-lock + `reindex(fill_value=0)`.
- `protonatedmass` is a feature when mass is ON; for `Non-glycan` rows it is forced to `0.0` and remaining NaNs → `0.0` (ported `_prep_ms1_block`, v14:6563). `--no-mass` drops it.
- Metadata columns are never features: `unique_ID`, `WURCS`, `GlyToucan ID`, `Glycanannotation2`, `IUPACname(optional)`, plus `UID`/`Origin_*` provenance (ported `feature_exclude` v14:6662 + `_drop_meta_and_get_X` v14:6507).
- `random_state=42` throughout (matches the hard-coded MAS workflow seed).

## 4. Feature and label handling

- **Label column:** `Structure`, normalized via the ported `normalize_structure_for_training` (v14:6555). MAS labels here are already manual strings (e.g. `F1H5N4S2`), so normalization is effectively passthrough; the real `pretrain_normalizer` is imported read-only if available.
- **Row identifier:** `unique_ID` (GlycoMSP's MAS key), falling back to `UID` then a generated stable id. Carried through to per-row outputs for scan-level traceability.
- **Features:** the 283-column canonical set; `feature_columns_used.json` records the exact order used.

## 5. Evaluation design

- **A — within-tissue stratified hold-out** (per tissue): ported `balance_and_split` Mode B — 80/20 stratified split, `min_count=5` global rare-class drop, then `Non-glycan` capped to 3× max-minority **on the training fold only** (test uncapped). This is the formal in-sample ML evaluation and the validation gate.
- **B — leave-one-tissue-out (LOTO)**: three folds — train on two tissues, predict the fully held-out third. This is the **independent test set**. *Framing note:* this is **leave-one-GROUP-out** across 3 biological tissues, NOT sample-level LOOCV (which Bioinformatics editorially rejects). State the special circumstance: only three expert-annotated tissues exist; each held-out tissue is a genuinely independent biological sample.
- **C — pooled coverage + stratified k-fold CV** (k=5 default): the "cross-validation sets" the policy requires, plus the coverage demonstration (pooling raises the count of trainable compositions). Train-only `Non-glycan` capping is applied per fold.

Bioinformatics asks for both cross-validation sets **and** an independent test set, and says CV averages alone are insufficient. This design delivers both; the LOTO fold table is the headline independent-test result.

## 6. Seen / unseen class logic

For each LOTO fold, the model's class set = compositions surviving `min_count` in the **training** pool. Held-out test rows are partitioned into:

- **seen glycan** — true `Structure` in the model class set, ≠ `Non-glycan` → **primary metrics**;
- **unseen glycan** — true `Structure` absent from the model class set → reported **separately** as out-of-label-space, **never** counted as ordinary false negatives (supervised RF cannot predict an unseen label);
- **non-glycan** — handled in the secondary metric only.

`loto_seen_unseen_classes.csv` gives per-class train/test counts and seen-status per fold.

## 7. Metrics to report

- **Primary (manuscript):** macro precision / recall / F1 and accuracy on **seen glycan classes only**. The macro label set is the union of true and predicted labels on the evaluated rows, **excluding `Non-glycan`**. A true-glycan row predicted `Non-glycan` is kept and counts as a miss against that glycan class's recall (`Non-glycan` is never a scored class in the primary).
- **Secondary (operational, `*_seen_all`):** macro precision / recall / F1 and accuracy over seen-glycan + `Non-glycan` rows, label set = union of true and predicted (may include `Non-glycan`). `Non-glycan` here is an estimated negative generated by GlycoMSP (low ion-score / no-hit), **not** ground truth — reported separately, never folded into the primary.
- **Diagnostic (`*_training_label_space_diagnostic`):** the earlier macro over the full training label space; kept for transparency, **not** a manuscript metric.
- Per-class confusion matrices for seen-glycan classes; k-fold mean ± SD.

### Results obtained (this run)

**A — within-tissue hold-out (glycan-only):**

| tissue | n_test | accuracy | macro-F1 |
|---|---|---|---|
| brain | 219 | 0.957 | **0.931** |
| intestine | 144 | 0.936 | 0.935 |
| ovary | 196 | 0.978 | 0.952 |

Brain reproduces the frozen GUI run (`acc_all 0.9909`, glycan macro-F1 `0.9306`, 46 glycan test rows) → port faithful.

**B — LOTO independent test (seen-glycan, primary):**

| held-out tissue | seen-glycan rows | unseen rows | accuracy | macro-F1 |
|---|---|---|---|---|
| brain | 173 | 115 | 0.607 | 0.635 |
| intestine | 152 | 37 | 0.743 | 0.599 |
| ovary | 98 | 152 | 0.908 | 0.832 |

Macro metrics are averaged over the union of true and predicted seen-glycan labels in each held-out evaluation, excluding Non-glycan (the same label convention as the within-tissue and pooled-CV sections). The earlier all-trainable-label-space values (brain 0.46 / intestine 0.33 / ovary 0.49) are retained as `*_training_label_space_diagnostic` columns for transparency and are **not** manuscript metrics.

**C — pooled 5-fold CV (seen-glycan):** macro-F1 **0.904 ± 0.038**, accuracy 0.928 ± 0.029. Pooled coverage: 2908 rows, 59 glycan labels, **44** meeting min support (vs. ~20 per single tissue).

**Interpretation (for drafting, not asserted by the script):** in-sample MAS-ML is strong (macro-F1 ~0.93); cross-tissue transfer is lower and tissue-dependent (0.60–0.83), reflecting non-identical glycomes and a partially shared label space rather than classifier failure. Ovary recovers best because the brain+intestine training pool covers most of its shared compositions. Pooling restores coverage. This supports the thesis claim that MAS-ML performance depends on accumulated annotated spectra and label balance.

## 8. Output files / tables / figures

`input_audit.csv`, `feature_columns_used.json`, `within_dataset_holdout_summary.csv`, `loto_fold_summary.csv`, `loto_seen_unseen_classes.csv`, `loto_predictions_with_truth.csv` (per-row, with `max_proba` + margin for traceability), `split_manifest.csv` (every row's train/test role per evaluation), `pooled_coverage_summary.csv`, `pooled_class_counts.csv`, `pooled_stratified_5fold_cv.csv`, `confusion_*` matrices, and figures `loto_macro_f1.png`, `loto_seen_unseen.png`. `README.md` records versions and settings.

## 9. Interpretation rules (manuscript-safe wording)

Use: "within-dataset held-out evaluation", "leave-one-tissue-out (leave-one-group-out) independent test", "seen glycan classes", "out-of-label-space held-out compositions", "pooled class coverage summary". Do **not** call pooled-model prediction on the same pooled rows independent validation. Do not present spectrum counts as abundance. Avoid "universal prediction" / "de novo structure discovery". Disclose AI assistance in the cover letter / Methods per Bioinformatics policy.

## 10. Explicitly out of scope

CGA/CGA-ML validation (candidate labels, not ground truth); per-spectrum manual TP/FP/FN re-validation; hyperparameter tuning / model search; algorithmic novelty claims; any edit to GlycoMSP v1.10 source. If a reviewer demands exhaustive CGA-ML validation, that is a separate, scoped effort.

---

### Reproduce

```bash
python mas_ml_loto_manuscript_analysis.py \
  --brain   6_1_brain/zf_sPerMeNG_brain_merged_20260519_224028.csv \
  --intestine 6_1_intestine/zf_sPerMeNG_intestine_merged_20260519_225627.csv \
  --ovary   6_1_ovary/zf_sPerMeNG_ovary_merged_20260519_230200.csv \
  --features 6_1_brain/training_run.json \
  --outdir  MAS_ML_LOTO_results \
  --validate-brain 6_1_brain/zf_sPerMeNG_brain_merged_20260519_224028_withMass_rf_performance.txt
```

Sandbox used sklearn 1.7.2; the frozen reference was sklearn 1.2.0 (cross-version reproducibility was established byte-identical during the v1.10 freeze sprint). For exact manuscript numbers, re-run in your frozen environment (sklearn 1.2.0).
