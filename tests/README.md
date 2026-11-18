# GlycoMSP automated tests

The pytest suite in this folder runs from a clean clone with **no data files**:

```
pip install -r src/requirements.txt      # Python 3.10-3.12 (3.13: src/requirementspy313.txt)
pip install pytest
python3 -m pytest tests -q -rs
```

Developed and run on Python 3.12 / pytest 9 on macOS; runs on Windows, Linux or CI have not been recorded here. Every test builds its own input in a temporary folder (synthetic mzML, a fake Thermo reader, small tables) and each file finishes in seconds. Nothing is downloaded; files are written to temporary folders.

A test that cannot run on the machine **skips with a reason** (`-rs` prints it); none should fail. One test is a strict expected failure (`xfail`): `test_ml_features_train_predict.py::test_cli_loads_a_model_saved_with_the_pinned_skops` records a known defect of the prediction command line with current skops; it will fail loudly once the defect is fixed, so that the marker is removed on purpose. Two tests in the same file are labelled CHARACTERISATION: they record what the command line does today and are not a claim that it is desired.

| Skip | Reason |
|---|---|
| a few GUI tests ("real Tk") | no display (headless server); the same logic is tested without Tk elsewhere |
| `test_reader_diag_mzml.py::test_pilot_excerpt_rebuilt` | optional local excerpt of a real file, not distributed |
| `test_s1_glytoucan_resolver.py::test_writer_default_path_equals_v110_baseline_from_git` | compares the MAS writer with the v1.10 development commit `f27a872`, which exists only in the development repository |
| SIGTERM test in `test_db1_manage_references.py` | Windows only |

## What is covered

| Workflow stage | Files |
|---|---|
| RAW and mzML conversion, charge and polarity handling, diagnostics, exit codes | `test_reader_diag_*.py`, `_synth_mzml.py` (synthetic mzML builder) |
| GUI wiring for conversion (options dialog, reader result) | `test_v14_conversion_options.py`, `test_v14_reader_result.py` |
| Composition library enumeration (N-glycan constraints, O-glycan cores) | `test_compnewv4_enumeration.py` |
| Features, Random Forest, `.skops` round trip, command-line prediction, tau rule, prediction report | `test_ml_features_train_predict.py` |
| Precursor mass calculation (Native, PerMe, SO3/PO3H) | `test_native_precursor_mass.py`, `test_cli_precursor_mass.py`, `test_b05_kdn_mass_smoke.py` |
| CGA: Score A, tolerances in ppm or Da, recording of settings | `test_b24_score_counter.py`, `test_cga_tolerance.py`, `test_cga_tolerance_record.py` |
| Score B composition tokenizer | `test_hygiene_scoreb_tokens.py` |
| Trainable CSV command line, code hygiene guards | `test_hygiene_cli_import.py`, `test_hygiene_shadowed_defs.py` |
| GlyTouCan / WURCS reference database and Manage References | `test_s1_*.py`, `test_db1_manage_references.py`, `fixtures/db1_stub_cli.py` (stub command line) |
| Random Forest timing script | `test_rf_timing.py` (see `tools/validation/`) |
| Controlled validation scripts: CGA tolerance decoy count, shared-peak count | `test_cga_tolerance_decoy_count.py` (31 tests; the scripts are in `tools/validation/`) |

`conftest.py` records environment overrides so one test file cannot leak into another. `_v14_extract.py` executes `src/mspfileloaderv14.py` up to the creation of the Tk root so that GUI helpers can be tested headless. Several files therefore read the v14 source text and look for fixed anchors; if a refactor moves an anchor the test **skips with "coverage LOST"** rather than passing silently.

## Script-style checks (not collected by pytest)

Five older checks run as scripts, need no data files, and finish in under a second each:

```
python3 tests/test_b01_kdn_smoke.py
python3 tests/test_b03_apply_to_exp_smoke.py
python3 tests/test_b04_excepthook_smoke.py
python3 tests/test_b06_status_label_smoke.py      # needs a display
python3 tests/test_b17_glytoucan_columns_smoke.py
```

Some of them mirror small GUI helpers instead of importing them ("keep in lockstep"), so they can drift from the product; treat them as regression guards for the fixed bug, not as proof of current behaviour.

## Not part of the public set

These files exist in the development repository and are listed in `tests/public_exclude.txt`.

* `test_smoke_02_raw_to_csv.py`, `test_smoke_02a_mzml_to_csv.py`, `test_smoke_03_metadata_json.py`, `test_smoke_04_mas_to_trainable.py`, `test_smoke_06_ml_training.py`, `test_smoke_08_cga_insilico.py`, `test_smoke_09_cga_composition_fitting.py`: smoke checks against reference outputs of the manuscript runs, about 8 GB, not distributed. Four of them (02, 06, 08, 09) compare stored reference outputs and do not execute GlycoMSP code; 02a, 03 and 04 re-run the converter or the validator on the reference inputs.
* `test_json_upgrade.py`: tests of the JSON format upgrade for versions before v1.10; it does not apply to the released code.

**Regenerating the withheld reference outputs.** They were derived from the public zebrafish deposit (GlycoPOST GPST000224). Converting the deposited RAW files with v1.14 (`python src/msp_extract_cli.py convert --raw <file.raw> --out <folder> --alias <name>`; RAW needs Windows and Thermo MSFileReader, mzML works anywhere with `--mzml`) regenerates the converter outputs; the later stages (annotation workbook, trainable CSV, Random Forest, CGA) are GUI runs on those files with the settings listed in `docs/TEST_REGISTER.md`, section "Regenerating the withheld reference outputs", where the commands and the expected values are given. The smoke scripts then compare a regenerated output with the stored one; they do not create it. One reference, the O-glycan run used in smoke 09, comes from an unpublished dataset and cannot be regenerated from the deposit.

## Controlled validation scripts

`tools/validation/` holds scripts that a person or an AI assistant runs on chosen inputs (decoy-library chance baseline for CGA tolerances, shared-peak counts, Random Forest timing). They are validation tools, not application code; the tests of the decoy-count script (`tests/test_cga_tolerance_decoy_count.py`) are part of the normal `python3 -m pytest tests` run (see `tools/validation/README.md`).
