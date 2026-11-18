# Known Issues

**DRAFT — internal reference**
**Generated:** 2026-03-18 | **Updated:** 2026-03-31; Issues 16–20 added 2026-10-20 (READER-DIAG CP-E); Issue 22 added 2026-10-27 (V14-FOLLOWUP CP-A); Issues 17, 18, 21 status and Issues 23–26 added 2026-10-27 (V14-FOLLOWUP CP-D); Issue 17 rationale and Issue 27 added 2026-10-30 (V14-FOLLOWUP manual acceptance); Issue 27 status and Issue 28 added 2026-11-01 (CGA-TOLERANCE CP-D); Issues 30–34 added 2026-11-18 (TESTS-PUBLIC)
**Sources:** ARCHITECTURE_VALIDATION.md (2026-03-18), CGA_SCORING_ANALYSIS.md (2026-02-06), session history, v14 analysis (2026-03-31)

---

## Issue 1: score_counter Anchor Shadowing (Active Production Bug)

**Module:** `msp_insilicomarker_withGPT.py`
**Lines:** 293 (simple version), 1241 (anchor-gated version)
**Severity:** Behavioral — affects all Route A annotation scores

### Description

`score_counter` is defined twice in the module:

```python
# Line 293 — simple hit ratio
def score_counter(hits, ionlist):
    if len(ionlist) == 0:
        return 0.0
    return len(hits) / len(ionlist)

# Line 1241 — anchor-gated (active due to Python late binding)
def score_counter(hits_mz, ionlist_mz):
    if not has_anchors(hits_mz):
        return 0.0
    if len(ionlist_mz) == 0:
        return 0.0
    return len(hits_mz) / len(ionlist_mz)
```

Python's module-level late binding means the second definition at line 1241 **shadows** the first. Every call to `score_counter()` at runtime uses the anchor-gated version.

The anchor gate (`has_anchors` at line 1238) returns a count of how many anchor ions from `ANCHORS = [204.087, 366.140, 512.197]` (line 1235) appear within `ANCHOR_TOL = 0.02` Da in the matched hits. `score_counter` (line 1243) returns 0.0 if this count is `< 2`. If fewer than 2 anchor ions match, the score is 0.0 regardless of other hits.

The anchor values `[204.087, 366.140, 512.197]` correspond approximately to HexNAc oxonium ions common in permethylated glycan MS2 spectra. The comment in the source notes these are placeholder values that should be adjusted for the specific derivatization and adduct mode.

### Impact

All Route A annotation scores are anchor-gated. Scans lacking the required anchor ions receive a score of 0.0 and will not be annotated as glycan candidates. This may produce conservative/restrictive annotation in conditions where anchor ions are absent.

ARCHITECTURE.md incorrectly tagged the line-1241 function as `[DEAD]` logic — it is actively enforced.

### Identifier Correction

ARCHITECTURE.md referred to `ANCHOR_IONS` (dict). The actual name is `ANCHORS` (list) at line 1235. `ANCHOR_TOL` at line 1236 is the associated tolerance. There is no `anchor_required` parameter anywhere in the module.

**Status:** Known, working as designed (anchor enforcement is intentional), but was undocumented and misclassified.

---

## Issue 2: log10 Formula — Status as of 20260122

**Modules:** `ml_ng_utils.py`, `mspvalidator_merger.py`, `msp_insilicomarker_withGPT.py`, `mspfileloaderv14.py`
**Severity:** RESOLVED — formula is correct in all active code as of 20260122

### Description

The correct formula for intensity feature transformation is:
```
feature = log10(I) + 1       # CORRECT
```
Not:
```
feature = log10(I + 1)       # WRONG — different mathematical meaning
```

These differ when I is small: `log10(10) + 1 = 2.0` vs `log10(11) ≈ 1.041`. For I = 0, the correct formula should return baseline 1.0 (code guards I ≤ 0 with baseline).

### Current Status

All four production locations now use the correct formula:

| File | Line | Formula | Status |
|------|------|---------|--------|
| `ml_ng_utils.py` | 51 | `math.log10(I) + 1.0` | Correct |
| `mspvalidator_merger.py` | 345 | `np.log10(intensity) + 1` | Correct |
| `msp_insilicomarker_withGPT.py` | 490 | `log10(I) + 1` | Correct (fixed 20260122) |
| `mspfileloaderv14.py` | ~9880 | `log10(I) + 1` | Correct (carried from v13 fix 20260122) |

### Remaining Inconsistency

A **comment** at marker line 396 still says `log10(I+1)` but the actual code implements the correct formula. This comment is misleading but does not affect behavior.

ARCHITECTURE.md Section 10 claimed lines 9424 and 9525 of v13 contained incorrect code — this is false. The fix was applied before the validation report was written. Lines 9424/9525 do not contain log10 operations.

---

## Issue 3: pretrain_normalizer.py Duplicate COMP_ORDER

**Module:** `pretrain_normalizer.py`
**Lines:** 13 (first definition), 76 (second definition)
**Severity:** Low — affects label canonicalization order; risk of unexpected composition string format

### Description

`COMP_ORDER` is defined twice with different orderings:

```python
# Line 13:
COMP_ORDER = ["H", "N", "S", "G", "K", "F"]  # F = Fuc is last

# Line 76 (shadows line 13 at runtime):
COMP_ORDER = ["H", "N", "F", "S", "G", "K"]  # F = Fuc is third
```

At runtime, the line-76 definition `["H", "N", "F", "S", "G", "K"]` is active. Labels produced by `_canon_structure_label` will use Fuc (F) as the third component, yielding strings like `"H5N4F1S1"` rather than `"H5N4S1G0K0F1"`.

Similarly, `_canon_structure_label` is defined twice (lines 34 and 78). The line-78 definition is active.

**Update 2026-10-01 (hygiene H3):** the shadowed first `_canon_structure_label` was removed; the effective definition is unchanged. The duplicate `COMP_ORDER` remains (the first is now read only by the unused helper `_compose_label_from_counts`).

### Impact

If any caller relies on the line-13 `COMP_ORDER` order, it will receive unexpected label formats. The two versions are internally consistent within each definition block, but the duplicate is confusing and a maintenance hazard.

**Recommendation:** Remove the line-13 definitions and keep only the line-76 block.

---

## Issue 4: mspvalidator_merger.py Duplicate sample_real_negatives

**Module:** `mspvalidator_merger.py`
**Lines:** 749 (first), 846 (second)
**Severity:** Low — second definition is active, but first exists and is dead

### Description

`sample_real_negatives` is defined twice. Python late binding means the line-846 version is the active one. The line-749 version is unreachable and should be removed.

---

## Issue 5: compnewv4.py Duplicate decide_functional_groups

**Module:** `compnewv4.py`
**Lines:** 201 (first), 258 (second)
**Severity:** Low — second definition is active

### Description

`decide_functional_groups` (N-glycan variant) is defined at line 201 and again at line 258. Python uses the line-258 version. The line-201 version is shadowed.

`compcheck` is similarly defined at lines 915 and 979; line-979 version is active.

---

## Issue 6: mspfileloaderv14.py Duplicate _ensure_method_stub

**Module:** `mspfileloaderv14.py` (also present in v13)
**Lines:** 2415 (first), 2970 (second) [v13: 2320, 2875]
**Severity:** Low — may cause confusion during maintenance

### Description

`_ensure_method_stub` is a nested function defined twice within `open_prepare_dataset_window()`. The second definition at line 2970 (v14) is in the active inner scope and is the version used by downstream logic (`_ensure_methods_before_saving_exp_v1`, etc.). The first definition at line 2415 (v14) is also present.

Also, `build_final_path` is defined at both lines 1818 and 1830 (v14) at module level; the second definition is active.

---

## Issue 7: mspfileloaderv14.py Active Debug Print Statements

**Module:** `mspfileloaderv14.py` (carried from v13)
**Lines:** ~2785, ~2796, ~2830 [v13: 2690, 2701, 2735]
**Severity:** Low — informational noise in production

### Description

Three debug print statements remain active inside `import_experiment_json`:

| String | Line |
|--------|------|
| `"[LEGACY DEBUG v1]"` | 2690 |
| `"[ML-1 DEBUG]"` | 2701 |
| `"[LEGACY DEBUG]"` | 2735 |

These will print to console on every experiment import operation. They should be gated with a `debug=False` parameter or removed before public release.

---

## Issue 8: msp_CLI_make_trainable_with_negatives.py Legacy Import

**Module:** `msp_CLI_make_trainable_with_negatives.py`
**Line:** 48
**Severity:** RESOLVED 2026-10-01 (hygiene H2) — was Medium, ImportError at runtime if `--add-decoys` is used. The import had since moved to `mspfileloaderv14` (still not defining the name); it now reads `from mspvalidator_merger import append_feature_space_decoys`. Test: `tests/test_hygiene_cli_import.py`.

### Description

The `--add-decoys` branch imports `append_feature_space_decoys` from `mspfileloaderv10` (a legacy module):

```python
from mspfileloaderv10 import append_feature_space_decoys
```

`mspfileloaderv10` may not exist in the current working directory. The correct import should be:

```python
from mspvalidator_merger import append_feature_space_decoys
```

`mspvalidator_merger.append_feature_space_decoys` exists at line 679.

---

## Issue 9: tau Default Inconsistency

**Modules:** `ml_ng_utils.py`, `prediction_report.py`, `msp_CLI_predict_skops.py`
**Severity:** Low — but callers must specify tau explicitly

### Description

The confidence threshold tau has different defaults in different contexts:

| Location | Default tau | Notes |
|----------|------------|-------|
| `predict_with_threshold()` in `ml_ng_utils.py` line 178 | 0.55 | Applied during inference |
| `ReportParams.tau` in `prediction_report.py` line 35 | 0.60 | Applied during acceptance filtering |
| Coverage display in `msp_CLI_predict_skops.py` line 76 | 0.65 | Informational only, not applied to output |

These are intentionally different use-cases (inference floor vs report acceptance threshold vs coverage metric), but callers must specify tau explicitly rather than relying on defaults for consistent behavior.

---

## Issue 10: PPM Asymmetry in peaks_ppm

**Module:** `mspvalidator_merger.py`
**Line:** ~353–400
**Severity:** Low — minor systematic error, consistent across all comparisons

### Description

The PPM error calculation in `peaks_ppm()` uses:

```python
denominator = max(observed_peak, reference_ion)
```

The symmetric formula would use:
```python
denominator = reference_ion
```

The asymmetric formula causes slightly narrower effective tolerance when the observed peak is higher than the reference. Since this applies uniformly to all comparisons, it affects scores consistently rather than selectively. See [CGA_SCORING_ANALYSIS.md](../CGA_SCORING_ANALYSIS.md) Section 7 for full analysis.

---

## Issue 11: msp_insilicomarker_withGPT.py Duplicate extract_ionmasslist

**Module:** `msp_insilicomarker_withGPT.py`
**Lines:** 27 (first), 1278 (second)
**Severity:** RESOLVED 2026-10-01 (hygiene H3) — first definition removed together with the hard-coded U937 `if __name__ == "__main__":` dev block (its only module-load-time caller); `extract_ionmasslist(ionmass_sheet, mass_col="mass")` is the single definition. Test: `tests/test_hygiene_shadowed_defs.py`.

### Description

`extract_ionmasslist` is defined at line 27 (no `mass_col` parameter) and again at line 1278 (with `mass_col="mass"` parameter). The line-1278 version is active at runtime.

---

## Issue 12: mspfileloaderv14.py CGA_glycotope_score Dead Body

**Module:** `mspfileloaderv14.py`
**Line:** 2163
**Severity:** Low — function exists but body is dead code

### Description

`CGA_glycotope_score()` is defined at line 2163 as a new structural scoring stub. The function sets up local state (lines 2168–2171), but the entire scoring loop body is wrapped in a triple-quoted string (lines 2174–2196). The function returns `None` implicitly.

This appears to be planned work that was superseded by the `msp_CGA_structscore.py` integration. The function is not called anywhere in the active codebase.

---

## Issue 13: mspfileloaderv14.py Duplicate os/traceback Imports

**Module:** `mspfileloaderv14.py`
**Lines:** 1, 23, 27 (`os`); 17, 27, 34 (`traceback`)
**Severity:** Lint — harmless duplicates

### Description

`os` is imported three times at module top level. `traceback` is imported three times. These are harmless but messy.

---

## Issue 14: score_b_skeleton_v3.py run_candidate_case Missing Return

**Module:** `score_b_skeleton_v3.py`
**Line:** 1868
**Severity:** Low — type annotation bug in experimental module

### Description

`run_candidate_case()` has return type `ScoreBResult` but no `return` statement — returns `None`. Callers expecting a return value will get `None`.

---

## Issue 15: MAS-ML Feature Writer Keeps Last In-Window Peak (B-38)

**Module:** `mspvalidator_merger.py`
**Function:** `createnormailzedionlistcsv` (`data_dict` loop after `findingions(...)`). Grep to confirm the line.
**Severity:** Very low. Latent; 0 occurrences on the published zebrafish MAS-ML fixtures. Filed 2026-10-12.

### Description

`peaks_ppm()` returns one entry per in-window peak. The writer's `data_dict` loop overwrites earlier entries, so when two or more peaks fall within one target's ppm window, the MAS-ML feature is `log10(I)+1` of the **last** (highest-m/z) peak, not the most intense one. CGA-ML training (`ml_ng_utils.build_features_from_peaks_log10_plus1`) and prediction (`extract_ion_intensities`) both use the max. Tracker: `handoff/B38_mas_peaks_ppm_last_peak_wins_stub.md`.

---

## Issue 16: RAW Conversion Gave Empty Peak Lists for Ion-Trap Scans (READER-DIAG)

**Module:** `msprawextractor.py` (`convert_raw_to_csv`; now `read_scan_peaks`). Grep to confirm the line.
**Severity:** High for ion-trap (LTQ, ITMS) RAW files; none for FTMS data. Filed 2026-10-13.
**Status:** **Diagnosed and fixed on branch `reader-diag` (CP-C.2), not yet merged into `main`.** `main` still has the defect.

### Description

The RAW path read every MS2 scan's peaks with `GetLabelData`. Label data exist only for FT scans; for ion-trap scans the call succeeds and returns 0 points, so an external LTQ ion-trap RAW file (the pilot RAW) converted to rows with empty peak lists. The read-only probe (`msp_extract_cli inspect`, Lab PC 2026-10-13) showed the mass list holds the points and that the vendor-centroided mass list is bit-identical to MSConvert's vendor peak picking. Fix: FTMS scans keep label data; other scans read the mass list (vendor-centroided when the scan is profile). MS3 rows now read their own peaks (they repeated the previous MS2's arrays). Lab PC re-runs 2026-10-16 and 2026-10-20: peaks equal to MSConvert's vendor peak picking. Record: `handoff/READER_DIAG_handoff_20261013.md` ("Unresolved RAW items", CP-C.2).

---

## Issue 17: Rows With Unknown Charge Reach Downstream Steps as Mass 0 (READER-DIAG, v14 follow-up)

**Modules:** `mspfileloaderv14.py` (CGA precursor match, CGA→Train precursor gate, mass feature, training, prediction, unlabeled build), `msp_CLI_predict_skops.py`. Grep to confirm the lines.
**Severity:** Medium for files without reported charges; none for files whose charges are all reported. Filed 2026-10-13 (Codex READER-DIAG F1/B8).
**Status:** **Closed by decision (V14-FOLLOWUP, 2026-10-27; branch `v14-reader-followup`, not yet merged).** D-1: no unknown-charge gate and no per-row guards; downstream behaviour is deliberately unchanged (the rows still reach downstream steps as mass 0), and any gate would be a separate, re-run-and-re-reported change. What changed: the GUI now shows the reader's warning with the count of MS2 without a charge (log pane, orange status), and the comb / dataset charge are available in the GUI ("Conversion options" dialog, both off by default), so fewer rows are affected when they are used. Converter side (READER-DIAG): unknown charges stay `chargeState 0`, `protonatedmass 0` (CSV unchanged), and are marked in `<stem>.scan_status.tsv` (`charge_status unknown`, `mass_status unavailable`) and counted in the diagnostics (exit 3).

**Why closed (Henry, 2026-10-29):** the issue is left unfixed on purpose; the effect is assumed negligible and not to affect the contents of the datasets currently expected. Measured at the V14-FOLLOWUP manual acceptance (2026-10-30): on the LTQ native files 20–29 % of exported MS2 rows have unknown charge, a bounded loss that is not the cause of the near-zero CGA result on those files (that is Issue 27).

### Description

No consumer checks `mass > 0` or `charge != 0`: CGA keeps such a scan as "unmatched" with no reason, and a mass-trained model sees mass 0 as a Non-glycan signature. Required behaviour per site: `handoff/READER_DIAG_handoff_20261013.md` B8; list: `handoff/READER_DIAG_deferred_items_20261013.md`. Henry, 2026-10-15: v14 `mass > 0` guards are not planned.

---

## Issue 18: GUI Conversion Leaves the Diagnostics Files Next to the Input File (READER-DIAG, v14 follow-up)

**Module:** `mspfileloaderv14.py` (`background_conversion` → finalize at save). Grep to confirm the lines.
**Severity:** Low (findability). Filed 2026-10-20 (Lab PC CP-D run).
**Status:** **Fixed on branch `v14-reader-followup` (V14-FOLLOWUP CP-B, 2026-10-27), not yet merged.** At save, v14 moves `.reader_diag.json/.txt` and `.scan_status.tsv` with the CSVs to the chosen folder under the save name (`<savename>.reader_diag.json` …) and rewrites their output references to the final names; the outcome and warnings are shown in the log pane, with a dialog for failed or partial results. The diagnostics are used only when they belong to this conversion (D-4). Later files of a batch are already converted into the chosen folder.

### Description

v14 converts before it asks for the output folder, so the readers write the temporary CSVs and the three CP-D files (`<input stem>.reader_diag.json/.txt`, `.scan_status.tsv`) into the RAW/mzML file's folder. At save, v14 promotes only the two CSVs to the chosen folder. The GUI also does not show the outcome or warnings yet.

---

## Issue 19: Profile MS2 in an mzML Is Written Point by Point (READER-DIAG)

**Module:** `mzmlreader.py`. **Severity:** Medium when it happens (rows hold hundreds of profile points, zeros included, per scan). Filed 2026-10-20.
**Status:** Warned, not changed. The reader does not centroid (no new centroiding path); since CP-D the diagnostics warn ("profile data") and the run exits 3.

### Description

Seen when MSConvert ran `zeroSamples` **before** `peakPicking vendor`: MSConvert then skips peak picking entirely (it prints "peakPicking is not the first filter"). Re-convert with `peakPicking vendor msLevel=2-` as the first filter. With that, the mzML route reproduces the RAW route's rows, peaks and charges (Lab PC 2026-10-20).

---

## Issue 20: MS1 Isotope Comb Limits (READER-DIAG, by design)

**Module:** `msp_charge_comb.py`. **Severity:** informational. Filed 2026-10-20.
**Status:** By design; stated in the user guide (`handoff/READER_DIAG_charge_supply_summary_20261015.md` §6, §6a).

### Description

The comb tests z = 1–3 on profile, non-FT MS1 only; a mixture of charge states is reported as the stronger one without a flag (3 of 8 hand-checked z = 3 calls were mixtures); z ≥ 4 is never called; 24–27% of MS2 scans in the LTQ files tried stayed `below_margin`. It is a feasibility-checked heuristic, not a validation; calls are checked downstream by scoring.

---

## Issue 21: Consumers Must Check Array Integrity and the Reader Outcome (READER-DIAG R2, v14 follow-up)

**Modules:** `mspfileloaderv14.py` (`background_conversion` → `on_conversion_complete`; scoring that indexes `peakintensity`), `mspvalidator_merger.py`. Grep to confirm the lines.
**Severity:** Low–Medium (malformed input only). Filed 2026-10-23 (Codex R2-F15).
**Status:** **Partly done (V14-FOLLOWUP, 2026-10-27; branch `v14-reader-followup`, not yet merged).** The reader outcome part is done: a `failed` result stops promotion in the GUI (CP-B, D-2). The array-integrity part is **unchanged and open**: WP-C (the check in `mspvalidator_merger.expandpeaklist`, D-3) was not triggered, because no CSV written by the current readers (pilot RAW, U-937, zebrafish, LTQ B/C; RAW and mzML) has a `()` row or an m/z-vs-intensity length mismatch; only two CSVs written before READER-DIAG CP-C.3 crash (`handoff/V14_FOLLOWUP_wpc_check_20261027.md`). Readers unchanged by design (they write what the input holds).

### Description

A row whose m/z and intensity arrays differ in length is exported with `arrays_status length_mismatch` in the status file; a consumer that matches the second m/z and reads `peakintensity[1]` raises `IndexError`. Consumers should check array integrity before scoring. Where native IDs repeat (several controllers), sidecar joins should use `native_id`, and charge / mass-source provenance should be kept in derived files. A reader result with `outcome: failed` must stop promotion in the GUI, not only be shown (Issue 18).

---

## Issue 22: Route B (Manual Annotation) Uses PerMe Masses for Native Glycans (deferred until a usable dataset exists)

**Modules:** `mspvalidator_merger.py` (`directassign_files` → `adding_masses_mode_aware` → `calc_neutral_mass`; the Native branch of `import_sugarinfo` only prints), called from `mspfileloaderv14.py` (Route B, `Derivatization Type` passed as `derivatization_type`). Grep to confirm the lines.
**Severity:** Medium for a Native project annotated through Route B (the theoretical masses are wrong); none for PerMe. Filed 2026-10-27 (V14-FOLLOWUP CP-A trace).
**Status:** Deferred by design (Henry, 2026-10-27).

### Description

`calc_neutral_mass` always sums the permethylated residue table (`monosaccharide_masses`), whatever the derivatization, so a Native or Native(Reduced) project going through Route B gets PerMe theoretical masses. V14-FOLLOWUP CP-A gave the CGA route (Route A, `compnewv4.precursormassv4`) native and native-reduced masses; Route B was left unchanged.

**Disclaimer:** native glycans are currently supported through the CGA (constraint-based) route only. The manual-annotation route for native glycans is deferred: the patch lands when a usable, expert-annotated native dataset is available, so that the change is validated on real annotations rather than on annotations made for the purpose.

---

## Issue 23: A Batch Ends With "All metadata have been completed" After a Failed Finalize (V14-FOLLOWUP, pre-existing)

**Module:** `mspfileloaderv14.py` (`MetadataEditorWindow.generate`: CSV promotion → batch advance). Grep to confirm the lines.
**Severity:** Low (the log has the `[WARNING] MS2/MS3 finalize failed` line; the message is misleading). Filed 2026-10-27 (Codex V14F-F5).
**Status:** Open; deferred (Henry, 2026-10-27). Predates the V14-FOLLOWUP branch.

### Description

When `promote_temp_to_final` fails for a file, `generate` sets `success = False`, skips the log save and `finalize_extraction`, and still advances the batch; the last file then shows "All metadata have been completed." even though one file's CSVs were not finalized. The V14-FOLLOWUP branch removed the most common macOS/Linux cause (EXDEV across volumes now falls back to copy-then-delete) but did not change the message. A fix would report the failed files in the final message (a GUI text change, so a separate decision).

---

## Issue 24: Two Simultaneous Conversions of the Same Stem Into One Folder Can Mix Their Diagnostics (V14-FOLLOWUP, by design)

**Module:** `mspfileloaderv14.py` (`reader_diag_fingerprint`, `load_reader_diagnostics`, `summarize_reader_result`). Grep to confirm the lines.
**Severity:** Very low (needs two GUI windows, or the GUI and `msp_extract_cli`, converting inputs with the same stem into the same folder at the same moment). Filed 2026-10-27 (Codex V14F-F2 remainder).
**Status:** Accepted for now (Henry, 2026-10-27); a run-owned staging folder is a possible future optimisation.

### Description

v14 decides which `<stem>.reader_diag.json` belongs to a conversion from the reader's return (a current reader names its files), the pre-call fingerprint (a file unchanged by the call, or one whose state could not be read, is not used), the timestamp and the input identity. Another process writing the same path during the call changes the file, so its diagnostics could be taken as this run's. Only the shown explanation and the moved diagnostics are affected: whether the CSVs are promoted comes from the reader's own return. A per-run temporary output folder would remove the case, at the cost of changing where temporary files live during conversion.

---

## Issue 25: Conversion Log Lines Are Written to the Tk Text Widget From the Worker Thread (pre-existing pattern)

**Module:** `mspfileloaderv14.py` (`AppLogger.log` → `write_to_gui`; called from `background_conversion`). Grep to confirm the lines.
**Severity:** Low (timing-dependent; not reproduced). Filed 2026-10-27 (Codex V14F-F6).
**Status:** Open; recorded, no change (Henry, 2026-10-27). A logger queue drained on the Tk thread would be a separate ticket.

### Description

`AppLogger.log` calls the GUI writer directly, which inserts into the Tk text widget. `background_conversion` runs in a worker thread and logs from it on every conversion (the temp-file checks, the `[FATAL]` line, and since V14-FOLLOWUP the `[reader] …` lines on the failure path). Dialogs and status updates are marshalled with `_safe_after`; log lines are not. If the window closes at the wrong moment, a Tk error could interrupt the failure recovery.

---

## Issue 26: PerMe Masses for Mixed SO3 + PO3H, or Three or More of Either, Are Not Modelled (compnewv4, v1.10 behaviour)

**Module:** `compnewv4.py` (`precursormassv4`, PerMe branch: the `so3 == 1/2`, `po3h == 1/2` rules). Grep to confirm the lines.
**Severity:** None for the real analyses (these compositions are not observed); the values are wrong for them. Filed 2026-10-27 (V14-FOLLOWUP, Codex round 2 check; Henry).
**Status:** Open; documented, no change to the frozen arithmetic. `msp_CLI_precursor_mass` 0.3 refuses these compositions (exit 2); the CGA path still computes them as in v1.10.

### Description

The PerMe sulfate/phosphate patch was fitted to observed species, one group type at a time. Henry's charge rule (2026-10-27): each SO3 / PO3H carries one charge; two PO3H is the observed species that loses one charge by elimination, so it is singly charged. The code models SO3 1 (z = 1), SO3 2 (halved, z = 2), PO3H 1 and PO3H 2 (z = 1). For mixed compositions it has no combined charge state: s1p1 (charge 2), s2p1 and s1p2 (charge 3) give singly or doubly charged values, and for s2p2 the PO3H 2 correction is applied after the halving (a full-mass correction on a half-mass value). SO3 3 or more gets the flat SO3 only, and PO3H 3 or more gets no phosphate term at all. Native masses are not affected (flat additions, always singly charged).

---

## Issue 27: CGA Mass Tolerances Are Fixed in ppm, So Ion-Trap Data Get Almost No Annotations (v1.10 design)

**Modules:** `mspfileloaderv14.py` (`CGASetupWindow` → `on_start` → `run_pseudolabeling`; CGA→Train gate and eligibility; features and prediction), `msp_insilicomarker_withGPT.py` (`find_compositions_for_mass`, `attach_ion_score_on_matched`), `mspvalidator_merger.py` (`findingions` / `peaks_ppm`), `ml_ng_utils.py` (`build_features_from_peaks_log10_plus1`). Grep to confirm the lines.
**Severity:** High for ion-trap (LTQ) data whose mass is the isolation m/z; none for Orbitrap-class data with a monoisotopic m/z. Filed 2026-10-30 (V14-FOLLOWUP manual acceptance step 5).
**Status:** **Fixed by option on branch `cga-tolerance` (GlycoMSP v14 1.14, CGA-TOLERANCE, 2026-11-01), off by default, not yet merged.** `main` still has the fixed tolerances. CGA Setup has a *Mass tolerances (CGA matching)* frame: precursor match (T1) and Score A fragment match (T2), each in ppm or Da (`_parse_tolerance` → `_cga_tolerance_kwargs` → `run_pseudolabeling` `tol_da` / `ion_tol_da` → `find_compositions_for_mass` `tol_da`, `attach_ion_score_on_matched` / `findingions` / `peaks_ppm` `tol_da` → `_peaks_da`). The CGA → Trainable precursor gate, which had no widget since v1.10 (on, 10 ppm, theoretical-mass relative), now has a checkbox and ppm entry (CP-C.4, `_trainable_gate_ppm`). The CGA → Trainable fields, the negatives and the prediction feature window stay ppm-only. Tolerances are recorded (T3: method JSON `parameters.cga`, `training_run.json` `cga_tolerances`, `run_summary.json` `model_tolerances`) and non-default values prefill prediction. Defaults unchanged: the default CGA TSV is byte-identical to `7910068` on the brain fixture. Not done: T4 (deferred); Score B enrichment keeps its fixed 20 ppm. The real-data check on the ion-trap files is Henry's manual test. Ticket `handoff/TICKET_cga_tolerance_iontrap_20261030.md`; reviews `handoff/codex_CGA_TOLERANCE_20261031*` (rounds 1–5).

### Description

The CGA precursor match uses a fixed 20 ppm (`run_pseudolabeling` default; `on_start` passes no tolerance and CGA Setup has no field for it) and Score A matches fragments at a fixed 10 ppm. On LTQ native files (mono m/z 0, so the mass is the isolation m/z) the precursor error is −20 … −90 ppm at z = 1 and +150 … +250 ppm at z = 2, and fragment errors are 0.06–0.14 Da. 20 ppm kept about 5 % of the real precursor matches, and Score A at 10 ppm removed nearly all of those, so almost no trainable CSV could be built; about 300 ppm / 0.3 Da would give tens of annotations per file against a low decoy rate. The error changes sign with charge, so a wider window (ppm or Da) is needed, not an offset. The user manual is wrong on this point: it gives a 10 ppm default and says the tolerance is configurable in CGA Setup (corrected in the manual on 2026-11-01: §2.4 (d), §2.5, §5.5). The offset was recorded during READER-DIAG (`handoff/READER_DIAG_comb_design_20261015.md` §13.4), but the CGA precursor tolerance was assumed to be exposed.

Found during the fix (2026-11-01): CGA → Trainable also applied a second, hidden precursor filter, `build_trainable_from_CGA(precursor_gate_ppm=10)`, |1e6 · (protonatedmass − theoretical) / theoretical| ≤ 10, with no widget. It is stricter than the visible *max |ppm_error|* (20), so in v1.10 the effective precursor cut at CGA → Trainable was 10 ppm, although CGA analysis matched at 20 ppm. The v1.10 manual (§2.6) did not mention it.

---

## Issue 28: Score A Counts Peak–Ion Pairs, So One Ion or One Peak Can Count More Than Once (B-39, v1.10 behaviour)

**Modules:** `mspvalidator_merger.py` (`peaks_ppm`, `_peaks_da`), `msp_insilicomarker_withGPT.py` (`attach_ion_score_on_matched`), `mspfileloaderv14.py` (`_fallback_simple_ion_scoring`; `score_counter`). Grep to confirm the lines.
**Severity:** Low at the default 10 ppm; Medium at wide fragment windows (Score A and `ion hit count` rise without added evidence). Filed 2026-10-31 (CGA-TOLERANCE CP-B).
**Status:** Open; kept on purpose (Henry, 2026-10-31): GlycoMSP 1.14 keeps the v1.10 counting rule in ppm and in Da, so existing results stay reproducible. A fix changes Score A and must ship as its own option, or be re-run and re-reported; C-001 and C-003 must then be re-checked. Tracker `handoff/B39_scoreA_counts_each_in_window_peak_stub.md`. Explained to users in the user manual §4.1 (*How hits are counted*).

### Description

For each ion, every peak inside the ion's window appends one hit; nothing marks a peak or an ion as used. Case 1: several peaks in one ion's window count that ion several times (`ion hits m/z` repeats it). Case 2: one peak inside the windows of two ions closer than twice the window counts for both. Brain fixture (zebrafish PerMe, Orbitrap, 282 ions): case 1 0 % of matched rows at 20 ppm / 10 ppm, 31.8 % at 300 ppm / 0.3 Da; case 2 in 323 scans at 10 ppm, 1,285 at 0.3 Da, 8,844 at 0.5 Da (of 23,193). Same loop as Issue 15 (B-38).

---

## Issue 30: The Prediction Command Line Cannot Load a Model With the Pinned skops (TESTS-PUBLIC, v1.10 code)

**Modules:** `msp_CLI_predict_skops.py`, `msp_CLI_joblibtoskops_compat.py` (both call `skops.io.load(..., trusted=True)`). Grep to confirm the lines. Filed 2026-11-18 (TESTS-PUBLIC CP-B).
**Severity:** Medium: the command-line prediction path does not run at all with the pinned library; the GUI is not affected.
**Status:** Open. Not fixed before the manuscript submission (Henry, 2026-11-18); planned for v1.15. Tripwire: `tests/test_ml_features_train_predict.py::test_cli_loads_a_model_saved_with_the_pinned_skops` is a strict expected failure and will fail loudly once the CLI is fixed, so that the marker is removed on purpose.
**Workaround:** predict in the GUI (Predict tab), which loads the file with the types it lists; or install skops older than 0.10 for the command line (`pip install "skops<0.10"`).

### Description

skops 0.10 and later reject `trusted=True` with `TypeError: trusted must be a list of strings`. `src/requirements.txt` pins skops 0.13.0, so a clean install following the documented requirements cannot run `msp_CLI_predict_skops.py` or the compatibility converter. The GUI works because it passes the untrusted types reported by `skops.io.get_untrusted_types` for the file.

---

## Issue 31: The Prediction Command Line Drops the Mass Feature by Default (TESTS-PUBLIC, v1.10 code)

**Modules:** `msp_CLI_predict_skops.py` (`--drop-cols`, default `MS2scan_no protonatedmass`). Filed 2026-11-18 (TESTS-PUBLIC CP-B).
**Severity:** Medium for a model trained with the mass feature (which production training requires).
**Status:** Open. Not fixed before submission; planned for v1.15. Recorded as a CHARACTERISATION test (what the command line does today, not a claim that it is desired).
**Workaround:** predict in the GUI; or pass `--drop-cols` with an explicit list that does not name the mass column the model was trained with.

### Description

The default drop list removes `protonatedmass`. A model trained with the mass feature then misses that column, and the missing column is filled with zero, which is also the value of the mass feature for Non-glycan rows, so the predictions are made on a changed input without any message.

---

## Issue 32: The Prediction Command Line Writes Integer Class Labels (TESTS-PUBLIC, v1.10 code)

**Modules:** `msp_CLI_predict_skops.py`; GUI training (`_labelencoder.joblib` sidecar). Filed 2026-11-18 (TESTS-PUBLIC CP-B).
**Severity:** Low to Medium: output is not human-readable and can be misread.
**Status:** Open. Not fixed before submission; planned for v1.15. Recorded as a CHARACTERISATION test.
**Workaround:** predict in the GUI; or map the integers back with the saved `_labelencoder.joblib` (`LabelEncoder.inverse_transform`).

### Description

The GUI trains the Random Forest on LabelEncoder integers and saves the encoder next to the model. The command line never reads the sidecar, so its output column holds the integer codes instead of the glycan composition labels.

---

## Issue 33: A Zero or Negative Peak Intensity Stops Feature Building (TESTS-PUBLIC, v1.10 code)

**Modules:** `ml_ng_utils.py` (`build_features_from_peaks_log10_plus1`). Filed 2026-11-18 (TESTS-PUBLIC CP-B).
**Severity:** Low as far as known: a RAW or mzML conversion writes the intensities of detected peaks; whether any dataset carries a zero or negative value has not been checked.
**Status:** Open. Not fixed before submission; planned for v1.15.
**Workaround:** remove peaks with intensity ≤ 0 from the converted CSV before training or prediction.

### Description

The guard `if I < 0: I = 0.0` is followed by `math.log10(I)`, so a matched peak with intensity 0 or a negative one raises `ValueError: math domain error` and the whole run stops; the guard cannot take effect. A miss correctly gives the baseline 1.0.

---

## Issue 34: Four Different Default Confidence Thresholds (tau) in the Code (TESTS-PUBLIC, v1.10 code)

**Modules:** `mspfileloaderv14.py` (training default 0.65 with margin 0.05; Predict-tab fallback 0.60), `prediction_report.py` (`ReportParams` 0.60), `ml_ng_utils.py` (`predict_with_threshold` 0.55). Grep to confirm the lines. Filed 2026-11-18 (TESTS-PUBLIC CP-B).
**Severity:** Low: the rule that is applied is the margin-aware demotion in v14 (`demote = (y_pred != majority) & (maxp < tau) & (p_major >= maxp - margin)`); the other values are defaults of paths that are not reached or are overridden.
**Status:** Open. Not changed before submission; planned for v1.15 (one documented default). Both copies of the v14 demotion expression are tested by reading the source text.
**Workaround:** set tau and margin explicitly in the GUI and record them (they are written to the run summary).

### Description

`CLAUDE.md` and the GUI state tau = 0.65. The Predict tab falls back to 0.60 when the entry cannot be read, the prediction report object defaults to 0.60, and `predict_with_threshold` defaults to 0.55 but is only called by `tau_sweep_summary`, which the GUI does not call. A script that uses one of the library defaults gets a different cut-off than the GUI.

---

## Resolved Issues (carried forward for record)

| Issue | Resolution date | Notes |
|-------|----------------|-------|
| log10 formula wrong in marker and v13 | 20260122 | Fixed in all active code; misleading comment remains at marker L396 |
| ARCHITECTURE.md v13 line numbers all off by 300–600 | 2026-03-18 | Corrected in this wiki; ARCHITECTURE.md preserved as-is |
| `make_trainable_dataframe` documented as existing in ml_ng_utils | 2026-03-18 | Confirmed does not exist; removed from wiki |
| `ANCHOR_IONS` / `anchor_required` documented as existing | 2026-03-18 | Actual names: `ANCHORS` (list), no `anchor_required` parameter |
| `score_counter` line 1241 tagged [DEAD] | 2026-03-18 | Confirmed active via Python late binding |
| Score B "not imported by any active module" | 2026-03-31 | Score B now integrated via `msp_CGA_structscore.py`, imported by `mspfileloaderv14.py` |
