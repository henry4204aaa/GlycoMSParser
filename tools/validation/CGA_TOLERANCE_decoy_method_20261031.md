# CGA-TOLERANCE — decoy count method (chance-match baseline for widened CGA tolerances)

**Provenance:** method introduced on 2026-10-30 in the V14-FOLLOWUP step-5 investigation (session `945bdee9-3428-43c1-8170-9565deba5a67`, read-only; ticket `handoff/TICKET_cga_tolerance_iontrap_20261030.md`); adopted as the acceptance criterion by Henry on 2026-10-31 (kickoff addendum 8, option b); implemented and tested in the CGA-TOLERANCE session `73787ea0-bc0c-4655-a304-41288139a488` on branch `cga-tolerance` (GlycoMSP v14 1.14).
**Script:** `handoff/cga_tolerance_decoy_count.py` (standalone; Python 3.9+, pandas, numpy; openpyxl only for an `.xlsx` ion list). **Tests:** `handoff/test_cga_tolerance_decoy_count.py` (synthetic, 6 tests). Not in `src/`, not part of the application.
**Naming rule:** outputs carry only the `--alias` given on the command line (e.g. `pilot`, `run1`); no input file name, path or scan list is printed or written. Error messages show only the exception type for anything that could carry a path.

## Question it answers

When a CGA tolerance is widened (precursor 20 ppm → e.g. 300 ppm; Score A fragments 10 ppm → e.g. 0.3 Da), how many of the extra matches are chance?

## Method

1. **Decoy library.** Every in silico library mass is shifted by **+3.37 Da**, an arbitrary offset chosen to avoid residue, isotope (~1.003 Da) and charge spacings. No decoy mass is a plausible glycan of the same sample, so every decoy match is chance.
2. **Same scans, same tolerance.** The converted MS2 CSV is matched against the real and the decoy library with the same tolerance. Count per tolerance: **real scans** = scans with ≥ 1 real candidate, **decoy scans** = scans with ≥ 1 decoy candidate (and candidate pairs, capped at the CGA top-3 per scan). The decoy count estimates the chance matches among the real ones; `decoy / real` is the chance fraction.
3. **Fragments (optional).** The ion list is shifted by the same +3.37 Da. Score A is computed per scan against the real and the decoy ion list. A scan passes when Score A ≥ 0.07 (`--min-score`) or, to compare ion lists of different size at equal hit counts, hits ≥ N (`--min-hits`). Four counts per (precursor, fragment) pair: real precursor + real ions (what CGA would keep), decoy precursor + decoy ions (full chance baseline), and the two mixed combinations (which gate carries the specificity).
4. **Several tolerances in one run** (`--precursor 20ppm,300ppm,0.5Da --fragment 10ppm,0.3Da`); every precursor × fragment pair is reported.

## Matching rules (identical to GlycoMSP 1.14)

| Step | Rule |
|---|---|
| precursor, ppm | \|obs − theo\| ≤ obs · ppm / 10⁶ (window centred on the observed `protonatedmass`, inclusive) |
| precursor, Da | \|obs − theo\| ≤ Da |
| scans counted | rows with `protonatedmass` > 0 (charge-unknown rows have mass 0 and never match; known issue 17) |
| fragment, ppm | peak ≤ ion: `1000000 * (ion − peak) / ion <= ppm`; peak > ion: `1000000 * (peak − ion) / peak < ppm` (the v1.10 asymmetry, known issue 10); the scan of an ion stops at the first peak above it that is outside the window |
| fragment, Da | `ion − peak <= Da` / `peak − ion <= Da` (same stop rule) |
| fragment hit | `log10(I) + 1 > 1` with I < 0 taken as 0, computed as GlycoMSP computes it; a repeated m/z takes the intensity of its first occurrence |
| exactness | every final hit decision uses GlycoMSP's float operations above; binary search only finds candidates (with a 1e-6 Da margin). Binary search is used only where it is provably exact (ppm ≤ 500,000, Da ≤ 1,000, masses finite and ≤ 10⁶ Da, and for a ppm window every ion mass > 0); everything else, and every unsorted or repeated-m/z peak list (count reported as "irregular peak lists"), runs the literal GlycoMSP loop, so the counts equal GlycoMSP's for every finite input. The script reproduces GlycoMSP, it does not check the inputs make physical sense: a zero or negative ion mass gives GlycoMSP's own result (or its division-by-zero error), not a rejection. Tested against `mspvalidator_merger.peaks_ppm` / `_peaks_da` on random spectra, on peaks at the window edges and one float either side, at extreme windows, and on mixed-sign masses (Codex rounds 1–4) |
| Score A | hits / ion-list size, capped at 1.0, per scan (independent of the candidate composition) |
| counting | `glycomsp` = one hit per in-window peak (what GlycoMSP writes, B-39); `unique` = each ion at most once. Neither removes B-39 case 2 (one peak counted for two close ions); `handoff/B39_shared_peaks.py` measures that |

**Checked against GlycoMSP itself** (Mac, 2026-10-31, in-house fixture `tests/fixtures/smokev1.20/1c_brain/CGA`, library `…011456`, `zf_ionlist.csv` 282 ions): the script's real counts equal the CGA TSVs written by `run_pseudolabeling` on the same inputs: scans matched 4,782 (20 ppm) and 8,767 (300 ppm); candidate pairs 8,917 and 17,431; scans with Score A ≥ 0.07 1,571 (20 ppm / 10 ppm) and 2,844 (300 ppm / 0.3 Da). The same run gives 15 decoy scans at 20 ppm but **11,664 decoy vs 8,767 real at 300 ppm** on this Orbitrap PerMe data (dense library): a wide window is a per-dataset choice, which is what this check is for.

## Difference from the 2026-10-30 investigation numbers

The ticket's reference numbers (ten native LTQ files, precursor only: 20 ppm 65 real / 8 decoy, 300 ppm 1,315 / 57, 0.5 Da 1,192 / 47) came from a scratch script with the same decoy library and scan definition (rows with mass > 0, nearest library mass within the window) but: ppm relative to the **theoretical** mass instead of the observed mass (negligible at these errors); fragment Score A with **each ion once** (`unique`), symmetric ppm, and **no intensity filter**. Expect the script's precursor counts to match those numbers closely; fragment counts can differ slightly and should be compared with the `unique` row.

## Run (Lab PC, per file)

```
python cga_tolerance_decoy_count.py --alias pilot --converted <converted MS2 CSV> --library <in silico CSV> ^
    --ions <ion list CSV/XLSX> --precursor 20ppm,300ppm,0.5Da --fragment 10ppm,0.3Da --out pilot_decoy.csv
```

Repeat with the 97-ion candidate list (`--ions …`, and `--min-hits 1`/`2` to compare at equal hit counts: one hit = 0.071 with 14 ions, 0.010 with 97). Only the alias-labelled table / CSV comes back from the Lab PC.

## Reading the result

Use the decoy count as the chance estimate at the same tolerance: real − decoy ≈ matches above chance; `decoy / real` ≈ fraction expected to be chance. Prefer the tolerance where real − decoy stops growing while decoy keeps growing. The script reports counts only; it does not choose a tolerance and does not touch any GlycoMSP output.
