# tools/validation

Controlled validation scripts for GlycoMSP. They are driven by a person or an AI assistant on purpose-chosen inputs, print counts or timings, and answer one question each. **They are not application code**: nothing in `src/` imports them, they change no GlycoMSP output, and they are not part of the GUI or the command-line tools. They need Python 3.9+, numpy and pandas (`openpyxl` only for an `.xlsx` ion list); the scripts that reuse product logic import it from `src/`.

Every script reports under a user-chosen `--alias` and never prints or writes an input file name, path or scan list.

| Script | Question it answers | Tests |
|---|---|---|
| `rf_timing.py` | How long does the GUI's Random Forest train step take on a given trainable CSV (fit alone, and CSV load to fitted model)? | `tests/test_rf_timing.py` (synthetic 30-row trainable) |
| `cga_tolerance_decoy_count.py` | When a CGA mass tolerance is widened, how many of the extra matches are chance? (decoy-library baseline, below) | `tests/test_cga_tolerance_decoy_count.py` (31 tests, synthetic) |
| `B39_dup_hits.py` | How many Score A rows list the same ion more than once (several peaks inside one ion's window)? | — (counts a CGA TSV you supply) |
| `B39_shared_peaks.py` | Which ion pairs are close enough to share one peak, and how many peaks does GlycoMSP count as a hit for two or more ions? | `tests/test_cga_tolerance_decoy_count.py` (shared matching code) |

The tests of these scripts live in `tests/` and run with the rest of the suite (`python3 -m pytest tests`); the scripts stay here because they are tools, not tests.

## Decoy count method (`cga_tolerance_decoy_count.py`)

**Question.** When a CGA tolerance is widened (precursor 20 ppm to, say, 300 ppm; Score A fragments 10 ppm to, say, 0.3 Da), how many of the extra matches are chance?

1. **Decoy library.** Every in silico library mass is shifted by **+3.37 Da**, an arbitrary offset chosen to avoid residue, isotope (about 1.003 Da) and charge spacings. No decoy mass is a plausible glycan of the same sample, so every decoy match is chance.
2. **Same scans, same tolerance.** The converted MS2 CSV is matched against the real and the decoy library with the same tolerance. Per tolerance the script counts *real scans* (at least one real candidate), *decoy scans* (at least one decoy candidate) and candidate pairs (capped at the CGA top 3 per scan). The decoy count estimates the chance matches among the real ones; `decoy / real` is the chance fraction.
3. **Fragments (optional).** The ion list is shifted by the same +3.37 Da and Score A is computed against the real and the decoy ion list. A scan passes when Score A >= 0.07 (`--min-score`) or, to compare ion lists of different size at equal hit counts, hits >= N (`--min-hits`). Four counts per (precursor, fragment) pair: real precursor + real ions (what CGA keeps), decoy + decoy (full chance baseline), and the two mixed combinations.
4. **Several tolerances in one run:** `--precursor 20ppm,300ppm,0.5Da --fragment 10ppm,0.3Da`; every precursor x fragment pair is reported.

**Matching rules (identical to GlycoMSP 1.14):**

| Step | Rule |
|---|---|
| precursor, ppm | `abs(obs - theo) <= obs * ppm / 1e6` (window centred on the observed `protonatedmass`, inclusive) |
| precursor, Da | `abs(obs - theo) <= Da` |
| scans counted | rows with `protonatedmass` > 0 (charge-unknown rows have mass 0 and never match) |
| fragment, ppm | peak <= ion: `1000000 * (ion - peak) / ion <= ppm`; peak > ion: `1000000 * (peak - ion) / peak < ppm` (the v1.10 asymmetry); an ion's scan stops at the first peak above it that is outside the window |
| fragment, Da | `ion - peak <= Da` / `peak - ion <= Da` (same stop rule) |
| fragment hit | `log10(I) + 1 > 1` with I < 0 taken as 0; a repeated m/z takes the intensity of its first occurrence |
| exactness | every final hit decision uses GlycoMSP's float operations; binary search only finds candidates and only where it is provably exact; everything else runs a literal copy of GlycoMSP's loop. The tests compare the script with `mspvalidator_merger.peaks_ppm` / `_peaks_da` on random spectra, at the window edges (and one float either side), at extreme windows and on mixed-sign masses |
| Score A | hits / ion-list size, capped at 1.0, per scan |
| counting | `glycomsp` = one hit per in-window peak (what GlycoMSP writes); `unique` = each ion at most once. Neither removes the case of one peak counted for two close ions, which `B39_shared_peaks.py` measures |

**Checked against GlycoMSP itself:** on one PerMe Orbitrap MS2 run (alias `perme run A`), with its in silico library and a 282-ion list, the script's real counts equalled the CGA TSVs written by `run_pseudolabeling` on the same inputs, at 20 ppm / 10 ppm and at 300 ppm / 0.3 Da. On the same run the decoy count at 300 ppm exceeded the real count: a wide window is a per-dataset choice, which is what this check is for.

**Run:**

```
python tools/validation/cga_tolerance_decoy_count.py --alias run1 --converted <converted MS2 CSV> --library <in silico CSV> \
    --ions <ion list CSV/XLSX> --precursor 20ppm,300ppm,0.5Da --fragment 10ppm,0.3Da --out run1_decoy.csv
```

**Reading the result.** Use the decoy count as the chance estimate at the same tolerance: real minus decoy is roughly the matches above chance; `decoy / real` is roughly the fraction expected to be chance. Prefer the tolerance where real minus decoy stops growing while decoy keeps growing. The script reports counts only; it does not choose a tolerance and does not touch any GlycoMSP output.
