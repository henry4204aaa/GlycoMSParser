# Changelog

All notable changes to GlycoMSP are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

GlycoMSP is developed **additively**: the default behaviour of the published
v1.10 is preserved, and new behaviour ships behind explicit options, so results
produced with v1.10 remain reproducible.

## [1.14] – 2026-11-01

Publication build for the manuscript revision: the first version that can
annotate low-resolution (ion-trap) MS2 data. The header and window titles
report v1.14 and no longer say "preview". With every new field left at its
default, results produced with v1.10 are unchanged (verified on a v1.10 project
end to end).

### Added
- **Configurable CGA mass tolerances.** *CGA Setup* has a *Mass tolerances*
  frame with two fields — *Precursor match* (default 20 ppm) and *Score A
  fragment match* (default 10 ppm) — each accepting a value in ppm or Da. A Da
  window is symmetric around the target mass. Invalid or empty values stop the
  run with a message before anything is written.
- **Tolerances are recorded.** A run with non-default tolerances is marked in
  the log pane, the run log and the CGA method JSON (`parameters.cga`), and the
  value is carried forward: *CGA → Trainable* and the prediction dialogs are
  prefilled from the method or model record, and a different value entered by
  hand produces a notice. The data files themselves (CGA TSV, trainable CSV,
  prediction tables) are unchanged in layout and carry no marker.
- **Precursor gate in *CGA → Trainable*.** The 10 ppm precursor gate that
  v1.10 applied silently when building a training set is now a visible,
  adjustable field (checkbox + ppm) in box *4) Thresholds / selection*; the
  default reproduces v1.10.
- **Unlabeled-dataset dialog** accepts decimal values (1–500 ppm) and is
  prefilled from the selected model's tolerance record.
- **Controlled validation scripts** (`tools/validation/`, not part of the
  application): a decoy-library count that reports how many scans pass the
  precursor and Score A criteria against the real library and a +3.37 Da
  shifted copy at chosen tolerances, and two Score A peak-counting checks.
  They document the method used to choose tolerances for ion-trap data.
- **Negative-mode native N-glycan ion lists** in `templates/`:
  `GlycoMSP_ionlist_native_neg_14.xlsx` (primary set) and
  `GlycoMSP_ionlist_native_neg_candidate_97.xlsx` (extended set), for use as the
  Score A ion list and the feature ion list with *Native(Reduced)*, negative mode.

### Changed
- The CGA Setup window re-packs its footer buttons so the new frame does not
  push *Start CGA analysis* / *Close* off a small screen.
- The user manual (`docs/glycomsp_user_manual_v1.0.md`) describes the new
  fields, the ion-trap workflow (§7.2), the defaults table (§5.5, corrected
  against the code) and the method JSON provenance block (§5.2).

### Known limitations
- Score A counts peak–ion pairs: with a wide fragment window (0.3 Da or more)
  and a dense ion list, one peak can be counted for two neighbouring ions. With
  the 14-ion native list no scan is affected; with a 97-ion list about one scan
  in eight is. The counts are reported by the validation scripts.
- The *ion list ppm* feature window in *CGA → Trainable* and in prediction is
  ppm-only and capped at 500 ppm (0.13 Da at m/z 260), so low-m/z ion-trap
  fragments get a narrower window than the CGA step's Da setting. A ppm-or-Da
  selector for feature extraction is planned.
- Prediction run summaries report *Accepted 0*, τ 0 and zero probability
  statistics when thresholding is off (the default); the predicted labels and
  per-class counts are correct. The summary will be corrected.
- The CGA result TSV carries every converted scan, so its size follows the raw
  file (about 1.7 GB for a run of 100,000 MS2 scans). A matched-rows-only
  output is planned.
- Native masses apply to the CGA route only; rows with unknown charge are kept
  with mass `0` and are not yet excluded downstream (unchanged from 1.13).


## [1.13] – 2026-10-31 (preview)

Preview of the second group of v1.20 changes: data readers that report what
they read and native (underivatized) glycan masses in the CGA route. With all
new options left at their defaults, results produced with v1.10 are unchanged
except for the corrections listed under *Fixed*. The GUI header reports v1.13;
the user manual still describes v1.10.

### Added
- **Conversion diagnostics.** Every RAW or mzML conversion now writes, next to
  the two CSV files, a per-scan status table (`<name>.scan_status.tsv`) and a
  diagnostics report (`<name>.reader_diag.json` / `.txt`) with the outcome
  (complete, complete with warnings, partial, failed), counts of exported and
  skipped spectra, known and unknown precursor charges, and the reason for any
  failure. The GUI shows the outcome in its log; a failed conversion is never
  saved as if it had succeeded, and a partial one is saved with a warning.
- **Spectra without a precursor charge are kept.** When the instrument did not
  assign a charge (common for ion-trap data), the spectrum is exported with its
  peaks and precursor m/z, charge `0` and mass `0`, and marked as "unknown
  charge" in the status table, instead of stopping the conversion.
- **Conversion options** (dialog before each GUI conversion; both off by default):
  estimate the precursor charge (1–3) from the isotope pattern of the parent
  MS1 scan, or supply one charge for the whole dataset. A charge reported by the
  instrument is never overwritten. The choice is recorded in the project
  metadata.
- **Native glycan masses.** Choosing derivatization *Native* or *Native(Reduced)*
  now gives correct underivatized precursor masses in the CGA library (previously
  the masses were zero or silently permethylated). Monosaccharide, reducing-end,
  sulfate and phosphate masses are documented in the module page.
- **`src/msp_CLI_precursor_mass.py`** — prints the theoretical precursor mass of
  a composition exactly as the CGA library computes it (`--composition`,
  `--deri perme|native`, `--reduced`, `--mode +|-`, `--json`, `--table`,
  `--selftest`).
- **`src/msp_extract_cli.py`** — command-line conversion (`convert`) and a
  read-only probe of the Thermo backend for selected scans (`inspect`), using the
  same readers and options as the GUI; exit codes 0/3/4/5 mirror the outcomes.
  General-purpose scripting support beyond these two commands is still in progress.

### Changed
- Readers (`msprawextractor`, `mzmlreader`) share one set of rules for missing
  metadata and one CSV writer; the 9-column MS2 / 7-column MS3 format is
  unchanged.

### Fixed
- Negative-mode precursor masses used the positive-mode formula.
- RAW files from ion-trap instruments exported empty peak lists; peaks are now
  read from the vendor mass list when no label data exist.
- MS3 rows carried the parent's peaks and wrong header names; they now carry
  their own peaks.
- Spectra with no peaks were written as empty rows; they are no longer exported.
- mzML files without a polarity term, or with the precursor on a second
  `selectedIon`, were misread; polarity now comes from the filter string and
  precursor fields from the first selected ion of the first precursor.
- A conversion saved across different volumes failed on some systems.

### Known limitations
- The CGA precursor match (20 ppm) and Score A fragment match (10 ppm) are fixed
  in this version; ion-trap data with isolation-m/z precursors therefore receive
  almost no annotations. (Resolved in 1.14 by the configurable tolerances.)
- Native masses apply to the CGA route only; manual (MAS) annotation still uses
  permethylated residue masses.
- Rows with unknown charge are kept in the CSV with mass `0`; downstream steps do
  not yet exclude them.



## [1.12] – 2026-10-01 (preview)

Preview build on top of 1.11; GUI panels and the user manual remain at the v1.10 level except where noted. With the reference-database option off, all outputs are identical to v1.10.

### Added
- **Manage References window** (Prepare Dataset → *Manage References…*): shows the active reference database (path, name/version, entry counts by status, last update, content hash, log size), lets you select another database file for the session, run *Prebuild / Update from source…* on a composition table or MAS workbook with a live log and Cancel, rebuild the transaction log, and set the chosen database as default. Prebuild runs the existing command-line tool in a child process, so the window stays responsive; batches already committed are kept when a run is cancelled and a rerun resumes from them.
- **Active database label** in Prepare Dataset, next to the fill option.
- **Reference database info** on a sample's right-click menu: which database filled the sample's CGA/MAS results, how many compositions were filled, and when.
- **Settings file** (`settings.json` in the GlycoMSP data folder) that remembers the fill option and the chosen database — written only when you press *Set as default* or tick *Remember these settings*; never created by browsing.
- **`references/` folder** with a showcase reference-database snapshot and its data attribution (GlyTouCan / GlyCosmos, CC BY 4.0). Select it through *Manage References…*; see `references/README.md`.
- Placeholders (disabled) for the next steps: *Check again* for unresolved rows, *Resolve missing references*, and *Update GlyToucan ID / WURCS in this file…*.

### Changed
- The reference database now records the source, time and status of its last prebuild, and its identity (hash, entry count, version) is read in a single consistent snapshot.
- Score B composition parsing is case-sensitive: lowercase `s` (sulfate) and `p` (phosphate) are no longer read as NeuAc / another residue. Tokens that Score B does not know are left out of scoring, recorded per row, and reported once at the end of the run instead of stopping it.
- The MAS KDN residue mass is aligned with the CGA value (320.1471).

### Fixed
- `msp_CLI_make_trainable_with_negatives.py` imported its decoy helper from the wrong module.
- Two functions in the in-silico marker and pre-training normaliser modules were defined twice; the shadowed copies are removed (no behaviour change).

### Notes
- The reference-database option is still off by default. Prebuild needs network access and the `certifi` package (listed in `src/requirements.txt`); macOS Python builds without a system certificate store fail with `CERTIFICATE_VERIFY_FAILED` until it is installed.
- Shipped database snapshots carry no compatibility promise; merging a shipped snapshot into your own database is planned for a later preview.

## [1.11] – 2026-09-23 (preview)

Preview of the first v1.20 feature. The application header and window title
report v1.11; GUI panels and the user manual still describe v1.10 behaviour.
v1.10 remains the release cited by the manuscript and archived on Zenodo.

### Added
- **GlyTouCan / WURCS reference database.** GlycoMSP can annotate glycan
  compositions with their GlyTouCan accession and WURCS string from a local
  SQLite reference database (default `glycomsp_reference.db` under
  `~/.glycomsp/references/` on macOS/Linux and `%LOCALAPPDATA%\GlycoMSP\references\`
  on Windows; override the database and log locations with the
  `GLYCOMSP_REFERENCE_DB` and `GLYCOMSP_GLYTOUCAN_LOG` environment variables).
  Every successful fetch is also appended to a JSONL transaction log
  (`glytoucan_cache.jsonl` in the same data directory).
- **`src/msp_CLI_glytoucan_prebuild.py`** — the only component that contacts the
  network. It reads the compositions in a CGA TSV, trainable CSV, MAS workbook
  (`MSlist` sheet by default) or in-silico library (`--source`, `--column`,
  `--sheet`), fetches only those not yet in the database through the GlyCosmos
  composition API (≤ 30 compositions per request, 10 s pause between requests),
  and stores the results batch by batch. A network failure stops further requests
  but keeps everything already fetched; re-running resumes from what is missing,
  and a second run over an unchanged source makes zero API calls. Options:
  `--db`, `--log`, `--no-network` (report what would be fetched), `--batch-size`,
  `--pause`, `--timeout`, `--summary`, `--json`, `--rebuild-log` (regenerate the
  transaction log from an existing database, keeping the previous log as `.bak`;
  it never creates or modifies a database). Exit code 0 covers completion and
  network fail-soft; exit code 2 signals a storage problem (database write
  failed, transaction log incomplete, or `--rebuild-log` refused). See `--help`.
- **"Fill GlyToucan ID / WURCS from reference DB"** checkbox in *Prepare Dataset*
  (off by default). When ticked, empty GlyToucan ID / WURCS cells are filled from
  the local database — never from the network — in: CGA annotation tables, the
  MAS merge (values already in your sheet always win), the CGA → trainable
  conversion, and prediction output (two additional columns,
  `Predicted_GlyToucan_ID` and `Predicted_WURCS`; `Predicted_Label` is unchanged).
  A row that already carries your own accession keeps its WURCS cell untouched,
  and a row with your own WURCS keeps its accession cell untouched, so the
  database never pairs a generic value with a structure-specific one.
  Each fill writes a one-line summary to the main-window log.
- **Method JSON provenance.** When a fill happened, the CGA or MAS method JSON
  gains an `artifacts.glytoucan_resolution` block (counts, database identity,
  timestamp). Files produced with the option off are identical to v1.10;
  `schema_version` stays `1.0.0`.
- Reference-database rows carry a status: `resolved` (WURCS + accession),
  `wurcs_only` (WURCS available, composition not registered in GlyTouCan — seen
  for some KDN-, HexA- and sulfate-containing compositions), or `unresolvable`.
  Values you entered yourself are never overwritten.
- Third-party data attribution for GlyTouCan / GlyCosmos (CC BY 4.0) in
  `LICENSE.md` and `README.md`.

### Changed
- `mspvalidator_merger` 0.94 → 0.95: the MAS writer accepts an optional
  enrichment hook. With the option off its output is byte-identical to v1.10.
- `certifi` added to `src/requirements.txt` and `src/requirementspy313.txt`
  (HTTPS certificate bundle for the prebuild CLI; certificate verification is
  never disabled).
- An obsolete test that referenced a removed module was moved to
  `archive/tests_retired/test_p1_validation.py`.

### Known issues
- The fill option is not remembered between launches (by design until the
  *Manage References* panel arrives in a later v1.20 preview).
- No reference database is shipped yet; run the prebuild CLI once per dataset
  to populate your local copy. The final database location may still change.
- `wurcs_only` rows keep an empty accession until the composition is
  registered in GlyTouCan; a "check again" action is planned.
- Pre-existing (not introduced in 1.11): a loaded method whose family is
  recognised as UNKNOWN hides its `artifacts` until *Export Method v1 (MAS)* is
  used once.
- The user manual (v1.0) does not yet cover the reference-database features.

## [1.10] – 2026-06-26

Published, frozen release cited by the manuscript. Reproducibility archive:
<https://doi.org/10.5281/zenodo.20823042>. Earlier history is recorded in the
module headers under `src/`.
