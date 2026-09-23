# Changelog

All notable changes to GlycoMSP are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

GlycoMSP is developed **additively**: the default behaviour of the published
v1.10 is preserved, and new behaviour ships behind explicit options, so results
produced with v1.10 remain reproducible.

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
