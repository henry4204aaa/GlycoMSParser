# references/ — shipped GlyTouCan / WURCS reference databases

This folder mirrors the per-user runtime folder `<data_dir>/references/` that GlycoMSP uses for its local
reference database (`<data_dir>` = `%LOCALAPPDATA%\GlycoMSP` on Windows, `~/.glycomsp` elsewhere).

## What belongs here

| File | Meaning |
|---|---|
| `glycomsp_reference_<db_version>.db` | A reference database snapshot, e.g. `glycomsp_reference_v0.02.20261001.db`. SQLite; schema `1.0.0`. |
| `ATTRIBUTION.md` | Data source and licence (CC BY 4.0) for the composition data inside the databases. |
| `README.md` | This file. |

The v1.12 preview ships one **showcase snapshot** (`glycomsp_reference_v0.01.20260921.db`, built from a zebrafish
brain permethylated N-glycan CGA table; entries are either `resolved` with a GlyTouCan accession or `wurcs_only`). A
published database is a snapshot of a database in use at the time. Later design or version changes may alter its
structure and behaviour; there is **no compatibility promise** for shipped snapshots.

## How to use a shipped database

1. Open **Prepare Dataset → Manage References…**.
2. Press **Select database…** and choose the `.db` file. GlycoMSP opens it read-only and checks that it is a
   GlycoMSP reference database; anything else is refused and the previous selection is kept.
3. The chosen database is used for this session. Press **Set as default** and tick **Remember these settings**
   to use it again at the next launch.
4. Tick **Fill GlyToucan ID / WURCS from reference DB on run** to fill empty `GlyToucan ID` / `WURCS` cells in
   CGA results, MAS sheets, trainable CSVs and predictions. With the option off, outputs are identical to v1.10.

A selected database that is not the default one keeps its own transaction log beside it
(`<db_folder>/<db_name>.jsonl`). Merging a shipped database into your own local database is not available yet
(ticket T1-f).

## Annotation without registration

GlycoMSP does not ask you to register structures on GlyTouCan. A row with a blank accession and your own WURCS
is a complete local annotation, not a gap.

## Attribution

GlyTouCan / GlyCosmos composition data, CC BY 4.0 (https://glytoucan.org, https://glycosmos.org). See
`ATTRIBUTION.md`.
