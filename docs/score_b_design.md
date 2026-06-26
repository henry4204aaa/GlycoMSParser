# Score B — Motif-Aware Re-ranking Design

**Module:** `src/msp_CGA_structscore.py` (version 0.7, last_update 20260411)
**Status:** optional, experimental — **not used in the main analyses**
**Scope of this document:** describes the *current* Score B implementation as it
exists in source. No source logic is changed by this document.

> **Manuscript alignment (Supplementary Data 2, verbatim).** "Score B is an optional
> experimental motif-aware re-ranking module in GlycoMSP and was not used in the main
> analyses. It evaluates candidate compositions using user-defined glyco-feature rules
> encoded in a curated Excel workbook. Conceptually, Score B combines a motif-support
> term with penalty terms for composition-inconsistent or unexpected fragment
> evidence, and the final score is bounded between 0 and 1. When a valid Score B
> workbook is supplied, GlycoMSP can calculate both the fragment-ion glycan-fitting
> score (Score A) and the motif-aware Score B, then use Score B for candidate
> re-ranking when sufficient motif evidence is available. If Score B evidence is
> insufficient, CGA ranking remains based on Score A and the ion-score threshold.
> Because Score B depends on user-defined motif rules and weights, it is provided as
> an optional advanced feature rather than a validated structural interpretation
> module. The reference workbook for positive-ion-mode permethylated N-glycans is
> provided with GlycoMSP and in the Zenodo package."

This document follows the manuscript's terminology: **glyco-feature** (the
user-defined motif rules). Where the source code uses a different literal identifier
for the same concept, that is noted explicitly (see
[§7 Implementation ↔ manuscript notes](#7-implementation--manuscript-notes)).

---

## 1. Overview

Score B is a per-candidate, motif-aware score in `[0, 1]`. For a given MS2 scan, each
CGA candidate composition is scored from observed fragment-ion evidence interpreted
through user-defined glyco-feature rules in a curated Excel workbook. The candidate
with the highest Score B is preferred for that scan; if Score B evidence is not
available, the scan falls back to the existing Score A (ion-score) selection.

The pipeline has three layers:

1. **Workbook load & validation** — `ScoreBLoader` parses six sheets into typed rules
   and compiles a target hierarchy index (`ScoreBConfig`).
2. **Per-spectrum motif evidence** — `ScoreBEngine` evaluates interpretation rules in
   two passes against the observed ions, yielding per-target evidence strengths.
3. **Per-candidate scoring** — `CandidateScoreBScorer` combines a motif-support term
   with two penalty terms and clamps the result to `[0, 1]`.

A thin runtime layer (`enrich_cga_tsv_with_score_b`) applies all three layers to an
existing CGA pseudolabel TSV and writes Score B columns back in place.

---

## 2. The Score B formula

The final candidate score is computed in `CandidateScoreBScorer.evaluate_candidate`:

```python
final_score = max(0.0, min(1.0,
    motif_support - composition_penalty - unexpected_penalty))
```

i.e.

> **b = max(0, min(1, support − penalty_comp − penalty_unexpected))**

The three terms are described below.

### 2.1 Support term (`motif_support`)

The motif-support term is the weighted mean of the **policy-selected scoring units**.

- Candidate flags are matched against `MotifPolicy` rows
  (`build_candidate_target_mask_trace`) to decide which targets (glyco-features) are
  relevant to this candidate. Each active policy row optionally expands to its
  hierarchy descendants (`include_descendants`).
- Each active policy row becomes a `ScoringUnit`:
  - a **leaf** unit (single target) takes that target's `partial_score`
    (= normalized per-target evidence strength), or
  - a **group** unit (multiple expanded targets) takes the **max** `partial_score`
    over its members.
- Support is the weight-weighted average of the unit scores
  (`_compute_motif_support_score_default`; all unit weights are currently `1.0`):

  ```
  support = Σ(unit.score · unit.weight) / Σ(unit.weight)
  ```

- **Fallback when no MotifPolicy row is active:** if there are no scoring units,
  `_compute_motif_support_score` averages `partial_score` over *all* targets. (When
  `candidate_flags` is `None`/empty, no policy rows activate, so this fallback path
  is what runs.)

Per-target evidence strength comes from the two-pass engine (§4); `partial_score` is
the target's `normalized_score = total_raw_score / total_max_score`, in `[0, 1]`.

### 2.2 Composition-inconsistency penalty (`composition_penalty`)

Computed in `_compute_composition_penalty` from the `CompositionConsistency` sheet.
For each rule whose `when_value` condition matches the candidate composition:

- `composition_field` selects one count from the candidate tuple
  (order `(H, N, S, G, K, F)`; see `COMPOSITION_INDEX`).
- `when_value` is a comparator string (`0`, `==0`, `>=1`, `>0`, `<=2`, `!=1`, …)
  evaluated against that count (`_matches_when_value`).
- `evidence_target_id` is resolved to concrete targets (a single target or a
  hierarchy group), and the **strongest** target evidence among them is taken as
  `max_evidence`.
- Rule mode (`_infer_comp_rule_mode`, `mode="auto"` by default):
  - `penalize_present` (e.g. `when_value` `0`/`==0`): `evidence_term = max_evidence`
    — penalize when forbidden evidence *is present*.
  - `penalize_absent` (e.g. `when_value` `>=1`/`>0`):
    `evidence_term = 1 − max_evidence` — penalize when expected evidence is *absent*.
- `rule_penalty = penalty · evidence_term`; penalties are **summed across rules and
  clamped to `1.0`**.

### 2.3 Unexpected-evidence penalty (`unexpected_penalty`)

Computed in `_compute_unexpected_penalty` from the `UnexpectedEvidencePenalty` sheet.

- A rule is active only when its `flag_id` is set true in `candidate_flags`
  (so with no flags, this penalty is `0`).
- `target_id` (optionally expanded to descendants) is resolved to concrete targets;
  `max_evidence` is the strongest evidence among them.
- `rule_penalty = penalty · max_evidence`; penalties are **summed and clamped to
  `1.0`**.

This penalizes a candidate when fragment evidence for a *forbidden* glyco-feature is
observed for a flag the candidate asserts.

---

## 3. Workbook sheet schema

The workbook is read by `ScoreBLoader._read_required_sheets`. All six sheets are
**required**; column headers are normalized through `COLUMN_ALIASES` (e.g.
`"Target ID"` → `target_id`), so human-friendly headers are accepted. Required
columns per sheet (`REQUIRED_SHEETS` / `REQUIRED_COLUMNS`):

| Sheet | Required columns | Purpose |
|---|---|---|
| `ionlist` | `glycotope_name`, `ion_type`, `fragmentation_mass`, `structural_identifier`, `charge_mode`, `derivatization` | Defines diagnostic fragment ions. `structural_identifier` must be unique. |
| `interpretation_v2` | `display_name`, `target_id`, `target_level`, `target_group`, `parent_target_id`, `ion_struct_id`, `ion_struct_target_id`, `ion_role`, `logic`, `min_hits`, `weight`, `charge_mode_required`, `derivatization_required` | Maps observed ions to glyco-feature targets. `min_hits ≥ 1`, `weight ≥ 0`, `logic ∈ {AND, OR}` (one logic mode per `target_id`). |
| `TargetHierarchy` | `parent_target_id`, `child_target_id` | Parent→child edges between targets/groups; compiled into a descendants index (cycles rejected). |
| `MotifPolicy` | `flag_id`, `target_id`, `activation_mode`, `include_descendants` | Selects which targets are relevant per candidate flag. `activation_mode ∈ {always, when_true/enable, when_false/disable}`. |
| `CompositionConsistency` | `rule_id`, `composition_field`, `when_value`, `evidence_target_id`, `penalty` | Composition-inconsistency penalties. `composition_field ∈ {H, N, S, G, K, F}`, `0 ≤ penalty ≤ 1`. Optional `mode` column. |
| `UnexpectedEvidencePenalty` | `flag_id`, `target_id`, `include_descendants`, `penalty` | Unexpected-evidence penalties. `0 ≤ penalty ≤ 1`. |

**Cross-sheet validation** (`_validate_cross_references`): every rule's
`ion_struct_id` must exist in `ionlist`; every `parent_target_id` /
`evidence_target_id` / penalty `target_id` must be a known concrete target or group
node. Optional `notes` columns are accepted on every sheet.

### Reference workbook & intended Score B path

The reference workbook example currently ships at:

```
templates/GlycoMSP_scoring_update_example_v7.xlsx
```

**Intended renamed path** (rename performed separately by Henry — *not* applied by
this document): a Score-B-named copy under `templates/`, e.g.

```
templates/GlycoMSP_ScoreB_reference_Nglycan_posmode_permethylated_v7.xlsx
```

The exact final filename is at Henry's discretion; the README and this doc should be
updated to point at it once the rename lands. The reference workbook is for
**positive-ion-mode permethylated N-glycans** and is also included in the Zenodo
package.

---

## 4. Per-spectrum evidence engine (two passes)

`ScoreBEngine.evaluate_global_motif_evidence` evaluates interpretation rules against
the observed ions of one spectrum:

- **Pass 1** — rows with empty `parent_target_id`; always evaluated (direct
  evidence). A row is *supported* when `charge_mode_required` and
  `derivatization_required` match the spectrum **and** the ion hit count
  ≥ `min_hits`; a supported row contributes its `weight`.
- **Pass 2** — rows with a `parent_target_id`; evaluated only if that parent acts as
  an **open context gate**. The gate opens when the parent target (or any expanded
  group descendant) has Pass-1 normalized evidence `> gate_threshold` (default
  `0.0`). `parent_target_id` is a *context gate*, not hierarchy inheritance.

Per target, evidence accumulates into `pass1_score/pass1_max` and
`pass2_score/pass2_max`; `normalized_score = total_raw / total_max ∈ [0, 1]`.
Target-level finalization (`_finalize_target_evidence_logic`) sets `logic_satisfied`:
**OR** targets need ≥1 supported row; **AND** targets need *all* active rows
supported. `partial_score` carries the continuous strength used by the support term.

Ion matching (`extract_observed_ion_hits_from_peaks`) keeps, per `ion_struct_id`, the
single peak with the smallest absolute ppm error within tolerance (default 20 ppm).

---

## 5. "Score B preferred, Score A fallback" selection rule

The selection rule lives in the CGA→Trainable build
(`src/mspfileloaderv14.py`, ~lines 4160–4255), consuming the Score B columns written
by `enrich_cga_tsv_with_score_b`:

1. **Score A eligibility** is computed first (legacy selection): rows passing the
   ion-score / ppm thresholds are ranked by `ion score` (or `ppm_error`) and the
   per-scan Top-N are taken as `legacy_sel`.
2. **If `score_b_selected` is present**, the rows with `score_b_selected == True`
   (one per scan, lowest `score_b_rank`) form `scoreb_sel`, tagged
   `selection_source = "score_b"`.
3. **Score A fallback** (`score_a_fallback`) supplies a selection only for scans
   **not** covered by a valid Score B-selected row.
4. The final selection is `concat(scoreb_sel, fallback_sel)`. Log line:
   *"ranking method change to: Score B preferred, Score A fallback."*
5. **If `score_b_selected` is absent, or no valid Score B rows exist**, the build uses
   the legacy Score A selection unchanged (*"Legacy Score A based"*).

Within `enrich_cga_tsv_with_score_b`, per-scan ranking sorts candidates by
`(final_score_b, motif_support_score, −composition_penalty, −unexpected_penalty)`
descending; rank 1 sets `score_b_selected = True` and `selected_composition`.

---

## 6. Non-fatal fallback behavior

Score B is designed to be additive and **non-destructive**:

- The selection rule prefers Score B *only per scan that has a valid Score B
  selection*; every other scan keeps its Score-A-only result. A scan with no valid
  Score B row, or a missing/empty `selected_composition`, simply falls back to Score
  A — the Score-A selection is never discarded.
- If the `score_b_selected` column is entirely absent (workbook not supplied / Score
  B not run), the trainable build silently uses the legacy Score A path. No exception
  is raised by the *selection* layer for missing Score B.

> **Caveat (loader is strict).** The *selection* layer degrades gracefully, but
> `ScoreBLoader` itself raises (`WorkbookValidationError` / `HierarchyValidationError`
> / `FileNotFoundError`) on a malformed or missing workbook. A supplied-but-invalid
> workbook therefore surfaces an error at enrichment time rather than silently
> producing a Score-A-only result; the graceful "leaves Score A intact" behavior
> applies once Score B columns are (or are not) present on the TSV the trainable build
> reads.

---

## 7. Implementation ↔ manuscript notes

Discrepancies and wording nuances to reconcile with Supplementary Data 2 (flagged for
Henry; **no code change implied**):

1. **"glyco-feature" vs `glycotope_name`.** The manuscript standardizes on
   *glyco-feature*. The workbook schema still uses the literal column name
   `glycotope_name` in the `ionlist` sheet (and `IonDef.glycotope_name`), and several
   identifiers/debug strings use *glycotope*. The concept is identical; this doc uses
   "glyco-feature" in prose but documents the real column name as-is. If the workbook
   header is to match the manuscript term, that is a workbook + alias change for a
   later version, not part of the v1.10 freeze.

2. **"sufficient motif evidence" — engine vs selection.** The manuscript says Score B
   is used "when sufficient motif evidence is available." In the implementation there
   are two distinct notions of sufficiency: (a) the *engine* gate
   (`gate_threshold`, Pass-1→Pass-2) controls per-rule evidence flow, and (b) the
   *selection* layer falls back to Score A based on the **presence** of a valid Score
   B-selected row (non-empty `selected_composition`), not on a Score B magnitude
   threshold. There is no numeric "Score B ≥ X" cutoff governing the
   prefer/fallback decision. The manuscript wording is consistent at the conceptual
   level but does not map to a single tunable threshold — worth a one-line
   clarification if reviewers ask.

3. **Two support-score code paths.** `_compute_motif_support_score_default` (weighted
   mean of policy-selected units) is the active path; `_compute_motif_support_score`
   provides the all-targets-average fallback when no MotifPolicy row activates, and
   `_compute_motif_support_score_experimental` is an inert placeholder that currently
   delegates to the default. `aggregation_debug` (strength/coverage/hybrid) is a
   debug payload and does **not** feed the final score. These are internal details,
   consistent with the manuscript's "combines a motif-support term"; noted so the
   description is not read as a single fixed formula.

4. **Score A is computed elsewhere.** This module computes Score B only; Score A (the
   fragment-ion glycan-fitting / ion score) is produced upstream and read from the
   `ion score` / `score_a` column. Consistent with the manuscript ("calculate both …
   Score A and … Score B").

No contradictions were found between the implementation and the bounded-[0,1],
support-minus-penalties description in Supplementary Data 2. The clamp, the two
penalty terms, and the prefer/fallback rule all match the manuscript text.
