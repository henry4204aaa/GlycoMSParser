"""Per-spectrum record and missing-metadata rules shared by the RAW and mzML readers.

version 0.42, 20261023 (READER-DIAG CP-A; v0.2 CP-C.3: readable empty spectra are not exported;
v0.3 CP-C.5: charge order of sources, `resolve_charge`, `validate_dataset_charge`;
v0.4 CP-D: per-scan status file, reader diagnostics files, conversion outcome and exit codes;
v0.41 CP-D.1: `Redactor` is public and replaces in one pass;
v0.42 R2 (Codex round 2): Redactor covers repr/escaped paths, short folders and case, stem as a whole word;
every string redacted when the diagnostics are written; `finish` never raises; failing scan recorded; rows of a
CSV that was never written marked `not_written`; one mass-eligibility rule; polarity from the filter string)

Both readers (`msprawextractor`, `mzmlreader`) build one `SpectrumRecord` per MSn spectrum
and apply the rules below, so the GUI and the extraction CLI see the same behaviour.

Rules:
- Peak arrays are read and validated independently of precursor metadata. A spectrum whose
  charge, precursor m/z or parent reference is missing keeps its peaks.
- Missing metadata is recorded as None plus an issue flag. Charge is never substituted silently:
  not 1, not -1, not the instrument-method default. A reported charge of 0 means "not determined"
  (Thermo convention) and is treated as unknown.
- Charge order (CP-C.5, Henry 2026-10-15): a valid reported charge is never overwritten; otherwise the
  comb's charge (RAW, when the user turns the comb on); otherwise the user's single dataset charge;
  otherwise unknown. `charge_source` records which one was used.
- The derived precursor mass is computed only when both precursor m/z and charge are known.
- A spectrum whose arrays cannot be decoded is not exported; it is counted and reported.
- A readable spectrum with no peaks is not exported either (as MSConvert drops it); it is counted
  as `empty` and listed under the `arrays_empty` issue. Henry's decision D2, 2026-10-13.

CSV serialisation of unknown values is interim (see `legacy_ms2_row`): it reuses the 0 sentinel
the RAW path already writes, until the CP-B export contract is agreed.
"""
import json
import math
import ntpath
import os
import platform
import re
import sys
import time
import traceback
from dataclasses import dataclass, field

version = "0.42"
last_update = 20261023

# Issue flags (stable strings; they appear in diagnostics output)
CHARGE_UNKNOWN = "charge_unknown"
PRECURSOR_MZ_MISSING = "precursor_mz_missing"
PARENT_REF_MISSING = "parent_ref_missing"
PRECURSOR_METADATA_ERROR = "precursor_metadata_error"
DERIVED_MASS_UNAVAILABLE = "derived_mass_unavailable"
ARRAYS_EMPTY = "arrays_empty"
ARRAYS_LENGTH_MISMATCH = "arrays_length_mismatch"
DECODE_ERROR = "decode_error"
# Status file only (R2, Henry D3 2026-10-23): read, but the conversion stopped before its CSV was written.
NOT_WRITTEN = "not_written"

# Where the charge used for the mass came from (CP-C.5)
CHARGE_SOURCE_REPORTED = "reported"
CHARGE_SOURCE_COMB = "comb"
CHARGE_SOURCE_DATASET = "dataset_value"
CHARGE_SOURCE_NONE = "none"
CHARGE_SOURCES = (CHARGE_SOURCE_REPORTED, CHARGE_SOURCE_COMB, CHARGE_SOURCE_DATASET, CHARGE_SOURCE_NONE)

ALL_ISSUES = (
    CHARGE_UNKNOWN, PRECURSOR_MZ_MISSING, PARENT_REF_MISSING, PRECURSOR_METADATA_ERROR,
    DERIVED_MASS_UNAVAILABLE, ARRAYS_EMPTY, ARRAYS_LENGTH_MISMATCH, DECODE_ERROR,
)


@dataclass
class SpectrumRecord:
    """One MSn spectrum as read from a backend, before CSV serialisation.

    `scan` and `parent_scan` keep the value the backend reports (the legacy CSV writes them
    with str()), `native_id` keeps the full native identifier when the backend has one.
    `mz`/`intensity` are None only when decoding failed; a readable empty spectrum is ().
    """
    ms_level: int
    scan: object
    native_id: object = None
    parent_scan: object = None
    precursor_mz: object = None
    monoiso_mz: object = None
    precursor_intensity: object = None
    charge: object = None
    charge_source: object = None     # one of CHARGE_SOURCES once resolved
    derived_mass: object = None
    polarity: object = None          # "+", "-" or None
    polarity_source: object = None   # "cv", "filter" or None (R2 D5)
    parent_ref: object = None        # the precursor's full spectrumRef, when the backend has one (R2-F8)
    representation: object = None    # "centroid", "profile" or None
    mz: object = None
    intensity: object = None
    issues: list = field(default_factory=list)
    errors: list = field(default_factory=list)   # [(stage, exception class, message)]

    def flag(self, issue):
        if issue not in self.issues:
            self.issues.append(issue)

    def note_error(self, stage, exc):
        self.errors.append((stage, type(exc).__name__, str(exc)))

    @property
    def decoded(self):
        return DECODE_ERROR not in self.issues

    @property
    def exportable(self):
        return self.decoded and ARRAYS_EMPTY not in self.issues

    @property
    def mass_eligible(self):
        """True only when a known charge gave a usable (finite, non-zero) precursor mass."""
        return self.decoded and mass_is_usable(self.charge, self.derived_mass)


def _known(value):
    """A number that is set, finite and not 0 (the CSV's unknown value)."""
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return False
    return math.isfinite(number) and number != 0.0


def mass_is_usable(charge, mass):
    """The one mass-eligibility rule (R2-F11): a known charge and a finite, non-zero mass."""
    return charge is not None and _known(mass)


# R2 D5 (Henry 2026-10-23): the polarity token of a Thermo filter string ("ITMS - p ESI ...", "FTMS + c NSI ...").
_POLARITY_RE = re.compile(r"(?:^|\s)([+-])\s+[pc](?:\s|$)")


def polarity_from_filter(scan_filter):
    """'+', '-' or None from a Thermo scan filter string."""
    match = _POLARITY_RE.search(scan_filter or "")
    return match.group(1) if match else None


def normalize_charge(value):
    """Return the reported charge as int, or None when it is absent or not determined.

    0, None, blanks and unparsable values all mean "unknown". Negative values are kept
    as reported; the sign convention belongs to the backend.
    """
    if value is None or isinstance(value, bool):
        return None
    try:
        as_float = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(as_float) or as_float != int(as_float):   # NaN, +-inf (R2-F10) or non-integer
        return None
    charge = int(as_float)
    return None if charge == 0 else charge


def validate_dataset_charge(value):
    """The user's single dataset charge as a positive int, or None when not given.

    Raises ValueError for 0, a bool, a non-integer or a non-number. A negative value is read as its
    magnitude: the converters take the polarity from the spectrum.
    """
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"dataset charge must be a non-zero integer, got {value!r}")
    if value == 0:
        raise ValueError("dataset charge must be a non-zero integer, got 0")
    return abs(value)


def resolve_charge(reported, comb_z=None, dataset_charge=None):
    """(charge, source) in the agreed order: reported, comb, dataset value, none.

    `reported` is the output of `normalize_charge` (None when absent or 0). The charge is None
    when no source has one.
    """
    if reported is not None:
        return reported, CHARGE_SOURCE_REPORTED
    if comb_z is not None:
        return comb_z, CHARGE_SOURCE_COMB
    if dataset_charge is not None:
        return dataset_charge, CHARGE_SOURCE_DATASET
    return None, CHARGE_SOURCE_NONE


def derive_mass(precursor_mz, charge, calc):
    """Derived precursor mass via the reader's own `calc(mz, charge)`, or None.

    `calc` is passed in so each reader keeps its legacy formula byte-for-byte. None is returned
    when either input is unknown or the formula has no value for this charge.
    """
    if precursor_mz is None or charge is None:
        return None
    return calc(precursor_mz, charge)


def finish_metadata_flags(rec):
    """Set the metadata issue flags from the record's fields. Call after filling them."""
    if rec.charge is None:
        rec.flag(CHARGE_UNKNOWN)
    if rec.precursor_mz is None:
        rec.flag(PRECURSOR_MZ_MISSING)
    if rec.parent_scan is None:
        rec.flag(PARENT_REF_MISSING)
    if rec.derived_mass is None:
        rec.flag(DERIVED_MASS_UNAVAILABLE)


def check_arrays(rec):
    """Flag empty or length-mismatched arrays. Decode errors are flagged by the reader."""
    if not rec.decoded:
        return
    if len(rec.mz) != len(rec.intensity):
        rec.flag(ARRAYS_LENGTH_MISMATCH)
    if len(rec.mz) == 0 and len(rec.intensity) == 0:
        rec.flag(ARRAYS_EMPTY)


def _sentinel(value):
    # Interim (pre CP-B): unknown numeric metadata is written as 0, as the RAW path does today.
    return 0 if value is None else value


def legacy_ms2_row(entry_no, rec):
    """The 9-column MS2 CSV row. Identical to the pre-READER-DIAG output when nothing is missing."""
    # With m/z and charge known, the mass cell is whatever the reader's formula returned;
    # otherwise it is the interim sentinel.
    if rec.precursor_mz is not None and rec.charge is not None:
        mass_cell = rec.derived_mass
    else:
        mass_cell = 0
    return (
        entry_no,
        _sentinel(rec.parent_scan),
        _sentinel(rec.precursor_mz),
        _sentinel(rec.monoiso_mz),
        _sentinel(rec.charge),
        mass_cell,
        rec.scan,
        tuple(rec.mz),
        tuple(rec.intensity),
    )


def legacy_ms3_row(entry_no, rec):
    """The 7-column MS3 CSV row: entry, parent MS2 scan, MS2 isolation m/z, MS2 mono m/z, MS3 scan, peaks."""
    return (
        entry_no,
        _sentinel(rec.parent_scan),
        _sentinel(rec.precursor_mz),
        _sentinel(rec.monoiso_mz),
        rec.scan,
        tuple(rec.mz),
        tuple(rec.intensity),
    )


def new_summary():
    """Counters every reader fills in the same shape."""
    return {
        "spectra_by_ms_level": {},
        "ms2": {"seen": 0, "exported": 0, "empty": 0, "failed": 0,
                "charge_known": 0, "charge_unknown": 0, "mass_eligible": 0,
                "charge_sources": {k: 0 for k in CHARGE_SOURCES},
                "issues": {k: [] for k in ALL_ISSUES},   # issue -> native scan numbers
                "errors": []},                           # [{"scan", "stage", "exception", "message"}]
        "ms3": {"seen": 0, "exported": 0, "empty": 0, "failed": 0,
                "issues": {k: [] for k in ALL_ISSUES},
                "errors": []},
    }


def count_level(summary, ms_level):
    lvl = summary["spectra_by_ms_level"]
    lvl[str(ms_level)] = lvl.get(str(ms_level), 0) + 1


def tally(summary, rec):
    """Add one MSn record to a summary made by `new_summary` (also counts its MS level)."""
    count_level(summary, rec.ms_level)
    key = {2: "ms2", 3: "ms3"}.get(rec.ms_level)
    if key is None:
        return
    block = summary[key]
    for issue in rec.issues:
        block["issues"][issue].append(rec.scan)
    for stage, exc_name, message in rec.errors:
        block["errors"].append({"scan": rec.scan, "stage": stage,
                                "exception": exc_name, "message": message})
    block["seen"] += 1
    if not rec.decoded:
        block["failed"] += 1
        return
    if not rec.exportable:          # readable but empty: counted, not exported
        block["empty"] += 1
        return
    block["exported"] += 1
    if key == "ms2":
        block["charge_known" if rec.charge is not None else "charge_unknown"] += 1
        if rec.charge_source is not None:
            block["charge_sources"][rec.charge_source] += 1
        if rec.mass_eligible:
            block["mass_eligible"] += 1


# ---------------------------------------------------------------------------------------------------------------
# CP-D (20261020): diagnostics persistence. Schema from the READER-DIAG handoff: B3 (per-scan status file),
# B6 (outcome and exit codes), B7 (diagnostics files). Henry, 2026-10-20: both readers write the files, so GUI
# conversions get them too (D3); exit codes as B6 (D4); spectra skipped without peaks are a note, not a warning.
# The CSVs are not touched by any of this.

OUTCOME_COMPLETE, OUTCOME_PARTIAL, OUTCOME_FAILED = "complete", "partial", "failed"
EXIT_OK, EXIT_USAGE, EXIT_WARNINGS, EXIT_PARTIAL, EXIT_FAILED = 0, 2, 3, 4, 5
EXIT_CODE_MEANING = {
    EXIT_OK: "complete: converted, nothing to check",
    EXIT_USAGE: "usage error (CLI only): bad arguments, nothing converted",
    EXIT_WARNINGS: "complete with warnings: converted; read the warnings (MS2 charges unknown, "
                   "no usable precursor mass, or profile MS2)",
    EXIT_PARTIAL: "partial: converted, but some spectra could not be decoded and were not exported",
    EXIT_FAILED: "failed: the input could not be read, the conversion stopped, or no MS2 spectrum with peaks "
                 "was exported",
}
NO_MASS_ELIGIBLE_EXPLANATION = ("No MS2 spectrum has a known precursor charge; charge-dependent matching "
                                "(CGA precursor search, precursor gates, mass feature) cannot use this file.")
# R2-F11: the same consequence when charges are known but no precursor m/z gives a usable mass.
NO_MASS_ELIGIBLE_MZ_EXPLANATION = ("No MS2 spectrum has a usable precursor mass: charges are known, but the precursor "
                                   "m/z is missing or 0; charge-dependent matching (CGA precursor search, precursor "
                                   "gates, mass feature) cannot use this file.")

STATUS_COLUMNS = ("ms_level", "scan", "native_id", "exported", "parent_scan", "precursor_mz", "charge",
                  "charge_source", "charge_status", "mass_status", "arrays_status", "n_peaks", "polarity",
                  "polarity_source", "representation", "peak_source", "comb_status", "comb_z", "comb_best_z", "comb_score",
                  "comb_margin", "comb_anchor", "mass_mz_source", "issues")
_COMB_FIELDS = ("comb_status", "comb_z", "comb_best_z", "comb_score", "comb_margin", "comb_anchor", "mass_mz_source")


def status_row(rec, comb_entry=None):
    """Status-file row (a dict) for one SpectrumRecord; `comb_entry` is its summary["charge"]["scans"] entry."""
    if not rec.decoded:
        arrays = "decode_error"
    elif ARRAYS_EMPTY in rec.issues:
        arrays = "empty"
    elif ARRAYS_LENGTH_MISMATCH in rec.issues:
        arrays = "length_mismatch"
    else:
        arrays = "ok"
    row = {"ms_level": rec.ms_level, "scan": rec.scan, "native_id": rec.native_id, "exported": rec.exportable,
           "parent_scan": rec.parent_scan, "precursor_mz": rec.precursor_mz, "charge": rec.charge,
           "charge_source": rec.charge_source, "mass": rec.derived_mass, "arrays_status": arrays,
           "n_peaks": len(rec.mz) if rec.decoded else None, "polarity": rec.polarity,
           "polarity_source": rec.polarity_source, "representation": rec.representation, "issues": list(rec.issues)}
    if comb_entry:
        row.update({k: comb_entry[k] for k in _COMB_FIELDS if k in comb_entry})
    return row


def _complete_row(row):
    """Fill the derived columns. Charge and mass status are set for exported MS2 rows only."""
    row = dict(row)
    exported_ms2 = row.get("ms_level") == 2 and row.get("exported")
    charge = normalize_charge(row.get("charge"))
    row["charge"] = charge
    if exported_ms2:
        row["charge_status"] = "known" if charge is not None else "unknown"
        row["mass_status"] = "eligible" if mass_is_usable(charge, row.get("mass")) else "unavailable"
        row.setdefault("charge_source", CHARGE_SOURCE_REPORTED if charge is not None else CHARGE_SOURCE_NONE)
        if row.get("comb_status") is None:
            row["comb_status"] = "not_run"
    if row.get("issues") is None:    # RAW rows: the same flag names, derived here
        issues = []
        if row.get("arrays_status") == "empty":
            issues.append(ARRAYS_EMPTY)
        if row.get("arrays_status") == "length_mismatch":
            issues.append(ARRAYS_LENGTH_MISMATCH)
        if row.get("exported"):
            if charge is None and row.get("ms_level") == 2:
                issues.append(CHARGE_UNKNOWN)
            if not _known(row.get("precursor_mz")):
                issues.append(PRECURSOR_MZ_MISSING)
            if not _known(row.get("parent_scan")):
                issues.append(PARENT_REF_MISSING)
            if exported_ms2 and row["mass_status"] == "unavailable":
                issues.append(DERIVED_MASS_UNAVAILABLE)
        row["issues"] = issues
    return row


def _cell(value):
    if value is None:
        return ""
    if value is True:
        return "yes"
    if value is False:
        return "no"
    if isinstance(value, (list, tuple)):
        return ",".join(str(v) for v in value)
    return str(value)


def write_status_tsv(path, rows, redact=None):
    """The per-scan status file: tab-separated, header STATUS_COLUMNS, unknown values blank.

    `redact` (a Redactor) is applied to every cell (R2-F3: a native ID can hold a file name)."""
    redact = redact or str
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(STATUS_COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(redact(_cell(row.get(col))) for col in STATUS_COLUMNS) + "\n")


def _count_rows(rows):
    """Counts, representation, polarity and affected scans from completed status rows."""
    counts = {"ms2": {"seen": 0, "exported": 0, "empty": 0, "failed": 0, "not_written": 0, "charge_known": 0,
                      "charge_unknown": 0, "mass_eligible": 0, "charge_sources": {k: 0 for k in CHARGE_SOURCES}},
              "ms3": {"seen": 0, "exported": 0, "empty": 0, "failed": 0, "not_written": 0}}
    representation = {lvl: {"centroid": 0, "profile": 0, "unknown": 0} for lvl in ("ms2", "ms3")}
    polarity = {lvl: {"+": 0, "-": 0, "unknown": 0} for lvl in ("ms2", "ms3")}
    affected = {"ms2": {}, "ms3": {}}
    for row in rows:
        lvl = {2: "ms2", 3: "ms3"}.get(row.get("ms_level"))
        if lvl is None:
            continue
        block = counts[lvl]
        block["seen"] += 1
        for issue in row["issues"]:
            affected[lvl].setdefault(issue, []).append(row.get("scan"))
        if row.get("arrays_status") == "decode_error":
            block["failed"] += 1
            continue
        if NOT_WRITTEN in row["issues"]:
            block["not_written"] += 1
            continue
        if not row.get("exported"):
            block["empty"] += 1
            continue
        block["exported"] += 1
        rep = row.get("representation")
        representation[lvl][rep if rep in ("centroid", "profile") else "unknown"] += 1
        pol = row.get("polarity")
        polarity[lvl][pol if pol in ("+", "-") else "unknown"] += 1
        if lvl == "ms2":
            block["charge_known" if row["charge_status"] == "known" else "charge_unknown"] += 1
            block["mass_eligible"] += row["mass_status"] == "eligible"
            block["charge_sources"][row["charge_source"]] = block["charge_sources"].get(row["charge_source"], 0) + 1
    return counts, representation, polarity, affected


def classify_outcome(counts, representation, failure=None, polarity=None):
    """(outcome, exit code, warnings, notes, explanation) by the B6 rules (Henry's 2026-10-20 warning set)."""
    ms2, ms3 = counts["ms2"], counts["ms3"]
    warnings, notes, explanation = [], [], None
    if failure is not None:
        outcome = OUTCOME_FAILED
        explanation = (f"The input could not be opened ({failure['exception']})." if failure["stage"] == "open" else
                       f"The conversion stopped at the {failure['stage']} stage ({failure['exception']}); "
                       f"outputs are missing or incomplete.")
        if ms2.get("not_written") or ms3.get("not_written"):
            warnings.append(f"{ms2.get('not_written', 0)} MS2 and {ms3.get('not_written', 0)} MS3 spectra were read, "
                            f"but the conversion stopped before their CSV was written (status file: issue "
                            f"{NOT_WRITTEN}).")
    elif ms2["exported"] == 0:
        outcome = OUTCOME_FAILED
        explanation = "No MS2 spectrum with peaks was exported."
    elif ms2["failed"] + ms3["failed"]:
        outcome = OUTCOME_PARTIAL
        warnings.append(f"{ms2['failed']} MS2 and {ms3['failed']} MS3 spectra could not be decoded and were not "
                        f"exported (status file: arrays_status decode_error).")
    else:
        outcome = OUTCOME_COMPLETE
    if ms2["exported"]:
        if ms2["charge_unknown"]:
            warnings.append(f"{ms2['charge_unknown']} of {ms2['exported']} exported MS2 spectra have no charge "
                            f"(written as chargeState 0, protonatedmass 0); charge-dependent steps cannot use them.")
        if ms2["mass_eligible"] == 0:
            reason = NO_MASS_ELIGIBLE_MZ_EXPLANATION if ms2["charge_known"] else NO_MASS_ELIGIBLE_EXPLANATION
            warnings.append(reason)
            explanation = explanation or reason
        if representation["ms2"]["profile"]:
            warnings.append(f"{representation['ms2']['profile']} exported MS2 spectra are profile data: their rows "
                            f"hold every profile point, not peaks. Re-convert with peak picking on MS2 (MSConvert: "
                            f"'peakPicking vendor msLevel=2-' as the first filter).")
    if ms2["empty"] or ms3["empty"]:
        notes.append(f"{ms2['empty']} MS2 and {ms3['empty']} MS3 spectra without peaks were not exported "
                     f"(as MSConvert does).")
    if polarity and polarity["ms2"].get("unknown"):   # R2-F13
        notes.append(f"{polarity['ms2']['unknown']} exported MS2 spectra state no polarity; their mass uses the "
                     f"positive-mode formula (status file: polarity blank).")
    for source, label in ((CHARGE_SOURCE_COMB, "the MS1 isotope comb"), (CHARGE_SOURCE_DATASET, "the dataset charge")):
        if ms2["charge_sources"].get(source):
            notes.append(f"{ms2['charge_sources'][source]} MS2 charges came from {label} "
                         f"(status file: charge_source {source}).")
    if outcome == OUTCOME_FAILED:
        code = EXIT_FAILED
    elif outcome == OUTCOME_PARTIAL:
        code = EXIT_PARTIAL
    else:
        code = EXIT_WARNINGS if warnings else EXIT_OK
    return outcome, code, warnings, notes, explanation


_WINDOWS_PATH_RE = re.compile(r"^(?:[A-Za-z]:[\\/]|\\\\)")
_WORD_EDGE = "A-Za-z0-9"


def _is_root(path, pathmod):
    return not path or pathmod.dirname(path) == path or bool(re.fullmatch(r"[A-Za-z]:[\\/]?|[\\/]+", path))


class Redactor:
    """Replace the input path, its folder, file name and stem with the output stem in any string.

    One regex pass over all tokens, longest first, so text already replaced is never scanned again
    (CP-D.1: replacing tokens one after another turned "bad.raw" with alias "bad_raw" into "bad_raw_raw").
    R2 (Codex round 2): every path also in its other-separator and repr/JSON-escaped forms (doubled backslashes),
    matched ignoring case; the parent folder at any length except a drive or root; the stem only as a whole word
    (stem "a" must not touch "Charge"). A stem that is a common word still replaces that word: over-redaction.
    """

    def __init__(self, input_path, stem):
        self.stem = stem
        self._pattern = None
        if not input_path:
            return
        input_path = os.fspath(input_path)
        pathmod = ntpath if (os.name == "nt" or _WINDOWS_PATH_RE.match(input_path)) else os.path
        full = pathmod.abspath(input_path) if pathmod is os.path else pathmod.normpath(input_path)
        base = pathmod.basename(input_path)
        paths = {full, input_path}
        parent = pathmod.dirname(full)
        if not _is_root(parent, pathmod):
            paths.add(parent)
        tokens = set()
        for path in paths:
            for form in (path, path.replace("\\", "/"), path.replace("/", "\\")):
                tokens |= {form, form.replace("\\", "\\\\")}
        tokens.add(base)
        words = {pathmod.splitext(base)[0]}
        tokens -= {"", stem}
        words -= {"", stem}
        if not tokens and not words:
            return
        # The output stem is an alternative too (replaced by itself), so text that was already redacted keeps its
        # alias when it is redacted again (alias "bad_raw" must not become "bad_raw_raw"). Longest first.
        edge = f"(?<![{_WORD_EDGE}])", f"(?![{_WORD_EDGE}])"
        parts = {t: edge[0] + re.escape(t) + edge[1] for t in tokens | words}
        if stem:
            parts[stem] = "(?-i:" + re.escape(stem) + ")"   # exact case: an identity replacement
        self._pattern = re.compile("|".join(parts[t] for t in sorted(parts, key=len, reverse=True)), re.IGNORECASE)

    def __call__(self, text):
        text = str(text)
        return self._pattern.sub(lambda _: self.stem, text) if self._pattern else text


def redact_tree(value, redact):
    """Redact every string inside dicts/lists/tuples (dict keys are kept: they are schema names). R2-F3."""
    if isinstance(value, str):
        return redact(value)
    if isinstance(value, dict):
        return {k: redact_tree(v, redact) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [redact_tree(v, redact) for v in value]
    return value


class ReaderDiagnostics:
    """Collects one conversion's status rows and context; `finish` writes the status file and diagnostics files.

    The reader calls `step(scan, ms_level)` for every spectrum it reads, `add(row)` for every MS2/MS3 spectrum, sets
    `stage` ("open", then "iterate" once the input is open, then "write"), `declared` and `instrument`, calls
    `wrote(level, path)` after each CSV is written, and calls `finish` once, also when the conversion raised. Only
    the stem (alias or input stem) and basenames reach the files; every string is redacted when they are written.
    """

    def __init__(self, reader, reader_version, reader_last_update, input_path, stem, fmt, entry=None,
                 settings=None, backend=None):
        self.started = time.time()
        self.reader, self.reader_version, self.reader_last_update = reader, reader_version, reader_last_update
        self.input_path, self.stem, self.fmt = input_path, stem, fmt
        self.entry = entry or "direct"
        self.settings = settings or {}
        self.backend = backend or {}
        self.redact = Redactor(input_path, stem)
        self.rows = []
        self.declared, self.processed, self.last_scan, self.last_level = None, 0, None, None
        self.levels = {}
        self.stage = "open"
        self.instrument = {}
        self.outputs = {}
        self.charge = None
        self.failure = None

    def step(self, scan, ms_level):
        self.stage = "iterate"
        self.processed += 1
        self.last_scan = scan
        if isinstance(ms_level, float) and ms_level.is_integer():   # pymzml may give 2.0
            ms_level = int(ms_level)
        self.last_level = ms_level
        self.levels[str(ms_level)] = self.levels.get(str(ms_level), 0) + 1

    def add(self, row):
        self.rows.append(row)

    def wrote(self, level, path):
        """Record one CSV after it was written ("ms2" or "ms3"). R2-F2: only written files are named."""
        self.outputs[f"{level}_csv"] = path

    def fail(self, exc):
        """Record the exception that stopped the conversion at the current stage, and the last scan started."""
        self.failure = {"stage": self.stage, "exception": type(exc).__name__, "message": self.redact(exc),
                        "last_scan_started": self.last_scan if self.stage != "open" else None,
                        "last_ms_level": self.last_level if self.stage != "open" else None,
                        "traceback": [f"{os.path.basename(f.filename)}:{f.lineno}"
                                      for f in traceback.extract_tb(exc.__traceback__)]}

    def paths(self, outdir):
        return {"json": os.path.join(outdir, f"{self.stem}.reader_diag.json"),
                "txt": os.path.join(outdir, f"{self.stem}.reader_diag.txt"),
                "status_tsv": os.path.join(outdir, f"{self.stem}.scan_status.tsv")}

    def _mark_not_written(self, rows):
        """R2 D3: after a failure, rows whose CSV was never written are not exported (issue not_written)."""
        if self.failure is None:
            return rows
        for row in rows:
            level = {2: "ms2", 3: "ms3"}.get(row.get("ms_level"))
            if level and row.get("exported") and f"{level}_csv" not in self.outputs:
                row["exported"] = False
                row["issues"] = list(row["issues"]) + [NOT_WRITTEN]
        return rows

    def build(self):
        rows = self._mark_not_written([_complete_row(r) for r in self.rows])
        counts, representation, polarity, affected = _count_rows(rows)
        outcome, code, warnings, notes, explanation = classify_outcome(counts, representation, self.failure, polarity)
        try:
            size = os.path.getsize(self.input_path)
        except (OSError, TypeError, ValueError):
            size = None
        finished = time.time()
        charge = self.charge or {}
        diag = {
            "json_type": "glycomsp.reader_diagnostics", "schema_version": "1.1.0",
            "outcome": outcome, "exit_code": code, "exit_code_meaning": EXIT_CODE_MEANING[code],
            "warnings": warnings, "notes": notes, "explanation": explanation,
            "input": {"alias": self.stem, "format": self.fmt, "size_bytes": size},
            "tool": {"reader": self.reader, "reader_version": str(self.reader_version),
                     "reader_last_update": self.reader_last_update,
                     "shared_rules": f"msp_spectrum_record {version}", "entry": self.entry},
            "environment": {"python": sys.version.split()[0], "platform": platform.platform(),
                            "backend": self.backend},
            "settings": self.settings,
            "instrument": self.instrument,
            "representation": representation, "polarity": polarity,
            "counts": {"spectra_by_ms_level": dict(sorted(self.levels.items())), **counts},
            "charge": {k: charge[k] for k in ("sources", "comb_status", "comb_vs_dataset_disagree", "mass_mz_source")
                       if k in charge},
            "affected_scans": affected,
            "coverage": {"spectra_declared": self.declared, "spectra_processed": self.processed,
                         "last_scan_processed": self.last_scan,
                         "aborted": self.failure is not None and self.failure["stage"] != "open"},
            "failure": self.failure,
            "outputs": {**{k: os.path.basename(v) if v else None for k, v in self.outputs.items()},
                        "partial": outcome == OUTCOME_PARTIAL or (self.failure is not None and bool(self.outputs))},
            "timing": {"started_utc": _utc(self.started), "finished_utc": _utc(finished),
                       "duration_s": round(finished - self.started, 1)},
        }
        return diag, rows

    def finish(self, outdir, exc=None):
        """Write the three files and return {outcome, exit_code, warnings, diagnostics: paths}. Never raises.

        R2-F6: everything is guarded, so a diagnostics error never replaces the conversion's own exception. When the
        files cannot be built or written, the return also holds `diagnostics_error` (the CLI then exits 5, Henry D2).
        """
        try:
            if exc is not None and self.failure is None:
                self.fail(exc)
            diag, rows = self.build()
        except Exception as build_exc:
            return self._report(None, {}, build_exc)
        paths = self.paths(outdir)
        diag["outputs"]["status_tsv"] = os.path.basename(paths["status_tsv"])
        error = None
        try:
            os.makedirs(outdir, exist_ok=True)
            write_status_tsv(paths["status_tsv"], rows, self.redact)
            write_diagnostics(paths["json"], redact_tree(diag, self.redact))
        except Exception as write_exc:   # the diagnostics must never hide the conversion's own result
            error, paths = write_exc, {}
        return self._report(diag, paths, error)

    def _report(self, diag, paths, error):
        if diag is None:
            diag = {"outcome": OUTCOME_FAILED, "exit_code": EXIT_FAILED, "warnings": [], "notes": []}
        result = {"outcome": diag["outcome"], "exit_code": diag["exit_code"], "warnings": diag["warnings"],
                  "diagnostics": paths}
        try:
            if error is not None:
                result["diagnostics_error"] = f"{type(error).__name__}: {self.redact(error)}"
                print(f"[diag] could not write the diagnostics files: {result['diagnostics_error']}")
            label = diag["outcome"] + (" with warnings" if diag["exit_code"] == EXIT_WARNINGS else "")
            print(f"[diag] outcome: {label} (exit code {diag['exit_code']}); warnings {len(diag['warnings'])}, "
                  f"notes {len(diag['notes'])}; {os.path.basename(paths.get('txt', '')) or 'no diagnostics file'}")
        except Exception:
            pass
        return result


def _utc(stamp):
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(stamp))


def write_diagnostics(json_path, diag):
    """Write `<stem>.reader_diag.json` and the text summary next to it (`.txt`, same stem)."""
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(diag, fh, indent=2, default=str)
    with open(json_path[:-len(".json")] + ".txt", "w", encoding="utf-8") as fh:
        fh.write(diagnostics_text(diag))


def update_diagnostics_outputs(json_path, failure=None, **outputs):
    """Rename entries in "outputs" (the CLI renames the CSVs after the reader wrote them) and rewrite both files.

    `failure` (R2-F2: a CLI step after the reader failed, e.g. moving a CSV) turns the outcome into failed, so the
    reader's diagnostics and the CLI's record agree.
    """
    with open(json_path, encoding="utf-8") as fh:
        diag = json.load(fh)
    diag["outputs"].update(outputs)
    if failure is not None:
        diag.update(outcome=OUTCOME_FAILED, exit_code=EXIT_FAILED, exit_code_meaning=EXIT_CODE_MEANING[EXIT_FAILED],
                    explanation=f"The conversion stopped at the {failure['stage']} stage ({failure['exception']}); "
                                f"outputs are missing or incomplete.")
        diag["failure"] = diag.get("failure") or failure
    write_diagnostics(json_path, diag)
    return diag


def diagnostics_text(diag):
    """Human summary of a diagnostics dict, including what the exit code means."""
    c2, c3 = diag["counts"]["ms2"], diag["counts"]["ms3"]
    label = diag["outcome"] + (" with warnings" if diag["exit_code"] == EXIT_WARNINGS else "")
    tool, inst, cov = diag["tool"], diag.get("instrument") or {}, diag["coverage"]
    lines = [f"GlycoMSP reader diagnostics - {diag['input']['alias']} ({diag['input']['format']}; "
             f"{tool['reader']} {tool['reader_version']}, {tool['shared_rules']}; entry {tool['entry']})",
             f"Outcome: {label} - exit code {diag['exit_code']}: {diag['exit_code_meaning']}",
             "Exit codes: " + "; ".join(f"{k} {v.split(':')[0]}" for k, v in sorted(EXIT_CODE_MEANING.items()))]
    if diag.get("explanation"):
        lines.append(f"Explanation: {diag['explanation']}")
    for title, items in (("Warnings", diag["warnings"]), ("Notes", diag["notes"])):
        lines.append(f"{title}:" + ("" if items else " none"))
        lines.extend(f"  - {item}" for item in items)
    lines += [
        f"MS2: seen {c2['seen']}, exported {c2['exported']}, without peaks {c2['empty']}, not decodable "
        f"{c2['failed']}, CSV not written {c2.get('not_written', 0)}; charge known {c2['charge_known']}, unknown {c2['charge_unknown']}; usable precursor mass "
        f"{c2['mass_eligible']}",
        f"MS2 charge sources: {c2['charge_sources']}",
        f"MS3: seen {c3['seen']}, exported {c3['exported']}, without peaks {c3['empty']}, not decodable "
        f"{c3['failed']}, CSV not written {c3.get('not_written', 0)}",
        f"Spectra by MS level: {diag['counts']['spectra_by_ms_level']}",
        f"Exported MS2 representation: {diag['representation']['ms2']}; polarity: {diag['polarity']['ms2']}",
        f"Coverage: processed {cov['spectra_processed']} of {cov['spectra_declared']} declared spectra; last "
        f"scan {cov['last_scan_processed']}; stopped early: {'yes' if cov['aborted'] else 'no'}",
        f"Instrument: {inst or 'not recorded'}",
        f"Settings: {diag['settings']}",
        f"Backend: {diag['environment']['backend']}; Python {diag['environment']['python']}",
        f"Files: {diag['outputs']}",
    ]
    if diag.get("failure"):
        f = diag["failure"]
        lines.append(f"Failure: stage {f['stage']}, {f['exception']}: {f['message']} "
                     f"[{' > '.join(f.get('traceback') or [])}]")
        if f.get("last_scan_started") is not None:
            lines.append(f"Stopped at or just after scan {f['last_scan_started']} (MS{f.get('last_ms_level')}); "
                         f"that scan has no status row when it was the failing one.")
    lines.append(f"Time: {diag['timing']['started_utc']}, {diag['timing']['duration_s']} s")
    return "\n".join(lines) + "\n"
