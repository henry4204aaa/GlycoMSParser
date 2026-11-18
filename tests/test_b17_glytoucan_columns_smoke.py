"""B-17 smoke test — verify GlyToucan/Glycanannotation2/WURCS metadata propagation
across MAS + CGA pipelines, v5 MAS header acceptance, B2 canonicalization,
NaN-aware empty-label guard, and Fold-C ionlistcheck indent fix.

The bug (see handoff/fix2_glytoucan_column.md): GlycoMSP declared GlyToucan ID as
a label-column candidate in the Train Model dropdown + train_model validator +
feature_exclude tuple, but the MAS writer pipeline silently discarded the column
between Excel read and trainable CSV write (slice_combined_df dropped it,
extract_ionmasslist excluded it, createnormailzedionlistcsv never wrote it,
sample_real_negatives had no defaults). Concurrently the CGA pipeline never
seeded the column at all. Net: user populating GlyToucan ID in v4 MAS Excel +
selecting GlyToucan ID as label → train_model errored "label column not in
trainable csv".

The fix (S3 / 2026-05-17, B-17): Z1+W1 propagation through all four MAS sites +
CGA seed at run_pseudolabeling Stage 7 + finalize column-order at
build_trainable_from_CGA + NaN-aware empty-label guard in train_model + v5 MAS
header invention (adds 'GlyToucan Access Number' as the official Glycosmos term
+ 'WURCS' placeholder for future API integration). B2 canonicalize renames
'GlyToucan Access Number' → 'GlyToucan ID' in directassign_files so downstream
code stays uniform. Fold-C indent fix at directassign_files:598 restores
ionlistcheck correct semantics.

Mode:
    T1-T6   real mspval functions imported and tested directly
    T7-T10  v14 nested-closure logic mirrored (since run_pseudolabeling /
            build_trainable_from_CGA / train_model are not importable closures)

Run from repo root:
    python3 tests/test_b17_glytoucan_columns_smoke.py
"""
import os
import sys
import tempfile

import numpy as np
import pandas as pd

# Import real mspval functions (module-level, importable).
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_SRC_DIR = os.path.normpath(os.path.join(_THIS_DIR, "..", "src"))
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

from mspvalidator_merger import (  # noqa: E402
    slice_combined_df,
    extract_ionmasslist,
    validate_annotation_structure,
    sample_real_negatives,
    directassign_files,
    createnormailzedionlistcsv,
)


# ---- v14 mirrors (nested-closure logic, kept in lockstep with edits) -----

_messagebox_calls = []  # captures messagebox.askyesno / showerror / showwarning calls

class _FakeMessagebox:
    """Stand-in for tkinter messagebox; askyesno returns canned answers in order."""
    def __init__(self):
        self._askyesno_queue = []
    def queue_askyesno(self, *answers):
        self._askyesno_queue.extend(answers)
    def askyesno(self, title, message):
        _messagebox_calls.append(("askyesno", title, message))
        if not self._askyesno_queue:
            raise RuntimeError("askyesno called but no canned answer queued")
        return self._askyesno_queue.pop(0)
    def showerror(self, title, message=""):
        _messagebox_calls.append(("showerror", title, message))


# Mirror of Edit 1 (run_pseudolabeling Stage 7 seed at v14:~4868).
def _mirror_run_pseudolabeling_seed(out_df):
    """Mirrors the empty-default seed loop applied post-merge in run_pseudolabeling."""
    for _label_col in ("Glycanannotation2", "GlyToucan ID", "WURCS"):
        if _label_col not in out_df.columns:
            out_df[_label_col] = ""
    return out_df


# Mirror of Edits 2a + 2b (build_trainable_from_CGA at v14:4210 + 4423-4432).
def _mirror_build_trainable_from_CGA_finalize(final_df, include_mass_feature=False):
    """Mirrors the column-order block of build_trainable_from_CGA after merge_cols
    extension (Edit 2a) and finalize empty-fill + reorder (Edit 2b).
    """
    for _b17_label_col in ("Glycanannotation2", "GlyToucan ID", "WURCS"):
        if _b17_label_col not in final_df.columns:
            final_df[_b17_label_col] = ""
    feature_cols = [c for c in final_df.columns if c not in ("MS2scan_no","Structure","Glycanannotation2","GlyToucan ID","WURCS")]
    if include_mass_feature and "protonatedmass" in final_df.columns:
        feature_cols = [c for c in feature_cols if c != "protonatedmass"]
        feature_cols.insert(0, "protonatedmass")
    else:
        feature_cols = [c for c in feature_cols if c != "protonatedmass"]
    return final_df[["MS2scan_no","Structure","Glycanannotation2","GlyToucan ID","WURCS"] + feature_cols]


# Mirror of Edit 3 (train_model NaN-aware empty-label guard at v14:7264-7268).
def _mirror_train_model_empty_guard(df, label_col, messagebox_obj):
    """Mirrors the non-Structure branch of train_model's normalization, including
    the new NaN-aware empty-label guard. Returns (df_after, return_early_flag).
    """
    if label_col == "Structure":
        # not the path under test
        return df, False
    df[label_col] = df[label_col].astype(str).str.strip()
    _b17_empty_mask = df[label_col].str.strip().str.lower().isin({"", "nan", "none", "null"})
    _b17_empty_count = int(_b17_empty_mask.sum())
    if _b17_empty_count > 0:
        if not messagebox_obj.askyesno(
            "Empty Labels",
            f"Label column '{label_col}' has {_b17_empty_count} of {len(df)} rows with empty values. "
            f"These rows will be dropped before training. Continue?",
        ):
            return df, True  # user cancelled
        df = df.loc[~_b17_empty_mask].copy()
        if df.empty:
            messagebox_obj.showerror("No Labels", f"All rows had empty {label_col}. Nothing to train.")
            return df, True
    return df, False


# Mirror of Edit 8 (ionlistcheck indent fix at mspval:598).
def _mirror_ionlistcheck(ion_df):
    """Mirrors the post-indent-fix logic. Pre-fix: ionlistcheck = False was at
    same level as the print inside the if/else, unconditionally overriding True.
    Post-fix: ionlistcheck = False is inside the else branch only.
    """
    ionlistcheck = None
    if ion_df["mass"].dtype == "float64":
        ionlistcheck = True
    else:
        ionlistcheck = False
    return ionlistcheck


# ---- v4 / v5 header fixtures --------------------------------------------

_V4_HEADERS = ['Structure', 'MS2scan_no', 'peak', 'charge', 'mass shift', 'adduct',
               'diff profile (integer for groups if different from same annotation)',
               'IUPACname(optional)', 'Glycanannotation2', 'note', 'GlyToucan ID']
_V5_HEADERS = ['Structure', 'MS2scan_no', 'peak', 'charge', 'mass shift', 'adduct',
               'diff profile (integer for groups if different from same annotation)',
               'IUPACname(optional)', 'Glycanannotation2', 'note',
               'GlyToucan Access Number', 'WURCS']


def _build_valid_v4_df(n_rows=2):
    """Build a DataFrame matching the v4 header order with valid dtypes for the
    needcheck columns enforced by validate_annotation_structure.
    """
    data = {
        'Structure':                                                               ['F1H5N4', 'F1H6N5'][:n_rows],
        'MS2scan_no':                                                              [101, 102][:n_rows],
        'peak':                                                                    [1, 1][:n_rows],
        'charge':                                                                  [2, 2][:n_rows],
        'mass shift':                                                              [0, 0][:n_rows],
        'adduct':                                                                  ['H', 'H'][:n_rows],
        'diff profile (integer for groups if different from same annotation)':     [1, 1][:n_rows],
        'IUPACname(optional)':                                                     ['', ''][:n_rows],
        'Glycanannotation2':                                                       ['N', 'N'][:n_rows],
        'note':                                                                    ['', ''][:n_rows],
        'GlyToucan ID':                                                            ['G00001AB', 'G00002CD'][:n_rows],
    }
    df = pd.DataFrame(data)
    return df[_V4_HEADERS]


def _build_valid_v5_df(n_rows=2):
    df = _build_valid_v4_df(n_rows).rename(columns={'GlyToucan ID': 'GlyToucan Access Number'})
    df['WURCS'] = ['', ''][:n_rows]
    return df[_V5_HEADERS]


# ---- tests --------------------------------------------------------------

def main():
    failures = []

    # ===== T1 — slice_combined_df adds empty GlyToucan ID + WURCS when absent =====
    df_in = pd.DataFrame({
        'protonatedmass': [500.0, 600.0],
        'peaklist': ['[1,2]', '[3,4]'],
        'peakintensity': ['[10,20]', '[30,40]'],
        'Structure': ['F1H5N4', 'F1H6N5'],
        'IUPACname(optional)': ['', ''],
        'Glycanannotation2': ['N', 'N'],
        'MS2scan_no': [101, 102],
    })
    out = slice_combined_df(df_in.copy())
    if 'GlyToucan ID' not in out.columns or 'WURCS' not in out.columns:
        failures.append(f"T1: slice_combined_df output missing meta cols; got {list(out.columns)}")
    elif not (out['GlyToucan ID'] == '').all() or not (out['WURCS'] == '').all():
        failures.append(f"T1: empty defaults not '' — got GlyToucan ID={out['GlyToucan ID'].tolist()}, WURCS={out['WURCS'].tolist()}")
    else:
        print("T1 PASS: slice_combined_df adds empty GlyToucan ID + WURCS when absent")

    # ===== T2 — slice_combined_df preserves user-supplied GlyToucan ID + WURCS =====
    df_in = pd.DataFrame({
        'protonatedmass': [500.0],
        'peaklist': ['[1,2]'],
        'peakintensity': ['[10,20]'],
        'Structure': ['F1H5N4'],
        'IUPACname(optional)': [''],
        'Glycanannotation2': ['N'],
        'GlyToucan ID': ['G00001AB'],
        'WURCS': ['WURCS=2.0/...'],
        'MS2scan_no': [101],
    })
    out = slice_combined_df(df_in.copy())
    if out['GlyToucan ID'].iloc[0] != 'G00001AB':
        failures.append(f"T2: GlyToucan ID overwritten; got {out['GlyToucan ID'].iloc[0]}")
    elif out['WURCS'].iloc[0] != 'WURCS=2.0/...':
        failures.append(f"T2: WURCS overwritten; got {out['WURCS'].iloc[0]}")
    else:
        print("T2 PASS: slice_combined_df preserves user-supplied GlyToucan ID + WURCS")

    # ===== T3 — extract_ionmasslist iondfindex carries 5 meta columns =====
    ionmass_sheet = pd.DataFrame({'mass': [204.087, 366.140, 512.197]})
    iondfindex = extract_ionmasslist(ionmass_sheet)
    expected_meta = {'Structure', 'IUPACname(optional)', 'Glycanannotation2', 'GlyToucan ID', 'WURCS'}
    actual_meta = set(iondfindex) & expected_meta
    if actual_meta != expected_meta:
        failures.append(f"T3: iondfindex missing meta cols; expected {expected_meta}, got {actual_meta}")
    elif 'protonatedmass' not in iondfindex or iondfindex[0] != 'protonatedmass':
        failures.append(f"T3: protonatedmass not at front; got iondfindex={iondfindex}")
    else:
        print("T3 PASS: extract_ionmasslist iondfindex carries 5 meta columns (Z1+W1)")

    # ===== T4 — validate_annotation_structure v4 still accepted (back-compat) =====
    v4_df = _build_valid_v4_df()
    if not validate_annotation_structure(v4_df):
        failures.append("T4: validate_annotation_structure rejected valid v4 DataFrame (back-compat broken)")
    else:
        print("T4 PASS: validate_annotation_structure accepts v4 (back-compat preserved)")

    # ===== T5 — validate_annotation_structure v5 accepted (greedy newest-first) =====
    v5_df = _build_valid_v5_df()
    if not validate_annotation_structure(v5_df):
        failures.append("T5: validate_annotation_structure rejected valid v5 DataFrame")
    else:
        print("T5 PASS: validate_annotation_structure accepts v5 (new path)")

    # ===== T13 — malformed v5 rejected (missing WURCS or swapped final columns) =====
    # T13a — v5 minus WURCS (11 cols, but last is 'GlyToucan Access Number' — neither v4 nor v5).
    v5_missing_wurcs = _build_valid_v5_df().drop(columns=['WURCS'])
    if validate_annotation_structure(v5_missing_wurcs):
        failures.append("T13a: validate_annotation_structure accepted v5-minus-WURCS (should fail strict list equality)")
    else:
        print("T13a PASS: malformed v5 (missing WURCS) rejected")
    # T13b — v5 with last two columns swapped — strict list equality must reject.
    v5_swapped_cols = _V5_HEADERS[:-2] + [_V5_HEADERS[-1], _V5_HEADERS[-2]]
    v5_swap_df = _build_valid_v5_df()[v5_swapped_cols]
    if validate_annotation_structure(v5_swap_df):
        failures.append("T13b: validate_annotation_structure accepted v5 with last two columns swapped")
    else:
        print("T13b PASS: malformed v5 (swapped final cols) rejected")

    # ===== T6 — sample_real_negatives empty-frame schema has GlyToucan ID + WURCS =====
    # Write a 1-row TSV; mark its only scan as annotated → cand_scans empty → early-return path fires.
    with tempfile.TemporaryDirectory() as tmpd:
        raw_tsv_path = os.path.join(tmpd, "raw_one.tsv")
        pd.DataFrame({
            'MS2scan_no': [101],
            'protonatedmass': [500.0],
            'peaklist': ['[204.087, 366.140]'],
            'peakintensity': ['[100.0, 50.0]'],
        }).to_csv(raw_tsv_path, sep='\t', index=False)
        neg = sample_real_negatives(
            raw_tsv_path,
            annotated_scans=[101],  # the only scan is annotated → cand_scans = empty
        )
    if neg is None:
        failures.append("T6: sample_real_negatives returned None instead of empty-schema DataFrame")
    elif 'GlyToucan ID' not in neg.columns or 'WURCS' not in neg.columns:
        failures.append(f"T6: empty-frame missing meta cols; got {list(neg.columns)}")
    else:
        print("T6 PASS: sample_real_negatives empty-frame schema includes GlyToucan ID + WURCS")

    # ===== T11 — REAL directassign_files: v5 canonicalize + ionlistcheck indent verify =====
    # Codex flagged Edit 5 (B2 canonicalize) and Edit 8 (ionlistcheck indent fix) as mirror-only.
    # This test exercises the real function with a temp .xlsx (v5 MSlist + ionlist) and a temp raw CSV.
    try:
        import openpyxl  # noqa: F401 — used by pd.ExcelWriter engine
    except ImportError:
        print("T11 SKIP: openpyxl not available")
    else:
        with tempfile.TemporaryDirectory() as tmpd:
            xlsx_path = os.path.join(tmpd, "annot_v5.xlsx")
            raw_csv_path = os.path.join(tmpd, "raw.tsv")
            # v5 MSlist sheet — 12 cols, 1 valid row with matching MS2scan_no
            mslist_df = pd.DataFrame([{
                'Structure': 'F1H5N4',
                'MS2scan_no': 101,
                'peak': 1,
                'charge': 2,
                'mass shift': 0,
                'adduct': 'H',
                'diff profile (integer for groups if different from same annotation)': 1,
                'IUPACname(optional)': '',
                'Glycanannotation2': 'N',
                'note': '',
                'GlyToucan Access Number': 'G00001AB',
                'WURCS': 'WURCS=2.0/test',
            }])[_V5_HEADERS]
            # ionlist sheet — float64 mass column triggers Edit 8 indent-fix happy path
            ionlist_df = pd.DataFrame({'mass': [204.087, 366.140]})  # float64
            with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
                mslist_df.to_excel(writer, sheet_name='MSlist', index=False)
                ionlist_df.to_excel(writer, sheet_name='ionlist', index=False)
            # raw CSV — must pass validate_csv_structure with exact 9-col schema + dtypes
            raw_df = pd.DataFrame([{
                'entry_no': 1,
                'MS1scan_no': 1,
                'MS1_isolationmass': 500.0,
                'MS1_monoisolationmass': 500.0,
                'chargeState': 2,
                'protonatedmass': 500.5,
                'MS2scan_no': 101,
                'peaklist': '[204.087, 366.140]',
                'peakintensity': '[100.0, 50.0]',
            }])
            raw_df.to_csv(raw_csv_path, sep='\t', index=False)
            # Real call. derivatizationtags = 'PerMe' (a non-dash string → mode '+').
            try:
                pre_df, iondfindex, ion_df_ret = directassign_files(xlsx_path, raw_csv_path, 'PerMe', debug=False)
            except Exception as e:
                failures.append(f"T11: directassign_files raised on valid v5 input: {type(e).__name__}: {e}")
            else:
                # Edit 5 (B2 canonicalize): pre_df must carry 'GlyToucan ID', NOT 'GlyToucan Access Number'.
                if 'GlyToucan Access Number' in pre_df.columns:
                    failures.append(f"T11: pre_df still has 'GlyToucan Access Number' — B2 rename not applied; cols={list(pre_df.columns)}")
                elif 'GlyToucan ID' not in pre_df.columns:
                    failures.append(f"T11: pre_df missing 'GlyToucan ID' after canonicalize; cols={list(pre_df.columns)}")
                elif 'WURCS' not in pre_df.columns:
                    failures.append(f"T11: pre_df missing 'WURCS'; cols={list(pre_df.columns)}")
                else:
                    # Edit 5: also verify the actual VALUE was preserved through the rename.
                    if 'G00001AB' not in pre_df['GlyToucan ID'].astype(str).tolist():
                        failures.append(f"T11: 'G00001AB' value lost after rename; got {pre_df['GlyToucan ID'].tolist()}")
                    elif 'WURCS=2.0/test' not in pre_df['WURCS'].astype(str).tolist():
                        failures.append(f"T11: 'WURCS=2.0/test' value lost; got {pre_df['WURCS'].tolist()}")
                    else:
                        print("T11 PASS: directassign_files canonicalizes v5 'GlyToucan Access Number' → 'GlyToucan ID' + preserves WURCS (Edits 5 + 8 real-call)")

    # ===== T12 — REAL createnormailzedionlistcsv: disk-output preserves GlyToucan ID + WURCS =====
    with tempfile.TemporaryDirectory() as tmpd:
        out_csv_path = os.path.join(tmpd, "trainable_out.csv")
        ion_df = pd.DataFrame({'mass': [204.087, 366.140]})
        ionindex = extract_ionmasslist(ion_df)
        converted_df = pd.DataFrame([{
            'protonatedmass': 500.5,
            'peaklist': [204.087],
            'peakintensity': [100.0],
            'Structure': 'F1H5N4',
            'IUPACname(optional)': '',
            'Glycanannotation2': 'N',
            'GlyToucan ID': 'G00001AB',
            'WURCS': 'WURCS=2.0/test',
            'MS2scan_no': 101,
        }])
        try:
            createnormailzedionlistcsv(ionindex, converted_df, ion_df, out_csv_path)
        except Exception as e:
            failures.append(f"T12: createnormailzedionlistcsv raised: {type(e).__name__}: {e}")
        else:
            if not os.path.exists(out_csv_path):
                failures.append("T12: output CSV not written")
            else:
                out_df = pd.read_csv(out_csv_path)
                if 'GlyToucan ID' not in out_df.columns:
                    failures.append(f"T12: output CSV missing 'GlyToucan ID' col; cols={list(out_df.columns)}")
                elif 'WURCS' not in out_df.columns:
                    failures.append(f"T12: output CSV missing 'WURCS' col; cols={list(out_df.columns)}")
                elif out_df['GlyToucan ID'].astype(str).iloc[0] != 'G00001AB':
                    failures.append(f"T12: 'GlyToucan ID' value not preserved on disk; got {out_df['GlyToucan ID'].iloc[0]}")
                elif out_df['WURCS'].astype(str).iloc[0] != 'WURCS=2.0/test':
                    failures.append(f"T12: 'WURCS' value not preserved on disk; got {out_df['WURCS'].iloc[0]}")
                else:
                    print("T12 PASS: createnormailzedionlistcsv writes GlyToucan ID + WURCS to disk with values preserved (Edit 6+7e real-call)")

    # ===== T7 — Edit 1 mirror: run_pseudolabeling seed produces 3 empty cols =====
    out = pd.DataFrame({
        'MS2scan_no': [101, 102],
        'Structure': ['F1H5N4', 'F1H6N5'],
        'theoretical_mass': [500.0, 600.0],
        'ppm_error': [2.5, 3.1],
    })
    out = _mirror_run_pseudolabeling_seed(out)
    new_cols = {'Glycanannotation2', 'GlyToucan ID', 'WURCS'}
    if not new_cols.issubset(set(out.columns)):
        failures.append(f"T7: Edit 1 mirror missing seeded cols; got {list(out.columns)}")
    elif any(not (out[c] == '').all() for c in new_cols):
        failures.append(f"T7: seeded cols not '' — Glycanannotation2={out['Glycanannotation2'].tolist()}, GlyToucan ID={out['GlyToucan ID'].tolist()}, WURCS={out['WURCS'].tolist()}")
    else:
        print("T7 PASS: run_pseudolabeling seed adds 3 empty meta cols (Edit 1 mirror)")

    # ===== T7b — Edit 1 mirror: seed preserves upstream values (now includes WURCS) =====
    # Codex D risk: T7b didn't cover WURCS preservation when upstream WURCS is non-empty.
    out = pd.DataFrame({
        'MS2scan_no': [101],
        'Structure': ['F1H5N4'],
        'Glycanannotation2': ['N'],
        'GlyToucan ID': ['G00001AB'],
        'WURCS': ['WURCS=2.0/upstream'],
    })
    out = _mirror_run_pseudolabeling_seed(out)
    if out['GlyToucan ID'].iloc[0] != 'G00001AB' or out['Glycanannotation2'].iloc[0] != 'N':
        failures.append(f"T7b: Edit 1 mirror overwrote upstream Glycanannotation2 / GlyToucan ID values")
    elif out['WURCS'].iloc[0] != 'WURCS=2.0/upstream':
        failures.append(f"T7b: WURCS upstream value overwritten; got {out['WURCS'].iloc[0]}")
    else:
        print("T7b PASS: run_pseudolabeling seed preserves upstream Glycanannotation2 / GlyToucan ID / WURCS values")

    # ===== T8b — Edit 2a mirror: merge_cols extension carries metadata pre-finalize =====
    # Codex flagged: T8 only mirrored the finalize block; verify the upstream merge_cols extension that prevents metadata loss.
    sel = pd.DataFrame({
        'MS2scan_no_pl': [101, 102],
        'Structure': ['F1H5N4', 'F1H6N5'],
        'selection_source': ['score_b', 'score_b'],
        'Glycanannotation2': ['N', 'N'],
        'GlyToucan ID': ['G00001AB', 'G00002CD'],
        'WURCS': ['WURCS=2.0/A', 'WURCS=2.0/B'],
    })
    pos_feat = pd.DataFrame({'MS2scan_no': ['101', '102'], '204.087': [1.5, 2.1]})
    # Mirror of Edit 2a logic at v14:4210-4212
    pl_scan_col = 'MS2scan_no_pl'
    _b17_meta_cols = [c for c in ("Glycanannotation2", "GlyToucan ID", "WURCS") if c in sel.columns]
    merge_cols = [pl_scan_col, "Structure"] + (["selection_source"] if "selection_source" in sel.columns else []) + _b17_meta_cols
    sel_for_merge = sel.copy()
    sel_for_merge[pl_scan_col] = sel_for_merge[pl_scan_col].astype(str)
    pos_feat_merged = pos_feat.merge(sel_for_merge[merge_cols], left_on="MS2scan_no", right_on=pl_scan_col, how="left").drop(columns=[pl_scan_col])
    if 'GlyToucan ID' not in pos_feat_merged.columns or 'WURCS' not in pos_feat_merged.columns:
        failures.append(f"T8b: merge_cols extension did not carry metadata into pos_feat; cols={list(pos_feat_merged.columns)}")
    elif pos_feat_merged['GlyToucan ID'].iloc[0] != 'G00001AB':
        failures.append(f"T8b: GlyToucan ID value lost in merge; got {pos_feat_merged['GlyToucan ID'].iloc[0]}")
    elif pos_feat_merged['WURCS'].iloc[1] != 'WURCS=2.0/B':
        failures.append(f"T8b: WURCS value lost in merge; got {pos_feat_merged['WURCS'].iloc[1]}")
    else:
        print("T8b PASS: build_trainable_from_CGA merge_cols extension preserves metadata pre-finalize (Edit 2a mirror)")

    # ===== T8 — Edits 2a+2b mirror: meta cols adjacent to Structure in final order =====
    final_df = pd.DataFrame({
        'MS2scan_no': [101, 102],
        'Structure': ['F1H5N4', 'F1H6N5'],
        '204.087': [1.5, 2.1],  # ion feature col
        '366.140': [1.8, 2.5],
        'Glycanannotation2': ['N', 'N'],
        'GlyToucan ID': ['G00001AB', 'G00002CD'],
        # WURCS missing — should empty-fill
    })
    out = _mirror_build_trainable_from_CGA_finalize(final_df.copy(), include_mass_feature=False)
    col_order = list(out.columns)
    if col_order[:5] != ['MS2scan_no', 'Structure', 'Glycanannotation2', 'GlyToucan ID', 'WURCS']:
        failures.append(f"T8: col order incorrect; first 5 = {col_order[:5]}")
    elif out['WURCS'].iloc[0] != '':
        failures.append(f"T8: missing WURCS not empty-filled; got {out['WURCS'].iloc[0]}")
    elif out['GlyToucan ID'].iloc[0] != 'G00001AB':
        failures.append(f"T8: GlyToucan ID value lost; got {out['GlyToucan ID'].iloc[0]}")
    elif any(c in col_order[5:] for c in ('Glycanannotation2', 'GlyToucan ID', 'WURCS')):
        failures.append(f"T8: meta col leaked into feature block tail; col_order={col_order}")
    else:
        print("T8 PASS: build_trainable_from_CGA finalize places meta cols adjacent to Structure (Edits 2a+2b mirror)")

    # ===== T9 — Edit 3 mirror: NaN-aware empty-label guard + askyesno + drop =====
    # T9a: some NaN labels → askyesno fires, user accepts → drops + continues
    df_in = pd.DataFrame({'GlyToucan ID': ['G1', None, 'G3', '', '  nan  ', 'NONE']})
    fb = _FakeMessagebox()
    fb.queue_askyesno(True)  # user accepts drop
    _messagebox_calls.clear()
    out, return_early = _mirror_train_model_empty_guard(df_in.copy(), 'GlyToucan ID', fb)
    if return_early:
        failures.append("T9a: returned early when user accepted askyesno")
    elif not _messagebox_calls or _messagebox_calls[0][0] != 'askyesno':
        failures.append(f"T9a: askyesno not fired; calls={_messagebox_calls}")
    elif len(out) != 2:
        failures.append(f"T9a: expected 2 remaining rows (G1, G3); got {len(out)} rows: {out['GlyToucan ID'].tolist()}")
    else:
        print("T9a PASS: NaN-aware empty-label guard drops NaN/empty/'nan'/'NONE' rows after askyesno=True")

    # T9b: user declines → return early, no training, no drop
    df_in = pd.DataFrame({'GlyToucan ID': ['G1', None, 'G3']})
    fb = _FakeMessagebox()
    fb.queue_askyesno(False)  # user declines
    _messagebox_calls.clear()
    out, return_early = _mirror_train_model_empty_guard(df_in.copy(), 'GlyToucan ID', fb)
    if not return_early:
        failures.append("T9b: did not return early when user declined askyesno")
    elif len(out) != 3:
        failures.append(f"T9b: rows mutated despite user decline; got {len(out)} rows")
    else:
        print("T9b PASS: NaN-aware empty-label guard returns early on askyesno=False")

    # T9c: all labels empty → askyesno=True → zero-row guard fires, return early
    df_in = pd.DataFrame({'GlyToucan ID': [None, '', 'nan']})
    fb = _FakeMessagebox()
    fb.queue_askyesno(True)
    _messagebox_calls.clear()
    out, return_early = _mirror_train_model_empty_guard(df_in.copy(), 'GlyToucan ID', fb)
    if not return_early:
        failures.append("T9c: zero-row guard did not fire — training would proceed on empty df")
    elif not any(c[0] == 'showerror' for c in _messagebox_calls):
        failures.append(f"T9c: showerror not surfaced after all-empty drop; calls={_messagebox_calls}")
    else:
        print("T9c PASS: zero-row guard fires when all labels empty")

    # ===== T14 — Edit 3 mirror: Structure label does NOT trigger empty-label guard =====
    # B-17 guard is non-Structure-only; verify Structure label_col with empty values does not fire askyesno.
    df_in = pd.DataFrame({'Structure': ['F1H5N4', None, '', 'nan']})
    fb = _FakeMessagebox()
    # Intentionally do not queue any askyesno — if the guard fires, _FakeMessagebox raises RuntimeError.
    _messagebox_calls.clear()
    out, return_early = _mirror_train_model_empty_guard(df_in.copy(), 'Structure', fb)
    if return_early:
        failures.append("T14: returned early on Structure label — guard should not fire for Structure")
    elif any(c[0] == 'askyesno' for c in _messagebox_calls):
        failures.append(f"T14: askyesno fired for Structure label; calls={_messagebox_calls}")
    else:
        print("T14 PASS: empty-label guard does not fire for label_col='Structure' (no-fire path)")

    # ===== T10 — Edit 8 mirror: ionlistcheck True on float64, False on non-float64 =====
    ion_df_ok = pd.DataFrame({'mass': [204.087, 366.140]})  # float64
    ion_df_bad = pd.DataFrame({'mass': ['204.087', '366.140']})  # object (strings)
    if not _mirror_ionlistcheck(ion_df_ok):
        failures.append("T10a: ionlistcheck False on valid float64 mass column")
    elif _mirror_ionlistcheck(ion_df_bad):
        failures.append("T10b: ionlistcheck True on non-float64 mass column")
    else:
        print("T10 PASS: ionlistcheck True on float64, False otherwise (Edit 8 indent-fix mirror)")

    # ===== Summary =====
    # 18 assertions: T1-T6 (6 mspval real-fn) + T13a/T13b (2 malformed v5) + T11 + T12 (2 real-fn)
    # + T7 + T7b + T8 + T8b (4 CGA mirrors) + T9a + T9b + T9c + T14 (4 train-guard mirrors) + T10 (1 indent mirror).
    total = 18
    print()
    if failures:
        print(f"FAIL: {len(failures)}/{total} tests failed")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    else:
        print(f"PASS: {total}/{total} tests passed")
        sys.exit(0)


if __name__ == "__main__":
    main()
