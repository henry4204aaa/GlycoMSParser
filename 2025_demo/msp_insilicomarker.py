#part of it was copied from GPT5 answer based on the discussion


import numpy as np
import csv
import mspvalidator_merger as validator
import pandas as pd

#fast converting fixed path
converted_csv = r"C:\Users\Sakazuki\Desktop\Khoolab_2025data\U937cells_NG_new\ms2_U937_20250729_U937_ST1OE_NGneu.csv"
composition_list = r"G:\其他電腦\My Computer\GlycoMSParser\U937NG.csv" #this is fault file, mass wrong  #fix wrong mass before adding comp list file
ion_sheet = r"C:\Users\Sakazuki\Downloads\ionlist_from_zebrafish_NG.csv"
#v20250815 we have old version to use
#need to restrict numbers, currently set to float, could be better to limit 4 dec
def getppm(refmass, comparemass, ppm):
    if refmass == 0 or comparemass == 0:
        return np.nan
    elif refmass > comparemass:
        return ((refmass-comparemass)/refmass)* 1e6
    elif comparemass > refmass:
        return ((comparemass-refmass)/comparemass)* 1e6
    else:
        raise ValueError("Invalid input of mass or ppm value provided")

#extract_ionmasslist in mspvalidator
def extract_ionmasslist(ionmass_sheet):
    ion_df = ionmass_sheet[["mass"]]
    iondfindex = ion_df.values.flatten().tolist()
    return iondfindex

#sort mass for later searching
def parse_composition(compfile):
    compfile1 = pd.read_csv(compfile)
    print(compfile1.head())
    sort = compfile1.sort_values('Mass')
    print(sort.head())
    #print("sort composition list by mass ascend")
    #sortcomp = compfile1.sort_values('mass').reset_index(drop=True)
    #print(sortcomp.head(10))

#test line
parse_composition(composition_list)

#v20250817 we have old version, not sure where it was
def readunlabeledcsv(csv, ionlist, compositionlist):
    #validation of data format
    loadcsv = pd.read_csv(csv, sep='\t')
    valid_csv = validator.validate_csv_structure(loadcsv)
    if valid_csv:
        print("Call mspvalidator function expandpeaklist to deal with peaklist and intensity")
        peak_flattern_csv = validator.expandpeaklist(loadcsv)
    elif not valid_csv:
        print("probably the csv is truncated or not supported format")
        raise ValueError("the unlabel csv file is invalid for GlycoMSP!")
    else:
        raise Exception("The validation process failed for unknown reason")
    loadionlist = pd.read_csv(ionlist)
    flatternions = extract_ionmasslist(loadionlist)
    
    #complist = pd.read_csv(compositionlist)
    
    peak_flattern_csv.extend(['PseudoComposition'])
    print("haha")




#my requirements
"""
part A - Get peaklist and peakintensity while keep the index (peaklist[i] and peakintensity[i]) 
bound and the peaklist and peakintensity pair should still remember it's original location (row[k] has peaklist row[k_i] and peakintensity[k_i] ) 
so even if we use numpy to do the math, we can attach the results back to proper place - 
Probably need to get a ion list, suppose it's ion_filter_for_pseudolabel(mode, input=None, debug=True) 
while mode == "legacy" we will directly import an existing list and calculate if the peaklist contains the m/z value matching the list 
using the peaks_ppm from mspvalidator_merger.py - Matched ions peaklist paired to the intensity and a score should be appended to the original df 
(could be a copy to prevent warning of editing original file) 

part B - compare the protonatedmass in unlabeled csv file with loaded in_silico_composition using that binary search we've mentioned, 
also we may need to sort the mass value in in_silico csv file because it's sorted by number of hexoses when we generate the list 
- record all composition matching that row (equals to certain ms2 scan) - do the part A only if the row has composition 
(means probably a glycan spectra) Please also estimate the priority of these 2 task and decide which first to avoid waste of calculation. 
I still feel I should grant you the access because we already touched the old files. Let me know if you need it to work better.

Part A

Keep peaklist/peakintensity paired and row-aligned

There isn’t a special “alignment” helper, but all current code keeps them paired per row by design.

Parsing/unwrapping: expandpeaklist(df) converts the tuple-strings to Python sequences so you can index peaklist[i] ↔ peakintensity[i] safely (mspvalidator_merger.py).

Matching is done per row, so alignment is preserved.

Get an ion list & ppm-match to peaks

“ion_filter_for_pseudolabel(..)” does not exist.

Closest pieces you already have:

extract_fragment_masses(excel_path, sheet_name="ionlist") → list of target m/z from the annotation Excel (mspfileloaderv10.py; also in v09).

extract_ionmasslist(ionmass_sheet) → builds a header/index that includes masses (mspvalidator_merger.py).

peaks_ppm(inputpeaklist, peakintensity, reflist, setppm) and findingions(dfrow, iondf, ppm) → do the per-spectrum ppm matching (mspvalidator_merger.py).

Append matched-ion intensities/score back to data

You have a “wide” exporter: createnormailzedionlistcsv(ionindex, converted_df, ion_df, filename) which calls findingions(...), collects best hits per ion, log-scales intensities, and writes a normalized table (mspvalidator_merger.py).

There’s no single function that returns a long/tidy per-scan×ion frame or that merges a score back into the original DataFrame in-memory; you’d add a thin wrapper around findingions(...) for that.

Status summary (Part A):

Alignment handling: ✅ (implicit via row-wise logic; expandpeaklist exists)

Get ions (“legacy” import): ⚠️ partial (use extract_fragment_masses, no ion_filter_for_pseudolabel)

ppm matching to peaks: ✅ (peaks_ppm, findingions)

Append hits/score to df: ⚠️ partial (wide CSV via createnormailzedionlistcsv; tidy/merge helper missing)

Part B

Compare protonatedmass to in-silico compositions with binary search; keep all matches

No binary-search (or equivalent) composition matcher is present in the project files you shared. The only related code is in your old notebook (CompositionFinder2.ipynb).

None of the current modules implement “sorted theoretical mass array → (multi-match) lookup” against protonatedmass.

Record all compositions per row (per MS2 scan)

Not implemented in the shared files.

Run Part A only if a composition match exists

Not implemented; would be the orchestration glue around a new composition-match step.

Status summary (Part B):

Precursor↔composition match: ❌ missing

Collect all matches per scan: ❌ missing

Gate Part A on Part B hits: ❌ missing

What I recommend doing first (to save compute)

Do Part B first, then Part A on just the matched scans.

Complexity: Part A is ~O(N×M) (N scans × M ions). Part B is ~O(N log K + neighbors) (K theoretical comps). Filtering first usually shrinks N a lot (non-glycan spectra drop out), so you avoid a big fragment-matching pass on obvious non-hits.

Your own requirement says “do Part A only if the row has composition” — that’s exactly this order.

Minimal glue you’re missing (so you can slot it in)

ion_filter_for_pseudolabel(mode, input=None, debug=True) → add this (tiny):

mode=="legacy": if input is a path to CSV/XLSX, reuse extract_fragment_masses(...); if it’s a Python list, return it; else load your canned list.

match_precursor_to_compositions(unlabeled_df, insilico_df, ppm) → new:

Sort theoretical masses once; use np.searchsorted (or your old binary search) to collect all comps within ±ppm at each protonatedmass.

Return a long table: [MS2scan_no, observed_mass, theoretical_mass, composition, ppm_error] plus a grouped list per scan if you want.

attach_ion_hits(df_subset, ion_list, ppm) → new:

For the scan subset from Part B, call findingions(...) per row (or reuse extract_ion_intensities(...)), compute a score (e.g., sum of log-intensities of matched ions, or hit-count), and merge back on MS2scan_no.

If you want, I can draft these three helpers so you can drop them into mspvalidator_merger.py (or a new pseudolabel_utils.py) with your preferred column names.

"""




"""
# Write the requested Python module to a file so the user can download it.
from textwrap import dedent

code = dedent('''
    
    composition_finder.py
    ---------------------
    Utilities to (1) compare observed precursor masses in a converted MS2 CSV
    to theoretical composition masses, using ppm tolerance; and
    (2) (optional) extract tidy fragment-ion intensities from peaklists.

    Designed to be flexible with your current GlycoMSP merged/trainable CSV.
    If a function named `ppm` exists (e.g., provided elsewhere in your project),
    this module will use `ppm(20)` when you pass ppm_value=20. If not, it falls
    back to the standard tolerance: abs_delta = target_mass * ppm_value / 1e6.

    Output is long/tidy format ready for manual validation or ML.
    

    from __future__ import annotations

    import ast
    import json
    import math
    import os
    from typing import Callable, Iterable, List, Optional, Sequence, Tuple, Union

    import numpy as np
    import pandas as pd


    # -----------------------------
    # Helpers
    # -----------------------------

    def _maybe_call_ppm_function(ppm_value: float) -> float:
        \"\"\"
        If a global function named `ppm` exists, call it with `ppm_value`
        and use the returned number as the effective ppm. Otherwise return
        ppm_value unchanged.

        This matches your request: \"assume calling function ppm(20)\".
        \"\"\"
        try:
            # Look up ppm in global namespace (if user defines it elsewhere)
            fn = globals().get("ppm", None)
            if callable(fn):
                return float(fn(ppm_value))
        except Exception:
            pass
        return float(ppm_value)


    def _ppm_abs_tol(target_mass: float, ppm_value: float) -> float:
        \"\"\"Return the absolute m/z tolerance for a target mass at `ppm_value`.\"\"\"
        effective_ppm = _maybe_call_ppm_function(ppm_value)
        return abs(target_mass) * (effective_ppm / 1e6)


    def compute_ppm_error(observed: float, theoretical: float) -> float:
        \"\"\"Compute signed ppm error: (obs - theo) / theo * 1e6.\"\"\"
        if theoretical == 0:
            return np.nan
        return (float(observed) - float(theoretical)) / float(theoretical) * 1e6


    def _safe_parse_list(cell) -> List[float]:
        \"\"\"Parse peaklist/peakintensity cell into a list of floats (robust).\"\"\"
        if isinstance(cell, (list, tuple, np.ndarray)):
            return [float(x) for x in cell]
        if isinstance(cell, str):
            s = cell.strip()
            # Try JSON first
            try:
                val = json.loads(s)
                if isinstance(val, (list, tuple)):
                    return [float(x) for x in val]
            except Exception:
                pass
            # Try Python literal
            try:
                val = ast.literal_eval(s)
                if isinstance(val, (list, tuple)):
                    return [float(x) for x in val]
            except Exception:
                pass
        # Fallback: empty
        return []


    def _resolve_sheet_or_path(path: Union[str, os.PathLike], sheet: Optional[str] = None) -> pd.DataFrame:
        \"\"\"Load CSV or Excel into a DataFrame. If Excel, optionally choose sheet.\"\"\"
        path = str(path)
        if path.lower().endswith(('.xlsx', '.xls')):
            return pd.read_excel(path, sheet_name=sheet)
        return pd.read_csv(path)


    # -----------------------------
    # Core: precursor ↔ composition matching
    # -----------------------------

    def match_precursor_to_compositions(
        merged_csv: Union[str, os.PathLike, pd.DataFrame],
        compositions: Union[str, os.PathLike, pd.DataFrame],
        ppm_value: float = 20.0,
        observed_mass_col: Optional[str] = None,
        composition_mass_col: Optional[str] = None,
        composition_id_cols: Optional[Sequence[str]] = None,
        keep_top_n_per_scan: Optional[int] = None,
    ) -> pd.DataFrame:
        \"\"\"
        Compare observed precursor masses in a converted CSV to theoretical masses.

        Parameters
        ----------
        merged_csv
            Path to your merged/trainable CSV (or a DataFrame). Must contain
            at least one of these precursor mass columns (in this priority):
            ['protonatedmass', 'MS1_monoisolationmass', 'MS1_isolationmass'].
        compositions
            Path to a CSV/Excel (or DataFrame) with theoretical masses. Expected
            columns:
              • either a mass column like ['theoretical_mass', 'mass', 'mz']
              • OR component columns that you already pre-summed into a 'mass' column
                (this function does not try to calculate from monosaccharide counts).
        ppm_value
            PPM tolerance (e.g., 20). If a global `ppm` function exists, we call
            `ppm(ppm_value)` first and use that number as the effective ppm.
        observed_mass_col
            If provided, use this column as the observed precursor mass.
            Otherwise we auto-pick from the known candidates.
        composition_mass_col
            If provided, use this column as the theoretical mass column.
            Otherwise we auto-pick from ['theoretical_mass', 'mass', 'mz'].
        composition_id_cols
            Columns in the composition table that identify a composition (e.g. ['name','composition']).
            These will be carried through to the tidy output.
        keep_top_n_per_scan
            If set, keep only the top-N (lowest |ppm_error|) composition matches per scan.

        Returns
        -------
        pd.DataFrame
            Long-format tidy table with columns including:
            ['MS2scan_no','observed_mass','theoretical_mass','ppm_error','abs_ppm',
             'match_within_ppm', <composition_id_cols...>, <*other composition fields>]

        Notes
        -----
        • This function does not attempt adduct handling or derivatization math; it assumes
          your composition table already contains the correct theoretical (protonated) masses.
        • If needed later, we can add an optional callback to compute masses on the fly.
        \"\"\"
        # Load inputs
        if isinstance(merged_csv, (str, os.PathLike)):
            obs_df = pd.read_csv(merged_csv)
        else:
            obs_df = merged_csv.copy()

        if isinstance(compositions, (str, os.PathLike)):
            comp_df = _resolve_sheet_or_path(compositions)
        else:
            comp_df = compositions.copy()

        # Pick observed mass column
        if observed_mass_col is None:
            for cand in ['protonatedmass', 'MS1_monoisolationmass', 'MS1_isolationmass']:
                if cand in obs_df.columns:
                    observed_mass_col = cand
                    break
        if observed_mass_col is None:
            raise ValueError(\"No precursor mass column found. Please pass `observed_mass_col`.\" )

        # Pick composition mass column
        if composition_mass_col is None:
            for cand in ['theoretical_mass', 'mass', 'mz']:
                if cand in comp_df.columns:
                    composition_mass_col = cand
                    break
        if composition_mass_col is None:
            raise ValueError(\"No theoretical mass column found in composition table. Provide `composition_mass_col`.\")

        # Identify composition ID columns to retain
        if composition_id_cols is None:
            # sensible defaults if present
            auto_cols = [c for c in ['name', 'composition', 'glycan', 'id'] if c in comp_df.columns]
            composition_id_cols = auto_cols or []

        # Streamlined cross-match:
        # For efficiency, we vectorize by expanding compositions once and joining via a condition.
        obs_small = obs_df[['MS2scan_no', observed_mass_col]].dropna().copy()
        obs_small.rename(columns={observed_mass_col: 'observed_mass'}, inplace=True)

        comp_small = comp_df[composition_id_cols + [composition_mass_col]].copy() if composition_id_cols else comp_df[[composition_mass_col]].copy()
        comp_small.rename(columns={composition_mass_col: 'theoretical_mass'}, inplace=True)

        # Cartesian join (small → medium-sized tables expected). If very large, switch to interval indexing.
        obs_small['key'] = 1
        comp_small['key'] = 1
        cross = obs_small.merge(comp_small, on='key').drop(columns='key')

        # Compute ppm error and match flag
        cross['ppm_error'] = compute_ppm_error(cross['observed_mass'].values, cross['theoretical_mass'].values)
        cross['abs_ppm'] = cross['ppm_error'].abs()

        # Filter by tolerance (mass-dependent absolute window)
        # Use vectorized absolute tolerance per-theoretical mass
        eff_ppm = _maybe_call_ppm_function(ppm_value)
        abs_tol = cross['theoretical_mass'].abs() * (eff_ppm / 1e6)
        cross['match_within_ppm'] = (cross['observed_mass'] - cross['theoretical_mass']).abs() <= abs_tol

        matched = cross[cross['match_within_ppm']].copy()

        # Optional: keep top-N per scan
        if keep_top_n_per_scan is not None and keep_top_n_per_scan > 0:
            matched = (
                matched.sort_values(['MS2scan_no', 'abs_ppm'], ascending=[True, True])
                       .groupby('MS2scan_no', as_index=False)
                       .head(keep_top_n_per_scan)
            )

        # Tidy column order
        lead_cols = ['MS2scan_no', 'observed_mass', 'theoretical_mass', 'ppm_error', 'abs_ppm', 'match_within_ppm']
        id_cols = [c for c in composition_id_cols if c in matched.columns]
        other_cols = [c for c in matched.columns if c not in lead_cols + id_cols]
        ordered = matched[lead_cols + id_cols + other_cols]

        return ordered.reset_index(drop=True)


    # -----------------------------
    # Optional: fragment-ion feature extraction (tidy long format)
    # -----------------------------

    def match_fragment_ions(
        merged_csv: Union[str, os.PathLike, pd.DataFrame],
        ion_table: Union[str, os.PathLike, pd.DataFrame],
        ppm_value: float = 20.0,
        ion_mass_col: Optional[str] = None,
        ion_id_cols: Optional[Sequence[str]] = None,
    ) -> pd.DataFrame:
        \"\"\"
        Extract tidy fragment-ion intensities per spectrum by ppm matching.

        Parameters
        ----------
        merged_csv
            Path or DataFrame of your converted MS2 CSV containing columns:
            ['MS2scan_no', 'peaklist', 'peakintensity'] (strings or lists).
        ion_table
            Path or DataFrame containing fragment ion m/z values. It should have
            a mass column (auto-detected from ['mass','mz']). If additional
            identifying columns exist (e.g., 'ion', 'name', 'type'), include them
            via `ion_id_cols`.
        ppm_value
            PPM tolerance for peak ↔ ion matching.
        ion_mass_col
            Override the mass column name if auto-detection is wrong.
        ion_id_cols
            Columns to carry through for identification.

        Returns
        -------
        pd.DataFrame
            Long-format table with rows per (MS2scan_no × ion). Columns include:
            ['MS2scan_no','ion_mass','found','best_peak_mz','ppm_error','intensity_log10','intensity_raw', <ion_id_cols...>]
        \"\"\"
        # Load inputs
        if isinstance(merged_csv, (str, os.PathLike)):
            ms_df = pd.read_csv(merged_csv)
        else:
            ms_df = merged_csv.copy()

        if isinstance(ion_table, (str, os.PathLike)):
            ion_df = _resolve_sheet_or_path(ion_table)
        else:
            ion_df = ion_table.copy()

        # Detect columns
        if ion_mass_col is None:
            for cand in ['mass', 'mz', 'ion_mz']:
                if cand in ion_df.columns:
                    ion_mass_col = cand
                    break
        if ion_mass_col is None:
            raise ValueError(\"No ion mass column found in ion table. Provide `ion_mass_col`.\" )

        if ion_id_cols is None:
            ion_id_cols = [c for c in ['ion', 'name', 'label', 'type'] if c in ion_df.columns]

        ions = ion_df[ion_id_cols + [ion_mass_col]].copy() if ion_id_cols else ion_df[[ion_mass_col]].copy()
        ions.rename(columns={ion_mass_col: 'ion_mass'}, inplace=True)

        out_rows = []
        eff_ppm = _maybe_call_ppm_function(ppm_value)

        for _, row in ms_df.iterrows():
            scan = int(row['MS2scan_no']) if 'MS2scan_no' in row else None
            peaks = np.array(_safe_parse_list(row.get('peaklist', [])), dtype=float)
            intens = np.array(_safe_parse_list(row.get('peakintensity', [])), dtype=float)

            # Skip if no peaks
            if peaks.size == 0 or intens.size == 0 or len(peaks) != len(intens):
                for _, irow in ions.iterrows():
                    rec = {'MS2scan_no': scan, 'ion_mass': float(irow['ion_mass']), 'found': False,
                           'best_peak_mz': np.nan, 'ppm_error': np.nan, 'intensity_log10': 0.0, 'intensity_raw': 0.0}
                    for c in ion_id_cols:
                        rec[c] = irow[c]
                    out_rows.append(rec)
                continue

            for _, irow in ions.iterrows():
                target = float(irow['ion_mass'])
                abs_tol = abs(target) * (eff_ppm / 1e6)
                diffs = np.abs(peaks - target)
                mask = diffs <= abs_tol
                if np.any(mask):
                    # Choose the strongest peak in window
                    idx = np.argmax(intens[mask])
                    matched_peaks = peaks[mask]
                    matched_ints = intens[mask]
                    best_mz = float(matched_peaks[idx])
                    best_int = float(matched_ints[idx])
                    ppm_err = compute_ppm_error(best_mz, target)
                    rec = {'MS2scan_no': scan, 'ion_mass': target, 'found': True,
                           'best_peak_mz': best_mz, 'ppm_error': ppm_err,
                           'intensity_log10': float(np.log10(best_int + 1.0)), 'intensity_raw': best_int}
                else:
                    rec = {'MS2scan_no': scan, 'ion_mass': target, 'found': False,
                           'best_peak_mz': np.nan, 'ppm_error': np.nan,
                           'intensity_log10': 0.0, 'intensity_raw': 0.0}
                for c in ion_id_cols:
                    rec[c] = irow[c]
                out_rows.append(rec)

        return pd.DataFrame(out_rows)


    # -----------------------------
    # Convenience: file IO
    # -----------------------------

    def save_dataframe(df: pd.DataFrame, out_path: str) -> str:
        \"\"\"Save DataFrame as CSV (UTF-8) and return the path.\"\"\"
        out_path = str(out_path)
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
        df.to_csv(out_path, index=False, encoding='utf-8')
        return out_path


    # -----------------------------
    # CLI
    # -----------------------------

    def _cli():
        import argparse

        p = argparse.ArgumentParser(description=\"Compare observed precursor masses to theoretical compositions; optional ion features.\")
        p.add_argument('--merged', required=True, help='Path to merged/trainable CSV')
        p.add_argument('--compositions', required=True, help='Path to composition table (CSV/XLSX)')
        p.add_argument('--ppm', dest='ppm_value', type=float, default=20.0, help='PPM tolerance (default: 20)')
        p.add_argument('--out-matches', default='matches_long.csv', help='Output CSV for composition matches')
        p.add_argument('--ions', default=None, help='Optional ion/fragment table (CSV/XLSX)')
        p.add_argument('--out-ions', default='ion_features_long.csv', help='Output CSV for ion features (if --ions used)')
        p.add_argument('--topn', type=int, default=None, help='Keep top-N matches per scan (by |ppm|)')

        args = p.parse_args()

        matches = match_precursor_to_compositions(
            merged_csv=args.merged,
            compositions=args.compositions,
            ppm_value=args.ppm_value,
            keep_top_n_per_scan=args.topn,
        )
        save_dataframe(matches, args.out_matches)

        if args.ions:
            ion_feats = match_fragment_ions(args.merged, args.ions, ppm_value=args.ppm_value)
            save_dataframe(ion_feats, args.out_ions)

    if __name__ == '__main__':
        _cli()
''')

path = '/mnt/data/composition_finder.py'
with open(path, 'w', encoding='utf-8') as f:
    f.write(code)

path


"""