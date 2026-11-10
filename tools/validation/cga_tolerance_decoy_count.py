#!/usr/bin/env python3
"""CGA-TOLERANCE decoy count — chance-match baseline for widened CGA mass tolerances (standalone, counts only).

Method (handoff/CGA_TOLERANCE_decoy_method_20261031.md): match the same MS2 scans against a DECOY library made by
shifting every library mass by +3.37 Da (an arbitrary offset that avoids residue, isotope and charge spacings),
with the same tolerance. Every decoy hit is chance, so the decoy count estimates the chance matches among the
real-library ones. Optionally the same for fragments: the ion list is shifted by the same offset and Score A is
computed on the decoy ions.

Matching rules mirror GlycoMSP 1.14 (v14 run_pseudolabeling, CGA-TOLERANCE):
  precursor  |obs - theo| <= obs * ppm / 1e6   (ppm)   or   |obs - theo| <= Da   (Da)   — obs = `protonatedmass`;
             rows with no mass (0 / empty) are not scans with a precursor mass and are never matched.
  fragments  peaks_ppm: peak <= ion: 1000000 * (ion - peak) / ion <= ppm;  peak > ion: 1000000 * (peak - ion) / peak
             < ppm (the v1.10 asymmetry, known issue 10);  Da: ion - peak <= Da / peak - ion <= Da.  A peak is a hit
             only if log10(I) + 1 > 1 (I < 0 taken as 0), as GlycoMSP computes it.  Every final hit decision uses
             exactly these float operations; a binary search is used only to find candidates, with a margin.
  Score A    hits / ion-list size, capped at 1.0, computed per scan (it does not depend on the candidate).
             Two counts (B-39): "glycomsp" = one hit per in-window peak (what GlycoMSP writes today),
             "unique" = each ion at most once.
  pass       Score A >= --min-score (default 0.07), or hits >= --min-hits when given (to compare ion lists of
             different size at equal hit counts).
  irregular  a scan whose peak list is not strictly ascending (unsorted or a repeated m/z), or any scan outside the
             fast path's proven domain (ppm > 500,000, Da > 1,000, or a mass above 1e6 Da / not finite), is scored
             with a literal copy of GlycoMSP's loop (first-occurrence intensity, early stop), so the result equals
             GlycoMSP for every finite input; converter output is ascending and in range, so this path is normally
             unused (the irregular-list count is reported).

Output: counts only, under --alias. No input file name, path or scan list is printed or written.
Needs Python 3.9+, pandas, numpy (openpyxl only for an .xlsx ion list).

Example (Lab PC):
  python cga_tolerance_decoy_count.py --alias pilot --converted <converted.csv> --library <insilico.csv> \
      --ions <ionlist.csv> --precursor 20ppm,300ppm,0.5Da --fragment 10ppm,0.3Da
"""
import argparse
import csv
import json
import re
import sys

import numpy as np
import pandas as pd

DECOY_SHIFT_DA = 3.37
_TOL_RE = re.compile(r"^\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*(ppm|da)\s*$", re.IGNORECASE)


def parse_tolerances(text):
    """'20ppm,0.3Da' -> [(20.0, 'ppm'), (0.3, 'Da')]."""
    out = []
    for part in (text or "").split(","):
        if not part.strip():
            continue
        m = _TOL_RE.match(part)
        if not m or float(m.group(1)) <= 0:
            raise SafeError("bad tolerance (use e.g. 20ppm or 0.3Da)")
        out.append((float(m.group(1)), "ppm" if m.group(2).lower() == "ppm" else "Da"))
    return out


def tol_label(t):
    return f"{t[0]:g} {t[1]}"


# ---------------------------------------------------------------- readers (column names as GlycoMSP uses them)

def _parse_list(x):
    if isinstance(x, (list, tuple)):
        return [float(v) for v in x]
    s = str(x).strip().strip("()[]")
    return [float(v) for v in s.split(",") if v.strip()] if s else []


def read_converted(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if "protonatedmass" not in df.columns or len(df.columns) < 2:
        df = pd.read_csv(path, low_memory=False)
    for c in ("protonatedmass", "peaklist", "peakintensity"):
        if c not in df.columns:
            raise SafeError(f"converted CSV has no '{c}' column")
    mass = pd.to_numeric(df["protonatedmass"], errors="coerce").fillna(0.0).to_numpy(float)
    scans = []
    for m, pk, it in zip(mass, df["peaklist"], df["peakintensity"]):
        p, i = _parse_list(pk), _parse_list(it)
        if len(p) != len(i):
            p, i = [], []
        scans.append(make_scan(float(m), p, i))
    return scans


def _intensity_ok(i):
    """GlycoMSP's hit rule on one intensity: intensity_normalization (I < 0 -> 0), then log10(I) + 1 > 1."""
    i = 0.0 if i < 0 else float(i)
    with np.errstate(divide="ignore"):
        return bool(np.log10(i) + 1 > 1)


def make_scan(mass, peaks, intensities):
    """(mass, Spectrum). All peaks are kept (GlycoMSP's early stop can be triggered by a low-intensity peak)."""
    p = [float(x) for x in peaks]
    i = [float(x) for x in intensities]
    regular = not (len(p) > 1 and not all(b > a for a, b in zip(p, p[1:])))
    return (float(mass), (Spectrum if regular else Irregular)(p, i))


class Spectrum:
    """Strictly ascending peak list: candidates found by binary search, decided with GlycoMSP's operations."""
    __slots__ = ("peaks", "ints", "arr", "ok")

    def __init__(self, peaks, ints):
        self.peaks, self.ints = peaks, ints
        self.arr = np.asarray(peaks, dtype=float)
        self.ok = [_intensity_ok(x) for x in ints]


class Irregular:
    """Peak list that is not strictly ascending (unsorted or a repeated m/z): the literal GlycoMSP loop."""
    __slots__ = ("peaks", "ints")

    def __init__(self, peaks, ints):
        self.peaks, self.ints = peaks, ints


def read_library(path):
    lib = pd.read_csv(path)
    col = next((c for c in lib.columns if c.strip().lower() == "mass"), None)
    if col is None:
        raise SafeError("library has no 'Mass' column")
    return np.sort(pd.to_numeric(lib[col], errors="coerce").dropna().to_numpy(float))


def read_ions(path):
    p = path.lower()
    if p.endswith((".xlsx", ".xls")):
        sheets = pd.read_excel(path, sheet_name=None)
        df = next((d for n, d in sheets.items() if n.strip().lower() in {"ionlist", "ions", "ion_list"}), None)
        if df is None:
            df = next((d for d in sheets.values()
                       if {c.strip().lower() for c in d.columns} & {"mass", "mz", "ion_mz", "m/z"}), None)
        if df is None:
            raise SafeError("ion workbook has no ion-list sheet with a mass column")
    else:
        df = pd.read_csv(path)
    col = next((c for c in df.columns if c.strip().lower() in ("mass", "mz", "ion_mz", "m/z")), None)
    if col is None:
        raise SafeError("ion list has no mass column")
    return pd.to_numeric(df[col], errors="coerce").dropna().to_numpy(float)


# ---------------------------------------------------------------- matching

def precursor_window(obs, tol):
    v, u = tol
    return v if u == "Da" else abs(obs) * v / 1e6


def count_candidates(obs, lib_sorted, tol):
    w = precursor_window(obs, tol)
    lo = np.searchsorted(lib_sorted, obs - w, side="left")
    hi = np.searchsorted(lib_sorted, obs + w, side="right")
    return int(hi - lo)


def _in_window(pk, ref, tol):
    """GlycoMSP's window test, with its exact float operations (peaks_ppm / _peaks_da)."""
    v, u = tol
    if pk <= ref:
        return (ref - pk <= v) if u == "Da" else (1000000 * (ref - pk) / ref <= v)
    return (pk - ref <= v) if u == "Da" else (1000000 * (pk - ref) / pk < v)


# Fast path (binary search for candidates) only where its candidate range provably contains every peak GlycoMSP
# would count (CP-C.3 R3-F3): ppm <= 500,000 (1 - f >= 0.5, so ref / (1 - f) is well conditioned), Da <= 1,000, all
# ion and peak masses finite and <= 1e6 Da. There the float error of the bounds is below ~1e-9 Da, far inside the
# 1e-6 Da margin. A ppm window also needs every ion mass > 0 (the bounds ion * (1 - f) .. ion / (1 - f) assume a
# positive reference; Codex CGAT-R4-F1). Anything else (and every unsorted / repeated-m/z list) runs GlycoMSP's literal
# loop, so the count equals GlycoMSP for every finite input, including GlycoMSP's own result (or error) for a zero or
# negative ion mass; the helper reproduces such inputs, it does not reject them.
_CANDIDATE_MARGIN = 1e-6
_FAST_MAX_PPM, _FAST_MAX_DA, _FAST_MAX_MASS = 5e5, 1e3, 1e6


def _fast_ok(spec, ions, tol):
    v, u = tol
    if not isinstance(spec, Spectrum):
        return False
    if (u == "ppm" and not v <= _FAST_MAX_PPM) or (u == "Da" and not v <= _FAST_MAX_DA):
        return False
    if ions.size and not (np.all(np.isfinite(ions)) and np.max(np.abs(ions)) <= _FAST_MAX_MASS):
        return False
    if u == "ppm" and ions.size and not np.all(ions > 0):
        return False
    arr = spec.arr
    return not arr.size or (np.all(np.isfinite(arr)) and np.max(np.abs(arr)) <= _FAST_MAX_MASS)


def _literal_hit_indices(peaklist, ints, ions, tol):
    """Literal copy of the mspvalidator_merger.peaks_ppm / _peaks_da loop: per ion, the positions of the peaks GlycoMSP
    counts as hits (window test with its float operations; intensity of the FIRST peak with the same m/z via
    peaklist.index; stop at the first peak above the ion that is outside the window)."""
    out = []
    for ref in (float(x) for x in ions):
        hits = []
        for k, pk in enumerate(peaklist):
            inwin = _in_window(pk, ref, tol)
            if pk > ref and not inwin:
                break
            if inwin and _intensity_ok(ints[peaklist.index(pk)]):
                hits.append(k)
        out.append(hits)
    return out


def _fast_hit_indices(spec, ions, tol):
    """Same result as _literal_hit_indices for a Spectrum inside the fast domain: the literal loop restricted to the
    candidate range (peaks outside it are outside the window, so skipping them changes nothing)."""
    peaks = spec.arr
    if peaks.size == 0 or ions.size == 0:
        return [[] for _ in range(ions.size)]
    v, u = tol
    if u == "Da":
        lo_b, hi_b = ions - v, ions + v
    else:
        f = v / 1e6
        lo_b, hi_b = ions * (1.0 - f), ions / (1.0 - f)
    lo = np.searchsorted(peaks, lo_b - _CANDIDATE_MARGIN, side="left")
    hi = np.searchsorted(peaks, hi_b + _CANDIDATE_MARGIN, side="right")
    pl, ok = spec.peaks, spec.ok
    out = [[] for _ in range(ions.size)]
    for j in np.nonzero(hi > lo)[0]:
        ref = float(ions[j])
        for k in range(int(lo[j]), int(hi[j])):
            pk = pl[k]
            inwin = _in_window(pk, ref, tol)
            if pk > ref and not inwin:
                break
            if inwin and ok[k]:
                out[j].append(k)
    return out


def fragment_hit_indices(spec, ions, tol):
    """Per ion, the positions of the peaks GlycoMSP counts (shared by the decoy counts and B39_shared_peaks.py)."""
    ions = np.asarray(ions, dtype=float)
    if _fast_ok(spec, ions, tol):
        return _fast_hit_indices(spec, ions, tol)
    return _literal_hit_indices(spec.peaks, spec.ints, ions, tol)


def fragment_hits(spec, ions, tol):
    """Per ion: GlycoMSP hit count (one hit per in-window peak, B-39)."""
    return np.asarray([len(h) for h in fragment_hit_indices(spec, ions, tol)], dtype=int)


def score_a(peaks, ions, tol):
    """(glycomsp_hits, unique_hits) for one scan."""
    per_ion = fragment_hits(peaks, ions, tol)
    return int(per_ion.sum()), int((per_ion > 0).sum())


# ---------------------------------------------------------------- driver

def run(scans, lib, ions, prec_tols, frag_tols, *, shift=DECOY_SHIFT_DA, top_n=3, min_score=0.07, min_hits=None):
    lib_d = lib + shift
    ions_d = None if ions is None else ions + shift
    with_mass = [(m, p) for m, p in scans if m > 0]
    res = {"n_rows": len(scans), "n_scans_with_mass": len(with_mass), "decoy_shift_da": shift,
           "top_n": top_n, "precursor": [], "fragments": [],
           "n_irregular_peak_lists": sum(1 for _, p in scans if isinstance(p, Irregular))}
    n_ions = 0 if ions is None else int(ions.size)
    res["ion_list_size"] = n_ions
    rule = (f"hits >= {min_hits}" if min_hits is not None else f"Score A >= {min_score:g}")
    res["pass_rule"] = rule

    def passes(h):
        if min_hits is not None:
            return h >= min_hits
        return n_ions > 0 and min(1.0, h / n_ions) >= min_score

    for pt in prec_tols:
        real_c = np.array([count_candidates(m, lib, pt) for m, _ in with_mass], dtype=int)
        dec_c = np.array([count_candidates(m, lib_d, pt) for m, _ in with_mass], dtype=int)
        res["precursor"].append({
            "precursor": tol_label(pt),
            "real_scans": int((real_c > 0).sum()), "decoy_scans": int((dec_c > 0).sum()),
            "real_pairs_topN": int(np.minimum(real_c, top_n).sum()) if top_n else int(real_c.sum()),
            "decoy_pairs_topN": int(np.minimum(dec_c, top_n).sum()) if top_n else int(dec_c.sum()),
        })
        if ions is None:
            continue
        for ft in frag_tols:
            row = {"precursor": tol_label(pt), "fragment": tol_label(ft)}
            for rule_name, idx in (("glycomsp", 0), ("unique", 1)):
                cnt = {"real_prec+real_ions": 0, "decoy_prec+decoy_ions": 0,
                       "real_prec+decoy_ions": 0, "decoy_prec+real_ions": 0}
                for (m, peaks), rc, dc in zip(with_mass, real_c, dec_c):
                    if rc == 0 and dc == 0:
                        continue
                    pr = passes(score_a(peaks, ions, ft)[idx])
                    pd_ = passes(score_a(peaks, ions_d, ft)[idx])
                    cnt["real_prec+real_ions"] += int(rc > 0 and pr)
                    cnt["decoy_prec+decoy_ions"] += int(dc > 0 and pd_)
                    cnt["real_prec+decoy_ions"] += int(rc > 0 and pd_)
                    cnt["decoy_prec+real_ions"] += int(dc > 0 and pr)
                row[rule_name] = cnt
            res["fragments"].append(row)
    return res


def format_report(alias, res):
    L = [f"alias: {alias}",
         f"rows: {res['n_rows']}   scans with a precursor mass: {res['n_scans_with_mass']}   "
         f"decoy shift: +{res['decoy_shift_da']:g} Da   candidates capped at top {res['top_n']} per scan   "
         f"irregular peak lists (exact GlycoMSP loop): {res['n_irregular_peak_lists']}",
         "", "precursor only (scans with >= 1 candidate)",
         f"{'precursor':>12} {'real':>8} {'decoy':>8} {'decoy/real':>10} {'real pairs':>11} {'decoy pairs':>12}"]
    for r in res["precursor"]:
        ratio = f"{r['decoy_scans'] / r['real_scans']:.3f}" if r["real_scans"] else "-"
        L.append(f"{r['precursor']:>12} {r['real_scans']:>8} {r['decoy_scans']:>8} {ratio:>10} "
                 f"{r['real_pairs_topN']:>11} {r['decoy_pairs_topN']:>12}")
    if res["fragments"]:
        L += ["", f"precursor + fragments (ion list: {res['ion_list_size']} ions; pass: {res['pass_rule']}; "
                  "Score A per scan)",
              f"{'precursor':>12} {'fragment':>10} {'count':>9} {'real+real':>10} {'decoy+decoy':>12} "
              f"{'real+decoyIons':>15} {'decoy+realIons':>15}"]
        for r in res["fragments"]:
            for rule in ("glycomsp", "unique"):
                c = r[rule]
                L.append(f"{r['precursor']:>12} {r['fragment']:>10} {rule:>9} {c['real_prec+real_ions']:>10} "
                         f"{c['decoy_prec+decoy_ions']:>12} {c['real_prec+decoy_ions']:>15} "
                         f"{c['decoy_prec+real_ions']:>15}")
    return "\n".join(L)


def write_csv(path, alias, res):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["alias", "section", "precursor", "fragment", "count_rule", "pass_rule", "metric", "value"])
        for r in res["precursor"]:
            for k in ("real_scans", "decoy_scans", "real_pairs_topN", "decoy_pairs_topN"):
                w.writerow([alias, "precursor", r["precursor"], "", "", "", k, r[k]])
        for r in res["fragments"]:
            for rule in ("glycomsp", "unique"):
                for k, v in r[rule].items():
                    w.writerow([alias, "fragments", r["precursor"], r["fragment"], rule, res["pass_rule"], k, v])


class SafeError(ValueError):
    """A message built by this script from fixed text only (never from input content): safe to print."""


def _err(e):
    """Never echo input content or paths (CGAT-F7): only SafeError text, else the exception type."""
    return str(e) if isinstance(e, SafeError) else type(e).__name__


class _QuietParser(argparse.ArgumentParser):
    """argparse echoes offending arguments (which can be file names); print a fixed message instead."""
    def error(self, message):
        self.print_usage(sys.stderr)
        self.exit(2, "error: invalid or missing options (see --help); arguments are not echoed\n")


def main(argv=None):
    ap = _QuietParser(description="CGA decoy count (counts only; see module docstring).")
    ap.add_argument("--alias", required=True, help="name printed instead of any file name")
    ap.add_argument("--converted", required=True, help="converted MS2 CSV (tab-separated)")
    ap.add_argument("--library", required=True, help="in silico library CSV (column 'Mass')")
    ap.add_argument("--ions", help="ion list (CSV or XLSX with a 'mass' column); enables the fragment section")
    ap.add_argument("--precursor", default="20ppm", help="comma list, e.g. 20ppm,300ppm,0.5Da")
    ap.add_argument("--fragment", default="10ppm", help="comma list, e.g. 10ppm,0.3Da (with --ions)")
    ap.add_argument("--min-score", type=float, default=0.07)
    ap.add_argument("--min-hits", type=int, default=None, help="pass on hit count instead of Score A")
    ap.add_argument("--top-n", type=int, default=3, help="candidate cap per scan for the pair counts (CGA: 3)")
    ap.add_argument("--shift", type=float, default=DECOY_SHIFT_DA, help="decoy offset in Da (default 3.37)")
    ap.add_argument("--out", help="also write the counts as CSV (alias-labelled, no paths)")
    ap.add_argument("--json", action="store_true", help="print the counts as JSON instead of the table")
    a = ap.parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9_.\-]{1,64}", a.alias):
        sys.exit("error: --alias may use letters, digits, '_', '-', '.' only (max 64)")
    try:
        prec, frag = parse_tolerances(a.precursor), parse_tolerances(a.fragment)
    except Exception as e:
        sys.exit(f"[{a.alias}] {_err(e)}")
    try:
        scans = read_converted(a.converted)
    except Exception as e:
        sys.exit(f"[{a.alias}] cannot read the converted CSV ({_err(e)})")
    try:
        lib = read_library(a.library)
    except Exception as e:
        sys.exit(f"[{a.alias}] cannot read the library ({_err(e)})")
    ions = None
    if a.ions:
        try:
            ions = read_ions(a.ions)
        except Exception as e:
            sys.exit(f"[{a.alias}] cannot read the ion list ({_err(e)})")
    try:
        res = run(scans, lib, ions, prec, frag, shift=a.shift, top_n=a.top_n, min_score=a.min_score,
                  min_hits=a.min_hits)
    except Exception as e:
        sys.exit(f"[{a.alias}] counting failed ({_err(e)})")
    res["alias"] = a.alias
    print(json.dumps(res, indent=2) if a.json else format_report(a.alias, res))
    if a.out:
        try:
            write_csv(a.out, a.alias, res)
        except Exception as e:
            sys.exit(f"[{a.alias}] cannot write the --out file ({_err(e)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
