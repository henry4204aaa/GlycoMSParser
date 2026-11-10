"""B-39 case 2 check: (1) ion pairs in an ion list closer than twice the fragment window (a property of the ion list
only: such ions CAN share one peak); (2) with a converted CSV, how many scans have a peak that GlycoMSP counts as a
hit for two or more ions, and how many such peaks.

Part (2) uses the same exact matching as cga_tolerance_decoy_count.py (fragment_hit_indices: GlycoMSP's window
float operations, intensity rule log10(I) + 1 > 1, first-occurrence intensity for a repeated m/z, early stop), so
its counts equal what GlycoMSP's peaks_ppm / _peaks_da would count (CP-C.3, Codex CGAT-R3-F4).
Counts only, no file names in the output. Keep next to cga_tolerance_decoy_count.py (it reuses its readers).

usage: python B39_shared_peaks.py <ion list>[,<ion list>...] [<converted MS2 CSV>]
"""
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cga_tolerance_decoy_count as D  # noqa: E402

TOLS = ((10.0, "ppm"), (0.3, "Da"), (0.5, "Da"))


def close_ion_pairs(ions, tol):
    """Ion pairs closer than twice the window (ion-list property; ppm window taken at the heavier ion)."""
    s = np.sort(np.asarray(ions, dtype=float))
    out = 0
    for i in range(len(s)):
        for j in range(i + 1, len(s)):
            gap = s[j] - s[i]
            lim = 2 * tol[0] if tol[1] == "Da" else 2 * s[j] * tol[0] / 1e6
            if gap > lim:
                break
            out += 1
    return out


def shared_in_scan(spec, ions, tol):
    """Number of peak positions counted as a hit for two or more ions in one scan."""
    per_ion = D.fragment_hit_indices(spec, ions, tol)
    c = Counter(k for hits in per_ion for k in set(hits))
    return sum(1 for n in c.values() if n >= 2)


def shared(scans, ions, tol):
    """(scans with >= 1 shared peak, shared peaks in total)."""
    n_sc = n_pk = 0
    for _, spec in scans:
        k = shared_in_scan(spec, ions, tol)
        n_pk += k
        n_sc += int(k > 0)
    return n_sc, n_pk


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        sys.exit(__doc__)
    ion_files = argv[0].split(",")
    try:
        scans = D.read_converted(argv[1]) if len(argv) > 1 else None
    except Exception as e:
        sys.exit(f"cannot read the converted CSV ({D._err(e)})")
    for idx, f in enumerate(ion_files, 1):
        try:
            ions = D.read_ions(f)
        except Exception as e:
            sys.exit(f"cannot read ion list {idx} ({D._err(e)})")
        print(f"ion list {idx}: {ions.size} ions, {np.unique(ions).size} unique masses")
        for tol in TOLS:
            line = f"   {D.tol_label(tol):>7}: ion pairs closer than 2x window = {close_ion_pairs(ions, tol)}"
            if scans is not None:
                sc, pk = shared(scans, ions, tol)
                line += f" | scans with a peak counted for >=2 ions = {sc} of {len(scans)} ({pk} peaks)"
            print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
