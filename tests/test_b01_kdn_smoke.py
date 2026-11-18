"""B-01 smoke test — verify the KDN token-convention fixes round-trip cleanly.

The bug (see handoff/kdnbug.md): four closure-local helpers inside `run_prediction`
plus the sister `_comp_tuple_to_label` outside it used three different KDN-token
conventions, breaking the precursor gate for any KDN-containing prediction.

The fix (S2 / 2026-05-17, B-01): unify all helpers on the 3-char `KDN` token,
with a back-compat alias on the canonicalizer that folds legacy bare `K` into
`KDN` (so MAS-trained models that emit `K1`-form predictions still resolve
against the in-silico library).

This test mirrors the post-fix source of each helper and exercises:
    T1 - _canon_comp canonical round-trip
    T2 - _canon_comp legacy K-alias to KDN
    T3 - _parse_counts no HexNAc inflation, no KDN drop
    T4 - _counts_to_comp emits 3-char KDN
    T5 - end-to-end gate trace (KDN-form AND legacy K-form) hits fake lib_counts
    T6 - _comp_tuple_to_label (sister site outside run_prediction)

If src/mspfileloaderv14.py is edited so this test diverges from source behavior,
update the mirrored helpers below in lockstep.

Run from repo root:
    python3 tests/test_b01_kdn_smoke.py
"""
import re
import sys


# ---- mirrors of the post-fix v14 helpers --------------------------------
# Keep in lockstep with src/mspfileloaderv14.py edits 1-5 (lines 7767-8011, 3811).


# A2: _canon_comp (post-Edit 1 — strict case, KDN-first, K aliased to KDN)
_comp_pat = re.compile(r'(KDN|K|F|H|N|S|G|A|s|p)\s*([0-9]+)')


def _canon_comp(s):
    if not isinstance(s, str) or not s.strip():
        return ""
    counts = {"F": 0, "H": 0, "N": 0, "S": 0, "G": 0, "KDN": 0,
              "A": 0, "s": 0, "p": 0}
    for k, v in _comp_pat.findall(s):
        if k == "K":
            k = "KDN"
        counts[k] = counts.get(k, 0) + int(v)
    parts = []
    for key in ("F", "H", "N", "S", "G", "KDN", "A", "s", "p"):
        n = counts.get(key, 0)
        if n > 0:
            parts.append(f"{key}{n}")
    return "".join(parts)


# A4: _counts_to_comp (post-Edit 2 — emits 3-char KDN)
def _counts_to_comp(row):
    F = row.get("Fuc", 0)
    H = row.get("Hex", 0)
    N = row.get("HexNAc", 0)
    NeuAc = row.get("NeuAc", 0)
    NeuGc = row.get("NeuGc", 0)
    KDN = row.get("KDN", 0)
    HexA = row.get("HexA", 0)
    SO3 = row.get("SO3", 0) + row.get("s", 0)
    PO3H = row.get("PO3H", 0) + row.get("p", 0)
    parts = []
    if F:     parts.append(f"F{F}")
    if H:     parts.append(f"H{H}")
    if N:     parts.append(f"N{N}")
    if NeuAc: parts.append(f"S{NeuAc}")
    if NeuGc: parts.append(f"G{NeuGc}")
    if KDN:   parts.append(f"KDN{KDN}")
    if HexA:  parts.append(f"A{HexA}")
    if SO3:   parts.append(f"s{SO3}")
    if PO3H:  parts.append(f"p{PO3H}")
    return "".join(parts)


# A5: _parse_counts (post-Edit 3 — KDN-aware regex, no bare K)
def _parse_counts(label):
    d = {"F": 0, "H": 0, "N": 0, "S": 0, "G": 0, "KDN": 0,
         "A": 0, "s": 0, "p": 0}
    for m in re.finditer(r'(KDN|[FHNSGAsp])(\d+)', str(label)):
        ch, num = m.group(1), int(m.group(2))
        if ch in d:
            d[ch] += num
    return d


# A6: _masses_by_counts (post-Edit 4 — filter key K -> KDN)
def _masses_by_counts(lib_counts_rows, lbl):
    """lib_counts_rows is a list of dicts mirroring the in-silico library rows."""
    cnt = _parse_counts(lbl)
    matches = []
    for r in lib_counts_rows:
        if (r["Fuc"] == cnt["F"]
                and r["Hex"] == cnt["H"]
                and r["HexNAc"] == cnt["N"]
                and r["NeuAc"] == cnt["S"]
                and r["NeuGc"] == cnt["G"]
                and r["KDN"] == cnt["KDN"]
                and r["HexA"] == cnt["A"]
                and r["SO3"] == cnt["s"]
                and r["PO3H"] == cnt["p"]):
            matches.append(r["__Mass__"])
    return sorted(set(matches))


# B1: _comp_tuple_to_label (post-Edit 5 — emits 3-char KDN)
def _comp_tuple_to_label(x):
    if x is None:
        return None
    if not isinstance(x, (list, tuple)) or len(x) < 6:
        return None
    hex_, hexc, neuac, neugc, kdn, fuc = [int(v) for v in x[:6]]
    parts = []
    if fuc:   parts.append(f"F{fuc}")
    if hex_:  parts.append(f"H{hex_}")
    if hexc:  parts.append(f"N{hexc}")
    if neuac: parts.append(f"S{neuac}")
    if neugc: parts.append(f"G{neugc}")
    if kdn:   parts.append(f"KDN{kdn}")
    return "".join(parts) if parts else None


# ---- tests --------------------------------------------------------------


def main():
    failures = []

    # T1 - canonical round-trip
    got = _canon_comp("F1H4N3KDN1")
    expected = "F1H4N3KDN1"
    if got != expected:
        failures.append(f"T1: _canon_comp('F1H4N3KDN1') == {got!r}, expected {expected!r}")
    else:
        print("T1 PASS: _canon_comp('F1H4N3KDN1') round-trips to canonical KDN form")

    # T2 - legacy K-alias to KDN (concern 4: MAS-trained model back-compat)
    got = _canon_comp("H4N3K1")
    expected = "H4N3KDN1"
    if got != expected:
        failures.append(f"T2: _canon_comp('H4N3K1') == {got!r}, expected {expected!r}")
    else:
        print("T2 PASS: _canon_comp('H4N3K1') aliases legacy K to KDN")

    # T3 - _parse_counts no HexNAc inflation, no KDN drop (failure mode 2)
    got = _parse_counts("H4N3KDN1")
    expected = {"F": 0, "H": 4, "N": 3, "S": 0, "G": 0, "KDN": 1,
                "A": 0, "s": 0, "p": 0}
    if got != expected:
        failures.append(f"T3: _parse_counts('H4N3KDN1') == {got}, expected {expected}")
    else:
        print("T3 PASS: _parse_counts('H4N3KDN1') correctly assigns N=3, KDN=1 (no inflation)")

    # T4 - _counts_to_comp emits 3-char KDN
    row = {"Fuc": 0, "Hex": 4, "HexNAc": 3, "NeuAc": 0, "NeuGc": 0,
           "KDN": 1, "HexA": 0, "SO3": 0, "PO3H": 0}
    got = _counts_to_comp(row)
    expected = "H4N3KDN1"
    if got != expected:
        failures.append(f"T4: _counts_to_comp(row) == {got!r}, expected {expected!r}")
    else:
        print("T4 PASS: _counts_to_comp emits 3-char KDN token")

    # T5 - end-to-end gate trace: failure mode 3 fixed for BOTH label forms
    fake_lib = [
        {"Fuc": 0, "Hex": 4, "HexNAc": 3, "NeuAc": 0, "NeuGc": 0,
         "KDN": 1, "HexA": 0, "SO3": 0, "PO3H": 0, "__Mass__": 1234.5},
        # Decoy row WITHOUT KDN at same H4N3 — pre-fix code would pool them.
        {"Fuc": 0, "Hex": 4, "HexNAc": 3, "NeuAc": 0, "NeuGc": 0,
         "KDN": 0, "HexA": 0, "SO3": 0, "PO3H": 0, "__Mass__": 999.9},
    ]
    # T5a - prediction label in canonical 3-char form
    lbl_canon = _canon_comp("H4N3KDN1")
    masses_canon = _masses_by_counts(fake_lib, lbl_canon)
    if masses_canon != [1234.5]:
        failures.append(f"T5a: gate trace 'H4N3KDN1' -> {masses_canon}, expected [1234.5]")
    else:
        print("T5a PASS: gate trace 'H4N3KDN1' resolves to KDN library row only (no decoy contamination)")

    # T5b - prediction label in legacy 1-char K form (MAS-trained model output)
    lbl_legacy = _canon_comp("H4N3K1")
    masses_legacy = _masses_by_counts(fake_lib, lbl_legacy)
    if masses_legacy != [1234.5]:
        failures.append(f"T5b: gate trace 'H4N3K1' -> {masses_legacy}, expected [1234.5] (K-alias path)")
    else:
        print("T5b PASS: gate trace 'H4N3K1' resolves via K-alias to KDN library row")

    # T6 - sister site B1
    got = _comp_tuple_to_label((4, 3, 0, 0, 1, 0))
    expected = "H4N3KDN1"
    if got != expected:
        failures.append(f"T6: _comp_tuple_to_label((4,3,0,0,1,0)) == {got!r}, expected {expected!r}")
    else:
        print("T6 PASS: _comp_tuple_to_label emits 3-char KDN (B1 sister fix)")

    # Summary
    total = 7  # T1, T2, T3, T4, T5a, T5b, T6
    if failures:
        print()
        print(f"FAIL: {len(failures)}/{total} tests failed")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    else:
        print()
        print(f"PASS: {total}/{total} tests passed")
        sys.exit(0)


if __name__ == "__main__":
    main()
