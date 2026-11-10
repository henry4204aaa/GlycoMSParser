"""GlycoMSP precursor mass CLI (V14-FOLLOWUP CP-A.2) - how GlycoMSP forms a glycan precursor mass.

version 0.3, 20261027. Every number comes from compnewv4.precursormassv4, the function the CGA in-silico library
uses, so this tool cannot drift from CGA: it has no mass table of its own. stdlib + compnewv4 only (no GUI import).

  python src/msp_CLI_precursor_mass.py --composition F1H5N4S2 [--deri perme|native] [--reduced] [--mode +|-]
         [--hexa N] [--so3 N] [--po3h N] [--json]
  python src/msp_CLI_precursor_mass.py --table OUT.csv|OUT.md [--composition C ...] [--deri ...] [--reduced] [--mode ...]
  python src/msp_CLI_precursor_mass.py --selftest

Composition labels use the CGA short form: F (Fuc), H (Hex), N (HexNAc), S (NeuAc), G (NeuGc), KDN (or K), A (HexA),
s (SO3), p (PO3H), each followed by a count, in any order (H5N4F1S2 = F1H5N4S2).

Masses: precursormassv4 returns an m/z. Native, and PerMe without SO3/PO3H, give the singly charged ion ([M+H]+ for
mode +, [M-H]- for mode -); the proton term is symmetric, so the neutral mass printed here is the midpoint of the
two ions. For PerMe with SO3/PO3H the function applies fixed charge-state corrections, so the value is labelled
"m/z (z=+1)" / "m/z (z=-1)", or [M+2H]2+ / [M-2H]2- for two SO3 (the function halves the mass), and no neutral
mass is derived. Charge rule (Henry 2026-10-27): each SO3 / PO3H carries one charge; two PO3H is the observed
species that loses one by elimination, so it is singly charged. compnewv4 models only one group type at a time
(SO3 1-2 or PO3H 1-2), so PerMe with both SO3 and PO3H, or more than two of either, is refused (exit 2).

--table without --composition writes the CP-A.1 acceptance rows; with --composition it writes those compositions
for the chosen options, --hexa/--so3/--po3h included (both polarities unless --mode is given). The format follows
the extension (.csv or .md). --selftest checks the acceptance rows: PerMe rows must equal the base output exactly,
Native rows must be within 0.002 Da of their reference; a row without a reference fails. Exit codes: 0 ok,
1 selftest mismatch, 2 usage error.
"""
import argparse
import csv
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import compnewv4  # noqa: E402

version = "0.3"   # 0.3 (Codex V14F-R2, Henry): PerMe mixed SO3 + PO3H refused
# 0.2 (Codex V14F-F3/F4/F8): table group flags, PerMe SO3/PO3H labels, rows 7-9 references
last_update = 20261027

TOL_NATIVE = 0.002  # Da, kickoff pass rule
_TOKEN = re.compile(r"(KDN|[FHNSGKAsp])(\d+)")
_KEYS = {"F": "Fuc", "H": "Hex", "N": "HexNAc", "S": "NeuAc", "G": "NeuGc", "KDN": "KDN", "K": "KDN",
         "A": "HexA", "s": "SO3", "p": "PO3H"}
_DERI = {"perme": "PerMe", "native": "Native"}
_ION = {"+": "[M+H]+", "-": "[M-H]-"}


def parse_composition(label):
    """'H5N4F1S2' -> dict of counts (Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc, HexA, SO3, PO3H). Raises ValueError."""
    text = (label or "").strip()
    if not text:
        raise ValueError("empty composition")
    counts = dict.fromkeys(("Hex", "HexNAc", "NeuAc", "NeuGc", "KDN", "Fuc", "HexA", "SO3", "PO3H"), 0)
    seen, pos = set(), 0
    for m in _TOKEN.finditer(text):
        if m.start() != pos:
            break
        key = _KEYS[m.group(1)]
        if key in seen:
            raise ValueError(f"{key} given twice in {label!r}")
        seen.add(key)
        counts[key] = int(m.group(2))
        pos = m.end()
    if pos != len(text):
        raise ValueError(f"cannot read {label!r} at {text[pos:]!r} (letters: F H N S G KDN A s p, each with a count)")
    return counts


def label_of(counts):
    """Canonical short label, CGA order: F, H, N, S, G, KDN, A, s, p (zero counts left out)."""
    order = (("F", "Fuc"), ("H", "Hex"), ("N", "HexNAc"), ("S", "NeuAc"), ("G", "NeuGc"), ("KDN", "KDN"),
             ("A", "HexA"), ("s", "SO3"), ("p", "PO3H"))
    return "".join(f"{tag}{counts[key]}" for tag, key in order if counts[key]) or "(empty)"


PERME_GROUP_MAX = 2   # compnewv4 has PerMe SO3/PO3H rules for counts 1 and 2 only


def ion_label(counts, deri, mode):
    """What precursormassv4 returns for these inputs (Codex V14F-F4, Henry 2026-10-27).

    Native and plain PerMe: the singly charged ion. PerMe with SO3/PO3H: an m/z with fixed charge-state corrections,
    singly charged ("m/z (z=+1)" / "m/z (z=-1)"): one SO3, one PO3H, or two PO3H (observed with one charge after
    an elimination, Henry 2026-10-27). Two SO3, where the function halves the mass: [M-2H]2- / [M+2H]2+. Mixed
    SO3 + PO3H never reaches here (check_supported).
    """
    if deri == "PerMe" and (counts["SO3"] or counts["PO3H"]):
        if counts["SO3"] == 2:
            return "[M-2H]2-" if mode == "-" else "[M+2H]2+"
        return "m/z (z=-1)" if mode == "-" else "m/z (z=+1)"
    return _ION[mode]


def check_supported(counts, deri):
    """Raise ValueError for PerMe inputs whose charge state compnewv4 does not model: more than two SO3 or PO3H, or
    SO3 and PO3H together (additive charge 2 or more; Henry 2026-10-27: not observed, refused)."""
    if deri == "PerMe":
        if counts["SO3"] and counts["PO3H"]:
            raise ValueError("PerMe with both SO3 and PO3H is not supported: compnewv4 has no rule for the combined "
                             "charge state")
        for key in ("SO3", "PO3H"):
            if counts[key] > PERME_GROUP_MAX:
                raise ValueError(f"PerMe with {counts[key]} {key} is not supported: compnewv4 has PerMe {key} "
                                 f"rules for 1 and 2 only")


def masses(counts, deri="PerMe", reduced=False):
    """The m/z for both polarities and the neutral mass, all from compnewv4.precursormassv4."""
    comp = tuple(counts[k] for k in ("Hex", "HexNAc", "NeuAc", "NeuGc", "KDN", "Fuc"))
    kw = dict(hexA=counts["HexA"], so3=counts["SO3"], po3h=counts["PO3H"])
    ions = {mode: compnewv4.precursormassv4(comp, deri=deri, reduced=reduced, mode=mode, **kw) for mode in "+-"}
    has_groups = counts["SO3"] or counts["PO3H"]
    neutral = None if (deri == "PerMe" and has_groups) else (ions["+"] + ions["-"]) / 2
    return neutral, ions


# CP-A.1 acceptance rows: (label, deri, reduced, mode, reference neutral, reference ion, source).
# PerMe references are the base output at 53f3b9b and must match exactly (repr). Native references are the
# kickoff's figures, residue constants checked against GlycoMod by Henry 2026-10-27. Rows 7-9: the values Henry
# confirmed against GlycoMod and GlycoWorkBench 2026-10-27 (row 8 neutral = GlycoMod 1947.6824546, ion = that
# + 1.00728); tests/test_cli_precursor_mass.py also recomputes them from elemental formulas.
ACCEPTANCE_ROWS = [
    ("H5N4", "PerMe", False, "+", None, 2048.0532900000003, "base output 53f3b9b (kickoff 2048.053)"),
    ("H5N4", "PerMe", True, "+", None, 2064.0845900000004, "base output 53f3b9b"),
    ("H5N4", "Native", False, "+", 1640.592, 1641.599, "kickoff / GlycoMod"),
    ("H5N4", "Native", True, "+", 1642.608, 1643.615, "kickoff / GlycoMod"),
    ("F1H5N4S2", "Native", False, "+", 2368.841, 2369.848, "kickoff / GlycoMod"),
    ("H3N2", "Native", False, "+", 910.328, 911.335, "kickoff / GlycoMod"),
    ("H5N4KDN1", "Native", False, "+", 1890.6610, 1891.6683, "Henry 2026-10-27: GlycoMod / GlycoWorkBench"),
    ("H5N4G1", "Native", False, "+", 1947.6825, 1948.6897, "Henry 2026-10-27: GlycoMod 1947.6824546 / GlycoWorkBench"),
    ("H5N4A1", "Native", False, "+", 1816.6242, 1817.6315, "Henry 2026-10-27: GlycoMod / GlycoWorkBench"),
    ("H5N4", "Native", False, "-", 1640.592, 1639.585, "kickoff / GlycoMod"),
    ("H5N4", "Native", True, "-", 1642.608, 1641.6005, "kickoff / GlycoMod (pilot conditions)"),
    ("F1H5N4S2", "Native", True, "-", 2370.857, 2369.8493, "kickoff / GlycoMod (pilot conditions)"),
]


def evaluate(row):
    """Compute one acceptance row through precursormassv4 and judge it against its reference."""
    label, deri, reduced, mode, ref_neutral, ref_ion, source = row
    counts = parse_composition(label)
    neutral, ions = masses(counts, deri, reduced)
    ion = ions[mode]
    if deri == "PerMe":
        ok = ref_ion is not None and repr(ion) == repr(ref_ion)
        result = "identical" if ok else "MISMATCH"
    elif ref_ion is None:
        ok, result = False, "NO REFERENCE"
    else:
        ok = abs(ion - ref_ion) <= TOL_NATIVE and (ref_neutral is None or abs(neutral - ref_neutral) <= TOL_NATIVE)
        result = "within 0.002 Da" if ok else "MISMATCH"
    return {"composition": label, "derivatization": deri, "reduced": reduced, "mode": mode,
            "neutral_mass": neutral, "ion": ion_label(counts, deri, mode), "ion_mz": ion,
            "reference_neutral": ref_neutral, "reference_ion": ref_ion,
            "delta_ion": None if ref_ion is None else ion - ref_ion,
            "source": source, "result": result, "ok": ok}


def custom_rows(compositions, deri, reduced, modes):
    """Table rows for count dicts that already went through parse_composition and the group flags."""
    rows = []
    for counts in compositions:
        neutral, ions = masses(counts, deri, reduced)
        for mode in modes:
            rows.append({"composition": label_of(counts), "derivatization": deri, "reduced": reduced, "mode": mode,
                         "neutral_mass": neutral, "ion": ion_label(counts, deri, mode), "ion_mz": ions[mode],
                         "reference_neutral": None, "reference_ion": None, "delta_ion": None,
                         "source": "compnewv4.precursormassv4", "result": "", "ok": True})
    return rows


_COLUMNS = ["composition", "derivatization", "reduced", "mode", "neutral_mass", "ion", "ion_mz",
            "reference_neutral", "reference_ion", "delta_ion", "source", "result"]


def _fmt(value, digits=5):
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def write_table(path, rows):
    if path.lower().endswith(".md"):
        lines = ["| " + " | ".join(_COLUMNS) + " |", "|" + "---|" * len(_COLUMNS)]
        for r in rows:
            lines.append("| " + " | ".join(_fmt(r[c]) for c in _COLUMNS) + " |")
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(lines) + "\n")
    else:
        with open(path, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(_COLUMNS)
            for r in rows:
                w.writerow([_fmt(r[c]) for c in _COLUMNS])


def _parse_mode(text):
    s = str(text).strip().lower()
    if s in ("+", "pos", "positive"):
        return "+"
    if s in ("-", "−", "neg", "negative"):
        return "-"
    raise argparse.ArgumentTypeError(f"mode must be + or - (got {text!r})")


def build_parser():
    p = argparse.ArgumentParser(prog="msp_CLI_precursor_mass",
                                description="Glycan precursor masses from compnewv4.precursormassv4 (the CGA path).")
    p.add_argument("--composition", action="append", default=[], metavar="LABEL",
                   help="CGA short label, e.g. H5N4F1S2 (repeat for --table)")
    p.add_argument("--deri", default="perme", type=str.lower, choices=sorted(_DERI), help="derivatization (default perme)")
    p.add_argument("--reduced", action="store_true", help="reduced (alditol) reducing end; default free end")
    p.add_argument("--mode", type=_parse_mode, default=None, help="+ or - (default: print both ions)")
    p.add_argument("--hexa", type=int, default=None, metavar="N", help="HexA count (same as A<N> in the label)")
    p.add_argument("--so3", type=int, default=None, metavar="N", help="SO3 count (same as s<N> in the label)")
    p.add_argument("--po3h", type=int, default=None, metavar="N", help="PO3H count (same as p<N> in the label)")
    p.add_argument("--json", action="store_true", help="print JSON")
    p.add_argument("--table", metavar="FILE", help="write a table (.csv or .md)")
    p.add_argument("--selftest", action="store_true", help="check the CP-A.1 acceptance rows; exit 1 on a mismatch")
    return p


def _apply_group_flags(counts, args, parser):
    for flag, key in (("hexa", "HexA"), ("so3", "SO3"), ("po3h", "PO3H")):
        value = getattr(args, flag)
        if value is None:
            continue
        if value < 0:
            parser.error(f"--{flag} must be >= 0")
        if counts[key] and counts[key] != value:
            parser.error(f"--{flag} {value} conflicts with {key} {counts[key]} in the composition")
        counts[key] = value
    return counts


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    deri = _DERI[args.deri]
    group_flags = [f"--{f}" for f in ("hexa", "so3", "po3h") if getattr(args, f) is not None]
    if group_flags and (args.selftest or (args.table and not args.composition)):
        parser.error(f"{' '.join(group_flags)}: only with --composition (the acceptance rows are fixed)")

    if args.selftest:
        rows = [evaluate(r) for r in ACCEPTANCE_ROWS]
        for r in rows:
            print(f"{r['composition']:<10} {r['derivatization']:<6} reduced={_fmt(r['reduced']):<3} {r['ion']:<7} "
                  f"{r['ion_mz']:.5f}  ref {_fmt(r['reference_ion']) or '-':<11} {r['result']}")
        bad = [r for r in rows if not r["ok"]]
        print(f"selftest: {len(rows) - len(bad)} of {len(rows)} rows ok"
              f" (compnewv4 {compnewv4.version}, {compnewv4.last_update})")
        if args.table:
            write_table(args.table, rows)
        return 1 if bad else 0

    if args.table:
        if args.composition:
            compositions = []
            for label in args.composition:   # the same parsing, group flags and checks as one --composition (V14F-F3)
                try:
                    counts = parse_composition(label)
                except ValueError as exc:
                    parser.error(str(exc))
                counts = _apply_group_flags(counts, args, parser)
                try:
                    check_supported(counts, deri)
                except ValueError as exc:
                    parser.error(str(exc))
                compositions.append(counts)
            modes = [args.mode] if args.mode else ["+", "-"]
            rows = custom_rows(compositions, deri, args.reduced, modes)
        else:
            rows = [evaluate(r) for r in ACCEPTANCE_ROWS]
        write_table(args.table, rows)
        print(f"wrote {len(rows)} rows to {os.path.basename(args.table)}")
        return 0

    if len(args.composition) != 1:
        parser.error("give one --composition (or use --table / --selftest)")
    try:
        counts = parse_composition(args.composition[0])
    except ValueError as exc:
        parser.error(str(exc))
    counts = _apply_group_flags(counts, args, parser)
    try:
        check_supported(counts, deri)
    except ValueError as exc:
        parser.error(str(exc))
    neutral, ions = masses(counts, deri, args.reduced)
    modes = [args.mode] if args.mode else ["+", "-"]
    labels = {m: ion_label(counts, deri, m) for m in modes}
    record = {"composition": label_of(counts), "counts": counts, "derivatization": deri,
              "reduced": args.reduced, "neutral_mass": neutral,
              "ions": {labels[m]: ions[m] for m in modes},
              "function": "compnewv4.precursormassv4",
              "compnewv4": {"version": compnewv4.version, "last_update": compnewv4.last_update}}
    if args.json:
        print(json.dumps(record, indent=2))
        return 0
    print(f"composition     {record['composition']}  ("
          + ", ".join(f"{k} {v}" for k, v in counts.items()) + ")")
    print(f"derivatization  {deri}, reducing end {'reduced' if args.reduced else 'free'}")
    print(f"neutral M       {_fmt(neutral) or 'n/a (PerMe SO3/PO3H values carry fixed charge corrections)'}")
    for m in modes:
        print(f"{labels[m]:<15} {ions[m]:.5f}")
    print(f"source          compnewv4.precursormassv4 (compnewv4 {compnewv4.version}, {compnewv4.last_update})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
