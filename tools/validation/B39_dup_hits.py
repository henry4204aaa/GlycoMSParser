"""Count Score A rows whose hit list repeats an ion (several peaks inside one ion's window).

usage: python B39_dup_hits.py <alias>=<CGA TSV> [<alias>=<CGA TSV> ...]
Only the alias is printed (naming rule, 2026-11-01); an argument without "<alias>=" is printed as file1, file2, ..."""
import sys
import pandas as pd

for _i, _arg in enumerate(sys.argv[1:], 1):
    alias, sep, p = _arg.partition("=")
    if not sep:
        alias, p = f"file{_i}", _arg
    t = pd.read_csv(p, sep="\t", usecols=lambda c: c in ("MS2scan_no", "composition", "ion score", "ion hit count", "ion hits m/z"))
    t = t[t["composition"].notna()]
    hits = t["ion hits m/z"].fillna("").astype(str)
    n_all = hits.map(lambda s: len([x for x in s.split(";") if x]))
    n_uni = hits.map(lambda s: len(set(x for x in s.split(";") if x)))
    dup = n_all > n_uni
    print(alias, f"rows={len(t)} rows_with_repeated_ion={int(dup.sum())} "
          f"({100 * dup.mean():.1f}%) extra_hits={int((n_all - n_uni).sum())} "
          f"score>=0.07 counted={int((t['ion score'] >= 0.07).sum())} "
          f"max_hits_per_row={int(n_all.max()) if len(t) else 0} max_unique={int(n_uni.max()) if len(t) else 0}")
    if "ion score" in t:
        # Score A if each ion counted once (denominator = list size, recovered as hits/score)
        denom = (t["ion hit count"] / t["ion score"]).where(t["ion score"] > 0)
        uniq_score = (n_uni / denom).fillna(0)
        print("   score>=0.07 if each ion counted once:", int((uniq_score >= 0.07).sum()))
