"""Hygiene H1 (LF1): Score B composition tokenizer is case-sensitive.

Locks the 2026-10-01 contract of `parse_compact_composition_label`:
- only uppercase F/H/N/S/G and KDN are counted; sulfate `s`, phosphate `p`,
  HexA `A` and other leftovers no longer alias to NeuAc or vanish silently
- labels with leftovers still parse (no raise) and are recorded in
  `_UNPARSED_TOKEN_LABELS`; a label with no recognized token still raises
- `enrich_cga_tsv_with_score_b` prints exactly one summary warning per run

The enrich tests stub the workbook loader and scorer so no Score B workbook is
needed; they exercise the TSV loop, the per-run reset and the warning only.
See handoff/hygiene_LF1_kickoff_20260926.md (H1).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pandas as pd
import pytest

import msp_CGA_structscore as sb

WARN_TAG = "[score_b][WARN]"


@pytest.fixture(autouse=True)
def _clear_unparsed():
    sb._UNPARSED_TOKEN_LABELS.clear()
    yield
    sb._UNPARSED_TOKEN_LABELS.clear()


# ---------------------------------------------------------------- tokenizer

def test_valid_uppercase_label_unchanged():
    assert sb.parse_compact_composition_label("H3N2S2") == (3, 2, 2, 0, 0, 0)
    assert sb.parse_compact_composition_label("F1H3N3KDN1") == (3, 3, 0, 0, 1, 1)
    assert sb.parse_compact_composition_label("H5N4S1") == (5, 4, 1, 0, 0, 0)
    assert sb._UNPARSED_TOKEN_LABELS == []


def test_kdn_leading_still_parses():
    assert sb.parse_compact_composition_label("KDN1H3N2") == (3, 2, 0, 0, 1, 0)
    assert sb._UNPARSED_TOKEN_LABELS == []


def test_sulfate_not_counted_as_neuac():
    # pre-fix: (3, 2, 1, 0, 0, 0)
    assert sb.parse_compact_composition_label("H3N2s1") == (3, 2, 0, 0, 0, 0)
    assert sb._UNPARSED_TOKEN_LABELS == ["H3N2s1"]
    assert sb.unparsed_composition_tokens("H3N2s1") == ["s1"]


def test_sulfate_alongside_neuac_keeps_neuac_count():
    # pre-fix: (3, 2, 3, 0, 0, 0)
    assert sb.parse_compact_composition_label("H3N2S2s1") == (3, 2, 2, 0, 0, 0)
    assert sb._UNPARSED_TOKEN_LABELS == ["H3N2S2s1"]


def test_hexa_and_phosphate_reported_as_residuals():
    assert sb.parse_compact_composition_label("H3N2A1p1") == (3, 2, 0, 0, 0, 0)
    assert sb.unparsed_composition_tokens("H3N2A1p1") == ["A1", "p1"]
    assert sb._UNPARSED_TOKEN_LABELS == ["H3N2A1p1"]


def test_live_cga_producer_labels():
    # v14 builds the CGA `composition` column with marker.canonical_label_from_row(style="short"),
    # which emits A/s/p for HexA/SO3/PO3H (Codex H1 review, packet correction).
    import msp_insilicomarker_withGPT as marker
    lab = marker.canonical_label_from_row({"Hex": 3, "HexNAc": 2, "NeuAc": 1, "SO3": 1}, style="short")
    assert lab == "H3N2S1s1"
    assert sb.parse_compact_composition_label(lab) == (3, 2, 1, 0, 0, 0)
    assert sb.unparsed_composition_tokens(lab) == ["s1"]
    plain = marker.canonical_label_from_row({"Hex": 5, "HexNAc": 4, "KDN": 1, "Fuc": 1}, style="short")
    assert sb.parse_compact_composition_label(plain) == (5, 4, 0, 0, 1, 1)
    assert sb._UNPARSED_TOKEN_LABELS == [lab]


def test_lowercase_only_label_raises():
    with pytest.raises(ValueError, match="no recognized tokens"):
        sb.parse_compact_composition_label("h3n2s1")
    assert sb._UNPARSED_TOKEN_LABELS == []


def test_separators_are_not_residuals():
    assert sb.parse_compact_composition_label("H5 N4") == (5, 4, 0, 0, 0, 0)
    assert sb.parse_compact_composition_label("H5-N4,F1") == (5, 4, 0, 0, 0, 1)
    assert sb._UNPARSED_TOKEN_LABELS == []


def test_unicode_residuals_are_recorded():
    # Codex H1 review F1: non-ASCII leftovers must not escape the warning
    assert sb.parse_compact_composition_label("H5N4Ω") == (5, 4, 0, 0, 0, 0)
    assert sb.unparsed_composition_tokens("H5N4Ω") == ["Ω"]
    assert sb.unparsed_composition_tokens("H5N4Ｓ¹") == ["Ｓ¹"]
    assert sb._UNPARSED_TOKEN_LABELS == ["H5N4Ω"]


# ---------------------------------------------------------------- run-level warning

def _stub_scoring(monkeypatch, seen):
    """Replace workbook loading and scoring with deterministic stubs."""
    monkeypatch.setattr(sb, "load_scoreb_workbook", lambda path: object())
    monkeypatch.setattr(sb, "extract_observed_ion_hits_from_peaks", lambda *a, **k: [])
    monkeypatch.setattr(sb, "evaluate_global_motif_evidence", lambda cfg, spectrum: None)
    monkeypatch.setattr(sb, "_build_score_b_motif_summary", lambda result: "")

    def _fake_eval(cfg, global_evidence, *, candidate_composition, candidate_flags=None):
        seen.append(candidate_composition)
        total = float(sum(candidate_composition))
        return sb.ScoreBResult(candidate_composition, total, 0.0, 0.0, total, {})

    monkeypatch.setattr(sb, "evaluate_candidate_score_b", _fake_eval)


def _write_tsv(path, compositions):
    rows = [
        {"MS2scan_no": i // 2 + 1, "composition": c,
         "peaklist": "[100.0, 200.0]", "peakintensity": "[10.0, 20.0]"}
        for i, c in enumerate(compositions)
    ]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_enrich_warns_once_per_run(tmp_path, monkeypatch, capsys):
    seen = []
    _stub_scoring(monkeypatch, seen)
    tsv = tmp_path / "cga.tsv"
    _write_tsv(tsv, ["H3N2s1", "H5N4S1", "H3N2s1", "H3N2A1p1"])

    sb.enrich_cga_tsv_with_score_b(tsv, tmp_path / "unused.xlsx")
    out = capsys.readouterr().out
    warn_lines = [ln for ln in out.splitlines() if WARN_TAG in ln]
    assert len(warn_lines) == 1
    assert "3 row(s) / 2 label(s)" in warn_lines[0]
    assert "H3N2s1 [s1]" in warn_lines[0] and "H3N2A1p1 [A1 p1]" in warn_lines[0]

    # sulfated rows were scored on (H, N) only, not as NeuAc
    assert (3, 2, 0, 0, 0, 0) in seen and (3, 2, 1, 0, 0, 0) not in seen
    df = pd.read_csv(tsv, sep="\t")
    assert df["score_b"].notna().all()


def test_warning_stays_on_one_line(capsys):
    # Codex H1 review F2: a label with an embedded newline must not split the warning
    sb.parse_compact_composition_label("H3N2\ns1")
    sb._warn_unparsed_composition_tokens()
    out = capsys.readouterr().out
    assert out.count("\n") == 1 and WARN_TAG in out
    assert "H3N2\\ns1 [s1]" in out
    assert sb._UNPARSED_TOKEN_LABELS == ["H3N2\ns1"]  # recorded label is not altered


def test_enrich_resets_between_runs(tmp_path, monkeypatch, capsys):
    seen = []
    _stub_scoring(monkeypatch, seen)
    dirty = tmp_path / "dirty.tsv"
    clean = tmp_path / "clean.tsv"
    _write_tsv(dirty, ["H3N2s1", "H5N4"])
    _write_tsv(clean, ["H5N4S1", "F1H3N3KDN1"])

    sb.enrich_cga_tsv_with_score_b(dirty, tmp_path / "unused.xlsx")
    first = capsys.readouterr().out
    assert sum(WARN_TAG in ln for ln in first.splitlines()) == 1

    # a label recorded outside a run (e.g. curated-library scoring) must not leak in
    sb.parse_compact_composition_label("H3N2p1")
    sb.enrich_cga_tsv_with_score_b(clean, tmp_path / "unused.xlsx")
    second = capsys.readouterr().out
    assert not any(WARN_TAG in ln for ln in second.splitlines())
