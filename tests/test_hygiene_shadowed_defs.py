"""Hygiene H3: shadowed duplicate definitions removed (second definition kept).

- msp_insilicomarker_withGPT.extract_ionmasslist: only the
  (ionmass_sheet, mass_col="mass") version remains; the dead hard-coded
  `if __name__ == "__main__":` dev block that called it was removed with it.
- pretrain_normalizer._canon_structure_label: only the effective second
  definition remains.

See handoff/hygiene_LF1_kickoff_20260926.md (H3).
"""
import ast
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
SRC = os.path.join(os.path.dirname(__file__), '..', 'src')

import pandas as pd
import pytest

import msp_insilicomarker_withGPT as marker
import pretrain_normalizer as normalizer


def _top_level_def_count(filename, name):
    with open(os.path.join(SRC, filename), encoding="utf-8") as fh:
        tree = ast.parse(fh.read())
    return sum(isinstance(n, ast.FunctionDef) and n.name == name for n in tree.body)


@pytest.mark.parametrize("filename,name", [
    ("msp_insilicomarker_withGPT.py", "extract_ionmasslist"),
    ("pretrain_normalizer.py", "_canon_structure_label"),
])
def test_single_definition(filename, name):
    assert _top_level_def_count(filename, name) == 1


def test_extract_ionmasslist_mass_col():
    sheet = pd.DataFrame({"mz": [204.0867, "bad", 366.1395], "mass": [1.0, 2.0, 3.0]})
    assert marker.extract_ionmasslist(sheet, mass_col="mz") == [204.0867, 366.1395]
    assert marker.extract_ionmasslist(sheet) == [1.0, 2.0, 3.0]


def test_extract_ionmasslist_fallback_column():
    sheet = pd.DataFrame({"m/z": [100.5, 200.5]})
    assert marker.extract_ionmasslist(sheet) == [100.5, 200.5]


def test_canon_structure_label_effective_behaviour():
    assert normalizer._canon_structure_label("H5N4F1") == "H5N4F1"
    assert normalizer._canon_structure_label("Hex5 HexNAc4 NeuAc1") == "H5N4S1"
    assert normalizer._canon_structure_label("Non-glycan") == "Non-glycan"
