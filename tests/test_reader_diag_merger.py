"""READER-DIAG: MAS merge stops with a clear error on an invalid converted CSV (2026-10-13).

Before: validate_csv_structure returned False, directassign_files carried on and raised
UnboundLocalError on df2. Trigger used here is the one Codex F2 found: an all-unknown-mass
file whose protonatedmass column is integer text 0, so pandas infers int64.
"""
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import mspvalidator_merger as M  # noqa: E402

HEADER = ['entry_no', 'MS1scan_no', 'MS1_isolationmass', 'MS1_monoisolationmass', 'chargeState',
          'protonatedmass', 'MS2scan_no', 'peaklist', 'peakintensity']


def _csv(path, mass_cell):
    path.write_text('\t'.join(HEADER) + '\n'
                    + f'1\t10\t500.25\t500.25\t0\t{mass_cell}\t11\t(150.0, 200.5)\t(1.0, 2.0)\n')
    return str(path)


def test_default_return_is_still_bool(tmp_path):
    assert M.validate_csv_structure(_csv(tmp_path / 'ok.csv', '999.5')) is True
    assert M.validate_csv_structure(_csv(tmp_path / 'bad.csv', '0')) is False


def test_return_reason(tmp_path):
    assert M.validate_csv_structure(_csv(tmp_path / 'ok.csv', '999.5'), return_reason=True) == (True, None)
    ok, reason = M.validate_csv_structure(_csv(tmp_path / 'bad.csv', '0'), return_reason=True)
    assert ok is False and "protonatedmass" in reason and "int64" in reason
    ok, reason = M.validate_csv_structure(str(tmp_path / 'missing.csv'), return_reason=True)
    assert ok is False and reason.startswith("could not be read")


def test_directassign_files_raises_clear_error(tmp_path):
    pytest.importorskip("openpyxl")
    xlsx = tmp_path / 'anno.xlsx'
    with pd.ExcelWriter(xlsx) as writer:
        pd.DataFrame({'Structure': ['H5N4'], 'MS2scan_no': [11], 'charge': [1], 'mass shift': [0.0],
                      'adduct': ['H'], 'note': ['']}).to_excel(writer, sheet_name='MSlist', index=False)
        pd.DataFrame({'mass': [204.087]}).to_excel(writer, sheet_name='ionlist', index=False)
    with pytest.raises(ValueError, match=r"bad\.csv' failed structure validation: column 'protonatedmass'"):
        M.directassign_files(str(xlsx), _csv(tmp_path / 'bad.csv', '0'), 'PerMe', debug=False)
