"""Codex S1 R3-F3 — collected AFTER the S1 modules (by filename; `test_smoke_*` sort later) and asserts
the S1 modules restored the environment. Assumes default alphabetical collection; explicit selection or
parallel runners can bypass it, in which case the fixture's own MonkeyPatch cleanup is the guarantee.

`tests/test_s1_v14_headless.py` scopes its GLYCOMSP_REFERENCE_DB override to a MonkeyPatch context; this
module checks the exact prior state (including 'unset') survived, so a leak would fail the suite rather than
silently redirect later in-process CLI/GUI database selection.
"""
import os

import conftest


def test_reference_db_env_restored_after_s1_modules():
    for k, prior in conftest.GLYCOMSP_ENV_AT_SESSION_START.items():
        assert os.environ.get(k) == prior, f"{k} leaked: {os.environ.get(k)!r} != {prior!r}"


def test_headless_fixture_has_no_broad_skip():
    """R3-F4: the headless fixture must not convert arbitrary exec/compile errors into skips."""
    src = open(os.path.join(os.path.dirname(__file__), "test_s1_v14_headless.py"), encoding="utf-8").read()
    fixture = src[src.index("def v14ns("):src.index("def _sample(")]
    assert "except Exception" not in fixture and "pytest.skip" not in fixture
    assert "pytest.skip" in src[src.index("def _find_anchors("):src.index("@pytest.fixture")]   # only anchors skip
